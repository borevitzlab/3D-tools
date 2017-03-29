"""
Tools to read and write georeferenced pointclouds such as output from Pix4D,
based on the plyfile module.
"""

import collections
import itertools
import json
import tempfile
import warnings

import numpy as np
import plyfile


UTM_COORD = collections.namedtuple(
    'UTMCoord', ['easting', 'northing', 'zone', 'northern'])


def get_tmpfile():
    """Create a temporary file, for easy use of np.memmap"""
    # TODO:  allow user to specify where this is, eg. $JOBFS
    return tempfile.SpooledTemporaryFile(max_size=2**20)


class GeoPly(plyfile.PlyData):
    """A pointcloud, with the UTM georeference for the origin coordinate.

    Spatial units are meters east (X) and north (Y) of the origin.
    Pointcloud data typically has a substantial XY offset from the UTM
    origin, which is retained due to the precision limits of 32-bit floats.
    """
    # Marker for special UTM coord comments
    _COORD_MARKER = 'UTM_COORD='

    #pylint:disable=too-many-arguments
    def __init__(self, elements=None, text=False, byte_order='=',
                 comments=None, obj_info=None, *, utm_coord, memmap=None):
        """Create a GeoPly instance.  utm_coord is a required keyword arg."""
        # Validate utm_coord a little
        self.utm_coord = utm_coord
        if not isinstance(self.utm_coord, UTM_COORD):
            raise ValueError('Must include the UTM coords of the local origin')
        # Handle the more flexible argument types allowed here
        if isinstance(elements, np.ndarray):
            elements = [plyfile.PlyElement.describe(elements, 'vertex')]
        # Call parent __init__, avoiding mutable default arguments
        elements = elements or []
        comments = comments or []
        obj_info = obj_info or []
        super().__init__(elements, text, byte_order, comments, obj_info)
        # Memmap if requested, or autodetecting and many vertices
        if memmap is None:
            memmap = self['vertex'].data.size >= 10**7
        if memmap and not isinstance(self['vertex'].data, np.memmap):
            mmap = np.memmap(get_tmpfile(), dtype=self['vertex'].data.dtype,
                             shape=self['vertex'].data.shape)
            mmap[:] = self['vertex'].data[:]
            self['vertex'].data = mmap


    @staticmethod
    def read(stream):
        """Reads vertices from ``stream``, with UTM offset and data cleaning.
        stream may be a filename, or a file-like object.

        The UTM coordinate (self.utm_coord) is read from
        - comments in the file header, if the pointcloud was created
          by this class.
        - the corresponding '_ply_offset.xyz' file, if the pointcloud was
          created by Pix4D.  In this case, the Z-offset is added to vertex
          coordinates.

        Data cleaning consists of:
        - discarding non-"vertex" elements (if present)
        - removing the marker comment if written by Meshlab, and if a uniform
          alpha channel was added removing that too
        """
        data = plyfile.PlyData.read(stream)
        verts = data['vertex']

        # Remove meshlab cruft
        if 'VCGLIB generated' in data.comments:
            names = verts.data.dtype.names  # field names of each vertex
            if 'alpha' in names and len(np.unique(verts['alpha'])) == 1:
                # properties of the PlyElement instance are manually updated
                verts.properties = [p for p in verts.properties
                                    if p.name != 'alpha']
                # removal of a vertex field is via fancy indexing
                verts.data = verts.data[[n for n in names if n != 'alpha']]
            data.comments.remove('VCGLIB generated')

        # Add UTM coordinates if known or discoverable
        utm_coord = None
        coords = []
        for c in data.comments:
            if c.startswith(GeoPly._COORD_MARKER):
                data.comments.remove(c)
                serialised = c.lstrip(GeoPly._COORD_MARKER)
                coords.append(UTM_COORD(**json.loads(serialised)))
        if coords:
            utm_coord = coords[0]
            if len(coords) > 1:
                msg = RuntimeWarning('Found multiple coordinates in comments:'
                                     '{}, using first...'.format(coords))
                warnings.warn(msg)
        else:
            # Try to find and apply the Pix4D offset, which may raise...
            z_offset, utm_coord = GeoPly._offset_from_pix4d(stream)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', FutureWarning)
                # Numpy wories about writing to multiple columns here
                verts['z'] += z_offset

        # Return as GeoPly instance with only vertex elements
        return GeoPly([verts], data.text, data.byte_order,
                      data.comments, data.obj_info, utm_coord=utm_coord)


    def write(self, stream):
        """Write to a file, serialising utm_coord as a special comment."""
        assert not any(c.startswith(self._COORD_MARKER) for c in self.comments)
        # Serialise as JSON dict following the marker string
        serialised = self._COORD_MARKER + json.dumps(self.utm_coord._asdict())
        # Insert, write, pop - keeps comments in correct state
        self.comments.insert(0, serialised)
        super().write(stream)
        self.comments.pop(0)


    @staticmethod
    def _offset_from_pix4d(ply_filename, utm_zone=55,
                           in_northern_hemisphere=False):
        """Return a (offset_z, UTM coord) for a .ply file from the
        corresponding Pix4D offset file and provided zone.
        Raises FileNotFoundError if the offset file is invalid or nonexistent.
        """
        assert ply_filename.endswith('.ply')
        offset_file = ply_filename[:-4] + '_ply_offset.xyz'
        try:
            with open(offset_file) as f:
                line = f.readline().strip()
            x, y, z = (float(n) for n in line.split(' '))
        except Exception:
            raise FileNotFoundError
        return z, UTM_COORD(x, y, utm_zone, in_northern_hemisphere)


    @property
    def vertices(self):
        """Return a read-only view of the vertex data."""
        vertices = np.ndarray.view(self['vertex'].data)
        vertices.flags.writeable = False
        return vertices


    @staticmethod
    def from_iterable(iterable, utm_coord, **kwargs):
        """Create a GeoPly from an iterable of vertices and a UTM offset.
        The iterable must contain numpy scalars with a consistent dtype.
        """
        it = iter(iterable)
        first = next(it)
        assert isinstance(first, np.void)
        array = np.fromiter(itertools.chain([first], it), first.dtype)
        return GeoPly(array, utm_coord=utm_coord, **kwargs)


    @staticmethod
    def from_array(array, utm_coord, **kwargs):
        """Create a GeoPly from a Numpy array of vertices and a UTM offset."""
        assert isinstance(array, np.ndarray)
        assert all(dim in array.dtype.names for dim in 'xyz')
        return GeoPly(array, utm_coord=utm_coord, **kwargs)


    @classmethod
    def from_geoplys(cls, *geoplys):
        """Create a new geoply by combining two or more GeoPly instances.

        All inputs must have compatible georeferences and datatypes.
        The output GeoPly uses the base georeference and comcatenates all
        input vertices, applying relative offsets.  If any of the inputs
        stored vertices in a np.memmap, so will the output.
        """
        assert len(geoplys) >= 2
        assert all(isinstance(p, cls) for p in geoplys)
        assert all(p.utm_coord is not None for p in geoplys)
        assert len(set(p.utm_coord.zone for p in geoplys)) == 1
        assert len(set(p.utm_coord.northern for p in geoplys)) == 1
        assert len(set(p['vertex'].data.dtype for p in geoplys)) == 1

        # flatten and deduplicate comments
        comments = [c for pf in geoplys for c in pf.comments]
        comments = sorted(set(comments), key=comments.index)

        # paste arrays into single memmap, handling UTM offsets
        using_memmap = any(isinstance(p['vertex'].data, np.memmap)
                           for p in geoplys)
        to_arr = (np.memmap if using_memmap else np.array)(
            get_tmpfile(), dtype=geoplys[0]['vertex'].data.dtype,
            shape=(sum([p['vertex'].data.size for p in geoplys]),))
        base, *other_files = geoplys
        start = base['vertex'].data.size
        to_arr[:start] = base['vertex'].data
        for pf in other_files:
            arr = np.copy(pf['vertex'].data)
            arr['x'] += (pf.utm_coord.easting - base.utm_coord.easting)
            arr['y'] += (pf.utm_coord.northing - base.utm_coord.northing)
            to_arr[start:start+arr.size] = arr
            start += arr.size

        # Load data back into the complete structure and return
        return cls(to_arr, comments=comments, utm_coord=base.utm_coord)

