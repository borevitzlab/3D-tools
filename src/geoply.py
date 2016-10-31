"""
Tools to read and write georeferenced pointclouds such as output from Pix4D,
based on the plyfile module.
"""

import collections
import contextlib
import tempfile

import numpy as np
import plyfile


UTM_COORD = collections.namedtuple(
    'UTMCoord', ['easting', 'northing', 'zone', 'northern'])


def get_tmpfile():
    """Create a temporary file, for easy use of np.memmap"""
    # TODO:  allow user to specify where this is, eg. $JOBFS
    return tempfile.SpooledTemporaryFile(max_size=2**20)


class GeoPly(plyfile.PlyData):
    """A pointcloud, with the UTM georeference for it's origin coordinates.

    Spatial units are meters east (X) and north (Y) of the origin.
    Pointcloud data typically hsa a substantial XY offset from the UTM
    origin point, which is retained due to the precision limits of 32-bit
    floats.  Z-offsets are
    """

    def __init__(self, elements=None, text=False, byte_order='=',
                 comments=None, obj_info=None, *, utm_coord=None, memmap=True):
        #
        self.utm_coord = utm_coord
        if not isinstance(self.utm_coord, UTM_COORD):
            raise ValueError('Must include the UTM coords of the local origin')
        # Call parent __init__, avoiding mutable default arguments
        elements = elements or []
        comments = comments or []
        obj_info = obj_info or []
        super().__init__(elements, text, byte_order, comments, obj_info)
        # Use a memmap for vertex data by default
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
          created by Pix4D.

        Data cleaning consists of:
        - discarding non-"vertex" elements (if present)
        -
        Additionally:
            - applies Pix4D altitude offset
              (unlike XY, Z is small enougth that precision is not an issue)
            - removes Meshlab cruft if present (comment, empty alpha channel)
        """
        data = plyfile.PlyData.read(stream)
        verts = data['vertex']

        # Remove meshlab cruft
        if 'VCGLIB generated' in data.comments:
            names = verts.data.dtype.names  # field names of each vertex
            if 'alpha' in names and 1 == len(np.unique(verts['alpha'])):
                # properties of the PlyElement instance are manually updated
                verts.properties = [p for p in verts.properties
                                    if p.name != 'alpha']
                # removal of a vertex field is via fancy indexing
                verts.data = verts.data[[n for n in names if n != 'alpha']]
            data.comments.remove('VCGLIB generated')

        # Add UTM coordinates if known or discoverable
        utm_coord = UTM_COORD(0, 0, 55, False)  # FIXME:  better defaults
        if isinstance(stream, str):  # FIXME:  or Path???
            with contextlib.suppress(Exception):
                # First, let's try for a Pix4D-style offset file
                with open(stream[:-4] + '_ply_offset.xyz') as f:
                    x, y, z = tuple(
                        float(n) for n in f.readline().strip().split(' '))
                utm_coord = utm_coord._replace(easting=x, northing=y)
                verts['z'] = verts['z'] + z
        for c in data.comments:
            if c.startswith('utm_'):
                types = dict(zip(UTM_COORD._fields,
                                 [float, float, int, lambda b: b == 'True']))
                attr, val = c.replace('utm_', '').split()
                setattr(utm_coord, attr, types[attr](val))
                data.comments.remove(c)
        if not (utm_coord.easting and utm_coord.northing):
            utm_coord = None

        # Return as GeoPly instance with only vertex elements
        return GeoPly([verts], data.text, data.byte_order,
                      data.comments, data.obj_info, utm_coord=utm_coord)


    def write(self, stream):
        """Write to a file, serialising utm_coord as a special comment."""
        if self.utm_coord is not None:
            if any(c.startswith('utm_') for c in self.comments):
                # TODO:  log warning
                self.comments = [c for c in self.comments if c[:4] != 'utm_']
            for attr, val in zip(UTM_COORD._fields, self.utm_coord):
                self.comments.append('utm_{0} {1}'.format(attr, val))
        super().write(stream)


    @property
    def vertices(self):
        """Return a read-only view of the vertex data."""
        vertices = np.ndarray.view(self['vertex'].data)
        vertices.flags.writeable = False
        return vertices


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
        return GeoPly([plyfile.PlyElement.describe(to_arr, 'vertex')],
                      comments=comments, utm_coord=base.utm_coord)

