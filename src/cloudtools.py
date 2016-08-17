#!/usr/bin/env python3
"""
`cloudtools` filters and optionally combines .ply pointclouds.

Units are assumed to be meters in the UTM coordinate system, with some
offset to avoid float32 imprecision.  This matches the output of
software such as Pix4D.
"""

import argparse
import collections
import contextlib
import functools
import glob
import os
import pkg_resources
import tempfile
import time
import warnings

import numpy as np
import plyfile


# I actually *want* the new numpy behaviour...
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

UTM_COORD = collections.namedtuple(
    'UTMCoord', ['easting', 'northing', 'zone', 'northern'])


def get_tmpfile():
    return tempfile.SpooledTemporaryFile(max_size=2**20)


class reify:
    """Use as a class method decorator.  It operates almost exactly like the
    Python ``@property`` decorator, but it puts the result of the method it
    decorates into the instance dict after the first call, effectively
    replacing the function it decorates with an instance variable."""
    # From Pyramid.  This allows us to write code as if the terrain array
    # was an attribute calculated in advance, while only paying for this
    # for instances where it is accessed.
    def __init__(self, wrapped):
        self.wrapped = wrapped
        functools.update_wrapper(self, wrapped)

    def __get__(self, inst, objtype=None):
        if inst is None:
            return self
        val = self.wrapped(inst)
        setattr(inst, self.wrapped.__name__, val)
        return val


class GeoPly(plyfile.PlyData):
    """Represents a pointcloud in the UTM coordinate system including zone,
    hemisphere, and offset (due to f32 imprecision)."""

    def __init__(self, elements=None, text=False, byte_order='=',
                 comments=None, obj_info=None, utm_coord=None):
        elements = elements or []
        comments = comments or []
        obj_info = obj_info or []
        super().__init__(elements, text, byte_order, comments, obj_info)
        self.utm_coord = utm_coord
        # replace the vertex data with a memmapped version of itself
        if not isinstance(self['vertex'].data, np.memmap):
            mmap = np.memmap(get_tmpfile(), dtype=self['vertex'].data.dtype,
                             shape=self['vertex'].data.shape)
            mmap[:] = self['vertex'].data[:]
            self['vertex'].data = mmap

    @staticmethod
    def read(filename):
        """Reads a vertex set from the given file.
        Additionally:
            - drops non-"vertex" elements if present
            - adds a self.utm_coord attribute
            - applies Pix4D altitude offset
              (unlike XY, Z is small enougth that precision is not an issue)
            - removes Meshlab cruft if present (comment, empty alpha channel)
        """
        # Read file as GeoPly instance
        with open(filename, 'rb') as stream:
            data = GeoPly._parse_header(stream)
            for elt in data:
                elt._read(stream, data.text, data.byte_order)
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
        utm_coord = UTM_COORD(0, 0, args.utmzone, args.northern)
        with contextlib.suppress(Exception):
            # First, let's try for a Pix4D-style offset file
            with open(filename[:-4] + '_ply_offset.xyz') as f:
                x, y, z = tuple(
                    float(n) for n in f.readline().strip().split(' '))
            utm_coord = utm_coord._replace(easting=x, northing=y)
            verts['z'] += z
        for c in data.comments:
            fword, *rest = c.strip().split()
            if fword == 'utm_offset':
                utm_coord = utm_coord._replace(
                    easting=float(rest[0]), northing=float(rest[1]))
            elif fword == 'utm_zone':
                utm_coord = utm_coord._replace(zone=int(rest[0]))
            elif fword == 'utm_northern':
                utm_coord = utm_coord._replace(northern=(rest[0] == 'True'))
        if not (utm_coord.easting and utm_coord.northing):
            utm_coord = None

        # Return as GeoPly instance with only vertex elements
        return GeoPly([verts], data.text, data.byte_order,
                      data.comments, data.obj_info, utm_coord)

    def write(self, stream):
        """Write to a file, serialising utm_coord as a special comment."""
        if self.utm_coord is not None:
            self.comments = [c for c in self.comments if c[:4] != 'utm_']
            self.comments.extend(
                'utm_offset {} {}\nutm_zone {}\nutm_northern {}'.format(
                    *self.utm_coord).split('\n'))
        super().write(stream)

    @reify
    def terrain(self):
        """An array describing the terrain and trees in the pointcloud."""
        verts = self['vertex'].data
        min_x, min_y = (np.floor(verts[c].min()) for c in 'xy')
        other_names = tuple(n for n in verts.dtype.names if n not in 'xyz')

        shape = tuple(np.ceil((verts[c].max() - verts[c].min() + 1) /
                              args.cellsize) for c in 'xy')
        dtype = [('points', 'i4'), ('ground', 'f4'), ('canopy', 'f4')] + [
            (name, 'f8' if 'f' in dt else ('u4' if 'u' in dt else 'i4'))
            for name, dt in verts.dtype.descr if name not in 'xyz']
        terrain = np.zeros(shape, dtype=dtype)
        terrain['ground'] = np.inf
        terrain['canopy'] = -np.inf

        # TODO:  move XY->idx conversion to array operations for speed

        for vert in verts:
            tx = terrain[(np.floor((vert['x'] - min_x) / args.cellsize),
                          np.floor((vert['y'] - min_y) / args.cellsize))]
            tx['points'] += 1
            if tx['ground'] > vert['z']:
                tx['ground'] = vert['z']
            if tx['canopy'] < vert['z']:
                tx['canopy'] = vert['z']
            for n in other_names:
                tx[n] += vert[n]

        return terrain


def concat(*geoplys):
    """Return a GeoPly instance with the vertex elements from all inputs.
    """
    if not geoplys:
        return GeoPly()
    # do some sanity-checking
    assert all(isinstance(p, GeoPly) for p in geoplys)
    assert all(p.utm_coord is not None for p in geoplys)
    assert len(set(p.utm_coord.zone for p in geoplys)) == 1
    assert len(set(p.utm_coord.northern for p in geoplys)) == 1
    assert len(set(p['vertex'].data.dtype for p in geoplys)) == 1

    # flatten and deduplicate comments
    comments = [c for pf in geoplys for c in pf.comments]
    comments = sorted(set(comments), key=lambda k: comments.index(k))

    # paste arrays into single memmap, handling UTM offsets
    ply_1, *ply_rest = geoplys
    to_arr = np.memmap(
        get_tmpfile(), dtype=geoplys[0]['vertex'].data.dtype,
        shape=(sum([p['vertex'].data.size for p in geoplys]),))
    start = ply_1['vertex'].data.size
    to_arr[:start] = ply_1['vertex'].data
    for pf in ply_rest:
        arr = pf['vertex'].data
        arr['x'] += (pf.utm_coord.easting - ply_1.utm_coord.easting)
        arr['y'] += (pf.utm_coord.northing - ply_1.utm_coord.northing)
        to_arr[start:start+arr.size] = arr
        start += arr.size
    verts = plyfile.PlyElement.describe(to_arr, 'vertex')

    # Load data back into the complete structure and return
    return GeoPly([verts], comments=comments, utm_coord=ply_1.utm_coord)


def add_clouds(args):
    """A tester function to add things togther."""
    print('Concatenating input files...')
    output_file = os.path.join(args.output_dir, 'out.ply')
    t1 = time.time()
    out = concat(*[GeoPly.read(f) for f in args.input_glob])
    out.write(output_file)
    t2 = time.time()
    print('Saved to {} in {} seconds'.format(output_file, int(t2-t1)))
    print(out.terrain['points'])
    print('Took {} seconds to generate terrain'.format(int(time.time()-t2)))


def get_args():
    """Handle command-line arguments, including default values.

    The goal is to make argparse do as much of the work as possible, eg with
    type validation functions, to 'fail fast' and simplify other code.
    """
    def glob_arg_type(val):
        if not val.endswith('.ply'):
            raise argparse.ArgumentTypeError('glob pattern must end ".ply"')
        files = tuple(glob.glob(val))
        if not files:
            raise argparse.ArgumentTypeError('glob pattern matchs no files')
        return files

    def dir_arg_type(val):
        if not os.path.exists(val):
            raise argparse.ArgumentTypeError('path does not exist')
        if not os.path.isdir(val):
            raise argparse.ArgumentTypeError('path is not a directory')
        return val

    df = ' default=%(default)s'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-V', '--version', action='version',
        version=pkg_resources.get_distribution('forestutils').version)

    io_g = parser.add_argument_group('Input/Output Args')
    io_g.add_argument(
        'input_glob', type=glob_arg_type,
        help='glob pattern for the file(s) to process, ending ".ply"')
    io_g.add_argument(
        'output_dir', type=dir_arg_type, help='directory for output files')

    filt_g = parser.add_argument_group('Filtering Args')
    filt_g.add_argument(
        '--cellsize', default=0.1, type=float,
        help=('edge length of the analysis pixels in meters.  Around 10 times'
              ' point spacing seems to give good results.' + df))
    filt_g.add_argument(
        '--grassdepth', metavar='DEPTH', default=0.2, type=float,
        help=('filter out points less than this height in meters above lowest'
              ' point per cell.  n <= 0 disables this filter.' + df))

    geo_g = parser.add_argument_group('Geolocation Args')
    geo_g.add_argument(
        '--utmzone', default=55, type=int, choices=range(1, 61),
        metavar='[1..60]', help='the UTM coordinate zone of the data,' + df)
    ns_g = geo_g.add_mutually_exclusive_group()
    ns_g.add_argument('--northern', action='store_true',
                      help='for data in the northern hemisphere,' + df)
    ns_g.add_argument('--southern', action='store_true', default=True,
                      help='for data in the southern hemisphere,' + df)

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    add_clouds(args)
