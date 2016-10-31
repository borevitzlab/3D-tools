#!/usr/bin/env python3
"""
`cloudtools` filters and optionally combines .ply pointclouds.

Units are assumed to be meters in the UTM coordinate system, with some
offset to avoid float32 imprecision.  This matches the output of
software such as Pix4D.
"""

import argparse
import glob
import math
import os
import time

import numpy as np

from geoply import GeoPly


def terrain(verts, cellsize):
    """An array describing the terrain and trees in the pointcloud."""
    # FIXME: replace with faster implementation based on forestutils
    other_names = tuple(n for n in verts.dtype.names if n not in 'xyz')
    # set up terrain array
    shape = tuple(math.ceil((verts[c].max() - verts[c].min() + 1) /
                            cellsize) for c in 'xy')
    dtype = [('points', 'i4'), ('ground', 'f4'), ('canopy', 'f4')] + [
        (name, 'f8' if 'f' in dt else ('u4' if 'u' in dt else 'i4'))
        for name, dt in verts.dtype.descr if name not in 'xyz']
    arr = np.zeros(shape, dtype=dtype)
    arr['ground'] = np.inf
    arr['canopy'] = -np.inf
    for idx, vert in zip(indices_for(verts, cellsize), verts):
        tx = arr[idx['x index'], idx['y index']]
        tx['points'] += 1
        z = vert['z']
        if z < tx['ground']:
            tx['ground'] = z
        if z > tx['canopy']:
            tx['canopy'] = z
        for n in other_names:
            tx[n] += vert[n]
    return arr


def indices_for(verts, cellsize):
    """Return index array for the given verticies."""
    # XY->idx conversion is done with array operations for speed
    indicies = np.zeros((verts.size,),
                        dtype=[('x index', 'u2'), ('y index', 'u2')])
    indicies['x index'] = (verts['x'] - verts['x'].min()) // cellsize
    indicies['y index'] = (verts['y'] - verts['y'].min()) // cellsize
    return indicies


def filter_ply(geoply, cellsize, grassdepth):
    """Filter undesirable points out of the geoply."""
    # TODO:  filter on polygon boundary (GeoJSON) by marking terrain
    if grassdepth <= 0:
        return geoply
    verts = geoply['vertex'].data
    ground = np.array(
        [terrain(verts, cellsize)[idx['x index'], idx['y index']]['ground']
         for idx in indices_for(verts, cellsize)])
    geoply['vertex'].data = verts[np.equal(verts['z'], ground) |
                                  ((verts['z'] - ground) > grassdepth)]
    return geoply


def add_clouds(args):
    """A tester function to add things togther."""
    print('Concatenating input files...')
    output_file = os.path.join(args.output_dir, 'out.ply')
    t1 = time.time()
    in_ = [GeoPly.read(f) for f in args.input_glob]
    out = GeoPly.from_geoplys(*in_) if len(in_) > 1 else in_[0]
    t2 = time.time()
    print('Concat took {:.3f} seconds'.format(t2-t1))
    filter_ply(out, args.cellsize, args.grassdepth).write(output_file)
    print('{:.3f} seconds for terrain and write'.format(time.time()-t2))


def get_args():
    """Handle command-line arguments, including default values.

    The goal is to make argparse do as much of the work as possible, eg with
    type validation functions, to 'fail fast' and simplify other code.
    """
    def glob_arg_type(val):
        """Custom type converter for argparse."""
        if not val.endswith('.ply'):
            raise argparse.ArgumentTypeError('glob pattern must end ".ply"')
        files = tuple(glob.glob(val))
        if not files:
            raise argparse.ArgumentTypeError('glob pattern matchs no files')
        return files

    def dir_arg_type(val):
        """Custom type converter for argparse."""
        if not os.path.exists(val):
            raise argparse.ArgumentTypeError('path does not exist')
        if not os.path.isdir(val):
            raise argparse.ArgumentTypeError('path is not a directory')
        return val

    df = ' default=%(default)s'
    parser = argparse.ArgumentParser(description=__doc__)

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
    add_clouds(get_args())
