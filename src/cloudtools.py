#!/usr/bin/env python3
"""
`cloudtools` filters and optionally combines .ply pointclouds.

Units are assumed to be meters in the UTM coordinate system, with some
offset to avoid float32 imprecision.  This matches the output of
software such as Pix4D.
"""

import argparse
from collections import namedtuple
import contextlib
import glob
import os
import pkg_resources

import numpy as np
import plyfile


UTM_COORD = namedtuple('UTMCoord', ['easting', 'northing', 'zone', 'northern'])


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
        # Read file, reconstruct as GeoPly instance, drop non-vertex elements
        with open(filename, 'rb') as stream:
            data = GeoPly._parse_header(stream)
            for elt in data:
                elt._read(stream, data.text, data.byte_order)
        verts = data['vertex']
        data.elements = [verts]

        # Remove meshlab cruft
        if 'VCGLIB generated' in data.comments:
            names = verts.data.dtype.names  # field names of each vertex
            if 'alpha' in names and 1 == len(np.unique(verts['alpha'])):
                # removal of a vertex field is via fancy indexing
                verts.data = verts.data[[n for n in names if n != 'alpha']]
                # properties of the PlyElement instance are manually updated
                verts.properties = [p for p in verts.properties
                                    if p.name != 'alpha']
            data.comments.remove('VCGLIB generated')

        # Add UTM coordinates if known or discoverable
        utm_coord = UTM_COORD(0, 0, args.utmzone, args.northern)
        with contextlib.suppress(Exception):
            pass
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

        return GeoPly(data.elements, data.text, data.byte_order,
                      data.comments, data.obj_info, utm_coord)

    def write(self, stream):
        """Write to a file, serialising utm_coord as a special comment."""
        if self.utm_coord is not None:
            self.comments = [c for c in self.comments if c[:4] != 'utm_']
            self.comments.extend(
                'utm_offset {} {}\nutm_zone {}\nutm_northern {}'.format(
                    *self.utm_coord).split('\n'))
        super().write(stream)


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

    # handle the relative UTM offsets
    prime_offset = geoplys[0].utm_coord
    all_verts = []
    for pf in geoplys:
        arr = pf['vertex'].data
        if pf.utm_coord is not None:
            arr['x'] += (pf.utm_coord.easting - prime_offset.easting)
            arr['y'] += (pf.utm_coord.northing - prime_offset.northing)
        all_verts.append(arr)
    verts = plyfile.PlyElement.describe(np.concatenate(all_verts), 'vertex')

    # Load data back into the complete structure and return
    return GeoPly([verts], comments=comments, utm_coord=prime_offset)


def add_clouds(args):
    """A tester function to add things togther."""
    print('Concatenating input files...')
    output_file = os.path.join(args.output_dir, 'out.ply')
    concat(*[GeoPly.read(f) for f in args.input_glob]).write(output_file)
    print('Saved to ' + output_file)


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
        help='size of the grid; optimal at ~10x point spacing')
    filt_g.add_argument(
        '--grounddepth', default=0.2, type=float,
        help='depth above lowest point per cell to omit from output')

    geo_g = parser.add_argument_group('Geolocation Args')
    geo_g.add_argument(
        '--utmzone', default=55, type=int, choices=range(1, 61),
        metavar='[1..60]', help='the UTM coordinate zone of the data')
    ns_g = geo_g.add_mutually_exclusive_group()
    ns_g.add_argument('--northern', action='store_true',
                      help='set if the data is in the northern hemisphere')
    ns_g.add_argument('--southern', action='store_true',
                      default=True,
                      help='set if the data is in the southern hemisphere')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    add_clouds(args)
