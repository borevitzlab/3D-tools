#!/usr/bin/env python3
"""Tools for analysing forest point clouds.

Inputs: a coloured pointcloud in ``.ply`` format (XYZRGB vertices), which
can be obtained by putting images from drone photography through
structure-from-motion software.

Outputs (most are optional):

* A 'sparse' point cloud, with most ground points discarded.  This eases
  further analysis, storage, etc without compromising coverage of vegetation.
* A ``.csv`` file listing identified trees, with traits including location,
  height, canopy area, colour, and point count.
* Individual files containing the point cloud for each tree.

Extensive use of mutable coordinate-property mappings and streamed input
ensure that even files to large to load in memory can be processed.  In extreme
cases, the resolution can be decreased to trade accuracy for memory.

Example outputs (from an older version):
`a map <https://www.google.com/maps/d/viewer?mid=z1pH7HaTWL9Q.kzQflQGYVRIU>`_,
and `pointclouds <http://phenocam.org.au/pointclouds>`_.
"""
# pylint:disable=unsubscriptable-object

import argparse
import collections
import csv
import math
import os
from typing import MutableMapping, NamedTuple, Tuple, Set

import utm  # type: ignore

import geoply


# User-defined types
XY_Coord = NamedTuple('XY_Coord', [('x', int), ('y', int)])
Coord_Labels = MutableMapping[XY_Coord, int]


def coords(pos, cellsize):
    """Return a tuple of integer coordinates as keys for the dict/map.
    * pos can be a full point tuple, or just (x, y)
    * use floor() to avoid imprecise float issues"""
    return XY_Coord(
        math.floor(pos['x'] / cellsize), math.floor(pos['y'] / cellsize))


def neighbors(key: XY_Coord) -> Tuple[XY_Coord, ...]:
    """Return the adjacent keys, whether they exist or not."""
    return tuple(XY_Coord(key.x + a, key.y + b)
                 for a in (-1, 0, 1) for b in (-1, 0, 1) if a or b)


def connected_components(input_dict: Coord_Labels) -> None:
    """Connected components in a dict of coordinates.

    Uses depth-first search.  Non-component cells are absent from the input.
    """
    def expand(old_key: XY_Coord, com: MutableMapping) -> None:
        """Implement depth-first search."""
        for key in neighbors(old_key):
            if com.get(key) is None:
                return
            if com[key] == com[old_key]:
                continue
            elif com[key] < com[old_key]:
                com[old_key] = com[key]
            else:
                com[key] = com[old_key]
                expand(key, com)

    for key in tuple(input_dict):
        try:
            expand(key, input_dict)
        except RuntimeError:
            # Recursion depth; finish on next pass.
            continue


class MapObj:
    """Stores a maximum and minimum height map of the cloud, in GRID_SIZE
    cells.  Hides data structure and accessed through coordinates."""
    # pylint:disable=too-many-instance-attributes

    def __init__(self, input_file, *, args, colours=True):
        """
        Args:
            input_file (path): the ``.ply`` file to process.  If dealing with
                Pix4D outputs, ``*_part_1.ply``.
            args (namespace): parsed command-line arguments.
            colours (bool): whether to read colours from the file.  Set to
                False for eg. LIDAR data where mean colour is not useful.
        """
        self.__file = input_file
        self.canopy = dict()
        self.density = dict()
        self.ground = dict()
        self.colours = collections.defaultdict(
            lambda: collections.defaultdict(int))
        self.trees = dict()

        self.geoply = geoply.GeoPly.read(input_file)

        self.cellsize = args.cellsize
        self.slicedepth = args.slicedepth
        self.grounddepth = args.grassdepth
        self.joinedcells = args.joinedcells

        self.update_spatial()
        if colours:
            self.update_colours()

    def update_spatial(self):
        """Expand, correct, or maintain map with a new observed point."""
        # Fill out the spatial info in the file
        for p in self.geoply.vertices:
            idx = coords(p, self.cellsize)
            if self.density.get(idx) is None:
                self.density[idx] = 1
                self.canopy[idx] = p['z']
                self.ground[idx] = p['z']
                continue
            self.density[idx] += 1
            if self.ground[idx] > p['z']:
                self.ground[idx] = p['z']
            elif self.canopy[idx] < p['z']:
                self.canopy[idx] = p['z']
        self.trees = self._tree_components()

    def update_colours(self):
        """Expand, correct, or maintain map with a new observed point."""
        # We assume that vertex attributes not named "x", "y" or "z"
        # are colours, and thus accumulate a total to get the mean
        for p in self.geoply.vertices:
            if self.is_ground(p):
                continue
            idx = coords(p, self.cellsize)
            for name, value in zip(p.dtype.names, p):
                if name not in 'xyz':
                    self.colours[idx][name] += value.item()

    def is_ground(self, point) -> bool:
        """Returns boolean whether the point is not classified as ground - ie
        True if within GROUND_DEPTH of the lowest point in the cell.
        If not lossy, also true for lowest ground point in a cell."""
        return point['z'] - self.ground[coords(point, self.cellsize)] \
            < self.grounddepth

    def is_lowest(self, point) -> bool:
        """Returns boolean whether the point is lowest in that grid cell."""
        return point['z'] == self.ground[coords(point, self.cellsize)]

    def __len__(self) -> int:
        """Total observed points."""
        return sum(self.density.values())

    def _tree_components(self) -> Coord_Labels:
        """Returns a dict where keys refer to connected components.
        NB: Not all keys in other dicts exist in this output."""
        # Set up a boolean array of larger keys to search
        key_scale_record = {}  # type: Dict[XY_Coord, Set[XY_Coord]]
        for key in self.density:
            if self.canopy[key] - self.ground[key] > self.slicedepth:
                cc_key = XY_Coord(int(math.floor(key.x / self.joinedcells)),
                                  int(math.floor(key.y / self.joinedcells)))
                if cc_key not in key_scale_record:
                    key_scale_record[cc_key] = {key}
                else:
                    key_scale_record[cc_key].add(key)
        # Assign a unique integer value to each large key, then search
        # Final labels are positive ints, but not ordered or consecutive
        trees = {k: i for i, k in enumerate(tuple(key_scale_record))}
        connected_components(trees)
        # Copy labels to grid of original scale
        return {s: trees[k] for k, v in key_scale_record.items() for s in v}

    def tree_data(self, keys: Set[XY_Coord]) -> dict:
        """Return a dictionary of data about the tree in the given keys."""
        # Calculate positional information
        utmx, utmy, zone, northern = self.geoply.utm_coord
        x = utmx + args.cellsize * sum(k.x for k in keys) / len(keys)
        y = utmy + args.cellsize * sum(k.y for k in keys) / len(keys)
        lat, lon = utm.to_latlon(x, y, zone, northern=northern)
        out = {
            'latitude': lat,
            'longitude': lon,
            'UTM_X': x,
            'UTM_Y': y,
            'UTM_zone': args.utmzone,
            'height': 0,
            'area': len(keys) * args.cellsize**2,
            'base_altitude': sum(self.ground[k] for k in keys) / len(keys),
            'point_count': 0,
            }
        colours = collections.defaultdict(list)
        for k in keys:
            out['height'] = max(out['height'], self.canopy[k] - self.ground[k])
            out['point_count'] += self.density[k]
            for colour, total in self.colours[k].items():
                colours[colour].append(total / self.density[k])
        for name, values in colours.items():
            out[name] = sum(values) / len(values)
        return out

    def all_trees(self):
        """Yield the characteristics of each tree."""
        # Group coord-keys by tree ID; pass each set of coords to tree_data()
        tree_coords = collections.defaultdict(set)
        for coord, tree_id in self.trees.items():
            if tree_id is not None:
                tree_coords[tree_id].add(coord)
        for coord_group in tree_coords.values():
            data = self.tree_data(coord_group)
            if data['height'] > 1.5 * self.slicedepth:
                yield data

    def save_sparse_cloud(self, new_fname, lowest=True, canopy=True):
        """Yield points for a sparse point cloud, eliminating ~3/4 of all
        points without affecting analysis."""
        assert lowest or canopy, 'Must keep at least one of canopy or lowest'
        newpoints = (point for point in self.geoply.vertices
                     if canopy and not self.is_ground(point) or
                     lowest and self.is_lowest(point))
        output = geoply.GeoPly.from_iterable(newpoints, self.geoply.utm_coord)
        output.write(new_fname)
        if lowest and canopy:
            self.geoply = output

    def save_individual_trees(self, trees_dir):
        """Save single trees to files."""
        # Validate the directory to be written to
        if not trees_dir:
            print('ERROR: cannot save trees to nonexistent path')
            return
        if os.path.isfile(trees_dir):
            raise IOError('Output dir for trees is already a file')
        if not os.path.isdir(trees_dir):
            os.makedirs(trees_dir)
        # Group points into dict of ID number to list of points
        # This is relatively memory-intensive.
        tree_clouds = collections.defaultdict(list)
        for point in self.geoply.vertices:
            tree = self.trees.get(coords(point, self.cellsize))
            if tree is not None:
                tree_clouds[tree].append(point)
        # Write them to disk
        for tree, points in tree_clouds.items():
            out = geoply.GeoPly.from_iterable(points, self.geoply.utm_coord)
            fname = 'tree_{}.ply'.format(tree)
            out.write(os.path.join(trees_dir, fname))


    def stream_analysis(self, out: str) -> None:
        """Saves the list of trees with attributes to the file 'out'."""
        header = ('latitude', 'longitude', 'UTM_X', 'UTM_Y', 'UTM_zone',
                  'height', 'area', 'base_altitude', 'point_count') + tuple(
                      a for a in self.geoply.vertices.dtype.names
                      if a not in 'xyz')
        with open(out, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=header)
            writer.writeheader()
            for data in self.all_trees():
                writer.writerow(data)


def get_args():
    """Handle command-line arguments, including default values."""
    df = ' default=%(default)s'
    parser = argparse.ArgumentParser(
        description=__doc__.strip(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Contact:  zac.hatfield.dodds@anu.edu.au',
        )

    io_g = parser.add_argument_group('Input/Output Args')
    io_g.add_argument(
        'file', help='the file or files to process, ending ".ply"')
    io_g.add_argument(
        '--out', default='.', help='directory for output files')
    parser.add_argument(
        '--savetrees', metavar='DIR',
        help='where to save individual trees (default: not saved)')

    filt_g = parser.add_argument_group('Filtering Args')
    filt_g.add_argument(
        '--cellsize', default=0.1, type=float,
        help=('edge length of the analysis pixels in meters.  Around 10 times'
              ' point spacing seems to give good results.' + df))
    filt_g.add_argument(
        '--grassdepth', metavar='DEPTH', default=0.2, type=float,
        help=('filter out points less than this height in meters above lowest'
              ' point per cell.  n <= 0 disables this filter.' + df))
    filt_g.add_argument(
        '--joinedcells', metavar='DEPTH', default=3, type=float,
        help='use cells X times larger to detect gaps between trees')
    filt_g.add_argument(
        '--slicedepth', default=0.6, type=float,
        help='slice depth for canopy area and feature extraction')

    geo_g = parser.add_argument_group('Geolocation Args')
    geo_g.add_argument(
        '--utmzone', default=55, type=int, choices=range(1, 61),
        metavar='[1..60]', help='the UTM coordinate zone of the data,' + df)
    ns_g = geo_g.add_mutually_exclusive_group()
    ns_g.add_argument('--northern', action='store_true',
                      help='for data in the northern hemisphere,' + df)
    ns_g.add_argument('--southern', action='store_true', default=True,
                      help='for data in the southern hemisphere,' + df)

    args = parser.parse_args()

    if not os.path.isfile(args.file):
        parser.error('Input file not found, ' + args.file)
    if not args.file.endswith('.ply'):
        parser.error('Input file requires ".ply" extension')
    if not os.path.isdir(args.out):
        parser.error('Output path does not exist: ' + args.out)

    return args


def main():
    """Logic on which functions to call, and efficient order."""
    args = get_args()
    print('Reading from "{}" ...'.format(args.file))
    sparse = os.path.join(args.out, os.path.basename(args.file))
    if not args.file.endswith('_sparse.ply'):
        sparse = os.path.join(
            args.out, os.path.basename(args.file)[:-4] + '_sparse.ply')
    sparse = sparse.replace('_part_1', '')
    if os.path.isfile(sparse):
        attr_map = MapObj(sparse, args=args)
        print('Read {} points into {} cells'.format(
            len(attr_map), len(attr_map.canopy)))
    else:
        attr_map = MapObj(args.file, args=args, colours=False)
        print('Read {} points into {} cells, writing "{}" ...'.format(
            len(attr_map), len(attr_map.canopy), sparse))
        attr_map.save_sparse_cloud(sparse)
        print('Reading colours from ' + sparse)
        attr_map.update_colours()
    print('File IO complete, starting analysis...')
    table = '{}_analysis.csv'.format(sparse[:-4].replace('_sparse', ''))
    attr_map.stream_analysis(table)
    if args.savetrees:
        print('Saving individual trees...')
        attr_map.save_individual_trees(args.savetrees)
    print('Done.')


if __name__ == '__main__':
    main()
