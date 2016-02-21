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
# TODO:  index point attributess by name only, not position

import argparse
import csv
import math
import os
from typing import Tuple, Set

from . import pointcloudfile, utm_convert


# User-defined types
XY_Coord = Tuple[float, float]  # float is duck-type compatible with int too


def coords(pos):
    """Return a tuple of integer coordinates as keys for the dict/map.
    * pos can be a full point tuple, or just (x, y)
    * use floor() to avoid imprecise float issues"""
    x = math.floor(pos.x / args.cellsize)
    y = math.floor(pos.y / args.cellsize)
    return x, y


class MapObj:
    """Stores a maximum and minimum height map of the cloud, in GRID_SIZE
    cells.  Hides data structure and accessed through coordinates."""
    # pylint:disable=too-many-instance-attributes

    def __init__(self, input_file, *, colours=True, zone=55, south=True):
        """
        Args:
            input_file (path): the ``.ply`` file to process.  If dealing with
                Pix4D outputs, ``*_part_1.ply``.
            colours (bool): whether to read colours from the file.  Set to
                False for eg. LIDAR data where mean colour is not useful.
            prev_csv (path): path to a csv file which associates a name
                with coordinates, to correctly name detected trees.
            zone (int): the UTM zone of the site.
            south (bool): if the site is in the southern hemisphere.
        """
        self.input_file = input_file
        self.file = input_file
        self.canopy = dict()
        self.density = dict()
        self.ground = dict()
        self.colours = dict()
        self.trees = dict()

        x, y, _ = pointcloudfile.offset_for(self.file)
        self.offset = utm_convert.UTM_coords(x, y, zone, south)

        self.update_spatial()
        if colours:
            self.update_colours()

    def update_spatial(self):
        """Expand, correct, or maintain map with a new observed point."""
        # Fill out the spatial info in the file
        for p in pointcloudfile.read(self.file):
            idx = coords(p)
            if self.density.get(idx) is None:
                self.density[idx] = 1
                self.canopy[idx] = p.z
                self.ground[idx] = p.z
                continue
            self.density[idx] += 1
            if self.ground[idx] > p.z:
                self.ground[idx] = p.z
            elif self.canopy[idx] < p.z:
                self.canopy[idx] = p.z
        self.__smooth_ground()
        self.trees = self._tree_components()

    def update_colours(self):
        """Expand, correct, or maintain map with a new observed point."""
        # We assume that vertex attributes not named "x", "y" or "z"
        # are colours, and thus accumulate a total to get the mean
        for p in pointcloudfile.read(self.file):
            if self.is_ground(p):
                continue
            p_cols = {k: v for k, v in p._asdict().items() if k not in 'xyz'}
            idx = coords(p)
            if idx not in self.colours:
                self.colours[idx] = p_cols
            else:
                for k, v in p_cols.items():
                    self.colours[idx][k] += v

    def is_ground(self, point) -> bool:
        """Returns boolean whether the point is not classified as ground - ie
        True if within GROUND_DEPTH of the lowest point in the cell.
        If not lossy, also true for lowest ground point in a cell."""
        return point[2] - self.ground[coords(point)] < args.grounddepth

    def is_lowest(self, point) -> bool:
        """Returns boolean whether the point is lowest in that grid cell."""
        return point[2] == self.ground[coords(point)]

    def __len__(self) -> int:
        """Total observed points."""
        return sum(self.density.values())

    @staticmethod
    def _neighbors(key):
        """Return the adjacent keys, whether they exist or not."""
        return ((key[0]+a, key[1]+b) for a, b in
                [(1, 1), (1, 0), (1, -1),
                 (0, 1), (0, -1),
                 (-1, 1), (-1, 0), (-1, -1)])

    def __problematic(self, prior: set=None) -> Set[XY_Coord]:
        """Identifies cells with more than 2:1 slope to 3+ adjacent cells."""
        problematic = set()
        if prior is None:
            prior = self.ground.keys()
        for k in prior:
            adjacent = {self.ground.get(n) for n in self._neighbors(k)}
            adjacent.discard(None)
            if len(adjacent) < 6:
                continue
            height = self.ground[k]
            # Number of cells at more than 2:1 slope - suspiciously steep.
            # 3+ usually indicates a misclassified cell or data artefact.
            probs = sum(abs(height - n) > 2*args.cellsize for n in adjacent)
            if probs >= 3:
                problematic.add(k)
        return problematic

    def __smooth_ground(self):
        """Smooths the ground map, to reduce the impact of spurious points, eg.
        points far underground or misclassification of canopy as ground."""
        problematic = self.__problematic()
        for _ in range(10):
            for key in problematic:
                adjacent = {self.ground.get(n) for n in self._neighbors(key)
                            if n not in problematic}
                adjacent.discard(None)
                if not adjacent:
                    continue
                self.ground[key] = min(adjacent) + 2*args.cellsize
            prev = len(problematic)
            problematic = self.__problematic(problematic)
            if prev == 0 or len(problematic) / prev > 0.9:
                break

    def _tree_components(self):
        """Returns a dict where keys refer to connected components.
        NB: Not all keys in other dicts exist in this output."""
        # pylint:disable=too-many-branches
        def expand(old_key, com):
            """Implement depth-first search."""
            for key in self._neighbors(old_key):
                if com.get(key) is None:
                    return
                if com[key] == com[old_key]:
                    continue
                elif com[key] < com[old_key]:
                    com[old_key] = com[key]
                else:
                    com[key] = com[old_key]
                    expand(key, com)

        # Set up a boolean array of larger keys to search
        key_scale_record = dict()
        for key in self.density.keys():
            if self.canopy[key] - self.ground[key] > args.slicedepth:
                cc_key = tuple(math.floor(k/args.joinedcells) for k in key)
                if cc_key not in key_scale_record:
                    key_scale_record[cc_key] = {key}
                else:
                    key_scale_record[cc_key].add(key)
        # Assign a unique integer value to each large key, then search
        # Final labels are positive ints, but not ordered or consecutive
        com = dict()
        for idx, key in enumerate(tuple(key_scale_record)):
            com[key] = idx
        for key in tuple(com):
            try:
                expand(key, com)
            except RuntimeError:
                # Recursion depth; finish on next pass.
                continue
        # Copy labels to grid of original scale
        ret = dict()
        for key, values in key_scale_record.items():
            for skey in values:
                ret[skey] = com[key]
        return ret

    def tree_data(self, keys):
        """Return a dictionary of data about the tree in the given keys."""
        out = dict()
        # Calculate positional information
        out['UTM_X'] = x = self.offset.y + args.cellsize * (
            max(k[0] for k in keys) + min(k[0] for k in keys)) / 2
        out['UTM_Y'] = y = self.offset.y + args.cellsize * (
            max(k[1] for k in keys) + min(k[1] for k in keys)) / 2
        out['UTM_zone'] = args.utmzone
        out['latitude'], out['longitude'] = utm_convert.UTM_to_LatLon(
            x, y, args.utmzone)
        out['area'] = len(keys) * args.cellsize**2
        out['base_altitude'] = sum(self.ground[k] for k in keys) / len(keys)
        # Loop over grid cells to accumulate other information
        for k in keys:
            out['height'] = max(out.get('height', 0),
                                self.canopy[k] - self.ground[k])
            out['point_count'] = out.get('point_count', 0) + self.density[k]
        # TODO:  handle more flexible colour channels
        # Finally, normalise colour values and round position to confidence
        for k in ['height', 'area', 'base_altitude']:
            out[k] = round(out[k], 2)
        for k in ['UTM_X', 'UTM_Y']:
            out[k] = round(out[k], 1)
        return out

    def all_trees(self):
        """Yield the characteristics of each tree."""
        ids = list(set(self.trees.values()))
        keys = {v: set() for v in ids}
        for k, v in self.trees.items():
            if v is None:
                continue
            keys[v].add(k)
        for v in ids:
            data = self.tree_data(keys[v])
            if data['height'] > 1.5 * args.slicedepth:
                # Filter trees by height
                yield data

    def save_sparse_cloud(self, new_fname, lowest=True, canopy=True):
        """Yield points for a sparse point cloud, eliminating ~3/4 of all
        points without affecting analysis."""
        newpoints = (point for point in pointcloudfile.read(self.file)
                     if canopy and not self.is_ground(point) or
                     lowest and self.is_lowest(point))
        pointcloudfile.write(newpoints, new_fname, self.input_file)
        if lowest and canopy:
            self.file = new_fname

    def save_individual_trees(self):
        """Save single trees to files."""
        if not args.savetrees:
            return
        if os.path.isfile(args.savetrees):
            raise IOError('Output dir for trees is already a file')
        if not os.path.isdir(args.savetrees):
            os.makedirs(args.savetrees)
        # Maps tree ID numbers to a incremental writer for that tree
        tree_to_file = {tree_ID: pointcloudfile.IncrementalWriter(
            os.path.join(args.savetrees, 'tree_{}.ply'.format(tree_ID)),
            self.input_file) for tree_ID in set(self.trees.values())}
        # For non-ground, find the appropriate writer and call with the point
        for point in pointcloudfile.read(self.file):
            val = self.trees.get(coords(point))
            if val is not None:
                tree_to_file[val](point)


def stream_analysis(attr, out: str) -> None:
    """Saves the list of trees with attributes to the file 'out'."""
    header = ('latitude', 'longitude', 'UTM_X', 'UTM_Y', 'UTM_zone',
              'height', 'area', 'base_altitude', 'point_count')
    with open(out, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        for data in attr.all_trees():
            writer.writerow(data)


def get_args():
    """Handle command-line arguments, including default values."""
    parser = argparse.ArgumentParser(
        description=('Takes a .ply forest  point cloud; outputs a sparse point'
                     'cloud and a .csv file of attributes for each tree.'))
    parser.add_argument(
        'file', help='name of the file to process', type=str)
    parser.add_argument(
        'out', default='.', nargs='?', type=str,
        help='directory for output files (optional)')
    parser.add_argument(
        '--savetrees', default='', nargs='?', type=str,
        help='where to save individual trees (default "", not saved)')
    parser.add_argument(  # analysis scale
        '--cellsize', default=0.1, nargs='?', type=float,
        help='grid scale; optimal at ~10x point spacing')
    parser.add_argument(  # georeferenced location
        '--utmzone', default=55, type=int,
        help='the UTM coordinate zone for georeferencing')
    parser.add_argument(  # feature extraction
        '--joinedcells', default=3, type=float,
        help='use cells X times larger to detect gaps between trees')
    parser.add_argument(  # feature extraction
        '--slicedepth', default=0.6, type=float,
        help='slice depth for canopy area and feature extraction')
    parser.add_argument(  # feature classification
        '--grounddepth', default=0.2, type=float,
        help='depth to omit from sparse point cloud')
    return parser.parse_args()


def main_processing():
    """Logic on which functions to call, and efficient order."""
    print('Reading from "{}" ...'.format(args.file))
    sparse = os.path.join(args.out, os.path.basename(args.file))
    if not args.file.endswith('_sparse.ply'):
        sparse = os.path.join(
            args.out, os.path.basename(args.file)[:-4] + '_sparse.ply')
    sparse = sparse.replace('_part_1', '')
    if os.path.isfile(sparse):
        attr_map = MapObj(sparse)
        print('Read {} points into {} cells'.format(
            len(attr_map), len(attr_map.canopy)))
    else:
        attr_map = MapObj(args.file, colours=False)
        print('Read {} points into {} cells, writing "{}" ...'.format(
            len(attr_map), len(attr_map.canopy), sparse))
        attr_map.save_sparse_cloud(sparse)
        print('Reading colours from ' + sparse)
        attr_map.update_colours()
    print('File IO complete, starting analysis...')
    table = '{}_analysis.csv'.format(sparse[:-4].replace('_sparse', ''))
    stream_analysis(attr_map, table)
    if args.savetrees:
        print('Saving individual trees...')
        attr_map.save_individual_trees()
    print('Done.')


def main():
    """Interface to call from outside the package."""
    # pylint:disable=global-statement
    global args
    args = get_args()
    if not os.path.isfile(args.file):
        raise IOError('Input file not found, ' + args.file)
    main_processing()

if __name__ == '__main__':
    # Call to get_args is duplicated to work in static analysis, from
    # command line, and when installed as package (calls main directly)
    args = get_args()
    main()
