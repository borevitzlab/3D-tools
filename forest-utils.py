#!/usr/bin/env python3
"""Tools for analysing forest point clouds."""

import argparse
import math
import os.path
import shutil

import pointcloudfile
from utm_convert import UTM_to_LatLon

# Location, defaults to UTM zone 55
UTM_ZONE = 55

# Default to 10cm squares for analysis; optimum depends on cloud density
CELL_SIZE = 0.1

# Classifies anything in the lowest GROUND_DEPTH in CELL_SIZE as ground
GROUND_DEPTH = 0.2

# SLICE_DEPTH is used in the connected component checks; max width of canopy
# will be >= this.
SLICE_DEPTH = 0.6

# The component analysis uses cells JOINED_CELLS * CELL_SIZE on each side
JOINED_CELLS = 4


class MapObj(object):
    """Stores a maximum and minimum height map of the cloud, in GRID_SIZE
    cells.  Hides data structure and accessed through coordinates."""
    def __init__(self, input_file, colours=True):
        """Initialises data"""
        self.canopy = dict()
        self.density = dict()
        self.ground = dict()
        self.colours = dict()
        self.trees = dict()
        self.file = input_file
        self.offset = pointcloudfile.UTM_offset_for(self.file)
        self.update_spatial()
        if colours:
            self.update_colours()

    def update_spatial(self):
        """Expand, correct, or maintain map with a new observed point."""
        # Fill out the spatial info in the file
        for point in pointcloudfile.read(self.file):
            z = point[2]
            idx = self.coords(point)
            if self.density.get(idx) is None:
                self.density[idx] = 1
                self.canopy[idx] = z
                self.ground[idx] = z
                continue
            self.density[idx] += 1
            if self.ground[idx] > z:
                self.ground[idx] = z
            elif self.canopy[idx] < z:
                self.canopy[idx] = z
        self.__smooth_ground()
        self.trees = self._tree_components()

    def update_colours(self):
        """Expand, correct, or maintain map with a new observed point."""
        for x, y, z, r, g, b in pointcloudfile.read(self.file):
            if self.is_ground([x, y, z]):
                continue
            idx = self.coords([x, y])
            R, G, B = self.colours.get(idx, (0, 0, 0))
            self.colours[idx] = (R+r, G+g, B+b)

    @staticmethod
    def coords(pos):
        """Return a tuple of integer coordinates as keys for the dict/map.
        * pos can be a full point tuple, or just (x, y)
        * use floor() to avoid imprecise float issues"""
        return math.floor(pos[0] / CELL_SIZE), math.floor(pos[1] / CELL_SIZE)

    def is_ground(self, point):
        """Returns boolean whether the point is not classified as ground - ie
        True if within GROUND_DEPTH of the lowest point in the cell.
        If not lossy, also true for lowest ground point in a cell."""
        return point[2] - self.ground[self.coords(point)] < GROUND_DEPTH

    def is_lowest(self, point):
        """Returns boolean whether the point is lowest in that grid cell."""
        return point[2] == self.ground[self.coords(point)]

    def __len__(self):
        """Total observed points."""
        if not self.density:
            return 0
        return sum(self.density.values())

    @staticmethod
    def _neighbors(key):
        """Return the adjacent keys, whether they exist or not."""
        return ((key[0]+a, key[1]+b) for a, b in
                [(1, 1), (1, 0), (1, -1),
                 (0, 1), (0, -1),
                 (-1, 1), (-1, 0), (-1, -1)])

    def __problematic(self, prior=None):
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
            OK = 3 > sum(abs(height - n) > 2*CELL_SIZE for n in adjacent)
            if not OK:
                problematic.add(k)
        return problematic

    def __smooth_ground(self):
        """Smooths the ground map, to reduce misclassification of canopy as
        ground, and reduce the impact of spurious points."""
        problematic = self.__problematic()
        p1 = len(problematic)
        for _ in range(10):
            # Smooth problematic cells where practical
            for key in problematic:
                adjacent = {self.ground.get(n) for n in self._neighbors(key)
                            if not n in problematic}
                adjacent.discard(None)
                if not adjacent:
                    continue
                self.ground[key] = min(adjacent) + 2*CELL_SIZE
            p2 = len(problematic)
            problematic = self.__problematic(problematic)
            if p2/p1 < 0.2 or p2/len(problematic) > 0.8:
                break

    def _tree_components(self):
        """Returns a dict where keys refer to connected components.
        NB: Not all keys in other dicts exist in this output."""
        #pylint:disable=too-many-branches
        def expand(old_key, com):
            """Implement depth-first search."""
            for key in self._neighbors(old_key):
                if com.get(key) is None:
                    # No above-ground point in that joined cell
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
            if self.canopy[key] - self.ground[key] > SLICE_DEPTH:
                cc_key = tuple(math.floor(k/JOINED_CELLS) for k in key)
                if not cc_key in key_scale_record:
                    key_scale_record[cc_key] = {key}
                else:
                    key_scale_record[cc_key].add(key)
        # Assign a unique integer value to each large key, then search
        com = dict()
        for idx, key in enumerate(tuple(key_scale_record)):
            com[key] = idx
        for key in tuple(com):
            try:
                expand(key, com)
            except RuntimeError:
                # Recursion depth; finish on next pass.
                continue
        # Reduce component labels to consecutive integers from zero
        out = dict()
        translate = dict((v, k) for k, v in enumerate(set(com.values())))
        for k, v in com.items():
            out[k] = translate[v]
        # Copy labels to grid of original scale
        ret = dict()
        for key, values in key_scale_record.items():
            for skey in values:
                ret[skey] = out[key]
        return ret

    def tree_data(self, keys):
        """Lat, Lon UTM X, UTM Y, height, area, red, green, blue, points."""
        x = self.offset.x + (max(k[0] for k in keys) +
                             min(k[0] for k in keys)) * CELL_SIZE / 2
        y = self.offset.y + (max(k[1] for k in keys) +
                             min(k[1] for k in keys)) * CELL_SIZE / 2
        out = list(UTM_to_LatLon(x, y, UTM_ZONE))
        out.extend([x, y, UTM_ZONE])
        height, r, g, b, points = 0, 0, 0, 0, 0
        for k in keys:
            height = max(height, self.canopy[k] - self.ground[k])
            r, g, b = (a+b for a, b in zip(
                [r, g, b], self.colours.get(k, [0, 0, 0])))
            points += self.density[k]
        out.extend([height, len(keys) * CELL_SIZE**2,
                    r // points, g // points, b // points, points])
        return tuple(out)

    def all_trees(self):
        """Yield the characteristics of each tree."""
        ids = list(set(self.trees.values()))
        keys = [set() for _ in ids]
        for k, v in self.trees.items():
            if v is None:
                continue
            keys[v].add(k)
        for v in set(self.trees.values()):
            yield self.tree_data(keys[v])

def make_sparse(filename, attr=None, *, keep_lowest=True, keep_canopy=True):
    """Yield points for a sparse point cloud.
    By default, it keeps both canopy points and the equal lowest points in each
    grid cell; eliminating ~3/4 of all points without affecting analysis.  It
    can also keep canopy only or lowest only."""
    if not keep_lowest and not keep_canopy:
        raise StopIteration
    if attr is None:
        attr = MapObj(filename, colours=False)
    if attr.file != filename:
        raise ValueError('Supplied attributes are for a different file.')
    for point in pointcloudfile.read(filename):
        if keep_canopy and not attr.is_ground(point):
            yield point
        elif keep_lowest and attr.is_lowest(point):
            yield point

def match_across_takes(trees, previous_trees=None):
    """Applies ID strings to trees based on a previous dataset."""
    if previous_trees is None:
        return [['unknown'] + list(t) for t in trees]
    # TODO: implement this once a second dataset is available.
    raise NotImplementedError('Cannot import names from nonexistent data.')

def stream_analysis(attr, out):
    """Saves the list of trees with attributes to the file 'out'."""
    form_str = '"{}",{},{},{:.1f},{:.1f},{},{:.2f},{:.2f},{},{},{},{},\n'
    trees = match_across_takes(attr.all_trees(), None)
    with open(out, 'w') as f:
        f.write('name (ID), latitude, longitude, UTM_X, UTM_Y, UTM_zone, '
                'height, area, red, green, blue, point_count,\n')
        for data in trees:
            f.write(form_str.format(*data))

def get_args():
    """Handle setup, parsing, and validation.  Return args object."""
    parser = argparse.ArgumentParser(
        description='Takes a .ply forest  point cloud; ' +
        'outputs attributes of each tree to .csv')
    parser.add_argument('file', help='name of the file to process')
    parser.add_argument('out', nargs='?', default='.',
                        help='directory for output files (optional)')
    args = parser.parse_args()
    if not os.path.isfile(args.file):
        raise FileNotFoundError(args.file)
    if not args.file.endswith('.ply'):
        raise ValueError('file must be a point cloud in .ply format')
    return args

def main_processing(args):
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
        pointcloudfile.write(make_sparse(args.file, attr_map), sparse)
        if os.path.isfile(args.file[:-4] + '_ply_offset.xyz'):
            shutil.copy(args.file[:-4] + '_ply_offset.xyz',
                        sparse[:-4] + '_ply_offset.xyz')
        attr_map.file = sparse
        print('Reading colours from ' + sparse)
        attr_map.update_colours()
    print('File IO complete, starting analysis...')
    table = '{}_analysis.csv'.format(sparse[:-4].replace('_sparse', ''))
    stream_analysis(attr_map, table)
    print('Done.')

if __name__ == '__main__':
    main_processing(get_args())

