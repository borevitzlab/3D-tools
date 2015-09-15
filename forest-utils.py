#!/usr/bin/env python3
"""Tools for analysing forest point clouds."""

import argparse
import math
import os.path
import shutil

import pointcloudfile
from utm_convert import UTM_to_LatLon


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
        self.offset = pointcloudfile.UTM_offset_for(self.file, args.utmzone)
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
        x = math.floor(pos[0] / args.cellsize)
        y = math.floor(pos[1] / args.cellsize)
        return x, y

    def is_ground(self, point):
        """Returns boolean whether the point is not classified as ground - ie
        True if within GROUND_DEPTH of the lowest point in the cell.
        If not lossy, also true for lowest ground point in a cell."""
        return point[2] - self.ground[self.coords(point)] < args.grounddepth

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
            OK = 3 > sum(abs(height - n) > 2*args.cellsize for n in adjacent)
            if not OK:
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
            if len(problematic) / prev > 0.9:
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
        """Lat, Lon, UTM X, UTM Y, height, area, red, green, blue, points."""
        x = self.offset.x + (max(k[0] for k in keys) +
                             min(k[0] for k in keys)) * args.cellsize / 2
        y = self.offset.y + (max(k[1] for k in keys) +
                             min(k[1] for k in keys)) * args.cellsize / 2
        out = list(UTM_to_LatLon(x, y, args.utmzone))
        out.extend([x, y, args.utmzone])
        height, r, g, b, points = 0, 0, 0, 0, 0
        for k in keys:
            height = max(height, self.canopy[k] - self.ground[k])
            r, g, b = (a+b for a, b in zip(
                [r, g, b], self.colours.get(k, [0, 0, 0])))
            points += self.density[k]
        out.extend([height, len(keys) * args.cellsize**2,
                    r // points, g // points, b // points, points])
        return tuple(out)

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
            if data[4] > 2*args.slicedepth and data[-1] > 100:
                # Filter trees by height and point count
                yield data

    def save_sparse_cloud(self, new_fname, lowest=True, canopy=True):
        """Yield points for a sparse point cloud, eliminating ~3/4 of all
        points without affecting analysis."""
        newpoints = (point for point in pointcloudfile.read(self.file)
                     if canopy and not self.is_ground(point) or
                     lowest and self.is_lowest(point))
        pointcloudfile.write(newpoints, new_fname, self.offset)
        if os.path.isfile(self.file[:-4] + '_ply_offset.xyz'):
            shutil.copy(self.file[:-4] + '_ply_offset.xyz',
                        new_fname[:-4] + '_ply_offset.xyz')
        self.file = new_fname

    def save_individual_trees(self):
        """Prototype function to save single trees to files."""
        if not args.savetrees:
            return
        if os.path.isfile(args.savetrees):
            raise IOError('Output dir for trees is already a file')
        if not os.path.isdir(args.savetrees):
            os.makedirs(args.savetrees)
        tree_to_file = {v: os.path.join(args.savetrees, 'tree_'+str(v)+'.ply')
                        for v in set(self.trees.values())}
        newpoints = (point for point in pointcloudfile.read(self.file)
                     if not self.is_ground(point))
        # TODO: rewrite with one-to-many generators, to save memory
        points_saved = {v: set() for v in tree_to_file}
        for p in newpoints:
            c = self.coords(p)
            v = self.trees.get(c)
            if v is not None:
                points_saved[v].add(p)
        for v, points in points_saved.items():
            if len(points) > 100:
                pointcloudfile.write(points, tree_to_file[v])


def match_across_takes(trees, previous_trees=None):
    """Applies ID strings to trees based on a previous dataset."""
    if previous_trees is None:
        return [['unknown'] + list(t) for t in trees]
    # TODO: implement this once a second dataset is available.
    raise NotImplementedError('Cannot import names from nonexistent data.')
    # inputs are lists of (name, utmx, utmy)
    for _, x1, y1 in trees:
        dist = 2 * 10**2
        name = 'unknown'
        for t in trees:
            n, x2, y2 = t
            dxy = (x1-x2)**2 + (y1-y2)**2
            if dist > dxy:
                dist = dxy
                t[0] = n
    return trees


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
    """Handle command-line arguments, including default values."""
    parser = argparse.ArgumentParser(
        description=('Takes a .ply forest  point cloud; outputs a sparse point'
                     'cloud and a .csv file of attributes for each tree.'))
    parser.add_argument(
        'file', help='name of the file to process')
    parser.add_argument(
        'out', default='.', nargs='?',
        help='directory for output files (optional)')
    parser.add_argument(
        '--savetrees', default='', nargs='?',
        help='where to save individual trees (default "", not saved)')
    parser.add_argument(  # analysis scale
        '--cellsize', default=0.1, nargs='?',
        help='grid scale; optimal at ~10x point spacing')
    parser.add_argument(  # georeferenced location
        '--utmzone', default=55, nargs='?',
        help='the UTM coordinate zone for georeferencing')
    parser.add_argument(  # feature extraction
        '--joinedcells', default=3, nargs='?',
        help='use cells X times larger to detect gaps between trees')
    parser.add_argument(  # feature extraction
        '--slicedepth', default=0.6, nargs='?',
        help='max tree width >= this for feature extraction')
    parser.add_argument(  # feature classification
        '--grounddepth', default=0.2, nargs='?',
        help='depth to omit from sparse point cloud')
    return parser.parse_args()

# args are globally available
args = get_args()
if not os.path.isfile(args.file):
    raise IOError('Input file not found, ' + args.file)


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

if __name__ == '__main__':
    main_processing()
