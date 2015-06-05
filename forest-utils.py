#!/usr/bin/env python3
"""Tools for analysing forest point clouds.
"""

import argparse
import math
import os

import pointcloudfile

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
    # Until I think of something better, maps are stored as dictionaries
    # with an (x, y) tuple as keys.  x and y values are rounded down to
    # the nearest int after division by CELL_SIZE.
    # This means we can build the map before knowing it's extent
    def __init__(self, initial=None):
        """Initialises data"""
        self.canopy = dict()
        self.density = dict()
        self.ground = dict()
        self.colours = dict()
        self._components_counted_at = None
        self._components_map = None
        self.count = 0
        if initial is not None:
            for point in initial:
                self.update(point)

    def update(self, point):
        """Expand, correct, or maintain map with a new observed point."""
        self.count += 1
        x, y, z, r, g, b = point
        idx = self.coords((x, y))
        if not self.density.get(idx):
            self.density[idx] = 1
            self.canopy[idx] = z
            self.ground[idx] = z
            return
        self.density[idx] += 1
        if self.canopy[idx] < z:
            self.canopy[idx] = z
        if self.ground[idx] > z:
            self.ground[idx] = z
        R, G, B = self.colours.get(idx, (0, 0, 0))
        self.colours[idx] = (R+r, G+g, B+b)

    @staticmethod
    def coords(pos):
        """Return a tuple of integer coordinates as keys for the dict/map.
        * pos can be a full point tuple, or just (x, y)
        * use floor() to avoid nasty float issues
        * integer labels so we can do index-based stuff like finding neighbors
        """
        x, y, *_ = pos
        return math.floor(x / CELL_SIZE), math.floor(y / CELL_SIZE)

    def point_not_ground(self, point, keep_lowest=False):
        """Returns boolean whether the point is not classified as ground - ie
        True if within GROUND_DEPTH of the lowest point in the cell.
        If not lossy, also true for lowest ground point in a cell."""
        x, y, z, *_ = point
        height = z - self.ground[self.coords((x, y))]
        return (keep_lowest and height == 0) or height > GROUND_DEPTH

    def _get_corners(self):
        """Return the co-ordinates of the corners of the XY bounding box."""
        loc = set(self.density.keys())
        x, y = set(x[0] for x in loc), set(y[1] for y in loc)
        return (min(x), min(y)), (max(x), max(y))
    corners = property(_get_corners)

    def _tree_components(self):
        """Returns a dict where keys refer to connected components.
        NB:  use out.get(pos) to avoid KeyError for not-in-tree points."""
        #pylint:disable=bad-whitespace,missing-docstring,too-many-branches
        def neighbors(key):
            return [(key[0]+a, key[1]+b) for a, b in {
                (1,1), (1,0), (1,-1), (1,0), (0,1), (1,-1), (0,-1), (-1,-1)}]
        def expand(old_key):
            for key in neighbors(old_key):
                if com.get(key) is None:
                    # We never observed a point in that cell above GROUND_DEPTH
                    return
                if com[key] == com[old_key]:
                    continue
                elif com[key] < com[old_key]:
                    com[old_key] = com[key]
                else:
                    com[key] = com[old_key]
                    expand(key)

        non_ground = set()
        for key in self.density.keys():
            if self.canopy[key] - self.ground[key] > SLICE_DEPTH:
                cc_key = tuple(math.floor(k/JOINED_CELLS) for k in key)
                non_ground.add(cc_key)
        com = dict()
        for idx, key in enumerate(non_ground):
            com[key] = idx
        for key in com.keys():
            try:
                expand(key)
            except RuntimeError:
                # Recursion depth; finish on next pass.
                pass
        out = dict()
        translate = dict((v, k) for k, v in enumerate(set(com.values())))
        for k, v in com.items():
            out[k] = translate[v]
        return out

    def _get_components(self):
        """Accessor to avoid recalculating the connected components unless
        new points have been observed."""
        # Memorised map; recalculated if required
        if self._components_counted_at != self.count:
            self._components_counted_at = self.count
            self._components_map = self._tree_components()
        return self._components_map
    trees = property(_get_components)

    def get_component_at(self, point):
        """Get the component ID at a particular point; int>=0 or None"""
        x, y, *_ = point
        key = self.coords((x, y))
        key = tuple(math.floor(k/JOINED_CELLS) for k in key)
        return self.trees.get(key)

    def _get_len_components(self):
        """Returns the number of trees in the study area."""
        return max(set(self.trees.values())) + 1 # zero-indexed
    len_components = property(_get_len_components)

    def tree_data(self, tree_id):
        """Return (ID, X, Y, height, area, mean_red, mean_green, mean_blue,
                   point_count)"""
        #pylint:disable=too-many-locals
        assert isinstance(tree_id, int)
        if tree_id > self.len_components:
            raise ValueError('Given tree ID is out of range')
        retval = []
        keys = [k for k, v in self.trees.items() if v == tree_id]
        retval.append((max(k[0] for k in keys) + min(k[0] for k in keys)) / 2)
        retval.append((max(k[1] for k in keys) + min(k[1] for k in keys)) / 2)
        height = 0
        r, g, b = 0, 0, 0
        points = 0
        for x, y in keys:
            small_keys = ((JOINED_CELLS*x+a, JOINED_CELLS*y+b)
                          for a in range(JOINED_CELLS)
                          for b in range(JOINED_CELLS))
            for k in small_keys:
                height = max(height,
                             self.canopy.get(k, 0) - self.ground.get(k, 0))
                R, G, B = self.colours.get(k, (0, 0, 0))
                r += R
                g += G
                b += B
                points += self.density.get(k, 0)
        retval.append(height)
        retval.append(len(keys) * (CELL_SIZE * JOINED_CELLS) ** 2)
        for c in (r, g, b):
            retval.append(c/points)
        retval.append(points)
        return tuple(retval)


def remove_ground(filename, keep_lowest=False):
    """Precise ground removal, without hurting tree height (much).
    Operates by dividing the cloud into square columns col_w metres wide,
    and removes the bottom depth meters of each."""
    attr = MapObj(pointcloudfile.read(filename))
    for point in pointcloudfile.read(filename):
        if attr.point_not_ground(point, keep_lowest):
            yield point

def bin_trees(filename):
    """Takes a forest, and returns a generator of point clouds.
    Each such subcloud of the forest is a tree.
    Holds all above-ground points in memory at maximum..."""
    attr = MapObj(pointcloudfile.read(filename))
    features = [[] for _ in range(attr.len_components)]
    for point in pointcloudfile.read(filename):
        idx = attr.get_component_at(point)
        if idx is not None and attr.point_not_ground(point):
            features[idx].append(point)
    return heuristic_tree_filter(features)

def heuristic_tree_filter(cloud_list):
    """A heuristic: which clouds in this list represent a tree?
    Generator yields those that pass the heuristic.
    """
    while cloud_list:
        t = list(cloud_list.pop())
        height = [m[2] for m in t]
        if len(t) > 50 and max(height) - min(height) > SLICE_DEPTH:
            yield t

def offset_for(filename):
    """Return the (x, y, z) offset for a Pix4D .ply cloud."""
    offset = filename[:-4] + '_ply_offset.xyz'
    if not os.path.isfile(offset):
        return 0, 0, 0
    with open(offset) as f:
        return [float(n) for n in f.readline().split(' ')]

def count_trees(fname):
    """Count the number of trees in a pointcloud."""
    attr = MapObj(pointcloudfile.read(fname))
    print('There are {} trees'.format(attr.len_components))
    return attr.len_components

def stream_spatial(fname):
    """Analysis of the forest, without allowing for point export or color
    analysis."""
    # NB: colours include ground points.
    # Until that's fixed, first reduce density to minimise impact.
    lines = ['GPS?_X, GPS?_Y, X, Y, height, area, '
             'mean_red, mean_green, mean_blue, point_count,\n']
    form_str = ('{:.1f}, {:.1f}, {:.1f}, {:.1f}, {:.2f}, {:.2f}, '
                '{:.3f}, {:.3f}, {:.3f}, {:.0f},\n')
    x_, y_, _ = offset_for(fname)
    attr = MapObj(pointcloudfile.read(fname))
    print('Got points, analysing...')
    for ID in range(attr.len_components):
        x, y, *rest = attr.tree_data(ID)
        lines.append(form_str.format(x+x_, y+y_, x, y, *rest))
    with open('analysis_spatial_'+fname+'.csv', 'w') as f:
        f.writelines(lines)

def get_args():
    """Return CLI arguments to determine functions to run."""
    parser = argparse.ArgumentParser(
        description='Takes a .ply forest  point cloud; ' +
        'outputs attributes of each tree to .csv')
    parser.add_argument('file', help='name of the file to process')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('-a', '--analyse',
                      help='save .csv list of attributes for each tree',
                      action='store_true')
    mode.add_argument('-r', '--reduce_density',
                      help='discard all but lowest ground points in each cell',
                      action='store_true')
    return parser.parse_args()

if __name__ == '__main__':
    # pylint doesn't deal with argparse very well
    #pylint:disable=undefined-variable
    args = get_args()
    if not os.path.isfile(args.file):
        raise FileNotFoundError
    if not args.file.endswith('.ply'):
        raise ValueError('file must be a point cloud in .ply format')
    if args.analyse:
        stream_spatial(args.file)
    elif args.reduce_density:
        pointcloudfile.write(remove_ground(args.file, True),
                             args.file[:-4] + '_groundless.ply')
    else:
        print('A mode argument is required.  Try "--help"" for more info')

