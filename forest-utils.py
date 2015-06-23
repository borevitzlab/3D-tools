#!/usr/bin/env python3
"""Tools for analysing forest point clouds."""

# TODO: fix occasional classification of canopy as ground
# TODO: refactor to split accumulation of data from analysis
# TODO: improve colour accuracy
# TODO: output UTM location
# TODO: support persistent unique identifiers for each tree

import argparse
import math
import os

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
        R, G, B = self.colours.get(idx, (0, 0, 0))
        self.colours[idx] = (R+r, G+g, B+b)
        if not self.density.get(idx):
            self.density[idx] = 1
            self.canopy[idx] = z
            self.ground[idx] = z
            return
        self.density[idx] += 1
        if self.ground[idx] > z:
            self.ground[idx] = z
        elif self.canopy[idx] < z:
            self.canopy[idx] = z

    @staticmethod
    def coords(pos, invert=False):
        """Return a tuple of integer coordinates as keys for the dict/map.
        * pos can be a full point tuple, or just (x, y)
        * use floor() to avoid imprecise float issues
        invert goes from stored to true coordinates
        """
        if invert:
            return pos[0] * CELL_SIZE, pos[1] * CELL_SIZE
        return math.floor(pos[0] / CELL_SIZE), math.floor(pos[1] / CELL_SIZE)

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
        #pylint:disable=too-many-branches
        def expand(old_key):
            """Implement depth-first search."""
            neighbors = [(old_key[0]+a, old_key[1]+b) for a, b in
                         {(1, 1), (1, 0), (1, -1), (1, 0),
                          (0, 1), (1, -1), (0, -1), (-1, -1)}]
            for key in neighbors:
                if com.get(key) is None:
                    # No above-ground point in that joined cell
                    return
                if com[key] == com[old_key]:
                    continue
                elif com[key] < com[old_key]:
                    com[old_key] = com[key]
                else:
                    com[key] = com[old_key]
                    expand(key)

        # Set up a boolean array of larger keys to search
        non_ground = set()
        key_scale_record = dict()
        for key in self.density.keys():
            if self.canopy[key] - self.ground[key] > SLICE_DEPTH:
                cc_key = tuple(math.floor(k/JOINED_CELLS) for k in key)
                if not cc_key in key_scale_record:
                    key_scale_record[cc_key] = {key}
                else:
                    key_scale_record[cc_key].add(key)
                non_ground.add(cc_key)
        # Assign a unique integer value to each large key, then search
        com = dict()
        for idx, key in enumerate(non_ground):
            com[key] = idx
        for key in com.keys():
            try:
                expand(key)
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
        return self.trees.get(self.coords(point))

    def _get_len_components(self):
        """Returns the number of trees in the study area."""
        return max(set(self.trees.values())) + 1 # zero-indexed
    len_components = property(_get_len_components)

    def tree_data(self, tree_id):
        """Yield X, Y, height, area, mean_red, mean_green, mean_blue, points"""
        assert isinstance(tree_id, int)
        if tree_id > self.len_components:
            raise ValueError('Given tree ID is out of range')
        keys = [k for k, v in self.trees.items() if v == tree_id]
        location = [self.coords(pos, True) for pos in keys]
        yield (max(k[0] for k in location) + min(k[0] for k in location)) / 2
        yield (max(k[1] for k in location) + min(k[1] for k in location)) / 2
        height, r, g, b, points = 0, 0, 0, 0, 0
        for k in keys:
            height = max(height, self.canopy[k] - self.ground[k])
            R, G, B = self.colours[k]
            r += R
            g += G
            b += B
            points += self.density[k]
        yield height
        yield len(keys) * CELL_SIZE**2
        for c in (r, g, b):
            yield c / points
        yield points


def remove_ground(filename, attr=None, keep_lowest=True):
    """Precise ground removal, without hurting tree height (much).
    Operates by dividing the cloud into square columns col_w metres wide,
    and removes the bottom depth meters of each."""
    if attr is None:
        attr = MapObj(pointcloudfile.read(filename))
    for point in pointcloudfile.read(filename):
        if attr.point_not_ground(point, keep_lowest):
            yield point

def stream_analysis(fname, attr, out):
    """Analysis of the forest, without allowing for point export or color
    analysis."""
    #pylint:disable=too-many-locals
    lines = ['latitude, longitude, UTM_X, UTM_Y, UTM_zone, height, area, '
             'mean_red, mean_green, mean_blue, point_count,\n']
    form_str = ('{}, {}, {:.1f}, {:.1f}, {:.2f}, {:.2f}, '
                '{:.3f}, {:.3f}, {:.3f}, {:.0f},\n')
    x_, y_, _ = pointcloudfile.offset_for(fname)
    for ID in range(attr.len_components):
        x, y, *rest = attr.tree_data(ID)
        utm_x, utm_y = x+x_, y+y_
        lat, lon = UTM_to_LatLon(utm_x, utm_y, UTM_ZONE)
        lines.append(form_str.format(lat, lon, utm_x, utm_y, UTM_ZONE, *rest))
    with open(out, 'w') as f:
        f.writelines(lines)

def get_args():
    """Handle setup, parsing, and validation.  Return args object."""
    parser = argparse.ArgumentParser(
        description='Takes a .ply forest  point cloud; ' +
        'outputs attributes of each tree to .csv')
    parser.add_argument('file', help='name of the file to process')
    parser.add_argument('out', nargs='?', default='.',
                        help='directory for output files (optional)')
    args = parser.parse_args()
    # TODO:  remove as file locator becomes smarter?
    if not os.path.isfile(args.file):
        raise FileNotFoundError(args.file)
    if not args.file.endswith('.ply'):
        raise ValueError('file must be a point cloud in .ply format')
    return args

def main_processing(args):
    """Logic on which functions to call, and efficient order."""
    print('Working...')
    # Set up name of output file, when input may be groundless or partial...
    groundless = os.path.join(args.out, os.path.basename(args.file))
    if not args.file.endswith('_groundless.ply'):
        groundless = os.path.join(
            args.out, os.path.basename(args.file)[:-4] + '_groundless.ply')
    groundless = groundless.replace('_part_1', '')
    # Execute main logic
    if os.path.isfile(groundless):
        attr_map = MapObj(pointcloudfile.read(groundless))
    else:
        attr_map = MapObj(pointcloudfile.read(args.file))
        pointcloudfile.write(remove_ground(args.file, attr_map), groundless)
    table = '{}_analysis.csv'.format(
        groundless[:-4].replace('_groundless', ''))
    stream_analysis(groundless, attr_map, table)
    print('Done.')

if __name__ == '__main__':
    main_processing(get_args())

