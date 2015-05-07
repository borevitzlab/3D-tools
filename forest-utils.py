#!/usr/bin/env python3
"""
Testing the idea of using objects that are updated by each point.

Distance units are assumed to be meters.
"""

import pointcloudfile
import math

# Default to 10cm squares for analysis; optimum depends on cloud density
CELL_SIZE = 0.1

# Classifies anything in the lowest GROUND_DEPTH in CELL_SIZE as ground
GROUND_DEPTH = 0.2

# SLICE_DEPTH is used in the connected component checks; max width of canopy
# will be >= this.
SLICE_DEPTH = 0.6

# The component analysis uses cells JOINED_CELLS * CELL_SIZE on each side
JOINED_CELLS = 4

class CloudAttributes(object):
    """Parent class for ExtentObj and ColorObj"""
    def __init__(self, initial=None):
        """Initialises the accumulator object"""
        self.count = 0
        if initial is not None:
            self.start(initial)

    def start(self, initial):
        """Updates the object with all points in the cloud argument"""
        for point in initial:
            self.update(point)

    def update(self, point):
        """The update method; override in children."""
        #pylint:disable=unused-argument
        self.count += 1

    @staticmethod
    def reformat(point):
        """Formats points as an (x, y, z, r, g, b) tuple.  This method allows
        for easy changes to the format in future; just change this interface.
        """
        return point


class ExtentObj(CloudAttributes):
    """Stores the extent of a cloud; max and min in x, y, and z dimensions."""
    def __init__(self, initial=None):
        """Initialises data"""
        self.x_max = self.x_min = self.y_max = None
        self.y_min = self.z_max = self.z_min = None
        super().__init__(initial)

    def update(self, point):
        """Updates stored data based on a point."""
        self.count += 1
        x, y, z, *_ = self.reformat(point)
        if self.x_max is None:
            self.x_max = self.x_min = x
            self.y_max = self.y_min = y
            self.z_max = self.z_min = z
            return
        if self.x_max < x:
            self.x_max = x
        elif self.x_min > x:
            self.x_min = x
        if self.y_max < y:
            self.y_max = y
        elif self.y_min > y:
            self.y_min = y
        if self.z_max < z:
            self.z_max = z
        elif self.z_min > z:
            self.z_min = z

    def _get_corners(self):
        """Return the co-ordinates of the corners of the XY bounding box."""
        return ((self.x_min, self.y_min), (self.x_max, self.y_max))
    corners = property(_get_corners)

    def _get_height(self):
        """The height of the cloud"""
        return self.z_max - self.z_min
    height = property(_get_height)

    def _get_centre(self):
        """The middle of the cloud; tuple of (x, y)"""
        return (self.x_max + self.x_min)/2, (self.y_max + self.y_min)/2
    centre = property(_get_centre)

    def __str__(self):
        """String method for printing, etc."""
        x = 'X extent: {:.1f} to {:.1f}'.format(self.x_min, self.x_max)
        y = 'Y extent: {:.1f} to {:.1f}'.format(self.y_min, self.y_max)
        z = 'Z extent: {:.1f} to {:.1f}'.format(self.z_min, self.z_max)
        return '{}\n{}\n{}\n'.format(x, y, z)


class ColorsObj(CloudAttributes):
    """Stores a histogram of encountered colors"""
    def __init__(self, initial=None):
        """Initialises data"""
        self.red = [0]*256
        self.green = [0]*256
        self.blue = [0]*256
        super().__init__(initial)

    def update(self, point):
        """Updates stored data based on a point.
        The point is a tuple of (x, y, z, r, g, b)
        """
        self.count += 1
        _, _, _, r, g, b = self.reformat(point)
        self.red[r] += 1
        self.green[g] += 1
        self.blue[b] += 1

    def _get_sum_rgb(self):
        """Returns an RGB tuple of the sum of observed values."""
        red = green = blue = 0
        for idx, r in enumerate(self.red):
            red += idx*r
        for idx, g in enumerate(self.green):
            green += idx*g
        for idx, b in enumerate(self.blue):
            blue += idx*b
        return (red, green, blue)
    rgb_sum = property(_get_sum_rgb)

    def _get_mean_rgb(self):
        """Returns an RGB tuple of the mean observed values."""
        return tuple(int(c/self.count) for c in self.rgb_sum)
    rgb = property(_get_mean_rgb)

    def _get_GCC(self):
        """Returns a float value for Green Chromatic Co-ordinate."""
        return self.rgb_sum[1] / sum(self.rgb_sum)
    GCC = property(_get_GCC)


class MapObj(CloudAttributes):
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
        self._components_counted_at = None
        self._components_map = None
        super().__init__(initial)

    def update(self, point):
        """Expand, correct, or maintain map with a new observed point."""
        self.count += 1
        x, y, z, *_ = self.reformat(point)
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

    @staticmethod
    def coords(pos):
        """Return a tuple of integer coordinates as keys for the dict/map.
        * pos can be a full point tuple, or just (x, y)
        * use floor() to avoid nasty float issues
        * integer labels so we can do index-based stuff like finding neighbors
        """
        x, y, *_ = pos
        return math.floor(x / CELL_SIZE), math.floor(y / CELL_SIZE)

    def point_not_ground(self, point):
        """Returns boolean wether the point is not classified as ground - ie
        True if within GROUND_DEPTH of the lowest point in the cell (30cm)."""
        x, y, z, *_ = self.reformat(point)
        idx = self.coords((x, y))
        return z - self.ground[idx] > GROUND_DEPTH

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
        x, y, *_ = self.reformat(point)
        key = self.coords((x, y))
        key = tuple(math.floor(k/JOINED_CELLS) for k in key)
        return self.trees.get(key)

    def _get_len_components(self):
        """Returns the number of trees in the study area."""
        return max(set(self.trees.values()))
    len_components = property(_get_len_components)



def remove_ground(filename):
    """Precise ground removal, without hurting tree height (much).
    Operates by dividing the cloud into square columns col_w metres wide,
    and removes the bottom depth meters of each."""
    attr = MapObj(pointcloudfile.read(filename))
    for point in pointcloudfile.read(filename):
        if attr.point_not_ground(point):
            yield point

def bin_trees(filename):
    """Takes a forest, and returns a generator of point clouds.
    Each such subcloud of the forest is a tree.
    Holds all above-ground points in memory at maximum..."""
    attr = MapObj(pointcloudfile.read(filename))
    print('Finished first pass')
    features = [[] for _ in range(attr.len_components + 1)]
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
        if len(t) > 50 and max(m[2] for m in t) - min(m[2] for m in t) > 0.5:
            yield t

def test_demo():
    """A test case to demonstrate the power of the module"""
    fname = 'M00002_ascii_forest_crop.ply'
    with open('analysis_'+fname+'.csv', 'w') as f:
        f.write('X, Y, height, GCC,\n')
        for i, tree in enumerate(bin_trees(fname)):
            E, C = ExtentObj(), ColorsObj()
            for p in tree:
                E.update(p)
                C.update(p)
            f.write('{:.1f}, {:.1f}, {:.2f}, {:.4f},\n'.format(
                E.centre[0], E.centre[1], E.height, C.GCC))
            pointcloudfile.write(tree, 'groundless/tree_'+str(i)+'.ply')

if __name__ == '__main__':
    test_demo()
    print('Done!')

