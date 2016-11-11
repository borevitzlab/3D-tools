#! python3
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
"""
# pylint:disable=unsubscriptable-object

import argparse
from collections import defaultdict, namedtuple, OrderedDict
import csv
import math
import os
import statistics
from typing import MutableMapping, NamedTuple, Tuple, Set

import utm  # type: ignore

try:
    import geoply
except ImportError:
    from . import geoply


# User-defined types
XY_Coord = NamedTuple('XY_Coord', [('x', int), ('y', int)])
Coord_Labels = MutableMapping[XY_Coord, int]


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


class MapDict(defaultdict):
    """Calling this defaultdict returns the grid cell for the point"""
    def __init__(self, *args, cellsize, **kwargs):
        self._cellsize = cellsize
        super().__init__(*args, **kwargs)

    def __call__(self, pos=None):
        try:
            key = self.to_key(pos)
        except KeyError:
            if pos is None:
                return self
            raise
        return self[key]

    def to_key(self, pos):
        """Convert a vertex to the corresponding hashable key"""
        return XY_Coord(math.floor(pos['x'] / self._cellsize),
                        math.floor(pos['y'] / self._cellsize))


class ForestMap:
    """Easy-to-use data storage; basically a defaultdict collection."""

    def __init__(self, cellsize):
        self.density = MapDict(int, cellsize=cellsize)
        self.ground = MapDict(lambda: math.inf, cellsize=cellsize)
        self.canopy = MapDict(lambda: -math.inf, cellsize=cellsize)
        self.colors = MapDict(lambda: defaultdict(int), cellsize=cellsize)

    def add_vertices(self, vertices):
        for point in vertices:
            idx = self.density.to_key(point)
            for name, value in zip(point.dtype.names, point):
                if name not in 'xyz':
                    self.colors[idx][name] += value.item()
            if self.ground[idx] > point['z']:
                self.ground[idx] = point['z']
            elif self.canopy[idx] < point['z']:
                self.canopy[idx] = point['z']
            self.density[idx] += 1

    @classmethod
    def from_parts(cls, *forestmaps):
        assert len(forestmaps) >= 2
        assert all(isinstance(m, cls) for m in forestmaps)
        self, *others = forestmaps
        for other in forestmaps:
            for key in other.density:
                for name in self.colors[key]:
                    self.colors[key][name] += other.colors[key][name]
                if self.canopy[key] < other.canopy[key]:
                    self.canopy[key] == other.canopy[key]
                if self.ground[key] < other.ground[key]:
                    self.canopy[key] == other.ground[key]
                self.density[key] += other.density[key]
        return self

    @property
    def as_array(self):
        """Return the map as a 2D structured array, offset, and cellsize."""
        keys = self.density.keys()
        offset = XY_Coord(min(k.x for k in self.canopy),
                          min(k.y for k in self.canopy))
        keys = {k._replace(x=k.x-offset.x, y=k.y-offset.y) for k in keys}
        shape = max(k.x for k in keys), max(k.y for k in keys)
        raise NotImplementedError('TODO:  make 2D structured array')
        return namedtuple('Map', 'array offset cellsize')(
            None, offset, self.density._cellsize)


class MapObj:
    """Stores a maximum and minimum height map of the cloud, in GRID_SIZE
    cells.  Hides data structure and accessed through coordinates."""
    # pylint:disable=too-many-instance-attributes

    def __init__(self, input_file, *, args):
        """
        Args:
            input_file (path): the ``.ply`` file to process.  If dealing with
                Pix4D outputs, ``*_part_1.ply``.
            args (namespace): parsed command-line arguments.
            colours (bool): whether to read colours from the file.  Set to
                False for eg. LIDAR data where mean colour is not useful.
        """
        self.cellsize = args.cellsize
        self.geoply = geoply.GeoPly.read(input_file)
        self.map = ForestMap(self.geoply.vertices, self.cellsize)

        self.slicedepth = args.slicedepth
        self.grounddepth = args.grassdepth
        self.joinedcells = args.joinedcells


    def __len__(self) -> int:
        """Total observed points."""
        return sum(self.map.density.values())


    def _tree_components(self) -> Coord_Labels:
        """Returns a dict where keys refer to connected components.
        NB: Not all keys in other dicts exist in this output."""
        # Set up a boolean array of larger keys to search
        key_scale_record = {}  # type: Dict[XY_Coord, Set[XY_Coord]]
        for key in self.map.canopy:
            if self.map.canopy[key] - self.map.ground[key] > self.slicedepth:
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
        return MapDict(None, (
            (s, trees[k]) for k, v in key_scale_record.items() for s in v),
                       cellsize=self.cellsize)


    def tree_data(self, keys: Set[XY_Coord]) -> dict:
        """Return a dictionary of data about the tree in the given keys."""
        assert keys
        utmx, utmy, zone, northern = self.geoply.utm_coord
        x = utmx + self.cellsize * sum(k.x for k in keys) / len(keys)
        y = utmy + self.cellsize * sum(k.y for k in keys) / len(keys)
        lat, lon = utm.to_latlon(x, y, zone, northern=northern)
        return OrderedDict((
            ('latitude', lat),
            ('longitude', lon),
            ('UTM_easting', x),
            ('UTM_northing', y),
            ('height', max(self.map.canopy[k] - self.map.ground[k]
                           for k in keys)),
            ('area', len(keys) * self.cellsize**2),
            ('base_altitude', sum(self.map.ground[k]
                                  for k in keys) / len(keys)),
            ('point_count', sum(self.map.density[k] for k in keys)),
            *((band, statistics.mean(
                [self.map.colors[k][band] / self.map.density[k]
                 for k in keys]))
              for band in self.map.colors[next(iter(keys))].keys()),
            ))


    def all_trees(self):
        """Yield the characteristics of each tree."""
        # Group coord-keys by tree ID; pass each set of coords to tree_data()
        tree_coords = defaultdict(set)
        for coord, tree_id in self._tree_components().items():
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
        newpoints = (
            point for point in self.geoply.vertices if
            canopy and point['z'] - self.map.ground(point) > self.grounddepth
            or lowest and point['z'] == self.map.ground(point))
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
        tree_clouds = defaultdict(list)
        trees = self._tree_components()
        for point in self.geoply.vertices:
            tree = trees(point)
            if tree is not None:
                tree_clouds[tree].append(point)
        # Write them to disk
        for tree, points in tree_clouds.items():
            geoply.GeoPly.from_iterable(points, self.geoply.utm_coord).write(
                os.path.join(trees_dir, 'tree_{}.ply'.format(tree)))


    def stream_analysis(self, out: str) -> None:
        """Saves the list of trees with attributes to the file 'out'."""
        rows = list(self.all_trees())
        with open(out, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)


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
        sparse = os.path.join(args.out, os.path.basename(args.file)
                              .rstrip('.ply') + '_sparse.ply')
    if not os.path.isfile(sparse):
        attr_map = MapObj(args.file, args=args)
        print('Read {} points into {} cells, writing "{}" ...'.format(
            len(attr_map), len(attr_map.map.canopy), sparse))
        attr_map.save_sparse_cloud(sparse)
        print('Reading colours from ' + sparse)
    attr_map = MapObj(sparse, args=args)
    print('Read {} points into {} cells'.format(
        len(attr_map), len(attr_map.map.canopy)))
    print('File IO complete, starting analysis...')
    table = '{}_analysis.csv'.format(sparse[:-4].replace('_sparse', ''))
    attr_map.stream_analysis(table)
    if args.savetrees:
        print('Saving individual trees...')
        attr_map.save_individual_trees(args.savetrees)
    print('Done.')


if __name__ == '__main__':
    main()
