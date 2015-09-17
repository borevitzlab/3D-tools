#!/usr/bin/env python3
"""Tools for adjusting trees from forest-utils.py based on other datasets.

The required input to this module is a .csv file associating the name and
location of each tree in the dataset.  A prior analysis run from forest-utils
works well for this purpose if the names have been manually added at some
point, as other columns are ignored.

Data which can be corrected includes later analysis runs of forest-utils,
and the names of individual tree files output.

It may be worth moving to an iterative closest point algorithm at some point,
and matching the whole datasets at once.  More flexible input, but much more
complex code.
"""

# TODO:  enforce same UTM zone for all inputs in some elegant way.

import argparse
import csv
import os

import pointcloudfile
from utm_convert import LatLon_to_UTM


def name_from_location(in_list, new_loc, tolerance=1):
    """Return the name of the closest location in the list to new_loc,
    if the seperation is within tolerance.  Otherwise return None.

    in_list:    a sequence of (name, x, y) tuples
    new_loc:    an (x, y) tuple
    """
    return None


def name_loc_list_from_csv(fname):
    """Return a list of (name_str, UTM_x, UTM_Y) tuples from a csv file.

    The name_str is the first column of the file.
    Coordinates are taken from columns with names in (latitude, longitude,
    UTM_X, UTM_Y), case-insensitive.
    """
    raise NotImplementedError


def paste_names_across_analyses(old_file, new_file, tolerance_metres=1):
    """Copy the names of trees (locations) from one csv file to another.

    The matching algorithm is very simple; get the XY location of each name
    from latitude and longitude colums (or directly from UTM, if available).

    Each entry in the new_file takes the name of the closest in the old_file
    by XY distance.  If the nearest is above tolerance_metres away, the new
    name is kept.
    """
    raise NotImplementedError


def rename_tree_clouds(from_file, tree_dir):
    """Name files containing individual trees based on an input file.
    """
    for f in os.listdir(tree_dir):
        pass
    raise NotImplementedError


def get_args():
    """Commandline controls, for standalone usage."""
    parser = argparse.ArgumentParser(description=('Unfinished module; avoid.'))
    parser.add_argument('dummyarg', help='no-op command')
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    raise NotImplementedError
