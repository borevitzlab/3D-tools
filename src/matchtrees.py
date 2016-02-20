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

import argparse
import csv

from . import utm_convert


def name_from_location(in_list, new_loc, tolerance=1):
    """Return the name of the closest location in the list to new_loc.

    Args:
        in_list: a sequence of (name, x, y) tuples
        new_loc: an (x, y) tuple

    Returns:
        The name of the closest location (if any) within tolerance, or None.
    """
    min_dist = 2 * tolerance
    name = None
    for name, x, y in in_list:
        dist_square = abs(new_loc[0] - x)**2 + abs(new_loc[1] - y)**2
        min_dist = min(min_dist, dist_square)
    return name if min_dist**0.5 <= tolerance else None


def name_loc_list_from_csv(fname):
    """Yield (name_str, UTM_x, UTM_Y) tuples from a csv file.

    The name_str is the first column of the file.
    Coordinates are taken from columns with names in (latitude, longitude,
    UTM_X, UTM_Y), case-insensitive.

    Args:
        fname (path): the filename of the .csv input data file

    Yields:
        Tuples of a name and it's UTM coordinates (name, x, y) from the file.
    """
    with open(fname) as csvfile:
        firstcol = next(csv.reader(csvfile))[0].lower()
    with open(fname) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            row = {k.strip().lower(): v for k, v in row.items()}
            if 'utm_x' in row and 'utm_y' in row:
                yield row[firstcol], float(row['utm_x']), float(row['utm_y'])
            elif 'latitude' in row and 'longitude' in row:
                x, y, _ = utm_convert.LatLon_to_UTM(
                    float(row['latitude']), float(row['longitude']))
                yield row[firstcol], x, y


def paste_names_across_analyses(old_file, new_file, tolerance=1):
    """Copy the names of trees (locations) from one csv file to another.

    The matching algorithm is very simple; get the XY location of each name
    from latitude and longitude colums (or directly from UTM, if available).

    Each entry in the new_file takes the name of the closest in the old_file
    by XY distance.  If the nearest is above tolerance_metres away, the new
    name is kept.

    NOT YET IMPLEMENTED
    """
    #pylint:disable=unused-argument
    raise NotImplementedError


def get_args():
    """Commandline controls, for standalone usage."""
    parser = argparse.ArgumentParser(description=('Unfinished module; avoid.'))
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
