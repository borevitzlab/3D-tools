#!/usr/bin/env python3
"""A specialised io module for .ply files containing XYZRGB points.

This is useful for all the example forest data I've seen so far, and is
easy enough to extend later for other attributes.

Returns generator objects, which yield points, due to out-of-memory issues.
"""

import struct
import os.path
from tempfile import SpooledTemporaryFile

from utm_convert import UTM_coords


def _offset_for(filename):
    """Return the (x, y, z) offset for a Pix4D .ply cloud."""
    offset = filename[:-4] + '_ply_offset.xyz'
    try:
        with open(offset) as f:
            x, y, z = tuple(float(n) for n in f.readline().strip().split(' '))
            return x, y, z
    except FileNotFoundError:
        return 0, 0, 0


def _check_input(fname, ending=''):
    """Checks that the file exists and has the right ending"""
    if not os.path.isfile(fname):
        raise FileNotFoundError('Cannot read points from a nonexistent file')
    if not fname.endswith(ending):
        raise ValueError('Tried to read file type {}, expected {}.'.format(
            fname[-4:], ending))


def read(fname):
    """Passes the file to a read function for that format."""
    if fname.endswith('_point_cloud_part_1.ply'):
        parts, p = [fname], 1
        stub = fname.replace('_point_cloud_part_1.ply', '')
        while True:
            p += 1
            part = stub + '_point_cloud_part_{}.ply'.format(p)
            if os.path.isfile(part):
                parts.append(part)
            else:
                return _read_pix4d_ply_parts(parts)
    return _read_ply(fname)


def _read_pix4d_ply_parts(fname_list):
    """Yield points from a list of Pix4D ply files as if they were one file.

    Pix4D usually exports point clouds in parts, with an xyz offset for the
    origin.  This means that we can yield the points from each, correcting
    for the offset in the origin coordinate in each.

    We can further move the altitude information into the points without loss
    of precision (to any significant degree).  However UTM XY coordinates
    can't be added; we don't know the UTM zone and loss of precision may
    be noticible if we did.
    """
    for f in fname_list:
        _check_input(f, '.ply')
    first = fname_list.pop(0)
    ox, oy, oz = _offset_for(first)
    for point in _read_ply(first):
        x, y, z, r, g, b = point
        yield x, y, z+oz, r, g, b
    for f in fname_list:
        dx, dy, dz = [b - a for a, b in zip([ox, oy, 0], _offset_for(f))]
        for x, y, z, r, g, b in _read_ply(f):
            yield x+dx, y+dy, z+dz, r, g, b


def process_header(file_handle):
    """Return key information from a list of bytestrings (the raw header);
    a Struct format string (empty == ascii mode), and the index order tuple."""
    head = []
    while True:
        head.append(next(file_handle).decode())
        if head[-1].strip() == 'end_header':
            break

    head = [line.strip().split(' ') for line in head]
    is_big_endian = ['format', 'binary_big_endian', '1.0'] in head
    while not head[-1] == ['end_header']:
        head.pop()
        if not head:
            raise ValueError('Invalid header for .ply file.')

    offset = None
    for line in head:
        if line[:6] == ['comment', 'UTM', 'x', 'y', 'zone', 'south']:
            x, y, zone, S = line[-4:]
            offset = UTM_coords(float(x), float(y), int(zone), bool(int(S)))

    typeorder = [line[1] for line in head if line[0] == 'property']
    form_str = '>' if is_big_endian else '<'
    ply_types = {'float': 'f', 'double': 'd',
                 'uchar': 'B', 'char': 'b',
                 'ushort': 'H', 'short': 'h',
                 'uint': 'I', 'int': 'i'}
    form_str += ''.join(ply_types[t] for t in typeorder)

    head = [line[2] for line in head if line[0] == 'property']
    cols = {'r': 'red', 'g': 'green', 'b': 'blue'}
    if any(l.startswith('diffuse_') for l in head):
        cols = {k: 'diffuse_'+v for k, v in cols.items()}
    ind = {c: head.index(c) for c in 'xyz'}
    ind.update({c: head.index(cols[c]) for c in 'rgb'})
    ind = tuple(ind[c] for c in 'xyzrgb')

    return struct.Struct(form_str), ind, offset


def _read_ply(fname):
    """Opens the specified file, and returns a point set in the format required
    by attributes_from_cloud.  Only handles xyzrgb point clouds, but that's
    a fine subset of the format.  See http://paulbourke.net/dataformats/ply/"""
    _check_input(fname, '.ply')
    with open(fname, 'rb') as f:
        point, ind, _ = process_header(f)
        raw = f.read(point.size)
        if ind == tuple(range(6)):  # Faster in the most common case
            while raw:
                yield point.unpack(raw)
                raw = f.read(point.size)
        else:
            while raw:
                p_tup = point.unpack(raw)
                yield tuple(p_tup[i] for i in ind)
                raw = f.read(point.size)


class IncrementalWriter(object):
    """A streaming file writer for point clouds.

    Using the IncrementalWriter with spooled temporary files, which are
    only flushed to disk if they go above the given size, allows for
    streaming points to disk even when the header is unknown in advance.
    This allows some nice tricks, including splitting a point cloud into
    multiple files in a single pass, without memory issues.

    Args:
        filename: final place to save the file on disk.  Parent directory
            must exist.  This file will not be created until the
            :py:meth:`save_to_disk` method is called.
        utm_coords: The (x, y, z, zone) offset added to find the
            UTM coordinates of each point.
        buffer (int): The number of bytes to hold in RAM before flushing
            the temporary file to disk.  Default 1MB, which holds ~8300
            points - enough for most objects but still practical to hold
            thousands in memory.  Set a smaller buffer for large forests.
    """

    def __init__(self, filename, utm_coords=None, buffer=2**22):
        """Set up the object."""
        if not os.path.isdir(os.path.dirname(filename)):
            raise FileNotFoundError('Parent directory of the given filename '
                                    'must exist.')
        self.filename = filename
        self.temp_storage = SpooledTemporaryFile(max_size=buffer, mode='w+b')
        self.count = 0
        self.binary = struct.Struct('<fffBBB')
        try:
            self.utm_coords = UTM_coords(*utm_coords)
        except AssertionError:
            self.utm_coords = None

    def __call__(self, point):
        """Add a single point to this pointcloud, saving in binary format.

        Args:
            point: a six-element tuple of (x, y, x, red, green, blue) values,
                in that order.  (x, y, z) values are coerced to 32-bit floats.
                (red, green, blue) values are coerced to 8-bit unsigned ints.
        """
        self.temp_storage.write(self.binary.pack(*point))
        self.count += 1

    def save_to_disk(self):
        """Flush data to disk and clean up."""
        head = ['ply',
                'format binary_little_endian 1.0',
                'element vertex {}'.format(self.count),
                'property float x',
                'property float y',
                'property float z',
                'property uchar red',
                'property uchar green',
                'property uchar blue',
                'end_header']
        if self.utm_coords is not None:
            x, y, zone, south = self.utm_coords
            head.insert(-2, 'comment UTM x y zone south {:.2f} {:.2f} {:d} {}'
                        .format(x, y, zone, '1' if south else '0'))
        with open(self.filename, 'wb') as f:
            f.write(('\n'.join(head) + '\n').encode('ascii'))
            self.temp_storage.seek(0)
            chunk = self.temp_storage.read(8192)
            while chunk:
                f.write(chunk)
                chunk = self.temp_storage.read(8192)
        self.temp_storage.close()


def write(cloud, fname, utm_coords=None):
    """Write the given cloud to disk."""
    writer = IncrementalWriter(fname, utm_coords)
    for p in cloud:
        writer(p)
    writer.save_to_disk()
