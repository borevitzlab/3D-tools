#!/usr/bin/env python3
"""A specialised io module for .ply files containing XYZRGB points.

This is useful for all the example forest data I've seen so far, and is
easy enough to extend later for other attributes.

Returns generator objects, which yield points, due to out-of-memory issues.
"""

import struct
import os.path
from tempfile import SpooledTemporaryFile
from collections import namedtuple

UTM_offset = namedtuple('UTM_offset', ['x', 'y', 'z', 'zone'])


def UTM_offset_for(filename, zone=55):
    """Return the (x, y, z) offset for a Pix4D .ply cloud."""
    with open(filename, 'rb') as f:
        raw_header = []
        while True:
            raw_header.append(next(f))
            if raw_header[-1].startswith(b'end_header'):
                # ply files can have basically any line ending...
                break
        offset = _process_header(raw_header)[-1]
    if isinstance(offset, UTM_offset) and offset is not None:
        return offset
    offset = filename[:-4] + '_ply_offset.xyz'
    if not os.path.isfile(offset):
        return UTM_offset(0, 0, 0, zone)
    with open(offset) as f:
        x, y, z = (float(n) for n in f.readline().strip().split(' '))
        return UTM_offset(x, y, z, zone)


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
    """Read from a collection of ply files as if they were one file.
    Pix4D usually exports point clouds in parts, with a UTM offset for the
    origin.  We can thus read from all (though only write to one file."""
    for f in fname_list:
        _check_input(f, '.ply')
    first = fname_list.pop(0)
    for point in _read_ply(first):
        yield point
    for f in fname_list:
        dx, dy, dz, _ = [b - a for a, b in
                         zip(UTM_offset_for(first), UTM_offset_for(f))]
        for x, y, z, r, g, b in _read_ply(f):
            yield x+dx, y+dy, z+dz, r, g, b


def _process_header(head):
    """Return key information from a list of bytestrings (the raw header);
    a Struct format string (empty == ascii mode), and the index order tuple."""
    head = [line.strip().split(b' ') for line in head]
    is_big_endian = ['format', 'binary_big_endian', '1.0'] in head
    while not head[-1] == ['end_header']:
        head.pop()
        if not head:
            raise ValueError('Invalid header for .ply file.')

    offset = None
    for line in head:
        if line[:4] == ['comment', 'UTM', 'easting', 'northing']:
            x, y, z, zone = line[-4:]
            offset = UTM_offset(float(x), float(y), float(z), int(zone))

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
        raw_header = []
        while True:
            raw_header.append(next(f).decode('ascii'))
            if raw_header[-1].strip() == 'end_header':
                break
        point, ind, _ = _process_header(raw_header)
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
        self.utm_coords = utm_coords

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
        if self.utm_coords is not None and len(self.utm_coords) == 4:
            head.insert(-2, 'comment UTM_easting UTM_northing altitude zone ' +
                        '{:.2f} {:.2f} {:.2f} {:d}'.format(*self.utm_coords))
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
