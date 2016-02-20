#!/usr/bin/env python3
"""A specialised io module for binary ``.ply`` files containing XYZRGB points.

Most uses of this module should go through :py:func:`read` to iterate over
points in the file, or :py:func:`write` to save an iterable of points.
Neither function accumulates much data in memory.

:py:class:`IncrementalWriter` is useful when accumulating data in memory to
write many files is impractical.  :py:func:`offset_for` and
:py:func:`read_header` provide location metadata if possible.

In all cases a "point" is tuple of (x, y, z, r, g, b).  XYZ are floats denoting
spatial coordinates.  RGB is the color, each an unsigned 8-bit integer.
While intentionally limited in scope, most data can be converted to this
format easily enough.
"""
# Pylint can freak out about mypy type notation; it's fine at runtime
# pylint:disable=unsubscriptable-object,invalid-sequence-index

from collections import namedtuple
import itertools
import struct
import os.path
from tempfile import SpooledTemporaryFile
from typing import Any, Iterator, List, Tuple

from . import utm_convert

# User-defined types:
Complicated = Tuple[int, Any, struct.Struct, List[str]]


def offset_for(filename: str) -> Tuple[float, float, float]:
    """Return the (x, y, z) UTM offset for a Pix4D or forestutils .ply file."""
    offset = filename[:-4] + '_ply_offset.xyz'
    if os.path.isfile(offset):
        with open(offset) as f:
            x, y, z = tuple(float(n) for n in f.readline().strip().split(' '))
            return x, y, z
    for com in parse_ply_header(ply_header_text(filename))[3]:
        if com.startswith('comment UTM x y zone south'):
            return float(com.split(' ')[-4]), float(com.split(' ')[-3]), 0
    return 0, 0, 0


def _check_input(fname, ending='.ply'):
    """Checks that the file exists and has the right ending"""
    if not os.path.isfile(fname):
        raise FileNotFoundError('Cannot read points from a nonexistent file')
    if not fname.endswith(ending):
        raise ValueError('Tried to read file type {}, expected {}.'.format(
            fname[-4:], ending))


def read(fname: str) -> Iterator:
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


def _read_pix4d_ply_parts(fname_list: List[str]) -> Iterator:
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
        _check_input(f)
    f = fname_list.pop(0)
    ox, oy, oz = offset_for(f)
    for p in _read_ply(f):
        yield p._replace(z=p.z+oz)
    for f in fname_list:
        dx, dy, dz = [b - a for a, b in zip([ox, oy, 0], offset_for(f))]
        for p in _read_ply(f):
            yield p._replace(x=p.x+dx)._replace(y=p.y+dy)._replace(z=p.z+dz)


def ply_header_text(filename: str) -> bytes:
    """Return the exact text of the header of the given .ply file, as bytes.

    Using bytes to allow len(header) to give index to start of data; it's
    trivial to decode in the parsing function.
    """
    _check_input(filename)
    header = b''
    with open(filename, 'rb') as f:
        while b'end_header' not in header:
            header += next(f)  # type: ignore
    return header


def parse_ply_header(header_text: bytes) -> Complicated:
    """Parse the bytes of a .ply header to useful data about the vertices.

    Deliberately discards the non-vertex data - this is a pointcloud module!
    """
    # Get file lines, do some basic validation
    lines = [l.strip() for l in header_text.decode('ascii').split('\n')]
    magic_num, data_format, *lines = lines
    if magic_num != 'ply':
        raise ValueError('Not a valid .ply file (wrong magic number).')
    if not data_format.startswith('format'):
        raise ValueError(
            'Unknown data format "{}" for .ply file.'.format(data_format))
    if 'ascii' in data_format:
        raise ValueError('ASCII format .ply files not supported at this time.')

    # Extract comments from lines
    comments = [c for c in lines if c.startswith('comment ')]
    lines = [l for l in lines if not l.startswith('comment ')]

    # Get vertex count
    element, _, vertex_count = lines.pop(0).rpartition(' ')
    if element != 'element vertex':
        raise ValueError('File must begin with vertex data!')

    # Get list of (type, name) pairs from the list of vertex properties
    lines = [l.split(' ') for l in lines]  # type: ignore
    properties = [(t, n) for _, t, n in
                  itertools.takewhile(lambda l: l[0] == 'property', lines)]

    # Get Struct format from list of property types
    PLY_TYPES = {'float': 'f', 'double': 'd', 'uchar': 'B', 'char': 'b',
                 'ushort': 'H', 'short': 'h', 'uint': 'I', 'int': 'i'}
    form_str = '>' if 'binary_big_endian' in data_format else '<'
    form_str += ''.join(PLY_TYPES[t] for t, n in properties)

    # Get Namedtuple instance from property names
    names = [n for t, n in properties]
    if not all(p in names for p in ('x', 'y', 'z')):
        raise ValueError('Pointcloud verticies must have x, y, z attributes!')
    point_tuple = namedtuple('Point', names)  # type: ignore

    # Finally, return our values
    return int(vertex_count), point_tuple, struct.Struct(form_str), comments


def _read_ply(fname: str) -> Iterator:
    """Opens the specified file, and returns a point set in the format required
    by attributes_from_cloud.  Only handles xyzrgb point clouds, but that's
    a fine subset of the format.  See http://paulbourke.net/dataformats/ply/"""
    header = ply_header_text(fname)
    vertex_count, Point, fmt, _ = parse_ply_header(header)
    with open(fname, 'rb') as f:
        f.seek(len(header))
        for _ in range(vertex_count):  # type: ignore
            yield Point(*fmt.unpack(f.read(fmt.size)))


class IncrementalWriter:
    """A streaming file writer for point clouds.

    Using the IncrementalWriter with spooled temporary files, which are
    only flushed to disk if they go above the given size, allows for
    streaming points to disk even when the header is unknown in advance.
    This allows some nice tricks, including splitting a point cloud into
    multiple files in a single pass, without memory issues.
    """
    # pylint:disable=too-few-public-methods
    # TODO:  add methods to allow use as context manager ("with _ as _: ...")

    def __init__(self, filename, utm_coords=None, buffer=2**22) -> None:
        """
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
        self.filename = filename
        self.temp_storage = SpooledTemporaryFile(max_size=buffer, mode='w+b')
        self.count = 0
        self.binary = struct.Struct('<fffBBB')
        try:
            self.utm_coords = utm_convert.UTM_coords(*utm_coords)
        except AssertionError:
            self.utm_coords = None

    def __call__(self, point) -> None:
        """Add a single point to this pointcloud, saving in binary format.

        Args:
            point: a six-element tuple of (x, y, x, red, green, blue) values,
                in that order.  (x, y, z) values are coerced to 32-bit floats.
                (red, green, blue) values are coerced to 8-bit unsigned ints.
        """
        # TODO:  URGENT - make writer method flexible enough to match reader
        self.temp_storage.write(self.binary.pack(*point[:6]))
        self.count += 1

    def __del__(self):
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
        if not os.path.isdir(os.path.dirname(self.filename)):
            os.makedirs(os.path.dirname(self.filename))
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
