#!/usr/bin/env python3
"""A specialised io module for .ply files containing XYZRGB points.

This is useful for all the example forest data I've seen so far, and is
easy enough to extend later for other attributes.

Returns generator objects, which yield points, due to out-of-memory issues.
"""

import struct
import os.path

def read(fname):
    """Passes the file to a read function for that format."""
    if not isinstance(fname, str):
        raise TypeError('The filename to read must be a string!')
    if not os.path.isfile(fname):
        raise FileNotFoundError('Cannot read points from a nonexistent file')
    if fname.endswith('.ply'):
        return _read_ply(fname)
    raise ValueError('Tried to read unsupported file type: ' + fname[-4:])

def _read_ply(fname):
    """Opens the specified file, and returns a point set in the format required
    by attributes_from_cloud.  Currently a special case."""
    #pylint:disable=too-many-branches
    f = open(fname, 'rb')
    raw_header = []
    while True:
        raw_header.append(next(f))
        if raw_header[-1] == b'end_header\n':
            break
    header = [line.strip().split(b' ') for line in raw_header]
    header = [(line[2], line[1]) for line in header if line[0] == b'property']
    if b'format ascii 1.0\n' in raw_header:
        head = {}
        for i, line in enumerate(header):
            head[line[0]] = (i, b'float' in line[1] or line[1] == b'double')
        for line in f:
            point = line.strip().split(b' ')
            out = tuple([float(point[head[s][0]]) if head[s][1]
                         else int(point[head[s][0]])
                         for s in (b'x', b'y', b'z',
                                   b'red', b'green', b'blue')])
            yield out
    else:
        head = [h[0] for h in header]
        idx = {'x':head.index(b'x'), 'y':head.index(b'y'), 'z':head.index(b'z'),
               'r':head.index(b'red'), 'g':head.index(b'green'),
               'b':head.index(b'blue')}
        form_str = '<'
        if b'format binary_big_endian 1.0\n' in raw_header:
            form_str = '>'
        for _, type_ in header:
            if type_ == b'float':
                form_str += 'f'
            elif type_ == b'uchar':
                form_str += 'B'
            elif type_ == b'char':
                form_str += 'b'
            elif type_ == b'short':
                form_str += 'h'
            elif type_ == b'ushort':
                form_str += 'H'
            elif type_ == b'int':
                form_str += 'i'
            elif type_ == b'uint':
                form_str += 'I'
            elif type_ == b'double':
                form_str += 'd'
        point = struct.Struct(form_str)
        raw = f.read(point.size)
        while raw:
            p_tup = point.unpack(raw)
            yield (p_tup[idx['x']], p_tup[idx['y']], p_tup[idx['z']],
                   p_tup[idx['r']], p_tup[idx['g']], p_tup[idx['b']])
            raw = f.read(point.size)
    f.close()

def write(cloud, fname, count=None):
    """Write the given cloud to disk."""
    if count is None:
        # If count is not given, load all to find out and hope for the best
        cloud = list(cloud)
        count = len(cloud)
    header = """ply
format binary_little_endian 1.0
comment written by read_tree.py
element vertex {}
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
element face 0
property list uchar int vertex_indices
end_header
""".format(count)
    binary = struct.Struct('<fffBBB')
    with open(fname, 'wb') as f:
        f.write(header.encode('ascii'))
        for point in cloud:
            #pylint:disable=star-args
            f.write(binary.pack(*point))

