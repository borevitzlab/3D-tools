========
3D Tools
========
A set of tools for dealing with point clouds representing a forest. Basic
object recognition and extraction, calculation of spatial and spectral
attributes, automation of data processing, etc.


forest-utils.py
===============
Analyses a point cloud representing a forest; outputs an attribute table
with xy coordinates, height, colour indices, and (in future) various other
data.

A key design principle is robustness over large files; this is achieved by
making minimal assumptions and processing a stream instead of loading too much
into memory where various accumulators will do.  Both are a work in progress,
but reliable on the ~33M point datasets tested so far.

Performance is desired only to a degree that does not decrease readability;
processing images to points takes several hours and this tool only a few
minutes.


pointcloudfile.py
=================
A simple compatibility layer for dealign with point cloud files; this provides
two general methods.

In all cases a "point" is tuple of (x, y, z, r, g, b).  XYZ are floats denoting
spatial coordinates.  RGB is the color, each an int between 0 and 255 inclusive.

``read(filename)`` returns a generator of points.  ``*.ply`` files are currently
supported, but other formats are planned.

``write(cloud, filename, count=None)`` writes the supplied points to
``filename`` in binary ``.ply`` format.  ``cloud`` must be an iterable of
points, and if ``count`` (the number of points) is not supplied will be
converted to a list as the format requires point count in the header.
