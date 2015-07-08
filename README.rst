========
3D Tools
========
A set of tools for dealing with point clouds representing a forest. Basic
object recognition and extraction, calculation of spatial and spectral
attributes, automation of data processing, etc.


forest-utils.py
===============
Analyses a point cloud representing a forest; outputs an attribute table
with location, size, colour, and point coverage.  A "sparse" point cloud is
also created, which preserves all trees but reduces ground point density to
only the lowest point in each (10cm) grid cell.

Extensive use of mutable coordinate-property mappings and streamed input
ensure that even files to large to load in memory can be processed.  In extreme
cases, the resolution can be decreased to trade accuracy for memory.

Here's a map made from example output: 
https://www.google.com/maps/d/viewer?mid=z1pH7HaTWL9Q.kzQflQGYVRIU

pointcloudfile.py
=================
A simple compatibility layer for dealing with point cloud files; this provides
two general methods.

In all cases a "point" is tuple of (x, y, z, r, g, b).  XYZ are floats denoting
spatial coordinates.  RGB is the color, each an int between 0 and 255 inclusive.

``read(filename)`` returns a generator of points.  ``*.ply`` files are currently
supported, but other formats are planned.  ``.ply`` files exported from Pix4D
have some special handling to combine multi-file point clouds, and preserve
correct georeferencing.

``write(cloud, filename)`` writes the supplied points to ``filename`` in
binary ``.ply`` format.  ``cloud`` must be an iterable of points.

utm_convert.py
==============
A library to convert co-ordinates between the Universial Tranverse Mercator
projection, and WGSM84 (GPS-projection) latitude and longitude in decimal
degrees - in either direction.

The UTM system is nice in that it uses eastings, northings, and altitude in
metres - which is really good for local use.  Pix4D exports geolocation in
this format.  Lat/Lon degrees are a common standard for georeferencing, and
often an easier interface for other software.

This module is a translation of `Chuck Taylor's Javascript version 
<http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html>`_ to idiomatic
Python.  Floating-point imprecision is negligible; a few microns per cycle.
