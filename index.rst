Welcome to 3D-tools's documentation!
====================================
``3D-tools`` is a set of tools for analysing pointclouds of a forest.

``forestutils.py`` is the main interface, with the other files
providing functionality that may be indepedently useful.

Features:

- Detect, map, and extract trees from the full pointcloud
- Calculate location, height, canopy area, and colour of each tree
- Losslessly reduce pointcloud size by discarding ground points

It is written in pure Python (3.4+), available under the GPL3 license,
and can analyse multi-gigabyte datasets in surprisingly little memory.

Useful links:

- `documentation <https://3d-tools.rtfd.org>`_
- `processed pointclouds <https://pointclouds.phenocam.org.au>`_
- `maps (from output) <https://www.google.com/maps/d/viewer?mid=zH6cETTFeU6A.kawPmizn8suk>`_
- `source code <https://github.com/borevitzlab/3D-tools>`_

Contact: *zac.hatfield.dodds at gmail.com*

Contents
========

.. toctree::
   :maxdepth: 2
   :glob:

   docs/*

* :ref:`genindex`
