3D Forest Tools
###############

|license_badge| |travis_badge| |docs_badge|

.. |license_badge| image:: https://img.shields.io/badge/license-GPL3+-blue.svg
   :target: https://en.wikipedia.org/wiki/GNU_General_Public_License

.. |travis_badge| image:: https://travis-ci.org/borevitzlab/3D-tools.svg?branch=master
   :target: https://travis-ci.org/borevitzlab/3D-tools

.. |docs_badge| image:: https://readthedocs.org/projects/3d-tools/badge
   :target: https://3d-tools.readthedocs.org

.. END_BADGES_TAG

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
