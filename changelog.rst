.. _changelog:

Changelog
#########

Please note that ``forestutils`` is considered to be in alpha.
`Semantic versions <http://semver.org/>`_ below ``1.x.y``
indicate that anything may change without notice.

::

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

0.2.0
=====
Upgraded to handle pointclouds with any vertex attributes - no longer
just XYZRGB.  Useful for FLIR, multispectral, or even hyperspectral
data!

Consequential changes also made .ply parsing more robust
towards unexpected data or unknown formats.

More tests again, adding type analysis with ``mypy``.


0.1.2
=====
Improved testing structure and declared dependencies for testing.


0.1.1
=====
Improved documentation.


0.1.0
=====
First public release on PyPi.

Earlier versions were developed for a remote sensing project,
and internal use in the `Borevitz Lab <http://borevitzlab.anu.edu.au>`_.
