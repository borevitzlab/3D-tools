Collecting Input Data
#####################

.. note:

    Parts of this document have been adapted from academic reports,
    derived from blog posts, copied from internal memos, and so on.
    If some parts seem disjointed, file an issue and forgive me if
    it was useful anyway.

Wanted: XYZRGB Pointclouds
==========================
In ``.ply`` format.  See `usage tips <usage>` for more details.


LIDAR
=====
If you have access to LIDAR hardware expensive enough to do full-colour
scans, you shouldn't need my help operating it.  Have fun, fortunate user!


Photogrammetry
==============
Photogrammetry, sometimes called Structure-from-Motion (SFM) is the technique
of reconstructing a 3D scene from multiple 2D images.

With a quadcopter to make landscape-scale scanning practical, this is a
viable and reasonably-priced option for research, forestry, and citizen
forest monitoring programmes.

The tutorials on flightriot.com_ are a good overview of the whole process.

.. _flightriot.com: http://flightriot.com/visualsfm-cmvs-post-processing-tutorial


Flight paths
------------
The key property of a flight path for photogrammetry is that every part of
the scene must be visible in multiple images, from multiple angles.
Images should have a roughly one-third overlap in each direction.

Note that while it is tempting to use a fisheye lens for a wider field of
view, the avoiding image distortion is much more important.  The less
reprojection must be done in software, the better!  Also ensure that
the shutter speed is fast enough to avoid motion blur, and that geolocation
metadata is of reasonably high precision.

- The classic option is a simple grid pattern with the camera facing directly
  downwards.  This generally works well at low altitude; coverage of
  horizontal surfaces is good but the sides of trees are often missed.

- Grid patterns with an angled camera (eg 45 degrees) work very well,
  but require many more photographs - essentially flying the pattern
  multiple times in different directions.  One is similar to a nadir
  grid; two would be better, and I suspect three would give very good
  results.

- An 'orbit' pattern works well if only a single landscape feature is of
  much interest.  Ensure a variety of angles (including altitude) and
  distances - without this the reconstruction will likely fail.


Software options
----------------
There is a wide variety of photogrammetry software out there; this section
only mentions those useful at landscape-scale.  Feel free to contribute
examples!

I do not know of any photogrammetry packages with a license that meets
the `Open Source definition <https://opensource.org/osd>`_, which is
frustrating for would-be developers and non-research users.  I encourage
anyone interested to start or support a project to fill this gap.

=============== =========================== =================================
Software        License                     Comments
=============== =========================== =================================
VisualSFM_      Free for noncommercial use, Works well with limited compute
                unmodified redistribution   power available; low cost.

CMPMVS_         Research use only (free).   Requires unusually powerful
                                            hardware (eg 4GB+ video memory).

Pix4D_          Sells expensive commercial  Very expensive, but excellent
                and costly educational      results.  Trial version works,
                licenses.                   except exporting files!
=============== =========================== =================================

.. _VisualSFM: http://ccwu.me/vsfm
.. _CMPMVS: http://ptak.felk.cvut.cz/sfmservice/websfm.pl?menu=cmpmvs
.. _Pix4D: https://pix4d.com


Suggested settings
------------------
The various software options are quite distinct, but a few tips carry over.

- Performance

  3D reconstruction is more computationally intensive than any common desktop
  application.  Look at all the options; some shave hours off the running time.

  - Ensure that GPU acceleration is enabled
  - Test a subset of your data, and run the full set overnight
  - Reserve 1GB RAM and a CPU core, if you want to run anything else
  
- Completeness of reconstruction

  Trees look very different from slightly different angles, so some settings
  tweaks are useful to ensure they are in the scene as well as the ground.
  Remember that ``forestutils`` will drop spurious ground points later.

  - Generally a point is constructed if it is recognised in 3 images.
    Turn this down to 2 images.
  - Use reduced image scale for the keypoint reconstruction stage.
  - Use full image scale (or even more) for point densification.
  - Turn desired point density to maximum.
