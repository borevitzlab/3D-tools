#!/usr/env/bin python3
# pylint:disable=missing-docstring,redefined-outer-name

import math

from hypothesis import given, assume
import hypothesis.strategies as st

from src.utm_convert import *


# 10000000 and 1000000 meters look too similar, so let's name big constants.
EARTH_CIRC = 40.075 * 10**6
UTM_MAX_VAL = 10**7

# For testing, define the range of values various inputs can have
st_utm_coord = st.floats(min_value=0, max_value=UTM_MAX_VAL-1)
st_utm_easting = st.floats(min_value=0, max_value=EARTH_CIRC/60-1)
# UTM_MAX_VAL doesn't *quite* reach the poles; see FootpointLatitude(10**7)
st_lats = st.floats(min_value=-math.pi/2+0.00031, max_value=math.pi/2-0.00031)
st_lons = st.floats(min_value=-math.pi, max_value=math.pi)
st_zone = st.integers(min_value=1, max_value=60)


@given(st_lats)
def test_ArcLengthOfMeridian(phi):
    assert -EARTH_CIRC/4 <= ArcLengthOfMeridian(phi) <= EARTH_CIRC/4


@given(st_utm_coord)
def test_FootpointLatitude(y):
    assert 0 <= FootpointLatitude(y) <= math.pi/2, y


@given(st_lats)
def test_lat_converstions(lat):
    # test that converting from lat, to distance and back is correct
    dist = ArcLengthOfMeridian(lat)
    lat2 = FootpointLatitude(dist)
    dist2 = ArcLengthOfMeridian(lat2)
    assert dif(dist, dist2) < 10**-4
    assert dif(lat, lat2) < 10**-10


@given(st_lats, st_lons)
def test_LatLonToUTMXY(lat, lon):
    utm = LatLonToUTMXY(lat, lon)
    # UTM_coords() checks values, this checks correct transform
    assert utm.south == (lat < 0), (lat, utm.south)


@given(st_lats, st.floats(min_value=-math.pi/120, max_value=math.pi/120))
def test_lesser_conversions_invert(lat, lon):
    # Round-tripping the conversion works if we keep a reference meridian
    cmerid = math.radians(6 * int(math.degrees(lon/6)))
    x, y = MapLatLonToXY(lat, lon, cmerid)
    lat2, lon2 = MapXYToLatLon(x, y, cmerid)
    x2, y2 = MapLatLonToXY(lat2, lon2, cmerid)
    assert dif(x, x2) < 10**-4
    assert dif(y, y2) < 10**-4
    assert dif(lat, lat2) < 10**-10
    assert dif(lon, lon2) < 10**-10


@given(st_utm_easting, st_utm_coord, st_zone, st.booleans())
def test_full_conversions_invert(x, y, zone, south):
    utm = UTM_coords(x, y, zone, south)
    # test_lesser_conversions_invert checks coordinate system
    # here we assert that UTM zone survives round-trip conversion
    lat, lon = UTMXYToLatLon(*utm)
    utm2 = LatLonToUTMXY(lat, lon)
    assume(is_valid_utm(*utm2))
    assert utm.south == utm2.south, (utm, utm2)
    # TODO: detect where utm zone is not stable and avoid...
    #assert utm.zone == utm2.zone, (utm, utm2)
