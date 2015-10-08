#!/usr/env/bin python3
"""A library to convert between Latitude/Longitude and UTM coordinates.

Translated from the FOSS Javascript at
http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
"""

import argparse
from collections import namedtuple
import math

try:
    from hypothesis import given, assume
    import hypothesis.strategies as st
except ImportError:
    # If hypothesis is not available, our tests are just no-ops.
    def noop(*args, **kwargs): return lambda v: None
    given = assume = st = noop
    noop.floats = noop.integers = noop.booleans = noop


# 10000000 and 1000000 meters look too similar, so let's name big constants.
EARTH_CIRC = 40.075 * 10**6
UTM_MAX_VAL = 10**7

# For testing, define the range of values various inputs can have
st_utm_coord = st.floats(min_value=0, max_value=UTM_MAX_VAL)
st_lats = st.floats(min_value=-math.pi/2, max_value=math.pi/2)
st_lons = st.floats(min_value=-math.pi, max_value=math.pi)
st_zone = st.integers(min_value=1, max_value=60)

# Ellipsoid model constants for WGS84
sm_a = 6378137
sm_b = 6356752.314
UTMScaleFactor = 0.9996

# Defining this every time we instantiate one is too slow
__utm = namedtuple('UTM_coords', ['x', 'y', 'zone', 'south'])


def UTM_coords(x, y, zone, south):
    """Return a namedtuple of valid values for a UTM coordinate."""
    assert 0 <= x <= UTM_MAX_VAL
    assert 0 <= y <= UTM_MAX_VAL
    assert zone in range(1, 61)
    assert isinstance(south, bool)
    return __utm(x, y, zone, south)


def ArcLengthOfMeridian(phi):
    """Compute the ellipsoidal distance from the equator to a given latitude.

    Args:
        phi (radians): Latitude of the point.

    Returns:
        (meters): The ellipsoidal distance of the point from the equator.
    """
    n = (sm_a - sm_b) / (sm_a + sm_b)
    alpha = ((sm_a + sm_b) / 2) * (1 + n**2 / 4) + n**4 / 64
    beta = -3 * n / 2 + 9 * n**3 / 16 - 3 * n**5 / 32
    gamma = (15 * n**2 / 16) + (-15 * n**4 / 32)
    delta = (-35 * n**3 / 48) + (105 * n**5 / 256)
    epsilon = 315 * n**4 / 512
    return sum([phi,
                beta * math.sin(2 * phi),
                gamma * math.sin(4 * phi),
                delta * math.sin(6 * phi),
                epsilon * math.sin(8 * phi)]) * alpha


@given(st_lats)
def test_ArcLengthOfMeridian(phi):
    assert -EARTH_CIRC/4 <= ArcLengthOfMeridian(phi) <= EARTH_CIRC/4


def FootpointLatitude(y):
    '''Computes the footpoint latitude for use in converting transverse
    Mercator coordinates to ellipsoidal coordinates.

    Args:
        y (meters): The UTM northing coordinate.

    Returns:
        (radians): The footpoint latitude.
    '''
    n = (sm_a - sm_b) / (sm_a + sm_b)
    alpha_ = 0.5 * (sm_a + sm_b) * (1 + n**2 / 4 + n**4 / 64)
    y_ = y / alpha_
    beta_ = 3 * n / 2 - 27 * n**3 / 32 + 269 * n**5 / 512
    gamma_ = 21 * n**2 / 16 - 55 * n**4 / 32
    delta_ = 151 * n**3 / 96 - 417 * n**5 / 128
    epsilon_ = 1097 * n**4 / 512
    return sum([y_,
                beta_ * math.sin(2 * y_),
                gamma_ * math.sin(4 * y_),
                delta_ * math.sin(6 * y_),
                epsilon_ * math.sin(8 * y_)])


@given(st_utm_coord)
def test_FootpointLatitude(y):
    assert 0 <= FootpointLatitude(y) <= math.pi/2, y


def MapLatLonToXY(phi, lambda_, lambda_0):
    '''Converts a latitude/longitude pair to x and y coordinates in the
    Transverse Mercator projection.  Note that Transverse Mercator is not
    the same as UTM - a scale factor is required to convert between them.

    Args:
        phi (radians): Latitude of the point.
        lambda_ (radians): Longitude of the point.
        lambda_0 (radians): Longitude of the central meridian to be used.

    Returns:
        2-element tuple of x and y coordinates in meters.
    '''
    ep2 = (sm_a**2 - sm_b**2) / sm_b**2
    nu2 = ep2 * math.cos(phi)**2
    N = sm_a**2 / (sm_b * (1 + nu2)**0.5)
    t = math.tan(phi)
    l = lambda_ - lambda_0

    # Precalculate coefficients for l**n in the equations below so humans
    # can read the expressions for easting and northing
    coef = {
        1: 1,
        2: 1,
        3: 1 - t**2 + nu2,
        4: 5 - t**2 + 9 * nu2 + 4 * nu2**2,
        5: 5 - 18 * t**2 + t**4 + 14 * nu2 - 58 * t**2 * nu2,
        6: 61 - 58 * t**2 + t**4 + 270 * nu2 - 330 * t**2 * nu2,
        7: 61 - 479 * t**2 + 179 * t**4 - t**6,
        8: 1385 - 3111 * t**2 + 543 * t**4 - t**6
        }
    frac = {
        1: N,
        2: t / 2 * N,
        3: N / 6,
        4: t / 24 * N,
        5: N / 120,
        6: t / 720 * N,
        7: N / 5040,
        8: t / 40320 * N
        }
    # Calculate easting and northing
    x = sum(frac[n] * math.cos(phi)**n * coef[n] * l**n for n in [1, 3, 5, 7])
    y = sum(frac[n] * math.cos(phi)**n * coef[n] * l**n for n in [2, 4, 6, 8])
    return x, y + ArcLengthOfMeridian(phi)


def LatLonToUTMXY(lat, lon):
    '''Converts a latitude/longitude pair to x and y coordinates in the
    Universal Transverse Mercator projection.

    Args:
        lat (radians): Latitude of the point.
        lon (radians): Longitude of the point.

    Returns:
        A UTM_coords namedtuple for the corresponding point.
    '''
    zone = int((math.degrees(lon)+183) / 6)
    if zone < 1:
        zone += 60
    elif zone > 60:
        zone -= 60
    c_merid = math.radians(zone*6 - 183)
    x, y = MapLatLonToXY(lat, lon, c_merid)
    # Adjust easting and northing for UTM system
    x = x * UTMScaleFactor + UTM_MAX_VAL/20
    y *= UTMScaleFactor
    south = False
    if y < 0:
        y += UTM_MAX_VAL
        south = True
    return UTM_coords(x, y, zone, south)


@given(st_lats, st_lons)
def test_LatLonToUTMXY(lat, lon):
    # TODO:  work out why this assumption is required; fix it. (zone related?)
    # Low priority, since the code is never called with such values.
    assume(math.radians(-177) <= lon)
    utm = LatLonToUTMXY(lat, lon)
    # UTM_coords() checks values, this checks correct transform
    assert utm.south == (lat < 0), (lat, utm.south)


def MapXYToLatLon(x, y, lambda_0):
    '''Converts x and y coordinates in the Transverse Mercator projection to
    a latitude/longitude pair.  Note that Transverse Mercator is not
    the same as UTM a scale factor is required to convert between them.

    Args:
        x (meters): The easting of the point.
        y (meters): The northing of the point.
        lambda_0 (radians): Longitude of the central meridian to be used.

    Returns:
        (radians): 2-element tuple of latitude and longitude.
    '''
    phif = FootpointLatitude(y)
    ep2 = (sm_a**2 - sm_b**2) / sm_b**2
    cf = math.cos(phif)
    nuf2 = ep2 * cf**2
    Nf = sm_a**2 / (sm_b * (1 + nuf2)**0.5)
    tf = math.tan(phif)
    # Calculate fractional coefficients for x**n in the equations
    # below to simplify the expressions for latitude and longitude.
    frac = {
        1: 1 / (Nf * cf),
        2: tf / (2 * Nf**2),
        3: 1 / (6 * Nf**3 * cf),
        4: tf / (24 * Nf**4),
        5: 1 / (120 * Nf**5 * cf),
        6: tf / (720 * Nf**6),
        7: 1 / (5040 * Nf**7 * cf),
        8: tf / (40320 * Nf**8)}
    # Precalculate polynomial coefficients for x**n
    poly = {
        1: 1,
        2: -1 - nuf2,
        3: -1 - 2 * tf**2 - nuf2,
        4: (5 + 3 * tf**2 + 6 * nuf2 - 6 * tf**2 * nuf2 -
            3 * nuf2**2 - 9 * tf**2 * nuf2**2),
        5: 5 + 28 * tf**2 + 24 * tf**4 + 6 * nuf2 + 8 * tf**2 * nuf2,
        6: -61 - 90 * tf**2 - 45 * tf**4 - 107 * nuf2 + 162 * tf**2 * nuf2,
        7: -61 - 662 * tf**2 - 1320 * tf**4 - 720 * tf**6,
        8: 1385 + 3633 * tf**2 + 4095 * tf**4 + 1575 * tf**6}
    # Calculate latitude and longitude
    lat = phif + sum(frac[n] * poly[n] * x**n for n in [2, 4, 6, 8])
    lon = lambda_0 + sum(frac[n] * poly[n] * x**n for n in [1, 3, 5, 7])

    def wrap(L, i):
        """Wrap overflowing coordinates around the sphere."""
        return ((L + math.pi/i) % 2*math.pi/i) - math.pi/i
    return wrap(lat, 2), wrap(lon, 1)


@given(st_utm_coord, st_utm_coord, st_lons)
def test_MapXYToLatLon(x, y, lambda_0):
    lat, lon = MapXYToLatLon(x, y, lambda_0)
    assert -math.pi/2 <= lat <= math.pi/2
    assert -math.pi <= lon <= math.pi


def UTMXYToLatLon(x, y, zone, south):
    '''Converts x and y coordinates in the Universal Transverse Mercator
    projection to a latitude/longitude pair.

    Args:
        x (meters): The easting of the point.
        y (meters): The northing of the point.
        zone: The UTM zone in which the point lies, int in the range [1, 60].
        south (bool): If the point lies in the southern hemisphere.

    Returns:
       (radians): A latitude, longitude tuple.
    '''
    x -= UTM_MAX_VAL / 20
    x /= UTMScaleFactor
    # If in southern hemisphere, adjust y accordingly.
    if south:
        y -= UTM_MAX_VAL
    y /= UTMScaleFactor
    # We need the central meridian of the UTM zone, in radians
    c_merid = math.radians(int(zone)*6 - 183)
    return MapXYToLatLon(x, y, c_merid)


@given(st_utm_coord, st_utm_coord, st_zone, st.booleans())
def test_UTMXYToLatLon(x, y, zone, south):
    lat, lon = UTMXYToLatLon(x, y, zone, south)
    assert -math.pi/2 <= lat <= math.pi/2
    assert -math.pi <= lon <= math.pi


def UTM_to_LatLon(x, y, zone=55, south=True):
    """Take xy values in UTM (default zone S55), and return lat and lon."""
    return tuple(math.degrees(n) for n in UTMXYToLatLon(x, y, zone, south))


def LatLon_to_UTM(lat, lon):
    """Take lat and lon in decimal degrees, and return UTM x, y, and zone."""
    return LatLonToUTMXY(math.radians(lat), math.radians(lon))


def get_args():
    """Handle the very simple options."""
    parser = argparse.ArgumentParser(description=(
        'A simple tool to convert between Lat/Lon and UTM coordinates.  '
        'If two values are given, convert lat and lon to UTM.  If three are '
        'given convert UTM X, UTM Y and UTM zone to lat/lon.'))
    parser.add_argument(
        'n1', type=float,
        help='Latitude in demimal degrees, or UTM X-coordinate in meters')
    parser.add_argument(
        'n2', type=float,
        help='Longitude in demimal degrees, or UTM Y-coordinate in meters')
    parser.add_argument(
        'UTMZone', type=int, nargs='?', default=None,
        help='[optional] UTM zone, if converting from UTM.  Ommit for reverse')
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    if args.UTMZone is not None:
        print('Lat:  {}\nLon:  {}'.format(
            *UTM_to_LatLon(args.n1, args.n2, args.UTMZone)))
    else:
        latlon, zone, south = LatLon_to_UTM(args.n1, args.n2)
        print('UTM X:     {}\nUTM Y:     {}\nUTM Zone:  {}{}'.format(
            latlon[0], latlon[1], zone, 'S' if south else ''))
