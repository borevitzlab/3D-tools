#!/usr/env/bin python3
"""A library to convert between Latitude/Longitude and UTM coordinates.

Translated from the FOSS Javascript at
http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
"""

import math

# Ellipsoid model constants for WGS84
sm_a = 6378137
sm_b = 6356752.314
UTMScaleFactor = 0.9996

def ArcLengthOfMeridian(phi):
    '''Computes the ellipsoidal distance from the equator to a point at a
    given latitude.
    Inputs:
        phi - Latitude of the point, in radians.
    Returns:
        The ellipsoidal distance of the point from the equator, in meters.
    '''
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

def UTMCentralMeridian(zone):
    '''Determines the central meridian for the given UTM zone.'''
    return math.radians(int(zone)*6 - 183)

def FootpointLatitude(y):
    '''Computes the footpoint latitude for use in converting transverse
    Mercator coordinates to ellipsoidal coordinates.
    Inputs:
        y - The UTM northing coordinate, in meters.
    Returns:
        The footpoint latitude, in radians.
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

def MapLatLonToXY(phi, lambda_, lambda_0):
    '''Converts a latitude/longitude pair to x and y coordinates in the
    Transverse Mercator projection.  Note that Transverse Mercator is not
    the same as UTM - a scale factor is required to convert between them.
    Inputs:
        phi - Latitude of the point, in radians.
        lambda_ - Longitude of the point, in radians.
        lambda_0 - Longitude of the central meridian to be used, in radians.
    Outputs:
        xy - A 2-element array containing the x and y coordinates
             of the computed point.
    '''
    ep2 = (sm_a**2 - sm_b**2) / sm_b**2
    nu2 = ep2 * math.cos(phi)**2
    N = sm_a**2 / (sm_b * (1 + nu2)**0.5)
    t = math.tan(phi)
    l = lambda_ - lambda_0

    #Precalculate coefficients for l**n in the equations below so humans
    #can read the expressions for easting and northing
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
    #Calculate easting and northing
    x = sum(frac[n] * math.cos(phi)**n * coef[n] * l**n for n in [1, 3, 5, 7])
    y = sum(frac[n] * math.cos(phi)**n * coef[n] * l**n for n in [2, 4, 6, 8])
    return x, y + ArcLengthOfMeridian(phi)

def LatLonToUTMXY(lat, lon):
    '''Converts a latitude/longitude pair to x and y coordinates in the
    Universal Transverse Mercator projection.
    Inputs:
        lat - Latitude of the point, in radians.
        lon - Longitude of the point, in radians.
    Outputs:
        xy - A 2-element array where the UTM x and y values will be stored.
        The UTM zone used for calculating the values of x and y.
    '''
    zone = int((math.degrees(lon)+183) / 6)
    x, y = MapLatLonToXY(lat, lon, UTMCentralMeridian(zone))
    # Adjust easting and northing for UTM system
    x = x * UTMScaleFactor + 500000
    y *= UTMScaleFactor
    if y < 0:
        y += 10000000
    return (x, y), zone

def MapXYToLatLon(x, y, lambda_0):
    '''Converts x and y coordinates in the Transverse Mercator projection to
    a latitude/longitude pair.  Note that Transverse Mercator is not
    the same as UTM a scale factor is required to convert between them.
    Inputs:
        x - The easting of the point, in meters.
        y - The northing of the point, in meters.
        lambda_0 - Longitude of the central meridian to be used, in radians.
    Outputs:
        a tuple of latitude and longitude
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
    #Precalculate polynomial coefficients for x**n
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
    #Calculate latitude and longitude
    return (phif + sum(frac[n] * poly[n] * x**n for n in [2, 4, 6, 8]),
            lambda_0 + sum(frac[n] * poly[n] * x**n for n in [1, 3, 5, 7]))

def UTMXYToLatLon(x, y, zone, south):
    '''Converts x and y coordinates in the Universal Transverse Mercator
    projection to a latitude/longitude pair.
    Inputs:
        x - The easting of the point, in meters.
        y - The northing of the point, in meters.
        zone - The UTM zone in which the point lies.
        south - If the point is in the southern hemisphere
    Outputs:
       a tuple of latitude and longitude of the point, in radians
    '''
    x -= 500000
    x /= UTMScaleFactor
    #If in southern hemisphere, adjust y accordingly.
    if south:
        y -= 10000000
    y /= UTMScaleFactor
    # We need the central meridian of the UTM zone, in radians
    c_merid = UTMCentralMeridian(zone)
    return MapXYToLatLon(x, y, c_merid)

def UTM_to_LatLon(x, y, zone=55, south=True):
    """Take xy values in UTM (default zone S55), and return lat and lon."""
    return tuple([math.degrees(n) for n in UTMXYToLatLon(x, y, zone, south)])

def LatLon_to_UTM(lat, lon):
    """Take lat and lon in decimal degrees, and return UTM x, y, and zone."""
    return LatLonToUTMXY(math.radians(lat), math.radians(lon))

if __name__ == '__main__':
    a, b = (float(n) for n in input('What XY coords?  ').split())
    print(UTM_to_LatLon(a, b))

