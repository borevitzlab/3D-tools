#!/usr/env/bin python3
'''A library to convert UTM values to lat/lon values.
Translated from the FOSS Javascript at
http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
'''

import math

# Ellipsoid model constants for WGS84
sm_a = 6378137
sm_b = 6356752.314
UTMScaleFactor = 0.9996

def DegToRad(deg):
    '''Converts degrees to radians.'''
    return deg / 180 * math.pi

def RadToDeg(rad):
    '''Converts  radians to degrees.'''
    return rad / math.pi * 180

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
    return DegToRad(int(zone)*6 - 183)

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
    gamma_ = 1 * n**2 / 16 - 55 * n**4 / 32
    delta_ = 151 * n**3 / 96 - 417 * n**5 / 128
    epsilon_ = 1097 * n**4 / 512
    return sum([y_,
                beta_ * math.sin(2 * y_),
                gamma_ * math.sin(4 * y_),
                delta_ * math.sin(6 * y_),
                epsilon_ * math.sin(8 * y_)])

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
    return (phif + sum((frac[n] * poly[n] * x**n for n in [2, 4, 6, 8])),
            lambda_0 + sum((frac[n] * poly[n] * x**n for n in [1, 3, 5, 7])))

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
    c_merid = DegToRad(UTMCentralMeridian(zone))
    return MapXYToLatLon(x, y, c_merid)

def UTM_to_LatLon(x, y, zone=55, south=True):
    """Take xy values in UTM (default zone S55), and return lat and lon."""
    return tuple([RadToDeg(n) for n in UTMXYToLatLon(x, y, zone, south)])

if __name__ == '__main__':
    a, b = (float(n) for n in input('What XY coords?  ').split())
    print(UTM_to_LatLon(a, b))

