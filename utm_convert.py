#!/usr/env/bin python3
''' Ellipsoid model constants (actual values here are for WGS84) '''
import math

# Translated from http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
# Thanks!

#pylint:disable=bad-whitespace

sm_a = 6378137
sm_b = 6356752.314
sm_EccSquared = 6.69437999013e-03

UTMScaleFactor = 0.9996

def DegToRad(deg):
    '''Converts degrees to radians.'''
    return deg / 180 * math.pi

def RadToDeg (rad):
    '''Converts  radians to degrees.'''
    return rad / math.pi * 180

def ArcLengthOfMeridian (phi):
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
    gamma = (15 * n**2 / 16) + (-15 * math.pow(n, 4) / 32)
    delta = (-35 * n**3 / 48) + (105 * math.pow(n, 5) / 256)
    epsilon = 315 * n**4 / 512
    return alpha * sum([phi,
                        beta * math.sin(2 * phi),
                        gamma * math.sin(4 * phi),
                        delta * math.sin(6 * phi),
                        epsilon * math.sin(8 * phi)])


def UTMCentralMeridian (zone):
    '''Determines the central meridian for the given UTM zone.
    '''
    return DegToRad(int(zone)*6 - 183)


def FootpointLatitude (y):
    '''Computes the footpoint latitude for use in converting transverse
    * Mercator coordinates to ellipsoidal coordinates.
    *
    * Inputs:
    *   y - The UTM northing coordinate, in meters.
    *
    * Returns:
    *   The footpoint latitude, in radians.
    '''
    n = (sm_a - sm_b) / (sm_a + sm_b)
    alpha_ = ((sm_a + sm_b) / 2) * (1 + (n**2 / 4) + n**4 / 64)
    y_ = y / alpha_
    beta_ = (3 * n / 2) + (-27 * n**3 / 32) + (269 * n**5 / 512)
    gamma_ = (21 * math.pow (n, 2) / 16) + (-55 * math.pow (n, 4) / 32)
    delta_ = (151 * math.pow (n, 3) / 96) + (-417 * math.pow (n, 5) / 128)
    epsilon_ = 1097 * math.pow (n, 4) / 512
    return sum([y_,
                beta_ * math.sin (2 * y_),
                gamma_ * math.sin (4 * y_),
                delta_ * math.sin (6 * y_),
                epsilon_ * math.sin (8 * y_)])



def MapLatLonToXY (phi, lambda_, lambda_0, xy):
    '''
    * Converts a latitude/longitude pair to x and y coordinates in the
    * Transverse Mercator projection.  Note that Transverse Mercator is not
    * the same as UTM a scale factor is required to convert between them.
    *
    * Inputs:
    *    phi - Latitude of the point, in radians.
    *    lambda_ - Longitude of the point, in radians.
    *    lambda_0 - Longitude of the central meridian to be used, in radians.
    *
    * Outputs:
    *    xy - A 2-element array containing the x and y coordinates
    *         of the computed point.
    '''
    #pylint:disable=too-many-locals
    ep2 = (sm_a**2 - sm_b**2) / sm_b**2
    nu2 = ep2 * math.cos(phi)**2
    N = sm_a**2 / (sm_b * (1 + nu2)**0.5)
    t = math.tan(phi)
    t2 = t**2
    l = lambda_ - lambda_0

    #Precalculate coefficients for l**n in the equations below
    #so a normal human being can read the expressions for easting
    #and northing -- l**1 and l**2 have coefficients of 1
    l3coef = 1 - t2 + nu2
    l4coef = 5 - t2 + 9 * nu2 + 4 * (nu2 * nu2)
    l5coef = 5 - 18 * t2 + t**4 + 14 * nu2 - 58 * t2 * nu2
    l6coef = 61 - 58 * t2 + t**4 + 270 * nu2 - 330 * t2 * nu2
    l7coef = 61 - 479 * t2 + 179 * t**4 - t**6
    l8coef = 1385 - 3111 * t2 + 543 * t**4 - t**6

    #Calculate easting (x)
    xy[0] = sum([N * math.cos (phi) * l,
                 N / 6 * math.cos(phi)**3 * l3coef,
                 N / 120 * math.cos(phi)**5 * l5coef,
                 N / 5040 * math.cos(phi)**7 * l7coef])
    #Calculate northing (y)
    xy[1] = sum([ArcLengthOfMeridian(phi),
                 t / 2 * N * math.cos(phi)**2,
                 t / 24 * N * math.cos(phi)**4 * l4coef,
                 t / 720 * N * math.cos(phi)**6 * l6coef,
                 t / 40320 * N * math.cos(phi)**8 * l8coef])
    return xy



def MapXYToLatLon (x, y, lambda_0):
    '''Converts x and y coordinates in the Transverse Mercator projection to
    * a latitude/longitude pair.  Note that Transverse Mercator is not
    * the same as UTM a scale factor is required to convert between them.
    *
    * Inputs:
    *   x - The easting of the point, in meters.
    *   y - The northing of the point, in meters.
    *   lambda_0 - Longitude of the central meridian to be used, in radians.
    *
    * Outputs:
    *   philambda_ - A 2-element containing the latitude and longitude
    *               in radians.
    * Remarks:
    *   The local variables Nf, nuf2, tf, and tf2 serve the same purpose as
    *   N, nu2, t, and t2 in MapLatLonToXY, but they are computed with respect
    *   to the footpoint latitude phif.
    *
    *   x1frac, x2frac, x2poly, x3poly, etc. are to enhance readability and
    *   to optimize computations.
    '''
    #pylint:disable=too-many-locals

    #Get the value of phif, the footpoint latitude
    phif = FootpointLatitude (y)
    ep2 = (sm_a**2 - sm_b**2) / sm_b**2
    cf = math.cos(phif)
    nuf2 = ep2 * cf**2
    Nf = sm_a**2 / (sm_b * (1 + nuf2)**0.5)
    Nfpow = Nf
    tf = math.tan(phif)
    tf2 = tf * tf
    tf4 = tf2 * tf2
    # Calculate fractional coefficients for x**n in the equations
    # below to simplify the expressions for latitude and longitude. '''
    x1frac = 1 / (Nfpow * cf)

    Nfpow *= Nf   #now equals Nf**2)
    x2frac = tf / (2 * Nfpow)

    Nfpow *= Nf   #now equals Nf**3)
    x3frac = 1 / (6 * Nfpow * cf)

    Nfpow *= Nf   #now equals Nf**4)
    x4frac = tf / (24 * Nfpow)

    Nfpow *= Nf   #now equals Nf**5)
    x5frac = 1 / (120 * Nfpow * cf)

    Nfpow *= Nf   #now equals Nf**6)
    x6frac = tf / (720 * Nfpow)

    Nfpow *= Nf   #now equals Nf**7)
    x7frac = 1 / (5040 * Nfpow * cf)

    Nfpow *= Nf   #now equals Nf**8)
    x8frac = tf / (40320 * Nfpow)

    #Precalculate polynomial coefficients for x**n
    x2poly = -1 - nuf2
    x3poly = -1 - 2 * tf2 - nuf2
    x4poly = 5 + 3 * tf2 + 6 * nuf2 - 6 * tf2 * nuf2 \
        - 3 * nuf2**2 - 9 * tf2 * nuf2**2
    x5poly = 5 + 28 * tf2 + 24 * tf4 + 6 * nuf2 + 8 * tf2 * nuf2
    x6poly = -61 - 90 * tf2 - 45 * tf4 - 107 * nuf2 + 162 * tf2 * nuf2
    x7poly = -61 - 662 * tf2 - 1320 * tf4 - 720 * tf4 * tf2
    x8poly = 1385 + 3633 * tf2 + 4095 * tf4 + 1575 * tf4 * tf2

    #Calculate latitude
    philambda_ = [0, 0]
    philambda_[0] = sum([phif,
                         x2frac * x2poly * x**2,
                         x4frac * x4poly * x**4,
                         x6frac * x6poly * x**6,
                         x8frac * x8poly * x**8])
    #Calculate longitude
    philambda_[1] = sum([lambda_0,
                         x1frac * x,
                         x3frac * x3poly * x**3,
                         x5frac * x5poly * x**5,
                         x7frac * x7poly * x**7])
    return philambda_


def LatLonToUTMXY (lat, lon, zone):
    '''Converts a latitude/longitude pair to x and y coordinates in the
    * Universal Transverse Mercator projection.
    *
    * Inputs:
    *   lat - Latitude of the point, in radians.
    *   lon - Longitude of the point, in radians.
    *   zone - UTM zone to be used for calculating values for x and y.
    *          If zone is less than 1 or greater than 60, the routine
    *          will determine the appropriate zone from the value of lon.
    *
    * Outputs:
    *   xy - A 2-element array where the UTM x and y values will be stored.
    *
    * Returns:
    *   The UTM zone used for calculating the values of x and y.
    '''
    xy = MapLatLonToXY (lat, lon, UTMCentralMeridian (zone), xy)
    # Adjust easting and northing for UTM system
    xy[0] = xy[0] * UTMScaleFactor + 500000
    xy[1] = xy[1] * UTMScaleFactor
    if xy[1] < 0:
        xy[1] = xy[1] + 10000000
    return xy, zone


def UTMXYToLatLon (x, y):
    '''Converts x and y coordinates in the Universal Transverse Mercator
    * projection to a latitude/longitude pair.
    *
    * Inputs:
    *   x - The easting of the point, in meters.
    *   y - The northing of the point, in meters.
    *   zone - The UTM zone in which the point lies.
    *   southhemi - True if the point is in the southern hemisphere
    *               false otherwise.
    *
    * Outputs:
    *   latlon - A 2-element array containing the latitude and
    *            longitude of the point, in radians.
    '''
    x -= 500000
    x /= UTMScaleFactor
    #If in southern hemisphere, adjust y accordingly.
    if True:
        y -= 10000000
    y /= UTMScaleFactor
    return MapXYToLatLon (x, y, DegToRad(147))

if __name__ == '__main__':
    x, y = (float(n) for n in input('What XY coords?    ').split())
    print([RadToDeg(n) for n in UTMXYToLatLon (x, y)])






