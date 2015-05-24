/*
 
  sunrise_sunset.c
 
  Transcription in C of the algorithms described in
  Meeus, Astronomical algorithms, 1998
  ISBN 0-943396-61-1
 
 The MIT License (MIT)

 Copyright (c) 2015-2016 Fabrice Ferino

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 
*/

#include <math.h>

#define RADEG   ( 180.0 / M_PI )
#define DEGRAD  ( M_PI / 180.0 )

/* 360 deg = 24 hours = 24 * 60 minutes = 1440 minutes */
/* 1 deg = 1440/360 minutes = 4 minutes */
#define MIN_PER_DEG (4.0)
#define MIN_PER_DAY (1440.0)

const double kRiseOrSet = 90.833; /* 90°50' for atmospheric refraction at sunrise/sunset */
const double kCivilTwilight = 96.00; /* civil twilight */
const double kNauticalTwilight = 102.00; /* nautical twilight */
const double kAstronomicalTwilight = 108.00; /* astronomical twilight */

/* Trig functions in degrees */
#define sind(x)      sin(x*DEGRAD)
#define cosd(x)      cos(x*DEGRAD)
#define tand(x)      tan(x*DEGRAD)

#define asind(x)     (RADEG*asin(x))
#define acosd(x)     (RADEG*acos(d))
#define atand(x)     (RADEG*atan(x))
#define atan2d(y,x)  (RADEG*atan2(y,x))

#define DEG2RAD(x)  ((x)*DEGRAD)
#define RAD2DEG(x)  ((x)*RADEG)

#define DEG2MIN(x)  ((x)*MIN_PER_DEG)

/* p.163 */
#define DAYS_IN_CENTURY (36525.0)
#define JAN_1_2000_JD   (2451545.0)

double JulianCenturyFromJulianDay( double jd)
{
    return (jd - JAN_1_2000_JD) / DAYS_IN_CENTURY;  /* p.163, 25.1 */
}

double JulianDayFromJulianCentury( double centuryTime)
{
    return centuryTime * DAYS_IN_CENTURY + JAN_1_2000_JD;
}

int IsLeapYear(int year)
{
    return (0 == year % 4 && 0 != year % 100) ||
        0 == year % 400;
}


double DayOfYearFromJulianDay(double jd)
{
    /* Meeus, Astronomical algorithms, p. 63 */
    double z = floor(jd + 0.5); /* integral part */
    double f = (jd + 0.5) - z;  /* fractional part */
    double b = z + 1524.0;
    if (z >= 2299161.0)
    {
        double alpha = floor(( z - 1867216.25)/ 36524.25);
        b += 1 + alpha - floor(alpha/4.0);
    }
    double c = floor((b - 122.1)/365.25);
    double d = floor(365.25 *c);
    double e = floor( (b-d)/30.6001);
    double day = b - d - floor(30.6001 * e) + f;
    double month = (e < 14.0) ? e - 1.0 : e - 13.0;
    double year = ( month > 2) ? c - 4716.0 : c - 4715.0;
    /* ibid. p. 65 */

    double n = floor( (275.0*month)/9.0) -
        (IsLeapYear(year)? 1: 2) * floor((month+9.0)/12.0) + day - 30;
    return n;
}

/* p.147, 22.2 */
double MeanObliquityEcliptic(double centuryTime)
{
    double arcSeconds = 21.448 - centuryTime *
        (46.8150 + centuryTime *
         (0.00059 - centuryTime *  0.001813));
    double degrees = 23.0 + (26.0 + (arcSeconds/60.0))/ 60.0;
    return degrees;
}

double NormalizeDegrees(double degrees)
{
    while (degrees > 360.0)
    {
        degrees -= 360.0;
    }
    while (degrees < 0.0)
    {
        degrees += 360.0;
    }
    return degrees;
}


/* p. 163, 25.2 */
double GeometricMeanLongitudeSun( double centuryTime)
{
    double l0 = 280.46646 +  centuryTime * (36000.76983 + centuryTime * 0.0003032);
    return l0;
}

/* p. 163, 25.3 */
double GeometricMeanAnomalySun( double centuryTime)
{
    double m = 357.52911 + centuryTime * ( 35999.05029 - centuryTime * 0.0001537);
    return m;
}

/* p. 163, 25.4 */
double EccentricityEarth( double centuryTime)
{
    return 0.016708634 - centuryTime * (0.000042037 + centuryTime * 0.0000001267);
}

/* p. 164 */
double EquationOfCenterSunEx( double centuryTime, double meanAnomalySun)
{
    double mrad = DEG2RAD(meanAnomalySun);
    double sinm = sin(mrad);
    mrad += mrad;
    double sin2m = sin(mrad);
    mrad += mrad;
    double sin3m = sin(mrad);

    double c = (1.914602 - centuryTime * (0.004817 + centuryTime * 0.000014)) * sinm +
        (0.019993 - 0.000101 * centuryTime) * sin2m +
        0.000289 * sin3m;

    return c;
}

/* p. 164 */
double EquationOfCenterSun( double centuryTime)
{
    double meanAnomalySun = GeometricMeanAnomalySun(centuryTime);
    return EquationOfCenterSunEx( centuryTime, meanAnomalySun);
}


/* p. 164 */
double TrueLongitudeSun(double centuryTime)
{
    double l0 = GeometricMeanLongitudeSun( centuryTime);
    double c = EquationOfCenterSun( centuryTime);
    return l0 +c;
}

/* p. 164 */
double TrueAnomalySun( double centuryTime)
{
    double mdeg = GeometricMeanAnomalySun(centuryTime);
    double c = EquationOfCenterSunEx( centuryTime, mdeg);
    return mdeg + c;
}

/* p. 164 */
double Omega( double centuryTime)
{
    return 125.04 - 1934.136 * centuryTime;
}

double OmegaRad( double centuryTime)
{
    return DEG2RAD(Omega(centuryTime));
}

/* p. 164 */
double ApparentLongitudeSunEx( double centuryTime, double omegaRad)
{
    double trueLongitude = TrueLongitudeSun(centuryTime);

    return trueLongitude - 0.000569 - 0.00478 * sin(omegaRad);
}

double ApparentLongitudeSun( double centuryTime)
{
    return ApparentLongitudeSunEx( centuryTime, OmegaRad( centuryTime));
}

/* p. 165, 25.8 */
double ObliquityCorrectionEx( double centuryTime, double omegaRad)
{
    double moe = MeanObliquityEcliptic(centuryTime);
    return moe + 0.00256 * cos(omegaRad);
}

double ObliquityCorrection(double centuryTime)
{
    return ObliquityCorrectionEx( centuryTime, OmegaRad( centuryTime));
}

/* p. 165, 25.6 */
double SunRightAscensionRad( double centuryTime)
{
    double omegaRad = OmegaRad(centuryTime);
    double oc = ObliquityCorrectionEx( centuryTime, omegaRad);
    double al = ApparentLongitudeSunEx( centuryTime, omegaRad);
    return atan2(cosd(oc) * sind(al), cosd(al));
}

double SunRightAscension( double centuryTime)
{
    return RAD2DEG(SunRightAscensionRad(centuryTime));
}

/* p. 165, 25.7 */
double SunDeclinationRad( double centuryTime)
{
    double omegaRad = OmegaRad( centuryTime);
    double oc = ObliquityCorrectionEx( centuryTime, omegaRad);
    double al = ApparentLongitudeSunEx( centuryTime, omegaRad);

    return asin( sind(oc) * sind(al));
}

double SunDeclination( double centuryTime)
{
    return RAD2DEG( SunDeclinationRad(centuryTime));
}


/* p. 185, 28.3 */
double EquationOfTime( double centuryTime)
{
    double epsilon = ObliquityCorrection( centuryTime);
    double y = tand(epsilon/2.0);
    y *=y;

    double l0 = GeometricMeanLongitudeSun( centuryTime);
    double e = EccentricityEarth( centuryTime);
    double m = GeometricMeanAnomalySun( centuryTime);

    l0 = DEG2RAD(l0);
    double l0x2 = l0 + l0;
    double l0x4 = l0x2 + l0x2;

    m = DEG2RAD(m);
    double sinm = sin(m);
    double mx2 = m + m;

    double ex2 = e + e;

    double eRad = y * sin(l0x2) - ex2 * sinm + 2 * ex2 * y * sinm * cos(l0x2) -
        0.5 * y * y * sin(l0x4) - 1.25 * e * e * sin(mx2);

    /* convert radians to degrees to minutes of time */
    return DEG2MIN(RAD2DEG(eRad));
}


/* http://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF */
/* calculate the hour angle of the sun when it's at angle at the latitude */
double LocalHourAngleSunRad( double latitudeRad,
                             double declinationRad,
                             double angleRad)
{
    double cosHA = (cos(angleRad) - sin(latitudeRad) * sin(declinationRad)) /
        (cos(latitudeRad) * cos(declinationRad));
#if 0
    if (cosHA < -1.0 || cosHA > 1.0)
#endif

    return acos( cosHA); /* radians */
}

double UTCForSolarAngleAux( int rise, double jd,
                            double latitudeRad,
                            double angleRad)
{
    double centuryTime = JulianCenturyFromJulianDay( jd);
    double equationOfTime = EquationOfTime(centuryTime);
    double sunDeclinationRad = SunDeclinationRad(centuryTime);
    double hourAngle = LocalHourAngleSunRad( latitudeRad, sunDeclinationRad,
                                             angleRad);
    if (!rise) { hourAngle = -hourAngle; }

    return 720.0 - (4.0 * RAD2DEG(hourAngle)) - equationOfTime; /* minutes */
}

double UTCForSolarAngle( int rise, double jd, double latitude,
                         double angle)
{
    /* first call : approximation */
    double latitudeRad = DEG2RAD(latitude);
    double angleRad = DEG2RAD( angle);

    double firstTime = UTCForSolarAngleAux( rise, jd, latitudeRad, angleRad);
    double secondTime = UTCForSolarAngleAux( rise, jd + firstTime / MIN_PER_DAY,
                                             latitudeRad, angleRad);
    return secondTime; /* minutes */
}

/* p. 61, 7.1 */
double JulianDayEx( int y, int m, double dayFrac)
{
    int a, b;
    int JulianDay = 1; /* assume julian date */

    if ( y > 1582)
    {
        JulianDay = 0;
    }
    else if ( y == 1582)
    {
        if (m > 11)
        {
            JulianDay = 0;
        }
        else if ( m == 10)
        {
            if (dayFrac >= 15.0)
            {
                JulianDay = 0;
            }
        }
    }

    if (m <= 2)
    {
        --y;
        m += 12;
    }
    a = y/100;
    b = (JulianDay) ? 0 : 2 - a + a/4;

    return floor( 365.25 * ( y + 4716)) + floor(30.6001 * (m + 1)) +
        dayFrac + b - 1524.5;
}


double JulianDay( int year, int month, int day, int hour, int minute, int second)
{
    double dayFrac = ((double) day) + (hour / 24.0) + (minute / MIN_PER_DAY) +
        (second / (MIN_PER_DAY * 60.0));

    return JulianDayEx( year, month, dayFrac);
}


void D2DMS( double degreesFrac, int* degrees, int* minutes, double* seconds)
{
    *degrees = trunc(degreesFrac);
    double fraction = degreesFrac - *degrees;
    fraction *= 3600 * 1000; /* precision is milliseconds */
    long milliSeconds = (long) fraction;

    *minutes = (int) (milliSeconds / (60 * 1000));
    milliSeconds -= *minutes * 60L * 1000L;
    *seconds = milliSeconds / 1000.0;
}
