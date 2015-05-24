/*
  main.c

  SolarTimes
 
  Tests of the routines in sunset_sunrise.c/.h

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

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "sunrise_sunset.h"

/* struct and data used in the JulianDayTest */
typedef struct JulianDayTestInfo
{
    int y;
    int m;
    double dfrac;
    double jd;
} JulianDayTestInfo;

JulianDayTestInfo gJulianDayTestInfo[] =
{
    {   1957, 10,  4.81, 2436116.31 },
    {    333,  1,  27.5, 1842713.0 },
    {   2000,  1,   1.5, 2451545.0 },
    {   1999,  1,   1.0, 2451179.5 },
    {   1987,  1,  27.0, 2446822.5 },
    {   1987,  6,  19.5, 2446966.0 },
    {   1988,  1,  27.0, 2447187.5 },
    {   1988,  6,  19.5, 2447332.0 },
    {   1900,  1,     1, 2415020.5 },
    {   1600,  1,   1.0, 2305447.5 },
    {   1600, 12,  31.0, 2305812.5 },
    {    837,  4,  10.3, 2026871.8 },
    {  -1000,  7,  12.5, 1356001.0 },
    {  -4712,  1,   1.5,       0.0 }
};


typedef struct HourMinute
{
    int hour;
    int minute;
} HourMinute;

/* latitudes used in the Nautical Almanac */
const double kLatitudes[] =
{
    72.0,
    70.0,
    68.0,
    66.0,
    64.0,
    62.0,
    60.0,
    58.0,
    56.0,
    54.0,
    52.0,
    50.0,
    45.0,
    40.0,
    35.0,
    30.0,
    20.0,
    10.0,
     0.0,
    -10.0,
    -20.0,
    -30.0,
    -35.0,
    -40.0,
    -45.0,
    -50.0,
    -52.0,
    -54.0,
    -56.0,
    -58.0,
    -60.0
};

const int kNumLatitudes = sizeof(kLatitudes) / sizeof(kLatitudes[0]);

/* struct used in SunPosition test */
typedef struct SunPositionTest
{
    int y;
    int m;
    double dfrac;

    HourMinute nauticalRise[kNumLatitudes];
    HourMinute civilRise[kNumLatitudes];
    HourMinute sunRise[kNumLatitudes];
    HourMinute sunSet[kNumLatitudes];
    HourMinute civilSet[kNumLatitudes];
    HourMinute nauticalSet[kNumLatitudes];

} SunPositionTest;


int JulianDayTests()
{
    int i, retVal = 0;

    for (i = 0; i < sizeof(gJulianDayTestInfo)/sizeof(gJulianDayTestInfo[0]); ++i)
    {
        JulianDayTestInfo* t = gJulianDayTestInfo + i;
        double julianDay = JulianDayEx( t->y, t->m, t->dfrac);

        if (julianDay != t->jd)
        {
            printf("test %d, expected %lf instead of %lf\n",
                   i, t->jd, julianDay);
            ++retVal;
        }
    }
    return retVal;
}



int D2DMSTest()
{
    int deg, minutes;
    double seconds;

    D2DMS(121.1350000, &deg, &minutes, &seconds);

    return (deg == 121 && minutes == 8 && seconds == 6.0) ? 0 : 1;
}

int ObliquityTest()
{
    double jd = JulianDayEx(1987, 4, 10.0);
    double centuryTime = JulianCenturyFromJulianDay(jd);
    double e = MeanObliquityEcliptic(centuryTime);

    int degrees, minutes;
    double seconds;

    D2DMS(e, &degrees, &minutes, &seconds);

    return (degrees == 23 && minutes == 26 && seconds == 27.407)? 0 : 1;
}

int GeometricMeanLongitudeSunTest()
{
    double jd = JulianDayEx(1992, 10, 13.0);
    double centuryTime = JulianCenturyFromJulianDay(jd);
    double result = GeometricMeanLongitudeSun(centuryTime);
    result = NormalizeDegrees(result);
    /* result is 201.80720 deg */
    int resultInt = round( result * 100000.0);
    return (resultInt == 20180720 )? 0 : 1;
}

int GeometricMeanAnomalySunTest()
{
    double jd = JulianDayEx(1992, 10, 13.0);
    double centuryTime = JulianCenturyFromJulianDay(jd);
    double result = GeometricMeanAnomalySun(centuryTime);
    result = NormalizeDegrees(result);
    /* result is 278.99397 deg */
    int resultInt = round( result * 100000.0);
    return (resultInt == 27899397 )? 0 : 1;
}


void FormatMinutes( double minutesd, char s[/*6*/])
{
    int minutes, hh, mm;
    if (isnan(minutesd) || minutesd < 0)
    {
        strcpy(s, " N/A ");
        return;
    }
    minutes = (int) round(minutesd);
    hh = minutes / 60;
    mm = minutes % 60;
    snprintf(s, 6, "%02d %02d", hh, mm);
}

typedef char TimeString[6];

int SunRiseTest(double jd)
{
    int i;

    for ( i = 0; i < kNumLatitudes; ++i)
    {
        double nauticalTwilightBegin = UTCForSolarAngle(1, jd, kLatitudes[i], kNauticalTwilight);
        double civilTwilightBegin = UTCForSolarAngle(1, jd, kLatitudes[i], kCivilTwilight);
        double sunRise = UTCForSolarAngle(1, jd, kLatitudes[i], kRiseOrSet);
        double sunSet = UTCForSolarAngle(0, jd, kLatitudes[i], kRiseOrSet);
        double civilTwilightEnd = UTCForSolarAngle(0, jd, kLatitudes[i], kCivilTwilight);
        double nauticalTwilightEnd = UTCForSolarAngle(0, jd, kLatitudes[i], kNauticalTwilight);


        TimeString ntb, ctb, sr, ss, cte, nte;

        FormatMinutes(nauticalTwilightBegin, ntb);
        FormatMinutes(civilTwilightBegin, ctb);
        FormatMinutes(sunRise, sr);
        FormatMinutes(sunSet, ss);
        FormatMinutes(civilTwilightEnd, cte);
        FormatMinutes(nauticalTwilightEnd, nte);


        printf("| %+2.0lf | %s | %s | %s | %s | %s | %s |\n",
               kLatitudes[i],
               ntb, ctb, sr, ss, cte, nte);

    }

    return 0;
}


int SunRiseTests()
{
    SunRiseTest(JulianDayEx(1994, 5, 8.0)); /* match results of Nautical Almanac */
    printf("\n\n");

    SunRiseTest(JulianDayEx(1994, 12, 9.0));
    printf("\n\n");

    SunRiseTest(JulianDayEx(1994, 12, 10.0));
    printf("\n\n");

    SunRiseTest(JulianDayEx(1994, 12, 10.5));
    printf("\n\n");
    
    SunRiseTest(JulianDayEx(2012, 1, 2.0));
    printf("\n\n");

    SunRiseTest(JulianDayEx(2012, 11, 21.0));
    printf("\n\n");


    return 0;
}

int main( int argc, char* argv[])
{
    int retVal = 0;

    if (argc == 1)
    {
        /* run built-in tests */
        retVal += JulianDayTests();
        retVal += D2DMSTest();
        retVal += ObliquityTest();
        retVal += GeometricMeanLongitudeSunTest();
        retVal += GeometricMeanAnomalySunTest();

        SunRiseTests();
    }

    return retVal;
}



