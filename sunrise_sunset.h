/* 
 sunrise_sunset.h

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
#ifndef SUNRISE_SUNSET_HEADER
#define SUNRISE_SUNSET_HEADER


extern const double kRiseOrSet;
extern const double kCivilTwilight;
extern const double kNauticalTwilight;
extern const double kAstronomicalTwilight;

double JulianCenturyFromJulianDay(double jd);
double JulianDayFromJulianCentury(double centuryTime);
int IsLeapYear(int year);
double DayOfYearFromJulianDay(double jd);
double MeanObliquityEcliptic(double centuryTime);
double GeometricMeanLongitudeSun(double centuryTime);
double GeometricMeanAnomalySun(double centuryTime);
double EccentricityEarth(double centuryTime);
double EquationOfCenterSunEx(double centuryTime, double meanAnomalySun);
double EquationOfCenterSun(double centuryTime);
double TrueLongitudeSun(double centuryTime);
double TrueAnomalySun(double centuryTime);
double Omega(double centuryTime);
double OmegaRad(double centuryTime);
double ApparentLongitudeSunEx(double centuryTime, double omegaRad);
double ApparentLongitudeSun(double centuryTime);
double ObliquityCorrectionEx(double centuryTime, double omegaRad);
double ObliquityCorrection(double centuryTime);
double SunRightAscensionRad(double centuryTime);
double SunRightAscension(double centuryTime);
double SunDeclinationRad(double centuryTime);
double SunDeclination(double centuryTime);
double EquationOfTime(double centuryTime);
double LocalHourAngleSunRad(double latitudeRad, double declinationRad, double angleRad);
double UTCForSolarAngleAux(int rise, double jd, double latitudeRad, double angleRad);
double UTCForSolarAngle(int rise, double jd, double latitude, double angle);
double JulianDayEx(int y, int m, double dayFrac);
double JulianDay(int year, int month, int day, int hour, int minute, int second);
void D2DMS( double degreesFrac, int* degrees, int* minutes, double* seconds);
double NormalizeDegrees(double degrees);


#endif