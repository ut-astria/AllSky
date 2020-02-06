//////////////////////////////
//AllSky project
//
//Daniel Kucharski, allsky@utexas.edu
//17 Jan 2020
//////////////////////////////




#include <iostream>
#include <limits>

#include <ctime>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <limits.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <time.h>

#include <math.h>
#include <errno.h>
#include <windows.h>//system time etc.

#include <sstream>      // std::stringstream
#include <fstream>      // std::ofstream

#include <fcntl.h>//Allows the user to perform a variety of operations on open files.

#include <sys/types.h>
#include "stdafx.h"

#include "AA+.h"
#include "SGP4.h"


#include "allsky.h"









/*
precision and significant digits
http://www.cplusplus.com/doc/tutorial/variables/
Name	                Description						Size and Range
char	                Character or small integer.		1B, signed: -128 to 127 unsigned: 0 to 255
short int (short)		Short Integer.	                2B, signed: -32768 to 32767 unsigned: 0 to 65535
int						Integer.						4B, signed: -2147483648 to 2147483647 unsigned: 0 to 4294967295
long int (long)	        Long integer.	                4B	signed: -2147483648 to 2147483647 unsigned: 0 to 4294967295
bool	                Boolean value.                  1B	true or false
float	                Floating point number.	        4B, 32b, 1.17*10^-38..+-3.4*10^38 (6 significant digits)
double	                Double precision                8B, 64b, 2.22*10^-308..+-1.79*10^308 (15 significant digits)
long double				Long double precision.	        80b	3.36*10^-4932..+-1.18*10^4932 (18 significant digits)
wchar_t	                Wide character.	                2 or 4 bytes	1 wide character
*/

const double pi = 3.14159265358979323846;
const double pihalf = pi / 2.;
const double pitwo = pi * 2.;
const double deg2rad = pi / 180.0;
const double rad2deg = 180.0 / pi;
double pi2 = 2.*pi;
const double AU_to_meter = 149597870700.;//1 AU [m]


///////////////////////////////////////
//Operators

v3 operator + (const v3& V1, const v3& V2)
{
	v3 Out;
	Out.x = V1.x + V2.x;
	Out.y = V1.y + V2.y;
	Out.z = V1.z + V2.z;
	return Out;
}

v3 operator - (const v3& V1, const v3& V2)
{
	v3 Out;
	Out.x = V1.x - V2.x;
	Out.y = V1.y - V2.y;
	Out.z = V1.z - V2.z;
	return Out;
}

v3 operator - (const v3& V1)
{
	v3 Out;
	Out.x = -V1.x;
	Out.y = -V1.y;
	Out.z = -V1.z;
	return Out;
}

v3 operator * (const v3 &V1, const v3 &V2)
{
	v3 Out;
	Out.x = V1.x*V2.x;
	Out.y = V1.y*V2.y;
	Out.z = V1.z*V2.z;
	return Out;
}

v3 operator * (const v3 &V, const double &C)
{
	v3 out;
	out.x = V.x*C;
	out.y = V.y*C;
	out.z = V.z*C;
	return out;
}

v3 operator += (v3& V1, const v3& V2)
{
	V1.x += V2.x;
	V1.y += V2.y;
	V1.z += V2.z;
	return V1;
}

v3 operator -= (v3& V1, const v3& V2)
{
	V1.x -= V2.x;
	V1.y -= V2.y;
	V1.z -= V2.z;
	return V1;
}

v3 operator * (const double &C, v3 &V)
{
	v3 out;
	out.x = V.x*C;
	out.y = V.y*C;
	out.z = V.z*C;
	return out;
}

v3 operator *= (v3 &V, const double &C)
{
	V.x *= C;
	V.y *= C;
	V.z *= C;
	return V;
}

v3 operator /= (v3 &V, const double &C)
{
	V.x /= C;
	V.y /= C;
	V.z /= C;
	return V;
}

v3 operator * (const v33 & M, const v3& V)//matrix vector product
{
	v3 out;
	out.x = M.v[0][0] * V.x + M.v[0][1] * V.y + M.v[0][2] * V.z;
	out.y = M.v[1][0] * V.x + M.v[1][1] * V.y + M.v[1][2] * V.z;
	out.z = M.v[2][0] * V.x + M.v[2][1] * V.y + M.v[2][2] * V.z;
	return out;
}

v33 operator * (const v33 &L, const v33 &R)//left matrix, right matrix product
{
	v33 M;//out matrix

	M.v[0][0] = L.v[0][0] * R.v[0][0] + L.v[0][1] * R.v[1][0] + L.v[0][2] * R.v[2][0];
	M.v[0][1] = L.v[0][0] * R.v[0][1] + L.v[0][1] * R.v[1][1] + L.v[0][2] * R.v[2][1];
	M.v[0][2] = L.v[0][0] * R.v[0][2] + L.v[0][1] * R.v[1][2] + L.v[0][2] * R.v[2][2];

	M.v[1][0] = L.v[1][0] * R.v[0][0] + L.v[1][1] * R.v[1][0] + L.v[1][2] * R.v[2][0];
	M.v[1][1] = L.v[1][0] * R.v[0][1] + L.v[1][1] * R.v[1][1] + L.v[1][2] * R.v[2][1];
	M.v[1][2] = L.v[1][0] * R.v[0][2] + L.v[1][1] * R.v[1][2] + L.v[1][2] * R.v[2][2];

	M.v[2][0] = L.v[2][0] * R.v[0][0] + L.v[2][1] * R.v[1][0] + L.v[2][2] * R.v[2][0];
	M.v[2][1] = L.v[2][0] * R.v[0][1] + L.v[2][1] * R.v[1][1] + L.v[2][2] * R.v[2][1];
	M.v[2][2] = L.v[2][0] * R.v[0][2] + L.v[2][1] * R.v[1][2] + L.v[2][2] * R.v[2][2];

	return M;
}

v3 operator / (const v3 &V, const double &C)
{
	v3 out;
	out.x = V.x / C;
	out.y = V.y / C;
	out.z = V.z / C;
	return out;
}

v2 operator + (const v2 &V1, const v2 &V2)
{
	v2 Out;
	Out.x = V1.x + V2.x;
	Out.y = V1.y + V2.y;
	return Out;
}

v2 operator + (const v2 &V1, const double &C)
{
	v2 Out;
	Out.x = V1.x + C;
	Out.y = V1.y + C;
	return Out;
}

v2 operator - (const v2 &V1, const v2 &V2)
{
	v2 Out;
	Out.x = V1.x - V2.x;
	Out.y = V1.y - V2.y;
	return Out;
}

v2 operator - (const v2 &V1, const double &C)
{
	v2 Out;
	Out.x = V1.x - C;
	Out.y = V1.y - C;
	return Out;
}
///////////////////////////////////////




///////////////////////////////////////
//SGP4
//source: http://www.celestrak.com/software/vallado-sw.asp
char typerun, typeinput, opsmode;
gravconsttype  whichconst;
int whichcon;

int calculate_sat_positions_from_TLE(
	double MJD,
	char TLE_Line_1[130],
	char TLE_Line_2[130],
	double ICSm[3]//out
);
///////////////////////////////////////

///////////////////////////////////////
//Astronomical Algorithms
//source: http://www.naughter.com/aa.html
//example functions in AATest.cpp

double Distance_to_the_SUN_AU(double JD);
///////////////////////////////////////

///////////////////////////////////////
//Date time MJD
void give_present_time(
	unsigned short *YearNow,
	unsigned short *MonthNow,
	unsigned short *DayNow,
	unsigned short *HourNow,
	unsigned short *MinNow,
	unsigned short *SecNow,
	unsigned short *mSecNow
);

double present_MJD_full();

double MJD_full(
	int year,//2006
	int month,
	int day,
	int hour,
	int minute,
	double second
);

void print_date(unsigned short Year, unsigned short MonthNow, unsigned short DayNow);
///////////////////////////////////////

///////////////////////////////////////
//Geodesy and coordinate systems
v3 Give_ITRFpoint_for_WGS84_ellipse_Lon_Lat_height(//out: ITRF_m
	ae Ellipse_rad,
	double Ellipse_height_m
);

double GAST_rad(double MJD);
///////////////////////////////////////


int Round_value(long double Value);
void ry(double Angle_rad, double M[3][3]);
void rz(double Angle_rad, double M[3][3]);
void inv(double in[3][3], double out[3][3]);
void mul(double left[3][3], double right[3][3], double out[3][3]);
void mul(double M[3][3], double V[3], double out[3]);
double det(double in[3][3]);

float RSO_apparent_magnitude(
	double MJD,
	double Sun_ICS_m[3],
	double Satellite_ICS_m[3],
	double Observer_ICS_m[3]
);

double norm(double V[3]);

void topo_AzEl_rad(//topocentric
	double Observer_TCS_m[3],//TCS = Terrestrial Coordinate System (ITRF)
	double Point_TCS_m[3],
	double *out_az,
	double *out_el
);

double shadow_function_Hubaux(//shadow function
	double MJD,
	v3 Sat_ICSm,
	ae Sun_ICS_rad,
	double *Out_Distance_Sun_Earth_m,
	v3 &Out_Sun_ICS_m
);

double norm(v3 v);//returns vector length
v3 get_v3(ae Vrad);//get unit vector at this orientation
v3 get_v3(ae Vrad, double Length);//get unit vector at this orientation
ae get_ae(v3 v);
double IArad(v3 A, v3 B);
v3 normalize(v3 v);//returns unit vector
double dot(v3 left, v3 right);
void pass_prediction();


//endd

//main function starts here
int main(int argc, const char *argv[])
{
	printf("\nAll-Sky\n");

	/*
	//////////////////////////////////////
	//TLE lines for TOPEX/Poseidon
	//https://www.n2yo.com/satellite/?s=22076
	char TLE_line1[70];//69+1
	char TLE_line2[70];//69+1
	sprintf_s(TLE_line1, "%s", "");
	sprintf_s(TLE_line2, "%s", "");
	sprintf_s(TLE_line1, "%s", "1 22076U 92052A   20020.32091875 -.00000058 +00000-0 +14373-4 0  9990");
	sprintf_s(TLE_line2, "%s", "2 22076 066.0418 313.3236 0008080 269.6122 268.9326 12.81024647283823");

	printf("\nTopex TLE\n");
	printf("%s\n", TLE_line1);
	printf("%s\n", TLE_line2);


	double jd;
	double jdFrac;
	SGP4Funcs::jday(2020, 1, 21, 2, 57, 7.2, jd, jdFrac);//test epoch
	double JD = jd + jdFrac;
	double MJD_UTC = JD - 2400000.5;
	double Sat_ICSm[3];//out
	int err = calculate_sat_positions_from_TLE(MJD_UTC, TLE_line1, TLE_line2, Sat_ICSm);

	printf("\nRSO position (J2000)\n");
	printf("x = %.1lf m, y = %.1lf m, z = %.1lf m\n", Sat_ICSm[0], Sat_ICSm[1], Sat_ICSm[2]);

	//get Sun position	
	double Epsilon = CAANutation::TrueObliquityOfEcliptic(JD);
	//sun apparent equatorial coordinates
	double Apparent_ecliptic_Lon_deg = CAASun::ApparentEclipticLongitude(JD, false);
	double Apparent_ecliptic_Lat_deg = CAASun::ApparentEclipticLatitude(JD, false);
	CAA2DCoordinate Sun1 = CAACoordinateTransformation::Ecliptic2Equatorial(Apparent_ecliptic_Lon_deg, Apparent_ecliptic_Lat_deg, Epsilon);
	double Sun_RA_rad = CAACoordinateTransformation::HoursToRadians(Sun1.X);
	double Sun_Dec_rad = CAACoordinateTransformation::DegreesToRadians(Sun1.Y);
	double Earth_Sun_m = Distance_to_the_SUN_AU(JD)*AU_to_meter;
	double Sun_ICS_m[3];
	double rXY = fabs(cos(Sun_Dec_rad));
	Sun_ICS_m[0] = Earth_Sun_m * rXY * cos(Sun_RA_rad);
	Sun_ICS_m[1] = Earth_Sun_m * rXY * sin(Sun_RA_rad);
	Sun_ICS_m[2] = sin(Sun_Dec_rad);



	//Observer ITRF2000
	//McDonald Observatory
	double Observer_ITRF_m[3] = { -1330021.0 , -5328401.8, 3236480.7 };
	double Observer_elev_m = 2006.2210;
	//get observer inertial coordinates at given epoch; no polar motion	
	double ICStoTCS[3][3];
	double TCStoICS[3][3];
	rz(GAST_rad(MJD_UTC), ICStoTCS);
	inv(ICStoTCS, TCStoICS);
	double Observer_ICS_m[3];

	mul(TCStoICS, Observer_ITRF_m, Observer_ICS_m);

	//RSO topocentric elevation
	/*
	void topo_AzEl_rad(//topocentric
	double Observer_TCS_m[3],//TCS = Terrestrial Coordinate System (ITRF)
	double Point_TCS_m[3],
	double *out_az,
	double *out_el
	*/

	/*
	//Use shadow function
	v3 SatICSm, SunICS_m;
	SatICSm.x = Sat_ICSm[0];
	SatICSm.y = Sat_ICSm[1];
	SatICSm.z = Sat_ICSm[2];
	ae Sun_ICS_rad;
	Sun_ICS_rad.az = Sun_RA_rad;
	Sun_ICS_rad.el = Sun_Dec_rad;

	//shadow function
	float ShadowFunction = shadow_function_Hubaux(MJD_UTC, SatICSm, Sun_ICS_rad, &Earth_Sun_m, SunICS_m);

	//apparent magnitude
	float RSO_app_mag = 999.;
	
	if (ShadowFunction > 0.)
	RSO_app_mag = RSO_apparent_magnitude(MJD_UTC, Sun_ICS_m, Sat_ICSm, Observer_ICS_m);

	printf("\nRSO apparent magnitude = %.1lf\n", RSO_app_mag);
	*/

	//SGP4Funcs::invjday(JD, 0., year, month, day, hour, minute, second);




	pass_prediction();




	printf("\npress any key to exit\n");
	_getch();

	return 0;



}

int Round_value(long double Value)
{
	int Coeff1 = int(fabs(Value));
	int Result = 0;
	if (fabs(fabs(Value) - double(Coeff1)) >= 0.5)
		Result = Coeff1 + 1;
	else
		Result = Coeff1;
	if (Value < 0.)
		Result *= -1.;
	return Result;
}


int calculate_sat_positions_from_TLE(
	double MJD,
	char TLE_Line_1[130],
	char TLE_Line_2[130],
	double ICSm[3]//out
)
{
	int ReturnError = 0;
	//0 - no error
	//1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
	//2 - mean motion less than 0.0
	//3 - pert elements, ecc < 0.0  or  ecc > 1.0
	//4 - semi-latus rectum < 0.0
	//5 - epoch elements are sub-orbital
	//6 - satellite has decayed

	ICSm[0] = 0.;
	ICSm[1] = 0.;
	ICSm[2] = 0.;

	char str[2];
	double ro[3];
	double vo[3];
	errno_t err;
	// ----------------------------  locals  -------------------------------
	double p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper;
	double sec, jd, jdFrac, rad, startmfe, stopmfe, deltamin;
	int  year; int mon; int day; int hr; int min;
	char longstr1[130];
	typedef char str3[4];
	str3 monstr[13];
	char outname[64];
	char longstr2[130];
	elsetrec satrec;

	sprintf_s(longstr1, "%s", "");
	sprintf_s(longstr2, "%s", "");
	sprintf_s(longstr1, "%s", TLE_Line_1);
	sprintf_s(longstr2, "%s", TLE_Line_2);

	rad = 180.0 / pi;

	//opsmode = 'a' best understanding of how afspc code works
	//opsmode = 'i' improved sgp4 resulting in smoother behavior
	opsmode = 'a';

	//typerun = c compare 1 year of full satcat data
	//typerun = v verification run, requires modified elm file with start, stop, and delta times
	//typerun = m manual operation- either mfe, epoch, or day of yr
	char typerun = 'c';
	typeinput = 'e';//input start stop ymd hms

	whichcon = 84;
	whichconst = wgs84;

	// includes initialization of sgp4	
	SGP4Funcs::twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst,
		startmfe, stopmfe, deltamin,
		satrec);

	startmfe = 0.;
	stopmfe = 0.;
	deltamin = 0.;

	double ReferenceEpochJD = satrec.jdsatepoch + satrec.jdsatepochF;

	// call the propagator to get the initial state vector value
	// no longer need gravconst since it is assigned in sgp4init
	SGP4Funcs::sgp4(satrec, 0.0, ro, vo);

	double tsince = (MJD + 2400000.5 - ReferenceEpochJD)*1440.;//[minute]
	SGP4Funcs::sgp4(satrec, tsince, ro, vo);

	for (int i = 0; i < 3; i++)
		ICSm[i] = ro[i] * 1.e3;//km->m

	return ReturnError;
}

double Distance_to_the_SUN_AU(double JD)
{
	//astronomical algorithms, chapter 24, solar coordinates
	double Distance = 0.;
	double T = (JD - 2451545.) / 36525.;
	double T2 = T * T;
	double T3 = T * T*T;
	// double L0=280.46645+36000.76983*T+0.0003032*T2;//deg
	double M_deg = 357.52910 + 35999.05030*T - 0.0001559*T2 - 0.00000048*T3;//deg
	double M_rad = M_deg * deg2rad;
	double e = 0.016708617 - 0.000042037*T - 0.0000001236*T2;
	double C_deg = (1.914600 - 0.004817*T - 0.000014*T2)*sin(M_rad) + (0.019993 - 0.000101*T)*sin(2.*M_rad) + 0.000290*sin(3.*M_rad);
	double C_rad = C_deg * deg2rad;
	double v_rad = M_rad + C_rad;
	double Numerator = 1.000001018*(1. - e * e);
	double Denominator = 1 + e * cos(v_rad);
	if (Denominator != 0.)
		Distance = Numerator / Denominator;
	return Distance;
}

double GAST_rad(double MJD)
{
	//GAST/sideral time, increases over time
	double Result = 0.;
	long double Diff_time = MJD - 33282.;
	long double Parameter_1 = (100.075542 + 360.*(Diff_time - double(int(Diff_time))) + 0.985647346*Diff_time + 2.9e-13*std::pow(Diff_time, 2.)) / 360.;
	long double Parameter_2;
	Result = modfl(Parameter_1, &Parameter_2)*pi2;
	return Result;//rad
}

//rotation matrix RY
void ry(double Angle_rad, double M[3][3])
{
	//Seeber
	double s = sin(Angle_rad);
	double c = cos(Angle_rad);
	M[0][0] = c;	M[0][1] = 0.0;	M[0][2] = -s;
	M[1][0] = 0.0;	M[1][1] = 1.0;	M[1][2] = 0.0;
	M[2][0] = s;	M[2][1] = 0.0;	M[2][2] = c;
}

//rotation matrix RZ
void rz(double Angle_rad, double M[3][3])
{
	//Seeber
	double s = sin(Angle_rad);
	double c = cos(Angle_rad);
	M[0][0] = c;	M[0][1] = s;	M[0][2] = 0.0;
	M[1][0] = -s;	M[1][1] = c;	M[1][2] = 0.0;
	M[2][0] = 0.0;	M[2][1] = 0.0;	M[2][2] = 1.0;
}

void inv(double in[3][3], double out[3][3])
{
	//inverse of 3x3 matrix
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			out[i][j] = 0.;

	double D = det(in);//get determinant

	if (D != 0.)
	{
		double a11 = in[0][0];
		double a12 = in[0][1];
		double a13 = in[0][2];
		double a21 = in[1][0];
		double a22 = in[1][1];
		double a23 = in[1][2];
		double a31 = in[2][0];
		double a32 = in[2][1];
		double a33 = in[2][2];

		out[0][0] = (a22*a33 - a23 * a32) / D;
		out[0][1] = (a13*a32 - a12 * a33) / D;
		out[0][2] = (a12*a23 - a13 * a22) / D;

		out[1][0] = (a23*a31 - a21 * a33) / D;
		out[1][1] = (a11*a33 - a13 * a31) / D;
		out[1][2] = (a13*a21 - a11 * a23) / D;

		out[2][0] = (a21*a32 - a22 * a31) / D;
		out[2][1] = (a12*a31 - a11 * a32) / D;
		out[2][2] = (a11*a22 - a12 * a21) / D;
	}
}

double det(double in[3][3])
{
	//calculate determinant of 3x3 matrix	
	double D =
		in[0][0] * in[1][1] * in[2][2] -
		in[0][0] * in[1][2] * in[2][1] -
		in[0][1] * in[1][0] * in[2][2] +
		in[0][1] * in[1][2] * in[2][0] +
		in[0][2] * in[1][0] * in[2][1] -
		in[0][2] * in[1][1] * in[2][0];

	return D;
}


float RSO_apparent_magnitude(
	double MJD,
	double Sun_ICS_m[3],
	double Satellite_ICS_m[3],
	double Observer_ICS_m[3]
)
{
	double apparent_magnitude = 0.;

	//get vectors 
	double SatSun_ICSuv[3];
	double SatObserver_ICSuv[3];
	for (int i = 0; i < 3; i++)
	{
		SatSun_ICSuv[i] = Sun_ICS_m[i] - Satellite_ICS_m[i];
		SatObserver_ICSuv[i] = Observer_ICS_m[i] - Satellite_ICS_m[i];
	}
	//normalize vectors
	double lengthA = norm(SatSun_ICSuv);
	double sat_observer_distance_m = norm(SatObserver_ICSuv);
	bool okay = true;//okay flag
	if ((lengthA > 0.) && (sat_observer_distance_m > 0.))
	{
		for (int i = 0; i < 3; i++)
		{
			SatSun_ICSuv[i] /= lengthA;
			SatObserver_ICSuv[i] /= sat_observer_distance_m;
		}
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			SatSun_ICSuv[i] = 0.;
			SatObserver_ICSuv[i] = 0.;
		}
		okay = false;
	}

	if (okay)
	{
		//get phase angle between the unit vectors
		double Phase_rad = acos(
			SatSun_ICSuv[0] * SatObserver_ICSuv[0] +
			SatSun_ICSuv[1] * SatObserver_ICSuv[1] +
			SatSun_ICSuv[2] * SatObserver_ICSuv[2]);

		//calculate RSO magnitude 
		//assume 1-m diameter Lambertian sphere
		double RSO_sphere_radius_m = 0.5;
		double Surface_area_m2 = pi * RSO_sphere_radius_m * RSO_sphere_radius_m;

		//phase reflection function for a spherical object
		double Cd = 0.8;//diffuse reflection coefficient
		double Psi_sphere = (Cd / pi)*(sin(Phase_rad) + (pi - Phase_rad)*cos(Phase_rad));

		//Irradiation at the ground observing point
		double Irradiation_sun = 1368.0;// [W/m^2], solar constant; solar irradiation at 1AU
		double Irradiation_topo = Irradiation_sun * Surface_area_m2*Psi_sphere / (sat_observer_distance_m*sat_observer_distance_m);

		double Sun_apparent_mag = -26.74;//V, apparent magnitude of Sun
		double RSO_apparent_mag = Sun_apparent_mag - 2.5*log10(Irradiation_topo / Irradiation_sun);//absolute magnitude of RSO
		//Apparent magnitudes, based on measurements, differ from the absolute magnitudes. For absolute magnitudes, reference values are needed.

		apparent_magnitude = RSO_apparent_mag;
	}

	return apparent_magnitude;
}

double norm(double V[3])
{
	return std::pow(V[0] * V[0] + V[1] * V[1] + V[2] * V[2], 0.5);
}

void topo_AzEl_rad(//topocentric
	double Observer_TCS_m[3],//TCS = Terrestrial Coordinate System (ITRF)
	double Point_TCS_m[3],
	double *out_az,
	double *out_el
)
{
	*out_az = 0.;
	*out_el = 0.;

	double Observer_TCS_uv[3];
	double length = norm(Observer_TCS_m);
	bool okay = true;//okay flag
	if (length > 0.)
	{
		for (int i = 0; i < 3; i++)
			Observer_TCS_uv[i] /= length;

		double Observer_lon_rad = 0., Observer_lat_rad = 0.;

		Observer_lon_rad = atan2(Observer_TCS_uv[1], Observer_TCS_uv[0]);
		if (Observer_lon_rad < 0.) Observer_lon_rad += pitwo;

		Observer_lat_rad = pihalf - acos(Observer_TCS_uv[2]);

		//Azimuth - from North (true North) CW towards East
		//North Az=0, East Az=90, South Az=180, West Az=270
		double toptrt[3];
		for (int i = 0; i < 3; i++)
			toptrt[i] = Point_TCS_m[i] - Observer_TCS_m[i];

		double Rz[3][3], Ry[3][3], M[3][3];
		rz(Observer_lon_rad, Rz);
		ry(pihalf - Observer_lat_rad, Ry);
		mul(Ry, Rz, M);

		double out[3];
		mul(M, toptrt, out);

		out[0] *= -1.;//because azimuth is counted from North vector towards East: https://commons.wikimedia.org/wiki/File:ECEF_ENU_Longitude_Latitude_relationships.svg

		*out_az = atan2(out[1], out[0]);

		if ((*out_az) < 0.) (*out_az) += pitwo;
		(*out_el) = pihalf - acos(out[2]);
	}
}

void mul(double left[3][3], double right[3][3], double out[3][3])
{
	//matrix multiplication
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			out[i][j] = 0.0;
			for (int k = 0; k < 3; k++)
				out[i][j] += left[i][k] * right[k][j];
		}

}

void mul(double M[3][3], double V[3], double out[3])
{
	for (int i = 0; i < 3; i++)
	{
		out[i] = 0.;
		for (int j = 0; j < 3; j++)
			out[i] += M[i][j] * V[j];
	}
}

double shadow_function_Hubaux(//shadow function
	double MJD,
	v3 Sat_ICSm,
	ae Sun_ICS_rad,
	double *Out_Distance_Sun_Earth_m,
	v3 &Out_Sun_ICS_m
)
{
	*Out_Distance_Sun_Earth_m = 0.;
	Out_Sun_ICS_m.x = 0.;
	Out_Sun_ICS_m.y = 0.;
	Out_Sun_ICS_m.z = 0.;

	double Out_ShadowFunction = 0.;
	//1=sun, penumbra, 0=umbra/shadow
	//this is a solution from Hubaux paper, 2012
	//"Sympletic integration of space debris motion considering several Earth's shadowing models"

	double Earth_radius_km = 6378.137;
	double Sun_radius_km = 696000.;

	v3 Sat_CCS_km = Sat_ICSm / 1000.;
	double r_km = norm(Sat_CCS_km);

	if (r_km > Earth_radius_km)
	{
		//Sun CCS
		double Distance_Sun_Earth_AU = Distance_to_the_SUN_AU(MJD);
		double AU_2_km = 149597870.7;//1 AU = km
		double Distance_Sun_Earth_km = Distance_Sun_Earth_AU * AU_2_km;

		v3 Sun_CCS_km = get_v3(Sun_ICS_rad, Distance_Sun_Earth_km);

		double valueIArad = IArad(Sun_CCS_km, Sat_CCS_km);//always positive, result is 0..PI

		double Sc = r_km * cos(valueIArad) + std::pow(r_km*r_km - Earth_radius_km * Earth_radius_km, 0.5);//eq5

		//in the paper there is a mistake, no need for 'r' in the denominator
		double Alpha_rad = atan((Sun_radius_km - Earth_radius_km) / Distance_Sun_Earth_km);
		double Beta_rad = atan((Sun_radius_km + Earth_radius_km) / Distance_Sun_Earth_km);

		double argument = std::pow(std::pow(r_km, 2.) - std::pow(Earth_radius_km*cos(Alpha_rad), 2.), 0.5);
		double Su = r_km * cos(valueIArad) + cos(Alpha_rad)*(argument + Earth_radius_km * sin(Alpha_rad));//eq6

		argument = std::pow(std::pow(r_km, 2.) - std::pow(Earth_radius_km*cos(Beta_rad), 2.), 0.5);
		double Sp = r_km * cos(valueIArad) + cos(Beta_rad)*(argument - Earth_radius_km * sin(Beta_rad));//eq7

		double Delta_h = Su - Sp;

		//penumbra function
		double delta = 8;//was 8
	   //  double vp=0.5*(1.+tanh(delta*PI2*Earth_radius_km*Sc/Delta_h));//original eq13
		double vp2 = 0.5*(1. + tanh(delta*Sc / Delta_h));//without constant factor of 'PI2*Earth_radius_km'

		if (Su <= 0.) Out_ShadowFunction = 0.; else
			if (Sp <= 0.) Out_ShadowFunction = vp2; else
				Out_ShadowFunction = 1.;

		*Out_Distance_Sun_Earth_m = Distance_Sun_Earth_AU * AU_2_km*1000.;
		Out_Sun_ICS_m= Sun_CCS_km*1000.;
	}

	return Out_ShadowFunction;
}
//------------------------------------------------------------------------------

double norm(v3 v)//returns vector length
{
	return std::pow(v.x*v.x + v.y*v.y + v.z*v.z, 0.5);
}

v3 get_v3(ae Vrad)//get unit vector at this orientation
{
	double Radius_XY = fabs(cos(Vrad.el));
	v3 Out;
	Out.x = Radius_XY * cos(Vrad.az);
	Out.y = Radius_XY * sin(Vrad.az);
	Out.z = sin(Vrad.el);
	return Out;
}

v3 get_v3(ae Vrad, double Length)//get unit vector at this orientation
{
	double Radius_XY = fabs(cos(Vrad.el));
	v3 Out;
	Out.x = Length * Radius_XY*cos(Vrad.az);
	Out.y = Length * Radius_XY*sin(Vrad.az);
	Out.z = Length * sin(Vrad.el);
	return Out;
}

double IArad(v3 A, v3 B)
{
	//the angle between two vectors is given by acos of the dot product of the two (normalised) vectors
	v3 Auv = normalize(A);
	v3 Buv = normalize(B);
	return acos(dot(Auv, Buv));
}

v3 normalize(v3 v)//returns unit vector
{
	v3 Out = v;
	double Length = norm(v);
	if (Length > 0.)
	{
		Out.x /= Length;
		Out.y /= Length;
		Out.z /= Length;
	}
	return Out;
}

double dot(v3 left, v3 right)
{
	double Sum = left.x*right.x + left.y*right.y + left.z*right.z;
	return Sum;
}

ae get_ae(v3 v)
{
	ae Vrad;
	Vrad.az = 0.;
	Vrad.el = 0.;

	v3 uv = normalize(v);
	Vrad.az = atan2(uv.y, uv.x);
	if (Vrad.az < 0.) Vrad.az += pi2;
	Vrad.el = pi / 2. - acos(uv.z);
	return Vrad;
}

void give_present_time(
	unsigned short *YearNow,
	unsigned short *MonthNow,
	unsigned short *DayNow,
	unsigned short *HourNow,
	unsigned short *MinNow,
	unsigned short *SecNow,
	unsigned short *mSecNow
)
{
	SYSTEMTIME st;
	//FILETIME ft;
	//LARGE_INTEGER li;
	GetSystemTime(&st);
	//SystemTimeToFileTime(&st, &ft);
	//li.LowPart = ft.dwLowDateTime;
	//li.HighPart = ft.dwHighDateTime;

	/*
	WORD wYear;
	WORD wMonth;
	WORD wDayOfWeek;
	WORD wDay;
	WORD wHour;
	WORD wMinute;
	WORD wSecond;
	WORD wMilliseconds;
	*/

	*YearNow = st.wYear;
	*MonthNow = st.wMonth;
	*DayNow = st.wDay;
	*HourNow = st.wHour;
	*MinNow = st.wMinute;
	*SecNow = st.wSecond;
	*mSecNow = st.wMilliseconds;
}

double present_MJD_full()
{
	unsigned short YearNow, MonthNow, DayNow;
	unsigned short HourNow, MinNow, SecNow, mSecNow;
	give_present_time(&YearNow, &MonthNow, &DayNow, &HourNow, &MinNow, &SecNow, &mSecNow);
	return MJD_full(YearNow, MonthNow, DayNow, HourNow, MinNow, double(SecNow));
}

double MJD_full(
	int year,//2006
	int month,
	int day,
	int hour,
	int minute,
	double second
)
{
	long double result = 0.;
	static int k[] = { 0,31,59,90,120,151,181,212,243,273,304,334 };
	if (year % 4 == 0 && month > 2) ++day;
	result =
		double((k[--month] + day + (year - 1972) * 365 + (year - 1969) / 4) + (long)41316.) +
		(double(hour * 3600 + minute * 60) + second) / 86400.;
	return result;
}

void pass_prediction()
{
	printf("\nPass prediction\n");
	//printf("\nPass prediction\n", RSO_app_mag);

	unsigned short YearNow, MonthNow, DayNow;
	unsigned short HourNow, MinNow, SecNow, mSecNow;
	give_present_time(&YearNow, &MonthNow, &DayNow, &HourNow, &MinNow, &SecNow, &mSecNow);

	printf("\nPredictions for date: ");
	print_date(YearNow, MonthNow, DayNow);
	printf("\n");

	double MJDnow = MJD_full(YearNow, MonthNow, DayNow, HourNow, MinNow, double(SecNow));


	char Directory[MAX_PATH];
	char TLE_filename[MAX_PATH];
	char out_filename[MAX_PATH];
	sprintf(Directory, "%s", "");
	sprintf(TLE_filename, "%s", "");
	sprintf(out_filename, "%s", "");

	sprintf(Directory, "%s", "C:\\Predictions\\tle");
	sprintf(TLE_filename, "%s\\%s", Directory, "20200129.3le");
	sprintf(out_filename, "%s\\%s", Directory, "pass.txt");

	//observer location
	v3 Observer_ITRFm;
	double Observer_elevation_m = 0.;
	ae Observer_LonLat_rad;

	/*
	//McDonald
	//https://ilrs.cddis.eosdis.nasa.gov/network/stations/active/MDOL_general.html
	Observer_ITRFm.x = -1330021.0;
	Observer_ITRFm.y = -5328401.8;
	Observer_ITRFm.z = 3236480.7;
	Observer_elevation_m = 2006.2210;
	Observer_LonLat_rad.az = 255.9848*deg2rad;
	Observer_LonLat_rad.el = 30.6802*deg2rad;
	*/

	//New Mexico Skies
	//https://www.nmskies.com/
	Observer_elevation_m =2225.04;
   // Observer_LonLat_rad.az=-(105.+31./60.+46.4/3600.)*deg2rad;//-105.5295556
	Observer_LonLat_rad.az = 254.4704444*deg2rad;
	Observer_LonLat_rad.el=(32.+54./60.+11./3600.)*deg2rad;//32.90305556
	Observer_ITRFm =Give_ITRFpoint_for_WGS84_ellipse_Lon_Lat_height(Observer_LonLat_rad, Observer_elevation_m);
	////////////////////////////////////////////////

	//read full TLE catalogue









}

void print_date(unsigned short Year, unsigned short Month, unsigned short Day)
{
	std::string MthName = "";
	switch (Month)
	{
	case 1: MthName = "Jan"; break;
	case 2: MthName = "Feb"; break;
	case 3: MthName = "Mar"; break;
	case 4: MthName = "Apr"; break;
	case 5: MthName = "May"; break;
	case 6: MthName = "June"; break;
	case 7: MthName = "July"; break;
	case 8: MthName = "Aug"; break;
	case 9: MthName = "Sept"; break;
	case 10: MthName = "Oct"; break;
	case 11: MthName = "Nov"; break;
	case 12: MthName = "Dec"; break;
	}

	printf("%s-%d-%d", MthName.c_str(), Day, Year);
}

v3 Give_ITRFpoint_for_WGS84_ellipse_Lon_Lat_height(//out: ITRF_m
	ae Ellipse_rad,
	double Ellipse_height_m
)
{
	//geodetic to cartesian
	//World Geodetic System
	//http://en.wikipedia.org/wiki/WGS84
	//this function gives the lon lat coordintates on the ellipse WGS84 of the site
	//http://www.ga.gov.au/geodesy/datums/calcs.jsp#trans
	v3 OutPos_m;
	OutPos_m.x = 0.;
	OutPos_m.y = 0.;
	OutPos_m.z = 0.;

	double Semi_major_axis = 6378160.00;//a
	double Inverse_flattening = 298.2500;//(1/f)
	double f = 1. / Inverse_flattening;
	double e2 = 2.*f - f * f;

	double Denominator = 0.;
	double sqrt_argument = 1. - e2 * sin(Ellipse_rad.el)*sin(Ellipse_rad.el);
	if (sqrt_argument >= 0.)
		Denominator = sqrt(sqrt_argument);

	double n = 0.;
	if (Denominator != 0.)
		n = Semi_major_axis / Denominator;

	OutPos_m.x = (n + Ellipse_height_m)*cos(Ellipse_rad.el)*cos(Ellipse_rad.az);
	OutPos_m.y = (n + Ellipse_height_m)*cos(Ellipse_rad.el)*sin(Ellipse_rad.az);
	OutPos_m.z = ((1. - e2)*n + Ellipse_height_m)*sin(Ellipse_rad.el);

	return OutPos_m;
}