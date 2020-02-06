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
#include <vector>

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
#define MAXBUF 200

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

int RSO_ephemeris_from_TLE(
	double EphSat_start_MJD,
	double EphSat_step_s,
	int EphSat_steps,
	std::string TLE_Line_1,
	std::string TLE_Line_2,
	//out
	v3 **EphSat_ICS_m//must be of size [EphSat_steps]
);
///////////////////////////////////////

///////////////////////////////////////
//Astronomical Algorithms
//source: http://www.naughter.com/aa.html
//example functions in AATest.cpp

double Distance_to_the_SUN_AU(double JD);

ae Give_Sun_topocentric_Az_El_rad(
	double MJD,
	v3 Observer_TCS_m,
	ae Observer_TCS_rad
);

double Give_Sun_topocentric_El_rad(
	double MJD,
	v3 Observer_TCS_m,
	ae Observer_TCS_rad
);

ae Sun_ICS_rad(double MJD);
///////////////////////////////////////

///////////////////////////////////////
//Date time MJD
void give_present_time_UTC(
	unsigned short *YearNow,
	unsigned short *MonthNow,
	unsigned short *DayNow,
	unsigned short *HourNow,
	unsigned short *MinNow,
	unsigned short *SecNow,
	unsigned short *mSecNow
);

double present_MJD_UTC_full();

double MJD_full(
	int year,//2006
	int month,
	int day,
	int hour,
	int minute,
	double second
);

void date_and_time_from_MJD(
	double MJD,
	int *year,
	int *month,
	int *day,
	int *hour,
	int *minute,
	double *sec
);

double MJD_full_DoY1_d(
	int year,//2006
	double doy1//day of year starting from 1, Jan 1=1doy
);

void date_from_year_and_doy1(//doy1=day of year starting from 1, Jan 1=1doy
	int year,
	int doy1,
	int *month,
	int *day
);

void print_date_time(double MJD);
///////////////////////////////////////

///////////////////////////////////////
//Geodesy and coordinate systems
v3 Give_ITRFpoint_for_WGS84_ellipse_Lon_Lat_height(//out: ITRF_m
	ae Ellipse_rad,
	double Ellipse_height_m
);

double GAST_rad(double MJD);

v3 eph_point(
	double MJD,
	//ephemeris
	double Eph_start_MJD,
	double Eph_step_s,
	int Eph_steps,
	v3 **Eph_m
);

ae topo_AzEl_rad(//topocentric
	v3 StaTCSm,
	ae StaTCSrad,
	v3 PointTCSm
);

v33 rx(double Angle_rad);
v33 ry(double Angle_rad);
v33 rz(double Angle_rad);
v3 ICS_to_TCS(double MJD, v3 ICS);//Inertial to Terrestrial

///////////////////////////////////////

///////////////////////////////////////
//I/O

inline bool file_exists(char *name);
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

void give_present_time_UTC(
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
	//Retrieves the current system date and time.The system time is expressed in Coordinated Universal Time(UTC).
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

double present_MJD_UTC_full()
{
	unsigned short YearNow, MonthNow, DayNow;
	unsigned short HourNow, MinNow, SecNow, mSecNow;
	give_present_time_UTC(&YearNow, &MonthNow, &DayNow, &HourNow, &MinNow, &SecNow, &mSecNow);
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
	FILE *FFout;
	errno_t err;
	char OutFilename[MAX_PATH];




	printf("\nNight-time pass prediction\n");
	//printf("\nPass prediction\n", RSO_app_mag);

	double MJD_UTC = present_MJD_UTC_full();

	printf("\nPredictions for UTC date: ");
	print_date_time(MJD_UTC);
	printf("\n");

	char Directory[MAX_PATH];
	char TLE_filename[MAX_PATH];
	char out_filename[MAX_PATH];
	sprintf_s(Directory, "%s", "");
	sprintf_s(TLE_filename, "%s", "");
	sprintf_s(out_filename, "%s", "");

	sprintf_s(Directory, "%s", "C:\\Predictions\\tle");
	sprintf_s(TLE_filename, "%s\\%s", Directory, "full.3le");//www.space-track.org
	sprintf_s(out_filename, "%s\\%s", Directory, "pass.txt");

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
	Observer_elevation_m = 2225.04;
	// Observer_LonLat_rad.az=-(105.+31./60.+46.4/3600.)*deg2rad;//-105.5295556
	Observer_LonLat_rad.az = 254.4704444*deg2rad;
	Observer_LonLat_rad.el = (32. + 54. / 60. + 11. / 3600.)*deg2rad;//32.90305556
	Observer_ITRFm = Give_ITRFpoint_for_WGS84_ellipse_Lon_Lat_height(Observer_LonLat_rad, Observer_elevation_m);
	////////////////////////////////////////////////


	//identify time-span of the current/upcoming night 
	double Sun_elevation_max_rad = -10.*deg2rad;
	double Night_start_MJD = -1.;
	double Night_stop_MJD = -1.;
	double MJD_step = 10. / 86400.;

	double Sun_elevation_rad = 0.;
	double MJD = MJD_UTC;
	double MJD_search_max = MJD_UTC + 1.;
	do
	{
		Sun_elevation_rad = Give_Sun_topocentric_El_rad(MJD, Observer_ITRFm, Observer_LonLat_rad);

		if (Sun_elevation_rad <= Sun_elevation_max_rad)
		{
			if (Night_start_MJD < 30000.)
				Night_start_MJD = MJD;

			Night_stop_MJD = MJD;
		}

		MJD += MJD_step;

	} while ((MJD < MJD_search_max) &&
		((Night_start_MJD < MJD_UTC) || (Sun_elevation_rad <= Sun_elevation_max_rad))
		);

	printf("\n");
	printf("Night start: ");
	print_date_time(Night_start_MJD);
	printf("\n");
	printf("Night stop: ");
	print_date_time(Night_stop_MJD);
	printf("\n");
	printf("Night duration: %.1lf h\n", (Night_stop_MJD - Night_start_MJD)*24.);

	/*
	//test sun topocentric angles	
	sprintf_s(OutFilename, "%s", "C:\\work\\test\\test1.txt");
	err = fopen_s(&FFout, OutFilename, "w");
	for (double MJD = Night_start_MJD; MJD <= Night_stop_MJD; MJD += MJD_step)
	{
		double dt = (MJD - Night_start_MJD)*86400.;		
		ae Sun_topo_rad = Give_Sun_topocentric_Az_El_rad(MJD, Observer_ITRFm, Observer_LonLat_rad);
		fprintf(FFout, "%.6lf\t%.6lf\t%.6lf\n", dt, Sun_topo_rad.az*rad2deg, Sun_topo_rad.el*rad2deg);
	}
	fclose(FFout);
	*/

	//Ephemeris parameters
	double Eph_step_s = 5.;
	double Eph_start_MJD = Night_start_MJD;
	int Eph_steps = 0;
	if (Eph_step_s > 0.)
		Eph_steps = int((Night_stop_MJD - Eph_start_MJD)*86400. / Eph_step_s) + 1;

	//Sun ephemeris
	v3 *Eph_Sun_ICS_m = new v3[Eph_steps];

	for (int iStep = 0; iStep < Eph_steps; iStep++)
	{
		double MJD = Eph_start_MJD + (double(iStep)*Eph_step_s) / 86400.;
		ae ICS_rad = Sun_ICS_rad(MJD);
		double Earth_Sun_m = Distance_to_the_SUN_AU(MJD)*AU_to_meter;
		Eph_Sun_ICS_m[iStep]= get_v3(ICS_rad, Earth_Sun_m);
	}


	/*
	//test sun ephemeris
	sprintf_s(OutFilename, "%s", "C:\\work\\test\\test2.txt");
	err = fopen_s(&FFout, OutFilename, "w");
	double test_step_MJD = 1. / 86400.;
	for (double MJD = Night_start_MJD; MJD < Night_stop_MJD; MJD += test_step_MJD)
	{
		double dt = (MJD - Night_start_MJD)*86400.;
		v3 ICS = eph_point(MJD, Eph_start_MJD, Eph_step_s, Eph_steps, &Eph_Sun_ICS_m);			
		v3 TCS=ICS_to_TCS(MJD, ICS);//Inertial to Terrestrial
		ae topo = topo_AzEl_rad(Observer_ITRFm, Observer_LonLat_rad, TCS);
		fprintf(FFout, "%.6lf\t%.6lf\t%.6lf\n", dt, topo.az*rad2deg, topo.el*rad2deg);
	}
	fclose(FFout);
	*/

	//read full TLE catalogue	
	std::vector <int> TLEcat_NORAD;
	std::vector <double> TLEcat_MJD;
	std::vector <std::string> TLEcat_Name;
	std::vector <std::string> TLEcat_L1;
	std::vector <std::string> TLEcat_L2;

	TLEcat_NORAD.resize(0);
	TLEcat_MJD.resize(0);
	TLEcat_Name.resize(0);
	TLEcat_L1.resize(0);
	TLEcat_L2.resize(0);


	if (file_exists(TLE_filename))
	{
		printf("Reading TLE catalogue\n");

		FILE *FFin;
		errno_t err;
		char buf[MAXBUF];
		char Line[MAX_PATH], L0[MAX_PATH], L1[MAX_PATH], L2[MAX_PATH];//3LE
		double record_MJD = 0.;
		char record_name[MAXBUF];
		char record_L1[MAXBUF];
		char record_L2[MAXBUF];
		
		//read TL3, 3 lines		
		sprintf_s(Line, "%s", "");
		sprintf_s(L1, "%s", "");
		sprintf_s(L2, "%s", "");

		err = fopen_s(&FFin, TLE_filename, "r");		

		int lines = 0;

		while (fgets(Line, sizeof(Line), FFin) != NULL)
		{
			// strip trailing '\n' if it exists
			int length = strlen(Line) - 1;
			if (Line[length] == '\n')
				Line[length] = 0;

			if (length > 1)
			{
				if (!memcmp(Line, "0 ", 2))
				{
					lines++;
					sprintf_s(L0, "%s", Line);
				}
				else
					if (!memcmp(Line, "1 ", 2))
					{
						lines++;
						sprintf_s(L1, "%s", Line);
					}
					else
						if (!memcmp(Line, "2 ", 2))
						{
							lines++;
							sprintf_s(L2, "%s", Line);

							int L0_length = strlen(L0);
							int L1_length = strlen(L1);
							int L2_length = strlen(L2);

							if ((L1_length >= 68) && (L2_length >= 68))
							{
								record_MJD = 0.;
								sprintf_s(record_name, "%s", "");
								sprintf_s(record_L1, "%s", "");
								sprintf_s(record_L2, "%s", "");
								sprintf_s(record_L1, "%s", L1);
								sprintf_s(record_L2, "%s", L2);

								if (L0_length > 2)
								{
									int name_length = L0_length - 2;
									memcpy(buf, L0 + 2, name_length);
									buf[name_length] = 0;
									sprintf_s(record_name, "%s", buf);
								}

								//get MJD
								//000000000011111111112222222222333333333344444444445555555555666666666
								//012345678901234567890123456789012345678901234567890123456789012345678
								//1 16908U 86061A   19305.80007987 -.00000099 +00000-0 -86868-5 0  9993
								int year = 0;//19-20 Epoch Year (Last two digits of year)
								double doy1 = 0.;//21-32 Epoch (Day of the year and fractional portion of the day)
								int Norad = 0;
								//sscanf(L1, "%*18c%2d%12lf", &year, &doy1);
								sscanf_s(L1, "%*2c%5d%*11c%2d%12lf", &Norad, &year, &doy1);

								if (year < 50.)
									year += 2000.;
								else
									year += 1900.;

								record_MJD = MJD_full_DoY1_d(year, doy1);

								if (record_MJD > 30000.)
								{
									TLEcat_NORAD.push_back(Norad);
									TLEcat_MJD.push_back(record_MJD);
									TLEcat_Name.push_back(std::string(record_name));
									TLEcat_L1.push_back(std::string(record_L1));
									TLEcat_L2.push_back(std::string(record_L2));
								}
							}
						}
			}
		}


		fclose(FFin);
	}
	
	/*
	//test
	printf("Objects in catalogue: %d\n", TLEcat_NORAD.size());
	for (int i = 0; i < TLEcat_NORAD.size(); i++)
	//for (int i = 0; i < 10; i++)
	{
		printf("%d ", TLEcat_NORAD[i]);
		print_date_time(TLEcat_MJD[i]);
		printf("\n");
		printf("%s\n", TLEcat_Name[i].c_str());
		printf("%s\n", TLEcat_L1[i].c_str());
		printf("%s\n", TLEcat_L2[i].c_str());
		printf("------------------------------------\n");
	}
	*/


	//get TLE for selected RSO

	//int NORAD = 44966;//Starlink-1088
	int NORAD = 44933;//Starlink-1144
	double TLE_MJD = 0.;
	std::string TLE_name = "";
	std::string TLE_L1="";
	std::string TLE_L2="";
	bool RSO_ok = false;
	for (int i = 0; i < TLEcat_NORAD.size(); i++)
		if (TLEcat_NORAD[i] == NORAD)
		{
			TLE_MJD = TLEcat_MJD[i];
			TLE_name = TLEcat_Name[i];
			TLE_L1 = TLEcat_L1[i];
			TLE_L2 = TLEcat_L2[i];
			RSO_ok = true;
		}

	if (RSO_ok)
	{
		printf("\nSelected RSO: %d ", NORAD);
		print_date_time(TLE_MJD);
		printf("\n");
		printf("%s\n", TLE_name.c_str());
		printf("%s\n", TLE_L1.c_str());
		printf("%s\n", TLE_L2.c_str());
		printf("------------------------------------\n");
	}
	else
	{
		printf("\nNo TLE for selected RSO.\n");
		printf("------------------------------------\n");
	}

	//ephemeris for selected satellite
	double EphSat_step_s = 1.;
	double EphSat_start_MJD = Night_start_MJD;
	int EphSat_steps = 0;
	if (EphSat_step_s > 0.)
		EphSat_steps = int((Night_stop_MJD - EphSat_start_MJD)*86400. / EphSat_step_s) + 1;
	
	v3 *EphSat_ICS_m = new v3[EphSat_steps];

	err = RSO_ephemeris_from_TLE(EphSat_start_MJD, EphSat_step_s, EphSat_steps, TLE_L1, TLE_L2, &EphSat_ICS_m);

	/*
	//test save sat ephemeris
	sprintf_s(OutFilename, "%s", "C:\\work\\test\\eph_sat.txt");
	err = fopen_s(&FFout, OutFilename, "w");
	for (int i = 0; i < EphSat_steps; i++)
	{
		double MJD = EphSat_start_MJD + (double(i)*EphSat_step_s) / 86400.;
		double dt = (MJD - Night_start_MJD)*86400.;
		v3 TCS = ICS_to_TCS(MJD, EphSat_ICS_m[i]);//Inertial to Terrestrial
		ae topo = topo_AzEl_rad(Observer_ITRFm, Observer_LonLat_rad, TCS);
		fprintf(FFout, "%.6lf\t%.6lf\t%.6lf\n", dt, topo.az*rad2deg, topo.el*rad2deg);
		//fprintf(FFout, "%d\t%lf\t%lf\t%lf\n", i, EphSat_ICS_m[i].x, EphSat_ICS_m[i].y, EphSat_ICS_m[i].z);
	}
	fclose(FFout);
	*/

	double pass_duration_min_s = 10.;
	double RSO_topo_elevation_min_rad = 10.*deg2rad;

	double pass_search_step_day = pass_duration_min_s / 2. / 86400.;
	double pass_fine_search_step_day = 1. / 86400.;

	std::vector <double> Pass_start_MJD;
	std::vector <double> Pass_stop_MJD;
	Pass_start_MJD.resize(0);
	Pass_stop_MJD.resize(0);

	double temp_pass_start_MJD = -1.;
	double temp_pass_stop_MJD = -1.;
	for (double MJD = Night_start_MJD; MJD < Night_stop_MJD; MJD += pass_search_step_day)
	{
		v3 Sat_ICS = eph_point(MJD, EphSat_start_MJD, EphSat_step_s, EphSat_steps, &EphSat_ICS_m);
		v3 Sat_TCS = ICS_to_TCS(MJD, Sat_ICS);//Inertial to Terrestrial
		ae Sat_topo_rad = topo_AzEl_rad(Observer_ITRFm, Observer_LonLat_rad, Sat_TCS);

		if (Sat_topo_rad.el >= RSO_topo_elevation_min_rad)
		{
			if (temp_pass_start_MJD < Night_start_MJD)
			temp_pass_start_MJD = MJD;

			temp_pass_stop_MJD = MJD;
		}
		else
		{
			//check current pass
			if ((temp_pass_start_MJD >= Night_start_MJD) && (temp_pass_stop_MJD >= temp_pass_start_MJD))
			{
				//fine search to determine exact start/stop epoch times
				double dd = 2.*pass_search_step_day;
				bool search = true;
				
				for (double dt = -dd; dt<=dd; dt+= pass_fine_search_step_day)
					if (search)
					{
						double MJD2 = temp_pass_start_MJD + dt;
						Sat_ICS = eph_point(MJD2, EphSat_start_MJD, EphSat_step_s, EphSat_steps, &EphSat_ICS_m);
						Sat_TCS = ICS_to_TCS(MJD2, Sat_ICS);//Inertial to Terrestrial
						Sat_topo_rad = topo_AzEl_rad(Observer_ITRFm, Observer_LonLat_rad, Sat_TCS);

						if (Sat_topo_rad.el >= RSO_topo_elevation_min_rad)
						{
							temp_pass_start_MJD = MJD2;
							search = false;
						}
					}

				search = true;
				for (double dt = -dd; dt <= dd; dt += pass_fine_search_step_day)
					if (search)
					{
						double MJD2 = temp_pass_stop_MJD + dt;
						Sat_ICS = eph_point(MJD2, EphSat_start_MJD, EphSat_step_s, EphSat_steps, &EphSat_ICS_m);
						Sat_TCS = ICS_to_TCS(MJD2, Sat_ICS);//Inertial to Terrestrial
						Sat_topo_rad = topo_AzEl_rad(Observer_ITRFm, Observer_LonLat_rad, Sat_TCS);

						if (Sat_topo_rad.el >= RSO_topo_elevation_min_rad)
						{
							temp_pass_stop_MJD = MJD2;							
						}
						else
						{
							search = false;
						}
					}

				Pass_start_MJD.push_back(temp_pass_start_MJD);
				Pass_stop_MJD.push_back(temp_pass_stop_MJD);
			}

			temp_pass_start_MJD = -1.;
			temp_pass_stop_MJD = -1.;
		}

//		v3 Sun_ICS = eph_point(MJD, Eph_start_MJD, Eph_step_s, Eph_steps, &Eph_Sun_ICS_m);
//		v3 TCS = ICS_to_TCS(MJD, ICS);//Inertial to Terrestrial
//		ae topo = topo_AzEl_rad(Observer_ITRFm, Observer_LonLat_rad, TCS);
//		fprintf(FFout, "%.6lf\t%.6lf\t%.6lf\n", dt, topo.az*rad2deg, topo.el*rad2deg);
	}

	//list passes
	if (Pass_start_MJD.size()>0)
	{ 
		for (int iPass = 0; iPass < Pass_start_MJD.size(); iPass++)
		{
			double pass_duration_s = (Pass_stop_MJD[iPass] - Pass_start_MJD[iPass])*86400.;

			printf("\nPass %d, start ", iPass+1);
			print_date_time(Pass_start_MJD[iPass]);
			printf(", stop ");
			print_date_time(Pass_stop_MJD[iPass]);
			printf(", duration %6.1lf s\n", pass_duration_s);
		}
	}
	





	if (Eph_steps > 0)
	{
		delete[] Eph_Sun_ICS_m;
	}

	if (EphSat_steps > 0)
		delete[] EphSat_ICS_m;


	TLEcat_NORAD.resize(0);
	TLEcat_MJD.resize(0);
	TLEcat_Name.resize(0);
	TLEcat_L1.resize(0);
	TLEcat_L2.resize(0);
	Pass_start_MJD.resize(0);
	Pass_stop_MJD.resize(0);
}

void print_date_time(double MJD)
{
	int year, month, day, hour, minute;
	double sec;

	date_and_time_from_MJD(MJD,&year, &month, &day, &hour, &minute, &sec);
	
	std::string MthName = "";
	switch (month)
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

	printf("%s-%d-%d %d:%02d:%04.1lf", MthName.c_str(), day, year, hour, minute, sec);
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

void date_and_time_from_MJD(
	double MJD,
	int *year,
	int *month,
	int *day,
	int *hour,
	int *minute,
	double *sec
)
{
	long int st = 678881L;
	long int st1 = 146097L;
	long int a = 4L * (int(MJD) + st) + 3L;
	long int b = 0.;
	if (st1 != 0.)
		b = a / st1;
	else
		b = 0.;//error
	long int c = a - b * st1;
	long int f = 4 * (c / 4) + 3;
	c = f / 1461;
	f = (f - c * 1461 + 4) / 4;
	long int e = f * 5 - 3;
	*year = (int)(b * 100 + c);
	*day = (int)((e - (e / 153) * 153 + 5) / 5);
	*month = (int)(e / 153 + 3);
	if (*month > 12)
	{
		*month -= 12;
		++*year;
	}
	double time = (MJD - int(MJD))*86400.;
	*hour = int(time) / 3600;
	*minute = (int(time) - *hour * 3600) / 60;
	*sec = time - double(*hour * 3600) - double(*minute * 60);
}

ae Give_Sun_topocentric_Az_El_rad(
	double MJD,
	v3 Observer_TCS_m,
	ae Observer_TCS_rad	
)
{
	ae ICS_rad = Sun_ICS_rad(MJD);
	double Earth_Sun_m = Distance_to_the_SUN_AU(MJD)*AU_to_meter;
	v3 ICS_m = get_v3(ICS_rad, Earth_Sun_m);
	v3 TCS_m = ICS_to_TCS(MJD, ICS_m);//Inertial to Terrestrial
	ae Topo_rad = topo_AzEl_rad(Observer_TCS_m, Observer_TCS_rad, TCS_m);
	return Topo_rad;	
}

double Give_Sun_topocentric_El_rad(
	double MJD,
	v3 Observer_TCS_m,
	ae Observer_TCS_rad
)
{
	ae Sun_Topo_rad = Give_Sun_topocentric_Az_El_rad(MJD, Observer_TCS_m, Observer_TCS_rad);
	return Sun_Topo_rad.el;
}

v3 eph_point(
	double MJD,
	//ephemeris
	double Eph_start_MJD,
	double Eph_step_s,
	int Eph_steps,
	v3 **Eph_m
)
{
	v3 Out_m;
	Out_m.x = 0.;
	Out_m.y = 0.;
	Out_m.z = 0.;

	if ((MJD >= Eph_start_MJD) && (Eph_step_s > 0.))
	{
		double sec = (MJD - Eph_start_MJD)*86400.;
		int i0 = int(sec / Eph_step_s);
		int i1 = i0 + 1;
		double dt = (sec - double(i0)*Eph_step_s) / Eph_step_s;

		if ((i0 >= 0) && (i0 < Eph_steps) &&
			(i1 >= 0) && (i1 < Eph_steps) &&
			(i1 >= i0))
		{
			v3 dd = (*Eph_m)[i1] - (*Eph_m)[i0];
			Out_m = dd * dt + (*Eph_m)[i0];
		}
	}
	return Out_m;
}

ae topo_AzEl_rad(//topocentric
	v3 StaTCSm,
	ae StaTCSrad,
	v3 PointTCSm
)
{
	//Azimuth - from North (true North) CW towards East
	//North Az=0, East Az=90, South Az=180, West Az=270
	v3 toptrt = PointTCSm - StaTCSm;
	v33 m3 = rz(StaTCSrad.az);
	v33 m2 = ry(pi / 2. - StaTCSrad.el);
	v3 out = m2 * m3*toptrt;
	out.x *= -1.;//because azimuth is counted from North vector towards East: https://commons.wikimedia.org/wiki/File:ECEF_ENU_Longitude_Latitude_relationships.svg
	ae Out_rad = get_ae(out);
	return Out_rad;
}

v33 rx(double Angle_rad)
{
	//Seeber
	double s = sin(Angle_rad);
	double c = cos(Angle_rad);
	v33 out;
	out.v[0][0] = 1.0;  out.v[0][1] = 0.0;  out.v[0][2] = 0.0;
	out.v[1][0] = 0.0;  out.v[1][1] = c;  out.v[1][2] = s;
	out.v[2][0] = 0.0;  out.v[2][1] = -s;  out.v[2][2] = c;
	return out;
}

v33 ry(double Angle_rad)
{
	//Seeber
	double s = sin(Angle_rad);
	double c = cos(Angle_rad);
	v33 out;
	out.v[0][0] = c;  out.v[0][1] = 0.0;  out.v[0][2] = -s;
	out.v[1][0] = 0.0;  out.v[1][1] = 1.0;  out.v[1][2] = 0.0;
	out.v[2][0] = s;  out.v[2][1] = 0.0;  out.v[2][2] = c;
	return out;
}

v33 rz(double Angle_rad)
{
	//Seeber
	double s = sin(Angle_rad);
	double c = cos(Angle_rad);
	v33 out;
	out.v[0][0] = c;  out.v[0][1] = s;  out.v[0][2] = 0.0;
	out.v[1][0] = -s;  out.v[1][1] = c;  out.v[1][2] = 0.0;
	out.v[2][0] = 0.0;  out.v[2][1] = 0.0;  out.v[2][2] = 1.0;
	return out;
}


v3 ICS_to_TCS(double MJD, v3 ICS)//Inertial to Terrestrial
{
	//no polar motion
	//   Matrix Rot_CCStoTCS(R_y(-PM_x_rad)*R_x(-PM_y_rad)*R_z(GAST_rad));
	//-> polar motion changes station coordinates by ~10 meters; polar motion is below 0.5 arcsec
	//polar motion will affect orientation of staton vector on Ajisai below 0.4 mdeg.
	double gast_rad = GAST_rad(MJD);
	v3 TCS = rz(gast_rad) * ICS;
	return TCS;
}

ae Sun_ICS_rad(double MJD)
{
	double JD = MJD + 2400000.5;
	double Epsilon = CAANutation::TrueObliquityOfEcliptic(JD);

	//sun apparent equatorial coordinates
	double Apparent_ecliptic_Lon_deg = CAASun::ApparentEclipticLongitude(JD, false);
	double Apparent_ecliptic_Lat_deg = CAASun::ApparentEclipticLatitude(JD, false);
	CAA2DCoordinate Sun1 = CAACoordinateTransformation::Ecliptic2Equatorial(Apparent_ecliptic_Lon_deg, Apparent_ecliptic_Lat_deg, Epsilon);

	ae CCSrad;
	CCSrad.az = CAACoordinateTransformation::HoursToRadians(Sun1.X);
	CCSrad.el = CAACoordinateTransformation::DegreesToRadians(Sun1.Y);

	return CCSrad;
}
//------------------------------------------------------------------------------

inline bool file_exists(char *name)
{
	struct stat buffer;
	return (stat(name, &buffer) == 0);
}

double MJD_full_DoY1_d(
	int year,//2006
	double doy1//day of year starting from 1, Jan 1=1doy
)
{
	int month, day, Hour, Minute;
	date_from_year_and_doy1(year, int(doy1), &month, &day);
	double SOD = (doy1 - double(int(doy1)))*86400.;
	Hour = int(SOD) / 3600;
	Minute = int(SOD) - Hour * 3600;
	double Second = SOD - double(Hour * 3600) - double(Minute * 60);
	double Full_MJD = MJD_full(year, month, day, Hour, Minute, Second);
	return Full_MJD;
}

void date_from_year_and_doy1(//doy1=day of year starting from 1, Jan 1=1doy
	int year,
	int doy1,
	int *month,
	int *day
)
{
	int d[12] = { 31,59,90,120,151,181,212,243,273,304,334,365 };
	//is leap year? XXI century
	if (year % 4 == 0)
	{
		for (unsigned int i = 1; i < 12; i++)
			d[i] = d[i] + 1;
	}

	*month = 1;
	for (unsigned int i = 0; i < 12; i++)
	{
		if (doy1 > d[i])
			*month = i + 2;
	}
	if (*month > 1)
		*day = doy1 - d[*month - 2];
	else
		*day = doy1;
}

int RSO_ephemeris_from_TLE(
	double EphSat_start_MJD,
	double EphSat_step_s,
	int EphSat_steps,
	std::string TLE_Line_1,
	std::string TLE_Line_2,
	//out
	v3 **EphSat_ICS_m//must be of size [EphSat_steps]
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

	//reset output
	for (int i = 0; i < EphSat_steps; i++)
	{
		(*EphSat_ICS_m)[i].x = 0.;
		(*EphSat_ICS_m)[i].y = 0.;
		(*EphSat_ICS_m)[i].z = 0.;
	}

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
	sprintf_s(longstr1, "%s", TLE_Line_1.c_str());
	sprintf_s(longstr2, "%s", TLE_Line_2.c_str());

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

	for (int iStep = 0; iStep < EphSat_steps; iStep++)
	{
		double MJD = EphSat_start_MJD + (double(iStep)*EphSat_step_s) / 86400.;
		double tsince = (MJD + 2400000.5 - ReferenceEpochJD)*1440.;//[minute]
		SGP4Funcs::sgp4(satrec, tsince, ro, vo);

		(*EphSat_ICS_m)[iStep].x = ro[0] * 1.e3;//km->m
		(*EphSat_ICS_m)[iStep].y = ro[1] * 1.e3;//km->m
		(*EphSat_ICS_m)[iStep].z = ro[2] * 1.e3;//km->m
	}

	return ReturnError;
}