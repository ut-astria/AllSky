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

#include <math.h>
#include <errno.h>

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

int Round_value(long double Value);
double GAST_rad(double MJD);
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



//endd

//main function starts here
int main(int argc, const char *argv[])
{
	printf("\nAll-Sky\n");


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

	//Use shadow function
	float Shadow = 0.;

	float RSO_app_mag = RSO_apparent_magnitude(MJD_UTC, Sun_ICS_m, Sat_ICSm, Observer_ICS_m);

	printf("\nRSO apparent magnitude = %.1lf\n", RSO_app_mag);

	//SGP4Funcs::invjday(JD, 0., year, month, day, hour, minute, second);

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