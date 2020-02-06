#pragma once



////////////////////////////////////////////////////////////////////////////////
//vector data type

//2-D
typedef struct { double x, y; } v2;
typedef struct { double az, el; } ae;
typedef struct { float az, el; } aef;
//3-D
typedef struct { double x, y, z; } v3;

typedef struct
{
	double v[3][3];//v[row][column], v[i=row][j=col]
} v33;

//user defined data types and operators
//vector operations

//v3
v3 operator + (const v3& V1, const v3& V2);
v3 operator - (const v3& V1, const v3& V2);
v3 operator - (const v3& V1);
v3 operator * (const v3 &V1, const v3 &V2);
v3 operator * (const v3 &V, const double &C);

v3 operator += (v3& V1, const v3& V2);
v3 operator -= (v3& V1, const v3& V2);
// v3 operator =() (v3& V, double& x, double& y, double& z);
v3 operator * (const double &C, v3 &V);
v3 operator *= (v3 &V, const double &C);
v3 operator /= (v3 &V, const double &C);

v3 operator * (const v33 &M, const v3& V);//matrix vector product

v33 operator * (const v33 &L, const v33 &R);//left matrix, right matrix product
v3 operator / (const v3 &V, const double &C);

//v2
v2 operator + (const v2 &V1, const v2 &V2);
v2 operator + (const v2 &V1, const double &C);
v2 operator - (const v2 &V1, const v2 &V2);
v2 operator - (const v2 &V1, const double &C);

////////////////////////////////////////////////////////////////////////////////
const COLORREF clAqua = RGB(0, 255, 255);
const COLORREF clBlack = RGB(0, 0, 0);
const COLORREF clBlue = RGB(0, 0, 255);
const COLORREF clCream = RGB(255, 251, 240);
const COLORREF clGray = RGB(128, 128, 128);
const COLORREF clFuchsia = RGB(255, 0, 255);
const COLORREF clGreen = RGB(0, 128, 0);
const COLORREF clLime = RGB(0, 255, 0);
const COLORREF clMaroon = RGB(128, 0, 0);
const COLORREF clNavy = RGB(0, 0, 128);
const COLORREF clOlive = RGB(128, 128, 0);
const COLORREF clPurple = RGB(255, 0, 255);
const COLORREF clRed = RGB(255, 0, 0);
const COLORREF clSilver = RGB(192, 192, 192);
const COLORREF clSilverLight = RGB(200, 200, 200);
const COLORREF clTeal = RGB(0, 128, 128);
const COLORREF clWhite = RGB(255, 255, 255);
const COLORREF clYellow = RGB(255, 255, 0);
const COLORREF clSun = RGB(255, 220, 100);


