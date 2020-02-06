#pragma once



////////////////////////////////////////////////////////////////////////////////
//vector data type

//2-D
typedef struct { double x, y; } v2;
typedef struct { double az, el; } ae;
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


