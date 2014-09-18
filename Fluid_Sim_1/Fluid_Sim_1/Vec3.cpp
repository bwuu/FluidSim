#include "stdafx.h"
#include "Vec3.h"

Vec3::Vec3() {}

Vec3 Vec3::operator - ()
{
	Vec3 res;
	res.x[0] = -x[0];
	res.x[1] = -x[1];
	res.x[2] = -x[2];
	return res;
}

Vec3::Vec3 (float a, float b, float c) 
{
	x[0] = a;
	x[1] = b;
	x[2] = c;
}

Vec3 Vec3::operator + (Vec3 arg)
{
	Vec3 res;
	res.x[0] = x[0] + arg.x[0];
	res.x[1] = x[1] + arg.x[1];
	res.x[2] = x[2] + arg.x[2];
	return (res);
}


Vec3 Vec3::operator - (Vec3 arg)
{
	Vec3 res;
	res.x[0] = x[0] - arg.x[0];
	res.x[1] = x[1] - arg.x[1];
	res.x[2] = x[2] - arg.x[2];
	return (res);
}
	

Vec3 Vec3::operator * (float arg)
{
	Vec3 res;
	res.x[0] = x[0] * arg;
	res.x[1] = x[1] * arg;
	res.x[2] = x[2] * arg;
	return (res);
}

Vec3 Vec3::operator += (Vec3 arg)
{
	x[0] += arg.x[0];
	x[1] += arg.x[1];
	x[2] += arg.x[2];
	return *this;
}

float Vec3::sqr_mag ()
{
	return (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]); 
}

float Vec3::mag() 
{
	return sqrt(sqr_mag());
}

Vec3 Vec3::unit_vec () 
{
	Vec3 res;
	const float norm = mag();
	assert(norm > 0.0);
	res.x[0] = x[0] / norm;
	res.x[1] = x[1] / norm;
	res.x[2] = x[2] / norm;
	return (res);
}

float dot(Vec3 v1, Vec3 v2)
{
	return (v1.x[0] * v2.x[0] + v1.x[1] * v2.x[1] + v1.x[2] * v2.x[2]);
}

Vec3 cross(Vec3 v1, Vec3 v2)
{
	Vec3 res;
	res.x[0] = v1.x[1] * v2.x[2] - v1.x[2] * v2.x[1];
	res.x[1] = v1.x[2] * v2.x[0] - v1.x[0] * v2.x[2];
	res.x[2] = v1.x[0] * v2.x[1] - v1.x[1] * v2.x[0];
	return (res);
}