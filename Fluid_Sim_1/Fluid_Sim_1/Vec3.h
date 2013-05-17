#ifndef VECTOR_MATHS_LIB
#define VECTOR_MATHS_LIB

class Vec3 {
public:
	double x[3];

	Vec3 ();
	Vec3 (double, double, double);

	Vec3 operator + (Vec3);
	Vec3 operator - (Vec3);
	Vec3 operator * (double);
	Vec3 operator += (Vec3);
	Vec3 operator - ();
	
	double sqr_mag ();
	double mag();
	Vec3 unit_vec();
};

double dot(Vec3 v1, Vec3 v2);

Vec3 cross(Vec3 v1, Vec3 v2);

typedef Vec3 RGB;

#endif