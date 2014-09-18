#ifndef VECTOR_MATHS_LIB
#define VECTOR_MATHS_LIB

class Vec3 {
public:
	float x[3];

	Vec3 ();
	Vec3 (float, float, float);

	Vec3 operator + (Vec3);
	Vec3 operator - (Vec3);
	Vec3 operator * (float);
	Vec3 operator += (Vec3);
	Vec3 operator - ();
	
	float sqr_mag ();
	float mag();
	Vec3 unit_vec();
};

float dot(Vec3 v1, Vec3 v2);

Vec3 cross(Vec3 v1, Vec3 v2);

typedef Vec3 RGB;

#endif