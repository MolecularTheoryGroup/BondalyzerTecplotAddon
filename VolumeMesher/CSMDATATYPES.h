#ifndef CSMDATATYPES_H_
#define CSMDATATYPES_H_

#define GTADEBUG 0

#define   DEG2RAD 0.017453292519943295769236907684886127134428718885417

typedef double CSM_ResultType_t;

enum CSM_CompDir_e{
	LessThan = 1,
	GreaterThan
};

/*
 *	Quick declaration so that the CSM_Vec3_s constructor from
 *	a vec4 can work
 */
struct CSM_Vec4_s;

struct CSM_Vec3_s{
	double X, Y, Z;

	CSM_Vec3_s(double x = 0, double y = 0, double z = 0);
	CSM_Vec3_s(CSM_Vec4_s const &a);
	CSM_Vec3_s(CSM_Vec3_s const &a);

	CSM_Vec3_s & operator=(CSM_Vec3_s const &a);
	CSM_Vec3_s & operator=(double a);
	const bool		operator==(CSM_Vec3_s const &a) const;
	const bool		operator!=(CSM_Vec3_s const &a) const;
	double &	operator[](unsigned int i);
	double const at(unsigned int i);

	CSM_Vec3_s const operator-();
	CSM_Vec3_s & operator-=(CSM_Vec3_s const &a);
	CSM_Vec3_s & operator+=(CSM_Vec3_s const &a);
	CSM_Vec3_s & operator*=(double const &a);
	CSM_Vec3_s & operator*=(CSM_Vec3_s const &a);
	const CSM_Vec3_s	operator-(CSM_Vec3_s const &a) const;
	const CSM_Vec3_s	operator+(CSM_Vec3_s const &a) const;
	const CSM_Vec3_s	operator*(double a) const;
	const CSM_Vec3_s	operator*(CSM_Vec3_s const &a) const;

	const double		Dot(CSM_Vec3_s const &a) const;
	const CSM_Vec3_s	Cross(CSM_Vec3_s const &a) const;
	const CSM_Vec3_s	Normalize() const;
	const double		Magnitude() const;
	const double		MagSqr() const;
	const CSM_Vec3_s	RotationAngles(CSM_Vec3_s const &a) const;
};

const CSM_Vec3_s		Normalize(CSM_Vec3_s const &a);
const double			CSM_Distance(CSM_Vec3_s const &a, CSM_Vec3_s const &b);


struct CSM_Mat3_s{
	CSM_Vec3_s R1, R2, R3;

	CSM_Mat3_s(CSM_Vec3_s r1 = 0, CSM_Vec3_s r2 = 0, CSM_Vec3_s r3 = 0);
	CSM_Mat3_s(CSM_Mat3_s const & a);

	CSM_Vec3_s & operator[](unsigned int i);
	CSM_Vec3_s & at(unsigned int i);

	CSM_Mat3_s & operator*=(CSM_Mat3_s &a);
	CSM_Vec3_s const operator*(CSM_Vec3_s const &a) const;
	CSM_Mat3_s const operator*(CSM_Mat3_s &a) const;
};

struct CSM_Vec4_s{
	double X, Y, Z, W;

	CSM_Vec4_s(double x = 0, double y = 0, double z = 0, double w = 0);
	CSM_Vec4_s(CSM_Vec3_s const &a, double w = 0);

	CSM_Vec4_s & operator=(CSM_Vec4_s const &a);
	CSM_Vec4_s & operator=(double a);
	const bool		operator==(CSM_Vec4_s const &a) const;
	const bool		operator!=(CSM_Vec4_s const &a) const;
	double &	operator[](unsigned int i);

	CSM_Vec4_s const operator-();
	CSM_Vec4_s & operator-=(CSM_Vec4_s const &a);
	CSM_Vec4_s & operator+=(CSM_Vec4_s const &a);
	CSM_Vec4_s & operator*=(double const &a);
	CSM_Vec4_s & operator*=(CSM_Vec4_s const &a);
	const CSM_Vec4_s	operator-(CSM_Vec4_s const &a) const;
	const CSM_Vec4_s	operator+(CSM_Vec4_s const &a) const;
	const CSM_Vec4_s	operator*(double a) const;
	const CSM_Vec4_s	operator*(CSM_Vec4_s const &a) const;

	const double		Dot(CSM_Vec4_s const &a) const;
	const CSM_Vec4_s	Cross(CSM_Vec4_s const &a) const;
	const CSM_Vec4_s	Normalize() const;
	const double		Magnitude() const;
	const CSM_Vec4_s	RotationAngles(CSM_Vec4_s const &a) const;
};

const CSM_Vec4_s		Normalize(CSM_Vec4_s const &a);
const double			CSM_Distance(CSM_Vec4_s const &a, CSM_Vec4_s const &b);


struct CSM_Mat4_s{
	CSM_Vec4_s R1, R2, R3, R4;

	CSM_Mat4_s(CSM_Vec4_s r1 = 0, CSM_Vec4_s r2 = 0, CSM_Vec4_s r3 = 0, CSM_Vec4_s r4 = 0);
	CSM_Mat4_s(CSM_Mat4_s const & a);

	CSM_Vec4_s & operator[](unsigned int i);
	CSM_Vec4_s & at(unsigned int i);

	CSM_Mat4_s & operator*=(CSM_Mat4_s &a);
	CSM_Vec4_s const operator*(CSM_Vec4_s const &a) const;
	CSM_Mat4_s const operator*(CSM_Mat4_s &a) const;
};

const CSM_Mat4_s		CSM_Rotate(double Angle, CSM_Vec3_s const &a);



#endif
