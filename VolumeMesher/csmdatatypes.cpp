
#include <cmath>
#include <vector>

#include "CSMDATATYPES.h"


/*
	Begin methods for CSM_Vec3_s
*/
CSM_Vec3_s::CSM_Vec3_s(double x, double y, double z){
	X = x;
	Y = y;
	Z = z;
}
CSM_Vec3_s::CSM_Vec3_s(CSM_Vec4_s const &a){
	X = a.X;
	Y = a.Y;
	Z = a.Z;
}
CSM_Vec3_s::CSM_Vec3_s(CSM_Vec3_s const &a){
	X = a.X;
	Y = a.Y;
	Z = a.Z;
}

CSM_Vec3_s & CSM_Vec3_s::operator=(CSM_Vec3_s const &a){
	X = a.X;
	Y = a.Y;
	Z = a.Z;
	return *this;
}
CSM_Vec3_s & CSM_Vec3_s::operator=(double a){
	X = Y = Z = a;
	return *this;
}
bool const CSM_Vec3_s::operator==(CSM_Vec3_s const &a) const{
	return (
		X == a.X &&
		Y == a.Y &&
		Z == a.Z);
}
bool const CSM_Vec3_s::operator!=(CSM_Vec3_s const &a) const{
	return !(*this == a);
}
double &	CSM_Vec3_s::operator[](unsigned int i){
	switch (i){
		case 0:
			return X;
			break;
		case 1:
			return Y;
			break;
		case 2:
			return Z;
			break;
		default:
			throw - 1;
	}
}

double const CSM_Vec3_s::at(unsigned int i)
{
	switch (i){
		case 0:
			return X;
			break;
		case 1:
			return Y;
			break;
		case 2:
			return Z;
			break;
		default:
			throw - 1;
	}
}

CSM_Vec3_s const CSM_Vec3_s::operator-(){
	return CSM_Vec3_s(-X, -Y, -Z);
}
CSM_Vec3_s & CSM_Vec3_s::operator-=(CSM_Vec3_s const &a){
	X -= a.X;
	Y -= a.Y;
	Z -= a.Z;

	return *this;
}
CSM_Vec3_s & CSM_Vec3_s::operator+=(CSM_Vec3_s const &a){
		X += a.X;
		Y += a.Y;
		Z += a.Z;

		return *this;
}
CSM_Vec3_s & CSM_Vec3_s::operator*=(double const &a){
	X *= a;
	Y *= a;
	Z *= a;

	return *this;
}
CSM_Vec3_s & CSM_Vec3_s::operator*=(CSM_Vec3_s const &a){
	X *= a.X;
	Y *= a.Y;
	Z *= a.Z;

	return *this;
}
CSM_Vec3_s const CSM_Vec3_s::operator-(CSM_Vec3_s const &a) const{
	return CSM_Vec3_s(*this) -= a;
}
CSM_Vec3_s const CSM_Vec3_s::operator+(CSM_Vec3_s const &a) const{
	return CSM_Vec3_s(*this) += a;
}
const CSM_Vec3_s	CSM_Vec3_s::operator*(double a) const{
	return CSM_Vec3_s(*this) *= a;
}
const CSM_Vec3_s	CSM_Vec3_s::operator*(CSM_Vec3_s const &a) const{
	return CSM_Vec3_s(*this) *= a;
}

double const CSM_Vec3_s::Dot(CSM_Vec3_s const &a) const{
	return X*a.X + Y*a.Y + Z*a.Z;
}
CSM_Vec3_s const CSM_Vec3_s::Cross(CSM_Vec3_s const &a) const{
	CSM_Vec3_s Tmp;
	Tmp.X = -Z*a.Y + Y*a.Z;
	Tmp.Y = Z*a.X - X*a.Z;
	Tmp.Z = -Y*a.X + X*a.Y;
	return Tmp;
}
const CSM_Vec3_s	CSM_Vec3_s::Normalize() const{
	CSM_Vec3_s a(*this);
	double Mag = a.Magnitude();
	if (Mag > 0){
		a.X /= Mag;
		a.Y /= Mag;
		a.Z /= Mag;
	}
	return a;
}
double const CSM_Vec3_s::Magnitude() const
{
	return sqrt(X*X + Y*Y + Z*Z);
}
double const CSM_Vec3_s::MagSqr() const
{
	return X*X + Y*Y + Z*Z;
}
const CSM_Vec3_s	CSM_Vec3_s::RotationAngles(CSM_Vec3_s const &a) const{
	CSM_Vec3_s RotAngles;
	for (int i = 0; i < 3; ++i){
		CSM_Vec3_s Tmp1 = *this;
		CSM_Vec3_s Tmp2 = a;
		Tmp1[i] = Tmp2[i] = 0;
		RotAngles[i] = acos(Tmp1.Dot(Tmp2) / (Tmp1.Magnitude() * Tmp2.Magnitude()));
	}
	return RotAngles;
}
/*
	End methods for CSM_Vec3_s
*/

const CSM_Vec3_s		Normalize(CSM_Vec3_s const &a){
	return CSM_Vec3_s(a).Normalize();
}
double const CSM_Distance(CSM_Vec3_s const &a, CSM_Vec3_s const &b){
	return CSM_Vec3_s(a - b).Magnitude();
}


/*
	Begin methods for CSM_Mat3_s
*/

CSM_Mat3_s::CSM_Mat3_s(CSM_Vec3_s r1, CSM_Vec3_s r2, CSM_Vec3_s r3){
	R1 = r1;
	R2 = r2;
	R3 = r3;
}
CSM_Mat3_s::CSM_Mat3_s(CSM_Mat3_s const & a){
	R1 = a.R1;
	R2 = a.R2;
	R3 = a.R3;
}

CSM_Vec3_s & CSM_Mat3_s::operator[](unsigned int i){
	switch (i){
		case 0:
			return R1;
			break;
		case 1:
			return R2;
			break;
		case 2:
			return R3;
			break;
		default:
			throw - 1;
	}
}
CSM_Vec3_s & CSM_Mat3_s::at(unsigned int i){
	return this->operator [](i);
}

CSM_Mat3_s & CSM_Mat3_s::operator*=(CSM_Mat3_s &a){
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
				this->at(i)[j] += this->at(i)[k] * a[j][k];

	return *this;
}
CSM_Vec3_s const CSM_Mat3_s::operator*(CSM_Vec3_s const &a) const{
	CSM_Vec3_s Pdt;
	CSM_Mat3_s aa(*this);
	for (int i = 0; i < 3; ++i)
		Pdt[i] = aa[i].Dot(a);
	return Pdt;
}
CSM_Mat3_s const CSM_Mat3_s::operator*(CSM_Mat3_s &a) const{
	return CSM_Mat3_s(*this) *= a;
}

/*
Begin methods for CSM_Vec4_s
*/
CSM_Vec4_s::CSM_Vec4_s(double x, double y, double z, double w){
	X = x;
	Y = y;
	Z = z;
	W = w;
}
CSM_Vec4_s::CSM_Vec4_s(CSM_Vec3_s const &a, double w){
	X = a.X;
	Y = a.Y;
	Z = a.Z;
	W = w;
}

CSM_Vec4_s & CSM_Vec4_s::operator=(CSM_Vec4_s const &a){
	X = a.X;
	Y = a.Y;
	Z = a.Z;
	W = a.W;
	return *this;
}
CSM_Vec4_s & CSM_Vec4_s::operator=(double a){
	X = Y = Z = W = a;
	return *this;
}
bool const CSM_Vec4_s::operator==(CSM_Vec4_s const &a) const{
	return (
		X == a.X &&
		Y == a.Y &&
		Z == a.Z &&
		W == a.W);
}
bool const CSM_Vec4_s::operator!=(CSM_Vec4_s const &a) const{
	return !(*this == a);
}
double &	CSM_Vec4_s::operator[](unsigned int i){
	switch (i){
		case 0:
			return X;
			break;
		case 1:
			return Y;
			break;
		case 2:
			return Z;
			break;
		case 3:
			return W;
			break;
		default:
			throw - 1;
	}
}

CSM_Vec4_s const CSM_Vec4_s::operator-(){
	return CSM_Vec4_s(-X, -Y, -Z, -W);
}
CSM_Vec4_s & CSM_Vec4_s::operator-=(CSM_Vec4_s const &a){
	X -= a.X;
	Y -= a.Y;
	Z -= a.Z;
	W -= a.W;

	return *this;
}
CSM_Vec4_s & CSM_Vec4_s::operator+=(CSM_Vec4_s const &a){
	X += a.X;
	Y += a.Y;
	Z += a.Z;
	W += a.W;

	return *this;
}
CSM_Vec4_s & CSM_Vec4_s::operator*=(double const &a){
	X *= a;
	Y *= a;
	Z *= a;
	W *= a;

	return *this;
}
CSM_Vec4_s & CSM_Vec4_s::operator*=(CSM_Vec4_s const &a){
	X *= a.X;
	Y *= a.Y;
	Z *= a.Z;
	W *= a.W;

	return *this;
}
CSM_Vec4_s const CSM_Vec4_s::operator-(CSM_Vec4_s const &a) const{
	return CSM_Vec4_s(*this) -= a;
}
CSM_Vec4_s const CSM_Vec4_s::operator+(CSM_Vec4_s const &a) const{
	return CSM_Vec4_s(*this) += a;
}
const CSM_Vec4_s	CSM_Vec4_s::operator*(double a) const{
	return CSM_Vec4_s(*this) *= a;
}
const CSM_Vec4_s	CSM_Vec4_s::operator*(CSM_Vec4_s const &a) const{
	return CSM_Vec4_s(*this) *= a;
}

double const CSM_Vec4_s::Dot(CSM_Vec4_s const &a) const{
	return X*a.X + Y*a.Y + Z*a.Z + W*a.W;
}
// CSM_Vec4_s CSM_Vec4_s::Cross(CSM_Vec4_s const &a){
// 	CSM_Vec4_s Tmp;
// 	Tmp.X = -Z*a.Y + Y*a.Z;
// 	Tmp.Y = Z*a.X - X*a.Z;
// 	Tmp.Z = -Y*a.X + X*a.Y;
// 	return Tmp;
// }
const CSM_Vec4_s	CSM_Vec4_s::Normalize() const{
	double Mag = this->Magnitude();
	CSM_Vec4_s a = *this;
	if (Mag > 0){
		a.X /= Mag;
		a.Y /= Mag;
		a.Z /= Mag;
		a.W /= Mag;
	}
	return a;
}
double const CSM_Vec4_s::Magnitude() const
{
	return sqrt(X*X + Y*Y + Z*Z + W*W);
}
// CSM_Vec4_s	CSM_Vec4_s::RotationAngles(CSM_Vec4_s const &a){
// 	CSM_Vec4_s RotAngles;
// 	for (int i = 0; i < 3; ++i){
// 		CSM_Vec4_s Tmp1 = *this;
// 		CSM_Vec4_s Tmp2 = a;
// 		Tmp1[i] = Tmp2[i] = 0;
// 		RotAngles[i] = acos(Tmp1.Dot(Tmp2) / (Tmp1.Magnitude() * Tmp2.Magnitude()));
// 	}
// 	return RotAngles;
// }
/*
End methods for CSM_Vec4_s
*/

const CSM_Vec4_s		Normalize(CSM_Vec4_s const &a){
	return CSM_Vec4_s(a).Normalize();
}
double const CSM_Distance(CSM_Vec4_s const &a, CSM_Vec4_s const &b)
{
	return CSM_Vec4_s(a - b).Magnitude();
}


/*
Begin methods for CSM_Mat4_s
*/

CSM_Mat4_s::CSM_Mat4_s(CSM_Vec4_s r1, CSM_Vec4_s r2, CSM_Vec4_s r3, CSM_Vec4_s r4){
	R1 = r1;
	R2 = r2;
	R3 = r3;
	R4 = r4;
}
CSM_Mat4_s::CSM_Mat4_s(CSM_Mat4_s const & a){
	R1 = a.R1;
	R2 = a.R2;
	R3 = a.R3;
	R4 = a.R4;
}

CSM_Vec4_s & CSM_Mat4_s::operator[](unsigned int i){
	switch (i){
		case 0:
			return R1;
			break;
		case 1:
			return R2;
			break;
		case 2:
			return R3;
			break;
		case 3:
			return R4;
			break;
		default:
			throw - 1;
	}
}
CSM_Vec4_s & CSM_Mat4_s::at(unsigned int i){
	return this->operator [](i);
}

CSM_Mat4_s & CSM_Mat4_s::operator*=(CSM_Mat4_s &a){
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			for (int k = 0; k < 4; ++k)
				this->at(i)[j] += this->at(i)[k] * a[j][k];

	return *this;
}
CSM_Vec4_s const CSM_Mat4_s::operator*(CSM_Vec4_s const &a) const{
	CSM_Vec4_s Pdt;
	CSM_Mat4_s aa(*this);
	for (int i = 0; i < 4; ++i)
		Pdt[i] = aa[i].Dot(a);
	return Pdt;
}
CSM_Mat4_s const CSM_Mat4_s::operator*(CSM_Mat4_s &a) const{
	return CSM_Mat4_s(*this) *= a;
}

/*
	End methods for CSM_Mat4_s
*/

/*
 *	Takes Angle in radians
 *	
 *	This code taken from
 *		http://www.programming-techniques.com/2012/03/3d-rotation-algorithm-about-arbitrary.html
 *	
 *	Modified the code written by Bibek Subedi 3/31/2012
 */
const CSM_Mat4_s	CSM_Rotate(double Angle, CSM_Vec3_s const &a){
	double L = CSM_Vec3_s(a).Magnitude();
	double LSqr = L * L;
	CSM_Vec3_s a2 = a * a;

	return CSM_Mat4_s(
		CSM_Vec4_s((a2.X + (a2.Y + a2.Z) * cos(Angle)) / LSqr,
		(a.X * a.Y * (1 - cos(Angle)) - a.Z * L * sin(Angle)) / LSqr,
		(a.X * a.Z * (1 - cos(Angle)) + a.Y * L * sin(Angle)) / LSqr,
		0.0),

		CSM_Vec4_s((a.X * a.Y * (1 - cos(Angle)) + a.Z * L * sin(Angle)) / LSqr,
		(a2.Y + (a2.X + a2.Z) * cos(Angle)) / LSqr,
		(a.Y * a.Z * (1 - cos(Angle)) - a.X * L * sin(Angle)) / LSqr,
		0.0),

		CSM_Vec4_s((a.X * a.Z * (1 - cos(Angle)) - a.Y * L * sin(Angle)) / LSqr,
		(a.Y * a.Z * (1 - cos(Angle)) + a.X * L * sin(Angle)) / LSqr,
		(a2.Z + (a2.X + a2.Y) * cos(Angle)) / LSqr,
		0.0),

		CSM_Vec4_s(0.0,
		0.0,
		0.0,
		1.0)
	);
}