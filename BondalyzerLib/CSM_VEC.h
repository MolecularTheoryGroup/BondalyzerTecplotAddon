/*
#pragma once
#ifndef CSMVEC_H_
#define CSMVEC_H_

#include "SmallVector.h"
#include <gsl/gsl_vector.h>

// #define USE_SMALL_VECTOR

using std::vector;

#ifdef USE_SMALL_VECTOR
typedef llvm::SmallVector<double, 3> VecData;
#else
typedef vector<double> VecData;
#endif // USE_SMALL_VECTOR

namespace CSM{

	class vec
	{
	public:
		/ *
		 *	Constructors
		 * /
		vec();
		//	Sets size of vec only
		vec(const unsigned int & N);
		//	Copy constructor
		vec(const vec &rhs);
		//	Set size and value of all elements
		vec(const unsigned int & N, const double & rhs);
		//	Construct from a vector of doubles (also takes care of C-style arrays)
		vec(const vector<double> & rhs);
		//	Construct from GSL vector
		vec(const gsl_vector * rhs);

		//	Destructor
		~vec();

		/ *
		 *	Element access
		 * /
		double &			operator[](const unsigned int & i);
		double &			X();
		double &			Y();
		double &			Z();
		double &			W();

		/ *
		 *	Assignment operators
		 * /
		vec &			operator=(const vec &rhs);
		vec &			operator=(const double & rhs);
		vec &			operator=(const vector<double> & rhs);
		vec &			operator=(const gsl_vector * GSLVec3);

		/ *
		 *	Comparison operators
		 * /
		const Boolean_t		operator==(const vec &rhs) const;
		const Boolean_t		operator!=(const vec &rhs) const;
		/ *
		*	Note: >, <, >=, <= operators tell if ANY of the
		*	values are >,<,>=,<= the other's.
		* /
		const Boolean_t		operator<=(const vec &rhs) const;
		const Boolean_t		operator>=(const vec &rhs) const;
		const Boolean_t		operator>(const vec &rhs) const;
		const Boolean_t		operator<(const vec &rhs) const;

		const Boolean_t		IsPos() const;
		const Boolean_t		IsNeg() const;

		/ *
		 *	Arithmetic operators
		 * /
		const vec		operator-() const;

		vec &			operator+=(const vec &rhs);
		vec &			operator+=(const double &rhs);
		vec &			operator-=(const vec &rhs);
		vec &			operator-=(const double &rhs);
		vec &			operator*=(const vec &rhs);
		vec &			operator*=(const double &rhs);
		vec &			operator/=(const vec &rhs);
		vec &			operator/=(const double &rhs);

		const vec		operator+(const vec &rhs) const;
		const vec		operator+(const double &rhs) const;
		const vec		operator-(const vec &rhs) const;
		const vec		operator-(const double &rhs) const;
		const vec		operator*(const vec &rhs) const;
		const vec		operator*(const double & rhs) const;
		const vec		operator/(const vec &rhs) const;
		const vec		operator/(const double & rhs) const;

		/ *
		 *	Other methods
		 * /
		const unsigned int	Size() const;
		const double		At(const unsigned int & i) const;
		const double		Dot(const vec &rhs) const;
		vec &			NormalizeSelf();
		const vec		Normalize() const;
		const double		Magnitude() const;
		const double		MagSqr() const;
		const double		Distance(const vec & rhs) const;
		const double		DistSqr(const vec & rhs) const;

	protected:
		VecData m_Data;
#ifdef _DEBUG
		double *x, *y, *z, *w;
#endif
		void SetPtrs();
	};

	class vec4;
	class vec3;

	class vec2 : public vec
	{
	public:
		using vec::operator=;

		vec2(const double & x = 0, const double & y = 0);
		vec2(const vec &rhs);
		vec2(const double & rhs);
		vec2(const vector<double> & rhs);
		vec2(const vec3 &rhs);
		vec2(const vec4 &rhs);
		vec2(const gsl_vector * rhs);
		~vec2();
	};

	class vec3 : public vec
	{
	public:
		using vec::operator=;

		vec3(const double & x = 0, const double & y = 0, const double & z = 0);
		vec3(const vec &rhs);
		vec3(const double & rhs);
		vec3(const double * rhs);
		vec3(const vector<double> & rhs);
		vec3(const vec2 &rhs, const double & z = 0);
		vec3(const vec4 &rhs);
		vec3(const gsl_vector * rhs);
		~vec3();

		vec3 &			operator=(const double * rhs);

		const vec3 Cross(const vec3 &rhs) const;
		const vec3 RotationAngles(const vec3 &rhs) const;
	};

	const vec3 Transform2dTo3d(const vec2 & TwoPt,
		const vector<vec3> & BasisVectors,
		const vec3 & Origin);

	class vec4 : public vec
	{
	public:
		using vec::operator=;

		vec4(const double & x = 0, const double & y = 0, const double & z = 0, const double & w = 0);
		vec4(const vec &rhs);
		vec4(const double & rhs);
		vec4(const vector<double> & rhs);
		vec4(const vec2 &rhs, const double & z = 0, const double & w = 0);
		vec4(const vec3 &rhs, const double & w = 0);
		vec4(const gsl_vector * rhs);
		~vec4();
	};

	const vec		Normalize(const vec &rhs);
	const double		Distance(const vec &lhs, const vec &rhs);

}


#endif*/