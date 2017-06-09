/*
#pragma once
#ifndef CSMMAT_H_
#define CSMMAT_H_

#include "SmallVector.h"

using std::vector;

namespace CSM{

#ifdef USE_SMALL_VECTOR
	typedef llvm::SmallVector<vec, 3> MatData;
#else
	typedef vector<vec> MatData;
#endif // USE_SMALL_VECTOR

	void func();

	class mat
	{
		friend class mat33;
		friend class mat44;

	public:
		/ *
		*	Constructors
		* /
		mat();
		//	Sets shape only, where I and J are rows and columns
		mat(const unsigned int & I, const unsigned int & J);
		//	Copy constructor
		mat(const mat & rhs);
		//	Sets shape and value of all elements
		mat(const unsigned int & I, const unsigned int & J, const double & rhs);
		//	Construct from vector of vec's (also takes care of C-style 2-d arrays)
		mat(const vector<vec> & rhs);
		//	Construct from GSL matrix
		mat(const gsl_matrix * rhs);

		//	Destructor
		~mat();

		/ *
		*	Element access
		* /
		vec &			operator[](const unsigned int & i);

		/ *
		*	Assignment operators
		* /
		mat &			operator=(const double & rhs);
		mat &			operator=(const vector<vec> & rhs);
		mat &			operator=(const mat & rhs);
		mat &			operator=(const gsl_matrix * rhs);

		/ *
		*	Comparison operators
		* /
		const Boolean_t		operator==(const mat & rhs) const;
		const Boolean_t		operator!=(const mat & rhs) const;
		/ *
		*	Note: >, <, >=, <= operators tell if ANY of the
		*	values are >,<,>=,<= the other's.
		* /
		const Boolean_t		operator<=(const mat &rhs) const;
		const Boolean_t		operator>=(const mat &rhs) const;
		const Boolean_t		operator>(const mat &rhs) const;
		const Boolean_t		operator<(const mat &rhs) const;

		const Boolean_t		IsPos() const;
		const Boolean_t		IsNeg() const;

		/ *
		*	Arithmetic operators
		* /
		const mat		operator-() const;

		/ *
		 *	All binary operators are element-wise
		 *	EXCEPT multiplication, which is proper
		 *	matrix multiplication, hence the extra
		 *	methods for performing element-wise
		 *	multiplication.
		 * /
		mat &			operator+=(const mat &rhs);
		mat &			operator+=(const double &rhs);
		mat &			operator-=(const mat &rhs);
		mat &			operator-=(const double &rhs);
		mat &			operator*=(const mat &rhs);
		mat &			ElemMultAssign(const mat & rhs);
		mat &			operator*=(const double &rhs);
		mat &			operator/=(const mat &rhs);
		mat &			operator/=(const double &rhs);

		const mat		operator+(const mat &rhs) const;
		const mat		operator+(const double &rhs) const;
		const mat		operator-(const mat &rhs) const;
		const mat		operator-(const double &rhs) const;
		/ *
		 *	operator* is real matrix multiplication
		 * /
		const mat		operator*(const mat &rhs) const;
		const mat		ElemMult(const mat & rhs) const;
		const mat		operator*(const double &rhs) const;
		const mat		operator/(const mat &rhs) const;
		const mat		operator/(const double &rhs) const;
		//	Matrix-vector multiplication
		const vec		operator*(const vec &rhs) const;

		/ *
		*	Other methods
		* /
		const Boolean_t		SameDims(const mat & rhs) const;
		const int			Rank() const;
		const vec		At(const unsigned int & i) const;
		const vec		C(const unsigned int & i) const;
		mat &			TransposeSelf();
		const mat		Transpose() const;
		const Boolean_t		EigenSystem(mat & EigVecs, vec & EigVals) const;

	protected:
		MatData m_Vecs;
#ifdef _DEBUG
		vec *r1, *r2, *r3, *r4;
#endif
		void SetPtrs();
	};

	class mat44;
	class mat33;

	class mat22 : public mat
	{
	public:
		using mat::operator=;

		mat22();
		mat22(const mat & rhs);
		mat22(const vec2 & R1, const vec2 & R2);
		mat22(const double & rhs);
		mat22(const vector<vec> & rhs);
		mat22(const gsl_matrix * rhs);
		~mat22();
	};


	class mat33 : public mat
	{
	public:
		using mat::operator=;

		mat33();
		mat33(const mat & rhs);
		mat33(const vec3 & R1, const vec3 & R2, const vec3 & R3);
		mat33(const double & rhs);
		mat33(const vector<vec> & rhs);
		mat33(const gsl_matrix * rhs);
		~mat33();
	};

	class mat44 : public mat
	{
	public:
		using mat::operator=;

		mat44();
		mat44(const mat & rhs);
		mat44(const vec4 & R1, const vec4 & R2, const vec4 & R3, const vec4 & R4);
		mat44(const double & rhs);
		mat44(const vector<vec> & rhs);
		mat44(const gsl_matrix * rhs);
		~mat44();
	};

	const mat44		Rotate(const double & Angle, vec3 &Axis);
	const Boolean_t			EigenSystem(const mat & Hessian, mat & EigVecs, vec & EigVals);

}

#endif // !CSMMAT_H_
*/
