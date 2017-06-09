/*

#include <vector>
#include "SmallVector.h"

#include <gsl/gsl_vector.h>
#include "TECADDON.h"

#include "CSM_VEC.h"


namespace CSM{

	/ *
	 *	Begin vec base class methods
	 * /

	/ *
	*	Constructors
	* /
	vec::vec(){}
	//	Sets size of vec only
	vec::vec(const unsigned int & N){
		REQUIRE(N > 0);
		m_Data.resize(N);
		SetPtrs();
	}
	//	Copy constructor
	vec::vec(const vec &rhs){ *this = rhs; }
	//	Set size and value of all elements
	vec::vec(const unsigned int & N, const double & rhs) : vec(N) { *this = rhs; }
	//	Construct from a vector of doubles (also takes care of C-style arrays)
	vec::vec(const vector<double> & rhs){ *this = rhs; }
	//	Construct from GSL vector
	vec::vec(const gsl_vector * rhs){ *this = rhs; }

	//	Destructor
	vec::~vec(){}

	/ *
	*	Element access
	* /
	double &	vec::operator[](const unsigned int & i){
		// 	REQUIRE(i >= 0 && i < m_Data.size());

		return m_Data[i];
	}
	double & vec::X(){
		REQUIRE(m_Data.size() >= 1);
		return m_Data[0];
	}
	double & vec::Y(){
		REQUIRE(m_Data.size() >= 2);
		return m_Data[1];
	}
	double & vec::Z(){
		REQUIRE(m_Data.size() >= 3);
		return m_Data[2];
	}
	double & vec::W(){
		REQUIRE(m_Data.size() >= 4);
		return m_Data[3];
	}

	/ *
	*	Assignment operators
	* /
	vec & vec::operator=(const vec &rhs){
		if (this == &rhs)
			return *this;

		REQUIRE(m_Data.size() == rhs.m_Data.size() || m_Data.size() == 0);

		m_Data = rhs.m_Data;
		SetPtrs();
		return *this;
	}
	vec & vec::operator=(const double & rhs){
		m_Data.resize(m_Data.size(), rhs);
		SetPtrs();
		return *this;
	}
	vec & vec::operator=(const vector<double> & rhs){
#ifdef USE_SMALL_VECTOR
		m_Data.resize(rhs.size());
		for (unsigned int i = 0; i < rhs.size(); ++i)
			m_Data[i] = rhs[i];
#else
		m_Data = rhs;
#endif
		SetPtrs();
		return *this;
	}
	vec & vec::operator=(const gsl_vector * rhs){
		REQUIRE(rhs->size > 0);

		m_Data.resize(rhs->size);
		for (unsigned int i = 0; i < m_Data.size(); ++i)
			m_Data[i] = gsl_vector_get(rhs, i);

		SetPtrs();
		return *this;
	}

	/ *
	*	Comparison operators
	* /
	const Boolean_t vec::operator==(const vec &rhs) const{
		return (m_Data == rhs.m_Data);
	}
	const Boolean_t vec::operator!=(const vec &rhs) const{
		return !(*this == rhs);
	}

	const Boolean_t		vec::operator<=(const vec &rhs) const{
		REQUIRE(m_Data.size() == rhs.m_Data.size());

		for (unsigned int i = 0; i < m_Data.size(); ++i)
			if (!(m_Data[i] <= rhs.m_Data[i]))
				return FALSE;

		return TRUE;
	}
	const Boolean_t		vec::operator>=(const vec &rhs) const{
		REQUIRE(m_Data.size() == rhs.m_Data.size());

		for (unsigned int i = 0; i < m_Data.size(); ++i)
			if (!(m_Data[i] >= rhs.m_Data[i]))
				return FALSE;

		return TRUE;
	}
	const Boolean_t		vec::operator>(const vec &rhs) const{
		return !(*this <= rhs);
	}
	const Boolean_t		vec::operator<(const vec &rhs) const{
		return !(*this >= rhs);
	}

	const Boolean_t		vec::IsPos() const{
		for (const auto & i : m_Data)
			if (i < 0)
				return FALSE;

		return TRUE;
	}
	const Boolean_t		vec::IsNeg() const{
		return !IsPos();
	}

	/ *
	*	Arithmetic operators
	* /
	const vec vec::operator-() const{
		vec Tmp = *this;
		for (auto & i : Tmp.m_Data)
			i *= -1.0;

		return Tmp;
	}

	vec & vec::operator+=(const vec &rhs){
		REQUIRE(m_Data.size() == rhs.m_Data.size());

		for (unsigned int i = 0; i < m_Data.size(); ++i)
			m_Data[i] += rhs.m_Data[i];

		return *this;
	}
	vec & vec::operator+=(const double &rhs){
		for (auto & i : m_Data)
			i += rhs;

		return *this;
	}
	vec & vec::operator-=(const vec &rhs){
		REQUIRE(m_Data.size() == rhs.m_Data.size());

		for (unsigned int i = 0; i < m_Data.size(); ++i)
			m_Data[i] -= rhs.m_Data[i];

		return *this;
	}
	vec & vec::operator-=(const double &rhs){
		for (auto & i : m_Data)
			i -= rhs;

		return *this;
	}
	vec & vec::operator*=(const vec &rhs){
		REQUIRE(m_Data.size() == rhs.m_Data.size());

		for (unsigned int i = 0; i < m_Data.size(); ++i)
			m_Data[i] *= rhs.m_Data[i];

		return *this;
	}
	vec & vec::operator*=(const double &rhs){
		for (auto & i : m_Data)
			i *= rhs;

		return *this;
	}
	vec & vec::operator/=(const vec &rhs){
		REQUIRE(m_Data.size() == rhs.m_Data.size());

		for (unsigned int i = 0; i < m_Data.size(); ++i){
			REQUIRE(rhs.m_Data[i] != 0.0);
			m_Data[i] /= rhs.m_Data[i];
		}

		return *this;
	}
	vec & vec::operator/=(const double &rhs){
		REQUIRE(rhs != 0.0);

		for (auto & i : m_Data)
			i /= rhs;

		return *this;
	}

	const vec vec::operator+(const vec &rhs) const{
		return vec(*this) += rhs;
	}
	const vec	vec::operator+(const double & rhs) const{
		return vec(*this) += rhs;
	}
	const vec vec::operator-(const vec &rhs) const{
		return vec(*this) -= rhs;
	}
	const vec	vec::operator-(const double & rhs) const{
		return vec(*this) -= rhs;
	}
	const vec	vec::operator*(const vec &rhs) const{
		return vec(*this) *= rhs;
	}
	const vec	vec::operator*(const double & rhs) const{
		return vec(*this) *= rhs;
	}
	const vec	vec::operator/(const vec &rhs) const{
		return vec(*this) /= rhs;
	}
	const vec	vec::operator/(const double & rhs) const{
		return vec(*this) /= rhs;
	}

	/ *
	*	Other methods
	* /
	const unsigned int  vec::Size() const{
		return static_cast<unsigned int>(m_Data.size());
	}
	const double		vec::At(const unsigned int & i) const{
		// 	REQUIRE(i >= 0 && i < m_Data.size());

		return m_Data[i];
	}
	const double vec::Dot(const vec &rhs) const{
		REQUIRE(m_Data.size() == rhs.m_Data.size());

		double Value = 0.0;
		for (unsigned int i = 0; i < m_Data.size(); ++i)
			Value += m_Data[i] * rhs.m_Data[i];

		return Value;
	}
	vec &		 vec::NormalizeSelf(){
		double Mag = Magnitude();
		if (Mag > 0)
			operator/=(Mag);

		return *this;
	}
	const vec	vec::Normalize() const{
		return vec(*this).NormalizeSelf();
	}
	const double vec::Magnitude() const
	{
		double Mag = 0.0;
		for (const auto & i : m_Data)
			Mag += pow(i, 2);

		return sqrt(Mag);
	}
	const double vec::MagSqr() const
	{
		double Mag = 0.0;
		for (const auto & i : m_Data)
			Mag += pow(i, 2);

		return Mag;
	}
	const double		vec::Distance(const vec & rhs) const{
		return sqrt(DistSqr(rhs));
	}
	const double		vec::DistSqr(const vec & rhs) const{
		REQUIRE(m_Data.size() == rhs.m_Data.size());
		double SqrDist = 0.0;
		for (unsigned int i = 0; i < m_Data.size(); ++i)
			SqrDist += pow(rhs.m_Data[i] - m_Data[i], 2);

		return SqrDist;
	}

	void vec::SetPtrs(){
#ifdef _DEBUG
		x = y = z = w = NULL;
		if (m_Data.size() >= 1) x = &m_Data[0];
		if (m_Data.size() >= 2) y = &m_Data[1];
		if (m_Data.size() >= 3) z = &m_Data[2];
		if (m_Data.size() >= 4) w = &m_Data[3];
#endif
	}

	/ *
	*	Begin vec2 methods
	* /

	vec2::vec2(const double & x, const double & y) : vec(2){
		*this = vector < double > { x, y };
	}
	vec2::vec2(const vec &rhs) : vec(rhs){
		REQUIRE(m_Data.size() == 2);
	}
	vec2::vec2(const double & rhs) : vec(2, rhs){

	}
	vec2::vec2(const vector<double> & rhs) : vec(rhs){
		REQUIRE(m_Data.size() == 2);
	}
	//	Construct from Vec3, using only x,y (first 2) values
	vec2::vec2(const vec3 &rhs) : vec(2){
		for (unsigned int i = 0; i < 2; ++i)
			m_Data[i] = rhs.At(i);
	}
	//	Construct from Vec4, using only x,y (first 2) values
	vec2::vec2(const vec4 &rhs) : vec(2){
		for (unsigned int i = 0; i < 2; ++i)
			m_Data[i] = rhs.At(i);
		SetPtrs();
	}
	vec2::vec2(const gsl_vector * rhs) : vec(rhs){
		REQUIRE(m_Data.size() == 2);
	}
	vec2::~vec2(){}


	/ *
	 *	Begin vec3 methods
	 * /

	vec3::vec3(const double & x, const double & y, const double & z) : vec(3){
		*this = vector < double > { x, y, z };
	}
	vec3::vec3(const vec &rhs) : vec(rhs){ REQUIRE(m_Data.size() == 3); }
	vec3::vec3(const double & rhs) : vec(3, rhs){}
	vec3::vec3(const double * rhs) : vec(3){ *this = rhs; }
	vec3::vec3(const vector<double> & rhs) : vec(rhs){ REQUIRE(m_Data.size() == 3); }
	//	Construct from Vec2 and a third value
	vec3::vec3(const vec2 &rhs, const double & z) : vec(rhs){
		REQUIRE(m_Data.size() == 2);
		m_Data.resize(3);
		m_Data[2] = z;
		SetPtrs();
	}
	//	Construct from Vec4, using only x,y,z (first 3) values
	vec3::vec3(const vec4 &rhs) : vec(3){
		for (unsigned int i = 0; i < 3; ++i)
			m_Data[i] = rhs.At(i);
		SetPtrs();
	}
	vec3::vec3(const gsl_vector * rhs) : vec(rhs){ REQUIRE(m_Data.size() == 3); }
	vec3::~vec3(){}

	vec3 & vec3::operator=(const double * rhs){
#ifdef USE_SMALL_VECTOR
		m_Data.resize(3);
		for (unsigned int i = 0; i < 3; ++i)
			m_Data[i] = rhs[i];
#else
		m_Data = { rhs[0], rhs[1], rhs[2] };
#endif // USE_SMALL_VECTOR

		SetPtrs();
		return *this;
	}

	const vec3 vec3::Cross(const vec3 &rhs) const{
		vec3 Tmp;

		Tmp.m_Data[0] = -m_Data[2] * rhs.m_Data[1] + m_Data[1] * rhs.m_Data[2];
		Tmp.m_Data[1] = m_Data[2] * rhs.m_Data[0] - m_Data[0] * rhs.m_Data[2];
		Tmp.m_Data[2] = -m_Data[1] * rhs.m_Data[0] + m_Data[0] * rhs.m_Data[1];
		// 	Tmp.X = -Z*rhs.Y + Y*rhs.Z;
		// 	Tmp.Y = Z*rhs.X - X*rhs.Z;
		// 	Tmp.Z = -Y*rhs.X + X*rhs.Y;
		return Tmp;
	}
	const vec3 vec3::RotationAngles(const vec3 &rhs) const{
		vec3 RotAngles;
		for (int i = 0; i < 3; ++i){
			vec3 Tmp1 = *this;
			vec3 Tmp2 = rhs;
			Tmp1[i] = Tmp2[i] = 0;
			RotAngles[i] = acos(Tmp1.Dot(Tmp2) / (Tmp1.Magnitude() * Tmp2.Magnitude()));
		}
		return RotAngles;
	}

	const vec3 Transform2dTo3d(const vec2 & TwoPt,
		const vector<vec3> & BasisVectors,
		const vec3 & Origin)
	{
		vec3 ThreePt(Origin);
		for (int i = 0; i < 2; ++i){
			ThreePt += BasisVectors[i] * TwoPt.At(i);
		}

		return ThreePt;
	}

	/ *
	 *	Begin vec4 methods
	 * /

	vec4::vec4(const double & x, const double & y, const double & z, const double & w) : vec(4){
		*this = vector < double > { x, y, z, w };
	}
	vec4::vec4(const vec &rhs) : vec(rhs){ REQUIRE(m_Data.size() == 4); }
	vec4::vec4(const double & rhs) : vec(4, rhs){}
	vec4::vec4(const vector<double> & rhs) : vec(rhs){ REQUIRE(m_Data.size() == 4); }
	//	Construct from Vec3 and a fourth value
	vec4::vec4(const vec2 &rhs, const double & z, const double & w) : vec(rhs){
		REQUIRE(m_Data.size() == 2);
		m_Data.resize(4);
		m_Data[2] = z;
		m_Data[3] = w;
		SetPtrs();
	}
	//	Construct from Vec3 and a fourth value
	vec4::vec4(const vec3 &rhs, const double & w) : vec(rhs){
		REQUIRE(m_Data.size() == 3);
		m_Data.resize(4);
		m_Data[3] = w;
		SetPtrs();
	}
	vec4::vec4(const gsl_vector * rhs) : vec(rhs){ REQUIRE(m_Data.size() == 4); }
	vec4::~vec4(){}

	/ *
	 *	Begin global functions for vec
	 * /

	const vec	Normalize(const vec &rhs){
		return rhs.Normalize();
	}
	const double Distance(const vec &rhs, const vec &b){
		return rhs.Distance(b);
	}

}*/