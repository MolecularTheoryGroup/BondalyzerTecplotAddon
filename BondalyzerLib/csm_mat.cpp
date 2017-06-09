/*
#include <vector>
#include "SmallVector.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "TECADDON.h"

#include "CSM_VEC.h"
#include "CSM_MAT.h"

using std::vector;

namespace CSM{

	/ *
	* Public domain code for finding eigensystems of 3x3 real matrices
	* taken from: http://barnesc.blogspot.com/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
	* /

	/ * Eigen-decomposition for symmetric 3x3 real matrices.
	Public domain, copied from the public domain Java library JAMA. * /
	/ * Symmetric matrix A => eigenvectors in columns of V, corresponding
	eigenvalues in d. * /
	const Boolean_t eigen_decomposition(const mat & A, mat & V, vec & d);

	/ *
	 *	Begin mat methods
	 * /

	/ *
	*	Constructors
	* /
	mat::mat(){}
	//	Sets shape only, where I and J are rows and columns
	mat::mat(const unsigned int & I, const unsigned int & J){
		REQUIRE(I > 0 && J > 0);
		m_Vecs.resize(I, vec(J));
		SetPtrs();
	}
	//	Copy constructor
	mat::mat(const mat & rhs){ *this = rhs; }
	//	Sets shape and value of all elements
	mat::mat(const unsigned int & I, const unsigned int & J, const double & rhs) : mat(I, J){ *this = rhs; }
	//	Construct from vector of vec's (also takes care of C-style 2-d arrays)
	mat::mat(const vector<vec> & rhs){ *this = rhs; }
	//	Construct from GSL matrix
	mat::mat(const gsl_matrix * rhs){ *this = rhs; }

	//	Destructor
	mat::~mat(){}

	/ *
	*	Element access
	* /
	vec &			mat::operator[](const unsigned int & i){
		// 	REQUIRE(i >= 0 && i < m_Vecs.size());

		return m_Vecs[i];
	}

	/ *
	*	Assignment operators
	* /
	mat &			mat::operator=(const double & rhs){
		m_Vecs.resize(m_Vecs.size(), vec(m_Vecs[0].Size(), rhs));
		SetPtrs();
		return *this;
	}
	mat &			mat::operator=(const vector<vec> & rhs){
#ifdef USE_SMALL_VECTOR
		m_Vecs.resize(rhs.size());
		for (unsigned int i = 0; i < rhs.size(); ++i)
			m_Vecs[i] = rhs[i];
#else
		m_Vecs = rhs;
#endif
		SetPtrs();
		return *this;
	}
	mat &			mat::operator=(const mat & rhs){
		if (this == &rhs)
			return *this;

		REQUIRE(SameDims(rhs) || m_Vecs.size() == 0);

		m_Vecs = rhs.m_Vecs;
		SetPtrs();
		return *this;
	}
	mat &			mat::operator=(const gsl_matrix * rhs){
		unsigned int
			N = static_cast<unsigned int>(rhs->size1),
			M = static_cast<unsigned int>(rhs->size2);

		REQUIRE(N > 0 && M > 0);

		m_Vecs.resize(N, vec(M));
		for (unsigned int i = 0; i < N; ++i)
			for (unsigned int j = 0; j < M; ++j)
				m_Vecs[i][j] = gsl_matrix_get(rhs, i, j);

		SetPtrs();
		return *this;
	}

	/ *
	*	Comparison operators
	* /
	const Boolean_t		mat::operator==(const mat & rhs) const{
		return (m_Vecs == rhs.m_Vecs);
	}
	const Boolean_t		mat::operator!=(const mat & rhs) const{
		return !(*this == rhs);
	}
	/ *
	*	Note: >, <, >=, <= operators tell if ANY of the
	*	values are >,<,>=,<= the other's.
	* /
	const Boolean_t		mat::operator<=(const mat &rhs) const{
		REQUIRE(SameDims(rhs));

		for (unsigned int i = 0; i < m_Vecs.size(); ++i)
			if (!(m_Vecs[i] <= rhs.m_Vecs[i]))
				return FALSE;

		return TRUE;
	}
	const Boolean_t		mat::operator>=(const mat &rhs) const{
		REQUIRE(SameDims(rhs));

		for (unsigned int i = 0; i < m_Vecs.size(); ++i)
			if (!(m_Vecs[i] >= rhs.m_Vecs[i]))
				return FALSE;

		return TRUE;
	}
	const Boolean_t		mat::operator>(const mat &rhs) const{
		return !(*this <= rhs);
	}
	const Boolean_t		mat::operator<(const mat &rhs) const{
		return !(*this >= rhs);
	}

	const Boolean_t		mat::IsPos() const{
		for (const auto & i : m_Vecs)
			if (!i.IsPos())
				return FALSE;

		return TRUE;
	}
	const Boolean_t		mat::IsNeg() const{
		return !IsPos();
	}

	/ *
	*	Arithmetic operators
	* /
	const mat		mat::operator-() const{
		mat Tmp = *this;
		for (auto & i : Tmp.m_Vecs)
			i *= -1.0;

		return Tmp;
	}

	mat &			mat::operator+=(const mat &rhs){
		REQUIRE(SameDims(rhs));

		for (unsigned int i = 0; i < m_Vecs.size(); ++i)
			m_Vecs[i] += rhs.m_Vecs[i];

		return *this;
	}
	mat &			mat::operator+=(const double &rhs){
		for (auto & i : m_Vecs)
			i += rhs;

		return *this;
	}
	mat &			mat::operator-=(const mat &rhs){
		REQUIRE(SameDims(rhs));

		for (unsigned int i = 0; i < m_Vecs.size(); ++i)
			m_Vecs[i] -= rhs.m_Vecs[i];

		return *this;
	}
	mat &			mat::operator-=(const double &rhs){
		for (auto & i : m_Vecs)
			i -= rhs;

		return *this;
	}
	/ *
	 *	This is real matrix multiplication
	 * /
	mat &			mat::operator*=(const mat &rhs){
		unsigned int
			N1 = static_cast<unsigned int>(m_Vecs.size()),
			N2 = static_cast<unsigned int>(rhs.m_Vecs.size()),
			M1 = m_Vecs[0].Size(),
			M2 = rhs.m_Vecs[0].Size();

		REQUIRE(M1 == N2);

		mat Pdt(N1, M2);
		for (unsigned int i = 0; i < N1; ++i)
			for (unsigned int j = 0; j < M2; ++j)
				for (unsigned int k = 0; k < M1; ++k)
					Pdt[i][j] += At(i).At(k) * rhs.At(j).At(k);

		*this = Pdt;
		return *this;
	}
	mat &			mat::ElemMultAssign(const mat &rhs){
		REQUIRE(SameDims(rhs));

		for (unsigned int i = 0; i < m_Vecs.size(); ++i)
			m_Vecs[i] += rhs.m_Vecs[i];

		return *this;
	}
	mat &			mat::operator*=(const double &rhs){
		for (auto & i : m_Vecs)
			i *= rhs;

		return *this;
	}
	mat &			mat::operator/=(const mat &rhs){
		REQUIRE(SameDims(rhs));

		for (unsigned int i = 0; i < m_Vecs.size(); ++i)
			m_Vecs[i] /= rhs.m_Vecs[i];

		return *this;
	}
	mat &			mat::operator/=(const double &rhs){
		REQUIRE(rhs != 0.0);

		for (auto & i : m_Vecs)
			i /= rhs;

		return *this;
	}

	const mat		mat::operator+(const mat &rhs) const{
		return mat(*this) += rhs;
	}
	const mat		mat::operator+(const double &rhs) const{
		return mat(*this) += rhs;
	}
	const mat		mat::operator-(const mat &rhs) const{
		return mat(*this) -= rhs;
	}
	const mat		mat::operator-(const double &rhs) const{
		return mat(*this) -= rhs;
	}
	const mat		mat::operator*(const mat &rhs) const{
		return mat(*this) *= rhs;
	}
	const mat		mat::ElemMult(const mat &rhs) const{
		return mat(*this).ElemMultAssign(rhs);
	}
	const mat		mat::operator*(const double &rhs) const{
		return mat(*this) *= rhs;
	}
	const mat		mat::operator/(const mat &rhs) const{
		return mat(*this) /= rhs;
	}
	const mat		mat::operator/(const double &rhs) const{
		return mat(*this) /= rhs;
	}

	const vec		mat::operator*(const vec &rhs) const{
		REQUIRE(m_Vecs[0].Size() == rhs.Size());

		unsigned int N = static_cast<unsigned int>(m_Vecs.size());
		vec Pdt(N);
		for (unsigned int i = 0; i < N; ++i)
			Pdt[i] = At(i).Dot(rhs);

		return Pdt;
	}

	/ *
	*	Other methods
	* /
	const Boolean_t		mat::SameDims(const mat & rhs) const{
		Boolean_t Same = (m_Vecs.size() == rhs.m_Vecs.size());

		unsigned int N = static_cast<unsigned int>(m_Vecs.size());
		for (unsigned int i = 0; i < N && Same; ++i)
			Same = (m_Vecs[i].Size() == rhs.m_Vecs[i].Size());

		return Same;
	}
	const int			mat::Rank() const{
		int Dim = static_cast<int>(m_Vecs.size());

		for (const auto & i : m_Vecs)
			if (Dim != i.Size()){
				Dim = -1;
				break;
			}

		return Dim;
	}
	const vec		mat::At(const unsigned int & i) const{
		// 	REQUIRE(i >= 0 && i < m_Vecs.size());

		return m_Vecs[i];
	}
	const vec		mat::C(const unsigned int & i) const{
		// 	REQUIRE(i >= 0 && i < m_Vecs[0].Size());

		vector<double> Col;
		Col.reserve(m_Vecs.size());

		for (const auto & j : m_Vecs)
			Col.push_back(j.At(i));

		return vec(Col);
	}
	mat &			mat::TransposeSelf(){
		MatData Vecs;
		Vecs.reserve(m_Vecs[0].Size());

		for (unsigned int i = 0; i < m_Vecs[0].Size(); ++i)
			Vecs.push_back(C(i));

		m_Vecs = Vecs;

		return *this;
	}
	const mat		mat::Transpose() const{
		return mat(*this).TransposeSelf();
	}
	const Boolean_t		mat::EigenSystem(mat & EigVecs, vec & EigVals) const{
		int MyRank = Rank();
		REQUIRE(MyRank > 0);

		Boolean_t IsOk = (MyRank > 0);

		if (IsOk){
			/ *
			 *	Method using public domain code. Info at the top of the file.
			 * /
			IsOk = eigen_decomposition(*this, EigVecs, EigVals);

			// 		/ *
			// 		 *	Method using GSL's eigensystem solver. Pretty slow even though it's for symmetric matrices...
			// 		 * /
			// 		/ *
			// 		*	Prepare gsl data structures for the Hessian and eigen vector matrices
			// 		*	and eigen value vector.
			// 		* /
			// 		gsl_matrix * Hess = gsl_matrix_alloc(MyRank, MyRank);
			// 		gsl_matrix * EigenVectors = gsl_matrix_alloc(MyRank, MyRank);
			// 		gsl_vector * EigenValues = gsl_vector_alloc(MyRank);
			// 
			// 		for (int i = 0; i < MyRank; ++i)
			// 			for (int j = 0; j < MyRank; ++j)
			// 				gsl_matrix_set(Hess, i, j, m_Vecs[i].At(j));
			// 
			// 
			// 		// Setup the GSL eigensystem workspace
			// 		gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(MyRank);
			// 
			// 		// Solve the eigensystem
			// 		gsl_eigen_symmv(Hess, EigenValues, EigenVectors, w);
			// 		/ *
			// 		*	Simultaneously sort the eigen values and vectors in ascending order according
			// 		*	to the eigenvalues.
			// 		* /
			// 		gsl_eigen_symmv_sort(EigenValues, EigenVectors, GSL_EIGEN_SORT_VAL_ASC);
			// 
			// 		/ *
			// 		*	Store the results in the CSM vector and matrix provided
			// 		* /
			// 		EigVals = EigenValues;
			// 		EigVecs = EigenVectors;
			// 
			// 		IsOk = (EigVals.Size() == MyRank && EigVecs.Rank() == MyRank);
			// 
			// 		/ *
			// 		*	Clear the workspace
			// 		* /
			// 		gsl_eigen_symmv_free(w);
			// 		gsl_matrix_free(Hess);
			// 		gsl_matrix_free(EigenVectors);
			// 		gsl_vector_free(EigenValues);
		}

		return IsOk;
	}

	void mat::SetPtrs(){
#ifdef _DEBUG
		r1 = r2 = r3 = r4 = NULL;
		if (m_Vecs.size() >= 1) r1 = &m_Vecs[0];
		if (m_Vecs.size() >= 2) r2 = &m_Vecs[1];
		if (m_Vecs.size() >= 3) r3 = &m_Vecs[2];
		if (m_Vecs.size() >= 4) r4 = &m_Vecs[3];
#endif
	}

	/ *
	*	Begin mat22 methods
	* /

	mat22::mat22() : mat(2, 2, 0.0){}
	mat22::mat22(const mat & rhs) : mat(rhs){ REQUIRE(Rank() == 2); }
	mat22::mat22(const vec2 & R1, const vec2 & R2) : mat(vector < vec > {R1, R2}){}
	mat22::mat22(const double & rhs) : mat(2, 2, rhs){}
	mat22::mat22(const vector<vec> & rhs) : mat(rhs){ REQUIRE(Rank() == 2); }
	mat22::mat22(const gsl_matrix * rhs) : mat(rhs){ REQUIRE(Rank() == 2); }

	mat22::~mat22(){}

	/ *
	 *	Begin mat33 methods
	 * /

	mat33::mat33() : mat(3, 3, 0.0){}
	mat33::mat33(const mat & rhs) : mat(rhs){ REQUIRE(Rank() == 3); }
	mat33::mat33(const vec3 & R1, const vec3 & R2, const vec3 & R3) : mat(vector < vec > {R1, R2, R3}){}
	mat33::mat33(const double & rhs) : mat(3, 3, rhs){}
	mat33::mat33(const vector<vec> & rhs) : mat(rhs){ REQUIRE(Rank() == 3); }
	mat33::mat33(const gsl_matrix * rhs) : mat(rhs){ REQUIRE(Rank() == 3); }

	mat33::~mat33(){}

	/ *
	 *	Begin mat44 methods
	 * /

	mat44::mat44() : mat(4, 4, 0.0){}
	mat44::mat44(const mat & rhs) : mat(rhs){ REQUIRE(Rank() == 4); }
	mat44::mat44(const vec4 & R1, const vec4 & R2, const vec4 & R3, const vec4 & R4) : mat(vector < vec > {R1, R2, R3, R4}){}
	mat44::mat44(const double & rhs) : mat(4, 4, rhs){}
	mat44::mat44(const vector<vec> & rhs) : mat(rhs){ REQUIRE(Rank() == 4); }
	mat44::mat44(const gsl_matrix * rhs) : mat(rhs){ REQUIRE(Rank() == 4); }

	mat44::~mat44(){}

	/ *
	 *	Begin global mat methods
	 * /

	const mat44		Rotate(const double & Angle, vec3 &Axis){
		double L = Axis.Magnitude();
		double LSqr = L * L;
		vec3 AxisSqr = Axis * Axis;

		return mat44(
			vec4((AxisSqr(0) + (AxisSqr(1) + AxisSqr(2)) * cos(Angle)) / LSqr,
			(Axis(0) * Axis(1) * (1 - cos(Angle)) - Axis(2) * L * sin(Angle)) / LSqr,
			(Axis(0) * Axis(2) * (1 - cos(Angle)) + Axis(1) * L * sin(Angle)) / LSqr,
			0.0),

			vec4((Axis(0) * Axis(1) * (1 - cos(Angle)) + Axis(2) * L * sin(Angle)) / LSqr,
			(AxisSqr(1) + (AxisSqr(0) + AxisSqr(2)) * cos(Angle)) / LSqr,
			(Axis(1) * Axis(2) * (1 - cos(Angle)) - Axis(0) * L * sin(Angle)) / LSqr,
			0.0),

			vec4((Axis(0) * Axis(2) * (1 - cos(Angle)) - Axis(1) * L * sin(Angle)) / LSqr,
			(Axis(1) * Axis(2) * (1 - cos(Angle)) + Axis(0) * L * sin(Angle)) / LSqr,
			(AxisSqr(2) + (AxisSqr(0) + AxisSqr(1)) * cos(Angle)) / LSqr,
			0.0),

			vec4(0.0,
			0.0,
			0.0,
			1.0)
			);
	}

	const Boolean_t			EigenSystem(const mat & Hessian, mat & EigVecs, vec & EigVals){
		Boolean_t IsOk = Hessian.EigenSystem(EigVecs, EigVals);

		return IsOk;
	}

	/ *
	* Below is come public domain code for finding the eigensystem of 3x3 real symmetric matrices
	* taken from: http://barnesc.blogspot.com/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
	* /

	/ * Eigen decomposition code for symmetric 3x3 matrices, copied from the public
	domain Java Matrix library JAMA. * /

	/ *
	* C++ module 'eig3' by Connelly Barnes
	------------------------------------

	License: public domain.

	The source files in this directory have been copied from the public domain
	Java matrix library JAMA.  The derived source code is in the public domain
	as well.

	Usage:

	Symmetric matrix A => eigenvectors in columns of V, corresponding
	eigenvalues in d.
	void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);
	* /


	// #ifdef MAX
	// #undef MAX
	// #endif
	// 
	// #define MAX(a, b) ((a)>(b)?(a):(b))
	// 
	// #define n 3

	static double hypot2(double x, double y) {
		return sqrt(x*x + y*y);
	}

	// Symmetric Householder reduction to tridiagonal form.

	static void tred2(mat & V, vec & d, vec & e) {
		int n = d.Size();

		//  This is derived from the Algol procedures tred2 by
		//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		//  Fortran subroutine in EISPACK.

		for (int j = 0; j < n; j++) {
			d[j] = V[n - 1][j];
		}

		// Householder reduction to tridiagonal form.

		for (int i = n - 1; i > 0; i--) {

			// Scale to avoid under/overflow.

			double scale = 0.0;
			double h = 0.0;
			for (int k = 0; k < i; k++) {
				scale = scale + fabs(d[k]);
			}
			if (scale == 0.0) {
				e[i] = d[i - 1];
				for (int j = 0; j < i; j++) {
					d[j] = V[i - 1][j];
					V[i][j] = 0.0;
					V[j][i] = 0.0;
				}
			}
			else {

				// Generate Householder vector.

				for (int k = 0; k < i; k++) {
					d[k] /= scale;
					h += d[k] * d[k];
				}
				double f = d[i - 1];
				double g = sqrt(h);
				if (f > 0) {
					g = -g;
				}
				e[i] = scale * g;
				h = h - f * g;
				d[i - 1] = f - g;
				for (int j = 0; j < i; j++) {
					e[j] = 0.0;
				}

				// Apply similarity transformation to remaining columns.

				for (int j = 0; j < i; j++) {
					f = d[j];
					V[j][i] = f;
					g = e[j] + V[j][j] * f;
					for (int k = j + 1; k <= i - 1; k++) {
						g += V[k][j] * d[k];
						e[k] += V[k][j] * f;
					}
					e[j] = g;
				}
				f = 0.0;
				for (int j = 0; j < i; j++) {
					e[j] /= h;
					f += e[j] * d[j];
				}
				double hh = f / (h + h);
				for (int j = 0; j < i; j++) {
					e[j] -= hh * d[j];
				}
				for (int j = 0; j < i; j++) {
					f = d[j];
					g = e[j];
					for (int k = j; k <= i - 1; k++) {
						V[k][j] -= (f * e[k] + g * d[k]);
					}
					d[j] = V[i - 1][j];
					V[i][j] = 0.0;
				}
			}
			d[i] = h;
		}

		// Accumulate transformations.

		for (int i = 0; i < n - 1; i++) {
			V[n - 1][i] = V[i][i];
			V[i][i] = 1.0;
			double h = d[i + 1];
			if (h != 0.0) {
				for (int k = 0; k <= i; k++) {
					d[k] = V[k][i + 1] / h;
				}
				for (int j = 0; j <= i; j++) {
					double g = 0.0;
					for (int k = 0; k <= i; k++) {
						g += V[k][i + 1] * V[k][j];
					}
					for (int k = 0; k <= i; k++) {
						V[k][j] -= g * d[k];
					}
				}
			}
			for (int k = 0; k <= i; k++) {
				V[k][i + 1] = 0.0;
			}
		}
		for (int j = 0; j < n; j++) {
			d[j] = V[n - 1][j];
			V[n - 1][j] = 0.0;
		}
		V[n - 1][n - 1] = 1.0;
		e[0] = 0.0;
	}

	// Symmetric tridiagonal QL algorithm.

	static void tql2(mat & V, vec & d, vec & e) {
		int n = d.Size();

		//  This is derived from the Algol procedures tql2, by
		//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		//  Fortran subroutine in EISPACK.

		for (int i = 1; i < n; i++) {
			e[i - 1] = e[i];
		}
		e[n - 1] = 0.0;

		double f = 0.0;
		double tst1 = 0.0;
		double eps = pow(2.0, -52.0);
		for (int l = 0; l < n; l++) {

			// Find small subdiagonal element

			tst1 = MAX(tst1, fabs(d[l]) + fabs(e[l]));
			int m = l;
			while (m < n) {
				if (fabs(e[m]) <= eps*tst1) {
					break;
				}
				m++;
			}

			// If m == l, d[l] is an eigenvalue,
			// otherwise, iterate.

			if (m > l) {
				int iter = 0;
				do {
					iter = iter + 1;  // (Could check iteration count here.)

					// Compute implicit shift

					double g = d[l];
					double p = (d[l + 1] - g) / (2.0 * e[l]);
					double r = hypot2(p, 1.0);
					if (p < 0) {
						r = -r;
					}
					d[l] = e[l] / (p + r);
					d[l + 1] = e[l] * (p + r);
					double dl1 = d[l + 1];
					double h = g - d[l];
					for (int i = l + 2; i < n; i++) {
						d[i] -= h;
					}
					f = f + h;

					// Implicit QL transformation.

					p = d[m];
					double c = 1.0;
					double c2 = c;
					double c3 = c;
					double el1 = e[l + 1];
					double s = 0.0;
					double s2 = 0.0;
					for (int i = m - 1; i >= l; i--) {
						c3 = c2;
						c2 = c;
						s2 = s;
						g = c * e[i];
						h = c * p;
						r = hypot2(p, e[i]);
						e[i + 1] = s * r;
						s = e[i] / r;
						c = p / r;
						p = c * d[i] - s * g;
						d[i + 1] = h + s * (c * g + s * d[i]);

						// Accumulate transformation.

						for (int k = 0; k < n; k++) {
							h = V[k][i + 1];
							V[k][i + 1] = s * V[k][i] + c * h;
							V[k][i] = c * V[k][i] - s * h;
						}
					}
					p = -s * s2 * c3 * el1 * e[l] / dl1;
					e[l] = s * p;
					d[l] = c * p;

					// Check for convergence.

				} while (fabs(e[l]) > eps*tst1);
			}
			d[l] = d[l] + f;
			e[l] = 0.0;
		}

		// Sort eigenvalues and corresponding vectors.

		for (int i = 0; i < n - 1; i++) {
			int k = i;
			double p = d[i];
			for (int j = i + 1; j < n; j++) {
				if (d[j] < p) {
					k = j;
					p = d[j];
				}
			}
			if (k != i) {
				d[k] = d[i];
				d[i] = p;
				for (int j = 0; j < n; j++) {
					p = V[j][i];
					V[j][i] = V[j][k];
					V[j][k] = p;
				}
			}
		}
	}

	const Boolean_t eigen_decomposition(const mat & A, mat & V, vec & d) {
		int n = A.Rank();

		Boolean_t IsOk = (n > 0 && n == V.Rank() && n == d.Size());

		vec e(n);
		V = A;

		tred2(V, d, e);
		tql2(V, d, e);

		V.TransposeSelf();

		return IsOk;
	}

	//	Testing Mat and Vec for compiler errors
	void func(){
		vec3 v3(1, 2, 3), v3a(2, 3, 4), EigVals;
		vec4 v4(1, 2, 3, 4), EVal4;
		mat33 m3 = { { 1, 2, 3 }, { 2, 4, 5 }, { 3, 5, 6 } }, m3a, EigVecs;
		mat44 m4 = {
			{ 1, 2, 3, 4 },
			{ 2, 4, 5, 6 },
			{ 3, 5, 6, 7 },
			{ 4, 6, 7, 8 } }, EVec4;

		vec3 v3b = v3 + v3a;
		v3b += v3;
		v3 = v3.Normalize();
		v3.NormalizeSelf();
		v3 = Normalize(v3);

		const vec3 what = { 1, 2, 3 };

		// 	double a = what(0);

		v3 = m3.C(0);
		v3 = m3.At(0);

		v3 = m3[0];

		Boolean_t IsOk = m3.EigenSystem(EigVecs, EigVals);

		IsOk = m4.EigenSystem(EVec4, EVal4);

		m3a = m3;
		m3a = m3 + m3a;
		m3a = m3 * m3a;

		m3a = m3.ElemMult(m3a);
		m3a.ElemMultAssign(m3);

		m3a = -m3a;
		m3a.TransposeSelf();
		m3a = m3.Transpose();
		m3[0] = m3a[1].Normalize();
		double dist = m3[1].Distance(m3a[2]);
	}

}*/