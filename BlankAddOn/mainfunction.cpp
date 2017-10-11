
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

//#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include <armadillo>
using namespace arma;

using std::vector;
using std::stringstream;
using std::string;
using std::to_string;


#include "MAINFUNCTION.h"




Boolean_t StatusUpdate(unsigned int CurrentNum, unsigned int TotalNum, const char* ProgresssText){
	unsigned int Percent = (int)((double)CurrentNum / (double)TotalNum * 100.);

	std::stringstream ss;
	ss << ProgresssText << "  (" << Percent << "% Complete)";

	TecUtilDialogSetPercentDoneText(ss.str().c_str());
	Boolean_t IsOk = TecUtilDialogCheckPercentDone(Percent);

	return IsOk;
}

vector<double> split(const string &s, char delim) {
	stringstream ss(s);
	string item;
	vector<double> tokens;
	while (getline(ss, item, delim)) {
		tokens.push_back(std::stod(item));
	}
	return tokens;
}

const mat LoadFile(const string & Path, const char & Delim){
	mat Points(1, 3);
	std::ifstream File(Path.c_str());

	if (File.is_open()){
		string Line;
		while (!File.eof()){
			std::getline(File, Line);
			vector<double> tmp = split(Line, Delim);
			if (tmp.size() > 0)
				Points = join_cols(Points, mat(tmp).t());
		}
		File.close();
	}
	else{
		TecUtilDialogErrMsg("Failed to open file");
	}

	return Points.tail_rows(Points.n_rows - 1);
}

const mat getPoints(const mat & Edge, const int & NumPts){
	mat Pts(NumPts, 3);

	mat EdgeDiff = Edge.tail_rows(Edge.n_rows - 1) - Edge.head_rows(Edge.n_rows - 1);

	vec Dist = join_cols(mat(vector<double>({ 0. })), sqrt(square(EdgeDiff.col(0)) + square(EdgeDiff.col(1)) + square(EdgeDiff.col(2))));

	for (int i = 1; i < Dist.n_elem; ++i)
		Dist(i) += Dist(i - 1);

	vec PtDist = linspace(0, Dist(Dist.n_elem - 1), NumPts);

	Pts.row(0) = Edge.row(0);
	Pts.row(NumPts - 1) = Edge.row(Edge.n_rows - 1);

	int Ind = 0;
	for (int i = 1; i < NumPts - 1; ++i){
		double d = PtDist(i);
		for (int j = Ind; j < Dist.n_elem - 1; ++j){
			if (Dist(j) <= d && Dist(j + 1) > d){
				Ind = j;
				break;
			}
		}
		double Rem = d - Dist(Ind);
		double Tot = Dist(Ind + 1) - Dist(Ind);
		double Rat = Rem / Tot;
		Pts.row(i) = Edge.row(Ind) + (Edge.row(Ind + 1) - Edge.row(Ind)) * Rat;
	}

	return Pts;
}

mat cubTrans(const mat & e1, const mat & e2, const mat & e3){

	mat Rpt;
	Rpt << 0 << 0 << 0 << endr
		<< 1. / 3. << 0 << 0 << endr
		<< 2. / 3. << 0 << 0 << endr

		<< 1 << 0 << 0 << endr
		<< 0 << 1. / 3. << 0 << endr
		<< 0 << 2. / 3. << 0 << endr

		<< 0 << 1 << 0 << endr
		<< 0 << 0 << 1. / 3. << endr
		<< 0 << 0 << 2. / 3. << endr

		<< 0 << 0 << 1 << endr
		<< 2. / 3. << 1. / 3. << 0 << endr
		<< 1. / 3. << 2. / 3. << 0 << endr

		<< 2. / 3. << 0 << 1. / 3. << endr
		<< 1. / 3. << 0 << 2. / 3. << endr
		<< 0 << 2. / 3. << 1. / 3. << endr

		<< 0 << 1. / 3. << 2. / 3. << endr
		<< 1. / 3. << 1. / 3. << 0 << endr
		<< 1. / 3. << 0 << 1. / 3. << endr

		<< 0 << 1. / 3. << 1. / 3. << endr
		<< 1. / 3. << 1. / 3. << 1. / 3.;

	// 		mat Rpt({
	// 			{ 0, 0, 0 }, { 1. / 3., 0, 0 }, { 2. / 3., 0, 0 },
	// 			{ 1, 0, 0 }, { 0, 1. / 3., 0 }, { 0, 2. / 3., 0 },
	// 			{ 0, 1, 0 }, { 0, 0, 1. / 3. }, { 0, 0, 2. / 3. },
	// 			{ 0, 0, 1 }, { 2. / 3., 1. / 3., 0 }, { 1. / 3., 2. / 3., 0 },
	// 			{ 2. / 3., 0, 1. / 3. }, { 1. / 3., 0, 2. / 3. }, { 0, 2. / 3., 1. / 3. },
	// 			{ 0, 1. / 3., 2. / 3. }, { 1. / 3., 1. / 3., 0 }, { 1. / 3., 0, 1. / 3. },
	// 			{ 0, 1. / 3., 1. / 3. }, { 1. / 3., 1. / 3., 1. / 3. }
	// 		});

	mat pt1 = getPoints(e1, 4),
		pt2 = getPoints(e2, 4),
		pt3 = getPoints(e3, 4);

	rowvec xx = getPoints(e1, 3).row(1),
		yy = getPoints(e2, 3).row(1),
		zz = getPoints(e3, 3).row(1);

	mat pt(20, 3);
	int pti = 0;
	for (int i = 0; i < 4; ++i)
		pt.row(pti++) = pt1.row(i);
	for (int i = 1; i < 4; ++i)
		pt.row(pti++) = pt2.row(i);
	for (int i = 1; i < 4; ++i)
		pt.row(pti++) = pt3.row(i);
	pt.row(pti++) = (pt1.row(3) * (2. / 3.)) + (pt2.row(3) * (1. / 3.));
	pt.row(pti++) = (pt1.row(3) * (1. / 3.)) + (pt2.row(3) * (2. / 3.));

	pt.row(pti++) = (pt1.row(3) * (2. / 3.)) + (pt3.row(3) * (1. / 3.));
	pt.row(pti++) = (pt1.row(3) * (1. / 3.)) + (pt3.row(3) * (2. / 3.));

	pt.row(pti++) = (pt2.row(3) * (2. / 3.)) + (pt3.row(3) * (1. / 3.));
	pt.row(pti++) = (pt2.row(3) * (1. / 3.)) + (pt3.row(3) * (2. / 3.));

	pt.row(pti++) = (xx * (1. / 2.)) + (yy * (1. / 2.));
	pt.row(pti++) = (xx * (1. / 2.)) + (zz * (1. / 2.));
	pt.row(pti++) = (yy * (1. / 2.)) + (zz * (1. / 2.));

	pt.row(pti++) = (pt1.row(3) * (1. / 3.)) + (pt2.row(3) * (1. / 3.)) + (pt3.row(3) * (1. / 3.));

	vec C1 = Rpt.col(0), C2 = Rpt.col(1), C3 = Rpt.col(2);
	vec C12 = square(C1), C22 = square(C2), C32 = square(C3);

	mat A(Rpt.n_rows, Rpt.n_rows);
	pti = 0;

	A.col(pti++) = vec(Rpt.n_rows, fill::ones);
	A.col(pti++) = C1;
	A.col(pti++) = C2;
	A.col(pti++) = C3;
	A.col(pti++) = C12;
	A.col(pti++) = C22;
	A.col(pti++) = C32;
	A.col(pti++) = C1%C2;
	A.col(pti++) = C1%C3;
	A.col(pti++) = C2%C3;
	A.col(pti++) = pow(C1, 3);
	A.col(pti++) = pow(C2, 3);
	A.col(pti++) = pow(C3, 3);
	A.col(pti++) = C12%C2;
	A.col(pti++) = C12%C3;
	A.col(pti++) = C22%C1;
	A.col(pti++) = C22%C3;
	A.col(pti++) = C32%C1;
	A.col(pti++) = C32%C2;
	A.col(pti++) = C1%C2%C3;

	mat abc(pt.n_rows, 3);

	for (int i = 0; i < 3; ++i){
		abc.col(i) = solve(A, pt.col(i));
	}

	return abc.t();
}

const mat cubTrans_func(mat & abc,
	vec & s,
	vec & t,
	vec & u)
{
	mat xyz(s.n_elem, 3);
	for (int i = 0; i < 3; ++i){
		xyz.col(i) =
			(s*abc(i, 1) + abc(i, 0))
			+ t*abc(i, 2)
			+ u*abc(i, 3)
			+ square(s)*abc(i, 4)
			+ square(t)*abc(i, 5)
			+ square(u)*abc(i, 6)
			+ s%t*abc(i, 7)
			+ s%u*abc(i, 8)
			+ t%u*abc(i, 9)
			+ pow(s, 3)*abc(i, 10)
			+ pow(t, 3)*abc(i, 11)
			+ pow(u, 3)*abc(i, 12)
			+ square(s) % t*abc(i, 13)
			+ square(s) % u*abc(i, 14)
			+ square(t) % s*abc(i, 15)
			+ square(t) % u*abc(i, 16)
			+ square(u) % s*abc(i, 17)
			+ square(u) % t*abc(i, 18)
			+ s%t%u*abc(i, 19);
	}
	return xyz;
}

const double cubJacobian(mat & abc, const double & s, const double & t, const double & u){
	vector<vec> dstu({ vector<double>({ 0., 1., 0., 0., s * 2, 0., 0., t, u, 0.,
		(s*s) * 3, 0., 0., s*t * 2, s*u * 2, t*t, 0., u*u, 0., t*u }),
		vector<double>({ 0., 0., 1., 0., 0., t * 2, 0., s, 0., u,
		0., (t*t) * 3, 0., s*s, 0., t*s * 2, t*u * 2, 0., u*u, s*u }),
		vector<double>({ 0., 0., 0., 1., 0., 0., u * 2, 0., s, t,
		0., 0., (u*u) * 2, 0., s*s, 0., t*t, u*s * 2, u*t * 2, s*t }) });

	mat33 J;
	for (unsigned int i = 0; i < 3; ++i){
		for (unsigned int j = 0; j < 3; ++j){
			J.at(i, j) = dot(abc.row(i), dstu[j]);
		}
	}
	return det(J);
}

void rquad(const int & N, const double & k, vec & x, vec & w)
{
	double k1 = k + 1, k2 = k + 2;

	vec n = linspace<vec>(1, N, N);

	vec nnk = 2 * n + k;

	rowvec A = join_rows(mat(vector<double>({ k / k2 })), ((ones<vec>(N) * (k*k)) / (nnk % (nnk + 2))).t());

	n = n.tail(N - 1);
	nnk = nnk.tail(N - 1);

	double B1 = 4. * k1 / (k2*k2*(k + 3));

	vec nk = n + k, nnk2 = square(nnk);
	vec B = 4. * square(n%nk) / (square(nnk2) - nnk2);

	mat ab = join_rows(A.t(), join_cols(vec(vector<double>({ pow(2, k1) / k1, B1 })), B));

	vec s = sqrt(ab(span(1, N - 1), 1));

	mat V;
	vec X;

	eig_sym(X, V, diagmat(ab(span(0, N - 1), 0)) + diagmat(s, -1) + diagmat(s, 1));

	x = (X + 1) / 2;

	w = pow(0.5, k1) * ab(0, 1) * square(V.row(0).t());

	return;
}

const vector<vec> tetraquad(const int & N, const mat & Verts)
{

	vector<vec> q(3), w123(3);

	for (int i = 0; i < 3; ++i)
		rquad(N, 2 - i, q[i], w123[i]);

	int N2 = N*N;
	int N3 = N2*N;
	vec q1(N3), q2(N3), q3(N3);

	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N; ++j){
			for (int k = 0; k < N; ++k){
				int Ind = i*N2 + j*N + k;
				q1[Ind] = q[0][j];
				q2[Ind] = q[1][k];
				q3[Ind] = q[2][i];
			}
		}
	}

	vec x = 1 - q1,
		y = (1 - q2) % q1,
		z = q1 % q2 % q3;

	vec w = reshape(reshape(w123[1] * w123[0].t(), N2, 1) * w123[2].t(), N3, 1);
	// 		mat c = mat({ { 1, 0, 0, 0 }, { -1, 1, 0, 0 }, { -1, 0, 1, 0 }, { -1, 0, 0, 1 } }) * Verts;
	mat c = reshape(mat(vector<double>({ 1, 0, 0, 0, -1, 1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 1 })), 4, 4).t() * Verts;
	vec W = fabs(det(c.rows(1, 3))) * w;

	mat XYZ = join_rows(ones<vec>(N3), join_rows(x, join_rows(y, z))) * c;

	return vector<vec>({ XYZ.col(0), XYZ.col(1), XYZ.col(2), W });
}



int GetWeightsPoints(){
	int N = 20;
	vector<int> NumElemsList = { 80, 320, 1280, 5120 };

	mat Verts;
	Verts << 0 << 0 << 0 << endr
		<< 1 << 0 << 0 << endr
		<< 0 << 1 << 0 << endr
		<< 0 << 0 << 1;
	// 		mat Verts = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

	vector<vec> stuW = tetraquad(N, Verts);

	for (auto NumElems : NumElemsList){
		int NumW = stuW[3].n_elem;
		vec WW(NumW*NumElems, fill::zeros);
		mat PP(NumW*NumElems, 3, fill::zeros);

		for (int ElemNum = 0; ElemNum < NumElems; ++ElemNum){
			vector<mat> Points(3);
			for (int i = 0; i < 3; ++i){
				Points[i] = LoadFile(string("P:\\Archive\\" + to_string(NumElems) + "\\Elem_" + to_string(ElemNum) + "_" + to_string(i) + ".csv"), ' ');
			}

			mat abc = cubTrans(Points[0], Points[1], Points[2]);
			mat xyz = cubTrans_func(abc, stuW[0], stuW[1], stuW[2]);
			vec WTemp(NumW);

			for (int i = 0; i < NumW; ++i){
				WTemp(i) = stuW[3](i) * cubJacobian(abc, stuW[0](i), stuW[1](i), stuW[2](i));
			}

			WW.subvec(ElemNum * NumW, (ElemNum + 1) * NumW - 1) = WTemp;
			PP.rows(ElemNum * NumW, (ElemNum + 1) * NumW - 1) = xyz;
		}

		TecUtilDialogMessageBox(string("n = " + to_string(N) + ", " + to_string(NumElems) + " elements, sum = " + to_string(sum(WW))).c_str(), MessageBoxType_Information);
	}

	return 0;
}

int Number = 1;

const string AuxDataMakeStringValidName(string Str){
	if (Str.length() <= 0) return "Blank_Name";

	if (!(Str[0] == 95 || (Str[0] >= 65 && Str[0] <= 90) || (Str[0] >= 97 && Str[0] <= 122)))
		Str = "_" + Str;

	for (char & Chr : Str){
		if (!(Chr == 95
			|| Chr == 46
			|| (Chr >= 48 && Chr <= 57)
			|| (Chr >= 65 && Chr <= 90)
			|| (Chr >= 97 && Chr <= 122)))
		{
			Chr = '_';
		}
	}

	return Str;
}

const Boolean_t AuxDataDataSetSetItem(const string & AuxDataName, const string & AuxDataValue)
{
	AuxData_pa TempAuxData = TecUtilAuxDataDataSetGetRef();
	Boolean_t IsOk = VALID_REF(TempAuxData);
	if (IsOk){
		IsOk = TecUtilAuxDataSetStrItem(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), AuxDataValue.c_str(), TRUE);
	}
	return IsOk;
}

const Boolean_t AuxDataZoneSetItem(const int & ZoneNum, const string & AuxDataName, const string & AuxDataValue)
{
	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	Boolean_t IsOk = VALID_REF(TempAuxData);
	if (IsOk){
		IsOk = TecUtilAuxDataSetStrItem(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), AuxDataValue.c_str(), TRUE);
	}
	return IsOk;
}

void MainFunction(){
	TecUtilDialogMessageBox("Just your run-of-the-mill add-on skeleton!", MessageBox_Warning);
	AuxDataDataSetSetItem("TempAuxData", to_string(Number));
	AuxDataZoneSetItem(1, "TempAuxData", to_string(Number++));
// 	return IsOk;
// 	GetWeightsPoints();
}