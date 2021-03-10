
#include <iostream>
#include <cmath>
#include <vector>

#include <string>
#include <stdio.h>
#include <set>
#include <map>
#include <list>

#include <armadillo>

using namespace std;
using namespace arma;
// 
// /*
//        v2
//       /|\
//      / | \
//   e1/  |  \e3
//    /   |e2 \
// v1/____|____\v3
//   \ e6 |    /
//    \   |   /
//   e4\  |  /e5
//      \ | /
//       \|/
//        v4
// 
// * V{1-4} are the vertices and
// * e{1-6} are the edge lengths
// */
// 
// /*
// *	Calculate tetrahedron volume using squared edge lengths.
// *	Taken from    http://keisan.casio.com/exec/system/1329962711
// */
// double const TetVolume(vector<double> const & e2){
// 
// 	double v = 0.0069444444444444 * (
// 		e2[0] * e2[4] * (e2[1] + e2[2] + e2[3] + e2[5] - e2[0] - e2[4])
// 		+ e2[1] * e2[5] * (e2[0] + e2[2] + e2[3] + e2[4] - e2[1] - e2[5])
// 		+ e2[2] * e2[3] * (e2[0] + e2[1] + e2[4] + e2[5] - e2[2] - e2[3])
// 		- e2[0] * e2[1] * e2[3]
// 		- e2[1] * e2[2] * e2[4]
// 		- e2[0] * e2[2] * e2[5]
// 		- e2[3] * e2[4] * e2[5]
// 		);
// 
// 	return sqrt(v);
// }
// 
// /*
// *	Calculate tetrahedron volume given vertices
// *	Taken from https://en.wikipedia.org/wiki/Tetrahedron#Volume
// */
// 
// /*
// *	Calculate total volume of all tets formed between
// *	the vertices of a hexahedron and an internal point.
// *	V is list of vertices, where 0 (i) is the internal point.
// */
// double const TetVolume(vec3 const & a, vec3 const & b, vec3 const & c, vec3 const & d){
// 	return abs(dot(a - d, cross(b - d, c - d))) / 6.0;
// }
// 
// 
// /*
//    5-------6	    e-------f
//   /|      /|	   /|      /|
//  / |  *0 / |	  / |  *i / |
// 1-------2  |	 a-------b  |
// |  7----|--8	 |  g----|--h
// | /     | /	 | /     | /
// 3-------4	     c-------d
// 
// 
//      5--------------6
//     /|             /|
//    / |            / |
//   /  |           /  |
//  /   |   *0     /   |
// 1--------------2    |
// |    |         |    |
// |    7---------|----8
// |   /          |   /
// |  /           |  /
// | /            | /
// 3--------------4
// 
// */
// // 26 edges:
// //	12 for the edges
// //	6 for face crossing
// //	8 for interior point to vertices
// vector<vector<int> > EdgeInds = {
// 	{ 1, 2 },	// 0
// 	{ 1, 3 },	// 1
// 	{ 1, 5 },	// 2
// 	{ 2, 4 }, 	// 3
// 	{ 2, 6 },	// 4
// 	{ 3, 4 }, 	// 5
// 	{ 3, 7 },	// 6
// 	{ 5, 6 },	// 7
// 	{ 5, 7 },	// 8
// 	{ 4, 8 }, 	// 9
// 	{ 6, 8 }, 	// 10
// 	{ 7, 8 },	// 11
// 
// 	{ 1, 4 },	// 12
// 	{ 1, 6 },	// 13
// 	{ 1, 7 },	// 14
// 	{ 2, 8 },	// 15
// 	{ 3, 8 },	// 16
// 	{ 5, 8 },	// 17
// 
// 	{ 0, 1 },	// 18
// 	{ 0, 2 },	// 19
// 	{ 0, 3 },	// 20
// 	{ 0, 4 },	// 21
// 	{ 0, 5 },	// 22
// 	{ 0, 6 },	// 23
// 	{ 0, 7 },	// 24
// 	{ 0, 8 }	// 25
// };
// // indices of all the tets with 0 as vertex
// vector<vector<int> > TetInds = {
// 	{ 0, 1, 2, 3 }, { 0, 2, 3, 4 },
// 	{ 0, 1, 3, 5 }, { 0, 3, 5, 7 },
// 	{ 0, 1, 2, 5 }, { 0, 2, 5, 6 },
// 	{ 0, 2, 4, 6 }, { 0, 4, 6, 8 },
// 	{ 0, 3, 4, 7 }, { 0, 4, 7, 8 },
// 	{ 0, 5, 6, 7 }, { 0, 6, 7, 8 }
// };
// double const HexahedronInternalPointTetVolume(vector<vec3> const & V){
// 
// 	double vol = 0.0;
// 	for (auto & vi : TetInds) vol += TetVolume(V[vi[0]], V[vi[1]], V[vi[2]], V[vi[3]]);
// 
// 	return vol;
// }
// 
// /*
// *	Calculate volume of parallelpiped using the lattice vector that defines it
// */
// double const ParallepipedVolume(vector<vec3> const & LV){
// 
// 	return  abs(dot(LV[0], cross(LV[1], LV[2])));
// }
// 
// double const ParallepipedVolume(mat33 const & LV){
// 	return det(LV);
// }
// 
// bool const ParallelpidedPointIsInternal(mat33 const & LV, vec3 const & Origin, vec3 const & Pt){
// 
// 	vector<vec3> V(9);
// 	V[0] = Pt;
// 	V[1] = Origin;
// 	V[2] = Origin + LV.col(0);
// 	V[3] = Origin + LV.col(1);
// 	V[4] = V[2] + LV.col(1);
// 	V[5] = Origin + LV.col(2);
// 	V[6] = V[5] + LV.col(0);
// 	V[7] = V[5] + LV.col(1);
// 	V[8] = V[6] + LV.col(1);
// 
// 	double PtTetVol = HexahedronInternalPointTetVolume(V);
// 	double FullVol = ParallepipedVolume(LV);
// 
// 	return (PtTetVol <= FullVol);
// }
// 
// void TestParallelpipedInteriorCheck(){
// 	
// 	mat33 lv;
// 	lv <<
// 		2 << 0 << 0 << endr <<
// 		0 << 2 << 0 << endr <<
// 		0 << 0 << 2 << endr;
// 
// 	vec3 origin, pt;
// 	origin << 0 << 0 << 0;
// 	pt << 1 << 1 << 1;
// 
// 	vector<vec3> V(9);
// 	V[0] = pt;
// 	V[1] = origin;
// 	V[2] = origin + lv.col(0);
// 	V[3] = origin + lv.col(1);
// 	V[4] = V[2] + lv.col(1);
// 	V[5] = origin + lv.col(2);
// 	V[6] = V[5] + lv.col(0);
// 	V[7] = V[5] + lv.col(1);
// 	V[8] = V[6] + lv.col(1);
// 
// 	cout << "lattice vector:" << endl << lv << endl << "volume by direct calculation: " << ParallepipedVolume(lv) << endl << endl
// 		<< "origin:" << endl << origin << endl << "internal point: " << endl << pt << endl << "volume by sum of tets: " << HexahedronInternalPointTetVolume(V) << endl << endl;
// 
// 	pt << -.5 << .5 << .5;
// 	bool IsInternal = ParallelpidedPointIsInternal(lv, origin, pt);
// 	cout << "internal point: " << endl << pt << endl << "Is internal? " << boolalpha << IsInternal << endl;
// 
// 	pt << 2.5 << .5 << .5;
// 	IsInternal = ParallelpidedPointIsInternal(lv, origin, pt);
// 	cout << "internal point: " << endl << pt << endl << "Is internal? " << boolalpha << IsInternal << endl;
// 
// 	pt << 2 << 2 << 2;
// 	IsInternal = ParallelpidedPointIsInternal(lv, origin, pt);
// 	cout << "internal point: " << endl << pt << endl << "Is internal? " << boolalpha << IsInternal << endl;
// }
// 
// void TestTetVolume(){
// 	mat33 lv;
// 	lv <<
// 		1 << 0 << 0 << endr <<
// 		0 << 1 << 0 << endr <<
// 		0 << 0 << 1 << endr;
// 
// 	vec3 origin;
// 	origin << 0 << 0 << 0;
// 	cout << TetVolume(origin, lv.col(0), lv.col(1), lv.col(2)) << endl;
// }
// 
// void TransformationTest(){
// 	vec3 origin;
// 	origin << 2 << 2 << 2;
// 	vector<vec3> a;
// 	for (double z = 0; z <= 1; z += 0.5){
// 		for (double y = 0; y <= 1; y += 0.5){
// 			for (double x = 0; x <= 1; x += 0.5){
// 				a.push_back(vec3());
// 				a.back() << x << y << z;
// 			}
// 		}
// 	}
// 
// 	mat b(3, a.size());
// 	for (int i = 0; i < a.size(); ++i) b.col(i) = a[i] + origin;
// 
// // 	for (auto & i : A) cout << endl << i << endl;
// 
// 	cout << "original:" << endl;
// 	cout << endl << b << endl;
// 
// 	mat33 lv;
// 	lv <<
// 		0		<< 6.99925	<< 6.99925	<< endr <<
// 		6.99925 << 0		<< 6.99925	<< endr <<
// 		6.99925 << 6.99925	<< 0		<< endr;
// 
// 	lv = mat33(normalise(lv));
// 
// 	mat c = lv * b;
// 
// 	cout << endl << c << endl;
// 
// 	c = lv.i() * c;
// 
// 	cout << endl << c << endl;
// }

void SetBinarySearchInterpolationTesting(bool flip){
	double val = 4.5;
	std::vector<double> a = { 1,3,5,7,9 };
	if (flip) a = std::vector<double>(a.rbegin(), a.rend());

	std::vector<double> b = { 10,20,40,70,110 };

	std::set<double> s(a.begin(), a.end());

	auto lBound = s.lower_bound(val);

	int ind = std::distance(s.begin(), lBound) - 1;

	if (a.back() < a[0])
		ind = a.size() - 2 - ind;

	int ind2 = ind + 1;

	double weight = (val - a[ind]) / (a[ind2] - a[ind]);

	if (a.back() < a[0])
		ind2 = ind - 1;

	double outVal = b[ind] + weight * (b[ind2] - b[ind]);


	return;
}

void SetIteratingWhileModifyingTesting(){
	std::set<int> a = { 1,3,5,7,9 };
	auto it = a.begin();
	
	cout << *it << endl;

	it++;

	cout << *it << endl;

	a.insert(2);

	cout << *it << endl;

	it--;

	cout << *it << endl;

	a.insert(4);

	it++;
	it++;

	cout << *it;

	return;
}

void ListIteratingWhileModifyingTesting() {
	std::list<int> a = { 1,3,5,7,9 };
	auto it = a.begin();

	cout << *it << endl;

	it++;

	cout << *it << endl;

	a.insert(it,2);

	cout << *it << endl;

	it--;

	cout << *it << endl;

	a.insert(it,4);

	for (auto i : a)
		cout << i << ",";

	return;
}

vector<int> IntVecConcat(vector<int> const & a, vector<int> const & b) {
	int MinSqrDist = abs(a.back() - b[0]);
	int MinPair = 1;
	int TempSqrDist = abs(a.back() - b.back());
	if (MinSqrDist > TempSqrDist) {
		MinSqrDist = TempSqrDist;
		MinPair = 2;
	}
	TempSqrDist = abs(a[0] - b.back());
	if (MinSqrDist > TempSqrDist) {
		MinSqrDist = TempSqrDist;
		MinPair = 3;
	}
	TempSqrDist = abs(a[0] - b[0]);
	if (MinSqrDist > TempSqrDist) {
		MinSqrDist = TempSqrDist;
		MinPair = 4;
	}

	vector<int> out = a;

	switch (MinPair) {
	case 1:
		out.insert(out.end(), b.cbegin(), b.cend());
		break;
	case 2:
		out.insert(out.end(), b.crbegin(), b.crend());
		break;
	case 3:
		out.insert(out.begin(), b.cbegin(), b.cend());
		break;
	case 4:
		out.insert(out.begin(), b.crbegin(), b.crend()-1);
		break;
	}

	return out;
}

void VectorConcatenateTesting(){
	vector<int> aForward = { 1,2,3,4,5 },
		bForward = { 11,12,13,14,15 };
	vector<int> aBackward(aForward.rbegin(), aForward.rend()),
		bBackward(bForward.rbegin(), bForward.rend());

	auto c1 = IntVecConcat(aForward, bForward),
		c2 = IntVecConcat(aForward, bBackward),
		c3 = IntVecConcat(aBackward, bBackward),
		c4 = IntVecConcat(aBackward, bForward);

	return;
}

void PairSortTesting(){
	vector<std::pair<double, int> > a;
	a.push_back(std::make_pair(2.3, 1));
	a.push_back(std::make_pair(5.4, 2));
	a.push_back(std::make_pair(1.2, 3));
	a.push_back(std::make_pair(1.3, 4));

	std::set<std::pair<double, int> > b(a.begin(), a.end());

	std::make_heap(a.begin(), a.end());



	return;
}

double TriArea(vec3 const & p1, vec3 const & p2, vec3 const & p3) {
	vec3 ABC;

	ABC[0] = p1[2] * (p2[1] - p3[1])
		+ p2[2] * (p3[1] - p1[1])
		+ p3[2] * (p1[1] - p2[1]);

	ABC[1] = p1[0] * (p2[2] - p3[2])
		+ p2[0] * (p3[2] - p1[2])
		+ p3[0] * (p1[2] - p2[2]);

	ABC[2] = p1[1] * (p2[0] - p3[0])
		+ p2[1] * (p3[0] - p1[0])
		+ p3[1] * (p1[0] - p2[0]);

	return 0.5 * norm(ABC);
}

vector<string> SplitString(string const &s, string const & delim = " ", bool RemoveAllBlanks = true, bool RemoveBlankAtEnd = true) {
	vector<string> tokens;
	auto start = 0U;
	auto end = s.find(delim);
	if (RemoveAllBlanks) {
		RemoveBlankAtEnd = true;
		while (end != std::string::npos) {
			if (end - start > 0) {
				tokens.push_back(s.substr(start, end - start));
			}
			start = end + delim.length();
			end = s.find(delim, start);
		}
	}
	else {
		while (end != std::string::npos) {
			tokens.push_back(s.substr(start, end - start));
			start = end + delim.length();
			end = s.find(delim, start);
		}
	}
	if (!RemoveBlankAtEnd || s.length() - start > 0)
		tokens.push_back(s.substr(start, s.length() - start));
	// 	if (tokens.size() == 0) tokens.push_back(s);
	return tokens;
}

bool StringIsInt(string const & s) {
	bool IsInt = true;
	for (int i = 0; i < s.length() && IsInt; ++i)
		IsInt = isdigit(s[i]);

	return IsInt;
}

bool StringIsFloat(string const & s) {
	for (auto const & i : s) if (!isdigit(i) && i != '.') return false;

	return true;
}





int main(char *vargv, int argc){
// 	SetBinarySearchInterpolationTesting(false);
// 	SetBinarySearchInterpolationTesting(true);
// 	ListIteratingWhileModifyingTesting();

// 	VectorConcatenateTesting();

// 	int a = 5, b = 10;
// 	int c = (a + b) * 0.5;

// 	PairSortTesting();

	vec3 p1, p2, p3;
	p1 << 0 << 0 << 0;
	p2 << 1 << 0 << 0;
	p3 << 0 << 1 << 0;

	double a = TriArea(p1,p2,p3);

	cout << endl << endl;
	system("pause");

	return 0;
}