#pragma once
#include<vector>

#include "CSMDATATYPES.h"

using namespace std;
class Point3d
{
public:
	Point3d(double _x, double _y, double _z);
	Point3d(CSM_Vec3_s point);
	Point3d();
	~Point3d();

	double x, y, z;
	vector<int> adjacentPoints;


	void addConnection(int connection);

	double operator[](int i) const{
		switch(i){
			case 0: return x; break;
			case 1: return y; break;
			case 2: return z; break;
			default: throw -1;
		}
	}
};

