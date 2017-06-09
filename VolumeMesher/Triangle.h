#pragma once
#include <vector>
#include "TECADDON.h"
#include "Point3d.h"
#include "CSMDATATYPES.h"
class Triangle
{
public:
	Triangle();
	Triangle(CSM_Vec3_s _a, CSM_Vec3_s _b, CSM_Vec3_s _c); 
	//Triangle(Point3d* _a, Point3d* _b, Point3d* _c);
	Triangle(int _a, int _b, int _c, std::vector<Point3d> pointList);
	~Triangle();

	bool isFront;
	void findCentroid();
	CSM_Vec3_s findNormal();
	int intersectWithRay(Point3d R_P1, Point3d R_P2) const;


	int aIdx = -1;
	int bIdx = -1;
	int cIdx = -1;
	Point3d* a;
	Point3d* b;
	Point3d* c;

	CSM_Vec3_s T_P0; //a
	CSM_Vec3_s T_P1; //b
	CSM_Vec3_s T_P2; //c
	CSM_Vec3_s centroid;
	CSM_Vec3_s normal;

	int operator[](const int & i) const{
		switch (i){
			case 0: return aIdx; break;
			case 1: return bIdx; break;
			case 2: return cIdx; break;
			default: throw - 1;
		}
	}

	const Boolean_t debug_save_tecplot_zone(const string & ZoneName, const int & ColorNum) const;

};

