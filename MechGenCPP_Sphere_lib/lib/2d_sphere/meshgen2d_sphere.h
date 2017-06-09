/*
 * meshgen2d_sphere.h
 *
 *  Created on: Dec 13, 2014
 *      Author: Haiiro
 */

#ifndef LIB_2D_SPHERE_MESHGEN2D_SPHERE_H_
#define LIB_2D_SPHERE_MESHGEN2D_SPHERE_H_

#include <vector>
#include "mesh.h"

using std::vector;

enum MeshStatus_e{
	SUCCESS = 0,
	FAIL_NEEDREFINEMENT,
	FAIL_INVALIDCONSTRAINT,
	SUCCESS_NOTCONVERGED
};

MeshStatus_e meshgen2D_sphere(double Radius,
	int Level,
	vector<point> & ConstrainedVertices,
	vector<int> & MovedPointNums,
	point *& pIn,
	triangle *& tIn,
	int **& e,
	int &NumPts,
	int &NumTri,
	int &NumEdges);


#endif /* LIB_2D_SPHERE_MESHGEN2D_SPHERE_H_ */
