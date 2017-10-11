#include "TECADDON.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif
#include "GLOBAL.h"

#include <vector>
#include <armadillo>

#include "GSM_GEOMETRY.h"
#include "CSM_DATA_TYPES.h"

using namespace arma;

using std::vector;

static const int ArrowRotationSteps = 10;

const bool CSMArrow(const vec3 & Origin, 
	const vec3 & Dir, 
	const double & Length, 
	const double & Radius, 
	const double & ArrowheadLengthRatio,
	const double & ArrowheadRadiusRatio,
	vector<vec3> & Nodes,
	vector<vector<int> > & ElemList)
{
	int NumNodes = ArrowRotationSteps * 7 + 1; 
	// RotSteps for bottom of cylinder, then minor and major arrowhead, then 1 for center of bottom and RotSteps for tip
	// The bottom, minor, major, and tip nodes are doubled in order to produce proper shading on arrow
	int NumElems = ArrowRotationSteps * 6; // RotSteps triangles for bottom and top and RotSteps quads (two triangles) for cylinder sides and bottom of cap
	Nodes.resize(NumNodes); 
	ElemList.resize(NumElems, vector<int>(3));

	double InnerRadius = Radius * (1.0 - ArrowheadRadiusRatio);
	double RotStep = 2.0 * PI / static_cast<double>(ArrowRotationSteps);

	// First get normal vectors
	vec3 v1 = normalise(Dir);
	vec3 v2 = v1;
	v2[2] += 1;
	vec3 v3 = normalise(cross(v1, v2));

	v2 = normalise(cross(v1, v3));

	// Cylinder bottom points
	int NodeInd = 0;

	for (int l = 0; l < 2; ++l){
		for (int i = 0; i < ArrowRotationSteps; ++i){
			Nodes[NodeInd++] = Origin + Rotate(v2 * InnerRadius, static_cast<double>(i)* RotStep, v1);
		}
	}

	// Arrowhead inner points
	vec3 ArrowheadBottom = Origin + v1 * Length * (1.0 - ArrowheadLengthRatio);

	for (int l = 0; l < 2; ++l){
		for (int i = 0; i < ArrowRotationSteps; ++i){
			Nodes[NodeInd++] = ArrowheadBottom + Rotate(v2 * InnerRadius, static_cast<double>(i)* RotStep, v1);
		}
	}

	// Arrowhead outer points
	for (int l = 0; l < 2; ++l){
		for (int i = 0; i < ArrowRotationSteps; ++i){
			Nodes[NodeInd++] = ArrowheadBottom + Rotate(v2 * Radius, static_cast<double>(i)* RotStep, v1);
		}
	}

	// Arrow tip points
	for (int i = 0; i < ArrowRotationSteps; ++i){
		Nodes[NodeInd++] = Origin + v1 * Length;
	}

	//Cylinder Bottom center
	int BottomCenter = NodeInd++;
	Nodes[BottomCenter] = Origin;

	// Now make the connectivity list.
	// First the cylinder bottom
	int ElemInd = 0, ei = 0;

	for (int i = 0; i < ArrowRotationSteps; ++i){
		ElemList[ElemInd++] = { BottomCenter, ei, i < ArrowRotationSteps - 1 ? ei + 1 : 0 };
		ei++;
	}

// 	ei = 0;

	// Now the cylinder sides, with each side as a pair of triangles
	// This happens twice, first for the cylinder, then for the bottom of the cap
	for (int l = 0; l < 2; ++l){
		for (int i = 0; i < ArrowRotationSteps - 1; ++i){
			ElemList[ElemInd++] = { ei, ei + 1, ei + ArrowRotationSteps };
			ei++;
			ElemList[ElemInd++] = { ei, ei + ArrowRotationSteps, ei + ArrowRotationSteps - 1 };
		}
		ElemList[ElemInd++] = { ei, ei - ArrowRotationSteps + 1, ei + ArrowRotationSteps };
		ei++;
		ElemList[ElemInd++] = { ei - ArrowRotationSteps, ei, ei + ArrowRotationSteps - 1 };
		ei += ArrowRotationSteps;
	}

// 	ei -= ArrowRotationSteps;

// 	ei += ArrowRotationSteps;

	// Now the top of the arrowhead
	for (int i = 0; i < ArrowRotationSteps - 1; ++i){
		// 		ElemList[ElemInd++] = { Tip, ei, i < ArrowRotationSteps - 1 ? ei + 1 : ei - ArrowRotationSteps + 1 };
		ElemList[ElemInd++] = { ei, ei + 1, ei + ArrowRotationSteps};
		ei++;
	}
	ElemList[ElemInd++] = { ei, ei + 1 - ArrowRotationSteps, ei + ArrowRotationSteps };

	return true;
}