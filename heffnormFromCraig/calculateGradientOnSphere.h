#pragma once

#include <vector>

#include "SimpleTypes.h"
#include "Vec3.h"
#include "TriNodes.h"

namespace tpcsm {

typedef enum _GradientCalcMethod
{
    GradientCalcMethod_ThreeNodesAndCellCenter,
    GradientCalcMethod_FourCellCenters,
    GradientCalcMethod_AllConnectedCellCenters
} GradientCalcMethod_e;

// returns true if successful, false otherwise
    bool calculateGradientOnSphere(
    std::vector<Vec3> const&     nodeXYZs, // in: nodal
    std::vector<double> const&   scalarValues, // in: cell-centered
    std::vector<TriNodes> const& surfaceTriangles, // in: cell-based
    Vec3 const&                  sphereCenter,// in: single XYZ
    double                       sphereRadius, // in: single scalar
    GradientCalcMethod_e         gradientCalcMethod, // in: number of pts to consider for gradient calculation
    std::vector<Vec3>&           gradients, // out: cell-centered dx,dy,dz
    char const*&                 statusMessage); // OUT: contains description of error or "Calculation Successful" if successful

}