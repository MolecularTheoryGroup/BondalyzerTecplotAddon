#pragma once

#include <vector>

#include "SimpleTypes.h"
#include "Vec3.h"
#include "Edge.h"
#include "TriNodes.h"

namespace tpcsm {

// returns true if successful, false otherwise
bool updateSphericalTriangulation(
    std::vector<Vec3> const&     initialNodeXYZs, // IN: location of the initial triangulation's nodes
    std::vector<TriNodes> const& initialTriangleNodes, // IN: connectivity of each triangle, defined in terms of points above
    Vec3 const&                  sphereCenter,
    double                       radius,
    std::vector<Vec3> const&     constraintNodes, // IN:
    std::vector<Edge> const&     constraintSegments, // IN: defined in terms of constraintNodes above
    bool                         includeConstraintSegmentsOnly, // IN: false means include intersections 
    std::vector<Vec3>&           nodeXYZs, // OUT: location of the resulting triangulation's nodes (includes all of those in initialNodeXYZs)
    std::vector<Index_t>&        nodeEdge, // OUT: for each node in nodeXYZ which constraintSegments it lines on or -1 if none
    std::vector<TriNodes>&       triangleNodes, // OUT: connectivity of each triangle, defined in terms of points above
    char const*&                 statusMessage); // OUT: contains description of error or "Triangulation Successful" if successful

}