// heffNorm.h
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include "TASSERT.h"

namespace tpcsm {
    /*
    * Calculate the "norm" from the Hessian Eigensystem to the Frenet-Frame. This norm is the minimum
    * rotation angle (in radians) required to rotate the Frenet into alignment with the Eigensystem, including
    * all variations of that system using vectors in the opposite direction.
    * The eigenvectors can be specified in any order, and the same for the Frenet-Frame vectors
    * (The parenthetical comments below are just suggestions).
    */
    double calculateHeffNorm(
        double he1x, double he1y, double he1z, // Hessian eigenvector 1
        double he2x, double he2y, double he2z, // Hessian eigenvector 2
        double he3x, double he3y, double he3z, // Hessian eigenvector 3
        double ff1x, double ff1y, double ff1z, // Frenet-Frame vector 1 (tangent)
        double ff2x, double ff2y, double ff2z, // Frenet-Frame vector 2 (normal)
        double ff3x, double ff3y, double ff3z); // Frenet-Frame vector 3 (binormal)
}
