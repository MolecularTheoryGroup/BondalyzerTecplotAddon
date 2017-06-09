#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

const double pi    = 3.1415926535897932385; // M_PI; // 4.0*atan(1.0);

// Common square roots encountered in the code:
const double sq2   = 1.4142135623730950488;
const double sq3   = 1.7320508075688772935;
const double sq5   = 2.2360679774997896964;
const double sq6   = 2.4494897427831780982;
const double sq7   = 2.6457513110645905905;
const double sq10  = 3.1622776601683793320;
const double sq11  = 3.3166247903553998491;
const double sq13  = 3.6055512754639892931;
const double sq15  = 3.8729833462074168852;
const double sq19  = 4.3588989435406735522;
const double sq23  = 4.7958315233127195416;
const double sq30  = 5.4772255750516611346;
const double sq35  = 5.9160797830996160426;
const double sq71  = 8.4261497731763586306;
const double sq105 = 1.0246950765959598383e+01;

// Common recipricols of the above square roots enountered in the code:
const double osq2  = 7.0710678118654752440e-01;
const double osq3  = 5.7735026918962576449e-01;
const double osq5  = 4.4721359549995793928e-01;
const double osq7  = 3.7796447300922722721e-01;
const double osq13 = 2.7735009811261456101e-01;
const double osq19 = 2.2941573387056176591e-01;

// Common fractions where we would prefer multiplication over division:
const double onehalf    = 0.5;
const double onethird   = 3.3333333333333333333e-01;
const double oneseventh = 1.4285714285714285714e-01;
const double oneninth   = 1.1111111111111111111e-01;

// Used to test for machine zero:
const double EPSILON = 0.5e-15;

// Faster lookup for small factorials.  See also: dog_math.h
const double  factorial[11] = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0 };
const int    ifactorial[11] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 };

#endif
