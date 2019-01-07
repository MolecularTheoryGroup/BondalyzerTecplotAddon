#pragma once
#ifndef CSMDATATYPES_H_
#define CSMDATATYPES_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

#include <vector>

// #include "CSM_VEC.h"
// #include "CSM_MAT.h"
#include "CSM_FIELD_DATA_POINTER.h"

#include <armadillo>
using namespace arma;

#define GTADEBUG 0

#define DEG2RAD 0.017453292519943295769236907684886127134428718885417
#define Ang3PerBohr3 0.1481847435

using std::vector;
using std::string;

#define SIGN(x) (x < 0 ? -1 : (x > 0 ? 1 : 0))

double dot2(const vec & v);

vector<string> const ElementSymbolList = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo"};
vector<string> const ElementNameList = {
"Hydrogen",
"Helium",
"Lithium",
"Beryllium",
"Boron",
"Carbon",
"Nitrogen",
"Oxygen",
"Fluorine",
"Neon",
"Sodium",
"Magnesium",
"Aluminum",
"Silicon",
"Phosphorus",
"Sulfur",
"Chlorine",
"Argon",
"Potassium",
"Calcium",
"Scandium",
"Titanium",
"Vanadium",
"Chromium",
"Manganese",
"Iron",
"Cobalt",
"Nickel",
"Copper",
"Zinc",
"Gallium",
"Germanium",
"Arsenic",
"Selenium",
"Bromine",
"Krypton",
"Rubidium",
"Strontium",
"Yttrium",
"Zirconium",
"Niobium",
"Molybdenum",
"Technetium",
"Ruthenium",
"Rhodium",
"Palladium",
"Silver",
"Cadmium",
"Indium",
"Tin",
"Antimony",
"Tellurium",
"Iodine",
"Xenon",
"Cesium",
"Barium",
"Lanthanum",
"Cerium",
"Praseodymium",
"Neodymium",
"Promethium",
"Samarium",
"Europium",
"Gadolinium",
"Terbium",
"Dysprosium",
"Holmium",
"Erbium",
"Thulium",
"Ytterbium",
"Lutetium",
"Hafnium",
"Tantalum",
"Tungsten",
"Rhenium",
"Osmium",
"Iridium",
"Platinum",
"Gold",
"Mercury",
"Thallium",
"Lead",
"Bismuth",
"Polonium",
"Astatine",
"Radon",
"Francium",
"Radium",
"Actinium",
"Thorium",
"Protactinium",
"Uranium",
"Neptunium",
"Plutonium",
"Americium",
"Curium",
"Berkelium",
"Californium",
"Einsteinium",
"Fermium",
"Mendelevium",
"Nobelium",
"Lawrencium",
"Rutherfordium",
"Dubnium",
"Seaborgium",
"Bohrium",
"Hassium",
"Meitnerium",
"Darmstadtium",
"Roentgenium",
"Copernicium",
"Ununtrium",
"Flerovium",
"Ununpentium",
"Livermorium",
"Ununseptium",
"Ununoctium",
"Atom "
};

vec const LogSpace(double const & low, double const & high, int n);


double DistSqr(vec const & A, vec const & B);
double Distance(vec const & A, vec const & B);
double VectorAngle(vec3 const & A, vec3 const & B);
vec3 const SphericalToCartesian(double const & r, double const & theta, double const & phi);
const mat44	RotationMatrix(double const & Angle, vec3 Axis);
vec3 const Rotate(vec3 const & Point, double const & Angle, vec3 Axis);
double const TriangleArea(vec3 const & A, vec3 const & B, vec3 const & C);
double const TetVolume(vector<vec3> const & V); 
double const HexahedronInternalPointTetVolume(vector<vec3> const & V);
double const ParallepipedVolume(vector<vec3> const & BV);
double const ParallepipedVolume(mat33 const & BV);

bool const ParallelpidedPointIsInternal(mat33 const & LV, vec3 const & Origin, vec3 const & Pt);

typedef double ImportType_t;

enum GPType_e
{
	GPType_Classic = 0,
	GPType_NormalPlaneRhoCP,
	GPType_NormalPlaneEberlyCP,
	GPType_NEB,

	GPType_Invalid = -1
};

/*
*	Struct for parameters for GSL multidimensional root finder
*/
struct VolExtentIndexWeights_s;

struct MultiRootParams_s{
	GPType_e CalcType = GPType_Invalid;
	VolExtentIndexWeights_s * VolInfo = nullptr;
	Boolean_t IsPeriodic;
	Boolean_t HasGrad;
	Boolean_t HasHess;
	FieldDataPointer_c const * RhoPtr = nullptr;
	vector<FieldDataPointer_c> const * GradPtrs = nullptr;
	vector<FieldDataPointer_c> const * HessPtrs = nullptr;
	mat33 const * BasisVectors = nullptr;
	vec3 * Origin = nullptr;
	vec3 * EquilPos = nullptr;
	int Index = -1;

	double KDisp = 0.0, KGrad = 0.0;
};

struct MultiRootObjects_s{
	gsl_multiroot_function_fdf Func;
	gsl_vector * pos;
	gsl_multiroot_fdfsolver_type const * T;
	gsl_multiroot_fdfsolver * s;
};

enum CompDir_e{
	LessThan = 1,
	GreaterThan
};

string GetEdgeString(int ei, int ej);


#endif
