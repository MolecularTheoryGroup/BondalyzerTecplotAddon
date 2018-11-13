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

const vector<string> ElementSymbolList = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo"};
const vector<string> ElementNameList = {
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

const vec LogSpace(const double & low, const double & high, const int & n);

double DistSqr(const vec & A, const vec & B);
double Distance(const vec & A, const vec & B);
double VectorAngle(const vec3 & A, const vec3 & B);
const vec3 SphericalToCartesian(const double & r, const double & theta, const double & phi);
const mat44	RotationMatrix(const double & Angle, vec3 Axis);
const vec3 Rotate(const vec3 & Point, const double & Angle, vec3 Axis);
const double TriangleArea(const vec3 & A, const vec3 & B, const vec3 & C);
const double TetVolume(const vector<vec3> & V); 
const double HexahedronInternalPointTetVolume(const vector<vec3> & V);
const double ParallepipedVolume(const vector<vec3> & BV);
const double ParallepipedVolume(const mat33 & BV);

const bool ParallelpidedPointIsInternal(const mat33 & LV, const vec3 & Origin, const vec3 & Pt);

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
	VolExtentIndexWeights_s * VolInfo = NULL;
	Boolean_t IsPeriodic;
	Boolean_t HasGrad;
	Boolean_t HasHess;
	const FieldDataPointer_c * RhoPtr = NULL;
	const vector<FieldDataPointer_c> * GradPtrs = NULL;
	const vector<FieldDataPointer_c> * HessPtrs = NULL;
	const mat33 * BasisVectors = NULL;
	vec3 * Origin = NULL;
	vec3 * EquilPos = NULL;
	int Index = -1;

	double KDisp = 0.0, KGrad = 0.0;
};

struct MultiRootObjects_s{
	gsl_multiroot_function_fdf Func;
	gsl_vector * pos;
	const gsl_multiroot_fdfsolver_type * T;
	gsl_multiroot_fdfsolver * s;
};

enum CompDir_e{
	LessThan = 1,
	GreaterThan
};


#endif
