
#ifndef CSMDATALOADERFUNCTIONS_H_
#define CSMDATALOADERFUNCTIONS_H_

#include <string>
#include <vector>

#include "CSM_DATA_TYPES.h"

#include <armadillo>
using namespace arma;



using std::string;
using std::vector;

/*
 *	Used to know whether data was already loaded or not
 */
enum LoaderStatus_e{
	LoaderFailed = 0,
	Loading,
	Loaded
};

/*
* I use my own data types for easy
* changing of precision
*/
typedef int lg_int_t;
typedef double lg_double_t;

//	Enumerator for variable types and accompanying strings
enum VarType_e{
	SCF = 0,
	FRAG,
	ORTHO,
	TRANS,
	GEOMETRY,
	NOTYPE
};
//	Enumerator for A/B (spin polarization) types and accompanying strings
enum ab_e{
	A = 0,
	B,
	NOAB
};

/*
* Datatype structure: Contains all the necessary information and methods for
* a given variable for it's identification and the calculation of any extra
* variables for which it's used
*/
struct T41Var_s{
	string T41NameStr;
	string T41FullNameStr;
	string NameStr;
	string BaseNameStr;
	VarType_e VarType;
	ab_e DTAB;
	std::streamoff ByteOffset;
	double FirstValue;
	Boolean_t UseForIsosurface;

	T41Var_s();
	T41Var_s(const char* iName);
	T41Var_s(const char* iT41Name, const char* iName, const char* iBaseName, Boolean_t UseForIsoSurfaceDisplay);
	~T41Var_s();

	void Set(const char* iT41Name, const char* iName, const char* iBaseName, Boolean_t UseForIsoSurfaceDisplay);
	void SetAB(const string &TempVarName);
	void SetAB(ab_e In);
	void SetType(const string &TypeStr);
	void SetType(VarType_e In);
	const Boolean_t operator== (const T41Var_s & T2) const;
	const Boolean_t operator!= (const T41Var_s & T2) const;
};
/*
*	Atom group structure for storing the atomic positions and charges of
*	all atoms. Atoms are stored in the T41 in groups, so the process of
*	storing their information is somewhat simplified since you know you'll
*	find all the Carbons etc. contiguously.
*/

struct AtomColor_s{
	vector<string> Names;
	ColorIndex_t Color;

	AtomColor_s(vector<string> TheNames, ColorIndex_t TheColor){
		Names = TheNames;
		Color = TheColor;
	}
};


struct AtomGroup_s{
	string Name;
	ColorIndex_t AtomColor;
	int Count;
	vector<double> Positions[3];
	vector<double> Charges;

	AtomGroup_s(){}
	AtomGroup_s(const string & AtomTypeName);
	void AddPosition(vec3 & Position);
};

const string GetAtomStr(const int & AtomNum);

//	Checks if a line in the input text file is the start of a new block
int IsVarLine(const string &CheckLine);
//	Searching the datatype list for a datatype
Boolean_t IsInDataTypeList(const vector<T41Var_s> &SearchVec, const T41Var_s &ToFind);
//	Searching a list of strings for a string
Boolean_t VectorContainsString(const vector<string> &SearchVec, const string &ToFind);
//	Populate the list that associates atoms with their respective colors.
//	Uses ADF as the basis for atom colors.
void PopulateAtomColorList(vector<AtomColor_s> &AtomColorList);
//	Provides the correct color based on an atom's name.
const ColorIndex_t GetAtomColor(const vector<AtomColor_s> &AtomColorList, const string &AtomName);
//	Creates the list of types outputable by densf and inputs their properties
int PopulateTypeList(vector<T41Var_s> &TypeList);

const Boolean_t CreateAtomZonesFromAtomGroupList(const vector<AtomGroup_s> & AtomGroupList,
	const vector<string> & XYZVarNames,
	vector<FieldDataType_e> & VarDataTypes,
	const int & GroupIndex);

#endif