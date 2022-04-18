#pragma once

// #include <gsl/gsl_vector.h>
// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_multiroots.h>

#include "CSM_DATA_SET_INFO.h"
#include "CSM_DATA_TYPES.h"
#include "Edge.h"

#include <vector>
#include <map>
#include <set>

#include <armadillo>
using namespace arma;



#define MaxCPIter 100
#define CheckPosIter 1
#define DefaultCellSpacing 0.1

using std::vector;

// 0.27 bohr is the value used in the BAND cp search, but then I went higher!
static double const SpuriousCPCheckDistance = 0.2;
static double const SpuriousCPDistanceRatioOfSearchGrid = 0.5;

enum CPType_e{
	CPType_Nuclear = -3,
	CPType_Bond = -1,
	CPType_Ring = 1,
	CPType_Cage = 3,
	CPType_RingFF = 11,
	CPType_CageFF = 13,

	CPType_Invalid = -99
};
enum CPTypeNum_e {
	CPTypeNum_Invalid = -1,
	CPTypeNum_Nuclear,
	CPTypeNum_Bond,
	CPTypeNum_Ring,
	CPTypeNum_Cage,
	CPTypeNum_RingFF,
	CPTypeNum_CageFF
};

// static char const CPTypeList[] = { -3, -1, 1, 3, 11, 13 };
static vector<CPType_e> const CPTypeList = { CPType_Nuclear, CPType_Bond, CPType_Ring, CPType_Cage, CPType_RingFF, CPType_CageFF };

static vector<string> const CPNameList = { "Nuclear", "Bond", "Ring", "Cage", "Ring FF", "Cage FF" };

static vector<int> const CPPrincDirInds = {
	0, // atoms, most negative direction
	2, // bonds, most (only) positive direction
	0, // rings, most (only) negative direction
	2, // cages, most positive direction
	-1, // ringFF n/a
	-1 // cageFF n/a
};

static vector<CPTypeNum_e> const CPSaddleTypeNums = { CPTypeNum_Bond, CPTypeNum_Ring };
static vector<CPTypeNum_e> const CPNearFieldTypes = { CPTypeNum_Nuclear, CPTypeNum_Bond, CPTypeNum_Ring, CPTypeNum_Cage };

static vector<ColorIndex_t> const CPColorList = { White_C, Red_C, Green_C, Cyan_C, Custom5_C, Custom6_C };

/*
	*	Group of critical points
	*/
class CritPoints_c
{
public:
	CritPoints_c();
	// Specify a cutoff value during construction
	CritPoints_c(double const & RhoCutoff, int NumDimensions);
	// Construct from a set of other CritPoints_c's
	CritPoints_c(vector<CritPoints_c> const & CPLists);
	// Construct from existing CP zone
	CritPoints_c(int CPZoneNum, 
		vector<int> const & XYZVarNums,
		int CPTypeVarNum,
		int RhoVarNum = -1, 
		MultiRootParams_s *MR = nullptr);
	~CritPoints_c();

	/*
		*	Operator overloads
		*/
	CritPoints_c & operator+=(CritPoints_c const & rhs);

	/*
		*	Getter methods
		*/
	double GetRhoCutoff() const { return m_RhoCutoff; }

	int NumCPs() const { return m_TotNumCPs; }
	int NumCPs(int TypeNum) const { return m_NumCPs[TypeNum]; }
	int NumAtoms() const { return m_NumCPs[CPTypeNum_Nuclear]; }
	int NumBonds() const { return m_NumCPs[CPTypeNum_Bond]; }
	int NumRings() const { return m_NumCPs[CPTypeNum_Ring]; }
	int NumCages() const { return m_NumCPs[CPTypeNum_Cage]; }
	int NumFFRings() const { return m_NumCPs[CPTypeNum_RingFF]; }
	int NumFFCages() const { return m_NumCPs[CPTypeNum_CageFF]; }
	int NumDimensions() const { return m_Dimensions; }

	vector<int> GetTypeNumOffsetFromTotOffset(int TotOffset) const;
	CPType_e GetTypeFromTotOffset(int TotOffset) const;
	int GetTotOffsetFromTypeNumOffset(int TypeNum, int TypeOffset) const;

	double GetMinCPDist(vector<CPType_e> const & CPTypes = CPTypeList);
	double GetMinCPDist(int CPTypeInd, int CPOffset, vector<CPType_e> const & CPTypes = CPTypeList);
	double GetMinCPDist(int CPTotOffset, vector<CPType_e> const & CPTypes = CPTypeList);

	double GetRho(int TypeNum, int Offset) const { return m_Rho[TypeNum][Offset]; }
	double GetRho(int TotOffset) const;
	vec3 GetXYZ(int TypeNum, int Offset) const { return m_XYZ[TypeNum][Offset]; }
	vec3 GetXYZ(int TotOffset) const;
	vec3 GetPrincDir(int TypeNum, int Offset) const { return m_PrincDir[TypeNum][Offset]; }
	vec3 GetPrincDir(int TotOffset) const;

	vec3 GetEigVals(int TypeNum, int Offset) const { return m_EigVals[TypeNum][Offset]; }
	vec3 GetEigVals(int TotOffset) const;

	mat33 GetEigVecs(int TypeNum, int Offset) const { return m_EigVecs[TypeNum][Offset]; }
	mat33 GetEigVecs(int TotOffset) const;

	bool HasEdge(Edge const & e, int * zoneNum = nullptr) const;
	bool HasEdge(int e1, int e2, int * zoneNum = nullptr) const { return HasEdge(MakeEdge(e1, e2), zoneNum); }
	bool HasEdgeNoPathRecorded(Edge const & e, int searchDepthLimit) const;
	bool HasEdge(Edge const & e, int searchDepthLimit, vector<int> * Path = nullptr) const;
	bool HasEdge(int e1, int e2, int searchDepth, vector<int> * Path = nullptr) const { return HasEdge(MakeEdge(e1, e2), searchDepth, Path); }

	vec3 ClosestPoint(vec3 const & Pt, int & TotCPOffset, double & MinDist) const;

	Boolean_t IsValid() const;

	/*
		*	Setter methods
		*/
	void SetMinCPDist(double const & MinCPDist){ m_MinCPDist = MinCPDist; }
	Boolean_t AddPoint(double const & Rho,
		vec3 const & Pos,
		vec3 const & PrincDir,
		char Type);
	void RemPoint(int TypeNum, int PointIndex);
	void Append(CritPoints_c const & rhs);

	/*
		*	Mutators and other methods
		*/

	Boolean_t FindMinCPDist(vector<CPType_e> const & CPTypes);
	void RemoveSpuriousCPs(double const & CheckDist = SpuriousCPCheckDistance);
	void GenerateCPGraph(vector<string> const & AuxDataSubTypeList = CSMAuxData.CC.OneSkeletonGPSubTypes);
	
	vector<int> SaveAsOrderedZone(vector<int> const & XYZVarNum, int RhoVarNum = -1, Boolean_t SaveCPTypeZones = FALSE, int VolZoneNum = -1);

private:


	/*
		*	m_Rho, m_XYZ, m_PrincDir, and m_NumCPs are length 6 so they store
		*	information for the 6 types of critical points.
		*/
	vector<double> m_Rho[6];
	vector<vec3> m_XYZ[6], m_PrincDir[6], m_EigVals[6];
	vector<mat33> m_EigVecs[6];
	int m_Dimensions;
	int m_TotNumCPs, m_NumCPs[6];
	double m_MinCPDist;
	Boolean_t m_MinCPDistFound;
	double m_RhoCutoff;
	int m_ZoneNum = -1;

	std::map<Edge,int> m_CPEdgeToZoneNumMap;
	std::map<int,std::set<int> > m_CPAdjacencyList;

	vector<CPType_e> m_MinDistCPTypes;
};

void SetCPZone(int ZoneNum);

Boolean_t FindCPs(CritPoints_c & CPs,
	VolExtentIndexWeights_s VolInfo,
	Boolean_t IsPeriodic,
	vector<int> const & StartIJK,
	vector<int> const & EndIJK,
	FieldDataPointer_c const & RhoPtr,
	vector<FieldDataPointer_c> const & GradPtrs,
	vector<FieldDataPointer_c> const & HessPtrs);

Boolean_t FindCPs(CritPoints_c & CPs,
	VolExtentIndexWeights_s const & VolInfo,
	double const & CellSpacing,
	double & RhoCutoff,
	Boolean_t IsPeriodic,
	FieldDataPointer_c & RhoPtr,
	vector<FieldDataPointer_c> & GradXYZPtrs,
	vector<FieldDataPointer_c> & HessPtrs,
	double AgitationFactor = -1.0,
	int AgitationMaxNumIter = 1000);

Boolean_t CritPointInCell(vector<int> const & IJK,
	vec3 & Point,
	vec3 & PrincDir,
	double & RhoValue,
	double const & RhoCutoff,
	char & Type,
	MultiRootParams_s & RootParams,
	MultiRootObjects_s & MR);

Boolean_t CritPointInCell(
	vec3 const & CellMinXYZ,
	vec3 const & CellMaxXYZ,
	vec3 & Point,
	vec3 & PrincDir,
	double & RhoValue,
	double const & RhoCutoff,
	char & Type,
	MultiRootParams_s & RootParams,
	MultiRootObjects_s & MR);
