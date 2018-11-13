#pragma once
#ifndef CSMGRADPATH_H_
#define CSMGRADPATH_H_

#include <vector>

#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"

#include <armadillo>
using namespace arma;



using std::vector;

#define GP_NumPointsBufferFactor	0.2
#define GP_StallPointCount			50
#define GP_StallNumPointsToCheck	50
#define GP_StallPointDistTol		1e-3
#define GP_MaxNumPoints				10000
#define GP_PlaneCPStallCount		30
#define GP_PlaneCPMaxIter			100

class CritPoints_c;

enum GPTerminate_e
{
	GPTerminate_AtBoundary = 0,
	GPTerminate_AtPoint,
	GPTerminate_AtPointRadius,
	GPTerminate_AtCP,
	GPTerminate_AtCPRadius,
	GPTerminate_AtRhoValue,

	GPTerminate_Invalid = -1
};

struct GradPathParams_s{
	/*
	*	Raw pointers to Tec360 field data.
	*	These behave like C-style arrays.
	*/

	vector<FieldDataPointer_c> HessPtrs;
	Boolean_t HasHess;
	vector<FieldDataPointer_c> GradPtrs;
	Boolean_t HasGrad;
	FieldDataPointer_c RhoPtr;

	VolExtentIndexWeights_s VolZoneInfo;

	StreamDir_e Direction;

	GradPathParams_s & operator=(const GradPathParams_s & rhs);
	const Boolean_t operator==(const GradPathParams_s & rhs) const;
};





class GradPathBase_c
{
public:
	GradPathBase_c();
	~GradPathBase_c();

	/*
	*	Constructor for loading an existing 1-D ordered
	*	zone as a GradPath_c
	*/
	GradPathBase_c(EntIndex_t ZoneNum,
		const vector<EntIndex_t> & XYZRhoVarNums,
		const AddOn_pa & AddOnID);


	/*
	*	Operator declarations
	*/
	GradPathBase_c & operator=(const GradPathBase_c & rhs);
	const Boolean_t IsSame(const GradPathBase_c & rhs) const;
	const Boolean_t operator==(const GradPathBase_c & rhs) const;
	GradPathBase_c & operator+=(const GradPathBase_c & rhs);
	const GradPathBase_c operator+(const GradPathBase_c & rhs) const;
	const vec3 operator[](const int & i) const;


	const Boolean_t SetStartEndCPNum(const int & CPNum, const int & StartEnd) {
		if (CPNum >= 0){
			m_StartEndCPNum[StartEnd] = CPNum;
			return TRUE;
		}
		else return FALSE;
	}
	const Boolean_t SetStartEndCPNum(int * CPNums){
		for (int i = 0; i < 2; ++i){
			if (CPNums[i] >= 0)
				m_StartEndCPNum[i] = CPNums[i];
			else
				return FALSE;
		}
		return TRUE;
	}

	const Boolean_t Resample(const int & NumPoints);
	const Boolean_t Reverse();
	GradPathBase_c & ConcatenateResample(GradPathBase_c & rhs, const int & NumPoints);
	GradPathBase_c & ConcatenateResample(GradPathBase_c & rhs, const int & NumPoints, int & BrigePtNum);
	GradPathBase_c & Concatenate(const GradPathBase_c & rhs) { return *this += rhs; }
	const Boolean_t Trim(const vec3 & Point, const double & Radius);
	void PointAppend(const vec3 & Point, const double & Rho);
	void PointPrepend(const vec3 & Point, const double & Rho);

	const LgIndex_t GetZoneNum() const { return m_ZoneNum; }
	const EntIndex_t SaveAsOrderedZone(const string & ZoneName = "Gradient Path", const ColorIndex_t MeshColor = Black_C);
	const EntIndex_t SaveAsOrderedZone(const string & ZoneName,
		vector<FieldDataType_e> & VarDataTypes,
		const vector<int> & XYZVarNums,
		const int & RhoVarNum,
		const Boolean_t DoActivate = FALSE,
		const ColorIndex_t MeshColor = Black_C);
	const Boolean_t SaveAsCSV(const string & PathToFile, const Boolean_t & IncludeVars = FALSE);

	/*
	*	Getter methods
	*/
	const Boolean_t IsMade() const { return m_GradPathMade; }
	const Boolean_t IsReady() const { return m_GradPathReady; }
	const double GetLength();
	const int GetCount() const {
		return static_cast<const int>(m_XYZList.size());
	}
	vector<int> GetStartEndCPNum() const{
		vector<int> Nums = { m_StartEndCPNum[0], m_StartEndCPNum[1] };
		return Nums;
	}
	const int GetStartEndCPNum(const unsigned int & i) const{
		REQUIRE(i < 2);
		return m_StartEndCPNum[i];
	}

	const vec3 XYZAt(const int & i) const { return operator[](i); }
	const double RhoAt(const int & i) const;

	const GradPathBase_c SubGP(int BegPt, int EndPt) const;

	const vec3 ClosestPoint(const vec3 & rhs) const;
	const vec3 ClosestPoint(const vec3 & rhs, int & PtNum) const;

private:

	const unsigned int GetInd(const int & i) const{
		REQUIRE(abs(i) < m_XYZList.size());
		return (i >= 0) ? i : m_XYZList.size() + i;
	}

protected:
	/*
	*	Containers for GradPath values.
	*/
	vector<vec3> m_XYZList;
	vector<double> m_RhoList;

	// Overall 0-based index of start-end CPs
	int m_StartEndCPNum[2];

	int m_NumGPPoints;

	double m_Length = -1;

	Boolean_t m_GradPathReady = FALSE;
	Boolean_t m_GradPathMade = FALSE;
	LgIndex_t m_ZoneNum = 0;
};


class GradPath_c : public GradPathBase_c
{
	friend class FESurface_c;
public:
	/*
		*	Constructors and destructors
		*/

	/*
		*	Default constructor
		*/
	GradPath_c();

	/*
		*	Constructor for grad path that will make itself
		*/
	GradPath_c(const vec3 & StartPoint,
		const StreamDir_e & Direction,
		const int & NumGPPoints,
		const GPTerminate_e & HowTerminate,
		vec3 * TermPoint,
		const vector<FieldDataPointer_c> & CPXYZPtrs,
		int * NumCPs,
		double * TermPointRadius,
		double * TermValue,
		const vector<int> & MaxIJK,
		const vec3 & MaxXYZ,
		const vec3 & MinXYZ,
		const vector<FieldDataPointer_c> & GradPtrs,
		const FieldDataPointer_c & RhoPtr);

	GradPath_c(const vec3 & StartPoint,
		const StreamDir_e & Direction,
		const int & NumGPPoints,
		const GPType_e & GPType,
		const GPTerminate_e & HowTerminate,
		vec3 * TermPoint,
		CritPoints_c * CPs,
		double * TermPointRadius,
		double * TermValue,
		VolExtentIndexWeights_s & VolInfo,
		const vector<FieldDataPointer_c> & HessPtrs,
		const vector<FieldDataPointer_c> & GradPtrs,
		const FieldDataPointer_c & RhoPtr);

	GradPath_c::GradPath_c(EntIndex_t ZoneNum,
		const vector<EntIndex_t> & XYZRhoVarNums,
		const AddOn_pa & AddOnID);

	/*
		* Copy constructor
		*/
	GradPath_c(const GradPath_c & rhs);
	GradPath_c(const GradPathBase_c & rhs);
	~GradPath_c();

	/*
	*	Operator declarations
	*/
	GradPath_c & operator=(const GradPath_c & rhs);
	GradPath_c & operator=(const GradPathBase_c & rhs);
	const Boolean_t operator==(const GradPath_c & rhs) const;


	/*
		*	Getter methods
		*/

	/*
		*	Setter methods
		*/
	const Boolean_t SetupGradPath(const vec3 & StartPoint,
		const StreamDir_e & Direction,
		const int & NumGPPoints,
		const GPTerminate_e & HowTerminate,
		vec3 * TermPoint,
		const vector<FieldDataPointer_c> & CPXYZPtrs,
		int * NumCPs,
		double * TermPointRadius,
		double * TermValue,
		const vector<int> & MaxIJK,
		const vec3 & MaxXYZ,
		const vec3 & MinXYZ,
		const vector<FieldDataPointer_c> & GradPtrs,
		const FieldDataPointer_c & RhoPtr);

	const Boolean_t SetupGradPath(const vec3 & StartPoint,
		const StreamDir_e & Direction,
		const int & NumGPPoints,
		const GPType_e & GPType,
		const GPTerminate_e & HowTerminate,
		vec3 * TermPoint,
		CritPoints_c * CPs,
		double * TermPointRadius,
		double * TermValue,
		VolExtentIndexWeights_s & VolInfo,
		const vector<FieldDataPointer_c> & HessPtrs,
		const vector<FieldDataPointer_c> & GradPtrs,
		const FieldDataPointer_c & RhoPtr);

	const Boolean_t Seed(const bool DoResample = true);
	const Boolean_t SetMixingFactor(const double & MixFactor){
		if (MixFactor >= 0.0 && MixFactor <= 1.0)
			m_DirMixFactor = MixFactor;
		else
			return FALSE;

		return TRUE;
	}

private:
	//const Boolean_t SetIndexAndWeightsForPoint(vec3 & Point);
	const double RhoByCurrentIndexAndWeights();
	const Boolean_t SeedInDirection(const StreamDir_e & Direction);

	/*
		* Special gradient path algorithm functions
		*/
		
	GPType_e m_GPType = GPType_Classic;
	Boolean_t m_SGPMade = FALSE;
	char m_RidgeRank;
	double m_DirMixFactor;

	size_t m_ODE_NumDims = 3;


	/*
		*	Data that the GSL ODE solver must access.
		*	Kept in a struct because the grad path itself
		*	doesn't need them except to pass them to the
		*	function provided to the ODE solver.
		*/
	GradPathParams_s m_ODE_Data;


	// If terminating at:
	// Point
	vec3 m_TermPoint;
	// Any CP using tecplot zones
	vector<FieldDataPointer_c> m_CPXYZPtrs;
	int m_NumCPs;
	// ...using CritPoints_c
	CritPoints_c * m_CPs;
	// For point or CP:
	double m_TermPointRadiusSqr;
	// Variable value
	// Which item in dep var list used for terminal value?
	double m_TermValue;

	GPTerminate_e m_HowTerminate;
	vec3 m_StartPoint;
};



const Boolean_t GPsStraddleIB(const GradPath_c & GP1,
	const GradPath_c & GP2,
	const double & IBCheckAngle,
	const double & IBCheckDistRatio);

const Boolean_t CPInNormalPlane(vec3 & StartPt, const vec3 & PlaneBasis, MultiRootObjects_s & MR);


class NEBGradPath_c : public GradPathBase_c
{
public:
	NEBGradPath_c();
	~NEBGradPath_c();

	NEBGradPath_c(const vec3 & StartPt, 
		const vec3 & EndPt, 
		const unsigned int & NumPts);

	const Boolean_t Relax(const double & StepRatio,
		const double & Tol, 
		const unsigned int MaxIter,
		MultiRootParams_s & Params);

private:

};

/*
*	Redo of path stitching algorithm, formally TryPolyLines
*/
void StitchPaths(
	const vector<int> &     L,       // indices of points in P
	const vector<int> &     R,
	const vector<vec3> &     P,
	vector<vector<int> > &     T       // triplets of integers specifying nodes of triangles
	);
void StitchCapPaths(
	const vector<int> &     L,       // indices of points in P
	const vector<int> &     R,
	const vector<int> &		C,       // indices of points in the cap, C
	const vector<vec3> &     P,
	vector<vector<int> > &     T       // triplets of integers specifying nodes of triangles
	);


#endif