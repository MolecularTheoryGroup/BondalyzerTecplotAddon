#pragma once

#include <vector>
#include <set>

#include "CSM_DATA_TYPES.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"

#include <armadillo>
using namespace arma;



using std::vector;

#define GP_NumPointsBufferFactor	10
#define GP_StallPointCount			20
#define GP_StallNumPointsToCheck	10
#define GP_StallPointDistTol		1e-6
#define GP_MaxNumPoints				10000
#define GP_PlaneCPStallCount		30
#define GP_PlaneCPMaxIter			100
#define GP_MaxStepSize				0.01

static double const GP_DeviationAngleMaxCutoff = 170. / 180. * PI;


static double const DefaultGPDeviationCheckAngleFactor = 0.15;
static double const DefaultGPDeviationCheckDistFactor = 0.05;
static double const DefaultGPDeviationStartPtLengthFactor = 0.1;


static double const GB_DefaultKinkCheckAngle = 0.3 * PI;

class CritPoints_c;
class FESurface_c;

enum GPResampleMethod_e
{
	GPResampleMethod_Linear = 0,
	GPResampleMethod_Adaptive,

	GPResampleMethod_Invalid = -1
};

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

	GradPathParams_s & operator=(GradPathParams_s const & rhs);
	Boolean_t const operator==(GradPathParams_s const & rhs) const;
};





class GradPathBase_c
{
	friend class FESurface_c;
public:
	GradPathBase_c();
	~GradPathBase_c();

	/*
	*	Constructor for loading an existing 1-D ordered
	*	zone as a GradPath_c
	*/
	GradPathBase_c(EntIndex_t ZoneNum,
		vector<EntIndex_t> const & XYZRhoVarNums,
		AddOn_pa const & AddOnID);


	/*
	*	Operator declarations
	*/
	GradPathBase_c & operator=(GradPathBase_c const & rhs);
	Boolean_t IsSame(GradPathBase_c const & rhs) const;
	Boolean_t operator==(GradPathBase_c const & rhs) const;
	GradPathBase_c & operator+=(GradPathBase_c const & rhs);
	GradPathBase_c operator+(GradPathBase_c const & rhs) const;
	vec3 operator[](int i) const;


	Boolean_t SetStartEndCPNum(int CPNum, int StartEnd) {
		if (CPNum >= 0){
			m_StartEndCPNum[StartEnd] = CPNum;
			return TRUE;
		}
		else return FALSE;
	}
	Boolean_t SetStartEndCPNum(int * CPNums){
		for (int i = 0; i < 2; ++i){
			if (CPNums[i] >= 0)
				m_StartEndCPNum[i] = CPNums[i];
			else
				return FALSE;
		}
		return TRUE;
	}

	void Clear();


	void RemoveKinks(double const & AngleCutoff = GB_DefaultKinkCheckAngle);


	vector<int> GradPathBase_c::GetCPCoincidentPoints(CritPoints_c const * CPs, std::set<CPType_e> const & CPTypes = {}, double tol = 1e-8) const;

	Boolean_t Resample(int NumPoints, vector<int> & ProtectedPoints, GPResampleMethod_e Method = GPResampleMethod_Adaptive);
	Boolean_t Resample(int NumPoints, GPResampleMethod_e Method = GPResampleMethod_Adaptive) { vector<int> v; return Resample(NumPoints, v, Method); }
	Boolean_t ResampleByRhoVals(vector<double> const & RhoVals, vector<int> const & ProtectedPoints);
	Boolean_t Reverse();
	GradPathBase_c & ConcatenateResample(GradPathBase_c & rhs, int NumPoints, GPResampleMethod_e Method = GPResampleMethod_Adaptive);
	GradPathBase_c & ConcatenateResample(GradPathBase_c & rhs, int NumPoints, int & BrigePtNum, GPResampleMethod_e Method = GPResampleMethod_Adaptive);
	GradPathBase_c & Concatenate(GradPathBase_c const & rhs) { return *this += rhs; }
	Boolean_t Trim(vec3 const & Point, double const & Radius);
	void PointAppend(vec3 const & Point, double const & Rho);
	void PointPrepend(vec3 const & Point, double const & Rho);
	bool ReplaceRhoList(vector<double> const & RhoList) { if (m_RhoList.size() == RhoList.size()) { m_RhoList = RhoList; return true; } else return false; }

	LgIndex_t GetZoneNum() const { return m_ZoneNum; }
	EntIndex_t SaveAsOrderedZone(string const & ZoneName = "Gradient Path", ColorIndex_t const MeshColor = Black_C);
	EntIndex_t SaveAsOrderedZone(string const & ZoneName,
		vector<FieldDataType_e> & VarDataTypes,
		vector<int> const & XYZVarNums,
		int RhoVarNum,
		Boolean_t const DoActivate = FALSE,
		ColorIndex_t const MeshColor = Black_C);
	Boolean_t SaveAsCSV(string const & PathToFile, Boolean_t IncludeVars = FALSE);

	Boolean_t TruncateAtRhoValue(double const & RhoVal);

	/*
	*	Getter methods
	*/
	Boolean_t IsMade() const { return m_GradPathMade; }
	Boolean_t IsReady() const { return m_GradPathReady; }
	double GetLength();
	double GetLength(int AtInd = -1) const;
	int GetIndAtLength(double const & Length) const;
	int GetCount() const {
		return static_cast<const int>(m_XYZList.size());
	}
	vector<int> GetStartEndCPNum() const{
		vector<int> Nums = { m_StartEndCPNum[0], m_StartEndCPNum[1] };
		return Nums;
	}
	int GetStartEndCPNum(unsigned int i) const{
		REQUIRE(i < 2);
		return m_StartEndCPNum[i];
	}

	vec3 ClosestMaxCurvaturePoint(vec3 const & CheckPt) const { int i; return ClosestMaxCurvaturePoint(CheckPt, i); }
	vec3 ClosestMaxCurvaturePoint(vec3 const & CheckPt, int & PtInd) const; 

	bool GetPointAtRhoValue(double const & RhoValue, vec3 & OutVec, int & Ind, double & Weight) const;

	bool GetDeviationMidpointRhoBasedFromOtherGP(GradPathBase_c const & GP, double const & CheckAngle, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2) const;
	bool GetDeviationMidpointClosestPointBasedFromOtherGP(GradPathBase_c const & GP, double const & CheckAngle, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2) const;
	bool GetDeviationMidpointAsTerminalAngleAndDistanceFactorFromOtherGP(GradPathBase_c const & GP, 
		double const & CheckTermAngleFactor, 
		double const & CheckTermDistFactor, 
		double const & StartPtAsFactorOfLength, 
		vec3 & OutPt, int & Ind1, 
		int & Ind2, vec3 & Pt2,
		double const & MinCheckDist = 0.0,
		double const & MinCheckAngle = PIOVER2) const;
	bool GetMaxSeparationMidpointFromOtherGP(GradPathBase_c const & GP, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2, double & MaxDist, int step=1) const;
	bool GetMaxSeparationMidpointFromOtherGPRhoBased(GradPathBase_c const & GP, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2, double & MaxDist, int step=1) const;
	bool GetSeparationMidpointAtDistFromOtherGP(GradPathBase_c const & GP, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2, const double & MaxDist, const double & StartPtAsFactorOfLength) const;
	double GetMaxSeparationFromNeighboringGPs(vector<GradPathBase_c const *> GPs, vector<int> & IndList, int & Ind) const;

	vector<vec3> const * GetXYZListPtr() const { return &m_XYZList; }
	vector<double> const * GetRhoListPtr() const { return &m_RhoList; }

	vec3 XYZAt(int i) const { return operator[](i); }
	double RhoAt(int i) const;

	GradPathBase_c SubGP(int BegPt, int EndPt) const;

	vec3 ClosestPoint(vec3 const & rhs) const;
	vec3 ClosestPoint(vec3 const & rhs, int & PtNum) const;
	int GetSphereIntersectionPoint(vec3 const & Center, double const & Radius, vec3 & IntersectionPoint) const;



private:

	void FixDeviationMidpoint(GradPathBase_c const & GP, vec3 & MidPt, int Ind1, int Ind2) const;
	Boolean_t ResampleAdaptive(int NumPoints, vector<int> & ProtectedPoints);

	unsigned int GetInd(int i) const{
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
	GradPath_c(vec3 const & StartPoint,
		StreamDir_e const & Direction,
		int NumGPPoints,
		GPTerminate_e const & HowTerminate,
		vec3 * TermPoint,
		vector<FieldDataPointer_c> const & CPXYZPtrs,
		int * NumCPs,
		double * TermPointRadius,
		double * TermValue,
		vector<int> const & MaxIJK,
		vec3 const & MaxXYZ,
		vec3 const & MinXYZ,
		vector<FieldDataPointer_c> const & GradPtrs,
		FieldDataPointer_c const & RhoPtr);

	GradPath_c(vec3 const & StartPoint,
		StreamDir_e const & Direction,
		int NumGPPoints,
		GPType_e const & GPType,
		GPTerminate_e const & HowTerminate,
		vec3 * TermPoint,
		CritPoints_c * CPs,
		double * TermPointRadius,
		double * TermValue,
		VolExtentIndexWeights_s & VolInfo,
		vector<FieldDataPointer_c> const & HessPtrs,
		vector<FieldDataPointer_c> const & GradPtrs,
		FieldDataPointer_c const & RhoPtr,
		FESurface_c const * Surf = nullptr);

	/* 
	 * Constructor for grad path from existing i-ordered zone
	 */
	GradPath_c(EntIndex_t ZoneNum,
		vector<EntIndex_t> const & XYZRhoVarNums,
		AddOn_pa const & AddOnID);


	// Constructor that makes a new GradPath using the midpoints
	// of the provided GradPaths.
	GradPath_c(vector<GradPathBase_c const *> const & GPs,
		double const & Path1Weight,
		double const & RhoVal,
		int GPStep);

	// Constructor that makes a surface from the provided 
	// GradPaths and then seeds a new grad path constrained to that
	// surface.
	GradPath_c(vector<GradPath_c const *> const & GPs,
		vec3 const & SeedPoint,
		StreamDir_e GPDir,
		vector<std::pair<int, int> > StartEndGPInds = {},
		vec3 const * TermPoint = nullptr,
		double const * TermPointRadiu = nullptr,
		GPTerminate_e GPTermType = GPTerminate_Invalid);

	// Construct from vector<vec3> points and optional vector<double> rho values
	GradPath_c(vector<vec3> const & XYZList, vector<double> const RhoList = vector<double>());

	/*
	* Copy constructor
	*/
	GradPath_c(GradPath_c const & rhs);
	GradPath_c(GradPathBase_c const & rhs);

	~GradPath_c();

	/*
	*	Operator declarations
	*/
	GradPath_c & operator=(GradPath_c const & rhs);
	GradPath_c & operator=(GradPathBase_c const & rhs);
	Boolean_t operator==(GradPath_c const & rhs) const;


	/*
	*	Getter methods
	*/

	GPTerminate_e GetGPTermType() { return m_HowTerminate; }
	vec3 GetStartPoint() { return m_StartPoint; };
	GPType_e GetGPType() { return m_GPType; }
	double GetTermPointRadius() { return sqrt(m_TermPointRadiusSqr); }
	double GetTermValue() { return m_TermValue; }
	VolExtentIndexWeights_s GetVolInfo() { return m_ODE_Data.VolZoneInfo; }
	FieldDataPointer_c GetRhoPtr() { return m_ODE_Data.RhoPtr; }
	vec3 GetTermPoint() { return m_TermPoint; }
	vector<CPTypeNum_e> GetTerminalCPTypeNums() { return m_TerminalCPTypeNums; }
	FESurface_c const * GetSurfPtr() { return m_Surface; }
	StreamDir_e GetDir() { return m_ODE_Data.Direction; }


	/*
	*	Setter methods
	*/
	Boolean_t SetupGradPath(vec3 const & StartPoint,
		StreamDir_e const & Direction,
		int NumGPPoints,
		GPTerminate_e const & HowTerminate,
		vec3 * TermPoint,
		vector<FieldDataPointer_c> const & CPXYZPtrs,
		int * NumCPs,
		double * TermPointRadius,
		double * TermValue,
		vector<int> const & MaxIJK,
		vec3 const & MaxXYZ,
		vec3 const & MinXYZ,
		vector<FieldDataPointer_c> const & GradPtrs,
		FieldDataPointer_c const & RhoPtr);

	Boolean_t SetupGradPath(vec3 const & StartPoint,
		StreamDir_e const & Direction,
		int NumGPPoints,
		GPType_e const & GPType,
		GPTerminate_e const & HowTerminate,
		vec3 * TermPoint,
		CritPoints_c * CPs,
		double * TermPointRadius,
		double * TermValue,
		VolExtentIndexWeights_s & VolInfo,
		vector<FieldDataPointer_c> const & HessPtrs,
		vector<FieldDataPointer_c> const & GradPtrs,
		FieldDataPointer_c const & RhoPtr,
		FESurface_c const * Surf = nullptr);

	Boolean_t Seed(bool const DoResample = true);
	Boolean_t SetMixingFactor(double const & MixFactor){
		if (MixFactor >= 0.0 && MixFactor <= 1.0)
			m_DirMixFactor = MixFactor;
		else
			return FALSE;

		return TRUE;
	}

	void Clear();

	void SetNumGPPoints(int NumPts) { m_NumGPPoints = NumPts; }

	void SetGPTermType(GPTerminate_e Type) { m_HowTerminate = Type; }
	void SetStartPoint(vec3 const & pt);
	void SetGPType(GPType_e Type) { m_GPType = Type; }
	void SetTermPointRadius(double Rad) { m_TermPointRadiusSqr = Rad * Rad; }
	void SetTermValue(double Val) { m_TermValue = Val; }
	void SetVolInfo(VolExtentIndexWeights_s VolInfo) { m_ODE_Data.VolZoneInfo = VolInfo; }
	void SetRhoPtr(FieldDataPointer_c & RhoPtr) { m_ODE_Data.RhoPtr = RhoPtr; }
	void SetTermPoint(vec3 const & pt) { m_TermPoint = pt; }
	void SetTerminalCPTypeNum(CPTypeNum_e CPTypeNum);
	void SetTerminalCPTypeNums(vector<CPTypeNum_e> const & CPTypeNums);
	void SetSurfPtr(FESurface_c const * SurfPtr, int SurfProjectionFrequency = -1) { m_Surface = SurfPtr; m_SurfGPProjectionFrequency = SurfProjectionFrequency; }
	void SetDir(StreamDir_e Dir) { m_ODE_Data.Direction = Dir; }

	void MakeRhoValuesMonotomic(VolExtentIndexWeights_s * VolInfo = nullptr, FieldDataPointer_c * RhoPtr = nullptr);


	Boolean_t ReinterpolateRhoValuesFromVolume(VolExtentIndexWeights_s * VolInfo = nullptr, FieldDataPointer_c * RhoPtr = nullptr);
	Boolean_t AlignToOtherPath(GradPath_c const & rhs, vector<int> AlignmentPointInds = {}, vector<int> ProtectedPointInds = {});

	Boolean_t ProjectPathToSurface(FESurface_c const * SurfPtr = nullptr);

private:
	//Boolean_t SetIndexAndWeightsForPoint(vec3 & Point);
	double RhoByCurrentIndexAndWeights();
	Boolean_t SeedInDirection(StreamDir_e const & Direction);

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

	vector<int> m_CPCoincidentPoints;

	// If terminating at:
	// Point
	vec3 m_TermPoint;
	// Any CP using tecplot zones
	vector<FieldDataPointer_c> m_CPXYZPtrs;
	int m_NumCPs;
	// ...using CritPoints_c
	CritPoints_c * m_CPs = nullptr;
	// if surface grad path
	FESurface_c const * m_Surface;
	int m_SurfGPProjectionFrequency = -1;
	// For point or CP:
	double m_TermPointRadiusSqr = 0;
	// Variable value
	// Which item in dep var list used for terminal value?
	double m_TermValue = -1;

	vector<CPTypeNum_e> m_TerminalCPTypeNums;

	GPTerminate_e m_HowTerminate;
	vec3 m_StartPoint;
};



Boolean_t GPsStraddleIB(GradPath_c const & GP1,
	GradPath_c const & GP2,
	double const & IBCheckAngle,
	double const & IBCheckDistRatio);

Boolean_t CPInNormalPlane(vec3 & StartPt, vec3 const & PlaneBasis, MultiRootObjects_s & MR);


class NEBGradPath_c : public GradPathBase_c
{
public:
	NEBGradPath_c();
	~NEBGradPath_c();

	NEBGradPath_c(vec3 const & StartPt, 
		vec3 const & EndPt, 
		unsigned int NumPts);

	Boolean_t Relax(double const & StepRatio,
		double const & Tol, 
		unsigned int const MaxIter,
		MultiRootParams_s & Params);

private:

};


GradPath_c ConcatenateResample(vector<GradPath_c> GPList, int NumPoints, vector<int> GPNumPointsList = {}, GPResampleMethod_e Method = GPResampleMethod_Adaptive);

/*
*	Redo of path stitching algorithm, formally TriPolyLines
*/
void StitchPaths(
	vector<int> const &     L,       // indices of points in P
	vector<int> const &     R,
	vector<vec3> const &     P,
	vector<vector<int> > &     T       // triplets of integers specifying nodes of triangles
	);
void StitchCapPaths(
	vector<int> const &     L,       // indices of points in P
	vector<int> const &     R,
	vector<int> &			C,       // indices of points in the cap, C
	vector<vec3> const &     P,
	vector<vector<int> > &     T,       // triplets of integers specifying nodes of triangles
	int LeftPathDeviationPointInd = -1,
	int RightPathDeviationPointInd = -1,
	int CapMidPointInd = -1
	);

