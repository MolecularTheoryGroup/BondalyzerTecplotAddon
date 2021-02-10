

#include "TECADDON.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

#include <list>
#include <deque>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <queue>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>

#include <Set.h>

#include "CSM_DATA_TYPES.h"
#include "CSM_CALC_VARS.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GRAD_PATH.h"

#include <armadillo>
using namespace arma;
using namespace tecplot::toolbox;

using std::list;
using std::vector;
using std::string;
using std::to_string;



int F2DGrad(gsl_vector const * pos, void * params, gsl_vector * GradValues);
int DF2DGrad(gsl_vector const * pos, void * params, gsl_matrix * Jacobian);
int FDF2DGrad(gsl_vector const * pos, void * params, gsl_vector * GradValues, gsl_matrix * Jacobian);

/*
*	Functions for GSL ODE solver.
*	These are not a member functions because it's a pain to do
*	member function with GSL.
*/
int GP_ODE_Gradient(double t, double const pos[], double dydt[], void* params);
int GP_ODE_Jacobian(double t, double const pos[], double *dfdy, double dydt[], void* params);

/*
*	GradPathParams_s methods
*/
GradPathParams_s & GradPathParams_s::operator=(GradPathParams_s const & rhs)
{
	if (this == &rhs)
		return *this;

	GradPtrs = rhs.GradPtrs;
	HasGrad = rhs.HasGrad;

	HessPtrs = rhs.HessPtrs;
	HasHess = rhs.HasHess;

	RhoPtr = rhs.RhoPtr;

	VolZoneInfo = rhs.VolZoneInfo;

	Direction = rhs.Direction;

	return *this;
}//	GradPathParams_s & GradPathParams_s::operator=(GradPathParams_s const & rhs)
Boolean_t const GradPathParams_s::operator==(GradPathParams_s const & rhs) const
{
	Boolean_t AreEqual = (

		GradPtrs == rhs.GradPtrs &&
		HasGrad == rhs.HasGrad &&

		HessPtrs == rhs.HessPtrs &&
		HasHess == rhs.HasHess &&

		RhoPtr == rhs.RhoPtr &&

		VolZoneInfo == rhs.VolZoneInfo &&

		Direction == rhs.Direction
		);

	return AreEqual;
}// Boolean_t GradPathParams_s::operator==(GradPathParams_s const & rhs) const


/*
*	GradPathBase_c methods
*/


GradPathBase_c::GradPathBase_c()
{
	m_StartEndCPNum[0] = m_StartEndCPNum[1] = -1;
	m_GradPathReady = FALSE;
	m_GradPathMade = FALSE;
	m_NumGPPoints = -1;
	m_Length = -1;
}

GradPathBase_c::~GradPathBase_c()
{
}



/*
*	Constructor for loading an existing 1-D ordered
*	zone as a GradPath_c
*/
GradPathBase_c::GradPathBase_c(EntIndex_t ZoneNum,
	vector<EntIndex_t> const & XYZRhoVarNums,
	AddOn_pa const & AddOnID)
{
	m_GradPathMade = FALSE;
	m_ZoneNum = ZoneNum;
	Boolean_t IsOk = TecUtilZoneIsOrdered(ZoneNum);
	if (!IsOk){
		TecUtilDialogErrMsg("Can't import grad path: Zone is not ordered.");
		TecUtilLockFinish(AddOnID);
		return;
	}
// 	if (!TecUtilZoneIsActive(ZoneNum)){
// 		Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// 		TecUtilSetAddMember(TmpSet, ZoneNum, FALSE);
// 		TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);
// 		TecUtilSetDealloc(&TmpSet);
// 	}

	LgIndex_t MaxIJK[3];
	TecUtilZoneGetIJK(ZoneNum, &MaxIJK[0], &MaxIJK[1], &MaxIJK[2]);
	IsOk = (MaxIJK[0] > 1 && MaxIJK[1] == 1 && MaxIJK[2] == 1);
	if (!IsOk){
		TecUtilDialogErrMsg("Can't import grad path: Zone is not 1-D.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	vector<FieldDataPointer_c> RawPtrs(4);
	for (int i = 0; i < XYZRhoVarNums.size(); ++i){
		RawPtrs[i].GetReadPtr(ZoneNum, XYZRhoVarNums[i]);
	}

	m_XYZList.resize(MaxIJK[0]);
	for (int i = 0; i < 4; ++i){
		if (i < 3){
			for (int j = 0; j < MaxIJK[0]; ++j)
				m_XYZList[j][i] = RawPtrs[i][j];
		}
		else{
			m_RhoList.resize(MaxIJK[0]);
			if (RawPtrs[i].IsReady()){
				for (int j = 0; j < MaxIJK[0]; ++j)
					m_RhoList[j] = RawPtrs[i][j];
			}
		}
	}

	m_StartEndCPNum[0] = m_StartEndCPNum[1] = -1;
	string tmpStr;
	for (int i = 0; i < 2; ++i){
		if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.CC.GPEndNumStrs[i], tmpStr))
			m_StartEndCPNum[i] = stoi(tmpStr) - 1;
	}

	m_GradPathMade = TRUE;
	m_GradPathReady = TRUE;
	return;
}



/*
*	Operator declarations
*/
GradPathBase_c & GradPathBase_c::operator=(GradPathBase_c const & rhs)
{
	if (this == &rhs)
		return *this;

	m_XYZList = rhs.m_XYZList;

	m_RhoList = rhs.m_RhoList;

	if (m_XYZList.size() == 0){
		m_XYZList.reserve(rhs.m_XYZList.capacity());
		m_RhoList.reserve(rhs.m_RhoList.capacity());
	}

	for (int i = 0; i < 2; ++i)
		m_StartEndCPNum[i] = rhs.m_StartEndCPNum[i];

	m_GradPathReady = rhs.m_GradPathReady;
	m_GradPathMade = rhs.m_GradPathMade;

	m_ZoneNum = rhs.m_ZoneNum;

	m_NumGPPoints = rhs.m_NumGPPoints;

	m_Length = rhs.m_Length;

	return *this;
}// GradPathBase_c & GradPathBase_c::operator=(GradPathBase_c const & rhs)

Boolean_t GradPathBase_c::IsSame(GradPathBase_c const & rhs) const
{
	Boolean_t AreSame = (
		m_RhoList == rhs.m_RhoList &&

		m_ZoneNum == rhs.m_ZoneNum &&

		m_NumGPPoints == rhs.m_NumGPPoints &&

		m_Length == rhs.m_Length &&

		m_GradPathReady == rhs.m_GradPathReady &&
		m_GradPathMade == rhs.m_GradPathMade
		);

	if (AreSame){
		for (int i = 0; i < m_XYZList.size() && AreSame; ++i)
			AreSame = sum(m_XYZList[i] == rhs.m_XYZList[i]) == 3;
	}

	return AreSame;
}
Boolean_t GradPathBase_c::operator==(GradPathBase_c const & rhs) const
{
	return IsSame(rhs);
}
GradPathBase_c & GradPathBase_c::operator+=(GradPathBase_c const & rhs)
{
	if (m_XYZList.size() > 0 && rhs.GetCount() > 0){
		double MinSqrDist = DistSqr(m_XYZList.back(), rhs.m_XYZList[0]);
		int MinPair = 1;
		double TempSqrDist = DistSqr(m_XYZList.back(), rhs.m_XYZList.back());
		if (MinSqrDist > TempSqrDist){
			MinSqrDist = TempSqrDist;
			MinPair = 2;
		}
		TempSqrDist = DistSqr(m_XYZList[0], rhs.m_XYZList.back());
		if (MinSqrDist > TempSqrDist){
			MinSqrDist = TempSqrDist;
			MinPair = 3;
		}
		TempSqrDist = DistSqr(m_XYZList[0], rhs.m_XYZList[0]);
		if (MinSqrDist > TempSqrDist){
			MinSqrDist = TempSqrDist;
			MinPair = 4;
		}

		int offset = 0; // used to 

		switch (MinPair){
			case 1:
				if (approx_equal(this->XYZAt(-1), rhs.XYZAt(0), "absdiff", 0.001)) 
					offset = 1;
				m_XYZList.insert(m_XYZList.end(), rhs.m_XYZList.cbegin() + offset, rhs.m_XYZList.cend());
				m_RhoList.insert(m_RhoList.end(), rhs.m_RhoList.cbegin() + offset, rhs.m_RhoList.cend());
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[1];
				break;
			case 2:
				if (approx_equal(this->XYZAt(-1), rhs.XYZAt(-1), "absdiff", 0.001))
					offset = -1;
				m_XYZList.insert(m_XYZList.end(), rhs.m_XYZList.crbegin(), rhs.m_XYZList.crend() + offset);
				m_RhoList.insert(m_RhoList.end(), rhs.m_RhoList.crbegin(), rhs.m_RhoList.crend() + offset);
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[0];
				break;
			case 3:
				if (approx_equal(this->XYZAt(0), rhs.XYZAt(-1), "absdiff", 0.001))
					offset = -1;
				m_XYZList.insert(m_XYZList.begin(), rhs.m_XYZList.cbegin(), rhs.m_XYZList.cend() + offset);
				m_RhoList.insert(m_RhoList.begin(), rhs.m_RhoList.cbegin(), rhs.m_RhoList.cend() + offset);
				m_StartEndCPNum[0] = rhs.m_StartEndCPNum[0];
				break;
			case 4:
				if (approx_equal(this->XYZAt(0), rhs.XYZAt(0), "absdiff", 0.001))
					offset = -1;
				m_XYZList.insert(m_XYZList.begin(), rhs.m_XYZList.crbegin(), rhs.m_XYZList.crend() + offset);
				m_RhoList.insert(m_RhoList.begin(), rhs.m_RhoList.crbegin(), rhs.m_RhoList.crend() + offset);
				m_StartEndCPNum[0] = rhs.m_StartEndCPNum[1];
				break;
		}
	}
	else if (rhs.GetCount() > 0){
		m_XYZList.insert(m_XYZList.end(), rhs.m_XYZList.cbegin(), rhs.m_XYZList.cend());
		m_RhoList.insert(m_RhoList.end(), rhs.m_RhoList.cbegin(), rhs.m_RhoList.cend());
		for (int i = 0; i < 2; ++i)
			m_StartEndCPNum[i] = rhs.m_StartEndCPNum[i];
	}
	return *this;
}
GradPathBase_c GradPathBase_c::operator+(GradPathBase_c const & rhs) const
{
	return GradPathBase_c(*this) += rhs;
}


/*
*	Getter methods
*/

double GradPathBase_c::GetLength() {
// 	if (m_GradPathMade && m_Length >= 0)
// 		return m_Length;

	m_Length = 0.0;

	for (int i = 0; i < GetCount() - 1; ++i)
		m_Length += Distance(m_XYZList[i], m_XYZList[i + 1]);

	return m_Length;
}

double GradPathBase_c::GetLength(int AtInd) const {

	double Length = 0.0;


	if (AtInd < 0) {
		AtInd += GetCount();
	}
	AtInd = CLAMP(AtInd, 0, GetCount() - 1);

	for (int i = 0; i < AtInd; ++i)
		Length += Distance(m_XYZList[i], m_XYZList[i + 1]);

	return Length;
}

int GradPathBase_c::GetIndAtLength(double const & Length) const{
	int i = 0;
	if (Length > 0.0){
		double checkLength = 0.0;
		for (i = 1; i < GetCount() - 1; ++i) {
			checkLength += Distance(m_XYZList[i-1], m_XYZList[i]);
			if (checkLength > Length){
				i--;
				break;
			}
		}
	}

	return i;
}

GradPathBase_c GradPathBase_c::SubGP(int BegPt, int EndPt) const{
	BegPt = GetInd(BegPt);
	EndPt = GetInd(EndPt);

	GradPathBase_c Out = *this;
	Out.m_XYZList.assign(m_XYZList.begin() + BegPt, m_XYZList.begin() + EndPt + 1);
	Out.m_RhoList.assign(m_RhoList.begin() + BegPt, m_RhoList.begin() + EndPt + 1);
	Out.GradPathBase_c::m_NumGPPoints = Out.m_XYZList.size();

	Out.m_GradPathMade = true;

	return Out;
}

vec3 GradPathBase_c::ClosestPoint(vec3 const & CheckPt) const
{
	int iJunk;

	return ClosestPoint(CheckPt, iJunk);
}

vec3 GradPathBase_c::ClosestPoint(vec3 const & CheckPt, int & PtNum) const
{
	return ClosestPointToPath(m_XYZList, CheckPt, PtNum);
}

vec3 GradPathBase_c::ClosestMaxCurvaturePoint(vec3 const & CheckPt, int & PtInd) const{
	vector<double> Curvatures(m_XYZList.size() - 2);
	vector<vec3> Segments(m_XYZList.size() - 1);

#pragma omp parallel for
	for (int i = 0; i < Segments.size(); ++i)
		Segments[i] = m_XYZList[i + 1] - m_XYZList[i];

#pragma omp parallel for
	for (int i = 0; i < Curvatures.size(); ++i)
		Curvatures[i] = VectorAngle(Segments[i], Segments[i+1]) / (norm(Segments[i]) + norm(Segments[i+1]));

	double MinDistSqr = DBL_MAX;
	PtInd = 0;
	for (int i = 1; i < Curvatures.size() - 1; ++i){
		if (Curvatures[i] > Curvatures[i-1] && Curvatures[i] > Curvatures[i+1]){
			double TmpDistSqr = DistSqr(CheckPt, m_XYZList[i + 1]);
			if (TmpDistSqr < MinDistSqr){
				MinDistSqr = TmpDistSqr;
				PtInd = i + 1;
			}
		}
	}

	return m_XYZList[PtInd];
}


/*
 *	Adaptive resampling of gradien path based on line simplification
 *	https://bost.ocks.org/mike/simplify/
 */
struct Tri{
	int Ind[3];
	double Area;
 	double MidLength;

	bool operator==(Tri const & rhs){
		return (Ind[0] == rhs.Ind[0]
			&& Ind[1] == rhs.Ind[1]
			&& Ind[2] == rhs.Ind[2]);
	}

	Tri operator=(Tri const & rhs) {
		if (*this == rhs)
			return *this;

		for (int i = 0; i < 3; ++i) Ind[i] = rhs.Ind[i];
		Area = rhs.Area;
 		MidLength = rhs.MidLength;

		return *this;
	}

	double GetMidLength(vector<vec3> const & PointVec){
		vector<double> Lengths = { DistSqr(PointVec[Ind[0]], PointVec[Ind[1]]), DistSqr(PointVec[Ind[0]], PointVec[Ind[2]]), DistSqr(PointVec[Ind[1]], PointVec[Ind[2]]) };
		std::sort(Lengths.begin(), Lengths.end());
		double len = sqrt(Lengths[1]);
		return len;
	}

	void Initialize(vector<vec3> const & PointVec) {
		Area = TriArea(PointVec[Ind[0]], PointVec[Ind[1]], PointVec[Ind[2]]);
		MidLength = GetMidLength(PointVec);
	}

	void Initialize(int i1, int i2, int i3, vector<vec3> const & PointVec) {
		Ind[0] = i1;
		Ind[1] = i2;
		Ind[2] = i3;
		Initialize(PointVec);
	}

// 	bool operator>(Tri const & rhs){
// 		return Area > rhs.Area;
// 	}
// 	bool operator<(Tri const & rhs){
// 		return !(*this > rhs);
// 	}
};
using ResampleAdaptive_MyPair_t = Tri;
using ResampleAdaptive_MyContainer_t = vector<ResampleAdaptive_MyPair_t>;
// bool ResampleAdaptive_CompFunc(ResampleAdaptive_MyPair_t const & e1, ResampleAdaptive_MyPair_t const & e2) { return e1.first < e2.first; }

Boolean_t GradPathBase_c::ResampleAdaptive(int TargetNumPoints, vector<int> & ProtectedPoints) {
	this->Resample(TargetNumPoints * 2, ProtectedPoints, GPResampleMethod_Linear);

	int NumPoints = m_XYZList.size();
	int NumTris = TargetNumPoints - 2;
	ResampleAdaptive_MyContainer_t TriVec(NumPoints - 2);
	vector<bool> NodeIsIncluded(NumPoints, true),
		NodeIsProtected(NumPoints, false);

	for (int i : ProtectedPoints)
		NodeIsProtected[i] = true;

	ProtectedPoints.clear();

	int iEnd = m_XYZList.size() - 1;
#pragma omp parallel for
	for (int i = 1; i < iEnd; ++i){
		int im1 = i - 1;
		TriVec[im1].Initialize(im1, i, i + 1, m_XYZList);
	}

	double MaxSegLength = this->GetLength() / (double)TargetNumPoints * 2.0;

	auto ResampleAdaptive_AreaCompFunc = [](ResampleAdaptive_MyPair_t const & e1, ResampleAdaptive_MyPair_t const & e2) { return e1.Area > e2.Area; };
//   	auto ResampleAdaptive_LengthCompFunc = [](ResampleAdaptive_MyPair_t const & e1, ResampleAdaptive_MyPair_t const & e2) { return e1.MidLength < e2.MidLength; };

	std::priority_queue<ResampleAdaptive_MyPair_t, ResampleAdaptive_MyContainer_t, decltype(ResampleAdaptive_AreaCompFunc)> AreaQueue(ResampleAdaptive_AreaCompFunc, TriVec);
// 	std::priority_queue<ResampleAdaptive_MyPair_t, ResampleAdaptive_MyContainer_t, decltype(ResampleAdaptive_LengthCompFunc)> LengthQueue(ResampleAdaptive_LengthCompFunc, TriVec);
	while (NumPoints > TargetNumPoints && !AreaQueue.empty()){
		while (NodeIsProtected[AreaQueue.top().Ind[1]]) 
			AreaQueue.pop();

		auto MinAreaTri = AreaQueue.top();// ,
// 			MaxLengthTri = LengthQueue.top();
		AreaQueue.pop();

// 		if (MinAreaTri == MaxLengthTri){
// 			LengthQueue.pop();
// 			continue;
// 		}
// 		
		if (MinAreaTri.MidLength > MaxSegLength && AreaQueue.size() > TargetNumPoints)
			continue;

		// See if MinAreaTri is a valid triangle (i.e. all three of its vertices are 'true' in
		// NodeIsIncluded).
		bool TriIsValid = true;
		for (int i = 0; i < 3 && TriIsValid; ++i)
			TriIsValid = NodeIsIncluded[MinAreaTri.Ind[i]];

		if (!TriIsValid)
			continue;

		// MinAreaTri is valid, so remove it and add two (or one) new triangles to the queue
		NodeIsIncluded[MinAreaTri.Ind[1]] = false;
		NumPoints--;

		// Get indices of closest valid triangles to the left and right.
		int LeftInd = MinAreaTri.Ind[0],
			RightInd = MinAreaTri.Ind[2];
		while (LeftInd > 0 && !NodeIsIncluded[LeftInd]){
			LeftInd--;
		}
		while (RightInd < m_XYZList.size() - 1 && !NodeIsIncluded[RightInd]){
			RightInd++;
		}

		if (LeftInd > 0){
			auto * Tri = &TriVec[LeftInd - 1];
			Tri->Ind[2] = RightInd;
			Tri->Initialize(m_XYZList);
			AreaQueue.push(*Tri);
// 			LengthQueue.push(*Tri);
		}
		if (RightInd < m_XYZList.size() - 1) {
			auto * Tri = &TriVec[RightInd - 1];
			Tri->Ind[0] = LeftInd;
			Tri->Initialize(m_XYZList);
			AreaQueue.push(*Tri);
// 			LengthQueue.push(*Tri);
		}
	}

	// Now the true elements of NodeIsIncluded form the path
	vector<vec3> NewXYZ;
	vector<double> NewRho;

	NewXYZ.reserve(TargetNumPoints);
	NewRho.reserve(TargetNumPoints);

	int NumTrue = 0;
	for (auto i : NodeIsIncluded)
		NumTrue += (int)i;



	for (int i = 0; i < NodeIsIncluded.size(); ++i){
		if (NodeIsIncluded[i]) {
			if (NodeIsProtected[i])
				ProtectedPoints.push_back(NewXYZ.size());
			NewXYZ.push_back(m_XYZList[i]);
			NewRho.push_back(m_RhoList[i]);
		}
	}

	m_XYZList = NewXYZ;
	m_RhoList = NewRho;

	m_XYZList.shrink_to_fit();
	m_RhoList.shrink_to_fit();

	return TRUE;
}

/*
 *	Redefines m_XYZList according to the nodes of rhs, but along the current m_XYZlist.
 *	That is, for each point i on rhs, define this-> point i as the closest point
 *	along the path to point i.
 */
Boolean_t GradPath_c::AlignToOtherPath(GradPath_c const & rhs, vector<int> AlignmentPointInds, vector<int> ProtectedPointInds){
	REQUIRE(this->IsMade() && rhs.IsMade() && !AlignmentPointInds.empty());
	GradPath_c thiscopy = GradPath_c(*this);
	if (this->GetCount() != rhs.GetCount()) {
		this->m_XYZList.resize(rhs.GetCount());
		this->m_RhoList.resize(rhs.GetCount());
		this->m_XYZList[0] = thiscopy[0];
		this->m_XYZList.back() = thiscopy[-1];
		this->m_RhoList[0] = thiscopy.RhoAt(0);
		this->m_RhoList.back() = thiscopy.RhoAt(-1);
	}

	GradPath_c NewGP;

	vector<vec3> ProtectedPoints;
	for (int i : ProtectedPointInds){
		ProtectedPoints.push_back(this->XYZAt(i));
	}

	// 1. Find the closest points on *this to the constrained points of rhs.
	// 2. Reconstruct *this as the concatenated set of gradient paths starting/ending
	//    at the closest points (or at *this's original starting/ending points).
	// At least one AlignmentPointInds is required, so at least two paths will be
	// made and concatenated. The first and last paths will be made explicitly, and 
	// the rest made in a loop.
	// 
	// 
	
	AlignmentPointInds.insert(AlignmentPointInds.begin(), 0);
	AlignmentPointInds.push_back(rhs.GetCount() - 1);

// 	{
// 		// Do the first path.
// 		// First, get *this's closest point to constraint point
// 		int ClosestInd;
// 		vec3 ClosestPt = this->ClosestPoint(rhs[AlignmentPointInds.front()], ClosestInd);
// 		// Get length of new path
// 		double tmpLength = this->GetLength(ClosestInd) + Distance(this->XYZAt(ClosestInd), ClosestPt);
// 		int NumPts = AlignmentPointInds.front();
// 		// Step size for getting new points
// 		double ArcLenStep = tmpLength / double(NumPts - 1);
// 		// Prepare workspace
// 		double CurLen = 0.0;
// 		vector<vec3> tmpPts(NumPts);
// 		// Add first point
// 		tmpPts.front() = this->XYZAt(0);
// 		// Add all but last point
// 		for (int i = 1; i < NumPts - 1; ++i) {
// 			CurLen += ArcLenStep; // Incremented arc length at which to get new point
// 			int CurInd = this->GetIndAtLength(CurLen); // index of point BEFORE CurLen is met
// 			double IndLen = this->GetLength(CurInd); // exact length at CurInd
// 			double LenDiff = CurLen - IndLen; // length discrepancy
// 			vec3 NewPt;
// 			if (LenDiff > 0.0) {
// 				// Linearly interpolate down current segment of path to achieve desired length
// 				double SegLen = Distance(this->XYZAt(CurInd), this->XYZAt(CurInd + 1));
// 				double LenDiffRatio = LenDiff / SegLen; // SegLen cannot be zero because of how CurInd was determined
// 				vec3 Seg = this->XYZAt(CurInd + 1) - this->XYZAt(CurInd);
// 				NewPt = this->XYZAt(CurInd) + Seg * LenDiffRatio;
// 			}
// 			else {
// 				// CurLen corresponds perfectly to one of the existing points, so just use that
// 				NewPt = this->XYZAt(CurInd);
// 			}
// 		}
// 		// Add last point
// 		tmpPts.push_back(NewPt);
// 		// Add to NewGP
// 		NewGP += GradPath_c(tmpPts);
// 	}
// 	
	
	// Do each path segment
	for (int i = 0; i < AlignmentPointInds.size() - 1; ++i){
		int j = i + 1;
		// Do the first path.
		// First, get *this's closest points to constraint points
		int ClosestIndI, ClosestIndJ;
		vec3 ClosestPtI, ClosestPtJ;
		if (i == 0) {
			ClosestPtI = this->XYZAt(0);
			ClosestIndI = 0;
		}
		else {
			ClosestPtI = this->ClosestPoint(rhs[AlignmentPointInds[i]], ClosestIndI);
		}
		if (j == AlignmentPointInds.size() - 1){
			ClosestPtJ = this->XYZAt(-1);
			ClosestIndJ = this->GetCount() - 1;
		}
		else{
			ClosestPtJ = this->ClosestPoint(rhs[AlignmentPointInds[j]], ClosestIndJ);
		}
		// Get length of new path
		double LenI = this->GetLength(ClosestIndI) + Distance(this->XYZAt(ClosestIndI), ClosestPtI);
		double LenJ = this->GetLength(ClosestIndJ) + Distance(this->XYZAt(ClosestIndJ), ClosestPtJ);
		double tmpLength = LenJ - LenI;
		int NumPts = AlignmentPointInds[j] - AlignmentPointInds[i] + 1;
		// Step size for getting new points
		double ArcLenStep = tmpLength / double(NumPts - 1);
		// Prepare workspace
		vector<vec3> tmpPts;
		tmpPts.reserve(NumPts);
		// Add first point
		tmpPts.push_back(ClosestPtI);
		// Add all but last point
		double CurLen = LenI;
		for (int i = 1; i < NumPts - 1; ++i) {
			CurLen += ArcLenStep; // Incremented arc length at which to get new point
			int CurInd = this->GetIndAtLength(CurLen); // index of point BEFORE CurLen is met
			double IndLen = this->GetLength(CurInd); // exact length at CurInd
			double LenDiff = CurLen - IndLen; // length discrepancy
			vec3 NewPt;
			if (LenDiff > 0.0) {
				// Linearly interpolate down current segment of path to achieve desired length
				double SegLen = Distance(this->XYZAt(CurInd), this->XYZAt(CurInd + 1));
				double LenDiffRatio = LenDiff / SegLen; // SegLen cannot be zero because of how CurInd was determined
				vec3 Seg = this->XYZAt(CurInd + 1) - this->XYZAt(CurInd);
				NewPt = this->XYZAt(CurInd) + Seg * LenDiffRatio;
			}
			else {
				// CurLen corresponds perfectly to one of the existing points, so just use that
				NewPt = this->XYZAt(CurInd);
			}
			tmpPts.push_back(NewPt);
		}
		// Add last point
		tmpPts.push_back(ClosestPtJ);
		// Add to NewGP
		NewGP += GradPath_c(tmpPts);
	}

	this->m_XYZList = NewGP.m_XYZList;

	for (vec3 p : ProtectedPoints){
		int MinDistInd;
		this->ClosestPoint(p, MinDistInd);
		this->m_XYZList[MinDistInd] = p;
	}


// 	AlignmentPointInds.insert(AlignmentPointInds.begin(), 0);
// 	AlignmentPointInds.push_back(rhs.GetCount() - 1);
// 
// 	for (int i = 1; i < AlignmentPointInds.size(); ++i) {
// 		auto tmpGP = this->SubGP(AlignmentPointInds[i - 1], AlignmentPointInds[i]);
// 		tmpGP.Resample(AlignmentPointInds[i] - AlignmentPointInds[i - 1] + 1, GPResampleMethod_Linear);
// 		NewGP += tmpGP;
// 	}
// 
// 	this->m_XYZList = NewGP.m_XYZList;

// 		for (int i = 1; i < this->GetCount() - 1; ++i) {
// 			this->m_XYZList[i] = thiscopy.ClosestPoint(rhs[i]);
// 		}

	return TRUE;
}

// Removes occurrances of having neighboring segments with > 90 degree turns
// (gradient paths shouldn't do that)
void GradPathBase_c::RemoveKinks(double const & AngleCutoff){
	int NumPoints = m_XYZList.size();
	vector<vec3> Segments(NumPoints - 1);
#pragma omp parallel for
	for (int i = 1; i < NumPoints; ++i)
		Segments[i - 1] = m_XYZList[i] - m_XYZList[i - 1];

	for (int i = 2; i < NumPoints - 1; ++i){
		if (VectorAngle(Segments[i], Segments[i-1]) > AngleCutoff
			&& VectorAngle(Segments[i-1],Segments[i-2]) > AngleCutoff)
		{
			vec3 MidPt = (m_XYZList[i - 1] + m_XYZList[i]) * 0.5;
			m_XYZList[i - 1] = MidPt;
			m_XYZList[i] = MidPt;
		}
	}
}


Boolean_t GradPathBase_c::Resample(int NumPoints, vector<int> & ProtectedPoints, GPResampleMethod_e Method){
	Boolean_t IsOk = m_GradPathMade && NumPoints > 3;

	if (Method == GPResampleMethod_Adaptive && NumPoints < m_XYZList.size())
		return ResampleAdaptive(NumPoints, ProtectedPoints);

	vector<int> NewProtectedPoints;
	vector<vec3> NewXYZList;
	vector<double> NewRhoList;

	int OldCount = GetCount();

	// 	if (IsOk && NumPoints < OldCount){
	if (IsOk) {
		NewXYZList.resize(NumPoints);

		NewRhoList.resize(NumPoints);

		double DelLength = GetLength() / static_cast<double>(NumPoints - 1);

		double ArcLength = 0.0,
			ArcLengthI = 0.0,
			ArcLengthIm1 = 0.0;

		vec3 PtI, PtIm1;
		double RhoI, RhoIm1;

		PtI = m_XYZList[0];
		RhoI = m_RhoList[0];

		NewXYZList[0] = PtI;

		NewRhoList[0] = RhoI;

		int OldI = 0;

		for (int NewI = 1; NewI < NumPoints - 1; ++NewI) {
			ArcLength += DelLength;

			while (OldI < OldCount - 1 && ArcLengthI < ArcLength) {
				++OldI;

				ArcLengthIm1 = ArcLengthI;
				PtIm1 = PtI;
				RhoIm1 = RhoI;

				PtI = m_XYZList[OldI];
				RhoI = m_RhoList[OldI];

				ArcLengthI += Distance(PtI, PtIm1);

				if (std::find(ProtectedPoints.begin(), ProtectedPoints.end(), OldI) != ProtectedPoints.end()) {
					NewProtectedPoints.push_back(NewI);
					break;
				}
			}

			double Ratio;

			if (std::find(NewProtectedPoints.begin(), NewProtectedPoints.end(), OldI) != NewProtectedPoints.end()) {
				NewXYZList[NewI] = m_XYZList[OldI];
				NewRhoList[NewI] = m_RhoList[OldI];
			}
			else {
				Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
				NewXYZList[NewI] = PtIm1 + (PtI - PtIm1) * Ratio;
				NewRhoList[NewI] = RhoIm1 + Ratio * (RhoI - RhoIm1);
			}


			if (OldI >= OldCount) {
				while (NewI < NumPoints) {
					NewI++;
					if (NewI < NumPoints) {
						NewXYZList[NewI] = PtIm1 + (PtI - PtIm1) * Ratio;

						NewRhoList[NewI] = RhoIm1 + Ratio * (RhoI - RhoIm1);
					}
				}
			}
		}

		/*
		*	Add last point
		*/

		NewXYZList[NumPoints - 1] = m_XYZList[OldCount - 1];

		NewRhoList[NumPoints - 1] = m_RhoList[OldCount - 1];


		IsOk = NewXYZList.size() == NumPoints;

		if (IsOk) {
			IsOk = (NewRhoList.size() == NumPoints);
		}

		if (IsOk) {
			m_XYZList.swap(NewXYZList);

			m_RhoList.swap(NewRhoList);
		}

		m_Length = -1;
	}
	else IsOk = FALSE;

	m_NumGPPoints = m_XYZList.size();

	m_XYZList.shrink_to_fit();
	m_RhoList.shrink_to_fit();

	ProtectedPoints = NewProtectedPoints;

	return IsOk;
}

/*
 *	Resample path so that all but the terminal points are those determined by RhoVals[1] through RhoVals[-2]
 */
Boolean_t GradPathBase_c::ResampleByRhoVals(vector<double> const & RhoVals, vector<int> const & ProtectedPoints){

	vector<vec3> ProtectedPts;
	ProtectedPts.reserve(ProtectedPoints.size());
	for (auto i : ProtectedPoints){
		if (i > 0 && i < this->GetCount() - 1) {
			ProtectedPts.push_back(this->XYZAt(i));
		}
	}
	Boolean_t IsOk = TRUE;
	vector<vec3> NewXYZList(RhoVals.size());
	int ind;
	double weight;
	for (int i = 1; i < RhoVals.size() - 1; ++i){
		bool isfound = false;
		for (int j = 0; j < this->m_RhoList.size() && !isfound; ++j){
			if (abs(RhoVals[i] - this->RhoAt(j)) < 1e-10){
				NewXYZList[i] = this->XYZAt(j);
				isfound = true;
			}
		}
		if (!isfound) {
			this->GetPointAtRhoValue(RhoVals[i], NewXYZList[i], ind, weight);
		}
	}

	vector<double> NewRhoList = RhoVals;
	NewRhoList.front() = this->RhoAt(0);
	NewRhoList.back() = this->RhoAt(-1);

	if (abs(NewRhoList.front() - NewRhoList[1]) > abs(NewRhoList.front() - NewRhoList[NewRhoList.size()-1])) {
		// Ordering of rho values was opposite that of the gradient path, so flip the end points to match.
		NewXYZList.front() = this->XYZAt(-1);
		NewRhoList.front() = this->RhoAt(-1);
		NewXYZList.back() = this->XYZAt(0);
		NewRhoList.back() = this->RhoAt(0);
	}
	else{
		NewXYZList.front() = this->XYZAt(0);
		NewXYZList.back() = this->XYZAt(-1);
	}

	for (int i = 0; i < NewXYZList.size() - 2; ++i) {
		if (sum(NewXYZList[i] == NewXYZList[i + 1]) == 3){
			NewXYZList[i + 1] = (NewXYZList[i] + NewXYZList[i + 2]) * 0.5;
		}
	}

	this->m_XYZList = NewXYZList;
	this->m_RhoList = NewRhoList;

	for (auto const & p : ProtectedPts){
		int ind;
		this->ClosestPoint(p, ind);
		this->m_XYZList[ind] = p;
	}

	return IsOk;
}

Boolean_t GradPathBase_c::Reverse(){
	Boolean_t IsOk = m_GradPathMade;
	if (IsOk){
// 		vector<vec3> NewXYZList;
// 		vector<double> NewRhoList;
// 
// 		NewXYZList.insert(NewXYZList.begin(), m_XYZList.crbegin(), m_XYZList.crend());
// 		m_XYZList = NewXYZList;
		std::reverse(m_XYZList.begin(), m_XYZList.end());

// 		NewRhoList.insert(NewRhoList.begin(), m_RhoList.crbegin(), m_RhoList.crend());
// 		m_RhoList = NewRhoList;
		std::reverse(m_RhoList.begin(), m_RhoList.end());

		int TmpInt = m_StartEndCPNum[0];
		m_StartEndCPNum[0] = m_StartEndCPNum[1];
		m_StartEndCPNum[1] = TmpInt;
	}
	return IsOk;
}

vector<int> GradPathBase_c::GetCPCoincidentPoints(CritPoints_c const * CPs, std::set<CPType_e> const & CPTypes, double tol) const {
	vector<int> outPoints;
	// Check distance to saddle CPs, don't care about max/min
	// keep track of the closest point to each saddle point, then 
	// check that distance against the tolerance.
	if (CPs != nullptr) {
		std::map<int, std::pair<int, double> > MinSaddlePointDist;
		for (int i = 0; i < m_XYZList.size(); ++i) {
			int CPOffSet;
			double MinDist;
			vec ClosestPt = CPs->ClosestPoint(m_XYZList[i], CPOffSet, MinDist);
			if (CPOffSet >= 0) {
				auto CPType = CPs->GetTypeFromTotOffset(CPOffSet);
				if (CPTypes.empty() || CPTypes.count(CPType)) {
					if (MinSaddlePointDist.count(CPOffSet)) {
						if (MinDist < MinSaddlePointDist[CPOffSet].second) {
							MinSaddlePointDist[CPOffSet] = std::make_pair(i, MinDist);
						}
					}
					else {
						MinSaddlePointDist[CPOffSet] = std::make_pair(i, MinDist);
					}
				}
			}
		}

		for (auto p : MinSaddlePointDist) {
			if (p.second.second <= tol) {
				outPoints.push_back(p.second.first);
			}
		}
	}
	return outPoints;

	// Find bends close to 90 deg, which only* happens when a GP passes through a saddle CP

// 	for (int i = 1; i < m_XYZList.size() - 1; ++i){
// 		if (VectorAngle(m_XYZList[i-1] - m_XYZList[i], m_XYZList[i+1] - m_XYZList[i]) < PIOVER2 * 1.1){
// 			outPoints.push_back(i);
// 		}
// 	}
// 
// 	return outPoints;
}

/*
*	Concatenates and resamples at once.
*	Resampling a concatenated grad path that has a
*	sharp bend can cut the corner at the bend, so
*	this method avoids that.
*/
GradPathBase_c & GradPathBase_c::ConcatenateResample(GradPathBase_c & rhs, int NumPoints, GPResampleMethod_e Method)
{
	int iJunk;

	ConcatenateResample(rhs, NumPoints, iJunk, Method);

	return *this;
}

GradPathBase_c & GradPathBase_c::ConcatenateResample(GradPathBase_c & rhs, int NumPoints, int & BrigePtNum, GPResampleMethod_e Method)
{
// 	double MyLength = GetLength();
// 	double rhsLength = rhs.GetLength();
// 	int MyNumPtsInit = GetCount();
// 	int rhsNumPtsInit = rhs.GetCount();
// 	double TotalLength = MyLength + rhsLength;
// 
// 	int MyNumPts = static_cast<int>((MyLength / TotalLength) * static_cast<double>(NumPoints));
// 	MyNumPts = MIN(MyNumPts, MyNumPtsInit);
// 	BrigePtNum = MyNumPts - 1;
// 	int rhsNumPts = MIN(NumPoints - MyNumPts, rhsNumPtsInit);
// 
// 	Resample(MyNumPts);
// 	GradPathBase_c NewRHS(rhs);
// 	NewRHS.Resample(rhsNumPts);
// 
// 	Concatenate(NewRHS);



	GradPathBase_c TmpGP = *this;
	TmpGP += rhs;

	TmpGP.RemoveKinks();

	vector<int> ProtectedPoints;
	if (sum(TmpGP[0] == this->XYZAt(0)) == 3)
		ProtectedPoints = { GetCount() - 1 };
	else
		ProtectedPoints = { rhs.GetCount() - 1 };

	TmpGP.Resample(NumPoints, ProtectedPoints, Method);

	BrigePtNum = ProtectedPoints.front();

	*this = TmpGP;

	return *this;
}

GradPath_c ConcatenateResample(vector<GradPath_c> GPList, int NumPoints, vector<int> GPNumPointsList, GPResampleMethod_e Method)
{
	vector<double> LenghtList(GPList.size());
	vector<int> NumPointsList(GPList.size());
	if (GPNumPointsList.size() != GPList.size())
		GPNumPointsList = vector<int>(GPList.size(), -1);

	double TotalLength = 0.0;

	int NewNumPoints = NumPoints;

	for (int gpi = 0; gpi < GPList.size(); ++gpi){
		if (GPNumPointsList[gpi] < 0){
			LenghtList[gpi] = GPList[gpi].GetLength();
			TotalLength += LenghtList[gpi];
		}
		else{
			NewNumPoints -= GPNumPointsList[gpi];
		}
	}

	for (int gpi = 0; gpi < GPList.size(); ++gpi) {
		if (GPNumPointsList[gpi] < 0) {
			NumPointsList[gpi] = LenghtList[gpi] / TotalLength * double(NewNumPoints);
		}
		else{
			NumPointsList[gpi] = GPNumPointsList[gpi];
		}
	}

	int TotNumPoints = 0;
	for (auto n : NumPointsList) TotNumPoints += n;

	if (TotNumPoints != NumPoints){
		// Through roundoff there is a discrepancy between the specified
		// number of points and the calculated number of points,
		// So fix it by modifying the first path that doesn't
		// have it's number of points explicitly specified.
		
		for (int gpi = 0; gpi < GPList.size(); ++gpi) {
			if (GPNumPointsList[gpi] < 0) {
				NumPointsList[gpi] = NumPointsList[gpi] - (TotNumPoints - NumPoints);
				break;
			}
		}
	}

	TotNumPoints = 0;

	GradPath_c GP;

	while (TotNumPoints != NumPoints) {
		auto GPListCopy = GPList;

		GP = GPListCopy[0];
		GP.Resample(NumPointsList[0], Method);
		for (int i = 1; i < GPListCopy.size(); ++i) {
			GPListCopy[i].Resample(NumPointsList[i], Method);
			GP += GPListCopy[i];
		}

		TotNumPoints = GP.GetCount();
		if (TotNumPoints != NumPoints) {
			int maxPathNum = 0;
			for (int i = 1; i < NumPointsList.size(); ++i) {
				if (NumPointsList[i] > NumPointsList[maxPathNum]) {
					maxPathNum = i;
				}
			}
			if (TotNumPoints > NumPoints) {
				
				NumPointsList[maxPathNum]--;
			}
			else if (TotNumPoints < NumPoints) {
				NumPointsList[maxPathNum]++;
			}
		}
	}

	return GP;
}


/*
*	For making a grad path terminate at a particular point.
*	In this version, at the intersection of the grad path and
*	a sphere of r = Radius around Point.
*/
Boolean_t GradPathBase_c::Trim(vec3 const & Point, double const & Radius)
{
	Boolean_t IsOk = IsMade();

	if (IsOk){
		int Count = GetCount();

		vec3 IntPoint;
		int IntPointNum = this->GetSphereIntersectionPoint(Point, Radius, IntPoint);
		if (IntPointNum > 0) {
			double RhoVal = m_RhoList[IntPointNum] + (Distance(IntPoint, m_XYZList[IntPointNum]) / Distance(m_XYZList[IntPointNum + 1], m_XYZList[IntPointNum])) * (m_RhoList[IntPointNum + 1] - m_RhoList[IntPointNum]);
			GradPathBase_c GP = this->SubGP(0, IntPointNum);
			GP.PointAppend(IntPoint, RhoVal);
			*this = GP;
			return TRUE;
		}
		else
			return FALSE;

// 		vector<vec3> NewXYZList;
// 
// 		vector<double> NewRhoList;
// 
// 		if (DistSqr(Point, (m_XYZList[0])) < DistSqr(Point, (m_XYZList[Count - 1]))){
// 			/*
// 			*	Initial point is closer to Point, so start from beginning
// 			*/
// 			double RhoI, RhoIm1 = m_RhoList[0];
// 			vec3 PtI, PtIm1 = m_XYZList[0];
// 			double DistI, DistIm1 = Distance(PtIm1, Point);
// 
// 			for (int i = 1; i < Count; ++i){
// 				PtI = m_XYZList[i];
// 				RhoI = m_RhoList[i];
// 				DistI = Distance(PtI, Point);
// 
// 				if ((DistI < Radius && DistIm1 >= Radius) ||
// 					(DistI >= Radius && DistIm1 < Radius)){
// 					// Now PtI and PtIm1 are straddling the sphere.
// 					// Interpolate to find the intersection point.
// 					double DistanceRatio = (Radius - DistI) / (DistIm1 - DistI);
// 
// 					vec3 NewPoint = (PtI + ((PtIm1 - PtI) * DistanceRatio)) - Point;
// 					double NewRho = (RhoI + ((RhoIm1 - RhoI) * DistanceRatio));
// 
// 					/*
// 					* Project intersection point to sphere surface
// 					*/
// 					// Compute spherical coordinates
// 					double   rad = norm(NewPoint);
// 					double theta = acos(NewPoint[2] / rad);
// 					double   phi = atan2(NewPoint[1], NewPoint[0]);
// 
// 					// Project point onto a sphere of radius "radius"
// 					NewPoint[0] = Radius * sin(theta) * cos(phi);
// 					NewPoint[1] = Radius * sin(theta) * sin(phi);
// 					NewPoint[2] = Radius * cos(theta);
// 
// 					NewPoint += Point;
// 
// 					NewXYZList.reserve(m_XYZList.size() - i + 1);
// 					NewXYZList.push_back(NewPoint);
// 					NewXYZList.insert(NewXYZList.end(), m_XYZList.cbegin() + i, m_XYZList.cend());
// 
// 					m_XYZList[i] = NewXYZList[i];
// 
// 					NewRhoList.reserve(m_RhoList.size() - i + 1);
// 					NewRhoList.push_back(NewRho);
// 					NewRhoList.insert(NewRhoList.end(), m_RhoList.cbegin() + i, m_RhoList.cend());
// 					m_RhoList = NewRhoList;
// 
// 					break;
// 				}
// 				PtIm1 = PtI;
// 				RhoIm1 = RhoI;
// 				DistIm1 = DistI;
// 			}
// 		}
// 		else{
// 			/*
// 			*	Terminal point is closer to Point, so start from end
// 			*/
// 			double RhoI, RhoIm1 = m_RhoList[Count - 1];
// 			vec3 PtI, PtIm1 = m_XYZList[Count - 1];
// 			double DistI, DistIm1 = Distance(PtIm1, Point);
// 
// 			for (int i = Count - 2; i >= 0; --i){
// 				PtI = m_XYZList[i];
// 				RhoI = m_RhoList[i];
// 				DistI = Distance(PtI, Point);
// 
// 				if ((DistI < Radius && DistIm1 >= Radius) ||
// 					(DistI >= Radius && DistIm1 < Radius)){
// 					// Now PtI and PtIm1 are straddling the sphere.
// 					// Interpolate to find the intersection point.
// 					double DistanceRatio = (Radius - DistI) / (DistIm1 - DistI);
// 
// 					vec3 NewPoint = (PtI + ((PtIm1 - PtI) * DistanceRatio)) - Point;
// 					double NewRho = (RhoI + ((RhoIm1 - RhoI) * DistanceRatio));
// 
// 					/*
// 					* Project intersection point to sphere surface
// 					*/
// 					// Compute spherical coordinates
// 					double   rad = sqrt(pow(NewPoint[0], 2) + pow(NewPoint[1], 2) + pow(NewPoint[2], 2));
// 					double theta = acos(NewPoint[2] / rad);
// 					double   phi = atan2(NewPoint[1], NewPoint[0]);
// 
// 					// Project point onto a sphere of radius "radius"
// 					NewPoint[0] = Radius * sin(theta) * cos(phi);
// 					NewPoint[1] = Radius * sin(theta) * sin(phi);
// 					NewPoint[2] = Radius * cos(theta);
// 
// 					NewPoint += Point;
// 
// 					NewXYZList.reserve(i + 1);
// 					NewXYZList.insert(NewXYZList.end(), m_XYZList.cbegin(), m_XYZList.cbegin() + i);
// 					NewXYZList.push_back(NewPoint);
// 
// 					m_XYZList = NewXYZList;
// 
// 					NewRhoList.reserve(i + 1);
// 					NewRhoList.insert(NewRhoList.end(), m_RhoList.cbegin(), m_RhoList.cbegin() + i);
// 					NewRhoList.push_back(NewRho);
// 					m_RhoList = NewRhoList;
// 
// 					break;
// 				}
// 				PtIm1 = PtI;
// 				RhoIm1 = RhoI;
// 				DistIm1 = DistI;
// 			}
// 		}
// 
// 		IsOk = (NewXYZList.size() > 0 && NewRhoList.size() > 0);
	}

	return IsOk;
}

void GradPathBase_c::Clear(){
	m_XYZList.clear();
	m_RhoList.clear();
	m_GradPathMade = FALSE;
	m_Length = -1;
	m_StartEndCPNum[0] = -1;
	m_StartEndCPNum[1] = -1;
	m_ZoneNum = 0;
}

void GradPathBase_c::PointAppend(vec3 const & Point, double const & Rho){

	m_XYZList.push_back(Point);

	if (m_Length >= 0) m_Length += Distance(m_XYZList[-1], m_XYZList[-2]);

	m_RhoList.push_back(Rho);

	m_NumGPPoints++;
}


void GradPathBase_c::PointPrepend(vec3 const & Point, double const & Rho){

	m_XYZList.insert(m_XYZList.begin(), Point);

	if (m_Length >= 0) m_Length += Distance(m_XYZList[0], m_XYZList[1]);

	m_RhoList.insert(m_RhoList.begin(), Rho);

	m_NumGPPoints++;
}

EntIndex_t GradPathBase_c::SaveAsOrderedZone(string const & ZoneName, ColorIndex_t const MeshColor)
{
	Boolean_t IsOk = m_GradPathMade;

	if (IsOk){
		EntIndex_t VolZoneNum = ZoneNumByName(string("Full Volume"));
		EntIndex_t RhoVarNum = VarNumByName(string("Electron Density"));

		vector<EntIndex_t> XYZVarNums(3);
		TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);

		IsOk = (VolZoneNum > 0 &&
			XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);

		if (IsOk){
			EntIndex_t NumVars = TecUtilDataSetGetNumVars();
			vector<FieldDataType_e> VarTypes(NumVars);
			for (int i = 0; i < NumVars; ++i){
				VarTypes[i] = TecUtilDataValueGetType(VolZoneNum, i + 1);
			}

			for (int i = 0; i < 3; ++i)
				VarTypes[XYZVarNums[i] - 1] = FieldDataType_Double;

			if (RhoVarNum > 0) VarTypes[RhoVarNum - 1] = FieldDataType_Double;

			IsOk = (SaveAsOrderedZone(ZoneName, VarTypes, XYZVarNums, RhoVarNum, TRUE, MeshColor) > 0);
		}
	}

	if (IsOk)
		return m_ZoneNum;
	else return -1;
}

EntIndex_t GradPathBase_c::SaveAsOrderedZone(string const & ZoneName,
	vector<FieldDataType_e> & VarDataTypes,
	vector<int> const & XYZVarNums,
	int RhoVarNum,
	Boolean_t const DoActivate,
	ColorIndex_t const MeshColor)
{
	Boolean_t IsOk = (m_XYZList.size() > 0 && m_XYZList.size() == m_RhoList.size());

	if (IsOk){

		if (IsOk){

			if (VarDataTypes.size() == TecUtilDataSetGetNumVars())
				IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), GetCount(), 1, 1, ZoneType_Ordered, VarDataTypes.data());
			else IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), GetCount(), 1, 1, ZoneType_Ordered, nullptr);

			if (IsOk){
				m_ZoneNum = TecUtilDataSetGetNumZones();

				vector<vector<double> > TmpValues(3, vector<double>(m_XYZList.size()));
				for (int i = 0; i < m_XYZList.size(); ++i){
					for (int j = 0; j < 3; ++j)
						TmpValues[j][i] = m_XYZList[i][j];
				}

				TecUtilDataLoadBegin();

				for (int i = 0; i < 3; ++i){
					FieldData_pa GPXYZPtr = TecUtilDataValueGetWritableNativeRef(m_ZoneNum, XYZVarNums[i]);
					TecUtilDataValueArraySetByRef(GPXYZPtr, 1, GetCount(), TmpValues[i].data());
				}
				if (RhoVarNum > 0){
					FieldData_pa RhoFDPtr = TecUtilDataValueGetWritableNativeRef(m_ZoneNum, RhoVarNum);
					TecUtilDataValueArraySetByRef(RhoFDPtr, 1, GetCount(), m_RhoList.data());
				}

				TecUtilDataLoadEnd();

				Set_pa TempSet = TecUtilSetAlloc(FALSE);
				TecUtilSetAddMember(TempSet, m_ZoneNum, FALSE);
				if (DoActivate) TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
				else TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
				TecUtilZoneSetScatter(SV_SHOW, TempSet, 0.0, FALSE);
				TecUtilZoneSetContour(SV_SHOW, TempSet, 0.0, FALSE);
				TecUtilZoneSetMesh(SV_SHOW, TempSet, 0.0, TRUE);
				TecUtilZoneSetMesh(SV_COLOR, TempSet, 0.0, MeshColor);
				TecUtilSetDealloc(&TempSet);

				AuxDataZoneSetItem(m_ZoneNum, CSMAuxData.CC.ZoneType, CSMAuxData.CC.ZoneTypeGP);
			}
		}
	}

	return m_ZoneNum;
}

Boolean_t GradPathBase_c::SaveAsCSV(string const & PathToFile, Boolean_t IncludeVars){
	Boolean_t IsOk = TRUE;

	std::ofstream Out(PathToFile.c_str(), std::ios::out | std::ios::trunc);
	IsOk = Out.is_open();

	if (IsOk){
		for (int i = 0; i < m_XYZList.size(); ++i){
			for (int j = 0; j < 3; ++j){
				Out << std::setprecision(16) << std::scientific << m_XYZList[i][j] << ',';
			}
			if (IncludeVars)
				Out << std::setprecision(16) << std::scientific << m_RhoList[i];
			Out << '\n';
		}
	}

	Out.close();

	return IsOk;
}

Boolean_t GradPathBase_c::TruncateAtRhoValue(double const & RhoVal){
	REQUIRE(this->IsMade());
	bool RhoIsAscending = (m_RhoList.front() < m_RhoList.back());
	if (RhoIsAscending)
		REQUIRE(RhoVal > m_RhoList.front() && RhoVal < m_RhoList.back());

	int Ind;
	double Weight;
	vec3 Pt;

	if (this->GetPointAtRhoValue(RhoVal, Pt, Ind, Weight)){
		GradPathBase_c GP;
		if (RhoIsAscending){
			GP = this->SubGP(Ind, -1);
			GP.PointPrepend(Pt, RhoVal);
		}
		else{
			GP = this->SubGP(0, Ind);
			GP.PointAppend(Pt, RhoVal);
		}

		*this = GP;

		return TRUE;
	}

	return FALSE;
}

vec3 GradPathBase_c::operator[](int i) const{
	return m_XYZList[GetInd(i)];
}

/*
 *	Provide a rho value and get the point on the gradient path that
 *	has that value.
 *	I'm lazy, so this is all done using std::set
 */
bool GradPathBase_c::GetPointAtRhoValue(double const & RhoValue, vec3 & OutVec, int & Ind, double & Weight) const {
	REQUIRE(IsMade());
	REQUIRE(m_RhoList.size() >= 2);
	REQUIRE(RhoValue > 0.0);

	bool RhoIsAscending = (m_RhoList[0] < m_RhoList.back());

	// If RhoValue is out of bounds of the values in m_RhoList, then just return whichever end is closer
	double MinRho = MIN(m_RhoList[0], m_RhoList.back()),
		MaxRho = MAX(m_RhoList[0], m_RhoList.back());
	if (RhoValue <= MinRho || RhoValue >= MaxRho) {
		if ((RhoIsAscending && RhoValue <= MinRho) || (!RhoIsAscending && RhoValue >= MaxRho))
			Ind = 0;
		else
			Ind = m_RhoList.size() - 1;

		OutVec = m_XYZList[Ind];
		Weight = 0;

		return false;
	}

	// Get set of rho values in ascending order
	std::set<double> RhoSet;
	if (RhoIsAscending) RhoSet = std::set<double>(m_RhoList.cbegin(), m_RhoList.cend());
	else RhoSet = std::set<double>(m_RhoList.crbegin(), m_RhoList.crend());

	// Get iterator to lower bound of RhoSet for RhoValue.
	// That is, the value in RhoSet that is closest to but less than RhoValue
	auto lBound = RhoSet.lower_bound(RhoValue);

	bool IsFound = (lBound != RhoSet.end());

	Ind = -1;
	if (lBound != RhoSet.end()) {
		// Get index of lower bound in sorted rho list
		Ind = std::distance(RhoSet.begin(), lBound) - 1;

		if (!RhoIsAscending)
			Ind = m_RhoList.size() - 2 - Ind;

		int ind2 = Ind + 1;

		if (ind2 >= GetCount()){
			ind2 = GetCount() - 1;
			Ind = ind2 - 1;
		}

		Weight = (RhoValue - m_RhoList[Ind]) / (m_RhoList[ind2] - m_RhoList[Ind]);


// 		if (!RhoIsAscending)
// 			ind2 = Ind - 1;

		if (ind2 < 0){
			ind2 = 0;
			Ind = 1;
		}

		OutVec = m_XYZList[Ind] + (m_XYZList[ind2] - m_XYZList[Ind]) * Weight;
	}

	return IsFound;
}

bool GradPathBase_c::GetDeviationMidpointRhoBasedFromOtherGP(GradPathBase_c const & GP, double const & CheckAngle, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2) const{
	REQUIRE(IsMade() && m_XYZList.size() > 3);
	REQUIRE(GP.IsMade() && GP.GetCount() > 3);

	double Weight;

	vec3 v1, v2;

	for (Ind1 = 1; Ind1 < GetCount() - 1; ++Ind1) {
		if (GP.GetPointAtRhoValue(m_RhoList[Ind1], Pt2, Ind2, Weight)) {
			Ind2 = MIN(Ind2, GP.GetCount() - 2);
			v1 = m_XYZList[Ind1 + 1] - m_XYZList[Ind1];
			v2 = GP.XYZAt(Ind2 + 1) - GP.XYZAt(Ind2);

			if (VectorAngle(v1, v2) > CheckAngle) {
				OutPt = (m_XYZList[Ind1] + Pt2) * 0.5;
				return true;
			}
		}
		else
			break;
	}

	return false;
}

/*
 *	Function to guarantee that the midpoint between two GPs that deviate occurs between the GPs
 *	(not outside of the inter-GP region).
 *	Do this by checking the normals of the triangles formed by the deviating segments of each GP
 *	and comparing them to the normals of the triangles formed by the deviating segment of each GP
 *	with the segment from the deviating point of each GP with the current midpoint.
 *	Then do essentially a binary search until the midpoint is in the inter-GP region.
 */
void GradPathBase_c::FixDeviationMidpoint(GradPathBase_c const & GP, vec3 & MidPt, int Ind1, int Ind2) const
{
	// First get the normals for the deviating segments of the paths
	if (Ind1 < 0){
		Ind1++;
	}
	else if (Ind1 == this->GetCount()){
		Ind1--;
	}

	if (Ind2 < 0) {
		Ind2++;
	}
	else if (Ind2 == GP.GetCount()) {
		Ind2--;
	}

	int Ind1m = Ind1, Ind1p = Ind1,
		Ind2m = Ind2, Ind2p = Ind2;

	int NumBufferPoints = 5;

	for (int i = 0; i < NumBufferPoints; ++i){
		Ind1m = MAX(Ind1m - 1, 0);
		Ind2m = MAX(Ind2m - 1, 0);
		Ind1p = MIN(Ind1p + 1, this->GetCount() - 1);
		Ind2p = MIN(Ind2p + 1, GP.GetCount() - 1);
	}

	vec3 p1m = this->XYZAt(Ind1m), p1 = this->XYZAt(Ind1), p1p = this->XYZAt(Ind1p),
		p2m = GP.XYZAt(Ind2m), p2 = GP.XYZAt(Ind2), p2p = GP.XYZAt(Ind2p);

	vec3 v11 = p1 - p1m, v12 = p1p - p1m,
		v21 = p2 - p2m,	v22 =p2p - p2m;

	vec3 n1 = cross(v11, v12),
		n2 = cross(v21, v22);

	bool DoLoop;

	do 
	{
		DoLoop = false;
		if (dot(n1, cross(v11, MidPt - p1)) > 0){
			MidPt = (MidPt + p2) * 0.5;
			DoLoop = true;
		}
		else if (dot(n2, cross(v21, MidPt - p2)) > 0){
			MidPt = (MidPt + p1) * 0.5;
			DoLoop = true;
		}
	} while (DoLoop);
}

bool GradPathBase_c::GetDeviationMidpointClosestPointBasedFromOtherGP(GradPathBase_c const & GP, double const & CheckAngle, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2) const {
	REQUIRE(IsMade() && m_XYZList.size() > 3);
	REQUIRE(GP.IsMade() && GP.GetCount() > 3);

	double Weight;

	vec3 v1, v2;

	for (Ind1 = 1; Ind1 < GetCount() - 1; ++Ind1) {
		Pt2 = GP.ClosestPoint(this->XYZAt(Ind1), Ind2);
		Ind2 = MIN(Ind2, GP.GetCount() - 2);
		v1 = m_XYZList[Ind1 + 1] - m_XYZList[Ind1];
		v2 = GP.XYZAt(Ind2 + 1) - GP.XYZAt(Ind2);

		if (VectorAngle(v1, v2) > CheckAngle) {
			Ind1--;
			Pt2 = GP.ClosestPoint(this->XYZAt(Ind1), Ind2);
			Ind2 = MIN(Ind2, GP.GetCount() - 2);
			OutPt = (m_XYZList[Ind1] + Pt2) * 0.5;
// 			this->FixDeviationMidpoint(GP, OutPt, Ind1, Ind2);
			return true;
		}
	}

	return false;
}

/*
 *	Check for deviation with another gradient path.
 *	Assumes that both paths start from the same point (i.e. at m_XYZList[0]).
 *	Checks for deviation according to four criteria:
 *	 1.	Distance between points of equal rho value as a factor of the distance between terminal points
 *	 2. Angle between segments (starting at points of equal rho value) as a factor of the angle between terminal segments
 *	 3. Angle between segment on either path with the closest point on the other
 *	 4. Angle between segments (starting at points of equal rho value) compared to a defined maximum allowed angle
 *	Deviation is met when:
 *		(1) and (2) are met simultaneously, so that the distance and angle are both greater
 *			than the check distance and angles,
 *		(3) is met, so that the angle between a path segment and the closest point on the other path
 *			is greater than the maximum allowed angle, or
 *		(4) is met, so that the angle is greater than the maximum allowed angle.
 *	Deviation is checked by first finding the starting point for the search on Path1 (this path)
 *	and then by stepping down the paths checking for deviation at each point.
 *	Each iteration, a point from either Path1 or Path2 is chosen based on whichever 
 *	path has its next point closest to its previously used point.
 *	This makes sure that we're using which ever path is most dense locally.
 */
bool GradPathBase_c::GetDeviationMidpointAsTerminalAngleAndDistanceFactorFromOtherGP(GradPathBase_c const & GP, 
	double const & CheckTermAngleFactor, 
	double const & CheckTermDistFactor, 
	double const & StartPtAsFactorOfLength, 
	vec3 & OutPt, int & Ind1, 
	int & Ind2, vec3 & Pt2,
	double const & MinCheckDist,
	double const & MinCheckAngle) const 
{
	REQUIRE(IsMade() && m_XYZList.size() > 3);
	REQUIRE(GP.IsMade() && GP.GetCount() > 3);

	double TermDist = Distance(m_XYZList.back(), GP.XYZAt(-1));

	if (TermDist > MinCheckDist) {
		// Terminal points are not (nearly) coincident, so assume they terminate at
		// distinct CPs or at system boundary/cutoff.
		// Now get the angle of terminal segments and distance of terminal points
		// and define the deviation point as the point where
		// equivalent segments on the two paths reach some factor of the
		// terminal angle and the distance of equivalent points is greater then
		// some factor of the terminal distance.
		double Weight, TermAngle, CheckAngle, CheckDistSqr;
		vec3 v1 = XYZAt(-1) - XYZAt(-2), 
			v2 = GP.XYZAt(-1) - GP.XYZAt(-2);

		TermAngle = VectorAngle(v1, v2);
		CheckAngle = MAX(TermAngle, MinCheckAngle) * CheckTermAngleFactor;

		CheckDistSqr = MAX(TermDist, MinCheckDist) * CheckTermDistFactor;
		double SmallCheckDistFactor = 1.;
		double SmallCheckDistSqr = CheckDistSqr * SmallCheckDistFactor;
		CheckDistSqr *= CheckDistSqr;
		SmallCheckDistSqr *= SmallCheckDistSqr;


		GradPathBase_c const *Path1 = this,
			*Path2 = &GP;

		double StartGPLength = Path1->GetLength(),
			GP2Length = Path2->GetLength();
		if (GP2Length < StartGPLength){
			auto tmpPath = Path2;
			Path2 = Path1;
			Path1 = tmpPath;
			StartGPLength = GP2Length;
		}

		double StartLength = StartGPLength * StartPtAsFactorOfLength,
			CheckLength = 0.0;
		int StartInd;
		for (StartInd = 0; StartInd < Path1->GetCount() - 1; ++StartInd){
			CheckLength += Distance(Path1->XYZAt(StartInd), Path1->XYZAt(StartInd + 1));
			if (CheckLength > StartLength)
				break;
		}

		Ind1 = StartInd;
		Ind2 = 0;
		vec3 Pt1;

		// Pt1 and Pt2 go with Path1 and Path2, though once deviation is found,
		// Pt2 will be the point on GP where the deviation occurs.
		do {
			// Perform deviation checks for current points.
			// First check based on the segment angle and distance at the
			// current point on Path1.
			// Don't allow Ind2 to go to below what it was previously, as this
			// can result in an infinite loop.
			int OldInd2 = Ind2;
			if (Path2->GetPointAtRhoValue(Path1->RhoAt(Ind1), Pt2, Ind2, Weight)) {
				Ind2 = MIN(MAX(OldInd2, Ind2), Path2->GetCount() - 2);
				v1 = Path1->XYZAt(Ind1 + 1) - Path1->XYZAt(Ind1);
				v2 = Path2->XYZAt(Ind2 + 1) - Path2->XYZAt(Ind2);

				double PointDistSqr = DistSqr(Path1->XYZAt(Ind1), Pt2);
				double VecAngle = VectorAngle(v1, v2);

				bool IsDeviating = false;
				bool HardDeviation = false;
				if (VecAngle > GP_DeviationAngleMaxCutoff
					|| (PointDistSqr > CheckDistSqr && VecAngle > CheckAngle))
				{
					IsDeviating = true;
					HardDeviation = VecAngle > GP_DeviationAngleMaxCutoff;
// 					OutPt = (Path1->XYZAt(Ind1) + Pt2) * 0.5;
				}

				if (!IsDeviating) {
					// Second check: Path1 segment against closest point on Path2
					int Ind3;
					vec3 Pt3 = Path2->ClosestPoint(Path1->XYZAt(Ind1), Ind3);
					vec3 v3 = Pt3 - Path1->XYZAt(Ind1);
					VecAngle = VectorAngle(v1, v3);
					PointDistSqr = DistSqr(Pt3, Path1->XYZAt(Ind1));
					if (VecAngle > GP_DeviationAngleMaxCutoff && PointDistSqr > SmallCheckDistSqr){
// 						Ind2 = Ind3;
						IsDeviating = true;
						HardDeviation = true;
// 						OutPt = (Path1->XYZAt(Ind1) + Pt2) * 0.5;
					}
				}
				if (!IsDeviating) {
					// Third check: Path2 segment against closest point on Path1
					int Ind3;
					vec3 Pt3 = Path1->ClosestPoint(Path2->XYZAt(Ind2), Ind3);
					vec3 v3 = Pt3 - Path2->XYZAt(Ind2);
					VecAngle = VectorAngle(v2, v3);
					PointDistSqr = DistSqr(Pt3, Path2->XYZAt(Ind2));
					if (VecAngle > GP_DeviationAngleMaxCutoff && PointDistSqr > SmallCheckDistSqr){
						IsDeviating = true;
						HardDeviation = true;
// 						OutPt = (Path2->XYZAt(Ind2) + Pt2) * 0.5;
// 						Ind1 = Ind3;
					}
				}

// 				if (HardDeviation && Ind1 == StartInd && StartInd > 1)
// 				{
// 					StartInd--;
// 					Ind1 = StartInd;
// 					Ind2 = 0;
// 					continue;
// 				}

				if (IsDeviating) {
  					if (HardDeviation){
 						Path1->GetDeviationMidpointClosestPointBasedFromOtherGP(*Path2, GP_DeviationAngleMaxCutoff * 0.15, OutPt, Ind1, Ind2, Pt2);
  					}
 					else {
//  						double Len1 = Path1->GetLength(Ind1),
//  							Len2 = Path2->GetLength(Ind2),
//  							LenSum = Len1 + Len2,
//  							Len1Ratio = Len1 / LenSum;
// 						 					OutPt = Path1->XYZAt(Ind1) * (1.0 - Len1Ratio) + Pt2 * Len1Ratio;
 						
// 						Ind1--;
// 						Pt2 = Path2->ClosestPoint(Path1->XYZAt(Ind1), Ind2);
						OutPt = (Path1->XYZAt(Ind1) + Pt2) * 0.5;

 					}

					if (Path1 != this){
						int tmpInd = Ind2;
						Ind2 = Ind1;
						Ind1 = tmpInd;
						Pt2 = GP.XYZAt(Ind2);
					}

					return true;
				}

				if (dot(v1,v1) > dot(v2,v2)){
					// Path2 has the shorter length than Path1, so switch it to Path1 for the next iteration.
					auto tmpPath = Path2;
					Path2 = Path1;
					Path1 = tmpPath;

					int tmpInd = Ind2;
					Ind2 = Ind1;
					Ind1 = tmpInd;
				}
				Ind1++;
			}
			else
				break;

			// Set Path1 and Path2 based on distance
		} while (Ind1 < Path1->GetCount() - 2 && Ind2 < Path2->GetCount() - 2);
	}

	return false;
}

// Assumes only one max distance between paths, so quits after finding one.
struct MaxSeparationMidpointData_s {
	double DistSqr = -1;
	int Ind1;
	int Ind2;
	vec3 Pt;
};
bool GradPathBase_c::GetMaxSeparationMidpointFromOtherGP(GradPathBase_c const & GP, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2, double & MaxDist, int step) const{
	REQUIRE(IsMade() && m_XYZList.size() > 3);
	REQUIRE(GP.IsMade() && GP.GetCount() > 3);
// 	REQUIRE(DistSqr(m_XYZList[0], GP.XYZAt(0)) < 0.001);
// 	REQUIRE(DistSqr(m_XYZList.back(), GP.XYZAt(-1)) < 0.001);


// 	std::deque<double> SqrDistVals;

	double Weight;
// 	vec3 Pt2New;
	int StartInd = (double)GetCount() * 0.05,
		EndInd = (double)(GetCount()-1) * 0.95;

	vector<MaxSeparationMidpointData_s> DataVec(this->GetCount());

#pragma omp parallel for
	for (int i = StartInd; i < EndInd; i+=step){
		DataVec[i].Ind1 = i;
		DataVec[i].Pt = GP.ClosestPoint(m_XYZList[i], DataVec[i].Ind2);
		DataVec[i].DistSqr = DistSqr(m_XYZList[i], DataVec[i].Pt);
	}

	auto maxDistSqr = std::max_element(DataVec.begin(), DataVec.end(), [](MaxSeparationMidpointData_s const & e1, MaxSeparationMidpointData_s const & e2) {return e1.DistSqr < e2.DistSqr; });
	Ind1 = maxDistSqr->Ind1;
	Ind2 = maxDistSqr->Ind2;
	Pt2 = maxDistSqr->Pt;
// 		Pt2 = (this->XYZAt(Ind1) + GP.XYZAt(Ind2)) * 0.5;
	OutPt = (m_XYZList[maxDistSqr->Ind1] + Pt2) * 0.5;
	MaxDist = sqrt(maxDistSqr->DistSqr);
	return true;
}

bool GradPathBase_c::GetMaxSeparationMidpointFromOtherGPRhoBased(GradPathBase_c const & GP, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2, double & MaxDist, int step) const {
	REQUIRE(IsMade() && m_XYZList.size() > 3);
	REQUIRE(GP.IsMade() && GP.GetCount() > 3);
	// 	REQUIRE(DistSqr(m_XYZList[0], GP.XYZAt(0)) < 0.001);
	// 	REQUIRE(DistSqr(m_XYZList.back(), GP.XYZAt(-1)) < 0.001);


	// 	std::deque<double> SqrDistVals;

	double Weight;
	// 	vec3 Pt2New;
	int StartInd = (double)GetCount() * 0.05,
		EndInd = (double)(GetCount() - 1) * 0.95;

	vector<MaxSeparationMidpointData_s> DataVec(this->GetCount());

#pragma omp parallel for
	for (int i = StartInd; i < EndInd; i+=step) {
		DataVec[i].Ind1 = i;
		if (GP.GetPointAtRhoValue(RhoAt(i), DataVec[i].Pt, DataVec[i].Ind2, Weight))
			DataVec[i].DistSqr = DistSqr(m_XYZList[i], DataVec[i].Pt);
		else {
			DataVec[i].Pt = GP.ClosestPoint(m_XYZList[i], DataVec[i].Ind2);
			DataVec[i].DistSqr = DistSqr(m_XYZList[i], DataVec[i].Pt);
		}
	}

	auto maxDistSqr = std::max_element(DataVec.begin(), DataVec.end(), [](MaxSeparationMidpointData_s const & e1, MaxSeparationMidpointData_s const & e2) {return e1.DistSqr < e2.DistSqr; });
	Ind1 = maxDistSqr->Ind1;
	Ind2 = maxDistSqr->Ind2;
	Pt2 = maxDistSqr->Pt;
	// 		Pt2 = (this->XYZAt(Ind1) + GP.XYZAt(Ind2)) * 0.5;
	OutPt = (m_XYZList[maxDistSqr->Ind1] + Pt2) * 0.5;
	MaxDist = sqrt(maxDistSqr->DistSqr);
	return true;
}

// Assumes only one max distance between paths, so quits after finding one.
bool GradPathBase_c::GetSeparationMidpointAtDistFromOtherGP(GradPathBase_c const & GP, vec3 & OutPt, int & Ind1, int & Ind2, vec3 & Pt2, const double & MaxDist, const double & StartPtAsFactorOfLength) const {
	REQUIRE(IsMade() && m_XYZList.size() > 3);
	REQUIRE(GP.IsMade() && GP.GetCount() > 3);
	double distSqr = DistSqr(m_XYZList[0], GP.XYZAt(0));
	REQUIRE(distSqr < MaxDist * MaxDist);

	double MaxDistSqr = MaxDist * MaxDist;
	double Weight;

	GradPathBase_c const *Path1 = this,
		*Path2 = &GP;

	double StartGPLength = Path1->GetLength(),
		GP2Length = Path2->GetLength();
	if (GP2Length < StartGPLength) {
		auto tmpPath = Path2;
		Path2 = Path1;
		Path1 = tmpPath;
		StartGPLength = GP2Length;
	}

	double StartLength = StartGPLength * StartPtAsFactorOfLength,
		CheckLength = 0.0;
	int StartInd;
	for (StartInd = 0; StartInd < Path1->GetCount() - 1; ++StartInd) {
		CheckLength += Distance(Path1->XYZAt(StartInd), Path1->XYZAt(StartInd + 1));
		if (CheckLength > StartLength)
			break;
	}

	Ind1 = StartInd;
	Ind2 = 0;
	vec3 Pt1;

	bool DistanceFromStartPointReached = false;

	// Pt1 and Pt2 go with Path1 and Path2, though once deviation is found,
	// Pt2 will be the point on GP where the deviation occurs.
	do {
		// Perform deviation checks for current points.
		// First check based on the segment angle and distance at the
		// current point on Path1.
		// Don't allow Ind2 to go to below what it was previously, as this
		// can result in an infinite loop.
		int OldInd2 = Ind2;
		Path2->GetPointAtRhoValue(Path1->RhoAt(Ind1), Pt2, Ind2, Weight);
		Ind2 = CLAMP(Ind2, OldInd2, Path2->GetCount() - 2);

		double PointDistSqr = DistSqr(Path1->XYZAt(Ind1), Pt2);

		if (!DistanceFromStartPointReached){
			DistanceFromStartPointReached = VectorAngle(Path1->XYZAt(Ind1) - Path1->XYZAt(0), Path2->XYZAt(Ind1) - Path2->XYZAt(0)) < 3.0 * PIOVER2;
		}

		if (DistanceFromStartPointReached && PointDistSqr > MaxDistSqr){
			OutPt = (Path1->XYZAt(Ind1) + Pt2) * 0.5;
			if (Path1 != this) {
				int tmpInd = Ind2;
				Ind2 = Ind1;
				Ind1 = tmpInd;
				Pt2 = GP.XYZAt(Ind2);
			}
			return true;
		}

		vec3 v1 = Path1->XYZAt(Ind1 + 1) - Path1->XYZAt(Ind1);
		vec3 v2 = Path2->XYZAt(Ind2 + 1) - Path2->XYZAt(Ind2);

		if (dot(v1, v1) > dot(v2, v2)) {
			// Path2 has the shorter length than Path1, so switch it to Path1 for the next iteration.
			auto tmpPath = Path2;
			Path2 = Path1;
			Path1 = tmpPath;

			int tmpInd = Ind2;
			Ind2 = Ind1;
			Ind1 = tmpInd;
		}
		Ind1++;

		// Set Path1 and Path2 based on distance
	} while (Ind1 < Path1->GetCount() - 2 && Ind2 < Path2->GetCount() - 2);

	return false;

// 	double Weight;
// 	double MaxDistSqr = MaxDist * MaxDist;
// 	// 	vec3 Pt2New;
// 
// 	for (Ind1 = 1; Ind1 < GetCount() - 1; ++Ind1) {
// 		if (GP.GetPointAtRhoValue(m_RhoList[Ind1], Pt2, Ind2, Weight)) {
// 			double tmpDistSqr = DistSqr(m_XYZList[Ind1], Pt2);
// 			if (tmpDistSqr > MaxDistSqr) {
// 				OutPt = (m_XYZList[Ind1] + Pt2) * 0.5;
// 				return true;
// 			}
// 		}
// 		else
// 			break;
// 	}
// 
// 	return false;
}

using DistIndPair = std::pair<double, vector<int> >;
 
double GradPathBase_c::GetMaxSeparationFromNeighboringGPs(vector<GradPathBase_c const *> GPs, vector<int> & IndList, int & Ind) const
{
	REQUIRE(this->IsMade());
	for (auto i : GPs)
		REQUIRE(i->IsMade());

	int NumPoints = this->GetCount();
	vector<DistIndPair> DistIndList(NumPoints, std::make_pair(0.0, vector<int>(GPs.size()+1)));

#pragma omp parallel for
	for (int i = 0; i < NumPoints; ++i){
		DistIndList[i].second.back() = i;
		for (int igp = 0; igp < GPs.size(); ++igp) {
			auto gp = GPs[igp];
			vec3 OtherPt = gp->ClosestPoint(m_XYZList[i], DistIndList[i].second[igp]);
			DistIndList[i].first += Distance(m_XYZList[i], OtherPt);
		}
	}

	auto MaxDistInd = std::max_element(DistIndList.begin(), DistIndList.end(), [](DistIndPair const & e1, DistIndPair const & e2) {return e1.first < e2.first; });

	Ind = MaxDistInd->second.back();
	IndList = vector<int>(MaxDistInd->second.begin(), MaxDistInd->second.end()-1);
	return MaxDistInd->first;
}

double GradPathBase_c::RhoAt(int i) const { return m_RhoList[GetInd(i)]; }



/*
*	GradPath_c methods
*/

/*
*	Public Methods
*/

/*
*	Constructors and destructors
*/

/*
*	Default constructor
*/
GradPath_c::GradPath_c()
{
	m_StartPoint.fill(-1e50);
	m_ODE_Data.Direction = StreamDir_Invalid;
	m_HowTerminate = GPTerminate_Invalid;

	m_TermPoint.fill(-1e50);
	m_TermPointRadiusSqr = -1.0;

	m_NumCPs = -1;

	m_TermValue = -1.0;
}

/*
*	Constructor for grad path that will make itself
*/
GradPath_c::GradPath_c(vec3 const & StartPoint,
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
	FieldDataPointer_c const & RhoPtr)
{
	m_GradPathReady = SetupGradPath(StartPoint,
		Direction,
		NumGPPoints,
		HowTerminate,
		TermPoint,
		CPXYZPtrs,
		NumCPs,
		TermPointRadius,
		TermValue,
		MaxIJK,
		MaxXYZ,
		MinXYZ,
		GradPtrs,
		RhoPtr);
} //	GradPath_c::GradPath_c()


GradPath_c::GradPath_c(vec3 const & StartPoint,
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
	FESurface_c const * Surf)
{
	m_GradPathReady = SetupGradPath(StartPoint,
		Direction,
		NumGPPoints,
		GPType,
		HowTerminate,
		TermPoint,
		CPs,
		TermPointRadius,
		TermValue,
		VolInfo,
		HessPtrs,
		GradPtrs,
		RhoPtr,
		Surf);
} //	GradPath_c::GradPath_c()

GradPath_c::GradPath_c(EntIndex_t ZoneNum,
	vector<EntIndex_t> const & XYZRhoVarNums,
	AddOn_pa const & AddOnID) : GradPathBase_c(ZoneNum, XYZRhoVarNums, AddOnID){}

/*
	 *	Constructor for making a new GP down the midpoints of two neighboring
	 *	GPs. Specifically made for neighboring GPs in interatomic/ring surface
	 *	generation.
	 *	It is assumed that the GPs provided have the same origin or terminus that
	 *	will be used as the final point in the new path.
	 */
GradPath_c::GradPath_c(vector<GradPathBase_c const *> const & GPs,
	double const & Path1Weight,
	double const & RhoVal,
	int GPStep)
{
	REQUIRE(GPs.size() == 2);
	REQUIRE(Path1Weight >= 0 && Path1Weight <= 1);
	REQUIRE(RhoVal > 0);
	REQUIRE(GPStep == -1 || GPStep == 1);
	for (auto const & i : GPs) {
		REQUIRE(i->IsMade());
		REQUIRE(i->GetCount() > 3);
		REQUIRE((i->RhoAt(0) < RhoVal && RhoVal < i->RhoAt(-1))
			|| (i->RhoAt(-1) < RhoVal && RhoVal < i->RhoAt(0)));
	}

	int NumGPs = 2;

	// If LeftRight Factor is significantly closer to one path or the other,
	// then use the close path to define the point spacing of the new path.
	// Otherwise, who whichever path has closer spacing around the point
	// that corresponds to RhoVal.
	// Path1 is defined as the path whose spacing will be used.

	double Path2Weight = 1.0 - Path1Weight;
	int Path1, Path2;
	int Path2Ind;
	double Weight;
	vec3 Pt;

	int Inds[2];
	for (int i = 0; i < 2; ++i){
		GPs[i]->GetPointAtRhoValue(RhoVal, Pt, Inds[i], Weight);
	}
	

	if (GPStep < 0){
		if (Inds[0] > Inds[1]) {
			Path1 = 0;
			Path2 = 1;
		}
		else{
			Path1 = 1;
			Path2 = 0;
		}
	}
	else{
		if (GPs[0]->GetCount() - Inds[0] > GPs[1]->GetCount() - Inds[1]) {
			Path1 = 0;
			Path2 = 1;
		}
		else {
			Path1 = 1;
			Path2 = 0;
		}
	}

	// Start building the new path starting at the index of the closest
	// point to RhoVal on Path1
	int Path1Ind;

	if (GPs[Path1]->GetPointAtRhoValue(RhoVal, Pt, Path1Ind, Weight)) {
		m_RhoList.push_back(GPs[Path1]->RhoAt(Path1Ind));
		m_XYZList.push_back(GPs[Path1]->XYZAt(Path1Ind) * Path1Weight);
		GPs[Path2]->GetPointAtRhoValue(RhoVal, Pt, Path2Ind, Weight);
		m_XYZList.back() += Pt * Path2Weight;
	}


	// While loop terminates when we get to the first or last point of Path1,
	// or if we reach the end of Path2.
	while (Path1Ind > 0
		&& Path1Ind < GPs[Path1]->GetCount() - 1)
	{
		Path1Ind += GPStep;
		double Path1RhoVal = GPs[Path1]->RhoAt(Path1Ind);
		if (GPs[Path2]->GetPointAtRhoValue(Path1RhoVal, Pt, Path2Ind, Weight)) {
			if (Path2Ind < 0 || Path2Ind >= GPs[Path2]->GetCount())
				break;
			m_RhoList.push_back(Path1RhoVal);
			m_XYZList.push_back(GPs[Path1]->XYZAt(Path1Ind) * Path1Weight + Pt * Path2Weight);
		}
		else
			break;
	}

	// Add the final point
	m_XYZList.push_back(GPs[Path1]->XYZAt(Path1Ind));
	m_RhoList.push_back(GPs[Path1]->RhoAt(Path1Ind));

	if (GPStep > 0) {
		SetStartEndCPNum(GPs[0]->GetStartEndCPNum(1), 1);
	}
	else {
		SetStartEndCPNum(GPs[0]->GetStartEndCPNum(0), 1);
	}

	m_GradPathMade = m_GradPathReady = TRUE;
}

// Constructor that makes a surface from the provided 
	// GradPaths and then seeds a new grad path constrained to that
	// surface.
GradPath_c::GradPath_c(vector<GradPath_c const *> const & GPs,
	vec3 const & SeedPoint,
	StreamDir_e GPDir,
	vector<std::pair<int, int> > StartEndGPInds,
	vec3 const * TermPoint,
	double const * TermPointRadius,
	GPTerminate_e GPTermType)
{
	REQUIRE(GPs.size() > 1);
	REQUIRE(GPDir == StreamDir_Forward || GPDir == StreamDir_Reverse);
	for (int i = 0; i < StartEndGPInds.size(); ++i)
		REQUIRE(StartEndGPInds[i].first >= 0 && abs(StartEndGPInds[i].second) < GPs[i]->GetCount());

	StartEndGPInds.resize(GPs.size(), std::make_pair(0, -1));

	vector<GradPath_c> TmpGPs(GPs.size());
	auto GPPtrs = GPs;

	for (int i = 0; i < GPs.size(); ++i){
		if (StartEndGPInds[i] != std::make_pair(0, -1)) {
			TmpGPs[i] = GPs[i]->SubGP(StartEndGPInds[i].first, StartEndGPInds[i].second);
			GPPtrs[i] = &TmpGPs[i];
		}
// 		else
// 			TmpGPs[i] = *GPs[i];
	}

	FESurface_c Surf;
	Surf.MakeFromGPs(GPPtrs);
	Surf.GeneratePointElementDistanceCheckData();

	//DEBUG
	bool doSaveSurf = false;
	if (doSaveSurf) {
		Surf.SaveAsTriFEZone({ 1,2,3 }, "InterGP Surface");
		for (int i = 0; i < GPs.size(); ++i){
			GradPath_c tmpGP = *GPs[i];
			tmpGP.SaveAsOrderedZone("GP " + to_string(i + 1));
			TecUtilZoneSetActive(Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_PlusEquals);
		}
	}

	*this = *GPs[0];

	this->Clear();
	this->SetSurfPtr(&Surf, 5);
	this->SetStartPoint(SeedPoint);
	this->SetDir(GPDir);
	if (TermPoint != nullptr) {
		if (GPTermType == GPTerminate_Invalid){
			GPTermType = GPTerminate_AtPoint;
		}
		this->SetTermPoint(*TermPoint);
		this->SetGPTermType(GPTermType);
		this->SetTermPointRadius(TermPointRadius != nullptr ? *TermPointRadius : 1e-4);
	}
	else{
		this->SetGPTermType(GPTerminate_AtRhoValue);
	}
	this->SetTerminalCPTypeNums({});
	m_StartEndCPNum[0] = -1;

 	this->Seed(false);

// 	this->ProjectPathToSurface();
	this->RemoveKinks();
	this->MakeRhoValuesMonotomic();

	if (TermPoint != nullptr && GPTermType == GPTerminate_AtPoint){
		double TermDistSqr = DistSqr(this->XYZAt(-1), *TermPoint);
		if (TermDistSqr > 1e-2){
			if (GPs.size() == 2){
				// GP failed to reach it's designated terminus, so complete it with a midpoint
				// path between the two provided GPs.
				vec3 Pts[2];
				int Inds[2];
				double Weight;
				Pts[0] = GPs[0]->ClosestPoint(this->XYZAt(-1), Inds[0]);
				Pts[1] = GPs[1]->ClosestPoint(this->XYZAt(-1), Inds[1]);
	// 			while (DistSqr(*TermPoint, GPs[0]->XYZAt(Inds[0])) > DistSqr(*TermPoint, this->XYZAt(-1))
	// 				&& --Inds[0] > 5) {}
	// 			GPs[1]->GetPointAtRhoValue(GPs[0]->RhoAt(Inds[0]), Pts[1], Inds[1], Weight);
				if (TermDistSqr > 0.5
					|| (MIN(Inds[0], GPs[0]->GetCount() - Inds[0]) > 2 && MIN(Inds[1], GPs[1]->GetCount() - Inds[1]) > 2)) {
					double PathDist = Distance(Pts[0], Pts[1]);
					if (PathDist < 1e-6)
						Weight = 0.5;
					else
						Weight = Distance(Pts[1], this->XYZAt(-1)) / PathDist;
					TmpGPs.resize(2);
					// 			TmpGPs[0] = GPs[0]->SubGP(0, Inds[0]);
					// 			TmpGPs[1] = GPs[1]->SubGP(0, Inds[1]);
					vector<GradPathBase_c const *> Ptrs({ GPs[0], GPs[1] });
					Weight = CLAMP(Weight, 0.0, 1.0);
					GradPath_c EndGP(Ptrs, Weight, GPs[0]->RhoAt(Inds[0]), -1);

					*this += EndGP;
				}
			}
		}
		if (sum(*TermPoint == this->XYZAt(-1)) != 3){
			this->PointAppend(*TermPoint, GPs[0]->RhoAt(0));
			this->SetStartEndCPNum(GPs[0]->GetStartEndCPNum(0), 1);
		}
// 		TecUtilDialogErrMsg("GP failed to reach termpoint");
	}

	return;
}

	// Construct from vector<vec3> points and optional vector<double> rho values
	GradPath_c::GradPath_c(vector<vec3> const & XYZList, vector<double> const RhoList){
		REQUIRE(XYZList.size() > 1);
		REQUIRE(RhoList.size() == XYZList.size() || RhoList.size() == 0);

		this->m_XYZList = XYZList;
		if (RhoList.size() > 0) {
			this->m_RhoList = RhoList;
		}
		else{
			this->m_RhoList = vector<double>(XYZList.size(), 0.0);
		}
		this->m_GradPathMade = true;
		this->m_NumGPPoints = XYZList.size();
	}

/*
*	Copy constructor
*/
GradPath_c::GradPath_c(GradPath_c const & a) : GradPathBase_c()
{
	*this = a;
}
GradPath_c::GradPath_c(GradPathBase_c const & a) : GradPathBase_c()
{
	*this = a;
}

GradPath_c::~GradPath_c()
{
}

/*
*	Operator declarations
*/
GradPath_c & GradPath_c::operator=(GradPath_c const & rhs)
{
	if (this == &rhs)
		return *this;

	GradPathBase_c::operator=(rhs);

	m_ODE_Data = rhs.m_ODE_Data;

	m_TermPoint = rhs.m_TermPoint;
	m_CPXYZPtrs = rhs.m_CPXYZPtrs;
	m_NumCPs = rhs.m_NumCPs;
	m_CPs = rhs.m_CPs;
	m_Surface = rhs.m_Surface;

	m_TermPointRadiusSqr = rhs.m_TermPointRadiusSqr;
	m_TermValue = rhs.m_TermValue;
	m_TerminalCPTypeNums = rhs.m_TerminalCPTypeNums;

	m_HowTerminate = rhs.m_HowTerminate;
	m_StartPoint = rhs.m_StartPoint;

	return *this;
}// GradPath_c & GradPath_c::operator=(GradPath_c const & rhs)
GradPath_c & GradPath_c::operator=(GradPathBase_c const & rhs)
{
	if (this == &rhs)
		return *this;

	GradPathBase_c::operator=(rhs);

	return *this;
}// GradPath_c & GradPath_c::operator=(GradPath_c const & rhs)
Boolean_t GradPath_c::operator==(GradPath_c const & rhs) const
{
	Boolean_t AreSame = this->IsSame(rhs);

	if (AreSame){
		AreSame = (sum(m_TermPoint == rhs.m_TermPoint) == 3 &&
			m_CPXYZPtrs == rhs.m_CPXYZPtrs &&
			m_NumCPs == rhs.m_NumCPs &&
			m_CPs == rhs.m_CPs &&
			m_Surface == rhs.m_Surface &&
			m_TermPointRadiusSqr == rhs.m_TermPointRadiusSqr &&
			m_TermValue == rhs.m_TermValue &&

			m_HowTerminate == rhs.m_HowTerminate &&
			sum(m_StartPoint == rhs.m_StartPoint) == 3 &&

			m_TerminalCPTypeNums == rhs.m_TerminalCPTypeNums &&

			m_GradPathReady == rhs.m_GradPathReady &&
			m_GradPathMade == rhs.m_GradPathMade
			);
	}

	return AreSame;
}


/*
*	Getter methods
*/



/*
*	Setter methods
*/

/*
*	This is basically a constructor for a pre-constructed GP
*	to assign everything needed for it to make itself.
*/
Boolean_t GradPath_c::SetupGradPath(vec3 const & StartPoint,
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
	FieldDataPointer_c const & RhoPtr)
{
	m_XYZList.swap(vector<vec3>());
	m_RhoList.swap(vector<double>());

	m_NumGPPoints = NumGPPoints;

	m_StartPoint = StartPoint;
	m_ODE_Data.Direction = Direction;
	m_HowTerminate = HowTerminate;

	if (TermPoint != nullptr)
		m_TermPoint = *TermPoint;
	if (TermPointRadius != nullptr)
		m_TermPointRadiusSqr = (*TermPointRadius * *TermPointRadius);

	m_CPs = nullptr;

	m_CPXYZPtrs = CPXYZPtrs;

	if (NumCPs != nullptr)
		m_NumCPs = *NumCPs;

	if (TermValue != nullptr)
		m_TermValue = *TermValue;

	m_ODE_Data.VolZoneInfo.MaxIJK = MaxIJK;
	m_ODE_Data.VolZoneInfo.MaxXYZ = MaxXYZ;
	m_ODE_Data.VolZoneInfo.MinXYZ = MinXYZ;

	m_ODE_Data.GradPtrs = GradPtrs;

	m_ODE_Data.RhoPtr = RhoPtr;

	m_GradPathReady = (MaxIJK.size() == 3
		&& GradPtrs.size() == 3);

	if (m_GradPathReady){
		for (int i = 0; i < 3 && m_GradPathReady; ++i){
			m_GradPathReady = GradPtrs[i].IsReady();
		}
	}
	if (m_GradPathReady){
		m_GradPathReady = RhoPtr.IsReady();
	}

	if (m_GradPathReady){
		if (HowTerminate == GPTerminate_AtPoint || HowTerminate == GPTerminate_AtPointRadius)
			m_GradPathReady = (TermPoint != nullptr && TermPointRadius != nullptr);
		else if (HowTerminate == GPTerminate_AtRhoValue)
			m_GradPathReady = (TermValue != nullptr);
		else if (HowTerminate == GPTerminate_AtCP || m_HowTerminate == GPTerminate_AtCPRadius){
			m_GradPathReady = (CPXYZPtrs.size() == 3 && TermPointRadius != nullptr && NumCPs != nullptr);
			if (m_GradPathReady){
				for (int i = 0; i < 3 && m_GradPathReady; ++i){
					m_GradPathReady = CPXYZPtrs[i].IsReady();
				}
			}
		}
		int GPSize = GP_NumPointsBufferFactor *m_NumGPPoints;
		m_XYZList.reserve(GPSize);
		m_RhoList.reserve(GPSize);
	}

	m_StartEndCPNum[0] = m_StartEndCPNum[1] = -1;

	m_GradPathMade = FALSE;

	return m_GradPathReady;
}

Boolean_t GradPath_c::SetupGradPath(vec3 const & StartPoint,
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
	FESurface_c const * Surf)
{
	m_XYZList.swap(vector<vec3>());
	m_RhoList.swap(vector<double>());

	m_NumGPPoints = NumGPPoints;

	m_StartPoint = StartPoint;
	m_ODE_Data.Direction = Direction;
	m_HowTerminate = HowTerminate;
	m_GPType = GPType;

	if (TermPoint != nullptr)
		m_TermPoint = *TermPoint;
	if (TermPointRadius != nullptr)
		m_TermPointRadiusSqr = (*TermPointRadius * *TermPointRadius);

	m_NumCPs = -1;

	m_CPs = CPs;
	m_Surface = Surf;

	if (TermValue != nullptr)
		m_TermValue = *TermValue;

	m_ODE_Data.VolZoneInfo = VolInfo;

	m_ODE_Data.HessPtrs = HessPtrs;
	m_ODE_Data.HasHess = HessPtrs.size() == 6;
	for (int i = 0; i < HessPtrs.size() && m_ODE_Data.HasHess; ++i){
		m_ODE_Data.HasHess = HessPtrs[i].IsReady();
	}

	m_ODE_Data.GradPtrs = GradPtrs;
	m_ODE_Data.HasGrad = GradPtrs.size() == 3;
	for (int i = 0; i < GradPtrs.size() && m_ODE_Data.HasGrad; ++i){
		m_ODE_Data.HasGrad = GradPtrs[i].IsReady();
	}

	m_ODE_Data.RhoPtr = RhoPtr;

	m_GradPathReady = (m_ODE_Data.VolZoneInfo.MaxIJK.size() == 3);

	if (m_GradPathReady){
		m_GradPathReady = RhoPtr.IsReady();
	}

	if (m_GradPathReady){
		if (HowTerminate == GPTerminate_AtPoint || HowTerminate == GPTerminate_AtPointRadius)
			m_GradPathReady = (TermPoint != nullptr && TermPointRadius != nullptr);
		else if (HowTerminate == GPTerminate_AtRhoValue)
			m_GradPathReady = (TermValue != nullptr);
		else if (HowTerminate == GPTerminate_AtCP || m_HowTerminate == GPTerminate_AtCPRadius){
			m_GradPathReady = (m_CPs != nullptr && TermPointRadius != nullptr && m_CPs->NumCPs() > 0);
			if (m_GradPathReady) m_NumCPs = m_CPs->NumCPs();
		}
	}

	m_StartEndCPNum[0] = m_StartEndCPNum[1] = -1;

	m_GradPathMade = FALSE;

	if (m_GradPathReady){
		int GPSize = GP_NumPointsBufferFactor * m_NumGPPoints;
		m_XYZList.reserve(GPSize);
		m_RhoList.reserve(GPSize);
	}

	return m_GradPathReady;
}

Boolean_t GradPath_c::SeedInDirection(StreamDir_e const & Direction){
	Boolean_t IsOk = m_GradPathReady && !m_GradPathMade;

	m_ODE_Data.HasGrad = (m_ODE_Data.GradPtrs.size() == 3);
	m_ODE_Data.HasHess = (m_ODE_Data.HessPtrs.size() == 6);

	StreamDir_e OldDir = m_ODE_Data.Direction;
	m_ODE_Data.Direction = Direction;

	gsl_odeiv2_system ODESys = { &GP_ODE_Gradient, &GP_ODE_Jacobian, m_ODE_NumDims, &m_ODE_Data };

// 	gsl_odeiv2_step_type const * T = gsl_odeiv2_step_rk4;
	gsl_odeiv2_step_type const * T = gsl_odeiv2_step_bsimp;
// 	gsl_odeiv2_step_type const * T = gsl_odeiv2_step_msbdf;
	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc(T, m_ODE_NumDims);
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new(1e-10, 1e-11);
	gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc(m_ODE_NumDims);



	// 	gsl_odeiv2_driver * ODEDriver;
	// 	ODEDriver = gsl_odeiv2_driver_alloc_yp_new(&ODESys, gsl_odeiv2_step_rk2, 1e-3, 1e-2, 0);
	// 	gsl_odeiv2_driver_set_hmin(ODEDriver, 1e-8);
	// 	gsl_odeiv2_driver_set_hmax(ODEDriver, 1e-1);

	if (IsOk){
		double tInit = 0.0;
		double tFinal = DBL_MAX;

		double h = 1e-10;

		double y[3] = { m_StartPoint[0], m_StartPoint[1], m_StartPoint[2] };

		int SurfLastProjectedElem = -1;
		int SurfProjectionIter = 0;
		bool SurfProjectionFound;
		vector<vec3> SurfXYZList;
		vector<double> SurfRhoList;

		if (m_Surface != nullptr && m_Surface->IsMade() && m_SurfGPProjectionFrequency <= 1){
			if (!m_Surface->ProjectPointToSurface(y, m_StartPoint, SurfLastProjectedElem, SurfProjectionFound)){
				TecUtilDialogErrMsg("Projection of point to FESurface failed");
			}
			for (int i = 0; i < 3; ++i) y[i] = m_StartPoint[i];
		}

		vec3 NewPoint, PtIm1 = m_StartPoint, PtI;

		IsOk = SetIndexAndWeightsForPoint(m_StartPoint, m_ODE_Data.VolZoneInfo);

		if (IsOk){
			m_XYZList.push_back(m_StartPoint);
			m_RhoList.push_back(RhoByCurrentIndexAndWeights());
		}

		m_StartEndCPNum[1] = -1;

		vector<CPTypeNum_e> CheckCPTypes;
		if (m_TerminalCPTypeNums.empty())
			CheckCPTypes = CPNearFieldTypes;
		else
			CheckCPTypes = m_TerminalCPTypeNums;

		int Status = GSL_SUCCESS;
		int Step = 1;

		unsigned int NumStalledPoints = 0;

		MultiRootObjects_s MR;
		MultiRootParams_s Params;
		vec3 StepDir, TmpPt, EigVals, DotPdts;
		mat33 EigVecs, Hess;
		vector<vec3> BV(3);
		BV[0] << 1 << 0 << 0;
		BV[1] << 0 << 1 << 0;
		BV[2] << 0 << 0 << 1;

		double StepSize, PlaneCPDist;
		Boolean_t SGPFound = FALSE;
		int PlaneCPFailIter = 0, PlaneCPIter = 0;
		mat33 I = eye<mat>(3, 3);
		if (m_GPType != GPType_Classic && m_GPType != GPType_Invalid){
			Params.CalcType = m_GPType;
			Params.HasHess = m_ODE_Data.HessPtrs.size() == 6;
			Params.HessPtrs = &m_ODE_Data.HessPtrs;
			for (int i = 0; i < 6 && Params.HasHess; ++i)
				Params.HasHess = Params.HessPtrs->at(i).IsReady();
			Params.IsPeriodic = m_ODE_Data.VolZoneInfo.IsPeriodic;
			Params.VolInfo = &m_ODE_Data.VolZoneInfo;
			Params.RhoPtr = &m_ODE_Data.RhoPtr;
			Params.HasGrad = m_ODE_Data.GradPtrs.size() == 3;
			Params.GradPtrs = &m_ODE_Data.GradPtrs;
			for (int i = 0; i < 3 && Params.HasGrad; ++i)
				Params.HasGrad = Params.GradPtrs->at(i).IsReady();
			Params.BasisVectors = &I;
			Params.VolInfo->BasisVectors = I;
			Params.Origin = &BV[0];

			MR.Func = { &F2DGrad, &DF2DGrad, &FDF2DGrad, 2, &Params };
			MR.pos = gsl_vector_alloc(2);
			MR.T = gsl_multiroot_fdfsolver_hybridsj;
			MR.s = gsl_multiroot_fdfsolver_alloc(MR.T, 2);
		}

		while (IsOk && Status == GSL_SUCCESS && Step < GP_MaxNumPoints){


			// 			Status = gsl_odeiv2_driver_apply(ODEDriver, &tInit, tInit + 1e-3, y);
			Status = gsl_odeiv2_evolve_apply(e, c, s, &ODESys, &tInit, tFinal, &h, y);

			h = MIN(h, GP_MaxStepSize);

			if (Status == GSL_SUCCESS || Status == GSL_EDOM){
				if (m_Surface != nullptr && m_Surface->IsMade()) {
					SurfProjectionIter++;
					SurfProjectionFound = false;
					if (NumStalledPoints < 2 && (m_SurfGPProjectionFrequency < 0 || SurfProjectionIter % m_SurfGPProjectionFrequency == 0)) {
						SurfProjectionIter = 0;
						if (!m_Surface->ProjectPointToSurface(y, PtI, SurfLastProjectedElem, SurfProjectionFound, 10, (m_SurfGPProjectionFrequency <= 0))) {
							TecUtilDialogErrMsg("Projection of point to FESurface failed");
						}
						for (int i = 0; i < 3; ++i) y[i] = PtI[i];
					}
				}
				
				PtI = y;

				if (m_GPType != GPType_Classic && m_GPType != GPType_Invalid){
					/*
					*	This allows a mixing of the last direction, as determined by the found 2d CP, and
					*	the current direction, as determined by the gradient at the last GP point.
					*/
					//   					if (!SGPFound){
					StepDir = normalise((PtI - PtIm1));

					/*
					*	To orient the plane used for the 2d CP search,
					*	need the eigenvector whose dot product with the
					*	gradient is farthest from zero.
					*/

					CalcEigenvecDotGradForPoint(PtI, DotPdts, EigVals, EigVecs, FALSE, Params);
					int MaxDir = 0;
					double MaxVal = 0.0;
					for (int i = 0; i < 3; ++i){
						if (abs(DotPdts[i]) > MaxVal){
							MaxVal = abs(DotPdts[i]);
							MaxDir = i;
						}
					}

					StepDir = normalise((EigVecs.row(MaxDir)));

					if (Step <= 1) StepSize = Distance(PtI, PtIm1);
					//   					}
					if (m_GPType != GPType_Classic && m_GPType != GPType_Invalid && CPInNormalPlane(PtI, StepDir, MR)){
						Params.BasisVectors = &I;
						PlaneCPIter = 0;
						SGPFound = FALSE;
						do
						{
							PlaneCPIter++;
							TmpPt = PtI;

							CalcEigenvecDotGradForPoint(PtI, DotPdts, EigVals, EigVecs, FALSE, Params);
							MaxDir = 0;
							MaxVal = 0.0;
							for (int i = 0; i < 3; ++i){
								if (abs(DotPdts[i]) > MaxVal){
									MaxVal = abs(DotPdts[i]);
									MaxDir = i;
								}
							}

							StepDir = normalise((EigVecs.row(MaxDir)));

							SGPFound = CPInNormalPlane(PtI, StepDir, MR);
							Params.BasisVectors = &I;
							PlaneCPDist = Distance(PtI, TmpPt);
						} while (SGPFound && PlaneCPDist >= 1e-4 && PlaneCPIter < 50);
						if (PlaneCPIter < 50 && PlaneCPDist < 1e-4){
							/*
							*	This allows a mixing of the last direction, as determined by the found 2d CP, and
							*	the current direction, as determined by the gradient at the last GP point.
							*/
							// 						StepDir = normalise((StepDir) * m_DirMixFactor + normalise( (PtI - PtIm1)) * normalise( (1.0 - m_DirMixFactor)));
							PtI = PtIm1 + normalise((PtI - PtIm1)) * StepSize;
							for (int i = 0; i < 3; ++i)
								y[i] = PtI[i];
						}
					}
					// 					if (!SGPFound){
					else{
						Params.BasisVectors = &I;
						PlaneCPFailIter++;
						if (PlaneCPFailIter > GP_PlaneCPMaxIter){
							Status = GSL_EFAILED;
							break;
						}
					}
				}

				double Rho = RhoByCurrentIndexAndWeights();

				if (m_HowTerminate == GPTerminate_AtRhoValue && Rho < m_TermValue){
					double OldRho;
					OldRho = m_RhoList[m_RhoList.size() - 1];

					NewPoint = PtIm1 + (PtI - PtIm1) * ((m_TermValue - OldRho) / (Rho - OldRho));

					m_XYZList.push_back(NewPoint);
					m_RhoList.push_back(m_TermValue);

					break;
				}
				else if (m_HowTerminate == GPTerminate_AtPoint || m_HowTerminate == GPTerminate_AtPointRadius){
					double PointRadiusSqr = DistSqr(PtI, (m_TermPoint));
					if (PointRadiusSqr <= m_TermPointRadiusSqr){
						if (m_HowTerminate == GPTerminate_AtPointRadius){
							double OldRadius = Distance(PtIm1, m_TermPoint);

							NewPoint = PtIm1 + (PtI - PtIm1) * ((sqrt(m_TermPointRadiusSqr) - OldRadius) / (sqrt(PointRadiusSqr) - OldRadius));
						}
						else
							NewPoint = m_TermPoint;

						IsOk = SetIndexAndWeightsForPoint(NewPoint, m_ODE_Data.VolZoneInfo);
						if (IsOk){
							Rho = RhoByCurrentIndexAndWeights();

							m_XYZList.push_back(NewPoint);
							m_RhoList.push_back(Rho);
						}

						break;
					}
				}
				else if (m_HowTerminate == GPTerminate_AtCP || m_HowTerminate == GPTerminate_AtCPRadius){
					Boolean_t PointFound = FALSE;
					if (m_CPs == nullptr){
						for (int CPNum = 0; CPNum < m_NumCPs && !PointFound; ++CPNum){
							if (CPNum != m_StartEndCPNum[0]){
								if (m_CPXYZPtrs.size() == 3) {
									for (int i = 0; i < 3; ++i) {
										NewPoint[i] = m_CPXYZPtrs[i][CPNum];
									}
								}
								else{
									NewPoint = m_CPs->GetXYZ(CPNum);
								}
								double PointRadiusSqr = DistSqr(PtI, NewPoint);
								if (PointRadiusSqr <= m_TermPointRadiusSqr){
									if (m_HowTerminate == GPTerminate_AtCPRadius){
										double OldRadius = Distance(PtIm1, NewPoint);

										NewPoint = PtIm1 + (PtI - PtIm1) * ((sqrt(m_TermPointRadiusSqr) - OldRadius) / (sqrt(PointRadiusSqr) - OldRadius));
									}

									IsOk = SetIndexAndWeightsForPoint(NewPoint, m_ODE_Data.VolZoneInfo);
									if (IsOk){
										Rho = RhoByCurrentIndexAndWeights();

										m_XYZList.push_back(NewPoint);
										m_RhoList.push_back(Rho);

										m_StartEndCPNum[1] = CPNum;
									}

									PointFound = TRUE;
								}
							}
						}
					}
					else{
// 						for (int TotCPNum = 0; TotCPNum < m_CPs->NumCPs() && !PointFound; ++TotCPNum){
// 							if (TotCPNum != m_StartEndCPNum[0]){
// 								NewPoint = m_CPs->GetXYZ(TotCPNum);
// 								double PointRadiusSqr = DistSqr(PtI, NewPoint);
// 								if (PointRadiusSqr <= m_TermPointRadiusSqr){
// 									if (m_HowTerminate == GPTerminate_AtCPRadius){
// 										double OldRadius = Distance(PtIm1, NewPoint);
// 
// 										NewPoint = PtIm1 + (PtI - PtIm1) * ((sqrt(m_TermPointRadiusSqr) - OldRadius) / (sqrt(PointRadiusSqr) - OldRadius));
// 									}
// 
// 									IsOk = SetIndexAndWeightsForPoint(NewPoint, m_ODE_Data.VolZoneInfo);
// 									if (IsOk){
// 										Rho = RhoByCurrentIndexAndWeights();
// 
// 										m_XYZList.push_back(NewPoint);
// 										m_RhoList.push_back(Rho);
// 
// 										m_StartEndCPNum[1] = TotCPNum;
// 									}
// 
// 									PointFound = TRUE;
// 								}
// 							}
// 						}
// 						
						for (auto CPTypeNum : CheckCPTypes) {
							double CheckRadiusSqr = m_TermPointRadiusSqr;
// 							if (m_HowTerminate != GPTerminate_AtCPRadius && CPTypeNum == CPTypeNum_Nuclear) {
// 								CheckRadiusSqr = 4 * sqrt(m_TermPointRadiusSqr);
// 								CheckRadiusSqr *= CheckRadiusSqr;
// 							}
							for (int CPNum = 0; CPNum < m_CPs->NumCPs(CPTypeNum) && !PointFound; ++CPNum){
								int TotCPNum = m_CPs->GetTotOffsetFromTypeNumOffset(CPTypeNum, CPNum);
								if (TotCPNum != m_StartEndCPNum[0]){
									NewPoint = m_CPs->GetXYZ(CPTypeNum, CPNum);
									double PointRadiusSqr = DistSqr(PtI, NewPoint);
									if (PointRadiusSqr <= CheckRadiusSqr){
										if (m_HowTerminate == GPTerminate_AtCPRadius){
											double OldRadius = Distance(PtIm1, NewPoint);
											NewPoint = PtIm1 + (PtI - PtIm1) * ((sqrt(m_TermPointRadiusSqr) - OldRadius) / (sqrt(PointRadiusSqr) - OldRadius));
										}

										IsOk = SetIndexAndWeightsForPoint(NewPoint, m_ODE_Data.VolZoneInfo);
										if (IsOk){
											Rho = RhoByCurrentIndexAndWeights();

											m_XYZList.push_back(NewPoint);
											m_RhoList.push_back(Rho);

											m_StartEndCPNum[1] = TotCPNum;
										}

										PointFound = TRUE;
									}
								}
							}
						}
					}
					if (PointFound)
						break;
				}
				if (m_TermValue > 0.0 && Rho < m_TermValue){
					double OldRho;
					OldRho = m_RhoList.back();

					NewPoint = PtIm1 + (PtI - PtIm1) * ((m_TermValue - OldRho) / (Rho - OldRho));

					m_XYZList.push_back(NewPoint);
					m_RhoList.push_back(m_TermValue);

					break;
				}

				if (IsOk){// && m_HowTerminate != GPTerminate_AtPoint && m_HowTerminate != GPTerminate_AtPointRadius){
					if (GP_StallNumPointsToCheck > 0 && m_XYZList.size() >= GP_StallNumPointsToCheck - 1){
						double NewPtDistSqr = Distance(PtI, m_XYZList[m_XYZList.size() - GP_StallNumPointsToCheck - NumStalledPoints + 1]);
						if (NewPtDistSqr < MIN(GP_StallPointDistTol, h) * GP_StallNumPointsToCheck){
							NumStalledPoints++;
							if (NumStalledPoints >= GP_StallPointCount && Status == GSL_SUCCESS){
								Status = GSL_ENOPROG;
							}
						}
						else
							NumStalledPoints = 0;
					}

					m_XYZList.push_back(PtI);
					m_RhoList.push_back(Rho);

					if (m_Surface != nullptr && SurfProjectionFound) {
						SurfXYZList.push_back(PtI);
						SurfRhoList.push_back(Rho);
					}

					PtIm1 = PtI;
				}
			}

			Step++;
		}

		IsOk = (IsOk && Status == GSL_SUCCESS || Status == GSL_EDOM || Status == GSL_ENOPROG);

 		if (Status == GSL_ENOPROG && m_XYZList.size() > GP_StallPointCount){
 		 	/*
 		 	*	Grad path stalled, so it was basically bouncing around the same point.
 		 	*	The last point can then be approximated as the midpoint of the stalled points.
 		 	*/
//  		 	vec3 Pt = zeros(3);
//  		 	double PtRho = 0.0;
//  		 	double PtCount = 0.0;
//  		 	for (int i = MAX(0, m_XYZList.size() - GP_StallPointCount); i < m_XYZList.size(); ++i){
//  		 		Pt += m_XYZList[i];
//  		 		PtRho += m_RhoList[i];
//  		 		PtCount += 1.0;
//  		 	}
//  		 
//  		 	Pt /= PtCount;
//  		 	PtRho /= PtCount;
//  		 
			vec3 Pt = m_XYZList.back();
			double PtRho = m_RhoList.back();
 		 
 		 	m_XYZList.resize(MAX(0, m_XYZList.size() - GP_StallPointCount));
 		 	m_XYZList.back() = Pt;
 		 	m_RhoList.resize(MAX(0, m_RhoList.size() - GP_StallPointCount));
 		 	m_RhoList.back() = PtRho;
//  		 	PtI = Pt;
 		}

// 		if (m_Surface != nullptr){
// 			m_XYZList = SurfXYZList;
// 			m_RhoList = SurfRhoList;
// 		}

		// 		gsl_odeiv2_driver_free(ODEDriver);

		gsl_odeiv2_evolve_free(e);
		gsl_odeiv2_control_free(c);
		gsl_odeiv2_step_free(s);

		if (m_GPType && MR.s != nullptr)
			gsl_multiroot_fdfsolver_free(MR.s);
		// 	if (m_GPType && MR.pos != nullptr)
		// 		gsl_vector_free(MR.pos);

		if (m_StartEndCPNum[1] < 0 && m_CPs != nullptr){
			/*
			*	Check to see if terminating point coincides with a CP
			*/
			Boolean_t PointFound = FALSE;
			for (auto CPTypeNum : CPNearFieldTypes) {
				for (int CPi = 0; CPi < m_CPs->NumCPs(CPTypeNum) && !PointFound; ++CPi) {
					int CPNum = m_CPs->GetTotOffsetFromTypeNumOffset(CPTypeNum, CPi);
					if (CPNum != m_StartEndCPNum[0]) {
						if (m_CPs == nullptr) {
							if (m_CPXYZPtrs.size() == 3) {
								for (int i = 0; i < 3; ++i) {
									NewPoint[i] = m_CPXYZPtrs[i][CPNum];
								}
							}
							else {
								NewPoint = m_CPs->GetXYZ(CPNum);
							}
						}
						else {
							NewPoint = m_CPs->GetXYZ(CPNum);
						}
						double PointRadiusSqr = DistSqr(PtI, NewPoint);
						if (PointRadiusSqr <= m_TermPointRadiusSqr) {
							m_StartEndCPNum[1] = CPNum;
							PointFound = TRUE;
						}
					}
				}
			}
		}

		m_GradPathMade = IsOk;
		if (IsOk && m_GPType)
			m_SGPMade = TRUE;
	}

	m_XYZList.shrink_to_fit();
	m_RhoList.shrink_to_fit();
	m_ODE_Data.Direction = OldDir;

	if (m_Surface != nullptr && m_Surface->IsMade() && m_SurfGPProjectionFrequency > 1) {
		this->ProjectPathToSurface(m_Surface);
	}

	return IsOk;
}

Boolean_t GradPath_c::Seed(bool const DoResample){
	Boolean_t IsOk = m_GradPathReady && !m_GradPathMade;


	if (m_ODE_Data.Direction != StreamDir_Both){
		IsOk = SeedInDirection(m_ODE_Data.Direction);
	}
	else{
		int StartCCPNum = m_StartEndCPNum[0];
		vector<GradPath_c> GPs(2, *this);
		vector<int> EndCPNums;
		m_GradPathMade = TRUE;
		for (int i = 0; i < 2 && IsOk; ++i){
			GPs[i].SetStartEndCPNum(-1, 0);
			IsOk = GPs[i].SeedInDirection((StreamDir_e)i);
			for (int j = 0; j < 2; ++j) {
				if (GPs[i].m_StartEndCPNum[j] >= 0 && VectorGetElementNum(EndCPNums, GPs[i].m_StartEndCPNum[j]) < 0)
					EndCPNums.push_back(GPs[i].m_StartEndCPNum[j]);
			}
// 			*this += GPs[i];
		}
		*this = GPs[0] + GPs[1];
		std::sort(EndCPNums.begin(), EndCPNums.end());
		for (int i = 0; i < 2 - EndCPNums.size(); ++i) EndCPNums.push_back(-1);
		m_StartEndCPNum[0] = EndCPNums.size() > 1 ? EndCPNums[1] : -1;
		m_StartEndCPNum[1] = EndCPNums[0];
	}

	if (IsOk && DoResample){
		Resample(m_NumGPPoints);
	}

	return IsOk;
}


Boolean_t GradPath_c::ReinterpolateRhoValuesFromVolume(VolExtentIndexWeights_s * VolInfo, FieldDataPointer_c * RhoPtr){
	FieldDataPointer_c * Ptr = &m_ODE_Data.RhoPtr;
	if (!Ptr->IsReady() && RhoPtr != nullptr)
		Ptr = RhoPtr;
	if (!Ptr->IsReady())
		return FALSE;

	VolExtentIndexWeights_s * vInfo = &m_ODE_Data.VolZoneInfo;
	if (!vInfo->IsReady() && VolInfo != nullptr)
		vInfo = VolInfo;
	if (!vInfo->IsReady())
		return FALSE;

	vector<VolExtentIndexWeights_s> ThVolInfo(omp_get_num_procs(), *vInfo);

#pragma omp parallel for
	for (int i = 0; i < m_XYZList.size(); ++i){
		m_RhoList[i] = ValAtPointByPtr(m_XYZList[i], ThVolInfo[omp_get_thread_num()], *Ptr);
	}

	return TRUE;
}

Boolean_t GradPath_c::ProjectPathToSurface(FESurface_c const * SurfPtr)
{
	REQUIRE(this->IsMade());
	REQUIRE((SurfPtr != nullptr && SurfPtr->IsMade()) || (m_Surface != nullptr && m_Surface->IsMade()));

	FESurface_c const * SPtr = (SurfPtr != nullptr ? SurfPtr : m_Surface);

	vector<Boolean_t> PointsInterior(omp_get_num_procs(), TRUE);
	vector<vec3> ProjectedPoint(omp_get_num_procs());

#pragma omp parallel for
	for (int i = 0; i < m_XYZList.size(); ++i){
		int ThNum = omp_get_thread_num();
		int ProjectedElemNum = -1;
		bool PointIsInterior;
		SPtr->ProjectPointToSurface(m_XYZList[i], ProjectedPoint[ThNum], ProjectedElemNum, PointIsInterior);
		m_XYZList[i] = ProjectedPoint[ThNum];
		PointsInterior[ThNum] = (PointsInterior[ThNum] && PointIsInterior);
	}

	Boolean_t AllPointsInterior = TRUE;
	for (auto i : PointsInterior)
		AllPointsInterior = (AllPointsInterior && i);

	return AllPointsInterior;
}

void GradPath_c::Clear(){
	GradPathBase_c::Clear();
}

void GradPath_c::SetStartPoint(vec3 const & pt){
	m_StartPoint = pt;
}

void GradPath_c::SetTerminalCPTypeNum(CPTypeNum_e CPTypeNum)
{
	m_TerminalCPTypeNums = {CPTypeNum};
}

void GradPath_c::SetTerminalCPTypeNums(vector<CPTypeNum_e> const & CPTypeNums)
{
	m_TerminalCPTypeNums = CPTypeNums;
}

void GradPath_c::MakeRhoValuesMonotomic(VolExtentIndexWeights_s * VolInfo, FieldDataPointer_c * RhoPtr)
{
	if (IsMade()){
		ReinterpolateRhoValuesFromVolume(VolInfo, RhoPtr);
		if (m_RhoList[0] < m_RhoList.back()) {
			std::sort(m_RhoList.begin(), m_RhoList.end());
		}
		else {
			std::sort(m_RhoList.rbegin(), m_RhoList.rend());
		}
		for (int i = 1; i < m_RhoList.size() - 1; ++i){
			if (m_RhoList[i-1] == m_RhoList[i] || m_RhoList[i] == m_RhoList[i+1]){
				m_RhoList[i] = (m_RhoList[i - 1] + m_RhoList[i + 1]) * 0.5;
			}
		}
	}
}

/*
*	Private Methods
*/

/*
*	Function for GSL ODE solver to use.
*	System of three first-order ODEs;
*	dX/dt = Grad(X),
*	dY/dt = Grad(Y),
*	dZ/dt = Grad(Z)
*
*	I need to provide the interpolation so
*	that arbitrary XYZ values can be queried, and to check
*	that they're in the bounds of the system.
*
*	This is NOT a member function, because it's a pain to
*	pass a member function to GSL.
*	Options are to make it static (unwanted because then it
*	would return the same value for all grad paths, I think...),
*	or to use std::bind and std::function to cast it properly, which
*	degrades performance, or to use template function-style casting,
*	which is scary and still uses static.
*/
int GP_ODE_Gradient(double t, double const pos[], double dydt[], void* params)
{
	GradPathParams_s *ODE_Data = reinterpret_cast<GradPathParams_s*>(params);

	int Status = GSL_SUCCESS;
	Boolean_t IsOk = TRUE;

	// 	vec3 TmpVec(pos);
	vec3 TmpVec = pos;

	/*
	*	Check that current position is in system bounds
	*/
// 	for (int i = 0; i < 3; ++i){
// 		if (TmpVec[i] < ODE_Data->VolZoneInfo.MinXYZ[i] || TmpVec[i] > ODE_Data->VolZoneInfo.MaxXYZ[i]){
// 			TmpVec[i] = MIN(ODE_Data->VolZoneInfo.MaxXYZ[i], MAX(TmpVec[i], ODE_Data->VolZoneInfo.MinXYZ[i]));
// 			Status = GSL_EDOM;
// 		}
// 	}

	if (!ODE_Data->VolZoneInfo.PointIsInterior(TmpVec))
		return GSL_EDOM;

	IsOk = SetIndexAndWeightsForPoint(TmpVec, ODE_Data->VolZoneInfo);

	/*
	*	Get gradient values at the actual position
	*/
	if (IsOk){
		vec3 TmpGrad;
		if (ODE_Data->HasGrad){
			for (int i = 0; i < 3; ++i){
				TmpGrad[i] = ValByCurrentIndexAndWeightsFromRawPtr(ODE_Data->VolZoneInfo, ODE_Data->GradPtrs[i]);
			}
		}
		else{
			CalcGradForPoint(TmpVec, ODE_Data->VolZoneInfo.PointSpacingV123, ODE_Data->VolZoneInfo, ODE_Data->VolZoneInfo.BasisNormalized, 0, ODE_Data->VolZoneInfo.IsPeriodic, TmpGrad, ODE_Data->RhoPtr, GPType_Invalid, nullptr);
		}

		if (ODE_Data->Direction == StreamDir_Reverse)
			TmpGrad *= -1.0;

// 		TmpVec = normalise(TmpVec);

		for (int i = 0; i < 3; ++i)
			dydt[i] = TmpGrad[i];
	}
	else
		Status = GSL_ESANITY;

	params = reinterpret_cast<void*>(ODE_Data);

	return Status;
} //	int GP_ODE_GradFunction()

/*
*	Jacobian for GSL ODE solver to use.
*/
int GP_ODE_Jacobian(double t, double const pos[], double *dfdy, double dydt[], void* params)
{
	GradPathParams_s *ODE_Data = reinterpret_cast<GradPathParams_s*>(params);

	int Status = GSL_SUCCESS;
	Boolean_t IsOk = TRUE;

	vec3 PosVec = pos;
	mat33 JacMat;

	/*
	*	Check that current position is in system bounds
	*/
// 	for (int i = 0; i < 3; ++i){
// 		if (PosVec[i] < ODE_Data->VolZoneInfo.MinXYZ[i] || PosVec[i] > ODE_Data->VolZoneInfo.MaxXYZ[i]){
// 			PosVec[i] = MIN(ODE_Data->VolZoneInfo.MaxXYZ[i], MAX(PosVec[i], ODE_Data->VolZoneInfo.MinXYZ[i]));
// 			Status = GSL_EDOM;
// 		}
// 	}
	if (!ODE_Data->VolZoneInfo.PointIsInterior(PosVec))
		return GSL_EDOM;

	IsOk = SetIndexAndWeightsForPoint(PosVec, ODE_Data->VolZoneInfo);

	/*
	*	Get gradient values at the actual position
	*/
	if (IsOk){
		vec3 TmpGrad;
		if (ODE_Data->HasGrad){
			for (int i = 0; i < 3; ++i){
				TmpGrad[i] = ValByCurrentIndexAndWeightsFromRawPtr(ODE_Data->VolZoneInfo, ODE_Data->GradPtrs[i]);
			}
		}
		else{
			CalcGradForPoint(PosVec, ODE_Data->VolZoneInfo.PointSpacingV123, ODE_Data->VolZoneInfo, ODE_Data->VolZoneInfo.BasisNormalized, 0, ODE_Data->VolZoneInfo.IsPeriodic, TmpGrad, ODE_Data->RhoPtr, GPType_Invalid, nullptr);
		}

		CalcHessForPoint(PosVec, ODE_Data->VolZoneInfo.PointSpacingV123, ODE_Data->VolZoneInfo, ODE_Data->VolZoneInfo.BasisNormalized, ODE_Data->VolZoneInfo.IsPeriodic, JacMat, ODE_Data->RhoPtr, GPType_Classic, nullptr);

		for (int i = 0; i < 3; ++i)
			dydt[i] = TmpGrad[i];

		gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 3, 3);
		gsl_matrix * m = &dfdy_mat.matrix;
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
			gsl_matrix_set(m, i, j, JacMat.at(i, j));
// 			dfdy[i * 3 + j] = JacMat.at(i, j);
	}
	else
		Status = GSL_ESANITY;

	params = reinterpret_cast<void*>(ODE_Data);

	return Status;
} //	int GP_ODE_Jacobian()


double GradPath_c::RhoByCurrentIndexAndWeights(){
	double Rho = 0;
	for (int i = 0; i < 8; ++i){
		Rho += m_ODE_Data.VolZoneInfo.Weights[i] * m_ODE_Data.RhoPtr[m_ODE_Data.VolZoneInfo.Index[i]];
	}

	return Rho;
}

/*
* Special gradient path algorithm functions
*/

/*
*	Function to return the actual function (gradient) value
*/
int F2DGrad(gsl_vector const * pos, void * params, gsl_vector * GradValues){
	MultiRootParams_s * RootParams = reinterpret_cast<MultiRootParams_s*>(params);

	vec3 ThreePoint = Transform2dTo3d(pos->data, *RootParams->BasisVectors, *RootParams->Origin);

	if (!SetIndexAndWeightsForPoint(ThreePoint, *RootParams->VolInfo))
		return GSL_EDOM;

	vec2 Grad;

	CalcGradForPoint(ThreePoint,
		RootParams->VolInfo->PointSpacingV123,
		*RootParams->VolInfo,
		*RootParams->BasisVectors,
		0,
		RootParams->IsPeriodic,
		Grad,
		*RootParams->RhoPtr,
		RootParams->CalcType, params);

	for (int i = 0; i < 2; ++i)
		gsl_vector_set(GradValues, i, Grad[i]);

	return GSL_SUCCESS;
}

/*
*	Function to return the derivatives (jacobian matrix of gradient, ie Hessian of rho)
*/
int DF2DGrad(gsl_vector const * pos, void * params, gsl_matrix * Jacobian){
	MultiRootParams_s * RootParams = reinterpret_cast<MultiRootParams_s*>(params);

	vec3 ThreePoint = Transform2dTo3d(pos->data, *RootParams->BasisVectors, *RootParams->Origin);

	if (!SetIndexAndWeightsForPoint(ThreePoint, *RootParams->VolInfo))
		return GSL_EDOM;

	mat22 J;

	CalcHessForPoint(ThreePoint,
		RootParams->VolInfo->PointSpacingV123,
		*RootParams->VolInfo,
		*RootParams->BasisVectors,
		RootParams->IsPeriodic, J,
		*RootParams->RhoPtr,
		RootParams->CalcType, params);

	for (int i = 0; i < 2; ++i)
	for (int j = 0; j < 2; ++j)
		gsl_matrix_set(Jacobian, i, j, J.at(i, j));

	return GSL_SUCCESS;
}

/*
*	Function to return both grad and dgrad (function and derivatives)
*/
int FDF2DGrad(gsl_vector const * pos, void * params, gsl_vector * GradValues, gsl_matrix * Jacobian){

	int Status = F2DGrad(pos, params, GradValues);

	if (Status == GSL_SUCCESS)
		Status = DF2DGrad(pos, params, Jacobian);

	return Status;
}

Boolean_t CPInNormalPlane(vec3 & StartPt, vec3 const & PlaneBasis, MultiRootObjects_s & MR){
	Boolean_t IsOk = TRUE;

	vector<vec3> BV(2);
	int Status = GSL_SUCCESS;
	int Iter = 0;

	MultiRootParams_s* Params = reinterpret_cast<MultiRootParams_s*>(MR.Func.params);

	/*
	*	Get the two orthonormal vectors to Dir, the grad path step direction
	*/

	Boolean_t DirFlipped = FALSE;

	BV[0] = PlaneBasis;
	for (int i = 0; i < 3 && !DirFlipped; ++i){
		if (abs(PlaneBasis[i]) > 1e-5){
			BV[0][i] *= -1.0;
			DirFlipped = TRUE;
		}
	}
	if (!DirFlipped) for (int i = 0; i < 3 && !DirFlipped; ++i){
		if (abs(PlaneBasis[i]) > 0){
			BV[0][i] *= -1.0;
			DirFlipped = TRUE;
		}
	}

	if (!DirFlipped) return FALSE;

	BV[0] = normalise(cross(PlaneBasis, BV[0]));
	BV[1] = normalise(cross(BV[0], PlaneBasis));

	mat33 BV3;
	BV3.col(0) = BV[0];
	BV3.col(1) = BV[1];

	Params->BasisVectors = &BV3;
	Params->Origin = &StartPt;

	for (int i = 0; i < 2; ++i)
		gsl_vector_set(MR.pos, i, 0.0);

	gsl_multiroot_fdfsolver_set(MR.s, &MR.Func, MR.pos);

#ifdef _DEBUG
	vec3 Pos = StartPt;
#endif
	double StepDist = 0.0;
	do
	{
		++Iter;

		Status = gsl_multiroot_fdfsolver_iterate(MR.s);

		if (Status != GSL_CONTINUE && Status != GSL_SUCCESS)
			break;

#ifdef _DEBUG
		StepDist = Distance(Pos, (Transform2dTo3d(MR.s->x->data, *Params->BasisVectors, StartPt)));
		Pos = Transform2dTo3d(MR.s->x->data, *Params->BasisVectors, StartPt);
#endif // _DEBUG

		Status = gsl_multiroot_test_residual(MR.s->f, 1e-9);

	} while (Status == GSL_CONTINUE && Iter < GP_PlaneCPMaxIter);



	if (Status == GSL_SUCCESS){
		vec3 EndPt = Transform2dTo3d(MR.s->x->data, *Params->BasisVectors, StartPt);

		if (Params->CalcType == GPType_NormalPlaneRhoCP){
			mat33 EigenVectors;
			vec3 EigenValues;
			CalcEigenSystemForPoint(EndPt,
				EigenValues,
				EigenVectors,
				*reinterpret_cast<MultiRootParams_s*>(MR.Func.params));
			char Rank = 0;
			for (int i = 0; i < 2; ++i){
				if (EigenValues[i] > 0)
					Rank++;
				else
					Rank--;
			}
			// 			m_RidgeRank = Rank;
		}

		StartPt = EndPt;

		IsOk = TRUE;
	}
	else {
		IsOk = FALSE;
	}
	MR.Func.params = reinterpret_cast<void*>(Params);

	return IsOk;
}

/*
*	Function for gradient bundle analysis that guesses whether or not
*	a pair of gradient paths straddle an irreducible bundle boundary
*	based on the angles of the last step of the paths.
*/
Boolean_t GPsStraddleIB(GradPath_c const & GP1,
	GradPath_c const & GP2,
	double const & IBCheckAngle,
	double const & IBCheckDistRatio)
{
	if (GP1.GetCount() > 1 && GP2.GetCount() > 1
		&& IBCheckAngle > 0 && IBCheckDistRatio > 0){
		vec3 Pts1[2] = { GP1[-1], GP1[-2] };
		vec3 Pts2[2] = { GP2[-1], GP2[-2] };

		double Dist = Distance(Pts1[0], Pts2[0]);

		vec3 V1 = Pts1[0] - Pts1[1];
		vec3 V2 = Pts2[0] - Pts2[1];

		double Angle = acos(dot(V1, V2) / (norm(V1) * norm(V2))) * 180.0 / PI;
		if (Angle > IBCheckAngle && Dist > IBCheckDistRatio)
			return TRUE;
	}

	return FALSE;
}


/*
*	Begin NEBGradPath_c methods
*/

NEBGradPath_c::NEBGradPath_c()
{
}

NEBGradPath_c::~NEBGradPath_c()
{
}

NEBGradPath_c::NEBGradPath_c(vec3 const & StartPt,
	vec3 const & EndPt,
	unsigned int NumPts)
{
	m_XYZList.reserve(NumPts);
	m_RhoList.resize(NumPts);

	m_NumGPPoints = NumPts;

	vec3 StepVec = (EndPt - StartPt) / (NumPts - 1);

	for (unsigned int StepNum = 0; StepNum < NumPts - 1; ++StepNum){
		m_XYZList.push_back(StartPt + (StepVec * StepNum));
	}

	// 	m_XYZList.back() = EndPt;

	m_XYZList.push_back(EndPt);
}


Boolean_t NEBGradPath_c::Relax(double const & StepRatio,
	double const & Tol,
	unsigned int const MaxIter,
	MultiRootParams_s & Params)
{
	REQUIRE(StepRatio >= 0.0 && StepRatio <= 1.0);
	REQUIRE(Tol > 0.0);
	REQUIRE(MaxIter > 1);

	Boolean_t IsOk;

	m_GradPathReady = m_GradPathMade = TRUE;

	vector<vec3> OldXYZ;
	vec3 OldPt;

	Params.EquilPos = &OldPt;

	double Cost;
	unsigned int Iter = 0;

	MultiRootObjects_s MR;

	MR.Func = { &F2DGrad, &DF2DGrad, &FDF2DGrad, 2, &Params };
	MR.pos = gsl_vector_alloc(2);
	MR.T = gsl_multiroot_fdfsolver_hybridsj;
	MR.s = gsl_multiroot_fdfsolver_alloc(MR.T, 2);

	do
	{
		Iter++;

		OldXYZ = m_XYZList;

		for (int Offset = 1; Offset <= 2; ++Offset){
			for (int i = Offset; i < m_XYZList.size() - 1; i += 2){
				OldPt = m_XYZList[i] * StepRatio + OldXYZ[i] * (1.0 - StepRatio);
				IsOk = CPInNormalPlane(m_XYZList[i], (m_XYZList[i + 1] - m_XYZList[i - 1]), MR);
			}
		}

#ifdef _DEBUG
		IsOk = SaveAsOrderedZone("NEB Iteration " + to_string(Iter) + " of " + to_string(MaxIter));
#endif // _DEBUG


		Cost = 0;
		for (int i = 1; i < m_XYZList.size() - 1; ++i)
			Cost += DistSqr(m_XYZList[i], OldXYZ[i]);

	} while (Cost > Tol && Iter <= MaxIter);

#ifdef _DEBUG
	IsOk = SaveAsOrderedZone("NEB Iteration " + to_string(Iter) + " of " + to_string(MaxIter));
#endif // _DEBUG

	IsOk = (Cost <= Tol && Iter <= MaxIter);

	return IsOk;
}


void StitchPaths(
	vector<int> const &     L,       // indices of points in P
	vector<int> const &     R,
	vector<vec3> const &     P,
	vector<vector<int> > &     T       // triplets of integers specifying nodes of triangles
	)
{
	int iL = 0, iR = 0;
	int nL = L.size() - 1, nR = R.size() - 1;
	// index of last element of L rather than count
	while (iL < nL || iR < nR)    // until exhaust both paths
	{
		int iL2 = MIN(iL + 1, nL),
			iR2 = MIN(iR + 1, nR);
		// next point along the path unless at end
		double dL = DistSqr(P[L[iL2]], P[R[iR]]),
			// length of next edge if we step down left path
			dR = DistSqr(P[L[iL]], P[R[iR2]]);
		// length of next edge if we step down right path
		// R and L are just indices, the actual points are in P
		if ((dL < dR && iL != iL2) || iR == iR2) {
			// if iR==iR2 then have to append left until finished
			T.push_back({ L[iL], L[iL2], R[iR] });
			iL = iL2;
		}
		else{// if ((dR < dL && iR != iR2) || iL == iL2) {
			T.push_back({ L[iL], R[iR2], R[iR] });
			iR = iR2;
		}
	}
}

void StitchCapPaths(
	vector<int> const &     L,       // indices of points in P
	vector<int> const &     R,
	vector<int> &			C,       // indices of points in the cap, C. that go from L[-1] to R[-1]
	vector<vec3> const &     P,
	vector<vector<int> > &     T,       // triplets of integers specifying nodes of triangles
	int LeftPathDeviationPointInd,
	int RightPathDeviationPointInd,
	int CapMidPointInd
	)
{
	int nC = C.size() - 1;

	if (nC < 0){
		// If there are no cap points (i.e. the endpoint are sufficiently close together)
		StitchPaths(L, R, P, T);
	}
	else {
		int jL = MIN(nC - 1, int(nC / 2) + 1);
		if (CapMidPointInd >= 0 && CapMidPointInd < C.size()) {
			if (DistSqr(P[C[0]], P[L.back()]) > DistSqr(P[C[0]], P[R.back()])){
				C = vector<int>(C.rbegin(), C.rend());
				CapMidPointInd = C.size() - 1 - CapMidPointInd;
			}
			jL = CapMidPointInd;
		}
		int jR = jL;
		if (LeftPathDeviationPointInd >= 0)
			LeftPathDeviationPointInd = CLAMP(LeftPathDeviationPointInd, 0, L.size() - 1);
		if (RightPathDeviationPointInd >= 0)
			RightPathDeviationPointInd = CLAMP(RightPathDeviationPointInd, 0, R.size() - 1);

		// midpoint of the cap, where we’ll transition to it
		int iL = 0, iR = 0;
		int nL = L.size() - 1, nR = R.size() - 1;
		while (iL < nL || iR < nR)
		{
			int iL2 = MIN(iL + 1, nL),
				iR2 = MIN(iR + 1, nR);
			double dL = DistSqr(P[L[iL2]], P[R[iR]]),
				dR = DistSqr(P[L[iL]], P[R[iR2]]),
				dLC = DistSqr(P[R[iR]], P[C[jR]]),
				dRC = DistSqr(P[L[iL]], P[C[jL]]);
			//  length of next edge from L/R if transition to cap
			if (dRC < dR / 2 
				|| dLC < dL / 2
				|| (LeftPathDeviationPointInd >= 0 && iL2 > LeftPathDeviationPointInd)
				|| (RightPathDeviationPointInd >= 0 && iR2 > RightPathDeviationPointInd))
			{
				// next point on paths is closer to cap than to other path, 
				// so switch to cap

				T.push_back({ L[iL], C[jL], R[iR] });
				// add transition triangle

				// Now stitch each half of the cap with the remaining legs of the paths

				vector<int> CVecL(&C[0], &C[jL] + 1);
				std::reverse(CVecL.begin(), CVecL.end());

				StitchPaths(vector<int>(&L[iL], &L[nL] + 1), CVecL, P, T);
				StitchPaths(vector<int>(&R[iR], &R[nR] + 1), vector<int>(&C[jR], &C[nC] + 1), P, T);

				break;
			}
			else if ((dL < dR && iL != iL2) || iR == iR2) {
				T.push_back({ L[iL], L[iL2], R[iR] });
				iL = iL2;
			}
			else{// if ((dR < dL && iR != iR2) || iL == iL2) {
				T.push_back({ L[iL], R[iR2], R[iR] });
				iR = iR2;
			}
		}
	}
}

int GradPathBase_c::GetSphereIntersectionPoint(vec3 const & Center, double const & Radius, vec3 & IntersectionPoint) const{
	if (!this->IsMade()) {
		TecUtilDialogErrMsg("Unmade gradient path being checked for sphere intersection");
		return -1;
	}
	/*
	 *	Since this is a check for intersection now with an arbitrary sphere, but
	 *	one with a nuclear CP at its center, any gradient path can intersect at most one time.
	 *	First verify that the path as a whole intersects the sphere (i.e. one end lies within
	 *	the radius).
	 */

	double RadSqr = Radius * Radius;
	double BegDistSqr = DistSqr(Center, this->XYZAt(0)),
		EndDistSqr = DistSqr(Center, this->XYZAt(-1));

	if ((BegDistSqr < RadSqr && EndDistSqr < RadSqr)
		|| (BegDistSqr >= RadSqr && EndDistSqr >= RadSqr))
		return -1;

	/*
	 *	Binary search on the grad path to find the point just inside the radius.
	 */
	int L = 0, R = this->m_XYZList.size() - 1, M;
	double TmpDistSqr;
	while (L != R){
		M = ceil(double(L + R) * 0.5);
		TmpDistSqr = DistSqr(Center, this->XYZAt(M));
		if (TmpDistSqr < RadSqr)
			R = M - 1;
		else
			L = M;
	}
	/*
	 *	Now m_XYZList[L] is just inside the radius, so interpolate
	 *	to get the intersection point.
	 */
	double LDist = Distance(Center, this->XYZAt(L)),
		RDist = Distance(Center, this->XYZAt(L + 1)),
		Ratio = (Radius - LDist) / (RDist - LDist);

	IntersectionPoint = this->XYZAt(L) + (this->XYZAt(L + 1) - this->XYZAt(L)) * Ratio;

	return L;
}