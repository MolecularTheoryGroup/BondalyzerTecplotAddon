

#include "TECADDON.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>

#include "CSM_DATA_TYPES.h"
#include "CSM_CALC_VARS.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"
#include "CSM_CRIT_POINTS.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GRAD_PATH.h"

#include <armadillo>
using namespace arma;

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
		double MinSqrDist = DistSqr(m_XYZList[GetCount() - 1], rhs[0]);
		int MinPair = 1;
		double TempSqrDist = DistSqr(m_XYZList[GetCount() - 1], rhs[-1]);
		if (MinSqrDist > TempSqrDist){
			MinSqrDist = TempSqrDist;
			MinPair = 2;
		}
		TempSqrDist = DistSqr(m_XYZList[0], rhs[0]);
		if (MinSqrDist > TempSqrDist){
			MinSqrDist = TempSqrDist;
			MinPair = 3;
		}
		TempSqrDist = DistSqr(m_XYZList[0], rhs[-1]);
		if (MinSqrDist > TempSqrDist){
			MinSqrDist = TempSqrDist;
			MinPair = 4;
		}

		int offset = 0; // used to 

		switch (MinPair){
			case 1:
				if (approx_equal(this->XYZAt(-1), rhs.XYZAt(0), "absdiff", 0.01)) 
					offset = 1;
				m_XYZList.insert(m_XYZList.end(), rhs.m_XYZList.cbegin() + offset, rhs.m_XYZList.cend());
				m_RhoList.insert(m_RhoList.end(), rhs.m_RhoList.cbegin() + offset, rhs.m_RhoList.cend());
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[1];
				break;
			case 2:
				if (approx_equal(this->XYZAt(-1), rhs.XYZAt(-1), "absdiff", 0.01))
					offset = 1;
				m_XYZList.insert(m_XYZList.end(), rhs.m_XYZList.crbegin() + offset, rhs.m_XYZList.crend());
				m_RhoList.insert(m_RhoList.end(), rhs.m_RhoList.crbegin() + offset, rhs.m_RhoList.crend());
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[0];
				break;
			case 3:
				if (approx_equal(this->XYZAt(0), rhs.XYZAt(-1), "absdiff", 0.01))
					offset = 1;
				m_XYZList.insert(m_XYZList.begin(), rhs.m_XYZList.crbegin() + offset, rhs.m_XYZList.crend());
				m_RhoList.insert(m_RhoList.begin(), rhs.m_RhoList.crbegin() + offset, rhs.m_RhoList.crend());
				m_StartEndCPNum[0] = m_StartEndCPNum[1];
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[1];
				break;
			case 4:
				if (approx_equal(this->XYZAt(0), rhs.XYZAt(0), "absdiff", 0.01))
					offset = 1;
				m_XYZList.insert(m_XYZList.begin(), rhs.m_XYZList.cbegin() + offset, rhs.m_XYZList.cend());
				m_RhoList.insert(m_RhoList.begin(), rhs.m_RhoList.cbegin() + offset, rhs.m_RhoList.cend());
				m_StartEndCPNum[0] = m_StartEndCPNum[1];
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[0];
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
	if (m_Length >= 0)
		return m_Length;

	m_Length = 0.0;

	if (m_GradPathMade){
		for (int i = 0; i < GetCount() - 1; ++i)
			m_Length += Distance(m_XYZList[i], m_XYZList[i + 1]);
	}

	return m_Length;
}

GradPathBase_c GradPathBase_c::SubGP(int BegPt, int EndPt) const{
	BegPt = GetInd(BegPt);
	EndPt = GetInd(EndPt);

	GradPathBase_c Out = *this;
	Out.m_XYZList.assign(m_XYZList.begin() + BegPt, m_XYZList.begin() + EndPt + 1);
	Out.m_RhoList.assign(m_RhoList.begin() + BegPt, m_RhoList.begin() + EndPt + 1);
	Out.GradPathBase_c::m_NumGPPoints = Out.m_XYZList.size();

	return Out;
}

vec3 GradPathBase_c::ClosestPoint(vec3 const & CheckPt) const
{
	int iJunk;

	return ClosestPoint(CheckPt, iJunk);
}

vec3 GradPathBase_c::ClosestPoint(vec3 const & CheckPt, int & PtNum) const
{
	vec3 ClosestPt;
	if (IsMade()){
		ClosestPt = m_XYZList[0];
		double MinSqrDist = DistSqr(m_XYZList[0], CheckPt);
		PtNum = 0;
		int Count = GetCount();
		for (int i = 1; i < Count; ++i){
			double TempSqrDist = DistSqr(m_XYZList[i], CheckPt);
			if (TempSqrDist < MinSqrDist){
				MinSqrDist = TempSqrDist;
				ClosestPt = m_XYZList[i];
				PtNum = i;
			}
		}
	}

	return ClosestPt;
}

Boolean_t GradPathBase_c::Resample(int NumPoints){
	Boolean_t IsOk = m_GradPathMade && NumPoints > 1;

	vector<vec3> NewXYZList;
	vector<double> NewRhoList;

	int OldCount = GetCount();

// 	if (IsOk && NumPoints < OldCount){
	if (IsOk){
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

		for (int NewI = 1; NewI < NumPoints - 1; ++NewI){
			ArcLength += DelLength;

			while (OldI < OldCount - 1 && ArcLengthI < ArcLength){
				++OldI;

				ArcLengthIm1 = ArcLengthI;
				PtIm1 = PtI;
				RhoIm1 = RhoI;

				PtI = m_XYZList[OldI];
				RhoI = m_RhoList[OldI];

				ArcLengthI += Distance(PtI, PtIm1);
			}

			double Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
			NewXYZList[NewI] = PtIm1 + (PtI - PtIm1) * Ratio;

			NewRhoList[NewI] = RhoIm1 + Ratio * (RhoI - RhoIm1);

			if (OldI >= OldCount){
				while (NewI < NumPoints){
					NewI++;
					if (NewI < NumPoints){
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

		if (IsOk){
			IsOk = (NewRhoList.size() == NumPoints);
		}

		if (IsOk){
			m_XYZList.swap(NewXYZList);

			m_RhoList.swap(NewRhoList);
		}

		m_Length = -1;
	}
	else IsOk = FALSE;

	m_NumGPPoints = m_XYZList.size();

	return IsOk;
}

Boolean_t GradPathBase_c::Reverse(){
	Boolean_t IsOk = m_GradPathMade;
	if (IsOk){
		vector<vec3> NewXYZList;
		vector<double> NewRhoList;

		NewXYZList.insert(NewXYZList.begin(), m_XYZList.crbegin(), m_XYZList.crend());
		m_XYZList = NewXYZList;

		NewRhoList.insert(NewRhoList.begin(), m_RhoList.crbegin(), m_RhoList.crend());
		m_RhoList = NewRhoList;

		int TmpInt = m_StartEndCPNum[0];
		m_StartEndCPNum[0] = m_StartEndCPNum[1];
		m_StartEndCPNum[1] = TmpInt;
	}
	return IsOk;
}

/*
*	Concatenates and resamples at once.
*	Resampling a concatenated grad path that has a
*	sharp bend can cut the corner at the bend, so
*	this method avoids that.
*/
GradPathBase_c & GradPathBase_c::ConcatenateResample(GradPathBase_c & rhs, int NumPoints)
{
	int iJunk;

	ConcatenateResample(rhs, NumPoints, iJunk);

	return *this;
}

GradPathBase_c & GradPathBase_c::ConcatenateResample(GradPathBase_c & rhs, int NumPoints, int & BrigePtNum)
{
	double MyLength = GetLength();
	double rhsLength = rhs.GetLength();
	int MyNumPtsInit = GetCount();
	int rhsNumPtsInit = rhs.GetCount();
	double TotalLength = MyLength + rhsLength;

	int MyNumPts = static_cast<int>((MyLength / TotalLength) * static_cast<double>(NumPoints));
	MyNumPts = MIN(MyNumPts, MyNumPtsInit);
	BrigePtNum = MyNumPts - 1;
	int rhsNumPts = MIN(NumPoints - MyNumPts, rhsNumPtsInit);

	Resample(MyNumPts);
	GradPathBase_c NewRHS(rhs);
	NewRHS.Resample(rhsNumPts);

	Concatenate(NewRHS);

	return *this;
}

GradPath_c ConcatenateResample(vector<GradPath_c> GPList, int NumPoints)
{
	vector<double> LenghtList;
	vector<int> NumPointsList;
	for (auto gp : GPList)	LenghtList.push_back(gp.GetLength());

	double TotalLength = 0.0;
	for (auto l : LenghtList) TotalLength += l;

	for (int li = 0; li < LenghtList.size() - 1; ++li)
		NumPointsList.push_back(LenghtList[li] / TotalLength * double(NumPoints));

	int TotNumPoints = 0;
	for (auto n : NumPointsList) TotNumPoints += n;
	NumPointsList.push_back(NumPoints - TotNumPoints);

	GradPath_c GP = GPList[0];
	GP.Resample(NumPointsList[0]);
	for (int i = 1; i < GPList.size(); ++i){
		GPList[i].Resample(NumPointsList[i]);
		GP += GPList[i];
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

		vector<vec3> NewXYZList;

		vector<double> NewRhoList;

		if (DistSqr(Point, (m_XYZList[0])) < DistSqr(Point, (m_XYZList[Count - 1]))){
			/*
			*	Initial point is closer to Point, so start from beginning
			*/
			double RhoI, RhoIm1 = m_RhoList[0];
			vec3 PtI, PtIm1 = m_XYZList[0];
			double DistI, DistIm1 = Distance(PtIm1, Point);

			for (int i = 1; i < Count; ++i){
				PtI = m_XYZList[i];
				RhoI = m_RhoList[i];
				DistI = Distance(PtI, Point);

				if ((DistI < Radius && DistIm1 >= Radius) ||
					(DistI >= Radius && DistIm1 < Radius)){
					// Now PtI and PtIm1 are straddling the sphere.
					// Interpolate to find the intersection point.
					double DistanceRatio = (Radius - DistI) / (DistIm1 - DistI);

					vec3 NewPoint = (PtI + ((PtIm1 - PtI) * DistanceRatio)) - Point;
					double NewRho = (RhoI + ((RhoIm1 - RhoI) * DistanceRatio));

					/*
					* Project intersection point to sphere surface
					*/
					// Compute spherical coordinates
					double   rad = norm(NewPoint);
					double theta = acos(NewPoint[2] / rad);
					double   phi = atan2(NewPoint[1], NewPoint[0]);

					// Project point onto a sphere of radius "radius"
					NewPoint[0] = Radius * sin(theta) * cos(phi);
					NewPoint[1] = Radius * sin(theta) * sin(phi);
					NewPoint[2] = Radius * cos(theta);

					NewPoint += Point;

					NewXYZList.reserve(m_XYZList.size() - i + 1);
					NewXYZList.push_back(NewPoint);
					NewXYZList.insert(NewXYZList.end(), m_XYZList.cbegin() + i, m_XYZList.cend());

					m_XYZList[i] = NewXYZList[i];

					NewRhoList.reserve(m_RhoList.size() - i + 1);
					NewRhoList.push_back(NewRho);
					NewRhoList.insert(NewRhoList.end(), m_RhoList.cbegin() + i, m_RhoList.cend());
					m_RhoList = NewRhoList;

					break;
				}
				PtIm1 = PtI;
				RhoIm1 = RhoI;
				DistIm1 = DistI;
			}
		}
		else{
			/*
			*	Terminal point is closer to Point, so start from end
			*/
			double RhoI, RhoIm1 = m_RhoList[Count - 1];
			vec3 PtI, PtIm1 = m_XYZList[Count - 1];
			double DistI, DistIm1 = Distance(PtIm1, Point);

			for (int i = Count - 2; i >= 0; --i){
				PtI = m_XYZList[i];
				RhoI = m_RhoList[i];
				DistI = Distance(PtI, Point);

				if ((DistI < Radius && DistIm1 >= Radius) ||
					(DistI >= Radius && DistIm1 < Radius)){
					// Now PtI and PtIm1 are straddling the sphere.
					// Interpolate to find the intersection point.
					double DistanceRatio = (Radius - DistI) / (DistIm1 - DistI);

					vec3 NewPoint = (PtI + ((PtIm1 - PtI) * DistanceRatio)) - Point;
					double NewRho = (RhoI + ((RhoIm1 - RhoI) * DistanceRatio));

					/*
					* Project intersection point to sphere surface
					*/
					// Compute spherical coordinates
					double   rad = sqrt(pow(NewPoint[0], 2) + pow(NewPoint[1], 2) + pow(NewPoint[2], 2));
					double theta = acos(NewPoint[2] / rad);
					double   phi = atan2(NewPoint[1], NewPoint[0]);

					// Project point onto a sphere of radius "radius"
					NewPoint[0] = Radius * sin(theta) * cos(phi);
					NewPoint[1] = Radius * sin(theta) * sin(phi);
					NewPoint[2] = Radius * cos(theta);

					NewPoint += Point;

					NewXYZList.reserve(i + 1);
					NewXYZList.insert(NewXYZList.end(), m_XYZList.cbegin(), m_XYZList.cbegin() + i);
					NewXYZList.push_back(NewPoint);

					m_XYZList = NewXYZList;

					NewRhoList.reserve(i + 1);
					NewRhoList.insert(NewRhoList.end(), m_RhoList.cbegin(), m_RhoList.cbegin() + i);
					NewRhoList.push_back(NewRho);
					m_RhoList = NewRhoList;

					break;
				}
				PtIm1 = PtI;
				RhoIm1 = RhoI;
				DistIm1 = DistI;
			}
		}

		IsOk = (NewXYZList.size() > 0 && NewRhoList.size() > 0);
	}

	return IsOk;
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

			IsOk = SaveAsOrderedZone(ZoneName, VarTypes, XYZVarNums, RhoVarNum, TRUE, MeshColor);
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

				for (int i = 0; i < 3; ++i){
					FieldData_pa GPXYZPtr = TecUtilDataValueGetWritableNativeRef(m_ZoneNum, XYZVarNums[i]);
					TecUtilDataValueArraySetByRef(GPXYZPtr, 1, GetCount(), TmpValues[i].data());
				}
				if (RhoVarNum > 0){
					FieldData_pa RhoFDPtr = TecUtilDataValueGetWritableNativeRef(m_ZoneNum, RhoVarNum);
					TecUtilDataValueArraySetByRef(RhoFDPtr, 1, GetCount(), m_RhoList.data());
				}

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

vec3 GradPathBase_c::operator[](int i) const{
	return m_XYZList[GetInd(i)];
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
		int GPSize = GP_NumPointsBufferFactor * GP_MaxNumPoints;
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
		int GPSize = GP_NumPointsBufferFactor * GP_MaxNumPoints;
		m_XYZList.reserve(GPSize);
		m_RhoList.reserve(GPSize);
	}

	return m_GradPathReady;
}

Boolean_t GradPath_c::SeedInDirection(StreamDir_e const & Direction){
	Boolean_t IsOk = m_GradPathReady && !m_GradPathMade;

	StreamDir_e OldDir = m_ODE_Data.Direction;
	m_ODE_Data.Direction = Direction;

	gsl_odeiv2_system ODESys = { &GP_ODE_Gradient, &GP_ODE_Jacobian, m_ODE_NumDims, &m_ODE_Data };

// 	gsl_odeiv2_step_type const * T = gsl_odeiv2_step_rk4;
	gsl_odeiv2_step_type const * T = gsl_odeiv2_step_bsimp;
// 	gsl_odeiv2_step_type const * T = gsl_odeiv2_step_msbdf;
	gsl_odeiv2_step * s = gsl_odeiv2_step_alloc(T, m_ODE_NumDims);
	gsl_odeiv2_control * c = gsl_odeiv2_control_y_new(1e-7, 1e-9);
	gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc(m_ODE_NumDims);



	// 	gsl_odeiv2_driver * ODEDriver;
	// 	ODEDriver = gsl_odeiv2_driver_alloc_yp_new(&ODESys, gsl_odeiv2_step_rk2, 1e-3, 1e-2, 0);
	// 	gsl_odeiv2_driver_set_hmin(ODEDriver, 1e-8);
	// 	gsl_odeiv2_driver_set_hmax(ODEDriver, 1e-1);

	if (IsOk){
		double tInit = 0.0;
		double tFinal = 1e12;

		double h = 0.01;

		double y[3] = { m_StartPoint[0], m_StartPoint[1], m_StartPoint[2] };

		int SurfLastProjectedElem = -1;
		bool SurfProjectionFound;

		if (m_Surface != nullptr && m_Surface->IsMade()){
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

			if (Status == GSL_SUCCESS || Status == GSL_EDOM){
				if (m_Surface != nullptr && m_Surface->IsMade()) {
					if (!m_Surface->ProjectPointToSurface(y, PtI, SurfLastProjectedElem, SurfProjectionFound)) {
						TecUtilDialogErrMsg("Projection of point to FESurface failed");
					}
					for (int i = 0; i < 3; ++i) y[i] = PtI[i];
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
						for (auto const & CPTypeNum : CPNearFieldTypes) {
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
				if (m_TermValue != -1.0 && Rho < m_TermValue){
					double OldRho;
					OldRho = m_RhoList.back();

					NewPoint = PtIm1 + (PtI - PtIm1) * ((m_TermValue - OldRho) / (Rho - OldRho));

					m_XYZList.push_back(NewPoint);
					m_RhoList.push_back(m_TermValue);

					break;
				}

				if (IsOk){
					if (m_XYZList.size() >= GP_StallNumPointsToCheck - 1){
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

					PtIm1 = PtI;
				}
			}

			Step++;
		}

		IsOk = (IsOk && Status == GSL_SUCCESS || Status == GSL_EDOM || Status == GSL_ENOPROG);

		// 		if (Status == GSL_ENOPROG && m_XYZList.size() > GP_StallPointCount){
		// 			/*
		// 			*	Grad path stalled, so it was basically bouncing around the same point.
		// 			*	The last point can then be approximated as the midpoint of the stalled points.
		// 			*/
		// 			vec3 Pt = zeros(3);
		// 			double PtRho = 0.0;
		// 			double PtCount = 0.0;
		// 			for (int i = MAX(0, m_XYZList.size() - GP_StallPointCount); i < m_XYZList.size(); ++i){
		// 				Pt += m_XYZList[i];
		// 				PtRho += m_RhoList[i];
		// 				PtCount += 1.0;
		// 			}
		// 
		// 			Pt /= PtCount;
		// 			PtRho /= PtCount;
		// 
		// 			m_XYZList.resize(MAX(0, m_XYZList.size() - GP_StallPointCount));
		// 			m_XYZList.back() = Pt;
		// 			m_RhoList.resize(MAX(0, m_RhoList.size() - GP_StallPointCount));
		// 			m_RhoList.back() = PtRho;
		// 			PtI = Pt;
		// 		}

		// 		gsl_odeiv2_driver_free(ODEDriver);

		gsl_odeiv2_evolve_free(e);
		gsl_odeiv2_control_free(c);
		gsl_odeiv2_step_free(s);

		if (m_GPType && MR.s != nullptr)
			gsl_multiroot_fdfsolver_free(MR.s);
		// 	if (m_GPType && MR.pos != nullptr)
		// 		gsl_vector_free(MR.pos);

		if (m_StartEndCPNum[1] < 0 && (m_HowTerminate == GPTerminate_AtCP || m_HowTerminate == GPTerminate_AtCPRadius)){
			/*
			*	Check to see if terminating point coincides with a CP
			*/
			Boolean_t PointFound = FALSE;
			for (int CPNum = 0; CPNum < m_NumCPs && !PointFound; ++CPNum){
				if (CPNum != m_StartEndCPNum[0]){
					if (m_CPs == nullptr){
						if (m_CPXYZPtrs.size() == 3) {
							for (int i = 0; i < 3; ++i) {
								NewPoint[i] = m_CPXYZPtrs[i][CPNum];
							}
						}
						else{
							NewPoint = m_CPs->GetXYZ(CPNum);
						}
					}
					else{
						NewPoint = m_CPs->GetXYZ(CPNum);
					}
					double PointRadiusSqr = DistSqr(PtI, NewPoint);
					if (PointRadiusSqr <= m_TermPointRadiusSqr){
						m_StartEndCPNum[1] = CPNum;
						PointFound = TRUE;
					}
				}
			}
		}

		m_GradPathMade = IsOk;
		if (IsOk && m_GPType)
			m_SGPMade = TRUE;
	}

	m_ODE_Data.Direction = OldDir;

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


Boolean_t GradPath_c::ReinterpolateRhoValuesFromVolume(VolExtentIndexWeights_s * VolInfo){
	if (!m_ODE_Data.RhoPtr.IsReady())
		return FALSE;

	if (VolInfo != nullptr)
		m_ODE_Data.VolZoneInfo = *VolInfo;
	for (int i = 0; i < m_XYZList.size(); ++i){
		m_RhoList[i] = ValAtPointByPtr(m_XYZList[i], m_ODE_Data.VolZoneInfo, m_ODE_Data.RhoPtr);
	}

	return TRUE;
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
	for (int i = 0; i < 3; ++i){
		if (TmpVec[i] < ODE_Data->VolZoneInfo.MinXYZ[i] || TmpVec[i] > ODE_Data->VolZoneInfo.MaxXYZ[i]){
			TmpVec[i] = MIN(ODE_Data->VolZoneInfo.MaxXYZ[i], MAX(TmpVec[i], ODE_Data->VolZoneInfo.MinXYZ[i]));
			Status = GSL_EDOM;
		}
	}

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
			CalcGradForPoint(TmpVec, ODE_Data->VolZoneInfo.DelXYZ, ODE_Data->VolZoneInfo, ODE_Data->VolZoneInfo.BasisNormalized, 0, ODE_Data->VolZoneInfo.IsPeriodic, TmpGrad, ODE_Data->RhoPtr, GPType_Invalid, nullptr);
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
	for (int i = 0; i < 3; ++i){
		if (PosVec[i] < ODE_Data->VolZoneInfo.MinXYZ[i] || PosVec[i] > ODE_Data->VolZoneInfo.MaxXYZ[i]){
			PosVec[i] = MIN(ODE_Data->VolZoneInfo.MaxXYZ[i], MAX(PosVec[i], ODE_Data->VolZoneInfo.MinXYZ[i]));
			Status = GSL_EDOM;
		}
	}

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
			CalcGradForPoint(PosVec, ODE_Data->VolZoneInfo.DelXYZ, ODE_Data->VolZoneInfo, ODE_Data->VolZoneInfo.BasisNormalized, 0, ODE_Data->VolZoneInfo.IsPeriodic, TmpGrad, ODE_Data->RhoPtr, GPType_Invalid, nullptr);
		}

		CalcHessForPoint(PosVec, ODE_Data->VolZoneInfo.DelXYZ, ODE_Data->VolZoneInfo, ODE_Data->VolZoneInfo.BasisNormalized, ODE_Data->VolZoneInfo.IsPeriodic, JacMat, ODE_Data->RhoPtr, GPType_Classic, nullptr);

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
		RootParams->VolInfo->DelXYZ,
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
		RootParams->VolInfo->DelXYZ,
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
	vector<int> const &		C,       // indices of points in the cap, C. that go from L[-1] to R[-1]
	vector<vec3> const &     P,
	vector<vector<int> > &     T       // triplets of integers specifying nodes of triangles
	)
{
	int nC = C.size() - 1;
	if (nC < 0){
		// If there are no cap points (i.e. the endpoint are sufficiently close together)
		StitchPaths(L, R, P, T);
	}
	else{
		int jL = MIN(nC - 1, int(nC / 2) + 1);
		int jR = jL;
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
			if (dRC < dR / 2 || dLC < dL / 2) {
				// next point on paths is closer to cap than to other path, 
				// so switch to cap

				T.push_back({ L[iL], C[jL], R[iR] });
				// add transition triangle

				// Now stitch each half of the cap with the remaining legs of the paths
				// 				vector<int> CVecL(&C[0], &C[jL]+1),
				// 					CVecR(&C[jR], &C[nC]+1),
				// 					LVec, RVec;
				// 				if (CVecL.size() == 0){
				// 					CVecL = C;
				// 					CVecR = C;
				// 				}
				// 				else std::reverse(CVecL.begin(), CVecL.end());
				// 
				// 				if (iL < nL) LVec.assign(&L[iL], &L[nL]+1);
				// 				else LVec.push_back(L[nL]);
				// 
				// 				if (iR < nR) RVec.assign(&R[iR], &R[nR]+1);
				// 				else RVec.push_back(R[nR]);
				// 
				// 				StitchPaths(LVec, CVecL, P, T);
				// 				StitchPaths(RVec, CVecR, P, T);

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

bool GradPathBase_c::GetSphereIntersectionPoint(vec3 const & Center, double const & Radius, vec3 & IntersectionPoint) const{
	if (!this->IsMade()) {
		TecUtilDialogErrMsg("Unmade gradient path being checked for sphere intersection");
		return false;
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
		return false;

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

	return true;
}