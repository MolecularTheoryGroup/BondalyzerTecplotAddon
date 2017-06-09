

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
#include "CSM_GRAD_PATH.h"

#include <armadillo>
using namespace arma;

using std::list;
using std::vector;
using std::string;
using std::to_string;


int F2DGrad(const gsl_vector * pos, void * params, gsl_vector * GradValues);
int DF2DGrad(const gsl_vector * pos, void * params, gsl_matrix * Jacobian);
int FDF2DGrad(const gsl_vector * pos, void * params, gsl_vector * GradValues, gsl_matrix * Jacobian);

/*
 *	GradPathParams_s methods
 */
GradPathParams_s & GradPathParams_s::operator=(const GradPathParams_s & rhs)
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
}//	GradPathParams_s & GradPathParams_s::operator=(const GradPathParams_s & rhs)
const Boolean_t GradPathParams_s::operator==(const GradPathParams_s & rhs) const
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
}// Boolean_t GradPathParams_s::operator==(const GradPathParams_s & rhs) const


/*
 *	GradPathBase_c methods
 */


GradPathBase_c::GradPathBase_c()
{
	m_StartEndCPNum[0] = m_StartEndCPNum[1] = -1;
	m_GradPathReady = FALSE;
	m_GradPathMade = FALSE;
	m_NumGPPoints = -1;
}

GradPathBase_c::~GradPathBase_c()
{
}



/*
*	Constructor for loading an existing 1-D ordered
*	zone as a GradPath_c
*/
GradPathBase_c::GradPathBase_c(EntIndex_t ZoneNum,
	const vector<EntIndex_t> & XYZRhoVarNums,
	const AddOn_pa & AddOnID)
{
	m_GradPathMade = FALSE;
	m_ZoneNum = ZoneNum;
	Boolean_t IsOk = TecUtilZoneIsOrdered(ZoneNum);
	if (!IsOk){
		TecUtilDialogErrMsg("Can't import grad path: Zone is not ordered.");
		TecUtilLockFinish(AddOnID);
		return;
	}
	if (!TecUtilZoneIsActive(ZoneNum)){
		Set_pa TmpSet = TecUtilSetAlloc(FALSE);
		TecUtilSetAddMember(TmpSet, ZoneNum, FALSE);
		TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);
		TecUtilSetDealloc(&TmpSet);
	}

	LgIndex_t MaxIJK[3];
	TecUtilZoneGetIJK(ZoneNum, &MaxIJK[0], &MaxIJK[1], &MaxIJK[2]);
	IsOk = (MaxIJK[0] > 1 && MaxIJK[1] == 1 && MaxIJK[2] == 1);
	if (!IsOk){
		TecUtilDialogErrMsg("Can't import grad path: Zone is not 1-D.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	vector<FieldDataPointer_c> RawPtrs(4);
	for (int i = 0; i < 4; ++i){
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
			for (int j = 0; j < MaxIJK[0]; ++j)
				m_RhoList[j] = RawPtrs[i][j];
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
GradPathBase_c & GradPathBase_c::operator=(const GradPathBase_c & rhs)
{
	if (this == &rhs)
		return *this;

	m_XYZList = rhs.m_XYZList;

	m_RhoList = rhs.m_RhoList;

	for (int i = 0; i < 2; ++i)
		m_StartEndCPNum[i] = rhs.m_StartEndCPNum[i];

	m_GradPathReady = rhs.m_GradPathReady;
	m_GradPathMade = rhs.m_GradPathMade;

	m_ZoneNum = rhs.m_ZoneNum;

	m_NumGPPoints = rhs.m_NumGPPoints;

	return *this;
}// GradPathBase_c & GradPathBase_c::operator=(const GradPathBase_c & rhs)
const Boolean_t GradPathBase_c::IsSame(const GradPathBase_c & rhs) const
{
	Boolean_t AreSame = (
		m_RhoList == rhs.m_RhoList &&

		m_GradPathReady == rhs.m_GradPathReady &&
		m_GradPathMade == rhs.m_GradPathMade
		);
		
	if (AreSame){
		for (int i = 0; i < m_XYZList.size() && AreSame; ++i)
			AreSame = sum(m_XYZList[i] == rhs.m_XYZList[i]) == 3;
	}

	return AreSame;
}
const Boolean_t GradPathBase_c::operator==(const GradPathBase_c & rhs) const
{
	return IsSame(rhs);
}
GradPathBase_c & GradPathBase_c::operator+=(const GradPathBase_c & rhs)
{
	Boolean_t IsOk = (m_GradPathMade && rhs.IsMade());
	if (IsOk){

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

		switch (MinPair){
			case 1:
				m_XYZList.insert(m_XYZList.end(), rhs.m_XYZList.cbegin(), rhs.m_XYZList.cend());
				m_RhoList.insert(m_RhoList.end(), rhs.m_RhoList.cbegin(), rhs.m_RhoList.cend());
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[1];
				break;
			case 2:
				m_XYZList.insert(m_XYZList.end(), rhs.m_XYZList.crbegin(), rhs.m_XYZList.crend());
				m_RhoList.insert(m_RhoList.end(), rhs.m_RhoList.crbegin(), rhs.m_RhoList.crend());
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[0];
				break;
			case 3:
				m_XYZList.insert(m_XYZList.begin(), rhs.m_XYZList.crbegin(), rhs.m_XYZList.crend());
				m_RhoList.insert(m_RhoList.begin(), rhs.m_RhoList.crbegin(), rhs.m_RhoList.crend());
				m_StartEndCPNum[0] = m_StartEndCPNum[1];
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[1];
				break;
			case 4:
				m_XYZList.insert(m_XYZList.begin(), rhs.m_XYZList.cbegin(), rhs.m_XYZList.cend());
				m_RhoList.insert(m_RhoList.begin(), rhs.m_RhoList.cbegin(), rhs.m_RhoList.cend());
				m_StartEndCPNum[0] = m_StartEndCPNum[1];
				m_StartEndCPNum[1] = rhs.m_StartEndCPNum[0];
				break;
		}
	}
	return *this;
}
const GradPathBase_c GradPathBase_c::operator+(const GradPathBase_c & rhs) const
{
	return GradPathBase_c(*this) += rhs;
}


/*
*	Getter methods
*/

const double GradPathBase_c::GetLength() const {
	double Length = 0.0;

	if (m_GradPathMade){
		for (int i = 0; i < GetCount() - 1; ++i)
			Length += Distance(m_XYZList[i], m_XYZList[i + 1]);
	}

	return Length;
}

const GradPathBase_c GradPathBase_c::SubGP(int BegPt, int EndPt) const{
	if (0 > EndPt || EndPt >= m_NumGPPoints)
		EndPt = m_NumGPPoints - 1;
	if (0 > BegPt)
		BegPt = 0;
	else if (BegPt > EndPt)
		BegPt = EndPt;
	REQUIRE(0 <= BegPt && BegPt <= EndPt && EndPt < m_NumGPPoints);
	
	GradPathBase_c Out = *this;
	Out.m_XYZList.assign(m_XYZList.begin() + BegPt, m_XYZList.begin() + EndPt + 1);
	Out.m_RhoList.assign(m_RhoList.begin() + BegPt, m_RhoList.begin() + EndPt + 1);
	Out.GradPathBase_c::m_NumGPPoints = Out.m_XYZList.size();

	return Out;
}

const vec3 GradPathBase_c::ClosestPoint(const vec3 & CheckPt) const
{
	int iJunk;

	return ClosestPoint(CheckPt, iJunk);
}

const vec3 GradPathBase_c::ClosestPoint(const vec3 & CheckPt, int & PtNum) const
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

const Boolean_t GradPathBase_c::Resample(const int & NumPoints){
	Boolean_t IsOk = m_GradPathMade && NumPoints > 1;

	vector<vec3> NewXYZList;
	vector<double> NewRhoList;

	int OldCount = GetCount();

	if (IsOk && NumPoints < OldCount){
		NewXYZList.resize(NumPoints);

		NewRhoList.resize(NumPoints);

		double DelLength = GetLength() / static_cast<double>(NumPoints);

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
	}
	else IsOk = FALSE;

	return IsOk;
}

const Boolean_t GradPathBase_c::Reverse(){
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
GradPathBase_c & GradPathBase_c::ConcatenateResample(const GradPathBase_c & rhs, const int & NumPoints)
{
	int iJunk;

	ConcatenateResample(rhs, NumPoints, iJunk);

	return *this;
}

GradPathBase_c & GradPathBase_c::ConcatenateResample(const GradPathBase_c & rhs, const int & NumPoints, int & BrigePtNum)
{
	double MyLength = GetLength();
	double rhsLength = rhs.GetLength();
	int MyNumPtsInit = GetCount();
	int rhsNumPtsInit = rhs.GetCount();
	double TotalLength = MyLength + rhsLength;

	int MyNumPts = static_cast<int>((MyLength / TotalLength) * static_cast<double>(NumPoints));
	MyNumPts = MIN(MyNumPts, MyNumPtsInit);
	BrigePtNum = MyNumPts-1;
	int rhsNumPts = MIN(NumPoints - MyNumPts, rhsNumPtsInit);

	Resample(MyNumPts);
	GradPathBase_c NewRHS(rhs);
	NewRHS.Resample(rhsNumPts);

	Concatenate(NewRHS);

	return *this;
}


/*
*	For making a grad path terminate at a particular point.
*	In this version, at the intersection of the grad path and
*	a sphere of r = Radius around Point.
*/
const Boolean_t GradPathBase_c::Trim(const vec3 & Point, const double & Radius)
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

void GradPathBase_c::PointAppend(const vec3 & Point, const double & Rho){

	m_XYZList.push_back(Point);

	m_RhoList.push_back(Rho);

	m_NumGPPoints++;
}


void GradPathBase_c::PointPrepend(const vec3 & Point, const double & Rho){

	m_XYZList.insert(m_XYZList.begin(), Point);

	m_RhoList.insert(m_RhoList.begin(), Rho);

	m_NumGPPoints++;
}

const Boolean_t GradPathBase_c::SaveAsOrderedZone(const string & ZoneName, const ColorIndex_t MeshColor)
{
	Boolean_t IsOk = m_GradPathMade;

	if (IsOk){
		EntIndex_t VolZoneNum = ZoneNumByName(string("Full Volume"));
		EntIndex_t RhoVarNum = VarNumByName(string("Electron Density"));

		vector<EntIndex_t> XYZVarNums(3);
		TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);

		IsOk = (VolZoneNum > 0 && RhoVarNum > 0 &&
			XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);

		if (IsOk){
			EntIndex_t NumVars = TecUtilDataSetGetNumVars();
			vector<FieldDataType_e> VarTypes(NumVars);
			for (int i = 0; i < NumVars; ++i){
				VarTypes[i] = TecUtilDataValueGetType(VolZoneNum, i + 1);
			}

			for (int i = 0; i < 3; ++i)
				VarTypes[XYZVarNums[i] - 1] = FieldDataType_Double;

			VarTypes[RhoVarNum - 1] = FieldDataType_Double;

			IsOk = SaveAsOrderedZone(ZoneName, VarTypes, XYZVarNums, RhoVarNum, MeshColor);
		}
	}

	return IsOk;
}

const EntIndex_t GradPathBase_c::SaveAsOrderedZone(const string & ZoneName, 
	vector<FieldDataType_e> & VarDataTypes,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum,
	const Boolean_t DoActivate,
	const ColorIndex_t MeshColor)
{
	Boolean_t IsOk = (m_XYZList.size() > 0 && m_XYZList.size() == m_RhoList.size());

	if (IsOk){

		if (IsOk){

			if (VarDataTypes.size() == TecUtilDataSetGetNumVars())
				IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), GetCount(), 1, 1, ZoneType_Ordered, VarDataTypes.data());
			else IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), GetCount(), 1, 1, ZoneType_Ordered, NULL);

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
				FieldData_pa RhoFDPtr = TecUtilDataValueGetWritableNativeRef(m_ZoneNum, RhoVarNum);
				TecUtilDataValueArraySetByRef(RhoFDPtr, 1, GetCount(), m_RhoList.data());

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

const Boolean_t GradPathBase_c::SaveAsCSV(const string & PathToFile, const Boolean_t & IncludeVars){
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
GradPath_c::GradPath_c() : GradPathBase_c()
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
GradPath_c::GradPath_c(const vec3 & StartPoint,
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
								const FieldDataPointer_c & RhoPtr)
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


GradPath_c::GradPath_c(const vec3 & StartPoint,
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
	const FieldDataPointer_c & RhoPtr)
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
		RhoPtr);
} //	GradPath_c::GradPath_c()

GradPath_c::GradPath_c(EntIndex_t ZoneNum,
	const vector<EntIndex_t> & XYZRhoVarNums,
	const AddOn_pa & AddOnID) : GradPathBase_c(ZoneNum, XYZRhoVarNums, AddOnID){}

/*
 *	Copy constructor
 */
GradPath_c::GradPath_c(const GradPath_c & a)
{
	*this = a;
}
GradPath_c::GradPath_c(const GradPathBase_c & a)
{
	*this = a;
}

GradPath_c::~GradPath_c()
{
}

/*
*	Operator declarations
*/
GradPath_c & GradPath_c::operator=(const GradPath_c & rhs)
{
	if (this == &rhs)
		return *this;

	GradPathBase_c::operator=(rhs);

	m_ODE_Data = rhs.m_ODE_Data;

	m_TermPoint = rhs.m_TermPoint;
	m_CPXYZPtrs = rhs.m_CPXYZPtrs;
	m_NumCPs = rhs.m_NumCPs;
	m_CPs = rhs.m_CPs;

	m_TermPointRadiusSqr = rhs.m_TermPointRadiusSqr;
	m_TermValue = rhs.m_TermValue;

	m_HowTerminate = rhs.m_HowTerminate;
	m_StartPoint = rhs.m_StartPoint;

	return *this;
}// GradPath_c & GradPath_c::operator=(const GradPath_c & rhs)
GradPath_c & GradPath_c::operator=(const GradPathBase_c & rhs)
{
	if (this == &rhs)
		return *this;

	GradPathBase_c::operator=(rhs);

	return *this;
}// GradPath_c & GradPath_c::operator=(const GradPath_c & rhs)
const Boolean_t GradPath_c::operator==(const GradPath_c & rhs) const
{
	Boolean_t AreSame = this->IsSame(rhs);

	if (AreSame){
		AreSame = (sum(m_TermPoint == rhs.m_TermPoint) == 3 &&
			m_CPXYZPtrs == rhs.m_CPXYZPtrs &&
			m_NumCPs == rhs.m_NumCPs &&
			m_CPs == rhs.m_CPs &&
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
const Boolean_t GradPath_c::SetupGradPath(const vec3 & StartPoint,
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
											const FieldDataPointer_c & RhoPtr)
{
	m_XYZList.swap(vector<vec3>());
	m_RhoList.swap(vector<double>());

	m_NumGPPoints = NumGPPoints;

	m_StartPoint = StartPoint;
	m_ODE_Data.Direction = Direction;
	m_HowTerminate = HowTerminate;

	if (TermPoint != NULL)
		m_TermPoint = *TermPoint;
	if (TermPointRadius != NULL)
		m_TermPointRadiusSqr = (*TermPointRadius * *TermPointRadius);

	m_CPs = NULL;

	m_CPXYZPtrs = CPXYZPtrs;

	if (NumCPs != NULL)
		m_NumCPs = *NumCPs;

	if (TermValue != NULL)
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
			m_GradPathReady = (TermPoint != NULL && TermPointRadius != NULL);
		else if (HowTerminate == GPTerminate_AtRhoValue)
			m_GradPathReady = (TermValue != NULL);
		else if (HowTerminate == GPTerminate_AtCP || m_HowTerminate == GPTerminate_AtCPRadius){
			m_GradPathReady = (CPXYZPtrs.size() == 3 && TermPointRadius != NULL && NumCPs != NULL);
			if (m_GradPathReady){
				for (int i = 0; i < 3 && m_GradPathReady; ++i){
					m_GradPathReady = CPXYZPtrs[i].IsReady();
				}
			}
		}
	}

	m_StartEndCPNum[0] = m_StartEndCPNum[1] = -1;

	m_GradPathMade = FALSE;

	return m_GradPathReady;
}

const Boolean_t GradPath_c::SetupGradPath(const vec3 & StartPoint,
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
	const FieldDataPointer_c & RhoPtr)
	{
		m_XYZList.swap(vector<vec3>());
		m_RhoList.swap(vector<double>());

		m_NumGPPoints = NumGPPoints;

		m_StartPoint = StartPoint;
		m_ODE_Data.Direction = Direction;
		m_HowTerminate = HowTerminate;
		m_GPType = GPType;

		if (TermPoint != NULL)
			m_TermPoint = *TermPoint;
		if (TermPointRadius != NULL)
			m_TermPointRadiusSqr = (*TermPointRadius * *TermPointRadius);

		m_NumCPs = -1;

		m_CPs = CPs;

		if (TermValue != NULL)
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
				m_GradPathReady = (TermPoint != NULL && TermPointRadius != NULL);
			else if (HowTerminate == GPTerminate_AtRhoValue)
				m_GradPathReady = (TermValue != NULL);
			else if (HowTerminate == GPTerminate_AtCP || m_HowTerminate == GPTerminate_AtCPRadius){
				m_GradPathReady = (m_CPs != NULL && TermPointRadius != NULL && m_CPs->NumCPs() > 0);
				if (m_GradPathReady) m_NumCPs = m_CPs->NumCPs();
			}
		}

		m_StartEndCPNum[0] = m_StartEndCPNum[1] = -1;

		m_GradPathMade = FALSE;

		return m_GradPathReady;
	}

const Boolean_t GradPath_c::Seed(const bool DoResample){
	Boolean_t IsOk = m_GradPathReady && !m_GradPathMade;

	gsl_odeiv2_system ODESys = { &GP_ODE_GradFunction, NULL, m_ODE_NumDims, &m_ODE_Data };
	gsl_odeiv2_driver * ODEDriver;

	ODEDriver = gsl_odeiv2_driver_alloc_yp_new(&ODESys, gsl_odeiv2_step_rk2, 1e-3, 1e-2, 0);
	gsl_odeiv2_driver_set_hmin(ODEDriver, 1e-8);
	gsl_odeiv2_driver_set_hmax(ODEDriver, 1e-1);

	if (IsOk){
		double tInit = 0.0;
		double tFinal = 10.0;

		double y[3] = { m_StartPoint[0], m_StartPoint[1], m_StartPoint[2] };

		vec3 NewPoint, PtIm1 = m_StartPoint, PtI;

		int GPSize = GP_NumPointsBufferFactor * m_NumGPPoints;

		IsOk = SetIndexAndWeightsForPoint(m_StartPoint, m_ODE_Data.VolZoneInfo);

		if (IsOk){
			m_XYZList.reserve(GPSize);
			m_XYZList.push_back(m_StartPoint);

			m_RhoList.reserve(GPSize);
			m_RhoList.push_back(RhoByCurrentIndexAndWeights());
		}


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
			Params.BasisVectors = &BV;
			Params.Origin = &BV[0];

			MR.Func = { &F2DGrad, &DF2DGrad, &FDF2DGrad, 2, &Params};
			MR.pos = gsl_vector_alloc(2);
			MR.T = gsl_multiroot_fdfsolver_hybridsj;
			MR.s = gsl_multiroot_fdfsolver_alloc(MR.T, 2);
		}

		while (IsOk && Status == GSL_SUCCESS && Step < GP_MaxNumPoints){


			Status = gsl_odeiv2_driver_apply(ODEDriver, &tInit, tInit + 1e-3, y);

			if (Status == GSL_SUCCESS || Status == GSL_EDOM){
				PtI = y;

				if (m_GPType != GPType_Classic && m_GPType != GPType_Invalid){
					/*
					*	This allows a mixing of the last direction, as determined by the found 2d CP, and
					*	the current direction, as determined by the gradient at the last GP point.
					*/
//   					if (!SGPFound){
					StepDir = normalise( (PtI - PtIm1));

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
						Params.BasisVectors = &BV;
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
							Params.BasisVectors = &BV;
							PlaneCPDist = Distance(PtI, TmpPt);
						} while (SGPFound && PlaneCPDist >= 1e-4 && PlaneCPIter < 50);
						if (PlaneCPIter < 50 && PlaneCPDist < 1e-4){
							/*
							 *	This allows a mixing of the last direction, as determined by the found 2d CP, and
							 *	the current direction, as determined by the gradient at the last GP point.
							 */
							// 						StepDir = normalise((StepDir) * m_DirMixFactor + normalise( (PtI - PtIm1)) * normalise( (1.0 - m_DirMixFactor)));
							PtI = PtIm1 + normalise( (PtI - PtIm1)) * StepSize;
							for (int i = 0; i < 3; ++i)
								y[i] = PtI[i];
						}
					}
// 					if (!SGPFound){
					else{
						Params.BasisVectors = &BV;
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
					for (int CPNum = 0; CPNum < m_NumCPs && !PointFound; ++CPNum){
						if (CPNum != m_StartEndCPNum[0]){
							if (m_CPs == NULL){
								for (int i = 0; i < 3; ++i){
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
					if (PointFound)
						break;
				}
				else if (m_TermValue != -1.0 && Rho < m_TermValue){
					double OldRho;
					OldRho = m_RhoList[m_RhoList.size() - 1];

					NewPoint = PtIm1 + (PtI - PtIm1) * ((m_TermValue - OldRho) / (Rho - OldRho));

					m_XYZList.push_back(NewPoint);
					m_RhoList.push_back(m_TermValue);

					break;
				}

				if (IsOk){
					double NewPtDistSqr = DistSqr(PtI, (PtIm1));
					if (NewPtDistSqr < GP_PointSqrDistTol){
						NumStalledPoints++;
						if (NumStalledPoints >= GP_StallPointCount && Status == GSL_SUCCESS){
							Status = GSL_ENOPROG;
						}
					}
					else
						NumStalledPoints = 0;

					m_XYZList.push_back(PtI);
					m_RhoList.push_back(Rho);

					PtIm1 = PtI;
				}
			}

			Step++;
		}

		IsOk = (IsOk && Status == GSL_SUCCESS || Status == GSL_EDOM || Status == GSL_ENOPROG);

		gsl_odeiv2_driver_free(ODEDriver);

		if (m_GPType && MR.s != NULL)
			gsl_multiroot_fdfsolver_free(MR.s);
		// 	if (m_GPType && MR.pos != NULL)
		// 		gsl_vector_free(MR.pos);

		IsOk = m_GradPathMade = (IsOk && m_RhoList.size() > 0);
		if (IsOk && m_GPType)
			m_SGPMade = TRUE;


		if (IsOk && DoResample){
			Resample(m_NumGPPoints);
		}
	}

	return IsOk;
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
 *	Of course, I need to provide the interpolation so
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
int GP_ODE_GradFunction(double t, const double pos[], double dydt[], void* params)
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
				TmpVec[i] = ValByCurrentIndexAndWeightsFromRawPtr(ODE_Data->VolZoneInfo, ODE_Data->GradPtrs[i]);
			}
		}
		else{
			vector<vec3> BV(3);
			BV[0] << 1 << 0 << 0;
			BV[1] << 0 << 1 << 0;
			BV[2] << 0 << 0 << 1;
			CalcGradForPoint(TmpVec, ODE_Data->VolZoneInfo.DelXYZ, ODE_Data->VolZoneInfo, BV, 0, ODE_Data->VolZoneInfo.IsPeriodic, TmpGrad, ODE_Data->RhoPtr, GPType_Invalid, NULL);
			TmpVec = TmpGrad;
		}

		if (ODE_Data->Direction == StreamDir_Reverse)
			TmpVec *= -1.0;

		TmpVec = normalise(TmpVec);

		for (int i = 0; i < 3; ++i)
			dydt[i] = TmpVec[i];
	}
	else
		Status = GSL_ESANITY;

	params = reinterpret_cast<void*>(ODE_Data);

	return Status;
} //	int GP_ODE_GradFunction()


const double GradPath_c::RhoByCurrentIndexAndWeights(){
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
int F2DGrad(const gsl_vector * pos, void * params, gsl_vector * GradValues){
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
int DF2DGrad(const gsl_vector * pos, void * params, gsl_matrix * Jacobian){
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
int FDF2DGrad(const gsl_vector * pos, void * params, gsl_vector * GradValues, gsl_matrix * Jacobian){

	int Status = F2DGrad(pos, params, GradValues);

	if (Status == GSL_SUCCESS)
		Status = DF2DGrad(pos, params, Jacobian);

	return Status;
}

const Boolean_t CPInNormalPlane(vec3 & StartPt, const vec3 & PlaneBasis, MultiRootObjects_s & MR){
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

	BV[0] = normalise( cross(PlaneBasis, BV[0]));
	BV[1] = normalise( cross(BV[0], PlaneBasis));

	Params->BasisVectors = &BV;
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
const Boolean_t GPsStraddleIB(const GradPath_c & GP1,
	const GradPath_c & GP2,
	const double & IBCheckAngle,
	const double & IBCheckDistRatio)
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

NEBGradPath_c::NEBGradPath_c(const vec3 & StartPt, 
	const vec3 & EndPt, 
	const unsigned int & NumPts)
{
	m_XYZList.reserve(NumPts);
	m_RhoList.resize(NumPts);

	m_NumGPPoints = NumPts;

	vec3 StepVec = (EndPt - StartPt) / (NumPts - 1);

	for (unsigned int StepNum = 0; StepNum < NumPts - 1; ++StepNum){
		m_XYZList.push_back(StartPt + (StepVec * StepNum));
	}

// 	m_XYZList[m_XYZList.size() - 1] = EndPt;

	m_XYZList.push_back(EndPt);
}


const Boolean_t NEBGradPath_c::Relax(const double & StepRatio,
	const double & Tol, 
	const unsigned int MaxIter, 
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