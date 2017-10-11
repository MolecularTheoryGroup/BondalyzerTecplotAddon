#include "TECADDON.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <armadillo>

#include "CSM_GRAD_PATH.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_FE_VOLUME.h"

// #define RECORD_INT_POINTS


using std::vector;
using std::string;
using std::to_string;
using std::stringstream;
using std::cout;
using std::endl;

using namespace arma;


/*
*	Public member functions
*/

/*
 *	Operators
 */

FESurface_c & FESurface_c::operator+=(const FESurface_c & rhs){
	int NodeCount = m_XYZList.size(),
		ElemCount = m_ElemList.size();

	// Copy nodes and elements from rhs
	m_XYZList.insert(m_XYZList.end(), rhs.m_XYZList.cbegin(), rhs.m_XYZList.cend());
	m_ElemList.insert(m_ElemList.end(), rhs.m_ElemList.cbegin(), rhs.m_ElemList.cend());

	// Fix new element numbers
	for (int e = ElemCount; e < m_ElemList.size(); ++e){
		for (int & i : m_ElemList[e]){
			i += NodeCount;
		}
	}

	m_NumElems = m_ElemList.size();
	m_NumNodes = m_XYZList.size();

	m_FEVolumeMade = (m_NumNodes > 0 && m_NumElems > 0);

	return *this;
}

const FESurface_c FESurface_c::operator+(const FESurface_c & rhs) const{
	return FESurface_c(*this) += rhs;
}

/*
*	Constructor/destructor
*/

FESurface_c::FESurface_c()
{
	m_NumGPs = -1;
	m_FEVolumeMade = FALSE;
	m_IntegrationResultsReady = FALSE;
}

FESurface_c::FESurface_c(const int & InZoneNum,
	const ZoneType_e & InZoneType,
	const vector<FieldDataPointer_c> & InXYZPtrs,
	const vec3 & InMaxXYZ,
	const vec3 & InMinXYZ,
	const vector<int> InMaxIJK,
	NodeMap_t* InConnectivityListPtr)
{
	m_FEVolumeMade = TRUE;
	m_IntegrationResultsReady = FALSE;

	m_XYZPtrs = InXYZPtrs;

	m_FEVolumeMade = m_XYZPtrs.size() == 3;

	for (int i = 0; i < 3 && m_FEVolumeMade; ++i){
		m_FEVolumeMade = m_XYZPtrs[i].IsReady();
	}

	if (m_FEVolumeMade)
		m_FEVolumeMade = (InMaxIJK.size() == 3);

	if (m_FEVolumeMade)
		m_FEVolumeMade = (InConnectivityListPtr != NULL);

	if (m_FEVolumeMade){
		TecUtilZoneGetIJK(InZoneNum, &m_NumNodes, &m_NumElems, &m_NumNodesPerElem);
		m_FEVolumeMade = ((InZoneType == ZoneType_FETriangle && m_NumNodesPerElem == 3)
			|| (InZoneType == ZoneType_FEQuad && m_NumNodesPerElem == 4));
		if (m_FEVolumeMade){
			m_VolZoneInfo.MaxIJK = InMaxIJK;
			m_VolZoneInfo.MaxXYZ = InMaxXYZ;
			m_VolZoneInfo.MinXYZ = InMinXYZ;
			m_ConnectivityListPtr = InConnectivityListPtr;
		}
	}
}

/*
*	Make FEVolume from lists of nodes (points) and elements (lists of point indices).
*/
FESurface_c::FESurface_c(const vector<vec3> & Nodes, vector<vector<int> > & Elements)
{
	m_FEVolumeMade = (Nodes.size() >= 3 && Elements.size() >= 1);
	if (m_FEVolumeMade){
		m_XYZList.clear();
		m_XYZList.insert(m_XYZList.begin(), Nodes.cbegin(), Nodes.cend());
		m_ElemList.insert(m_ElemList.begin(), Elements.cbegin(), Elements.cend());
	}
}

FESurface_c::FESurface_c(const int & InZoneNum,
	const int & VolZoneNum,
	const vector<int> & InXYZVarNums,
	const vector<int> & InIntVarNums)
{
	Setup(InZoneNum, VolZoneNum, InXYZVarNums, InIntVarNums);
}

FESurface_c::FESurface_c(const int & InZoneNum,
	const vector<int> & InXYZVarNums)
{
	Setup(InZoneNum, InXYZVarNums);
}


/*
	Make FEVolume from vector of pointers to GradPath_c's
*/
const Boolean_t FESurface_c::MakeGradientBundle(vector<GradPath_c*> GPs)
{
	Boolean_t IsOk = TRUE;
	for (int i = 0; i < GPs.size() && IsOk; ++i){
		IsOk = GPs[i]->IsMade();
		if (IsOk && i > 0){
			IsOk = (GPs[i]->GetCount() == GPs[i-1]->GetCount());
		}
	}
	
	if (IsOk){
		m_NumGPs = GPs.size();
		m_NumGPPts = GPs[0]->GetCount();
		m_NumNodes = m_NumGPs * m_NumGPPts;
		m_NumElems = 2 * (m_NumGPs - 2) + 2 * m_NumGPs * (m_NumGPPts - 1);

		vector<vector<vec3>::const_iterator> XYZIt(m_NumGPs);
		
		for (int i = 0; i < GPs.size(); ++i){
			XYZIt[i] = GPs[i]->m_XYZList.cbegin();
		}

		m_XYZList.resize(m_NumGPs * m_NumGPPts);

		for (int i = 0; i < m_NumGPs; ++i){
			for (int j = 0; j < m_NumGPPts; ++j){
				m_XYZList[i * m_NumGPPts + j] = *XYZIt[i];
				XYZIt[i]++;
			}
		}

		/*
		 * Generate triangle list.
		 */

		vector<int> V(m_NumGPs);
		for (int i = 0; i < m_NumGPs; ++i){
			V[i] = i * m_NumGPPts;
		}
 
		Domain_c D(V, this);
		D.Weight();
 
		TriPolyLines();
 
		V = vector<int>();
		V.reserve(m_NumGPs * 2);
		int OldNumNodes = m_XYZList.size() - m_NumGPs;
		for (int i = 0; i < m_NumGPs; ++i){
			V.push_back((i + 1) * m_NumGPPts - 1);
			V.push_back(OldNumNodes + i);
		}
		D.Setup(V, this);
		D.Weight();

		m_NumElems = m_ElemList.size();

		m_FEVolumeMade = TRUE;
	}

	return IsOk;
}

/*
Make FEVolume from vector of pointers to GradPath_c's
*/
const Boolean_t FESurface_c::MakeFromGPs(vector<GradPath_c*> GPs, const bool ConnectBeginningAndEndGPs)
{
	Boolean_t IsOk = TRUE;
	for (int i = 0; i < GPs.size() && IsOk; ++i){
		IsOk = GPs[i]->IsMade();
		if (IsOk && i > 0){
			IsOk = (GPs[i]->GetCount() == GPs[i - 1]->GetCount());
		}
	}

	if (IsOk){
		m_NumGPs = GPs.size();
		m_NumGPPts = GPs[0]->GetCount();
		m_NumNodes = m_NumGPs * m_NumGPPts;
		m_NumElems = 2 * (m_NumGPs - 2) + 2 * m_NumGPs * (m_NumGPPts - 1);

		vector<vector<vec3>::const_iterator> XYZIt(m_NumGPs);

		m_XYZList.clear();
		for (GradPath_c* GP : GPs)
			m_XYZList.insert(m_XYZList.begin(), GP->m_XYZList.cbegin(), GP->m_XYZList.cend());

// 		for (int i = 0; i < GPs.size(); ++i){
// 			XYZIt[i] = GPs[i]->m_XYZList.cbegin();
// 		}
// 
// 		m_XYZList.resize(m_NumGPs * m_NumGPPts);
// 
// 		for (int i = 0; i < m_NumGPs; ++i){
// 			for (int j = 0; j < m_NumGPPts; ++j){
// 				m_XYZList[i * m_NumGPPts + j] = *XYZIt[i];
// 				XYZIt[i]++;
// 			}
// 		}

		/*
		* Generate triangle list.
		*/

		TriPolyLines(ConnectBeginningAndEndGPs);

		m_NumElems = m_ElemList.size();

		m_FEVolumeMade = TRUE;
	}

	return IsOk;
}

/*
*	Two constructors for the 3- and 4-sided
*	FE volumes.
*/
const Boolean_t FESurface_c::MakeGradientBundle(const GradPath_c & GP1,
	const GradPath_c & GP2,
	const GradPath_c & GP3)
{
	Boolean_t IsOk = (GP1.IsMade()
		&& GP2.IsMade()
		&& GP3.IsMade()
		&& GP1.GetCount() == GP2.GetCount()
		&& GP1.GetCount() == GP3.GetCount());


	if (IsOk){
		m_GPList = { GP1, GP2, GP3 };

		m_NumGPs = 3;
		m_NumGPPts = GP1.GetCount();
		m_NumNodes = 3 * m_NumGPPts;
		m_NumElems = 3 * (m_NumGPPts - 1) + 2;

		vector<vec3>::const_iterator XYZIt[3];
		vector<double>::const_iterator RhoIt[3];

		int GPNum = 0;
		XYZIt[GPNum] = GP1.m_XYZList.cbegin();

		RhoIt[GPNum] = GP1.m_RhoList.cbegin();
		GPNum++;


		XYZIt[GPNum] = GP2.m_XYZList.cbegin();

		RhoIt[GPNum] = GP2.m_RhoList.cbegin();
		GPNum++;


		XYZIt[GPNum] = GP3.m_XYZList.cbegin();

		RhoIt[GPNum] = GP3.m_RhoList.cbegin();

		m_XYZList.resize(3 * m_NumGPPts);

		for (int j = 0; j < m_NumGPPts; ++j){
			for (int i = 0; i < 3; ++i){
				m_XYZList[i + 3 * j] = *XYZIt[i];
				XYZIt[i]++;
			}
		}

		m_RhoList.resize(3 * m_NumGPPts);

		for (int j = 0; j < m_NumGPPts; ++j){
			for (int i = 0; i < 3; ++i){
				m_RhoList[i + 3 * j] = *RhoIt[i];
				RhoIt[i]++;
			}
		}

		m_FEVolumeMade = TRUE;
	}

	return IsOk;
}
const Boolean_t FESurface_c::MakeGradientBundle(const GradPath_c & GP1,
	const GradPath_c & GP2,
	const GradPath_c & GP3,
	const GradPath_c & GP4)
{
	Boolean_t IsOk = (GP1.IsMade()
		&& GP2.IsMade()
		&& GP3.IsMade()
		&& GP4.IsMade()
		&& GP1.GetCount() == GP2.GetCount()
		&& GP1.GetCount() == GP3.GetCount()
		&& GP1.GetCount() == GP4.GetCount());


	if (IsOk){
		m_GPList = { GP1, GP2, GP3, GP4 };

		m_NumGPs = 4;
		m_NumGPPts = GP1.GetCount();
		m_NumNodes = 4 * m_NumGPPts;
		m_NumElems = 4 * (m_NumGPPts - 1) + 2;

		vector<vec3>::const_iterator XYZIt[4];
		vector<double>::const_iterator RhoIt[4];

		int GPNum = 0;
		XYZIt[GPNum] = GP1.m_XYZList.cbegin();

		RhoIt[GPNum] = GP1.m_RhoList.cbegin();
		GPNum++;

		XYZIt[GPNum] = GP2.m_XYZList.cbegin();

		RhoIt[GPNum] = GP2.m_RhoList.cbegin();
		GPNum++;

		XYZIt[GPNum] = GP3.m_XYZList.cbegin();

		RhoIt[GPNum] = GP3.m_RhoList.cbegin();
		GPNum++;

		XYZIt[GPNum] = GP4.m_XYZList.cbegin();

		RhoIt[GPNum] = GP4.m_RhoList.cbegin();

		m_XYZList.resize(4 * m_NumGPPts);

		for (int j = 0; j < m_NumGPPts; ++j){
			for (int i = 0; i < 4; ++i){
				m_XYZList[i + 4 * j] = *XYZIt[i];
				XYZIt[i]++;
			}
		}

		m_RhoList.resize(4 * m_NumGPPts);

		for (int j = 0; j < m_NumGPPts; ++j){
			for (int i = 0; i < 4; ++i){
				m_RhoList[i + 3 * j] = *RhoIt[i];
				RhoIt[i]++;
			}
		}

		m_FEVolumeMade = TRUE;
	}

	return IsOk;
}

FESurface_c::~FESurface_c()
{
}

/*
*	Getter methods
*/

const vector<double> FESurface_c::GetIntResults() const
{
	if (m_IntegrationResultsReady){
		return m_IntValues;
	}
	else
		return vector<double>(1, 0);
}


const int FESurface_c::SetZoneStyle(const int ZoneNum,
	const AssignOp_e ZoneActive,
	const Boolean_t ShowContour,
	const Boolean_t ShowMesh,
	const Boolean_t ShowScatter,
	const Boolean_t ShowShade)
{
	if (ZoneNum > 0 && ZoneNum <= TecUtilDataSetGetNumZones()){
		Set_pa ZoneSet = TecUtilSetAlloc(TRUE);
		TecUtilSetAddMember(ZoneSet, ZoneNum, TRUE);

		TecUtilZoneSetActive(ZoneSet, ZoneActive);
		TecUtilZoneSetMesh(SV_SHOW, ZoneSet, 0.0, ShowMesh);
		TecUtilZoneSetScatter(SV_SHOW, ZoneSet, 0.0, ShowScatter);
		TecUtilZoneSetContour(SV_SHOW, ZoneSet, 0.0, ShowContour);
		TecUtilZoneSetShade(SV_SHOW, ZoneSet, 0.0, ShowShade);
	}
	return ZoneNum;
}

/*
	Save an FEVolume_c that has it's elements and nodal structure
	already created.
	Only Saves the XYZ values.
*/
const int FESurface_c::SaveAsTriFEZone(const vector<int> & XYZVarNums, 
	string ZoneName)
{
	Boolean_t IsOk = IsMade();
	int ZoneNum = -1;

	if (ZoneName.length() == 0)
		ZoneName = "FE Volume";

	if (m_RefinedXYZList.size() <= 0) m_RefinedXYZList = m_XYZList;

	IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), m_RefinedXYZList.size(), m_ElemList.size(), 0, ZoneType_FETriangle, NULL);

	Set_pa TmpSet = TecUtilSetAlloc(FALSE);

	if (IsOk){
		ZoneNum = TecUtilDataSetGetNumZones();
		m_ZoneNum = ZoneNum;
		TecUtilSetAddMember(TmpSet, ZoneNum, FALSE);
		TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
		TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);

		ArgList_pa CurrentArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
		TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
		TecUtilStyleSetLowLevelX(CurrentArgList);
		TecUtilArgListDealloc(&CurrentArgList);

		TecUtilSetDealloc(&TmpSet);

		vector<vector<double> > TmpValues(3, vector<double>(m_RefinedXYZList.size()));
		for (int i = 0; i < m_RefinedXYZList.size(); ++i){
			for (int j = 0; j < 3; ++j)
				TmpValues[j][i] = m_RefinedXYZList[i][j];
		}

		for (int i = 0; i < 3 && IsOk; ++i){
			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
			IsOk = VALID_REF(SetFDPtr);
			if (IsOk){
				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_RefinedXYZList.size(), reinterpret_cast<void*>(const_cast<double*>(TmpValues[i].data())));
			}
		}

	}
	
	if (IsOk){
		NodeMap_pa NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
		IsOk = VALID_REF(NodeMap);
		for (int ei = 0; ei < m_ElemList.size(); ++ei){
			for (int i = 0; i < 3; ++i){
				TecUtilDataNodeSetByRef(NodeMap, ei + 1, i + 1, m_ElemList[ei][i] + 1);
			}
		}
	}

	TecUtilSetDealloc(&TmpSet);

	return ZoneNum;
}


/*
	Save a FEVolume_c that was made from multiple GPs
*/
const int FESurface_c::SaveAsTriFEZone(const string & ZoneName, 
	vector<FieldDataType_e> DataTypes,
	const vector<ValueLocation_e> & DataLocations,
	const vector<int> & XYZVarNums)
{
	Boolean_t IsOk = IsMade();
	int ZoneNum = -1;

	if (IsOk){
		for (int i = 0; i < 3; ++i)
			DataTypes[XYZVarNums[i] - 1] = FieldDataType_Double;

// 		IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), m_XYZList.size(), m_ElemList.size(), 0, ZoneType_FETriangle, DataTypes.data());
		ArgList_pa Args = TecUtilArgListAlloc();
		TecUtilArgListAppendString(Args, SV_NAME, ZoneName.c_str());
		TecUtilArgListAppendInt(Args, SV_ZONETYPE, ZoneType_FETriangle);
		TecUtilArgListAppendInt(Args, SV_IMAX, m_XYZList.size());
		TecUtilArgListAppendInt(Args, SV_JMAX, m_ElemList.size());
		TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataTypes.data());
		TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLocations.data());
		IsOk = TecUtilDataSetAddZoneX(Args);


		TecUtilArgListDealloc(&Args);
	}

	Set_pa TmpSet = TecUtilSetAlloc(FALSE);

	if (IsOk){
		ZoneNum = TecUtilDataSetGetNumZones();
		m_ZoneNum = ZoneNum;
		TecUtilSetAddMember(TmpSet, ZoneNum, FALSE);
		TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
		TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);

		ArgList_pa CurrentArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
		TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
		TecUtilStyleSetLowLevelX(CurrentArgList);
		TecUtilArgListDealloc(&CurrentArgList);

		TecUtilSetDealloc(&TmpSet);

		vector<vector<double> > TmpValues(3, vector<double>(m_XYZList.size()));
		for (int i = 0; i < m_XYZList.size(); ++i){
			for (int j = 0; j < 3; ++j)
				TmpValues[j][i] = m_XYZList[i][j];
		}

		for (int i = 0; i < 3 && IsOk; ++i){
			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
			IsOk = VALID_REF(SetFDPtr);
			if (IsOk){
				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_XYZList.size(), reinterpret_cast<void*>(const_cast<double*>(TmpValues[i].data())));
			}
		}
// 		if (IsOk){
// 			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, RhoVarNum);
// 			IsOk = VALID_REF(SetFDPtr);
// 			if (IsOk){
// 				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_NumNodes, reinterpret_cast<void*>(const_cast<double*>(m_RhoList.data())));
// 			}
// 		}

	}

	if (IsOk){
		/*
		*	Need to define the connectivity (nodal structure) of the FE volume.
		*	The order doesn't really matter here, just that the connectivity
		*	is correct between elements.
		*/
		NodeMap_pa NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
		IsOk = VALID_REF(NodeMap);
		if (IsOk){


			for (int ei = 0; ei < m_ElemList.size(); ++ei){
				for (int i = 0; i < 3; ++i){
					TecUtilDataNodeSetByRef(NodeMap, ei + 1, i + 1, m_ElemList[ei][i] + 1);
				}
			}
// 			int ei = 1;
// 			int ni = 1;
// 			vector<int> t;
// 			//	beginning cap
// 			for (int i = 1; i < m_NumGPs - 1; ++i){
// 				t = {ni, ni + i, ni + i + 1};
// 				for (int j = 0; j < 3; ++j){
// 					TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 				}
// 				ei++;
// 			}
// 			
// 			//	triangle's up sides of gradient bundle.
// 			//	Each while iteration is one more point moving
// 			//	down each GP in the GB.
// 			//	The ti for loop is the first m_NumGPs sides
// 			//	(2 triangular elements per side) of the GB
// 			//	and then after the loop the last side.
// 			while (ni < m_NumNodes - m_NumGPs && ei < m_NumElems){
// 				for (int ti = 0; ti < m_NumGPs - 1; ++ti){
// 					t = {ni, ni + m_NumGPs, ni + m_NumGPs + 1};
// 					for (int j = 0; j < 3; ++j){
// 						TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 					}
// 					ei++;
// 					t = {ni, ni + m_NumGPs + 1, ni + 1};
// 					for (int j = 0; j < 3; ++j){
// 						TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 					}
// 					ni++;
// 					ei++;
// 				}
// 				t = {ni, ni + m_NumGPs, ni + 1};
// 				for (int j = 0; j < 3; ++j){
// 					TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 				}
// 				ei++;
// 				t = {ni, ni + 1, ni - m_NumGPs + 1};
// 				for (int j = 0; j < 3; ++j){
// 					TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 				}
// 				ni++;
// 				ei++;
// 			}
// 			
// 			//	end cap
//   // 			ni -= m_NumGPs; 
// 			//	set ni back to the last node of the 
// 			//	"first" GP in the GB.
// 			for (int i = 1; i < m_NumGPs - 1; ++i){
// 				t = {ni, ni + i, ni + i + 1};
// 				for (int j = 0; j < 3; ++j){
// 					TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 				}
// 				ei++;
// 			}
		}
	}

	TecUtilSetDealloc(&TmpSet);

	return ZoneNum;
}


const int FESurface_c::SaveAsFEZone(vector<FieldDataType_e> DataTypes,
	const vector<int> & XYZVarNums,
	const int & RhoVarNum)
{
	Boolean_t IsOk = IsMade();
	int ZoneNum = -1;

	if (IsOk){
		for (int i = 0; i < 3; ++i)
			DataTypes[XYZVarNums[i] - 1] = FieldDataType_Double;

		DataTypes[RhoVarNum - 1] = FieldDataType_Double;

		IsOk = TecUtilDataSetAddZone(CSMZoneName.FESurface.c_str(), m_XYZList.size(), m_NumElems, 0, ZoneType_FEQuad, DataTypes.data());
	}

	Set_pa TmpSet = TecUtilSetAlloc(FALSE);

	if (IsOk){
		ZoneNum = TecUtilDataSetGetNumZones();
		m_ZoneNum = ZoneNum;
		TecUtilSetAddMember(TmpSet, ZoneNum, FALSE);
		TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
		TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);

		ArgList_pa CurrentArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
		TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
		TecUtilStyleSetLowLevelX(CurrentArgList);
		TecUtilArgListDealloc(&CurrentArgList);

		TecUtilSetDealloc(&TmpSet);

		vector<vector<double> > TmpValues(3, vector<double>(m_XYZList.size()));
		for (int i = 0; i < m_XYZList.size(); ++i){
			for (int j = 0; j < 3; ++j)
				TmpValues[j][i] = m_XYZList[i][j];
		}

		for (int i = 0; i < 3 && IsOk; ++i){
			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
			IsOk = VALID_REF(SetFDPtr);
			if (IsOk){
				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_XYZList.size(), reinterpret_cast<void*>(const_cast<double*>(TmpValues[i].data())));
			}
		}
		if (IsOk){
			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, RhoVarNum);
			IsOk = VALID_REF(SetFDPtr);
			if (IsOk){
				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_XYZList.size(), reinterpret_cast<void*>(const_cast<double*>(m_RhoList.data())));
			}
		}

	}

	if (IsOk){
		/*
		*	Need to define the connectivity (nodal structure) of the FE volume.
		*	The order doesn't really matter here, just that the connectivity
		*	is correct between elements.
		*/
		NodeMap_pa NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
		IsOk = VALID_REF(NodeMap);
		if (IsOk){
			if (m_NumGPs == 3){
				int ei = 1;
				//	Beginning triangular cap
				for (int i = 1; i <= 4; ++i){
					// 4th corner gets repeated
					TecUtilDataNodeSetByRef(NodeMap, ei, i, MIN(i, 3));
				}
				ei++;
				/*
				*	Now work down the FE volume, defining the elements starting
				*	at the sphere, wrapping around the FE volume as you work
				*	towards the FE volume's end.
				*/
				for (int i = 0; i < m_NumGPPts - 1 && ei < m_NumElems - 1; ++i){
					int ii = 1 + 3 * i;
					for (int j = 0; j < 2 && ei < m_NumElems - 1; ++j){
						/*
						*	First two elements of each iteration have
						*	the same relative nodal structure, offset by
						*	one for each element.
						*/
						int Vals[4] = { ii, ii + 1, ii + 4, ii + 3 };
						for (int k = 0; k < 4; ++k){
							TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
						}
						ii++;
						ei++;
					}
					/*
					*	The third element has a different nodal structure
					*/
					int Vals[4] = { ii, ii - 2, ii + 1, ii + 3 };
					for (int k = 0; k < 4; ++k){
						TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
					}
					ei++;
				}
				//	End triangular cap
				for (int i = 1; i <= 4; ++i){
					// 4th corner gets repeated
					TecUtilDataNodeSetByRef(NodeMap, ei, i, MIN(i + 3 * (m_NumGPPts - 1), m_NumNodes));
				}
			}
			else{
				/*
				*	Same thing for the 4-sided FE volumes, but the end cap
				*	is a quad and the "walls" of the volume have three with
				*	the same nodal structure instead of two.
				*/
				int ei = 1;
				for (int i = 1; i <= 4; ++i){
					TecUtilDataNodeSetByRef(NodeMap, ei, i, MAX(i % 4, 1));
				}
				ei++;
				for (int i = 0; i < m_NumGPPts - 1 && ei < m_NumElems - 1; ++i){
					int ii = 1 + 4 * i;
					for (int j = 0; j < 3 && ei < m_NumElems - 1; ++j){
						int Vals[4] = { ii, ii + 1, ii + 5, ii + 4 };
						for (int k = 0; k < 4; ++k){
							TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
						}
						ii++;
						ei++;
					}
					int Vals[4] = { ii, ii - 3, ii + 1, ii + 4 };
					for (int k = 0; k < 4; ++k){
						TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
					}
					ei++;
				}
				//	End quadragonal cap
				for (int i = 1; i <= 4; ++i){
					// 4th corner not repeated
					TecUtilDataNodeSetByRef(NodeMap, ei, i, i + 4 * (m_NumGPPts - 1));
				}
			}
		}
	}

	TecUtilSetDealloc(&TmpSet);

	return ZoneNum;
}



/*
*	Setter methods
*/

const Boolean_t FESurface_c::Setup(const int & InZoneNum,
	const int & VolZoneNum,
	const vector<int> & InXYZVarNums,
	const vector<int> & InIntVarNums,
	const bool CopyData)
{
	m_FEVolumeMade = (InZoneNum > 0
		&& VolZoneNum > 0
		&& InZoneNum != VolZoneNum
		&& TecUtilZoneIsFiniteElement(InZoneNum)
		&& InXYZVarNums.size() == 3);

	if (m_FEVolumeMade){
		m_ZoneNum = InZoneNum;
		ZoneType_e ZoneType = TecUtilZoneGetType(InZoneNum);
		TecUtilZoneGetIJK(InZoneNum, &m_NumNodes, &m_NumElems, &m_NumNodesPerElem);
		m_FEVolumeMade = (m_NumNodes > 0 && m_NumElems > 0
			&& (ZoneType == ZoneType_FETriangle && m_NumNodesPerElem == 3)
			|| (ZoneType == ZoneType_FEQuad && m_NumNodesPerElem == 4));
	}

	if (m_FEVolumeMade){
		m_VolZoneInfo.MaxIJK.resize(3);
		TecUtilZoneGetIJK(VolZoneNum, &m_VolZoneInfo.MaxIJK[0], &m_VolZoneInfo.MaxIJK[1], &m_VolZoneInfo.MaxIJK[2]);
		ZoneXYZVarGetMinMax_Ordered3DZone(InXYZVarNums, VolZoneNum, m_VolZoneInfo.MinXYZ, m_VolZoneInfo.MaxXYZ);
		m_VolZoneInfo.DelXYZ = GetDelXYZ_Ordered3DZone(InXYZVarNums, VolZoneNum);
	}

	/*
	*	Get pointers to necessary data
	*/
	if (m_FEVolumeMade){
		m_IntVarNums = InIntVarNums;
		m_XYZPtrs.resize(3);
		m_NumIntVars = static_cast<int>(InIntVarNums.size());
		m_IntVarPtrs.resize(m_NumIntVars);

		for (int i = 0; i < 3 && m_FEVolumeMade; ++i)
			m_FEVolumeMade = m_XYZPtrs[i].GetReadPtr(InZoneNum, InXYZVarNums[i]);

		for (int i = 0; i < m_NumIntVars && m_FEVolumeMade; ++i)
			m_FEVolumeMade = m_IntVarPtrs[i].GetReadPtr(VolZoneNum, m_IntVarNums[i]);
	}

	/*
	*	Get connectivity list for FE zone
	*/
	if (m_FEVolumeMade){
		TecUtilDataNodeGetReadableRawPtr(InZoneNum, &m_ConnectivityListPtr);
		m_FEVolumeMade = (m_ConnectivityListPtr != NULL);
	}
	if (m_FEVolumeMade){
		NodeToElemMap_pa NodeToElemMap = TecUtilDataNodeToElemMapGetReadableRef(InZoneNum);
		m_ConnectivityList.resize(m_NumNodes);
		ElemFaceOffset_t FaceOffset = 0;
		LgIndex_t NumUniqueNodes, UniqueNodesSize = 0, *UniqueNodes = NULL;
		for (LgIndex_t NodeNum = 0; NodeNum < m_NumNodes; ++NodeNum){
			m_ConnectivityList[NodeNum].reserve(8);
			LgIndex_t NumElems = TecUtilDataNodeToElemMapGetNumElems(NodeToElemMap, NodeNum+1);
			for (LgIndex_t ElemNum = 1; ElemNum <= NumElems; ++ElemNum){
				LgIndex_t ElemIndex = TecUtilDataNodeToElemMapGetElem(NodeToElemMap, NodeNum+1, ElemNum);
				TecUtilDataFECellGetUniqueNodes(InZoneNum, FaceOffset, ElemIndex, &NumUniqueNodes, &UniqueNodesSize, &UniqueNodes);
				for (LgIndex_t eNodeNum = 0; eNodeNum < UniqueNodesSize; ++eNodeNum){
					if (UniqueNodes[eNodeNum] > 0 && NodeNum+1 != UniqueNodes[eNodeNum] && std::find(m_ConnectivityList[NodeNum].begin(), m_ConnectivityList[NodeNum].end(), UniqueNodes[eNodeNum]-1) == m_ConnectivityList[NodeNum].end()){
						m_ConnectivityList[NodeNum].push_back(UniqueNodes[eNodeNum]-1);
					}
				}
			}
		}
		TecUtilArrayDealloc((void**)(&UniqueNodes));
	}

	if (m_FEVolumeMade && CopyData){
		m_XYZList.resize(m_NumNodes);
		for (int i = 0; i < m_NumNodes; ++i){
			for (int j = 0; j < 3; ++j)
				m_XYZList[i][j] = m_XYZPtrs[j][i];
		}
		m_ElemList.resize(m_NumElems);
		for (int i = 0; i < m_NumElems; ++i){
			m_ElemList[i].reserve(m_NumNodesPerElem);
			for (int j = 0; j < m_NumNodesPerElem; ++j)
				m_ElemList[i].push_back(m_ConnectivityListPtr[i * m_NumNodesPerElem + j]);
		}
	}

	m_IntegrationResultsReady = FALSE;

	return m_FEVolumeMade;
}

/*
 *	Make FEVolume from an existing zone.
 *	Can be IJ ordered zone or actual FE zone.
 *	If IJ ordered zone, need to convert to FE triangle zone representation
 *		with connectivity info and such. 
 *	This assumes the structure of the surfaces made by Bondalyzer
 */
const Boolean_t FESurface_c::Setup(const int InZoneNum,
								const vector<int> & InXYZVarNums){
	REQUIRE(1 <= InZoneNum && InZoneNum <= TecUtilDataSetGetNumZones());
	REQUIRE(InXYZVarNums.size() == 3);
	for (int i : InXYZVarNums) 
		REQUIRE(1 <= i && i <= TecUtilDataSetGetNumVars());

	m_FEVolumeMade = FALSE;

	if (TecUtilZoneIsActive(InZoneNum)){
		vector<int> IJK(3);
		m_ZoneNum = InZoneNum;

		TecUtilZoneGetIJK(m_ZoneNum, &IJK[0], &IJK[1], &IJK[2]);

		if (TecUtilZoneIsOrdered(m_ZoneNum) && (1 < IJK[0] && 1 < IJK[1] && 1 == IJK[2])){
			m_XYZList.resize(IJK[0] * IJK[1]);
			m_ElemList.reserve(m_XYZList.size() * 2);

			m_XYZPtrs.resize(3);
			for (int i = 0; i < 3; ++i)
				m_XYZPtrs[i].GetReadPtr(m_ZoneNum, InXYZVarNums[i]);

			for (int j = 0; j < IJK[1]; ++j){
				for (int i = 0; i < IJK[0]; ++i){
					int Index = IJK[0] * j + i;
					for (int Dir = 0; Dir < 3; ++Dir)
						m_XYZList[Index][Dir] = m_XYZPtrs[Dir][Index];

					if (i < IJK[0] - 1 && j < IJK[1] - 1){
						vector<int> Inds;
						Inds.reserve(4);
						for (int jj = 0; jj < 2; ++jj){
							for (int ii = 0; ii < 2; ++ii){
								Inds.push_back(IJK[0] * (j + jj) + (i + ii));
							}
						}
// 						Inds[0] = Index + 1;
// 						Inds[1] = Index + 1 + 1;
// 						Inds[2] = IJK[0] * (j + 1) + i + 1;
// 						Inds[3] = Inds[2] - 1 + 1;
						m_ElemList.push_back(vector<LgIndex_t>({ Inds[0], Inds[1], Inds[3] }));
						m_ElemList.push_back(vector<LgIndex_t>({ Inds[3], Inds[2], Inds[0] }));
					}
				}
			}
			m_NumNodes = m_XYZList.size();
			m_NumElems = m_ElemList.size();
			m_FEVolumeMade = TRUE;
			return TRUE;
		}
		else{
			return Setup(InZoneNum, ZoneNumByName("Full Volume"), InXYZVarNums, vector<int>(), true);
		}
	}
	return FALSE;
}

/*
 *	Remove duplicate nodes (nodes with same position)
 *	and revise element list to reflect the removal of
 *	the nodes.
 */
void FESurface_c::RemoveDupicateNodes(){
	vector<bool> NodeIsDuplicate(m_XYZList.size(), false);
	vector<int> NodeNums(m_XYZList.size());
	bool DupFound = false;

	for (int i = 0; i < m_XYZList.size(); ++i)
		NodeNums[i] = i;

	for (int i = 0; i < m_XYZList.size() - 1; ++i){
		for (int j = i + 1; j < m_XYZList.size(); ++j){
			if (!NodeIsDuplicate[j] && sum(m_XYZList[i] == m_XYZList[j]) == 3){
				NodeIsDuplicate[j] = true;
				NodeNums[j] = NodeNums[i];
				DupFound = true;
			}
		}
	}

	if (DupFound){
		vector<vec3> NewXYZ;
		NewXYZ.reserve(m_XYZList.size());

		for (int n = 0; n < NodeIsDuplicate.size(); ++n){
			if (!NodeIsDuplicate[n]){
				NewXYZ.push_back(m_XYZList[n]);
				NodeNums[n] = NewXYZ.size() - 1;
				for (int ni = n + 1; ni < NodeIsDuplicate.size(); ++ni){
					if (NodeIsDuplicate[ni] && NodeNums[ni] == n){
						NodeNums[ni] = NodeNums[n];
					}
				}
			}
		}

		m_XYZList = NewXYZ;

		vector<vector<int> > NewElems;
		NewElems.reserve(m_ElemList.size());
		for (auto & e : m_ElemList){
			for (int & ei : e){
				ei = NodeNums[ei];
			}
			if (e[0] != e[1] && e[0] != e[2] && e[1] != e[2]){
				NewElems.push_back(e);
			}
			else{
				int a = 1;
			}
		}

		m_ElemList = NewElems;
	}


}


const Boolean_t FESurface_c::Refine(){
	if (m_RefinedXYZList.size() == 0)
		m_RefinedXYZList = m_XYZList;

	m_RefinedXYZList.reserve(m_RefinedXYZList.size() * 2);
	m_ElemList.reserve(m_ElemList.size() * 4);

	for (int i = 0; i < m_NumElems; ++i)
		TriangleEdgeMidPointSubdivide(i);

	m_XYZList = m_RefinedXYZList;
	m_NumNodes = m_RefinedXYZList.size();
	m_NumElems = m_ElemList.size();

	RemoveDupicateNodes();

	m_RefinedXYZList = m_XYZList;
	m_NumNodes = m_RefinedXYZList.size();
	m_NumElems = m_ElemList.size();

	return TRUE;
}


const Boolean_t FESurface_c::DoIntegration(const Boolean_t & IntegrateVolume,
	const vector<vec> & stuW,
	const vector<int> & SplitPtNums,
	const vector<vec> & stuW2)
{
	Boolean_t IsOk = m_FEVolumeMade && m_NumNodesPerElem == 4;

	if (IsOk){
		Boolean_t IsBond = m_GPList.size() > 3;
		Boolean_t IsCage;
		if (IsBond){

		}
		else{
			/*
			*	Check to see if gradient bundle ends at a cage point.
			*	First check area of triandle made by last three points
			*	(should be 0 if the come together at a cage point).
			*/
		}


		vector<mat> xyzW = GetIntegrationPointsWeights(stuW);

// 		string TempStr = "Negative weights at i = {", TempStr2 = TempStr;
// 		for (int i = 0; i < xyzW[0].n_rows; ++i){
// 			if (xyzW[1](i,0) < 0){
// 				TempStr += to_string(i) + ",";
// 			}
// 		}
// 
// 		if (TempStr != TempStr2){
// 			TecUtilDialogMessageBox(TempStr.c_str(), MessageBoxType_Error);
// 			TempStr = "";
// 		}

// 		vector<vector<double> > vxyzW(xyzW[0].n_rows, vector<double>(3));
// 		vector<double> vW(xyzW[1].n_rows);
// 		for (int i = 0; i < 3; ++i){
// 			for (int j = 0; j < xyzW[0].n_rows; ++j){
// 				vxyzW[j][i] = xyzW[0](j, i);
// 			}
// 		}
// 		for (int i = 0; i < vW.size(); ++i){
// 			vW[i] = xyzW[1](i, 0);
// 		}

		m_IntValues.resize(m_NumIntVars + int(IntegrateVolume), 0.0);

		for (int i = 0; i < xyzW[0].n_rows; ++i){
			vec3 TmpVec = xyzW[0].submat(i,0,i,2).t();
			SetIndexAndWeightsForPoint(TmpVec, m_VolZoneInfo);
			if (IntegrateVolume){
				m_IntValues[m_NumIntVars] += xyzW[1](i, 0);
			}
			for (int j = 0; j < m_NumIntVars; ++j){
				m_IntValues[j] += xyzW[1](i, 0) * ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[j]);
			}
		}

		m_IntegrationResultsReady = TRUE;
	}

	return IsOk;
}

// const Boolean_t FEVolume_c::DoIntegration(const int & ResolutionScale, const Boolean_t & IntegrateVolume)
// {
// 	Boolean_t IsOk = m_FEVolumeMade;
// 
// 	if (IsOk && !m_IntegrationResultsReady){
// 		m_ZoneMinXYZ = m_VolZoneInfo.MaxXYZ;
// 		m_ZoneMaxXYZ = m_VolZoneInfo.MinXYZ;
// 		/*
// 		*	Get FE zone's max and min XYZ values
// 		*/
// 		for (int i = 0; i < 3; ++i){
// 			for (int j = 0; j < m_NumNodes; ++j){
// 				if (m_XYZPtrs[i][j] < m_ZoneMinXYZ[i])
// 					m_ZoneMinXYZ[i] = m_XYZPtrs[i][j];
// 				else if (m_XYZPtrs[i][j] > m_ZoneMaxXYZ[i])
// 					m_ZoneMaxXYZ[i] = m_XYZPtrs[i][j];
// 			}
// 		}
// 
// 		IsOk = sum(m_ZoneMaxXYZ > m_ZoneMinXYZ) == 3;
// 	}
// 
// 	/*
// 	if (IsOk){
// 	/ *
// 	*	Loop over every point in the system checking if it's interior
// 	*	or not. If interior, add it to the integral value.
// 	* /
// 	vec3 DelXYZ, TmpPoint;
// 
// 	for (int i = 0; i < 3; ++i){
// 	DelXYZ[i] = (m_VolZoneInfo.MaxXYZ[i] - m_VolZoneInfo.MinXYZ[i]) / static_cast<double>(m_VolZoneInfo.MaxIJK[i] - 1);
// 	}
// 	m_IntValues.resize(m_NumIntVars, 0.0);
// 	int IJK[3];
// 	for (IJK[2] = 0; IJK[2] < m_VolZoneInfo.MaxIJK[2]; ++IJK[2]){
// 	for (IJK[1] = 0; IJK[1] < m_VolZoneInfo.MaxIJK[1]; ++IJK[1]){
// 	for (IJK[0] = 0; IJK[0] < m_VolZoneInfo.MaxIJK[0]; ++IJK[0]){
// 	for (int i = 0; i < 3; ++i){
// 	TmpPoint[i] = m_VolZoneInfo.MinXYZ[i] + static_cast<double>(IJK[i]) * DelXYZ[i];
// 	}
// 	if (PointIsInterior(TmpPoint)){
// 	SetIndexAndWeightsForPoint(TmpPoint, m_VolZoneInfo);
// 	for (int i = 0; i < m_NumIntVars; ++i){
// 	if (IntVarFDTypes[i] == FieldDataType_Float)
// 	m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, IntVarPtrsFlt[i]);
// 	else
// 	m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, IntVarPtrsDbl[i]);
// 	}
// 	}
// 	}
// 	}
// 	}
// 	}
// 	*/
// 
// 
// 	if (IsOk){
// 		/*
// 		*	Loop over each cell in the full system volume.
// 		*	Maintain list of bools to see which nodes are in/out
// 		*	of the FE volume.
// 		*	If whole cell is in FE volume, add whole cell to integral.
// 		*	If part of cell in FE volume, subsample cell and add pieces
// 		*	of it to integral if they're in FE volume.
// 		*/
// 		vec3 FarPoint = m_VolZoneInfo.MinXYZ - 1.;
// 		
// 		vec3 DelXYZ = m_VolZoneInfo.DelXYZ, SubDelXYZ, TmpPoint, TmpSubPoint;
// 		double CellVolume = 1.0;
// 
// 		for (int i = 0; i < 3; ++i){
// // 			DelXYZ[i] = (m_VolZoneInfo.MaxXYZ[i] - m_VolZoneInfo.MinXYZ[i]) / static_cast<double>(m_VolZoneInfo.MaxIJK[i] - 1);
// 			CellVolume *= DelXYZ[i];
// 		}
// 
// 		SubDelXYZ = DelXYZ / static_cast<double>(ResolutionScale);
// 		double SubCellFactor = 1.0 / static_cast<double>(ResolutionScale * ResolutionScale * ResolutionScale);
// 		double SubCellVolume = CellVolume * SubCellFactor;
// 
// 		if (IntegrateVolume)
// 			m_IntValues = vector<double>(m_NumIntVars + 1, 0.0);
// 		else
// 			m_IntValues = vector<double>(m_NumIntVars, 0.0);
// 
// 		vector<Boolean_t> VolNodeInFEVolume(m_VolZoneInfo.MaxIJK[0] * m_VolZoneInfo.MaxIJK[1] * m_VolZoneInfo.MaxIJK[2]);
// 		vector<Boolean_t> VolNodeChecked(m_VolZoneInfo.MaxIJK[0] * m_VolZoneInfo.MaxIJK[1] * m_VolZoneInfo.MaxIJK[2], FALSE);
// 		int IJK[3], TmpIJK[3];
// 		int Index[8];
// 		int NumInteriorNodes;
// 		for (IJK[2] = 1; IJK[2] < m_VolZoneInfo.MaxIJK[2]; ++IJK[2]){
// 			for (IJK[1] = 1; IJK[1] < m_VolZoneInfo.MaxIJK[1]; ++IJK[1]){
// 				for (IJK[0] = 1; IJK[0] < m_VolZoneInfo.MaxIJK[0]; ++IJK[0]){
// 					/*
// 					*	Get I-D index numbers for the nodes of the cell.
// 					*/
// 					Index[0] = IndexFromIJK(IJK[0], IJK[1], IJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[1] = IndexFromIJK(IJK[0] + 1, IJK[1], IJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[2] = IndexFromIJK(IJK[0] + 1, IJK[1] + 1, IJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[3] = IndexFromIJK(IJK[0], IJK[1] + 1, IJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[4] = IndexFromIJK(IJK[0], IJK[1], IJK[2] + 1, m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[5] = IndexFromIJK(IJK[0] + 1, IJK[1], IJK[2] + 1, m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[6] = IndexFromIJK(IJK[0] + 1, IJK[1] + 1, IJK[2] + 1, m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[7] = IndexFromIJK(IJK[0], IJK[1] + 1, IJK[2] + 1, m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 
// 					// 					if (IJK[0] == 40 && IJK[1] == 40 && IJK[2] == 40){
// 					// 						int a = 1;
// 					// 					}
// 
// 
// 					/*
// 					*	Check the number of cell nodes that are interior to the volume.
// 					*/
// 					NumInteriorNodes = 0;
// 
// 					for (int i = 0; i < 8; ++i){
// 						if (!VolNodeChecked[Index[i]]){
// 							/*
// 							*	Node has not been check as interior or not, so get the XYZ
// 							*	coordinates and check if interior.
// 							*/
// 							// Get IJK for corner of cell (subtract 1 since IJK were 1-based
// 							for (int j = 0; j < 3; ++j)
// 								TmpIJK[j] = IJK[j] - 1;
// 
// 							GetCellCornerIndices(i, TmpIJK[0], TmpIJK[1], TmpIJK[2]);
// 
// 							for (int j = 0; j < 3; ++j)
// 								TmpPoint[j] = m_VolZoneInfo.MinXYZ[j] + static_cast<double>(TmpIJK[j]) * DelXYZ[j];
// 
// 							//Boolean_t IsInside = VolNodeInFEVolume[Index[i]] = PointIsInterior(TmpPoint);
// 							VolNodeInFEVolume[Index[i]] = PointIsInterior(TmpPoint, FarPoint);
// 							VolNodeChecked[Index[i]] = TRUE;
// 						}
// 
// 						if (VolNodeInFEVolume[Index[i]]) NumInteriorNodes++;
// 					}
// 
// 					if (NumInteriorNodes == 8){
// 						/*
// 						*	Entire cell is in FE volume, so get midpoint of the cell and add that
// 						*	value to the integral.
// 						*/
// 
// 						/*
// 						*	Get the midpoint of the cell by taking the 0,0,0 corner and adding 0.5*dXYZ to it.
// 						*/
// 						for (int i = 0; i < 3; ++i)
// 							TmpPoint[i] = m_VolZoneInfo.MinXYZ[i] + (static_cast<double>(IJK[i]) - 0.5) * DelXYZ[i];
// 
// 						/*
// 						*	Add the value of each integrating variable for the cell midpoint to the
// 						*	integral value for that variable.
// 						*	This assumes that the value is already in the correct units for the system.
// 						*/
// 						if (SetIndexAndWeightsForPoint(TmpPoint, m_VolZoneInfo)){
// 							if (IntegrateVolume)
// 								m_IntValues[m_NumIntVars] += CellVolume;
// 							for (int i = 0; i < m_NumIntVars; ++i){
// 								m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * CellVolume;
// 							}
// 						}
// 						else
// 							TecUtilDialogErrMsg("Error: Indices out of bounds");
// 
// 					}
// 					else if (NumInteriorNodes > 0){
// 						/*
// 						*	Cell is partially contained in FE volume, so need to
// 						*	subsample the cell, adding values of subcells to integral
// 						*	as they're found.
// 						*
// 						*	Since we're testing points INSIDE the cell, need to start at 1/2 of
// 						*	the sub-resolution into the cell so that the subcells of density we add or don't add
// 						*	are entirely contained in the cell and span the cell.
// 						*/
// 						for (int i = 0; i < 3; ++i)
// 							TmpPoint[i] = m_VolZoneInfo.MinXYZ[i] + static_cast<double>(IJK[i] - 1) * DelXYZ[i] + (SubDelXYZ[i] / 2.0);
// 
// 						for (TmpIJK[0] = 0; TmpIJK[0] < ResolutionScale; ++TmpIJK[0]){
// 							for (TmpIJK[1] = 0; TmpIJK[1] < ResolutionScale; ++TmpIJK[1]){
// 								for (TmpIJK[2] = 0; TmpIJK[2] < ResolutionScale; ++TmpIJK[2]){
// 									for (int i = 0; i < 3; ++i)
// 										TmpSubPoint[i] = TmpPoint[i] + SubDelXYZ[i] * static_cast<double>(TmpIJK[i]);
// 
// 									if (PointIsInterior(TmpSubPoint, FarPoint)){
// 										if (SetIndexAndWeightsForPoint(TmpSubPoint, m_VolZoneInfo)){
// 											if (IntegrateVolume)
// 												m_IntValues[m_NumIntVars] += SubCellVolume;
// 											for (int i = 0; i < m_NumIntVars; ++i){
// 												m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * SubCellVolume;
// 											}
// 										}
// 										else
// 											TecUtilDialogErrMsg("Error: Indices out of bounds");
// 									}
// 								}
// 							}
// 						}
// 					}
// 				}
// 
// 			}
// 		}
// 
// 		m_IntegrationResultsReady = TRUE;
// 	}
// 
// 	return IsOk;
// }

const Boolean_t FESurface_c::DoIntegrationNew(const int & ResolutionScale, const Boolean_t & IntegrateVolume)
{
	Boolean_t IsOk = m_FEVolumeMade;

	if (IsOk && !m_IntegrationResultsReady){
		m_ZoneMinXYZ = m_VolZoneInfo.MaxXYZ;
		m_ZoneMaxXYZ = m_VolZoneInfo.MinXYZ;
		/*
		*	Get FE zone's max and min XYZ values
		*/
		for (int i = 0; i < 3; ++i){
			for (int j = 0; j < m_NumNodes; ++j){
				if (m_XYZPtrs[i][j] < m_ZoneMinXYZ[i])
					m_ZoneMinXYZ[i] = m_XYZPtrs[i][j];
				else if (m_XYZPtrs[i][j] > m_ZoneMaxXYZ[i])
					m_ZoneMaxXYZ[i] = m_XYZPtrs[i][j];
			}
		}

		IsOk = sum(m_ZoneMaxXYZ > m_ZoneMinXYZ) == 3;
	}

	vector<int> ZoneMinIJK, ZoneMaxIJK, ZoneNumIJK(3);
	int ZoneNumPts = 1;

	/*
	*	Get the IJK indices that correspond to the min and max XYZ for the zone
	*/
	if (IsOk){
		ZoneMinIJK = GetIJKForPoint(m_ZoneMinXYZ, m_VolZoneInfo);
		ZoneMaxIJK = GetIJKForPoint(m_ZoneMaxXYZ, m_VolZoneInfo);
		for (int i = 0; i < 3; ++i){
			for (int j = 0; j < 2; ++j){
// 				if (ZoneMinIJK[i] > 1)
// 					ZoneMinIJK[i]--;
				if (ZoneMaxIJK[i] < m_VolZoneInfo.MaxIJK[i])
					ZoneMaxIJK[i]++;
			}
			ZoneNumIJK[i] = ZoneMaxIJK[i] - ZoneMinIJK[i] + 1;
			ZoneNumPts *= ZoneNumIJK[i];
		}
	}



	vector<Boolean_t> DoCheckNode(ZoneNumPts,FALSE), NodeIsInterior(ZoneNumPts, TRUE);
	vector<double> CheckNodeMinDistSqrToSurfaceNode(ZoneNumPts, 0.0);
	vector<int> CheckNodeMinDistNodeNum(ZoneNumPts, -1);

	m_IntValues = vector<double>(m_NumIntVars + static_cast<int>(IntegrateVolume), 0.0);

// 	if (IntegrateVolume)
// 		m_IntValues = vector<double>(m_NumIntVars + 1, 0.0);
// 	else
// 		m_IntValues = vector<double>(m_NumIntVars, 0.0);

	if (IsOk){
		CalcMaxNodeDistSqr();
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
		SaveAsTriFEZone(vector<int>({ 1, 2, 3 }));
// 		return IsOk;
#endif // _DEBUG

		vector<int> VolIJK(3), ZoneIJK(3);
		vec3 VolPt, NodePt, VolSubPt;
		vector<vec3> FarPoints(8), MinMaxXYZ = { m_VolZoneInfo.MinXYZ - 1, m_VolZoneInfo.MaxXYZ + 1 };
		for (int i = 0; i < 2; ++i){
			double x = MinMaxXYZ[i][0];
			for (int j = 0; j < 2; ++j){
				double y = MinMaxXYZ[j][1];
				for (int k = 0; k < 2; ++k){
					double z = MinMaxXYZ[k][2];
					FarPoints[4 * i + 2 * j + k] << x << y << z;
				}
			}
		}
		LgIndex_t VolIndex, ZoneIndex;
		ZoneIJK[2] = 1;
		for (VolIJK[2] = ZoneMinIJK[2]; VolIJK[2] <= ZoneMaxIJK[2]; ++VolIJK[2]){
			ZoneIJK[1] = 1;
			for (VolIJK[1] = ZoneMinIJK[1]; VolIJK[1] <= ZoneMaxIJK[1]; ++VolIJK[1]){
				ZoneIJK[0] = 1;
				for (VolIJK[0] = ZoneMinIJK[0]; VolIJK[0] <= ZoneMaxIJK[0]; ++VolIJK[0]){
// 					VolIndex = IndexFromIJK(VolIJK[0], VolIJK[1], VolIJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]);
					ZoneIndex = IndexFromIJK(ZoneIJK[0], ZoneIJK[1], ZoneIJK[2], ZoneNumIJK[0], ZoneNumIJK[1]) - 1;
					for (int ii = 0; ii < 3; ++ii){
						VolPt[ii] = m_VolZoneInfo.MinXYZ[ii] + static_cast<double>(VolIJK[ii]-1) * m_VolZoneInfo.DelXYZ[ii];
						NodeIsInterior[ZoneIndex] = NodeIsInterior[ZoneIndex] && (ZoneIJK[ii] > 1 && ZoneIJK[ii] < ZoneNumIJK[ii]);
					}
					if (ZoneIJK[0] == 11 && ZoneIJK[1] == 6 && ZoneIJK[2] == 6){
						int a = 1;
					}
					NodeIsInterior[ZoneIndex] = NodeIsInterior[ZoneIndex] && PointIsInterior(VolPt, FarPoints);
					if (!NodeIsInterior[ZoneIndex]){
						DoCheckNode[ZoneIndex] = DistSqrToSurfaceNodeWithinTolerance(VolPt, CheckNodeMinDistSqrToSurfaceNode[ZoneIndex], CheckNodeMinDistNodeNum[ZoneIndex]);
						/*for (int n = 0; n < m_NumNodes && !DoCheckNode[ZoneIndex]; ++n){
							for (int ii = 0; ii < 3; ++ii){
								NodePt[ii] = m_XYZPtrs[ii][n];
							}
							DoCheckNode[ZoneIndex] = (DistSqr(VolPt, NodePt) < m_MaxNeighborNodeDistSqr[n]);
							//DoCheckNode[ZoneIndex] = (DistSqr(VolPt, NodePt) < m_MaxNodeDistanceSqr);
						}*/
					}
					ZoneIJK[0]++;
				}
				ZoneIJK[1]++;
			}
			ZoneIJK[2]++;
		}


#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
		int NumCheckPoints = 0;
		for (const auto & i : DoCheckNode)
			NumCheckPoints += int(i);
		m_InteriorPts.reserve(NumCheckPoints);
		m_CheckPts.reserve(NumCheckPoints);
		m_InteriorSubPts.reserve(NumCheckPoints * (2 ^ ResolutionScale));
		m_CheckSubPts.reserve(NumCheckPoints * (2 ^ ResolutionScale));
#endif


		vec3 SubDelXYZ, TmpPoint, TmpSubPoint;
		double CellVolume = 1.0;

		for (int i = 0; i < 3; ++i){
			// 			DelXYZ[i] = (m_VolZoneInfo.MaxXYZ[i] - m_VolZoneInfo.MinXYZ[i]) / static_cast<double>(m_VolZoneInfo.MaxIJK[i] - 1);
			CellVolume *= m_VolZoneInfo.DelXYZ[i];
		}

		SubDelXYZ = m_VolZoneInfo.DelXYZ / static_cast<double>(ResolutionScale);
		double SubDivideFactorUser = 1.0 / static_cast<double>(ResolutionScale * ResolutionScale * ResolutionScale);
		double SubCellVolumeUser = CellVolume * SubDivideFactorUser;
		double SubDivideFactor, SubCellVolume;
		vector<LgIndex_t> NeighborIJK(3);
		vector<Boolean_t> NodeIsInteriorCopy = NodeIsInterior;
		LgIndex_t CheckInt;
		ZoneIJK[2] = 1;
		for (VolIJK[2] = ZoneMinIJK[2]; VolIJK[2] <= ZoneMaxIJK[2]; ++VolIJK[2]){
		//for (VolIJK[2] = ZoneMinIJK[2]; VolIJK[2] <= ZoneMaxIJK[2] - (ZoneNumIJK[2] * 3 / 4); ++VolIJK[2]){
			ZoneIJK[1] = 1;
			for (VolIJK[1] = ZoneMinIJK[1]; VolIJK[1] <= ZoneMaxIJK[1]; ++VolIJK[1]){
				ZoneIJK[0] = 1;
				for (VolIJK[0] = ZoneMinIJK[0]; VolIJK[0] <= ZoneMaxIJK[0]; ++VolIJK[0]){
					ZoneIndex = IndexFromIJK(ZoneIJK[0], ZoneIJK[1], ZoneIJK[2], ZoneNumIJK[0], ZoneNumIJK[1]) - 1; 
					for (int ii = 0; ii < 3; ++ii){
						VolPt[ii] = m_VolZoneInfo.MinXYZ[ii] + static_cast<double>(VolIJK[ii] - 1) * m_VolZoneInfo.DelXYZ[ii];
					}
					if (ZoneIJK[0] == 11 && ZoneIJK[1] == 6 && ZoneIJK[2] == 6){
						int a = 1;
					}
					if (NodeIsInterior[ZoneIndex]){
						NeighborIJK = ZoneIJK;
						for (NeighborIJK[0] = ZoneIJK[0] - 1; NeighborIJK[0] < ZoneIJK[0] + 2; ++NeighborIJK[0]){
							for (NeighborIJK[1] = ZoneIJK[1] - 1; NeighborIJK[1] < ZoneIJK[1] + 2; ++NeighborIJK[1]){
								for (NeighborIJK[2] = ZoneIJK[2] - 1; NodeIsInterior[ZoneIndex] && NeighborIJK[2] < ZoneIJK[2] + 2; ++NeighborIJK[2]){
									NodeIsInterior[ZoneIndex] = NodeIsInteriorCopy[IndexFromIJK(NeighborIJK[0], NeighborIJK[1], NeighborIJK[2], ZoneNumIJK[0], ZoneNumIJK[1]) - 1];
								}
							}

						}
// 						for (int iDir = 0; iDir < 3; ++iDir){
// 							for (int d = -1; d < 2 && NodeIsInterior[ZoneIndex]; d += 2){
// 								CheckInt = NeighborIJK[iDir] + d;
// 								NodeIsInterior[ZoneIndex] = (CheckInt >= 1 && CheckInt <= ZoneNumIJK[iDir]);
// 								if (NodeIsInterior[ZoneIndex]){
// 									NeighborIJK[iDir] = CheckInt;
// 									NodeIsInterior[ZoneIndex] = NodeIsInteriorCopy[IndexFromIJK(NeighborIJK[0], NeighborIJK[1], NeighborIJK[2], ZoneNumIJK[0], ZoneNumIJK[1]) - 1];
// 									NeighborIJK[iDir] -= d;
// 								}
// 							}
// 						}
						DoCheckNode[ZoneIndex] = !NodeIsInterior[ZoneIndex];
						if (NodeIsInterior[ZoneIndex]){
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
							m_InteriorPts.push_back(VolPt);
#endif
							VolIndex = IndexFromIJK(VolIJK[0], VolIJK[1], VolIJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
							SubDivideFactor = 1.;
							for (int dir = 0; dir < 3; ++dir){
								if (VolIJK[dir] == 1 || VolIJK[dir] == m_VolZoneInfo.MaxIJK[dir])
									SubDivideFactor *= 0.5;
							}
							SubCellVolume = CellVolume * SubDivideFactor;
							if (IntegrateVolume)
								m_IntValues[m_NumIntVars] += SubCellVolume;
							for (int i = 0; i < m_NumIntVars; ++i){
								m_IntValues[i] += m_IntVarPtrs[i][VolIndex] * SubCellVolume;
							}
						}
						else{
							DistSqrToSurfaceNodeWithinTolerance(VolPt, CheckNodeMinDistSqrToSurfaceNode[ZoneIndex], CheckNodeMinDistNodeNum[ZoneIndex]);
						}
					}
					if (!NodeIsInterior[ZoneIndex] && DoCheckNode[ZoneIndex]){
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
						m_CheckPts.push_back(VolPt);
#endif
						SubDivideIntegrateCellAtPoint(VolPt, FarPoints, m_VolZoneInfo.DelXYZ * 0.5, CheckNodeMinDistSqrToSurfaceNode[ZoneIndex], CheckNodeMinDistNodeNum[ZoneIndex], ResolutionScale - 1, IntegrateVolume);
// 						for (int ii = 0; ii < 3; ++ii){
// 							VolPt[ii] = m_VolZoneInfo.MinXYZ[ii] + static_cast<double>(VolIJK[ii]-1) * m_VolZoneInfo.DelXYZ[ii];
// 							VolPt[ii] = MAX(m_VolZoneInfo.MinXYZ[ii], MIN(m_VolZoneInfo.MaxXYZ[ii], VolPt[ii] - (m_VolZoneInfo.DelXYZ[ii] / 2.) + (SubDelXYZ[ii] / 2.)));
// 						}
// 						VolSubPt[0] = VolPt[0];
// 						for (int i = 0; i < ResolutionScale; ++i){
// 							VolSubPt[1] = VolPt[1];
// 							for (int j = 0; j < ResolutionScale; ++j){
// 								VolSubPt[2] = VolPt[2];
// 								for (int k = 0; k < ResolutionScale; ++k){
// #ifdef RECORD_INT_POINTS
// 									InteriorCheckPts.push_back(VolSubPt);
// #endif
// 									if (PointIsInterior(VolSubPt, FarPoint) && SetIndexAndWeightsForPoint(VolSubPt, m_VolZoneInfo)){
// #ifdef RECORD_INT_POINTS
// 										InteriorSubPts.push_back(VolSubPt);
// #endif
// 										if (IntegrateVolume)
// 											m_IntValues[m_NumIntVars] += SubCellVolumeUser;
// 										for (int i = 0; i < m_NumIntVars; ++i){
// 											m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * SubCellVolumeUser;
// 										}
// 									}
// 									VolSubPt[2] += SubDelXYZ[2];
// 								}
// 								VolSubPt[1] += SubDelXYZ[1];
// 							}
// 							VolSubPt[0] += SubDelXYZ[0];
// 						}
					}
					ZoneIJK[0]++;
				}
				ZoneIJK[1]++;
			}
			ZoneIJK[2]++;
		}
		m_IntegrationResultsReady = TRUE;
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
		vector<string> ZoneNames = { "Interior Points", "CheckPoints", "Interior Subpoints", "Sub Check Points" };
		for (auto & i : ZoneNames) i += " zone " + to_string(m_ZoneNum);
		vector<ColorIndex_t> Colors = { Green_C, Blue_C, Purple_C, Red_C };
		vector<vector<vec3> > Vecs = { m_InteriorPts, m_CheckPts, m_InteriorSubPts, m_CheckSubPts };
		vector<EntIndex_t> VarNums(3);
		for (int i = 0; i < 3; ++i)
			VarNums[i] = m_XYZPtrs[i].VarNum();
		for (int i = 0; i < Vecs.size(); ++i)
			SaveVec3VecAsScatterZone(Vecs[i], ZoneNames[i], Colors[i], VarNums);
#endif
	}

	m_ElemList = vector<vector<LgIndex_t> >();
	m_RefinedXYZList = vector<vec3>();
	
	return IsOk;
}



const Boolean_t FESurface_c::CalcMaxNodeDistSqr(){
	Boolean_t IsOk = m_XYZPtrs.size() == 3;

	for (int i = 0; i < m_XYZPtrs.size() && IsOk; ++i){
		IsOk = m_XYZPtrs[i].IsReady();
	}
	
	m_XYZList.resize(m_NumNodes);
	for (int n = 0; n < m_NumNodes; ++n){
		for (int i = 0; i < 3; ++i){
			m_XYZList[n][i] = m_XYZPtrs[i][n];
		}
	}

// 	if (IsOk){
// 		m_MaxNodeDistanceSqr = 0.0;
// 		m_MinNodeDistanceSqr = 1e150;
// 		m_MaxNeighborNodeDistSqr.resize(m_NumNodes, 0.0);
// 		m_MinNeighborNodeDistSqr.resize(m_NumNodes, 1e100);
// 		double TempDist;
// 		for (int i = 0; i < m_NumNodes; ++i){
// 			for (const auto & j : m_ConnectivityList[i]){
// 				TempDist = DistSqr(m_XYZList[i], m_XYZList[j]);
// 				m_MaxNeighborNodeDistSqr[i] = MAX(m_MaxNeighborNodeDistSqr[i], TempDist);
// 				if (TempDist > 1e-5){
// 					m_MinNeighborNodeDistSqr[i] = MIN(m_MinNeighborNodeDistSqr[i], TempDist);
// 				}
// 			}
// 			m_MaxNodeDistanceSqr = MAX(m_MaxNodeDistanceSqr, m_MaxNeighborNodeDistSqr[i]);
// 			m_MinNodeDistanceSqr = MIN(m_MinNodeDistanceSqr, m_MinNeighborNodeDistSqr[i]);
// 		}
// 	}
// 	
	m_MinNodeDistanceSqr = sum(square(m_VolZoneInfo.DelXYZ));
	
	m_RefinedXYZList = m_XYZList;
	m_ElemList.resize(m_NumElems, vector<LgIndex_t>(3, -1));
	
	for (int i = 0; i < m_NumElems; ++i){
		for (int j = 0; j < 3; ++j){
			m_ElemList[i][j] = m_ConnectivityListPtr[i * m_NumNodesPerElem + j];
		}
	}

	vector<int> ElemNumList(m_NumElems);
	for (int i = 0; i < m_NumElems; ++i){
		ElemNumList[i] = i;
	}
	RefineTriElems(ElemNumList);

	m_MinNodeDistanceSqr *= 0.49; // 70% of distance (70%^2 of square distance)

	return IsOk;
}

const vector<int> FESurface_c::TriangleEdgeMidPointSubdivide(const int & TriNum){
	vector<LgIndex_t> NewNodes(3);
	vector<int> NewTris(4);
	NewTris[0] = TriNum;
	int ni = m_RefinedXYZList.size();
	int ti = m_ElemList.size();
	for (int i = 0; i < 3; ++i){
		m_RefinedXYZList.push_back(
			(
				m_RefinedXYZList[m_ElemList[TriNum][i]]
				+ m_RefinedXYZList[m_ElemList[TriNum][(i + 1) % 3]]
			) / 2.
		);
		NewNodes[i] = ni++;
		NewTris[i+1] = ti++;
	}
	m_ElemList.push_back(vector<LgIndex_t>({m_ElemList[TriNum][1], NewNodes[1], NewNodes[0]}));
	m_ElemList.push_back(vector<LgIndex_t>({m_ElemList[TriNum][2], NewNodes[2], NewNodes[1]}));
	m_ElemList.push_back(NewNodes);
	m_ElemList[TriNum] = vector<LgIndex_t>({m_ElemList[TriNum][0], NewNodes[0], NewNodes[2]});
	return NewTris;
}

void FESurface_c::RefineTriElems(const vector<int> & TriNumList){
	for (const auto & t : TriNumList){
		double MaxNodeDistSqr = 0.0;
		for (int i = 0; i < 3; ++i){
			MaxNodeDistSqr = MAX(MaxNodeDistSqr, DistSqr(m_RefinedXYZList[m_ElemList[t][i]], m_RefinedXYZList[m_ElemList[t][(i+1) % 3]]));
		}
		if (MaxNodeDistSqr > m_MinNodeDistanceSqr){
			RefineTriElems(TriangleEdgeMidPointSubdivide(t));
		}
	}
}


void FESurface_c::TriPolyLines(const bool ConnectBeginningAndEndGPs)
{
	vec3 NewNode, EdgeVec;
	for (int iGP = 0; iGP < m_NumGPs - (ConnectBeginningAndEndGPs ? 0 : 1); ++iGP){
		int lInd = iGP,
			rInd = (iGP + 1) % m_NumGPs;
		int li = lInd * m_NumGPPts,
			ri = rInd * m_NumGPPts;
		vector<int> lr = { li, ri }, 
			lrMid = { li,ri };
		vector<bool> MidFound(2, false);
		int NumEdgePts;
		double EdgePtSpacing, EdgeLen, MinNodeScore, TmpNodeScore, MinNodeNum, lLen, rLen, TmpLen;
		EdgeVec = m_XYZList[ri + m_NumGPPts - 1] - m_XYZList[li + m_NumGPPts - 1];
		EdgeLen = norm(EdgeVec);
		EdgePtSpacing = norm(m_XYZList[li] - m_XYZList[li + 1]);
		NumEdgePts = static_cast<int>(EdgeLen / EdgePtSpacing);
		bool HasFarEdge = (NumEdgePts > 2);
		if (HasFarEdge){
			EdgePtSpacing = EdgeLen / static_cast<double>(NumEdgePts - 1);
			EdgeVec = normalise(EdgeVec) * EdgePtSpacing;
			MinNodeScore = DBL_MAX;
			for (int i = 1; i < NumEdgePts - 1; ++i){
				NewNode = m_XYZList[li + m_NumGPPts - 1] + EdgeVec * static_cast<double>(i);
				lLen = DBL_MAX;
				rLen = DBL_MAX;
				for (int j = 0; j < m_NumGPPts; ++j){
					lLen = MIN(lLen, DistSqr(NewNode, m_XYZList[li + j]));
					rLen = MIN(rLen, DistSqr(NewNode, m_XYZList[ri + j]));

				}
				TmpNodeScore = lLen + rLen;
				if (TmpNodeScore < MinNodeScore){
					MinNodeScore = TmpNodeScore;
					MinNodeNum = i;
				}
			}
			NewNode = m_XYZList[li + m_NumGPPts - 1] + EdgeVec * static_cast<double>(MinNodeNum);
		}
		else
			NewNode = (m_XYZList[li + m_NumGPPts - 1] + m_XYZList[ri + m_NumGPPts - 1]) / 2.0;
		if (HasFarEdge) m_XYZList.push_back(NewNode);
		int NewNodeNum = m_XYZList.size() - 1;
		while (lr[0] < (lInd + 1) * m_NumGPPts - 1 && lr[1] < (rInd + 1) * m_NumGPPts - 1){
			int FarPoint, MinSide, MinFarPoint;
			double MinLen = DBL_MAX, TmpLen;
			for (int i = 0; i < 2; ++i){
				for (int j = 0; j < 1 + int(HasFarEdge); ++j){
					if (j == 0)
						FarPoint = lr[(i + 1) % 2];
					else
						FarPoint = NewNodeNum;
					TmpLen = DistSqr(m_XYZList[lr[i] + 1], m_XYZList[FarPoint]);
					if (TmpLen < MinLen){
						MinSide = i;
						MinFarPoint = FarPoint;
						MinLen = TmpLen;
					}
				}
			}
			if (HasFarEdge && !MidFound[MinSide] && MinFarPoint == NewNodeNum){
				MidFound[MinSide] = true;
				lrMid[MinSide] = lr[MinSide];
			}
			m_ElemList.push_back({ lr[MinSide], lr[MinSide] + 1, MinFarPoint });
			lr[MinSide]++;
		}
		vector<int> lrInd = { lInd, rInd };
		for (int i = 0; i < 2; ++i){
			while (lr[i] < (lrInd[i] + 1) * m_NumGPPts - 1){
				if (DistSqr(m_XYZList[lr[i] + 1], m_XYZList[lr[(i + 1) % 2]]) < DistSqr(m_XYZList[lr[i] + 1], m_XYZList[NewNodeNum])){
					m_ElemList.push_back({ lr[i], lr[i] + 1, lr[(i + 1) % 2] });
				}
				else if (HasFarEdge){
					m_ElemList.push_back({ lr[i], lr[i] + 1, NewNodeNum });
					if (!MidFound[i]){
						MidFound[i] = true;
						lrMid[i] = lr[i];
					}
				}
				lr[i]++;
			}
		}
		if (MidFound[0] && MidFound[1]){
			m_ElemList.push_back({ lrMid[0], lrMid[1], NewNodeNum });
		}
	}
}

const Boolean_t FESurface_c::DistSqrToSurfaceNodeWithinTolerance(const vec3 & CheckPt,
																double & NewDistSqrUnderTol,
																int & CloseNodeNum,
																const double & DistSqrTol)
{
	Boolean_t IsFound = FALSE;
	// first check the old closest node to see if it's under tolerance
	if (CloseNodeNum >= 0){
		NewDistSqrUnderTol = DistSqr(CheckPt, m_RefinedXYZList[CloseNodeNum]);
		if (DistSqrTol > 0)
			IsFound = (NewDistSqrUnderTol <= DistSqrTol);
		else
			IsFound = (NewDistSqrUnderTol <= m_MinNodeDistanceSqr);
	}
	for (int CloseNodeNum = 0; CloseNodeNum < m_RefinedXYZList.size() && !IsFound; ++CloseNodeNum){
		NewDistSqrUnderTol = DistSqr(CheckPt, m_RefinedXYZList[CloseNodeNum]);
		if (DistSqrTol > 0)
			IsFound = (NewDistSqrUnderTol <= DistSqrTol);
		else
			IsFound = (NewDistSqrUnderTol <= m_MinNodeDistanceSqr);
		//DoCheckNode[ZoneIndex] = (DistSqr(VolPt, NodePt) < m_MaxNodeDistanceSqr);
	}
	return IsFound;
}

// const Boolean_t FEVolume_c::DistSqrToSurfaceEdgeWithinTolerance(const vec3 & CheckPt,
// 	double & NewDistSqrUnderTol,
// 	const double & DistSqrTol)
// {
// 	Boolean_t IsFound = FALSE;
// 	vec3 NodePt;
// 
// 	for (int n = 0; n < m_NumNodes && !IsFound; ++n){
// 
// 	}
// 
// 	return IsFound;
// }

const Boolean_t FESurface_c::SubDivideIntegrateCellAtPoint(const vec3 & Point,
														const vector<vec3> & FarPoints,
														const vec3 & DelXYZ,
														const double & MinDistSqrToSurfaceNode,
														int MinDistNodeNum,
														const int & SubDivideLevel,
														const Boolean_t & IntegrateVolume)
{
	Boolean_t IsOk = TRUE;
	vec3 TmpPoint;
	double TmpDistSqr;

	for (double xi = -1.0; xi < 1.1; xi += 2.0){
		TmpPoint[0] = Point[0] + DelXYZ[0] * xi * 0.5;
		for (double yi = -1.0; yi < 1.1; yi += 2.0){
			TmpPoint[1] = Point[1] + DelXYZ[1] * yi * 0.5;
			for (double zi = -1.0; zi < 1.1; zi += 2.0){
				TmpPoint[2] = Point[2] + DelXYZ[2] * zi * 0.5;
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
				m_CheckSubPts.push_back(TmpPoint);
				if (m_InteriorSubPts.size() == 3966){
					int a = 1;
				}
#endif // _DEBUG

				if (SubDivideLevel > 0){
					if (DistSqrToSurfaceNodeWithinTolerance(TmpPoint, TmpDistSqr, MinDistNodeNum, MinDistSqrToSurfaceNode)){
						IsOk = SubDivideIntegrateCellAtPoint(TmpPoint, FarPoints, DelXYZ * 0.5, TmpDistSqr, MinDistNodeNum, SubDivideLevel - 1, IntegrateVolume);
					}
					else if (PointIsInterior(TmpPoint, FarPoints) && SetIndexAndWeightsForPoint(TmpPoint, m_VolZoneInfo)){
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
						m_InteriorSubPts.push_back(TmpPoint);

#endif // _DEBUG

						double CellVolume = 1;
						for (int dir = 0; dir < 3; ++dir)
							CellVolume *= DelXYZ[dir];
						if (IntegrateVolume)
							m_IntValues[m_NumIntVars] += CellVolume;
						for (int i = 0; i < m_NumIntVars; ++i){
							m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * CellVolume;
						}
					}
				}
				else if (PointIsInterior(TmpPoint, FarPoints) && SetIndexAndWeightsForPoint(TmpPoint, m_VolZoneInfo)){
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
					m_InteriorSubPts.push_back(TmpPoint);
#endif // _DEBUG
					double CellVolume = 1;
					for (int dir = 0; dir < 3; ++dir)
						CellVolume *= DelXYZ[dir];
					if (IntegrateVolume)
						m_IntValues[m_NumIntVars] += CellVolume;
					for (int i = 0; i < m_NumIntVars; ++i){
						m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * CellVolume;
					}
				}
			}
		}
	}

	return IsOk;
}


/*
	*	Begin implementation of Charles' integration code
	*/


vector<double> split(const string &s, char delim) {
	stringstream ss(s);
	string item;
	vector<double> tokens;
	while (getline(ss, item, delim)) {
		tokens.push_back(std::stod(item));
	}
	return tokens;
}

const mat LoadFile(const string & Path, const char & Delim){
	mat Points(1, 3);
	std::ifstream File(Path.c_str());

	if (File.is_open()){
		string Line;
		while (!File.eof()){
			std::getline(File, Line);
			vector<double> tmp = split(Line, Delim);
			if (tmp.size() > 0)
				Points = join_cols(Points, mat(tmp).t());
		}
		File.close();
	}
	else{
		TecUtilDialogErrMsg("Failed to open file");
	}

	return Points.tail_rows(Points.n_rows - 1);
}

const mat getPoints(const mat & Edge, const int & NumPts){
	mat Pts(NumPts, 3);

	mat EdgeDiff = Edge.tail_rows(Edge.n_rows - 1) - Edge.head_rows(Edge.n_rows - 1);

	vec Dist = join_cols(mat(vector<double>({ 0. })), sqrt(square(EdgeDiff.col(0)) + square(EdgeDiff.col(1)) + square(EdgeDiff.col(2))));

	for (int i = 1; i < Dist.n_elem; ++i)
		Dist(i) += Dist(i - 1);

	vec PtDist = linspace(0, Dist(Dist.n_elem - 1), NumPts);

	Pts.row(0) = Edge.row(0);
	Pts.row(NumPts - 1) = Edge.row(Edge.n_rows - 1);

	int Ind = 0;
	for (int i = 1; i < NumPts - 1; ++i){
		double d = PtDist(i);
		for (int j = Ind; j < Dist.n_elem - 1; ++j){
			if (Dist(j) <= d && Dist(j + 1) > d){
				Ind = j;
				break;
			}
		}
		double Rem = d - Dist(Ind);
		double Tot = Dist(Ind + 1) - Dist(Ind);
		double Rat = Rem / Tot;
		Pts.row(i) = Edge.row(Ind) + (Edge.row(Ind + 1) - Edge.row(Ind)) * Rat;
	}

	return Pts;
}

mat cubTrans(const mat & e1, const mat & e2, const mat & e3){

	mat Rpt;
	Rpt << 0 << 0 << 0 << endr
		<< 1. / 3. << 0 << 0 << endr
		<< 2. / 3. << 0 << 0 << endr

		<< 1 << 0 << 0 << endr
		<< 0 << 1. / 3. << 0 << endr
		<< 0 << 2. / 3. << 0 << endr

		<< 0 << 1 << 0 << endr
		<< 0 << 0 << 1. / 3. << endr
		<< 0 << 0 << 2. / 3. << endr

		<< 0 << 0 << 1 << endr
		<< 2. / 3. << 1. / 3. << 0 << endr
		<< 1. / 3. << 2. / 3. << 0 << endr

		<< 2. / 3. << 0 << 1. / 3. << endr
		<< 1. / 3. << 0 << 2. / 3. << endr
		<< 0 << 2. / 3. << 1. / 3. << endr

		<< 0 << 1. / 3. << 2. / 3. << endr
		<< 1. / 3. << 1. / 3. << 0 << endr
		<< 1. / 3. << 0 << 1. / 3. << endr

		<< 0 << 1. / 3. << 1. / 3. << endr
		<< 1. / 3. << 1. / 3. << 1. / 3.;

	// 		mat Rpt({
	// 			{ 0, 0, 0 }, { 1. / 3., 0, 0 }, { 2. / 3., 0, 0 },
	// 			{ 1, 0, 0 }, { 0, 1. / 3., 0 }, { 0, 2. / 3., 0 },
	// 			{ 0, 1, 0 }, { 0, 0, 1. / 3. }, { 0, 0, 2. / 3. },
	// 			{ 0, 0, 1 }, { 2. / 3., 1. / 3., 0 }, { 1. / 3., 2. / 3., 0 },
	// 			{ 2. / 3., 0, 1. / 3. }, { 1. / 3., 0, 2. / 3. }, { 0, 2. / 3., 1. / 3. },
	// 			{ 0, 1. / 3., 2. / 3. }, { 1. / 3., 1. / 3., 0 }, { 1. / 3., 0, 1. / 3. },
	// 			{ 0, 1. / 3., 1. / 3. }, { 1. / 3., 1. / 3., 1. / 3. }
	// 		});

	mat pt1 = getPoints(e1, 4),
		pt2 = getPoints(e2, 4),
		pt3 = getPoints(e3, 4);

	rowvec xx = getPoints(e1, 3).row(1),
		yy = getPoints(e2, 3).row(1),
		zz = getPoints(e3, 3).row(1);

	mat pt(20, 3);
	int pti = 0;
	for (int i = 0; i < 4; ++i)
		pt.row(pti++) = pt1.row(i);
	for (int i = 1; i < 4; ++i)
		pt.row(pti++) = pt2.row(i);
	for (int i = 1; i < 4; ++i)
		pt.row(pti++) = pt3.row(i);
	pt.row(pti++) = (pt1.row(3) * (2. / 3.)) + (pt2.row(3) * (1. / 3.));
	pt.row(pti++) = (pt1.row(3) * (1. / 3.)) + (pt2.row(3) * (2. / 3.));

	pt.row(pti++) = (pt1.row(3) * (2. / 3.)) + (pt3.row(3) * (1. / 3.));
	pt.row(pti++) = (pt1.row(3) * (1. / 3.)) + (pt3.row(3) * (2. / 3.));

	pt.row(pti++) = (pt2.row(3) * (2. / 3.)) + (pt3.row(3) * (1. / 3.));
	pt.row(pti++) = (pt2.row(3) * (1. / 3.)) + (pt3.row(3) * (2. / 3.));

	pt.row(pti++) = (xx * (1. / 2.)) + (yy * (1. / 2.));
	pt.row(pti++) = (xx * (1. / 2.)) + (zz * (1. / 2.));
	pt.row(pti++) = (yy * (1. / 2.)) + (zz * (1. / 2.));

	pt.row(pti++) = (pt1.row(3) * (1. / 3.)) + (pt2.row(3) * (1. / 3.)) + (pt3.row(3) * (1. / 3.));

	vec C1 = Rpt.col(0), C2 = Rpt.col(1), C3 = Rpt.col(2);
	vec C12 = square(C1), C22 = square(C2), C32 = square(C3);

	mat A(Rpt.n_rows, Rpt.n_rows);
	pti = 0;

	A.col(pti++) = vec(Rpt.n_rows, fill::ones);
	A.col(pti++) = C1;
	A.col(pti++) = C2;
	A.col(pti++) = C3;
	A.col(pti++) = C12;
	A.col(pti++) = C22;
	A.col(pti++) = C32;
	A.col(pti++) = C1%C2;
	A.col(pti++) = C1%C3;
	A.col(pti++) = C2%C3;
	A.col(pti++) = pow(C1, 3);
	A.col(pti++) = pow(C2, 3);
	A.col(pti++) = pow(C3, 3);
	A.col(pti++) = C12%C2;
	A.col(pti++) = C12%C3;
	A.col(pti++) = C22%C1;
	A.col(pti++) = C22%C3;
	A.col(pti++) = C32%C1;
	A.col(pti++) = C32%C2;
	A.col(pti++) = C1%C2%C3;

	mat abc(pt.n_rows, 3);

	for (int i = 0; i < 3; ++i){
		abc.col(i) = solve(A, pt.col(i));
	}

	return abc.t();
}

const mat cubTrans_func(const mat & abc,
	const vec & s,
	const vec & t,
	const vec & u)
{
	mat xyz(s.n_elem, 3);
	for (int i = 0; i < 3; ++i){
		xyz.col(i) =
			(s*abc(i, 1) + abc(i, 0))
			+ t*abc(i, 2)
			+ u*abc(i, 3)
			+ square(s)*abc(i, 4)
			+ square(t)*abc(i, 5)
			+ square(u)*abc(i, 6)
			+ s%t*abc(i, 7)
			+ s%u*abc(i, 8)
			+ t%u*abc(i, 9)
			+ pow(s, 3)*abc(i, 10)
			+ pow(t, 3)*abc(i, 11)
			+ pow(u, 3)*abc(i, 12)
			+ square(s) % t*abc(i, 13)
			+ square(s) % u*abc(i, 14)
			+ square(t) % s*abc(i, 15)
			+ square(t) % u*abc(i, 16)
			+ square(u) % s*abc(i, 17)
			+ square(u) % t*abc(i, 18)
			+ s%t%u*abc(i, 19);
	}
	return xyz;
}

const double cubJacobian(mat & abc, const double & s, const double & t, const double & u){
	vector<vec> dstu({ vector<double>({ 0., 1., 0., 0., s * 2, 0., 0., t, u, 0.,
		(s*s) * 3, 0., 0., s*t * 2, s*u * 2, t*t, 0., u*u, 0., t*u }),
		vector<double>({ 0., 0., 1., 0., 0., t * 2, 0., s, 0., u,
		0., (t*t) * 3, 0., s*s, 0., t*s * 2, t*u * 2, 0., u*u, s*u }),
		vector<double>({ 0., 0., 0., 1., 0., 0., u * 2, 0., s, t,
		0., 0., (u*u) * 2, 0., s*s, 0., t*t, u*s * 2, u*t * 2, s*t }) });

	mat33 J;
	for (unsigned int i = 0; i < 3; ++i){
		for (unsigned int j = 0; j < 3; ++j){
			J.at(i, j) = dot(abc.row(i), dstu[j]);
		}
	}
	return det(J);
}

void rquad(const int & N, const double & k, vec & x, vec & w)
{
	double k1 = k + 1, k2 = k + 2;

	vec n = linspace<vec>(1, N, N);

	vec nnk = 2 * n + k;

	rowvec A = join_rows(mat(vector<double>({ k / k2 })), ((ones<vec>(N) * (k*k)) / (nnk % (nnk + 2))).t());

	n = n.tail(N - 1);
	nnk = nnk.tail(N - 1);

	double B1 = 4. * k1 / (k2*k2*(k + 3));

	vec nk = n + k, nnk2 = square(nnk);
	vec B = 4. * square(n%nk) / (square(nnk2) - nnk2);

	mat ab = join_rows(A.t(), join_cols(vec(vector<double>({ pow(2, k1) / k1, B1 })), B));

	vec s = sqrt(ab(span(1, N - 1), 1));

	mat V;
	vec X;

	eig_sym(X, V, diagmat(ab(span(0, N - 1), 0)) + diagmat(s, -1) + diagmat(s, 1));

	x = (X + 1) / 2;

	w = pow(0.5, k1) * ab(0, 1) * square(V.row(0).t());

	return;
}

const vector<vec> tetraquad(const int & N, const mat & Verts)
{

	vector<vec> q(3), w123(3);

	for (int i = 0; i < 3; ++i)
		rquad(N, 2 - i, q[i], w123[i]);

	int N2 = N*N;
	int N3 = N2*N;
	vec q1(N3), q2(N3), q3(N3);

	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N; ++j){
			for (int k = 0; k < N; ++k){
				int Ind = i*N2 + j*N + k;
				q1[Ind] = q[0][j];
				q2[Ind] = q[1][k];
				q3[Ind] = q[2][i];
			}
		}
	}

	vec x = ones<vec>(N3) - q1,
		y = (ones<vec>(N3) -q2) % q1,
		z = q1 % q2 % q3;

	vec w = reshape(reshape(w123[1] * w123[0].t(), N2, 1) * w123[2].t(), N3, 1);
	// 		mat c = mat({ { 1, 0, 0, 0 }, { -1, 1, 0, 0 }, { -1, 0, 1, 0 }, { -1, 0, 0, 1 } }) * Verts;
	mat c = reshape(mat(vector<double>({ 1, 0, 0, 0, -1, 1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 1 })), 4, 4).t() * Verts;
	vec W = w * fabs(det(c.rows(1, 3)));

	mat XYZ = join_rows(ones<vec>(N3), join_rows(x, join_rows(y, z))) * c;

// 	vector<vector<double> > vXYZ(XYZ.n_rows, vector<double>(XYZ.n_cols));
// 	vector<double> vW(W.n_rows);
// 
// 	for (int i = 0; i < XYZ.n_rows; ++i){
// 		for (int j = 0; j < XYZ.n_cols; ++j){
// 			vXYZ[i][j] = XYZ(i, j);
// 		}
// 	}
// 	for (int i = 0; i < W.n_rows; ++i){
// 		vW[i] = W(i);
// 	}

	return vector<vec>({ XYZ.col(0), XYZ.col(1), XYZ.col(2), W });
}



const vector<vec> GetWeightsPoints(const int & N){

	mat Verts;
	Verts << 0 << 0 << 0 << endr
		<< 1 << 0 << 0 << endr
		<< 0 << 1 << 0 << endr
		<< 0 << 0 << 1;
	// 		mat Verts = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

	return tetraquad(N, Verts);
}

const double FESurface_c::IntVolume(const int & N, const vec3 & StartPoint) const{
// 	vector<int> NumElemsList = { 80, 320, 1280, 5120 };

	mat Verts;
	Verts << 0 << 0 << 0 << endr
		<< 1 << 0 << 0 << endr
		<< 0 << 1 << 0 << endr
		<< 0 << 0 << 1;
	// 		mat Verts = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

	vector<vec> stuW = tetraquad(N, Verts);

	int NumW = stuW[3].n_elem;

	vector<mat> Points(3, mat(m_GPList[0].GetCount() + 1, 3));
	for (int i = 0; i < 3; ++i){
		Points[i].row(0) = StartPoint.t();
		for (int j = 0; j < m_GPList[i].GetCount(); ++j){
			Points[i].row(j + 1) = m_GPList[i][j].t();
		}
	}

	mat abc = cubTrans(Points[0], Points[1], Points[2]);
	mat xyz = cubTrans_func(abc, stuW[0], stuW[1], stuW[2]);
	vec WTemp(NumW);

	for (int i = 0; i < NumW; ++i){
		WTemp(i) = stuW[3](i) * cubJacobian(abc, stuW[0](i), stuW[1](i), stuW[2](i));
	}

	return sum(WTemp);

	return 0;
}

const vector<mat> FESurface_c::GetIntegrationPointsWeights(const vec3 & StartPoint, const vector<vec> & stuW) const{

	if (stuW.size() == 4){

		int NumW = stuW[3].n_elem;

		vector<mat> Points(3, mat(m_GPList[0].GetCount() + 1, 3));
		for (int i = 0; i < 3; ++i){
			Points[i].row(0) = StartPoint.t();
			for (int j = 0; j < m_GPList[i].GetCount(); ++j){
				Points[i].row(j + 1) = m_GPList[i][j].t();
			}
		}

		mat abc = cubTrans(Points[0], Points[1], Points[2]);
		mat xyz = cubTrans_func(abc, stuW[0], stuW[1], stuW[2]);
		vec WTemp(NumW);

		for (int i = 0; i < NumW; ++i){
			WTemp(i) = stuW[3](i) * cubJacobian(abc, stuW[0](i), stuW[1](i), stuW[2](i));
		}
		return vector<mat>({ xyz, WTemp });
	}
	return vector<mat>();
}

const vector<mat> FESurface_c::GetIntegrationPointsWeights(const vector<vec> & stuW) const{

	if (stuW.size() == 4){

		int NumW = stuW[3].n_elem;

		vector<mat> Points(3, mat(m_GPList[0].GetCount(), 3));
		for (int i = 0; i < 3; ++i){
			for (int j = 0; j < m_GPList[i].GetCount(); ++j){
				Points[i].row(j) = m_GPList[i][j].t();
			}
		}

		mat abc = cubTrans(Points[0], Points[1], Points[2]);
		mat xyz = cubTrans_func(abc, stuW[0], stuW[1], stuW[2]);
		vec WTemp(NumW);

		for (int i = 0; i < NumW; ++i){
			WTemp(i) = stuW[3](i) * cubJacobian(abc, stuW[0](i), stuW[1](i), stuW[2](i));
		}
		return vector<mat>({ xyz, WTemp });
	}
	return vector<mat>();
}




/*
*	Private member functions
*/

const Boolean_t FESurface_c::PointIsInterior(const vec3 & Point, const vector<vec3> & FarPoints) const{
	/*
	*	Check that Point is in bounding cube of FEVolume
	*/
// 	Boolean_t IsInterior = (sum(Point >= m_ZoneMinXYZ) == 3 && sum(Point <= m_ZoneMaxXYZ) == 3);
	Boolean_t IsInterior;
	/*
	*	Now loop over each element, checking the triangles in that element
	*	for intersections between the triangle
	*	and a line segment between the Point and the FarPoints.
	*/
	int NumIntersections;
	bool FarPointGood;
	for (const vec3 & FarPoint : FarPoints){
		NumIntersections = 0;
		FarPointGood = true;
		vec3 T[3];
		for (int ElemNum = 0; ElemNum < m_NumElems && FarPointGood; ++ElemNum){
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					T[j][i] = m_XYZPtrs[i][m_ConnectivityListPtr[ElemNum * m_NumNodesPerElem + j]];
				}
			}
			switch (TriangleIntersect(T[0], T[1], T[2], FarPoint, Point)){
				case 1:
					NumIntersections++;
					break;
				case 2:
					FarPointGood = false;
					break;
			}
		}
		if (FarPointGood) break;
	}

	IsInterior = FarPointGood && (NumIntersections % 2 != 0);

	return IsInterior;
}
const int FESurface_c::TriangleIntersect(const vec3 & T_P0,
	const vec3 & T_P1,
	const vec3 & T_P2,
	const vec3 & R_P1,
	const vec3 & R_P0) const
{
	if (sum(T_P0 == T_P1) == 3 || sum(T_P0 == T_P2) == 3 || sum(T_P1 == T_P2) == 3)
		return 0;

	vec3 u = T_P1 - T_P0;
	vec3 v = T_P2 - T_P0;
	vec3 n = cross(u, v);
	vec3 w0 = R_P0 - T_P0;
	vec3 D = R_P1 - R_P0;

	//TODO: check for degeneracy 

	//Check for ray-plane parallelism
	double a = abs(dot(n, D));
	if (abs(dot(n, D)) < 1e-10){
		return 0;
	}

	//Check if ray goes towards triangle plane
	double r = -dot(n, w0) / dot(n, D);
	if (r < 0)
		return 0;

	//Calc intersect or R and the plane
	vec3 P = R_P0 + (D * r);

	vec3 w = P - T_P0;
	double uv = dot(u, v);
	double wv = dot(w, v);
	double vv = dot(v, v);
	double wu = dot(w, u);
	double uu = dot(u, u);

	double Denom = ((uv * uv) - (uu * vv));

	double s = ((uv*wv) - (vv*wu)) / Denom;
	double t = ((uv*wu) - (uu*wv)) / Denom;

	//TODO : boundary cases
	if (s >= 0 && t >= 0 && s + t <= 1)
		return 1;
	else
		return 0;
}


const vector<vector<double> > FESurface_c::TriSphereIntValsByElem(vector<double> * SphereTriangleAreas) const
{
	Boolean_t IsOk = m_IntegrationResultsReady && m_NumNodesPerElem == 3;

	vector<vector<double> > IntVals;

	if (IsOk){
		double TotalArea = 0.0;
		IntVals.resize(m_NumElems, vector<double>(m_IntValues.size()));
		/*
		*	Get total surface area of sphere, storing element areas
		*	as they're found.
		*/

		vector<double> AreaFactors(m_NumElems);
		vec3 T[3];
		double A, B, C;

		for (int ElemNum = 0; ElemNum < m_NumElems; ++ElemNum){
			/*
			*	Get the corners of the triangle element.
			*/
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					T[j][i] = m_XYZPtrs[i][m_ConnectivityListPtr[ElemNum * m_NumNodesPerElem + j]];
				}
			}

			/*
			*	Get the area of the element
			*/
			A = T[0][2]*(T[1][1] - T[2][1])
				+ T[1][2]*(T[2][1] - T[0][1])
				+ T[2][2]*(T[0][1] - T[1][1]);

			B = T[0][0]*(T[1][2] - T[2][2])
				+ T[1][0]*(T[2][2] - T[0][2])
				+ T[2][0]*(T[0][2] - T[1][2]);

			C = T[0][1]*(T[1][0] - T[2][0])
				+ T[1][1]*(T[2][0] - T[0][0])
				+ T[2][1]*(T[0][0] - T[1][0]);

			AreaFactors[ElemNum] = 0.5 * sqrt(A*A + B*B + C*C);
			TotalArea += AreaFactors[ElemNum];
		}

		if (SphereTriangleAreas != NULL)
			*SphereTriangleAreas = AreaFactors;

		for (int ElemNum = 0; ElemNum < m_NumElems; ++ElemNum){
			AreaFactors[ElemNum] /= TotalArea;
			for (int VarNum = 0; VarNum < m_IntValues.size(); ++VarNum)
				IntVals[ElemNum][VarNum] = m_IntValues[VarNum] * AreaFactors[ElemNum];
		}
	}

	return IntVals;
}



/*
 *	Begin Domain_c functions
 */


void Domain_c::Setup(const vector<int> & V, FESurface_c *Vol)
{
	m_V = V;
	m_Vol = Vol;
	m_E.resize(m_V.size(), vector<int>(2));
	for (int i = 0; i < m_V.size(); ++i){
		for (int j = 0; j < 2; ++j){
			m_E[i][j] = (i + j) % m_V.size();
		}
	}
}

const double Domain_c::Weight() const
{
	if (m_V.size() <= 2)
		return 0.0;
	vector<double> MinTriAreas(m_E.size() - 2);
	vector<vector<int> > MinTris(m_E.size() - 2);
	vector<Domain_c> D1(m_E.size() - 2),
		D2(m_E.size() - 2);
	int tMin;
	double MinArea = DBL_MAX;
	for (int i = 0; i < MinTriAreas.size(); ++i){
		MinTriAreas[i] = Split(i, MinTris[i], D1[i], D2[i]);
		if (MinTriAreas[i] < MinArea){
			MinArea = MinTriAreas[i];
			tMin = i;
		}
	}
	m_Vol->m_ElemList.push_back(MinTris[tMin]);
	return MinTriAreas[tMin] + D1[tMin].Weight() + D2[tMin].Weight();
}

const double Domain_c::Split(const int & ei, vector<int> & t, Domain_c & D1, Domain_c & D2) const
{
	vector<vector<vector<int> > > TriList(m_V.size() - 2 - ei, vector<vector<int>  >(2, vector<int>(3)));
	for (int i = 0; i < TriList.size(); ++i){
		for (int j = 0; j < 2; ++j){
			TriList[i][0][j] = m_V[m_E[ei][j]];
			TriList[i][1][j] = m_E[ei][j];
		}
		TriList[i][0][2] = m_V[(m_E[ei][1] + 1 + i) % m_V.size()];
		TriList[i][1][2] = (m_E[ei][1] + 1 + i) % m_V.size();
	}
	int tMin;
	double TempArea, MinArea = 1e150;
	for (int i = 0; i < TriList.size(); ++i){
		TempArea = WeightFunc(TriList[i][0]);
		if (TempArea < MinArea){
			MinArea = TempArea;
			tMin = i;
		}
	}
	t = TriList[tMin][0];
	vector<int> ti = TriList[tMin][1];
	if (abs(ti[0] - ti[2]) < abs(ti[1] - ti[2]) && ti[0] < ti[1]){
		vector<int> V(ti[0] - ti[2] + 1);
		for (int i = 0; i < V.size(); ++i){
			V[i] = m_V[ti[2] + i];
		}
		D1.Setup(V, m_Vol);
		V = vector<int>(m_V.size() - ti[1] + ti[2] + 1);
		for (int i = 0; i < V.size(); ++i){
			V[i] = m_V[(ti[1] + i) % m_V.size()];
		}
		D2.Setup(V, m_Vol);
	}
	else{
		vector<int> V(m_V.size() + ti[0] - ti[2] + 1);
		for (int i = 0; i < V.size(); ++i){
			V[i] = m_V[(ti[2] + i) % m_V.size()];
		}
		D1.Setup(V, m_Vol);
		V = vector<int>(ti[2] - ti[1] + 1);
		for (int i = 0; i < V.size(); ++i){
			V[i] = m_V[ti[1] + i];
		}
		D2.Setup(V, m_Vol);
	}

	return MinArea;
}

const double Domain_c::WeightFunc(const vector<int> & t) const
{
	vec3 AB = m_Vol->m_XYZList[t[1]] - m_Vol->m_XYZList[t[0]],
		AC = m_Vol->m_XYZList[t[2]] - m_Vol->m_XYZList[t[0]];

	double Area = (pow(AB[1] * AC[2] - AB[2] * AC[1], 2)
		+ pow(AB[2] * AC[0] - AB[0] * AC[2], 2)
		+ pow(AB[0] * AC[1] - AB[1] * AC[0], 2));

	return Area;
}
