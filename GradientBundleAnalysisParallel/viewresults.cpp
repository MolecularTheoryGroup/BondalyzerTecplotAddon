
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

#include <cstring>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <ctime>

#include <armadillo>

#include "Set.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_DATA_TYPES.h"
#include "VIEWRESULTS.h"

using namespace arma;
using namespace tecplot::toolbox;

using std::string;
using std::to_string;
using std::stringstream;
using std::vector;
using std::stoi;
using std::ofstream;


string const GBADelim = " | ";
void GBAResultViewerPopulateGBs() {
	/*
	*	Populate list of gradient bundles
	*/
	TecGUIListDeleteAllItems(MLSelGB_MLST_T3_1);

	int NumZones = TecUtilDataSetGetNumZones();

	std::set<string> GBNameMap;
	string TmpStr;
	for (int ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum) {
		if ((AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasin)
			|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasin))
			&& AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinName, TmpStr)) {
			string intVar;
			if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.CondensedBasinDefiningVariable, intVar)){
				TmpStr += GBADelim + intVar;
			}
			
			GBNameMap.insert(TmpStr);
		}
	}
	for (auto const & i : GBNameMap)
		TecGUIListAppendItem(MLSelGB_MLST_T3_1, i.c_str());
}


void GBAResultViewerSelectSphere(){
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);


	Boolean_t IsOk = TRUE;

	/*
	*	Clear gradient bundle list
	*/
// 	TecGUIListDeleteAllItems(MLSelGB_MLST_T3_1);

	/*
	*	Get selected sphere name
	*/
	LgIndex_t SelSphereNum = TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1);
	char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_T3_1, SelSphereNum);
	string SphereNameStr = SphereNameCStr;
	TecUtilStringDealloc(&SphereNameCStr);
	/*
	*	Find sphere zone and
	*	1) check if it's active
	*	2) get it's gradient bundle volume CPs
	*/

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();

	EntIndex_t SphereZoneNum;

	if (IsOk){
// 		TecUtilDataLoadBegin();

		Boolean_t IsFound = FALSE;
		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && !IsFound; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)
				&& (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
					|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr)))
			{
				IsFound = TRUE;
				SphereZoneNum = ZoneNum;
			}
		}

// 		TecUtilDataLoadEnd();

		IsOk = IsFound;
	}



	if (IsOk)
		TecGUIToggleSet(TGLSphereVis_TOG_T3_1, TecUtilZoneIsActive(SphereZoneNum));


	int SelectedVarNum = 1;
	if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0)
		SelectedVarNum = TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1);

	TecGUIListDeleteAllItems(SLSelVar_SLST_T3_1);

	if (TecGUIToggleGet(TGLSphereVis_TOG_T3_1)) {
		/*
		 *	Get list of integrated variables for the selected sphere
		 */

		vector<string> IntVarNames;
		vector<string> IntCheckStrs = { "I: ", "IN: ", "INS: ", " Integration" };
		for (int i = 1; i <= TecUtilDataSetGetNumVars(); ++i) {
			char *VarName, *CheckStr;
			if (TecUtilVarGetName(i, &VarName)) {
				for (string const & Str : IntCheckStrs) {
					CheckStr = std::strstr(VarName, Str.c_str());
					if (CheckStr != nullptr) {
						// Integration variable found. Now make sure it's not bit type for sphere zone.
						if (TecUtilDataValueGetType(SphereZoneNum, i) != FieldDataType_Bit) {
							IntVarNames.push_back(VarName);
							break;
						}
					}
				}
				TecUtilStringDealloc(&VarName);
			}
		}

		if (IntVarNames.size() > 0) {
			for (string const & i : IntVarNames)
				TecGUIListAppendItem(SLSelVar_SLST_T3_1, i.c_str());
			TecGUIListSetSelectedItem(SLSelVar_SLST_T3_1, SelectedVarNum);
		}

	}

	TecUtilDrawGraphics(TRUE);
	TecUtilLockFinish(AddOnID);
	return;
}

void GBAResultViewerSelectIntVar(){
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	Boolean_t IsOk = TRUE;

	LgIndex_t *SelectedNums;
	LgIndex_t NumSelected;

	TecGUIListGetSelectedItems(SLSelVar_SLST_T3_1, &SelectedNums, &NumSelected);

	EntIndex_t VarNum = -1;
	EntIndex_t SphereZoneNum = -1;
	int NumContours = stoi(TecGUITextFieldGetString(TFNumContours_TF_T3_1));
	bool LogSpaceContours = (TecGUIRadioBoxGetToggle(RBLogLin_RADIO_T3_1) == 1);
	bool UseSelectedSphere = (TecGUIRadioBoxGetToggle(RBCntSrc_RADIO_T3_1) == 1);
	vec ContourLevelsVec;
	string SphereNameStr;
	if (NumSelected > 0){
		char *IntVarNameCStr = TecGUIListGetString(SLSelVar_SLST_T3_1, TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1));
		string IntVarNameStr = IntVarNameCStr;
		TecUtilStringDealloc(&IntVarNameCStr);

		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_T3_1, TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1));
		SphereNameStr = SphereNameCStr;
		TecUtilStringDealloc(&SphereNameCStr);

		// Get the variable number
		VarNum = TecUtilVarGetNumByName(IntVarNameStr.c_str());
// 		for (int i = 1; i <= TecUtilDataSetGetNumVars() && VarNum <= 0; ++i){
// 			char *ChkName;
// 			if (TecUtilVarGetName(i, &ChkName)){
// 				if (IntVarNameStr.compare(ChkName) == 0){
// 					VarNum = i;
// 				}
// 				TecUtilStringDealloc(&ChkName);
// 			}
// 		}

		// Find the selected sphere zone

		// get var min and max for sphere zone

		double VarMin = DBL_MAX, VarMax = -DBL_MAX, TmpMin, TmpMax;

		vector<int> SphereZoneNums;

		for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)
				&& (!UseSelectedSphere 
					|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
					|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr))){
				SphereZoneNums.push_back(ZoneNum);
				IsOk = TecUtilDataValueGetMinMaxByZoneVar(ZoneNum, VarNum, &TmpMin, &TmpMax);
				if (TmpMin != TmpMax){
					VarMin = MIN(VarMin, TmpMin);
					VarMax = MAX(VarMax, TmpMax);
				}
			}
		}

		if (VarMin == DBL_MAX || VarMax == -DBL_MAX){
			TecUtilDialogMessageBox("No values for selected variable for selected sphere(s)", MessageBoxType_Error);
			TecUtilDrawGraphics(TRUE);
			TecUtilLockFinish(AddOnID);
			return;
		}


		if (LogSpaceContours){
			if (VarMin <= 0){
				/*
				*	The variable spans (or stops/starts at) 0, so need to find the negative and positive numbers closest to 0
				*	and use their exponents as a range in the logspace call.
				*/
				double NegCloseToZero = -DBL_MAX,
					PosCloseToZero = DBL_MAX;
				for (int ZoneNum : SphereZoneNums){
					FieldData_pa TmpRef = TecUtilDataValueGetReadableNativeRef(ZoneNum, VarNum);
					if (VALID_REF(TmpRef)){
						vector<int> MaxIJK(3);
						TecUtilZoneGetIJK(ZoneNum, &MaxIJK[0], &MaxIJK[1], &MaxIJK[2]);
						vector<double> TmpVals;
						if (TecUtilDataValueGetLocationByRef(TmpRef) == ValueLocation_Nodal)
							TmpVals.resize(MaxIJK[0]);
						else
							TmpVals.resize(MaxIJK[1]);
						TecUtilDataValueArrayGetByRef(TmpRef, 1, TmpVals.size(), (void*)TmpVals.data());
						for (double const & Val : TmpVals){
							if (Val < 0.0 && Val > NegCloseToZero)
								NegCloseToZero = Val;
							else if (Val > 0 && Val < PosCloseToZero)
								PosCloseToZero = Val;
						}
					}
				}
				if (NegCloseToZero == -DBL_MAX)
					NegCloseToZero = -1e-3;
				if (PosCloseToZero == DBL_MAX)
					PosCloseToZero = 1e-3;

				if (VarMax <= 0){
					ContourLevelsVec = LogSpace(-VarMax, -VarMin, NumContours) * -1.0;
				}
				else if (VarMin >= 0.0){
					ContourLevelsVec = LogSpace(VarMin, VarMax, NumContours);
				}
				else{
					int NumNegContours = int(double(NumContours) * (-VarMin) / (VarMax - VarMin));

					ContourLevelsVec = join_cols(flipud(LogSpace(-NegCloseToZero, -VarMin, NumNegContours) * -1.0), LogSpace(PosCloseToZero, VarMax, NumContours - NumNegContours));
				}
			}
			else{
				ContourLevelsVec = LogSpace(VarMin, VarMax, NumContours);
			}
		}
		else{
			ContourLevelsVec = linspace(VarMin, VarMax, NumContours);
		}

		
	}

	if (IsOk && VarNum > 0){
		// Set the 8th contour setting to the new values
		ArgList_pa TempArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 8);
		TecUtilArgListAppendInt(TempArgList, SV_VAR, VarNum);
		TecUtilContourSetVariableX(TempArgList);
		TecUtilArgListClear(TempArgList);

		TecUtilArgListAppendInt(TempArgList, SV_CONTOURLABELACTION, ContourLabelAction_DeleteAll);
		TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 8);
		TecUtilContourLabelX(TempArgList);
		TecUtilArgListClear(TempArgList);

		TecUtilArgListAppendInt(TempArgList, SV_CONTOURLEVELACTION, ContourLevelAction_New);
		TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 8);
		TecUtilArgListAppendInt(TempArgList, SV_NUMVALUES, NumContours);
		TecUtilArgListAppendArray(TempArgList, SV_RAWDATA, ContourLevelsVec.memptr());
		TecUtilContourLevelX(TempArgList);
		TecUtilArgListClear(TempArgList);

		TecUtilArgListAppendString(TempArgList, SV_P1, SV_GLOBALCONTOUR);
		TecUtilArgListAppendString(TempArgList, SV_P2, SV_LEGEND);
		TecUtilArgListAppendString(TempArgList, SV_P3, SV_SHOW);
		TecUtilArgListAppendInt(TempArgList, SV_OFFSET1, 8);
		TecUtilArgListAppendArbParam(TempArgList, SV_IVALUE, FALSE);
		TecUtilStyleSetLowLevelX(TempArgList);
		TecUtilArgListClear(TempArgList);

		// Set all zones associated with the selected sphere to the right contour settings

		Set_pa ZoneSet = TecUtilSetAlloc(FALSE);

		for (int ZoneNum = 1; ZoneNum <= TecUtilDataSetGetNumZones(); ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr)){
				TecUtilSetAddMember(ZoneSet, ZoneNum, FALSE);
			}
		}

		TecUtilZoneSetContour(SV_SHOW, ZoneSet, 0.0, TRUE);
		TecUtilZoneSetContour(SV_LINECONTOURGROUP, ZoneSet, 0.0, ContourColoring_Group8);
		TecUtilZoneSetContour(SV_CONTOURTYPE, ZoneSet, 0.0, (ContourType_e)2);
		TecUtilZoneSetContour(SV_FLOODCOLORING, ZoneSet, 0.0, ContourColoring_Group8);
		TecUtilZoneSetMesh(SV_SHOW, ZoneSet, 0.0, FALSE);

		TecUtilArgListClear(TempArgList);
		TecUtilArgListAppendString(TempArgList, SV_P1, SV_GLOBALCONTOUR);
		TecUtilArgListAppendString(TempArgList, SV_P2, SV_LEGEND);
		TecUtilArgListAppendString(TempArgList, SV_P3, SV_SHOW);
		TecUtilArgListAppendInt(TempArgList, SV_OFFSET1, 8);
		TecUtilArgListAppendArbParam(TempArgList, SV_IVALUE, (Boolean_t)1);
		TecUtilStyleSetLowLevelX(TempArgList);

		TecUtilArgListClear(TempArgList);
		TecUtilArgListAppendString(TempArgList, SV_P1, SV_GLOBALCONTOUR);
		TecUtilArgListAppendString(TempArgList, SV_P2, SV_LABELS);
		TecUtilArgListAppendString(TempArgList, SV_P3, SV_AUTOLEVELSKIP);
		TecUtilArgListAppendInt(TempArgList, SV_OFFSET1, 8);
		TecUtilArgListAppendArbParam(TempArgList, SV_IVALUE, (SmInteger_t)(ContourLevelsVec.size() / 15));
		TecUtilStyleSetLowLevelX(TempArgList);

		TecUtilArgListDealloc(&TempArgList);
		TecUtilSetDealloc(&ZoneSet);
	}

	TecUtilDrawGraphics(TRUE);
	TecUtilLockFinish(AddOnID);
}

void GBAResultViewerSelectCondensedGBs() {
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	Boolean_t IsOk = TRUE;

	LgIndex_t *SelectedNums;
	LgIndex_t NumSelected;

	TecGUIListGetSelectedItems(MLSelGB_MLST_T3_1, &SelectedNums, &NumSelected);

	Set ActivateSet, DeactivateSet;
	std::set<int> ActivateSetSet, DeactivateSetSet;
	

	if (NumSelected > 0) {

		// 		TecUtilDataLoadBegin();

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		string ChkNameStr;

		for (int i = 0; i < NumSelected && IsOk; ++i) {

			char* ChkName = TecGUIListGetString(MLSelGB_MLST_T3_1, SelectedNums[i]);
			ChkNameStr = ChkName;
			TecUtilStringDealloc(&ChkName);

			string intVar;
			vector<string> strVec = SplitString(ChkNameStr, GBADelim);
			if (strVec.size() > 1){
				ChkNameStr = strVec[0];
				intVar = strVec[1];
			}

			bool useInfo = false;
			vector<string> tmpStrVec = SplitString(ChkNameStr, " basin");
			if (tmpStrVec.size() > 1) {
				ChkNameStr = tmpStrVec[0] + tmpStrVec[1];
				useInfo = true;
			}

			for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum) {
				if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasinWedge)
					|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasinWedge)
					|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasin)
					|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasin))
				{
					if ((useInfo && AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.CondensedBasinInfo, ChkNameStr))
						|| (!useInfo && AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.CondensedBasinName, ChkNameStr))
						&& (strVec.size() < 2 || AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.CondensedBasinDefiningVariable, intVar))) {
						ActivateSetSet.insert(ZoneNum);
					}
					else
						DeactivateSetSet.insert(ZoneNum);
				}
			}
		}

		if (IsOk) {
			for (auto const & i : ActivateSetSet)
				ActivateSet += i;
			for (auto const & i : DeactivateSetSet)
				if (ActivateSetSet.count(i) == 0)
					DeactivateSet += i;
			TecUtilZoneSetActive(ActivateSet.getRef(), AssignOp_PlusEquals);
			TecUtilZoneSetActive(DeactivateSet.getRef(), AssignOp_MinusEquals);
		}

		// 		TecUtilDataLoadEnd();
	}
	else if (TecGUIListGetItemCount(MLSelGB_MLST_T3_1) == 0) {
		// 		TecUtilDataLoadBegin();

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum) {
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasinWedge)
				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasinWedge))
			{
				DeactivateSet += ZoneNum;
			}
		}

		if (IsOk) {
			TecUtilZoneSetActive(DeactivateSet.getRef(), AssignOp_MinusEquals);
		}

		// 		TecUtilDataLoadEnd();
	}

	TecUtilArrayDealloc((void**)&SelectedNums);

	TecUtilDrawGraphics(TRUE);
	TecUtilLockFinish(AddOnID);
}

void GBAResultViewerSelectGB(){
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	Boolean_t IsOk = TRUE;

	LgIndex_t *SelectedNums;
	LgIndex_t NumSelected;

	TecGUIListGetSelectedItems(MLSelGB_MLST_T3_1, &SelectedNums, &NumSelected);

	if (NumSelected > 0){

// 		TecUtilDataLoadBegin();

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_T3_1, TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1));
		string SphereNameStr = SphereNameCStr;
		TecUtilStringDealloc(&SphereNameCStr);

		string ChkNameStr;

		Set_pa ActivateSet = TecUtilSetAlloc(FALSE);
		Set_pa DeactivateSet = TecUtilSetAlloc(FALSE);

		for (int i = 0; i < NumSelected && IsOk; ++i){

			char* ChkName = TecGUIListGetString(MLSelGB_MLST_T3_1, SelectedNums[i]);
			ChkNameStr = ChkName;
			TecUtilStringDealloc(&ChkName);

			for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
				if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
					|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr))
				{
					if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.VolumeCPName, ChkNameStr)){
						IsOk = TecUtilSetAddMember(ActivateSet, ZoneNum, FALSE);
					}
					else{
						Boolean_t IsFound = FALSE;
						string TempStr = AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.VolumeCPName);
						for (int j = 0; j < NumSelected && !IsFound; ++j){
							char *TempCStr2 = TecGUIListGetString(MLSelGB_MLST_T3_1, SelectedNums[j]);
							IsFound = (TempStr.compare(TempCStr2) == 0);
							TecUtilStringDealloc(&TempCStr2);
						}
						if (!IsFound)
							IsOk = TecUtilSetAddMember(DeactivateSet, ZoneNum, FALSE);
					}
				}
			}
		}

		if (IsOk){
			TecUtilZoneSetActive(ActivateSet, AssignOp_PlusEquals);
			TecUtilZoneSetActive(DeactivateSet, AssignOp_MinusEquals);
		}


		TecUtilSetDealloc(&ActivateSet);
		TecUtilSetDealloc(&DeactivateSet);

// 		TecUtilDataLoadEnd();
	}
	else if (TecGUIListGetItemCount(MLSelGB_MLST_T3_1) == 0){
// 		TecUtilDataLoadBegin();

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_T3_1, TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1));
		string SphereNameStr = SphereNameCStr;
		TecUtilStringDealloc(&SphereNameCStr);

		Set_pa DeactivateSet = TecUtilSetAlloc(FALSE);

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr))
			{
				IsOk = TecUtilSetAddMember(DeactivateSet, ZoneNum, FALSE);
			}
		}

		if (IsOk){
			TecUtilZoneSetActive(DeactivateSet, AssignOp_MinusEquals);
		}

		TecUtilSetDealloc(&DeactivateSet);

// 		TecUtilDataLoadEnd();
	}

	TecUtilArrayDealloc((void**)&SelectedNums);

	TecUtilDrawGraphics(TRUE);
	TecUtilLockFinish(AddOnID);
}

void GBAResultViewerToggleSphere(){
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	Boolean_t IsOk = (TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0);

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();

	Boolean_t ActivateZone = TecGUIToggleGet(TGLSphereVis_TOG_T3_1);

	if (IsOk)
	{
		/*
		 *	Find and activate or deactivate sphere zone
		 */
		char* ChkName = TecGUIListGetString(SLSelSphere_SLST_T3_1, TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1));
		string ChkNameStr = ChkName;
		TecUtilStringDealloc(&ChkName);

		TecUtilDataLoadBegin();

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)
				&& AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, ChkNameStr))
			{
				Set_pa TempSet = TecUtilSetAlloc(FALSE);
				TecUtilSetAddMember(TempSet, ZoneNum, FALSE);
				if (ActivateZone){
					TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
				}
				else{
					TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
					TecGUIListDeleteAllItems(MLSelGB_MLST_T3_1);
				}
				TecUtilSetDealloc(&TempSet);
				break;
			}
		}

		/*
		*	Find and activate or deactivate IB edge zone
		*/
		ChkName = TecGUIListGetString(SLSelSphere_SLST_T3_1, TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1));
		ChkNameStr = ChkName;
		TecUtilStringDealloc(&ChkName);

		TecUtilDataLoadBegin();

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)
				&& AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, ChkNameStr))
			{
				Set_pa TempSet = TecUtilSetAlloc(FALSE);
				TecUtilSetAddMember(TempSet, ZoneNum, FALSE);
				if (ActivateZone){
					TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
				}
				else{
					TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
					TecGUIListDeleteAllItems(MLSelGB_MLST_T3_1);
				}
				TecUtilSetDealloc(&TempSet);
				break;
			}
		}
	}
	else{
		TecGUIToggleSet(TGLSphereVis_TOG_T3_1, FALSE);
	}

	TecUtilDataLoadEnd();
	GBAResultViewerSelectSphere();

	TecUtilDrawGraphics(TRUE);
	TecUtilLockFinish(AddOnID);
}

void SortCPNameList(vector<string> & StrList){
	/*
	 *	These lists are never that big, so let's
	 *	bubble sort it up!
	 */


	/*	First get number of each CP type so that we can
	 *	simply sort by total CP index.
	 */
	Boolean_t IsOk = TRUE;
	EntIndex_t CPZoneNum = ZoneNumByName(string("Critical Points"));
	IsOk = (CPZoneNum > 0);
	int NumCPs[4];
	string TmpString;
	for (int i = 0; i < 4 && IsOk; ++i){
		IsOk = AuxDataZoneGetItem(CPZoneNum, CCDataNumCPs[i], TmpString);
		if (IsOk){
			NumCPs[i] = stoi(TmpString);
		}
	}

	/*
	 *	Now sort
	 */
	if (IsOk){
		string TmpStr;
		int TmpInt;
		for (int i = 0; i < StrList.size(); ++i){
			Boolean_t DidSwap = FALSE;
			for (int j = 0; j < StrList.size() - 1; ++j){
				stringstream ss;
				ss << StrList[j];
				ss >> TmpStr >> TmpInt;
				int CPOffset = 0;
				int CPj, CPjp1;
				for (int k = 0; k < 4; ++k){
					if (TmpStr == RankStrs[k]){
						CPj = TmpInt + CPOffset;
						break;
					}
					CPOffset += NumCPs[k];
				}
				ss.str(string());
				ss.clear();
				ss << StrList[j + 1];
				ss >> TmpStr >> TmpInt;
				CPOffset = 0;
				for (int k = 0; k < 4; ++k){
					if (TmpStr == RankStrs[k]){
						CPjp1 = TmpInt + CPOffset;
						break;
					}
					CPOffset += NumCPs[k];
				}
				if (CPjp1 < CPj){
					TmpStr = StrList[j];
					StrList[j] = StrList[j + 1];
					StrList[j + 1] = TmpStr;
					DidSwap = TRUE;
				}
			}
		}
	}
}

void GBAResultViewerDeleteSphere(){
	/*
	 *	Delete selected sphere zone and and all
	 *	GBA zones associated with it.
	 */
	TecUtilLockStart(AddOnID);

	if (TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0){
		TecUtilDrawGraphics(FALSE);
		/*
		*	Get selected sphere name
		*/
		LgIndex_t SelSphereNum = TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1);
		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_T3_1, SelSphereNum);
		string SphereNameStr = SphereNameCStr;
		string SphereNameForCondensedBasins = SphereNameStr + " Sphere";
		TecUtilStringDealloc(&SphereNameCStr);
		/*
		*	Find sphere zone and
		*	1) check if it's active
		*	2) get it's gradient bundle volume CPs
		*/

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		Set_pa DeleteSet = TecUtilSetAlloc(FALSE);

		TecUtilDataLoadBegin();
		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameForCondensedBasins)
				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr)){
				TecUtilSetAddMember(DeleteSet, ZoneNum, FALSE);
			}
		}

		if (!TecUtilSetIsEmpty(DeleteSet) && !TecUtilDataSetDeleteZone(DeleteSet)){
			TecUtilDialogErrMsg("Failed to delete zones.");
		}

		TecUtilSetDealloc(&DeleteSet);


		TecGUIListDeleteItemAtPos(SLSelSphere_SLST_T3_1, SelSphereNum);
		if (TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0){
			TecGUIListSetSelectedItem(SLSelSphere_SLST_T3_1, 1);
			GBAResultViewerSelectSphere();
		}
		else{
			TecGUIListDeleteAllItems(MLSelGB_MLST_T3_1);
			TecGUIToggleSet(TGLSphereVis_TOG_T3_1, FALSE);
		}

		if (TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) == 0){
			vector<string> VarNames(TecUtilDataSetGetNumVars());
			for (int i = 0; i < VarNames.size(); ++i){
				char *cStr;
				TecUtilVarGetName(i + 1,&cStr);
				VarNames[i] = cStr;
				TecUtilStringDealloc(&cStr);
			}
			Set DelVars;
			for (int vNum = 0; vNum < VarNames.size(); ++vNum){
				for (string s2 : vector<string>({ "I: ", "IN: ", "INS: " })){
					if (s2.length() <= VarNames[vNum].length() && VarNames[vNum].substr(0, s2.length()) == s2){
						DelVars += vNum + 1;
					}
				}
			}
			if (!DelVars.isEmpty()){
				TecUtilDataSetDeleteVar(DelVars.getRef());
			}
		}


		TecUtilDataLoadEnd();
		TecUtilDrawGraphics(TRUE);
	}
	TecUtilLockFinish(AddOnID);

}

void GBAResultViewerActivateAllGB(){
	TecUtilLockStart(AddOnID);

	if (TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0){
		TecUtilDrawGraphics(FALSE);
		/*
		*	Get selected sphere name
		*/
		LgIndex_t SelSphereNum = TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1);
		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_T3_1, SelSphereNum);
		string SphereNameStr = SphereNameCStr;

		/*
		*	Find sphere zone and
		*	1) check if it's active
		*	2) get it's gradient bundle volume CPs
		*/

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		Set_pa ActivateSet = TecUtilSetAlloc(FALSE);

		TecUtilDataLoadBegin();

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr))
			{
				TecUtilSetAddMember(ActivateSet, ZoneNum, FALSE);
			}
		}

		SetIndex_t NumGBZones = TecUtilSetGetMemberCount(ActivateSet);
		int NumActiveGB = 0;
		SetIndex_t ZoneNum;
		TecUtilSetForEachMember(ZoneNum, ActivateSet){
			NumActiveGB += TecUtilZoneIsActive((EntIndex_t)ZoneNum);
		}

		if ((double)NumActiveGB / (double)NumGBZones >= 0.5){
			TecUtilZoneSetActive(ActivateSet, AssignOp_MinusEquals);
			GBAResultViewerSelectSphere();
		}
		else{
			TecUtilZoneSetActive(ActivateSet, AssignOp_PlusEquals);
			int NumGBAtoms = TecGUIListGetItemCount(MLSelGB_MLST_T3_1);
			vector<int> SelNums(NumGBAtoms);
			for (int i = 1; i <= NumGBAtoms; ++i)
				SelNums[i - 1] = i;

			TecGUIListSetSelectedItems(MLSelGB_MLST_T3_1, SelNums.data(), NumGBAtoms);
		}

		TecUtilSetDealloc(&ActivateSet);

		TecUtilStringDealloc(&SphereNameCStr);

		TecUtilDataLoadEnd();
		TecUtilDrawGraphics(TRUE);
	}

	TecUtilLockFinish(AddOnID);
}

/*
 *	Volume zone toggle probe functions
 */

void ToggleFEVolumesProbeInstallCB(){
	ArgList_pa ProbeArgs = TecUtilArgListAlloc();
	TecUtilArgListAppendFunction(ProbeArgs, SV_CALLBACKFUNCTION, ToggleFEVolumesProbeCB);
	TecUtilArgListAppendString(ProbeArgs, SV_STATUSLINETEXT, "Select CPs to define a plane");
	TecUtilArgListAppendArbParam(ProbeArgs, SV_CLIENTDATA, ArbParam_t(NULL));
	if (!TecUtilProbeInstallCallbackX(ProbeArgs)){
		TecUtilDialogErrMsg("Failed to install probe callback.");
	}
	TecUtilArgListDealloc(&ProbeArgs);
}

void STDCALL ToggleFEVolumesProbeCB(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData)
{


	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	if (WasSuccessful){

		string TmpStr1, TmpStr2, TmpStr3;

		Boolean_t IsOk = TRUE;

		EntIndex_t ProbedZoneNum = TecUtilProbeFieldGetZone();

		string TmpStr;

		string CPName;

		TecUtilDataLoadBegin();

		if (IsOk){
			IsOk = AuxDataZoneGetItem(ProbedZoneNum, CSMAuxData.GBA.SphereCPName, CPName);
		}

		string ZoneType;
		if (IsOk){
			IsOk = AuxDataZoneGetItem(ProbedZoneNum, CSMAuxData.GBA.ZoneType, ZoneType);
		}

		if (IsOk){
			if (ZoneType == CSMAuxData.GBA.ZoneTypeSphereZone){

				EntIndex_t NumZones = TecUtilDataSetGetNumZones();
				if (isNearestPoint){
					/*
					 *	User selected a node, so activate all FE volumes
					 *	around that node.
					 */
					LgIndex_t NodeNum = TecUtilProbeGetPointIndex();
					int NumFound = 0;
					for (int CurZoneNum = 1; CurZoneNum < NumZones && NumFound < 6; ++CurZoneNum){
						if (TecUtilZoneGetType(CurZoneNum) == ZoneType_FETriangle){
							if (AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
								&& AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.SphereCPName, CPName))
							{
								for (int i = 0; i < 3; ++i){
									if (AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.NodeNums[i], to_string(NodeNum))){
										TecUtilZoneSetActive(Set(CurZoneNum).getRef(), AssignOp_PlusEquals);
										NumFound++;
										break;
									}
								}
							}
						}
					}
				}
				else{
					/*
					*	User selected an element, so activate its FE volume.
					*/
					LgIndex_t ElemNum = TecUtilProbeFieldGetCell();
					for (int CurZoneNum = 1; CurZoneNum < NumZones; ++CurZoneNum){
						if (TecUtilZoneGetType(CurZoneNum) == ZoneType_FETriangle){
							if (AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
								&& AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.SphereCPName, CPName)
								&& AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.ElemNum, to_string(ElemNum)))
							{
								TecUtilZoneSetActive(Set(CurZoneNum).getRef(), AssignOp_PlusEquals);
								break;
							}
						}
					}
				}
			}
			else if (ZoneType == CSMAuxData.GBA.ZoneTypeFEVolumeZone){
				/*
				 *	Two modes:
				 *	1. If user does nearest point probe, then the probed FE
				 *	volume and all volumes that share nodes with it are
				 *	deactivated.
				 *	2. If user does normal probe, only probed FE zone is deactivated.
				 */

				EntIndex_t FEZoneNum = TecUtilProbeFieldGetZone();

				if (isNearestPoint){
					/*
						*	Get node nums for this fe volume zone
						*/
					int NodeNums[3];
					for (int i = 0; i < 3 && IsOk; ++i){
						IsOk = AuxDataZoneGetItem(FEZoneNum, CSMAuxData.GBA.NodeNums[i], TmpStr);
						if (IsOk){
							NodeNums[i] = stoi(TmpStr);
						}
					}
					vector<int> NodeNeighborZones[3];
					for (int i = 0; i < 3; ++i)
						NodeNeighborZones[i].reserve(6);

					Set_pa ActiveZones;

					if (IsOk){
						IsOk = TecUtilZoneGetActive(&ActiveZones);
					}
					if (IsOk){
						IsOk = (TecUtilSetGetMemberCount(ActiveZones) > 0);
					}
					if (IsOk){
						SetIndex_t CurZoneNum = TecUtilSetGetNextMember(ActiveZones, TECUTILSETNOTMEMBER);
						while (IsOk && CurZoneNum != TECUTILSETNOTMEMBER){
							if (CurZoneNum != FEZoneNum && TecUtilZoneIsActive((EntIndex_t)CurZoneNum) && TecUtilZoneGetType((EntIndex_t)CurZoneNum) == ZoneType_FETriangle){
								if (AuxDataZoneItemMatches(static_cast<int>(CurZoneNum), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
									&& AuxDataZoneItemMatches(static_cast<int>(CurZoneNum), CSMAuxData.GBA.SphereCPName, CPName))
								{
									for (int i = 0; i < 3 && IsOk; ++i){
										IsOk = AuxDataZoneGetItem(static_cast<int>(CurZoneNum), CSMAuxData.GBA.NodeNums[i], TmpStr);
										if (IsOk){
											for (int j = 0; j < 3; ++j){
												if (NodeNums[j] == stoi(TmpStr)){
													NodeNeighborZones[j].push_back((int)CurZoneNum);
													break;
												}
											}
										}
									}
								}
							}

							CurZoneNum = TecUtilSetGetNextMember(ActiveZones, CurZoneNum);
						}
					}
					TecUtilSetDealloc(&ActiveZones);
					if (IsOk){
						Set_pa TempSet = TecUtilSetAlloc(FALSE);
						for (int i = 0; i < 3; ++i){
							for (int j = 0; j < NodeNeighborZones[i].size() && IsOk; ++j){
								IsOk = TecUtilSetAddMember(TempSet, NodeNeighborZones[i][j], FALSE);
							}
						}
						IsOk = TecUtilSetAddMember(TempSet, FEZoneNum, FALSE);
						TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
						TecUtilSetDealloc(&TempSet);
					}
				}
				else{
					Set_pa TempSet = TecUtilSetAlloc(FALSE);
					TecUtilSetAddMember(TempSet, FEZoneNum, FALSE);
					TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
					TecUtilSetDealloc(&TempSet);
				}
			}
		}

		TecUtilDataLoadEnd();
	}

	TecUtilDrawGraphics(TRUE);
	TecGUIDialogLaunch(Dialog1Manager);
	TecUtilLockFinish(AddOnID);
}

vec3 GetElemMidPoint(int ZoneNum, int ElemNum){
	int NodeNums[3];
	for (int i = 0; i < 3; ++i) NodeNums[i] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, i + 1);

	vec3 Nodes[3];
	for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) 
		Nodes[i][j] = TecUtilDataValueGetByZoneVar(ZoneNum, j + 1, NodeNums[i]);

	for (int i = 1; i < 3; ++i) Nodes[0] += Nodes[i];
	return Nodes[0] /= 3.0;
}

void SelectGBsInRegion(int const SphereZoneNum,
	int const InteriorElemNum,
	int const GroupNumberToWrite)
{
	/*
	 *	Check that provided zone number is correct and corresponds to a GBA sphere zone
	 */
	if (!AuxDataZoneItemMatches(SphereZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)){
		TecUtilDialogMessageBox("Failed to find source sphere zone", MessageBoxType_Error);
		return;
	}

	/*
	*	Get sphere zone info
	*/
	int SphereIJK[3];
	TecUtilZoneGetIJK(SphereZoneNum, &SphereIJK[0], &SphereIJK[1], &SphereIJK[2]);
	int NumElems = SphereIJK[1];
	string SphereName = AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SphereCPName);

	/*
	*	Get list of all active GBs for the sphere, which are assumed
	*	to be the GBs that define the region of interest.
	*/
	vector<int> ActiveElemNums;
	ActiveElemNums.reserve(NumElems);
	int ElemNum;
	for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z){
		if (TecUtilZoneIsActive(z)
			&& TecUtilZoneIsFiniteElement(z)
			&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
			&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.SphereCPName, SphereName))
		{
			ElemNum = stoi(AuxDataZoneGetItem(z, CSMAuxData.GBA.ElemNum));
			if (ElemNum != InteriorElemNum)
				ActiveElemNums.push_back(ElemNum);
		}
	}

	if (ActiveElemNums.size() <= 0){
		TecUtilDialogMessageBox("Failed to find any active GBs on sphere", MessageBoxType_Error);
		return;
	}

	/*
	 *	Get midpoints of the interior element and boundary elements
	 */
	vec3 InteriorMidPt = GetElemMidPoint(SphereZoneNum, InteriorElemNum);
	vector<vec3> BoundaryMidPts(ActiveElemNums.size());
	for (int i = 0; i < ActiveElemNums.size(); ++i)
		BoundaryMidPts[i] = GetElemMidPoint(SphereZoneNum, ActiveElemNums[i]);

	/*
	 *	Now for each element on the sphere:
	 *	1. find its minimum distance to any of the boundary nodes
	 *	2. find its distance to the interior node
	 *	3. find the vectors from it to the minimum distance boundary node
	 *		and to the interior element
	 *	4. if the distance to the interior node is less than to any of the 
	 *		boundary nodes, or if the dot product between the two vectors is
	 *		less than zero (i.e. they're pointing more than 180 degrees away
	 *		from each other) then consider the element to be interior to the 
	 *		region.
	 */
	vector<double> BoundaryDistSqrList(BoundaryMidPts.size());
	vector<int> ElemsToActivate;
	ElemsToActivate.reserve(NumElems);

	for (int e = 1; e <= NumElems; ++e){
		if (e != InteriorElemNum && std::find(ActiveElemNums.begin(), ActiveElemNums.end(), e) == ActiveElemNums.end()){
			/*
			 *	Get distance (squared) to all boundary elements
			 */
			vec3 ElemMidPt = GetElemMidPoint(SphereZoneNum, e);
			for (int i = 0; i < BoundaryMidPts.size(); ++i)
				BoundaryDistSqrList[i] = DistSqr(ElemMidPt, BoundaryMidPts[i]);

			/*
			 *	Find minimum distance (squared) to boundary nodes
			 */
			double MinBoundaryDistSqr = DBL_MAX;
			int MinBoundaryElemNum = 0;
			for (int i = 0; i < BoundaryDistSqrList.size(); ++i){
				if (BoundaryDistSqrList[i] < MinBoundaryDistSqr){
					MinBoundaryDistSqr = BoundaryDistSqrList[i];
					MinBoundaryElemNum = i;
				}
			}

			/*
			 *	Get vector from element to minimum distance boundary element
			 *	and to interior element.
			 *	Also get distance (squared) to interior element.
			 */
			vec3 BoundaryVec = BoundaryMidPts[MinBoundaryElemNum] - ElemMidPt,
				RefVec = InteriorMidPt - ElemMidPt;
			double RefDistSqr = DistSqr(ElemMidPt, InteriorMidPt);

			if (RefDistSqr < MinBoundaryDistSqr || dot(RefVec, BoundaryVec) < 0)
				ElemsToActivate.push_back(e);
		}	
	}

	ElemsToActivate.push_back(InteriorElemNum);
	ElemsToActivate.insert(ElemsToActivate.end(), ActiveElemNums.begin(), ActiveElemNums.end());

	/*
	 *	Now get set of zones to activate
	 */
	
	vector<int> ZonesToActivate;
	ZonesToActivate.reserve(ElemsToActivate.size());
	Set ZoneSet;

	for (auto const & e : ElemsToActivate){
		for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z){
			if (TecUtilZoneIsFiniteElement(z)
				&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
				&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.SphereCPName, SphereName)
				&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ElemNum, to_string(e)))
			{
				ZonesToActivate.push_back(z);
				ZoneSet += z;
			}
		}
	}

	TecUtilZoneSetActive(ZoneSet.getRef(), AssignOp_PlusEquals);
	
	/*
	 *	Change group number of zones
	 */
	if (GroupNumberToWrite > 0){
		string ZoneStr = "[" + to_string(ZonesToActivate[0]);
		for (int z = 1; z < ZonesToActivate.size(); ++z) ZoneStr += "," + to_string(ZonesToActivate[z]);
		ZoneStr += "]";
		TecUtilMacroExecuteCommand(string("$!FIELDMAP " + ZoneStr + " GROUP = " + to_string(GroupNumberToWrite)).c_str());
	}
}


void STDCALL SelectGBsInRegionProbeCB(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData)
{
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	if (WasSuccessful){
		Boolean_t IsOk = TRUE;

		EntIndex_t ProbedZoneNum = TecUtilProbeFieldGetZone();
		LgIndex_t ElemNum = TecUtilProbeFieldGetCell();
		LgIndex_t GroupNumberToWrite;
		TecGUITextFieldGetLgIndex(TFGrpNum_TF_T3_1, &GroupNumberToWrite);
		SelectGBsInRegion(ProbedZoneNum, ElemNum, GroupNumberToWrite);
	}

	TecUtilDrawGraphics(TRUE);
	TecGUIDialogLaunch(Dialog1Manager);
	TecUtilLockFinish(AddOnID);
}

void ExportGBAData(){
	if ((TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0 && TecGUIToggleGet(TGLExGBs_TOG_T3_1)) || TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0) {
		char* FolderNameCStr;
		if (TecUtilDialogGetFolderName("Select folder to save files", &FolderNameCStr)) {
			string FolderName = FolderNameCStr;
			TecUtilStringDealloc(&FolderNameCStr);
			if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0) {
				if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0) {
					Boolean_t ActiveZonesOnly = TecGUIToggleGet(TGLExGBs_TOG_T3_1);
					vector<int> IntVarNums;
					vector<string> IntVarNames;
					vector<string> IntCheckStrs = { "I: ", "IN: ", "INS: ", " Integration" };
					for (int VarNum = 1; VarNum <= TecUtilDataSetGetNumVars(); ++VarNum) {
						char *VarName, *CheckStr;
						if (TecUtilVarGetName(VarNum, &VarName)) {
							for (string const & Str : IntCheckStrs) {
								CheckStr = std::strstr(VarName, Str.c_str());
								if (CheckStr != nullptr) {
									string TmpStr = VarName;
									std::replace(TmpStr.begin(), TmpStr.end(), ',', '.');
									IntVarNames.push_back(TmpStr);
									IntVarNums.push_back(VarNum);
									break;
								}
								break; // Only want the "I:" integration values
							}
							TecUtilStringDealloc(&VarName);
						}
					}
					int NumZones = TecUtilDataSetGetNumZones();
					for (int ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum) {
						if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)) {
							char* ZoneName;
							if (TecUtilZoneGetName(ZoneNum, &ZoneName)) {
								vector<int> IJK(3);
								TecUtilZoneGetIJK(ZoneNum, &IJK[0], &IJK[1], &IJK[2]);

								vector<bool> ElemActive;
								int NumActive = IJK[1];
								if (ActiveZonesOnly) {
									string SphereName = AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.SphereCPName);
									/*
									*	Get list of elements for which the corresponding GB is active
									*/
									ElemActive.resize(IJK[1], false);
									NumActive = 0;
									for (int e = 0; e < IJK[1]; ++e) {
										for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z) {
											if (TecUtilZoneIsActive(z)
												&& TecUtilZoneIsFiniteElement(z)
												&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
												&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.SphereCPName, SphereName)
												&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ElemNum, to_string(e + 1)))
											{
												ElemActive[e] = true;
												NumActive++;
											}
										}
									}
								}
								else
									ElemActive.resize(IJK[1], true);

								if (NumActive > 0) {
									ofstream OutFile(string(FolderName + "/Zone_" + to_string(ZoneNum) + "_" + string(ZoneName) + "_IntegrationResults.csv").c_str(), std::ios::trunc);
									if (OutFile.is_open()) {
										OutFile << "Zone," << ZoneName << "\nNumber of gradient bundles (GBs)," << IJK[1];
										if (ActiveZonesOnly) {
											OutFile << "\nNumber of active GBs," << NumActive;
										}
										string TmpStr;
										if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.GPsPerGB, TmpStr))
											OutFile << "\nGradient paths (GPs) per GB," << TmpStr;
										if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.PointsPerGP, TmpStr))
											OutFile << "\nPoints per GP," << TmpStr;
										if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.IntPrecision, TmpStr))
											OutFile << "\nIntegration precision," << TmpStr;
										if (AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.IntWallTime, TmpStr))
											OutFile << "\nIntegration runtime total [seconds]," << TmpStr;
										OutFile << "\n\nIntegration totals\nVariable name,Variable number,Value\n";

										vector<FieldDataPointer_c> Ptrs(IntVarNums.size());
										for (int i = 0; i < IntVarNums.size(); ++i)
											Ptrs[i].GetReadPtr(ZoneNum, IntVarNums[i]);

										for (int i = 0; i < Ptrs.size(); ++i) {
											OutFile << IntVarNames[i] << "," << i + 1 << ",";
											double Total = 0.0;
											for (int j = 0; j < Ptrs[i].Size(); ++j) {
												if (ElemActive[j])
													Total += Ptrs[i][j];
											}
											OutFile << std::setprecision(16) << std::scientific << Total << '\n';
										}

										OutFile << "\nZone Name,Zone#,GB#";
										for (string const & i : IntVarNames)
											OutFile << "," << i;
										OutFile << '\n';

										for (int i = 0; i < Ptrs[0].Size(); ++i) {
											if (ElemActive[i]) {
												char* GBZoneName;
												TecUtilZoneGetName(ZoneNum + i + 1, &GBZoneName);
												OutFile << GBZoneName << "," << ZoneNum + i + 1 << "," << i + 1;
												TecUtilStringDealloc(&GBZoneName);
												for (auto const & j : Ptrs)
													OutFile << std::setprecision(16) << std::scientific << "," << j[i];
												OutFile << "\n";
											}
										}

										OutFile.close();
									}
									else
										TecUtilDialogMessageBox(string("Failed to open output file: " + string(FolderName + "/" + string(ZoneName) + "_IntegrationResults.csv")).c_str(), MessageBoxType_Error);
								}
								TecUtilStringDealloc(&ZoneName);
							}
						}
					}

					/*
					 *	Added this for exporting condensed basin and gradient bundle (e.g. bond bundle) information.
					 *	Loop through the zones again saving relevant information for three groups: Bader atoms, 
					 *	special gradient bundles, and condensed wedges.
					 */
					std::set<string> BaderAtomSet, SpecialGradientBundleSet, CondensedBasinSet, VarNameSet;
					std::map<string, vector<string> > CPNames;
					std::map<string, vector<string> > GBSphereZoneNums;
					std::map<string, int> CondensedBasinZoneNums;
					std::map<string, string> DefiningVars;
					std::map<string, std::map<string, double> > VarVals;

					vector<string> OmitStrs = { "IN: ", "INS: " };
					string delStr = "I: ";

					for (int z = 1; z <= NumZones; ++z) {
						string IntVars, IntVals;
						/*
						 *	Checking for bader atoms (GBA spheres)
						 */
						if (AuxDataZoneGetItem(z, CSMAuxData.GBA.AtomicBasinIntegrationValues, IntVals)
							&& AuxDataZoneGetItem(z, CSMAuxData.GBA.AtomicBasinIntegrationVariables, IntVars)){
							char* cstr;
							TecUtilZoneGetName(z, &cstr);
							string tmpStr = cstr;
							BaderAtomSet.insert(tmpStr);
							if (VarVals.count(tmpStr) > 0) {
								VarVals[tmpStr] = std::map<string, double>();
							}
							vector<string> varStrs = SplitString(IntVars, ",");
							vector<double> valVec = SplitStringDbl(IntVals, ",");
							if (varStrs.size() == valVec.size()){
								for (int i = 0; i < varStrs.size(); ++i){

									VarVals[tmpStr][varStrs[i]] = valVec[i];
									VarNameSet.insert(varStrs[i]);
								}
							}
							else {
								TecUtilDialogErrMsg("Discrepancy between number of integration variables and integration values in aux data for atomic basins.");
							}
							GBSphereZoneNums[tmpStr] = { to_string(z) };
							CPNames[tmpStr] = { AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceNucleusName) };
						}

						/*
						 *	Checking for special gradient bundles
						 */
						if (AuxDataZoneHasItem(z, CSMAuxData.GBA.CondensedBasinIsSpecialGradientBundle)
							&& (AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasin)
								|| AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasin))){
							IntVars = AuxDataZoneGetItem(z, CSMAuxData.GBA.IntVarNames);
							IntVals = AuxDataZoneGetItem(z, CSMAuxData.GBA.IntVarVals);

							string tmpStr = AuxDataZoneGetItem(z, CSMAuxData.GBA.CondensedBasinName);
							string defVar = AuxDataZoneGetItem(z, CSMAuxData.GBA.CondensedBasinDefiningVariable);
							tmpStr += GBADelim + defVar;
							DefiningVars[tmpStr] = AuxDataZoneGetItem(z, CSMAuxData.GBA.CondensedBasinDefiningVariable);
							SpecialGradientBundleSet.insert(tmpStr);
							if (CPNames.count(tmpStr) == 0)
								CPNames[tmpStr] = { AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceNucleusName) };
							else
								CPNames[tmpStr].push_back(AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceNucleusName));
							
							if (GBSphereZoneNums.count(tmpStr) == 0)
								GBSphereZoneNums[tmpStr] = { AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceZoneNum) };
							else
								GBSphereZoneNums[tmpStr].push_back(AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceZoneNum));


							vector<string> varStrs = SplitString(IntVars, ",");
							vector<double> valVec = SplitStringDbl(IntVals, ",");
							if (varStrs.size() == valVec.size()) {
								for (int i = 0; i < varStrs.size(); ++i){
									bool doOmit = false;
									for (int si = 0; si < OmitStrs.size() && !doOmit; ++si) {
										doOmit = (varStrs[i].find(OmitStrs[si]) != string::npos);
									}
									if (!doOmit) {
										varStrs[i] = StringReplaceSubString(varStrs[i], delStr, "");
										if (VarVals[tmpStr].count(varStrs[i]) == 0)
											VarVals[tmpStr][varStrs[i]] = valVec[i];
										else
											VarVals[tmpStr][varStrs[i]] += valVec[i];
										VarNameSet.insert(varStrs[i]);
									}
								}
							}
							else {
								TecUtilDialogErrMsg("Discrepancy between number of integration variables and integration values in aux data for special gradient bundles.");
							}

						}
					
						/*
						 *	Checking for condensed basins
						 */
						if (AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasin)
							|| AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasin))
						{
							IntVars = AuxDataZoneGetItem(z, CSMAuxData.GBA.IntVarNames);
							IntVals = AuxDataZoneGetItem(z, CSMAuxData.GBA.IntVarVals);

							string tmpStr = AuxDataZoneGetItem(z, CSMAuxData.GBA.CondensedBasinInfo);
							string defVar = AuxDataZoneGetItem(z, CSMAuxData.GBA.CondensedBasinDefiningVariable);
							tmpStr += GBADelim + defVar;
							CondensedBasinZoneNums[tmpStr] = z;
							CondensedBasinSet.insert(tmpStr);
							CPNames[tmpStr] = { AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceNucleusName) };
							GBSphereZoneNums[tmpStr] = { AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceZoneNum) };
							DefiningVars[tmpStr] = AuxDataZoneGetItem(z, CSMAuxData.GBA.CondensedBasinDefiningVariable);
							
							vector<string> varStrs = SplitString(IntVars, ",");
							vector<double> valVec = SplitStringDbl(IntVals, ",");
							if (varStrs.size() == valVec.size()) {
								for (int i = 0; i < varStrs.size(); ++i) {
									bool doOmit = false;
									for (int si = 0; si < OmitStrs.size() && !doOmit; ++si) {
										doOmit = (varStrs[i].find(OmitStrs[si]) != string::npos);
									}
									if (!doOmit) {
										varStrs[i] = StringReplaceSubString(varStrs[i], delStr, "");
										VarVals[tmpStr][varStrs[i]] = valVec[i];
										VarNameSet.insert(varStrs[i]);
									}
								}
							}
							else {
								TecUtilDialogErrMsg("Discrepancy between number of integration variables and integration values in aux data for condensed basins.");
							}
						}
					}

					/*
					 *	Now write it all out.
					 */
					char *cstr;
					TecUtilDataSetGetInfo(&cstr, nullptr, nullptr);
					string dataSetName = cstr;
					TecUtilStringDealloc(&cstr);
					ofstream OutFile(string(FolderName + "/" + dataSetName + "_Atomic_Basins.csv"));

					/*
					 *`Get current date
					 */
					auto t = std::time(nullptr);
					auto tm = *std::localtime(&t);
					stringstream dateSS;
					dateSS << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << std::endl;

					// First atomic basins
					if (OutFile.is_open()){
						/*
						 *	Print date and dataset name
						 */
						OutFile << dateSS.str() << endl << dataSetName << endl;

						/*
						 *	Print headings
						 */
						OutFile << "Atom,Zone name,Zone number,";
						// and the variable names
						for (string const & s : VarNameSet)
							OutFile << s << ",";
						OutFile << endl;

						/*
						 *	Now all the values
						 */
						for (string const & abStr : BaderAtomSet){
							OutFile << abStr << "," << CPNames[abStr][0] << "," << GBSphereZoneNums[abStr][0] << ",";
							for (string const & varStr : VarNameSet){
								if (VarVals[abStr].count(varStr) > 0)
									OutFile << std::setprecision(16) << VarVals[abStr][varStr];

								OutFile << ",";
							}
							OutFile << endl;
						}
						OutFile.close();
					}
					else {
						TecUtilDialogErrMsg("Failed to open output file for writing atomic basin data");
					}

					// Then special gradient bundles
					OutFile.open(string(FolderName + "/" + dataSetName + "_Special_Gradient_Bundles.csv"));
					if (OutFile.is_open()) {
						/*
							*	Print date and dataset name
							*/
						OutFile << dateSS.str() << endl << dataSetName << endl;

						/*
							*	Print headings
							*/
						OutFile << "Special gradient bundle,Contributing atoms,Atom sphere zone numbers,Defined using variable,";
						// and the variable names
						for (string const & s : VarNameSet)
							OutFile << s << ",";
						OutFile << endl;

						/*
							*	Now all the values
							*/
						for (string const & sgpStr : SpecialGradientBundleSet) {
							OutFile << SplitString(sgpStr, GBADelim)[0] << "," << VectorToString(CPNames[sgpStr], "-") << "," << VectorToString(GBSphereZoneNums[sgpStr],"-") << "," << DefiningVars[sgpStr] << ",";
							for (string const & varStr : VarNameSet) {
								if (VarVals[sgpStr].count(varStr) > 0)
									OutFile << std::setprecision(16) << VarVals[sgpStr][varStr];

								OutFile << ",";
							}
							OutFile << endl;
						}
						OutFile.close();
					}
					else {
						TecUtilDialogErrMsg("Failed to open output file for writing special gradient bundle data");
					}

					// Then the condensed basins, which includes any that were also counted in the special gradient bundles
					OutFile.open(string(FolderName + "/" + dataSetName + "_Condensed_Basins.csv"));
					if (OutFile.is_open()) {
						/*
							*	Print date and dataset name
							*/
						OutFile << dateSS.str() << endl << dataSetName << endl;

						/*
							*	Print headings
							*/
						OutFile << "Atom,Condensed basin information,Zone number,Atom sphere zone number,Defined using variable,";
						// and the variable names
						for (string const & s : VarNameSet)
							OutFile << s << ",";
						OutFile << endl;

						/*
							*	Now all the values
							*/
						for (string const & cbStr : CondensedBasinSet) {
							OutFile << CPNames[cbStr][0] << "," << SplitString(cbStr, GBADelim)[0] << "," << CondensedBasinZoneNums[cbStr] << GBSphereZoneNums[cbStr][0] << "," << DefiningVars[cbStr] << ",";
							for (string const & varStr : VarNameSet) {
								if (VarVals[cbStr].count(varStr) > 0)
									OutFile << std::setprecision(16) << VarVals[cbStr][varStr];

								OutFile << ",";
							}
							OutFile << endl;
						}
					}
					else {
						TecUtilDialogErrMsg("Failed to open output file for writing special gradient bundle data");
					}
				}
			}
		}
	}
}