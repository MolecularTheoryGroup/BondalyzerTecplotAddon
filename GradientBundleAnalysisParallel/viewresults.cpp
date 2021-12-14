
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
#include <queue>
#include <ctime>

#include <armadillo>

#include "Set.h"
#include "StyleValue.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GUI.h"
#include "VIEWRESULTS.h"
#include "CSM_GBAGUI.h"

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

	ResultsVarListReload();

	return;
}

void GBAResultViewerSelectIntVar(){

	Boolean_t IsOk = TRUE;

	LgIndex_t *SelectedNums;
	LgIndex_t NumSelected;

	TecGUIListGetSelectedItems(SLSelVar_SLST_T3_1, &SelectedNums, &NumSelected);

	EntIndex_t VarNum = -1;
	EntIndex_t SphereZoneNum = -1;
	int NumContours = stoi(TecGUITextFieldGetString(TFNumContours_TF_T3_1));
	bool LogSpaceContours = (TecGUIRadioBoxGetToggle(RBLogLin_RADIO_T3_1) == 1);
	bool UseSelectedSphere = (TecGUIRadioBoxGetToggle(RBCntSrc_RADIO_T3_1) == 1);
	bool INSOnly = TecGUIToggleGet(TGLINSV_TOG_T3_1);
	vec ContourLevelsVec;
	string SphereNameStr;
	if (NumSelected > 0){
		char *IntVarNameCStr = TecGUIListGetString(SLSelVar_SLST_T3_1, TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1));
		string IntVarNameStr = IntVarNameCStr;
		TecUtilStringDealloc(&IntVarNameCStr);

		bool HasINS = false;
		for (auto & s : { "I: ","IN: ", "INS: " }){
			if (IntVarNameStr.find(s) != string::npos){
				HasINS = true;
				break;
			}
		}
		if (!HasINS){
			if (IntVarNameStr.find("Average ") == string::npos) {
				IntVarNameStr = "INS: " + IntVarNameStr;
			}
			else{
				IntVarNameStr = "I: " + IntVarNameStr;
			}
		}

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
			return;
		}

		StyleValue styleValue;
		styleValue.set((Boolean_t)0, SV_GLOBALLINKING, SV_LINKCOLORMAPS);


		styleValue.set(8, VarNum, SV_GLOBALCONTOUR, SV_COLORMAPFILTER, SV_GROUP);

		if (VarMin >= 0. || VarMax <= 0.) {
			styleValue.set(28, 8, SV_GLOBALCOLORMAP, SV_CONTOURCOLORMAP);
			styleValue.set((Boolean_t)1, 8, SV_GLOBALCONTOUR, SV_COLORMAPFILTER, SV_REVERSECOLORMAP);
		}
		else {
			styleValue.set(8, 8, SV_GLOBALCOLORMAP, SV_CONTOURCOLORMAP);
			styleValue.set((Boolean_t)0, 8, SV_GLOBALCONTOUR, SV_COLORMAPFILTER, SV_REVERSECOLORMAP);
		}

		if (LogSpaceContours){
			if (VarMin <= 0.0){
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
							else if (Val > 0.0 && Val < PosCloseToZero)
								PosCloseToZero = Val;
						}
					}
				}
				if (NegCloseToZero == -DBL_MAX)
					NegCloseToZero = -1e-3;
				if (PosCloseToZero == DBL_MAX)
					PosCloseToZero = 1e-3;

				if (VarMax <= 0.){
					if (VarMax == 0.){
						ContourLevelsVec = join_cols(LogSpace(-0.95 * VarMin, -1e-2 * VarMin, NumContours - 1) * -1.0, vec({ 0. }));
					}
					else {
						ContourLevelsVec = LogSpace(-VarMax, -VarMin, NumContours) * -1.0;
					}
				}
				else if (VarMin >= 0.){
					if (VarMin == 0.){
						ContourLevelsVec = join_cols(vec({ 0. }), LogSpace(VarMax * 1e-5, VarMax, NumContours - 1));
					}
					else {
						ContourLevelsVec = LogSpace(VarMin, VarMax, NumContours);
					}
				}
				else{
					int NumNegContours = int(double(NumContours) * (-VarMin) / (VarMax - VarMin));

					ContourLevelsVec = join_cols(join_cols(flipud(LogSpace(-VarMin, -NegCloseToZero, 3) * -1.0), vec({ 0. })), flipud(LogSpace(VarMax, PosCloseToZero, 3)));
					NumContours = 7;
				}
			}
			else{
				ContourLevelsVec = LogSpace(VarMin, VarMax, NumContours);
			}
		}
		else{
			if (VarMin >= 0. || VarMax <= 0.) {
				ContourLevelsVec = linspace(VarMin, VarMax, NumContours);
			}
			else {
				ContourLevelsVec = join_cols(join_cols(flipud(linspace(-VarMin * 0.2, -VarMin * 0.8, 3) * -1.0), vec({ 0.0 })), linspace(VarMax * 0.2, VarMax * 0.8, 3));
				NumContours = 7;
			}
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
// 			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
// 				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr)
// 				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)) {
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)) {
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
}

void GBAResultViewerSelectCondensedGBs() {

	Boolean_t IsOk = TRUE;

	LgIndex_t *SelectedNums;
	LgIndex_t NumSelected;

	TecGUIListGetSelectedItems(MLSelGB_MLST_T3_1, &SelectedNums, &NumSelected);

	Set ActivateSet, DeactivateSet;
	std::set<int> ActivateSetSet, DeactivateSetSet;
	

	if (NumSelected > 0) {

		// 		TecUtilDataLoadBegin();

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		bool SurfacePresent = false;
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
						if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasinWedge)
							|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasinWedge)){
							SurfacePresent = true;
						}
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
		// turn on truncating isosurface if present
		if (SurfacePresent) {
			int isoZone = ZoneNumByName("gba-truncating-isosurface_", false, true);
			if (isoZone > 0) {
				TecUtilZoneSetActive(Set(isoZone).getRef(), AssignOp_PlusEquals);
			}
		}

		// 		TecUtilDataLoadEnd();
	}
	else{// if (TecGUIListGetItemCount(MLSelGB_MLST_T3_1) == 0) {
		// 		TecUtilDataLoadBegin();

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum) {
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasinWedge)
				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasinWedge)
				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasin)
				|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasin))
			{
				DeactivateSet += ZoneNum;
			}
		}

		if (IsOk) {
			TecUtilZoneSetActive(DeactivateSet.getRef(), AssignOp_MinusEquals);
		}

		// turn off truncating isosurface if present
		int isoZone = ZoneNumByName("gba-truncating-isosurface_", false, true);
		if (isoZone > 0) {
			TecUtilZoneSetActive(Set(isoZone).getRef(), AssignOp_MinusEquals);
		}

		// 		TecUtilDataLoadEnd();
	}

	TecUtilArrayDealloc((void**)&SelectedNums);
}

void GBAResultViewerSelectGB(){

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

}

void GBAResultViewerToggleSphere(){

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
	GBAResultViewerSelectSphere();
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

	if (TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0){
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

	}
}

void GBAResultViewerActivateAllGB(){

	if (TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0){
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

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB)
				&& (
				AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)
								|| AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SourceNucleusName, SphereNameStr))
				)
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
			LgIndex_t *selected, numselected;
			TecGUIListGetSelectedItems(MLSelGB_MLST_T3_1, &selected, &numselected);
			if (numselected == TecGUIListGetItemCount(MLSelGB_MLST_T3_1)){
				ListDeselect(MLSelGB_MLST_T3_1);
				GBAResultViewerSelectCondensedGBs();
			}
			else if (NumGBAtoms > 0) {
				vector<int> SelNums(NumGBAtoms);
				for (int i = 1; i <= NumGBAtoms; ++i)
					SelNums[i - 1] = i;

				TecGUIListSetSelectedItems(MLSelGB_MLST_T3_1, SelNums.data(), NumGBAtoms);
			}
		}

		TecUtilSetDealloc(&ActivateSet);

		TecUtilStringDealloc(&SphereNameCStr);
	}
}

/*
 *	Volume zone toggle probe functions
 */

std::map<unsigned int, bool> SphereZoneHasSavedGBs;

void SphereZoneCheckForSavedGBs(){
	SphereZoneHasSavedGBs.clear();
	int NumZones = TecUtilDataSetGetNumZones();
	string SphereName;
	for (int z = 1; z <= NumZones; ++z){
		if (AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)
			&& AuxDataZoneGetItem(z, CSMAuxData.GBA.SphereCPName, SphereName))
		{
			bool IsFound = false;
			for (int zgb = z + 1; zgb <= NumZones && !IsFound; ++zgb){
				if (AuxDataZoneItemMatches(zgb, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB)
					&& AuxDataZoneItemMatches(zgb, CSMAuxData.GBA.SphereCPName, SphereName)){
					IsFound = true;
				}
			}
			SphereZoneHasSavedGBs[z] = IsFound;
		}
	}
}

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

void MakeTriangularSphereElement(int ZoneNum, int NodeNum = -1, int ElemNum = -1){ // node and element base 0
	REQUIRE(ZoneNum > 0 && ZoneNum <= TecUtilDataSetGetNumZones());
	REQUIRE(NodeNum < 0 || ElemNum < 0);
	REQUIRE(AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone));

	FESurface_c Sphere = FESurface_c(ZoneNum,  { 1,2,3 });
	auto XYZPtr = Sphere.GetXYZListPtr();
	auto ElemPtr = Sphere.GetElemListPtr();

	vector<int> ElemNums;
	if (ElemNum >= 0) {
		REQUIRE(ElemNum < ElemPtr->size());
		ElemNums.push_back(ElemNum);
	}
	else{
		REQUIRE(NodeNum < XYZPtr->size());
		for (int ei = 0; ei < ElemPtr->size(); ++ei){
			for (auto c : ElemPtr->at(ei)){
				if (c == NodeNum){
					ElemNums.push_back(ei);
					break;
				}
			}
		}
	}

	// get info for new zone(s)
	string CPString = AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.SphereCPName), 
		NucleusName = AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.SourceNucleusName);
	vec3 SphereOrigin = vec(SplitStringDbl(AuxDataZoneGetItem(ZoneNum, CSMAuxData.GBA.SphereOrigin)));

	vector<vector<int> > Elems = { {0,1,2} };
	for (auto ei : ElemNums){
		vector<vec3> Nodes;
		for (auto ni : ElemPtr->at(ei)){
			Nodes.emplace_back(SphereOrigin + (XYZPtr->at(ni) - SphereOrigin) * 1.001);
		}
		FESurface_c GradientBundle(Nodes, Elems);

		GradientBundle.SaveAsTriFEZone({ 1,2,3 }, CPString + ": Gradient Bundle " + to_string(ei + 1));

		if (GradientBundle.IsMade() && GradientBundle.GetZoneNum() > 0) {
			TecUtilZoneSetActive(Set(GradientBundle.GetZoneNum()).getRef(), AssignOp_PlusEquals);

			TecUtilZoneSetShade(SV_SHOW, Set(GradientBundle.GetZoneNum()).getRef(), 0.0, TRUE);
			TecUtilZoneSetShade(SV_COLOR, Set(GradientBundle.GetZoneNum()).getRef(), 0.0, Cyan_C);

			AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.SourceZoneNum, to_string(ZoneNum));
			AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.SphereCPName, CPString);
			AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB);
			AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.ElemNum, to_string(ei + 1));
			AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.SourceNucleusName, NucleusName);

			for (int i = 0; i < 3; ++i)
				AuxDataZoneSetItem(GradientBundle.GetZoneNum(), CSMAuxData.GBA.NodeNums[i], to_string(ElemPtr->at(ei)[i] + 1));
		}
		else {
			TecUtilDialogErrMsg(string("Failed to save GB " + to_string(ei + 1)).c_str());
		}
	}
}

void STDCALL ToggleFEVolumesProbeCB(Boolean_t WasSuccessful,
	Boolean_t isNearestPoint,
	ArbParam_t ClientData)
{
	TecUtilLockStart(AddOnID);
	CSMGUILock();
	if (WasSuccessful){

		string TmpStr1, TmpStr2, TmpStr3;

		Boolean_t IsOk = TRUE;

		EntIndex_t ProbedZoneNum = TecUtilProbeFieldGetZone();

		string TmpStr;

		string CPName;

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
							if (AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB)
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

					// No zone found, so make one.
					MakeTriangularSphereElement(ProbedZoneNum, NodeNum - 1);
				}
				else{
					/*
					*	User selected an element, so activate its FE volume.
					*/
					LgIndex_t ElemNum = TecUtilProbeFieldGetCell();
					bool IsFound = false;
					for (int CurZoneNum = 1; CurZoneNum < NumZones; ++CurZoneNum){
						if (TecUtilZoneGetType(CurZoneNum) == ZoneType_FETriangle){
							if (AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB)
								&& AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.SphereCPName, CPName)
								&& AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.ElemNum, to_string(ElemNum)))
							{
								TecUtilZoneSetActive(Set(CurZoneNum).getRef(), AssignOp_PlusEquals);
								IsFound = true;
								break;
							}
						}
					}

					if (!IsFound){
						// No zone found, so make one.
						MakeTriangularSphereElement(ProbedZoneNum, -1, ElemNum - 1);
					}
				}
			}
			else if (ZoneType == CSMAuxData.GBA.ZoneTypeDGB){
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
								if (AuxDataZoneItemMatches(static_cast<int>(CurZoneNum), CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB)
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
	}

	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);
	TecGUIDialogLaunch(Dialog1Manager);
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
			&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB)
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
	 *	
	 *	Actually this method is terrible.
	 *	
	 *	Instead, start a breadth first search from the interior element.
	 *	When an element that has neighboring elements that are on the boundary,
	 *	only add its neighbors to the search queue if their distance to the 
	 *	interior element is less than the average distance of neighboring boundary
	 *	nodes to the interior element.
	 */

	FESurface_c Sphere(SphereZoneNum, { 1,2,3 });
	vec3 SphereOrigin = vec(SplitStringDbl(AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SphereOrigin)));
	Sphere.GenerateElemConnectivity();
	Sphere.GenerateElemMidpoints();
	auto ElemMidPoints = Sphere.GetElemMidpointsPtr();
	auto ElemConnectivity = Sphere.GetElemConnectivityListPtr();
	std::map<int, double> BoundaryElemNumToDistMap; // 0-based
	std::set<int> VisitedElems;
	for (auto i : ActiveElemNums){
// 		BoundaryElemNumToDistMap[i - 1] = Distance(ElemMidPoints->at(i - 1), ElemMidPoints->at(InteriorElemNum - 1));
		BoundaryElemNumToDistMap[i - 1] = VectorAngle(ElemMidPoints->at(i - 1) - SphereOrigin, ElemMidPoints->at(InteriorElemNum - 1) - SphereOrigin);
		VisitedElems.insert(i - 1);
	}

	std::queue<int> BFSQueue;
	BFSQueue.push(InteriorElemNum - 1);
	BFSQueue.push(-1);
	VisitedElems.insert(InteriorElemNum - 1);
	auto InteriorElems = VisitedElems;
	vector<double> NeighborBoundaryInteriorDists;
	NeighborBoundaryInteriorDists.reserve(20);

	while (!BFSQueue.empty()){
		auto ei = BFSQueue.front();
		BFSQueue.pop();

		if (ei < 0){
			BFSQueue.push(-1);
			if (BFSQueue.size() == 1){
				break;
			}
			continue;
		}

		for (auto ej : ElemConnectivity->at(ei)){
			if (!VisitedElems.count(ej)){ 
				VisitedElems.insert(ej);
				NeighborBoundaryInteriorDists.clear();
				for (auto ek : ElemConnectivity->at(ej)){
					if (ej != ek && BoundaryElemNumToDistMap.count(ek)){
						NeighborBoundaryInteriorDists.push_back(BoundaryElemNumToDistMap[ek]);
					}
				}
				if (NeighborBoundaryInteriorDists.size() > 0){
					// if (Distance(ElemMidPoints->at(ej), ElemMidPoints->at(InteriorElemNum - 1)) <= max(vec(NeighborBoundaryInteriorDists))) {
					if (VectorAngle(ElemMidPoints->at(ej) - SphereOrigin, ElemMidPoints->at(InteriorElemNum - 1) - SphereOrigin) <= max(vec(NeighborBoundaryInteriorDists))) {
						BFSQueue.push(ej);
						InteriorElems.insert(ej);
					}
				}
				else{
					BFSQueue.push(ej);
					InteriorElems.insert(ej);
				}
			}
		}
	}


	


// 
// 	/*
// 	 *	Get midpoints of the interior element and boundary elements
// 	 */
// 
// 	vec3 InteriorMidPt = GetElemMidPoint(SphereZoneNum, InteriorElemNum);
// 	vector<vec3> BoundaryMidPts(ActiveElemNums.size());
// 	for (int i = 0; i < ActiveElemNums.size(); ++i)
// 		BoundaryMidPts[i] = GetElemMidPoint(SphereZoneNum, ActiveElemNums[i]);
// 
// 	vector<double> BoundaryDistSqrList(BoundaryMidPts.size());
// 	vector<int> ElemsToActivate;
// 	ElemsToActivate.reserve(NumElems);
// 
// 	for (int e = 1; e <= NumElems; ++e){
// 		if (e != InteriorElemNum && std::find(ActiveElemNums.begin(), ActiveElemNums.end(), e) == ActiveElemNums.end()){
// 			/*
// 			 *	Get distance (squared) to all boundary elements
// 			 */
// 			vec3 ElemMidPt = GetElemMidPoint(SphereZoneNum, e);
// 			for (int i = 0; i < BoundaryMidPts.size(); ++i)
// 				BoundaryDistSqrList[i] = DistSqr(ElemMidPt, BoundaryMidPts[i]);
// 
// 			/*
// 			 *	Find minimum distance (squared) to boundary nodes
// 			 */
// 			double MinBoundaryDistSqr = DBL_MAX;
// 			int MinBoundaryElemNum = 0;
// 			for (int i = 0; i < BoundaryDistSqrList.size(); ++i){
// 				if (BoundaryDistSqrList[i] < MinBoundaryDistSqr){
// 					MinBoundaryDistSqr = BoundaryDistSqrList[i];
// 					MinBoundaryElemNum = i;
// 				}
// 			}
// 
// 			/*
// 			 *	Get vector from element to minimum distance boundary element
// 			 *	and to interior element.
// 			 *	Also get distance (squared) to interior element.
// 			 */
// 			vec3 BoundaryVec = BoundaryMidPts[MinBoundaryElemNum] - ElemMidPt,
// 				RefVec = InteriorMidPt - ElemMidPt;
// 			double RefDistSqr = DistSqr(ElemMidPt, InteriorMidPt);
// 
// 			if (RefDistSqr < MinBoundaryDistSqr || dot(RefVec, BoundaryVec) < 0)
// 				ElemsToActivate.push_back(e);
// 		}	
// 	}
// 
// 	ElemsToActivate.push_back(InteriorElemNum);
// 	ElemsToActivate.insert(ElemsToActivate.end(), ActiveElemNums.begin(), ActiveElemNums.end());

	vector<int> ElemsToActivate;
	for (auto ei : InteriorElems){
		ElemsToActivate.push_back(ei + 1);
	}

	/*
	 *	Now get set of zones to activate
	 */
	
	vector<int> ZonesToActivate;
	ZonesToActivate.reserve(ElemsToActivate.size());
	Set ZoneSet;

	for (auto const & e : ElemsToActivate){
		bool isFound = false;
		for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z){
			if (TecUtilZoneIsFiniteElement(z)
				&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB)
				&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.SphereCPName, SphereName)
				&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ElemNum, to_string(e)))
			{
				ZonesToActivate.push_back(z);
				isFound = true;
				ZoneSet += z;
				break;
			}
		}

		if (!isFound) {
			// Make new element zones
			MakeTriangularSphereElement(SphereZoneNum, -1, e - 1);
			ZonesToActivate.push_back(TecUtilDataSetGetNumZones());
			ZoneSet += TecUtilDataSetGetNumZones();
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
	CSMGUILock();

	if (WasSuccessful){
		Boolean_t IsOk = TRUE;

		EntIndex_t ProbedZoneNum = TecUtilProbeFieldGetZone();
		LgIndex_t ElemNum = TecUtilProbeFieldGetCell();
		LgIndex_t GroupNumberToWrite;
		TecGUITextFieldGetLgIndex(TFGrpNum_TF_T3_1, &GroupNumberToWrite);
		SelectGBsInRegion(ProbedZoneNum, ElemNum, GroupNumberToWrite);
	}

	CSMGUIUnlock();
	TecUtilLockFinish(AddOnID);
	TecGUIDialogLaunch(Dialog1Manager);
}

/*
 *	Export gradient bundle integration data for a variety of gradient bundles for every atom in the system for which GBA has been run:
 *	* Individual differential gradient bundles
 *		* if they or their representative elements are active, a separate CSV file with only active dGBs
 *		* a CSV with all dGBs
 *		* a CSV with dGBs for each condensed maximum/minimum basin
 *		* a CSV with dGBs for each special gradient bundle (ie. for a "bond" one or more bond wedges present and contributing to the bond)
 *		
 *	* Integration totals for for each condensed maximum/minimum basin for each sphere
 */
void ExportGBAData(){
	if (TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0 || TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0) {
		char* FolderNameCStr;
		vector<int> IntVarNums;
		vector<string> IntVarNames;
		if (TecUtilDialogGetFolderName("Select folder to save files", &FolderNameCStr)) {
			string FolderName = FolderNameCStr;
			TecUtilStringDealloc(&FolderNameCStr);
			if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0) {
				if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0) {
					Boolean_t IncludeAllGBs = FALSE; // !TecGUIToggleGet(TGLExGBs_TOG_T3_1);
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
// 								break; // Only want the "I:" integration values
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
								if (IncludeAllGBs) {
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
												&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeDGB)
												&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.SphereCPName, SphereName)
												&& AuxDataZoneItemMatches(z, CSMAuxData.GBA.ElemNum, to_string(e + 1)))
											{
												ElemActive[e] = true;
												NumActive++;
											}
										}
									}
								}
								else {
									ElemActive.resize(IJK[1], true);
								}

								if (!IncludeAllGBs || NumActive > 0) {
									ofstream OutFile(string(FolderName + "/Zone_" + to_string(ZoneNum) + "_" + string(ZoneName) + "_IntegrationResults.csv").c_str(), std::ios::trunc);
									if (OutFile.is_open()) {
										OutFile << "Zone," << ZoneName << "\nNumber of gradient bundles (GBs)," << IJK[1];
										if (IncludeAllGBs) {
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
											Ptrs[i].InitializeReadPtr(ZoneNum, IntVarNums[i]);

										for (int i = 0; i < Ptrs.size(); ++i) {
											OutFile << IntVarNames[i] << "," << i + 1 << ",";
											double Total = 0.0;
											for (int j = 0; j < Ptrs[i].Size(); ++j) {
												if (!IncludeAllGBs || ElemActive[j])
													Total += Ptrs[i][j];
											}
											OutFile << std::setprecision(16) << std::scientific << Total << '\n';
										}


										if (IncludeAllGBs) {
											OutFile << "\nGB number";
											for (string const & i : IntVarNames)
												OutFile << "," << i;
											OutFile << '\n';
											for (int i = 0; i < Ptrs[0].Size(); ++i) {
												if (ElemActive[i]) {
// 													int GBZoneNum = ZoneNum + i + 1;
// 													if (GBZoneNum >= 1 && GBZoneNum <= NumZones) {
// 														char* GBZoneName;
// 														TecUtilZoneGetName(ZoneNum + i + 1, &GBZoneName);
// 														OutFile << GBZoneName << "," << ZoneNum + i + 1 << "," << i + 1;
// 														TecUtilStringDealloc(&GBZoneName);
// 														for (auto const & j : Ptrs)
// 															OutFile << std::setprecision(16) << std::scientific << "," << j[i];
// 														OutFile << "\n";
// 													}
													OutFile << i + 1;
													for (auto const & j : Ptrs)
														OutFile << std::setprecision(16) << std::scientific << "," << j[i];
													OutFile << "\n";
												}
											}
										}
										else{
											OutFile << "\nGB number";
											for (string const & i : IntVarNames)
												OutFile << "," << i;
											OutFile << '\n';
											for (int i = 0; i < Ptrs[0].Size(); ++i) {
												OutFile << i + 1;
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
							OutFile << CPNames[cbStr][0] << "," << SplitString(cbStr, GBADelim)[0] << "," << CondensedBasinZoneNums[cbStr] << "," << GBSphereZoneNums[cbStr][0] << "," << DefiningVars[cbStr] << ",";
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

bool GetSphereOrigin(int SphereZoneNum, vec3 & Origin) {
	if (AuxDataZoneItemMatches(SphereZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone))
	{
		bool OriginFound = false;
		// If origin is saved to sphere, just use that
		string OriginString;
		if (AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SphereOrigin, OriginString)) {
			vector<double> OriginVec = SplitStringDbl(OriginString);
			if (OriginVec.size() == 3) {
				for (int i = 0; i < 3; ++i) {
					Origin[i] = OriginVec[i];
				}
				OriginFound = true;
			}
		}

		if (OriginFound)
			return OriginFound;

		vector<int> XYZVarNums(3);
		TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);

		// now see if origin CP zone is present and use that
		string CPZoneNumStr, CPNumStr;
		if (AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SourceZoneNum, CPZoneNumStr) && AuxDataZoneGetItem(SphereZoneNum, CSMAuxData.GBA.SphereCPNum, CPNumStr)) {
			int VolZoneNum = ZoneNumByName("Full Volume", false, true);
			int SourceCPZoneNum = stoi(CPZoneNumStr),
				SourceCPNum = stoi(CPNumStr);
			if (VolZoneNum <= 0){
				SourceCPZoneNum -= 1;
			}
			if (AuxDataZoneItemMatches(SourceCPZoneNum, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.CPSubTypes[0])) {
				for (int i = 0; i < 3; ++i) {
					Origin[i] = TecUtilDataValueGetByZoneVar(SourceCPZoneNum, XYZVarNums[i], SourceCPNum);
				}
				OriginFound = true;
			}
		}

		if (OriginFound)
			return OriginFound;


		// Just use nuclear CP closest to sphere node midpoint, since changes to zones can break the commented code below.
		FieldVecPointer_c SphereXYZ;
		SphereXYZ.InitializeReadPtr(SphereZoneNum, XYZVarNums);
		vec3 SphereMaxXYZ = vec3() * DBL_MIN,
			SphereMinXYZ = vec3() * DBL_MAX;
		vec3 SphereMidPoint = { 0,0,0 };
		for (int i = 0; i < SphereXYZ.Size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				SphereMaxXYZ[j] = MAX(SphereMaxXYZ[j], SphereXYZ[i][j]);
				SphereMinXYZ[j] = MAX(SphereMinXYZ[j], SphereXYZ[i][j]);
			}
		}
		SphereMidPoint = (SphereMaxXYZ + SphereMinXYZ) * 0.5;

		// 		SaveVec3VecAsScatterZone({ SphereMidPoint }, "sphere midpoint");

				// Now get closest nuclear CP
		double MinDistSqr = DBL_MAX;
		for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z) {
			if (AuxDataZoneItemMatches(z, CSMAuxData.CC.ZoneSubType, CSMAuxData.CC.CPSubTypes[0])) {
				FieldVecPointer_c CPsXYZ;
				CPsXYZ.InitializeReadPtr(z, XYZVarNums);
				for (int i = 0; i < CPsXYZ.Size(); ++i) {
					double TmpDistSqr = DistSqr(CPsXYZ[i], SphereMidPoint);
					if (TmpDistSqr < MinDistSqr) {
						MinDistSqr = TmpDistSqr;
						Origin = CPsXYZ[i];
					}
				}
			}
		}

		if (Distance(SphereMidPoint, Origin) > 0.1) {
			// too far from the sphere midpoint, so use the midpoint instead.
			Origin = SphereMidPoint;
		}

		return true;
	}
	return false;
}

void ResizeSpheres(double const & SizeFactor, Boolean_t AllSpheres, Boolean_t AbsoluteRadius){
	CSMGUILock();

	vector<EntIndex_t> XYZVarNums(3);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);

	string SelectedSphereName = TecGUIListGetString(SLSelSphere_SLST_T3_1, TecGUIListGetSelectedItem(SLSelSphere_SLST_T3_1));

	for (int z = 1; z <= TecUtilDataSetGetNumZones(); ++z){
		if (AuxDataZoneItemMatches(z, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeSphereZone)
			&& (AllSpheres || AuxDataZoneItemMatches(z, CSMAuxData.GBA.SphereCPName, SelectedSphereName)))
		{
			vec3 Origin;
			GetSphereOrigin(z, Origin);
// 			string CPZoneNumStr, CPNumStr;
// 			if (AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceZoneNum, CPZoneNumStr) && AuxDataZoneGetItem(z, CSMAuxData.GBA.SphereCPNum, CPNumStr)){
// 				int SourceCPZoneNum = stoi(CPZoneNumStr),
// 					SourceCPNum = stoi(CPNumStr);
// 				for (int i = 0; i < 3; ++i){
// 					Origin[i] = TecUtilDataValueGetByZoneVar(SourceCPZoneNum, XYZVarNums[i], SourceCPNum);
// 				}
// 			}
// 			else{
// 				TecUtilDialogErrMsg(string("Couldn't get source CP coordinates for sphere zone " + to_string(z)).c_str());
// 				continue;
// 			}

			// Get list of basin zone numbers
			string SphereName = AuxDataZoneGetItem(z, CSMAuxData.GBA.SourceNucleusName);
			vector<int> ZoneNums;
			for (int zi = z + 1; zi <= TecUtilDataSetGetNumZones(); ++zi){
				if ((AuxDataZoneItemMatches(zi, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedAttractiveBasin)
					|| AuxDataZoneItemMatches(zi, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeCondensedRepulsiveBasin)
					|| AuxDataZoneItemMatches(zi, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeTopoCageWedge))
					&& AuxDataZoneItemMatches(zi, CSMAuxData.GBA.SourceNucleusName, SphereName))
				{
					ZoneNums.push_back(zi);
				}
			}

			double NewRadius = SizeFactor;
			string OldRadStr;
			double OldRadius;
			if (!AbsoluteRadius && AuxDataZoneGetItem(z, CSMAuxData.GBA.RadialApprxRadius, OldRadStr)) {
				OldRadius = stof(OldRadStr);
				NewRadius *= OldRadius;
			}

			for (auto zi : ZoneNums) {
				FieldVecPointer_c XYZPtr;
				TecUtilDataLoadBegin();
				XYZPtr.InitializeWritePtr(zi, XYZVarNums);
				int NumNodes = XYZPtr.Size();

#pragma omp parallel for
				for (int i = 0; i < NumNodes; ++i) {
					XYZPtr.Write(i, Origin + normalise(XYZPtr[i] - Origin) * NewRadius * 1.01);
				}
				TecUtilDataLoadEnd();
			}

			FieldVecPointer_c XYZPtr;
			TecUtilDataLoadBegin();
			XYZPtr.InitializeWritePtr(z, XYZVarNums);
			int NumNodes = XYZPtr.Size();

#pragma omp parallel for
			for (int i = 0; i < NumNodes; ++i) {
				XYZPtr.Write(i, Origin + normalise(XYZPtr[i] - Origin) * NewRadius);
			}
			TecUtilDataLoadEnd();
		}
	}

	CSMGUIUnlock();
}


void ResultsVarListReload() {
	bool INSOnly = TecGUIToggleGet(TGLINSV_TOG_T3_1);
	string CurVal = "", CurValShort;
	if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0) {
		int CurSel = TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1);
		CurVal = TecGUIListGetString(SLSelVar_SLST_T3_1, CurSel);
		for (auto & s : { "I: ","IN: ","INS: " }) {
			CurValShort = StringRemoveSubString(CurVal, s);
		}
	}
	if (CurVal.find("Average ") == string::npos) {
		CurVal = "INS: " + CurValShort;
	}
	else{
		CurVal = "I: " + CurValShort;
	}

	TecGUIListDeleteAllItems(SLSelVar_SLST_T3_1);

	vector<string> SearchList;
	if (INSOnly) {
		SearchList = { "INS: " };
	}
	else {
		SearchList = { "I: ", "IN: ", "INS: " };
	}
	int SelVar = 1;
	for (int vi = 1; vi <= TecUtilDataSetGetNumVars(); ++vi) {
		bool IsMatch = false;
		char * TmpCStr;
		TecUtilVarGetName(vi, &TmpCStr);
		string VName = TmpCStr;
		string VNameShort;
		TecUtilStringDealloc(&TmpCStr);
		for (auto & s : SearchList) {
			if (VName.find(s) != string::npos) {
				IsMatch = true;
				VNameShort = StringRemoveSubString(VName, s);
			}
		}
		if (IsMatch) {
			if (INSOnly) {
				TecGUIListAppendItem(SLSelVar_SLST_T3_1, VNameShort.c_str());
				if (VNameShort == CurValShort) {
					SelVar = TecGUIListGetItemCount(SLSelVar_SLST_T3_1);
				}
			}
			else {
				TecGUIListAppendItem(SLSelVar_SLST_T3_1, VName.c_str());
				if (VName == CurVal) {
					SelVar = TecGUIListGetItemCount(SLSelVar_SLST_T3_1);
				}
			}
		}
	}

	if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0) {
		TecGUIListSetSelectedItem(SLSelVar_SLST_T3_1, SelVar);
	}
}
