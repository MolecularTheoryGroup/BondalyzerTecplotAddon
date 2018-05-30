
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include <armadillo>

#include "CSM_DATA_SET_INFO.h"
#include "VIEWRESULTS.h"

using namespace arma;

using std::string;
using std::to_string;
using std::stringstream;
using std::vector;
using std::stoi;





void GBAResultViewerSelectSphere(){
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);


	Boolean_t IsOk = TRUE;

	/*
	*	Clear gradient bundle list
	*/
	TecGUIListDeleteAllItems(MLSelGB_MLST_T3_1);

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

	string TmpStr1, TmpStr2, TmpStr3, TmpStr4;

	TmpStr1 = CSMAuxData.GBA.ZoneType;
	TmpStr2 = CSMAuxData.GBA.SphereCPName;

	if (IsOk){
// 		TecUtilDataLoadBegin();

		Boolean_t IsFound = FALSE;
		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && !IsFound; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, TmpStr1, CSMAuxData.GBA.ZoneTypeSphereZone)
				&& AuxDataZoneItemMatches(ZoneNum, TmpStr2, SphereNameStr))
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


	TecGUIListDeleteAllItems(MLSelGB_MLST_T3_1);
	int SelectedVarNum = 1;
	if (TecGUIListGetItemCount(SLSelVar_SLST_T3_1) > 0)
		SelectedVarNum = TecGUIListGetSelectedItem(SLSelVar_SLST_T3_1);

	TecGUIListDeleteAllItems(SLSelVar_SLST_T3_1);

	if (TecGUIToggleGet(TGLSphereVis_TOG_T3_1)){
		/*
		 *	Get list of integrated variables for the selected sphere
		 */

		vector<string> IntVarNames;
		vector<string> IntCheckStrs = { "I: ", "IN: ", "INS: ", " Integration" };
		for (int i = 1; i <= TecUtilDataSetGetNumVars(); ++i){
			char *VarName, *CheckStr;
			if (TecUtilVarGetName(i, &VarName)){
				for (const string & Str : IntCheckStrs){
					CheckStr = std::strstr(VarName, Str.c_str());
					if (CheckStr != NULL){
						// Integration variable found. Now make sure it's not bit type for sphere zone.
						if (TecUtilDataValueGetType(SphereZoneNum, i) != FieldDataType_Bit){
							IntVarNames.push_back(VarName);
							break;
						}
					}
				}
				TecUtilStringDealloc(&VarName);
			}
		}

		if (IntVarNames.size() > 0){
			for (const string & i : IntVarNames)
				TecGUIListAppendItem(SLSelVar_SLST_T3_1, i.c_str());
			TecGUIListSetSelectedItem(SLSelVar_SLST_T3_1, SelectedVarNum);
		}

		/*
		*	Get list of Gradient bundles for sphere zone
		*	(one's that correspond to a volume CP)
		*	and see if they're active.
		*	"Active" here means that ALL gradient bundles around
		*	a volume CP node are active.
		*/

		vector<string> GBFullCPNames;
		vector<Boolean_t> GBVolIsActive;
		vector<string> GBUniqueCPNames;
		vector<LgIndex_t> SelectNums;

		if (IsOk){
// 			TecUtilDataLoadBegin();
			/*
			*	Get full name list
			*/
			
			TmpStr3 = CSMAuxData.GBA.VolumeCPName;

			for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
				if (AuxDataZoneItemMatches(ZoneNum, TmpStr1, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
					&& AuxDataZoneItemMatches(ZoneNum, TmpStr2, SphereNameStr))
				{
					if (AuxDataZoneGetItem(ZoneNum, TmpStr3, TmpStr4)){
						GBFullCPNames.push_back(TmpStr4);
						GBVolIsActive.push_back(TecUtilZoneIsActive(ZoneNum));
					}
				}
			}

			/*
			*	Get unique name list
			*/
			for (const string & it1 : GBFullCPNames){
				if (VectorGetElementNum(GBUniqueCPNames, it1) < 0)
					GBUniqueCPNames.push_back(it1);
			}
			SortCPNameList(GBUniqueCPNames);

			/*
			*	Get list of active GB vols.
			*	Also add unique names to list while I'm at it.
			*/
			for (int i = 0; i < GBUniqueCPNames.size(); ++i){
				TecGUIListAppendItem(MLSelGB_MLST_T3_1, GBUniqueCPNames[i].c_str());
				int HitCount = 0;
				for (int j = 0; j < GBFullCPNames.size(); ++j){
					if (GBVolIsActive[j] && GBUniqueCPNames[i] == GBFullCPNames[j]){
						HitCount++;
					}
				}
				if (HitCount >= 5){
					SelectNums.push_back(i + 1);
				}
			}

			/*
			*	Select GB's in list that are active
			*/
// 			if (SelectNums.size() > 0){
// 				TecGUIListSetSelectedItems(MLSelGB_MLST_T3_1, SelectNums.data(), (LgIndex_t)SelectNums.size());
// 			}

// 			TecUtilDataLoadEnd();
		}
	}

// 	GBAResultViewerSelectGB();


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
				&& (!UseSelectedSphere || AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr))){
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
				for (const int & ZoneNum : SphereZoneNums){
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
						for (const double & Val : TmpVals){
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
				int NegExp = log10(-NegCloseToZero) - 1;
				int PosExp = log10(PosCloseToZero) - 1;

				if (VarMax < 0){
					ContourLevelsVec = logspace(log10(-VarMax)-1, log10(-VarMin)+1, NumContours) * -1.0;
				}
				else if (VarMax == 0.0){
					ContourLevelsVec = logspace(NegExp, log10(-VarMin) + 1, NumContours) * -1.0;
				}
				else if (VarMin == 0.0){
					ContourLevelsVec = logspace(PosExp, log10(VarMax) + 1, NumContours);
				}
				else{
					int NegExpUpper = log10(-VarMin) + 1;
					int PosExpUpper = log10(VarMax) + 1;
					int NumNegContours = static_cast<int>(static_cast<double>(NumContours)* (static_cast<double>(NegExpUpper - NegExp) / (static_cast<double>(NegExpUpper - NegExp) + static_cast<double>(PosExpUpper - PosExp))));
// 					int NumNegContours = static_cast<int>(static_cast<double>(NumContours)* (-VarMin / (VarMax - VarMin)));
					vec NegContours = flipud(logspace(NegExp, NegExpUpper, NumNegContours)) * -1.0;
					ContourLevelsVec = join_cols(NegContours, logspace(PosExp, PosExpUpper, NumContours - NumNegContours));
				}
			}
			else{
				ContourLevelsVec = logspace(log10(VarMin)-1, log10(VarMax)+1, NumContours);
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
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr)){
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
				if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
					&& AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr))
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
			if (AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.ZoneType, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
				&& AuxDataZoneItemMatches(ZoneNum, CSMAuxData.GBA.SphereCPName, SphereNameStr))
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

	string TmpStr1, TmpStr2;


	TmpStr1 = CSMAuxData.GBA.ZoneType;
	TmpStr2 = CSMAuxData.GBA.SphereCPName;

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
			if (AuxDataZoneItemMatches(ZoneNum, TmpStr1, CSMAuxData.GBA.ZoneTypeSphereZone)
				&& AuxDataZoneItemMatches(ZoneNum, TmpStr2, ChkNameStr))
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
			if (AuxDataZoneItemMatches(ZoneNum, TmpStr1, CSMAuxData.GBA.ZoneTypeSphereZone)
				&& AuxDataZoneItemMatches(ZoneNum, TmpStr2, ChkNameStr))
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

	string TmpStr1;

	TmpStr1 = CSMAuxData.GBA.SphereCPName;

	if (TecGUIListGetItemCount(SLSelSphere_SLST_T3_1) > 0){
		TecUtilDrawGraphics(FALSE);
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

		Set_pa DeleteSet = TecUtilSetAlloc(FALSE);

		TecUtilDataLoadBegin();
		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, TmpStr1, SphereNameStr)){
				TecUtilSetAddMember(DeleteSet, ZoneNum, FALSE);
			}
		}

		if (!TecUtilDataSetDeleteZone(DeleteSet)){
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


		TecUtilDataLoadEnd();
		TecUtilDrawGraphics(TRUE);
	}
	TecUtilLockFinish(AddOnID);

}

void GBAResultViewerActivateAllGB(){
	TecUtilLockStart(AddOnID);

	string TmpStr1, TmpStr2;

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

		TmpStr1 = CSMAuxData.GBA.ZoneType;
		TmpStr2 = CSMAuxData.GBA.SphereCPName;

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
			if (AuxDataZoneItemMatches(ZoneNum, TmpStr1, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
				&& AuxDataZoneItemMatches(ZoneNum, TmpStr2, SphereNameStr))
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

		TmpStr1 = CSMAuxData.GBA.ZoneType;
		TmpStr2 = CSMAuxData.GBA.SphereCPName;
		TmpStr3 = CSMAuxData.GBA.ElemNum;

		if (IsOk){
			IsOk = AuxDataZoneGetItem(ProbedZoneNum, TmpStr2, CPName);
		}

		string ZoneType;
		if (IsOk){
			IsOk = AuxDataZoneGetItem(ProbedZoneNum, TmpStr1, ZoneType);
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
					TmpStr = to_string(NodeNum);
					int NumFound = 0;
					for (int CurZoneNum = 1; CurZoneNum < NumZones && NumFound < 6; ++CurZoneNum){
						if (TecUtilZoneGetType(CurZoneNum) == ZoneType_FETriangle){
							if (AuxDataZoneItemMatches(CurZoneNum, TmpStr1, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
								&& AuxDataZoneItemMatches(CurZoneNum, TmpStr2, CPName))
							{
								for (int i = 0; i < 3; ++i){
									if (AuxDataZoneItemMatches(CurZoneNum, CSMAuxData.GBA.NodeNums[i], TmpStr)){
										Set_pa TempSet = TecUtilSetAlloc(FALSE);
										IsOk = TecUtilSetAddMember(TempSet, CurZoneNum, FALSE);
										TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
										TecUtilSetDealloc(&TempSet);
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
					TmpStr = to_string(ElemNum);
					for (int CurZoneNum = 1; CurZoneNum < NumZones; ++CurZoneNum){
						if (TecUtilZoneGetType(CurZoneNum) == ZoneType_FETriangle){
							if (AuxDataZoneItemMatches(CurZoneNum, TmpStr1, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
								&& AuxDataZoneItemMatches(CurZoneNum, TmpStr2, CPName)
								&& AuxDataZoneItemMatches(CurZoneNum, TmpStr3, TmpStr))
							{
								Set_pa TempSet = TecUtilSetAlloc(FALSE);
								IsOk = TecUtilSetAddMember(TempSet, CurZoneNum, FALSE);
								TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
								TecUtilSetDealloc(&TempSet);
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
								if (AuxDataZoneItemMatches(static_cast<int>(CurZoneNum), TmpStr1, CSMAuxData.GBA.ZoneTypeFEVolumeZone)
									&& AuxDataZoneItemMatches(static_cast<int>(CurZoneNum), TmpStr2, CPName))
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