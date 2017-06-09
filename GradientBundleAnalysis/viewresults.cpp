
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

//#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "ZONEVARINFO.h"
#include "VIEWRESULTS.h"

using std::string;
using std::to_string;
using std::stringstream;
using std::vector;
using std::stoi;

/*
*	Result viewer sidebar functions
*/
void GBAResultViewerPrepareGUI(){

	TecUtilLockStart(AddOnID);

	/*
	*	Clear lists and stuff
	*/

	TecGUIListDeleteAllItems(SLSelSphere_SLST_S1);
	TecGUIListDeleteAllItems(MLSelGB_MLST_S1);
	TecGUIToggleSet(TGLSphereVis_TOG_S1, FALSE);
	/*
	*	First, populate the list of spheres.
	*	Get a total list, then load them alphabetically
	*	with atoms first and cages second.
	*/

	Boolean_t IsOk = TRUE;

	vector<string> SphereCPNameList;

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();

	TecUtilDataLoadBegin();

	for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
		Boolean_t ZoneOK = TRUE;
		vector<string> AuxDataCheckStr = {
			GBAZoneType,
			GBASphereCPName
		};
		AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
		ZoneOK = VALID_REF(TempAuxData);
		for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK; ++i){
			char* TempCStr;
			AuxDataType_e ADTJunk;
			Boolean_t BJunk;
			ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
			if (ZoneOK){
				if (i == 0)
					ZoneOK = (string(TempCStr) == GBAZoneTypeSphereZone);
				else
					SphereCPNameList.push_back(string(TempCStr));
			}
			TecUtilStringDealloc(&TempCStr);
		}
	}

	IsOk = (SphereCPNameList.size() > 0);
	if (IsOk){
		/*
		*	Sort list of spheres and select first one
		*/
		SortCPNameList(SphereCPNameList);
		for (vector<string>::iterator it = SphereCPNameList.begin(); it != SphereCPNameList.end(); it++){
			TecGUIListAppendItem(SLSelSphere_SLST_S1, it->c_str());
		}
		TecGUIListSetSelectedItem(SLSelSphere_SLST_S1, 1);
		GBAResultViewerSelectSphere();
	}

	TecUtilDataLoadEnd();

	TecGUISidebarActivate(Sidebar1Manager);

	TecUtilLockFinish(AddOnID);
}


void GBAResultViewerSelectSphere(){
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	Boolean_t IsOk = TRUE;

	/*
	*	Clear gradient bundle list
	*/
	TecGUIListDeleteAllItems(MLSelGB_MLST_S1);

	/*
	*	Get selected sphere name
	*/
	LgIndex_t SelSphereNum = TecGUIListGetSelectedItem(SLSelSphere_SLST_S1);
	char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_S1, SelSphereNum);

	/*
	*	Find sphere zone and
	*	1) check if it's active
	*	2) get it's gradient bundle volume CPs
	*/

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();

	EntIndex_t SphereZoneNum;

	if (IsOk){
		TecUtilDataLoadBegin();

		Boolean_t IsFound = FALSE;
		vector<string> AuxDataCheckStr = {
			GBAZoneType,
			GBASphereCPName
		};
		vector<string> CheckStr = {
			GBAZoneTypeSphereZone,
			string(SphereNameCStr)
		};
		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && !IsFound; ++ZoneNum){
			Boolean_t ZoneOK = TRUE;

			AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
			ZoneOK = VALID_REF(TempAuxData);
			for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK; ++i){
				char* TempCStr;
				AuxDataType_e ADTJunk;
				Boolean_t BJunk;
				ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
				if (ZoneOK){
					ZoneOK = (string(TempCStr) == CheckStr[i]);
				}
				TecUtilStringDealloc(&TempCStr);
			}
			if (ZoneOK){
				IsFound = TRUE;
				SphereZoneNum = ZoneNum;
			}
		}

		TecUtilDataLoadEnd();

		IsOk = IsFound;
	}

	if (IsOk)
		TecGUIToggleSet(TGLSphereVis_TOG_S1, TecUtilZoneIsActive(SphereZoneNum));


	TecGUIListDeleteAllItems(MLSelGB_MLST_S1);

	if (TecGUIToggleGet(TGLSphereVis_TOG_S1)){
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
			TecUtilDataLoadBegin();
			/*
			*	Get full name list
			*/
			vector<string> AuxDataCheckStr = {
				GBAZoneType,
				GBASphereCPName,
				GBAVolumeCPName
			};
			vector<string> CheckStr = {
				GBAZoneTypeFEVolumeZone,
				string(SphereNameCStr)
			};

			for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
				Boolean_t ZoneOK = TRUE;

				AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
				ZoneOK = VALID_REF(TempAuxData);
				for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK; ++i){
					char* TempCStr;
					AuxDataType_e ADTJunk;
					Boolean_t BJunk;
					ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
					if (ZoneOK){
						if (i <= 1){
							ZoneOK = (string(TempCStr) == CheckStr[i]);
						}
						else{
							GBFullCPNames.push_back(string(TempCStr));
							GBVolIsActive.push_back(TecUtilZoneIsActive(ZoneNum));
						}

					}
					TecUtilStringDealloc(&TempCStr);
				}
			}

			/*
			*	Get unique name list
			*/
			for (vector<string>::iterator it1 = GBFullCPNames.begin(); it1 != GBFullCPNames.end(); it1++){
				Boolean_t IsFound = FALSE;
				for (vector<string>::iterator it2 = GBUniqueCPNames.begin(); it2 != GBUniqueCPNames.end() && !IsFound; it2++){
					IsFound = (*it1 == *it2);
				}
				if (!IsFound)
					GBUniqueCPNames.push_back(*it1);
			}
			SortCPNameList(GBUniqueCPNames);

			/*
			*	Get list of active GB vols.
			*	Also add unique names to list while I'm at it.
			*/
			for (int i = 0; i < GBUniqueCPNames.size(); ++i){
				TecGUIListAppendItem(MLSelGB_MLST_S1, GBUniqueCPNames[i].c_str());
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
			if (SelectNums.size() > 0){
				TecGUIListSetSelectedItems(MLSelGB_MLST_S1, SelectNums.data(), (LgIndex_t)SelectNums.size());
			}

			TecUtilDataLoadEnd();
		}
	}

	GBAResultViewerSelectGB();

	TecUtilStringDealloc(&SphereNameCStr);

	TecUtilDrawGraphics(TRUE);
	TecUtilLockFinish(AddOnID);
}

void GBAResultViewerSelectGB(){
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	Boolean_t IsOk = TRUE;

	LgIndex_t *SelectedNums;
	LgIndex_t NumSelected;

	TecGUIListGetSelectedItems(MLSelGB_MLST_S1, &SelectedNums, &NumSelected);

	if (NumSelected > 0){

		TecUtilDataLoadBegin();

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_S1, TecGUIListGetSelectedItem(SLSelSphere_SLST_S1));

		Set_pa ActivateSet = TecUtilSetAlloc(FALSE);
		Set_pa DeactivateSet = TecUtilSetAlloc(FALSE);

		for (int i = 0; i < NumSelected && IsOk; ++i){

			char* ChkName = TecGUIListGetString(MLSelGB_MLST_S1, SelectedNums[i]);
			vector<string> AuxDataCheckStr = {
				GBAZoneType,
				GBASphereCPName,
				GBAVolumeCPName
			};
			vector<string> CheckStr = {
				GBAZoneTypeFEVolumeZone,
				string(SphereNameCStr),
				string(ChkName)
			};
			TecUtilStringDealloc(&ChkName);

			for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
				Boolean_t ZoneOK = TRUE;

				AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
				ZoneOK = VALID_REF(TempAuxData);
				for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK && IsOk; ++i){
					char* TempCStr;
					string TempStr;
					AuxDataType_e ADTJunk;
					Boolean_t BJunk;
					ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
					if (ZoneOK){
						TempStr = TempCStr;
						ZoneOK = (TempStr == CheckStr[i]);
					}
					TecUtilStringDealloc(&TempCStr);
					if (!ZoneOK && i == 2){
						Boolean_t IsFound = FALSE;
						for (int j = 0; j < NumSelected && !IsFound; ++j){
							char *TempCStr2 = TecGUIListGetString(MLSelGB_MLST_S1, SelectedNums[j]);
							string TempStr2 = TempCStr2;
							IsFound = (TempStr == TempStr2);
							TecUtilStringDealloc(&TempCStr2);
						}
						if (!IsFound)
							IsOk = TecUtilSetAddMember(DeactivateSet, ZoneNum, FALSE);
					}
				}
				if (ZoneOK){
					IsOk = TecUtilSetAddMember(ActivateSet, ZoneNum, FALSE);
				}
			}
		}

		if (IsOk){
			TecUtilZoneSetActive(ActivateSet, AssignOp_PlusEquals);
			TecUtilZoneSetActive(DeactivateSet, AssignOp_MinusEquals);
		}


		TecUtilStringDealloc(&SphereNameCStr);
		TecUtilSetDealloc(&ActivateSet);
		TecUtilSetDealloc(&DeactivateSet);

		TecUtilDataLoadEnd();
	}
	else if (TecGUIListGetItemCount(MLSelGB_MLST_S1) == 0){
		TecUtilDataLoadBegin();

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_S1, TecGUIListGetSelectedItem(SLSelSphere_SLST_S1));

		Set_pa DeactivateSet = TecUtilSetAlloc(FALSE);

		vector<string> AuxDataCheckStr = {
			GBAZoneType,
			GBASphereCPName,
		};
		vector<string> CheckStr = {
			GBAZoneTypeFEVolumeZone,
			string(SphereNameCStr),
		};

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
			Boolean_t ZoneOK = TRUE;

			AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
			ZoneOK = VALID_REF(TempAuxData);
			for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK && IsOk; ++i){
				char* TempCStr;
				AuxDataType_e ADTJunk;
				Boolean_t BJunk;
				ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
				if (ZoneOK){
					ZoneOK = (string(TempCStr) == CheckStr[i]);
				}
				TecUtilStringDealloc(&TempCStr);
			}
			if (ZoneOK){
				IsOk = TecUtilSetAddMember(DeactivateSet, ZoneNum, FALSE);
			}
		}

		if (IsOk){
			TecUtilZoneSetActive(DeactivateSet, AssignOp_MinusEquals);
		}

		TecUtilStringDealloc(&SphereNameCStr);
		TecUtilSetDealloc(&DeactivateSet);

		TecUtilDataLoadEnd();
	}

	TecUtilArrayDealloc((void**)&SelectedNums);

	TecUtilDrawGraphics(TRUE);
	TecUtilLockFinish(AddOnID);
}

void GBAResultViewerToggleSphere(){
	TecUtilLockStart(AddOnID);
	TecUtilDrawGraphics(FALSE);

	Boolean_t IsOk = (TecGUIListGetItemCount(SLSelSphere_SLST_S1) > 0);

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();

	Boolean_t ActivateZone = TecGUIToggleGet(TGLSphereVis_TOG_S1);

	if (IsOk)
	{
		/*
		*	Find and activate or deactivate sphere zone
		*/
		char* ChkName = TecGUIListGetString(SLSelSphere_SLST_S1, TecGUIListGetSelectedItem(SLSelSphere_SLST_S1));
		vector<string> AuxDataCheckStr = {
			GBAZoneType,
			GBASphereCPName
		};
		vector<string> CheckStr = {
			GBAZoneTypeSphereZone,
			string(ChkName)
		};
		TecUtilStringDealloc(&ChkName);

		TecUtilDataLoadBegin();

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
			Boolean_t ZoneOK = TRUE;

			AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
			ZoneOK = VALID_REF(TempAuxData);
			for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK && IsOk; ++i){
				char* TempCStr;
				AuxDataType_e ADTJunk;
				Boolean_t BJunk;
				ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
				if (ZoneOK){
					ZoneOK = (string(TempCStr) == CheckStr[i]);
				}
				TecUtilStringDealloc(&TempCStr);
			}
			if (ZoneOK){
				Set_pa TempSet = TecUtilSetAlloc(FALSE);
				TecUtilSetAddMember(TempSet, ZoneNum, FALSE);
				if (ActivateZone){
					TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
				}
				else{
					TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
					TecGUIListDeleteAllItems(MLSelGB_MLST_S1);
				}
				TecUtilSetDealloc(&TempSet);
				break;
			}
		}

		/*
		*	Find and activate or deactivate IB edge zone
		*/
		ChkName = TecGUIListGetString(SLSelSphere_SLST_S1, TecGUIListGetSelectedItem(SLSelSphere_SLST_S1));
		CheckStr = {
			GBAZoneTypeIBEdgeZone,
			string(ChkName)
		};
		TecUtilStringDealloc(&ChkName);

		TecUtilDataLoadBegin();

		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && IsOk; ++ZoneNum){
			Boolean_t ZoneOK = TRUE;

			AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
			ZoneOK = VALID_REF(TempAuxData);
			for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK && IsOk; ++i){
				char* TempCStr;
				AuxDataType_e ADTJunk;
				Boolean_t BJunk;
				ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
				if (ZoneOK){
					ZoneOK = (string(TempCStr) == CheckStr[i]);
				}
				TecUtilStringDealloc(&TempCStr);
			}
			if (ZoneOK){
				Set_pa TempSet = TecUtilSetAlloc(FALSE);
				TecUtilSetAddMember(TempSet, ZoneNum, FALSE);
				if (ActivateZone){
					TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
				}
				else{
					TecUtilZoneSetActive(TempSet, AssignOp_MinusEquals);
					TecGUIListDeleteAllItems(MLSelGB_MLST_S1);
				}
				TecUtilSetDealloc(&TempSet);
				break;
			}
		}
	}
	else{
		TecGUIToggleSet(TGLSphereVis_TOG_S1, FALSE);
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
	if (IsOk){
		AuxData_pa CPAuxData = TecUtilAuxDataZoneGetRef(CPZoneNum);
		if (VALID_REF(CPAuxData)){
			for (int i = 0; i < 4; ++i){
				char* TempCStr;
				AuxDataType_e ADTJunk;
				Boolean_t BJunk;
				IsOk = TecUtilAuxDataGetItemByName(CPAuxData, CCDataNumCPs[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
				if (IsOk){
					NumCPs[i] = stoi(TempCStr);
				}
				TecUtilStringDealloc(&TempCStr);
			}
		}
	}

	/*
	*	Now sort
	*/
	if (IsOk){
		string TmpStr;
		int TmpInt;
		for (int i = 0; i < StrList.size(); ++i){
			bool DidSwap = false;
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
					DidSwap = true;
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

	if (TecGUIListGetItemCount(SLSelSphere_SLST_S1) > 0){
		TecUtilDrawGraphics(FALSE);
		/*
		*	Get selected sphere name
		*/
		LgIndex_t SelSphereNum = TecGUIListGetSelectedItem(SLSelSphere_SLST_S1);
		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_S1, SelSphereNum);

		/*
		*	Find sphere zone and
		*	1) check if it's active
		*	2) get it's gradient bundle volume CPs
		*/

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		Set_pa DeleteSet = TecUtilSetAlloc(FALSE);

		TecUtilDataLoadBegin();

		Boolean_t IsFound = FALSE;
		vector<string> AuxDataCheckStr = {
			GBASphereCPName
		};
		vector<string> CheckStr = {
			string(SphereNameCStr)
		};
		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && !IsFound; ++ZoneNum){
			Boolean_t ZoneOK = TRUE;

			AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
			ZoneOK = VALID_REF(TempAuxData);
			for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK; ++i){
				char* TempCStr;
				AuxDataType_e ADTJunk;
				Boolean_t BJunk;
				ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
				if (ZoneOK){
					ZoneOK = (string(TempCStr) == CheckStr[i]);
				}
				TecUtilStringDealloc(&TempCStr);
			}
			if (ZoneOK){
				TecUtilSetAddMember(DeleteSet, ZoneNum, FALSE);
			}
		}

		if (!TecUtilDataSetDeleteZone(DeleteSet)){
			TecUtilDialogErrMsg("Failed to delete zones.");
		}

		TecUtilSetDealloc(&DeleteSet);


		TecGUIListDeleteItemAtPos(SLSelSphere_SLST_S1, SelSphereNum);
		if (TecGUIListGetItemCount(SLSelSphere_SLST_S1) > 0){
			TecGUIListSetSelectedItem(SLSelSphere_SLST_S1, 1);
			GBAResultViewerSelectSphere();
		}
		else{
			TecGUIListDeleteAllItems(MLSelGB_MLST_S1);
			TecGUIToggleSet(TGLSphereVis_TOG_S1, FALSE);
		}

		TecUtilStringDealloc(&SphereNameCStr);

		TecUtilDataLoadEnd();
		TecUtilDrawGraphics(TRUE);
	}
	TecUtilLockFinish(AddOnID);

}

void GBAResultViewerActivateAllGB(){
	TecUtilLockStart(AddOnID);

	if (TecGUIListGetItemCount(SLSelSphere_SLST_S1) > 0){
		TecUtilDrawGraphics(FALSE);
		/*
		*	Get selected sphere name
		*/
		LgIndex_t SelSphereNum = TecGUIListGetSelectedItem(SLSelSphere_SLST_S1);
		char* SphereNameCStr = TecGUIListGetString(SLSelSphere_SLST_S1, SelSphereNum);

		/*
		*	Find sphere zone and
		*	1) check if it's active
		*	2) get it's gradient bundle volume CPs
		*/

		EntIndex_t NumZones = TecUtilDataSetGetNumZones();

		Set_pa ActivateSet = TecUtilSetAlloc(FALSE);

		TecUtilDataLoadBegin();

		Boolean_t IsFound = FALSE;
		vector<string> AuxDataCheckStr = {
			GBAZoneType,
			GBASphereCPName
		};
		vector<string> CheckStr = {
			GBAZoneTypeFEVolumeZone,
			string(SphereNameCStr)
		};
		for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones && !IsFound; ++ZoneNum){
			Boolean_t ZoneOK = TRUE;

			AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
			ZoneOK = VALID_REF(TempAuxData);
			for (int i = 0; i < AuxDataCheckStr.size() && ZoneOK; ++i){
				char* TempCStr;
				AuxDataType_e ADTJunk;
				Boolean_t BJunk;
				ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataCheckStr[i].c_str(), (ArbParam_t*)&TempCStr, &ADTJunk, &BJunk);
				if (ZoneOK){
					ZoneOK = (string(TempCStr) == CheckStr[i]);
				}
				TecUtilStringDealloc(&TempCStr);
			}
			if (ZoneOK){
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
			int NumGBAtoms = TecGUIListGetItemCount(MLSelGB_MLST_S1);
			vector<int> SelNums(NumGBAtoms);
			for (int i = 1; i <= NumGBAtoms; ++i)
				SelNums[i - 1] = i;

			TecGUIListSetSelectedItems(MLSelGB_MLST_S1, SelNums.data(), NumGBAtoms);
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

		Boolean_t IsOk = TRUE;

		EntIndex_t ProbedZoneNum = TecUtilProbeFieldGetZone();

		AuxData_pa ProbedZoneAuxDataPtr = TecUtilAuxDataZoneGetRef(ProbedZoneNum);
		IsOk = VALID_REF(ProbedZoneAuxDataPtr);

		string CPName;

		TecUtilDataLoadBegin();

		if (IsOk){
			char* TempCStr;
			AuxDataType_e AuxDataJunk;
			Boolean_t BooleanJunk;
			stringstream ss;
			ss << GBADataPrefix << GBASphereCPName;
			if (TecUtilAuxDataGetItemByName(ProbedZoneAuxDataPtr, ss.str().c_str(), (ArbParam_t*)&TempCStr, &AuxDataJunk, &BooleanJunk))
				CPName = TempCStr;
			else{
				IsOk = FALSE;
			}
			TecUtilStringDealloc(&TempCStr);
		}

		string ZoneType;
		if (IsOk){
			char* TempCStr;
			AuxDataType_e AuxDataJunk;
			Boolean_t BooleanJunk;
			stringstream ss;
			ss << GBADataPrefix << GBAZoneType;
			if (TecUtilAuxDataGetItemByName(ProbedZoneAuxDataPtr, ss.str().c_str(), (ArbParam_t*)&TempCStr, &AuxDataJunk, &BooleanJunk))
				ZoneType = TempCStr;
			else{
				IsOk = FALSE;
			}
			TecUtilStringDealloc(&TempCStr);
		}

		if (IsOk){
			if (ZoneType == GBAZoneTypeSphereZone){
				EntIndex_t NumZones = TecUtilDataSetGetNumZones();
				if (isNearestPoint){
					/*
					*	User selected a node, so activate all FE volumes
					*	around that node.
					*/
					LgIndex_t NodeNum = TecUtilProbeGetPointIndex();
					int NumFound = 0;
					for (int CurZoneNum = 1; CurZoneNum < NumZones && NumFound < 6; ++CurZoneNum){
						if (TecUtilZoneGetType(CurZoneNum) == ZoneType_FEQuad){
							AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(CurZoneNum);
							Boolean_t ZoneOK = TRUE;
							if (VALID_REF(TempAuxData)){
								if (ZoneOK){
									vector<string> CheckAuxNameStrings = {
										GBAZoneType,
										GBASphereCPName,
									};
									vector<string> CheckAuxValues = {
										GBAZoneTypeFEVolumeZone,
										CPName,
									};
									for (int i = 0; i < CheckAuxNameStrings.size() && ZoneOK; ++i){
										char* TempCStr;
										AuxDataType_e AuxDataJunk;
										Boolean_t BooleanJunk;
										ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, CheckAuxNameStrings[i].c_str(), (ArbParam_t*)&TempCStr, &AuxDataJunk, &BooleanJunk);
										if (ZoneOK){
											ZoneOK = (string(TempCStr) == CheckAuxValues[i]);
										}
										TecUtilStringDealloc(&TempCStr);
									}
								}
								if (ZoneOK){
									vector<string> CheckAuxNameStrings = {
										GBANodeNums[0],
										GBANodeNums[1],
										GBANodeNums[2]
									};
									Boolean_t NodeFound = FALSE;
									for (int i = 0; i < CheckAuxNameStrings.size() && ZoneOK && !NodeFound; ++i){
										char* TempCStr;
										AuxDataType_e AuxDataJunk;
										Boolean_t BooleanJunk;
										ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, CheckAuxNameStrings[i].c_str(), (ArbParam_t*)&TempCStr, &AuxDataJunk, &BooleanJunk);
										if (ZoneOK){
											NodeFound = (string(TempCStr) == to_string(NodeNum));
										}
										TecUtilStringDealloc(&TempCStr);
									}
									ZoneOK = NodeFound;
								}
								if (ZoneOK){
									Set_pa TempSet = TecUtilSetAlloc(FALSE);
									IsOk = TecUtilSetAddMember(TempSet, CurZoneNum, FALSE);
									TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
									TecUtilSetDealloc(&TempSet);
									NumFound++;
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
					Boolean_t FoundVolume = FALSE;
					for (int CurZoneNum = 1; CurZoneNum < NumZones && !FoundVolume; ++CurZoneNum){
						if (TecUtilZoneGetType(CurZoneNum) == ZoneType_FEQuad){
							AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(CurZoneNum);
							Boolean_t ZoneOK = TRUE;
							if (VALID_REF(TempAuxData)){
								vector<string> CheckAuxNameStrings = {
									GBAZoneType,
									GBASphereCPName,
									GBAElemNum
								};
								vector<string> CheckAuxValues = {
									GBAZoneTypeFEVolumeZone,
									CPName,
									to_string(ElemNum)
								};
								for (int i = 0; i < CheckAuxNameStrings.size() && ZoneOK; ++i){
									char* TempCStr;
									AuxDataType_e AuxDataJunk;
									Boolean_t BooleanJunk;
									ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, CheckAuxNameStrings[i].c_str(), (ArbParam_t*)&TempCStr, &AuxDataJunk, &BooleanJunk);
									if (ZoneOK){
										ZoneOK = (string(TempCStr) == CheckAuxValues[i]);
									}
									TecUtilStringDealloc(&TempCStr);
								}
								FoundVolume = ZoneOK;
							}
						}
						if (FoundVolume){
							Set_pa TempSet = TecUtilSetAlloc(FALSE);
							IsOk = TecUtilSetAddMember(TempSet, CurZoneNum, FALSE);
							TecUtilZoneSetActive(TempSet, AssignOp_PlusEquals);
							TecUtilSetDealloc(&TempSet);
						}
					}
				}
			}
			else if (ZoneType == GBAZoneTypeFEVolumeZone){
				/*
				*	Two modes:
				*	1. If user does nearest point probe, then the probed FE
				*	volume and all volumes that share nodes with it are
				*	deactivated.
				*	2. If user does normal probe, only probed FE zone is deactivated.
				*/

				EntIndex_t FEZoneNum = TecUtilProbeFieldGetZone();

				if (isNearestPoint){

					AuxData_pa FEAuxData = TecUtilAuxDataZoneGetRef(FEZoneNum);
					IsOk = VALID_REF(FEAuxData);
					/*
					*	Get node nums for this fe volume zone
					*/
					vector<string> CheckAuxNameStrings = {
						GBANodeNums[0],
						GBANodeNums[1],
						GBANodeNums[2]
					};
					int NodeNums[3];
					for (int i = 0; i < CheckAuxNameStrings.size() && IsOk; ++i){
						char* TempCStr;
						AuxDataType_e AuxDataJunk;
						Boolean_t BooleanJunk;
						IsOk = TecUtilAuxDataGetItemByName(FEAuxData, CheckAuxNameStrings[i].c_str(), (ArbParam_t*)&TempCStr, &AuxDataJunk, &BooleanJunk);
						if (IsOk){
							NodeNums[i] = stoi(string(TempCStr));
						}
						TecUtilStringDealloc(&TempCStr);
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
							if (CurZoneNum != FEZoneNum && TecUtilZoneIsActive((EntIndex_t)CurZoneNum) && TecUtilZoneGetType((EntIndex_t)CurZoneNum) == ZoneType_FEQuad){
								AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef((EntIndex_t)CurZoneNum);
								if (VALID_REF(TempAuxData)){
									Boolean_t ZoneOK = TRUE;
									if (ZoneOK){
										vector<string> CheckAuxNameStrings2 = {
											GBAZoneType,
											GBASphereCPName,
										};
										vector<string> CheckAuxValues = {
											GBAZoneTypeFEVolumeZone,
											CPName,
										};
										for (int i = 0; i < CheckAuxNameStrings2.size() && ZoneOK; ++i){
											char* TempCStr;
											AuxDataType_e AuxDataJunk;
											Boolean_t BooleanJunk;
											ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, CheckAuxNameStrings2[i].c_str(), (ArbParam_t*)&TempCStr, &AuxDataJunk, &BooleanJunk);
											if (ZoneOK){
												ZoneOK = (string(TempCStr) == CheckAuxValues[i]);
											}
											TecUtilStringDealloc(&TempCStr);
										}
									}
									if (ZoneOK){
										for (int i = 0; i < CheckAuxNameStrings.size() && IsOk; ++i){
											char* TempCStr;
											AuxDataType_e AuxDataJunk;
											Boolean_t BooleanJunk;
											IsOk = TecUtilAuxDataGetItemByName(TempAuxData, CheckAuxNameStrings[i].c_str(), (ArbParam_t*)&TempCStr, &AuxDataJunk, &BooleanJunk);
											if (IsOk){
												for (int j = 0; j < 3; ++j){
													if (NodeNums[j] == stoi(string(TempCStr))){
														NodeNeighborZones[j].push_back((int)CurZoneNum);
														break;
													}
												}
											}
											TecUtilStringDealloc(&TempCStr);
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
	TecUtilLockFinish(AddOnID);
}