
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

#include <omp.h>

#include <string>
#include <sstream>
#include <vector>
#include <fstream>

// for profiling
#include <chrono>
#include <ctime>
#include <ratio>

#include "CSM_DATA_SET_INFO.h"
#include "ZONEVARINFO.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_FE_VOLUME.h"
#include "VIEWRESULTS.h"
#include "CSM_GBAGUI.h"
#include "CSM_GEOMETRY.h"

#include "INTEGRATE.h"

#include <armadillo>
using namespace arma;

using std::string;
using std::to_string;
using std::stringstream;
using std::vector;
using std::ofstream;
using std::ios;
using std::endl;

//for profiling
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;




/*
*	Result viewer sidebar functions
*/
// void GBAIntegrationPrepareGUI(){
// 
// 	TecUtilLockStart(AddOnID);
// 
// 	/*
// 	*	Clear lists and stuff
// 	*/
// 
// 	TecGUIListDeleteAllItems(MLIntSelSph_MLST_T2_1);
// 	TecGUIListDeleteAllItems(MLIntSelVar_MLST_T2_1);
// 	TecGUIToggleSet(TGLIntVolInt_TOG_T2_1, TRUE);
// 	/*
// 	*	First, populate the list of spheres.
// 	*	Get a total list, then load them alphabetically
// 	*	with atoms first and cages second.
// 	*/
// 
// 	Boolean_t IsOk = TRUE;
// 
// 	vector<string> SphereCPNameList;
// 
// 	EntIndex_t NumZones = TecUtilDataSetGetNumZones();
// 
// 	TecUtilDataLoadBegin();
// 
// 	string TmpStr1, TmpStr2;
// 
// 	TmpStr1 = CSMAuxData.GBA.ZoneType;
// 	TmpStr2 = CSMAuxData.GBA.SphereCPName;
// 
// 	for (EntIndex_t ZoneNum = 1; ZoneNum <= NumZones; ++ZoneNum){
// 		if (AuxDataZoneItemMatches(ZoneNum, TmpStr1, CSMAuxData.GBA.ZoneTypeSphereZone)){
// 			SphereCPNameList.push_back(AuxDataZoneGetItem(ZoneNum, TmpStr2));
// 		}
// 	}
// 
// 	IsOk = (SphereCPNameList.size() > 0);
// 	if (IsOk){
// 		/*
// 		*	Sort list of spheres and select first one
// 		*/
// 		SortCPNameList(SphereCPNameList);
// 		for (string const & it : SphereCPNameList){
// 			TecGUIListAppendItem(MLIntSelSph_MLST_T2_1, it.c_str());
// 		}
// 		TecGUIListSetSelectedItem(MLIntSelSph_MLST_T2_1, 1);
// 	}
// 
// 	/*
// 	 *	Populate variable list
// 	 */
// 	ListPopulateWithVarNames(MLIntSelVar_MLST_T2_1);
// 
// 	TecGUIScaleShowNumericDisplay(SCIntPrecise_SC_T2_1, TRUE);
// 	TecGUIScaleSetLimits(SCIntPrecise_SC_T2_1, 0, 4, 0);
// 	TecGUIScaleSetValue(SCIntPrecise_SC_T2_1, IntPrecise);
// 	TecGUILabelSetText(LBLIntPrecis_LBL_T2_1, IntPrecisionLabels[IntPrecise].c_str());
// 
// 	TecUtilDataLoadEnd();
// 
// 	TecUtilLockFinish(AddOnID);
// }


Boolean_t PerformIntegration(vector<string> const & AtomNameList,
	vector<string> const & IntVarNameList,
	vector<int> const & IntVarNumList,
	Boolean_t IntegrateVolume,
	int IntResolution,
	Boolean_t ActiveGBsOnly)
{

	TecUtilLockStart(AddOnID);

	Boolean_t IsOk = TRUE;


	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	Boolean_t UserQuit = FALSE;

	char* DesktopPath = nullptr;
	DesktopPath = getenv("USERPROFILE");

	string OutFileName = string(DesktopPath) + string("\\Desktop\\Out.csv");

	ofstream OutFile;
	bool PrintOutput = false;

	if (PrintOutput){
		ofstream OutFile(OutFileName.c_str(), ios::out | ios::app);
		if (!OutFile.is_open()){
			IsOk = TecUtilDialogMessageBox("Filed to open output file. (are you accessing it with another program maybe?). Continue?", MessageBox_YesNo);
		}
		else
			OutFile.close();
	}

	if (!IsOk){
		TecUtilLockFinish(AddOnID);
		return IsOk;
	}

	Boolean_t DoTiming = TRUE;
	high_resolution_clock::time_point TimeStartLoop, TimeFinish, TimeStart;


	EntIndex_t VolZoneNum = ZoneNumByName("Full Volume");

	vector<EntIndex_t> XYZVarNums(3);
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	for (int i = 0; i < 3 && IsOk; ++i)
		IsOk = (XYZVarNums[i] > 0);

	vector<EntIndex_t> XYZRhoVarNums = XYZVarNums;
	XYZRhoVarNums.push_back(VarNumByName("Electron Density"));

	int NumZones = TecUtilDataSetGetNumZones();

	vector<string> NewVarNames;
	vector<int> NewVarNums;
	vector<FieldDataType_e> ZoneDataTypes(NumZones, FieldDataType_Bit);
	for (int i = 1; i <= NumZones; ++i){
		if (AuxDataZoneHasItem(i, CSMAuxData.GBA.ZoneType)){
			ZoneDataTypes[i - 1] = FieldDataType_Float;
		}
	}
	vector<ValueLocation_e> ZoneDataLocs(NumZones, ValueLocation_CellCentered);
// 	vector<ValueLocation_e> ZoneDataLocs(NumZones, ValueLocation_Nodal);
	vector<string> VarNameAppends = { "I", "N", "S" };
	for (int i = 0; i < IntVarNameList.size() && IsOk; ++i){
		for (int j = 0; j < 3; ++j){
			string TmpStr;
			for (int k = 0; k <= j; ++k)
				TmpStr += VarNameAppends[k];
			NewVarNames.push_back(TmpStr + ": " + IntVarNameList[i]);
			NewVarNums.push_back(VarNumByName(NewVarNames.back()));
			if (NewVarNums[i * 3 + j] < 0){
				ArgList_pa ArgList = TecUtilArgListAlloc();

				IsOk = TecUtilArgListAppendString(ArgList, SV_NAME, NewVarNames.back().c_str());
				if (IsOk)
					IsOk = TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, ZoneDataTypes.data());
				if (IsOk)
					IsOk = TecUtilArgListAppendArray(ArgList, SV_VALUELOCATION, ZoneDataLocs.data());

				if (IsOk)
					IsOk = TecUtilDataSetAddVarX(ArgList);

				if (IsOk)
					NewVarNums.back() = TecUtilDataSetGetNumVars();

				if (IsOk){
					Set_pa TmpSet = TecUtilSetAlloc(FALSE);
					IsOk = TecUtilSetAddMember(TmpSet, NewVarNums.back(), FALSE);
					if (IsOk)
						TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)TmpSet);
					TecUtilSetDealloc(&TmpSet);
				}

				TecUtilArgListDealloc(&ArgList);
			}
		}
	}
	if (IntegrateVolume && IsOk){
		for (int j = 0; j < 3; ++j){
			string TmpStr;
			for (int k = 0; k <= j; ++k)
				TmpStr += VarNameAppends[k];
			NewVarNames.push_back(TmpStr + ": Volume");
			NewVarNums.push_back(VarNumByName(NewVarNames.back()));
			if (NewVarNums.back() < 0){
				ArgList_pa ArgList = TecUtilArgListAlloc();

				IsOk = TecUtilArgListAppendString(ArgList, SV_NAME, NewVarNames.back().c_str());
				if (IsOk)
					IsOk = TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, ZoneDataTypes.data());
				if (IsOk)
					IsOk = TecUtilArgListAppendArray(ArgList, SV_VALUELOCATION, ZoneDataLocs.data());

				if (IsOk)
					IsOk = TecUtilDataSetAddVarX(ArgList);

				if (IsOk)
					NewVarNums.back() = TecUtilDataSetGetNumVars();

				if (IsOk){
					Set_pa TmpSet = TecUtilSetAlloc(FALSE);
					IsOk = TecUtilSetAddMember(TmpSet, NewVarNums.back(), FALSE);
					if (IsOk)
						TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)TmpSet);
					TecUtilSetDealloc(&TmpSet);
				}

				TecUtilArgListDealloc(&ArgList);
			}
		}
	}
	vector<string> AuxDataNames(IntVarNameList.size());
	for (int i = 0; i < IntVarNameList.size(); ++i)
		AuxDataNames[i] = RemoveStringChar(IntVarNameList[i], " ");
	if (IntegrateVolume)
		AuxDataNames.push_back(string("Volume"));

	if (!IsOk){
		TecUtilDialogErrMsg("Failed to identify or create variables for integration.");
	}

	vector<int> SphereZoneNums;
	double ContourMinMax[2] = { 1e50, -1e50 };

	EntIndex_t CPZoneNum = ZoneNumByName(string("Critical Points"));

// 	vector<vec> stuW = GetWeightsPoints(IntResolution);
// 
// 	for (int i = 0; i < stuW[3].n_cols; ++i){
// 		if (stuW[3][i] < 0){
// 			TecUtilDialogMessageBox(string("negative initial weight at i = " + to_string(i)).c_str(), MessageBoxType_Error);
// 		}
// 	}



// 	vector<vector<double> > vstuW(4, vector<double>(stuW[0].size()));
// 	for (int i = 0; i < 4; ++i){
// 		for (int j = 0; j < stuW[0].size(); ++j){
// 			vstuW[i][j] = stuW[i](j);
// 		}
// 	}


	Set_pa ZoneSet = TecUtilSetAlloc(FALSE);

	for (int AtomNum = 0; AtomNum < AtomNameList.size(); ++AtomNum){


		if (DoTiming){
			TimeStartLoop = high_resolution_clock::now();
			if (AtomNum == 0) TimeStart = TimeStartLoop;
		}

		stringstream ss;
		ss << "Integrating " << AtomNum + 1;
		if (AtomNameList.size() > 1)
			ss << " of " << AtomNameList.size();
		StatusLaunch(ss.str(), AddOnID, TRUE);

		TecUtilDataLoadBegin();

		vector<FESurface_c> VolumeList;
		vector<vector<int> > GPClosestPtNums;
		
		VolumeList.reserve(NumZones);
		for (int i = 1; i <= NumZones; ++i){
			if (TecUtilZoneIsFiniteElement(i)
				&& (!ActiveGBsOnly || TecUtilZoneIsActive(i))
				&& AuxDataZoneItemMatches(i, CSMAuxData.GBA.SphereCPName, AtomNameList[AtomNum]))
			{
				VolumeList.push_back(FESurface_c(i, VolZoneNum, XYZVarNums, IntVarNumList));
			}
		}

		SphereZoneNums.push_back(VolumeList[0].GetZoneNum());

// 		AuxDataZoneDeleteItemByName(SphereZoneNums[0], CSMAuxData.GBA.IntPrecision);
// 		AuxDataZoneDeleteItemByName(SphereZoneNums[0], CSMAuxData.GBA.IntWallTime);

		int CPNum = stoi(AuxDataZoneGetItem(VolumeList[0].GetZoneNum(), CSMAuxData.GBA.SphereCPNum));

		Boolean_t FreshIntegration = FALSE;
		for (int i = 0; i < AuxDataNames.size() && !FreshIntegration; ++i)
			FreshIntegration = !AuxDataZoneHasItem(VolumeList[0].GetZoneNum(), AuxDataNames[i]);

		if (FreshIntegration || TecUtilDialogMessageBox("Variables have already been integrated for this zone."
			" Integrate again?", MessageBox_YesNo)){
			int NumVolumes = static_cast<int>(VolumeList.size());

			int NumComplete = 0;

			high_resolution_clock::time_point StatusStartTime = high_resolution_clock::now();

// #ifndef _DEBUG
#pragma omp parallel for schedule(dynamic)
// #endif
			for (int i = 0; i < NumVolumes; ++i){
// 			for (int i = 0; i < NumVolumes; i += NumVolumes / 10){
// #else
// // 			vector<int> ind = { 105, 344 };
// // 			for (int & i : ind){
// // 				i -= (57 + 1);
// 			for (int i = 0; i < NumVolumes; i += NumVolumes / 10){
// #endif
// #pragma omp critical (UpdateIntegrationStatus)
				if (omp_get_thread_num() == 0)
				{
					UserQuit = !StatusUpdate(NumComplete, NumVolumes, ss.str(), AddOnID, StatusStartTime);
					#pragma omp flush (UserQuit)
				}
#pragma omp flush (UserQuit)
				if (!UserQuit){
// 					VolumeList[i].DoIntegration(IntResolution, IntegrateVolume);
					VolumeList[i].DoIntegrationNew(IntResolution, IntegrateVolume);
#pragma omp atomic
					NumComplete++;
// 					if (VolumeList[i].GetNumSides() == 3){
// 						VolumeList[i].DoIntegration(IntegrateVolume, stuW);
// 					}
// 					else if (VolumeList[i].GetNumSides() == 4){
// 
// 					}
					// 					TecUtilDataLoadEnd();

// 					TecUtilDrawGraphics(TRUE);
// 					TecUtilStatusSuspend(FALSE);
// 					TecUtilStatusFinishPercentDone();
// 					TecUtilLockFinish(AddOnID);
// 					return IsOk;
				}
// 				break;
			}

// #ifdef _DEBUG
// 			TecUtilDrawGraphics(TRUE);
// 			TecUtilStatusSuspend(FALSE);
// 			TecUtilStatusFinishPercentDone();
// 			TecUtilLockFinish(AddOnID);
// 			return IsOk;
// #endif // _DEBUG

			TecUtilDataLoadEnd();

			/*
			*	Now distribute the sphere integral among the GBs according to the
			*	percentage of surface area on the the sphere they occupy
			*/

			vector<double> SphereTriangleAreas;
			vector<int> ElemNums(VolumeList.size() - 1);

			vector<vector<double> > SphereElemIntVals = VolumeList[0].GetTriSphereIntValsByElem(&SphereTriangleAreas);
			vector<bool> ElemIncluded(SphereElemIntVals.size(), false);
			for (int i = 1; i < VolumeList.size() && VolumeList[i].IntResultsReady(); ++i){
				int ElemNum = stoi(AuxDataZoneGetItem(VolumeList[i].GetZoneNum(), CSMAuxData.GBA.ElemNum));
				ElemNums[i - 1] = ElemNum - 1;
				ElemIncluded[ElemNum - 1] = true;
				vector<double> IntVals = VolumeList[i].GetIntResults();
				for (int j = 0; j < SphereElemIntVals[ElemNum - 1].size(); ++j){
					SphereElemIntVals[ElemNum - 1][j] += IntVals[j];
				}
			}
			vector<vector<double> > NormalizedValues = SphereElemIntVals;
			vector<double> TotalList(SphereElemIntVals[0].size(), 0.0),
				TotalNormlizedList(SphereElemIntVals[0].size(), 0.0),
				IntScaleFactors(SphereElemIntVals[0].size());
// 			for (int i = 0; i < SphereElemIntVals.size(); ++i){
			for (int i : ElemNums){
				for (int j = 0; j < SphereElemIntVals[i].size(); ++j){
					NormalizedValues[i][j] /= SphereTriangleAreas[i];
					TotalNormlizedList[j] += NormalizedValues[i][j];
					TotalList[j] += SphereElemIntVals[i][j];
				}
			}
			for (int i = 0; i < SphereElemIntVals[0].size(); ++i){
				IntScaleFactors[i] = TotalList[i] / TotalNormlizedList[i];
			}
			// 			for (int i = 0; i < SphereElemIntVals.size(); ++i){
			for (int i : ElemNums){
				vector<double> TmpVec = SphereElemIntVals[i];
				SphereElemIntVals[i] = vector<double>();
				SphereElemIntVals[i].reserve(3 * TmpVec.size());
				for (int j = 0; j < TmpVec.size(); ++j){
					SphereElemIntVals[i].push_back(TmpVec[j]);
					SphereElemIntVals[i].push_back(NormalizedValues[i][j]);
					SphereElemIntVals[i].push_back(NormalizedValues[i][j] * IntScaleFactors[j]);
				}
			}


			TecUtilDataLoadBegin();
// 			vector<FieldData_pa> SphereResultVarRefs(NewVarNums.size());
// 			FieldDataType_e FDJunk;
// 			for (int i = 0; i < NewVarNums.size() && IsOk; ++i){
// 				SphereResultVarRefs[i] = TecUtilDataValueGetWritableNativeRef(VolumeList[0].GetZoneNum(), NewVarNums[i]);
// 				IsOk = VALID_REF(SphereResultVarRefs[i]);
// 			}
// 			if (!IsOk){
// 				TecUtilDialogErrMsg("Failed to get pointers to variables for integration.");
// 			}

			vector<vector<FieldDataPointer_c> > Ptrs(VolumeList.size());
			for (int i = 0; i < VolumeList.size(); ++i){
				if (VolumeList[i].IntResultsReady()){
					for (int j = 0; j < NewVarNums.size(); ++j){
						Ptrs[i].push_back(FieldDataPointer_c());
						Ptrs[i][j].InitializeWritePtr(VolumeList[i].GetZoneNum(), NewVarNums[j]);
					}
				}
			}

			// Write cell-centered integration values to sphere and FE zones
			for (int e = 0; e < ElemNums.size(); ++e){
				int i = ElemNums[e];
				for (int j = 0; j < SphereElemIntVals[i].size(); ++j){
					Ptrs[0][j].Write(i, SphereElemIntVals[i][j]);
					if (Ptrs.size() > e){
						for (int k = 0; k < Ptrs[e + 1][j].Size(); k++){
							Ptrs[e + 1][j].Write(k, SphereElemIntVals[i][j]);
						}
					}
				}
			}

			TecUtilDataLoadEnd();

// 			vector<vector<double> > AllIntValues = VolumeList[0].GetTriSphereIntValsByElem();
// 			vector<vector<double> > AllIntValues(VolumeList.size(), vector<double>(IntVarNumList.size() + 1 - int(IntegrateVolume), 0.0));


			if (PrintOutput){
				TotalList = vector<double>(NewVarNums.size(), 0.0);
				OutFile.open(OutFileName.c_str(), ios::out | ios::app);
				if (!OutFile.is_open()){
					TecUtilDialogErrMsg("Failed to open output file. (are you accessing it with another program maybe?)");
				}
				else{
					OutFile << endl << "Atom: " << AtomNameList[AtomNum] << endl <<
						"Elem #";
					for (auto it = IntVarNameList.cbegin(); it != IntVarNameList.cend(); it++)
						OutFile << "," << *it;
					if (IntegrateVolume)
						OutFile << ",Vol [Len.^3]";
					OutFile << endl;

					for (int i = 1; i < VolumeList.size(); ++i){
						if (VolumeList[i].IntResultsReady()){
							OutFile << VolumeList[i].GetZoneNum();
							// 						vector<double> IntResults = VolumeList[i].GetIntResults();

							for (int j = 0; j < SphereElemIntVals[i - 1].size(); ++j){
								// 							AllIntValues[i][j] += IntResults[j];
								OutFile << "," << SphereElemIntVals[i - 1][j];
								TotalList[j] += SphereElemIntVals[i - 1][j];

								// 							SphereResultVarPtrs[j].Write(i, IntResults[j]);
								// 							int IJK[3];
								// 							TecUtilZoneGetIJK(VolumeList[0].GetZoneNum(), &IJK[0], &IJK[1], &IJK[2]);
								// 							TecUtilDataValueSetByZoneVar(VolumeList[0].GetZoneNum(), NewVarNums[j], i, IntResults[j]);

								//AuxDataZoneSetItem(VolumeList[i].GetZoneNum(), AuxDataNames[j], to_string(AllIntValues[i][j]));
							}
							ContourMinMax[0] = MIN(ContourMinMax[0], SphereElemIntVals[i - 1][0]);
							ContourMinMax[1] = MAX(ContourMinMax[1], SphereElemIntVals[i - 1][0]);
							OutFile << endl;
						}
					}
					OutFile << endl << "ALL";
					for (int i = 0; i < IntVarNumList.size() + int(IntegrateVolume); ++i){
						OutFile << "," << TotalList[i];

						//AuxDataZoneSetItem(VolumeList[0].GetZoneNum(), AuxDataNames[i], to_string(TotalList[i]));
					}
					OutFile << endl << endl;
					// 				OutFile.close();
				}
			}


			/*
			 *	Now need to write the integration values to the sphere.
			 *	Since the sphere has node-centered data (could make it cell centered, but then I
			 *	don't know how to position the sphere...), need to get a volume/weighted average
			 *	of each nodes elements.
			 */
// 			for (int i = 1; i <= NumZones; ++i){
// 				if (TecUtilZoneIsOrdered(i)){
// 					string GPElemNums;
// 					if (AuxDataZoneGetItem(i, GBAGPElemNums, GPElemNums)){
// 						string GPNodeNum = AuxDataZoneGetItem(i, GBAGPNodeNum);
// 						vector<vector<double> > NodeIntVals;
// 						vector<double> NodeIntTotals(TotalList.size(), 0.0);
// 						for (int j = 0; j < VolumeList.size(); ++j){
// 							string FEElemNum;
// 							if (AuxDataZoneGetItem(VolumeList[j].GetZoneNum(), GBAElemNum, FEElemNum)){
// 								if (SearchVectorForString(SplitString(GPElemNums, ','), FEElemNum, false) >= 0){
// 									NodeIntVals.push_back(VolumeList[j].GetIntResults());
// 								}
// 							}
// 						}
// 						for (int j = 0; j < NodeIntTotals.size(); ++j){
// 							for (int k = 0; k < NodeIntVals.size(); ++k){
// 								NodeIntTotals[j] += NodeIntVals[k][j];
// 							}
// 						}
// 
// 						for (int j = 0; j < NodeIntTotals.size(); ++j){
// 							double WeightedVal = 0.0;
// 							for (int k = 0; k < NodeIntVals.size(); ++k){
// 								WeightedVal += NodeIntVals[k][j] * (NodeIntVals[k][NodeIntTotals.size() - 1] / NodeIntTotals[NodeIntTotals.size() - 1]);
// 							}
// 							SphereResultVarPtrs[j].Write(stoi(GPNodeNum) - 1, WeightedVal);
// 						}
// 					}
// 				}
// 			}
		}

		if (DoTiming){
			TimeFinish = high_resolution_clock::now();
			duration<double> TimeSpan = duration_cast<duration<double>>(TimeFinish - TimeStartLoop);
			TimeStart += TimeFinish - TimeStartLoop;
// 			TecUtilDialogMessageBox(string(string("Runtime: ") + to_string(TimeSpan.count()) + string(" seconds.")).c_str(), MessageBoxType_Information);
			if (PrintOutput)
				OutFile << "Precision," << IntResolution << '\n'
					<< "Wall time[s]," << TimeSpan.count();

// 			AuxDataZoneSetItem(VolumeList[0].GetZoneNum(), CSMAuxData.GBA.IntPrecision, to_string(IntResolution));
// 			AuxDataZoneSetItem(VolumeList[0].GetZoneNum(), CSMAuxData.GBA.IntWallTime, to_string(TimeSpan.count()));
		}

// 		Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// 		for (auto const & i : VolumeList){
// 			for (int j = 0; j < i.GetNumSides(); ++j){
// 				TecUtilSetAddMember(TmpSet, i.GetGPZoneNum(j), TRUE);
// 			}
// 		}
// 		TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
// 		TecUtilSetDealloc(&TmpSet);

// 		for (auto & i : VolumeList)
// 			IsOk = TecUtilSetAddMember(ZoneSet, i.GetZoneNum(), FALSE);
// 			
		TecUtilSetDealloc(&ZoneSet);

		StatusDrop(AddOnID);
	}

	if (DoTiming){
		TimeFinish = high_resolution_clock::now();
		duration<double> TimeSpan = duration_cast<duration<double>>(TimeFinish - TimeStart);
		TimeStart += TimeFinish - TimeStartLoop;
		// 			TecUtilDialogMessageBox(string(string("Runtime: ") + to_string(TimeSpan.count()) + string(" seconds.")).c_str(), MessageBoxType_Information);
		if (PrintOutput)
			OutFile << "\n\nTotal Wall time[s]," << TimeSpan.count();
	}
	if (PrintOutput)
		OutFile.close();

// 	if (IsOk){
// 		/*
// 		 *	Set 8th contour group with levels for coloring the 
// 		 *	sphere zones, then set their contour coloring to
// 		 *	the 8th group.
// 		 */
// 
// 		int NumContourLevels = 30;
// // 		vec ContourLevels = logspace(ContourMinMax[0], ContourMinMax[1], NumContourLevels);
// 		vec ContourLevels = linspace(ContourMinMax[0], ContourMinMax[1], NumContourLevels);
// // 		vector<double> ContourLevels(50);
// // 		for (int i = 0; i < NumContourLevels; ++i)
// // 			ContourLevels[i] = ContourMinMax[0] + (double)i / (double)(NumContourLevels - 1) * (ContourMinMax[1] - ContourMinMax[0]);
// 
// 		ArgList_pa TempArgList = TecUtilArgListAlloc();
// 		TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 8);
// 		TecUtilArgListAppendInt(TempArgList, SV_VAR, NewVarNums[0]);
// 		TecUtilContourSetVariableX(TempArgList);
// 		TecUtilArgListClear(TempArgList);
// 
// 		TecUtilArgListAppendInt(TempArgList, SV_CONTOURLABELACTION, ContourLabelAction_DeleteAll);
// 		TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 8);
// 		TecUtilContourLabelX(TempArgList);
// 		TecUtilArgListClear(TempArgList);
// 
// 		TecUtilArgListAppendInt(TempArgList, SV_CONTOURLEVELACTION, ContourLevelAction_New);
// 		TecUtilArgListAppendInt(TempArgList, SV_CONTOURGROUP, 8);
// 		TecUtilArgListAppendInt(TempArgList, SV_NUMVALUES, NumContourLevels);
// 		TecUtilArgListAppendArray(TempArgList, SV_RAWDATA, ContourLevels.memptr());
// 		TecUtilContourLevelX(TempArgList);
// 		TecUtilArgListClear(TempArgList);
// 
// 		TecUtilArgListAppendString(TempArgList, SV_P1, SV_GLOBALCONTOUR);
// 		TecUtilArgListAppendString(TempArgList, SV_P2, SV_LEGEND);
// 		TecUtilArgListAppendString(TempArgList, SV_P3, SV_SHOW);
// 		TecUtilArgListAppendInt(TempArgList, SV_OFFSET1, 8);
// 		TecUtilArgListAppendArbParam(TempArgList, SV_IVALUE, FALSE);
// 		TecUtilStyleSetLowLevelX(TempArgList);
// 		TecUtilArgListClear(TempArgList);
// 
// 		TecUtilZoneSetContour(SV_SHOW, ZoneSet, 0.0, TRUE);
// 		TecUtilZoneSetContour(SV_CONTOURGROUP, ZoneSet, 0.0, ContourColoring_Group8);
// 		TecUtilZoneSetContour(SV_CONTOURTYPE, ZoneSet, 0.0, ContourType_Overlay);
// 		TecUtilZoneSetContour(SV_FLOODCOLORING, ZoneSet, 0.0, ContourColoring_Group8);
// 
// 	}

	TecUtilLockFinish(AddOnID);

	return IsOk;
}

int GetRecommendedIntPrecision(){
	int RecPrecision = 2;

	if (TecUtilDataSetIsAvailable()){
		int ZoneNum = ZoneNumByName("Critical Points");
		if (ZoneNum > 0 && TecUtilZoneIsOrdered(ZoneNum)){
			double MinCPSpacing = DBL_MAX;
			int NumCPs;
			vec3 iCP, jCP;
			TecUtilZoneGetIJK(ZoneNum, &NumCPs, nullptr, nullptr);
			TecUtilDataLoadBegin();
			vector<int> XYZVarNums = { 1, 2, 3 };
			vector<FieldData_pa> CPXYZRefs(3);
			for (int i = 0; i < 3; ++i){
				CPXYZRefs[i] = TecUtilDataValueGetReadableNativeRef(ZoneNum, XYZVarNums[i]);
				if (!VALID_REF(CPXYZRefs[i]))
					return RecPrecision;
			}
			for (int i = 1; i <= NumCPs; ++i){
				for (int ii = 0; ii < 3; ++ii){
					iCP[ii] = TecUtilDataValueGetByRef(CPXYZRefs[ii], i);
				}
				for (int j = i + 1; j <= NumCPs; ++j){
					for (int ii = 0; ii < 3; ++ii){
						jCP[ii] = TecUtilDataValueGetByRef(CPXYZRefs[ii], j);
					}
					MinCPSpacing = MIN(MinCPSpacing, Distance(iCP, jCP));
				}
			}
			TecUtilDataLoadEnd();
			double DelXYZMag = sqrt(sum(square(GetDelXYZ_Ordered3DZone(XYZVarNums, ZoneNumByName("Full Volume")))));
			RecPrecision = MIN(16, MAX(0,static_cast<int>(std::log2(DelXYZMag / (MinCPSpacing * GBADefaultRadialSphereApprxRadius * PI * 2.0 / 40.0)))));
		}
	}

	return RecPrecision;
}

