
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

#include "CSM_DATA_TYPES.h"
#include "ZONEVARINFO.h"
#include "MACROFUNCTIONS.h"
#include "CSVFILE.h"

#include "GTAENGINE.h"

#include <armadillo>
using namespace arma;

// #define NumberOfSlices 100
// #define NumberOfStreamtraces 100
#define NumberOfContourLevels 20
// #define DistanceRadiusFactor 0.025
#define StreamtraceNumStepsFactor 25
#define MaxNumNewZones 1500


// #define KeepNewZones 1




Boolean_t GTARunGTA(std::vector<vec3> *Pts){
	/*
		Here calculate the three normal unit vectors
		for the three points chosen. The third vector,
		v3, is the vector normal to the plane we want,
		so use its x,y,z values as the normals in the 
		create slice call.

		v1 and v2 are the normal vectors that are used
		to position streamtraces.
	*/

	if (TecGUIDialogIsUp(Dialog1Manager))
		TecGUIDialogDrop(Dialog1Manager);

	if (Pts->size() < 2){
		TecUtilDialogErrMsg("Must select at least two points");
		return FALSE;
	}

	double CutoffValue,
		DistanceRadiusFactor,
		SliceBegAng,
		SliceEndAng,
		StreamtraceBegAng,
		StreamtraceEndAng;

	std::vector<double*> TFDoubleInputValues = {
		&CutoffValue,
		&DistanceRadiusFactor,
		&SliceBegAng,
		&SliceEndAng,
		&StreamtraceBegAng,
		&StreamtraceEndAng
	};
	std::vector<LgIndex_t> TFInputFieldNums = {
		TFCutoff_TF_S1,
		TFRadOffVal_TF_S1,
		TFSlcBegAng_TF_S1,
		TFSlcEndAng_TF_S1,
		TFStBegAng_TF_S1,
		TFStEndAng_TF_S1
	};
	std::vector<std::string> TFErrMsgs = { "Invalid cutoff value",
		"Invalid angle offset value",
		"Invalid starting slice angle",
		"Invalid ending slice angle",
		"Invalid starting streamtrace angle",
		"Invalid ending streamtrace angle" };

	for (int i = 0; i < TFDoubleInputValues.size(); ++i){
		if (!TecGUITextFieldGetDouble(TFInputFieldNums[i], TFDoubleInputValues[i])){
			TecUtilDialogErrMsg(TFErrMsgs[i].c_str());
			return FALSE;
		}
	}

	if (SliceEndAng - SliceBegAng <= 0){
		TecUtilDialogErrMsg("Slice angles invalid, ending must be greater than starting");
		return FALSE;
	}
	if (StreamtraceEndAng - StreamtraceBegAng <= 0){
		TecUtilDialogErrMsg("Streamtrace angles invalid, ending must be greater than starting");
		return FALSE;
	}

	TFDoubleInputValues.clear();
	TFInputFieldNums.clear();
	TFErrMsgs.clear();

	TecUtilLockStart(AddOnID);

	TecUtilStatusStartPercentDone("Initializing", TRUE, TRUE);
	Boolean_t IsOk = TecUtilStatusCheckPercentDone(0);
	TecUtilStatusSuspend(TRUE);

	Set_pa ActiveZones = TecUtilSetAlloc(TRUE);
	TecUtilZoneGetActive(&ActiveZones);

	if (Pts->size() < 3){
		Pts->push_back(Pts->at(0));
		Pts->at(2)[1] += 0.1;
	}

	vec3 v1 = Pts->at(1) - Pts->at(0);
	vec3 v2 = Pts->at(2) - Pts->at(1);
	vec3 v3 = cross(v1, v2);

	v1 = normalise( v1);
	v3 = normalise( v3);
	v2 = normalise(cross(v3, v1));

	if (IsOk){
		ArgList_pa argList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(argList, SV_P1, SV_INTERFACE);
		TecUtilArgListAppendString(argList, SV_P2, SV_AUTOREDRAWISACTIVE);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, FALSE);
		TecUtilStyleSetLowLevelX(argList);
	}

	/*
		Get the variable numbers of Rho and whatever variable the user
		wants to get data for.
	*/
	std::vector<std::string> RhoStrs = { "Electron Density", "Rho", "rho" };
	EntIndex_t RhoVarNum = VarNumByNameList(RhoStrs);
	EntIndex_t CutoffVarNum = TecGUIOptionMenuGet(OPTCutoffVar_OPT_S1);
	std::vector<EntIndex_t> FunctVarNums;
	LgIndex_t* SelectedVars;
	LgIndex_t NumSelected;
	TecGUIListGetSelectedItems(MLSelectVar_MLST_S1, &SelectedVars, &NumSelected);
	for (int i = 0; i < NumSelected; ++i)
		FunctVarNums.push_back(SelectedVars[i]);
	TecUtilArrayDealloc(reinterpret_cast<void**>(&SelectedVars));
	if (IsOk)
		IsOk = (RhoVarNum > 0 && FunctVarNums.size() > 0);
	EntIndex_t VolZoneNum = -1;

	Boolean_t KeepNewZones = !TecGUIToggleGet(TGLDeleteZon_TOG_S1);
	Boolean_t HideNewZones = TecGUIToggleGet(TGLHideZones_TOG_S1);
	Boolean_t ShowOnlyNewZones = TecGUIToggleGet(TGLShowOnly_TOG_S1);

	//TODO:need to update program so that data is stored and outputted correctly, since starting and ending slice/streamtrace angles can now be specified

	double SliceAngRange = SliceEndAng - SliceBegAng;
	double StreamtraceAndRange = StreamtraceEndAng - StreamtraceBegAng;
	unsigned int NumberOfSlices = (unsigned int)(180. / atof(TecGUIOptionMenuGetString(OptSliceAngl_OPT_S1, TecGUIOptionMenuGet(OptSliceAngl_OPT_S1))));
	unsigned int NumberOfStreamtraces = (unsigned int)(360. / atof(TecGUIOptionMenuGetString(OptStreamAng_OPT_S1, TecGUIOptionMenuGet(OptStreamAng_OPT_S1))));

	unsigned int NumNewZones = (NumberOfStreamtraces * 4 + 2) * NumberOfSlices;
	unsigned int PercentTotal = NumNewZones + NumNewZones * NumberOfSlices * (0.1 * (FunctVarNums.size() + 1));

	if (KeepNewZones && NumNewZones > MaxNumNewZones && NumberOfSlices > 1){
		IsOk = TecUtilDialogMessageBox("Because too many zones will be created, only the first slice will be preserved", MessageBoxType_Information);
	}

	/*
		Get the full volume zone number
	*/
	if (IsOk){
		//EntIndex_t VolZoneNum = ZoneNumByNameList(std::vector<std::string>(1, "Full Volume Zone"));
		VolZoneNum = ZoneNumByName(std::string("Full Volume Zone"));
		IsOk = (VolZoneNum > 0 && TecUtilZoneIsOrdered(VolZoneNum));
	}
	Set_pa VolZoneSet = TecUtilSetAlloc(TRUE);

	/*
		Activate ONLY full volume zone
	*/
	if (IsOk){
		if (VolZoneSet != NULL)
			TecUtilSetAddMember(VolZoneSet, VolZoneNum, TRUE);
		if (!TecUtilSetIsEmpty(VolZoneSet))
			TecUtilZoneSetActive(VolZoneSet, AssignOp_Equals);
		IsOk = (TecUtilZoneIsActive(VolZoneNum) && TecUtilZoneIsEnabled(VolZoneNum));
	}

	LgIndex_t IJK[3];
	LgIndex_t MaxIJK = 0;
	if (IsOk){
		TecUtilZoneGetIJK(VolZoneNum, &IJK[0], &IJK[1], &IJK[2]);
		for (int i = 0; i < 3; ++i)
			if (IJK[i] > MaxIJK)
				MaxIJK = IJK[i];
	}
	

	/*
		Get values for contour levels of the slice to be created.
		Based off volume zone rather than slice so that it's good
		for all slices.
	*/
	double RhoMax = -1, RhoMin = -1;
	std::vector<double> RhoContourLevels;
	if (IsOk){
		/*
			Get max and min rho in the volume
		*/
		TecUtilVarGetMinMax(RhoVarNum, &RhoMin, &RhoMax);
		RhoContourLevels.reserve(50);
// 		LogDataValues(RhoContourLevels, RhoMin, RhoMax, 0.1);
		X2DataValues(RhoContourLevels, RhoMin, RhoMax);
	}
	double SliceRadStep = 180. / (double)NumberOfSlices * DEG2RAD;
	double StreamtraceRadStep = 360. / (double)NumberOfStreamtraces * DEG2RAD;

	double Radius = Distance(v1, v2) * DistanceRadiusFactor;

	std::vector<std::vector<std::vector<ImportType_t> > > ResultsListFull;
	/*
	 *	The extra space is for the volume integration that's always included
	 */
	ResultsListFull.resize(FunctVarNums.size() + 1);

	Set_pa AllSlicesSet, AllStreamtracesSet, AllSurfacesSet;
	Set_pa FirstSliceSet, FirstStreamtracesSet, FirstSurfacesSet;
	AllSlicesSet = TecUtilSetAlloc(TRUE);
	AllStreamtracesSet = TecUtilSetAlloc(TRUE);
	AllSurfacesSet = TecUtilSetAlloc(TRUE);

	FirstSliceSet = TecUtilSetAlloc(TRUE);
	FirstStreamtracesSet = TecUtilSetAlloc(TRUE);
	FirstSurfacesSet = TecUtilSetAlloc(TRUE);

	EntIndex_t XYZVarNums[3];
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	char* AddOnFullPathWithNameCStr = TecUtilAddOnGetPath(AddOnID);
	char* AddOnFullPathCStr = TecUtilGetBasePath(AddOnFullPathWithNameCStr);
	TecUtilStringDealloc(&AddOnFullPathWithNameCStr);


	/*
	 *	Getting path of temporary file used to store integration results
	 */
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	static const std::string Slash = "\\";
	Boolean_t IsWindows = TRUE;
#else
	static const std::string Slash = "/";
	Boolean_t IsWindows = FALSE;
#endif
	std::stringstream ResultsFileNameSS;
	ResultsFileNameSS << AddOnFullPathCStr << Slash << "TempIntegrationResults.txt";

	if (IsWindows){
		std::string TempString = ResultsFileNameSS.str();
		ResultsFileNameSS.clear();
		ResultsFileNameSS.str(std::string());
		size_t StrPos = TempString.find_first_of(Slash.c_str());
		while (StrPos != std::string::npos){
			ResultsFileNameSS << TempString.substr(0, StrPos) << Slash << Slash;
			TempString = TempString.substr(StrPos + 1, TempString.length());
			StrPos = TempString.find_first_of(Slash.c_str());
		}
		ResultsFileNameSS << TempString;
	}

	/*
	 *	Make variables that will be used for integrations.
	 *	These are made here and are passive because the
	 *	"Data -> Alter -> Specify Equations" option that we
	 *	use doesn't let you modify a variable unless it
	 *	exists for all zones.
	 */

	ArgList_pa VarArgs = TecUtilArgListAlloc();
	Set_pa VarSet = TecUtilSetAlloc(TRUE);
	Set_pa RadiusVarSet = TecUtilSetAlloc(TRUE);
	Set_pa VolVarSet = TecUtilSetAlloc(TRUE);
	Set_pa IntVarSet = TecUtilSetAlloc(TRUE);

	std::vector<FieldDataType_e> DataTypes(TecUtilDataSetGetNumZones(), FieldDataType_Double);
	
	TecUtilArgListAppendString(VarArgs, SV_NAME, "Radius");
	TecUtilArgListAppendArray(VarArgs, SV_VARDATATYPE, DataTypes.data());
	TecUtilArgListAppendInt(VarArgs, SV_DEFERVARCREATION, TRUE);
	IsOk = TecUtilDataSetAddVarX(VarArgs);
	int RadVarNum;
	int VolVarNum;
	int IntVarNum;

	if (IsOk){
		RadVarNum = TecUtilDataSetGetNumVars();
		TecUtilSetAddMember(VarSet, RadVarNum, TRUE);
		TecUtilSetAddMember(RadiusVarSet, RadVarNum, TRUE);
		TecUtilArgListClear(VarArgs);
		TecUtilArgListAppendString(VarArgs, SV_NAME, "WedgeVolume");
		TecUtilArgListAppendArray(VarArgs, SV_VARDATATYPE, DataTypes.data());
		TecUtilArgListAppendInt(VarArgs, SV_DEFERVARCREATION, TRUE);
		IsOk = TecUtilDataSetAddVarX(VarArgs);
	}
	if (IsOk){
		VolVarNum = TecUtilDataSetGetNumVars();
		TecUtilSetAddMember(VarSet, VolVarNum, TRUE);
		TecUtilSetAddMember(VolVarSet, VolVarNum, TRUE);
		TecUtilArgListClear(VarArgs);
		TecUtilArgListAppendString(VarArgs, SV_NAME, "IntVar");
		TecUtilArgListAppendArray(VarArgs, SV_VARDATATYPE, DataTypes.data());
		TecUtilArgListAppendInt(VarArgs, SV_DEFERVARCREATION, TRUE);
		IsOk = TecUtilDataSetAddVarX(VarArgs);
	}
	TecUtilArgListDealloc(&VarArgs);
	if (IsOk){
		IntVarNum = TecUtilDataSetGetNumVars();
		TecUtilSetAddMember(VarSet, IntVarNum, TRUE);
		TecUtilSetAddMember(IntVarSet, IntVarNum, TRUE);
		ArgList_pa StateChangeArgs = TecUtilArgListAlloc();
		TecUtilArgListAppendInt(StateChangeArgs, SV_STATECHANGE, StateChange_VarsAdded);
		TecUtilStateChangedX(StateChangeArgs);
		TecUtilArgListDealloc(&StateChangeArgs);
	}

	unsigned int ZoneCount = 0;
	std::stringstream ss1, ss2;


	/*
		From here we'll iterate and do this stuff for each slice
		created.
	*/
	for (unsigned int SliceNum = 0; SliceNum < NumberOfSlices; ++SliceNum){

		ss1.str(std::string());
		ss1.clear();
		ss1 << "Slice " << SliceNum + 1 << " of " << NumberOfSlices << ": ";

		/*
		 *	Activate only the full volume zone
		 */
		if (IsOk){
			if (VolZoneSet != NULL)
				TecUtilSetAddMember(VolZoneSet, VolZoneNum, TRUE);
			if (!TecUtilSetIsEmpty(VolZoneSet))
				TecUtilZoneSetActive(VolZoneSet, AssignOp_Equals);
			IsOk = (TecUtilZoneIsActive(VolZoneNum) && TecUtilZoneIsEnabled(VolZoneNum));
		}

		/*
		 *	Make sure right frame is on top
		 */
		EntIndex_t FrameNum = 1;
		while (TecUtilFrameGetPlotType() != PlotType_Cartesian3D && IsOk){
			if (FrameNum <= TecUtilFrameGetCount()){
				IsOk = TecUtilFrameActivateByNumber(FrameNum);
				++FrameNum;
			}
		}

		char* TmpChar;
		TecUtilZoneGetName(1, &TmpChar);

		vec4 SliceVec4 = 
			RotationMatrix(SliceRadStep * SliceNum, v1) 
			* join_cols(v3, ones<vec>(1));
		vec3 PlaneNormVec;
		PlaneNormVec << SliceVec4[0] << SliceVec4[1] << SliceVec4[2];

		SliceVec4 =
			RotationMatrix(SliceRadStep * SliceNum, v1)
			* join_cols(v2, ones<vec>(1));
		vec3 PlaneParVec;
		PlaneParVec << SliceVec4[0] << SliceVec4[1] << SliceVec4[2];

		/*
			Create a new slice at point 0
		*/
		if (IsOk){
			IsOk = TecUtilCreateSliceZoneFromPlane(SliceSource_VolumeZones,
				Pts->at(0)[0],
				Pts->at(0)[1],
				Pts->at(0)[2],
				PlaneNormVec[0],
				PlaneNormVec[1],
				PlaneNormVec[2]);
		}		
// 		if (IsOk){
// 			IsOk = TecUtilCreateSliceZoneFromPlane(SliceSource_VolumeZones,
// 				Pts->at(0)[0],
// 				Pts->at(0)[1],
// 				Pts->at(0)[2],
// 				v3[0],
// 				v3[1],
// 				v3[2]);
// 		}
		/*
			Get number of new slice zone
		*/
		EntIndex_t SlizeZoneNum;
		if (IsOk)
			SlizeZoneNum = TecUtilDataSetGetNumZones();
		if (IsOk)
			IsOk = (SlizeZoneNum > 0);
		Set_pa SliceSet = TecUtilSetAlloc(TRUE);

		/*
			Activate ONLY new slice zone
		*/
		if (IsOk){
			ZoneCount++;
			if (SliceSet != NULL){
				IsOk = TecUtilSetAddMember(SliceSet, SlizeZoneNum, TRUE);
				if (IsOk)
					IsOk = TecUtilSetAddMember(AllSlicesSet, SlizeZoneNum, TRUE);
				if (SliceNum == 0 && IsOk)
					IsOk = TecUtilSetAddMember(FirstSliceSet, SlizeZoneNum, TRUE);
			}
			if (!TecUtilSetIsEmpty(SliceSet)){
				TecUtilZoneSetActive(SliceSet, AssignOp_Equals);
			}
			IsOk = (TecUtilZoneIsActive(SlizeZoneNum) && TecUtilZoneIsEnabled(SlizeZoneNum));
		}

		if (IsOk){
			//	Modify Contour Properties

			//	ContourShow
			TecUtilZoneSetContour(SV_SHOW, SliceSet, 0.0, TRUE);

			//	ContourType
			TecUtilZoneSetContour(SV_CONTOURTYPE, SliceSet, 0.0, ContourType_Lines);

			//	FloodColoring
			TecUtilZoneSetContour(SV_FLOODCOLORING, SliceSet, 0.0, ContourColoring_Group8);
			TecUtilZoneSetContour(SV_LINECONTOURGROUP, SliceSet, 0.0, ContourColoring_Group8);
			TecUtilZoneSetContour(SV_COLOR, SliceSet, 0.0, Blue_C);

			/*
			 *	Disable other styles (mesh, scatter, etc...)
			 */
			TecUtilZoneSetMesh(SV_SHOW, SliceSet, 0.0, FALSE);
			TecUtilZoneSetScatter(SV_SHOW, SliceSet, 0.0, FALSE);

			ArgList_pa ContourArgs = TecUtilArgListAlloc();
			TecUtilArgListAppendInt(ContourArgs, SV_CONTOURGROUP, 8);
			TecUtilArgListAppendInt(ContourArgs, SV_VAR, RhoVarNum);
			SetValueReturnCode_e Result = TecUtilContourSetVariableX(ContourArgs);

			TecUtilArgListClear(ContourArgs);

			IsOk = (Result == SetValueReturnCode_Ok || Result == SetValueReturnCode_DuplicateValue);

			if (IsOk){
				TecUtilArgListAppendInt(ContourArgs, SV_CONTOURLEVELACTION, ContourLevelAction_New);
				TecUtilArgListAppendInt(ContourArgs, SV_CONTOURGROUP, 8);
				TecUtilArgListAppendInt(ContourArgs, SV_NUMVALUES, (LgIndex_t)RhoContourLevels.size());
				TecUtilArgListAppendArray(ContourArgs, SV_RAWDATA, RhoContourLevels.data());
				IsOk = TecUtilContourLevelX(ContourArgs);
			}
			TecUtilArgListDealloc(&ContourArgs);
		}

		/*
		 *	Calculate the Radius variable for the slice.
		 *	Doing it here lets it carry through automatically
		 *	to streamtrace and polygrid zones.
		 */

		ZoneType_e SliceType = TecUtilZoneGetType(SlizeZoneNum);

		if (IsOk){
			IsOk = PopulateRadiusVar(RadVarNum, XYZVarNums, CutoffVarNum, CutoffValue, SlizeZoneNum, Pts->at(0), v1);
		}

		std::stringstream ss;
// 		if (IsOk){
// 			ss << "V" << RadVarNum << "=abs(" << Pts->at(0)[1] << "-Y" << ")";
// 			IsOk = TecUtilDataAlter(ss.str().c_str(), SliceSet,
// 				1, 0, 1,
// 				1, 0, 1,
// 				1, 0, 1,
// 				FieldDataType_Double);
// 		}
		if (IsOk){
			ArgList_pa StateChangeArgs = TecUtilArgListAlloc();
			TecUtilArgListAppendInt(StateChangeArgs, SV_STATECHANGE, StateChange_VarsAltered);
			TecUtilArgListAppendSet(StateChangeArgs, SV_VARLIST, RadiusVarSet);
			TecUtilArgListAppendSet(StateChangeArgs, SV_ZONELIST, SliceSet);
			TecUtilStateChangedX(StateChangeArgs);
			TecUtilArgListDealloc(&StateChangeArgs);
		}

		if (IsOk){
			ss.str(std::string());
			ss.clear();
			ss << "V" << VolVarNum << "=" << SliceRadStep << "*V" << RadVarNum;
			IsOk = TecUtilDataAlter(ss.str().c_str(), SliceSet,
				1, 0, 1,
				1, 0, 1,
				1, 0, 1,
				FieldDataType_Double);
		}
		if (IsOk){
			ArgList_pa StateChangeArgs = TecUtilArgListAlloc();
			TecUtilArgListAppendInt(StateChangeArgs, SV_STATECHANGE, StateChange_VarsAltered);
			TecUtilArgListAppendSet(StateChangeArgs, SV_VARLIST, VolVarSet);
			TecUtilArgListAppendSet(StateChangeArgs, SV_ZONELIST, SliceSet);
			TecUtilStateChangedX(StateChangeArgs);
			TecUtilArgListDealloc(&StateChangeArgs);
		}

		/*
			Now can start to seed streamtraces around point 0 on the slice.
		*/

		//	Set streamtrace options 
		if (IsOk){
			ArgList_pa argList = TecUtilArgListAlloc();
			TecUtilArgListAppendString(argList, SV_P1, SV_STREAMTRACELAYERS);
			TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
			SetValueReturnCode_e ReturnCode = TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListDealloc(&argList);
			IsOk = (ReturnCode == SetValueReturnCode_Ok || ReturnCode == SetValueReturnCode_DuplicateValue);
		}

		if (IsOk){
			ArgList_pa argList = TecUtilArgListAlloc();
			TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALTHREEDVECTOR);
			TecUtilArgListAppendString(argList, SV_P2, SV_UVAR);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, (EntIndex_t)VarNumByName(std::string("X Density Gradient")));
			SetValueReturnCode_e ReturnCode = TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListDealloc(&argList);
			IsOk = (ReturnCode == SetValueReturnCode_Ok || ReturnCode == SetValueReturnCode_DuplicateValue);
		}

		if (IsOk){
			ArgList_pa argList = TecUtilArgListAlloc();
			TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALTHREEDVECTOR);
			TecUtilArgListAppendString(argList, SV_P2, SV_VVAR);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, (EntIndex_t)VarNumByName(std::string("Y Density Gradient")));
			SetValueReturnCode_e ReturnCode = TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListDealloc(&argList);
			IsOk = (ReturnCode == SetValueReturnCode_Ok || ReturnCode == SetValueReturnCode_DuplicateValue);
		}

		if (IsOk){
			ArgList_pa argList = TecUtilArgListAlloc();
			TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALTHREEDVECTOR);
			TecUtilArgListAppendString(argList, SV_P2, SV_WVAR);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, (EntIndex_t)VarNumByName(std::string("Z Density Gradient")));
			SetValueReturnCode_e ReturnCode = TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListDealloc(&argList);
			IsOk = (ReturnCode == SetValueReturnCode_Ok || ReturnCode == SetValueReturnCode_DuplicateValue);
		}

		if (IsOk){
			ArgList_pa argList = TecUtilArgListAlloc();
			TecUtilArgListAppendString(argList, SV_P1, SV_STREAMATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_MAXSTEPS);
			TecUtilArgListAppendArbParam(argList, SV_IVALUE, (LgIndex_t)MAX(MaxIJK * StreamtraceNumStepsFactor, 1500));
			SetValueReturnCode_e ReturnCode = TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListDealloc(&argList);
			IsOk = (ReturnCode == SetValueReturnCode_Ok || ReturnCode == SetValueReturnCode_DuplicateValue);
		}
		if (IsOk){
			ArgList_pa argList = TecUtilArgListAlloc();
			TecUtilArgListAppendString(argList, SV_P1, SV_STREAMATTRIBUTES);
			TecUtilArgListAppendString(argList, SV_P2, SV_MINCELLFRACTION);
			TecUtilArgListAppendDouble(argList, SV_DVALUE, (double)0.00001);
			SetValueReturnCode_e ReturnCode = TecUtilStyleSetLowLevelX(argList);
			TecUtilArgListDealloc(&argList);
			IsOk = (ReturnCode == SetValueReturnCode_Ok || ReturnCode == SetValueReturnCode_DuplicateValue);
		}

		/*
		 *	Seed all the streamtraces
		 */
		/*
		 *	Do outer streamtraces first, then inner
		 */
		StreamDir_e StreamDirs[2] = { StreamDir_Reverse, StreamDir_Forward };
		Set_pa StreamtraceZoneSet;
		EntIndex_t StreamtraceZoneNumBegin[2];
		EntIndex_t StreamtraceZoneNumEnd[2];
		StreamtraceZoneSet = TecUtilSetAlloc(TRUE);

		ss2.str(std::string());
		ss2.clear();
		ss2 << ss1.str() << "Streamtraces";
// 		TecGUILabelSetText(Lab1Status_LBL_S1, ss2.str().c_str());

		TecUtilDataLoadBegin();

		for (int RunNum = 0; RunNum < 2 && IsOk; ++RunNum){
			for (unsigned int StreamTraceNum = 0; StreamTraceNum <= NumberOfStreamtraces && IsOk; ++StreamTraceNum)
			{
				vec3 SeedPoint;
				vec4 TmpVec4;
				double StreamtraceTheta;
				if (StreamTraceNum == 0 || StreamTraceNum == NumberOfStreamtraces){
					SeedPoint = v1 * (-Radius * 0.2);
					StreamtraceTheta = PI / 4.;
					if (StreamTraceNum == 0)
						StreamtraceTheta *= -1.;
					TmpVec4 = join_cols(SeedPoint, ones<vec>(1));
					TmpVec4 = RotationMatrix(StreamtraceTheta, PlaneNormVec)
						* TmpVec4;

					SeedPoint = vec(vector<double>({ TmpVec4[0], TmpVec4[1], TmpVec4[2] })) + Pts->at(1);
				}
				else{
					StreamtraceTheta = StreamtraceRadStep * (double)StreamTraceNum;
					SeedPoint = v1 * Radius;
					TmpVec4 = join_cols(SeedPoint, ones<vec>(1));
					TmpVec4 = RotationMatrix(StreamtraceTheta, PlaneNormVec) * TmpVec4;

					SeedPoint = vec(vector<double>({ TmpVec4[0], TmpVec4[1], TmpVec4[2] })) + Pts->at(0);
				}

				IsOk = TecUtilStreamtraceAdd(1, Streamtrace_SurfaceLine, StreamDirs[RunNum], 
					SeedPoint[0], SeedPoint[1], SeedPoint[2], 
					1, 1, 1);

				if (IsOk)
					ZoneCount++;

				if (IsOk)
					IsOk = SetPercent(ZoneCount, PercentTotal, ss2.str().c_str());
			}
			/*
			*	Extract all streamtraces as individual zones
			*/
			if (IsOk){
				TecUtilCreateStreamZones(FALSE);
				TecUtilStreamtraceDeleteAll();
			}
			if (IsOk){
				if (RunNum == 0)
					StreamtraceZoneNumBegin[RunNum] = SlizeZoneNum + 1;
				else
					StreamtraceZoneNumBegin[RunNum] = StreamtraceZoneNumEnd[RunNum - 1] + 1;

				StreamtraceZoneNumEnd[RunNum] = TecUtilDataSetGetNumZones();
			}

			REQUIRE(StreamtraceZoneNumEnd[RunNum] - StreamtraceZoneNumBegin[RunNum] == NumberOfStreamtraces);
		}

		TecUtilDataLoadEnd();

		ss2.str(std::string());
		ss2.clear();
		ss2 << ss1.str() << "Polyline Zones";
// 		TecGUILabelSetText(Lab1Status_LBL_S1, ss2.str().c_str());

		/*
			*	Only if the streamtraces stick around (for debugging and maybe if the user wants to see it),
			*	set their styles
			*/
		for (EntIndex_t StreamtraceZoneNum = StreamtraceZoneNumBegin[0]; StreamtraceZoneNum <= StreamtraceZoneNumEnd[1] && IsOk; ++StreamtraceZoneNum)
		{
			IsOk = TecUtilSetAddMember(StreamtraceZoneSet, StreamtraceZoneNum, TRUE);
			if (IsOk)
				IsOk = TecUtilSetAddMember(AllStreamtracesSet, StreamtraceZoneNum, TRUE);
			if (SliceNum == 0 && IsOk)
				IsOk = TecUtilSetAddMember(FirstStreamtracesSet, StreamtraceZoneNum, TRUE);
		}

		if (KeepNewZones && IsOk){
			TecUtilZoneSetScatter(SV_SHOW, StreamtraceZoneSet, 0.0, FALSE);
			TecUtilZoneSetContour(SV_SHOW, StreamtraceZoneSet, 0.0, FALSE);
			TecUtilZoneSetMesh(SV_SHOW, StreamtraceZoneSet, 0.0, TRUE);
			TecUtilZoneSetMesh(SV_COLOR, StreamtraceZoneSet, 0.0, Red_C);
		}

		/*
		 *	Create the surfaces zones containing the space between sequential streamtraces
		 */
		Set_pa SurfaceZoneSet;
		SurfaceZoneSet = TecUtilSetAlloc(TRUE);
		EntIndex_t SurfaceZoneNumBegin[2];
		EntIndex_t SurfaceZoneNumEnd[2];

		TecUtilDataLoadBegin();

		for (int RunNum = 0; RunNum < 2 && IsOk; ++RunNum){
			for (EntIndex_t SurfaceNum = StreamtraceZoneNumBegin[RunNum] + 1; SurfaceNum <= StreamtraceZoneNumEnd[RunNum] && IsOk; ++SurfaceNum)
			{
				std::vector<EntIndex_t> ZoneNums = { SurfaceNum - 1, SurfaceNum };
				EntIndex_t NewZoneNum = MacroCreateZoneFromPolylines(&ZoneNums, FALSE, AddOnFullPathCStr);
				if (NewZoneNum > 0){
					IsOk = TecUtilSetAddMember(SurfaceZoneSet, NewZoneNum, TRUE);
					if (IsOk){
						IsOk = TecUtilSetAddMember(AllSurfacesSet, NewZoneNum, TRUE);
						ZoneCount++;
					}
					if (SliceNum == 0 && IsOk)
						IsOk = TecUtilSetAddMember(FirstSurfacesSet, NewZoneNum, TRUE);
				}
				else{
					TecUtilDialogErrMsg("Failed to create zone from polylines");
					IsOk = FALSE;
				}

				if (IsOk)
					IsOk = SetPercent(ZoneCount, PercentTotal, ss2.str().c_str());
			}
			if (RunNum == 0)
				SurfaceZoneNumBegin[RunNum] = StreamtraceZoneNumEnd[1] + 1;
			else
				SurfaceZoneNumBegin[RunNum] = SurfaceZoneNumEnd[RunNum - 1] + 1;

			SurfaceZoneNumEnd[RunNum] = TecUtilDataSetGetNumZones();

			REQUIRE(SurfaceZoneNumEnd[RunNum] - SurfaceZoneNumBegin[RunNum] == NumberOfStreamtraces - 1);
		}

		TecUtilDataLoadEnd();

		TecUtilZoneSetScatter(SV_SHOW, SurfaceZoneSet, 0.0, FALSE);
		TecUtilZoneSetMesh(SV_SHOW, SurfaceZoneSet, 0.0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, SurfaceZoneSet, 0.0, FALSE);


		/*
		 *	Now to integrate all the surfaces
		 */

		ss2.str(std::string());
		ss2.clear();
		ss2 << ss1.str() << "Integrating";
// 		TecGUILabelSetText(Lab1Status_LBL_S1, ss2.str().c_str());

// 		EntIndex_t SurfaceZoneCount = SurfaceZoneNumEnd - SurfaceZoneNumBegin + 1;

		TecUtilZoneSetActive(StreamtraceZoneSet, AssignOp_MinusEquals);
		TecUtilZoneSetActive(SliceSet, AssignOp_MinusEquals);

		/*
		 *	Get the volume integration for the polygrid zones.
		 *	These are necessary in order to correctly distribute the
		 *	amount of integrated IntVar to the outer polygrid zones
		 *	(since the integral isn't trustworthy very close to the 
		 *	nucleus)
		 */
		std::vector<ImportType_t> InnerVolumes;
		if (IsOk){
			InnerVolumes.reserve(NumberOfStreamtraces);
			IsOk = MacroIntegrateVarByCellsOverZones(&InnerVolumes,
				SurfaceZoneNumBegin[1],
				SurfaceZoneNumEnd[1],
				VolVarNum,
				XYZVarNums,
				ResultsFileNameSS.str().c_str());
		}
		ImportType_t TotalInnerVolume = 0;
		for (int i = 0; i < InnerVolumes.size() && IsOk; ++i)
			TotalInnerVolume += InnerVolumes[i];

		std::vector<ImportType_t> VolumeRatios;
		if (IsOk)
			VolumeRatios = InnerVolumes;

		if (TotalInnerVolume > 0)
			for (int i = 0; i < VolumeRatios.size() && IsOk; ++i)
				VolumeRatios[i] /= TotalInnerVolume;

		FunctVarNums.push_back(VolVarNum);

		TecUtilDataLoadBegin();

		for (int i = 0; i < ResultsListFull.size() && IsOk; ++i){
			std::vector<ImportType_t> ResultsListSlice[2];
			EntIndex_t VarNum;
			if (FunctVarNums[i] != VolVarNum){
				VarNum = IntVarNum;
				ss.str(std::string());
				ss.clear();
				char* IntVarName = NULL;
				IsOk = TecUtilVarGetName(FunctVarNums[i], &IntVarName);

				if (IsOk && IntVarName != NULL){
					ss << "V" << IntVarNum << "=" << "{" << IntVarName << "}*V" << VolVarNum;
					TecUtilStringDealloc(&IntVarName);
					IsOk = TecUtilDataAlter(ss.str().c_str(), SurfaceZoneSet,
						1, 0, 1,
						1, 0, 1,
						1, 0, 1,
						FieldDataType_Double);
				}
				if (IsOk){
					ArgList_pa StateChangeArgs = TecUtilArgListAlloc();
					TecUtilArgListAppendInt(StateChangeArgs, SV_STATECHANGE, StateChange_VarsAltered);
					TecUtilArgListAppendSet(StateChangeArgs, SV_VARLIST, IntVarSet);
					TecUtilArgListAppendSet(StateChangeArgs, SV_ZONELIST, SurfaceZoneSet);
					TecUtilStateChangedX(StateChangeArgs);
					TecUtilArgListDealloc(&StateChangeArgs);
				}
			}
			else
				VarNum = VolVarNum;

			if (IsOk){
				for (int RunNum = 0; RunNum < 2 && IsOk; ++RunNum){
					ResultsListSlice[RunNum].reserve(NumberOfStreamtraces);
					IsOk = MacroIntegrateVarByCellsOverZones(&ResultsListSlice[RunNum],
						SurfaceZoneNumBegin[RunNum],
						SurfaceZoneNumEnd[RunNum],
						VarNum,
						XYZVarNums,
						ResultsFileNameSS.str().c_str());

					ZoneCount += NumNewZones * 0.1 * 0.5;

					if (IsOk)
						IsOk = SetPercent(ZoneCount, PercentTotal, ss2.str().c_str());
				}
			}

			if (IsOk){
				ImportType_t TotalInnerCharge = 0;
				for (int i = 0; i < ResultsListSlice[1].size() && IsOk; ++i)
					TotalInnerCharge += ResultsListSlice[1][i];
				if (TotalInnerCharge > 0){
					for (int i = 0; i < ResultsListSlice[0].size() && IsOk; ++i)
						ResultsListSlice[0][i] += TotalInnerCharge * VolumeRatios[i];
				}
			}

			if (IsOk)
				ResultsListFull[i].push_back(ResultsListSlice[0]);
		}

		TecUtilDataLoadEnd();



		/*
		 *	All finished with this slice, so delete surface and streamtrace zones if necessary
		 */

		if (IsOk){
			if (!KeepNewZones || (SliceNum != 0 && NumNewZones > MaxNumNewZones)){
				IsOk = TecUtilDataSetDeleteZone(SurfaceZoneSet);
				// 			if (IsOk){
				// 				ArgList_pa StateChangeArgs = TecUtilArgListAlloc();
				// 				TecUtilArgListAppendInt(StateChangeArgs, SV_STATECHANGE, StateChange_ZonesDeleted);
				// 				TecUtilArgListAppendSet(StateChangeArgs, SV_ZONELIST, SurfaceZoneSet);
				// 				TecUtilStateChangedX(StateChangeArgs);
				// 				TecUtilArgListDealloc(&StateChangeArgs);
				// 			}
				if (IsOk)
					IsOk = TecUtilDataSetDeleteZone(StreamtraceZoneSet);
				// 			if (IsOk){
				// 				ArgList_pa StateChangeArgs = TecUtilArgListAlloc();
				// 				TecUtilArgListAppendInt(StateChangeArgs, SV_STATECHANGE, StateChange_ZonesDeleted);
				// 				TecUtilArgListAppendSet(StateChangeArgs, SV_ZONELIST, StreamtraceZoneSet);
				// 				TecUtilStateChangedX(StateChangeArgs);
				// 				TecUtilArgListDealloc(&StateChangeArgs);
				// 			}
// 				if (IsOk)
// 					IsOk = TecUtilDataSetDeleteZone(SliceSet);
				// 			if (IsOk){
				// 				ArgList_pa StateChangeArgs = TecUtilArgListAlloc();
				// 				TecUtilArgListAppendInt(StateChangeArgs, SV_STATECHANGE, StateChange_ZonesDeleted);
				// 				TecUtilArgListAppendSet(StateChangeArgs, SV_ZONELIST, SliceSet);
				// 				TecUtilStateChangedX(StateChangeArgs);
				// 				TecUtilArgListDealloc(&StateChangeArgs);
				// 			}
			}
			else{
				TecUtilZoneSetScatter(SV_SHOW, SurfaceZoneSet, 0.0, FALSE);
				TecUtilZoneSetMesh(SV_SHOW, SurfaceZoneSet, 0.0, TRUE);
				TecUtilZoneSetContour(SV_SHOW, SurfaceZoneSet, 0.0, FALSE);

				TecUtilZoneSetActive(StreamtraceZoneSet, AssignOp_MinusEquals);
				TecUtilZoneSetActive(SurfaceZoneSet, AssignOp_MinusEquals);
				TecUtilZoneSetActive(SliceSet, AssignOp_MinusEquals);
			}
		}

		TecUtilSetDealloc(&StreamtraceZoneSet);
		TecUtilSetDealloc(&SurfaceZoneSet);
		TecUtilSetDealloc(&SliceSet);
	}

	/*
	 *	All slices done, so output data to file
	 */
	if (IsOk)
		IsOk = CSVCreateFileAndInsertData(AddOnFullPathCStr, &ResultsListFull, &FunctVarNums);

	/*
	 *	Finished! Delete 
	 */

// 	if (IsOk){
// 		if (!KeepNewZones){
// 			IsOk = TecUtilDataSetDeleteVar(VarSet);
			// 		if (IsOk){
			// 			ArgList_pa StateChangeArgs = TecUtilArgListAlloc();
			// 			TecUtilArgListAppendInt(StateChangeArgs, SV_STATECHANGE, StateChange_VarsDeleted);
			// 			TecUtilArgListAppendSet(StateChangeArgs, SV_VARLIST, VarSet);
			// 			TecUtilStateChangedX(StateChangeArgs);
			// 			TecUtilArgListDealloc(&StateChangeArgs);
			// 		}
// 		}

	if (ShowOnlyNewZones){
		TecUtilZoneSetActive(AllSlicesSet, AssignOp_Equals);
		TecUtilZoneSetActive(AllStreamtracesSet, AssignOp_PlusEquals);
	}
	else{
		if (HideNewZones)
			TecUtilZoneSetActive(ActiveZones, AssignOp_Equals);
		else
			TecUtilZoneSetActive(ActiveZones, AssignOp_PlusEquals);
	}
// 	}

	TecUtilSetDealloc(&VarSet);
	TecUtilSetDealloc(&IntVarSet);
	TecUtilSetDealloc(&RadiusVarSet);

	TecUtilSetDealloc(&AllStreamtracesSet);
	TecUtilSetDealloc(&AllSlicesSet);
	TecUtilSetDealloc(&AllSurfacesSet);

	TecUtilSetDealloc(&VolZoneSet);
	TecUtilStringDealloc(&AddOnFullPathCStr);

	/*
	 *	Disable 3d view
	 */
	ArgList_pa argList = TecUtilArgListAlloc();
	TecUtilArgListAppendString(argList, SV_P1, SV_THREEDVIEW);
	TecUtilArgListAppendString(argList, SV_P2, SV_DRAWINPERSPECTIVE);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, FALSE);
	TecUtilStyleSetLowLevelX(argList);
	TecUtilArgListDealloc(&argList);

	/*
	 *	Or just set the perspective to a reasonable number (10)
	 */
// 	ArgList_pa argList = TecUtilArgListAlloc();
// 	TecUtilArgListAppendString(argList, SV_P1, SV_THREEDVIEW);
// 	TecUtilArgListAppendString(argList, SV_P2, SV_FIELDOFVIEW);
// 	TecUtilArgListAppendDouble(argList, SV_DVALUE, (double)10);
// 	TecUtilStyleSetLowLevelX(argList);
// 
// 	TecUtilArgListClear(argList);
// 	TecUtilArgListAppendString(argList, SV_P1, SV_THREEDVIEW);
// 	TecUtilArgListAppendString(argList, SV_P2, SV_DRAWINPERSPECTIVE);
// 	TecUtilArgListAppendArbParam(argList, SV_IVALUE, FALSE);
// 	TecUtilStyleSetLowLevelX(argList);
// 	TecUtilArgListDealloc(&argList);

	TecUtilViewDataFit();

	TecUtilStatusSuspend(FALSE);

	TecUtilStatusFinishPercentDone();

	argList = TecUtilArgListAlloc();
	TecUtilArgListAppendString(argList, SV_P1, SV_INTERFACE);
	TecUtilArgListAppendString(argList, SV_P2, SV_AUTOREDRAWISACTIVE);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
	TecUtilStyleSetLowLevelX(argList);

	TecUtilLockFinish(AddOnID);

	if (IsOk)
		TecUtilDialogMessageBox("Finished!", MessageBoxType_Information);
	else
		TecUtilDialogErrMsg("Finished with errors!");

	return IsOk;
}

/*
 *	for every node on a slice, find the distance from the node
 *	to a specified axis using this thing I stole from
 *	http://stackoverflow.com/questions/5227373/minimal-perpendicular-vector-between-a-point-and-a-line
 */
// 	P - point
// 	D - direction of line(unit length)
// 	A - point in line
// 
// 	X - base of the perpendicular line
// 
// 	    P
// 	  / |
// 	 /  |
// 	/   v
// 	A---X----->D
// 
// 	(P - A).D == | X - A |
// 
// 	X == A + ((P - A).D)D
// 	Desired perpendicular : X - P
// 	Distance from P to X is |P - X|
// 	
Boolean_t PopulateRadiusVar(EntIndex_t RadVarNum, 
	EntIndex_t XYZVarNums[3],
	EntIndex_t CutoffVarNum,
	double     CutoffValue,
	EntIndex_t FEZoneNum,
	vec3 A,
	vec3 D)
{
	ZoneType_e SurfZoneType = TecUtilZoneGetType(FEZoneNum);
	Boolean_t IsOk = (SurfZoneType == ZoneType_FEQuad || SurfZoneType == ZoneType_FETriangle);

	LgIndex_t NumNodes = -1, NumElems = -1;
	TecUtilZoneGetInfo(FEZoneNum, &NumNodes, &NumElems, 
		NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	IsOk = (NumNodes > 0 && NumElems > 0);

	FieldData_pa XYZRef[3];
	FieldData_pa RadRef;
	FieldData_pa CutoffRef;

	for (int i = 0; i < 3 && IsOk; ++i){
		XYZRef[i] = TecUtilDataValueGetReadableRef(FEZoneNum, XYZVarNums[i]);
		IsOk = VALID_REF(XYZRef[i]);
	}

	if (IsOk){
		RadRef = TecUtilDataValueGetWritableNativeRef(FEZoneNum, RadVarNum);
		IsOk = VALID_REF(RadRef);
	}
	if (IsOk){
		CutoffRef = TecUtilDataValueGetWritableNativeRef(FEZoneNum, CutoffVarNum);
		IsOk = VALID_REF(CutoffRef);
	}

	FieldValueGetFunction_pf XYZGetFunc[3];
	FieldValueSetFunction_pf RadSetFunc;
	FieldValueGetFunction_pf CutoffGetFunc;

	for (int i = 0; i < 3 && IsOk; ++i){
		XYZGetFunc[i] = TecUtilDataValueRefGetGetFunc(XYZRef[i]);
		IsOk = VALID_REF(XYZGetFunc[i]);
	}

	if (IsOk){
		RadSetFunc = TecUtilDataValueRefGetSetFunc(RadRef);
		IsOk = VALID_REF(RadSetFunc);
	}

	if (IsOk){
		CutoffGetFunc = TecUtilDataValueRefGetGetFunc(CutoffRef);
		IsOk = VALID_REF(CutoffGetFunc);
	}

	CompDir_e Comp = (CompDir_e)TecGUIRadioBoxGetToggle(RBCutoffDi_RADIO_S1);

	for (LgIndex_t NodeNum = 1; NodeNum <= NumNodes && IsOk; ++NodeNum){
		double CutoffCheckValue = CutoffGetFunc(CutoffRef, NodeNum);
		if ((Comp == GreaterThan && CutoffCheckValue <= CutoffValue)
			|| (Comp == LessThan && CutoffCheckValue >= CutoffValue)){
			vec3 P;
			for (int i = 0; i < 3; ++i)
				P[i] = XYZGetFunc[i](XYZRef[i], NodeNum);

			vec3 X = A + (D * (dot((P - A),D)));
			RadSetFunc(RadRef, NodeNum,norm( (X - P)));
		}
		else
			RadSetFunc(RadRef, NodeNum, 0.0);
	}

	return IsOk;
}

Boolean_t SetPercent(unsigned int CurrentNum, unsigned int TotalNum, const char* ProgresssText){
	unsigned int Percent = (int)((double)CurrentNum / (double)TotalNum * 100.);

	std::stringstream ss;
	ss << ProgresssText << "  (" << Percent << "% Complete)";

	TecUtilStatusSuspend(FALSE);
	TecUtilStatusSetPercentDoneText(ss.str().c_str());
	Boolean_t IsOk = TecUtilStatusCheckPercentDone(Percent);
	TecUtilStatusSuspend(TRUE);

	return IsOk;
}

// void SetPercent(unsigned int CurrentNum, unsigned int TotalNum){
// 	unsigned int Percent = (int)((double)CurrentNum / (double)TotalNum * 100.);
// 
// 	std::stringstream ss;
// 	ss << Percent << "% Complete";
// 
// 	TecGUILabelSetText(Lab2Status_LBL_S1, ss.str().c_str());
// 	TecGUIScaleSetValue(SCProgress_SC_S1, Percent);
// }