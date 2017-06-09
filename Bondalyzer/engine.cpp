/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>

#include <omp.h>

#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "TASSERT.h"
#include "GUIDEFS.h"
#include "ENGINE.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "LINEARALGEBRA.h"
#include "ODERUNGEKUTTA.h"
#include "NORMALS.h"
#include "SURFELEMMAP.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
// #include "NORMALS.h"
#include "GRADPATH.h"
#include "ZEROFLUXSURFACE.h"
#include "BUNDLES.h"
#include "MIDPLANEGRADPATH.h"
#include "SPHEREGRADPATH.h"
#include "GEOMTOOLS.h"
#include "SURFTOPOSEG.h"
#include "ISOSURFGRADPATH.h"
#include "ORDEREDPAIRS.h"
#include "BONDBUNDLE.h"
#include "INTEGRATE.h"
#include "CREATEVOLUMEZONE.h"
#include "SETCHEMBONDSTYLE.h"

#include "CRITPOINTGRAPH.h"

using std::string;
using std::vector;



Boolean_t RecordMacroAndJournal(Strand_t StrandToAve,
								EntIndex_t VarToAve)
{
	Boolean_t IsOk = TRUE;

	EntIndex_t NumZones, NumVars;

	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(VarToAve    > 0 && VarToAve    <= NumVars);
	REQUIRE(StrandToAve > 0 && StrandToAve <= TecUtilDataSetGetMaxStrandID());

	if (IsOk)
	{
		char  S[1000];
		char *MacroString;

		sprintf_s(S,             "VarToAve = %d, ", VarToAve);
		sprintf(&S[strlen(S)], "StrandToAve = %d ", StrandToAve);

		S[strlen(S)] = '\0';
		MacroString = TecUtilStringAlloc((int)strlen(S), "MacroString");
		strcpy(MacroString, S);

		if (TecUtilDataSetJournalIsValid())
			IsOk = TecUtilDataSetAddJournalCommand(ADDON_NAME, MacroString, NULL);
		if (TecUtilMacroIsRecordingActive())
			IsOk = TecUtilMacroRecordAddOnCommand(ADDON_NAME, MacroString);
		TecUtilStringDealloc(&MacroString);
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}


/*
 * Utility to recompute the gradients at the boundaries for periodic
 * boundary conditions. Assumes that the X, Y, Z coordinates are the
 * I, J, K coordinates (cartesian, unit spacing).
 *
 * Input:
 *   ZoneNum, UVarNum, VVarNum, WVarNum, ChrgDensVarNum:
 *      Numbers of the zone, gradient variables, and charge-density variable.
 *   PeriodicBC: True if the data is periodic in all three coordinate directions.
 *
 * Return: TRUE if it works, FALSE otherwise.
 *
 */

Boolean_t FixGradForPeriodicBC(EntIndex_t ZoneNum,
							   EntIndex_t UVarNum,
							   EntIndex_t VVarNum,
							   EntIndex_t WVarNum,
							   EntIndex_t ChrgDensVarNum,
							   Boolean_t  PeriodicBC)
{
	Boolean_t  IsOk = TRUE;
	Boolean_t  IsStop = FALSE;
	Boolean_t  ShowStatusBar = TRUE;
	EntIndex_t NumZones, NumVars;
	LgIndex_t  IMax, JMax, KMax;
	LgIndex_t  ii, jj, kk;
	// Density information to be used in gradient calcs
	FieldData_pa ChrgDensFDPtr = NULL;
	EntIndex_t   XVarNum, YVarNum, ZVarNum;
	FieldData_pa XVarFDPtr = NULL;
	FieldData_pa YVarFDPtr = NULL;
	FieldData_pa ZVarFDPtr = NULL;

	REQUIRE(TecUtilDataSetIsAvailable());

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));
	REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
	REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
	REQUIRE(WVarNum > 0 && WVarNum <= NumVars);
	REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
	REQUIRE(VALID_BOOLEAN(PeriodicBC));

	// This is where the lengths of i,j,k are defined based on selected zone
	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */

	XVarNum = TecUtilVarGetNumByAssignment('X'); 
	YVarNum = TecUtilVarGetNumByAssignment('Y'); 
	ZVarNum = TecUtilVarGetNumByAssignment('Z');
	
	REQUIRE(XVarNum > 0 && XVarNum <= NumVars);
	REQUIRE(YVarNum > 0 && YVarNum <= NumVars);
	REQUIRE(ZVarNum > 0 && ZVarNum <= NumVars);

	XVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, XVarNum);
	YVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, YVarNum);
	ZVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, ZVarNum);

	// Giving values to the charge density field data (FD)
	// The FD is a "handle" or "reference" to the data, so it's like a read only pointer so you can pass 
	// its data to functions.
	ChrgDensFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, ChrgDensVarNum);

	// Recompute UVar grad on i-minus and i-plus faces
	if (IsOk)
	{
		LgIndex_t Index1, Index2, IndexMx, IndexMxM1;
		double UVar;
		double dx21, dxMxMxM1, dxtot, rdxtot;
		FieldData_pa UVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, UVarNum);

		// Good trick here; make for loop dependent on sentinal and incrementer value
		for (kk = 1; IsOk && kk <= KMax; kk++)
		{
			for (jj = 1; IsOk && jj <= JMax; jj++)
			{
				Index1 = 1 + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
				// rho1 = TecUtilDataValueGetByRef(ChrgDensFDPtr, Index1);

				Index2 = 2 + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
				// This returns the value from the FD reference at the specified index.
				double rho2 = TecUtilDataValueGetByRef(ChrgDensFDPtr, Index2);

				dx21 = TecUtilDataValueGetByRef(XVarFDPtr, Index2) - TecUtilDataValueGetByRef(XVarFDPtr, Index1);

				IndexMx = IMax + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
				// rhomx = TecUtilDataValueGetByRef(ChrgDensFDPtr, IndexMx);

				IndexMxM1 = IMax - 1 + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
				double rhomxm1 = TecUtilDataValueGetByRef(ChrgDensFDPtr, IndexMxM1);

				dxMxMxM1 = TecUtilDataValueGetByRef(XVarFDPtr, IndexMx) - TecUtilDataValueGetByRef(XVarFDPtr, IndexMxM1);

				dxtot = dx21 + dxMxMxM1;
				rdxtot = 0.0;
				if (ABS(dxtot) > SMALLFLOAT) rdxtot = 1.0/dxtot;

				// Compute and set i-gradient for 1, jj, kk using central diff
				// UVar = rho2 - rhomx;
				UVar = (rho2 - rhomxm1) * rdxtot;   // TODO: remove limitation that X aligned with I 
				TecUtilDataValueSetByRef(UVarFDPtr, Index1, UVar);

				// Compute and set i-gradient for IMax, jj, kk using central diff
				// UVar = rho1 - rhomxm1;
				TecUtilDataValueSetByRef(UVarFDPtr, IndexMx, UVar);
			}
		}
	}

	// Recompute VVar grad on j-minus and j-plus faces
	if (IsOk)
	{
		LgIndex_t Index1, Index2, IndexMx, IndexMxM1;
		double VVar;
		double dy21, dyMxMxM1, dytot, rdytot;
		FieldData_pa VVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, VVarNum);

		for (kk = 1; IsOk && kk <= KMax; kk++)
		{
			for (ii = 1; IsOk && ii <= IMax; ii++)
			{
				Index1 = ii + (kk - 1) * IMax * JMax;
				// rho1 = TecUtilDataValueGetByRef(ChrgDensFDPtr, Index1);

				Index2 = ii + IMax + (kk - 1) * IMax * JMax;
				double rho2 = TecUtilDataValueGetByRef(ChrgDensFDPtr, Index2);

				dy21 = TecUtilDataValueGetByRef(YVarFDPtr, Index2) - TecUtilDataValueGetByRef(YVarFDPtr, Index1);

				IndexMx = ii + (JMax - 1) * IMax + (kk - 1) * IMax * JMax;
				// rhomx = TecUtilDataValueGetByRef(ChrgDensFDPtr, IndexMx);

				IndexMxM1 = ii + (JMax - 2) * IMax + (kk - 1) * IMax * JMax;
				double rhomxm1 = TecUtilDataValueGetByRef(ChrgDensFDPtr, IndexMxM1);

				dyMxMxM1 = TecUtilDataValueGetByRef(YVarFDPtr, IndexMx) - TecUtilDataValueGetByRef(YVarFDPtr, IndexMxM1);

				dytot = dy21 + dyMxMxM1;
				rdytot = 0.0;
				if (ABS(dytot) > SMALLFLOAT) rdytot = 1.0/dytot;

				// Compute and set j-gradient for ii, 1, kk using central diff
				VVar = (rho2 - rhomxm1) * rdytot;    // TODO: remove limitation that Y aligned with J 
				TecUtilDataValueSetByRef(VVarFDPtr, Index1, VVar);

				// Compute and set j-gradient for ii, JMax, kk using central diff
				// VVar = rho1 - rhomxm1;
				TecUtilDataValueSetByRef(VVarFDPtr, IndexMx, VVar);
			}
		}
	}


	// Recompute WVar grad on k-minus and k-plus faces
	if (IsOk)
	{
		LgIndex_t Index1, Index2, IndexMx, IndexMxM1;
		double WVar;
		double dz21, dzMxMxM1, dztot, rdztot;
		FieldData_pa WVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, WVarNum);

		for (jj = 1; IsOk && jj <= JMax; jj++)
		{
			for (ii = 1; IsOk && ii <= IMax; ii++)
			{
				Index1 = ii + (jj - 1) * IMax;
				// rho1 = TecUtilDataValueGetByRef(ChrgDensFDPtr, Index1);

				Index2 = ii + (jj - 1) * IMax + IMax * JMax;
				double rho2 = TecUtilDataValueGetByRef(ChrgDensFDPtr, Index2);

				dz21 = TecUtilDataValueGetByRef(ZVarFDPtr, Index2) - TecUtilDataValueGetByRef(ZVarFDPtr, Index1);

				IndexMx = ii + (jj - 1) * IMax + (KMax - 1) * IMax * JMax;
				// rhomx = TecUtilDataValueGetByRef(ChrgDensFDPtr, IndexMx);

				IndexMxM1 = ii + (jj - 1) * IMax + (KMax - 2) * IMax * JMax;
				double rhomxm1 = TecUtilDataValueGetByRef(ChrgDensFDPtr, IndexMxM1);

				dzMxMxM1 = TecUtilDataValueGetByRef(ZVarFDPtr, IndexMx) - TecUtilDataValueGetByRef(ZVarFDPtr, IndexMxM1);

				dztot = dz21 + dzMxMxM1;
				rdztot = 0.0;
				if (ABS(dztot) > SMALLFLOAT) rdztot = 1.0/dztot;

				// Compute and set k-gradient for ii, jj, 1  using central diff
				WVar = (rho2 - rhomxm1) * rdztot;  // TODO: remove limitation that Y aligned with J 
				TecUtilDataValueSetByRef(WVarFDPtr, Index1, WVar);

				// Compute and set k-gradient for ii, jj, KMax  using central diff
				// WVar = rho1 - rhomxm1;
				TecUtilDataValueSetByRef(WVarFDPtr, IndexMx, WVar);
			}
		}
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}










int PrintTimeDate(std::ofstream & OutFile)
{
	std::time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
	std::string str(buffer);

	if (OutFile.good())
		OutFile << str;

	return 0;
}












Boolean_t ExtractTopology(EntIndex_t        ZoneNum,
						  Set_pa            refinedGridSet,
						  EntIndex_t        UVarNum,
						  EntIndex_t        VVarNum,
						  EntIndex_t        WVarNum,
						  EntIndex_t        GradMagVarNum,
						  EntIndex_t        ChrgDensVarNum,
						  CompletionLevel_e CompletionLevel,
						  Boolean_t         PeriodicBC)
{
	Boolean_t  IsOk = TRUE;
	Boolean_t  IsStop = FALSE;
	Boolean_t  ShowStatusBar = TRUE;
	EntIndex_t NumZones, NumVars;
	LgIndex_t  IMax, JMax, KMax;
	LgIndex_t  ii; //, jj, kk;

	//LgIndex_t NumOfDeletedZones = ClearChemBondZones(refinedGridSet);

	// TIM: All the pink stuff are #defines.
	// Some of them are simple definitions, and some are macros (basically preprocessed functions)(cool trick)
	//	NULL and TECUTILSETNOTMEMBER are set to 0

	// TIM: TecUtil are Tecplot utilities, part of 360, so I can't see all the code there.
	CritPoints_pa CritPoints = NULL;
	Bundles_pa    Bundles    = NULL;
	EntIndex_t    TypeVarNum = TECUTILSETNOTMEMBER;
	// EntIndex_t    GTotVarNum = TECUTILSETNOTMEMBER;
	EntIndex_t    SGradGTXVarNum = TECUTILSETNOTMEMBER;
	EntIndex_t    SGradGTYVarNum = TECUTILSETNOTMEMBER;
	EntIndex_t    SGradGTZVarNum = TECUTILSETNOTMEMBER;
	EntIndex_t    XVarNum = TECUTILSETNOTMEMBER;
	EntIndex_t    YVarNum = TECUTILSETNOTMEMBER;
	EntIndex_t    ZVarNum = TECUTILSETNOTMEMBER;
	FieldData_pa  XVarFDPtr = NULL;
	FieldData_pa  YVarFDPtr = NULL;
	FieldData_pa  ZVarFDPtr = NULL;
	double        DXYZ = 0.0;

	EntIndex_t    CPZoneNum = 0;
	Set_pa        originalZoneSet = NULL;

	OrderedPairs_pa AtomPairsForBonds = NULL;

	//CPG
	CritPointGraph CPGraph;
	//END CPG

	// TEMP
	CHECK(GeomToolsTest());

	TecUtilPleaseWait("Working on topology...", TRUE);

	// TIM: All the REQUIREs are just for debugging, verifying arguement data types/dimensions and return types/dimensions
	REQUIRE(TecUtilDataSetIsAvailable());

	IsOk = TecUtilZoneGetEnabled(&originalZoneSet);

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));
	REQUIRE(!TecUtilSetIsMember(refinedGridSet, ZoneNum));
	REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
	REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
	REQUIRE(WVarNum > 0 && WVarNum <= NumVars);
	REQUIRE(GradMagVarNum > 0 && GradMagVarNum <= NumVars);
	REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
	REQUIRE(VALID_BOOLEAN(PeriodicBC));

	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */

	XVarNum = TecUtilVarGetNumByAssignment('X');
	YVarNum = TecUtilVarGetNumByAssignment('Y'); 
	ZVarNum = TecUtilVarGetNumByAssignment('Z');
	
	REQUIRE(XVarNum > 0 && XVarNum <= NumVars);
	REQUIRE(YVarNum > 0 && YVarNum <= NumVars);
	REQUIRE(ZVarNum > 0 && ZVarNum <= NumVars);

	XVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, XVarNum);
	YVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, YVarNum);
	ZVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, ZVarNum);

	//	TIM:	temporary. run the style macro to set the contour/multicolor
	//			variables to match when isosurfaces are created

	ChemSysView();

	//	TIM: activate source zone
	Set_pa		TempObjectSet = TecUtilSetAlloc(TRUE);

	TecUtilSetAddMember(TempObjectSet, (SetIndex_t)1, TRUE);
	SetValueReturnCode_e SVRC = TecUtilZoneSetActive(TempObjectSet, AssignOp_Equals);
	TecUtilSetDealloc(&TempObjectSet);

	if (!(SVRC == SetValue_Ok && SVRC == SetValueReturnCode_Ok))
	{
		if(SVRC != SetValueReturnCode_DuplicateValue)
		{
			TecUtilDialogMessageBox("Failed to activate source zone",MessageBoxType_Information);
			IsOk = FALSE;
		}
	}

	//	TIM: Open log file and set log boolean.
	// LOG

	std::ofstream BondLog;
	Boolean_t LogBondInfo;

	LogBondInfo = TRUE;

	if (LogBondInfo){

		char* DesktopPath = NULL;
		DesktopPath = getenv("USERPROFILE");

		for (int i = 0 ; i < 100 ; ++i){
			std::stringstream ss;
			ss << DesktopPath << "\\Desktop\\Bondalyzer_Logfile_" << i << ".csv";
			//ss >> FileName;

			BondLog.open(ss.str().c_str(), std::ofstream::out | std::ofstream::app);
			if (BondLog.good())
				break;
			else
				BondLog.clear();
		}
		if (!BondLog.good()){
			std::string TempString = "Failed to initialize log file! Nothing will be logged this time...";
			TecUtilDialogMessageBox(TempString.c_str(),MessageBoxType_Error);
			LogBondInfo = FALSE;
		}
		else{
			PrintTimeDate(BondLog);
			BondLog << ",Bondalyzer initialized...\n\n";
		}
	}


	////	TIM:	for debugging: making sure there's an isosurface on source zone

	//ArgList_pa	CurrentArgList = TecUtilArgListAlloc();

	////	show isosurfaces
	//TecUtilArgListClear(CurrentArgList);
	//TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ISOSURFACELAYERS);
	//TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOW);
	//TecUtilArgListAppendInt(CurrentArgList, SV_OFFSET1, 1);
	//TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	//TecUtilStyleSetLowLevelX(CurrentArgList);

	//TecUtilArgListClear(CurrentArgList);

	////	set isosurface color 1 to variable 1 (rho)
	//TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	//TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_DEFINITIONCONTOURGROUP);
	//TecUtilArgListAppendInt(CurrentArgList, SV_OFFSET1, 1);
	//TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, (SmInteger_t)4);
	//TecUtilStyleSetLowLevelX(CurrentArgList);

	//TecUtilArgListClear(CurrentArgList);

	////	set value to create isosurface at to max bond rho
	//TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	//TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_ISOVALUE1);
	//TecUtilArgListAppendInt(CurrentArgList, SV_OFFSET1, 1);
	//TecUtilArgListAppendDouble(CurrentArgList, SV_DVALUE, (ATOMISOSURFFACTOR * 0.2721));
	//SVRC = TecUtilStyleSetLowLevelX(CurrentArgList);

	//TecUtilArgListClear(CurrentArgList);

	//TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
	//TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWGROUP);
	//TecUtilArgListAppendInt(CurrentArgList, SV_OFFSET1, 1);
	//TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
	//TecUtilStyleSetLowLevelX(CurrentArgList);

	//TecUtilArgListDealloc(&CurrentArgList);

	//TecUtilPleaseWait(NULL, FALSE);

	//return IsOk;


	// TIM: DXYZ is the (1,1,1) point, where 1 is like the (0,0,0) point.
	if (IsOk)
	{
		LgIndex_t IndexIJK2 = 2 + IMax + IMax * JMax;
		DXYZ = TecUtilDataValueGetByRef(XVarFDPtr, IndexIJK2) - TecUtilDataValueGetByRef(XVarFDPtr, 1);
		DXYZ = MIN(DXYZ, TecUtilDataValueGetByRef(YVarFDPtr, IndexIJK2) - TecUtilDataValueGetByRef(YVarFDPtr, 1));
		DXYZ = MIN(DXYZ, TecUtilDataValueGetByRef(ZVarFDPtr, IndexIJK2) - TecUtilDataValueGetByRef(ZVarFDPtr, 1));
	}

 

	/* Set DataSet aux data to completely transfer necessary information to Tecplot */
	if (IsOk)
	{
		AuxData_pa DSAuxDataRef = TecUtilAuxDataDataSetGetRef();
		if (DSAuxDataRef != NULL)
		{
			char valueString[200];
			if (IsOk)
			{
				sprintf_s(valueString, "%d", ZoneNum);
				IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.SourceZone", (ArbParam_t)valueString,
					AuxDataType_String, TRUE);
			}

			if (IsOk)
			{
				char *refinedGridString = NULL;
				Str_BuildSetString(refinedGridSet,
								   FALSE,/*IncludeSquareBrackets*/
								   &refinedGridString);
				IsOk = (refinedGridString!=NULL);
				if (IsOk)
					IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.RefinedGrids", (ArbParam_t)refinedGridString,
													AuxDataType_String, TRUE);
				if (refinedGridString != NULL)
					TecUtilStringDealloc(&refinedGridString);
			}

			if (IsOk)
			{
				sprintf_s(valueString, "%d", ChrgDensVarNum);
				IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.ChrgDensVarNum", (ArbParam_t)valueString,
					AuxDataType_String, TRUE);
			}

			if (IsOk)
			{
				sprintf_s(valueString, "%d", UVarNum);
				IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.GradXVarNum", (ArbParam_t)valueString,
											 AuxDataType_String, TRUE);
			}

			if (IsOk)
			{
				sprintf_s(valueString, "%d", VVarNum);
				IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.GradYVarNum", (ArbParam_t)valueString,
											 AuxDataType_String, TRUE);
			}

			if (IsOk)
			{
				sprintf_s(valueString, "%d", WVarNum);
				IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.GradZVarNum", (ArbParam_t)valueString,
											 AuxDataType_String, TRUE);
			}

			if (IsOk)
			{
				if (PeriodicBC)
					sprintf_s(valueString, "TRUE");
				else
					sprintf_s(valueString, "FALSE");
				IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.PeriodicBC", (ArbParam_t)valueString,
											 AuxDataType_String, TRUE);
			}
		}
		else IsOk = FALSE;
	}

	/*
	 * If it is a periodic structure, recompute the gradients on the boundaries.
	 */
	if (PeriodicBC)
		IsOk = FixGradForPeriodicBC(ZoneNum, UVarNum, VVarNum, WVarNum, ChrgDensVarNum, PeriodicBC);

	/*
	 * Add Variable for Critical Point Type (if it doesn't already exist).
	 * This variable is used in creation of critical points.
	 */
	if (IsOk)
	{
		char VarName[200];
		ArgList_pa ArgList;

		FieldDataType_e *ZoneVarDataType = NULL;

		sprintf_s(VarName, "CritPointType");
		TypeVarNum = TecUtilVarGetNumByName(VarName);
		if (TypeVarNum == TECUTILSETNOTMEMBER)
		{
			ZoneVarDataType = ALLOC_ARRAY(NumZones, FieldDataType_e, "ZoneVarDataType");
			IsOk = (ZoneVarDataType != NULL);

			for (EntIndex_t nz = 0; IsOk && nz < NumZones; nz++)
				ZoneVarDataType[nz] = FieldDataType_Int16;

			TecUtilDataLoadBegin();
			if (IsOk)
			{
				TecUtilLockStart(AddOnID);
				ArgList = TecUtilArgListAlloc();
				TecUtilArgListAppendString(ArgList, SV_NAME, VarName);
				TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)ZoneVarDataType);
				IsOk = TecUtilDataSetAddVarX(ArgList);
				if (IsOk)
				{
					NumVars++;
					TypeVarNum = NumVars;
				}
				TecUtilArgListDealloc(&ArgList);
				TecUtilLockFinish(AddOnID);

				// Inform Tecplot of new var
				if (IsOk)
				{
					Set_pa VSet = TecUtilSetAlloc(TRUE);
					TecUtilSetAddMember(VSet, NumVars, TRUE);
					TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)VSet);
					TecUtilSetDealloc(&VSet);
				}

				/* Set DataSet aux data for TypeVarNum */
				if (IsOk)
				{
					AuxData_pa DSAuxDataRef = TecUtilAuxDataDataSetGetRef();
					if (DSAuxDataRef != NULL)
					{
						char VarNumStr[200];

						sprintf_s(VarNumStr, "%d", TypeVarNum);
						IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.TypeVarNum", (ArbParam_t)VarNumStr,
													 AuxDataType_String, TRUE);

					}
					else IsOk = FALSE;
				}
			}
			TecUtilDataLoadEnd();
		}
		if (ZoneVarDataType != NULL)
			FREE_ARRAY(ZoneVarDataType, "ZoneVarDataType");
	}

	/*
	 * Add Variable for the gradient magnitude (if it doesn't already exist).
	 */
	// TEMP
	// GTotVarNum = 8;

	/*
	 * Add Variable for components of the isosurface (of ChargeDensity) tangent
	 * gradtot gradient components (if they don't already exist).
	 * These variables are used in the computation of far-field topology
	 * components (gradient lines and surfaces) as well as Atom-Cage-? surfaces.
	 */
	if (IsOk)
	{
		AuxData_pa AuxDataRef = TecUtilAuxDataDataSetGetRef();
		if (AuxDataRef != NULL)
		{
			 ArbParam_t    Value;
			 AuxDataType_e Type;
			 Boolean_t     Retain;

			 // Find or compute the X-Component of the gradient of total gradient
			 if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.SGradGTXVarNum",
										  &Value, &Type, &Retain))
			 {
				 char *ValueString = (char *)Value;
				 if (Type == AuxDataType_String)
				 {
					 sscanf_s(ValueString, "%d", &SGradGTXVarNum);
				 }
				 else  // Unsupported aux data type
					 IsOk = FALSE;

				 CHECK(SGradGTXVarNum > 0 && SGradGTXVarNum <= NumVars);
			 }
			 else
			 {
				 char VarName[200];
				 ArgList_pa ArgList;

				 FieldDataType_e *ZoneVarDataType = NULL;

				 sprintf_s(VarName, "GradXGTot");
				 SGradGTXVarNum = TecUtilVarGetNumByName(VarName);
				 if (SGradGTXVarNum == TECUTILSETNOTMEMBER)
				 {
					 ZoneVarDataType = ALLOC_ARRAY(NumZones, FieldDataType_e, "ZoneVarDataType");
					 IsOk = (ZoneVarDataType != NULL);

					 for (EntIndex_t nz = 0; IsOk && nz < NumZones; nz++)
						 ZoneVarDataType[nz] = FieldDataType_Float;

					 TecUtilDataLoadBegin();
					 if (IsOk)
					 {
						 TecUtilLockStart(AddOnID);
						 ArgList = TecUtilArgListAlloc();
						 TecUtilArgListAppendString(ArgList, SV_NAME, VarName);
						 TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)ZoneVarDataType);
						 IsOk = TecUtilDataSetAddVarX(ArgList);
						 if (IsOk)
						 {
							 NumVars++;
							 SGradGTXVarNum = NumVars;
						 }
						 TecUtilArgListDealloc(&ArgList);

						 TecUtilLockFinish(AddOnID);

						 // Inform Tecplot of new var
						 if (IsOk)
						 {
							 Set_pa VSet = TecUtilSetAlloc(TRUE);
							 TecUtilSetAddMember(VSet, SGradGTXVarNum, TRUE);
							 TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)VSet);
							 TecUtilSetDealloc(&VSet);
						 }

						 /* Set DataSet aux data for SGradGT[XYZ]VarNum variables */
						 if (IsOk)
						 {
							 AuxData_pa DSAuxDataRef = TecUtilAuxDataDataSetGetRef(); 
							 if (DSAuxDataRef != NULL)
							 {
								 char VarNumStr[200];

								 sprintf_s(VarNumStr, "%d", SGradGTXVarNum);
								 IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.SGradGTXVarNum", (ArbParam_t)VarNumStr,
															 AuxDataType_String, TRUE);
							 }
							 else IsOk = FALSE;
						 }
					 }
					 if (ZoneVarDataType != NULL)
						 FREE_ARRAY(ZoneVarDataType, "ZoneVarDataType");

					 // For the volume zone (ZoneNum), compute the components of the isosurface-tangent
					 // gradient of the total gradient of charge density.
					 if (IsOk)
					 {
						 double DGTotDX;
						 LgIndex_t ii, jj, kk;
						 LgIndex_t IndexI, IndexIP1, IndexIM1, IndexIMax, IndexIMaxM1, Index1, Index2;

						 FieldData_pa GTotVarFDPtr = NULL;
						 FieldData_pa GTXVarFDPtr = NULL;

						 GTotVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, GradMagVarNum);
						 GTXVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, SGradGTXVarNum);

						 // First compute the X,Y,Z partial derivatives of GTot and store in
						 // GT[XYZ]VarFDPtr
						 //
						 // X-derivatives (X=I direction)
						 for (kk = 1; kk <= KMax; kk++)
						 {
							 for (jj = 1; jj <= JMax; jj++)
							 {
								 for (ii = 2; ii < IMax; ii++)
								 {
									 IndexI   = (ii) + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
									 IndexIP1 = (ii + 1) + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
									 IndexIM1 = (ii - 1) + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
									 
									 DGTotDX = 0.5 * (TecUtilDataValueGetByRef(GTotVarFDPtr, IndexIP1)
													  - TecUtilDataValueGetByRef(GTotVarFDPtr, IndexIM1));
									 TecUtilDataValueSetByRef(GTXVarFDPtr, IndexI, DGTotDX);
								 }
								 // I=1 boundary
								 Index2 = 2 + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
								 Index1 = 1 + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
								 DGTotDX = (TecUtilDataValueGetByRef(GTotVarFDPtr, Index2)
											- TecUtilDataValueGetByRef(GTotVarFDPtr, Index1));
								 TecUtilDataValueSetByRef(GTXVarFDPtr, Index1, DGTotDX);

								 // I=IMax boundary
								 IndexIMax   = (IMax) + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
								 IndexIMaxM1 = (IMax - 1) + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
								 DGTotDX = (TecUtilDataValueGetByRef(GTotVarFDPtr, IndexIMax)
											- TecUtilDataValueGetByRef(GTotVarFDPtr, IndexIMaxM1));
								 TecUtilDataValueSetByRef(GTXVarFDPtr, IndexIMax, DGTotDX);
							 }
						 }
					 }
					 TecUtilDataLoadEnd();
				 }
			 }

			 // Find or compute the Y-Component of the gradient of total gradient
			 if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.SGradGTYVarNum",
										  &Value, &Type, &Retain))
			 {
				 char *ValueString = (char *)Value;
				 if (Type == AuxDataType_String)
				 {
					 sscanf_s(ValueString, "%d", &SGradGTYVarNum);
				 }
				 else  // Unsupported aux data type
					 IsOk = FALSE;

				 CHECK(SGradGTYVarNum > 0 && SGradGTYVarNum <= NumVars);
			 }
			 else
			 {
				 char VarName[200];
				 ArgList_pa ArgList;

				 FieldDataType_e *ZoneVarDataType = NULL;

				// Y-component
				 sprintf_s(VarName, "GradYGTot");
				 SGradGTYVarNum = TecUtilVarGetNumByName(VarName);
				 if (SGradGTYVarNum == TECUTILSETNOTMEMBER)
				 {
					 ZoneVarDataType = ALLOC_ARRAY(NumZones, FieldDataType_e, "ZoneVarDataType");
					 IsOk = (ZoneVarDataType != NULL);

					 for (EntIndex_t nz = 0; IsOk && nz < NumZones; nz++)
						 ZoneVarDataType[nz] = FieldDataType_Float;

					 TecUtilDataLoadBegin();
					 if (IsOk)
					 {
						 TecUtilLockStart(AddOnID);
						 ArgList = TecUtilArgListAlloc();
						 TecUtilArgListAppendString(ArgList, SV_NAME, VarName);
						 TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)ZoneVarDataType);
						 IsOk = TecUtilDataSetAddVarX(ArgList);
						 if (IsOk)
						 {
							 NumVars++;
							 SGradGTYVarNum = NumVars;
						 }
						 TecUtilArgListDealloc(&ArgList);

						 TecUtilLockFinish(AddOnID);

						 // Inform Tecplot of new var
						 if (IsOk)
						 {
							 Set_pa VSet = TecUtilSetAlloc(TRUE);
							 TecUtilSetAddMember(VSet, SGradGTYVarNum, TRUE);
							 TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)VSet);
							 TecUtilSetDealloc(&VSet);
						 }

						 /* Set DataSet aux data for SGradGTYVarNum variable */
						 if (IsOk)
						 {
							 AuxData_pa DSAuxDataRef = TecUtilAuxDataDataSetGetRef();
							 if (DSAuxDataRef != NULL)
							 {
								 char VarNumStr[200];

								 sprintf_s(VarNumStr, "%d", SGradGTYVarNum);
								 IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.SGradGTYVarNum", (ArbParam_t)VarNumStr,
															 AuxDataType_String, TRUE);
							 }
							 else IsOk = FALSE;
						 }
					 }
					 if (ZoneVarDataType != NULL)
						 FREE_ARRAY(ZoneVarDataType, "ZoneVarDataType");

					 // Compute the Y-Gradient of GTot
					if (IsOk)
					{
						double DGTotDY;
						LgIndex_t ii, jj, kk;
						LgIndex_t Index1, Index2;
						LgIndex_t IndexJ, IndexJP1, IndexJM1, IndexJMax, IndexJMaxM1;

						FieldData_pa GTotVarFDPtr = NULL;
						FieldData_pa GTYVarFDPtr = NULL;

						GTotVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, GradMagVarNum);
						GTYVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, SGradGTYVarNum);

						//
						// Y-derivatives (Y=J direction)
						for (kk = 1; kk <= KMax; kk++)
						{
							for (jj = 2; jj < JMax; jj++)
							{
								for (ii = 1; ii <= IMax; ii++)
								{
									IndexJ   = ii + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
									IndexJP1 = ii + (jj) * IMax + (kk - 1) * IMax * JMax;
									IndexJM1 = ii + (jj - 2) * IMax + (kk - 1) * IMax * JMax;
									DGTotDY = 0.5 * (TecUtilDataValueGetByRef(GTotVarFDPtr, IndexJP1)
													 - TecUtilDataValueGetByRef(GTotVarFDPtr, IndexJM1));
									TecUtilDataValueSetByRef(GTYVarFDPtr, IndexJ, DGTotDY);
								}
							}

							// J=1 and J=JMax boundaries
							for (ii = 1; ii <= IMax; ii++)
							{
								Index2 = ii + IMax + (kk - 1) * IMax * JMax;
								Index1 = ii + (kk - 1) * IMax * JMax;
								DGTotDY = (TecUtilDataValueGetByRef(GTotVarFDPtr, Index2)
										   - TecUtilDataValueGetByRef(GTotVarFDPtr, Index1));
								TecUtilDataValueSetByRef(GTYVarFDPtr, Index1, DGTotDY);

								IndexJMax   = ii + (JMax - 1) * IMax + (kk - 1) * IMax * JMax;
								IndexJMaxM1 = ii + (JMax - 2) * IMax + (kk - 1) * IMax * JMax;
								DGTotDY = (TecUtilDataValueGetByRef(GTotVarFDPtr, IndexJMax)
										   - TecUtilDataValueGetByRef(GTotVarFDPtr, IndexJMaxM1));
								TecUtilDataValueSetByRef(GTYVarFDPtr, IndexJMax, DGTotDY);
							}
						}
					 }
					 TecUtilDataLoadEnd();
				 }
			 }


			 // Find or compute the Z-Component of the gradient of total gradient
			 if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.SGradGTZVarNum",
										  &Value, &Type, &Retain))
			 {
				 char *ValueString = (char *)Value;
				 if (Type == AuxDataType_String)
				 {
					 sscanf_s(ValueString, "%d", &SGradGTZVarNum);
				 }
				 else  // Unsupported aux data type
					 IsOk = FALSE;

				 CHECK(SGradGTYVarNum > 0 && SGradGTZVarNum <= NumVars);
			 }
			 else
			 {
				 char VarName[200];
				 ArgList_pa ArgList;

				 FieldDataType_e *ZoneVarDataType = NULL;

				 // Z-component
				 sprintf_s(VarName, "GradZGTot");
				 SGradGTZVarNum = TecUtilVarGetNumByName(VarName);
				 if (SGradGTZVarNum == TECUTILSETNOTMEMBER)
				 {
					 ZoneVarDataType = ALLOC_ARRAY(NumZones, FieldDataType_e, "ZoneVarDataType");
					 IsOk = (ZoneVarDataType != NULL);

					 for (EntIndex_t nz = 0; IsOk && nz < NumZones; nz++)
						 ZoneVarDataType[nz] = FieldDataType_Float;

					 TecUtilDataLoadBegin();
					 if (IsOk)
					 {
						 TecUtilLockStart(AddOnID);
						 ArgList = TecUtilArgListAlloc();
						 TecUtilArgListAppendString(ArgList, SV_NAME, VarName);
						 TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)ZoneVarDataType);
						 IsOk = TecUtilDataSetAddVarX(ArgList);
						 if (IsOk)
						 {
							 NumVars++;
							 SGradGTZVarNum = NumVars;
						 }
						 TecUtilArgListDealloc(&ArgList);

						 TecUtilLockFinish(AddOnID);

						 // Inform Tecplot of new var
						 if (IsOk)
						 {
							 Set_pa VSet = TecUtilSetAlloc(TRUE);
							 TecUtilSetAddMember(VSet, SGradGTZVarNum, TRUE);
							 TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)VSet);
							 TecUtilSetDealloc(&VSet);
						 }

						 /* Set DataSet aux data for SGradGTZVarNum variable */
						 if (IsOk)
						 {
							 AuxData_pa DSAuxDataRef = TecUtilAuxDataDataSetGetRef();
							 if (DSAuxDataRef != NULL)
							 {
								 char VarNumStr[200];

								 sprintf_s(VarNumStr, "%d", SGradGTZVarNum);
								 IsOk = TecUtilAuxDataSetItem(DSAuxDataRef, "CompChem.SGradGTZVarNum", (ArbParam_t)VarNumStr,
															 AuxDataType_String, TRUE);
							 }
							 else IsOk = FALSE;
						 }
					 }
					 if (ZoneVarDataType != NULL)
						 FREE_ARRAY(ZoneVarDataType, "ZoneVarDataType");

					 // Compute the Z-Gradient of GTot
					if (IsOk)
					{
						double DGTotDZ;
						LgIndex_t ii, jj, kk;
						LgIndex_t Index1, Index2;

						LgIndex_t IndexK, IndexKP1, IndexKM1, IndexKMax, IndexKMaxM1;

						FieldData_pa GTotVarFDPtr = NULL;
						FieldData_pa GTZVarFDPtr = NULL;

						GTotVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, GradMagVarNum);
						GTZVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, SGradGTZVarNum);


						//
						// Z-derivatives (Z=K direction)
						for (kk = 2; kk < KMax; kk++)
						{
							for (jj = 1; jj <= JMax; jj++)
							{
								for (ii = 1; ii <= IMax; ii++)
								{
									IndexK   = ii + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
									IndexKP1 = ii + (jj - 1) * IMax + (kk) * IMax * JMax;
									IndexKM1 = ii + (jj - 1) * IMax + (kk - 2) * IMax * JMax;
									DGTotDZ = 0.5 * (TecUtilDataValueGetByRef(GTotVarFDPtr, IndexKP1)
													 - TecUtilDataValueGetByRef(GTotVarFDPtr, IndexKM1));
									TecUtilDataValueSetByRef(GTZVarFDPtr, IndexK, DGTotDZ);
								}
							}
						}

						// K=1 and K=KMax boundaries
						for (jj = 1; jj <= JMax; jj++)
						{
							for (ii = 1; ii <= IMax; ii++)
							{
								Index2 = ii + (jj - 1) * IMax +  IMax * JMax;
								Index1 = ii + (jj - 1) * IMax;
								DGTotDZ = (TecUtilDataValueGetByRef(GTotVarFDPtr, Index2)
										   - TecUtilDataValueGetByRef(GTotVarFDPtr, Index1));
								TecUtilDataValueSetByRef(GTZVarFDPtr, Index1, DGTotDZ);

								IndexKMax   = ii + (jj - 1) * IMax + (KMax - 1) * IMax * JMax;
								IndexKMaxM1 = ii + (jj - 1) * IMax + (KMax - 2) * IMax * JMax;
								DGTotDZ = (TecUtilDataValueGetByRef(GTotVarFDPtr, IndexKMax)
										   - TecUtilDataValueGetByRef(GTotVarFDPtr, IndexKMaxM1));
								TecUtilDataValueSetByRef(GTZVarFDPtr, IndexKMax, DGTotDZ);
							}
						}
					}
				}
			}
		}
	}

	// For the volume zone (ZoneNum), adjust the components of the isosurface-tangent
	// gradient of the total gradient of charge density.
	if (IsOk)
	{
		double Nx, Ny, Nz;
		LgIndex_t ii, jj, kk;
		LgIndex_t IndexI;

		FieldData_pa UVarFDPtr = NULL;
		FieldData_pa VVarFDPtr = NULL;
		FieldData_pa WVarFDPtr = NULL;
		FieldData_pa GTotVarFDPtr = NULL;
		FieldData_pa GTXVarFDPtr = NULL;
		FieldData_pa GTYVarFDPtr = NULL;
		FieldData_pa GTZVarFDPtr = NULL;

		UVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, UVarNum);
		VVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, VVarNum);
		WVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, WVarNum);
		GTotVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, GradMagVarNum);
		GTXVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, SGradGTXVarNum);
		GTYVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, SGradGTYVarNum);
		GTZVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, SGradGTZVarNum);


		// Now subtract the isosurface-normal component of the Grad of total charge
		// density Grad to compute the surface tangent components. Store in
		// GT[XYZ]VarFDPtr
		//
		for (kk = 1; kk <= KMax; kk++)
		{
			for (jj = 1; jj <= JMax; jj++)
			{
				for (ii = 1; ii <= IMax; ii++)
				{
					double RMag, GradGTX, GradGTY, GradGTZ, GradGTDotN;

					IndexI   = (ii) + (jj - 1) * IMax + (kk - 1) * IMax * JMax;

					Nx = TecUtilDataValueGetByRef(UVarFDPtr, IndexI);
					Ny = TecUtilDataValueGetByRef(VVarFDPtr, IndexI);
					Nz = TecUtilDataValueGetByRef(WVarFDPtr, IndexI);
					
					RMag = 1.0 / sqrt(Nx * Nx + Ny * Ny + Nz * Nz + SMALLFLOAT);
					// Normalizing x,y,z density values
					Nx *= RMag;
					Ny *= RMag;
					Nz *= RMag;

					GradGTX = TecUtilDataValueGetByRef(GTXVarFDPtr, IndexI);
					GradGTY = TecUtilDataValueGetByRef(GTYVarFDPtr, IndexI);
					GradGTZ = TecUtilDataValueGetByRef(GTZVarFDPtr, IndexI);
					// Orthonormal gradient at indexI
					GradGTDotN = GradGTX * Nx + GradGTY * Ny + GradGTZ * Nz;

					TecUtilDataValueSetByRef(GTXVarFDPtr, IndexI, (GradGTX - GradGTDotN * Nx));
					TecUtilDataValueSetByRef(GTYVarFDPtr, IndexI, (GradGTY - GradGTDotN * Ny));
					TecUtilDataValueSetByRef(GTZVarFDPtr, IndexI, (GradGTZ - GradGTDotN * Nz));
				}
			}
		}

	}

	

	// Temporarily add variables for X,Y,Z components of surface gradient of total gradient
	/*          if (IsOk && GradGradVarsAdded)
				{
				  FieldDataType_e *DataType = NULL;
				  char *dataset_title = NULL;
				  EntIndex_t NumZones;

				  if (IsOk)
					IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

				  DataType = ALLOC_ARRAY(NumZones, FieldDataType_e, "DataType array");
				  for (EntIndex_t nz = 0; nz < NumZones - 1; nz++)
					DataType[nz] = FieldDataType_Bit;

				  // Make the var's Double for IsoTopoZone
				  DataType[NumZones - 1] = FieldDataType_Double;

				  IsOk = TecUtilDataSetAddVar("GradXGTot", DataType);
				  if (IsOk)
					{
					  NumVars++;
					  GradXGTotVNum = NumVars;
					}

				  if (IsOk) IsOk = TecUtilDataSetAddVar("GradYGTot", DataType);
				  if (IsOk)
					{
					  NumVars++;
					  GradYGTotVNum = NumVars;
					}

				  if (IsOk) IsOk = TecUtilDataSetAddVar("GradZGTot", DataType);
				  if (IsOk)
					{
					  NumVars++;
					  GradZGTotVNum = NumVars;
					}

				  TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)NULL);

				  FREE_ARRAY(DataType, "DataType array");
				  GradGradVarsAdded = FALSE;
				}
	*/

	/* Set Critical Point Type Var to zero for source ZoneNum */
	if (IsOk)
	{
		FieldData_pa TypeVarFDPtr = NULL;
		LgIndex_t    ii;
		LgIndex_t    NumNodex = IMax * JMax * KMax;

		TecUtilDataLoadBegin();
		TypeVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, TypeVarNum);

		for (ii = 0; ii < NumNodex; ii++)
		{
			TecUtilDataValueSetByRef(TypeVarFDPtr, ii + 1, 0);
		}
		TecUtilDataLoadEnd();
	}

	// TEMP - Read critical point data from zone
	// CritPoints = CritPointsGetFromTPZone(ZoneNum+1, ChrgDensVarNum, TypeVarNum,
	//        UVarNum, VVarNum, WVarNum);



	/* Extract the critical points */
	if (IsOk)
	{
		// ZoneVarInfo_pa ZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "ZoneVarInfo");
		ZoneVarInfo_pa ZoneVarInfo = ZoneVarInfoAlloc();
		
		// a->b means the member 'b' of the objected pointed to by 'a'
		ZoneVarInfo->ZoneNum = ZoneNum;
		ZoneVarInfo->UVarNum = UVarNum;
		ZoneVarInfo->VVarNum = VVarNum;
		ZoneVarInfo->WVarNum = WVarNum;
		ZoneVarInfo->ChrgDensVarNum = ChrgDensVarNum;
		ZoneVarInfo->TypeVarNum = TypeVarNum;
		ZoneVarInfo->PeriodicBC = PeriodicBC;

		IsOk = ZoneVarInfoSetMinMax(ZoneVarInfo, refinedGridSet);

		CritPoints = CritPointsAlloc();
		// IsOk = ExtractCriticalPoints(ZoneNum, UVarNum, VVarNum, WVarNum, ChrgDensVarNum,
		//                              TypeVarNum, PeriodicBC, CritPoints);

		// LOG
		if (LogBondInfo){
			PrintTimeDate(BondLog);
			BondLog << ",Beginning critical point search...\n\n";
			PrintTimeDate(BondLog);
			BondLog << ",Critical point information: (IJK are actually XYZ in bohr maybe...)\n";
			BondLog << ",theta and phi (CP directionality) are defined such that where e1 e2 e3 are eigenvalues sorted in decreasing order\n";
			BondLog << ",,for bond (3 -1 CP): theta = atan{sqrt[abs(e1/e2)]},phi = atan{sqrt[abs(e1/e3)]}\n";
			BondLog << ",,for ring (3 1 CP): theta = atan{sqrt[abs(e3/e2)]},phi = atan{sqrt[abs(e3/e1)]}\n";
			BondLog << ",,for cage (3 3 CP): theta = atan{sqrt[abs(e1/e2)]},phi = atan{sqrt[abs(e1/e3)]}\n";
			PrintTimeDate(BondLog);
			BondLog << ",,CP Type,I (in [1 " << IMax << "]),J (in [1 " << JMax << "]),K (in [1 " << KMax << "])"
				<< ",Rho,theta [rad],phi [rad],theta [deg],phi [deg],eig. val. 1,eig. val. 2,eig. val. 3,"
				<< "eig. vec. 1.1,eig. vec. 1.2,eig. vec. 1.3,eig. vec. 2.1,"
				<< "eig. vec. 2.2,eig. vec. 2.3,eig. vec. 3.1,eig. vec. 3.2,eig. vec. 3.3"
				<< ",dRho/dx,dRho/dy,dRho/dz,||dRho||"
				<< ",hess r1c1,hess r1c2/r2c1,hess r2c2,hess r1c3/r3c1,hess r2c3/r3c2,hess r3c3\n";
		}
		ZoneVarInfo->BondLog = &BondLog;
		ZoneVarInfo->LogBondInfo = LogBondInfo;

		IsOk = ExtractCriticalPoints(ZoneVarInfo, NULL, FALSE, CritPoints);

		//CPG
		for (unsigned int CritPointID = 0; CritPointID < CritPoints->NumCrtPts; CritPointID++)
		{
			CPGraph.AddNode(CritPointID);
		}
		//END CPG

		if (LogBondInfo){
			BondLog << "\n";
			PrintTimeDate(BondLog);
			BondLog << ",Critical point search finished...\n\n\n";
			BondLog.close();
		}

		/* Clean up allocated structures/arrays */
		if (ZoneVarInfo != NULL)
			ZoneVarInfoDealloc(&ZoneVarInfo);
			// FREE_ITEM(ZoneVarInfo, "ZoneVarInfo");
	}

	/* Create Critical Point Zone */
	if (IsOk)
	{
		CPZoneNum = CritPointsCreateTPZone(CritPoints, ChrgDensVarNum, TypeVarNum,
										   UVarNum, VVarNum, WVarNum);
		if (CPZoneNum > 0) NumZones++;
	}



	// Allocate ordered-pairs array for saving the atom numbers associated with Bonds
	AtomPairsForBonds = OrderedPairsAlloc();


	/*
	 * To get Bond-Atom lines, seed streamtrace a small step in the
	 * PrincDir away from Bonds
	 */
	CHECK(GradPathTest());

	if (IsOk && CompletionLevel >= CompletionLevel_BondAtomLines)
	{
		GradPath_pa    GradPath = NULL;
		LgIndex_t      BegOffset = CritPointsGetBegOffset(CritPoints, -1);
		LgIndex_t      EndOffset = CritPointsGetEndOffset(CritPoints, -1);
		LgIndex_t      NumCrtPts = CritPointsGetCount(CritPoints);

		// creating new zone and giving it values from CP object
		ZoneVarInfo_pa ZoneVarInfo = ZoneVarInfoAlloc();
		ZoneVarInfo->ZoneNum = ZoneNum;
		ZoneVarInfo->UVarNum = UVarNum;
		ZoneVarInfo->VVarNum = VVarNum;
		ZoneVarInfo->WVarNum = WVarNum;
		ZoneVarInfo->ChrgDensVarNum = ChrgDensVarNum;
		ZoneVarInfo->PeriodicBC = PeriodicBC;

		// Setting dimensions of zone
		IsOk = ZoneVarInfoSetMinMax(ZoneVarInfo, refinedGridSet);

		// Creating new gradient path variable
		GradPath = GradPathAlloc();

		/* Set-up status bar */
		if (ShowStatusBar)
		{
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusStartPercentDone("Finding Bond Lines", TRUE, TRUE);
			TecUtilStatusSuspend(TRUE);
		}

		// Looping over CP array indices from first occurance of bond cp to last occurance of bond
		for (ii = BegOffset; IsOk && ii < EndOffset; ii++)
		{
			double XCrtPt, YCrtPt, ZCrtPt, PrincDirX, PrincDirY, PrincDirZ, dummy;
			char   TypeCrtPt;
			double XPos, YPos, ZPos;

			LgIndex_t EndCrtPtNum1, EndCrtPtNum2;
			
			// Fetches the position and principal directions (eigen vectors) of the ii-th critical point, 
			// for ii=BegOffset or ii=EndOffset this will always be a bond point
			IsOk = CritPointsGetPoint(CritPoints, ii, &XCrtPt, &YCrtPt, &ZCrtPt,
									  &dummy, &TypeCrtPt, &PrincDirX, &PrincDirY, &PrincDirZ);

			// X,Y,Z Pos are the coordinates of the seed point in the principal direction away from CP
			XPos = XCrtPt + DXYZ * PrincDirX;
			YPos = YCrtPt + DXYZ * PrincDirY;
			ZPos = ZCrtPt + DXYZ * PrincDirZ;
			// XPos = XCrtPt + 0.1 * PrincDirX;
			// YPos = YCrtPt + 0.1 * PrincDirY;
			// ZPos = ZCrtPt + 0.1 * PrincDirZ;

			/* Inform Tecplot that major data operation is beginning */
			TecUtilDataLoadBegin();
			TecUtilLockStart(AddOnID);

			/* Integrate to compute path lines */
			if (IsOk)
			{
				double    junk = 0.0;
				LgIndex_t NumPathPoints = 0;

				// Clear the temporary gradient path variable
				GradPathClear(GradPath);

				//	TIM:	GradPathSetPoint is the function to add individual points to gradient paths.
				//			This adds the seed point calculated earlier.
				//			The 0 for point offset means that the seed point is the first point of the grad path.
				//			So it doesn't "start" at the CP, but rather at the seed point.
				/* Seed first point */
				IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk);
				GradPath->BeginCrtPtNum = ii;

				/* Integrate gradient path line */
				IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Forward,
								   CritPoints, 1.0, &NumPathPoints, GradPath);

				/* Add a Tecplot zone for the Atom-Bond connector */
				if (IsOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts)
				{
					char       ZoneName[200];
					LgIndex_t  M1CrtPtNum = ii - BegOffset;
					EntIndex_t ConnectorZoneNum = 0;
					EntIndex_t ConnectorType = -4;  // Atom-Bond

					EndCrtPtNum1 = GradPath->EndCrtPtNum;

					sprintf_s(ZoneName, "Atom_%d_Bond_%d_Zone_%d", GradPath->EndCrtPtNum, M1CrtPtNum, ZoneNum);

					ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, GradPath);

					if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;

					//CPG
					CPGraph.AddEdge(ii, EndCrtPtNum1);
					//END CPG
				}
			}

			// IsOk = TecUtilStreamtraceAdd (1, Streamtrace_VolumeLine, StreamDir_Forward,
			//                               XPos, YPos, ZPos, 0.0, 0.0, 0.0);
			XPos = XCrtPt - DXYZ * PrincDirX;
			YPos = YCrtPt - DXYZ * PrincDirY;
			ZPos = ZCrtPt - DXYZ * PrincDirZ;
			// XPos = XCrtPt - 0.1 * PrincDirX;
			// YPos = YCrtPt - 0.1 * PrincDirY;
			// ZPos = ZCrtPt - 0.1 * PrincDirZ;
			// XPos = XCrtPt - 0.5 * PrincDirX;
			// YPos = YCrtPt - 0.5 * PrincDirY;
			// ZPos = ZCrtPt - 0.5 * PrincDirZ;

			/* Integrate to compute path lines */
			if (IsOk)
			{
				double    junk = 0.0;
				LgIndex_t NumPathPoints = 0;

				GradPathClear(GradPath);

				/* Seed first point */
				IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk);
				GradPath->BeginCrtPtNum = ii;

				/* Integrate gradient path line */
				IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Forward,
								   CritPoints, 1.0, &NumPathPoints, GradPath);

				/* Add a Tecplot zone for the Atom-Bond connector */
				if (IsOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts)
				{
					char       ZoneName[200];
					LgIndex_t  M1CrtPtNum = ii - BegOffset;
					EntIndex_t ConnectorZoneNum = 0;
					EntIndex_t ConnectorType = -4;  // Atom-Bond

					EndCrtPtNum2 = GradPath->EndCrtPtNum;

					sprintf_s(ZoneName, "Atom_%d_Bond_%d_Zone_%d", GradPath->EndCrtPtNum, M1CrtPtNum, ZoneNum);

					ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, GradPath);

					if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;

					//CPG
					CPGraph.AddEdge(ii, EndCrtPtNum2);
					//END CPG
				}
			}

			// Set the EndCrtPtNum pairs for the bond
			IsOk = OrderedPairsSetAtOffset(AtomPairsForBonds, ii - BegOffset, EndCrtPtNum1, EndCrtPtNum2);


			/* Update status bar */
			if (ShowStatusBar)
			{
				char PercentDoneText[200];
				int  PercentDone;
				PercentDone = (int)((100 * (ii - BegOffset)) / (EndOffset - BegOffset));
				sprintf_s(PercentDoneText, "Finding Bond Lines: %d Percent Done", PercentDone);
				TecUtilStatusSuspend(FALSE);
				TecUtilStatusSetPercentDoneText(PercentDoneText);
				if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
				TecUtilStatusSuspend(TRUE);
				if (IsStop) IsOk = FALSE;
			}

			/* Inform Tecplot that major data operation is ending */
			TecUtilLockFinish(AddOnID);
			TecUtilDataLoadEnd();
		}

		/* Clean up allocated structures/arrays */
		if (ZoneVarInfo != NULL)
			ZoneVarInfoDealloc(&ZoneVarInfo);
			// FREE_ITEM(ZoneVarInfo, "ZoneVarInfo");
		if (GradPath != NULL)
			GradPathDealloc(&GradPath);
		if (ShowStatusBar)
		{
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusFinishPercentDone();
			TecUtilStatusSuspend(TRUE);
		}
	}


	/*
	 * To get Ring-Cage lines, seed streamtrace a small step in the
	 * PrincDir away from Ring critical points.
	 */
	// Here principal direction is the direction in which the ring CP is a maximum, that way
	// as you move along it you'll decrease in electron density and head towards a cage.

	if (IsOk && CompletionLevel >= CompletionLevel_RingCageLines)
	{
		GradPath_pa    GradPath = NULL;
		LgIndex_t      BegOffset = CritPointsGetBegOffset(CritPoints, 1);
		LgIndex_t      EndOffset = CritPointsGetEndOffset(CritPoints, 1);
		LgIndex_t      NumCrtPts = CritPointsGetCount(CritPoints);

		// ZoneVarInfo_pa ZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "ZoneVarInfo");
		ZoneVarInfo_pa ZoneVarInfo = ZoneVarInfoAlloc();
		ZoneVarInfo->ZoneNum = ZoneNum;
		ZoneVarInfo->UVarNum = UVarNum;
		ZoneVarInfo->VVarNum = VVarNum;
		ZoneVarInfo->WVarNum = WVarNum;
		ZoneVarInfo->ChrgDensVarNum = ChrgDensVarNum;
		ZoneVarInfo->PeriodicBC = PeriodicBC;

		IsOk = ZoneVarInfoSetMinMax(ZoneVarInfo, refinedGridSet);

		GradPath = GradPathAlloc();

		/* Set-up status bar */
		if (ShowStatusBar)
		{
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusStartPercentDone("Finding Ring-Cage Lines", TRUE, TRUE);
			TecUtilStatusSuspend(TRUE);
		}

		for (ii = BegOffset; IsOk && ii < EndOffset; ii++)
		{
			double XCrtPt, YCrtPt, ZCrtPt, PrincDirX, PrincDirY, PrincDirZ, dummy;
			char   TypeCrtPt;
			double XPos, YPos, ZPos;

			IsOk = CritPointsGetPoint(CritPoints, ii, &XCrtPt, &YCrtPt, &ZCrtPt,
									  &dummy, &TypeCrtPt, &PrincDirX, &PrincDirY, &PrincDirZ);

			// XPos = XCrtPt + 0.1 * PrincDirX;
			// YPos = YCrtPt + 0.1 * PrincDirY;
			// ZPos = ZCrtPt + 0.1 * PrincDirZ;
			// XPos = XCrtPt + 0.5 * PrincDirX;
			// YPos = YCrtPt + 0.5 * PrincDirY;
			// ZPos = ZCrtPt + 0.5 * PrincDirZ;
			// XPos = XCrtPt + 1.0 * PrincDirX;
			// YPos = YCrtPt + 1.0 * PrincDirY;
			// ZPos = ZCrtPt + 1.0 * PrincDirZ;
			XPos = XCrtPt + DXYZ * PrincDirX;
			YPos = YCrtPt + DXYZ * PrincDirY;
			ZPos = ZCrtPt + DXYZ * PrincDirZ;

			/* Inform Tecplot that major data operation is beginning */
			TecUtilDataLoadBegin();
			TecUtilLockStart(AddOnID);

			/* Integrate to compute path lines */
			if (IsOk)
			{
				double    junk = 0.0;
				LgIndex_t NumPathPoints = 0;

				GradPathClear(GradPath);

				/* Seed first point */
				IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk);  // correct ChrgDens set later
				GradPath->BeginCrtPtNum = ii;

				/* Integrate gradient path line */
				IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Reverse,
								   CritPoints, 1.0, &NumPathPoints, GradPath);

				/* Add a Tecplot zone for the Ring-Cage connector */
				// if (IsOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts)
				if (IsOk && NumPathPoints > 1 && GradPath->EndCrtPtNum >= -1 && GradPath->EndCrtPtNum < NumCrtPts)
				{
					char       ZoneName[200];
					LgIndex_t  P1CrtPtNum = ii - BegOffset;
					EntIndex_t ConnectorZoneNum = 0;
					EntIndex_t ConnectorType = 4;  // Ring-Cage

					if (GradPath->EndCrtPtNum >= 0)
						sprintf_s(ZoneName, "Ring_%d_Cage_%d_Zone_%d", P1CrtPtNum, GradPath->EndCrtPtNum, ZoneNum);
					else
						sprintf_s(ZoneName, "Ring_%d_Cage_FF_Zone_%d", P1CrtPtNum, ZoneNum);

					ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, GradPath);
					// TODO: Store the ConnectorZoneNum in an array list

					if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;

					//CPG
					if (P1CrtPtNum >= 0 && GradPath->EndCrtPtNum >= 0)
						CPGraph.AddEdge(ii, GradPath->EndCrtPtNum);
					//END CPG
				}
			}

			// IsOk = TecUtilStreamtraceAdd (1, Streamtrace_VolumeLine, StreamDir_Reverse,
			//                               XPos, YPos, ZPos, 0.0, 0.0, 0.0);
			// XPos = XCrtPt - 0.1 * PrincDirX;
			// YPos = YCrtPt - 0.1 * PrincDirY;
			// ZPos = ZCrtPt - 0.1 * PrincDirZ;
			// XPos = XCrtPt - 0.5 * PrincDirX;
			// YPos = YCrtPt - 0.5 * PrincDirY;
			// ZPos = ZCrtPt - 0.5 * PrincDirZ;
			// XPos = XCrtPt - 1.0 * PrincDirX;
			// YPos = YCrtPt - 1.0 * PrincDirY;
			// ZPos = ZCrtPt - 1.0 * PrincDirZ;
			XPos = XCrtPt - DXYZ * PrincDirX;
			YPos = YCrtPt - DXYZ * PrincDirY;
			ZPos = ZCrtPt - DXYZ * PrincDirZ;

			/* Integrate to compute path lines */
			if (IsOk)
			{
				double    junk = 0.0;
				LgIndex_t NumPathPoints = 0;

				GradPathClear(GradPath);

				/* Seed first point */
				IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk); // correct ChrgDens set later
				GradPath->BeginCrtPtNum = ii;

				/* Integrate gradient path line */
				IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Reverse,
								   CritPoints, 1.0, &NumPathPoints, GradPath);

				/* Add a Tecplot zone for the Ring-Cage connector */
				// if (IsOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts)
				if (IsOk && NumPathPoints > 1 && GradPath->EndCrtPtNum >= -1 && GradPath->EndCrtPtNum < NumCrtPts)
				{
					char       ZoneName[200];
					LgIndex_t  P1CrtPtNum = ii - BegOffset;
					EntIndex_t ConnectorZoneNum = 0;
					EntIndex_t ConnectorType = 4;  // Atom-Bond

					if (GradPath->EndCrtPtNum >= 0)
						sprintf_s(ZoneName, "Ring_%d_Cage_%d_Zone_%d", P1CrtPtNum, GradPath->EndCrtPtNum, ZoneNum);
					else
						sprintf_s(ZoneName, "Ring_%d_Cage_FF_Zone_%d", P1CrtPtNum, ZoneNum);

					ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, GradPath);
					// TODO: Store the ConnectorZoneNum in an array list

					if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;

					//CPG
					if (P1CrtPtNum >= 0 && GradPath->EndCrtPtNum >= 0)
						CPGraph.AddEdge(ii, GradPath->EndCrtPtNum);
					//END CPG
				}
			}

			// IsOk = TecUtilStreamtraceAdd (1, Streamtrace_VolumeLine, StreamDir_Reverse,
			//                               XPos, YPos, ZPos, 0.0, 0.0, 0.0);

			/* Update status bar */
			if (ShowStatusBar)
			{
				char PercentDoneText[200];
				int  PercentDone;
				PercentDone = (int)((100 * (ii - BegOffset)) / (EndOffset - BegOffset));
				sprintf_s(PercentDoneText, "Finding Ring-Cage Lines: %d Percent Done", PercentDone);
				TecUtilStatusSuspend(FALSE);
				TecUtilStatusSetPercentDoneText(PercentDoneText);
				if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
				TecUtilStatusSuspend(TRUE);
				if (IsStop) IsOk = FALSE;
			}

			/* Inform Tecplot that major data operation is ending */
			TecUtilLockFinish(AddOnID);
			TecUtilDataLoadEnd();
		}

		/* Clean up allocated structures/arrays */
		if (ZoneVarInfo != NULL)
			ZoneVarInfoDealloc(&ZoneVarInfo);
			// FREE_ITEM(ZoneVarInfo, "ZoneVarInfo");
		if (GradPath != NULL)
			GradPathDealloc(&GradPath);

		if (ShowStatusBar)
		{
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusFinishPercentDone();
			TecUtilStatusSuspend(TRUE);
		}
	}

	/*
	 * To get Bond-Cage lines and Bond-Cage-Ring surface, seed streamtrace
	 * around a small circle surrounding a bond critical point. The circle
	 * must be perpendicular to principle directions.
	 *
	 * Note: This is also the point were the irreducible bundles are defined.
	 */
	// if (0)
	if (IsOk && CompletionLevel >= CompletionLevel_BRCSurfaces)
	{
		CircleGradPath_pa    CircleGradPath = NULL;
		LgIndex_t            BegOffset = CritPointsGetBegOffset(CritPoints, -1);
		LgIndex_t            EndOffset = CritPointsGetEndOffset(CritPoints, -1);
		LgIndex_t            NumCrtPts = CritPointsGetCount(CritPoints);
		LgIndex_t            P1BegOffset = CritPointsGetBegOffset(CritPoints, 1);
		LgIndex_t            P3BegOffset = CritPointsGetBegOffset(CritPoints, 3);


		// ZoneVarInfo_pa ZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "ZoneVarInfo");
		ZoneVarInfo_pa ZoneVarInfo = ZoneVarInfoAlloc();
		ZoneVarInfo->ZoneNum = ZoneNum;
		ZoneVarInfo->UVarNum = UVarNum;
		ZoneVarInfo->VVarNum = VVarNum;
		ZoneVarInfo->WVarNum = WVarNum;
		ZoneVarInfo->ChrgDensVarNum = ChrgDensVarNum;
		ZoneVarInfo->PeriodicBC = PeriodicBC;

		IsOk = ZoneVarInfoSetMinMax(ZoneVarInfo, refinedGridSet);

		CircleGradPath = CircleGradPathAlloc();

		/* Set-up status bar */
		if (ShowStatusBar)
		{
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusStartPercentDone("Finding Bond-Ring-Cage Surfaces", TRUE, TRUE);
			TecUtilStatusSuspend(TRUE);
		}

		for (ii = BegOffset; IsOk && ii < EndOffset; ii++)
		{
			/* Inform Tecplot that major data operation is beginning */
			TecUtilDataLoadBegin();
			TecUtilLockStart(AddOnID);

			/* Integrate to compute path lines */
			if (IsOk)
			{
				LgIndex_t nc;
				LgIndex_t NumConnectors = 0;

				CircleGradPathClear(CircleGradPath);

				if (ii == 20)
				{
					IsOk = IsOk;
				}

				/* Integrate gradient path lines from seeded circle */
				// IsOk = CircleGradPathAdd(ZoneVarInfo, ii, CritPoints,
				//                          StreamDir_Reverse, 1.0, CircleGradPath);
				IsOk = CircleGradPathAdd(ZoneVarInfo, ii, CritPoints,
										 StreamDir_Reverse, DXYZ, DXYZ, CircleGradPath);

				/* Decimate large GradPaths to conserve memory */
				// IsOk = CircleGradPathDecimate(CircleGradPath, MaxSize);

				/* Extract the Bond-Ring and Bond-Cage lines from the CircleGradPath */
				NumConnectors = CircleGradPathGetCPNCount(CircleGradPath);

				for (nc = 0; nc < NumConnectors; nc++)
				{
					LgIndex_t ConnectPathNum = CircleGradPathGetCPN(CircleGradPath, nc);

					GradPath_pa ConnectPath  = CircleGradPathGetGP(CircleGradPath, ConnectPathNum);
					LgIndex_t   EndCrtPtNum  = ConnectPath->EndCrtPtNum;
					LgIndex_t   BegCrtPtNum  = ConnectPath->BeginCrtPtNum;
					char        EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
					char        BegCrtPtType = CritPointsGetType(CritPoints, BegCrtPtNum);
					LgIndex_t   NumPathPoints = GradPathGetCount(ConnectPath);

					char        ConnectorType = EndCrtPtType + BegCrtPtType;
					if (ConnectorType == 0 && EndCrtPtType ==  1) ConnectorType =  1;
					if (ConnectorType == 0 && EndCrtPtType == -1) ConnectorType = -1;

					/* Add a Tecplot zone for the Bond-Ring or Bond-Cage connector */
					if (IsOk && NumPathPoints > 0 && EndCrtPtNum >= 0 &&
						EndCrtPtNum < CritPointsGetCount(CritPoints))
					{
						char       ZoneName[200];
						LgIndex_t  M1CrtPtNum = ii - BegOffset;
						EntIndex_t ConnectorZoneNum = 0;
						LgIndex_t  P1CrtPtNum = EndCrtPtNum - P1BegOffset;
						LgIndex_t  P3CrtPtNum = EndCrtPtNum - P3BegOffset;

						if (ConnectorType == 1)
							sprintf_s(ZoneName, "Bond_%d_Ring_%d_Zone_%d", M1CrtPtNum, P1CrtPtNum, ZoneNum);
						else if (ConnectorType == 2)
							sprintf_s(ZoneName, "Bond_%d_Cage_%d_Zone_%d", M1CrtPtNum, P3CrtPtNum, ZoneNum);

						// ConnectorZoneNum = CreateConnectorZone(ZoneName, ChrgDensVarNum, ConnectorType, TypeVarNum,
						//   NumPathPoints, ConnectPath);
						ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, ConnectPath);

						if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;

						//CPG
						CPGraph.AddEdge(ii, EndCrtPtNum);
						//END CPG
					}
				}

				/* The edges of the surface are made up of Bond-Cage, Bond-Ring, and
				 * Ring-Cage lines. To represent this "triangle" as an ordered grid
				 * we must append the Bond-Ring and Ring-Cage lines to create a
				 * single Bond-Ring-Cage line that represents either the J=1 or
				 * J=JMax line of the IJ-ordered Tecplot zone.
				 */
				/* In preparation for the creation of the bond-ring-cage surface,
				 * resample the gradient path lines in CircleGradPath to make them
				 * all have an equal number of points.
				 */
				// IsOk = CircleGradPathResample(CircleGradPath, brSurfaceIMax);

				NumConnectors = CircleGradPathGetCPNCount(CircleGradPath);

				if (NumConnectors > 0)
				{
					for (nc = 0; nc < NumConnectors; nc++)
						// for (nc=0; nc<2; nc++)
					{
						LgIndex_t  CGPNumBeg, CGPNumEnd;
						LgIndex_t  ncp1 = nc + 1;
						EntIndex_t SurfaceZoneNum = 0;

						if (nc == NumConnectors - 1) ncp1 = ncp1 - NumConnectors;

						CGPNumBeg = CircleGradPathGetCPN(CircleGradPath, nc);
						CGPNumEnd = CircleGradPathGetCPN(CircleGradPath, ncp1);
						if (CGPNumEnd <= CGPNumBeg) CGPNumEnd = CGPNumEnd + CircleGradPathGetCount(CircleGradPath);

						SurfaceZoneNum = CreateBRSurfaceZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints,
															 CircleGradPath, brSurfaceIMax, CGPNumBeg, CGPNumEnd);

						/* Create two bond bundles for this Bond-Ring-Cage combination */
						if (Bundles == NULL) Bundles = BundlesAlloc();
						IsOk =  CreateBundlesForBond(ZoneNum, CritPoints, CircleGradPath, CGPNumBeg, CGPNumEnd, Bundles);
					}
				}
				/* NumConnectors== 0, Create one big surface zone from the CircleGradPath */
				else
				{
					LgIndex_t  CGPNumBeg = 0;
					LgIndex_t  CGPNumEnd = CircleGradPathGetCount(CircleGradPath);
					EntIndex_t SurfaceZoneNum = 0;

					SurfaceZoneNum = CreateBRSurfaceZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints,
														 CircleGradPath, brSurfaceIMax, CGPNumBeg, CGPNumEnd);
				}


			}

			/* Update status bar */
			if (ShowStatusBar)
			{
				char PercentDoneText[200];
				int  PercentDone = (int)((100 * (ii - BegOffset)) / (EndOffset - BegOffset));
				sprintf_s(PercentDoneText, "Finding Bond-Ring-Cage Surfaces: %d Percent Done", PercentDone);
				TecUtilStatusSuspend(FALSE);
				TecUtilStatusSetPercentDoneText(PercentDoneText);
				if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
				TecUtilStatusSuspend(TRUE);
				if (IsStop) IsOk = FALSE;
			}

			/* Inform Tecplot that major data operation is ending */
			TecUtilLockFinish(AddOnID);
			TecUtilDataLoadEnd();
		}
		if (ZoneVarInfo != NULL)
			ZoneVarInfoDealloc(&ZoneVarInfo);
			// FREE_ITEM(ZoneVarInfo, "ZoneVarInfo");
		if (CircleGradPath != NULL)
			CircleGradPathDealloc(&CircleGradPath);
		if (ShowStatusBar)
		{
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusFinishPercentDone();
			TecUtilStatusSuspend(TRUE);
		}
	}


	/*
	 * To get Ring-Atom lines and Bond-Atom-Ring surface, seed streamtrace
	 * around a small circle surrounding a ring critical point. The circle
	 * must be perpendicular to principle directions.
	 */
	// if (0)
	if (IsOk && CompletionLevel >= CompletionLevel_RBASurfaces)
	{
		CircleGradPath_pa    CircleGradPath = NULL;
		LgIndex_t            BegOffset = CritPointsGetBegOffset(CritPoints, 1);
		LgIndex_t            EndOffset = CritPointsGetEndOffset(CritPoints, 1);
		LgIndex_t            NumCrtPts = CritPointsGetCount(CritPoints);
		LgIndex_t            M1BegOffset = CritPointsGetBegOffset(CritPoints, -1);
		LgIndex_t            M3BegOffset = CritPointsGetBegOffset(CritPoints, -3);


		// ZoneVarInfo_pa ZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "ZoneVarInfo");
		ZoneVarInfo_pa ZoneVarInfo = ZoneVarInfoAlloc();
		ZoneVarInfo->ZoneNum = ZoneNum;
		ZoneVarInfo->UVarNum = UVarNum;
		ZoneVarInfo->VVarNum = VVarNum;
		ZoneVarInfo->WVarNum = WVarNum;
		ZoneVarInfo->ChrgDensVarNum = ChrgDensVarNum;
		ZoneVarInfo->PeriodicBC = PeriodicBC;

		IsOk = ZoneVarInfoSetMinMax(ZoneVarInfo, refinedGridSet);

		CircleGradPath = CircleGradPathAlloc();

		/* Set-up status bar */
		if (ShowStatusBar)
		{
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusStartPercentDone("Finding Ring-Atom-Bond Surfaces", TRUE, TRUE);
			TecUtilStatusSuspend(TRUE);
		}

	   for (ii = BegOffset; IsOk && ii < EndOffset; ii++)
		{
			/* Inform Tecplot that major data operation is beginning */
			TecUtilDataLoadBegin();
			TecUtilLockStart(AddOnID);

			/* Integrate to compute path lines */
			if (IsOk)
			{
				CircleGradPathClear(CircleGradPath);

				/* Integrate gradient path lines from seeded circle */
				// IsOk = CircleGradPathAdd(ZoneVarInfo, ii, CritPoints,
				//                          StreamDir_Forward, 1.0, CircleGradPath);
				IsOk = CircleGradPathAdd(ZoneVarInfo, ii, CritPoints,
										 StreamDir_Forward, DXYZ, 2.0*DXYZ, CircleGradPath);

				if (IsOk == FALSE)
				{
					char Message[200];
					sprintf_s(Message, "Error computing CircleGradPath number %d", ii);
					TecUtilDialogMessageBox(Message, MessageBox_Warning);
				}
			}

			/* Decimate large GradPaths to conserve memory */
			// IsOk = CircleGradPathDecimate(CircleGradPath, MaxSize);

			if (IsOk)
			{
				LgIndex_t nc;
				LgIndex_t NumConnectors = 0;

				/* Extract the Ring-Atom lines from the CircleGradPath */
				NumConnectors = CircleGradPathGetCPNCount(CircleGradPath);

				for (nc = 0; nc < NumConnectors; nc++)
				{
					LgIndex_t ConnectPathNum = CircleGradPathGetCPN(CircleGradPath, nc);

					GradPath_pa ConnectPath  = CircleGradPathGetGP(CircleGradPath, ConnectPathNum);
					LgIndex_t   EndCrtPtNum  = ConnectPath->EndCrtPtNum;
					LgIndex_t   BegCrtPtNum  = ConnectPath->BeginCrtPtNum;
					if (EndCrtPtNum >= 0){
						char        EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
						char        BegCrtPtType = CritPointsGetType(CritPoints, BegCrtPtNum);
						LgIndex_t   NumPathPoints = GradPathGetCount(ConnectPath);

						char        ConnectorType = EndCrtPtType + BegCrtPtType;
						if (ConnectorType == 0 && EndCrtPtType == 1) ConnectorType = 1;
						if (ConnectorType == 0 && EndCrtPtType == -1) ConnectorType = -1;


						/* Add a Tecplot zone for the Ring-Atom connector */
						if (IsOk && NumPathPoints > 0 && EndCrtPtNum >= 0 &&
							EndCrtPtNum < CritPointsGetCount(CritPoints))
						{
							char       ZoneName[200];
							LgIndex_t  P1CrtPtNum = ii - BegOffset;
							EntIndex_t ConnectorZoneNum = 0;
							LgIndex_t  M1CrtPtNum = EndCrtPtNum - M1BegOffset;
							LgIndex_t  M3CrtPtNum = EndCrtPtNum - M3BegOffset;

							if (ConnectorType == -1)
								sprintf_s(ZoneName, "Ring_%d_Bond_%d_Zone_%d", P1CrtPtNum, M1CrtPtNum, ZoneNum);
							else if (ConnectorType == -2)
								sprintf_s(ZoneName, "Ring_%d_Atom_%d_Zone_%d", P1CrtPtNum, M3CrtPtNum, ZoneNum);

							/* Only do the Ring-Atom lines, Ring-Bond was already done. */
							if (ConnectorType == -2)
							{
								ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, ConnectPath);

								if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;

								//CPG
								CPGraph.AddEdge(ii, EndCrtPtNum);
								//END CPG
							}
						}
					}
				}

				/* The edges of the surface are made up of Ring-Atom, Ring-Bond, and
				 * Bond-Atom lines. To represent this "triangle" as an ordered grid
				 * we must append the Ring-Bond and Bond-Atom lines to create a
				 * single Ring-Bond-Atom line that represents either the J=1 or
				 * J=JMax line of the IJ-ordered Tecplot zone.
				 */
				/* In preparation for the creation of the Ring-Bond-Atom surface,
				 * resample the gradient path lines in CircleGradPath to make them
				 * all have an equal number of points.
				 */
				// IsOk = CircleGradPathResample(CircleGradPath, brSurfaceIMax);

				if (IsOk)
				{
					LgIndex_t NCMax;
					NumConnectors = CircleGradPathGetCPNCount(CircleGradPath);

					if (NumConnectors > 0)
					{
						NCMax = NumConnectors;
						if (!(CircleGradPath->CompleteCircle)) NCMax--;

						for (nc = 0; nc < NCMax; nc++)
							// for (nc=0; nc<2; nc++)
						{
							LgIndex_t  CGPNumBeg, CGPNumEnd;
							LgIndex_t  ncp1 = nc + 1;
							EntIndex_t SurfaceZoneNum = 0;

							if (nc == NumConnectors - 1) ncp1 = ncp1 - NumConnectors;

							CGPNumBeg = CircleGradPathGetCPN(CircleGradPath, nc);
							CGPNumEnd = CircleGradPathGetCPN(CircleGradPath, ncp1);
							if (CGPNumEnd <= CGPNumBeg) CGPNumEnd = CGPNumEnd + CircleGradPathGetCount(CircleGradPath);

							SurfaceZoneNum = CreateBRSurfaceZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints,
																 CircleGradPath, brSurfaceIMax, CGPNumBeg, CGPNumEnd);
						}
					}
					/* NumConnectors== 0, Create one big surface zone from the CircleGradPath */
					else
					{
						LgIndex_t  CGPNumBeg = 0;
						LgIndex_t  CGPNumEnd = CircleGradPathGetCount(CircleGradPath);
						EntIndex_t SurfaceZoneNum = 0;

						SurfaceZoneNum = CreateBRSurfaceZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints,
															 CircleGradPath, brSurfaceIMax, CGPNumBeg, CGPNumEnd);
					}
				}
			}

			/* Update status bar */
			if (ShowStatusBar)
			{
				char PercentDoneText[200];
				int  PercentDone;
				PercentDone = (int)((100 * (ii - BegOffset)) / (EndOffset - BegOffset));
				sprintf_s(PercentDoneText, "Finding Ring-Atom-Bond Surfaces: %d Percent Done", PercentDone);
				TecUtilStatusSuspend(FALSE);
				TecUtilStatusSetPercentDoneText(PercentDoneText);
				if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
				TecUtilStatusSuspend(TRUE);
				if (IsStop) IsOk = FALSE;
			}

			/* Inform Tecplot that major data operation is ending */
			TecUtilLockFinish(AddOnID);
			TecUtilDataLoadEnd();
		}
		if (ZoneVarInfo != NULL)
			ZoneVarInfoDealloc(&ZoneVarInfo);
			// FREE_ITEM(ZoneVarInfo, "ZoneVarInfo");
		if (CircleGradPath != NULL)
			CircleGradPathDealloc(&CircleGradPath);

		if (ShowStatusBar)
		{
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusFinishPercentDone();
			TecUtilStatusSuspend(TRUE);
		}
	}

	//CPG :: add file dump
	//CPGraph.WriteNodes();


	//
	// Compute the Atom-Cage and Atom-FF(ring & cage) lines,
	// and the Atom-Cage-Ring, Atom-Cage-Bond, and Atom-FF surfaces.
	// Do this by extracting an isosurface around each atom with rho
	// a little larger than rho at the largest connecting bond CP. Then
	// compute the topology (min, saddle, max critical points, and
	// connecting gradient paths) of magnitude(grad-rho) on the surface.
	//
	// if (0)  // <<<<<<<<<<<<< TEMPORARY FOR TESTING <<<<<<<<<<<<<<
	if (IsOk && CompletionLevel >= CompletionLevel_AtomIsoSurfTopo)
	{
		// IsoTopo_pa IsoTopo = NULL;
		Boolean_t      GradGradVarsAdded = TRUE;
		SurfTopoSeg_pa SurfTopoSeg = NULL;
		EntIndex_t     IsoSurfZone = -1;
		EntIndex_t     IsoTopoZone = -1;
		EntIndex_t     SurfCPZoneNum = 0;
		LgIndex_t      AtomNum = 1;
		LgIndex_t      NumAtoms = CritPoints->NumCrtPtsM3;

		CritPoints_pa  SurfCritPoints = NULL;
		Bundles_pa     SurfBundles = NULL;
		SurfElemMap_pa SurfElemMap = NULL;
		Normals_pa     Normals     = NULL;
		IsoSurfGradPath_pa IsoSurfGradPath = NULL;

		// Unit test normals
		CHECK(NormalsTest());

		// Unit test CritPoints
		CHECK(CritPointsTest());

		// ZoneVarInfo_pa VolZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "VolZoneVarInfo");
		ZoneVarInfo_pa VolZoneVarInfo = ZoneVarInfoAlloc();
		VolZoneVarInfo->ZoneNum = ZoneNum;
		VolZoneVarInfo->UVarNum = UVarNum;
		VolZoneVarInfo->VVarNum = VVarNum;
		VolZoneVarInfo->WVarNum = WVarNum;
		VolZoneVarInfo->XVarNum = XVarNum;
		VolZoneVarInfo->YVarNum = YVarNum;
		VolZoneVarInfo->ZVarNum = ZVarNum;
		VolZoneVarInfo->DGradXNum = SGradGTXVarNum;
		VolZoneVarInfo->DGradYNum = SGradGTYVarNum;
		VolZoneVarInfo->DGradZNum = SGradGTZVarNum;
		VolZoneVarInfo->GradMagVarNum = GradMagVarNum;
		VolZoneVarInfo->ChrgDensVarNum = ChrgDensVarNum;
		VolZoneVarInfo->TypeVarNum = TypeVarNum;
		VolZoneVarInfo->NumVars = NumVars;
		VolZoneVarInfo->PeriodicBC = PeriodicBC;

		IsOk = ZoneVarInfoSetMinMax(VolZoneVarInfo, refinedGridSet);

		//	TIM:	Show status bar
		if (ShowStatusBar)
		{
			char PercentDoneText[200];
			char AtomBar[] = "[                    ]";
			memset(AtomBar, '|', 1 + (int)(20*(AtomNum-1)/NumAtoms));
			memset(AtomBar, '[', 1);
			if (CompletionLevel > CompletionLevel_AtomIsoSurfTopo)
				sprintf_s(PercentDoneText, "Finding Atom-Ring-Cage Surfaces: Atom %d of %d - 0 Percent Done  %s", AtomNum, NumAtoms, AtomBar);
			else
				sprintf_s(PercentDoneText, "Finding Atom Isosurface Topology: Atom %d of %d - 0 Percent Done  %s", AtomNum, NumAtoms, AtomBar);
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusStartPercentDone(PercentDoneText, TRUE, TRUE);
			TecUtilStatusSuspend(TRUE);
		}

		//	Set streamtrace variables for use in gradpath generation
		GradPathSTSetVariables(VolZoneVarInfo);

		// Loop over atoms to compute the FF lines, FF surface, etc.
		// ZoneVarInfo_pa SurfZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "SurfZoneVarInfo");
		ZoneVarInfo_pa SurfZoneVarInfo = ZoneVarInfoAlloc();
		IsOk = (SurfZoneVarInfo != NULL);
		//	TIM: for debugging to isolate around a particular atom
		for (AtomNum = 1; IsOk && AtomNum <= NumAtoms; AtomNum++)
		//for (AtomNum = 6; IsOk && AtomNum <= 6; AtomNum++)
		{
			if (IsOk)
			{
				// SurfZoneVarInfo->ZoneNum = ZoneNum;
				SurfZoneVarInfo->UVarNum = SGradGTXVarNum;
				SurfZoneVarInfo->VVarNum = SGradGTYVarNum;
				SurfZoneVarInfo->WVarNum = SGradGTZVarNum;
				SurfZoneVarInfo->ChrgDensVarNum = ChrgDensVarNum;
				SurfZoneVarInfo->GradMagVarNum = GradMagVarNum;
				SurfZoneVarInfo->TypeVarNum = TypeVarNum;
				SurfZoneVarInfo->PeriodicBC = PeriodicBC;
			}

			// Extract the largest isosurface surrounding only AtomNum (no other critical points) and extract the
			// surface topology on the isosurface. Store in IsoSurfGradPath structure.
			if (IsOk)
			{
				IsoSurfGradPath = IsoSurfGradPathAlloc();

				IsOk = IsoSurfGradPathAdd(VolZoneVarInfo, SurfZoneVarInfo, AtomNum - 1, CritPoints, Bundles, 1.0, IsoSurfGradPath);
			}

			// For convenience, extract pointers to isosurface CritPoints and Bundles
			if (IsOk)
			{
				SurfCritPoints = IsoSurfGradPath->SurfCritPoints;
				SurfBundles   = IsoSurfGradPath->SurfBundles;
			}


			// For each min on the IsoTopo surface, find the corresponding volume
			// Bond point.
			// TEMP: turn off for debugging
			// if (0)
			if (IsOk && CompletionLevel >= CompletionLevel_ARCSurfaces)
			{
				double	   VolGPSurfCPTolerance = SurfCritPoints->MinCPDistance;
				ArrList_pa BondForMinList    = ArrListAlloc(10, ArrListType_Long);
				ArrList_pa RingForSaddleList = NULL;
				ArrList_pa RingForMaxList    = NULL;

				LgIndex_t  NumMinCPs    = CritPointsGetEndOffset(SurfCritPoints, (char)(2))
										  - CritPointsGetBegOffset(SurfCritPoints, (char)(2));
				LgIndex_t  NumSaddleCPs = CritPointsGetEndOffset(SurfCritPoints, (char)(0))
										  - CritPointsGetBegOffset(SurfCritPoints, (char)(0));
				LgIndex_t  NumMaxCPs    = CritPointsGetEndOffset(SurfCritPoints, (char)(-2))
										  - CritPointsGetBegOffset(SurfCritPoints, (char)(-2));
				LgIndex_t  NumBonds     = CritPointsGetEndOffset(CritPoints, (char)(-1))
										  - CritPointsGetBegOffset(CritPoints, (char)(-1));
				LgIndex_t  NumRings     = CritPointsGetEndOffset(CritPoints, (char)(1))
										  - CritPointsGetBegOffset(CritPoints, (char)(1));
				LgIndex_t  NumCages     = CritPointsGetEndOffset(CritPoints, (char)(3))
										  - CritPointsGetBegOffset(CritPoints, (char)(3));

				
				// Find the mapping of surface Min to volume Bond
				CHECK(NumBonds >= NumMinCPs); 
				if (NumMinCPs > 0 && NumBonds > 0) // TEMP disabled for debugging
				{
					// Loop over MinCP's determining the mapping between volume Bond's and surface Min's
					LgIndex_t MinOffset;
					for (MinOffset = 0; IsOk && MinOffset < NumMinCPs; MinOffset++)
					{
						/* Update status bar */
						//	Atom completion range: 0 to 8
						if (ShowStatusBar)
						{
							char PercentDoneText[200];
							char AtomBar[] = "[                    ]";
							memset(AtomBar, '|', 1 + (int)(20*(AtomNum-1)/NumAtoms));
							memset(AtomBar, '[', 1);
							int  PercentDone;
							PercentDone = (int)(8.0 * (MinOffset + 1) / NumMinCPs);
							if (CompletionLevel > CompletionLevel_AtomIsoSurfTopo)
								sprintf_s(PercentDoneText, "Finding Atom-Ring-Cage Surfaces: Atom %d of %d - %d Percent Done  %s", 
									AtomNum, NumAtoms, PercentDone, AtomBar);
							else
								sprintf_s(PercentDoneText, "Finding Atom Isosurface Topology: Atom %d of %d - %d Percent Done  %s", 
									AtomNum, NumAtoms, PercentDone, AtomBar);
							TecUtilStatusSuspend(FALSE);
							TecUtilStatusSetPercentDoneText(PercentDoneText);
							if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
							TecUtilStatusSuspend(TRUE);
						}

						// LgIndex_t BondOffset = BondCPFromAtomIsoTopoMin(CritPoints, ITZCritPoints,
						//                                             AtomNum-1, MinOffset, 1.0);
						LgIndex_t BondOffset = VolCPFromAtomIsoTopoSurfCP((char)(-1), CritPoints, SurfCritPoints,
																			AtomNum - 1, (char)(2), MinOffset, VolGPSurfCPTolerance);
						/*LgIndex_t BondOffset = VolCPFromAtomIsoTopoSurfCPClosest((char)(-1), CritPoints, ITZCritPoints,
																			AtomNum - 1, (char)(2), MinOffset, 1.0);*/
						if (BondOffset >= 0)
						{
							ArrListItem_u Item;
							Item.Long = BondOffset + 1;
							IsOk = ArrListSetItem(BondForMinList, MinOffset, Item);
						}
						/*
						else  // Must be a Bond for each Min
						  IsOk = FALSE;
						CHECK(BondOffset >= 0);
						*/

						if (IsStop) IsOk = FALSE;
					}
				}

				//	The below mapped rings are all near field rings, so
				//	record which saddles have near field to know which 
				//	surf grad paths need to have their ends appended later.

				// Find the mapping of surface Saddle to volume Ring
				if (NumSaddleCPs > 0 && NumRings > 0)
				{
					// Loop over SaddleCP's determining the mapping between volume Rings and surface Saddles
					LgIndex_t SaddleOffset;
					RingForSaddleList    = ArrListAlloc(10, ArrListType_Long);
					for (SaddleOffset = 0; IsOk && SaddleOffset < NumSaddleCPs; SaddleOffset++)
					{
						/* Update status bar */
						//	Atom completion range: 8 to 16
						if (ShowStatusBar)
						{
							char PercentDoneText[200];
							char AtomBar[] = "[                    ]";
							memset(AtomBar, '|', 1 + (int)(20*(AtomNum-1)/NumAtoms));
							memset(AtomBar, '[', 1);
							int  PercentDone;
							PercentDone = (int)(8.0 + 8.0 * (SaddleOffset + 1) / NumSaddleCPs);
							if (CompletionLevel > CompletionLevel_AtomIsoSurfTopo)
								sprintf_s(PercentDoneText, "Finding Atom-Ring-Cage Surfaces: Atom %d of %d - %d Percent Done  %s", 
									AtomNum, NumAtoms, PercentDone, AtomBar);
							else
								sprintf_s(PercentDoneText, "Finding Atom Isosurface Topology: Atom %d of %d - %d Percent Done  %s", 
									AtomNum, NumAtoms, PercentDone, AtomBar);
							TecUtilStatusSuspend(FALSE);
							TecUtilStatusSetPercentDoneText(PercentDoneText);
							if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
							TecUtilStatusSuspend(TRUE);
						}
						
						ArrListItem_u Item;
						/*LgIndex_t RingOffset = VolCPFromAtomIsoTopoSurfCP((char)(1), CritPoints, ITZCritPoints,
																		  AtomNum - 1, (char)(0), SaddleOffset, 1.0);*/
						LgIndex_t RingOffset = VolCPFromAtomIsoTopoSurfCPClosest((char)(1), CritPoints, SurfCritPoints,
																		  AtomNum - 1, (char)(0), SaddleOffset, VolGPSurfCPTolerance * 0.5);
						if (RingOffset >= 0)
						{
							Item.Long = RingOffset;
						}
						else
						{
							Item.Long = -1;
						}
						IsOk = ArrListSetItem(RingForSaddleList, SaddleOffset, Item);

						if (IsStop) IsOk = FALSE;
					}
				}

				// Noise sometimes causes a saddle-max-saddle sequence where a saddle should be.
				// Test max's to see if a Ring-Atom connector passes within Tolerance.
				// TODO: Do something with this!
				if (NumMaxCPs > 0 && NumRings > 0) // TEMP disabled for debugging
				{
					// Loop over MinCP's determining the mapping between volume cages and surface maxes
					LgIndex_t MaxOffset;
					RingForMaxList    = ArrListAlloc(10, ArrListType_Long);
					for (MaxOffset = 0; IsOk && MaxOffset < NumMaxCPs; MaxOffset++)
					{
						/* Update status bar */
						//	Atom completion range: 16 to 24
						if (ShowStatusBar)
						{
							char PercentDoneText[200];
							char AtomBar[] = "[                    ]";
							memset(AtomBar, '|', 1 + (int)(20*(AtomNum-1)/NumAtoms));
							memset(AtomBar, '[', 1);
							int  PercentDone;
							PercentDone = (int)(16.0 + 8.0 * (MaxOffset + 1) / NumMaxCPs);
							if (CompletionLevel > CompletionLevel_AtomIsoSurfTopo)
								sprintf_s(PercentDoneText, "Finding Atom-Ring-Cage Surfaces: Atom %d of %d - %d Percent Done  %s", 
								AtomNum, NumAtoms, PercentDone, AtomBar);
							else
								sprintf_s(PercentDoneText, "Finding Atom Isosurface Topology: Atom %d of %d - %d Percent Done  %s", 
								AtomNum, NumAtoms, PercentDone, AtomBar);
							TecUtilStatusSuspend(FALSE);
							TecUtilStatusSetPercentDoneText(PercentDoneText);
							if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
							TecUtilStatusSuspend(TRUE);
						}
						
						ArrListItem_u Item;
						LgIndex_t RingOffset = VolCPFromAtomIsoTopoSurfCP((char)(1), CritPoints, SurfCritPoints,
																		  AtomNum - 1, (char)(-2), MaxOffset, VolGPSurfCPTolerance);
						/*LgIndex_t RingOffset = VolCPFromAtomIsoTopoSurfCPClosest((char)(1), CritPoints, ITZCritPoints,
																		  AtomNum - 1, (char)(-2), MaxOffset, 1.0);*/
						if (RingOffset >= 0)
						{
							Item.Long = RingOffset;
						}
						else
						{
							Item.Long = -1;
						}

						IsOk = ArrListSetItem(RingForMaxList, MaxOffset, Item);

						if (IsStop) IsOk = FALSE;
					}
				}


				// Find the mapping of surface Max to volume Cage
				/*
				if (NumMaxCPs > 0 && NumCages > 0)
				  {
					// Loop over MinCP's determining the mapping between volume Cage's and surface Max's
					LgIndex_t MaxOffset;
					for (MaxOffset=0; IsOk && MaxOffset < NumMaxCPs; MaxOffset++)
					  {
						LgIndex_t CageOffset = CageCPFromAtomIsoTopoMax(CritPoints, ITZCritPoints,
																		AtomNum-1, CageOffset, 1.0);
						if (SaddleOffset >= 0)
						  {
							ArrListItem_u Item;
							Item.Long = SaddleOffset;
							IsOk = ArrListSetItem(SaddleForMinList, SaddleOffset, Item);
						  }
					  }
				  }
				  */

				// If Bundles (for volume zone) don't yet exist, create them
				if (IsOk)
				{
					LgIndex_t NumSurfBundles = BundlesGetCount(SurfBundles);

					if (Bundles == NULL)
					{
						Bundles = BundlesAlloc();
					}

					// Loop over surface bundles creating corresponding volume bundles
					// for (sboff=0; IsOk && sboff<NumSurfBundles; sboff++)
				}

				// Add a volume gradient path for each surface critical point that doesn't
				// already have a gradient path, and a volume surface for each new surface
				// line that doesn't already have a surface defined previously (saddle-max
				// lines may be represented by an bond-ring-cage defined by a bond
				// CircleGradPath.
				if (IsOk)
				{
					LgIndex_t MaxNum, SaddleNum;

					for (MaxNum = 0; IsOk && MaxNum < NumMaxCPs; MaxNum++)
					{
						/* Update status bar */
						//	Atom completion range: 24 to 32
						if (ShowStatusBar)
						{
							char PercentDoneText[200];
							char AtomBar[] = "[                    ]";
							memset(AtomBar, '|', 1 + (int)(20*(AtomNum-1)/NumAtoms));
							memset(AtomBar, '[', 1);
							int  PercentDone;
							PercentDone = (int)(24.0 + 8.0 * (MaxNum + 1) / NumMaxCPs);
							if (CompletionLevel > CompletionLevel_AtomIsoSurfTopo)
								sprintf_s(PercentDoneText, "Finding Atom-Ring-Cage Surfaces: Atom %d of %d - %d Percent Done  %s", 
								AtomNum, NumAtoms, PercentDone, AtomBar);
							else
								sprintf_s(PercentDoneText, "Finding Atom Isosurface Topology: Atom %d of %d - %d Percent Done  %s", 
								AtomNum, NumAtoms, PercentDone, AtomBar);
							TecUtilStatusSuspend(FALSE);
							TecUtilStatusSetPercentDoneText(PercentDoneText);
							if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
							TecUtilStatusSuspend(TRUE);
							if (IsStop) IsOk = FALSE;
						}

						GradPath_pa GradPath = NULL;
						XYZ_s SeedPos;
						double X, Y, Z, Rho, Px, Py, Pz;
						char Type;
						LgIndex_t NumPathPoints = 0;

						if (IsOk) IsOk = CritPointsGetPoint(SurfCritPoints, MaxNum, &X, &Y, &Z, &Rho, &Type, &Px, &Py, &Pz);

						if (IsOk)
						{
							EntIndex_t ConnectorZoneNum = 0;
							SeedPos.X = X;
							SeedPos.Y = Y;
							SeedPos.Z = Z;

							GradPath = GradPathAddMidway(VolZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
														 SeedPos, 1.0, &NumPathPoints);

							if (!GradPathIsValid(GradPath)) IsOk = FALSE;

							if (IsOk && GradPath->EndCrtPtNum == AtomNum)
								IsOk = GradPathReverse(&GradPath);

							if (IsOk && GradPath->EndCrtPtNum == -1)
							{
								LgIndex_t EndCrtPtNum;
								// TODO Extrapolate to the outer boundary and change last value in GradPath

								// Add a new far-field cage to CritPoints
								IsOk = GradPathGetPoint(GradPath, NumPathPoints - 1, &X, &Y, &Z, &Rho);
								if (IsOk)
									IsOk = CritPointsAppendPoint(CritPoints, X, Y, Z, Rho, (char)(13), 0.0, 0.0, 0.0);

								EndCrtPtNum = CritPointsGetEndOffset(CritPoints, (char)(13)) - 1;
								GradPath->EndCrtPtNum = EndCrtPtNum;
							}

							if (IsOk && GradPath->BeginCrtPtNum >= 0){
								ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, GradPath);

								if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
							}
						}

						/* Clean up allocated structures/arrays */
						if (GradPath != NULL) GradPathDealloc(&GradPath);
					}

					// Add the Atom-Ring lines
					for (SaddleNum = 0; IsOk && SaddleNum < NumSaddleCPs; SaddleNum++)
					{
						/* Update status bar */
						//	Atom completion range: 32 to 40
						if (ShowStatusBar)
						{
							char PercentDoneText[200];
							char AtomBar[] = "[                    ]";
							memset(AtomBar, '|', 1 + (int)(20*(AtomNum-1)/NumAtoms));
							memset(AtomBar, '[', 1);
							int  PercentDone;
							PercentDone = (int)(32.0 + 8.0 * (SaddleNum + 1) / NumSaddleCPs);
							if (CompletionLevel > CompletionLevel_AtomIsoSurfTopo)
								sprintf_s(PercentDoneText, "Finding Atom-Ring-Cage Surfaces: Atom %d of %d - %d Percent Done  %s", 
								AtomNum, NumAtoms, PercentDone, AtomBar);
							else
								sprintf_s(PercentDoneText, "Finding Atom Isosurface Topology: Atom %d of %d - %d Percent Done  %s", 
								AtomNum, NumAtoms, PercentDone, AtomBar);
							TecUtilStatusSuspend(FALSE);
							TecUtilStatusSetPercentDoneText(PercentDoneText);
							if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
							TecUtilStatusSuspend(TRUE);
							if (IsStop) IsOk = FALSE;
						}
						
						GradPath_pa GradPath = NULL;
						XYZ_s SeedPos;
						double X, Y, Z, Rho, Px, Py, Pz;
						char Type;
						LgIndex_t NumPathPoints = 0;
						LgIndex_t SurfCPNum = SaddleNum + CritPointsGetBegOffset(SurfCritPoints, (char)(0));

						ArrListItem_u Item;

						// If not found to already exist, add
						Item.Long = -1;
						if (RingForSaddleList != NULL && ArrListGetCount(RingForSaddleList) > 0)
							Item = ArrListGetItem(RingForSaddleList, SaddleNum);
						if (Item.Long < 0)
						{

							if (IsOk) IsOk = CritPointsGetPoint(SurfCritPoints, SurfCPNum, &X, &Y, &Z, &Rho, &Type, &Px, &Py, &Pz);

							if (IsOk)
							{
								EntIndex_t ConnectorZoneNum = 0;
								SeedPos.X = X;
								SeedPos.Y = Y;
								SeedPos.Z = Z;

								GradPath = GradPathAddMidway(VolZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
															 SeedPos, 1.0, &NumPathPoints);

								if (!GradPathIsValid(GradPath)) IsOk = FALSE;

								if (IsOk && GradPath->EndCrtPtNum == AtomNum)
									IsOk = GradPathReverse(&GradPath);

								if (IsOk && GradPath->EndCrtPtNum == -1)
								{
									LgIndex_t EndCrtPtNum;
									// TODO Extrapolate to the outer boundary and change last value in GradPath

									// Add a new far-field ring to CritPoints
									IsOk = GradPathGetPoint(GradPath, NumPathPoints - 1, &X, &Y, &Z, &Rho);
									if (IsOk)
									{
										LgIndex_t PointOffset = CritPoints->NumCrtPtsM3 + CritPoints->NumCrtPtsM1
																+ CritPoints->NumCrtPtsP1 + CritPoints->NumCrtPtsP3
																+ CritPoints->NumFFCrtPtsP1;
										IsOk =  CritPointsInsertPoint(CritPoints, PointOffset, X, Y, Z,
																	  Rho, (char)(11), 0.0, 0.0, 0.0);
									}

									EndCrtPtNum = CritPointsGetEndOffset(CritPoints, (char)(11)) - 1;
									GradPath->EndCrtPtNum = EndCrtPtNum;
								}

								if (IsOk && GradPath->BeginCrtPtNum >= 0){
									ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, GradPath);

									if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
								}
							}

							/* Clean up allocated structures/arrays */
							if (GradPath != NULL) GradPathDealloc(&GradPath);
						}
					}
				}

				//	TIM:	loop over all surface connectors and replace/append
				//			their end points with volGP-isosurf intersection points
				
				IsOk = IsoSurfGradPathSurfConnectSurfSaddleToVolRing(IsoSurfGradPath, SurfCritPoints, CritPoints, AtomNum-1);

				// Add a volume zero-flux surface for each surface critical-point
				// connector that doesn't already have a surface defined previously
				// (for example, saddle-max lines may be represented by an
				// ring-atom-cage surface defined by a ring CircleGradPath).

				//	TIM:	connector = line on the isosurface (ie. from max to min)

				// First do all of the atom-ring-cage zero-flux surfaces
				if (IsOk && NumMinCPs > 1)  // No surfaces for Hydrogen Atom
					// if (0)   // TEMP
				{
					LgIndex_t   Bundle;
					GradPath_pa SurfGradPath = NULL;

					OrderedPairs_pa SaddleMaxPairs = NULL;
					OrderedPairs_pa BundlePairs = NULL;
					LgIndex_t       PairOffset;

					EntIndex_t  ZeroFluxSurfZoneNum = TECUTILSETNOTMEMBER;

					LgIndex_t IMax = brSurfaceIMax;

					EntIndex_t SurfSourceZoneNum = SurfCritPoints->SourceZoneNum;

					/* Inform Tecplot that major data operation is beginning */
					TecUtilDataLoadBegin();
					TecUtilLockStart(AddOnID);

					// Allocate ordered-pairs array for avoiding duplicate surfaces
					SaddleMaxPairs = OrderedPairsAlloc();
					BundlePairs    = OrderedPairsAlloc();
					LgIndex_t NumOfSurfBundles = BundlesGetCount(SurfBundles);

					// Cycle through the SurfBundles to find the pairs of Bundles that have the
					// same Saddle and Max.
					for (Bundle = 0; IsOk && Bundle < NumOfSurfBundles; Bundle++)
					{
						// EntIndex_t  TPZoneForSurfLine;
						LgIndex_t   Min, Saddle, Max, Atom;
						ArrList_pa  GradPathList = ArrListAlloc(10, ArrListType_VoidPtr);

						IsOk = BundlesGetFromOffset(SurfBundles, Bundle, &Atom, &Min, &Saddle, &Max);

						PairOffset = OrderedPairsGetOffsetFromPair(SaddleMaxPairs, Saddle, Max);

						// If BundlePair exists, set the second bundle number, else create the BundlePair
						// and the SaddleMaxPair
						if (PairOffset >= 0)
						{
							LgIndex_t Bundle1, Bundle2;
							IsOk = OrderedPairsGetFromOffset(BundlePairs, PairOffset, &Bundle1, &Bundle2);
							Bundle2 = Bundle;
							IsOk = OrderedPairsSetAtOffset(BundlePairs, PairOffset, Bundle1, Bundle);
						}
						else
						{
							IsOk = OrderedPairsAppendAtEnd(SaddleMaxPairs, Saddle, Max);
							IsOk = OrderedPairsAppendAtEnd(BundlePairs, Bundle, 0);
						}
					}

					CHECK(OrderedPairsGetCount(SaddleMaxPairs) == OrderedPairsGetCount(BundlePairs));


					// Cycle through the BundlePairs (and SaddleMaxPairs) and add the
					// Atom-Ring-Cage Surfaces.
					LgIndex_t NumOfOrderedPairs = OrderedPairsGetCount(BundlePairs);
					for (PairOffset = 0; IsOk && PairOffset < NumOfOrderedPairs; PairOffset++)
					{
						/* Update status bar */
						//	Atom completion range: 40 to 100
						if (ShowStatusBar)
						{
							char PercentDoneText[200];
							char AtomBar[] = "[                    ]";
							memset(AtomBar, '|', 1 + (int)(20*(AtomNum-1)/NumAtoms));
							memset(AtomBar, '[', 1);
							int  PercentDone;
							PercentDone = (int)(40.0 + 60.0 * (PairOffset + 1) / NumOfOrderedPairs);
							if (CompletionLevel > CompletionLevel_AtomIsoSurfTopo)
								sprintf_s(PercentDoneText, "Finding Atom-Ring-Cage Surfaces: Atom %d of %d - %d Percent Done  %s", 
								AtomNum, NumAtoms, PercentDone, AtomBar);
							else
								sprintf_s(PercentDoneText, "Finding Atom Isosurface Topology: Atom %d of %d - %d Percent Done  %s", 
								AtomNum, NumAtoms, PercentDone, AtomBar);
							TecUtilStatusSuspend(FALSE);
							TecUtilStatusSetPercentDoneText(PercentDoneText);
							if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
							TecUtilStatusSuspend(TRUE);
							if (IsStop) IsOk = FALSE;
						}


						LgIndex_t   ii;
						// EntIndex_t  TPZoneForSurfLine;
						LgIndex_t   BundleOffset1, BundleOffset2;
						LgIndex_t   Min1, Saddle, Max, Atom;
						LgIndex_t   Min2, Saddle2, Max2, Atom2;
						LgIndex_t   Bond1 = -1;
						LgIndex_t   Bond2 = -1;

						ArrList_pa  GradPathList = ArrListAlloc(10, ArrListType_VoidPtr);

						if (IsOk) IsOk = OrderedPairsGetFromOffset(BundlePairs, PairOffset, &BundleOffset1, &BundleOffset2);

						if (IsOk) IsOk = BundlesGetFromOffset(SurfBundles, BundleOffset1, &Atom, &Min1, &Saddle, &Max);
						if (IsOk) IsOk = BundlesGetFromOffset(SurfBundles, BundleOffset2, &Atom2, &Min2, &Saddle2, &Max2);


						//	TIM:	This doesn't seem to do anything. 
 /*                       static Boolean_t keepChecking = TRUE;
						if ( keepChecking )
						{
							CHECK(Atom2 == Atom);
							CHECK(Saddle2 == Saddle);
							CHECK(Max2 == Max);
							keepChecking = (Atom2 == Atom && Saddle2 == Saddle && Max2 == Max);
						}*/

						// Extract the Bond numbers for the two Mins from the BondForMinList
						if (IsOk)
						{
							ArrListItem_u Item;
							LgIndex_t NumInBondList = ArrListGetCount(BondForMinList);

							if (Min1 <= NumInBondList)
							{
								Item = ArrListGetItem(BondForMinList, Min1 - 1);
								Bond1 = (LgIndex_t)Item.Long;
							}
							if (Min2 <= NumInBondList)
							{
								Item = ArrListGetItem(BondForMinList, Min2 - 1);
								Bond2 = (LgIndex_t)Item.Long;
							}
						}

						if (IsOk)
						{
							LgIndex_t SaddleCPOff = Saddle - 1 + CritPointsGetBegOffset(SurfCritPoints, (char)(-1));
							LgIndex_t MaxCPOff    = Max - 1;
							
							// TPZoneForSurfLine = GradPathTPZoneFromBegEndCP(Saddle, Max, SurfSourceZoneNum);
							// SurfGradPath = GradPathGetFromTPZone(TPZoneForSurfLine);
							SurfGradPath = IsoSurfGradPathGetGPByBegEndCP(IsoSurfGradPath, SaddleCPOff, MaxCPOff);
							LgIndex_t NumOfSurfGPPoints = GradPathGetCount(SurfGradPath);
							if (NumOfSurfGPPoints <= 0) IsOk = FALSE;

							//	For storing the last grad path created below, knowing it's the
							//	atom-cage line
							GradPath_pa AtomCageGP = NULL;

							// For each point on the SurfGradPath, seed a VolGradPath that
							// Becomes and I-line of the zero-flux surface zone.
							GradPath_pa RingAtomGP = NULL, RingCageGP = NULL, AtomRingCageGP = NULL;
							
							for (ii = 0; IsOk && ii < NumOfSurfGPPoints; ii++)
							{
								GradPath_pa GradPath = GradPathAlloc();
								XYZ_s SeedPos;
								double X, Y, Z, Rho;
								LgIndex_t NumPathPoints = 0;
													
								// If not found to already exist, add
								IsOk = GradPathGetPoint(SurfGradPath, ii, &X, &Y, &Z, &Rho);

								if (IsOk)
								{
									EntIndex_t ConnectorZoneNum = 0;
									SeedPos.X = X;
									SeedPos.Y = Y;
									SeedPos.Z = Z;

									GradPath = GradPathAddMidway(VolZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
																 SeedPos, 1.0, &NumPathPoints);

								}
								if (!GradPathIsValid(GradPath)) IsOk = FALSE;

								if (IsOk)
									IsOk = GradPathResample(&GradPath, IMax, atomISegmentRatio);

// 								if (IsOk && GradPath->BeginCrtPtNum != AtomNum - 1 && GradPath->EndCrtPtNum != AtomNum - 1){
// 									char BegType = 0, EndType = 0;
// 									if (GradPath->BeginCrtPtNum >= 0)
// 										BegType = CritPointsGetType(CritPoints, GradPath->BeginCrtPtNum);
// 									if (GradPath->EndCrtPtNum >= 0)
// 										EndType = CritPointsGetType(CritPoints, GradPath->EndCrtPtNum);
// 
// 									if (BegType == 0 && EndType != -3){
// 											GradPath->BeginCrtPtNum = AtomNum - 1;
// 									}
// 									if (EndType == 0 && BegType != -3 && EndType != BegType){
// 											GradPath->EndCrtPtNum = AtomNum - 1;
// 									}
// 								}

								if (IsOk)
								{
									ArrListItem_u Item;
									Item.VoidPtr = (void *)GradPath;
									IsOk = ArrListAppendItem(GradPathList, Item);
									if (ii == NumOfSurfGPPoints - 1) AtomCageGP = GradPath;
									if (ii == 0) RingAtomGP = GradPath;
								}


								// if (IsOk) ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, GradPath);
								// if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE
							}

							//	Since SurfGradPath starts at saddle (corresponds to atom-ring),
							//	assume that the first SurfGradPath point already lies on the
							//	atom-ring line. Fetch that line, combine it with the ring-cage
							//	line (assumed to exist for near-field ring-cage), and 
							//	resample the combined line.

							Boolean_t	HasRingCageGP = FALSE;
							if (RingAtomGP != NULL && AtomCageGP != NULL && RingAtomGP->BeginCrtPtNum >= 0 && AtomCageGP->BeginCrtPtNum >= 0){
								char TempTypeBeg = CritPointsGetType(CritPoints, RingAtomGP->BeginCrtPtNum);
								char TempTypeEnd;
								if (RingAtomGP->EndCrtPtNum == -1) TempTypeEnd = 11;
								else  TempTypeEnd = CritPointsGetType(CritPoints, RingAtomGP->EndCrtPtNum);
								char RingCPType = -1;

								LgIndex_t RingCPNum = -1;

								if (TempTypeBeg != (char)1 && TempTypeBeg != (char)11
									&& TempTypeEnd != (char)1 && TempTypeEnd != (char)11)
								{
									LgIndex_t RingOffset = VolCPFromAtomIsoTopoSurfCPClosest((char)1,
										CritPoints,
										SurfCritPoints,
										AtomNum - 1,
										(char)0,
										Saddle - 1,
										VolGPSurfCPTolerance*2.0);
									if (RingOffset >= 0)
									{
										if (IsOk) RingCPNum = CritPointsGetBegOffset(CritPoints, (char)1) + RingOffset;

										/*RingAtomGP = GradPathVolGPFromAtomIsoTopoSurfCPClosest((char)1, (char)11, CritPoints,
																								ITZCritPoints, AtomNum-1, (char)0,
																								Saddle-1, VolGPSurfCPTolerance*2.0);*/

										RingAtomGP = GradPathGetByBegEndCP(RingCPNum, AtomNum - 1);
									}

									ENSURE(GradPathIsValid(RingAtomGP));
									IsOk = GradPathIsValid(RingAtomGP);
								}
								else
								{
									if (TempTypeEnd == (char)11 || TempTypeEnd == (char)1)
									{
										RingCPNum = RingAtomGP->EndCrtPtNum;
										RingCPType = TempTypeEnd;
									}
									else if (TempTypeBeg == (char)11 || TempTypeBeg == (char)1)
									{
										RingCPNum = RingAtomGP->BeginCrtPtNum;
										RingCPType = TempTypeBeg;
									}
								}


								//	Get the cage number from the last grad path created above
								LgIndex_t	CageCPNum = AtomCageGP->EndCrtPtNum;
								HasRingCageGP = TRUE;
								int			ResampleWeight = 0;


								if (IsOk)
								{
									//if (CageOffset < 0) IsOk = FALSE;

									char	CageCPType;
									if (CageCPNum != -1)
										CageCPType = CritPointsGetType(CritPoints, CageCPNum);
									else
										CageCPType = (char)13;

									//if (IsOk) CageCPNum = CritPointsGetBegOffset(CritPoints, (char)3) + CageOffset;

									if (CageCPType == (char)13 && RingCPType == (char)11)
									{
										//	Both ring and cage are far field, so a ring-cage line isn't
										//	necessary
										HasRingCageGP = FALSE;
									}
									else if (CageCPType == (char)13 && RingCPType == (char)1)
									{
										//	Cage is farfield but ring is not, so ring-cage line will
										//	not terminate at the same cage cp as the atom-cage line.
										//	Instead, use the existing ring-ff line or create it if it
										//	can't be found.
										RingCageGP = GradPathGetByBegEndCP(RingCPNum, (LgIndex_t)-1);
										ResampleWeight = 1;
									}
									else if (CageCPType == (char)3 && RingCPType == 1)
									{
										//	Both ring and cage are nearfield, so use the existant line.
										RingCageGP = GradPathGetByBegEndCP(RingCPNum, CageCPNum);
										ResampleWeight = 0;
									}
									else
									{
										HasRingCageGP = FALSE;
									}

									//ENSURE(GradPathIsValid(RingCageGP));

									//	If no ring-cage grad path, create one.
									//	Since ring-cage lines were created earlier, this can be assumed 
									//	to be a far field ring.
									//if (!GradPathIsValid(RingCageGP))
									//{
									//	//	The beginning and end CPs of the desired gradpath are known.
									//	RingCageGP = GradPathAlloc();
									//	ENSURE(GradPathIsValid(RingCageGP));
									//	RingCageGP->BeginCrtPtNum = RingCPNum;
									//	RingCageGP->EndCrtPtNum = CageCPNum;

									//	IsOk = GradPathAdd2PtMinLen(VolZoneVarInfo, IMax, CritPoints, RingCageGP);
									//	
									//	//	TIM TODO:	Find out why the below code doesn't return the right cage point.
									//	////	Seed a small direction from the ring CP in the principal direction
									//	//double XPos, YPos, ZPos, Rho, 
									//	//		XCrtPt, YCrtPt, ZCrtPt,
									//	//		PrincDirX, PrincDirY, PrincDirZ;
									//	//char	junk;
									//	//LgIndex_t	NumOfPathPoints = 0;
									//	//IsOk = CritPointsGetPoint(CritPoints, RingCPNum, &XCrtPt, &YCrtPt, &ZCrtPt, 
									//	//						&Rho, &junk, &PrincDirX, &PrincDirY, &PrincDirZ);
									//	////	Attempt to add the grad path. If the seed point is out of bounds or
									//	////	the grad path doesn't terminate at the desired cage then switch direction.
									//	//Boolean_t Break = FALSE;
									//	//for (int Dir = 1 ; !Break && Dir > -2 ; Dir -= 2)
									//	//{
									//	//	if (IsOk)
									//	//	{
									//	//		XPos = XCrtPt + DXYZ * PrincDirX * Dir;
									//	//		YPos = YCrtPt + DXYZ * PrincDirY * Dir;
									//	//		ZPos = ZCrtPt + DXYZ * PrincDirZ * Dir;

									//	//		RingCageGP = GradPathAlloc();
									//	//		IsOk = GradPathSetPoint(RingCageGP, (LgIndex_t)0, XPos, YPos, ZPos, Rho);
									//	//		if (IsOk) IsOk = GradPathIsValid(RingCageGP);
									//	//		if (IsOk) RingCageGP->BeginCrtPtNum = RingCPNum;
									//	//	}
									//	//	if (IsOk)
									//	//	{
									//	//		IsOk = GradPathAdd(VolZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Reverse, 
									//	//							CritPoints, 1.0, &NumOfPathPoints, RingCageGP);
									//	//	}
									//	//	if (RingCageGP->EndCrtPtNum == CageCPNum) Break = TRUE;
									//	//}
									//	//if (!Break) IsOk = FALSE;
									//}
								}

								if (HasRingCageGP)
								{
									if (IsOk)
									{
										//	Both the ring-atom and ring-cage originate at the ring.
										//	Also, the grad paths that will make up the surface originate
										//	at the atom, so to conform with that, reverse the ring-atom
										//	line and then concatenate the ring-cage line to it, producing
										//	an atom-ring-cage line.
										if (TempTypeEnd != (char)1 && TempTypeEnd != (char)11)
											IsOk = GradPathReverse(&RingAtomGP);

										//	Need to resample the two lines to have a combined total of IMax
										//	points, but also need a point to fall at the intersection.
										double	AtomRingLen = 0.0, RingCageLen = 0.0;
										int		AtomRingNumPoints = 0;

										AtomRingLen = GradPathGetLength(RingAtomGP);
										if (GradPathIsValid(RingCageGP)) RingCageLen = GradPathGetLength(RingCageGP);
										else IsOk = FALSE;

										if (IsOk)
										{
											AtomRingNumPoints = (int)(AtomRingLen / (AtomRingLen + RingCageLen) * IMax) + ResampleWeight;
											IsOk = (AtomRingNumPoints > 0);

											if (IsOk) IsOk = GradPathResample(&RingAtomGP, AtomRingNumPoints + 1, atomISegmentRatio);
											if (IsOk) IsOk = GradPathResample(&RingCageGP, IMax - AtomRingNumPoints, atomISegmentRatio);

											//	Remove last point from the atom-ring line so that the j lines match up better.
											if (IsOk) IsOk = GradPathRemovePoint(RingAtomGP, AtomRingNumPoints);

											if (IsOk) GradPathAppend(RingAtomGP, RingCageGP);

											if (IsOk) AtomRingCageGP = RingAtomGP;
										}
									}

									IsOk = (AtomRingCageGP != NULL && AtomRingCageGP->BeginCrtPtNum == AtomNum - 1
										&& (AtomRingCageGP->EndCrtPtNum == CageCPNum || AtomRingCageGP->EndCrtPtNum == (char)(-1)));

									if (!GradPathIsValid(AtomRingCageGP)) IsOk = FALSE;

									if (IsOk)
									{
										ArrListItem_u Item;
										Item.VoidPtr = (void *)AtomRingCageGP;
										ArrListRemoveItem(GradPathList, (LgIndex_t)0);
										IsOk = ArrListInsertItem(GradPathList, (LgIndex_t)0, Item);
										ZeroFluxSurfZoneNum = CreateZFSurfZoneFrom3CPsGradpathList(ChrgDensVarNum,
											TypeVarNum, CritPoints,
											AtomNum - 1, RingCPNum, CageCPNum,
											GradPathList, IMax);
									}


									// Create the surface zone by using each GradPath, resampled,
									// as the I-line of a 2D ordered zone.
								}
							}
							if (!HasRingCageGP && RingAtomGP != NULL && AtomCageGP != NULL && RingAtomGP->BeginCrtPtNum >= 0 && AtomCageGP->BeginCrtPtNum >= 0)
							{
								ZeroFluxSurfZoneNum = CreateZFSurfZoneFromGPList(ChrgDensVarNum, TypeVarNum,
																			 CritPoints, GradPathList, IMax);
							}

// 							if (!IsOk)
// 							{
// 								char Message[200];
// 								sprintf_s(Message, "ZFS was not created for pair offset %d of atom number %d (atomnum = atom number - 1)",
// 									PairOffset, AtomNum);
// 								TecUtilDialogMessageBox(Message,MessageBox_Error);
// 							}

							// Add the bond aux data (bonds on either side of surface)
							if (IsOk && ZeroFluxSurfZoneNum != TECUTILSETNOTMEMBER)
							{
								char BondString[200];

								AuxData_pa ZoneAuxDataRef = TecUtilAuxDataZoneGetRef(ZeroFluxSurfZoneNum);
								IsOk = TecUtilAuxDataSetStrItem(ZoneAuxDataRef, "CompChem.OppositeCrtPtType", "Bond", TRUE);

								if (IsOk)
								{
									sprintf_s(BondString, "%d", Bond1);
									IsOk = TecUtilAuxDataSetStrItem(ZoneAuxDataRef, "CompChem.OppositeCrtPtNum1", BondString, TRUE);
								}

								if (IsOk)
								{
									sprintf_s(BondString, "%d", Bond2);
									IsOk = TecUtilAuxDataSetStrItem(ZoneAuxDataRef, "CompChem.OppositeCrtPtNum2", BondString, TRUE);
								}
							}
							else
							{
								IsOk = FALSE;
							}


							// Dealloc the volume GradPaths in the GradPathList and then
							// dealloc the list
							for (ii = 0; IsOk && ii < GradPathGetCount(SurfGradPath); ii++)
							{
								ArrListItem_u Item;
								GradPath_pa VolGradPath = NULL;
								Item = ArrListGetItem(GradPathList, ii);
								VolGradPath = (GradPath_pa)Item.VoidPtr;
								GradPathDealloc(&VolGradPath);
							}
							// RingAtomGP is deallocated when AtomRingCageGP is deallocated
							//GradPathDealloc(&RingAtomGP);
							if (IsOk)
								GradPathDealloc(&RingCageGP);
						}

						ArrListDealloc(&GradPathList);


						// GradPathDealloc(&SurfGradPath);
					}

					// Dealloc the ordered-pairs list used to avoid duplicate surfaces
					OrderedPairsDealloc(&SaddleMaxPairs);

					/* Inform Tecplot that major data operation is ending */
					TecUtilLockFinish(AddOnID);
					TecUtilDataLoadEnd();
				}

				// Now do all of the atom-cage-bond zero-flux surfaces
				// if (IsOk)
				#if(0)   // TEMP
				
					LgIndex_t   Bundle;
					GradPath_pa SurfGradPath = NULL;

					EntIndex_t  ZeroFluxSurfZoneNum = TECUTILSETNOTMEMBER;

					LgIndex_t IMax = brSurfaceIMax;

					EntIndex_t SurfSourceZoneNum = SurfCritPoints->SourceZoneNum;

					// Inform Tecplot that major data operation is beginning
					TecUtilDataLoadBegin();
					TecUtilLockStart(AddOnID);

					// Cycle through the SurfBundles and add the Atom-Ring-Cage Surfaces.
					for (Bundle = 0; IsOk && Bundle < BundlesGetCount(SurfBundles); Bundle++)
					{
						LgIndex_t   ii;
						// EntIndex_t  TPZoneForSurfLine;
						LgIndex_t   Min, Saddle, Max, Atom;
						ArrList_pa  GradPathList = ArrListAlloc(10, ArrListType_VoidPtr);

						IsOk = BundlesGetFromOffset(SurfBundles, Bundle, &Atom, &Min, &Saddle, &Max);

						// TPZoneForSurfLine = GradPathTPZoneFromBegEndCP(Max, Min, SurfSourceZoneNum);
						// SurfGradPath = GradPathGetFromTPZone(TPZoneForSurfLine);
						SurfGradPath = IsoSurfGradPathGetGPByBegEndCP(IsoSurfGradPath, Saddle, Max);

						// For each point on the SurfGradPath, seed a VolGradPath that
						// Becomes and I-line of the zero-flux surface zone.
						for (ii = 0; IsOk && ii < GradPathGetCount(SurfGradPath) - 1; ii++)
						{
							GradPath_pa GradPath = GradPathAlloc();
							XYZ_s SeedPos;
							double X, Y, Z, Rho;
							LgIndex_t NumPathPoints = 0;

							// If not found to already exist, add
							IsOk = GradPathGetPoint(SurfGradPath, ii, &X, &Y, &Z, &Rho);

							if (IsOk)
							{
								EntIndex_t ConnectorZoneNum = 0;
								SeedPos.X = X;
								SeedPos.Y = Y;
								SeedPos.Z = Z;

								GradPath = GradPathAddMidway(VolZoneVarInfo, MAXGRADPATHPOINTS, CritPoints,
															 SeedPos, 1.0, &NumPathPoints);

								if (!GradPathIsValid(GradPath)) IsOk = FALSE;

								// If the end critical point is a bond, substitute the line for the bond-atom
								// (plus the line for the bond-cage?)
								/*
								if (IsOk)
								  {
									LgIndex_t EndCrtPtNum = GradPath->EndCrtPtNum;
									if (EndCrtPtNum > 0 && CritPointsGetType(CritPoints, EndCrtPtNum) == (char)(-1))
									  {
										EntIndex_t BondAtomZone = GradPathTPZoneFromBegEndCP(EndCrtPtNum, Atom-1, ZoneNum);
										GradPathDealloc(&GradPath);
										GradPath = GradPathGetFromTPZone(BondAtomZone);
										IsOk = GradPathReverse(&GradPath);
									  }
								  }
								  */

								if (IsOk)
									IsOk = GradPathResample(&GradPath, IMax, atomISegmentRatio);

								if (IsOk)
								{
									ArrListItem_u Item;
									Item.VoidPtr = (void *)GradPath;
									IsOk = ArrListAppendItem(GradPathList, Item);
								}

							}
						}

						// End point of last GradPath should be a bond, substitute
						// the appropriate bond GradPath
						if (IsOk)
						{
							ArrListItem_u Item;
							LgIndex_t BondNum, EndCrtPtNum, MinOffset;
							GradPath_pa GradPath = GradPathAlloc();

							MinOffset = Min - CritPointsGetBegOffset(SurfCritPoints, (char)(2));
							Item = ArrListGetItem(BondForMinList, MinOffset);
							BondNum = Item.Long;
							EndCrtPtNum = BondNum + CritPointsGetBegOffset(CritPoints, (char)(-1));

							if (EndCrtPtNum >= 0)
							{
								// EntIndex_t BondAtomZone = GradPathTPZoneFromBegEndCP(EndCrtPtNum, Atom-1, ZoneNum);
								EntIndex_t BondAtomZone = GradPathTPZoneFromBegEndCP(EndCrtPtNum, Atom, ZoneNum);
								GradPathDealloc(&GradPath);
								GradPath = GradPathGetFromTPZone(BondAtomZone);
								IsOk = GradPathReverse(&GradPath);
							}

							if (IsOk)
								IsOk = GradPathResample(&GradPath, IMax, atomISegmentRatio);

							if (IsOk)
							{
								ArrListItem_u Item;
								Item.VoidPtr = (void *)GradPath;
								IsOk = ArrListAppendItem(GradPathList, Item);
							}
						}

						// Create the surface zone by using each GradPath, resampled,
						// as the I-line of a 2D ordered zone.
						ZeroFluxSurfZoneNum = CreateZFSurfZoneFromGPList(ChrgDensVarNum, TypeVarNum,
																		 CritPoints, GradPathList, IMax);

						if (ZeroFluxSurfZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;


						// Dealloc the volume GradPaths in the GradPathList and then
						// dealloc the list
						for (ii = 0; IsOk && ii < GradPathGetCount(SurfGradPath); ii++)
						{
							ArrListItem_u Item;
							GradPath_pa VolGradPath = NULL;
							Item = ArrListGetItem(GradPathList, ii);
							VolGradPath = (GradPath_pa)Item.VoidPtr;
							GradPathDealloc(&VolGradPath);
						}
						ArrListDealloc(&GradPathList);


						// GradPathDealloc(&SurfGradPath);
					}

					// Inform Tecplot that major data operation is ending
					TecUtilLockFinish(AddOnID);
					TecUtilDataLoadEnd();
				#endif



				ArrListDealloc(&BondForMinList);
				if (RingForSaddleList != NULL) ArrListDealloc(&RingForSaddleList);
				if (RingForMaxList != NULL) ArrListDealloc(&RingForMaxList);
			}


			/* Clean up allocated structures/arrays */
			// CritPointsDealloc(&ITZCritPoints);  Part of IsoSurfGradPath
			SurfTopoSegDealloc(&SurfTopoSeg);
			// BundlesDealloc(&SurfBundles);   Part of IsoSurfGradPath
			IsoSurfGradPathDealloc(&IsoSurfGradPath);

			/* Update status bar */
			if (ShowStatusBar)
			{
				char PercentDoneText[200];
				char AtomBar[] = "[                    ]";
				memset(AtomBar, '|', 1 + (int)(20*(AtomNum-1)/NumAtoms));
				memset(AtomBar, '[', 1);
				if (CompletionLevel > CompletionLevel_AtomIsoSurfTopo)
					sprintf_s(PercentDoneText, "Finding Atom-Ring-Cage Surfaces: Atom %d of %d - 0 Percent Done  %s", 
					AtomNum, NumAtoms, AtomBar);
				else
					sprintf_s(PercentDoneText, "Finding Atom Isosurface Topology: Atom %d of %d - 0 Percent Done  %s", 
					AtomNum, NumAtoms, AtomBar);
				TecUtilStatusSuspend(FALSE);
				TecUtilStatusSetPercentDoneText(PercentDoneText);
				if (!TecUtilStatusCheckPercentDone(0)) IsStop = TRUE;
				TecUtilStatusSuspend(TRUE);
				if (IsStop) IsOk = FALSE;
			}
			IsOk = TRUE;
		}
		if (SurfZoneVarInfo != NULL)
			ZoneVarInfoDealloc(&SurfZoneVarInfo);
			// FREE_ITEM(SurfZoneVarInfo, "SurfZoneVarInfo");
		if (VolZoneVarInfo != NULL)
			ZoneVarInfoDealloc(&VolZoneVarInfo);
			// FREE_ITEM(VolZoneVarInfo, "VolZoneVarInfo");

		if (ShowStatusBar)
		{
			TecUtilStatusSuspend(FALSE);
			TecUtilStatusFinishPercentDone();
			TecUtilStatusSuspend(TRUE);
		}
	}


	/*
	 * OLD METHOD.
	 * To get Atom-Cage lines, seed pathlines around a small sphere
	 * surrounding a atom critical point. Look for minimum length
	 * pathline between the atom and the cage points.
	 */
	// if (IsOk && 1)
	#if(0)  // <<<<<<<<<<<<< TEMPORARY FOR TESTING <<<<<<<<<<<<<<
	// if (IsOk && CompletionLevel >= CompletionLevel_AtomIsoSurfTopo)
	{
		SphereGradPath_pa    SphereGradPath = NULL;
		LgIndex_t            BegOffset = CritPointsGetBegOffset(CritPoints, -3);
		LgIndex_t            EndOffset = CritPointsGetEndOffset(CritPoints, -3);
		LgIndex_t            NumCrtPts = CritPointsGetCount(CritPoints);
		LgIndex_t            M1BegOffset = CritPointsGetBegOffset(CritPoints, -1);
		LgIndex_t            P3BegOffset = CritPointsGetBegOffset(CritPoints,  3);


		ZoneVarInfo_pa ZoneVarInfo = ZoneVarInfoAlloc();
		// ZoneVarInfo_pa ZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "ZoneVarInfo");
		ZoneVarInfo->ZoneNum = ZoneNum;
		ZoneVarInfo->UVarNum = UVarNum;
		ZoneVarInfo->VVarNum = VVarNum;
		ZoneVarInfo->WVarNum = WVarNum;
		ZoneVarInfo->ChrgDensVarNum = ChrgDensVarNum;
		ZoneVarInfo->TypeVarNum = TypeVarNum;
		ZoneVarInfo->PeriodicBC = PeriodicBC;

		IsOk = ZoneVarInfoSetMinMax(ZoneVarInfo, refinedGridSet);

		SphereGradPath = SphereGradPathAlloc();

		/* Set-up status bar */
		if (ShowStatusBar)
			TecUtilStatusStartPercentDone("Finding Atom-Cage lines, surfaces", TRUE, TRUE);

		// for (ii = BegOffset; IsOk && ii < EndOffset; ii++)
		for (ii=1; IsOk && ii<=2; ii++) //tmp
		{
			ArrList_pa CageList;

			/* Inform Tecplot that major data operation is beginning */
			TecUtilDataLoadBegin();
			TecUtilLockStart(AddOnID);

			/* Get the list of cages neighboring this atom, from the bundles */
			CageList = CreateAtomBondCageList(ii, 1, Bundles, CritPoints);  // TEMP: Bond 1 only

			/* test */
			if (0)
			{
				int jj;
				for (jj = 0; jj < ArrListGetCount(CageList); jj++)
				{
					ArrListItem_u Item;
					LgIndex_t CageNum;

					Item = ArrListGetItem(CageList, jj);
					CageNum = Item.Long;
					if (CageNum >= 0)
					{
						char Type;
						GradPath_pa GradPathTst = GradPathAlloc();

						Type = 3;
						GradPathTst->BeginCrtPtNum = ii;
						GradPathTst->EndCrtPtNum = CageNum +
												   CritPointsGetBegOffset(CritPoints, Type);

						IsOk = GradPathAdd2PtMinLen(ZoneVarInfo, brSurfaceIMax, CritPoints, GradPathTst);

						// Create the Atom/Cage line
						if (GradPathIsValid(GradPathTst))
						{
							EntIndex_t ConnectorZoneNum = 0;

							ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, GradPathTst);

							if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
						}
						else
							IsOk = FALSE;
					}
				}
			}

			/* Seed points on a segment of an intermdiate surface between the Atom and Cage.
			 * The segment is defined by the intersection with the surface of Atom-Ring and
			 * Atom-Bond lines that are part of the Atom/Cage bundle.
			 * Use least-square quadratic surface fit and iterate to get min-length.
			 */
			if (IsOk)
			{
				int jj;
				for (jj = 0; jj < ArrListGetCount(CageList); jj++)
				{
					ArrListItem_u Item;
					LgIndex_t CageNum;

					Item = ArrListGetItem(CageList, jj);
					CageNum = Item.Long;
					if (CageNum >= 0)
					{
						LgIndex_t EndCrtPtNum = 0;
						MidPlnGradPath_pa MidPlnGradPath = NULL;

						MidPlnGradPath = MidPlnGradPathAlloc();
						if (MidPlnGradPath == NULL) IsOk = FALSE;
						// old - delete
						//ArrList_pa  RingBondList = NULL;
						//RingBondList = CreateRingBondListForAtomCage(ii, CageNum, Bundles, CritPoints);

						// Define surface as normal to the straight line between the ring and the cage, and
						// half way between them.
						//
						// Find intersections of the Atom-Ring and Atom-Cage lines with the surface.
						//
						// Find the average of the intersection positions
						//
						// Initial set of seed points for gradpaths: the average intersection, and half-way
						// between the average and the Atom-Ring and Atom-Bond intesections.
						EndCrtPtNum = CageNum + CritPointsGetBegOffset(CritPoints, 3) - 1;
						IsOk = MidPlnGradPathAdd(ZoneVarInfo, ii, EndCrtPtNum, CritPoints, Bundles, 1.0, MidPlnGradPath);
						if (IsOk)
						{
							// LgIndex_t NumMinLenGP = MidPlnGradPathMinimizeLen(ZoneVarInfo, ii, EndCrtPtNum,
							//                                                  CritPoints, Bundles, 1.0, MidPlnGradPath);
						}

						// Iterate to find a progressively better approximation of the minimum-length gradpath

						// Create the Atom/Cage line
						if (IsOk)
						{
							GradPath_pa AtomCagePath = MidPlnGradPathGetGP(MidPlnGradPath, 0);

							if (AtomCagePath->BeginCrtPtNum != ii || AtomCagePath->EndCrtPtNum != EndCrtPtNum)
							{
								ENSURE(AtomCagePath->BeginCrtPtNum == ii && AtomCagePath->EndCrtPtNum == EndCrtPtNum);
							}

							if (GradPathIsValid(AtomCagePath))
							{
								EntIndex_t ConnectorZoneNum = 0;

								ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, AtomCagePath);

								if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
							}
							else
								IsOk = FALSE;
						}
					}
				}
			}

			/* Integrate to compute path lines */
			if (0 && IsOk)
			{
				SphereGradPathClear(SphereGradPath);


				/* Integrate gradient path lines from seeded circle */
				IsOk = SphereGradPathAdd(ZoneVarInfo, ii, CritPoints,
										 StreamDir_Reverse, 1.0, SphereGradPath);

				if (IsOk == FALSE)
				{
					char Message[200];
					sprintf_s(Message, "Error computing SphereGradPath number %d", ii);
					TecUtilDialogMessageBox(Message, MessageBox_Warning);
				}
			}

			/* Refine SphereGradPath triangles to find minimum pathlength connection to each cage */
			if (0 && IsOk)
			{
				int jj;
				for (jj = 0; jj < ArrListGetCount(CageList); jj++)
				{
					ArrListItem_u Item;
					LgIndex_t CageNum;
					GradPath_pa AtomCagePath = NULL;

					Item = ArrListGetItem(CageList, jj);
					CageNum = Item.Long;

					// Find min PathLength using shooting method
					if (CageNum >= 0)
					{
						LgIndex_t CageCPNum = CageNum + CritPointsGetBegOffset(CritPoints, 3);

						AtomCagePath = SphereGradPath2CPPath(ZoneVarInfo, ii, CageCPNum, CritPoints,
															 StreamDir_Reverse, 1.0, SphereGradPath);

						// Create the Atom/Cage line
						if (GradPathIsValid(AtomCagePath))
						{
							EntIndex_t ConnectorZoneNum = 0;

							ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, AtomCagePath);

							if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
						}
						else
							IsOk = FALSE;
					}
				}
			}

			/* Decimate large GradPaths to conserve memory */
			// IsOk = SphereGradPathDecimate(SphereGradPath, MaxSize);
			/*
					  if (IsOk)
						{
						  LgIndex_t nc;
						  LgIndex_t NumConnectors = 0;

						  // Extract the Ring-Atom lines from the CircleGradPath
						  NumConnectors = CircleGradPathGetCPNCount(CircleGradPath);

						  for (nc=0; nc<NumConnectors; nc++)
							{
							  LgIndex_t ConnectPathNum = CircleGradPathGetCPN(CircleGradPath, nc);

							  GradPath_pa ConnectPath  = CircleGradPathGetGP(CircleGradPath, ConnectPathNum);
							  LgIndex_t   EndCrtPtNum  = ConnectPath->EndCrtPtNum;
							  LgIndex_t   BegCrtPtNum  = ConnectPath->BeginCrtPtNum;
							  char        EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
							  char        BegCrtPtType = CritPointsGetType(CritPoints, BegCrtPtNum);
							  LgIndex_t   NumPathPoints = GradPathGetCount(ConnectPath);

							  char        ConnectorType = EndCrtPtType + BegCrtPtType;
							  if (ConnectorType == 0 && EndCrtPtType ==  1) ConnectorType =  1;
							  if (ConnectorType == 0 && EndCrtPtType == -1) ConnectorType = -1;

							  // Add a Tecplot zone for the Ring-Atom connector
							  if (IsOk && NumPathPoints > 0 && EndCrtPtNum >= 0 &&
										   EndCrtPtNum < CritPointsGetCount(CritPoints) )
								{
								  char       ZoneName[200];
								  LgIndex_t  P1CrtPtNum = ii - BegOffset;
								  EntIndex_t ConnectorZoneNum = 0;
								  LgIndex_t  M1CrtPtNum = EndCrtPtNum - M1BegOffset;
								  LgIndex_t  M3CrtPtNum = EndCrtPtNum - M3BegOffset;

								  if (ConnectorType == -1)
									sprintf_s(ZoneName, "Ring_%d_Bond_%d_Zone_%d", P1CrtPtNum, M1CrtPtNum, ZoneNum);
								  else if (ConnectorType == -2)
									sprintf_s(ZoneName, "Ring_%d_Atom_%d_Zone_%d", P1CrtPtNum, M3CrtPtNum, ZoneNum);

								  // Only do the Ring-Atom lines, Ring-Bond was already done.
								  if (ConnectorType == -2)
									{
									  ConnectorZoneNum = CreateConnectorZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints, ConnectPath);

									  if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
									}
								}
							}

						  // The edges of the surface are made up of Ring-Atom, Ring-Bond, and
						  // Bond-Atom lines. To represent this "triangle" as an ordered grid
						  // we must append the Ring-Bond and Bond-Atom lines to create a
						  // single Ring-Bond-Atom line that represents either the J=1 or
						  // J=JMax line of the IJ-ordered Tecplot zone.
						  //
						  // In preparation for the creation of the Ring-Bond-Atom surface,
						  // resample the gradient path lines in CircleGradPath to make them
						  // all have an equal number of points.
						  //
						  // IsOk = CircleGradPathResample(CircleGradPath, brSurfaceIMax);

						  NumConnectors = CircleGradPathGetCPNCount(CircleGradPath);

						  for (nc=0; nc<NumConnectors; nc++)
						  // for (nc=0; nc<2; nc++)
							{
							  LgIndex_t  CGPNumBeg, CGPNumEnd;
							  LgIndex_t  ncp1 = nc + 1;
							  EntIndex_t SurfaceZoneNum = 0;

							  if (nc == NumConnectors - 1) ncp1 = ncp1 - NumConnectors;

							  CGPNumBeg = CircleGradPathGetCPN(CircleGradPath, nc);
							  CGPNumEnd = CircleGradPathGetCPN(CircleGradPath, ncp1);
							  if (CGPNumEnd <= CGPNumBeg) CGPNumEnd = CGPNumEnd + CircleGradPathGetCount(CircleGradPath);

							  SurfaceZoneNum = CreateBRSurfaceZone(ZoneNum, ChrgDensVarNum, TypeVarNum, CritPoints,
								CircleGradPath, brSurfaceIMax, CGPNumBeg, CGPNumEnd);
							}

						}
					  */

			/* Update status bar */
			if (ShowStatusBar)
			{
				char PercentDoneText[200];
				int  PercentDone;
				PercentDone = (int)((100 * (ii - BegOffset)) / (EndOffset - BegOffset));
				sprintf_s(PercentDoneText, "Finding Atom-Cage lines, surfaces: %d Percent Done", PercentDone);
				TecUtilStatusSetPercentDoneText(PercentDoneText);
				if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
				if (IsStop) IsOk = FALSE;
			}

			/* Deallocate the CageList for Atom ii */
			ArrListDealloc(&CageList);

			/* Inform Tecplot that major data operation is ending */
			TecUtilLockFinish(AddOnID);
			TecUtilDataLoadEnd();
		}
		if ( ZoneVarInfo != NULL )
			FREE_ITEM(ZoneVarInfo, "ZoneVarInfo");
		if ( SphereGradPath != NULL )
			SphereGradPathDealloc(&SphereGradPath);

		if (ShowStatusBar)
			TecUtilStatusFinishPercentDone();
	}
	#endif


	/*
	 * Find the BondBundle volumes.
	 * Create an arraylist of pointers to BondBundles. Cycle through zones to fin
	 */
	// if (0)       //<<<<<<<<<<<<<<< TEMPORARY FOR TESTING <<<<<<<<<
	if (IsOk && CompletionLevel >= CompletionLevel_BondBundle)
	{
		LgIndex_t      NumBonds =  CritPointsGetEndOffset(CritPoints, (char)(-1)) -
								   CritPointsGetBegOffset(CritPoints, (char)(-1));

		BondBundles_pa BondBundles = BondBundlesAlloc(NumBonds);
		// ArrList_pa     BondBundlesList = ArrListAlloc(NumBonds, ArrListType_VoidPtr);

		IsOk = BondBundlesGetFromZoneAuxData(ZoneNum, CritPoints, BondBundles);
		// IsOk = BondBundlesGetFromZoneAuxData(ZoneNum, CritPoints, BondBundlesList);
		IsOk = BondBundlesCreateVolumeZones(ZoneNum, CritPoints, BondBundles,
											1.0, (double)IMax,
											1.0, (double)JMax,
											1.0, (double)KMax);
	}


	// Clean up
	if (AtomPairsForBonds != NULL) OrderedPairsDealloc(&AtomPairsForBonds);

	
	// Integrate
	if (IsOk && CompletionLevel >= CompletionLevel_BondBundle)
	{
		char const *auxDataName = "CompChem.BondIntegratedChargeDensity";

		// integrate over source zone
		if (IsOk)
		{
			Set_pa zoneSet = TecUtilSetAlloc(TRUE);
			IsOk = (zoneSet != NULL);
			if (IsOk)
			{
				IsOk = TecUtilSetAddMember(zoneSet, ZoneNum, TRUE);
				if (IsOk)
					IsOk = integrateScalarOverVolumeZones(ChrgDensVarNum, auxDataName, zoneSet);
				TecUtilSetDealloc(&zoneSet);
			}
		}

		// integrate over new zones
		if (IsOk)
		{
			// get current zone set and remove all original zones
			Set_pa zoneSet = NULL;
			IsOk = TecUtilZoneGetEnabled(&zoneSet);
			if (IsOk)
			{
				SetIndex_t zone;
				TecUtilSetForEachMember(zone, originalZoneSet)
					TecUtilSetRemoveMember(zoneSet, zone);

				IsOk = integrateScalarOverVolumeZones(ChrgDensVarNum, auxDataName, zoneSet);

				TecUtilSetDealloc(&zoneSet);
			}
		}
	}

	// Clean up from integration
	if (originalZoneSet != NULL)
		TecUtilSetDealloc(&originalZoneSet);

	TecUtilPleaseWait(NULL, FALSE);

	return(IsOk);
}


void CalcGradGradMag(Boolean_t IsPeriodic){

	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	int numCPU = sysinfo.dwNumberOfProcessors;
	omp_set_num_threads(numCPU);

	vector<string> GradXYZMagStr = {
		"X Density Gradient",
		"Y Density Gradient",
		"Z Density Gradient",
		"Density Gradient Magnitude"
	};
	vector<EntIndex_t> GradXYZMagVarNum(GradXYZMagStr.size(), -1);

	TecUtilLockStart(AddOnID);

	EntIndex_t RhoVarNum = VarNumByNameList(vector<string>({
		"Electron Density",
		"Rho",
		"rho"
	}));

	for (int i = 0; i < GradXYZMagStr.size(); ++i)
		GradXYZMagVarNum[i] = VarNumByName(GradXYZMagStr[i]);

	Boolean_t HasGrad = (GradXYZMagVarNum[0] > 0 && GradXYZMagVarNum[1] > 0 && GradXYZMagVarNum[2] > 0);
	Boolean_t HasGradMag = GradXYZMagVarNum[3] > 0;

	if (HasGrad && HasGradMag){
		TecUtilLockFinish(AddOnID);
		return;
	}

	if (RhoVarNum <= 0 && (!HasGrad || !HasGradMag)){
		TecUtilDialogErrMsg("Couldn't find electron density variable, so did not compute gradient or gradient magnitude!");
		TecUtilLockFinish(AddOnID);
		return;
	}

	EntIndex_t ZoneNum = ZoneNumByName(string("Full Volume"));
	vector<LgIndex_t> IJKMax(3, -1);
	if (ZoneNum > 0){
		TecUtilZoneGetIJK(ZoneNum, &IJKMax[0], &IJKMax[1], &IJKMax[2]);
	}
	if (ZoneNum <= 0 || IJKMax[0] <= 0 || IJKMax[1] <= 0 || IJKMax[2] <= 0){
		TecUtilDialogErrMsg("Failed to get zone information. Didn't calculation gradient or gradient magnitude.");
		TecUtilLockFinish(AddOnID);
		return;
	}

	vector<EntIndex_t> XYZNum(3, -1);
	TecUtilAxisGetVarAssignments(&XYZNum[0], &XYZNum[1], &XYZNum[2]);
	if (XYZNum[0] <= 0 || XYZNum[1] <= 0 || XYZNum[2] <= 0){
		vector<vector<string> > TmpStrs = {
			{ "X", "Y", "Z" }, { "I", "J", "K" }
		};
		Boolean_t VarsFound = FALSE;
		for (int i = 0; i < 2 && !VarsFound; ++i){
			for (int j = 0; j < 3; ++j)
				XYZNum[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
			VarsFound = (XYZNum[0] > 0 && XYZNum[1] > 0 && XYZNum[2] > 0);
		}
		if (!VarsFound){
			TecUtilDialogErrMsg("Failed to get x, y, z variable numbers. Didn't calculate gradient or gradient magnitude");
			TecUtilLockFinish(AddOnID);
			return;
		}
	}
	FieldDataType_e FDTJunk;

	vector<float*> XYZRawPtr(3, NULL);
	for (int i = 0; i < 3; ++i){
		TecUtilDataValueGetReadableRawPtr(ZoneNum, XYZNum[i], (void**)&XYZRawPtr[i], &FDTJunk);
	}

	vector<float*> GradXYZMagRawPtr(GradXYZMagVarNum.size(), NULL);

	float* RhoRawPtr;
	TecUtilDataValueGetReadableRawPtr(ZoneNum, RhoVarNum, (void**)&RhoRawPtr, &FDTJunk);

	EntIndex_t NumZones = TecUtilDataSetGetNumZones();
	vector<FieldDataType_e> DataType(NumZones, FieldDataType_Float);
	vector<ValueLocation_e> DataLoc(NumZones, ValueLocation_Nodal);

	int StatusKMax = IJKMax[2] / numCPU;

	bool TaskQuit = false;

	TecUtilDataLoadBegin();

	if (!HasGrad){
		vector<FieldValueSetFunction_pf> GradXYZPF(3, NULL);
		for (int i = 0; i < 3; ++i){
			if (GradXYZMagVarNum[i] <= 0){
				ArgList_pa Args = TecUtilArgListAlloc();
				TecUtilArgListAppendString(Args, SV_NAME, GradXYZMagStr[i].c_str());
				TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
				TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
				TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

				if (TecUtilDataSetAddVarX(Args)){
					Set_pa NewVar = TecUtilSetAlloc(FALSE);
					GradXYZMagVarNum[i] = TecUtilDataSetGetNumVars();
					TecUtilSetAddMember(NewVar, GradXYZMagVarNum[i], FALSE);
					TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)NewVar);
					TecUtilSetDealloc(&NewVar);
				}
				TecUtilArgListDealloc(&Args);
			}

			TecUtilDataValueGetRawPtr(ZoneNum, GradXYZMagVarNum[i], (void**)&GradXYZMagRawPtr[i], &FDTJunk);
		}

		TecUtilDialogLaunchPercentDone("Calculating gradient vector", TRUE);


#pragma omp parallel for
		for (LgIndex_t kk = 1; kk <= IJKMax[2]; ++kk){
			if (omp_get_thread_num() == 0 && !SetPercent(kk, StatusKMax, "Calculating gradient vector")){
				TaskQuit = true;
#pragma omp flush (TaskQuit)
			}
#pragma omp flush (TaskQuit)
			for (LgIndex_t jj = 1; jj <= IJKMax[1] && !TaskQuit; ++jj){
				for (LgIndex_t ii = 1; ii <= IJKMax[0]; ++ii){
					for (int Dir = 0; Dir < 3; ++Dir){
						LgIndex_t PM1[3], PP1[3], Pt[3];
						Pt[0] = PM1[0] = PP1[0] = ii;
						Pt[1] = PM1[1] = PP1[1] = jj;
						Pt[2] = PM1[2] = PP1[2] = kk;

						PM1[Dir]--;
						PP1[Dir]++;

						if (!IsPeriodic){
							LgIndex_t IJK[3] = { ii, jj, kk };
							for (int i = 0; i < 3; ++i){
								if (IJK[i] == 1)
									PM1[i] = Pt[i];
								else if (IJK[i] == IJKMax[i])
									PP1[i] = Pt[i];
							}
						}

						LgIndex_t Point = IndexFromIJK(Pt[0], Pt[1], Pt[2], IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1;
						LgIndex_t PtMinus = IndexFromIJK(PM1[0], PM1[1], PM1[2], IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1;
						LgIndex_t  PtPlus = IndexFromIJK(PP1[0], PP1[1], PP1[2], IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1;


						float PlusCoord = XYZRawPtr[Dir][PtPlus];
						float MinusCoord = XYZRawPtr[Dir][PtMinus];

						float Dist = PlusCoord - MinusCoord;

						float PlusRho = RhoRawPtr[PtPlus];
						float MinusRho = RhoRawPtr[PtMinus];
						float DRho = (PlusRho - MinusRho) / Dist;
						GradXYZMagRawPtr[Dir][Point] = (float)DRho;
					}
				}
			}
		}

		if (TaskQuit){
			TecUtilDialogDropPercentDone();
			TecUtilDataLoadEnd();
			Set_pa Vars = TecUtilSetAlloc(FALSE);
			for (int i = 0; i < 3; ++i)
				TecUtilSetAddMember(Vars, GradXYZMagVarNum[i], FALSE);
			TecUtilDataSetDeleteVar(Vars);
			TecUtilSetDealloc(&Vars);
			TecUtilLockFinish(AddOnID);
			return;
		}

		TecUtilDialogDropPercentDone();
	}

	if (!HasGradMag){
		vector<FieldValueGetFunction_pf> GradXYZPF(3, NULL);
		float* MagRawPtr;
		FieldData_pa MagPtr = NULL;
		FieldValueSetFunction_pf MagPF = NULL;
		for (int i = 0; i < 4; ++i){
			if (GradXYZMagVarNum[i] <= 0){
				ArgList_pa Args = TecUtilArgListAlloc();
				TecUtilArgListAppendString(Args, SV_NAME, GradXYZMagStr[i].c_str());
				TecUtilArgListAppendArray(Args, SV_VARDATATYPE, DataType.data());
				TecUtilArgListAppendArray(Args, SV_VALUELOCATION, DataLoc.data());
				TecUtilArgListAppendInt(Args, SV_SHAREVARWITHALLZONES, FALSE);

				if (TecUtilDataSetAddVarX(Args)){
					Set_pa NewVar = TecUtilSetAlloc(FALSE);
					GradXYZMagVarNum[i] = TecUtilDataSetGetNumVars();
					TecUtilSetAddMember(NewVar, GradXYZMagVarNum[i], FALSE);
					TecUtilStateChanged(StateChange_VarsAdded, (ArbParam_t)NewVar);
					TecUtilSetDealloc(&NewVar);
				}
				TecUtilArgListDealloc(&Args);
			}

			if (i < 3){
				TecUtilDataValueGetReadableRawPtr(ZoneNum, GradXYZMagVarNum[i], (void**)&GradXYZMagRawPtr[i], &FDTJunk);
			}
			else{
				TecUtilDataValueGetRawPtr(ZoneNum, GradXYZMagVarNum[i], (void**)&MagRawPtr, &FDTJunk);
			}
		}

		TecUtilDialogLaunchPercentDone("Calculating gradient magnitude", TRUE);

#pragma omp parallel for
		for (LgIndex_t kk = 1; kk <= IJKMax[2]; ++kk){
			if (omp_get_thread_num() == 0 && !SetPercent(kk, StatusKMax, "Calculating gradient magnitude")){
				TaskQuit = true;
#pragma omp flush (TaskQuit)
			}
#pragma omp flush (TaskQuit)
			for (LgIndex_t jj = 1; jj <= IJKMax[1] && !TaskQuit; ++jj){
				for (LgIndex_t ii = 1; ii <= IJKMax[0]; ++ii){

					LgIndex_t Point = IndexFromIJK(ii, jj, kk, IJKMax[0], IJKMax[1], IJKMax[2], IsPeriodic) - 1;

					float Magnitude = 0.0;
					for (int i = 0; i < 3; ++i)
						Magnitude += pow(GradXYZMagRawPtr[i][Point], 2);

					Magnitude = sqrt(Magnitude);

					MagRawPtr[Point] = Magnitude;
				}
			}
		}

		if (TaskQuit){
			TecUtilDialogDropPercentDone();
			TecUtilDataLoadEnd();
			Set_pa Vars = TecUtilSetAlloc(FALSE);
			TecUtilSetAddMember(Vars, GradXYZMagVarNum[3], FALSE);
			TecUtilDataSetDeleteVar(Vars);
			TecUtilSetDealloc(&Vars);
			TecUtilLockFinish(AddOnID);
			return;
		}

		TecUtilDialogDropPercentDone();
	}

	TecUtilDataLoadEnd();
	TecUtilLockFinish(AddOnID);
}

Boolean_t SetPercent(unsigned int CurrentNum, unsigned int TotalNum, const char* ProgresssText){
	unsigned int Percent = MIN(int((double)CurrentNum / (double)TotalNum * 100.), 100);

	std::stringstream ss;
	ss << ProgresssText << "  (" << Percent << "% Complete)";

	TecUtilLockStart(AddOnID);

	Boolean_t IsOk = TecUtilDialogCheckPercentDone(Percent);
	TecUtilDialogSetPercentDoneText(ss.str().c_str());

	TecUtilLockFinish(AddOnID);

	return IsOk;
}

