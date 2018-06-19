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
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "ODERUNGEKUTTA.h"
#include "NORMALS.h"
#include "SURFELEMMAP.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
#include "GRADPATH.h"
#include "BUNDLES.h"
#include "GEOMTOOLS.h"
#include "SURFTOPOSEG.h"
#include "ISOSURFGRADPATH.h"
#include "SURFACEFIT.h"
#include "ENGINE.h"



/**
 * Determine if the IsoSurfGradPath handle is sane.
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath structure in question.
 *
 * return
 *     TRUE if the IsoSurfGradPath structure is valid, otherwise FALSE.
 */
Boolean_t IsoSurfGradPathIsValid(IsoSurfGradPath_pa IsoSurfGradPath)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(IsoSurfGradPath) &&
               // VALID_REF(IsoSurfGradPath->VolBundles)     && BundlesIsValid(IsoSurfGradPath->VolBundles) &&
               VALID_REF(IsoSurfGradPath->SurfCritPoints) && CritPointsIsValid(IsoSurfGradPath->SurfCritPoints) &&
               VALID_REF(IsoSurfGradPath->SurfBundles)    && BundlesIsValid(IsoSurfGradPath->SurfBundles) &&
               VALID_REF(IsoSurfGradPath->GradPaths)      && ArrListIsValid(IsoSurfGradPath->GradPaths));


    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}





/**
 * Gets the Index'th GradPath handle.
 *
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath to get the GradPath from.
 *
 * param Offset
 *     Offset into the GradPath pointer ArrList.
 *
 * return
 *     GradPath handle available, otherewise a handle to NULL.
 */
GradPath_pa IsoSurfGradPathGetGP(IsoSurfGradPath_pa IsoSurfGradPath,
                                 LgIndex_t         Offset)
{
    GradPath_pa   Result = NULL;
    ArrListItem_u Item;

    REQUIRE(VALID_REF(IsoSurfGradPath));
    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));
    REQUIRE(Offset >= 0 && Offset < IsoSurfGradPathGetGPCount(IsoSurfGradPath));

    Item = ArrListGetItem(IsoSurfGradPath->GradPaths, Offset);
    Result = (GradPath_pa)Item.VoidPtr;

    ENSURE(VALID_REF(Result));
    return Result;
}









/**
 * Deallocates the IsoSurfGradPath handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a IsoSurfGradPath handle.
 */
void IsoSurfGradPathDealloc(IsoSurfGradPath_pa *IsoSurfGradPath)
{

    REQUIRE(VALID_REF(IsoSurfGradPath));
    // REQUIRE(IsoSurfGradPathIsValid(*IsoSurfGradPath) || *IsoSurfGradPath == NULL);

    if (*IsoSurfGradPath != NULL)
    {
        LgIndex_t ii;
        LgIndex_t Count = ArrListGetCount((*IsoSurfGradPath)->GradPaths);

        /* Dealloc the GradPaths */
        for (ii = 0; ii < Count; ii++)
        {
            GradPath_pa GradPath = IsoSurfGradPathGetGP((*IsoSurfGradPath), ii);
            GradPathDealloc(&GradPath);
        }


        /* release the ArrList's */
        if ((*IsoSurfGradPath)->GradPaths != NULL) ArrListDealloc(&((*IsoSurfGradPath)->GradPaths));
        if ((*IsoSurfGradPath)->SurfBundles != NULL) BundlesDealloc(&((*IsoSurfGradPath)->SurfBundles));
        if ((*IsoSurfGradPath)->SurfCritPoints != NULL) CritPointsDealloc(&((*IsoSurfGradPath)->SurfCritPoints));
        // TODO: Do we need to delete the IsoTopoZone?

        /* release the list structure itself */
        FREE_ITEM(*IsoSurfGradPath, "IsoSurfGradPath structure");
        *IsoSurfGradPath = NULL;
    }

    ENSURE(*IsoSurfGradPath == NULL);
}




/**
 * Gets the number of GradPaths currently in the IsoSurfGradPath
 * (maintained by the IsoSurfGradPath array lists).
 *
 * param
 *     IsoSurfGradPath structure in question.
 *
 * return
 *     Number of GradPaths/RelNormSeedVector's in the IsoSurfGradPath.
 */
LgIndex_t IsoSurfGradPathGetGPCount(IsoSurfGradPath_pa IsoSurfGradPath)
{
    LgIndex_t Result = 0;

    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));

    Result = ArrListGetCount(IsoSurfGradPath->GradPaths);

    ENSURE(Result >= 0);
    return Result;
}











/**
 * Empties the IsoSurfGradPath of all Seeds, GradPaths, and BBPIntersects
 * and resets the other variables.
 *
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath to clear.
 */
void IsoSurfGradPathClear(IsoSurfGradPath_pa IsoSurfGradPath)
{
    LgIndex_t ii;
    LgIndex_t Count = ArrListGetCount(IsoSurfGradPath->GradPaths);

    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));

    for (ii = 0; ii < Count; ii++)
    {
        GradPath_pa GradPath = IsoSurfGradPathGetGP(IsoSurfGradPath, ii);
        GradPathDealloc(&GradPath);
    }

    ArrListClear(IsoSurfGradPath->GradPaths);
    BundlesClear(IsoSurfGradPath->SurfBundles);
    CritPointsClear(IsoSurfGradPath->SurfCritPoints);

    IsoSurfGradPath->AtomNum        = -1;
    IsoSurfGradPath->VolBundles     = NULL;
    IsoSurfGradPath->IsoTopoZone    = 0;
    IsoSurfGradPath->SurfCPZoneNum  = 0;

    ENSURE(IsoSurfGradPathIsValid(IsoSurfGradPath) &&
           IsoSurfGradPathGetGPCount(IsoSurfGradPath) == 0 &&
           BundlesGetCount(IsoSurfGradPath->SurfBundles) == 0 &&
           CritPointsGetCount(IsoSurfGradPath->SurfCritPoints) == 0);
}






/**
 * Allocates a IsoSurfGradPath handle with a suitable default
 * capacity for the contained ArrList's.
 *
 *
 * return
 *     IsoSurfGradPath handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
IsoSurfGradPath_pa IsoSurfGradPathAlloc()
{
    IsoSurfGradPath_pa Result = NULL;

    Result = ALLOC_ITEM(IsoSurfGradPath_s, "IsoSurfGradPath structure");
    if (Result != NULL)
    {
        Result->AtomNum        = -1;
        Result->VolBundles     = NULL;  // Don't dealloc - this pointer is an input
        Result->IsoTopoZone    = 0;     // Don't delete zone - this input is a pointer
        Result->SurfCPZoneNum  = 0;     // TODO: should we dealloc this?
        Result->SurfCritPoints = CritPointsAlloc();
        Result->SurfBundles    = BundlesAlloc();

        Result->GradPaths      = ArrListAlloc(20, ArrListType_VoidPtr);

        /* If it failed to allocate any of the structures, clean-up and exit. */
        if (Result->SurfCritPoints == NULL || !CritPointsIsValid(Result->SurfCritPoints) ||
            Result->SurfBundles    == NULL || !BundlesIsValid(Result->SurfBundles) ||
            Result->GradPaths      == NULL || !ArrListIsValid(Result->GradPaths))
        {
            if (Result->SurfCritPoints != NULL) CritPointsDealloc(&(Result->SurfCritPoints));
            if (Result->SurfBundles    != NULL) BundlesDealloc(&(Result->SurfBundles));
            if (Result->GradPaths      != NULL) ArrListDealloc(&(Result->GradPaths));
            FREE_ITEM(Result, "IsoSurfGradPath structure");
            Result = NULL;
        }

    }

    ENSURE(IsoSurfGradPathIsValid(Result) || Result == NULL);
    return Result;
}






/**
 * Places GradPath pointer at the specified offset. If the
 * offset is beyond the end of the array lists, they are sized
 * accordingly and the intervening array values between the last
 * node of the original state and the last node of the new state are
 * assigned 0.
 *
 * note
 *     If GradPaths already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath target in which to set the Angle/GradPath.
 * param Offset
 *     Offset into array lists.
 * param GradPath
 *     GradPath handle to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t   IsoSurfGradPathSetGP(IsoSurfGradPath_pa IsoSurfGradPath,
                                 LgIndex_t   Offset,
                                 GradPath_pa GradPath)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));
    REQUIRE(Offset >= 0);

    if (IsOk)
    {
        Item.VoidPtr = (void *)GradPath;
        IsOk = ArrListSetItem(IsoSurfGradPath->GradPaths, Offset, Item);
    }

    ENSURE(IsoSurfGradPathIsValid(IsoSurfGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Inserts GradPath at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 * note
 *     If GradPath already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath target in which to set the coordinates.
 * param Offset
 *     Offset at which to insert the Angle/GradPath.
 * param GradPath
 *     GradPath handle to insert at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t IsoSurfGradPathInsertGP(IsoSurfGradPath_pa IsoSurfGradPath,
                                  LgIndex_t         Offset,
                                  GradPath_pa       GradPath)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));
    REQUIRE(0 <= Offset && Offset <= IsoSurfGradPathGetGPCount(IsoSurfGradPath));

    if (IsOk)
    {
        Item.VoidPtr = (void *)GradPath;
        IsOk = ArrListInsertItem(IsoSurfGradPath->GradPaths, Offset, Item);
    }

    ENSURE(IsoSurfGradPathIsValid(IsoSurfGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Deletes GradPath at the specified offset.
 *
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath target in which to set the coordinates.
 * param Offset
 *     Offset at which to insert the SeedVector/GradPath.
 *
 * return
 *     TRUE if successful operation, otherwise FALSE.
 */

Boolean_t IsoSurfGradPathRemoveGP(IsoSurfGradPath_pa IsoSurfGradPath,
                                  LgIndex_t         Offset)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));
    REQUIRE(0 <= Offset && Offset <= IsoSurfGradPathGetGPCount(IsoSurfGradPath));

    if (IsOk)
    {
        GradPath_pa GradPath;
        Item = ArrListRemoveItem(IsoSurfGradPath->GradPaths, Offset);
        GradPath = (GradPath_pa)Item.VoidPtr;
        GradPathDealloc(&GradPath);
    }

    ENSURE(IsoSurfGradPathIsValid(IsoSurfGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Appends the GradPath to the array lists in IsoSurfGradPath. The
 * array lists will be expanded to accommodate the additional items.
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath target to which the SeedVector/GradPath is to be appended.
 * param SeedPos
 *     Seed point position (on mid-plane) to append to the SeedVector/GradPath.
 * param GradPath
 *     GradPath handle to append to the IsoSurfGradPath.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t IsoSurfGradPathAppendGP(IsoSurfGradPath_pa  IsoSurfGradPath,
                                  GradPath_pa        GradPath)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));

    Count = IsoSurfGradPathGetGPCount(IsoSurfGradPath);

    IsOk = IsoSurfGradPathInsertGP(IsoSurfGradPath, Count, GradPath);

    ENSURE(IsoSurfGradPathIsValid(IsoSurfGradPath));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}










/**
 * Return the ending critical point number for the GradPath at Offset
 * in the array list.
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath target to extract the EndCrtPtNum from.
 * param Offset
 *     Offset into array list of GradPath from which to extract EndCrtPtNum.
 *
 * return
 *     EndCrtPtNum for the path.
 */
LgIndex_t IsoSurfGradPathGetEndCPNum(IsoSurfGradPath_pa  IsoSurfGradPath,
                                     LgIndex_t          Offset)
{
    LgIndex_t EndCrtPtNum = -2;

    GradPath_pa GradPath = NULL;
    ArrListItem_u Item;

    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));
    REQUIRE(0 <= Offset && Offset < IsoSurfGradPathGetGPCount(IsoSurfGradPath));

    Item = ArrListGetItem(IsoSurfGradPath->GradPaths, Offset);
    GradPath = (GradPath_pa)Item.VoidPtr;

    EndCrtPtNum = GradPath->EndCrtPtNum;

    ENSURE(-1 <= EndCrtPtNum);
    return EndCrtPtNum;
}




/**
 * Return the begining critical point number for the IsoSurfGradPath.
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath target to extract the BegCrtPtNum from.
 *
 * return
 *     BegCrtPtNum for the path.
 */
LgIndex_t IsoSurfGradPathGetBegCPNum(IsoSurfGradPath_pa  IsoSurfGradPath)
{
    LgIndex_t BegCrtPtNum = -2;
    LgIndex_t Offset = 0;

    GradPath_pa GradPath = NULL;
    ArrListItem_u Item;

    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));

    Item = ArrListGetItem(IsoSurfGradPath->GradPaths, Offset);
    GradPath = (GradPath_pa)Item.VoidPtr;

    BegCrtPtNum = GradPath->BeginCrtPtNum;

    ENSURE(-1 <= BegCrtPtNum);
    return BegCrtPtNum;
}









/*
 * IsoSurfGradPathGetGPByBegEndCP: Return a pointer to the GradPath with the
 *  given beginning and ending surface critical points.
 *
 * param IsoSurfGradPath
 *     IsoSurfGradPath structure.
 * param BegCPOffset
 *     Offset of the surface critical point at the beginning of the GradPath.
 * param EndCPOffset
 *     Offset of the surface critical point at the end of the GradPath.
 *
 *
 * return
 *     Pointer to GradPath structure if successful, NULL if there were errors.
 */
GradPath_pa  IsoSurfGradPathGetGPByBegEndCP(IsoSurfGradPath_pa IsoSurfGradPath,
                                            LgIndex_t          BegCPOffset,
                                            LgIndex_t          EndCPOffset)
{
    GradPath_pa  Result = NULL;
    LgIndex_t    NumCritPoints;
    LgIndex_t    NumGradPaths, GPOffset;

    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));

    NumCritPoints = CritPointsGetCount(IsoSurfGradPath->SurfCritPoints);
    REQUIRE(BegCPOffset >= 0 && BegCPOffset < NumCritPoints);
    REQUIRE(EndCPOffset >= 0 && EndCPOffset < NumCritPoints);

    // Loop over GradPaths in the IsoSurfGradPath structure, comparing the beginnning
    // and ending CPs with the desired offsets
    NumGradPaths = ArrListGetCount(IsoSurfGradPath->GradPaths);
    for (GPOffset = 0; !Result && GPOffset < NumGradPaths; GPOffset++)
    {
        GradPath_pa Candidate = NULL;
        ArrListItem_u Item;

        Item = ArrListGetItem(IsoSurfGradPath->GradPaths, GPOffset);
        Candidate = (GradPath_pa)(Item.VoidPtr);

        CHECK(GradPathIsValid(Candidate));

        if (Candidate->BeginCrtPtNum == BegCPOffset && Candidate->EndCrtPtNum == EndCPOffset)
            Result = Candidate;
    }


    ENSURE(Result == NULL || GradPathIsValid(Result));
    return Result;
}




/*
 * IsoTopoZoneForAtom: For a given max (atom) volume critical point, compute
 *  the "nearly largest" isosurface of the scalar variable (rho) surrounding the
 *  atom.
 *
 * param ZoneVarInfo
 *     Zone and variable information.
 * param AtomCrtPtNum
 *     Offset of the Atom volume critical point that is surrounded by the isosurface.
 * param VolCritPoints
 *     Critical point structure for the volume critical points.
 * param CPTolerance
 *     Tolerance - distance from critical point at which the pathline is said
 *     to terminate at the critical point.
 *
 * return
 *     IsoTopoZone number (1-based) if successful, -1 if there were errors.
 */
EntIndex_t IsoTopoZoneForAtom(const ZoneVarInfo_pa  ZoneVarInfo,
                              const LgIndex_t       AtomCrtPtNum,
                              const CritPoints_pa   VolCritPoints,
                              const float           CPTolerance)
{
    EntIndex_t     IsoTopoZone = -1;
    EntIndex_t     IsoSurfZone = -1;
    Boolean_t      IsOk = TRUE;
    double         XCrtPtAtom, YCrtPtAtom, ZCrtPtAtom;
    LgIndex_t      FirstElemNumOffset = -1;
    EntIndex_t     VolZoneNum = VolCritPoints->SourceZoneNum;
    SurfTopoSeg_pa SurfTopoSeg = NULL;

    REQUIRE(VALID_REF(ZoneVarInfo));
    REQUIRE(CritPointsIsValid(VolCritPoints));
    REQUIRE(AtomCrtPtNum >= 0 && AtomCrtPtNum < CritPointsGetCount(VolCritPoints));
    REQUIRE(CPTolerance > 0.0);
    // Atom critical points must be an Atom
    REQUIRE(CritPointsGetType(VolCritPoints, AtomCrtPtNum) == -3);


    // Set the isosurface value to just larger than the largest rho of
    // a connected bond critical point
    if (IsOk)
    {
        double MaxBondCPRho = 0.0;

		// Find the connected bond critical points
        if (IsOk)
        {
            LgIndex_t ic;
            LgIndex_t CPBondBeg = CritPointsGetBegOffset(VolCritPoints, (char)(-1));
            LgIndex_t CPBondEnd = CritPointsGetEndOffset(VolCritPoints, (char)(-1));
            for (ic = CPBondBeg; IsOk && ic < CPBondEnd; ic++)
            {
                EntIndex_t TPZoneGP = GradPathTPZoneFromBegEndCP(ic, AtomCrtPtNum, VolZoneNum);
                if (TPZoneGP > 0)
                {
                    double BondCPRho, dummy;
                    char   cdummy;
                    IsOk = CritPointsGetPoint(VolCritPoints, ic, &dummy, &dummy, &dummy,
                                              &BondCPRho, &cdummy, &dummy, &dummy, &dummy);
                    MaxBondCPRho = MAX(MaxBondCPRho, BondCPRho);
                }
            }
        }

		//	TEMP disabled for debugging
		REQUIRE(MaxBondCPRho > 0);
		//if (MaxBondCPRho <= 0) IsOk = FALSE;


        // Set the isosurface value
        if (IsOk)
        {
			//	show isosurfaces and assign rho as variable 1 using 
			//	TecUtilSetLowLevelX() instead.
			ArgList_pa CurrentArgList = TecUtilArgListAlloc();

			//	show isosurfaces
			TecUtilArgListClear(CurrentArgList);
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ISOSURFACELAYERS);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOW);
			TecUtilArgListAppendInt(CurrentArgList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
			SetValueReturnCode_e SVRC = TecUtilStyleSetLowLevelX(CurrentArgList);

			TecUtilArgListClear(CurrentArgList);

			if (!(SVRC == SetValue_Ok && SVRC == SetValueReturnCode_Ok))
			{
				if(SVRC != SetValueReturnCode_DuplicateValue)
				{
					TecUtilDialogMessageBox("Failed to activate isosurfaces",MessageBoxType_Information);
					IsOk = FALSE;
				}
			}
// 			std::vector<std::string> VarNameList(2);
// 			VarNameList[0] = "rho";
// 			VarNameList[1] = "Electron Density";
// 
// 			SmInteger_t RhoVarNum = VarNumByName(VarNameList);
// 			if (RhoVarNum < 0) RhoVarNum = 4;

			//	set isosurface color 1 to variable 1 (rho)
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_DEFINITIONCONTOURGROUP);
			TecUtilArgListAppendInt(CurrentArgList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, (SmInteger_t)4);
			SVRC = TecUtilStyleSetLowLevelX(CurrentArgList);

			TecUtilArgListClear(CurrentArgList);

			if (!(SVRC == SetValue_Ok && SVRC == SetValueReturnCode_Ok))
			{
				if(SVRC != SetValueReturnCode_DuplicateValue)
				{
					TecUtilDialogMessageBox("Failed to set isosurface color",MessageBoxType_Information);
					IsOk = FALSE;
				}
			}

			//	set variable of contour group 1 to rho
			TecUtilArgListAppendInt(CurrentArgList, SV_CONTOURGROUP, 4);
			TecUtilArgListAppendInt(CurrentArgList, SV_VAR, ZoneVarInfo->ChrgDensVarNum);
			SVRC = TecUtilContourSetVariableX(CurrentArgList);

			TecUtilArgListClear(CurrentArgList);

			if (!(SVRC == SetValue_Ok && SVRC == SetValueReturnCode_Ok))
			{
				if(SVRC != SetValueReturnCode_DuplicateValue)
				{
					TecUtilDialogMessageBox("Failed to set contour group variable",MessageBoxType_Information);
					IsOk = FALSE;
				}
			}

			//	set value to create isosurface at to max bond rho
			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_ISOVALUE1);
			TecUtilArgListAppendInt(CurrentArgList, SV_OFFSET1, 1);
			TecUtilArgListAppendDouble(CurrentArgList, SV_DVALUE, (ATOMISOSURFFACTOR * MaxBondCPRho));
			SVRC = TecUtilStyleSetLowLevelX(CurrentArgList);

			TecUtilArgListClear(CurrentArgList);

			if (!(SVRC == SetValue_Ok && SVRC == SetValueReturnCode_Ok))
			{
				if(SVRC != SetValueReturnCode_DuplicateValue)
				{
					TecUtilDialogMessageBox("Failed to set isosurface value",MessageBoxType_Information);
					IsOk = FALSE;
				}
			}

			TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_ISOSURFACEATTRIBUTES);
			TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_SHOWGROUP);
			TecUtilArgListAppendInt(CurrentArgList, SV_OFFSET1, 1);
			TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, TRUE);
			SVRC = TecUtilStyleSetLowLevelX(CurrentArgList);

			if (!(SVRC == SetValue_Ok && SVRC == SetValueReturnCode_Ok))
			{
				if(SVRC != SetValueReturnCode_DuplicateValue)
				{
					TecUtilDialogMessageBox("Failed to set isosurface color group",MessageBoxType_Information);
					IsOk = FALSE;
				}
			}

			TecUtilArgListDealloc(&CurrentArgList);

            //SetValueReturnCode_e SVRC = TecUtilStyleSetLowLevel((Widget)NULL,
            //                                                    // 1.01 * MaxBondCPRho,
            //                                                    // 1.02 * MaxBondCPRho,
            //                                                    ATOMISOSURFFACTOR * MaxBondCPRho,
            //                                                    (ArbParam_t)NULL,
            //                                                    (ArbParam_t)1,
            //                                                    AssignOp_Equals,
            //                                                    SV_ISOSURFACEATTRIBUTES,
            //                                                    SV_ISOVALUE1,
            //                                                    (char *)NULL,
            //                                                    (char *)NULL,
            //                                                    (char *)NULL,
            //                                                    (char *)NULL,
            //                                                    FALSE);

            //if (SVRC != SetValue_Ok && SVRC != SetValue_DuplicateValue) IsOk = FALSE;
        }

    }

    // Extract the isosurface to a zone
    if (IsOk)
        IsOk = TecUtilCreateIsoZones();
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &IsoSurfZone, NULL);

    // Create a zone that just has the topologically independent segment of the
    // isosurface around the atom critical point.
    if (IsOk)
    {
        SurfTopoSeg = SurfTopoSegAlloc();
        if (SurfTopoSeg == NULL) IsOk = FALSE;
    }
    if (IsOk)
        IsOk = SurfTopoSegInitialize(SurfTopoSeg, IsoSurfZone);

    // find Atom (max) critical point location
    if (IsOk)
    {
        double dummy;
        char   cdummy;

        IsOk = CritPointsGetPoint(VolCritPoints, AtomCrtPtNum, &XCrtPtAtom, &YCrtPtAtom, &ZCrtPtAtom,
                                  &dummy, &cdummy, &dummy, &dummy, &dummy);
    }

    // Find the number of the IsoSurfZone element closest to the Atom (Max) critical point
    if (IsOk)
    {
        double XJunk, YJunk, ZJunk;
        FirstElemNumOffset = GeomToolsClosestElemNum(IsoSurfZone, XCrtPtAtom, YCrtPtAtom, ZCrtPtAtom,
                                                     &XJunk, &YJunk, &ZJunk) - 1;
        if (FirstElemNumOffset < 0) IsOk = FALSE;
    }

    // Extract the topologically independent portion of the isosurface surrounding the Atom. Start
    // with the FirstElemNumOffset and find all connected cells.
    if (IsOk)
    {
        IsoTopoZone = SurfTopoSegComputeSeg(SurfTopoSeg, FirstElemNumOffset);
        if (IsoTopoZone <= 0) IsOk = FALSE;
    }

    // Clean up
    SurfTopoSegDealloc(&SurfTopoSeg);


    // TODO: Delete the full isosurface


    // TODO: consider grid improvement of the IsoTopo surface through small/skinny cell collapse


    // TODO: transfer code from engine.cpp to here.

    ENSURE(VALID_BOOLEAN(IsOk));
    ENSURE(IsoTopoZone >= -1);
    return IsoTopoZone;
}






/*
 * IsoSurfGradPathAdd: For a given max (atom) volume critical point, compute
 *  the "nearly largest" isosurface of the scalar variable surrounding the
 *  atom and no other volume critical point. Then compute the critical points
 *  of the scalar (charge density) gradient magrnitude on the isosurface and the
 *  surface gradient paths connecting the surface critical points. Finally,
 *  compute the surface bundles. These critical points, gradient paths, and
 *  surface "bundles" will be used to compute some of the volume gradient paths,
 *  zero-flux surfaces, and volume bundles.
 *
 * param ZoneVarInfo
 *     Zone and variable information.
 * param AtomCrtPtNum
 *     Offset of the Atom volume critical point that is surrounded by the isosurface.
 * param VolCritPoints
 *     Critical point structure for the volume critical points.
 * param VolBundles
 *     Lists containing the Atom, Bond, Ring, Cage combinations that make up the
 *     irreducible bundles for the volume. These will be reconciled with the
 *     Atom, min, saddle, max combinations created from the topology of the
 *     charge density gradient magrnitude on the isosurface.
 * param CPTolerance
 *     Tolerance - distance from critical point at which the pathline is said
 *     to terminate at the critical point.
 * param IsoSurfGradPath
 *     IsoSurfGradPath data structure to populate for an Atom volume critical point.
 *
 * return
 *     TRUE if successful, FALSE if there were errors.
 */
Boolean_t IsoSurfGradPathAdd(const ZoneVarInfo_pa  VolZoneVarInfo,
                             const ZoneVarInfo_pa  SurfZoneVarInfo,
                             const LgIndex_t       AtomCrtPtNum,
                             const CritPoints_pa   VolCritPoints,
                             const Bundles_pa      VolBundles,
                             const float           CPTolerance,
                             IsoSurfGradPath_pa    IsoSurfGradPath)
{
    Boolean_t      IsOk = TRUE;
    // double         XCrtPtAtom, YCrtPtAtom, ZCrtPtAtom, dummy;
    // char           cdummy;
    LgIndex_t      FirstElemNumOffset = -1;
    EntIndex_t     IsoTopoZone;
    EntIndex_t     VolZoneNum = VolCritPoints->SourceZoneNum;
    EntIndex_t     ChrgDensVarNum = VolZoneVarInfo->ChrgDensVarNum;
    EntIndex_t     UVarNum = VolZoneVarInfo->UVarNum;
    EntIndex_t     VVarNum = VolZoneVarInfo->VVarNum;
    EntIndex_t     WVarNum = VolZoneVarInfo->WVarNum;
    EntIndex_t     GradXGTotVNum = SurfZoneVarInfo->UVarNum;
    EntIndex_t     GradYGTotVNum = SurfZoneVarInfo->VVarNum;
    EntIndex_t     GradZGTotVNum = SurfZoneVarInfo->WVarNum;
    EntIndex_t     TypeVarNum    = SurfZoneVarInfo->TypeVarNum;
    Boolean_t      PeriodicBC    = FALSE;  // TODO: Is this correct?
    CritPoints_pa  SurfCritPoints = NULL;
    Bundles_pa     SurfBundles    = NULL;
    ZoneVarInfo_pa ZoneVarInfo    = NULL;
    SurfElemMap_pa SurfElemMap    = NULL;
	SurfElemMap_pa SurfElemMapPerturbed = NULL;
    Normals_pa     Normals        = NULL;
    EntIndex_t     SurfCPZoneNum  = -1;
    double         GridSpacing    = ZoneVarInfoGetDXFromZoneNum(VolZoneVarInfo, VolZoneVarInfo->ZoneNum);


    REQUIRE(VALID_REF(VolZoneVarInfo));
    REQUIRE(CritPointsIsValid(VolCritPoints));
    REQUIRE(AtomCrtPtNum >= 0 && AtomCrtPtNum < CritPointsGetCount(VolCritPoints));
    // REQUIRE(BundlesIsValid(VolBundles));
    REQUIRE(CPTolerance > 0.0);
    REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));
    // Atom critical points must be an Atom
    REQUIRE(CritPointsGetType(VolCritPoints, AtomCrtPtNum) == -3);

	if (AtomCrtPtNum == 0)
		IsOk = IsOk;
    // For a given max (atom) volume critical point, compute the "nearly largest" isosurface
    // of the scalar variable surrounding the  atom.
    IsoTopoZone = IsoTopoZoneForAtom(VolZoneVarInfo, AtomCrtPtNum, VolCritPoints, CPTolerance);
	ENSURE(IsoTopoZone > 0); // TEMP disabled for debuggin
	//IsOk = (IsoTopoZone > 0);

    IsoSurfGradPath->IsoTopoZone = IsoTopoZone;

    // if (0)
    {
        // For convenience, get pointers to SurfCritPoints and SurfBundles
        SurfCritPoints = IsoSurfGradPath->SurfCritPoints;
        SurfBundles   = IsoSurfGradPath->SurfBundles;

        if (IsOk && IsoTopoZone > 0)
        {
            // TEMPORARY <<<<
            LgIndex_t IMax, JMax, KMax;

            // Compute the surface element map for IsoTopoZone
            TecUtilZoneGetInfo(IsoTopoZone, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);

			//	A copy of the surface info 
            if (IsOk)
            {
                // Compute the element map
                SurfElemMap = SurfElemMapAlloc(IMax);
                IsOk = SurfElemMapCompute(SurfElemMap, IsoTopoZone);
            }

            // Compute the surface normals
            if (IsOk)
            {
                Normals = NormalsAlloc(IMax);
                // IsOk = NormalsComputeUsingSEM(Normals, IsoTopoZone, SurfElemMap);
                IsOk = NormalsComputeUsingIsoVarGrad(Normals, IsoTopoZone, VolZoneVarInfo->UVarNum, 
													VolZoneVarInfo->VVarNum, VolZoneVarInfo->WVarNum);
            }

            if (IsOk)
            {
                ZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "ZoneVarInfo");
                IsOk = ( ZoneVarInfo != NULL );
                if ( IsOk )
                {
                    ZoneVarInfo->ZoneNum = IsoTopoZone;
                    ZoneVarInfo->UVarNum = SurfZoneVarInfo->UVarNum;
                    ZoneVarInfo->VVarNum = SurfZoneVarInfo->VVarNum;
                    ZoneVarInfo->WVarNum = SurfZoneVarInfo->WVarNum;
                    ZoneVarInfo->ChrgDensVarNum = SurfZoneVarInfo->GradMagVarNum;
					ZoneVarInfo->GradMagVarNum = SurfZoneVarInfo->GradMagVarNum;
                    ZoneVarInfo->TypeVarNum = SurfZoneVarInfo->TypeVarNum;
                    ZoneVarInfo->PeriodicBC = PeriodicBC;
                    ZoneVarInfo->SurfElemMap = SurfElemMapPerturbed;
                    ZoneVarInfo->Normals = Normals;
                    ZoneVarInfo->PathDir = StreamDir_Forward;
                }
            }

			//	Get Dimensions of isosurface
			XYZ_s Center;
			double	djunk;
			char	cjunk;

			if (IsOk) IsOk = CritPointsGetPoint(VolCritPoints, AtomCrtPtNum, 
												&Center.X, &Center.Y, &Center.Z, 
												&djunk, &cjunk, &djunk, &djunk, &djunk);

			if (IsOk) IsOk = SurfTopoDimensions(IsoTopoZone, Center, &ZoneVarInfo->ElementCount, 
													&ZoneVarInfo->Radius, 
													&ZoneVarInfo->AvgElemWidth);

			if (IsOk)
			{
				IsoSurfGradPath->ElementCount = ZoneVarInfo->ElementCount;
				IsoSurfGradPath->Radius = VolZoneVarInfo->Radius = ZoneVarInfo->Radius;
				IsoSurfGradPath->AvgElemWidth = ZoneVarInfo->AvgElemWidth;
			}
			
			if (IsOk)
                IsOk = ExtractCriticalPoints(VolZoneVarInfo, ZoneVarInfo, TRUE, SurfCritPoints);
            // IsOk = ExtractCriticalPoints(IsoTopoZone, GradXGTotVNum, GradYGTotVNum, GradZGTotVNum,
            //   GradTotVNum, TypeVarNum, PeriodicBC, SurfCritPoints);

            /* TEMPORARY: Create Critical Point Zone */
            if (IsOk)
            {
                SurfCPZoneNum = CritPointsCreateTPZone(SurfCritPoints, ChrgDensVarNum, TypeVarNum,
                                                       UVarNum, VVarNum, WVarNum, TRUE);
                // if (SurfCPZoneNum > 0) NumZones++;
                CHECK(SurfCPZoneNum > 0);
            }

            // Compute the surface element map for IsoTopoZone
            // TecUtilZoneGetInfo(IsoTopoZone, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
            //                    NULL, NULL, NULL, NULL, NULL, NULL, NULL);

            // Interpolate variables from volume zone to new surface CritPoints zone
            if (IsOk && SurfCPZoneNum > 0)
            {
                Boolean_t IsOk;
                Set_pa Zones = TecUtilSetAlloc(TRUE);
                Set_pa Vars  = TecUtilSetAlloc(TRUE);

                TecUtilSetAddMember(Zones, VolZoneNum, TRUE);

                TecUtilSetAddMember(Vars, UVarNum, TRUE);
                TecUtilSetAddMember(Vars, VVarNum, TRUE);
                TecUtilSetAddMember(Vars, WVarNum, TRUE);

                IsOk = TecUtilLinearInterpolate(Zones, SurfCPZoneNum, Vars, 0.0,
                                                LinearInterpMode_SetToConst);
                TecUtilSetDealloc(&Zones);
                TecUtilSetDealloc(&Vars);
            }

            /*
             * To get surface Saddle-Max lines, seed streamtrace a small step in the
             * PrincDir away from Saddle
             */
            if (IsOk)
                // if (0)
            {
                GradPath_pa    GradPath = NULL;
				//	TIM:	beginning and end offset for saddle(ring) surfCPs
                LgIndex_t      BegOffset = CritPointsGetBegOffset(SurfCritPoints, -1);
                LgIndex_t      EndOffset = CritPointsGetEndOffset(SurfCritPoints, -1);
                LgIndex_t      NumCrtPts = CritPointsGetCount(SurfCritPoints);
                LgIndex_t      ii;

                /*
                ZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "ZoneVarInfo");
                ZoneVarInfo->ZoneNum = IsoTopoZone;
                ZoneVarInfo->UVarNum = SurfZoneVarInfo->UVarNum;
                ZoneVarInfo->VVarNum = SurfZoneVarInfo->VVarNum;
                ZoneVarInfo->WVarNum = SurfZoneVarInfo->WVarNum;
                ZoneVarInfo->ChrgDensVarNum = SurfZoneVarInfo->ChrgDensVarNum;
                ZoneVarInfo->PeriodicBC = PeriodicBC;
                ZoneVarInfo->SurfElemMap = SurfElemMap;
                ZoneVarInfo->Normals = Normals;
                ZoneVarInfo->PathDir = StreamDir_Forward;
                */

				//	Loop over all surf crit points to find minimum distance between any two
				//	of different type.
				//	Let DXYZ be one sixteenth the minimum distance between surf cps.
				double	DXYZFactor = 0.05, MinSaddleMaxDist = -1.0, MinSaddleMinDist = -1.0;
				const double	CPTolFactor = 0.5;

				if (IsOk)
				{
					if (SurfCritPoints->NumCrtPts > 1)
						IsOk = CritPointsMinDistance(SurfCritPoints, FALSE, 0, 0, &SurfCritPoints->MinCPDistance);
					else
						SurfCritPoints->MinCPDistance = ZoneVarInfo->AvgElemWidth;
				}
				
				if (SurfCritPoints->NumCrtPtsM3 > 0 && SurfCritPoints->NumCrtPtsM1 > 0 && SurfCritPoints->NumCrtPtsP3 > 0)
				{
					if (IsOk) IsOk = CritPointsMinDistance(SurfCritPoints, TRUE, (char)0, (char)-2, &MinSaddleMaxDist);
					if (IsOk) IsOk = CritPointsMinDistance(SurfCritPoints, TRUE, (char)0, (char)2, &MinSaddleMinDist);
					CHECK(MinSaddleMaxDist > 0 && MinSaddleMinDist > 0);
					IsOk = (MinSaddleMaxDist > 0 && MinSaddleMinDist > 0);
				}
				else
				{
					MinSaddleMaxDist = MinSaddleMinDist = SurfCritPoints->MinCPDistance;
				}
				
				
				

                /* Set-up status bar */
                // if (ShowStatusBar)
                //   TecUtilStatusStartPercentDone( "Finding Bond Lines", TRUE, TRUE);

				//	TIM:	The last saddle index is endoffset - 1, where ii stops
				//			This is looping over saddle surfCPs and building grad
				//			to other surfCPs, so begCP should always be between
				//			BegOffset and EndOffset, and endCP should not.

				XYZ_s		OldGPTestPt;
				const Boolean_t	MakeSurfGPZones = TRUE;

				LgIndex_t	MaxSteps = 1500;
				double		StepSize = 0.25;
				double		MinStepSize = 1e-5;

				IsOk = GradPathSTSetProperties(MaxSteps, StepSize, MinStepSize);

                for (ii = BegOffset; IsOk && ii < EndOffset; ii++)
                    // if (0)
                {
                    double XCrtPt, YCrtPt, ZCrtPt, PrincDirX, PrincDirY, PrincDirZ, dummy;
                    double Princ2DirX, Princ2DirY, Princ2DirZ;
                    char   TypeCrtPt;
                    double XPos, YPos, ZPos;

                    // Used in defining surface bundles
                    LgIndex_t MaxCPFirst = -1;
                    LgIndex_t MaxCPSecond = -1;

					/* Inform Tecplot that major data operation is beginning */
					TecUtilDataLoadBegin();
					TecUtilLockStart(AddOnID);
					
					IsOk = CritPointsGetPoint(SurfCritPoints, ii, &XCrtPt, &YCrtPt, &ZCrtPt,
                                              &dummy, &TypeCrtPt, &PrincDirX, &PrincDirY, &PrincDirZ);
					XYZ_s	CrtPtXYZ;
					CrtPtXYZ.X = XCrtPt;
					CrtPtXYZ.Y = YCrtPt;
					CrtPtXYZ.Z = ZCrtPt;

                    // Surface GradPath from saddle to max - first one (in PrincDir)

                    // XPos = XCrtPt + 0.5 * PrincDirX;
                    // YPos = YCrtPt + 0.5 * PrincDirY;
                    // ZPos = ZCrtPt + 0.5 * PrincDirZ;

                    //XPos = XCrtPt + 0.5 * GridSpacing * PrincDirX;
                    //YPos = YCrtPt + 0.5 * GridSpacing * PrincDirY;
                    //ZPos = ZCrtPt + 0.5 * GridSpacing * PrincDirZ;

					GradPath = GradPathAlloc();

					/* Integrate to compute path lines */
					if (IsOk)
					{
						double    junk = 0.0;
						LgIndex_t NumPathPoints = 0;

						Boolean_t	EndCPSuccess = FALSE;
						for (int GPRetry = 1 ; !EndCPSuccess && GPRetry <= 10 ; GPRetry++)
						{
							GradPathClear(GradPath);

							XPos = XCrtPt + (double)GPRetry * MinSaddleMaxDist * DXYZFactor * PrincDirX;
							YPos = YCrtPt + (double)GPRetry * MinSaddleMaxDist * DXYZFactor * PrincDirY;
							ZPos = ZCrtPt + (double)GPRetry * MinSaddleMaxDist * DXYZFactor * PrincDirZ;

							XYZ_s	TempPt1, TempPt2;

							TempPt1.X = XPos;
							TempPt1.Y = YPos;
							TempPt1.Z = ZPos;

							IsOk = GradPathProjectPointToSurf(ZoneVarInfo, 
								DistanceSquaredXYZ(CrtPtXYZ, TempPt1),
								TempPt1, &TempPt2);
							if (IsOk)
							{
								XPos = TempPt2.X;
								YPos = TempPt2.Y;
								ZPos = TempPt2.Z;
							}

							/* Seed first point */
							IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk);
							GradPath->BeginCrtPtNum = ii;

							/* Integrate gradient path line */
							ZoneVarInfo->PathDir = StreamDir_Forward;

							// TEMP
							//   {
							//     char Message[200];
							//     sprintf_s(Message, "Calling GradPathAddSurf for Saddle-Max %d",ii);
							//     TecUtilDialogMessageBox(Message, MessageBox_Warning);
							//     sprintf_s(Message, "XPos,YPos,ZPos=%g,%g,%g",XPos,YPos,ZPos);
							//     TecUtilDialogMessageBox(Message, MessageBox_Warning);
							//     sprintf_s(Message, "XCrtPt,YCrtPt,ZCrtPt=%g,%g,%g",XCrtPt,YCrtPt,ZCrtPt);
							//     TecUtilDialogMessageBox(Message, MessageBox_Warning);
							//     sprintf_s(Message, "PrincDirX,PrincDirY,PrincDirZ=%g,%g,%g",PrincDirX,PrincDirY,PrincDirZ);
							//     TecUtilDialogMessageBox(Message, MessageBox_Warning);
							//   }
							// IsOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Forward,
							//   SurfCritPoints, 1.0, &NumPathPoints, GradPath);

// 							IsOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Forward,
// 													SurfCritPoints, IncidentCPTol, &NumPathPoints, GradPath);
							IsOk = GradPathSTAddSurf(VolZoneVarInfo, StreamDir_Forward, SurfCritPoints, &NumPathPoints, MinSaddleMaxDist * CPTolFactor, GradPath);

							if (GradPath->EndCrtPtNum >= 0) 
							{
								if (IsOk) IsOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

								if (!GradPathCheckLenRatio(GradPath, SurfCritPoints, VolZoneVarInfo->Radius))
									GradPath->EndCrtPtNum = -1;

								CHECK(GradPathIsValid(GradPath));

							
								IsOk = GradPathGetPoint(GradPath, (LgIndex_t)(NUMRESAMPLEDPOINTSSURFACE / 4), &OldGPTestPt.X, &OldGPTestPt.Y, &OldGPTestPt.Z, &djunk);
								if (IsOk) EndCPSuccess = TRUE;
							}
							else
							{
								StepSize *= 0.9;
								IsOk = GradPathSTSetProperties(MaxSteps, StepSize, MinStepSize);
							}
						}

						//IsOk = EndCPSuccess;

						/* Add a Tecplot zone for the Surface Saddle(Ring)-Max(Cage) connector */
						if (IsOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts
							&& GradPath->EndCrtPtNum < BegOffset)
						{
							LgIndex_t  M1CrtPtNum = ii - BegOffset;
							EntIndex_t ConnectorZoneNum = 0;
							EntIndex_t ConnectorType = -4;  // Atom-Bond

							MaxCPFirst = GradPath->EndCrtPtNum;

							// Resample before writing the tecplot zone
							//if (IsOk) IsOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

							if (IsOk) IsOk = IsoSurfGradPathAppendGP(IsoSurfGradPath, GradPath);

							/*char	ZoneName[200];
                            
							sprintf_s(ZoneName, "Max_%d_Saddle_%d_Zone_%d", GradPath->EndCrtPtNum, M1CrtPtNum, ZoneNum);*/

							if (MakeSurfGPZones)
							{
								ConnectorZoneNum = CreateConnectorZone(IsoTopoZone, ChrgDensVarNum, TypeVarNum, SurfCritPoints, GradPath);

								if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
							}
						}
                    }

                    // Surface GradPath from saddle to max - second one (opposite PrincDir)

                    // IsOk = TecUtilStreamtraceAdd (1, Streamtrace_VolumeLine, StreamDir_Forward,
                    //                               XPos, YPos, ZPos, 0.0, 0.0, 0.0);
                    // XPos = XCrtPt - 0.5 * PrincDirX;
                    // YPos = YCrtPt - 0.5 * PrincDirY;
                    // ZPos = ZCrtPt - 0.5 * PrincDirZ;

					//XPos = XCrtPt - 0.5 * GridSpacing * PrincDirX;
					//YPos = YCrtPt - 0.5 * GridSpacing * PrincDirY;
					//ZPos = ZCrtPt - 0.5 * GridSpacing * PrincDirZ;

					GradPath = GradPathAlloc();

					StepSize = 0.25;

					/* Integrate to compute path lines */
					if (IsOk)
					{
						double    junk = 0.0;
						LgIndex_t NumPathPoints = 0;

						Boolean_t	EndCPSuccess = FALSE;
						for (int GPRetry = 1 ; !EndCPSuccess && GPRetry <= 10 ; GPRetry++)
						{	
							for (int ChangeSign = 1 ; ChangeSign > -2 ; ChangeSign -= 2)
							{
								XPos = XCrtPt - (double)GPRetry * (double)ChangeSign * MinSaddleMaxDist * DXYZFactor * PrincDirX;
								YPos = YCrtPt - (double)GPRetry * (double)ChangeSign * MinSaddleMaxDist * DXYZFactor * PrincDirY;
								ZPos = ZCrtPt - (double)GPRetry * (double)ChangeSign * MinSaddleMaxDist * DXYZFactor * PrincDirZ;

								XYZ_s	TempPt1, TempPt2;
								char	cjunk;

								TempPt1.X = XPos;
								TempPt1.Y = YPos;
								TempPt1.Z = ZPos;

								IsOk = GradPathProjectPointToSurf(ZoneVarInfo, 
									DistanceSquaredXYZ(CrtPtXYZ, TempPt1),
									TempPt1, &TempPt2);
								if (IsOk)
								{
									XPos = TempPt2.X;
									YPos = TempPt2.Y;
									ZPos = TempPt2.Z;
								}
							
								GradPathClear(GradPath);
								
								/* Seed first point */
								IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk);
								GradPath->BeginCrtPtNum = ii;

								/* Integrate gradient path line */
								ZoneVarInfo->PathDir = StreamDir_Forward;


								// TEMP
								//   {
								//     char Message[200];
								//     sprintf_s(Message, "Calling 2nd GradPathAddSurf for Saddle-Max %d",ii);
								//     TecUtilDialogMessageBox(Message, MessageBox_Warning);
								//   }
								//	TIM:	changed CPTolorance (5th arg) from 1.0 to 0.5 to match first
								//			grad path above
	// 							IsOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Forward,
	// 											SurfCritPoints, IncidentCPTol, &NumPathPoints, GradPath);

								IsOk = GradPathSTAddSurf(VolZoneVarInfo, StreamDir_Forward, SurfCritPoints, &NumPathPoints, MinSaddleMaxDist * CPTolFactor, GradPath);

								if (GradPath->EndCrtPtNum != -1)
								{
									if (IsOk) IsOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

									if (!GradPathCheckLenRatio(GradPath, SurfCritPoints, VolZoneVarInfo->Radius))
										GradPath->EndCrtPtNum = -1;
									
									CHECK(GradPathIsValid(GradPath));

									if (IsOk) IsOk = GradPathGetPoint(GradPath, (LgIndex_t)(NUMRESAMPLEDPOINTSSURFACE / 4), &TempPt1.X, &TempPt1.Y, &TempPt1.Z, &djunk);
									if (IsOk && DistanceSquaredXYZ(TempPt1, CrtPtXYZ) < DistanceSquaredXYZ(TempPt1, OldGPTestPt))
									{
										EndCPSuccess = TRUE;
										break;
									}
								}
								else
								{
									StepSize *= 0.9;
									IsOk = GradPathSTSetProperties(MaxSteps, StepSize, MinStepSize);
								}
							}
						}

						//IsOk = EndCPSuccess;

						/* Add a Tecplot zone for the Surface Saddle(Ring)-Max(Cage) connector */
						if (IsOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts
							&& GradPath->EndCrtPtNum < BegOffset)
						{
							// char       ZoneName[200];
							LgIndex_t  M1CrtPtNum = ii - BegOffset;
							EntIndex_t ConnectorZoneNum = 0;
							EntIndex_t ConnectorType = -4;  // Atom-Bond

							MaxCPSecond = GradPath->EndCrtPtNum;

							// Resample before writing the tecplot zone
							//if (IsOk) IsOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

							if (IsOk) IsOk = IsoSurfGradPathAppendGP(IsoSurfGradPath, GradPath);

							/*char	ZoneName[200];
                            
							sprintf_s(ZoneName, "Max_%d_Saddle_%d_Zone_%d", GradPath->EndCrtPtNum, M1CrtPtNum, ZoneNum);*/

							if (MakeSurfGPZones)
							{
								ConnectorZoneNum = CreateConnectorZone(IsoTopoZone, ChrgDensVarNum, TypeVarNum, SurfCritPoints, GradPath);

								if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
							}
						}
					}
                    
					// Find surface-tangent vector perpendicular to PrincDir
                    // Assume local ChargeDensityGradient is in direction of surface normal

                    if (IsOk)
                    {
                        FieldData_pa UVarFDPtr = TecUtilDataValueGetReadableRef(SurfCPZoneNum, UVarNum);
                        FieldData_pa VVarFDPtr = TecUtilDataValueGetReadableRef(SurfCPZoneNum, VVarNum);
                        FieldData_pa WVarFDPtr = TecUtilDataValueGetReadableRef(SurfCPZoneNum, WVarNum);

                        double Nx = TecUtilDataValueGetByRef(UVarFDPtr, ii + 1);
                        double Ny = TecUtilDataValueGetByRef(VVarFDPtr, ii + 1);
                        double Nz = TecUtilDataValueGetByRef(WVarFDPtr, ii + 1);

                        double RNTot = 1.0 / sqrt(Nx * Nx + Ny * Ny + Nz * Nz);

                        Nx *= RNTot;
                        Ny *= RNTot;
                        Nz *= RNTot;

                        Princ2DirX = Ny * PrincDirZ - Nz * PrincDirY;
                        Princ2DirY = Nz * PrincDirX - Nx * PrincDirZ;
                        Princ2DirZ = Nx * PrincDirY - Ny * PrincDirX;
                    }

                    // Surface GradPath from saddle to min - first one (aligned with Princ2Dir)

                    // IsOk = TecUtilStreamtraceAdd (1, Streamtrace_VolumeLine, StreamDir_Forward,
                    //                               XPos, YPos, ZPos, 0.0, 0.0, 0.0);
                    // XPos = XCrtPt + 0.1 * Princ2DirX;
                    // YPos = YCrtPt + 0.1 * Princ2DirY;
                    // ZPos = ZCrtPt + 0.1 * Princ2DirZ;
                    
					//XPos = XCrtPt + 0.5 * Princ2DirX;
     //               YPos = YCrtPt + 0.5 * Princ2DirY;
     //               ZPos = ZCrtPt + 0.5 * Princ2DirZ;

					//XPos = XCrtPt + 0.5 * GridSpacing * Princ2DirX;
					//YPos = YCrtPt + 0.5 * GridSpacing * Princ2DirY;
					//ZPos = ZCrtPt + 0.5 * GridSpacing * Princ2DirZ;

					//XPos = XCrtPt + DXYZ * Princ2DirX;
					//YPos = YCrtPt + DXYZ * Princ2DirY;
					//ZPos = ZCrtPt + DXYZ * Princ2DirZ;

                    GradPath = GradPathAlloc();

					StepSize = 0.25;

					/* Integrate to compute path lines */
					if (IsOk)
					{
						double    junk = 0.0;
						LgIndex_t NumPathPoints = 0;

						Boolean_t	EndCPSuccess = FALSE;
						for (int GPRetry = 1 ; !EndCPSuccess && GPRetry <= 10 ; GPRetry++)
						{
							GradPathClear(GradPath);

							XPos = XCrtPt + (double)GPRetry * MinSaddleMinDist * DXYZFactor * Princ2DirX;
							YPos = YCrtPt + (double)GPRetry * MinSaddleMinDist * DXYZFactor * Princ2DirY;
							ZPos = ZCrtPt + (double)GPRetry * MinSaddleMinDist * DXYZFactor * Princ2DirZ;

							XYZ_s	TempPt1, TempPt2;

							TempPt1.X = XPos;
							TempPt1.Y = YPos;
							TempPt1.Z = ZPos;

							IsOk = GradPathProjectPointToSurf(ZoneVarInfo, 
								DistanceSquaredXYZ(CrtPtXYZ, TempPt1),
								TempPt1, &TempPt2);
							if (IsOk)
							{
								XPos = TempPt2.X;
								YPos = TempPt2.Y;
								ZPos = TempPt2.Z;
							}

							/* Seed first point */
							IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk);
							GradPath->BeginCrtPtNum = ii;

							/* Integrate gradient path line */
							ZoneVarInfo->PathDir = StreamDir_Reverse;

							// TEMP
							//   {
							//     char Message[200];
							//     sprintf_s(Message, "Calling GradPathAddSurf for Saddle-Min %d",ii);
							//     TecUtilDialogMessageBox(Message, MessageBox_Warning);
							//   }
							// IsOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Reverse,
							//  SurfCritPoints, 1.0, &NumPathPoints, GradPath);
// 							IsOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Reverse,
// 												   SurfCritPoints, IncidentCPTol, &NumPathPoints, GradPath);

							IsOk = GradPathSTAddSurf(VolZoneVarInfo, StreamDir_Reverse, SurfCritPoints, &NumPathPoints, MinSaddleMinDist * CPTolFactor, GradPath);
							
							if (IsOk) 
							{
								if (GradPath->EndCrtPtNum >= 0)
								{
									IsOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

									if (!GradPathCheckLenRatio(GradPath, SurfCritPoints, VolZoneVarInfo->Radius))
										GradPath->EndCrtPtNum = -1;

									CHECK(GradPathIsValid(GradPath));
								
									IsOk = GradPathGetPoint(GradPath, (LgIndex_t)(NUMRESAMPLEDPOINTSSURFACE / 4), &OldGPTestPt.X, &OldGPTestPt.Y, &OldGPTestPt.Z, &djunk);
									EndCPSuccess = TRUE;
								}
								else
								{
									StepSize *= 0.9;
									IsOk = GradPathSTSetProperties(MaxSteps, StepSize, MinStepSize);
								}
							}
						}

						//IsOk = EndCPSuccess;

						/* Add a Tecplot zone and two SurfBundles for the saddle(ring)-min(bond) connector */
						if (IsOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts
							&& GradPath->EndCrtPtNum >= EndOffset)
						{
							// char       ZoneName[200];
							LgIndex_t  M1CrtPtNum = ii - BegOffset;
							// EntIndex_t ConnectorZoneNum = 0;
							EntIndex_t ConnectorType = -4;  // Atom-Bond

							LgIndex_t  MinNum = GradPath->EndCrtPtNum - CritPointsGetBegOffset(SurfCritPoints, (char)(3)) + 1;
							LgIndex_t  SaddleNum = ii - BegOffset + 1;

							// Define two surface bundles containing the saddle-min path
							if (MaxCPFirst >= 0)
							{
								LgIndex_t MaxNum = MaxCPFirst + 1;
								IsOk = BundlesAppendAtEnd(SurfBundles, AtomCrtPtNum, MinNum, SaddleNum, MaxNum);
							}
							if (MaxCPSecond >= 0)
							{
								LgIndex_t MaxNum = MaxCPSecond + 1;
								IsOk = BundlesAppendAtEnd(SurfBundles, AtomCrtPtNum, MinNum, SaddleNum, MaxNum);
							}

							// Resample before writing the tecplot zone
							//if (IsOk) IsOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

							if (IsOk) IsOk = IsoSurfGradPathAppendGP(IsoSurfGradPath, GradPath);

							/*char	ZoneName[200];
                            
							sprintf_s(ZoneName, "Max_%d_Saddle_%d_Zone_%d", GradPath->EndCrtPtNum, M1CrtPtNum, ZoneNum);*/

							if (MakeSurfGPZones)
							{
								EntIndex_t ConnectorZoneNum = CreateConnectorZone(IsoTopoZone, ChrgDensVarNum, TypeVarNum, SurfCritPoints, GradPath);

								if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
							}
						}
                    }

                    // Surface GradPath from saddle to min - Second one (opposite Princ2Dir)

                    // IsOk = TecUtilStreamtraceAdd (1, Streamtrace_VolumeLine, StreamDir_Forward,
                    //                               XPos, YPos, ZPos, 0.0, 0.0, 0.0);
                    // XPos = XCrtPt + 0.1 * Princ2DirX;
                    // YPos = YCrtPt + 0.1 * Princ2DirY;
                    // ZPos = ZCrtPt + 0.1 * Princ2DirZ;

                    //XPos = XCrtPt - 0.5 * Princ2DirX;
                    //YPos = YCrtPt - 0.5 * Princ2DirY;
                    //ZPos = ZCrtPt - 0.5 * Princ2DirZ;

					//XPos = XCrtPt - 0.5 * GridSpacing * Princ2DirX;
					//YPos = YCrtPt - 0.5 * GridSpacing * Princ2DirY;
					//ZPos = ZCrtPt - 0.5 * GridSpacing * Princ2DirZ;

					//XPos = XCrtPt - DXYZ * Princ2DirX;
					//YPos = YCrtPt - DXYZ * Princ2DirY;
					//ZPos = ZCrtPt - DXYZ * Princ2DirZ;
					
					GradPath = GradPathAlloc();

					StepSize = 0.25;

					/* Integrate to compute path lines */
					if (IsOk)
					{
						double    junk = 0.0;
						LgIndex_t NumPathPoints = 0;

						Boolean_t	EndCPSuccess = FALSE;
						for (int GPRetry = 1 ; !EndCPSuccess && GPRetry <= 10 ; GPRetry++)
						{
							//	This is checking that the current seed point isn't in between
							//	the beginning surf CP and the previously found end surf cp, in which
							//	case it would terminate at the same cp as previously found.
							for (int ChangeSign = 1 ; ChangeSign > -2 ; ChangeSign -= 2)
							{						
// 								for (int SeedPtShift = 1 ; SeedPtShift <= GPRetry ; SeedPtShift++)
// 								{
// 									if (SeedPtShift == 1)
// 									{
										XPos = XCrtPt - (double)GPRetry * (double)ChangeSign * DXYZFactor * Princ2DirX;
										YPos = YCrtPt - (double)GPRetry * (double)ChangeSign * DXYZFactor * Princ2DirY;
										ZPos = ZCrtPt - (double)GPRetry * (double)ChangeSign * DXYZFactor * Princ2DirZ;
// 									}
// 									else
// 									{
// 										XPos -= (double)ChangeSign * DXYZFactor * Princ2DirX;
// 										YPos -= (double)ChangeSign * DXYZFactor * Princ2DirY;
// 										ZPos -= (double)ChangeSign * DXYZFactor * Princ2DirZ;
// 									}

									XYZ_s	TempPt1, TempPt2;
									char	cjunk;

									TempPt1.X = XPos;
									TempPt1.Y = YPos;
									TempPt1.Z = ZPos;

									IsOk = GradPathProjectPointToSurf(ZoneVarInfo, 
														DistanceSquaredXYZ(CrtPtXYZ, TempPt1),
														TempPt1, &TempPt2);
									if (IsOk)
									{
										XPos = TempPt2.X;
										YPos = TempPt2.Y;
										ZPos = TempPt2.Z;
									}
//								}

								GradPathClear(GradPath);
									
								/* Seed first point */
								IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk);
								GradPath->BeginCrtPtNum = ii;

								/* Integrate gradient path line */
								ZoneVarInfo->PathDir = StreamDir_Reverse;


								// TEMP
								//   {
								//     char Message[200];
								//     sprintf_s(Message, "Calling GradPathAddSurf for Saddle-Min %d",ii);
								//     TecUtilDialogMessageBox(Message, MessageBox_Warning);
								//   }
	// 							IsOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Reverse,
	// 												SurfCritPoints, IncidentCPTol, &NumPathPoints, GradPath);

								IsOk = GradPathSTAddSurf(VolZoneVarInfo, StreamDir_Reverse, SurfCritPoints, &NumPathPoints, MinSaddleMinDist * CPTolFactor, GradPath);

								if (GradPath->EndCrtPtNum != -1)
								{
									if (IsOk) IsOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

									if (!GradPathCheckLenRatio(GradPath, SurfCritPoints, VolZoneVarInfo->Radius))
										GradPath->EndCrtPtNum = -1;

									CHECK(GradPathIsValid(GradPath));
								
									if (IsOk) IsOk = GradPathGetPoint(GradPath, (LgIndex_t)(NUMRESAMPLEDPOINTSSURFACE / 4), &TempPt1.X, &TempPt1.Y, &TempPt1.Z, &djunk);
									if (IsOk && DistanceSquaredXYZ(TempPt1, CrtPtXYZ) < DistanceSquaredXYZ(TempPt1, OldGPTestPt))
									{
										EndCPSuccess = TRUE;
										break;
									}
								}
								else
								{
									StepSize *= 0.9;
									IsOk = GradPathSTSetProperties(MaxSteps, StepSize, MinStepSize);
								}
							}
						}

						//IsOk = EndCPSuccess;

						/* Add a Tecplot zone and two SurfBundles for the saddle(ring)-min(bond) connector */
						if (IsOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts
							&& GradPath->EndCrtPtNum >= EndOffset)
						{
							// char       ZoneName[200];
							LgIndex_t  M1CrtPtNum = ii - BegOffset;
							// EntIndex_t ConnectorZoneNum = 0;
							EntIndex_t ConnectorType = -4;  // Atom-Bond

							LgIndex_t  MinNum = GradPath->EndCrtPtNum - CritPointsGetBegOffset(SurfCritPoints, (char)(3)) + 1;
							LgIndex_t  SaddleNum = ii - BegOffset + 1;

							// Define two surface bundles containing the saddle-min path
							if (MaxCPFirst >= 0)
							{
								LgIndex_t MaxNum = MaxCPFirst + 1;
								IsOk = BundlesAppendAtEnd(SurfBundles, AtomCrtPtNum, MinNum, SaddleNum, MaxNum);
							}
							if (MaxCPSecond >= 0)
							{
								LgIndex_t MaxNum = MaxCPSecond + 1;
								IsOk = BundlesAppendAtEnd(SurfBundles, AtomCrtPtNum, MinNum, SaddleNum, MaxNum);
							}

							// Resample before writing the tecplot zone
							//if (IsOk) IsOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

							if (IsOk) IsOk = IsoSurfGradPathAppendGP(IsoSurfGradPath, GradPath);

							/*char	ZoneName[200];
                            
							sprintf_s(ZoneName, "Max_%d_Saddle_%d_Zone_%d", GradPath->EndCrtPtNum, M1CrtPtNum, ZoneNum);*/

							if (MakeSurfGPZones)
							{
								EntIndex_t ConnectorZoneNum = CreateConnectorZone(IsoTopoZone, ChrgDensVarNum, TypeVarNum, SurfCritPoints, GradPath);

								if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
							}
						}
                    }

                    /* Update status bar */
                    /*
                    if (ShowStatusBar)
                      {
                        char PercentDoneText[200];
                        int  PercentDone;
                        PercentDone = (int)( (100 * (ii-BegOffset))/(EndOffset - BegOffset) );
                        sprintf_s(PercentDoneText,"Finding surface connect lines: %d Percent Done", PercentDone);
                        TecUtilStatusSetPercentDoneText(PercentDoneText);
                        if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop=TRUE;
                        if (IsStop) IsOk = FALSE;
                      } */

                    /* Inform Tecplot that major data operation is ending */
                    TecUtilLockFinish(AddOnID);
                    TecUtilDataLoadEnd();
                }

                /* Clean up allocated structures/arrays */
                if (!IsOk && GradPath != NULL) GradPathDealloc(&GradPath);
            }
            // if(ShowStatusBar) TecUtilStatusFinishPercentDone();
        }

        // Compute the surface min-max lines (seed atom-bond-cage surfaces in the volume topology)
        // TEMP: turn off for debugging
        //if (IsOk && BundlesGetCount(SurfBundles) > 0)
        if (0)
        {
            GradPath_pa  GradPath = NULL;
            LgIndex_t    NumBundles = BundlesGetCount(SurfBundles);
            LgIndex_t    ib;
            LgIndex_t    AtomVolCP, MinCP, SaddleCP, MaxCP;
            XYZ_s        SeedPos;

            for (ib = 0; IsOk && ib < NumBundles; ib++)
                // for (ib=0; IsOk && ib<1; ib++)
            {
                LgIndex_t FirstElemNumOffset;
                char      junkchar;
                double    junk1, junk2, junk3, junk4;
                double    XMin, YMin, ZMin, XMax, YMax, ZMax, XAve, YAve, ZAve, XSeed, YSeed, ZSeed;
                LgIndex_t NumPathPoints = 0;
                LgIndex_t NumCrtPts = CritPointsGetCount(SurfCritPoints);

                // Get the min and max critical point numbers
                IsOk = BundlesGetFromOffset(SurfBundles, ib, &AtomVolCP, &MinCP, &SaddleCP, &MaxCP);

                // Find the average mid-point on the line between MinCP and MaxCP
                if (MinCP >= 0 && MaxCP >= 0)
                {
                    IsOk = CritPointsGetPoint(SurfCritPoints, MinCP, &XMin, &YMin, &ZMin, &junk1, &junkchar, &junk2, &junk3, &junk4);
                    IsOk = CritPointsGetPoint(SurfCritPoints, MaxCP, &XMax, &YMax, &ZMax, &junk1, &junkchar, &junk2, &junk3, &junk4);
                    XAve = 0.5 * (XMin + XMax);
                    YAve = 0.5 * (YMin + YMax);
                    ZAve = 0.5 * (ZMin + ZMax);

                    // Find the nearest point on the IsoSurface to the mid-point
                    FirstElemNumOffset = GeomToolsClosestElemNum(IsoTopoZone, XAve, YAve, ZAve,
                                                                 &XSeed, &YSeed, &ZSeed) - 1;
                    SeedPos.X = XSeed;
                    SeedPos.Y = YSeed;
                    SeedPos.Z = ZSeed;

                    // Add a surface grad-path
                    GradPath = GradPathSurfAddMidway(ZoneVarInfo, MAXGRADPATHPOINTS, SurfCritPoints,
                                                     SeedPos, 1.0, &NumPathPoints);

                    // Create Tecplot zone for Min-Max connector
                    if (IsOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts)
                    {
                        // char       ZoneName[200];
                        // EntIndex_t ConnectorZoneNum = 0;
                        EntIndex_t ConnectorType = -4;  // Atom-Bond

                        // sprintf_s(ZoneName, "Max_%d_Saddle_%d_Zone_%d", GradPath->EndCrtPtNum, GradPath->BeginCrtPtNum, ZoneNum);

                        // Resample before writing the tecplot zone
                        IsOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

                        if (IsOk) IsOk = IsoSurfGradPathAppendGP(IsoSurfGradPath, GradPath);


                        /*
                        ConnectorZoneNum = CreateConnectorZone(IsoTopoZone, ChrgDensVarNum, TypeVarNum, SurfCritPoints, GradPath);

                        if (ConnectorZoneNum == TECUTILSETNOTMEMBER) IsOk = FALSE;
                        */
                    }

                    /* Clean up allocated structures/arrays */
                    if (!IsOk && GradPath != NULL) GradPathDealloc(&GradPath);
                }
            }
        }
        /* Clean up allocated structures/arrays */
        if (ZoneVarInfo != NULL)
            FREE_ITEM(ZoneVarInfo, "ZoneVarInfo");
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/** TIM
 * Take a given surface grad path, if its beginning and/or end
 * point(s) (surf CPs) are sufficiently close to the point of 
 * isosurface-volume grad paths intersection then replace it with 
 * the intersection point. If the point is far from the intersection
 * then append that point to the surface grad path
 *
 * NOTE:	Currently assumes that a volume grad path exists for each
 *			surf CP at ends of supplied surf grad path.
 *
 * param SurfGradPath
 *		Surface grad path to have its ends adjusted.
 * param SurfCritPoints
 *		Isosurface critical points
 * param VolCritPoints
 *		Volume critical points
 * param AtomNum
 *		Index of atom CP isosurface is based on (do I need this, or
 *		can I take it from the SurfCritPoints->SourceZone property?)
 * param GPFoundTolorance
 *		Tolerance to determine if a volume grad path corresponds
 *		to a surface CP
 * param ReplaceTolorance
 *		Tolerance to determine whether the surface grad path's end
 *		point should be replaced with the intersection point or if
 *		the intersection point should be appended to the surf grad
 *		path
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t	IsoSurfGradPathSurfConnectSurfSaddleToVolRing(IsoSurfGradPath_pa		IsoSurfGradPath,
														 const CritPoints_pa		SurfCritPoints,
														 const CritPoints_pa		VolCritPoints,
														 const LgIndex_t			AtomNum)
{
	Boolean_t	IsOk = TRUE;

	REQUIRE(IsoSurfGradPathIsValid(IsoSurfGradPath));
	REQUIRE(CritPointsIsValid(SurfCritPoints));
	REQUIRE(CritPointsIsValid(VolCritPoints));
	REQUIRE(SurfCritPoints->Dimensions == 2);
	REQUIRE(VolCritPoints->Dimensions == 3);
	REQUIRE(AtomNum >= 0);

	const char					RingCPType = 1;
	const char					SaddleCPType = 0, MaxType = -2, MinType = 2;
	std::vector <EntIndex_t>	AtomRingPathZoneNums;
	std::vector <XYZ_s>			IntersectionPTs;
	std::vector <double>		IntersectionRhoValues;
	LgIndex_t					RingBegOffset = CritPointsGetBegOffset(VolCritPoints, RingCPType);
	LgIndex_t					RingEndOffset = CritPointsGetEndOffset(VolCritPoints, RingCPType);
	LgIndex_t					SaddleBegOffset = CritPointsGetBegOffset(SurfCritPoints, SaddleCPType);
	LgIndex_t					SaddleEndOffset = CritPointsGetEndOffset(SurfCritPoints, SaddleCPType);

	//	First, loop over Ring CPs to find those that are connected to AtomNum.
	//	When found, find the intersection of the ring-atom GP with the isosurface.
	for (LgIndex_t RingCPNum = RingBegOffset ; IsOk && RingCPNum < RingEndOffset ; RingCPNum++)
	{
		GradPath_pa		RingAtomGP = NULL;
		XYZ_s			RingCP, SurfPT;
		double			Rho, djunk;
		char			cjunk;

		//	Check that path exists, and if it does get the GP object.
		EntIndex_t	RingAtomGPZoneNum = GradPathTPZoneFromBegEndCP(RingCPNum, AtomNum, 
														VolCritPoints->SourceZoneNum);
		if (RingAtomGPZoneNum > 0)
		{
			LgIndex_t I, J, K;
			TecUtilZoneGetIJK(RingAtomGPZoneNum, &I, &J, &K);
			if (I > 1 && J == 1 && K == 1){
				AtomRingPathZoneNums.push_back(RingAtomGPZoneNum);
				RingAtomGP = GradPathGetFromTPZone(RingAtomGPZoneNum);
				IsOk = GradPathIsValid(RingAtomGP);
			}
			else IsOk = FALSE;

			if (IsOk)
			{
				IsOk = CritPointsGetPoint(VolCritPoints, RingCPNum, 
									&RingCP.X, &RingCP.Y, &RingCP.Z, 
									&djunk, &cjunk, &djunk, &djunk, &djunk);

				//	Identify the closest element to the ring CP
				LgIndex_t	ClosestSurfElem = GeomToolsClosestElemNum(IsoSurfGradPath->IsoTopoZone, 
																	RingCP.X, RingCP.Y, RingCP.Z, 
																	&SurfPT.X, &SurfPT.Y, &SurfPT.Z);
				IsOk = (ClosestSurfElem > 0);
			}
			if (IsOk)
			{
				//	Now need to find the intersection point of the ring-atom GP and the surface.
				//	Do this by moving across the grad path, starting at the end (atom), and checking
				//	two adjacent points' distance to the previously found surface point. When the first
				//	point's distance becomes less than the second, or when the first's continues to
				//	decrease while the second's increases, then the two points are straddling the
				//	surface. Interpolate to find the intersection point.
				XYZ_s		LnBeg, LnEnd, Pt;
				double		Rho1, Rho2;
				LgIndex_t	NumOfVolGPPoints = GradPathGetCount(RingAtomGP);
				if (NumOfVolGPPoints > 0)
				{
					double	RSquare1, RSquare2, Min1, Min2;
					Boolean_t SurfFound = FALSE;
					for (LgIndex_t j = NumOfVolGPPoints - 1 ; j > 0 && !SurfFound ; j--)
					{
						GradPathGetPoint(RingAtomGP, j, 
							&LnBeg.X, 
							&LnBeg.Y, 
							&LnBeg.Z, 
							&Rho1);
						GradPathGetPoint(RingAtomGP, j-1, 
							&LnEnd.X, 
							&LnEnd.Y, 
							&LnEnd.Z, 
							&Rho2);
						RSquare1 = DistanceSquaredXYZ(LnBeg, SurfPT);
						RSquare2 = DistanceSquaredXYZ(LnEnd, SurfPT);

						if (j < NumOfVolGPPoints - 1 && RSquare2 > Min2 && RSquare1 < Min1)
							SurfFound = TRUE;
						else if (j < NumOfVolGPPoints - 1 && RSquare1 <= RSquare2)
						{
							j++;
							GradPathGetPoint(RingAtomGP, j, 
								&LnBeg.X, 
								&LnBeg.Y, 
								&LnBeg.Z, 
								&Rho1);
							GradPathGetPoint(RingAtomGP, j-1, 
								&LnEnd.X, 
								&LnEnd.Y, 
								&LnEnd.Z, 
								&Rho2);
							SurfFound = TRUE;
						}
						else
						{
							//	Still approaching surface, so reset Min variables
							//	and continue
							Min1 = RSquare1;
							Min2 = RSquare2;
						}

						if (SurfFound)
						{
							//	Interpolate to find intersection point.
							double DistanceRatio = RSquare1 / (RSquare1 + RSquare2);
							Pt.X = LnBeg.X + DistanceRatio * (LnEnd.X - LnBeg.X);
							Pt.Y = LnBeg.Y + DistanceRatio * (LnEnd.Y - LnBeg.Y);
							Pt.Z = LnBeg.Z + DistanceRatio * (LnEnd.Z - LnBeg.Z);

							//	Quadratic interpolation to find Rho at the point.
							double	Rho0, b1, b2;
							XYZ_s	Ln0;
							GradPathGetPoint(RingAtomGP, j+1,
								&Ln0.X,
								&Ln0.Y,
								&Ln0.Z,
								&Rho0);
							b1 = (Rho1 - Rho0) / sqrt(DistanceSquaredXYZ(Ln0, LnBeg));
							b2 = ((Rho2 - Rho1) / sqrt(DistanceSquaredXYZ(LnBeg, LnEnd)) - b1) 
												/ sqrt(DistanceSquaredXYZ(LnEnd, Ln0));
							double TempDist = sqrt(DistanceSquaredXYZ(Ln0, Pt));
							Rho = Rho0 + b1 * TempDist + b2 * TempDist
											* sqrt(DistanceSquaredXYZ(LnBeg, Pt));

							double SurfX, SurfY, SurfZ;
							LgIndex_t SurfElem;

							SurfElem = GeomToolsClosestElemNum(SurfCritPoints->SourceZoneNum, 
																Pt.X, Pt.Y, Pt.Z,
																&SurfX, &SurfY, &SurfZ);
							if (SurfElem > 0)
							{
								Pt.X = SurfX;
								Pt.Y = SurfY;
								Pt.Z = SurfZ;
							}
							else IsOk = FALSE;
							
							IntersectionPTs.push_back(Pt);
							IntersectionRhoValues.push_back(Rho);
						}
					}
					if (!SurfFound) IsOk = FALSE;
				}
				else IsOk = FALSE;
			}
		}
	}

	//	Now we have all the Rings that correspond to the given atom, and the intersection points
	//	of the ring-atom GP with the atom's isosurface. 
	//	Next, loop over the surface saddles to find the closest to an intersection point.
	//	When found, check the four (4) surf grad paths 
	size_t	NumOfIntPts = IntersectionPTs.size();
	for (LgIndex_t IntPt = 0 ; IsOk && IntPt < NumOfIntPts ; IntPt++)
	{
		LgIndex_t	ClosestSaddleNum = -1;
		double		MinDistSqr, CheckDistSqr, djunk;
		char		cjunk;
		XYZ_s		CheckPt;
		for (LgIndex_t SurfSaddle = SaddleBegOffset ; IsOk && SurfSaddle < SaddleEndOffset ; SurfSaddle++)
		{
			IsOk = CritPointsGetPoint(SurfCritPoints, SurfSaddle, 
									&CheckPt.X, &CheckPt.Y, &CheckPt.Z, 
									&djunk, &cjunk, &djunk, &djunk, &djunk);
			if (IsOk && SurfSaddle > SaddleBegOffset)
			{
				CheckDistSqr = DistanceSquaredXYZ(CheckPt, IntersectionPTs[IntPt]);
				if (CheckDistSqr < MinDistSqr)
				{
					MinDistSqr = CheckDistSqr;
					ClosestSaddleNum = SurfSaddle;
				}
			}
			else
			{
				MinDistSqr = DistanceSquaredXYZ(CheckPt, IntersectionPTs[IntPt]);
				ClosestSaddleNum = SurfSaddle;
			}
		}
		IsOk = (ClosestSaddleNum > 0);
		//CHECK(sqrt(MinDistSqr) < 0.5 * SurfCritPoints->MinCPDistance);
		
		//	Closest surf saddle found. Now get GPs starting with that saddle,
		//	eliminate outlying points and fix beginning point.
		if (IsOk)
		{	
			LgIndex_t	GPCount = 0;
			LgIndex_t	BegOffset = CritPointsGetBegOffset(SurfCritPoints, MaxType);
			LgIndex_t	EndOffset = CritPointsGetEndOffset(SurfCritPoints, MinType);
			CHECK(BegOffset >= 0 && EndOffset > 0);

			for (LgIndex_t SurfCPOffset = BegOffset ; IsOk && SurfCPOffset < EndOffset && GPCount < 4 ; SurfCPOffset++)
			{
				if (SurfCPOffset < SaddleBegOffset || SurfCPOffset >= SaddleEndOffset)
				{
					GradPath_pa	SurfGP = NULL;
					SurfGP = IsoSurfGradPathGetGPByBegEndCP(IsoSurfGradPath, ClosestSaddleNum, SurfCPOffset);

					if (SurfGP != NULL && GradPathIsValid(SurfGP))
					{
						//	Removing surf grad path points that are farther away from the end GP point than
						//	the intersection point found earlier
						double		TotalLenSqr = 0.0;
						XYZ_s		EndPt;
						int			SuccessNum = 0;

						IsOk = GradPathGetPoint(SurfGP, GradPathGetCount(SurfGP)-1, 
												&EndPt.X, &EndPt.Y, &EndPt.Z, &djunk);

						if (IsOk) TotalLenSqr = DistanceSquaredXYZ(IntersectionPTs[IntPt], EndPt);

						IsOk = (TotalLenSqr > 0);

						for (LgIndex_t CheckPtOffset = 0 ; IsOk && SuccessNum < 2 
							&& CheckPtOffset < GradPathGetCount(SurfGP) ; CheckPtOffset++)
						{
							GradPathGetPoint(SurfGP, CheckPtOffset, &CheckPt.X, &CheckPt.Y, &CheckPt.Z, &djunk);

							CheckDistSqr = DistanceSquaredXYZ(CheckPt, EndPt);

							if (CheckDistSqr < TotalLenSqr) SuccessNum++;
							else
							{
								GradPathRemovePoint(SurfGP, CheckPtOffset);
								CheckPtOffset--;
							}
						}

						//	Now the intersection point can be appended as the new first
						//	point of the surf grad path

						IsOk = GradPathInsertPoint(SurfGP, (LgIndex_t)0, IntersectionPTs[IntPt].X, 
															IntersectionPTs[IntPt].Y, 
															IntersectionPTs[IntPt].Z, 
															IntersectionRhoValues[IntPt]);
						GPCount++;
					}
				}
			}
		}
	}
	
	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}