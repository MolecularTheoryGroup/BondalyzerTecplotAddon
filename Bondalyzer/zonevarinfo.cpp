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
#include <string.h>
#include <fstream>
#include <vector>
#include <string>


SmInteger_t VarNumByNameList(std::vector<std::string> &VarNameList)
{
	SmInteger_t VarNum = -1;

	EntIndex_t NumberOfVars = TecUtilDataSetGetNumVars();

	for (EntIndex_t CurrentVar = 1; CurrentVar <= NumberOfVars && VarNum < 0; ++CurrentVar){
		VarName_t CurrentVarName;
		if (TecUtilVarGetName(CurrentVar, &CurrentVarName)){
			for (int CurrentNameCheck = 0; CurrentNameCheck < VarNameList.size() && VarNum < 0; ++CurrentNameCheck){
				if (!strncmp(CurrentVarName, VarNameList[CurrentNameCheck].c_str(), VarNameList[CurrentNameCheck].size())){
					VarNum = CurrentVar;
				}
			}
			TecUtilStringDealloc(&CurrentVarName);
		}
	}

	return VarNum;
}

SmInteger_t VarNumByName(std::string &VarName)
{
	SmInteger_t VarNum = -1;

	EntIndex_t NumberOfVars = TecUtilDataSetGetNumVars();

	for (EntIndex_t CurrentVar = 1; CurrentVar <= NumberOfVars && VarNum < 0; ++CurrentVar){
		VarName_t CurrentVarName;
		if (TecUtilVarGetName(CurrentVar, &CurrentVarName)){
			if (!strncmp(CurrentVarName, VarName.c_str(), VarName.size())){
				VarNum = CurrentVar;
			}
			TecUtilStringDealloc(&CurrentVarName);
		}
	}

	return VarNum;
}

SmInteger_t ZoneNumByNameList(std::vector<std::string> &ZoneNameList)
{
	SmInteger_t ZoneNum = -1;

	EntIndex_t NumberOfZones = TecUtilDataSetGetNumZones();

	for (EntIndex_t CurrentVar = 1; CurrentVar <= NumberOfZones && ZoneNum < 0; ++CurrentVar){
		VarName_t CurrentZoneName;
		if (TecUtilZoneGetName(CurrentVar, &CurrentZoneName)){
			for (int CurrentNameCheck = 0; CurrentNameCheck < ZoneNameList.size() && ZoneNum < 0; ++CurrentNameCheck){
				if (!strncmp(CurrentZoneName, ZoneNameList[CurrentNameCheck].c_str(), ZoneNameList[CurrentNameCheck].size())){
					ZoneNum = CurrentVar;
				}
			}
			TecUtilStringDealloc(&CurrentZoneName);
		}
	}

	return ZoneNum;
}

SmInteger_t ZoneNumByName(std::string &ZoneName)
{
	SmInteger_t ZoneNum = -1;

	EntIndex_t NumberOfZones = TecUtilDataSetGetNumZones();

	for (EntIndex_t CurrentVar = 1; CurrentVar <= NumberOfZones && ZoneNum < 0; ++CurrentVar){
		VarName_t CurrentZoneName;
		if (TecUtilZoneGetName(CurrentVar, &CurrentZoneName)){
			if (!strncmp(CurrentZoneName, ZoneName.c_str(), ZoneName.size())){
				ZoneNum = CurrentVar;
			}
			TecUtilStringDealloc(&CurrentZoneName);
		}
	}

	return ZoneNum;
}


/**
 * Verify that the ZoneVarInfo structure is valid.
 *
 * param ZoneVarInfo
 *     ZoneVarInfo structure being tested.
 *
 * return
 *     TRUE if valid, otherwise FALSE.
 */

Boolean_t ZoneVarInfoIsValid(ZoneVarInfo_pa ZoneVarInfo)
{
    Boolean_t IsOk = TRUE;
    EntIndex_t NumZones, NumVars;
    REQUIRE(VALID_REF(ZoneVarInfo));
    
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    IsOk = ( VALID_REF(ZoneVarInfo) && 
        ZoneVarInfo->ZoneNum >= 0 && ZoneVarInfo->ZoneNum <= NumZones &&
        ZoneVarInfo->ChrgDensVarNum >= 0 && ZoneVarInfo->ChrgDensVarNum <= NumVars &&
        ZoneVarInfo->UVarNum >= 0 && ZoneVarInfo->UVarNum <= NumVars &&
        ZoneVarInfo->VVarNum >= 0 && ZoneVarInfo->VVarNum <= NumVars &&
        ZoneVarInfo->WVarNum >= 0 && ZoneVarInfo->WVarNum <= NumVars &&
        ZoneVarInfo->TypeVarNum >= 0 && ZoneVarInfo->TypeVarNum <= NumVars &&
        (ZoneVarInfo->RefinedZoneNums == NULL || (VALID_REF(ZoneVarInfo->RefinedZoneNums) && TecUtilSetGetMemberCount(ZoneVarInfo->RefinedZoneNums) <= NumZones) ) &&
        (ZoneVarInfo->RefinedZoneXBeg == NULL || ArrListIsValid(ZoneVarInfo->RefinedZoneXBeg) ) &&
        (ZoneVarInfo->RefinedZoneXEnd == NULL || ArrListIsValid(ZoneVarInfo->RefinedZoneXEnd) ) &&
        (ZoneVarInfo->RefinedZoneYBeg == NULL || ArrListIsValid(ZoneVarInfo->RefinedZoneYBeg) ) &&
        (ZoneVarInfo->RefinedZoneYEnd == NULL || ArrListIsValid(ZoneVarInfo->RefinedZoneYEnd) ) &&
        (ZoneVarInfo->RefinedZoneZBeg == NULL || ArrListIsValid(ZoneVarInfo->RefinedZoneZBeg) ) &&
        (ZoneVarInfo->RefinedZoneZEnd == NULL || ArrListIsValid(ZoneVarInfo->RefinedZoneZEnd) ) &&
        (ZoneVarInfo->SurfElemMap == NULL || SurfElemMapIsValid(ZoneVarInfo->SurfElemMap) ) &&
        (ZoneVarInfo->Normals == NULL || NormalsIsValid(ZoneVarInfo->Normals) )
        );

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}

/**
 * Deallocates the ZoneVarInfo handle and set the handle to NULL.
 *
 *
 * param
 *     Reference to a ZoneVarInfo handle.
 */
void ZoneVarInfoDealloc(ZoneVarInfo_pa *ZoneVarInfo)
{
    REQUIRE(VALID_REF(ZoneVarInfo));
    REQUIRE(ZoneVarInfoIsValid(*ZoneVarInfo) || *ZoneVarInfo == NULL);

    if (*ZoneVarInfo != NULL)
    {
        if ((*ZoneVarInfo)->RefinedZoneNums != NULL) TecUtilSetDealloc(&((*ZoneVarInfo)->RefinedZoneNums));

        /* release the ArrList's */
        if ((*ZoneVarInfo)->RefinedZoneXBeg != NULL) ArrListDealloc(&((*ZoneVarInfo)->RefinedZoneXBeg));
        if ((*ZoneVarInfo)->RefinedZoneXEnd != NULL) ArrListDealloc(&((*ZoneVarInfo)->RefinedZoneXEnd));
        if ((*ZoneVarInfo)->RefinedZoneYBeg != NULL) ArrListDealloc(&((*ZoneVarInfo)->RefinedZoneYBeg));
        if ((*ZoneVarInfo)->RefinedZoneYEnd != NULL) ArrListDealloc(&((*ZoneVarInfo)->RefinedZoneYEnd));
        if ((*ZoneVarInfo)->RefinedZoneZBeg != NULL) ArrListDealloc(&((*ZoneVarInfo)->RefinedZoneZBeg));
        if ((*ZoneVarInfo)->RefinedZoneZEnd != NULL) ArrListDealloc(&((*ZoneVarInfo)->RefinedZoneZEnd));
        if ((*ZoneVarInfo)->RefinedZoneNums != NULL) TecUtilSetDealloc(&((*ZoneVarInfo)->RefinedZoneNums));

        if ((*ZoneVarInfo)->SurfElemMap != NULL) SurfElemMapDealloc(&((*ZoneVarInfo)->SurfElemMap));
        if ((*ZoneVarInfo)->Normals != NULL) NormalsDealloc(&((*ZoneVarInfo)->Normals));

        /* release the list structure itself */
        FREE_ITEM(*ZoneVarInfo, "ZoneVarInfo structure");
        *ZoneVarInfo = NULL;
    }

    ENSURE(*ZoneVarInfo == NULL);
}






/**
 * Empties the ZoneVarInfo structure and resets values to zero.
 *
 *
 * param ZoneVarInfo
 *     ZoneVarInfo structure to clear.
 */
void ZoneVarInfoClear(ZoneVarInfo_pa ZoneVarInfo)
{
    REQUIRE(ZoneVarInfoIsValid(ZoneVarInfo));

    ZoneVarInfo->ZoneNum   = 0;
    ZoneVarInfo->ChrgDensVarNum = 0;
    ZoneVarInfo->UVarNum = 0;
    ZoneVarInfo->VVarNum = 0;
    ZoneVarInfo->WVarNum = 0;
    ZoneVarInfo->TypeVarNum = 0;
    ZoneVarInfo->LastElemUsed = 0;

    TecUtilSetClear(ZoneVarInfo->RefinedZoneNums);

    ArrListClear(ZoneVarInfo->RefinedZoneXBeg);
    ArrListClear(ZoneVarInfo->RefinedZoneXEnd);
    ArrListClear(ZoneVarInfo->RefinedZoneYBeg);
    ArrListClear(ZoneVarInfo->RefinedZoneYEnd);
    ArrListClear(ZoneVarInfo->RefinedZoneZBeg);
    ArrListClear(ZoneVarInfo->RefinedZoneZEnd);

    if (SurfElemMapIsValid(ZoneVarInfo->SurfElemMap)) SurfElemMapClear(ZoneVarInfo->SurfElemMap);
    if (NormalsIsValid(ZoneVarInfo->Normals)) NormalsClear(ZoneVarInfo->Normals);

    ENSURE(ZoneVarInfoIsValid(ZoneVarInfo));
}


/**
 * Allocates a ZoneVarInfo handle with a suitable default
 * capacity for the contained ArrList's.
 *
 *
 * return
 *     ZoneVarInfo handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
ZoneVarInfo_pa ZoneVarInfoAlloc()
{
    ZoneVarInfo_pa Result = NULL;

    Result = ALLOC_ITEM(ZoneVarInfo_s, "ZoneVarInfo structure");
    if (Result != NULL)
    {
        Result->ZoneNum         = 0;
        Result->UVarNum         = 0;
        Result->VVarNum         = 0;
        Result->WVarNum         = 0;
        Result->ChrgDensVarNum  = 0;
        Result->TypeVarNum      = 0;
        Result->PeriodicBC      = FALSE;

        Result->XBeg            = 0.0;
        Result->XEnd            = 0.0;
        Result->YBeg            = 0.0;
        Result->YEnd            = 0.0;
        Result->ZBeg            = 0.0;
        Result->ZEnd            = 0.0;

        Result->RefinedZoneNums = NULL;
        Result->RefinedZoneXBeg = NULL;
        Result->RefinedZoneXEnd = NULL;
        Result->RefinedZoneYBeg = NULL;
        Result->RefinedZoneYEnd = NULL;
        Result->RefinedZoneZBeg = NULL;
        Result->RefinedZoneZEnd = NULL;

        Result->LastElemUsed    = 0;
        Result->SurfElemMap     = NULL;
        Result->Normals         = NULL;
    }

    ENSURE(ZoneVarInfoIsValid(Result) || Result == NULL);
    return Result;
}





/**
 * Load the coordinate ranges for the base and refined zones
 * into the ZoneVarInfo structure.
 *
 * param ZoneVarInfo
 *     ZoneVarInfo structure whose coordinate ranges are to be populated.
 *
 * param RefinedZoneNums
 *     Tecplot set of zone numbers defining refined-grid zones. This list
 *     should be copied into the ZoneVarInfo list.
 *
 *
 * NOTE: This function assumes that all zones are rectangular and aligned
 *       with the X, Y, Z coordinates. This greatly reduces the cost of
 *       verifying you are in a zone.
 *
 * return
 *     TRUE if successful, otherwise FALSE.
 */
Boolean_t ZoneVarInfoSetMinMax(ZoneVarInfo_pa ZoneVarInfo,
                               Set_pa         RefinedZoneNums)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t NumRefinedZones = 0;
    EntIndex_t NumZones, NumVars;
    EntIndex_t XVarNum, YVarNum, ZVarNum;

    REQUIRE(ZoneVarInfoIsValid(ZoneVarInfo));
    REQUIRE(ZoneVarInfo->ZoneNum > 0);

    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    NumRefinedZones = (LgIndex_t)TecUtilSetGetMemberCount(RefinedZoneNums);

    // Allocate the Set in the ZoneVarInfo structure for the refined-grid zone 
    // numbers and copy from the given set
    if (ZoneVarInfo->RefinedZoneNums != NULL) TecUtilSetDealloc(&(ZoneVarInfo->RefinedZoneNums));
    ZoneVarInfo->RefinedZoneNums = TecUtilSetAlloc(FALSE);
    IsOk = TecUtilSetCopy(ZoneVarInfo->RefinedZoneNums, RefinedZoneNums, FALSE);


    // Set up the array-lists for coordinate ranges of refined zones
    if (IsOk)
    {
        if (ArrListIsValid(ZoneVarInfo->RefinedZoneXBeg)) ArrListDealloc(&(ZoneVarInfo->RefinedZoneXBeg));
        ZoneVarInfo->RefinedZoneXBeg = ArrListAlloc(NumRefinedZones, ArrListType_Double);

        if (ArrListIsValid(ZoneVarInfo->RefinedZoneXEnd)) ArrListDealloc(&(ZoneVarInfo->RefinedZoneXEnd));
        ZoneVarInfo->RefinedZoneXEnd = ArrListAlloc(NumRefinedZones, ArrListType_Double);

        if (ArrListIsValid(ZoneVarInfo->RefinedZoneYBeg)) ArrListDealloc(&(ZoneVarInfo->RefinedZoneYBeg));
        ZoneVarInfo->RefinedZoneYBeg = ArrListAlloc(NumRefinedZones, ArrListType_Double);

        if (ArrListIsValid(ZoneVarInfo->RefinedZoneYEnd)) ArrListDealloc(&(ZoneVarInfo->RefinedZoneYEnd));
        ZoneVarInfo->RefinedZoneYEnd = ArrListAlloc(NumRefinedZones, ArrListType_Double);

        if (ArrListIsValid(ZoneVarInfo->RefinedZoneZBeg)) ArrListDealloc(&(ZoneVarInfo->RefinedZoneZBeg));
        ZoneVarInfo->RefinedZoneZBeg = ArrListAlloc(NumRefinedZones, ArrListType_Double);

        if (ArrListIsValid(ZoneVarInfo->RefinedZoneZEnd)) ArrListDealloc(&(ZoneVarInfo->RefinedZoneZEnd));
        ZoneVarInfo->RefinedZoneZEnd = ArrListAlloc(NumRefinedZones, ArrListType_Double);

        if (!ArrListIsValid(ZoneVarInfo->RefinedZoneXBeg) || !ArrListIsValid(ZoneVarInfo->RefinedZoneXEnd) ||
            !ArrListIsValid(ZoneVarInfo->RefinedZoneYBeg) || !ArrListIsValid(ZoneVarInfo->RefinedZoneYEnd) ||
            !ArrListIsValid(ZoneVarInfo->RefinedZoneZBeg) || !ArrListIsValid(ZoneVarInfo->RefinedZoneZEnd) ) IsOk = FALSE;
    }

    // Get X, Y, Z variable numbers
    if (IsOk)
    {
        XVarNum = TecUtilVarGetNumByAssignment('X');  
        YVarNum = TecUtilVarGetNumByAssignment('Y');  
        ZVarNum = TecUtilVarGetNumByAssignment('Z');
        if (XVarNum < 1 || XVarNum > NumVars ||
            YVarNum < 1 || YVarNum > NumVars || YVarNum == XVarNum ||
            ZVarNum < 1 || ZVarNum > NumVars || ZVarNum == XVarNum || ZVarNum == YVarNum) IsOk = FALSE;
    }

    // NOTE: Following assumes rectangular zones alligned with coordinates

    // Get X, Y, Z data ranges for base zone
    if (IsOk)
    {
        LgIndex_t IMax, JMax, KMax, IJKMax;

        FieldData_pa XFDPtr = TecUtilDataValueGetReadableNLRef(ZoneVarInfo->ZoneNum, XVarNum);
        FieldData_pa YFDPtr = TecUtilDataValueGetReadableNLRef(ZoneVarInfo->ZoneNum, YVarNum);
        FieldData_pa ZFDPtr = TecUtilDataValueGetReadableNLRef(ZoneVarInfo->ZoneNum, ZVarNum);

        TecUtilZoneGetInfo(ZoneVarInfo->ZoneNum, &IMax, &JMax, &KMax, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);


        IJKMax = IMax * JMax * KMax;

		
        ZoneVarInfo->XBeg = TecUtilDataValueGetByRef(XFDPtr, 1);
        ZoneVarInfo->XEnd = TecUtilDataValueGetByRef(XFDPtr, IJKMax);

        ZoneVarInfo->YBeg = TecUtilDataValueGetByRef(YFDPtr, 1);
        ZoneVarInfo->YEnd = TecUtilDataValueGetByRef(YFDPtr, IJKMax);

        ZoneVarInfo->ZBeg = TecUtilDataValueGetByRef(ZFDPtr, 1);
        ZoneVarInfo->ZEnd = TecUtilDataValueGetByRef(ZFDPtr, IJKMax);
    }

    // Get X, Y, Z data ranges for each refined zone
    if (IsOk)
    {
        LgIndex_t IMax, JMax, KMax, IJKMax;
        SetIndex_t NumRefinedZones = TecUtilSetGetMemberCount(ZoneVarInfo->RefinedZoneNums);
        
        if (NumRefinedZones > 0)
        {
            LgIndex_t sp;
            EntIndex_t ZoneNum;
            double XBeg, XEnd, YBeg, YEnd, ZBeg, ZEnd;
            ArrListItem_u Item;

            for (sp = 1; sp <= NumRefinedZones; sp++)
            {
                ZoneNum = (EntIndex_t)TecUtilSetGetMember(ZoneVarInfo->RefinedZoneNums, sp); 

                FieldData_pa XFDPtr = TecUtilDataValueGetReadableNLRef(ZoneNum, XVarNum);
                FieldData_pa YFDPtr = TecUtilDataValueGetReadableNLRef(ZoneNum, YVarNum);
                FieldData_pa ZFDPtr = TecUtilDataValueGetReadableNLRef(ZoneNum, ZVarNum);

                TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);


                IJKMax = IMax * JMax * KMax;

                XBeg = TecUtilDataValueGetByRef(XFDPtr, 1);
                XEnd = TecUtilDataValueGetByRef(XFDPtr, IJKMax);

                YBeg = TecUtilDataValueGetByRef(YFDPtr, 1);
                YEnd = TecUtilDataValueGetByRef(YFDPtr, IJKMax);

                ZBeg = TecUtilDataValueGetByRef(ZFDPtr, 1);
                ZEnd = TecUtilDataValueGetByRef(ZFDPtr, IJKMax);

                Item.Double = XBeg;
                IsOk = ArrListSetItem(ZoneVarInfo->RefinedZoneXBeg, sp-1, Item);

                if (IsOk)
                {
                    Item.Double = XEnd;
                    IsOk = ArrListSetItem(ZoneVarInfo->RefinedZoneXEnd, sp-1, Item);
                }

                if (IsOk)
                {
                    Item.Double = YBeg;
                    IsOk = ArrListSetItem(ZoneVarInfo->RefinedZoneYBeg, sp-1, Item);
                }

                if (IsOk)
                {
                    Item.Double = YEnd;
                    IsOk = ArrListSetItem(ZoneVarInfo->RefinedZoneYEnd, sp-1, Item);
                }

                if (IsOk)
                {
                    Item.Double = ZBeg;
                    IsOk = ArrListSetItem(ZoneVarInfo->RefinedZoneZBeg, sp-1, Item);
                }

                if (IsOk)
                {
                    Item.Double = ZEnd;
                    IsOk = ArrListSetItem(ZoneVarInfo->RefinedZoneZEnd, sp-1, Item);
                }
            }
        }
    }    

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}







/**
 * Compute and return the grid spacing for the base zone.
 *
 * param ZoneVarInfo
 *     ZoneVarInfo structure.
 *
 *
 * NOTE: This function assumes that all zones are rectangular and aligned
 *       with the X, Y, Z coordinates. This greatly reduces the cost of
 *       verifying you are in a zone.
 *
 * return
 *     Grid spacing if it worked, otherwise 1.0.
 */
double ZoneVarInfoGetDXForBaseZone(ZoneVarInfo_pa ZoneVarInfo)
{
    double GridSpacing = 1.0;

    EntIndex_t ZoneNum = ZoneVarInfo->ZoneNum;
    LgIndex_t  IMax, JMax, KMax;

    REQUIRE(ZoneVarInfoIsValid(ZoneVarInfo));
    REQUIRE(ZoneNum > 0);
    REQUIRE(TecUtilZoneIsOrdered(ZoneNum));

    // Use NULL for values we're not interested in
    TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, NULL,NULL,NULL,NULL,NULL,
                      NULL,NULL,NULL,NULL,NULL);

    REQUIRE(IMax > 1 && JMax > 1 && KMax > 1);

    GridSpacing = 0.333 * ( ABS(ZoneVarInfo->XEnd - ZoneVarInfo->XBeg) / (double)(IMax - 1) +
                            ABS(ZoneVarInfo->YEnd - ZoneVarInfo->YBeg) / (double)(JMax - 1) +
                            ABS(ZoneVarInfo->ZEnd - ZoneVarInfo->ZBeg) / (double)(KMax - 1) );

    if (GridSpacing < SMALLFLOAT) GridSpacing = 1.0;

    ENSURE(GridSpacing > 0.0);
    return GridSpacing;
}







/**
 * Given a zone number, determine the grid spacing.
 *
 * param ZoneVarInfo
 *     ZoneVarInfo structure.
 *
 * param ZoneNum
 *     The zone number from which to get the grid spacing.
 *
 *
 * NOTE: This function assumes that all zones are rectangular and aligned
 *       with the X, Y, Z coordinates. This greatly reduces the cost of
 *       verifying you are in a zone.
 *
 * return
 *     Grid spacing if it worked, otherwise 1.0.
 */
double ZoneVarInfoGetDXFromZoneNum(ZoneVarInfo_pa ZoneVarInfo,
                                   LgIndex_t      ZoneNum)
{
    double GridSpacing = 1.0;
    
    REQUIRE(ZoneVarInfoIsValid(ZoneVarInfo));
    REQUIRE(ZoneNum >= 1);

    // Find which zone (base zone, or offset to fine zone) is ZoneNum
    if (ZoneNum == ZoneVarInfo->ZoneNum)
    {
        GridSpacing = ZoneVarInfoGetDXForBaseZone(ZoneVarInfo);
    }
    else
    {
        // TODO: find spacing for refined grids
        CHECK("Refined Grid Spacing Not Implemented" == 0);
    }

    ENSURE(GridSpacing > 0.0);
    return GridSpacing;
}
