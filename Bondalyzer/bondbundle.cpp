/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "GEOMTOOLS.h"
#include "SURFELEMMAP.h"
#include "NORMALS.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
#include "GRADPATH.h"
#include "BUNDLES.h"
#include "BONDBUNDLE.h"
// #include "ENGINE.h"



/**
 * Determine if the bondBundle handle is sane.
 *
 * param bondBundle
 *     bondBundle structure in question.
 *
 * return
 *     TRUE if the bondBundle structure is valid, otherwise FALSE.
 */
Boolean_t BondBundleIsValid(const BondBundle_pa bondBundle)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(bondBundle) &&
               VALID_REF(bondBundle->BondRingCageSurfZnNums) && ArrListIsValid(bondBundle->BondRingCageSurfZnNums) &&
               VALID_REF(bondBundle->PlsAtomRingCageSurfZnNums) && ArrListIsValid(bondBundle->PlsAtomRingCageSurfZnNums) &&
               VALID_REF(bondBundle->MnsAtomRingCageSurfZnNums) && ArrListIsValid(bondBundle->MnsAtomRingCageSurfZnNums));


    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}












/**
 * Deallocates the bondBundle handle and set the handle to NULL.
 *
 * param
 *     Reference to a bondBundle handle.
 */
void BondBundleDealloc(BondBundle_pa *bondBundle)
{

    REQUIRE(VALID_REF(bondBundle));

    if (*bondBundle != NULL)
    {
        /* release the ArrList's */
        if ((*bondBundle)->BondRingCageSurfZnNums != NULL) ArrListDealloc(&((*bondBundle)->BondRingCageSurfZnNums));
        if ((*bondBundle)->PlsAtomRingCageSurfZnNums != NULL) ArrListDealloc(&((*bondBundle)->PlsAtomRingCageSurfZnNums));
        if ((*bondBundle)->MnsAtomRingCageSurfZnNums != NULL) ArrListDealloc(&((*bondBundle)->MnsAtomRingCageSurfZnNums));

        /* release the list structure itself */
        FREE_ITEM(*bondBundle, "bondBundle structure");
        *bondBundle = NULL;
    }

    ENSURE(*bondBundle == NULL);
}




/**
 * Gets the number of SurfZnNums currently in the specified ArrList of the bondBundle
 * (maintained by the bondBundle array lists).
 *
 * param
 *     bondBundle structure in question.
 * param surfType
 *     SurfType_BondRingCage, SurfType_PlsAtomRingCage, or SurfType_MnsAtomRingCage
 *
 * return
 *     Number of SurfZnNums in the surfType ArrList in the bondBundle.
 */
LgIndex_t BondBundleGetSurfCount(const BondBundle_pa bondBundle,
                                 const SurfType_e    surfType)
{
    LgIndex_t Result = 0;
    Boolean_t isOk = TRUE;

    REQUIRE(BondBundleIsValid(bondBundle));
    REQUIRE(VALID_ENUM(surfType, SurfType_e));

    // Get surfZoneNum count for specified surface type
    if (isOk)
    {
        switch (surfType)
        {
            case SurfType_BondRingCage:
                Result = ArrListGetCount(bondBundle->BondRingCageSurfZnNums);
                break;
            case SurfType_PlsAtomRingCage:
                Result = ArrListGetCount(bondBundle->PlsAtomRingCageSurfZnNums);
                break;
            case SurfType_MnsAtomRingCage:
                Result = ArrListGetCount(bondBundle->MnsAtomRingCageSurfZnNums);
                break;
            default:
                CHECK(0);
                break;
        }
    }

    ENSURE(Result >= 0);
    return Result;
}







/**
 * Empties the bondBundle of all surface, nullifies the reference
 * to volume bundles, and resets the other variables.
 *
 *
 * param bondBundle
 *     bondBundle to clear.
 */
void BondBundleClear(BondBundle_pa bondBundle)
{
    REQUIRE(BondBundleIsValid(bondBundle));

    ArrListClear(bondBundle->BondRingCageSurfZnNums);
    ArrListClear(bondBundle->PlsAtomRingCageSurfZnNums);
    ArrListClear(bondBundle->MnsAtomRingCageSurfZnNums);


    bondBundle->BondNum              = -1;

    bondBundle->PlsAtomNum           = -1;
    bondBundle->MnsAtomNum           = -1;

    bondBundle->BaseVolZoneNum = -1;

    bondBundle->VolZoneNum = -1;  // Or, should we delete the zone?

    ENSURE(BondBundleIsValid(bondBundle) &&
           BondBundleGetSurfCount(bondBundle, SurfType_BondRingCage) == 0 &&
           BondBundleGetSurfCount(bondBundle, SurfType_PlsAtomRingCage) == 0 &&
           BondBundleGetSurfCount(bondBundle, SurfType_MnsAtomRingCage) == 0);
}







/**
 * Allocates a bondBundle handle with a suitable default
 * capacity for the contained ArrList's.
 *
 *
 * return
 *     bondBundle handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
BondBundle_pa BondBundleAlloc()
{
    BondBundle_pa Result = NULL;

    Result = ALLOC_ITEM(BondBundle_s, "bondBundle structure");
    if (Result != NULL)
    {
        Result->BondNum        = -1;
        Result->PlsAtomNum     = -1;
        Result->MnsAtomNum     = -1;

        Result->BondRingCageSurfZnNums = ArrListAlloc(20, ArrListType_Long);
        Result->PlsAtomRingCageSurfZnNums = ArrListAlloc(20, ArrListType_Long);
        Result->MnsAtomRingCageSurfZnNums = ArrListAlloc(20, ArrListType_Long);

        Result->BaseVolZoneNum = -1;

        Result->VolZoneNum = -1;
    }

    ENSURE(BondBundleIsValid(Result) || Result == NULL);
    return Result;
}






/**
 * Gets the surface zone number at the specified offset for the
 * specified surface type ArrList.
 *
 * param bondBundle
 *     bondBundle to query for surface zone number.
 * param surfType
 *     SurfType_BondRingCage, SurfType_PlsAtomRingCage, or SurfType_MnsAtomRingCage
 * param Offset
 *     Offset into array list.
 * param SurfZnNum
 *     Surface zone number to set at the specified offset in ArrList.
 *
 * return
 *     Surface zone number if it works, otherwise 0 (zero).
 */
EntIndex_t   BondBundleGetSurfZnNum(const BondBundle_pa  bondBundle,
                                    const SurfType_e     surfType,
                                    const LgIndex_t      Offset)
{
    Boolean_t  isOk = TRUE;
    EntIndex_t SurfZnNum = 0;
    ArrListItem_u Item;


    REQUIRE(BondBundleIsValid(bondBundle));
    REQUIRE(VALID_ENUM(surfType, SurfType_e));
    REQUIRE(0 <= Offset && Offset <= BondBundleGetSurfCount(bondBundle, surfType) - 1);

    // Get surfZoneNum for specified surface type
    if (isOk)
    {
        switch (surfType)
        {
            case SurfType_BondRingCage:
                Item = ArrListGetItem(bondBundle->BondRingCageSurfZnNums, Offset);
                break;
            case SurfType_PlsAtomRingCage:
                Item = ArrListGetItem(bondBundle->PlsAtomRingCageSurfZnNums, Offset);
                break;
            case SurfType_MnsAtomRingCage:
                Item = ArrListGetItem(bondBundle->MnsAtomRingCageSurfZnNums, Offset);
                break;
            default:
                CHECK(0);
                break;
        }
        SurfZnNum = (EntIndex_t)Item.Long;
    }

    ENSURE(BondBundleIsValid(bondBundle));

    ENSURE(SurfZnNum > -1);
    return SurfZnNum;
}






/**
 * Sets ArrList surface zone number at the specified offset for the
 * specified surface type. If the
 * offset is beyond the end of the array list, it is sized
 * accordingly and the intervening array values between the last
 * node of the original state and the last node of the new state are
 * assigned 0.
 *
 * note
 *     If SurfZnNum already exists at the specified location it is
 *     replaced.
 *
 * param bondBundle
 *     bondBundle target in which to set the Angle/GradPath.
 * param surfType
 *     SurfType_BondRingCage, SurfType_PlsAtomRingCage, or SurfType_MnsAtomRingCage
 * param Offset
 *     Offset into array list.
 * param SurfZnNum
 *     Surface zone number to set at the specified offset in ArrList.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t   BondBundleSetSurfZnNum(BondBundle_pa     bondBundle,
                                   const SurfType_e  surfType,
                                   const LgIndex_t   Offset,
                                   const EntIndex_t  SurfZnNum)
{
    Boolean_t isOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(BondBundleIsValid(bondBundle));
    REQUIRE(VALID_ENUM(surfType, SurfType_e));
    REQUIRE(Offset >= 0);

    // Set surfZoneNum for specified surface type
    if (isOk)
    {
        Item.Long = SurfZnNum;
        switch (surfType)
        {
            case SurfType_BondRingCage:
                isOk = ArrListSetItem(bondBundle->BondRingCageSurfZnNums, Offset, Item);
                break;
            case SurfType_PlsAtomRingCage:
                isOk = ArrListSetItem(bondBundle->PlsAtomRingCageSurfZnNums, Offset, Item);
                break;
            case SurfType_MnsAtomRingCage:
                isOk = ArrListSetItem(bondBundle->MnsAtomRingCageSurfZnNums, Offset, Item);
                break;
            default:
                CHECK(0);
                break;
        }
    }

    ENSURE(BondBundleIsValid(bondBundle));

    ENSURE(VALID_BOOLEAN(isOk));
    return isOk;
}





/**
 * Inserts a surface zone number at the specified offset for the ArrList
 * of the specified surface type. The array will be expanded to accomodate
 * the additional value.
 *
 *
 * param bondBundle
 *     bondBundle target in which to set the coordinates.
 * param surfType
 *     SurfType_BondRingCage, SurfType_PlsAtomRingCage, or SurfType_MnsAtomRingCage
 * param Offset
 *     Offset into array list.
 * param SurfZnNum
 *     Surface zone number to insert at the specified offset in ArrList.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t BondBundleInsertSurfZnNum(BondBundle_pa     bondBundle,
                                    const SurfType_e  surfType,
                                    const LgIndex_t   Offset,
                                    const EntIndex_t  SurfZnNum)
{
    Boolean_t isOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(BondBundleIsValid(bondBundle));
    REQUIRE(VALID_ENUM(surfType, SurfType_e));
    REQUIRE(0 <= Offset && Offset <= BondBundleGetSurfCount(bondBundle, surfType));

    if (isOk)
    {
        Item.Long = SurfZnNum;
        switch (surfType)
        {
            case SurfType_BondRingCage:
                isOk = ArrListInsertItem(bondBundle->BondRingCageSurfZnNums, Offset, Item);
                break;
            case SurfType_PlsAtomRingCage:
                isOk = ArrListInsertItem(bondBundle->PlsAtomRingCageSurfZnNums, Offset, Item);
                break;
            case SurfType_MnsAtomRingCage:
                isOk = ArrListInsertItem(bondBundle->MnsAtomRingCageSurfZnNums, Offset, Item);
                break;
            default:
                CHECK(0);
                break;
        }
    }

    ENSURE(BondBundleIsValid(bondBundle));

    ENSURE(VALID_BOOLEAN(isOk));
    return isOk;
}







/**
 * Deletes surfZoneNum at the specified offset from the specified list..
 *
 *
 * param bondBundle
 *     bondBundle target in which to set the coordinates.
 * param surfType
 *     SurfType_BondRingCage, SurfType_PlsAtomRingCage, or SurfType_MnsAtomRingCage
 * param Offset
 *     Offset at which to remove the SurfZnNum.
 *
 * return
 *     TRUE if successful operation, otherwise FALSE.
 */
Boolean_t BondBundleRemoveSurfZnNum(BondBundle_pa     bondBundle,
                                    const SurfType_e  surfType,
                                    const LgIndex_t   Offset)
{
    Boolean_t isOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(BondBundleIsValid(bondBundle));
    REQUIRE(VALID_ENUM(surfType, SurfType_e));
    REQUIRE(0 <= Offset && Offset <= BondBundleGetSurfCount(bondBundle, surfType));

    if (isOk)
    {
        switch (surfType)
        {
            case SurfType_BondRingCage:
                Item = ArrListRemoveItem(bondBundle->BondRingCageSurfZnNums, Offset);
                break;
            case SurfType_PlsAtomRingCage:
                Item = ArrListRemoveItem(bondBundle->PlsAtomRingCageSurfZnNums, Offset);
                break;
            case SurfType_MnsAtomRingCage:
                Item = ArrListRemoveItem(bondBundle->MnsAtomRingCageSurfZnNums, Offset);
                break;
            default:
                CHECK(0);
                break;
        }
    }

    ENSURE(BondBundleIsValid(bondBundle));

    ENSURE(VALID_BOOLEAN(isOk));
    return isOk;
}






/**
 * Appends the SurfZnNum to the specified array lists in bondBundle. The
 * array list will be expanded to accommodate the additional items.
 *
 * param bondBundle
 *     bondBundle target to which the SurfZnNum is to be appended.
 * param surfType
 *     SurfType_BondRingCage, SurfType_PlsAtomRingCage, or SurfType_MnsAtomRingCage
 * param SurfZnNum
 *     Surface zone number to append to the specified ArrList.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t BondBundleAppendSurfZnNum(BondBundle_pa     bondBundle,
                                    const SurfType_e  surfType,
                                    const EntIndex_t  SurfZnNum)
{
    Boolean_t isOk = TRUE;
    LgIndex_t Count;

    REQUIRE(BondBundleIsValid(bondBundle));
    REQUIRE(VALID_ENUM(surfType, SurfType_e));

    Count = BondBundleGetSurfCount(bondBundle, surfType);

    isOk = BondBundleInsertSurfZnNum(bondBundle, surfType, Count, SurfZnNum);

    ENSURE(BondBundleIsValid(bondBundle));

    ENSURE(VALID_BOOLEAN(isOk));
    return isOk;
}











/*
 * BondBundleGetSurfs: For a given bond-saddle volume critical point, search
 *  for and add to bondBundle structure the zone numbers of the surfaces that
 *  make up the BondRingCage, PlsAtomRingCage, and MnsAtomRingCage surfaces of
 *  the bond bundle. For a bond PlsAtoms are in the direction ofthe principle
 *  eigenvector of the curvature tensor.
 *
 * param BondCrtPtNum
 *     Number (offset+1) of the Bond volume critical point that is surrounded by
 *     the isosurface.
 * param VolBundles
 *     Lists containing the Atom, Bond, Ring, Cage combinations that make up the
 *     irreducible bundles for the volume.
 * param bondBundle
 *     bondBundle data structure to populate for an Volume volume critical point.
 *
 * return
 *     TRUE if successful, FALSE if there were errors.
 */
/*
Boolean_t BondBundleGetSurfs(const LgIndex_t      BondCrtPtNum,
                             const Bundles_pa     VolBundles,
                             BondBundle_pa        bondBundle)
{
  Boolean_t      isOk = TRUE;

  LgIndex_t PlsAtom, MnsAtom;

  // The RingCage lists for the three surfaces. Unless the  bundle extends to the
  // outer boundary, these lists should all be the same.
  OrderedPairs_pa RingCageForBond = NULL;
  OrderedPairs_pa RingCageForPlsAtom = NULL;
  OrderedPairs_pa RingCageForMnsAtom = NULL;

  REQUIRE(BondCrtPtNum >= 0);
  REQUIRE(BundlesIsValid(VolBundles));

  // Allocate space for the RingCage lists
  RingCageForBond    = OrderedPairsAlloc();
  RingCageForPlsAtom = OrderedPairsAlloc();
  RingCageForMnsAtom = OrderedPairsAlloc();

  if (RingCageForBond == NULL || RingCageForPlsAtom == NULL || RingCageForMnsAtom == NULL) isOk = FALSE;


  // First, search the VolBundles structure to get the atom numbers for the two
  // atoms, and the ring-cage numbers for the outer edge of the surfaces.
  if (isOk)
    {
      LgIndex_t iboff;
      LgIndex_t NumIRBundles = BundlesGetCount(VolBundles);

      for (iboff = 0; iboff <NumIRBundles; iboff++)
        {
          LgIndex_t Atom, Bond, Ring, Cage;

          isOk = BundlesGetFromOffset(VolBundles, iboff, &Atom, &Bond, &Ring, &Cage);
          if (Bond == BondCrtPtNum)
            {
              if (MnsAtom == 0 && PlsAtom == 0) PlsAtom = Atom;
              if (PlsAtom >0 && PlsAtom != Atom && MnsAtom == 0) MnsAtom = Atom;

              // Add to the RingCage lists for Bond
              isOk = OrderedPairsAppendAtEnd(RingCageForBond, Ring, Cage);

              // Add to the RingCage lists for PlsAtom (if appropriate).
              if (PlsAtom == Atom)
                {
                  isOk = OrderedPairsAppendAtEnd(RingCageForPlsAtom, Ring, Cage);
                }

              // Add to the RingCage lists for MnsAtom (if appropriate0
              if (MnsAtom == Atom)
                {
                  isOk = OrderedPairsAppendAtEnd(RingCageForMnsAtom, Ring, Cage);
                }
            }
        }
    }

  // Now search for surfaces that contain the Bond-Ring-Cage numbers corresponding to
  // BondCrtPtNum-RingCageForBond (populate the bondBundle->BondRingCageSurfZnNums list),
  // BondCrtPtNum-RingCageForPlsAtom (populate the
  // bondBundle->PlsAtomRingCageSurfZnNums list), and
  // BondCrtPtNum-RingCageForMnsAtom (populate the
  // bondBundle->MnsAtomRingCageSurfZnNums list), and
  if (isOk)
    {
      EntIndex_t ZnNum;
      EntIndex_t NumZones, NumVars;

      if (isOk)
          isOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

      for(ZnNum=1; isOk && ZnNum <= NumZones; ZnNum++)
        {
          Boolean_t IsSurface = FALSE;
          // Is ZnNum a surface zone?
          if (TecUtilZoneIsOrdered(ZnNum))
            {
              LgIndex_t baseZoneIMax, baseZoneJMax, baseZoneKMax;
              // Use NULL for values we're not interested in
              TecUtilZoneGetInfo(ZnNum, &baseZoneIMax, &baseZoneJMax, &baseZoneKMax, NULL,NULL,NULL,NULL,NULL,
                                NULL,NULL,NULL,NULL,NULL);
              if (baseZoneIMax > 1 && baseZoneJMax > 1 && baseZoneKMax == 1) isSurface = TRUE;
            }

          // Read Aux Data and compare to Bond-Ring-Cage or Atom-Ring-Cage combinations
          if (isOk && IsSurface)
            {
              AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(ZnNum);
              if (AuxDataRef != NULL)
                {
                  char      *Value;
                  Boolean_t  Retain;
                  LgIndex_t  BeginCrtPtNum, EndingCrtPtNum;  // need something else

                  // Get BegCrtPtNum
                  if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
                                                     &Value, &Retain))
                    {
                      if (Value != NULL)
                        {
                          // Assign the parsed value to the structure
                          long BegCrtPtNum;
                          if (sscanf(Value, "%d", &BegCrtPtNum) == 1)
                            {
                              BeginCrtPtNum = (LgIndex_t)(BegCrtPtNum-1);
                            }
                          else
                            {
                              isOk = FALSE;
                            }
                            // release the allocated string copy
                          TecUtilStringDealloc(&Value);
                        }
                      else
                        {
                          isOk = FALSE;
                        }
                    }

                  // Get EndCrtPtNum
                  if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
                                                     &Value, &Retain))
                    {
                      if (Value != NULL)
                        {
                          // Assign the parsed value to the structure
                          long EndCrtPtNum;
                          if (sscanf(Value, "%d", &EndCrtPtNum) == 1)
                            {
                              EndingCrtPtNum = (LgIndex_t)(EndCrtPtNum-1);
                            }
                          else
                            {
                              isOk = FALSE;
                            }
                            // release the allocated string copy
                          TecUtilStringDealloc(&Value);
                        }
                      else
                        {
                          isOk = FALSE;
                        }
                    }
                }
              else  // ZoneAuxDataRef is Null
                isOk = FALSE;
            }
        }
    }

  // Cleanup
  OrderedPairsDealloc(&RingCageForBond);
  OrderedPairsDealloc(&RingCageForPlsAtom);
  OrderedPairsDealloc(&RingCageForMnsAtom);

  ENSURE(VALID_BOOLEAN(isOk));
  return isOk;
}
*/


/*
 * BondBundleProcessSurfaceLineSegIntersections: Call the specified callback
 * function (with the specified clientData) for each intersection between
 * a given line segment and the bond bundle surfaces (which must be fully defined).
 *
 * param bondBundle
 *     bondBundle with surfaces defined
 * param lineSegStart
 *     XYZ point of start of line segment
 * param lineSegEnd
 *     XYZ point of end of line segment
 * param intersectionCallback
 *     callback function
 * param clientData
 *     callback client data
 *
 */
void BondBundleProcessSurfaceLineSegIntersections(BondBundle_pa const     bondBundle,
                                                  XYZ_s const &           lineSegStart,
                                                  XYZ_s const &           lineSegEnd,
                                                  IntersectionCallback_pf intersectionCallback,
                                                  void *                  clientData)
{
    REQUIRE(BondBundleIsValid(bondBundle));
    REQUIRE(bondBundle->BaseVolZoneNum > 0);
    REQUIRE(lineSegStart.X != lineSegEnd.X || lineSegStart.Y != lineSegEnd.Y || lineSegStart.Z != lineSegEnd.Z);
    REQUIRE(VALID_FN_REF(intersectionCallback));

    // Inform Tecplot that major data operation is beginning
    TecUtilDataLoadBegin();

    // Lock Tecplot
    TecUtilLockStart(AddOnID);

    EntIndex_t baseZoneNum = bondBundle->BaseVolZoneNum;
    CHECK(TecUtilZoneIsOrdered(baseZoneNum));
        // Use NULL for values we're not interested in
    LgIndex_t baseZoneIMax;
    LgIndex_t baseZoneJMax;
    LgIndex_t baseZoneKMax;
    TecUtilZoneGetInfo(baseZoneNum, &baseZoneIMax, &baseZoneJMax, &baseZoneKMax,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL);

    double const baseZoneXMin = 1.0;
    double const baseZoneXMax = (double)baseZoneIMax;
    double const baseZoneYMin = 1.0;
    double const baseZoneYMax = (double)baseZoneJMax;
    double const baseZoneZMin = 1.0;
    double const baseZoneZMax = (double)baseZoneKMax;

    EntIndex_t const xVarNum = TecUtilVarGetNumByAssignment('X');
    EntIndex_t const yVarNum = TecUtilVarGetNumByAssignment('Y');
    EntIndex_t const zVarNum = TecUtilVarGetNumByAssignment('Z');

    CHECK(GeomToolsTest());

    // Check both + and - atom-ring-cage surfaces
    for ( int pass = 0; pass <= 1; pass++ )
    {
        SurfType_e surfType = (pass==0) ? SurfType_PlsAtomRingCage : SurfType_MnsAtomRingCage;
        LgIndex_t numSurfs = BondBundleGetSurfCount(bondBundle, surfType);
        CHECK(numSurfs >= 0);

        // For each cell of each surface, test for an intersectionPt until one is found
        for (LgIndex_t ns = 0; ns < numSurfs; ns++)
        {
            EntIndex_t surfZoneNum = BondBundleGetSurfZnNum(bondBundle, surfType, ns);
            GeomToolsLineSegSurfaceIntersection(surfZoneNum, lineSegStart, lineSegEnd, intersectionCallback, clientData);

        // Special code for cases where domain-spanning surfaces don't go into the corners of
        // the BaseZone. Search around the outer boundary of the zone, add triangles in the
        // empty corners, and search the new triangles.
            // Since we are looking for every intersection, we cannot just skip this
            FieldData_pa xFDPtr = TecUtilDataValueGetReadableRef(surfZoneNum, xVarNum);
            FieldData_pa yFDPtr = TecUtilDataValueGetReadableRef(surfZoneNum, yVarNum);
            FieldData_pa zFDPtr = TecUtilDataValueGetReadableRef(surfZoneNum, zVarNum);

            CHECK(TecUtilZoneIsOrdered(surfZoneNum));
            LgIndex_t iMax, jMax, kMax;
                    // Use NULL for values we're not interested in
            TecUtilZoneGetInfo(surfZoneNum, &iMax, &jMax, &kMax,
                               NULL, NULL, NULL, NULL, NULL,
                                       NULL, NULL, NULL, NULL, NULL);
            CHECK(iMax >= 2 && jMax >= 2 && kMax == 1);

            for (LgIndex_t jj = 1; jj < jMax; jj++)
                {
                XYZ_s cornerImaxJ;
                LgIndex_t index = jj * iMax; // iMax,jj
                cornerImaxJ.X = TecUtilDataValueGetByRef(xFDPtr, index);
                cornerImaxJ.Y = TecUtilDataValueGetByRef(yFDPtr, index);
                cornerImaxJ.Z = TecUtilDataValueGetByRef(zFDPtr, index);

                XYZ_s cornerBoxEdgePt1;
                cornerBoxEdgePt1.X = 0.0;
                cornerBoxEdgePt1.Y = 0.0;
                cornerBoxEdgePt1.Z = 0.0;

                Boolean_t isVolIMinSurf = FALSE;
                Boolean_t isVolIMaxSurf = FALSE;
                Boolean_t isVolJMinSurf = FALSE;
                Boolean_t isVolJMaxSurf = FALSE;
                Boolean_t isVolKMinSurf = FALSE;
                Boolean_t isVolKMaxSurf = FALSE;

                if (cornerImaxJ.X > 0.99 && cornerImaxJ.X < 1.01)
                    {
                    isVolIMinSurf = TRUE;
                    cornerBoxEdgePt1.X = baseZoneXMin;
                    }
                else if (cornerImaxJ.X > baseZoneXMax - 0.01 && cornerImaxJ.X < baseZoneXMax + 0.01)
                    {
                    isVolIMaxSurf = TRUE;
                    cornerBoxEdgePt1.X = baseZoneXMax;
                    }

                if (cornerImaxJ.Y > 0.99 && cornerImaxJ.Y < 1.01)
                    {
                    isVolJMinSurf = TRUE;
                    cornerBoxEdgePt1.Y = baseZoneYMin;
                    }
                else if (cornerImaxJ.Y > baseZoneYMax - 0.01 && cornerImaxJ.Y < baseZoneYMax + 0.01)
                    {
                    isVolJMaxSurf = TRUE;
                    cornerBoxEdgePt1.Y = baseZoneYMax;
                    }

                if (cornerImaxJ.Z > 0.99 && cornerImaxJ.Z < 1.01)
                    {
                    isVolKMinSurf = TRUE;
                    cornerBoxEdgePt1.Z = baseZoneZMin;
                    }
                else if (cornerImaxJ.Z > baseZoneZMax - 0.01 && cornerImaxJ.Z < baseZoneZMax + 0.01)
                    {
                    isVolKMaxSurf = TRUE;
                    cornerBoxEdgePt1.Z = baseZoneZMax;
                    }

                Boolean_t isCorner = FALSE;

                if (isVolIMinSurf || isVolIMaxSurf || isVolJMinSurf || isVolJMaxSurf || isVolKMinSurf || isVolKMaxSurf)
                    {
                    XYZ_s cornerImaxJp1;
                    index = (jj+1) * iMax; // iMax,jj+1
                    cornerImaxJp1.X = TecUtilDataValueGetByRef(xFDPtr, index);
                    cornerImaxJp1.Y = TecUtilDataValueGetByRef(yFDPtr, index);
                    cornerImaxJp1.Z = TecUtilDataValueGetByRef(zFDPtr, index);

                        // Vol zone IMin, ? corner
                    if (cornerImaxJp1.X > 0.99 && cornerImaxJp1.X < 1.01 && !isVolIMinSurf &&
                        (isVolJMinSurf || isVolJMaxSurf || isVolKMinSurf || isVolKMaxSurf))
                        {
                            isCorner = TRUE;
                        cornerBoxEdgePt1.X = baseZoneXMin;
                        }

                        // Vol zone IMax, ? corner
                    if (cornerImaxJp1.X > baseZoneXMax - 0.01 && cornerImaxJp1.X < baseZoneXMax + 0.01 &&  !isVolIMaxSurf &&
                        (isVolJMinSurf || isVolJMaxSurf || isVolKMinSurf || isVolKMaxSurf))
                        {
                            isCorner = TRUE;
                        cornerBoxEdgePt1.X = baseZoneXMax;
                        }

                        // Vol zone JMin, ? corner
                    if (cornerImaxJp1.Y > 0.99 && cornerImaxJp1.Y < 1.01 && !isVolJMinSurf &&
                        (isVolIMinSurf || isVolIMaxSurf || isVolKMinSurf || isVolKMaxSurf))
                        {
                            isCorner = TRUE;
                        cornerBoxEdgePt1.Y = baseZoneYMin;
                        }

                        // Vol zone JMax, ? corner
                    if (cornerImaxJp1.Y > baseZoneYMax - 0.01 && cornerImaxJp1.Y < baseZoneYMax + 0.01 &&  !isVolJMaxSurf &&
                        (isVolIMinSurf || isVolIMaxSurf || isVolKMinSurf || isVolKMaxSurf))
                        {
                            isCorner = TRUE;
                        cornerBoxEdgePt1.Y = baseZoneYMax;
                        }

                        // Vol zone KMin, ? corner
                    if (cornerImaxJp1.Z > 0.99 && cornerImaxJp1.Z < 1.01 && !isVolKMinSurf &&
                        (isVolIMinSurf || isVolIMaxSurf || isVolJMinSurf || isVolJMaxSurf))
                        {
                            isCorner = TRUE;
                        cornerBoxEdgePt1.Z = baseZoneZMin;
                        }

                        // Vol zone KMax, ? corner
                    if (cornerImaxJp1.Z > baseZoneZMax - 0.01 && cornerImaxJp1.Z < baseZoneZMax + 0.01 && !isVolKMaxSurf &&
                        (isVolIMinSurf || isVolIMaxSurf || isVolJMinSurf || isVolJMaxSurf))
                        {
                            isCorner = TRUE;
                        cornerBoxEdgePt1.Z = baseZoneZMax;
                        }

                        if (isCorner)
                        {
                        // okay, now we have an intersection in the corner,
                        // so we can now find a better intersection by expanding the
                        // quad near that corner. First get the two other
                        // ("interior") points

                        XYZ_s cornerImaxm1Jp1;
                        index = (jj+1)*iMax - 1; // iMax-1,jj+1
                        cornerImaxm1Jp1.X = TecUtilDataValueGetByRef(xFDPtr, index);
                        cornerImaxm1Jp1.Y = TecUtilDataValueGetByRef(yFDPtr, index);
                        cornerImaxm1Jp1.Z = TecUtilDataValueGetByRef(zFDPtr, index);

                        XYZ_s cornerImaxm1J;
                        index = jj*iMax - 1; // iMax-1,jj
                        cornerImaxm1J.X = TecUtilDataValueGetByRef(xFDPtr, index);
                        cornerImaxm1J.Y = TecUtilDataValueGetByRef(yFDPtr, index);
                        cornerImaxm1J.Z = TecUtilDataValueGetByRef(zFDPtr, index);

                        // now make new quad beyond the zone's outer boundary by extending the last
                        // (imax) cell out by some amount (factor) times the last cell
                        double maxFactor = 1000.0;
                        double factor;
                        for ( factor = 1.0; factor <= maxFactor; factor *= 2.0 )
                        {
                            XYZ_s cornerImaxp1J;
                            cornerImaxp1J.X = cornerImaxJ.X + factor*(cornerImaxJ.X-cornerImaxm1J.X);
                            cornerImaxp1J.Y = cornerImaxJ.Y + factor*(cornerImaxJ.Y-cornerImaxm1J.Y);
                            cornerImaxp1J.Z = cornerImaxJ.Z + factor*(cornerImaxJ.Z-cornerImaxm1J.Z);

                            XYZ_s cornerImaxp1Jp1;
                            cornerImaxp1Jp1.X = cornerImaxJp1.X + factor*(cornerImaxJp1.X-cornerImaxm1Jp1.X);
                            cornerImaxp1Jp1.Y = cornerImaxJp1.Y + factor*(cornerImaxJp1.Y-cornerImaxm1Jp1.Y);
                            cornerImaxp1Jp1.Z = cornerImaxJp1.Z + factor*(cornerImaxJp1.Z-cornerImaxm1Jp1.Z);

                            // first we have to make sure we have made the quad big enough
                            XYZ_s cornerBoxEdgePt2;
                            cornerBoxEdgePt2 = cornerBoxEdgePt1;
                            if (cornerBoxEdgePt1.X == 0.0)
                                cornerBoxEdgePt2.X = baseZoneXMax + 1.0;
                            else if (cornerBoxEdgePt1.Y == 0.0)
                                cornerBoxEdgePt2.Y = baseZoneYMax + 1.0;
                            else if (cornerBoxEdgePt1.Z == 0.0)
                                cornerBoxEdgePt2.Z = baseZoneZMax + 1.0;
                            else
                                CHECK(FALSE);

                            int numIntersections;
                            XYZ_s intersectionPts[2];
                            if ( GeomToolsLineSegQuadIntersections(cornerImaxJ, cornerImaxJp1, cornerImaxp1Jp1, cornerImaxp1J,
                                                                   cornerBoxEdgePt1, cornerBoxEdgePt2,
                                                                   &numIntersections, intersectionPts) )
                            {
                                // corner edge works, so try original segment and then break out of loop
                                if ( GeomToolsLineSegQuadIntersections(cornerImaxJ, cornerImaxJp1, cornerImaxp1Jp1, cornerImaxp1J,
                                                                       lineSegStart, lineSegEnd,
                                                                       &numIntersections, intersectionPts) )
                                {
                                    // make sure intersections are truly inside the zone
                                    for ( int pp = 0; pp < numIntersections; pp++ )
                                    {
                                        if ( baseZoneXMin-ijkEpsilon <= intersectionPts[pp].X && intersectionPts[pp].X <= baseZoneXMax+ijkEpsilon &&
                                             baseZoneYMin-ijkEpsilon <= intersectionPts[pp].Y && intersectionPts[pp].Y <= baseZoneYMax+ijkEpsilon &&
                                             baseZoneZMin-ijkEpsilon <= intersectionPts[pp].Z && intersectionPts[pp].Z <= baseZoneZMax+ijkEpsilon )
                                        {
                                            intersectionCallback(intersectionPts[pp], clientData);
                        }
                    }
                }
                                break;
            }
        }
                        if ( factor > maxFactor )
                        {
                            // just create a triangle out to the edge
                            if (cornerBoxEdgePt1.X == 0.0)
                                cornerBoxEdgePt1.X = 0.5*(cornerImaxJ.X + cornerImaxJp1.X);
                            else if (cornerBoxEdgePt1.Y == 0.0)
                                cornerBoxEdgePt1.Y = 0.5*(cornerImaxJ.Y + cornerImaxJp1.Y);
                            else if (cornerBoxEdgePt1.Z == 0.0)
                                cornerBoxEdgePt1.Z = 0.5*(cornerImaxJ.Z + cornerImaxJp1.Z);
                            else
                                CHECK(FALSE);
                            XYZ_s intersectionPt;
                            if ( GeomToolsLineSegTriangleIntersection(cornerImaxJ, cornerImaxJp1, cornerBoxEdgePt1,
                                                                      lineSegStart, lineSegEnd,
                                                                      NULL,
                                                                      &intersectionPt) )
                            {
                                intersectionCallback(intersectionPt, clientData);
    }
                        }
                    }
                }
            }
        }
    }

    // Finish lock for this function (may be nested)
    TecUtilLockFinish(AddOnID);

    // Inform Tecplot that major data operation is ending
    TecUtilDataLoadEnd();
    }








/**
 * Determine if the BondBundles handle is sane.
 *
 * param BondBundles
 *     BondBundles structure in question.
 *
 * return
 *     TRUE if the BondBundles structure is valid, otherwise FALSE.
 */
Boolean_t BondBundlesIsValid(const BondBundles_pa BondBundles)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(BondBundles) &&
               VALID_REF(BondBundles->BondBundlesList) && ArrListIsValid(BondBundles->BondBundlesList));


    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}












/**
 * Deallocates the BondBundles handle and set the handle to NULL.
 *
 * param
 *     Reference to a BondBundles handle.
 */
void BondBundlesDealloc(BondBundles_pa *BondBundles)
{

    REQUIRE(VALID_REF(BondBundles));

    if (*BondBundles != NULL)
    {
        /* release the ArrList's */
        if ((*BondBundles)->BondBundlesList != NULL)
        {
            LgIndex_t BondOffset;
            LgIndex_t Count = BondBundlesGetCount(*BondBundles);
            BondBundle_pa bondBundle = NULL;

            // Dealloc the bondBundle structures in the ArrList
            for (BondOffset = 0; BondOffset < Count; BondOffset++)
            {
                bondBundle = BondBundlesGetBondBundle((*BondBundles), BondOffset);
                if (bondBundle != NULL)
                    BondBundleDealloc(&bondBundle);
            }

            ArrListDealloc(&((*BondBundles)->BondBundlesList));
        }

        /* release the list structure itself */
        FREE_ITEM(*BondBundles, "BondBundles structure");
        *BondBundles = NULL;
    }

    ENSURE(*BondBundles == NULL);
}




/**
 * Gets the number of bondBundle structures currently in the ArrList of BondBundles
 *
 * param
 *     BondBundles structure in question.
 *
 * return
 *     Number of BondBundles in the ArrList.
 */
LgIndex_t BondBundlesGetCount(const BondBundles_pa BondBundles)
{
    LgIndex_t Result = 0;
    Boolean_t isOk = TRUE;

    REQUIRE(BondBundlesIsValid(BondBundles));

    // Get the BondBundles count
    Result = ArrListGetCount(BondBundles->BondBundlesList);

    ENSURE(Result >= 0);
    return Result;
}







/**
 * Empties the array list in the BondBundles structure of all BondBundles
 *
 *
 * param BondBundles
 *     BondBundles structure to clear.
 */
void BondBundlesClear(BondBundles_pa BondBundles)
{
    REQUIRE(BondBundlesIsValid(BondBundles));

    ArrListClear(BondBundles->BondBundlesList);


    ENSURE(BondBundlesIsValid(BondBundles) &&
           BondBundlesGetCount(BondBundles) == 0);
}







/**
 * Allocates a BondBundles handle with a suitable default
 * capacity for the contained ArrList's. The initial capacity
 * may, optionally, be specified.
 *
 *
 * return
 *     BondBundles handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
BondBundles_pa BondBundlesAlloc(const LgIndex_t ApproxNumBonds)
{
    BondBundles_pa Result = NULL;
    LgIndex_t InitialSize = 20;

    if (ApproxNumBonds > 0) InitialSize = ApproxNumBonds;

    Result = ALLOC_ITEM(BondBundles_s, "BondBundles structure");
    if (Result != NULL)
    {
        Result->BondBundlesList = ArrListAlloc(InitialSize, ArrListType_VoidPtr);
    }

    ENSURE(BondBundlesIsValid(Result) || Result == NULL);
    return Result;
}







/**
 * Gets the bondBundle pointer at the specified offset.
 *
 * param BondBundles
 *     BondBundles to query for bondBundle pointer.
 * param Offset
 *     Offset into array list.
 *
 * return
 *     Pointer to bondBundle structure if it works, otherwise NULL.
 */
BondBundle_pa   BondBundlesGetBondBundle(const BondBundles_pa  BondBundles,
                                         const LgIndex_t       Offset)
{
    Boolean_t  isOk = TRUE;
    BondBundle_pa Result = NULL;
    ArrListItem_u Item;


    REQUIRE(BondBundlesIsValid(BondBundles));
    REQUIRE(0 <= Offset && Offset <= BondBundlesGetCount(BondBundles) - 1);

    // Get bondBundle pointer stored as specific offest in ArrList
    if (isOk)
    {
        Item   = ArrListGetItem(BondBundles->BondBundlesList, Offset);
        Result = (BondBundle_pa)Item.VoidPtr;
    }

    ENSURE(BondBundleIsValid(Result));
    return Result;
}






/**
 * Sets ArrList bondBundle pointer at the specified offset. If the
 * offset is beyond the end of the array list, it is sized
 * accordingly and the intervening array values between the last
 * node of the original state and the last node of the new state are
 * assigned NULL.
 *
 * note
 *     If bondBundle already exists at the specified location it is
 *     deallocated.
 *
 * param BondBundles
 *     bondBundle target in which to set the Angle/GradPath.
 * param Offset
 *     Offset into array list.
 * param bondBundle
 *     bondBundle pointer to set at the specified offset in ArrList.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t   BondBundlesSetBondBundle(BondBundles_pa       BondBundles,
                                     const LgIndex_t      Offset,
                                     const BondBundle_pa  bondBundle)
{
    Boolean_t isOk = TRUE;
    ArrListItem_u Item;



    REQUIRE(BondBundlesIsValid(BondBundles));
    REQUIRE(Offset >= 0);
    REQUIRE(bondBundle == NULL || BondBundleIsValid(bondBundle));

    // If current bondBundle at the specified offset is valid, dealloc it
    if (BondBundlesGetCount(BondBundles) > Offset)
    {
        BondBundle_pa OldBondBundle = BondBundlesGetBondBundle(BondBundles, Offset);

        if (OldBondBundle != NULL && BondBundleIsValid(OldBondBundle))
        {
            BondBundleDealloc(&OldBondBundle);
        }
    }


    // Set bondBundle at specified offset
    if (isOk)
    {
        Item.VoidPtr = (void *)bondBundle;
        isOk = ArrListSetItem(BondBundles->BondBundlesList, Offset, Item);
    }

    ENSURE(BondBundlesIsValid(BondBundles));

    ENSURE(VALID_BOOLEAN(isOk));
    return isOk;
}







/*
 * BondBundlesGetFromZoneAuxData: For a given BaseZone, cycle through the zones
 *  to find the bondBundle surfaces for all bonds. Store these in a set of
 *  bondBundle structures referenced by the pointers in BondBundlesList.
 *
 * param baseZoneNum
 *     Number of the zone for which the topology is computed.
 * param VolCritPoints
 *     Critical point structure for the volume critical points.
 * param BondBundlesList
 *     ArrayList of pointers to bondBundle data structures for each bond.
 *
 * return
 *     TRUE if successful, FALSE if there were errors.
 */
Boolean_t      BondBundlesGetFromZoneAuxData(const EntIndex_t    baseZoneNum,
                                             const CritPoints_pa VolCritPoints,
                                             BondBundles_pa      BondBundles)
//                                             ArrList_pa          BondBundlesList)
{
    Boolean_t isOk = TRUE;

    EntIndex_t NumZones, ZoneNum;
    LgIndex_t  NumBonds, BondOffset;

    REQUIRE(baseZoneNum > 0);
    REQUIRE(CritPointsIsValid(VolCritPoints));
    REQUIRE(BondBundlesIsValid(BondBundles));
    // REQUIRE(ArrListIsValid(BondBundlesList));

    /* Get num zones in dataset */
    if (isOk)
        TecUtilDataSetGetInfo(NULL, &NumZones, NULL);

    REQUIRE(baseZoneNum <= NumZones);

    // Get the number of bonds
    NumBonds =  CritPointsGetEndOffset(VolCritPoints, (char)(-1)) -
                CritPointsGetBegOffset(VolCritPoints, (char)(-1));

    // Create the BondBundles and set pointers in the BondBundlesList
    for (BondOffset = 0; isOk && BondOffset < NumBonds; BondOffset++)
    {
        BondBundle_pa bondBundle = BondBundleAlloc();

        // ArrListItem_u Item;

        // Item.VoidPtr = bondBundle;

        bondBundle->BaseVolZoneNum = baseZoneNum;
        bondBundle->BondNum = BondOffset + 1;

        isOk = BondBundlesSetBondBundle(BondBundles, BondOffset, bondBundle);
        // isOk = ArrListSetItem(BondBundlesList, BondOffset, Item);
    }
    CHECK(BondBundlesIsValid(BondBundles));
    CHECK(BondBundlesGetCount(BondBundles) == NumBonds);
    // CHECK(ArrListIsValid(BondBundlesList));
    // CHECK(ArrListGetCount(BondBundlesList) == NumBonds);

    // Cycle over the zones of the dataset, looking for zero-flux atom-ring-cage surfaces
    // and set the zone numbers in the appropriate BondBundles.
    for (ZoneNum = 1; isOk && ZoneNum <= NumZones; ZoneNum++)
    {
        Boolean_t  IsSurface = FALSE;
        Boolean_t  IsLine = FALSE;
        Boolean_t  BegCrtPtIsAtom = FALSE;
        Boolean_t  BegCrtPtIsBond = FALSE;
        AuxData_pa AuxDataRef = NULL;
        LgIndex_t  BondCrtPtNum1, BondCrtPtNum2, AtomCrtPtNum;
        char      *Value;
        Boolean_t  Retain;

        // Is ZoneNum a surface zone?
        if (TecUtilZoneIsOrdered(ZoneNum))
        {
            LgIndex_t IMax, JMax, KMax;
            // Use NULL for values we're not interested in
            TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, NULL, NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL);
            if (IMax > 1 && JMax > 1 && KMax == 1) IsSurface = TRUE;
            if (IMax > 1 && JMax == 1 && KMax == 1) IsLine = TRUE;
        }

        // Read Aux Data - check for correct baseZoneNum
        if (isOk && (IsSurface || IsLine))
        {
            AuxDataRef = TecUtilAuxDataZoneGetRef(ZoneNum);
            if (AuxDataRef != NULL)
            {
                /* Check for correct baseZoneNum */
                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.baseZoneNum",
                                                   &Value, &Retain))
                {
                    if (Value != NULL)
                    {
                        /* Assign the parsed value to the structure */
                        long ZoneBaseZoneNum;
                        if (sscanf(Value, "%d", &ZoneBaseZoneNum) == 1)
                        {
                            if (ZoneBaseZoneNum != baseZoneNum)
                            {
                                IsSurface = FALSE;
                                IsLine    = FALSE;
                            }
                        }
                        else
                        {
                            isOk = FALSE;
                        }
                        // release the allocated string copy
                        TecUtilStringDealloc(&Value);
                    }
                    else
                    {
                        isOk = FALSE;
                    }
                }
            }
        }

        // Check BegCrtPtType against "Atom" or "Bond"
        if (isOk && (IsSurface || IsLine))
        {
            if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtType",
                                               &Value, &Retain))
            {
                if (Value != NULL)
                {
                    if (strcmp(Value, "Atom") == 0) BegCrtPtIsAtom = TRUE;
                    if (strcmp(Value, "Bond") == 0) BegCrtPtIsBond = TRUE;
                    TecUtilStringDealloc(&Value);
                }
                else
                    isOk = FALSE;
            }
        }

        if (isOk && IsSurface && BegCrtPtIsAtom)
        {
            // Get the atom number
            if (isOk && IsSurface)
            {
                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
                                                   &Value, &Retain))
                {
                    if (Value != NULL)
                    {
                        /* Assign the parsed value to the structure */
                        long AtomCrtPtNumber;
                        if (sscanf(Value, "%d", &AtomCrtPtNumber) == 1)
                        {
                            AtomCrtPtNum = AtomCrtPtNumber;
                        }
                        else
                        {
                            isOk = FALSE;
                        }
                        // release the allocated string copy
                        TecUtilStringDealloc(&Value);
                    }
                    else
                    {
                        isOk = FALSE;
                    }
                }
            }


            /* Check OppositeCrtPtType against "BOND" */
            if (isOk && IsSurface)
            {
                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.OppositeCrtPtType",
                                                   &Value, &Retain))
                {
                    if (Value != NULL)
                    {
                        if (strcmp(Value, "Bond") != 0) IsSurface = FALSE;
                    }
                    else
                        isOk = FALSE;
                }
            }

            /* Get OppositeCrtPtNum's */
            if (isOk && IsSurface)
            {
                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.OppositeCrtPtNum1",
                                                   &Value, &Retain))
                {
                    if (Value != NULL)
                    {
                        /* Assign the parsed value to the structure */
                        long BondCrtPtNum;
                        if (sscanf(Value, "%d", &BondCrtPtNum) == 1)
                        {
                            if (BondCrtPtNum >= 1 || BondCrtPtNum <= NumBonds)
                                BondCrtPtNum1 = BondCrtPtNum;
                            else
                                IsSurface = FALSE;
                        }
                        else
                        {
                            isOk = FALSE;
                        }
                        // release the allocated string copy
                        TecUtilStringDealloc(&Value);
                    }
                    else
                    {
                        isOk = FALSE;
                    }
                }
            }

            if (isOk && IsSurface)
            {
                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.OppositeCrtPtNum2",
                                                   &Value, &Retain))
                {
                    if (Value != NULL)
                    {
                        /* Assign the parsed value to the structure */
                        long BondCrtPtNum;
                        if (sscanf(Value, "%d", &BondCrtPtNum) == 1)
                        {
                            if (BondCrtPtNum >= 1 || BondCrtPtNum <= NumBonds)
                                BondCrtPtNum2 = BondCrtPtNum;
                            else
                                IsSurface = FALSE;
                        }
                        else
                        {
                            isOk = FALSE;
                        }
                        // release the allocated string copy
                        TecUtilStringDealloc(&Value);
                    }
                    else
                    {
                        isOk = FALSE;
                    }
                }
            }

            // Get the pointer for the BondCrtPtNum1 bundle and set the surface
            if (isOk && IsSurface)
            {
                BondBundle_pa BondBundle1 = NULL;
                BondBundle_pa BondBundle2 = NULL;

                // ArrListItem_u Item;

                // Set the appropriate AtomRingCage surface (Pls or Mns) for Bond 1
				if (BondCrtPtNum1 >= 1)
				{
					BondBundle1 = BondBundlesGetBondBundle(BondBundles, BondCrtPtNum1 - 1);
					// Item = ArrListGetItem(BondBundlesList, BondCrtPtNum1 - 1);
					// BondBundle1 = (BondBundle_pa)(Item.VoidPtr);

					if (BondBundle1->PlsAtomNum == AtomCrtPtNum || BondBundle1->PlsAtomNum == -1)
					{
						BondBundle1->PlsAtomNum = AtomCrtPtNum;
						isOk = BondBundleAppendSurfZnNum(BondBundle1, SurfType_PlsAtomRingCage, ZoneNum);
					}
					else if (BondBundle1->MnsAtomNum == AtomCrtPtNum || BondBundle1->MnsAtomNum == -1)
					{
						BondBundle1->MnsAtomNum = AtomCrtPtNum;
						isOk = BondBundleAppendSurfZnNum(BondBundle1, SurfType_MnsAtomRingCage, ZoneNum);
					}
					else
						isOk = FALSE;
				}

                // Set the appropriate AtomRingCage surface (Pls or Mns) for Bond 2
				if (BondCrtPtNum2 >= 1)
				{
					BondBundle2 = BondBundlesGetBondBundle(BondBundles, BondCrtPtNum2 - 1);
					// Item = ArrListGetItem(BondBundlesList, BondCrtPtNum2 - 1);
					// BondBundle2 = (BondBundle_pa)(Item.VoidPtr);

					if (BondBundle2->PlsAtomNum == AtomCrtPtNum || BondBundle2->PlsAtomNum == -1)
					{
						BondBundle2->PlsAtomNum = AtomCrtPtNum;
						isOk = BondBundleAppendSurfZnNum(BondBundle2, SurfType_PlsAtomRingCage, ZoneNum);
					}
					else if (BondBundle2->MnsAtomNum == AtomCrtPtNum || BondBundle2->MnsAtomNum == -1)
					{
						BondBundle2->MnsAtomNum = AtomCrtPtNum;
						isOk = BondBundleAppendSurfZnNum(BondBundle2, SurfType_MnsAtomRingCage, ZoneNum);
					}
					else
						isOk = FALSE;
				}
            }
        }
        else if (isOk && (IsSurface || IsLine) && BegCrtPtIsBond)
        {
            LgIndex_t BondNum = -1;

            // Get the Bond number
            if (isOk && (IsSurface || IsLine))
            {
                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
                                                   &Value, &Retain))
                {
                    if (Value != NULL)
                    {
                        /* Assign the parsed value to the structure */
                        long BondCrtPtNumber;
                        if (sscanf(Value, "%d", &BondCrtPtNumber) == 1)
                        {
                            BondNum = BondCrtPtNumber - CritPointsGetBegOffset(VolCritPoints, (char)(-1));
                        }
                        else
                        {
                            isOk = FALSE;
                        }
                        // release the allocated string copy
                        TecUtilStringDealloc(&Value);
                    }
                    else
                    {
                        isOk = FALSE;
                    }
                }
            }

            // Get the pointer for the BondCrtPtNum and set the surface
            if (isOk && IsSurface && BondNum > 0)
            {
                BondBundle_pa bondBundle = NULL;

                // ArrListItem_u Item;

                // Set the appropriate AtomRingCage surface (Pls or Mns) for Bond 1
                bondBundle = BondBundlesGetBondBundle(BondBundles, BondNum - 1);
                // Item = ArrListGetItem(BondBundlesList, BondNum - 1);
                // bondBundle = (BondBundle_pa)(Item.VoidPtr);

                if (bondBundle->BondNum == BondNum || bondBundle->BondNum == -1)
                {
                    bondBundle->BondNum = BondNum;
                    isOk = BondBundleAppendSurfZnNum(bondBundle, SurfType_BondRingCage, ZoneNum);
                }
                else
                    isOk = FALSE;
            }

            // If is is a line, check that it is a Bond-Atom line and set AtomNum and ZoneNum
            if (isOk && IsLine && BondNum > 0)
            {
                BondBundle_pa bondBundle = NULL;

                LgIndex_t AtomNum = 0;

                // Check EndCrtPtType against "Atom"
                if (isOk && IsLine)
                {
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtType",
                                                       &Value, &Retain))
                    {
                        if (Value != NULL)
                        {
                            if (strcmp(Value, "Atom") != 0) IsLine = FALSE;
                            TecUtilStringDealloc(&Value);
                        }
                        else
                            isOk = FALSE;
                    }
                }

                // Find EndCrtPtNum (Atom number)
                if (isOk && IsLine)
                {
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
                                                       &Value, &Retain))
                    {
                        if (Value != NULL)
                        {
                            /* Assign the parsed value to the structure */
                            long AtomCrtPtNumber;
                            if (sscanf(Value, "%d", &AtomCrtPtNumber) == 1)
                            {
                                AtomNum = AtomCrtPtNumber;
                            }
                            else
                            {
                                isOk = FALSE;
                            }
                            // release the allocated string copy
                            TecUtilStringDealloc(&Value);
                        }
                        else
                        {
                            isOk = FALSE;
                        }
                    }
                }

                bondBundle = BondBundlesGetBondBundle(BondBundles, BondNum - 1);

                // Set the Pls or Mns AtomNum and ZoneNum (if not already set)
                if (isOk && IsLine)
                {
                    if (bondBundle->PlsAtomNum <= 0 || bondBundle->PlsAtomNum == AtomNum)
                    {
                        bondBundle->PlsAtomNum = AtomNum;
                        bondBundle->BondPlsAtomLineZnNum = ZoneNum;
                    }
                    else if (bondBundle->MnsAtomNum <= 0 || bondBundle->MnsAtomNum == AtomNum)
                    {
                        bondBundle->MnsAtomNum = AtomNum;
                        bondBundle->BondMnsAtomLineZnNum = ZoneNum;
                    }
                }
            }

        }

    }





    ENSURE(VALID_BOOLEAN(isOk));
    return isOk;
}










/*
 * BondBundleAdd: For a given max (atom) volume critical point, compute
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
 * param bondBundle
 *     bondBundle data structure to populate for an Atom volume critical point.
 *
 * return
 *     TRUE if successful, FALSE if there were errors.
 */
Boolean_t BondBundleAdd(const ZoneVarInfo_pa VolZoneVarInfo,
                        const LgIndex_t      BondCrtPtNum,
                        const Bundles_pa     VolBundles,
                        BondBundle_pa        bondBundle)
{
    Boolean_t      isOk = TRUE;
    // double         XCrtPtAtom, YCrtPtAtom, ZCrtPtAtom, dummy;
    // char           cdummy;
    LgIndex_t      FirstElemNumOffset = -1;
    EntIndex_t     ChrgDensVarNum = VolZoneVarInfo->ChrgDensVarNum;
    EntIndex_t     UVarNum = VolZoneVarInfo->UVarNum;
    EntIndex_t     VVarNum = VolZoneVarInfo->VVarNum;
    EntIndex_t     WVarNum = VolZoneVarInfo->WVarNum;
    /*

    REQUIRE(VALID_REF(VolZoneVarInfo));
    REQUIRE(CritPointsIsValid(VolCritPoints));
    REQUIRE(AtomCrtPtNum >= 0 && AtomCrtPtNum < CritPointsGetCount(VolCritPoints) );
    // REQUIRE(BundlesIsValid(VolBundles));
    REQUIRE(CPTolerance > 0.0);
    REQUIRE(BondBundleIsValid(bondBundle));
    // Atom critical points must be an Atom
    REQUIRE(CritPointsGetType(VolCritPoints, AtomCrtPtNum) == -3);


    // For a given max (atom) volume critical point, compute the "nearly largest" isosurface
    // of the scalar variable surrounding the  atom.
    IsoTopoZone = IsoTopoZoneForAtom(VolZoneVarInfo, AtomCrtPtNum, VolCritPoints, CPTolerance);

    bondBundle->IsoTopoZone = IsoTopoZone;


    // if (0)
      {
        // For convenience, get pointers to SurfCritPoints and SurfBundles
        SurfCritPoints = bondBundle->SurfCritPoints;
        SurfBundles   = bondBundle->SurfBundles;

        if (isOk && IsoTopoZone > 0)
          {
            // TEMPORARY >>>>
            EntIndex_t GradTotVNum = 8;
            // TEMPORARY <<<<
            LgIndex_t baseZoneIMax, baseZoneJMax, baseZoneKMax;




            // Compute the surface element map for IsoTopoZone
            TecUtilZoneGetInfo(isoTopoZone, &baseZoneIMax, &baseZoneJMax, &baseZoneKMax,  NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            if (isOk)
              {
                // Compute the element map
                SurfElemMap = SurfElemMapAlloc(baseZoneIMax);
                isOk = SurfElemMapCompute(SurfElemMap, IsoTopoZone);
              }

            // Compute the surface normals
            if (isOk)
              {
                Normals = NormalsAlloc(baseZoneIMax);
                // isOk = NormalsComputeUsingSEM(Normals, IsoTopoZone, SurfElemMap);
                isOk = NormalsComputeUsingIsoVarGrad(Normals, IsoTopoZone, VolZoneVarInfo->UVarNum, VolZoneVarInfo->VVarNum, VolZoneVarInfo->WVarNum);
              }


            if (isOk)
              {
                ZoneVarInfo = ALLOC_ITEM(ZoneVarInfo_s, "ZoneVarInfo");
                ZoneVarInfo->ZoneNum = IsoTopoZone;
                ZoneVarInfo->UVarNum = SurfZoneVarInfo->UVarNum;
                ZoneVarInfo->VVarNum = SurfZoneVarInfo->VVarNum;
                ZoneVarInfo->WVarNum = SurfZoneVarInfo->WVarNum;
                ZoneVarInfo->ChrgDensVarNum = GradTotVNum;
                ZoneVarInfo->TypeVarNum = SurfZoneVarInfo->TypeVarNum;
                ZoneVarInfo->PeriodicBC = PeriodicBC;
                ZoneVarInfo->SurfElemMap = SurfElemMap;
                ZoneVarInfo->Normals = Normals;
                ZoneVarInfo->PathDir = StreamDir_Forward;
              }

            if (isOk)
              isOk = ExtractCriticalPoints(VolZoneVarInfo, ZoneVarInfo, TRUE, SurfCritPoints);
            // isOk = ExtractCriticalPoints(IsoTopoZone, GradXGTotVNum, GradYGTotVNum, GradZGTotVNum,
            //   GradTotVNum, TypeVarNum, PeriodicBC, SurfCritPoints);

            // TEMPORARY: Create Critical Point Zone
            if (isOk)
              {
                SurfCPZoneNum = CritPointsCreateTPZone(SurfCritPoints, ChrgDensVarNum, TypeVarNum,
                                                       UVarNum, VVarNum, WVarNum);
                // if (SurfCPZoneNum > 0) NumZones++;
                CHECK(SurfCPZoneNum > 0);
              }

            // Compute the surface element map for IsoTopoZone
            // TecUtilZoneGetInfo(isoTopoZone, &baseZoneIMax, &baseZoneJMax, &baseZoneKMax,  NULL, NULL, NULL,
            //                    NULL, NULL, NULL, NULL, NULL, NULL, NULL);

            // Interpolate variables from volume zone to new surface CritPoints zone
            if (isOk && SurfCPZoneNum > 0)
              {
                Boolean_t isOk;
                Set_pa Zones = TecUtilSetAlloc(TRUE);
                Set_pa Vars  = TecUtilSetAlloc(TRUE);

                TecUtilSetAddMember(Zones, VolZoneNum, TRUE);

                TecUtilSetAddMember(Vars, UVarNum, TRUE);
                TecUtilSetAddMember(Vars, VVarNum, TRUE);
                TecUtilSetAddMember(Vars, WVarNum, TRUE);

                isOk = TecUtilLinearInterpolate(Zones, SurfCPZoneNum, Vars, 0.0,
                                                LinearInterpMode_SetToConst);
                TecUtilSetDealloc(&Zones);
                TecUtilSetDealloc(&Vars);
              }

            //
            // To get surface Saddle-Max lines, seed streamtrace a small step in the
            // PrincDir away from Saddle
            //
            if (isOk)
            // if (0)
              {
                GradPath_pa    GradPath = NULL;
                LgIndex_t      BegOffset = CritPointsGetBegOffset(SurfCritPoints, -1);
                LgIndex_t      EndOffset = CritPointsGetEndOffset(SurfCritPoints, -1);
                LgIndex_t      NumCrtPts = CritPointsGetCount(SurfCritPoints);
                LgIndex_t      ii;



                // Set-up status bar
                // if (ShowStatusBar)
                //   TecUtilStatusStartPercentDone( "Finding Bond Lines", TRUE, TRUE);

                for (ii=BegOffset; isOk && ii<EndOffset; ii++)
                // if (0)
                  {
                    double XCrtPt, YCrtPt, ZCrtPt, PrincDirX, PrincDirY, PrincDirZ, dummy;
                    double Princ2DirX, Princ2DirY, Princ2DirZ;
                    char   TypeCrtPt;
                    double XPos, YPos, ZPos;

                    // Used in defining surface bundles
                    LgIndex_t MaxCPFirst = -1;
                    LgIndex_t MaxCPSecond = -1;

                    isOk = CritPointsGetPoint(SurfCritPoints, ii, &XCrtPt, &YCrtPt, &ZCrtPt,
                      &dummy, &TypeCrtPt, &PrincDirX, &PrincDirY, &PrincDirZ);

                    // Surface GradPath from saddle to max - first one (in PrincDir)

                    // XPos = XCrtPt + 0.1 * PrincDirX;
                    // YPos = YCrtPt + 0.1 * PrincDirY;
                    // ZPos = ZCrtPt + 0.1 * PrincDirZ;
                    XPos = XCrtPt + 0.5 * PrincDirX;
                    YPos = YCrtPt + 0.5 * PrincDirY;
                    ZPos = ZCrtPt + 0.5 * PrincDirZ;

                    GradPath = GradPathAlloc();

                    // Inform Tecplot that major data operation is beginning
                    TecUtilDataLoadBegin();
                    TecUtilLockStart(AddOnID);

                    // Integrate to compute path lines
                    if (isOk)
                      {
                        double    junk = 0.0;
                        LgIndex_t NumPathPoints = 0;

                        GradPathClear(GradPath);

                        // Seed first point
                        isOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk);
                        GradPath->BeginCrtPtNum = ii;

                        // Integrate gradient path line
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
                        // isOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Forward,
                        //   SurfCritPoints, 1.0, &NumPathPoints, GradPath);
                        isOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Forward,
                          SurfCritPoints, 0.5, &NumPathPoints, GradPath);

                        // Add a Tecplot zone for the Surface Saddle(Ring)-Max(Cage) connector
                        if (isOk && NumPathPoints > 0 && GradPath->EndCrtPtNum >= 0 && GradPath->EndCrtPtNum < NumCrtPts)
                          {
                            LgIndex_t  M1CrtPtNum = ii - BegOffset;
                            EntIndex_t ConnectorZoneNum = 0;
                            EntIndex_t ConnectorType = -4;  // Atom-Bond

                            MaxCPFirst = GradPath->EndCrtPtNum;

                            // Resample before writing the tecplot zone
                            isOk = GradPathResample(&GradPath, NUMRESAMPLEDPOINTSSURFACE);

                            if (isOk) isOk = BondBundleAppendGP(bondBundle, GradPath);

                          }
                      }

                    // Inform Tecplot that major data operation is ending
                    TecUtilLockFinish(AddOnID);
                    TecUtilDataLoadEnd();
                  }

                // Clean up allocated structures/arrays
                if ( ZoneVarInfo != NULL )
                  FREE_ITEM(ZoneVarInfo, "ZoneVarInfo");

                if (!isOk && GradPath != NULL) GradPathDealloc(&GradPath);
              }
            // if(ShowStatusBar) TecUtilStatusFinishPercentDone();
          }
      }
    */


    ENSURE(VALID_BOOLEAN(isOk));
    return isOk;
}
