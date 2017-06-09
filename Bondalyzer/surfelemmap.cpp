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
#include "ADDONVER.h"
#include "TASSERT.h"
#include "GUIDEFS.h"
#include "ENGINE.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "SURFELEMMAP.h"
#include <string.h>



/**
 * Determine if the SurfElemMap handle is sane.
 *
 * param SurfElemMap
 *     SurfElemMap structure in question.
 *
 * return
 *     TRUE if the SurfElemMap structure is valid, otherwise FALSE.
 */
Boolean_t SurfElemMapIsValid(SurfElemMap_pa SurfElemMap)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(SurfElemMap) &&
               VALID_REF(SurfElemMap->ElemMapListOffset) && ArrListIsValid(SurfElemMap->ElemMapListOffset) &&
               VALID_REF(SurfElemMap->ElemMapList) && ArrListIsValid(SurfElemMap->ElemMapList));

    // Unit test some functions
    // if (IsValid)
    // IsValid = SurfElemMapTestPointTriDist();

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the SurfElemMap handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a SurfElemMap handle.
 */
void SurfElemMapDealloc(SurfElemMap_pa *SurfElemMap)
{
    REQUIRE(VALID_REF(SurfElemMap));
    REQUIRE(SurfElemMapIsValid(*SurfElemMap) || *SurfElemMap == NULL);

    if (*SurfElemMap != NULL)
    {
        /* release the ArrList's */
        if ((*SurfElemMap)->ElemMapListOffset != NULL) ArrListDealloc(&((*SurfElemMap)->ElemMapListOffset));
        if ((*SurfElemMap)->ElemMapList != NULL) ArrListDealloc(&((*SurfElemMap)->ElemMapList));

        /* release the list structure itself */
        FREE_ITEM(*SurfElemMap, "SurfElemMap structure");
        *SurfElemMap = NULL;
    }

    ENSURE(*SurfElemMap == NULL);
}





/**
 * Empties the SurfElemMap structure.
 *
 * param SurfElemMap
 *     SurfElemMap to clear.
 */
void SurfElemMapClear(SurfElemMap_pa SurfElemMap)
{
    REQUIRE(SurfElemMapIsValid(SurfElemMap));

    ArrListClear(SurfElemMap->ElemMapListOffset);
    ArrListClear(SurfElemMap->ElemMapList);

    ENSURE(SurfElemMapIsValid(SurfElemMap) && ArrListGetCount(SurfElemMap->ElemMapListOffset) == 0 &&
           ArrListGetCount(SurfElemMap->ElemMapList) == 0);
}






/**
 * Allocates a SurfElemMap handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * param
 *     NumNodes for the surface zone.
 *
 * return
 *     SurfElemMap handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
SurfElemMap_pa SurfElemMapAlloc(const LgIndex_t NumNodes)
{
    SurfElemMap_pa Result = NULL;

    Result = ALLOC_ITEM(SurfElemMap_s, "SurfElemMap structure");
    if (Result != NULL)
    {
        Result->ElemMapListOffset  = ArrListAlloc(NumNodes + 1, ArrListType_Long);

        /* If it failed to allocate the array lists, clean-up and exit. */
        if (Result->ElemMapListOffset == NULL)
        {
            FREE_ITEM(Result, "SurfElemMap structure");
            Result = NULL;
        }
    }

    // Initialize offsets to zero (no element numbers)
    if (Result != NULL)
    {
        LgIndex_t nn;
        Boolean_t IsOk = TRUE;
        for (nn = 0; IsOk && nn <= NumNodes; nn++)
        {
            ArrListItem_u Item;
            Item.Long = 0;
            IsOk = ArrListSetItem(Result->ElemMapListOffset, nn, Item);
        }
        if (IsOk == FALSE)
        {
            ArrListDealloc(&(Result->ElemMapListOffset));

            FREE_ITEM(Result, "SurfElemMap structure");
            Result = NULL;
        }
    }

    if (Result != NULL)
    {
        Result->ElemMapList        = ArrListAlloc(120, ArrListType_Long);

        /* If it failed to allocate the array lists, clean-up and exit. */
        if (Result->ElemMapList == NULL)
        {
            ArrListDealloc(&(Result->ElemMapListOffset));

            FREE_ITEM(Result, "SurfElemMap structure");
            Result = NULL;
        }
    }


    ENSURE(SurfElemMapIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of nodes currently in the SurfElemMap structure.
 *
 * param
 *     SurfElemMap structure in question.
 *
 * return
 *     Number of nodes in the SurfElemMap structure.
 */
LgIndex_t SurfElemMapGetNodeCount(const SurfElemMap_pa SurfElemMap)
{
    LgIndex_t Result = 0;

    REQUIRE(SurfElemMapIsValid(SurfElemMap));

    // Subtract 1 because there is a final "offset" at the end to
    // conveniently handle the element count for the last node.
    Result = ArrListGetCount(SurfElemMap->ElemMapListOffset) - 1;

    ENSURE(Result >= 0);
    return Result;
}





/**
 * Gets the number of Elements currently specified for a Node
 * in the SurfElemMap structure.
 *
 * param
 *     SurfElemMap structure in question.
 * param
 *     NodeNum
 *
 * return
 *     Number of elements for a specified Node in the SurfElemMap
 *     structure.
 */
LgIndex_t SurfElemMapGetElemCountForNode(const SurfElemMap_pa SurfElemMap,
                                         const LgIndex_t      NodeNum)
{
    LgIndex_t Result = 0;
    ArrListItem_u Item;
    LgIndex_t BegOffset;

    REQUIRE(SurfElemMapIsValid(SurfElemMap));
    REQUIRE(NodeNum <= SurfElemMapGetNodeCount(SurfElemMap));

    Item = ArrListGetItem(SurfElemMap->ElemMapListOffset, NodeNum - 1);
    BegOffset = Item.Long;

    Item = ArrListGetItem(SurfElemMap->ElemMapListOffset, NodeNum);

    Result = Item.Long - BegOffset;

    ENSURE(Result >= 0);
    return Result;
}




/**
 * Places SurfElemMap point coordinates and variable value at
 * the specified offset. If the offset is beyond the end of the
 * array lists, they are sized accordingly and the intervening
 * array values between the last node of the original state and
 * the last node of the new state are assigned 0.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced.
 *
 * param SurfElemMap
 *     SurfElemMap target in which to set the coordinates and var value.
 * param PointOffset
 *     Offset of the point/node.
 * param Coord1, Coord2, Var
 *     Coordinates and var value to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
/*
Boolean_t SurfElemMapSetPointAtOffset(SurfElemMap_pa SurfElemMap,
                                  LgIndex_t  ElemOffset,
                                  LgIndex_t  IsoElemNum )
{
  Boolean_t IsOk = TRUE;
  ArrListItem_u Item;


  REQUIRE(SurfElemMapIsValid(SurfElemMap));
  REQUIRE(ElemOffset >= 0);

  Item.Long = IsoElemNum;
  IsOk = ArrListSetItem(SurfElemMap->IsoElemNums, ElemOffset, Item);

  ENSURE(SurfElemMapIsValid(SurfElemMap));

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}
*/




/**
 * Inserts SurfElemMap element number for the specified NodeNum at the specified
 * offset. The ElemMapList array will be expanded to accomodate the additional
 * value and the ElemMapListOffsets will be incremented for node numbers > NodeNum.
 *
 *
 * param SurfElemMap
 *     SurfElemMap target in which to set the coordinates.
 * param NodeNum
 *     Number of node for which the element number is added.
 * param ElemOffset
 *     Offset (for the node) at which to insert the ElemNum.
 * param ElemNum
 *     Number of an element in the Isosurface zone to set at the specified node
 *     offset. This must be an element that contains NodeNum.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t SurfElemMapInsertElemForNodeOffset(SurfElemMap_pa SurfElemMap,
                                             LgIndex_t      NodeNum,
                                             LgIndex_t      ElemOffset,
                                             LgIndex_t      ElemNum)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;
    LgIndex_t Offset;


    REQUIRE(SurfElemMapIsValid(SurfElemMap));
    REQUIRE(0 <= NodeNum  && NodeNum <= SurfElemMapGetNodeCount(SurfElemMap));
    REQUIRE(0 <= ElemOffset && ElemOffset <= SurfElemMapGetElemCountForNode(SurfElemMap, NodeNum));
    REQUIRE(0 < ElemNum);

    // Add Element number to ElemMapList
    Item = ArrListGetItem(SurfElemMap->ElemMapListOffset, NodeNum - 1);
    Offset = Item.Long + ElemOffset;

    Item.Long = ElemNum;
    IsOk = ArrListInsertItem(SurfElemMap->ElemMapList, Offset, Item);

    // Increment all ElemMapListOffset values for NodeNums > NodeNum
    // (TODO: This should be optimized using raw pointers.)
    if (IsOk)
    {
        LgIndex_t nn;
        LgIndex_t NumNodes = ArrListGetCount(SurfElemMap->ElemMapListOffset);
        ArrListItem_u Item;
        LgIndex_t NodeOffset;

        for (nn = NodeNum; IsOk && nn < NumNodes; nn++)
        {
            Item = ArrListGetItem(SurfElemMap->ElemMapListOffset, nn);
            NodeOffset = Item.Long;
            NodeOffset++;
            Item.Long = NodeOffset;
            IsOk = ArrListSetItem(SurfElemMap->ElemMapListOffset, nn, Item);
        }
    }

    ENSURE(SurfElemMapIsValid(SurfElemMap));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Appends ElemNum to desired NodeNum seciton of the SurfElemMap structure.
 * The array lists will be expanded to accommodate the additional items.
 *
 * param SurfElemMap
 *     SurfElemMap target to which the IsoElemNum is to be appended.
 * param NodeNum
 *     Number of node to which the element (ElemNum) contains.
 * param ElemNum
 *     Element number, in the isosurface zone, to be appended.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t SurfElemMapAppendElemForNode(SurfElemMap_pa  SurfElemMap,
                                       LgIndex_t       NodeNum,
                                       LgIndex_t       ElemNum)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;
    Boolean_t IsUnique = TRUE;

    REQUIRE(SurfElemMapIsValid(SurfElemMap));

    Count = SurfElemMapGetElemCountForNode(SurfElemMap, NodeNum);

    // Only add if it is a unique element number for this node
    if (Count > 0)
    {
        int Offset;
        for (Offset = 0; IsUnique && Offset < Count; Offset++)
        {
            LgIndex_t ExistingElem = SurfElemMapGetElemNumForNodeOffset(SurfElemMap, NodeNum, Offset);
            if (ExistingElem == ElemNum) IsUnique = FALSE;
        }
    }
    if (IsUnique)
        IsOk = SurfElemMapInsertElemForNodeOffset(SurfElemMap, NodeNum, Count, ElemNum);

    ENSURE(SurfElemMapIsValid(SurfElemMap));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Sets the ElemNum from the SurfElemMap structure for the
 * element associated with the specified NodeNum at the
 * specified offset.
 *
 * param SurfElemMap
 *     SurfElemMap structure containing the desired item.
 * param NodeNum
 *     Number of the node of interest.
 * param ElemOffset
 *     Offset to the ElemNum for the specified node.
 * param ElemNum
 *     Number of an element in the Isosurface zone to set at the specified node
 *     offset. This must be an element that contains NodeNum.
 *
 * return
 *     TRUE if successful, otherwise FALSE.
 */
LgIndex_t SurfElemMapSetElemNumForNodeOffset(const SurfElemMap_pa  SurfElemMap,
                                             const LgIndex_t       NodeNum,
                                             const LgIndex_t       ElemOffset,
                                             const LgIndex_t       ElemNum)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t NodeOffset = -1;
    LgIndex_t ElemCount = 0;
    ArrListItem_u Item;

    REQUIRE(SurfElemMapIsValid(SurfElemMap));
    REQUIRE(ElemOffset >= 0);

    ElemCount = SurfElemMapGetElemCountForNode(SurfElemMap, NodeNum);
    Item = ArrListGetItem(SurfElemMap->ElemMapListOffset, NodeNum - 1);
    NodeOffset = Item.Long;

    if (ElemOffset < ElemCount)
    {
        if (NodeOffset > -1)
        {
            Item.Long = ElemNum;
            IsOk = ArrListSetItem(SurfElemMap->ElemMapList, NodeOffset + ElemOffset, Item);
        }
    }
    else
    {
        IsOk = SurfElemMapInsertElemForNodeOffset(SurfElemMap, NodeNum, ElemOffset, ElemNum);
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Gets the ElemNum from the SurfElemMap structure for the
 * element associated with the specified NodeNum at the
 * specified offset.
 *
 * param SurfElemMap
 *     SurfElemMap structure containing the desired item.
 * param NodeNum
 *     Number of the node of interest.
 * param ElemOffset
 *     Offset to the ElemNum for the specified node.
 *
 * return
 *     ElemNum if it works, -1 otherwise.
 */

LgIndex_t SurfElemMapGetElemNumForNodeOffset(const SurfElemMap_pa  SurfElemMap,
                                             const LgIndex_t       NodeNum,
                                             const LgIndex_t       ElemOffset)
{
    LgIndex_t ElemNum = -1;
    LgIndex_t NodeOffset = -1;
    ArrListItem_u Item;

    REQUIRE(SurfElemMapIsValid(SurfElemMap));
    REQUIRE(ElemOffset >= 0 &&
            ElemOffset < SurfElemMapGetElemCountForNode(SurfElemMap, NodeNum));
    Item = ArrListGetItem(SurfElemMap->ElemMapListOffset, NodeNum - 1);
    NodeOffset = Item.Long;

    if (NodeOffset > -1)
    {
        Item = ArrListGetItem(SurfElemMap->ElemMapList, NodeOffset + ElemOffset);
        ElemNum = Item.Long;
    }

    ENSURE(ElemNum >= -1);
    return ElemNum;
}







/**
 * Create the element map for a surface (FEQuad or FETriangle) zone.
 *
 * param SurfElemMap
 *     SurfElemMap structure containing the input points for the curve fit.
 * param SurfZone
 *     Number of the surface zone.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t SurfElemMapCompute(SurfElemMap_pa SurfElemMap,
                             EntIndex_t     SurfZone)
{
    EntIndex_t NumZones, NumVars;
    ZoneType_e ZoneType;
    Boolean_t  IsOk = TRUE;
    LgIndex_t  NumNodes, NumElems;

    REQUIRE(SurfElemMapIsValid(SurfElemMap));

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    ZoneType = TecUtilZoneGetType(SurfZone);

    REQUIRE(SurfZone > 0 && SurfZone <= NumZones);
    REQUIRE(ZoneType == ZoneType_FETriangle || ZoneType == ZoneType_FEQuad);

    // Loop over elements, adding ElemNum to SurfElemMap for contained nodes
    if (IsOk)
    {
        LgIndex_t    ne;
        int          NodesPerElem = 3;
        Boolean_t    IsQuad = FALSE;

        // FETriangle or FEQuad?
        if (ZoneType == ZoneType_FEQuad)
        {
            IsQuad = TRUE;
            NodesPerElem = 4;
        }

        // Find the number of Elements in the IsoSurfZone
        TecUtilZoneGetInfo(SurfZone, &NumNodes, &NumElems,
                           NULL, NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL);

        // Loop over elements of SurfZone, adding element number for node numbers
        for (ne = 0; IsOk && ne < NumElems; ne++)
        {
            int       ii;
            LgIndex_t NodeNum;

            // Find unique node numbers
            for (ii = 0; IsOk && ii < NodesPerElem; ii++)
            {
                NodeNum = TecUtilDataNodeGetByZone(SurfZone, ne + 1, ii + 1);
                IsOk =  SurfElemMapAppendElemForNode(SurfElemMap, NodeNum, ne + 1);
            }
        }
    }

    REQUIRE(SurfElemMapIsValid(SurfElemMap));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}



/*
 * Unit test for SurfElemMap by creating an FEQuad zone (brick with
 * six Quads and 8 nodes) and see if it works correctly.
 *
 * return
 *     TRUE if successful, FALSE  if it fails.
 */
Boolean_t SurfElemMapTest()
{
    Boolean_t      IsOk = TRUE;
    EntIndex_t     NewZone;
    SurfElemMap_pa SurfElemMap = NULL;
    Set_pa         ZoneList = NULL;

    // Create a temporary quad zone
    IsOk = TecUtilDataSetAddZone("Temporary FEQuad Zone", 8, 6, 1,
                                 ZoneType_FEQuad, NULL);
    if (IsOk)
    {
        NodeMap_pa NodeMap = NULL;
        Set_pa zones_added = TecUtilSetAlloc(TRUE);
        // new zone is always last zone
        TecUtilDataSetGetInfo(NULL, &NewZone, NULL);

        // fill NodeMap
        //
        //         n8 -----------------  n7     f6 is in back
        //          /                / |
        //        /   |  f4        /   |
        //     n4 ---------------- n3  |
        //       | f5 |           |    |
        //       |       f1       | f3 |
        //       |  n5| _  _  _ _ | _ _| n6
        //       |   /            |   /
        //       | /      f2      | /
        //     n1 ---------------- n2
        //
        //
        NodeMap = TecUtilDataNodeGetWritableRef(NewZone);
        TecUtilDataNodeSetByRef(NodeMap, 1, 1, 1);
        TecUtilDataNodeSetByRef(NodeMap, 1, 2, 2);
        TecUtilDataNodeSetByRef(NodeMap, 1, 3, 3);
        TecUtilDataNodeSetByRef(NodeMap, 1, 4, 4);

        TecUtilDataNodeSetByRef(NodeMap, 2, 1, 1);
        TecUtilDataNodeSetByRef(NodeMap, 2, 2, 2);
        TecUtilDataNodeSetByRef(NodeMap, 2, 3, 6);
        TecUtilDataNodeSetByRef(NodeMap, 2, 4, 5);

        TecUtilDataNodeSetByRef(NodeMap, 3, 1, 2);
        TecUtilDataNodeSetByRef(NodeMap, 3, 2, 6);
        TecUtilDataNodeSetByRef(NodeMap, 3, 3, 7);
        TecUtilDataNodeSetByRef(NodeMap, 3, 4, 3);

        TecUtilDataNodeSetByRef(NodeMap, 4, 1, 4);
        TecUtilDataNodeSetByRef(NodeMap, 4, 2, 3);
        TecUtilDataNodeSetByRef(NodeMap, 4, 3, 7);
        TecUtilDataNodeSetByRef(NodeMap, 4, 4, 8);

        TecUtilDataNodeSetByRef(NodeMap, 5, 1, 1);
        TecUtilDataNodeSetByRef(NodeMap, 5, 2, 5);
        TecUtilDataNodeSetByRef(NodeMap, 5, 3, 8);
        TecUtilDataNodeSetByRef(NodeMap, 5, 4, 4);

        TecUtilDataNodeSetByRef(NodeMap, 6, 1, 5);
        TecUtilDataNodeSetByRef(NodeMap, 6, 2, 6);
        TecUtilDataNodeSetByRef(NodeMap, 6, 3, 7);
        TecUtilDataNodeSetByRef(NodeMap, 6, 4, 8);

        // inform Tecplot of new zone
        TecUtilSetAddMember(zones_added, NewZone, TRUE);
        TecUtilStateChanged(StateChange_ZonesAdded,
                            (ArbParam_t)zones_added);
        TecUtilSetDealloc(&zones_added);
    }


    // Compute the element map
    if (IsOk)
    {
        SurfElemMap = SurfElemMapAlloc(8);
        IsOk = SurfElemMapCompute(SurfElemMap, NewZone);
    }

    // Test the results: each node should have 3 elements in SurfElemMap
    if (IsOk) IsOk = (SurfElemMapGetElemCountForNode(SurfElemMap, 1) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemCountForNode(SurfElemMap, 2) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemCountForNode(SurfElemMap, 3) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemCountForNode(SurfElemMap, 4) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemCountForNode(SurfElemMap, 5) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemCountForNode(SurfElemMap, 6) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemCountForNode(SurfElemMap, 7) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemCountForNode(SurfElemMap, 8) == 3);

    // Test the results: node 1 should have elements 1, 2, 5
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 1, 0) == 1);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 1, 1) == 2);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 1, 2) == 5);

    // Test the results: node 2 should have elements 1, 2, 3
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 2, 0) == 1);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 2, 1) == 2);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 2, 2) == 3);

    // Test the results: node 3 should have elements 1, 3, 4
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 3, 0) == 1);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 3, 1) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 3, 2) == 4);

    // Test the results: node 4 should have elements 1, 4, 5
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 4, 0) == 1);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 4, 1) == 4);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 4, 2) == 5);

    // Test the results: node 5 should have elements 2, 5, 6
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 5, 0) == 2);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 5, 1) == 5);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 5, 2) == 6);

    // Test the results: node 6 should have elements 2, 3, 6
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 6, 0) == 2);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 6, 1) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 6, 2) == 6);

    // Test the results: node 7 should have elements 3, 4, 6
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 7, 0) == 3);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 7, 1) == 4);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 7, 2) == 6);

    // Test the results: node 8 should have elements 4, 5, 6
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 8, 0) == 4);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 8, 1) == 5);
    if (IsOk) IsOk = (SurfElemMapGetElemNumForNodeOffset(SurfElemMap, 8, 2) == 6);


    // Delete the temporary quad zone
    ZoneList = TecUtilSetAlloc(FALSE);
    if (ZoneList)
    {
        TecUtilSetAddMember(ZoneList, NewZone, FALSE);
        TecUtilDataSetDeleteZone(ZoneList);
        TecUtilSetDealloc(&ZoneList);
    }

    // Delete the SurfElemMap
    SurfElemMapDealloc(&SurfElemMap);

    REQUIRE(VALID_BOOLEAN(IsOk));
    return IsOk;
}