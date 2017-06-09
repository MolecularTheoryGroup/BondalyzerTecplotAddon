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
#include "ARRLISTTOOLS.h"
#include "SURFTOPOSEG.h"
#include "GEOMTOOLS.h"
#include <string.h>



/**
 * Determine if the SurfTopoSeg handle is sane.
 *
 * param SurfTopoSeg
 *     SurfTopoSeg structure in question.
 *
 * return
 *     TRUE if the SurfTopoSeg structure is valid, otherwise FALSE.
 */
Boolean_t SurfTopoSegIsValid(SurfTopoSeg_pa SurfTopoSeg)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(SurfTopoSeg) &&
               VALID_REF(SurfTopoSeg->SurfElemsRemaining) && ArrListIsValid(SurfTopoSeg->SurfElemsRemaining) &&
               VALID_REF(SurfTopoSeg->SurfSegZoneNums) && ArrListIsValid(SurfTopoSeg->SurfSegZoneNums));


    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the SurfTopoSeg handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a SurfTopoSeg handle.
 */
void SurfTopoSegDealloc(SurfTopoSeg_pa *SurfTopoSeg)
{
    REQUIRE(VALID_REF(SurfTopoSeg));
    REQUIRE(SurfTopoSegIsValid(*SurfTopoSeg) || *SurfTopoSeg == NULL);

    if (*SurfTopoSeg != NULL)
    {
        /* release the ArrList's */
        if ((*SurfTopoSeg)->SurfElemsRemaining != NULL) ArrListDealloc(&((*SurfTopoSeg)->SurfElemsRemaining));
        if ((*SurfTopoSeg)->SurfSegZoneNums != NULL) ArrListDealloc(&((*SurfTopoSeg)->SurfSegZoneNums));

        /* release the list structure itself */
        FREE_ITEM(*SurfTopoSeg, "SurfTopoSeg structure");
        *SurfTopoSeg = NULL;
    }

    ENSURE(*SurfTopoSeg == NULL);
}





/**
 * Empties the SurfTopoSeg structure.
 *
 * param SurfTopoSeg
 *     SurfTopoSeg to clear.
 */
void SurfTopoSegClear(SurfTopoSeg_pa SurfTopoSeg)
{
    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));

    ArrListClear(SurfTopoSeg->SurfElemsRemaining);
    ArrListClear(SurfTopoSeg->SurfSegZoneNums);

    ENSURE(SurfTopoSegIsValid(SurfTopoSeg) && SurfTopoSegGetElemCount(SurfTopoSeg) == 0);
    ENSURE(SurfTopoSegGetSegZoneCount(SurfTopoSeg) == 0);
}






/**
 * Allocates a SurfTopoSeg handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * return
 *     SurfTopoSeg handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
SurfTopoSeg_pa SurfTopoSegAlloc()
{
    SurfTopoSeg_pa Result = NULL;

    Result = ALLOC_ITEM(SurfTopoSeg_s, "SurfTopoSeg structure");
    if (Result != NULL)
    {
        Result->SurfElemsRemaining  = ArrListAlloc(120, ArrListType_Long);
        Result->SurfSegZoneNums     = ArrListAlloc(20, ArrListType_Long);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->SurfElemsRemaining == NULL || Result->SurfSegZoneNums == NULL)
        {
            ArrListDealloc(&(Result->SurfElemsRemaining));
            ArrListDealloc(&(Result->SurfSegZoneNums));
            FREE_ITEM(Result, "SurfTopoSeg structure");
            Result = NULL;
        }
    }

    ENSURE(SurfTopoSegIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of remaining Elements currently in the SurfTopoSeg
 * structure.
 *
 * param
 *     SurfTopoSeg structure in question.
 *
 * return
 *     Number of remaining elements in the SurfTopoSeg structure.
 */
LgIndex_t SurfTopoSegGetElemCount(const SurfTopoSeg_pa SurfTopoSeg)
{
    LgIndex_t Result = 0;

    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));

    Result = ArrListGetCount(SurfTopoSeg->SurfElemsRemaining);

    ENSURE(Result >= 0);
    return Result;
}



/**
 * Gets the number of segment zones currently in the SurfTopoSeg
 * structure.
 *
 * param
 *     SurfTopoSeg structure in question.
 *
 * return
 *     Number of segment zones in the SurfTopoSeg structure.
 */
LgIndex_t SurfTopoSegGetSegZoneCount(const SurfTopoSeg_pa SurfTopoSeg)
{
    LgIndex_t Result = 0;

    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));

    Result = ArrListGetCount(SurfTopoSeg->SurfSegZoneNums);

    ENSURE(Result >= 0);
    return Result;
}











/**
 * Inserts SurfTopoSeg remaining element number at the specified
 * offset. The arrays will be expanded to accomodate the additional value.
 *
 *
 * param SurfTopoSeg
 *     SurfTopoSeg target in which to set the coordinates.
 * param ElemOffset
 *     Offset at which to insert the ElemNum.
 * param ElemNum
 *     Number of an element in the Isosurface zone to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t SurfTopoSegInsertElemAtOffset(SurfTopoSeg_pa SurfTopoSeg,
                                        LgIndex_t      ElemOffset,
                                        LgIndex_t      ElemNum)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));
    REQUIRE(0 <= ElemOffset && ElemOffset <= SurfTopoSegGetElemCount(SurfTopoSeg));

    Item.Long = ElemNum;
    IsOk = ArrListInsertItem(SurfTopoSeg->SurfElemsRemaining, ElemOffset, Item);

    ENSURE(SurfTopoSegIsValid(SurfTopoSeg));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Appends ElemNum to SurfTopoSeg structure. The array lists
 * will be expanded to accommodate the additional items.
 *
 * param SurfTopoSeg
 *     SurfTopoSeg target to which the ElemNum is to be appended.
 * param ElemNum
 *     Element number, in the isosurface zone, to be appended to SurfTopoSeg.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t SurfTopoSegAppendElem(SurfTopoSeg_pa  SurfTopoSeg,
                                LgIndex_t       ElemNum)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));

    Count = SurfTopoSegGetElemCount(SurfTopoSeg);

    IsOk = SurfTopoSegInsertElemAtOffset(SurfTopoSeg, Count, ElemNum);

    ENSURE(SurfTopoSegIsValid(SurfTopoSeg));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}




/**
 * Remove element at the specified offset.
 *
 * param SurfTopoSeg
 *     SurfTopoSeg target from which to remove the element.
 * param Offset
 *     Offset of the element.
 *
 * return
 *     Element number removed if it worked, -1 otherwise.
 */

LgIndex_t SurfTopoSegRemoveElem(SurfTopoSeg_pa  SurfTopoSeg,
                                LgIndex_t       Offset)
{
    LgIndex_t Result = -1;
    ArrListItem_u Item;

    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));
    REQUIRE(Offset >= 0 && Offset < SurfTopoSegGetElemCount(SurfTopoSeg));

    Item = ArrListRemoveItem(SurfTopoSeg->SurfElemsRemaining, Offset);

    if (Item.Long > 0) Result = Item.Long;

    ENSURE(SurfTopoSegIsValid(SurfTopoSeg));

    ENSURE(Result > 0 || Result == -1);
    return Result;
}





/**
 * Remove element of specified number from SurfElemsRemaining. This
 * requires searching through the ArrList to find the correct offset.
 *
 * param SurfTopoSeg
 *     SurfTopoSeg target from which to remove the element.
 * param ElemNumToRemove
 *     The numver of the element to remove from SurfElemsRemaining.
 *
 * return
 *     Element number removed if it worked, -1 otherwise.
 */

LgIndex_t SurfTopoSegRemoveElemNum(SurfTopoSeg_pa  SurfTopoSeg,
                                   LgIndex_t       ElemNumToRemove)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Result = -1;
    LgIndex_t Offset = -1;
    // LgIndex_t nn;
    LgIndex_t NumElems = SurfTopoSegGetElemCount(SurfTopoSeg);
    Boolean_t IsFound = FALSE;
    ArrListItem_u Item;

    // Use native array for speed
    // LgIndex_t *SurfElemsRemaining = NULL;
    // LgIndex_t Count;
    // IsOk = ArrListToNative(SurfTopoSeg->SurfElemsRemaining, &SurfElemsRemaining, &Count);


    // Binary search parameters
    LgIndex_t Mid;
    LgIndex_t Low = 0;
    LgIndex_t High = NumElems;


    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));
    REQUIRE(ElemNumToRemove > 0);

    // Search to find the offset of the specified element
    // TODO: This is a sorted list - try binary search
    /*
    for (nn=0; !IsFound && nn<NumElems; nn++)
      {
        if (SurfElemsRemaining[nn] == ElemNumToRemove)
          {
            IsFound = TRUE;
            Offset = nn;
          }
        // Original
        // Item = ArrListGetItem(SurfTopoSeg->SurfElemsRemaining, nn);
        // if (Item.Long == ElemNumToRemove)
        //   {
        //    IsFound = TRUE;
        //    Offset = nn;
        //  }
        //
      }
    */

    // Binary Search to find the offset of the specified element
    while (Low < High)
    {
        Mid = Low + ((High - Low) / 2);
        Item = ArrListGetItem(SurfTopoSeg->SurfElemsRemaining, Mid);
        if (Item.Long < ElemNumToRemove)
            Low = Mid + 1;
        else
            //can't be High = Mid-1: here SurfElemsRemaining[Mid] >= ElemNumToRemove,
            //so High can't be < Mid if SurfElemsRemaining[Mid] == ElemNumToRemove
            High = Mid;
    }

    Item = ArrListGetItem(SurfTopoSeg->SurfElemsRemaining, Low);
    // if (Low < NumElems && SurfElemsRemaining[Low] == ElemNumToRemove)
    if (Low < NumElems && Item.Long == ElemNumToRemove)
        Offset = Low;


    // Remove the element from the list
    Item = ArrListRemoveItem(SurfTopoSeg->SurfElemsRemaining, Offset);

    if (Item.Long > 0) Result = Item.Long;

    ENSURE(SurfTopoSegIsValid(SurfTopoSeg));

    ENSURE(Result > 0 || Result == -1);
    return Result;
}







/**
 * Gets the ElemNum from the SurfTopoSeg structure for the
 * element at the specified offset.
 *
 * param SurfTopoSeg
 *     SurfTopoSeg structure containing the desired item.
 * param ElemOffset
 *     Offset to the IsoElemNum of the SurfTopoSeg.
 *
 * return
 *     ElemNum if it works, -1 otherwise.
 */

LgIndex_t SurfTopoSegGetElemNumFromOffset(const SurfTopoSeg_pa  SurfTopoSeg,
                                          const LgIndex_t       ElemOffset)
{
    LgIndex_t ElemNum = -1;
    ArrListItem_u Item;

    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));
    REQUIRE(ElemOffset >= 0 && ElemOffset < SurfTopoSegGetElemCount(SurfTopoSeg));

    Item = ArrListGetItem(SurfTopoSeg->SurfElemsRemaining, ElemOffset);
    ElemNum = Item.Long;

    ENSURE(ElemNum >= -1);
    return ElemNum;
}










/**
 * Given a surface zone, initialize the SurfTopoSeg structure.
 *
 * param SurfTopoSeg
 *     SurfTopoSeg structure containing the input points for the curve fit.
 * param SurfZone
 *     Number of the surface zone to be segmented.
 *
 * return
 *     TRUE if it worked, FALSE if failed.
 */
Boolean_t    SurfTopoSegInitialize(SurfTopoSeg_pa SurfTopoSeg,
                                   EntIndex_t     SurfZone)
{
    Boolean_t  IsOk = TRUE;
    EntIndex_t NumZones, NumVars;
    LgIndex_t  NumNodes, NumElems;
    LgIndex_t  ne;

    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);
    REQUIRE(SurfZone > 0 && SurfZone <= NumZones);

    // Set the SurfZone number in the structure
    SurfTopoSeg->SurfZone = SurfZone;

    // Find the number of Elements in the IsoSurfZone
    TecUtilZoneGetInfo(SurfZone, &NumNodes, &NumElems,
                       NULL, NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL);

    for (ne = 1; IsOk && ne <= NumElems; ne++)
        IsOk = SurfTopoSegAppendElem(SurfTopoSeg, ne);

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}







/**
 * Extract the topologically independent segment of the specified
 * extracted isosurface zone which is connected, topologically, to
 * the specified first element.
 *
 * param SurfTopoSeg
 *     SurfTopoSeg structure containing the input points for the curve fit.
 * param IsoSurfZone
 *     Number of the extraced isosurface zone.
 * param FirstElemNumOffset
 *     Offset of element number, within the SurfElemsRemaining, of the
 *     seed element.
 *
 * return
 *     Zone number of the new segment zone (if it worked), zero if it failed.
 */
EntIndex_t    SurfTopoSegComputeSeg(SurfTopoSeg_pa SurfTopoSeg,
                                    LgIndex_t      FirstElemNumOffset)
{
    EntIndex_t SurfTopoSegZone = 0;
    EntIndex_t NumZones, NumVars;
    EntIndex_t SurfZone;
    ZoneType_e SurfZoneType;
    Boolean_t  IsOk = TRUE;
    LgIndex_t  NumNodes, NumElems;
    Boolean_t  ShowStatusBar = TRUE;
    Boolean_t  IsStop = FALSE;

    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));

    SurfZone = SurfTopoSeg->SurfZone;

    /* Inform Tecplot that major data operation is beginning */
    TecUtilDataLoadBegin();

    /* Lock Tecplot */
    TecUtilLockStart(AddOnID);


    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    REQUIRE(SurfZone > 0 && SurfZone <= NumZones);

    SurfZoneType = TecUtilZoneGetType(SurfZone);
    REQUIRE(SurfZoneType == ZoneType_FETriangle || SurfZoneType == ZoneType_FEQuad);

    REQUIRE(FirstElemNumOffset >= 0 && FirstElemNumOffset < SurfTopoSegGetElemCount(SurfTopoSeg));


    // Find the number of Elements in the SurfZone
    TecUtilZoneGetInfo(SurfZone, &NumNodes, &NumElems,
                       NULL, NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL);


    // Start with the First cell and add all connected cells
    if (IsOk)
    {
        int nn;
        int NodesPerElem = 3;

        ArrList_pa ElemList = ArrListAlloc(60, ArrListType_Long);
        // Segment node numbers for segment elements
        ArrList_pa ElemNd1List = ArrListAlloc(60, ArrListType_Long);
        ArrList_pa ElemNd2List = ArrListAlloc(60, ArrListType_Long);
        ArrList_pa ElemNd3List = ArrListAlloc(60, ArrListType_Long);
        ArrList_pa ElemNd4List = ArrListAlloc(60, ArrListType_Long);

        ArrList_pa ElemBoundaryList = ArrListAlloc(60, ArrListType_Long);
        ArrList_pa NodeList = ArrListAlloc(60, ArrListType_Long);
        ArrListItem_u Item;

        // Create temporary reverse-map arrays for elements and nodes.
        // 0==unused, otherwise seg elem/node number (at the surf elem offset)
        LgIndex_t *ElemNumInSeg = NULL;
        LgIndex_t *NodeNumInSeg = NULL;
        ElemNumInSeg = ALLOC_ARRAY(NumElems + 1, LgIndex_t, "Element reverse map");
        NodeNumInSeg = ALLOC_ARRAY(NumNodes + 1, LgIndex_t, "Node reverse map");
        // Initialize reverse maps to 0 (not used yet)
        for (nn = 0; nn <= NumElems; nn++)
            ElemNumInSeg[nn] = 0;
        for (nn = 0; nn <= NumNodes; nn++)
            NodeNumInSeg[nn] = 0;

        if (SurfZoneType == ZoneType_FEQuad)
            NodesPerElem = 4;
        else if (SurfZoneType == ZoneType_FETriangle)
            NodesPerElem = 3;
        else
            IsOk = FALSE;


        // Add first element and initialize NodeBoundaryList with it's nodes
        if (IsOk)
        {
            // Get the first element number
            LgIndex_t FirstElem = SurfTopoSegGetElemNumFromOffset(SurfTopoSeg, FirstElemNumOffset);

            // Add first element to ElemList and ElemBoundaryList
            Item.Long = FirstElem;
            IsOk = ArrListAppendItem(ElemList, Item);
            ElemNumInSeg[FirstElem] = 1;
            IsOk = ArrListAppendItem(ElemBoundaryList, Item);

            // Add unique node numbers to the NodeList
            for (nn = 1; nn <= NodesPerElem; nn++)
            {
                LgIndex_t SegNodeNum, SurfNodeNum;
                SurfNodeNum = TecUtilDataNodeGetByZone(SurfZone, FirstElem, nn);

                // TEMPORARY try reverse map
                if (NodeNumInSeg[SurfNodeNum] == 0)
                {
                    Item.Long = SurfNodeNum;
                    IsOk = ArrListAppendItem(NodeList, Item);
                    SegNodeNum = ArrListGetCount(NodeList);
                    NodeNumInSeg[SurfNodeNum] = SegNodeNum;
                }
                else
                    SegNodeNum = NodeNumInSeg[SurfNodeNum];

                // SegNodeNum = ArrListAppendUniqueLongItem(NodeList, SurfNodeNum) + 1;
                if (SegNodeNum > 0)
                {
                    Item.Long = SegNodeNum;
                    switch (nn)
                    {
                        case 1:
                            IsOk = ArrListAppendItem(ElemNd1List, Item);
                            break;
                        case 2:
                            IsOk = ArrListAppendItem(ElemNd2List, Item);
                            break;
                        case 3:
                            IsOk = ArrListAppendItem(ElemNd3List, Item);
                            break;
                        case 4:
                            IsOk = ArrListAppendItem(ElemNd4List, Item);
                            break;
                    }
                }
            }
        }


        // Loop over Elements in ElemBoundaryList, each time adding the unique elements
        // to the ElemList and ElemBoundaryList and adding the unique nodes from
        // those elements to the NodeList. Continue until there are no more
        // nodes in the NodeBoundaryList. Only works for 2-manifold surfaces.
        while (IsOk && !IsStop && ArrListGetCount(ElemBoundaryList) > 0)
        {
            int Face;
            LgIndex_t ElemToRemove, SegElemNum;

            // Select first node in list
            Item = ArrListGetItem(ElemBoundaryList, 0);
            ElemToRemove = Item.Long;

            // Loop over faces, finding face neighbors and adding them
            for (Face = 1; Face <= NodesPerElem; Face++)
            {
                LgIndex_t  NeighborElem;
                FaceNeighbor_pa FNPtr = TecUtilDataFaceNbrGetRef(SurfZone);

                LgIndex_t NumElemsInList = ArrListGetCount(ElemList);
                LgIndex_t NewNumElemsInList;

                // Get the FaceNeighbor
                if (FNPtr != NULL)
                {
                    EntIndex_t NeighborZone;
                    LgIndex_t NumNeighbors = TecUtilDataFaceNbrGetNumNByRef(FNPtr, ElemToRemove, Face);

                    if (NumNeighbors > 0)
                        TecUtilDataFaceNbrGetNbrByRef(FNPtr, ElemToRemove, Face, 1,
                                                      &NeighborElem,
                                                      &NeighborZone);
                    else NeighborElem = -1;

                    // Add element to ElemList and ElemBoundaryList
                    if (NumNeighbors > 0 && (NeighborZone == SurfZone || NeighborZone == 0))
                        // SegElemNum = ArrListAppendUniqueLongItem(ElemList, NeighborElem) + 1;
                        // TEMPORARY try reverse map
                    {
                        if (ElemNumInSeg[NeighborElem] == 0)
                        {
                            Item.Long = NeighborElem;
                            IsOk = ArrListAppendItem(ElemList, Item);
                            SegElemNum = ArrListGetCount(ElemList);
                            ElemNumInSeg[NeighborElem] = SegElemNum;
                        }
                        else
                            SegElemNum = ElemNumInSeg[NeighborElem];
                    }

                    else
                        SegElemNum = -1;

                    // If the Elem is unique in the list, add Elem to ElemBoundaryList
                    // and add the elements nodes to the respective node lists.
                    NewNumElemsInList = ArrListGetCount(ElemList);
                    if (NewNumElemsInList > NumElemsInList)
                    {
                        LgIndex_t EBLElemNum = ArrListAppendUniqueLongItem(ElemBoundaryList, NeighborElem);

                        // Add unique node numbers from element to the NodeList
                        for (nn = 1; nn <= NodesPerElem; nn++)
                        {
                            LgIndex_t SegNodeNum, SurfNodeNum;
                            SurfNodeNum = TecUtilDataNodeGetByZone(SurfZone, NeighborElem, nn);
                            // TEMPORARY try reverse map
                            if (NodeNumInSeg[SurfNodeNum] == 0)
                            {
                                Item.Long = SurfNodeNum;
                                IsOk = ArrListAppendItem(NodeList, Item);
                                SegNodeNum = ArrListGetCount(NodeList);
                                NodeNumInSeg[SurfNodeNum] = SegNodeNum;
                            }
                            else
                                SegNodeNum = NodeNumInSeg[SurfNodeNum];

                            // SegNodeNum = ArrListAppendUniqueLongItem(NodeList, SurfNodeNum) + 1;
                            if (SegNodeNum > 0)
                            {
                                Item.Long = SegNodeNum;
                                switch (nn)
                                {
                                    case 1:
                                        IsOk = ArrListAppendItem(ElemNd1List, Item);
                                        break;
                                    case 2:
                                        IsOk = ArrListAppendItem(ElemNd2List, Item);
                                        break;
                                    case 3:
                                        IsOk = ArrListAppendItem(ElemNd3List, Item);
                                        break;
                                    case 4:
                                        IsOk = ArrListAppendItem(ElemNd4List, Item);
                                        break;
                                }
                            }
                        }
                    }
                }
            }

            // Remove the ElemToRemove from ElemBoundaryList
            Item = ArrListRemoveItem(ElemBoundaryList, 0);
            if (SurfTopoSegRemoveElemNum(SurfTopoSeg, ElemToRemove) != ElemToRemove) IsOk = FALSE;

            /* Update status bar */
            /*if (ShowStatusBar)
            {
                char PercentDoneText[200];
                int  PercentDone;
                sprintf_s(PercentDoneText, "Computing topo segmentation %d", SurfZone);
                TecUtilStatusSetPercentDoneText(PercentDoneText);
                PercentDone = (int)((100 * ArrListGetCount(ElemList)) / NumElems);
                IsStop = !TecUtilStatusCheckPercentDone(PercentDone);
            }*/
        }

        if (IsOk && !IsStop)
        {
            // Create the new Tecplot zone
            LgIndex_t  SegNumNodes = ArrListGetCount(NodeList);
            LgIndex_t  SegNumElems = ArrListGetCount(ElemList);
            FieldDataType_e *VarType = NULL;

            // Get the field data types of the source surface zone, using in new zone
            if (IsOk)
            {
                EntIndex_t nv;
                VarType = ALLOC_ARRAY(NumVars, FieldDataType_e, "VarType array");
                if (VarType != NULL)
                {
                    for (nv = 1; nv <= NumVars; nv++)
                    {
                        FieldDataType_e NVFDType = TecUtilDataValueGetType(SurfZone, nv);
                        VarType[nv-1] = NVFDType;
                    }
                }
            }

            // Create new zone for segment
            if (TecUtilDataSetAddZone("SurfTopoSegZone", SegNumNodes, SegNumElems, 1,
                                      SurfZoneType, VarType))
            {
                LgIndex_t ne, nn;
                EntIndex_t nv;
                NodeMap_pa NodeMapSurfZone = NULL;
                NodeMap_pa NodeMapSegZone = NULL;
                Set_pa ZonesAdded = TecUtilSetAlloc(TRUE);

                // new zone is always last zone
                TecUtilDataSetGetInfo(NULL, &SurfTopoSegZone, NULL);

                // Set the node map for the new sub-zone
                NodeMapSurfZone = TecUtilDataNodeGetWritableRef(SurfZone);
                NodeMapSegZone  = TecUtilDataNodeGetWritableRef(SurfTopoSegZone);
                for (ne = 1; ne <= SegNumElems; ne++)
                {
                    LgIndex_t SegElemOffset = ne - 1;
                    NodeMap_t Node;

                    // Set the SurfTopoSeg node map for element
                    // Node 1
                    Item = ArrListGetItem(ElemNd1List, SegElemOffset);
                    Node = Item.Long;
                    TecUtilDataNodeSetByRef(NodeMapSegZone, ne, 1, Node);
                    // Node 2
                    Item = ArrListGetItem(ElemNd2List, SegElemOffset);
                    Node = Item.Long;
                    TecUtilDataNodeSetByRef(NodeMapSegZone, ne, 2, Node);
                    // Node 3
                    Item = ArrListGetItem(ElemNd3List, SegElemOffset);
                    Node = Item.Long;
                    TecUtilDataNodeSetByRef(NodeMapSegZone, ne, 3, Node);
                    // Node 4
                    if (NodesPerElem > 3)
                    {
                        Item = ArrListGetItem(ElemNd4List, SegElemOffset);
                        Node = Item.Long;
                        TecUtilDataNodeSetByRef(NodeMapSegZone, ne, 4, Node);
                    }
                }

                // fill new zone with values for all variables from old zone
                for (nn = 1; nn <= SegNumNodes; nn++)
                {
                    LgIndex_t SurfNode;
                    LgIndex_t SegNodeOffset = nn - 1;
                    Item = ArrListGetItem(NodeList, SegNodeOffset);
                    SurfNode = Item.Long;
                    for (nv = 1; nv <= NumVars; nv++)
                    {
                        // TODO: replace with ByRef functions
                        double Value = TecUtilDataValueGetByZoneVar(SurfZone, nv, SurfNode);
                        IsOk = TecUtilDataValueSetByZoneVar(SurfTopoSegZone, nv, nn, Value);
                    }
                }


                // inform Tecplot of new zone
                TecUtilSetAddMember(ZonesAdded, SurfTopoSegZone, TRUE);
                TecUtilStateChanged(StateChange_ZonesAdded,
                                    (ArbParam_t)ZonesAdded);
                TecUtilSetDealloc(&ZonesAdded);
            }

            if (VarType != NULL) FREE_ARRAY(VarType, "VarType array");
        }

        // Clean up
        ArrListDealloc(&ElemList);
        ArrListDealloc(&NodeList);
        ArrListDealloc(&ElemNd1List);
        ArrListDealloc(&ElemNd2List);
        ArrListDealloc(&ElemNd3List);
        ArrListDealloc(&ElemNd4List);
        ArrListDealloc(&ElemBoundaryList);

        // Free temporary reverse-map arrays for elements and nodes
        FREE_ARRAY(ElemNumInSeg, "Element reverse map");
        FREE_ARRAY(NodeNumInSeg, "Node reverse map");
    }

    /* Finish lock for this function (may be nested) */
    TecUtilLockFinish(AddOnID);

    /* Inform Tecplot that major data operation is ending */
    TecUtilDataLoadEnd();

    REQUIRE(SurfTopoSegIsValid(SurfTopoSeg));
    ENSURE(VALID_BOOLEAN(IsOk));
    REQUIRE(SurfTopoSegZone == 0 || SurfTopoSegZone == NumZones + 1);
    return SurfTopoSegZone;
}




/**
 * Unit test for SurfTopoSegComputeSeg. Create a surface zone with two
 * blocks (six elements, eight nodes) and perform the surface segmentation.
 * The result should be two new zones containing the individual blocks.
 *
 * return
 *     TRUE if the test succeeds, FALSE if it fails.
 */

Boolean_t SurfTopoSegTest()
{
    Boolean_t  IsOk = TRUE;
    EntIndex_t NumZones, NumVars, NewZone;
    LgIndex_t  NumNodes, NumElems;
    double     XYZTol = 1.0e-6;

    /* Inform Tecplot that major data operation is beginning */
    TecUtilDataLoadBegin();

    /* Lock Tecplot */
    TecUtilLockStart(AddOnID);

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    // Create a temporary quad zone
    IsOk = TecUtilDataSetAddZone("Temporary FEQuad Zone", 16, 12, 1,
                                 ZoneType_FEQuad, NULL);
    if (IsOk)
    {
        NodeMap_pa NodeMap = NULL;
        FieldData_pa XFDPtr, YFDPtr, ZFDPtr;
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

        TecUtilDataNodeSetByRef(NodeMap, 7, 1, 9);
        TecUtilDataNodeSetByRef(NodeMap, 7, 2, 10);
        TecUtilDataNodeSetByRef(NodeMap, 7, 3, 11);
        TecUtilDataNodeSetByRef(NodeMap, 7, 4, 12);

        TecUtilDataNodeSetByRef(NodeMap, 8, 1, 9);
        TecUtilDataNodeSetByRef(NodeMap, 8, 2, 10);
        TecUtilDataNodeSetByRef(NodeMap, 8, 3, 14);
        TecUtilDataNodeSetByRef(NodeMap, 8, 4, 13);

        TecUtilDataNodeSetByRef(NodeMap, 9, 1, 10);
        TecUtilDataNodeSetByRef(NodeMap, 9, 2, 14);
        TecUtilDataNodeSetByRef(NodeMap, 9, 3, 15);
        TecUtilDataNodeSetByRef(NodeMap, 9, 4, 11);

        TecUtilDataNodeSetByRef(NodeMap, 10, 1, 12);
        TecUtilDataNodeSetByRef(NodeMap, 10, 2, 11);
        TecUtilDataNodeSetByRef(NodeMap, 10, 3, 15);
        TecUtilDataNodeSetByRef(NodeMap, 10, 4, 16);

        TecUtilDataNodeSetByRef(NodeMap, 11, 1, 9);
        TecUtilDataNodeSetByRef(NodeMap, 11, 2, 13);
        TecUtilDataNodeSetByRef(NodeMap, 11, 3, 16);
        TecUtilDataNodeSetByRef(NodeMap, 11, 4, 12);

        TecUtilDataNodeSetByRef(NodeMap, 12, 1, 13);
        TecUtilDataNodeSetByRef(NodeMap, 12, 2, 14);
        TecUtilDataNodeSetByRef(NodeMap, 12, 3, 15);
        TecUtilDataNodeSetByRef(NodeMap, 12, 4, 16);

        // Get FieldData pointers to X, Y, and Z vars.
        TecUtilZoneGetInfo(NewZone, &NumNodes, &NumElems, NULL, &XFDPtr, &YFDPtr, &ZFDPtr,
                           NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        //
        //
        //     0,1,1  -----------------  1,1,1
        //          /                / |
        //        /   |            /   |
        //  0,0,1 ----------------  1,1|0
        //       |    |           |    |
        //       |                |    |
        //      0|1,0 | _  _  _ _ | _ _| 1,1,0
        //       |   /            |   /
        //       | /              | /
        //  0,0,0 ---------------- 1,0,0
        //

        // Set X-Values for nodes
        TecUtilDataValueSetByRef(XFDPtr, 1, 0.0);
        TecUtilDataValueSetByRef(XFDPtr, 2, 1.0);
        TecUtilDataValueSetByRef(XFDPtr, 3, 1.0);
        TecUtilDataValueSetByRef(XFDPtr, 4, 0.0);
        TecUtilDataValueSetByRef(XFDPtr, 5, 0.0);
        TecUtilDataValueSetByRef(XFDPtr, 6, 1.0);
        TecUtilDataValueSetByRef(XFDPtr, 7, 1.0);
        TecUtilDataValueSetByRef(XFDPtr, 8, 0.0);
        TecUtilDataValueSetByRef(XFDPtr, 9, 2.0);
        TecUtilDataValueSetByRef(XFDPtr, 10, 3.0);
        TecUtilDataValueSetByRef(XFDPtr, 11, 3.0);
        TecUtilDataValueSetByRef(XFDPtr, 12, 2.0);
        TecUtilDataValueSetByRef(XFDPtr, 13, 2.0);
        TecUtilDataValueSetByRef(XFDPtr, 14, 3.0);
        TecUtilDataValueSetByRef(XFDPtr, 15, 3.0);
        TecUtilDataValueSetByRef(XFDPtr, 16, 2.0);

        // Set Y-Values for nodes
        TecUtilDataValueSetByRef(YFDPtr, 1, 0.0);
        TecUtilDataValueSetByRef(YFDPtr, 2, 0.0);
        TecUtilDataValueSetByRef(YFDPtr, 3, 1.0);
        TecUtilDataValueSetByRef(YFDPtr, 4, 1.0);
        TecUtilDataValueSetByRef(YFDPtr, 5, 0.0);
        TecUtilDataValueSetByRef(YFDPtr, 6, 0.0);
        TecUtilDataValueSetByRef(YFDPtr, 7, 1.0);
        TecUtilDataValueSetByRef(YFDPtr, 8, 1.0);
        TecUtilDataValueSetByRef(YFDPtr, 9, 0.0);
        TecUtilDataValueSetByRef(YFDPtr, 10, 0.0);
        TecUtilDataValueSetByRef(YFDPtr, 11, 1.0);
        TecUtilDataValueSetByRef(YFDPtr, 12, 1.0);
        TecUtilDataValueSetByRef(YFDPtr, 13, 0.0);
        TecUtilDataValueSetByRef(YFDPtr, 14, 0.0);
        TecUtilDataValueSetByRef(YFDPtr, 15, 1.0);
        TecUtilDataValueSetByRef(YFDPtr, 16, 1.0);

        // Set Z-Values for nodes
        TecUtilDataValueSetByRef(ZFDPtr, 1, 0.0);
        TecUtilDataValueSetByRef(ZFDPtr, 2, 0.0);
        TecUtilDataValueSetByRef(ZFDPtr, 3, 0.0);
        TecUtilDataValueSetByRef(ZFDPtr, 4, 0.0);
        TecUtilDataValueSetByRef(ZFDPtr, 5, 1.0);
        TecUtilDataValueSetByRef(ZFDPtr, 6, 1.0);
        TecUtilDataValueSetByRef(ZFDPtr, 7, 1.0);
        TecUtilDataValueSetByRef(ZFDPtr, 8, 1.0);
        TecUtilDataValueSetByRef(ZFDPtr, 9, 0.0);
        TecUtilDataValueSetByRef(ZFDPtr, 10, 0.0);
        TecUtilDataValueSetByRef(ZFDPtr, 11, 0.0);
        TecUtilDataValueSetByRef(ZFDPtr, 12, 0.0);
        TecUtilDataValueSetByRef(ZFDPtr, 13, 1.0);
        TecUtilDataValueSetByRef(ZFDPtr, 14, 1.0);
        TecUtilDataValueSetByRef(ZFDPtr, 15, 1.0);
        TecUtilDataValueSetByRef(ZFDPtr, 16, 1.0);



        // inform Tecplot of new zone
        TecUtilSetAddMember(zones_added, NewZone, TRUE);
        TecUtilStateChanged(StateChange_ZonesAdded,
                            (ArbParam_t)zones_added);
        TecUtilSetDealloc(&zones_added);
    }



    // Extract the first segment of SurfZone
    if (IsOk)
    {
        SurfTopoSeg_pa SurfTopoSeg = NULL;
        EntIndex_t  SurfZone = NumZones + 1;
        EntIndex_t  SegZone1 = -1;
        EntIndex_t  SegZone2 = -1;
        LgIndex_t   NumNodes, NumElems;
        Boolean_t ShowStatusBar = TRUE;
        Boolean_t IsStop = FALSE;

        SurfTopoSeg = SurfTopoSegAlloc();


        if (SurfTopoSeg != NULL)
        {
            IsOk = SurfTopoSegInitialize(SurfTopoSeg, SurfZone);
            if (IsOk) SegZone1 = SurfTopoSegComputeSeg(SurfTopoSeg, 0);
            if (SegZone1 == 0) IsOk = FALSE;
        }

        // Check NumElems and NumNodes
        if (IsOk)
        {
            TecUtilZoneGetInfo(SegZone1, &NumNodes, &NumElems, NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            if (NumNodes != 8 || NumElems != 6) IsOk = FALSE;
        }

        //
        //      (8)n8 -----------------  n7 (7)     f6 is in back
        //          /                / |
        //        /   |  f4        /   |
        //     n4 ---------------- n3  |
        //       | f5 |           |    |
        //       |       f1       | f3 |
        //      (|6)n5| _  _  _ _ | _ _| n6 (5)
        //       |   /            |   /
        //       | /      f2      | /
        //     n1 ---------------- n2
        //
        // Check the NodeMap
        if (IsOk)
        {
            NodeMap_pa nm;
            nm = TecUtilDataNodeGetWritableRef(SegZone1);
            if (nm)
            {
                if (TecUtilDataNodeGetByRef(nm, 1, 1) != 1) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 1, 2) != 2) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 1, 3) != 3) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 1, 4) != 4) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 2, 1) != 1) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 2, 2) != 2) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 2, 3) != 5) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 2, 4) != 6) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 3, 1) != 2) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 3, 2) != 5) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 3, 3) != 7) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 3, 4) != 3) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 4, 1) != 4) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 4, 2) != 3) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 4, 3) != 7) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 4, 4) != 8) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 5, 1) != 1) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 5, 2) != 6) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 5, 3) != 8) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 5, 4) != 4) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 6, 1) != 6) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 6, 2) != 5) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 6, 3) != 7) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 6, 4) != 8) IsOk = FALSE;

            }
        }
        if (IsOk)
        {
            FieldData_pa XFDPtr = NULL;
            FieldData_pa YFDPtr = NULL;
            FieldData_pa ZFDPtr = NULL;

            // Get FieldData pointers to X, Y, and Z vars.
            TecUtilZoneGetInfo(SegZone1, &NumNodes, &NumElems, NULL, &XFDPtr, &YFDPtr, &ZFDPtr,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            // Test X-Values for nodes of SegZone1
            if (XFDPtr != NULL)
            {
                if (ABS(TecUtilDataValueGetByRef(XFDPtr, 1)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 2) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 3) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 4)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 5) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 6)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 7) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 8)) > XYZTol) IsOk = FALSE;
            }
            // Test Y-Values for nodes of SegZone1
            if (IsOk && YFDPtr != NULL)
            {
                if (ABS(TecUtilDataValueGetByRef(YFDPtr, 1)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 2)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 3) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 4) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 5)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 6)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 7) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 8) - 1.0) > XYZTol) IsOk = FALSE;
            }
            // Test Z-Values for nodes of SegZone1
            if (IsOk && ZFDPtr != NULL)
            {
                if (ABS(TecUtilDataValueGetByRef(ZFDPtr, 1)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 2)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 3)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 4)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 5) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 6) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 7) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 8) - 1.0) > XYZTol) IsOk = FALSE;
            }
        }


        // Second segment
        if (SurfTopoSeg != NULL)
        {
            SegZone2 = SurfTopoSegComputeSeg(SurfTopoSeg, 0);
        }


        // Check NumElems and NumNodes
        if (IsOk)
        {
            TecUtilZoneGetInfo(SegZone2, &NumNodes, &NumElems, NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            if (NumNodes != 8 || NumElems != 6) IsOk = FALSE;
        }

        // Check the NodeMap
        if (IsOk)
        {
            NodeMap_pa nm;
            nm = TecUtilDataNodeGetWritableRef(SegZone2);
            if (nm)
            {
                if (TecUtilDataNodeGetByRef(nm, 1, 1) != 1) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 1, 2) != 2) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 1, 3) != 3) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 1, 4) != 4) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 2, 1) != 1) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 2, 2) != 2) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 2, 3) != 5) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 2, 4) != 6) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 3, 1) != 2) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 3, 2) != 5) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 3, 3) != 7) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 3, 4) != 3) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 4, 1) != 4) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 4, 2) != 3) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 4, 3) != 7) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 4, 4) != 8) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 5, 1) != 1) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 5, 2) != 6) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 5, 3) != 8) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 5, 4) != 4) IsOk = FALSE;

                if (IsOk && TecUtilDataNodeGetByRef(nm, 6, 1) != 6) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 6, 2) != 5) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 6, 3) != 7) IsOk = FALSE;
                if (IsOk && TecUtilDataNodeGetByRef(nm, 6, 4) != 8) IsOk = FALSE;

            }
        }

        if (IsOk)
        {
            FieldData_pa XFDPtr = NULL;
            FieldData_pa YFDPtr = NULL;
            FieldData_pa ZFDPtr = NULL;

            // Get FieldData pointers to X, Y, and Z vars.
            TecUtilZoneGetInfo(SegZone2, &NumNodes, &NumElems, NULL, &XFDPtr, &YFDPtr, &ZFDPtr,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            // Test X-Values for nodes of SegZone1
            if (XFDPtr != NULL)
            {
                if (ABS(TecUtilDataValueGetByRef(XFDPtr, 1) - 2.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 2) - 3.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 3) - 3.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 4) - 2.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 5) - 3.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 6) - 2.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 7) - 3.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(XFDPtr, 8) - 2.0) > XYZTol) IsOk = FALSE;
            }
            // Test Y-Values for nodes of SegZone1
            if (IsOk && YFDPtr != NULL)
            {
                if (ABS(TecUtilDataValueGetByRef(YFDPtr, 1)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 2)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 3) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 4) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 5)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 6)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 7) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(YFDPtr, 8) - 1.0) > XYZTol) IsOk = FALSE;
            }
            // Test Z-Values for nodes of SegZone1
            if (IsOk && ZFDPtr != NULL)
            {
                if (ABS(TecUtilDataValueGetByRef(ZFDPtr, 1)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 2)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 3)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 4)) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 5) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 6) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 7) - 1.0) > XYZTol) IsOk = FALSE;
                if (IsOk && ABS(TecUtilDataValueGetByRef(ZFDPtr, 8) - 1.0) > XYZTol) IsOk = FALSE;
            }
        }


        SurfTopoSegDealloc(&SurfTopoSeg);


        // Delete the zones
        if (SurfZone > 0 && SegZone1 > 0 && SegZone2 > 0)
        {
            Set_pa ZoneList = TecUtilSetAlloc(FALSE);
            if (ZoneList)
            {
                TecUtilSetAddMember(ZoneList, SurfZone, FALSE);
                TecUtilSetAddMember(ZoneList, SegZone1, FALSE);
                TecUtilSetAddMember(ZoneList, SegZone2, FALSE);
                if (TecUtilDataSetDeleteZone(ZoneList))
                    TecUtilSetDealloc(&ZoneList);
            }
        }
    }

    /* Finish lock for this function (may be nested) */
    TecUtilLockFinish(AddOnID);

    /* Inform Tecplot that major data operation is ending */
    TecUtilDataLoadEnd();

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


/*
	Loops over all elements of a surface zone and returns
	the average of cell max extent (max distance from any
	node to any other node of a single cell), the number of
	elements in the surface, and the max distance from one
	element to another (diameter of the sphere).
*/
Boolean_t	SurfTopoDimensions(EntIndex_t	SurfZone,
							   const XYZ_s	Center,
							LgIndex_t	*ElemCount,
							double		*SurfRadius,
							double		*AvgMaxExtent)
{
	Boolean_t	IsOk = TRUE;

	REQUIRE(VALID_REF(ElemCount));
	REQUIRE(VALID_REF(SurfRadius));
	REQUIRE(VALID_REF(AvgMaxExtent));
	
	LgIndex_t   NumNodes, ne;
	XYZ_s		Nodes[4];
	double		MaxCellDist;
	int			AvgCountCellExtent = 0;

	*SurfRadius = -1.0;

	/* Inform Tecplot that major data operation is beginning */
	TecUtilDataLoadBegin();

	/* Lock Tecplot */
	TecUtilLockStart(AddOnID);

	REQUIRE(1 <= SurfZone && SurfZone <= TecUtilDataSetGetNumZones());

	// Find the number of Elements in the SurfZone
	TecUtilZoneGetInfo(SurfZone, &NumNodes, &*ElemCount,
		NULL, NULL, NULL, NULL, NULL, NULL,
		NULL, NULL, NULL, NULL, NULL);

	IsOk = (*ElemCount > 0);
	
	*AvgMaxExtent = 0.0;
	
	// Loop over elements of SurfZone to find the one closest to X,Y,Z
	for (ne = 1; IsOk && ne < *ElemCount; ne++)
	{
		//	First get information about the extent of current cell
		SurfTopoGetElemNodeXYZByElemNum(SurfZone, ne, Nodes);
	
		MaxCellDist = -1.0;
		for (int i = 0 ; i < 3 ; i++)
		{
			*SurfRadius = MAX(*SurfRadius, DistanceSquaredXYZ(Nodes[i], Center));
			for (int j = i+1 ; j < 4 ; j++)
				MaxCellDist = MAX(MaxCellDist, DistanceSquaredXYZ(Nodes[i], Nodes[j]));
		}

		*SurfRadius = MAX(*SurfRadius, DistanceSquaredXYZ(Nodes[3], Center));

		IsOk = (MaxCellDist >= 0);

		if (MaxCellDist > 0)
		{
			*AvgMaxExtent += sqrt(MaxCellDist);
			AvgCountCellExtent++;
		}
	}

	IsOk = (AvgCountCellExtent > 0);

	if (IsOk)
	{
		*AvgMaxExtent /= (double)AvgCountCellExtent;
		*SurfRadius = sqrt(*SurfRadius);
	}


	/* Finish lock for this function (may be nested) */
	TecUtilLockFinish(AddOnID);

	/* Inform Tecplot that major data operation is ending */
	TecUtilDataLoadEnd();

	REQUIRE(VALID_BOOLEAN(IsOk) && IsOk);
	return IsOk;
}


/*
	Sets XYZ variables to the xyz coordinates of the nodes
	for a specified surface element.
*/
ZoneType_e	SurfTopoGetElemNodeXYZByElemNum(EntIndex_t	SurfZone,
											LgIndex_t	ElemNum,
											XYZ_s*		Nodes)
{

	REQUIRE(1 <= SurfZone && SurfZone <= TecUtilDataSetGetNumZones());


	// Find the closest element to the point
	EntIndex_t   XVarNum, YVarNum, ZVarNum;
	FieldData_pa XVarFDPtr, YVarFDPtr, ZVarFDPtr;
	ZoneType_e   ZoneType;

	// FETriangle or FEQuad?
	ZoneType = TecUtilZoneGetType(SurfZone);

	// Get X, Y, Z variable numbers
	XVarNum = TecUtilVarGetNumByAssignment('X');
	YVarNum = TecUtilVarGetNumByAssignment('Y');
	ZVarNum = TecUtilVarGetNumByAssignment('Z');

	XVarFDPtr = TecUtilDataValueGetReadableRef(SurfZone, XVarNum);
	YVarFDPtr = TecUtilDataValueGetReadableRef(SurfZone, YVarNum);
	ZVarFDPtr = TecUtilDataValueGetReadableRef(SurfZone, ZVarNum);

	LgIndex_t NdNums[4];

	// Find unique node numbers
	for (int ii = 0; ii < 4; ii++)
		NdNums[ii] = TecUtilDataNodeGetByZone(SurfZone, ElemNum, ii + 1);

	for (int ii = 0 ; ii < 4 ; ii++)
	{
	Nodes[ii].X = TecUtilDataValueGetByRef(XVarFDPtr, NdNums[ii]);
	Nodes[ii].Y = TecUtilDataValueGetByRef(YVarFDPtr, NdNums[ii]);
	Nodes[ii].Z = TecUtilDataValueGetByRef(ZVarFDPtr, NdNums[ii]);
	}

	REQUIRE(ZoneType == ZoneType_FETriangle || ZoneType == ZoneType_FEQuad);
	return ZoneType;
}