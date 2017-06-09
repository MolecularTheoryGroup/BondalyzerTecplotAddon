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
#include "NORMALS.h"
#include <string.h>


/**
 * Determine if the Normals handle is sane.
 *
 * param Normals
 *     Normals structure in question.
 *
 * return
 *     TRUE if the Normals structure is valid, otherwise FALSE.
 */
Boolean_t NormalsIsValid(Normals_pa Normals)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(Normals) &&
               VALID_REF(Normals->Nx) && ArrListIsValid(Normals->Nx) &&
               VALID_REF(Normals->Ny) && ArrListIsValid(Normals->Ny) &&
               VALID_REF(Normals->Nz) && ArrListIsValid(Normals->Nz));

    // Unit test some functions
    // if (IsValid)
    // IsValid = NormalsTest();

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the Normals handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a Normals handle.
 */
void NormalsDealloc(Normals_pa *Normals)
{
    REQUIRE(VALID_REF(Normals));
    REQUIRE(NormalsIsValid(*Normals) || *Normals == NULL);

    if (*Normals != NULL)
    {
        /* release the ArrList's */
        if ((*Normals)->Nx != NULL) ArrListDealloc(&((*Normals)->Nx));
        if ((*Normals)->Ny != NULL) ArrListDealloc(&((*Normals)->Ny));
        if ((*Normals)->Nz != NULL) ArrListDealloc(&((*Normals)->Nz));

        /* release the list structure itself */
        FREE_ITEM(*Normals, "Normals structure");
        *Normals = NULL;
    }

    ENSURE(*Normals == NULL);
}





/**
 * Empties the Normals structure.
 *
 * param Normals
 *     Normals to clear.
 */
void NormalsClear(Normals_pa Normals)
{
    REQUIRE(NormalsIsValid(Normals));

    ArrListClear(Normals->Nx);
    ArrListClear(Normals->Ny);
    ArrListClear(Normals->Nz);

    ENSURE(NormalsIsValid(Normals) && ArrListGetCount(Normals->Nx) == 0 &&
           ArrListGetCount(Normals->Ny) == 0 && ArrListGetCount(Normals->Nz) == 0);
}






/**
 * Allocates a Normals handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * param
 *     NumNodes for the surface zone (if zero, it will still
 *     work suboptimally).
 *
 * return
 *     Normals handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
Normals_pa NormalsAlloc(const LgIndex_t NumNodes)
{
    Boolean_t  IsOk = TRUE;
    Normals_pa Result = NULL;
    LgIndex_t InitialSize = NumNodes + 1;

    if (InitialSize < 2) InitialSize = 60;

    Result = ALLOC_ITEM(Normals_s, "Normals structure");
    if (Result != NULL)
    {
        Result->Nx  = ArrListAlloc(InitialSize, ArrListType_Double);
        Result->Ny  = ArrListAlloc(InitialSize, ArrListType_Double);
        Result->Nz  = ArrListAlloc(InitialSize, ArrListType_Double);

        /* If it failed to allocate the array lists, clean-up and exit. */
        if (IsOk) IsOk = ArrListIsValid(Result->Nx);
        if (IsOk) IsOk = ArrListIsValid(Result->Ny);
        if (IsOk) IsOk = ArrListIsValid(Result->Nz);

        if (!IsOk)
        {
            ArrListDealloc(&(Result->Nx));
            ArrListDealloc(&(Result->Ny));
            ArrListDealloc(&(Result->Nz));
            FREE_ITEM(Result, "Normals structure");
            Result = NULL;
        }
    }

    ENSURE(NormalsIsValid(Result) || Result == NULL);
    return Result;
}





/**
 * Gets the number of nodes currently in the Normals structure.
 *
 * param
 *     Normals structure in question.
 *
 * return
 *     Number of nodes in the Normals structure.
 */
LgIndex_t NormalsGetNodeCount(const Normals_pa Normals)
{
    LgIndex_t Result = 0;

    REQUIRE(NormalsIsValid(Normals));

    // Subtract 1 because there is a final "offset" at the end to
    // conveniently handle the element count for the last node.
    Result = ArrListGetCount(Normals->Nx);

    ENSURE(Result >= 0);
    return Result;
}








/**
 * Places Normals components at
 * the specified offset. If the offset is beyond the end of the
 * array lists, they are sized accordingly and the intervening
 * array values between the last node of the original state and
 * the last node of the new state are assigned 0.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced.
 *
 * param Normals
 *     Normals target in which to set the coordinates and var value.
 * param PointOffset
 *     Offset of the point/node.
 * param Nx, Ny, Nz
 *     Components of the normal vector.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
/*
Boolean_t NormalsSetPointAtOffset(Normals_pa Normals,
                                  LgIndex_t  PointOffset,
                                  double     Nx,
                                  double     Ny,
                                  double     Nz)
{
  Boolean_t IsOk = TRUE;
  ArrListItem_u Item;


  REQUIRE(NormalsIsValid(Normals));
  REQUIRE(ElemOffset >= 0);

  Item.Long = IsoElemNum;
  IsOk = ArrListSetItem(Normals->IsoElemNums, ElemOffset, Item);

  ENSURE(NormalsIsValid(Normals));

  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}
*/




/**
 * Insert a normal vector at the specified node number.
 *
 * param Normals
 *     Normals target in which to set the coordinates.
 * param NodeNum
 *     Number of node for which the element number is added.
 * param ElemNum
 *     Components of normal vector
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t NormalsInsertForNode(Normals_pa Normals,
                               LgIndex_t  NodeNum,
                               XYZ_s      Normal)
{
    Boolean_t     IsOk = TRUE;
    ArrListItem_u Item;
    LgIndex_t     Offset;


    REQUIRE(NormalsIsValid(Normals));
    REQUIRE(0 <= NodeNum  && NodeNum <= NormalsGetNodeCount(Normals));

    Offset = NodeNum - 1;

    Item.Double = Normal.X;
    IsOk = ArrListInsertItem(Normals->Nx, Offset, Item);

    Item.Double = Normal.Y;
    IsOk = ArrListInsertItem(Normals->Ny, Offset, Item);

    Item.Double = Normal.Z;
    IsOk = ArrListInsertItem(Normals->Nz, Offset, Item);

    ENSURE(NormalsIsValid(Normals));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}








/**
 * Sets Normal vector to Normals structure.
 * The array lists may be expanded to accommodate the additional items.
 *
 * param Normals
 *     Normals target to which the Normal vector is to be appended.
 * param Normal
 *     Components of normal vector at node.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t NormalsSetForNode(Normals_pa  Normals,
                            LgIndex_t   NodeNum,
                            XYZ_s       Normal)
{
    Boolean_t     IsOk = TRUE;
    ArrListItem_u Item;
    LgIndex_t     Offset;


    REQUIRE(NormalsIsValid(Normals));
    REQUIRE(0 <= NodeNum);

    Offset = NodeNum - 1;

    Item.Double = Normal.X;
    IsOk = ArrListSetItem(Normals->Nx, Offset, Item);

    Item.Double = Normal.Y;
    IsOk = ArrListSetItem(Normals->Ny, Offset, Item);

    Item.Double = Normal.Z;
    IsOk = ArrListSetItem(Normals->Nz, Offset, Item);

    ENSURE(NormalsIsValid(Normals));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}




/**
 * Appends Normal vector to Normals structure.
 * The array lists will be expanded to accommodate the additional items.
 *
 * param Normals
 *     Normals target to which the Normal vector is to be appended.
 * param Nx, Ny, Nz
 *     Components of normal vector at node.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t NormalsAppendForNode(Normals_pa  Normals,
                               XYZ_s       Normal)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(NormalsIsValid(Normals));

    Count = NormalsGetNodeCount(Normals);

    IsOk = NormalsInsertForNode(Normals, Count + 1, Normal);

    ENSURE(NormalsIsValid(Normals));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Gets the components of the normal vector from the
 * Normals structure for the specified NodeNum.
 *
 * param Normals
 *     Normals structure containing the desired item.
 * param NodeNum
 *     Number of the node of interest.
 * param *Nx, *Ny, *Nz
 *     Pointers to components of normal vector.
 *
 * return
 *     Normal vector.
 */
XYZ_s NormalsGetNormalForNode(const Normals_pa  Normals,
                              const LgIndex_t   NodeNum)
{
    XYZ_s         Normal;
    ArrListItem_u Item;

    REQUIRE(NormalsIsValid(Normals));
    REQUIRE(NodeNum > 0 && NodeNum <= NormalsGetNodeCount(Normals));

    Item = ArrListGetItem(Normals->Nx, NodeNum - 1);
    Normal.X = Item.Double;

    Item = ArrListGetItem(Normals->Ny, NodeNum - 1);
    Normal.Y = Item.Double;

    Item = ArrListGetItem(Normals->Nz, NodeNum - 1);
    Normal.Z = Item.Double;

    return Normal;
}








/**
 * Gets the components of the normal vector for the specified
 * element of the specified FETriangle zone.
 *
 * param SurfZone
 *     Zone number.
 * param ElemNum
 *     Element number.
 *
 * return
 *     Normal vector.
 */
XYZ_s NormalsComputeTriNorm(EntIndex_t SurfZone,
                            LgIndex_t  ElemNum)
{
    XYZ_s         Normal;
    NodeMap_pa    NodeMap = NULL;
    FieldData_pa  XFDPtr = NULL;
    FieldData_pa  YFDPtr = NULL;
    FieldData_pa  ZFDPtr = NULL;
    EntIndex_t    XVarNum, YVarNum, ZVarNum;
    NodeMap_t     n1, n2, n3;

    REQUIRE(SurfZone > 0);
    REQUIRE(ElemNum > 0);

    // Normal Calculation
    NodeMap = TecUtilDataNodeGetWritableRef(SurfZone);
    if (NodeMap != NULL)
    {
        // TEMP
        {
            LgIndex_t IMax, JMax, KMax;
            TecUtilZoneGetInfo(SurfZone, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);

            if (ElemNum < 1 || ElemNum > JMax)
                TecUtilDialogMessageBox("Elem error in NormalsComputeTriNorm",
                                        MessageBox_Warning);
        }
        n1 = TecUtilDataNodeGetByRef(NodeMap, ElemNum, 1);
        n2 = TecUtilDataNodeGetByRef(NodeMap, ElemNum, 2);
        n3 = TecUtilDataNodeGetByRef(NodeMap, ElemNum, 3);
    }

    // Compute the cross-product of two edge 1-2 and edge 1-3
    XVarNum = TecUtilVarGetNumByAssignment('X');
    YVarNum = TecUtilVarGetNumByAssignment('Y');
    ZVarNum = TecUtilVarGetNumByAssignment('Z');
    XFDPtr  = TecUtilDataValueGetReadableRef(SurfZone, XVarNum);
    YFDPtr  = TecUtilDataValueGetReadableRef(SurfZone, YVarNum);
    ZFDPtr  = TecUtilDataValueGetReadableRef(SurfZone, ZVarNum);

    if (XFDPtr != NULL && YFDPtr != NULL && ZFDPtr != NULL)
    {
        double x1 = TecUtilDataValueGetByRef(XFDPtr, n1);
        double y1 = TecUtilDataValueGetByRef(YFDPtr, n1);
        double z1 = TecUtilDataValueGetByRef(ZFDPtr, n1);

        // Edge 1-2 components
        double dx12 = TecUtilDataValueGetByRef(XFDPtr, n2) - x1;
        double dy12 = TecUtilDataValueGetByRef(YFDPtr, n2) - y1;
        double dz12 = TecUtilDataValueGetByRef(ZFDPtr, n2) - z1;

        // Edge 1-3 components
        double dx13 = TecUtilDataValueGetByRef(XFDPtr, n3) - x1;
        double dy13 = TecUtilDataValueGetByRef(YFDPtr, n3) - y1;
        double dz13 = TecUtilDataValueGetByRef(ZFDPtr, n3) - z1;

        // Cross product
        double cpx  = dy12 * dz13 - dy13 * dz12;
        double cpy  = dz12 * dx13 - dz13 * dx12;
        double cpz  = dx12 * dy13 - dx13 * dy12;
        double cptot = sqrt(MAX((cpx * cpx + cpy * cpy + cpz * cpz), SMALLFLOAT));

        // Normal is the normalized cross-product
        Normal.X = cpx / cptot;
        Normal.Y = cpy / cptot;
        Normal.Z = cpz / cptot;
    }
    else
    {
        Normal.X = 0.0;
        Normal.Y = 0.0;
        Normal.Z = 0.0;
    }

    return Normal;
}








/**
 * Gets the components of the normal vector for the specified
 * corner of the specified element of the specified FEQuad zone.
 *
 * param SurfZone
 *     Zone number.
 * param ElemNum
 *     Element number.
 * param CornerNum
 *     Corner number (first, second, third, or fourth node of elem).
 *
 * return
 *     Normal vector.
 */
XYZ_s NormalsComputeQuadNorm(EntIndex_t SurfZone,
                             LgIndex_t  ElemNum,
                             LgIndex_t  CornerNum)
{
    XYZ_s         Normal;
    NodeMap_pa    NodeMap = NULL;
    FieldData_pa  XFDPtr = NULL;
    FieldData_pa  YFDPtr = NULL;
    FieldData_pa  ZFDPtr = NULL;
    EntIndex_t    XVarNum, YVarNum, ZVarNum;
    NodeMap_t     n1, n2, n3;
    LgIndex_t     CornerM1, CornerP1;

    REQUIRE(SurfZone > 0);
    REQUIRE(ElemNum > 0);
    REQUIRE(CornerNum > 0 && CornerNum <= 4);

    NodeMap = TecUtilDataNodeGetWritableRef(SurfZone);

    // Find nodes on either side of CornerNum
    if (NodeMap != NULL)
    {
        n1 = TecUtilDataNodeGetByRef(NodeMap, ElemNum, CornerNum);

        CornerP1 = CornerNum + 1;
        if (CornerP1 == 5) CornerP1 -= 4;
        n2 = TecUtilDataNodeGetByRef(NodeMap, ElemNum, CornerP1);

        if (n2 == n1)
        {
            CornerP1++;
            if (CornerP1 == 0) CornerP1 -= 4;
            n2 = TecUtilDataNodeGetByRef(NodeMap, ElemNum, CornerP1);
        }

        CornerM1 = CornerNum - 1;
        if (CornerM1 == 0) CornerM1 += 4;
        n3 = TecUtilDataNodeGetByRef(NodeMap, ElemNum, CornerM1);

        if (n3 == n1)
        {
            CornerM1--;
            if (CornerM1 == 0) CornerM1 += 4;
            n3 = TecUtilDataNodeGetByRef(NodeMap, ElemNum, CornerM1);
        }

        CHECK(n3 != n2);
    }

    // Compute the cross-product of two edge 1-2 and edge 1-3
    XVarNum = TecUtilVarGetNumByAssignment('X');
    YVarNum = TecUtilVarGetNumByAssignment('Y');
    ZVarNum = TecUtilVarGetNumByAssignment('Z');
    XFDPtr  = TecUtilDataValueGetReadableRef(SurfZone, XVarNum);
    YFDPtr  = TecUtilDataValueGetReadableRef(SurfZone, YVarNum);
    ZFDPtr  = TecUtilDataValueGetReadableRef(SurfZone, ZVarNum);

    if (XFDPtr != NULL && YFDPtr != NULL && ZFDPtr != NULL)
    {
        double x1 = TecUtilDataValueGetByRef(XFDPtr, n1);
        double y1 = TecUtilDataValueGetByRef(YFDPtr, n1);
        double z1 = TecUtilDataValueGetByRef(ZFDPtr, n1);

        // Edge 1-2 components
        double dx12 = TecUtilDataValueGetByRef(XFDPtr, n2) - x1;
        double dy12 = TecUtilDataValueGetByRef(YFDPtr, n2) - y1;
        double dz12 = TecUtilDataValueGetByRef(ZFDPtr, n2) - z1;

        // Edge 1-3 components
        double dx13 = TecUtilDataValueGetByRef(XFDPtr, n3) - x1;
        double dy13 = TecUtilDataValueGetByRef(YFDPtr, n3) - y1;
        double dz13 = TecUtilDataValueGetByRef(ZFDPtr, n3) - z1;

        // Cross product
        double cpx  = dy12 * dz13 - dy13 * dz12;
        double cpy  = dz12 * dx13 - dz13 * dx12;
        double cpz  = dx12 * dy13 - dx13 * dy12;
        double cptot = sqrt(MAX((cpx * cpx + cpy * cpy + cpz * cpz), SMALLFLOAT));

        // Normal is the normalized cross-product
        Normal.X = cpx / cptot;
        Normal.Y = cpy / cptot;
        Normal.Z = cpz / cptot;
    }
    else
    {
        Normal.X = 0.0;
        Normal.Y = 0.0;
        Normal.Z = 0.0;
    }

    return Normal;
}


/**
 * Compute the node-normals for a surface (FEQuad or FETriangle) zone.
 *
 * param Normals
 *     Normals structure containing array lists of normal vector components.
 * param SurfZone
 *     Number of the surface zone.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t NormalsCompute(Normals_pa Normals,
                         EntIndex_t SurfZone)
{
    EntIndex_t NumZones, NumVars;
    ZoneType_e ZoneType;
    Boolean_t  IsOk = TRUE;
    LgIndex_t  NumNodes, NumElems;

    REQUIRE(NormalsIsValid(Normals));

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    ZoneType = TecUtilZoneGetType(SurfZone);

    REQUIRE(SurfZone > 0 && SurfZone <= NumZones);
    REQUIRE(ZoneType == ZoneType_FETriangle || ZoneType == ZoneType_FEQuad);

    // Loop over elements, adding normal for ElemNum to Normals for contained nodes
    if (IsOk)
    {
        LgIndex_t    ElemNum, Offset;
        LgIndex_t    NodeNum;
        int          NodesPerElem = 3;
        Boolean_t    IsQuad = FALSE;
        XYZ_s        Normal;
        ArrListItem_u Item;
        ArrList_pa   NumElemsForNode = NULL;


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

        // Array list to track number of elements contributing to each node
        NumElemsForNode = ArrListAlloc(NumNodes + 1, ArrListType_UnsignedChar);
        for (Offset = 0; IsOk && Offset < NumNodes + 1; Offset++)
        {
            Item.UnsignedChar = 0;
            IsOk = ArrListSetItem(NumElemsForNode, Offset, Item);
        }

        // Loop over elements of SurfZone, adding element number for node numbers
        for (ElemNum = 1; IsOk && ElemNum <= NumElems; ElemNum++)
        {
            int       ii;

            // Compute triangle normal
            if (!IsQuad)
            {
                Normal = NormalsComputeTriNorm(SurfZone, ElemNum);
            }

            // Find unique node numbers
            for (ii = 0; IsOk && ii < NodesPerElem; ii++)  // ii is corner number
            {
                XYZ_s TmpNormal;
                int   NumContribElems;
                NodeNum = TecUtilDataNodeGetByZone(SurfZone, ElemNum, ii + 1);

                // Compute quad normal
                if (IsQuad)
                {
                    Normal = NormalsComputeQuadNorm(SurfZone, ElemNum, ii + 1);
                }

                Item = ArrListGetItem(NumElemsForNode, NodeNum - 1);
                NumContribElems = Item.UnsignedChar;

                // sum if other elements have already contributed to node, otherwise set
                if (NumContribElems > 0)
                {
                    TmpNormal = NormalsGetNormalForNode(Normals, NodeNum);
                    Normal.X = Normal.X + TmpNormal.X;
                    Normal.Y = Normal.Y + TmpNormal.Y;
                    Normal.Z = Normal.Z + TmpNormal.Z;
                }
                IsOk =  NormalsSetForNode(Normals, NodeNum, Normal);

                // Increment count of elements contributing to node
                Item.UnsignedChar = NumContribElems + 1;
                IsOk = ArrListSetItem(NumElemsForNode, NodeNum - 1, Item);
            }
        }

        // Divide by the number of elements contributing to the nodes
        for (NodeNum = 1; IsOk && NodeNum <= NumNodes; NodeNum++)
        {
            int    NumContribElems;

            Item = ArrListGetItem(NumElemsForNode, NodeNum - 1);
            NumContribElems = Item.UnsignedChar;

            if (NumContribElems > 0)
            {
                double Magnitude, RMagnitude;
                Normal = NormalsGetNormalForNode(Normals, NodeNum);
                Magnitude = sqrt(Normal.X * Normal.X + Normal.Y * Normal.Y + Normal.Z * Normal.Z);
                RMagnitude = 1.0 / (Magnitude + SMALLFLOAT);
                Normal.X *= RMagnitude;
                Normal.Y *= RMagnitude;
                Normal.Z *= RMagnitude;
            }
            IsOk =  NormalsSetForNode(Normals, NodeNum, Normal);
        }


        // Clean up
        ArrListDealloc(&NumElemsForNode);
    }

    ENSURE(NormalsIsValid(Normals));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}



/**
 * Compute the node-normals for a surface (FEQuad or FETriangle) zone using
 * a surface element map.
 *
 * param Normals
 *     Normals structure containing array lists of normal vector components.
 * param SurfZone
 *     Number of the surface zone.
 * param SurfElemMap
 *     Surface element map structure
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t NormalsComputeUsingSEM(Normals_pa     Normals,
                                 const EntIndex_t     SurfZone,
                                 const SurfElemMap_pa SurfElemMap)
{
    EntIndex_t NumZones, NumVars;
    ZoneType_e ZoneType;
    Boolean_t  IsOk = TRUE;
    LgIndex_t  NumNodes, NumElems;

    REQUIRE(NormalsIsValid(Normals));
    REQUIRE(SurfElemMapIsValid(SurfElemMap));

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    ZoneType = TecUtilZoneGetType(SurfZone);

    REQUIRE(SurfZone > 0 && SurfZone <= NumZones);
    REQUIRE(ZoneType == ZoneType_FETriangle || ZoneType == ZoneType_FEQuad);

    // Loop over nodes, adding normal for SurfElemMap ElemNum to Normals for the nodes
    if (IsOk)
    {
        LgIndex_t    ElemNum;
        LgIndex_t    NodeNum;
        int          NodesPerElem = 3;
        Boolean_t    IsQuad = FALSE;
        XYZ_s        Normal;


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

        for (NodeNum = 1; IsOk && NodeNum <= NumNodes; NodeNum++)
        {
            int ElemOffset, ii;
            LgIndex_t NumElemsForNode = SurfElemMapGetElemCountForNode(SurfElemMap, NodeNum);

            for (ElemOffset = 0; IsOk && ElemOffset < NumElemsForNode; ElemOffset++)
            {
                LgIndex_t NodeNumTry;
                XYZ_s TmpNormal;

                ElemNum = SurfElemMapGetElemNumForNodeOffset(SurfElemMap, NodeNum, ElemOffset);

                // Compute triangle normal
                if (!IsQuad)
                {
                    Normal = NormalsComputeTriNorm(SurfZone, ElemNum);
                }


                // Compute quad normal
                if (IsQuad)
                {
                    int ElemNodeNum;
                    Boolean_t IsFound = FALSE;

                    // Find element node numbers
                    for (ii = 0; !IsFound && ii < NodesPerElem; ii++)  // ii is corner number
                    {
                        NodeNumTry = TecUtilDataNodeGetByZone(SurfZone, ElemNum, ii + 1);
                        if (NodeNumTry = NodeNum)
                        {
                            IsFound = TRUE;
                            ElemNodeNum = ii + 1;
                        }
                    }

                    Normal = NormalsComputeQuadNorm(SurfZone, ElemNum, ElemNodeNum);
                }

                // sum if other elements have already contributed to node, otherwise set
                if (ElemOffset > 0)
                {
                    TmpNormal = NormalsGetNormalForNode(Normals, NodeNum);
                    Normal.X = Normal.X + TmpNormal.X;
                    Normal.Y = Normal.Y + TmpNormal.Y;
                    Normal.Z = Normal.Z + TmpNormal.Z;
                }
                IsOk =  NormalsSetForNode(Normals, NodeNum, Normal);

            }

        }


        // Normalize the vectors (unit length)
        for (NodeNum = 1; IsOk && NodeNum <= NumNodes; NodeNum++)
        {
            LgIndex_t NumElemsForNode = SurfElemMapGetElemCountForNode(SurfElemMap, NodeNum);

            if (NumElemsForNode > 0)
            {
                double Magnitude, RMagnitude;
                Normal = NormalsGetNormalForNode(Normals, NodeNum);
                Magnitude = sqrt(Normal.X * Normal.X + Normal.Y * Normal.Y + Normal.Z * Normal.Z);
                RMagnitude = 1.0 / (Magnitude + SMALLFLOAT);
                Normal.X *= RMagnitude;
                Normal.Y *= RMagnitude;
                Normal.Z *= RMagnitude;
            }
            IsOk =  NormalsSetForNode(Normals, NodeNum, Normal);
        }

    }



    ENSURE(NormalsIsValid(Normals));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}




/**
 * Compute the node-normals for a isosurface zone using the gradient of the
 * isosurface variable at the nodes. Assume that this gradient has been
 * precomputed for the volume zone and interpolated to the isosurface zone.
 *
 * param Normals
 *     Normals structure containing array lists of normal vector components.
 * param SurfZone
 *     Number of the surface zone.
 * param GradXVar, GradYVar, GradZVar
 *     Variables for components of the isosurface variable.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t NormalsComputeUsingIsoVarGrad(Normals_pa        Normals,
                                        const EntIndex_t  SurfZone,
                                        const EntIndex_t  GradXVar,
                                        const EntIndex_t  GradYVar,
                                        const EntIndex_t  GradZVar)
{
    EntIndex_t NumZones, NumVars;
    ZoneType_e ZoneType;
    Boolean_t  IsOk = TRUE;
    LgIndex_t  NumNodes, NumElems;
    FieldData_pa XGradFDPtr = NULL;
    FieldData_pa YGradFDPtr = NULL;
    FieldData_pa ZGradFDPtr = NULL;


    REQUIRE(NormalsIsValid(Normals));

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    ZoneType = TecUtilZoneGetType(SurfZone);

    REQUIRE(SurfZone > 0 && SurfZone <= NumZones);
    REQUIRE(ZoneType == ZoneType_FETriangle || ZoneType == ZoneType_FEQuad);
    REQUIRE(GradXVar > 0 && GradXVar <= NumVars);
    REQUIRE(GradYVar > 0 && GradYVar <= NumVars);
    REQUIRE(GradZVar > 0 && GradZVar <= NumVars);
    REQUIRE(GradXVar != GradYVar);
    REQUIRE(GradXVar != GradZVar);
    REQUIRE(GradYVar != GradZVar);

    // Get the field data pointers for the components of the gradient
    XGradFDPtr = TecUtilDataValueGetReadableRef(SurfZone, GradXVar);
    YGradFDPtr = TecUtilDataValueGetReadableRef(SurfZone, GradYVar);
    ZGradFDPtr = TecUtilDataValueGetReadableRef(SurfZone, GradZVar);
    if (XGradFDPtr == NULL || YGradFDPtr == NULL || ZGradFDPtr == NULL) IsOk = FALSE;

    // Loop over nodes, adding normal for SurfElemMap ElemNum to Normals for the nodes
    if (IsOk)
    {
        LgIndex_t    NodeNum;
        XYZ_s        Normal;


        // Find the number of Elements in the IsoSurfZone
        TecUtilZoneGetInfo(SurfZone, &NumNodes, &NumElems,
                           NULL, NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL);

        // Set the normal components to the negative of the gradient components and
        // nrmalize the normal vectors (unit length)
        for (NodeNum = 1; IsOk && NodeNum <= NumNodes; NodeNum++)
        {
            double MagnitudeSqr, RMagnitude;
            Normal.X = TecUtilDataValueGetByRef(XGradFDPtr, NodeNum);
            Normal.Y = TecUtilDataValueGetByRef(YGradFDPtr, NodeNum);
            Normal.Z = TecUtilDataValueGetByRef(ZGradFDPtr, NodeNum);
            MagnitudeSqr = Normal.X * Normal.X + Normal.Y * Normal.Y + Normal.Z * Normal.Z;

            if (MagnitudeSqr > SMALLFLOAT)
            {
                RMagnitude = 1.0 / sqrt(MagnitudeSqr);
                Normal.X *= RMagnitude;
                Normal.Y *= RMagnitude;
                Normal.Z *= RMagnitude;
            }
            else
            {
                IsOk = FALSE;
            }

            IsOk =  NormalsSetForNode(Normals, NodeNum, Normal);
        }

    }

    ENSURE(NormalsIsValid(Normals));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


double    NormalsCompare(Normals_pa   Normals1,
                         Normals_pa   Normals2)
{
    double MaxAngleDiff = 0.0;
    XYZ_s Normal1;
    XYZ_s Normal2;
    LgIndex_t Count1, Count2, CountMin, NodeOffset;

    REQUIRE(NormalsIsValid(Normals1));
    REQUIRE(NormalsIsValid(Normals2));

    Count1 = NormalsGetNodeCount(Normals1);
    Count2 = NormalsGetNodeCount(Normals2);
    CountMin = MIN(Count1, Count2);

    for (NodeOffset = 1; NodeOffset <= CountMin; NodeOffset++)
    {
        double DotProd, Normal1Mag, Normal2Mag, AngleDiff;

        Normal1 = NormalsGetNormalForNode(Normals1, NodeOffset);
        Normal2 = NormalsGetNormalForNode(Normals2, NodeOffset);

        DotProd = Normal1.X * Normal2.X + Normal1.Y * Normal2.Y + Normal1.Z * Normal2.Z;
        Normal1Mag = sqrt(Normal1.X * Normal1.X + Normal1.Y * Normal1.Y + Normal1.Z * Normal1.Z);
        Normal2Mag = sqrt(Normal2.X * Normal2.X + Normal2.Y * Normal2.Y + Normal2.Z * Normal2.Z);

        AngleDiff = acos(DotProd / (Normal1Mag * Normal2Mag));

        MaxAngleDiff = MAX(MaxAngleDiff, AngleDiff);
    }


    ENSURE(MaxAngleDiff >= 0.0);
    return MaxAngleDiff;
}



/*
 * Verify the unit test Normals for a regular cube FEQuad zone (cube with
 * six Quads and 8 nodes).
 *
 * return
 *     TRUE if successful, FALSE  if it fails.
 */
Boolean_t NormalsTestVerifyQuadBrick(Normals_pa Normals)
{
    Boolean_t      IsOk = TRUE;
    double         SqrRtOneThrd = sqrt(1.0 / 3.0);
    double         SmallDiff = 1.0e-6;
    XYZ_s          Normal;

    REQUIRE(NormalsIsValid(Normals));

    if (NormalsGetNodeCount(Normals) != 8) IsOk = FALSE;

    // Test the results: node 1 should have normal = (-sqrt(1/3), -sqrt(1/3), -sqrt(1/3))
    if (IsOk) Normal = NormalsGetNormalForNode(Normals, 1);
    if (IsOk) IsOk = (ABS(Normal.X + SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Y + SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Z + SqrRtOneThrd) < SmallDiff);

    // Test the results: node 2 should have normal = (sqrt(1/3), -sqrt(1/3), -sqrt(1/3))
    if (IsOk) Normal = NormalsGetNormalForNode(Normals, 2);
    if (IsOk) IsOk = (ABS(Normal.X - SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Y + SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Z + SqrRtOneThrd) < SmallDiff);

    // Test the results: node 3 should have normal = (sqrt(1/3), sqrt(1/3), -sqrt(1/3))
    if (IsOk) Normal = NormalsGetNormalForNode(Normals, 3);
    if (IsOk) IsOk = (ABS(Normal.X - SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Y - SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Z + SqrRtOneThrd) < SmallDiff);

    // Test the results: node 4 should have normal = (-sqrt(1/3),  sqrt(1/3), -sqrt(1/3))
    if (IsOk) Normal = NormalsGetNormalForNode(Normals, 4);
    if (IsOk) IsOk = (ABS(Normal.X + SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Y - SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Z + SqrRtOneThrd) < SmallDiff);

    // Test the results: node 5 should have normal = (-sqrt(1/3), -sqrt(1/3), sqrt(1/3))
    if (IsOk) Normal = NormalsGetNormalForNode(Normals, 5);
    if (IsOk) IsOk = (ABS(Normal.X + SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Y + SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Z - SqrRtOneThrd) < SmallDiff);

    // Test the results: node 6 should have normal = (sqrt(1/3), -sqrt(1/3), sqrt(1/3))
    if (IsOk) Normal = NormalsGetNormalForNode(Normals, 6);
    if (IsOk) IsOk = (ABS(Normal.X - SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Y + SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Z - SqrRtOneThrd) < SmallDiff);

    // Test the results: node 7 should have normal = (sqrt(1/3), sqrt(1/3), sqrt(1/3))
    if (IsOk) Normal = NormalsGetNormalForNode(Normals, 7);
    if (IsOk) IsOk = (ABS(Normal.X - SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Y - SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Z - SqrRtOneThrd) < SmallDiff);

    // Test the results: node 8 should have normal = (-sqrt(1/3), sqrt(1/3), sqrt(1/3))
    if (IsOk) Normal = NormalsGetNormalForNode(Normals, 8);
    if (IsOk) IsOk = (ABS(Normal.X + SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Y - SqrRtOneThrd) < SmallDiff);
    if (IsOk) IsOk = (ABS(Normal.Z - SqrRtOneThrd) < SmallDiff);

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/*
 * Unit test for Normals by creating an FEQuad zone (brick with
 * six Quads and 8 nodes) and see if it works correctly.
 *
 * return
 *     TRUE if successful, FALSE  if it fails.
 */
Boolean_t NormalsTest()
{
    Boolean_t      IsOk = TRUE;
    EntIndex_t     NewZone;
    Normals_pa     Normals = NULL;
    SurfElemMap_pa SurfElemMap = NULL;
    double         SqrRtOneThrd = sqrt(1.0 / 3.0);
    double         SmallDiff = 1.0e-6;
    // XYZ_s          Normal;
    Set_pa         ZoneList = NULL;

    // Create a temporary quad zone
    IsOk = TecUtilDataSetAddZone("Temporary FEQuad Zone", 8, 6, 1,
                                 ZoneType_FEQuad, NULL);
    if (IsOk)
    {
        NodeMap_pa NodeMap = NULL;
        Set_pa zones_added = TecUtilSetAlloc(TRUE);
        FieldData_pa XVarFD = NULL;
        FieldData_pa YVarFD = NULL;
        FieldData_pa ZVarFD = NULL;

        // new zone is always last zone
        TecUtilDataSetGetInfo(NULL, &NewZone, NULL);

        // fill NodeMap
        //      z                                 f4 is on back
        //      /|\n8 -----------------  n7
        //       |  /                / |
        //       |/   |    f6      /   |
        //     n5 ---------------- n6  |
        //       | f5 |  _ y      |    |
        //       |       /|  f2   | f3 |
        //       |  n4|/_  _  _ _ | _ _| n3
        //       |   /            |   /
        //       | /     f1       | /
        //     n1 ---------------- n2   ------->x
        //
        //
        NodeMap = TecUtilDataNodeGetWritableRef(NewZone);
        TecUtilDataNodeSetByRef(NodeMap, 1, 1, 1);
        TecUtilDataNodeSetByRef(NodeMap, 1, 2, 4);
        TecUtilDataNodeSetByRef(NodeMap, 1, 3, 3);
        TecUtilDataNodeSetByRef(NodeMap, 1, 4, 2);

        TecUtilDataNodeSetByRef(NodeMap, 2, 1, 1);
        TecUtilDataNodeSetByRef(NodeMap, 2, 2, 2);
        TecUtilDataNodeSetByRef(NodeMap, 2, 3, 6);
        TecUtilDataNodeSetByRef(NodeMap, 2, 4, 5);

        TecUtilDataNodeSetByRef(NodeMap, 3, 1, 2);
        TecUtilDataNodeSetByRef(NodeMap, 3, 2, 3);
        TecUtilDataNodeSetByRef(NodeMap, 3, 3, 7);
        TecUtilDataNodeSetByRef(NodeMap, 3, 4, 6);

        TecUtilDataNodeSetByRef(NodeMap, 4, 1, 3);
        TecUtilDataNodeSetByRef(NodeMap, 4, 2, 4);
        TecUtilDataNodeSetByRef(NodeMap, 4, 3, 8);
        TecUtilDataNodeSetByRef(NodeMap, 4, 4, 7);

        TecUtilDataNodeSetByRef(NodeMap, 5, 1, 1);
        TecUtilDataNodeSetByRef(NodeMap, 5, 2, 5);
        TecUtilDataNodeSetByRef(NodeMap, 5, 3, 8);
        TecUtilDataNodeSetByRef(NodeMap, 5, 4, 4);

        TecUtilDataNodeSetByRef(NodeMap, 6, 1, 5);
        TecUtilDataNodeSetByRef(NodeMap, 6, 2, 6);
        TecUtilDataNodeSetByRef(NodeMap, 6, 3, 7);
        TecUtilDataNodeSetByRef(NodeMap, 6, 4, 8);

        // Set the coordinates of the nodes (corners)
        TecUtilZoneGetInfo(NewZone, NULL, NULL, NULL, &XVarFD, &YVarFD, &ZVarFD, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL);

        if (XVarFD)
        {
            TecUtilDataValueSetByRef(XVarFD, 1, 0.0);
            TecUtilDataValueSetByRef(XVarFD, 2, 1.0);
            TecUtilDataValueSetByRef(XVarFD, 3, 1.0);
            TecUtilDataValueSetByRef(XVarFD, 4, 0.0);
            TecUtilDataValueSetByRef(XVarFD, 5, 0.0);
            TecUtilDataValueSetByRef(XVarFD, 6, 1.0);
            TecUtilDataValueSetByRef(XVarFD, 7, 1.0);
            TecUtilDataValueSetByRef(XVarFD, 8, 0.0);
        }

        if (YVarFD)
        {
            TecUtilDataValueSetByRef(YVarFD, 1, 0.0);
            TecUtilDataValueSetByRef(YVarFD, 2, 0.0);
            TecUtilDataValueSetByRef(YVarFD, 3, 1.0);
            TecUtilDataValueSetByRef(YVarFD, 4, 1.0);
            TecUtilDataValueSetByRef(YVarFD, 5, 0.0);
            TecUtilDataValueSetByRef(YVarFD, 6, 0.0);
            TecUtilDataValueSetByRef(YVarFD, 7, 1.0);
            TecUtilDataValueSetByRef(YVarFD, 8, 1.0);
        }

        if (ZVarFD)
        {
            TecUtilDataValueSetByRef(ZVarFD, 1, 0.0);
            TecUtilDataValueSetByRef(ZVarFD, 2, 0.0);
            TecUtilDataValueSetByRef(ZVarFD, 3, 0.0);
            TecUtilDataValueSetByRef(ZVarFD, 4, 0.0);
            TecUtilDataValueSetByRef(ZVarFD, 5, 1.0);
            TecUtilDataValueSetByRef(ZVarFD, 6, 1.0);
            TecUtilDataValueSetByRef(ZVarFD, 7, 1.0);
            TecUtilDataValueSetByRef(ZVarFD, 8, 1.0);
        }


        // inform Tecplot of new zone
        TecUtilSetAddMember(zones_added, NewZone, TRUE);
        TecUtilStateChanged(StateChange_ZonesAdded,
                            (ArbParam_t)zones_added);
        TecUtilSetDealloc(&zones_added);
    }


    // Compute the element map
    if (IsOk)
    {
        Normals = NormalsAlloc(8);
        IsOk = NormalsCompute(Normals, NewZone);
    }

    if (IsOk) IsOk = NormalsTestVerifyQuadBrick(Normals);

    // Delete the Normals
    NormalsDealloc(&Normals);


    // Try the second technique for computing the normal
    if (IsOk)
    {
        // Compute the element map
        SurfElemMap = SurfElemMapAlloc(8);
        IsOk = SurfElemMapCompute(SurfElemMap, NewZone);
    }

    if (IsOk)
    {
        Normals = NormalsAlloc(8);
        IsOk = NormalsComputeUsingSEM(Normals, NewZone, SurfElemMap);
    }

    if (IsOk) IsOk = NormalsTestVerifyQuadBrick(Normals);

    // Delete the Normals
    NormalsDealloc(&Normals);

    // Delete the temporary quad zone
    ZoneList = TecUtilSetAlloc(FALSE);
    if (ZoneList)
    {
        TecUtilSetAddMember(ZoneList, NewZone, FALSE);
        TecUtilDataSetDeleteZone(ZoneList);
        TecUtilSetDealloc(&ZoneList);
    }


    REQUIRE(VALID_BOOLEAN(IsOk));
    return IsOk;
}


