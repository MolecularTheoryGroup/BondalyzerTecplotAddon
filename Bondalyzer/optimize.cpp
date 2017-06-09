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
#include "OPTIMIZE.h"
#include <string.h>



/**
 * Determine if the OptPoint handle is sane.
 *
 * param OptPoint
 *     OptPoint structure in question.
 *
 * return
 *     TRUE if the OptPoint structure is valid, otherwise FALSE.
 */
Boolean_t OptPointIsValid(OptPoint_pa OptPoint)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(OptPoint) && VALID_REF(OptPoint->Coords));

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}






/**
 * Deallocates the OptPoint handle and set the handle to NULL.
 *
 * param
 *     Reference to a OptPoint handle.
 */
void OptPointDealloc(OptPoint_pa *OptPoint)
{
    REQUIRE(VALID_REF(OptPoint));
    REQUIRE(OptPointIsValid(*OptPoint) || *OptPoint == NULL);

    if (*OptPoint != NULL)
    {
        /* release the ArrList's */
        if ((*OptPoint)->Coords != NULL)
            FREE_ARRAY(((*OptPoint)->Coords), "OptPoint-Coords array");

        /* release the list structure itself */
        FREE_ITEM(*OptPoint, "OptPoint structure");
        *OptPoint = NULL;
    }

    ENSURE(*OptPoint == NULL);
}






/**
 * Allocates a OptPoint handle with requested capacity (NDim)
 * for the contained Coords array.
 *
 * param
 *     NDim: Dimension of the parameter space.
 *
 * return
 *     OptPoint handle if sufficient memory was available,
 *     otherwise a handle to NULL.
 */
OptPoint_pa OptPointAlloc(int NDim)
{
    OptPoint_pa Result = NULL;

    Result = ALLOC_ITEM(OptPoint_s, "OptPoint structure");
    if (Result != NULL)
        Result->Coords = ALLOC_ARRAY(NDim, double, "OptPoint->Coords array");

    ENSURE(OptPointIsValid(Result) || Result == NULL);
    return Result;
}








/**
 * Returns a (double *) pointer to the Coords array of the OptPoint structure.
 *
 * param
 *     OptPoint: OptPoint structure of interest.
 *
 * return
 *     Pointer to a double array.
 */
double*       OptPointGetCoordPointer(const OptPoint_pa OptPoint)
{
    double *Result = NULL;

    REQUIRE(OptPointIsValid(OptPoint));

    Result = OptPoint->Coords;

    ENSURE(VALID_REF(Result));
    return Result;
}






/**
 * Determine if the Optimize handle is sane.
 *
 * param Optimize
 *     Optimize structure in question.
 *
 * return
 *     TRUE if the Optimize structure is valid, otherwise FALSE.
 */
Boolean_t OptimizeIsValid(Optimize_pa Optimize)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(Optimize) &&
               VALID_REF(Optimize->OptPoints) && ArrListIsValid(Optimize->OptPoints));

    if (IsValid)
        IsValid = VALID_ENUM(Optimize->OptimizeMethod, OptimizeMethod_e);

    /* Require that each Optimize->OptPoints points to a valid OptPoint structure. */
    if (IsValid)
    {
        int ii;
        LgIndex_t Count = ArrListGetCount(Optimize->OptPoints);
        for (ii = 0; IsValid && ii < Count; ii++)
        {
            ArrListItem_u Item;
            OptPoint_pa OptPoint = NULL;

            Item = ArrListGetItem(Optimize->OptPoints, ii);
            OptPoint = (OptPoint_pa)Item.VoidPtr;
            IsValid = (VALID_REF(OptPoint) && OptPointIsValid(OptPoint));
        }
    }

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the Bundles handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a Optimize handle.
 */
void OptimizeDealloc(Optimize_pa *Optimize)
{
    REQUIRE(VALID_REF(Optimize));
    REQUIRE(OptimizeIsValid(*Optimize) || *Optimize == NULL);

    if (*Optimize != NULL)
    {
        /* Dealloc each OptPoint structured pointed to by the Optimize->OptPoints array list. */
        if ((*Optimize)->OptPoints != NULL)
        {
            int ii;
            LgIndex_t Count = ArrListGetCount((*Optimize)->OptPoints);
            for (ii = 0; ii < Count; ii++)
            {
                ArrListItem_u Item;
                OptPoint_pa OptPoint = NULL;

                Item = ArrListGetItem((*Optimize)->OptPoints, ii);
                OptPoint = (OptPoint_pa)Item.VoidPtr;
                if (OptPoint != NULL) OptPointDealloc(&OptPoint);
            }

            /* release the ArrList's */
            ArrListDealloc(&((*Optimize)->OptPoints));
        }

        /* release the list structure itself */
        FREE_ITEM(*Optimize, "Optimize structure");
        *Optimize = NULL;
    }

    ENSURE(*Optimize == NULL);
}





/**
 * Empties the Optimize structure.
 *
 * param Optimize
 *     Optimize to clear.
 */
void OptimizeClear(Optimize_pa Optimize)
{
    REQUIRE(OptimizeIsValid(Optimize));

    Optimize->OptimizeMethod = OptimizeMethod_NoOpt;

    /* Dealloc each OptPoint structured pointed to by the Optimize->OptPoints array list. */
    if (Optimize->OptPoints != NULL)
    {
        int ii;
        LgIndex_t Count = ArrListGetCount(Optimize->OptPoints);
        for (ii = 0; ii < Count; ii++)
        {
            ArrListItem_u Item;
            OptPoint_pa OptPoint = NULL;

            Item = ArrListGetItem(Optimize->OptPoints, ii);
            OptPoint = (OptPoint_pa)Item.VoidPtr;
            if (OptPoint != NULL) OptPointDealloc(&OptPoint);
        }

        /* Clear the OptPoints ArrList */
        ArrListClear(Optimize->OptPoints);
    }


    ENSURE(OptimizeIsValid(Optimize) && OptimizeGetPointCount(Optimize) == 0);
}




/**
 * Allocates a Optimize handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * return
 *     Optimize handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
Optimize_pa OptimizeAlloc(int NDim)
{
    Optimize_pa Result = NULL;

    REQUIRE(NDim > 0);

    Result = ALLOC_ITEM(Optimize_s, "Optimize structure");
    if (Result != NULL)
    {
        Result->OptimizeMethod = OptimizeMethod_DownhillSimplex;
        Result->NDim = NDim;

        Result->OptPoints   = ArrListAlloc(NDim + 1, ArrListType_VoidPtr);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->OptPoints == NULL)
        {
            FREE_ITEM(Result, "Optimize structure");
            Result = NULL;
        }
    }

    ENSURE(OptimizeIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of Points currently in the Optimize structure.
 *
 * param
 *     Optimize structure in question.
 *
 * return
 *     Number of points in the Optimize structure.
 */
LgIndex_t OptimizeGetPointCount(const Optimize_pa Optimize)
{
    LgIndex_t Result = 0;

    REQUIRE(OptimizeIsValid(Optimize));

    Result = ArrListGetCount(Optimize->OptPoints);

    ENSURE(Result >= 0);
    return Result;
}




/**
 * Gets the number of parameter dimensions currently in the Optimize structure.
 *
 * param
 *     Optimize structure in question.
 *
 * return
 *     Number of paramenter in the Optimize structure.
 */
LgIndex_t OptimizeGetDimensions(const Optimize_pa Optimize)
{
    LgIndex_t Result = 0;

    REQUIRE(OptimizeIsValid(Optimize));

    Result = Optimize->NDim;

    ENSURE(Result >= 0);
    return Result;
}





/**
 * Places OptPoint structure (point coordinates and variable value)
 * at the specified offset. If the offset is beyond the end of the
 * array lists, they are sized accordingly and the intervening
 * array values between the last node of the original state and
 * the last node of the new state are assigned NULL.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced.
 *
 * param Optimize
 *     Optimize target in which to set the coordinates and var value.
 * param PointOffset
 *     Offset of the point/node.
 * param OptPoint
 *     OptPoint structure point to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t OptimizeSetOptPointAtOffset(Optimize_pa Optimize,
                                      LgIndex_t   PointOffset,
                                      OptPoint_pa OptPoint)
{
    Boolean_t     IsOk = TRUE;
    ArrListItem_u Item;
    LgIndex_t     NDim;
    LgIndex_t     PointCount;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(OptPointIsValid(OptPoint));
    REQUIRE(PointOffset >= 0);

    PointCount = OptimizeGetPointCount(Optimize);
    NDim = OptimizeGetDimensions(Optimize);

    // If OptPoint already exists at PointOffset, Dealloc
    if (PointOffset < PointCount)
    {
        OptPoint_pa OldOptPoint = NULL;
        Item = ArrListGetItem(Optimize->OptPoints, PointOffset);
        OldOptPoint = (OptPoint_pa)Item.VoidPtr;
        OptPointDealloc(&OldOptPoint);
    }

    // Set pointer to OptPoint structure in OptPoints array list
    Item.VoidPtr = (void *)OptPoint;
    IsOk = ArrListSetItem(Optimize->OptPoints, PointOffset, Item);

    ENSURE(OptimizeIsValid(Optimize));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Places Optimize point coordinates and variable value at
 * the specified offset. If the offset is beyond the end of the
 * array lists, they are sized accordingly and the intervening
 * array values between the last node of the original state and
 * the last node of the new state are assigned 0.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced.
 *
 * param Optimize
 *     Optimize target in which to set the coordinates and var value.
 * param PointOffset
 *     Offset of the point/node.
 * param Coord1, Coord2, Var
 *     Coordinates and var value to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t OptimizeSetPointAtOffset(Optimize_pa Optimize,
                                   LgIndex_t   PointOffset,
                                   double     *Coords,
                                   double      Value)
{
    Boolean_t     IsOk = TRUE;
    ArrListItem_u Item;
    OptPoint_pa   OptPoint = NULL;
    double       *OptPointCoords = NULL;
    LgIndex_t     NDim;
    LgIndex_t     PointCount;
    int           ii;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(PointOffset >= 0);

    PointCount = OptimizeGetPointCount(Optimize);
    NDim = OptimizeGetDimensions(Optimize);

    if (PointOffset < PointCount)
    {
        Item = ArrListGetItem(Optimize->OptPoints, PointOffset);
        OptPoint = (OptPoint_pa)Item.VoidPtr;
    }
    else
    {
        OptPoint = OptPointAlloc(NDim);
        Item.VoidPtr = (void *)OptPoint;
        IsOk = ArrListSetItem(Optimize->OptPoints, PointOffset, Item);
    }

    OptPointCoords = OptPointGetCoordPointer(OptPoint);
    for (ii = 0; ii < NDim; ii++)
        OptPointCoords[ii] = Coords[ii];

    OptPoint->Value = Value;

    ENSURE(OptimizeIsValid(Optimize));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Gets the OptPoint structured from the Optimize structure for the
 * point at the specified offset.
 *
 * param Optimize
 *     Optimize structure containing the desired item.
 * param PointOffset
 *     Offset to the components of the Optimize.
 *
 * return
 *     OptPoint_pa point if it works, NULL otherwise.
 */

OptPoint_pa OptimizeGetPointFromOffset(const Optimize_pa  Optimize,
                                       const LgIndex_t    PointOffset)
{
    OptPoint_pa Result = NULL;
    ArrListItem_u Item;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(PointOffset >= 0 && PointOffset < OptimizeGetPointCount(Optimize));

    Item = ArrListGetItem(Optimize->OptPoints, PointOffset);
    Result = (OptPoint_pa)Item.VoidPtr;

    ENSURE(Result == NULL || VALID_REF(Result));

    return Result;
}






/**
 * Gets the Value for the point at the specified offset
 * from the Optimize structure.
 *
 * param Optimize
 *     Optimize structure containing the desired item.
 * param PointOffset
 *     Offset to the components of the Optimize.
 *
 * return
 *     double value.
 */

double OptimizeGetValueAtOffset(const Optimize_pa  Optimize,
                                const LgIndex_t    PointOffset)
{
    double Result;
    OptPoint_pa Point = NULL;
    ArrListItem_u Item;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(PointOffset >= 0 && PointOffset < OptimizeGetPointCount(Optimize));

    Item = ArrListGetItem(Optimize->OptPoints, PointOffset);
    Point = (OptPoint_pa)Item.VoidPtr;

    CHECK(OptPointIsValid(Point));

    Result = Point->Value;

    return Result;
}







/**
 * Gets the OptPoints in the Optimize structure at the two specified offsets.
 *
 * param Optimize
 *     Optimize structure containing the desired item.
 * param PointOffset1, PointOffset2
 *     Offsets to two point of the Optimize structure that are to be swapped.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t     OptimizeSwapPointsAtOffsets(Optimize_pa Optimize,
                                          LgIndex_t   PointOffset1,
                                          LgIndex_t   PointOffset2)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;
    OptPoint_pa Point1 = NULL;
    OptPoint_pa Point2 = NULL;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(PointOffset1 >= 0 && PointOffset1 < OptimizeGetPointCount(Optimize));
    REQUIRE(PointOffset2 >= 0 && PointOffset2 < OptimizeGetPointCount(Optimize));

    // Swap pointers
    Item = ArrListGetItem(Optimize->OptPoints, PointOffset1);
    Point1 = (OptPoint_pa)Item.VoidPtr;

    Item = ArrListGetItem(Optimize->OptPoints, PointOffset2);
    Point2 = (OptPoint_pa)Item.VoidPtr;

    if (Point1 == NULL || Point2 == NULL) IsOk = FALSE;

    if (IsOk)
    {
        Item.VoidPtr = (void *)Point2;
        IsOk = ArrListSetItem(Optimize->OptPoints, PointOffset1, Item);
    }
    if (IsOk)
    {
        Item.VoidPtr = (void *)Point1;
        IsOk = ArrListSetItem(Optimize->OptPoints, PointOffset2, Item);
    }

    ENSURE(OptimizeIsValid(Optimize));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Gets Optimize method (NotOpt, DownhillSimplex, etc.).
 *
 * param Optimize
 *     Optimize structure containing the desired item.
 *
 * return
 *     OptimizeMethod if the surface fit has been computed, OptimizeMethod_NoOpt otherwise.
 */

OptimizeMethod_e OptimizeGetMethod(const Optimize_pa  Optimize)
{
    OptimizeMethod_e Method = OptimizeMethod_NoOpt;

    REQUIRE(OptimizeIsValid(Optimize));

    Method = Optimize->OptimizeMethod;

    ENSURE(VALID_ENUM(Method, OptimizeMethod_e));
    return Method;
}






/**
 * Computes the average coordinates for the DownhillSimplex Optimize
 * Simplex, and stores them in the NDim+2 position of the Optimize->
 * OptPoints array list.
 *
 * param Optimize
 *     Optimize structure containing the simplex coordinates to average.
 *
 * returns
 *     TRUE if it worked, FALSE if there was an error.
 */

Boolean_t OptimizeAveSimplexPoints(Optimize_pa  Optimize)
{
    Boolean_t IsOk  = TRUE;
    LgIndex_t Count = OptimizeGetPointCount(Optimize);
    LgIndex_t NDim  = OptimizeGetDimensions(Optimize);
    OptPoint_pa  Ave = NULL;
    double      *AveCoords = NULL;
    int ii;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(OptimizeGetMethod(Optimize) == OptimizeMethod_DownhillSimplex);
    REQUIRE(Count >= NDim + 1);

    if (Count == NDim + 1)
    {
        Ave = OptPointAlloc(NDim);
    }
    else if (Count == NDim + 2)
    {
        Ave = OptimizeGetPointFromOffset(Optimize, NDim + 1);
        if (Ave == NULL) Ave = OptPointAlloc(NDim);
    }
    else
        IsOk = FALSE;

    // Get AveCoords array and initialize it and Ave->Value to zero
    if (IsOk)
    {
        int jj;
        AveCoords = OptPointGetCoordPointer(Ave);
        CHECK(VALID_REF(AveCoords));
        for (jj = 0; jj < NDim; jj++)
            AveCoords[jj] = 0.0;

        Ave->Value = 0.0;
    }

    // Sum the coordinates/Value of all NDim+1 points of the Simplex
    for (ii = 0; IsOk && ii < NDim + 1; ii++)
    {
        int jj;
        OptPoint_pa  Point = OptimizeGetPointFromOffset(Optimize, ii);
        double      *PointCoords = OptPointGetCoordPointer(Point);
        CHECK(VALID_REF(PointCoords));
        for (jj = 0; IsOk && jj < NDim; jj++)
            AveCoords[jj] += PointCoords[jj];

        Ave->Value += Point->Value;
    }

    // Divide by NDim+1 (number of points in simplex) to get averages
    if (IsOk)
    {
        for (ii = 0; ii < NDim; ii++)
            AveCoords[ii] /= (double)(NDim + 1);

        Ave->Value /= (double)(NDim + 1);

        CHECK(OptPointIsValid(Ave));
        IsOk = OptimizeSetOptPointAtOffset(Optimize, NDim + 2, Ave);
    }


    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}








/*
 * Extrapolate the Highest point and the average of the points on the opposite
 * "face" of the simplex. Accept the extrapolation if it makes the Highest value
 * smaller.
 *
 * Input:
 *   Optimize:    Structure containing data used in the optimization
 *                process. Note, this must contain at least one point
 *                that can be used as a starting point for the
 *                optimization
 *   Highest      Number, in the Optimize->OptPoints array list, of
 *                the point that currently has the largest Value.
 *   Factor       -1.0 == extrapolation of Highest through the opposite face
 *                        (mirror image)
 *                 2.0 == extrapolation (elongation) in the direction of Highest
 *                 0.5 == Contraction of Highest to a point 1/2 way between
 *                        Highest and opposite face.
 *   *ClientData: Pointer to structure containing various need info
 *                needed by Function.
 *   *Function    Pointer to a function that computes Value given Coords.
 * Return:
 *   New value of Highest..
 */
double OptimizeSimplexExtrap(Optimize_pa  Optimize,
                             LgIndex_t    Highest,
                             double       Factor,
                             void        *ClientData,
                             Boolean_t (*Function)(const void *, const double *, double *))
{
    double    NewValue = 0.0;
    Boolean_t IsOk = TRUE;
    LgIndex_t NDim;
    OptPoint_pa TryPoint = NULL;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(Highest >= 0 && Highest <= OptimizeGetDimensions(Optimize));
    REQUIRE(Factor > -2.1 && Factor < 2.1);  // Abritrary
    REQUIRE(VALID_REF(ClientData));
    REQUIRE(VALID_REF(Function));

    NDim = OptimizeGetDimensions(Optimize);
    TryPoint = OptPointAlloc(NDim);

    if (NDim <= 0 || TryPoint == NULL) IsOk = FALSE;

    CHECK(OptPointIsValid(TryPoint));

    // Extrapolate using the Highest point and the average of the opposite face
    if (IsOk)
    {
        OptPoint_pa SumOfPoints  = OptimizeGetPointFromOffset(Optimize, NDim + 2);
        OptPoint_pa HighestPoint = OptimizeGetPointFromOffset(Optimize, Highest);
        double *SumOfPointsCoords  = OptPointGetCoordPointer(SumOfPoints);
        double *HighestPointCoords = OptPointGetCoordPointer(HighestPoint);
        double *TryPointCoords     = OptPointGetCoordPointer(TryPoint);
        double  TryValue;

        // The following extrapolation formula is derived from
        //   P_face = ((NDim + 1) * P_sum - P_highest) / NDim
        double Fac1 = (NDim + 1) * (1.0 - Factor) / NDim;
        double Fac2 = Fac1 / (NDim + 1) - Factor;
        int jj;
        for (jj = 0; jj < NDim; jj++)
            TryPointCoords[jj] = Fac1 * SumOfPointsCoords[jj] - Fac2 * HighestPointCoords[jj];

        IsOk = Function(ClientData, TryPointCoords, &TryValue);
        if (IsOk)
        {
            TryPoint->Value = TryValue;

            if (TryValue < HighestPoint->Value)
            {
                NewValue = TryValue;
                // Adjust the SumOfPoints coordinates
                for (jj = 0; jj < NDim; jj++)
                    SumOfPointsCoords[jj] += (TryPointCoords[jj] - HighestPointCoords[jj]) / NDim;

                // Replace the HighestPoint with TryPoint
                OptPointDealloc(&HighestPoint);
                OptimizeSetOptPointAtOffset(Optimize, Highest, TryPoint);
            }
            else
                NewValue = HighestPoint->Value;
        }
    }

    ENSURE(VALID_BOOLEAN(IsOk) && IsOk == TRUE);
    return NewValue;
}









/*
 * Contract the edges toward the Lowest point.
 *
 * Input:
 *   Optimize:    Structure containing data used in the optimization
 *                process. Note, this must contain at least one point
 *                that can be used as a starting point for the
 *                optimization
 *   Lowest       Number, in the Optimize->OptPoints array list, of
 *                the point that currently has the lowest Value.
 *   Factor       Interpolation factor by which edges are contracted
 *                (0.5 == Average, 1.0 == No contraction)
 *   *ClientData: Pointer to structure containing various need info
 *                needed by Function.
 *   *Function    Pointer to a function that computes Value given Coords.
 * Return:
 *   New value of Highest..
 */
Boolean_t OptimizeSimplexContract(Optimize_pa  Optimize,
                                  LgIndex_t    Lowest,
                                  double       Factor,
                                  void        *ClientData,
                                  Boolean_t (*Function)(const void *, const double *, double *))
{
    Boolean_t IsOk = TRUE;
    LgIndex_t NDim;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(Lowest >= 0 && Lowest <= OptimizeGetDimensions(Optimize));
    REQUIRE(Factor > 0.0 && Factor < 2.1);  // Abritrary
    REQUIRE(VALID_REF(ClientData));
    REQUIRE(VALID_REF(Function));

    NDim = OptimizeGetDimensions(Optimize);

    if (NDim <= 0) IsOk = FALSE;

    // Contract all edges toward the Lowest point
    if (IsOk)
    {
        OptPoint_pa SumOfPoints   = OptimizeGetPointFromOffset(Optimize, NDim + 2);
        OptPoint_pa LowestPoint   = OptimizeGetPointFromOffset(Optimize, Lowest);
        double *SumOfPointsCoords = OptPointGetCoordPointer(SumOfPoints);
        double *LowestPointCoords = OptPointGetCoordPointer(LowestPoint);

        double Fac1 = 1.0 - Factor;
        int ii, jj;

        for (ii = 0; ii < NDim + 1; ii++)
        {
            if (ii != Lowest)
            {
                OptPoint_pa OtherEdgePoint = OptimizeGetPointFromOffset(Optimize, ii);
                double *OtherEdgePntCoords = OptPointGetCoordPointer(OtherEdgePoint);
                double OtherEdgePntValue;
                for (jj = 0; jj < NDim; jj++)
                {
                    double NewCoord = Factor * OtherEdgePntCoords[jj] + Fac1 * LowestPointCoords[jj];
                    double DeltaCoord = NewCoord - OtherEdgePntCoords[jj];
                    SumOfPointsCoords[jj] += DeltaCoord / (NDim + 1);
                    OtherEdgePntCoords[jj] = NewCoord;
                }
                IsOk = Function(ClientData, OtherEdgePntCoords, &OtherEdgePntValue);
                OtherEdgePoint->Value = OtherEdgePntValue;
            }
        }
    }


    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/*
 * Compute the point (Optimize->OptPoint[]) that minimizes Function
 * using the Downhill Simplex method.
 *
 * Input:
 *   Optimize:    Structure containing data used in the optimization
 *                process. Note, this must contain at least one point
 *                that can be used as a starting point for the
 *                optimization.
 *   *ClientData: Pointer to structure containing various need info
 *                needed by Function.
 *   *Function    Pointer to a function that computes Value given Coords.
 * Return:
 *   Index of the optimal point in Optimize if successful, -1 if failed.
 */
LgIndex_t OptimizeComputeDownhillSimplex(Optimize_pa  Optimize,
                                         void        *ClientData,
                                         Boolean_t (*Function)(const void *, const double *, double *))
{
    LgIndex_t OptPtIndex = -1;
    Boolean_t IsOk = TRUE;
    LgIndex_t NDim;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(OptimizeGetPointCount(Optimize) > 0);
    REQUIRE(VALID_ENUM(Optimize->OptimizeMethod, OptimizeMethod_e));
    REQUIRE(VALID_REF(ClientData));
    REQUIRE(VALID_REF(Function));

    NDim = OptimizeGetDimensions(Optimize);
    if (NDim < 2) IsOk = FALSE;


    // Compute other corners of the initial simplex using simple perturbation
    // in each of the Coords.
    if (IsOk)
    {
        int ii;
        double *Coords;
        double *Coords0;
        double DCoords = 0.1; // Temporary, need better method to estimate
        OptPoint_pa Point0 = OptimizeGetPointFromOffset(Optimize, 0);

        Coords0 = OptPointGetCoordPointer(Point0);
        Coords = ALLOC_ARRAY(NDim, double, "Temporary Coords array");

        for (ii = 0; IsOk && ii < NDim; ii++)
            Coords[ii] = Coords0[ii];


        for (ii = 1; IsOk && ii <= NDim; ii++)
        {
            OptPoint_pa Point = OptPointAlloc(NDim);
            double Value;

            Coords[ii-1] = Coords0[ii-1] + DCoords;

            IsOk = Function(ClientData, Coords, &Value);

            if (IsOk)
            {
                IsOk = OptimizeSetPointAtOffset(Optimize, ii, Coords, Value);

                Coords[ii-1] = Coords0[ii-1];
            }
        }

        FREE_ARRAY(Coords, "Temporary Coords array");
    }

    // Perform a series of reflections or contractions of the simplex to move
    // it toward a minimum
    if (IsOk)
    {
        Boolean_t Done = FALSE;
        LgIndex_t Highest = 0;
        LgIndex_t NextHighest = 0;
        LgIndex_t Lowest = 0;
        int iter;

        IsOk = OptimizeAveSimplexPoints(Optimize);

        for (iter = 0; !Done && IsOk && iter < Optimize->MaxIter; iter++)
        {
            int ii;
            double ConvergRatio;
            double ValHighest = OptimizeGetValueAtOffset(Optimize, 0);
            double ValNextHighest = ValHighest;
            double ValLowest = ValHighest;

            // First determine which points are the highest (worst), next
            // highest, and lowest (best)
            for (ii = 1; IsOk && ii <= NDim; ii++)
            {
                double ValTry = OptimizeGetValueAtOffset(Optimize, ii);
                if (ValTry <= ValLowest)
                {
                    Lowest = ii;
                    ValLowest = ValTry;
                }
                if (ValTry > ValHighest)
                {
                    Highest = ii;
                    ValHighest = ValTry;
                }
            }

            ValNextHighest = ValLowest - 0.1;
            for (ii = 1; IsOk && ii <= NDim; ii++)
            {
                if (ii != Highest && ii != Lowest)
                {
                    double ValTry = OptimizeGetValueAtOffset(Optimize, ii);
                    if (ValTry > ValNextHighest)
                    {
                        NextHighest = ii;
                        ValNextHighest = ValTry;
                    }
                }
            }

            // Check to see if minimization is adaquately converged
            ConvergRatio = 2.0 * ABS(ValHighest - ValLowest) / (ABS(ValHighest) + ABS(ValLowest) + SMALLFLOAT);
            if (ConvergRatio < Optimize->Tolerance)
            {
                IsOk = OptimizeSwapPointsAtOffsets(Optimize, 0, Lowest); // Put lowest in first spot
                Done = TRUE;
            }

            // Begin manipulation of the simplex
            if (!Done && IsOk)
            {
                double ValTry;

                // Extrapolate by a factor of -1 through the centroid - reflect from high point
                // If extrapolated point is better it replaces Highest point.
                ValTry = OptimizeSimplexExtrap(Optimize, Highest, -1, ClientData, Function);

                // If gives a better result than the previous lowest,
                // try an additional extrapolation - factor of 2
                if (ValTry <= ValLowest)
                    ValTry = OptimizeSimplexExtrap(Optimize, Highest, 2, ClientData, Function);

                // If the reflected point is worse than the previous next worst, try a
                // one-dimensional contraction (move highest point toward opposite face)
                else if (ValTry >= ValNextHighest)
                {
                    double ValSave = OptimizeGetValueAtOffset(Optimize, Highest);
                    ValTry = OptimizeSimplexExtrap(Optimize, Highest, 0.5, ClientData, Function);

                    // If it is still hasn't reduced the highest value, contract around the
                    // lowest point.
                    if (ValTry >= ValSave)
                    {
                        IsOk = OptimizeSimplexContract(Optimize, Lowest, 0.5, ClientData, Function);
                    }
                }
            }
        }
    }

    if (IsOk)
        OptPtIndex = 0;  // Above logic sets the optimial point to the zero location.
    else
        OptPtIndex = -1;

    ENSURE(VALID_BOOLEAN(IsOk));

    ENSURE(OptPtIndex >= -1 && OptPtIndex < OptimizeGetPointCount(Optimize));
    return OptPtIndex;
}





/*
 * Compute the point (Optimize->OptPoint[]) that minimizes Function.
 *
 * Input:
 *   Optimize:    Structure containing data used in the optimization
 *                process. Note, this must contain at least one point
 *                that can be used as a starting point for the
 *                optimization.
 *   OptimizeMethod: Method to use for optimization.
 *   *ClientData: Pointer to structure containing various need info
 *                needed by Function.
 *   *Function    Pointer to a function that computes Value given Coords.
 * Return:
 *   Index of the optimal point in Optimize if successful, -1 if failed.
 */
LgIndex_t OptimizeCompute(Optimize_pa      Optimize,
                          OptimizeMethod_e OptimizeMethod,
                          void        *ClientData,
                          Boolean_t (*Function)(const void *, const double *, double *))
{
    LgIndex_t OptPtIndex = -1;
    Boolean_t IsOk = TRUE;
    LgIndex_t NDim;

    REQUIRE(OptimizeIsValid(Optimize));
    REQUIRE(OptimizeGetPointCount(Optimize) > 0);
    REQUIRE(VALID_ENUM(OptimizeMethod, OptimizeMethod_e));
    REQUIRE(VALID_REF(ClientData));
    REQUIRE(VALID_REF(Function));

    NDim = OptimizeGetDimensions(Optimize);
    if (NDim < 1) IsOk = FALSE;

    switch (OptimizeMethod)
    {
        case OptimizeMethod_DownhillSimplex:
            Optimize->OptimizeMethod = OptimizeMethod;
            OptPtIndex = OptimizeComputeDownhillSimplex(Optimize, ClientData, Function);
            break;
        case OptimizeMethod_PowellsConjugateGrad:
        default:
            CHECK(FALSE);
            break;
    }

    if (IsOk == FALSE) OptPtIndex = -1;

    ENSURE(VALID_BOOLEAN(IsOk));

    ENSURE(OptPtIndex >= -1 && OptPtIndex < OptimizeGetPointCount(Optimize));
    return OptPtIndex;
}






/*
 * Test callback from OptimizeCompute used by the test case.
 * The test is a parabolic surface. The function computes
 * a value given the coordinates.
 *
 * Input:
 *   *ClientData: Pointer to structure containing various need info
 *                (in this case, )
 *   *Coords:     Independent variables of the multi-dimensional function.
 * Output:
 *   *Value:      Function value at the *Coords.
 *
 */
Boolean_t TestFunction(const void       *ClientData,
                       const double     *Coords,
                       double           *Value)
{
    Boolean_t  IsOk = TRUE;
    int NDims, ii;
    double *Origin;

    OptTestInfo_pa OptTestInfo = (OptTestInfo_pa)(ClientData);

    REQUIRE(VALID_REF(ClientData));
    REQUIRE(VALID_REF(Coords));
    REQUIRE(VALID_REF(Value));

    NDims  = OptTestInfo->NDims;
    Origin = OptTestInfo->Origin;

    // Compute the value of the quadratic function
    *Value = 0.0;
    for (ii = 0; ii < NDims; ii++)
    {
        double CoordsDiff = Coords[ii] - Origin[ii];
        *Value += CoordsDiff * CoordsDiff;
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}








/**
 * Compute the surface fit, of requested type, based on the input Coord1,
 * Coord2, and Var data stored in the Optimize structure. Use a singular
 * value decomposition to compute the coefficients of the surface fit
 * equations.
 *
 * param Optimize
 *     Optimize structure containing the input points for the curve fit.
 * param SurfFitType
 *     Type of surface fit (NoFit, Planar, Quadratic, etc.)..
 *
 * return
 *     SurfFitType if the surface fit has been computed, SurfFitType_NoFit otherwise.
 */
/*
Boolean_t OptimizeCompute(Optimize_pa Optimize,
                            SurfFitType_e SurfFitType)
{
  Boolean_t IsOk = TRUE;

  REQUIRE(OptimizeIsValid(Optimize));
  REQUIRE(VALID_ENUM(SurfFitType, SurfFitType_e));

  if (SurfFitType == SurfFitType_NoFit) IsOk = FALSE;

  // TODO: Add the surface fit
  if (IsOk)
    {
      int i;
      double wmax = 0.0;
      double wmin = 0.0;
      double *b, *x, *w, **a, **v;
      double r2sum = 0.0;
      double r2ave = 0.0;
      double r2min = LARGEFLOAT;

      switch(SurfFitType)
        {
          case SurfFitType_Quadratic:
            {
              // Allocate memory for temporary matrices and vectors
              LgIndex_t NumPoints = OptimizeGetPointCount(Optimize);

              a = SVDAllocMatrix(1, NumPoints, 1, 6);
              v = SVDAllocMatrix(1, 6, 1, 6);
              b = SVDAllocVector(1, NumPoints);
              w = SVDAllocVector(1, NumPoints);
              x = SVDAllocVector(1, NumPoints);

              if (a == NULL || b == NULL) IsOk = FALSE;
              if (!IsOk)
                break;

              // Build up the coefficient matrix and RHS vector
              for (i = 0; i < NumPoints; i++)
                {
                  double Coord1, Coord2, Var, dc1, dc2;
                  double r2;

                  IsOk = OptimizeGetPointFromOffset(Optimize, i, &Coord1, &Coord2, &Var);

                  if (!IsOk)
                    break;

                  dc1 = Coord1;  // TODO: make xn - x0;
                  dc2 = Coord2;  // TODO: make yn - y0;
                  a[i+1][1] = 1;
                  a[i+1][2] = dc1;
                  a[i+1][3] = dc2;
                  a[i+1][4] = 0.5*dc1*dc1;
                  a[i+1][5] = 0.5*dc2*dc2;
                  a[i+1][6] = dc1*dc2;

                  b[i+1]    = Var;

                  r2 = dc1*dc1 + dc2*dc2;
                  r2sum += r2;
                  if (i == 0)
                    r2min = r2;
                  else
                    r2min = MIN(r2, r2min);
                }

              if (!IsOk)
                break;

              //
              //  Scale the rows to by the weighting parameter 1/(r**2/r2ave + 0.1).
              //  Use the fact that dx and dy are stored in the a matrix.
              //
              r2ave = r2sum / (double)NumPoints;
              for (i = 0; i < NumPoints; i++)
                {
                  double weight, r2;
                  double dc1 = a[i+1][2];
                  double dc2 = a[i+1][3];

                  r2 = dc1*dc1 + dc2*dc2;
                  weight = r2ave / (r2 + 0.1*r2ave);

                  a[i+1][1] = weight*a[i+1][1];
                  a[i+1][2] = weight*a[i+1][2];
                  a[i+1][3] = weight*a[i+1][3];
                  a[i+1][4] = weight*a[i+1][4];
                  a[i+1][5] = weight*a[i+1][5];
                  a[i+1][6] = weight*a[i+1][6];

                  b[i+1]    = weight*b[i+1];
                }

              // Compute the singular value decomposition.
              ComputeSVD(a, NumPoints, 6, w, v);

              // If singular values are less than a certain level, set them to zero.
              for (i = 1; i <= 6; i++) if (wmax < w[i]) wmax = w[i];
              wmin = 1.0e-6 * wmax;
              for (i = 1; i <= 6; i++)
                {
                  if (w[i] < wmin)
                    {
                      w[i] = 0.0;
                    }
                }

              //
              // Solve the system of equations for the coefficients of
              // the curve fit.
              //
              BackSubstituteSVD(a, w, v, NumPoints, 6, b, x);

              // Store coefficients in the Optimize structure
              for (i = 0; IsOk && i < 6; i++)
                IsOk = OptimizeSetCoefAtOffset(Optimize, i, x[i+1]);

              // Free temporary matrices
              SVDFreeMatrix(a, 1, NumPoints, 1, 6);
              SVDFreeMatrix(v, 1, 6, 1, 6);
              SVDFreeVector(b, 1, NumPoints);
              SVDFreeVector(w, 1, NumPoints);
              SVDFreeVector(x, 1, NumPoints);

            }
            break;
          case SurfFitType_Planar:
            {
              // TODO: planar surface fit
              // Allocate memory for temporary matrices and vectors

              // Build up the coefficient matrix and RHS vector
            }
            break;
        }

      if (IsOk) Optimize->SurfFitType = SurfFitType;


    }

  ENSURE(VALID_ENUM(Optimize->SurfFitType, SurfFitType_e));
  ENSURE(VALID_BOOLEAN(IsOk));
  return IsOk;
}
*/





