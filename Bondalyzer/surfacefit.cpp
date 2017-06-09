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
#include "SURFACEFIT.h"
#include "SVD.h"
#include <string.h>



/**
 * Determine if the SurfaceFit handle is sane.
 *
 * param SurfaceFit
 *     SurfaceFit structure in question.
 *
 * return
 *     TRUE if the SurfaceFit structure is valid, otherwise FALSE.
 */
Boolean_t SurfaceFitIsValid(SurfaceFit_pa SurfaceFit)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(SurfaceFit) &&
               VALID_REF(SurfaceFit->Coord1) && ArrListIsValid(SurfaceFit->Coord1) &&
               VALID_REF(SurfaceFit->Coord2) && ArrListIsValid(SurfaceFit->Coord2) &&
               VALID_REF(SurfaceFit->Var) && ArrListIsValid(SurfaceFit->Var));

    if (IsValid)
        IsValid = VALID_ENUM(SurfaceFit->SurfFitType, SurfFitType_e);

    /* Require the same count for each array list in SurfaceFit structure. */
    if (IsValid)
    {
        LgIndex_t Count = ArrListGetCount(SurfaceFit->Coord1);
        IsValid = (ArrListGetCount(SurfaceFit->Coord2) == Count);
        if (IsValid) IsValid = (ArrListGetCount(SurfaceFit->Var) == Count);
    }

    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the SurfaceFit handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a SurfaceFit handle.
 */
void SurfaceFitDealloc(SurfaceFit_pa *SurfaceFit)
{
    REQUIRE(VALID_REF(SurfaceFit));
    REQUIRE(SurfaceFitIsValid(*SurfaceFit) || *SurfaceFit == NULL);

    if (*SurfaceFit != NULL)
    {
        /* release the ArrList's */
        if ((*SurfaceFit)->Coord1 != NULL) ArrListDealloc(&((*SurfaceFit)->Coord1));
        if ((*SurfaceFit)->Coord2 != NULL) ArrListDealloc(&((*SurfaceFit)->Coord2));
        if ((*SurfaceFit)->Var != NULL) ArrListDealloc(&((*SurfaceFit)->Var));

        /* release the list structure itself */
        FREE_ITEM(*SurfaceFit, "SurfaceFit structure");
        *SurfaceFit = NULL;
    }

    ENSURE(*SurfaceFit == NULL);
}





/**
 * Empties the SurfaceFit structure.
 *
 * param SurfaceFit
 *     SurfaceFit to clear.
 */
void SurfaceFitClear(SurfaceFit_pa SurfaceFit)
{
    REQUIRE(SurfaceFitIsValid(SurfaceFit));

    SurfaceFit->SurfFitType = SurfFitType_NoFit;

    ArrListClear(SurfaceFit->Coord1);
    ArrListClear(SurfaceFit->Coord2);
    ArrListClear(SurfaceFit->Var);
    ArrListClear(SurfaceFit->FitCoefs);

    ENSURE(SurfaceFitIsValid(SurfaceFit) && SurfaceFitGetPointCount(SurfaceFit) == 0);
}






/**
 * Allocates a SurfaceFit handle with a suitable default
 * capacity for the contained ArrList's.
 *
 * return
 *     SurfaceFit handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
SurfaceFit_pa SurfaceFitAlloc()
{
    SurfaceFit_pa Result = NULL;

    Result = ALLOC_ITEM(SurfaceFit_s, "SurfaceFit structure");
    if (Result != NULL)
    {
        Result->SurfFitType = SurfFitType_NoFit;

        Result->Coord1   = ArrListAlloc(20, ArrListType_Double);
        Result->Coord2   = ArrListAlloc(20, ArrListType_Double);
        Result->Var      = ArrListAlloc(20, ArrListType_Double);
        Result->FitCoefs = ArrListAlloc(20, ArrListType_Double);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->Coord1 == NULL || Result->Coord2 == NULL || Result->Var == NULL ||
            Result->FitCoefs == NULL)
        {
            if (Result->Coord1   != NULL) ArrListDealloc(&(Result->Coord1));
            if (Result->Coord2   != NULL) ArrListDealloc(&(Result->Coord2));
            if (Result->Var      != NULL) ArrListDealloc(&(Result->Var));
            if (Result->FitCoefs != NULL) ArrListDealloc(&(Result->FitCoefs));
            FREE_ITEM(Result, "SurfaceFit structure");
            Result = NULL;
        }
    }

    ENSURE(SurfaceFitIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of Points currently in the SurfaceFit structure.
 *
 * param
 *     SurfaceFit structure in question.
 *
 * return
 *     Number of points in the SurfaceFit structure.
 */
LgIndex_t SurfaceFitGetPointCount(const SurfaceFit_pa SurfaceFit)
{
    LgIndex_t Result = 0;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));

    Result = ArrListGetCount(SurfaceFit->Coord1);

    ENSURE(Result >= 0);
    return Result;
}






/**
 * Places SurfaceFit point coordinates and variable value at
 * the specified offset. If the offset is beyond the end of the
 * array lists, they are sized accordingly and the intervening
 * array values between the last node of the original state and
 * the last node of the new state are assigned 0.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced.
 *
 * param SurfaceFit
 *     SurfaceFit target in which to set the coordinates and var value.
 * param PointOffset
 *     Offset of the point/node.
 * param Coord1, Coord2, Var
 *     Coordinates and var value to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t SurfaceFitSetPointAtOffset(SurfaceFit_pa SurfaceFit,
                                     LgIndex_t  PointOffset,
                                     double     Coord1,
                                     double     Coord2,
                                     double     Var)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(SurfaceFitIsValid(SurfaceFit));
    REQUIRE(PointOffset >= 0);

    Item.Double = Coord1;
    IsOk = ArrListSetItem(SurfaceFit->Coord1, PointOffset, Item);

    if (IsOk)
    {
        Item.Double = Coord2;
        IsOk = ArrListSetItem(SurfaceFit->Coord2, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = Var;
        IsOk = ArrListSetItem(SurfaceFit->Var, PointOffset, Item);
    }

    ENSURE(SurfaceFitIsValid(SurfaceFit));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Inserts SurfaceFit point coordinates and variable value at the specified
 * offset. The arrays will be expanded to accomodate the additional value.
 *
 *
 * param SurfaceFit
 *     SurfaceFit target in which to set the coordinates.
 * param PointOffset
 *     Offset at which to insert the point/node.
 * param Coord1, Coord2, Var
 *     Coordinates and var value to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t SurfaceFitInsertPointAtOffset(SurfaceFit_pa SurfaceFit,
                                        LgIndex_t  PointOffset,
                                        double     Coord1,
                                        double     Coord2,
                                        double     Var)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(SurfaceFitIsValid(SurfaceFit));
    REQUIRE(0 <= PointOffset && PointOffset <= SurfaceFitGetPointCount(SurfaceFit));

    Item.Double = Coord1;
    IsOk = ArrListInsertItem(SurfaceFit->Coord1, PointOffset, Item);

    if (IsOk)
    {
        Item.Double = Coord2;
        IsOk = ArrListInsertItem(SurfaceFit->Coord2, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = Var;
        IsOk = ArrListInsertItem(SurfaceFit->Var, PointOffset, Item);
    }

    ENSURE(SurfaceFitIsValid(SurfaceFit));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Appends the point components to SurfaceFit structure. The array lists
 * will be expanded to accommodate the additional items.
 *
 * param SurfaceFit
 *     SurfaceFit target to which the point components are to be appended.
 * param Coord1, Coord2, Var
 *     Coordinates and var value to append.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t SurfaceFitAppendPointAtEnd(SurfaceFit_pa  SurfaceFit,
                                     double     Coord1,
                                     double     Coord2,
                                     double     Var)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));

    Count = SurfaceFitGetPointCount(SurfaceFit);

    IsOk = SurfaceFitInsertPointAtOffset(SurfaceFit, Count, Coord1, Coord2, Var);

    ENSURE(SurfaceFitIsValid(SurfaceFit));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Gets the components of the SurfaceFit structure for the
 * point at the specified offset.
 *
 * param SurfaceFit
 *     SurfaceFit structure containing the desired item.
 * param PointOffset
 *     Offset to the components of the SurfaceFit.
 * param *Coord1, *Coord2, *Var
 *     Pointers to components of the point
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */

Boolean_t SurfaceFitGetPointFromOffset(const SurfaceFit_pa  SurfaceFit,
                                       const LgIndex_t   PointOffset,
                                       double     *Coord1,
                                       double     *Coord2,
                                       double     *Var)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));
    REQUIRE(PointOffset >= 0 && PointOffset < SurfaceFitGetPointCount(SurfaceFit));
    REQUIRE(VALID_REF(Coord1));
    REQUIRE(VALID_REF(Coord2));
    REQUIRE(VALID_REF(Var));

    Item = ArrListGetItem(SurfaceFit->Coord1, PointOffset);
    *Coord1 = Item.Double;

    Item = ArrListGetItem(SurfaceFit->Coord2, PointOffset);
    *Coord2 = Item.Double;

    Item = ArrListGetItem(SurfaceFit->Var, PointOffset);
    *Var = Item.Double;

    ENSURE(VALID_REF(Coord1));
    ENSURE(VALID_REF(Coord2));
    ENSURE(VALID_REF(Var));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}








/**
 * Gets SurfaceFit type (NotFit, Planar, Quadratic, etc.).
 *
 * param SurfaceFit
 *     SurfaceFit structure containing the desired item.
 *
 * return
 *     SurfFitType if the surface fit has been computed, SurfFitType_NoFit otherwise.
 */

SurfFitType_e SurfaceFitGetType(const SurfaceFit_pa  SurfaceFit)
{
    SurfFitType_e Type = SurfFitType_NoFit;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));

    Type = SurfaceFit->SurfFitType;

    ENSURE(VALID_ENUM(Type, SurfFitType_e));
    return Type;
}




/**
 * Gets the number of coefficients currently in the SurfaceFit structure.
 *
 * param
 *     SurfaceFit structure in question.
 *
 * return
 *     Number of coefficients in the SurfaceFit structure.
 */
LgIndex_t SurfaceFitGetCoefCount(const SurfaceFit_pa SurfaceFit)
{
    LgIndex_t Result = 0;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));

    Result = ArrListGetCount(SurfaceFit->FitCoefs);

    ENSURE(Result >= 0);
    return Result;
}







/**
 * Places SurfaceFit coefficients at
 * the specified offset. If the offset is beyond the end of the
 * array lists, they are sized accordingly and the intervening
 * array values between the last node of the original state and
 * the last node of the new state are assigned 0.
 *
 * note
 *     If coefficients already exists at the specified location they are
 *     replaced.
 *
 * param SurfaceFit
 *     SurfaceFit target in which to set the coordinates and var value.
 * param CoefOffset
 *     Offset of the coefficient.
 * param Coef
 *     Coordinates and var value to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t SurfaceFitSetCoefAtOffset(SurfaceFit_pa SurfaceFit,
                                    LgIndex_t  CoefOffset,
                                    double     Coef)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(SurfaceFitIsValid(SurfaceFit));
    REQUIRE(CoefOffset >= 0);

    Item.Double = Coef;
    IsOk = ArrListSetItem(SurfaceFit->FitCoefs, CoefOffset, Item);


    ENSURE(SurfaceFitIsValid(SurfaceFit));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Inserts SurfaceFit coefficient at the specified
 * offset. The arrays will be expanded to accomodate the additional value.
 *
 *
 * param SurfaceFit
 *     SurfaceFit target in which to set the coordinates.
 * param CoefOffset
 *     Offset at which to insert the coefficient.
 * param Coef
 *     Coefficient to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t SurfaceFitInsertCoefAtOffset(SurfaceFit_pa SurfaceFit,
                                       LgIndex_t  CoefOffset,
                                       double     Coef)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(SurfaceFitIsValid(SurfaceFit));
    REQUIRE(0 <= CoefOffset && CoefOffset <= SurfaceFitGetCoefCount(SurfaceFit));

    Item.Double = Coef;
    IsOk = ArrListInsertItem(SurfaceFit->FitCoefs, CoefOffset, Item);


    ENSURE(SurfaceFitIsValid(SurfaceFit));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Appends the coefficient to SurfaceFit structure. The array lists
 * will be expanded to accommodate the additional items.
 *
 * param SurfaceFit
 *     SurfaceFit target to which the point components are to be appended.
 * param Coef
 *     Coefficient to append.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t SurfaceFitAppendCoefAtEnd(SurfaceFit_pa  SurfaceFit,
                                    double     Coef)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));

    Count = SurfaceFitGetCoefCount(SurfaceFit);

    IsOk = SurfaceFitInsertCoefAtOffset(SurfaceFit, Count, Coef);

    ENSURE(SurfaceFitIsValid(SurfaceFit));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Gets the coefficient of the SurfaceFit structure for the
 * point at the specified offset.
 *
 * param SurfaceFit
 *     SurfaceFit structure containing the desired item.
 * param CoefOffset
 *     Offset to the coefficient of the SurfaceFit.
 * param *Coef
 *     Pointers to surface fit coefficient
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */

double SurfaceFitGetCoefFromOffset(const SurfaceFit_pa  SurfaceFit,
                                   const LgIndex_t   CoefOffset)
{
    Boolean_t IsOk = TRUE;
    double    Coef;
    ArrListItem_u Item;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));
    REQUIRE(CoefOffset >= 0 && CoefOffset < SurfaceFitGetCoefCount(SurfaceFit));

    Item = ArrListGetItem(SurfaceFit->FitCoefs, CoefOffset);
    Coef = Item.Double;

    ENSURE(VALID_BOOLEAN(IsOk));
    return Coef;
}






/**
 * Compute the surface fit, of requested type, based on the input Coord1,
 * Coord2, and Var data stored in the SurfaceFit structure. Use a singular
 * value decomposition to compute the coefficients of the surface fit
 * equations.
 *
 * param SurfaceFit
 *     SurfaceFit structure containing the input points for the curve fit.
 * param SurfFitType
 *     Type of surface fit (NoFit, Planar, Quadratic, etc.)..
 *
 * return
 *     SurfFitType if the surface fit has been computed, SurfFitType_NoFit otherwise.
 */
Boolean_t SurfaceFitCompute(SurfaceFit_pa SurfaceFit,
                            SurfFitType_e SurfFitType)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));
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

        switch (SurfFitType)
        {
            case SurfFitType_Quadratic:
            {
                // Allocate memory for temporary matrices and vectors
                LgIndex_t NumPoints = SurfaceFitGetPointCount(SurfaceFit);

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

                    IsOk = SurfaceFitGetPointFromOffset(SurfaceFit, i, &Coord1, &Coord2, &Var);

                    if (!IsOk)
                        break;

                    dc1 = Coord1;  // TODO: make xn - x0;
                    dc2 = Coord2;  // TODO: make yn - y0;
                    a[i+1][1] = 1;
                    a[i+1][2] = dc1;
                    a[i+1][3] = dc2;
                    a[i+1][4] = dc1 * dc1;
                    a[i+1][5] = dc1 * dc2;
                    a[i+1][6] = dc2 * dc2;

                    b[i+1]    = Var;

                    r2 = dc1 * dc1 + dc2 * dc2;
                    r2sum += r2;
                    if (i == 0)
                        r2min = r2;
                    else
                        r2min = MIN(r2, r2min);
                }

                if (!IsOk)
                    break;

                /*
                 *  Scale the rows to by the weighting parameter 1/(r**2/r2ave + 0.1).
                 *  Use the fact that dx and dy are stored in the a matrix.
                 */
                r2ave = r2sum / (double)NumPoints;
                for (i = 0; i < NumPoints; i++)
                {
                    double weight, r2;
                    double dc1 = a[i+1][2];
                    double dc2 = a[i+1][3];

                    r2 = dc1 * dc1 + dc2 * dc2;
                    weight = r2ave / (r2 + 0.1 * r2ave);

                    a[i+1][1] = weight * a[i+1][1];
                    a[i+1][2] = weight * a[i+1][2];
                    a[i+1][3] = weight * a[i+1][3];
                    a[i+1][4] = weight * a[i+1][4];
                    a[i+1][5] = weight * a[i+1][5];
                    a[i+1][6] = weight * a[i+1][6];

                    b[i+1]    = weight * b[i+1];
                }

                /* Compute the singular value decomposition. */
                ComputeSVD(a, NumPoints, 6, w, v);

                /* If singular values are less than a certain level, set them to zero. */
                for (i = 1; i <= 6; i++) if (wmax < w[i]) wmax = w[i];
                wmin = 1.0e-6 * wmax;
                for (i = 1; i <= 6; i++)
                {
                    if (w[i] < wmin)
                    {
                        w[i] = 0.0;
                    }
                }

                /*
                 * Solve the system of equations for the coefficients of
                 * the curve fit.
                 */
                BackSubstituteSVD(a, w, v, NumPoints, 6, b, x);

                // Store coefficients in the SurfaceFit structure
                for (i = 0; IsOk && i < 6; i++)
                    IsOk = SurfaceFitSetCoefAtOffset(SurfaceFit, i, x[i+1]);

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
                // Allocate memory for temporary matrices and vectors
                LgIndex_t NumPoints = SurfaceFitGetPointCount(SurfaceFit);

                a = SVDAllocMatrix(1, NumPoints, 1, 3);
                v = SVDAllocMatrix(1, 3, 1, 3);
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

                    IsOk = SurfaceFitGetPointFromOffset(SurfaceFit, i, &Coord1, &Coord2, &Var);

                    if (!IsOk)
                        break;

                    dc1 = Coord1;  // TODO: make xn - x0;
                    dc2 = Coord2;  // TODO: make yn - y0;
                    a[i+1][1] = 1;
                    a[i+1][2] = dc1;
                    a[i+1][3] = dc2;

                    b[i+1]    = Var;

                    r2 = dc1 * dc1 + dc2 * dc2;
                    r2sum += r2;
                    if (i == 0)
                        r2min = r2;
                    else
                        r2min = MIN(r2, r2min);
                }

                if (!IsOk)
                    break;

                /*
                 *  Scale the rows to by the weighting parameter 1/(r**2/r2ave + 0.1).
                 *  Use the fact that dx and dy are stored in the a matrix.
                 */
                r2ave = r2sum / (double)NumPoints;
                for (i = 0; i < NumPoints; i++)
                {
                    double weight, r2;
                    double dc1 = a[i+1][2];
                    double dc2 = a[i+1][3];

                    r2 = dc1 * dc1 + dc2 * dc2;
                    weight = r2ave / (r2 + 0.1 * r2ave);

                    a[i+1][1] = weight * a[i+1][1];
                    a[i+1][2] = weight * a[i+1][2];
                    a[i+1][3] = weight * a[i+1][3];

                    b[i+1]    = weight * b[i+1];
                }

                /* Compute the singular value decomposition. */
                ComputeSVD(a, NumPoints, 3, w, v);

                /* If singular values are less than a certain level, set them to zero. */
                for (i = 1; i <= 3; i++) if (wmax < w[i]) wmax = w[i];
                wmin = 1.0e-6 * wmax;
                for (i = 1; i <= 3; i++)
                {
                    if (w[i] < wmin)
                    {
                        w[i] = 0.0;
                    }
                }

                /*
                 * Solve the system of equations for the coefficients of
                 * the curve fit.
                 */
                BackSubstituteSVD(a, w, v, NumPoints, 3, b, x);

                // Store coefficients in the SurfaceFit structure
                for (i = 0; IsOk && i < 3; i++)
                    IsOk = SurfaceFitSetCoefAtOffset(SurfaceFit, i, x[i+1]);

                // Free temporary matrices
                SVDFreeMatrix(a, 1, NumPoints, 1, 3);
                SVDFreeMatrix(v, 1, 3, 1, 3);
                SVDFreeVector(b, 1, NumPoints);
                SVDFreeVector(w, 1, NumPoints);
                SVDFreeVector(x, 1, NumPoints);
            }
            break;
        }

        if (IsOk) SurfaceFit->SurfFitType = SurfFitType;


    }

    ENSURE(VALID_ENUM(SurfaceFit->SurfFitType, SurfFitType_e));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * For an existing SurfaceFit and independent coordinates, compute the
 * interpolated value of Var.
 *
 * param SurfaceFit
 *     SurfaceFit structure containing the input points for the curve fit.
 * param SurfFitType
 *     Type of surface fit (NoFit, Planar, Quadratic, etc.)..
 *
 * return
 *     Var at Coord1,Coord2 if everything is OK, -LARGEDOUBLE otherwise.
 */
double    SurfaceFitInterpolateVar(SurfaceFit_pa SurfaceFit,
                                   double        Coord1,
                                   double        Coord2)
{
    double Result = -LARGEDOUBLE;
    Boolean_t IsOk = TRUE;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));
    ENSURE(VALID_ENUM(SurfaceFit->SurfFitType, SurfFitType_e));

    switch (SurfaceFit->SurfFitType)
    {
        case SurfFitType_Quadratic:
        {
            double a, b, c, d, e, f;

            CHECK(SurfaceFitGetCoefCount(SurfaceFit) == 6);

            if (SurfaceFitGetCoefCount(SurfaceFit) == 6)
            {
                a = SurfaceFitGetCoefFromOffset(SurfaceFit, 0);
                b = SurfaceFitGetCoefFromOffset(SurfaceFit, 1);
                c = SurfaceFitGetCoefFromOffset(SurfaceFit, 2);
                d = SurfaceFitGetCoefFromOffset(SurfaceFit, 3);
                e = SurfaceFitGetCoefFromOffset(SurfaceFit, 4);
                f = SurfaceFitGetCoefFromOffset(SurfaceFit, 5);
                Result = a + (b + d * Coord1 + e * Coord2) * Coord1
                         + (c              + f * Coord2) * Coord2;
            }
        }
        break;
        case SurfFitType_Planar:
        {
            double a, b, c;

            CHECK(SurfaceFitGetCoefCount(SurfaceFit) == 3);

            if (SurfaceFitGetCoefCount(SurfaceFit) == 3)
            {
                a = SurfaceFitGetCoefFromOffset(SurfaceFit, 0);
                b = SurfaceFitGetCoefFromOffset(SurfaceFit, 1);
                c = SurfaceFitGetCoefFromOffset(SurfaceFit, 2);
                Result = a + b * Coord1 + c * Coord2;
            }
        }
        break;
        case SurfFitType_NoFit:
            IsOk = FALSE;
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return Result;
}



/**
 * For an existing SurfaceFit, compute the location of the critical point
 * (level point in Var) and the eigenvalues of the curvature matrix.
 *
 * param SurfaceFit
 *     SurfaceFit structure containing the input points for the curve fit.
 *
 * return
 *     TRUE if it worked (valid non-planar SurfaceFit). FALSE otherwise.
 */
Boolean_t SurfaceFitCompCritPoint(SurfaceFit_pa SurfaceFit)
{
    Boolean_t IsOk = TRUE;

    REQUIRE(SurfaceFitIsValid(SurfaceFit));

    switch (SurfaceFit->SurfFitType)
    {
        case SurfFitType_Quadratic:
        {
            if (SurfaceFitGetCoefCount(SurfaceFit) == 6)
            {
                double a = SurfaceFitGetCoefFromOffset(SurfaceFit, 0);
                double b = SurfaceFitGetCoefFromOffset(SurfaceFit, 1);
                double c = SurfaceFitGetCoefFromOffset(SurfaceFit, 2);
                double d = SurfaceFitGetCoefFromOffset(SurfaceFit, 3);
                double e = SurfaceFitGetCoefFromOffset(SurfaceFit, 4);
                double f = SurfaceFitGetCoefFromOffset(SurfaceFit, 5);
                //
                // Quadratic surface fit:
                //      v = a + b x + c y + d x^2 + e x y + f y^2
                //
                // Solve the following 2x2 system for dv/dx = dv/dy = 0.
                //  _        _   _  _       _  _
                // | 2d     e | | c1 |     | -b |
                // |          | |    |  =  |    |
                // |_ e    2f_| |_c2_|     |_-c_|
                //
                double det = 4 * d * f - e * e;
                if (ABS(det) > 1.0e-10)
                {
                    SurfaceFit->c1_crt = (-2 * f * b + e * c) / det;
                    SurfaceFit->c2_crt = (e * b - 2 * d * c) / det;
                }
                else  // System singular, no single critical point exists
                    IsOk = FALSE;

                //
                // Find the eigenvalues and eigenvectors of the curvature matrix
                //   _        _
                //  | 2d     e |
                //  |          |
                //  |_ e    2f_|
                //
                // Use a single Jacobi rotation to zero out the off-diagonal elements
                //   _       _     _        _   _        _   _        _
                //  | l1   0  |   |  c    s  | | 2d    e  | |  c   -s  |
                //  |         | = |          | |          | |          |
                //  |_ 0  l2 _|   |_-s    c _| |_ e   2f _| |_ s    c _|
                //
                // In this equation, l1 and l2 are the eigenvalues and the columns of
                // the rightmost transformation matrix are the eigenvectors.
                //
                if (IsOk)
                {
                    double Cos = 1.0;
                    double Sin = 0.0;
                    double Tan = 0.0;

                    // Compute the rotation angle for the eigenvectors (default is zero)
                    if (ABS(e) > 0.000001)
                    {
                        double Theta = (d - f) / e;
                        double SgnThta;

                        SgnThta = 1.0;
                        if (Theta < 0.0) SgnThta = -1;

                        Tan = SgnThta / (ABS(Theta) + sqrt(Theta * Theta + 1.0));
                        Cos = 1.0 / sqrt(Tan * Tan + 1.0);
                        Sin = Tan * Cos;
                        // Save in structure somehow
                        SurfaceFit->RotationAngle = asin(Sin);
                    }
                    else
                        SurfaceFit->RotationAngle = 0.0;


                    // Compute eigenvalues
                    if (IsOk)
                    {
                        SurfaceFit->lambda1 = 2.0 * d + Tan * e;
                        SurfaceFit->lambda2 = 2.0 * f - Tan * e;
                    }

                }

            }
            else
                IsOk = FALSE;  // Incorrect Num Coefs for quadratic surface fit.
        }
        break;
        case SurfFitType_Planar:
        {
            IsOk = FALSE;  //Doesn't work with planar surfaces
        }
        break;
        case SurfFitType_NoFit:
            IsOk = FALSE;
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}




/**
 * Test the SurfaceFit module. Create a set of points using a known quadratic
 * function and see if the SurfaceFit module accurately duplicates the
 * quadratic surface.
 *
 * TODO: do the same thing for linear surface.
 *
 *
 * return
 *     TRUE if it worked (valid non-planar SurfaceFit). FALSE otherwise.
 */
Boolean_t SurfaceFitTest()
{
    Boolean_t IsOk = TRUE;
    SurfaceFit_pa TestSurfaceFit = NULL;

    // Allocate the surface fit structure
    TestSurfaceFit = SurfaceFitAlloc();
    if (TestSurfaceFit == NULL) IsOk = FALSE;

    CHECK(SurfaceFitIsValid(TestSurfaceFit));

    //  Simple Test:
    //
    //                   ^ y = Coord2
    //                   |
    //                   |
    //                   |
    //         ______________________ (1, 1)
    //        |          |           |
    //        |          |           |
    //        |          |           |
    //        |          |           |
    //        |----------------------|  ---------> x = Coord1
    //        |          |           |
    //        |          |           |
    //        |          |           |
    //        |__________|___________|
    //      (-1,-1)
    //
    //  Simple Test 1:  p = x**2 + 0.5 * y**2
    //
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 0, -1.0, -1.0, 1.5);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 1,  0.0, -1.0, 0.5);
    if (IsOk) IsOk = SurfaceFitAppendPointAtEnd(TestSurfaceFit,     1.0, -1.0, 1.5);
    if (IsOk) IsOk = SurfaceFitAppendPointAtEnd(TestSurfaceFit,    -1.0,  0.0, 1.0);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 4,  0.0,  0.0, 0.0);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 5,  1.0,  0.0, 1.0);
    if (IsOk) IsOk = SurfaceFitAppendPointAtEnd(TestSurfaceFit,    -1.0,  1.0, 1.5);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 7,  0.0,  1.0, 0.5);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 8,  1.0,  1.0, 1.5);
    if (IsOk) IsOk = SurfaceFitCompute(TestSurfaceFit, SurfFitType_Quadratic);

    // Verify the result a=b=c=e=0.0, d=f=1.0
    if (IsOk)
    {
        LgIndex_t Count = SurfaceFitGetCoefCount(TestSurfaceFit);
        double a, b, c, d, e, f;

        if (Count != 6) IsOk = FALSE;

        if (IsOk)
        {
            a = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 0);
            b = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 1);
            c = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 2);
            d = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 3);
            e = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 4);
            f = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 5);

            if (a > 1.0e-5 || a < -1.0e-5) IsOk = FALSE;
            if (IsOk)
                if (b > 1.0e-5 || b < -1.0e-5) IsOk = FALSE;
            if (IsOk)
                if (c > 1.0e-5 || c < -1.0e-5) IsOk = FALSE;
            if (IsOk)
                if (e > 1.0e-5 || e < -1.0e-5) IsOk = FALSE;
            if (IsOk)
                if (d < 0.9999999 || d > 1.000001) IsOk = FALSE;
            if (IsOk)
                if (f < 0.4999999 || f > 0.500001) IsOk = FALSE;
        }
    }

    // Compute critical points and verify
    if (IsOk)
        IsOk = SurfaceFitCompCritPoint(TestSurfaceFit);





    //
    //  Simple Test 2:  p = xp**2 - 0.5 * yp**2
    //     where xp,yp are translated/rotated (angle alpha) coordinates
    //     xp =   (x - x0) * cos(alpha) + (y - y0) * sin(alpha)
    //     yp = - (x - x0) * sin(alpha) + (y - y0) * cos(alpha)
    //       choose (x0, y0) = (0.6, 0.3)  and  alpha = pi/6
    //
    SurfaceFitClear(TestSurfaceFit);
    if (IsOk)
    {
        double x0 = 0.6;
        double y0 = 0.3;
        double alpha = 3.1415927 / 6.0;
        double sna = sin(alpha);
        double csa = cos(alpha);

        if (IsOk)
        {
            double xx = -1.0;
            double yy = -1.0;
            double xp = (xx - x0) * csa + (yy - y0) * sna;
            double yp = - (xx - x0) * sna + (yy - y0) * csa;
            double pp = xp * xp - 0.5 * yp * yp;
            IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 0, xx, yy, pp);
        }
        if (IsOk)
        {
            double xx =  0.0;
            double yy = -1.0;
            double xp = (xx - x0) * csa + (yy - y0) * sna;
            double yp = - (xx - x0) * sna + (yy - y0) * csa;
            double pp = xp * xp - 0.5 * yp * yp;
            IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 1,  xx, yy, pp);
        }
        if (IsOk)
        {
            double xx =  1.0;
            double yy = -1.0;
            double xp = (xx - x0) * csa + (yy - y0) * sna;
            double yp = - (xx - x0) * sna + (yy - y0) * csa;
            double pp = xp * xp - 0.5 * yp * yp;
            IsOk = SurfaceFitAppendPointAtEnd(TestSurfaceFit,     xx, yy, pp);
        }
        if (IsOk)
        {
            double xx = -1.0;
            double yy =  0.0;
            double xp = (xx - x0) * csa + (yy - y0) * sna;
            double yp = - (xx - x0) * sna + (yy - y0) * csa;
            double pp = xp * xp - 0.5 * yp * yp;
            IsOk = SurfaceFitAppendPointAtEnd(TestSurfaceFit,    xx,  yy, pp);
        }
        if (IsOk)
        {
            double xx =  0.0;
            double yy =  0.0;
            double xp = (xx - x0) * csa + (yy - y0) * sna;
            double yp = - (xx - x0) * sna + (yy - y0) * csa;
            double pp = xp * xp - 0.5 * yp * yp;
            IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 4,  xx,  yy, pp);
        }
        if (IsOk)
        {
            double xx =  1.0;
            double yy =  0.0;
            double xp = (xx - x0) * csa + (yy - y0) * sna;
            double yp = - (xx - x0) * sna + (yy - y0) * csa;
            double pp = xp * xp - 0.5 * yp * yp;
            IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 5,  xx,  yy, pp);
        }
        if (IsOk)
        {
            double xx = -1.0;
            double yy =  1.0;
            double xp = (xx - x0) * csa + (yy - y0) * sna;
            double yp = - (xx - x0) * sna + (yy - y0) * csa;
            double pp = xp * xp - 0.5 * yp * yp;
            IsOk = SurfaceFitAppendPointAtEnd(TestSurfaceFit,    xx,  yy, pp);
        }
        if (IsOk)
        {
            double xx =  0.0;
            double yy =  1.0;
            double xp = (xx - x0) * csa + (yy - y0) * sna;
            double yp = - (xx - x0) * sna + (yy - y0) * csa;
            double pp = xp * xp - 0.5 * yp * yp;
            IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 7,  xx,  yy, pp);
        }
        if (IsOk)
        {
            double xx =  1.0;
            double yy =  1.0;
            double xp = (xx - x0) * csa + (yy - y0) * sna;
            double yp = - (xx - x0) * sna + (yy - y0) * csa;
            double pp = xp * xp - 0.5 * yp * yp;
            IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 8,  xx,  yy, pp);
        }

        if (IsOk) IsOk = SurfaceFitCompute(TestSurfaceFit, SurfFitType_Quadratic);


        // Verify the result a=b=c=e=0.0, d=f=1.0
        if (IsOk)
        {
            LgIndex_t Count = SurfaceFitGetCoefCount(TestSurfaceFit);
            double a, b, c, d, e, f;

            double ac = x0 * x0 * (csa * csa - 0.5 * sna * sna) +
                        3.0 * x0 * y0 * csa * sna +
                        y0 * y0 * (sna * sna - 0.5 * csa * csa);
            double bc = x0 * (-2.0 * csa * csa + sna * sna) -
                        3.0 * y0 * csa * sna;
            double cc = y0 * (-2.0 * sna * sna + csa * csa) -
                        3.0 * x0 * csa * sna;
            double dc = csa * csa - 0.5 * sna * sna;
            double ec = 3.0 * csa * sna;
            double fc = sna * sna - 0.5 * csa * csa;

            // TEST
            double vc11 = ac + bc + cc + dc + ec + fc;
            double v11  = SurfaceFitInterpolateVar(TestSurfaceFit, 1.0, 1.0);

            if (Count != 6) IsOk = FALSE;

            // TODO: Test temporarily disabled. Isn't working. Find out why
            if (IsOk)
            {
                a = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 0);
                b = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 1);
                c = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 2);
                d = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 3);
                e = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 4);
                f = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 5);

                if (a >= 0.0)
                {
                    if (a < 0.9999999 * ac || a > 1.000001 * ac) IsOk = FALSE;
                }
                else
                {
                    if (a > 0.9999999 * ac || a < 1.000001 * ac) IsOk = FALSE;
                }

                if (IsOk)
                {
                    if (b >= 0.0)
                    {
                        if (b < 0.9999999 * bc || b > 1.000001 * bc) IsOk = FALSE;
                    }
                    else
                    {
                        if (b > 0.9999999 * bc || b < 1.000001 * bc) IsOk = FALSE;
                    }
                }

                if (IsOk)
                {
                    if (c >= 0.0)
                    {
                        if (c < 0.9999999 * cc || c > 1.000001 * cc) IsOk = FALSE;
                    }
                    else
                    {
                        if (c > 0.9999999 * cc || c < 1.000001 * cc) IsOk = FALSE;
                    }
                }

                if (IsOk)
                {
                    if (d >= 0.0)
                    {
                        if (d < 0.9999999 * dc || d > 1.000001 * dc) IsOk = FALSE;
                    }
                    else
                    {
                        if (d > 0.9999999 * dc || d < 1.000001 * dc) IsOk = FALSE;
                    }
                }

                if (IsOk)
                {
                    if (e >= 0.0)
                    {
                        if (e < 0.9999999 * ec || e > 1.000001 * ec) IsOk = FALSE;
                    }
                    else
                    {
                        if (e > 0.9999999 * ec || e < 1.000001 * ec) IsOk = FALSE;
                    }
                }

                if (IsOk)
                {
                    if (f >= 0.0)
                    {
                        if (f < 0.9999999 * fc || f > 1.000001 * fc) IsOk = FALSE;
                    }
                    else
                    {
                        if (f > 0.9999999 * fc || f < 1.000001 * fc) IsOk = FALSE;
                    }
                }
            }
        }
    }

    // Compute critical points and verify
    if (IsOk)
        IsOk = SurfaceFitCompCritPoint(TestSurfaceFit);

    if (IsOk)
    {
        double XCrt = TestSurfaceFit->c1_crt;
        double YCrt = TestSurfaceFit->c2_crt;

        if (XCrt < 0.59999 || XCrt > 0.60001) IsOk = FALSE;
        if (IsOk && YCrt < 0.29999 || YCrt > 0.30001) IsOk = FALSE;
        if (IsOk)
        {
            double RotationAngle = TestSurfaceFit->RotationAngle;
            double Lambda1 = TestSurfaceFit->lambda1;
            double Lambda2 = TestSurfaceFit->lambda2;

            // Rotation angle (eigenvectors)
            if (RotationAngle >= 0.0)
            {
                if (RotationAngle < 0.5235978 || RotationAngle > 0.5235998) IsOk = FALSE;
            }

            // First Eigenvalue (should be +1.0)
            if (Lambda1 < 1.999999 || Lambda1 > 2.00001) IsOk = FALSE;

            // Second Eigenvalue (should be -0.5)
            if (Lambda2 > -0.999999 || Lambda2 < -1.00001) IsOk = FALSE;
        }
    }



    //  Simple Test:  Linear surface fit
    //
    //                   ^ y = Coord2
    //                   |
    //                   |
    //                   |
    //         ______________________ (1, 1)
    //        |          |           |
    //        |          |           |
    //        |          |           |
    //        |          |           |
    //        |----------------------|  ---------> x = Coord1
    //        |          |           |
    //        |          |           |
    //        |          |           |
    //        |__________|___________|
    //      (-1,-1)
    //
    //  Simple Test 1:  p = (x - 0.5) + 0.5 * (y + 0.5)
    //
    SurfaceFitClear(TestSurfaceFit);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 0, -1.0, -1.0, -1.75);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 1,  0.0, -1.0, -0.75);
    if (IsOk) IsOk = SurfaceFitAppendPointAtEnd(TestSurfaceFit,     1.0, -1.0,  0.25);
    if (IsOk) IsOk = SurfaceFitAppendPointAtEnd(TestSurfaceFit,    -1.0,  0.0, -1.25);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 4,  0.0,  0.0, -0.25);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 5,  1.0,  0.0,  0.75);
    if (IsOk) IsOk = SurfaceFitAppendPointAtEnd(TestSurfaceFit,    -1.0,  1.0, -0.75);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 7,  0.0,  1.0,  0.25);
    if (IsOk) IsOk = SurfaceFitSetPointAtOffset(TestSurfaceFit, 8,  1.0,  1.0,  1.25);
    if (IsOk) IsOk = SurfaceFitCompute(TestSurfaceFit, SurfFitType_Planar);

    // Verify the result a=-0.25, b=1.0, c=0.5
    if (IsOk)
    {
        LgIndex_t Count = SurfaceFitGetCoefCount(TestSurfaceFit);
        double a, b, c;

        if (Count != 3) IsOk = FALSE;

        if (IsOk)
        {
            a = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 0);
            b = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 1);
            c = SurfaceFitGetCoefFromOffset(TestSurfaceFit, 2);

            if (a > -0.249999 || a < -0.250001) IsOk = FALSE;
            if (IsOk)
                if (b > 1.00001 || b < 0.99999) IsOk = FALSE;
            if (IsOk)
                if (c > 0.50001 || c < 0.49999) IsOk = FALSE;
        }
    }



    // Clean up
    SurfaceFitDealloc(&TestSurfaceFit);

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}
