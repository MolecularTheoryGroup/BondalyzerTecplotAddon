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
#include "TASSERT.h"
#include "ALLOC.h"
#include "LINEARALGEBRA.h"

#define MAX_STATUS_MSG_LEN 127


void NullMatrix(Matrix_s *matrixptr)
{
    matrixptr->val = NULL;
    matrixptr->raw_data = NULL;
}

void FreeMatrix(Matrix_s *matrixptr)
{
    if (matrixptr)
    {
        if (matrixptr->val)
        {
            FREE_ARRAY(matrixptr->val, "Matrix pointers");
            matrixptr->val = NULL;
        }
        if (matrixptr->raw_data)
        {
            FREE_ARRAY(matrixptr->raw_data, "Matrix data");
            matrixptr->raw_data = NULL;
        }
    }
} /* FreeMatrix() */


Matrix_s CreateMatrix(LgIndex_t dim1, LgIndex_t dim2)
{
    Matrix_s matrix;
    matrix.val = ALLOC_ARRAY(dim1, double *, "Matrix pointers");
    matrix.raw_data = ALLOC_ARRAY(dim1 * dim2, double, "Matrix data");

    if (matrix.val == NULL || matrix.raw_data == NULL)
    {
        FreeMatrix(&matrix);
    }
    else
    {
        LgIndex_t ii;

        for (ii = 0; ii < dim1; ii++)
            matrix.val[ii] = matrix.raw_data + ii * dim2;
    }
    return matrix;

} /* CreateMatrix() */


/*
 * cmludcmp: perform LU decomposition of matrix
 *  (Adapted from "Numerical Recipes in C", Page 41)
 *   matrix       (ndim)x(ndim)
 *   ndim         dimension of matrix
 *   indx         (ndim) array of LgIndex_t to hold row permutation
 *   TmpSpace  temp space for ndim float's
 */
Boolean_t cmludcmp(Matrix_s   matrix,
                   LgIndex_t  ndim,
                   LgIndex_t *indx,
                   double    *TmpSpace)
{
    LgIndex_t ii = 0;
    LgIndex_t imax = 0;
    LgIndex_t jj = 0;
    LgIndex_t pp = 0;
    double    big = 0.0;
    double    sum = 0.0;

    /* loop over rows to get scaling info */
    for (ii = 0; ii < ndim; ii++)
    {
        double *row = matrix.val[ii];
        big = 0;
        for (jj = 0; jj < ndim; jj++)
        {
            double tmp = fabs(row[jj]);
            if (tmp > big)
                big = tmp;
        }
        if (big == 0)
        {
            TecUtilDialogErrMsg("Singular matrix.");
            return(FALSE);
        }
        TmpSpace[ii] = 1 / big;
    }

    for (jj = 0; jj < ndim; jj++)
    {
        for (ii = 0; ii < jj; ii++)
        {
            double *row = matrix.val[ii];
            sum = row[jj];
            for (pp = 0; pp < ii; pp++)
                sum -= row[pp] * matrix.val[pp][jj];
            row[jj] = sum;
        }
        big = 0.0;
        for (; ii < ndim; ii++)
        {
            double tmp;
            double *row = matrix.val[ii];
            sum = row[jj];
            for (pp = 0; pp < jj; pp++)
                sum -= row[pp] * matrix.val[pp][jj];
            row[jj] = sum;
            tmp = TmpSpace[ii] * fabs(sum);
            if (tmp >= big)
            {
                big = tmp;
                imax = ii;
            }
        }

        if (jj != imax)
        {
            /* swap rows */
            double *tmp = matrix.val[imax];
            matrix.val[imax] = matrix.val[jj];
            matrix.val[jj] = tmp;
            TmpSpace[imax] = TmpSpace[jj];
        }
        indx[jj] = imax;

        if (-SMALLDOUBLE < matrix.val[jj][jj] && matrix.val[jj][jj] < SMALLDOUBLE)
            matrix.val[jj][jj] = SMALLDOUBLE;
        if (jj != ndim - 1)
        {
            double tmp = 1 / matrix.val[jj][jj];
            for (ii = jj + 1; ii < ndim; ii++)
                matrix.val[ii][jj] *= tmp;
        }
    }

    return(TRUE);

} /* cmludcmp() */


/*
 * cmlubksb: solves AX=B, bb is both input (B) and return (X) matrix
 *  A is LU decomposed matrix (from cmludcmp)
 *  (Adapted from "Numerical Recipes in C", page 44)
 *  matrix  LU decomp matrix (from cmludcmp)
 *  ndim    dimension of aa, indx, and bb
 *  indx    permutation matrix of aa (from cmludcmp)
 *  bb      input (B) and return (X) matrix
 */
Boolean_t cmlubksb(Matrix_s   matrix,
                   LgIndex_t  ndim,
                   LgIndex_t *indx,
                   double    *bb)
{
    LgIndex_t ii, i2, ip, jj;
    double sum;

    i2 = -1;
    for (ii = 0; ii < ndim; ii++)
    {
        ip = indx[ii];
        sum = bb[ip];
        bb[ip] = bb[ii];
        if (i2 != -1)
        {
            double *row = matrix.val[ii];
            for (jj = i2; jj < ii; jj++)
                sum -= row[jj] * bb[jj];
        }
        else
        {
            if (sum != 0)
                i2 = ii;
        }
        bb[ii] = sum;
    }
    for (ii = ndim - 1; ii > -1; ii--)
    {
        double *row = matrix.val[ii];
        double tmp;
        sum = bb[ii];
        for (jj = ii + 1; jj < ndim; jj++)
            sum -= row[jj] * bb[jj];
        tmp = row[ii];
        if (-SMALLDOUBLE < tmp && tmp < SMALLDOUBLE)
        {
            //TecUtilDialogErrMsg("Singular matrix.");
            return FALSE;
        }
        else
        {
            double tmp2 = sum / tmp;
            if ((tmp2 > LARGEDOUBLE) || (tmp2 < -LARGEDOUBLE))
            {
                //TecUtilDialogErrMsg("Value overflow.");
                return FALSE;
            }
            else
                bb[ii] = tmp2;
        }
    }

    return(TRUE);

} /* cmlubksb() */


