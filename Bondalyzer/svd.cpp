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
#include "SVD.h"
#include <string.h>


/* * LOCAL CONSTANTS * */
#define VerySmall SMALLDOUBLE
#define VeryLarge LARGEDOUBLE


#define SIGN(a,b) ((b)<0.0 ? -fabs(a): fabs(a))

#define COMPUTE_PATHAG(a,b) ((at = fabs(a)) > (bt = fabs(b)) ? \
(ct = bt/at, at*sqrt(1.0+ct*ct)): (bt ? (ct = at/bt, bt*sqrt(1.0+ct*ct)): 0.0))




/**
 * Allocate memory for a vector (function from Numerical Recipes in C).
 */
double *SVDAllocVector(int nl, int nh)
{
    double *v;
    double *Result;

    REQUIRE(nl <= nh);

    v = ALLOC_ARRAY((unsigned)(nh - nl + 1), double, "vector");
    Result = v - nl;

    ENSURE(VALID_REF(v));
    ENSURE(VALID_REF(Result));
    return Result;
}


/**
 * Free memory allocated for SVDAllocVector (function from Numerical Recipes in C).
 */
void SVDFreeVector(double *v, int nl, int nh)
{
    double *Vector = v + nl;

    REQUIRE(VALID_REF(v));
    REQUIRE(VALID_REF(Vector));
    REQUIRE(nl <= nh);

    FREE_ARRAY(Vector, "vector");
}


/**
 * Allocate memory for a matrix (function from Numerical Recipes in C).
 */
double **SVDAllocMatrix(int nrl, int nrh, int ncl, int nch)
/* Allocates a double matrix with range [nrl..nrh][ncl..nch]. */
{
    int i;
    double **m;
    double **Result;

    REQUIRE(nrl <= nrh);
    REQUIRE(ncl <= nch);

    /* Allocate pointers to rows. */
    m = ALLOC_ARRAY((unsigned)(nrh - nrl + 1), double*, "matrix rows");
    ENSURE(VALID_REF(m));
    m -= nrl;
    Result = m;
    ENSURE(VALID_REF(Result));

    /* Allocate rows and set pointers to them. */
    for (i = nrl; i <= nrh; i++)
    {
        m[i] = ALLOC_ARRAY((unsigned)(nch - ncl + 1), double, "matrix row");
        // ENSURE(VALID_REF(m[i]));
        // ENSURE(VALID_REF(m[i]-ncl));
        m[i] -= ncl;
    }
    return Result;
}


/**
 * Free memory allocated for SVDAllocMatrix (function from Numerical Recipes in C).
 */
void SVDFreeMatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **Matrix;

    REQUIRE(VALID_REF(m));
    // REQUIRE(VALID_REF(m+nrl));

    for (i = nrh; i >= nrl; i--)
    {
        double *MatrixRow = m[i] + ncl;
        REQUIRE(VALID_REF(m[i]));
        REQUIRE(VALID_REF(MatrixRow));
        FREE_ARRAY(MatrixRow, "matrix row");
    }
    Matrix = m + nrl;
    FREE_ARRAY(Matrix, "matrix rows");
}



/**
 * Function to perform the singular value decomposition (function from
 * Numerical Recipes in C).
 *
 * Given a matrix a[1..m][1..n], this routine computes its singular value
 * decomposition, A = U?W?V T. Thematrix U replaces a on output. The diagonal
 * matrix of singular values W is output as a vector w[1..n]. The matrix V (not
 * the transpose V T ) is output as v[1..n][1..n].
 */
void ComputeSVD(double **a, int m, int n, double w[], double **v)
{
    int flag, i, its, j, jj, k;
    int l  = 0; /* ... quiet compiler */
    int nm = 0; /* ... quiet compiler */
    double anorm, c, f, g, h, s, scale, x, y, z, *rv1;
    double at, bt, ct;
    rv1 = SVDAllocVector(1, n);
    g = scale = anorm = 0.0;

    /* Householder reduction to bidiagonal form.  */
    for (i = 1; i <= n; i++)
    {
        l = i + 1;
        CHECK(1 <= i && i <= n);
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i <= m)
        {
            for (k = i; k <= m; k++)
            {
                CHECK((1 <= k && k <= m) && (1 <= i && i <= n));
                scale += ABS(a[k][i]);
            }
            if (scale)
            {
                for (k = i; k <= m; k++)
                {
                    CHECK((1 <= k && k <= m) && (1 <= i && i <= n));
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                CHECK((1 <= i && i <= m) && (1 <= i && i <= n));
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                CHECK((1 <= i && i <= m) && (1 <= i && i <= n));
                a[i][i] = f - g;
                for (j = l; j <= n; j++)
                {
                    for (s = 0.0, k = i; k <= m; k++)
                    {
                        CHECK((1 <= k && k <= m) && (1 <= i && i <= n) && (1 <= j && j <= n));
                        s += a[k][i] * a[k][j];
                    }
                    f = s / h;
                    for (k = i; k <= m; k++)
                    {
                        CHECK((1 <= k && k <= m) && (1 <= i && i <= n) && (1 <= j && j <= n));
                        a[k][j] += f * a[k][i];
                    }
                }
                for (k = i; k <= m; k++)
                {
                    CHECK((1 <= k && k <= m) && (1 <= i && i <= n));
                    a[k][i] *= scale;
                }
            }
        }
        CHECK(1 <= i && i <= n);
        w[i] = scale * g;
        g = s = scale = 0.0;
        if (i <= m && i != n)
        {
            for (k = l ; k <= n; k++)
            {
                CHECK((1 <= i && i <= m) && (1 <= k && k <= n));
                scale += ABS(a[i][k]);
            }
            if (scale)
            {
                for (k = l; k <= n; k++)
                {
                    CHECK((1 <= i && i <= m) && (1 <= k && k <= n));
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                CHECK((1 <= i && i <= m) && (1 <= l && l <= n));
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = f - g;
                for (k = l; k <= n; k++)
                {
                    CHECK((1 <= i && i <= m) && (1 <= k && k <= n));
                    rv1[k] = a[i][k] / h;
                }
                if (i != m)
                {
                    for (j = l; j <= m; j++)
                    {
                        for (s = 0.0, k = l; k <= n; k++)
                        {
                            CHECK((1 <= i && i <= m) && (1 <= i && i <= m) && (1 <= k && k <= n));
                            s += a[j][k] * a[i][k];
                        }
                        for (k = l; k <= n; k++)
                        {
                            CHECK((1 <= j && j <= m) && (1 <= k && k <= n));
                            a[j][k] += s * rv1[k];
                        }
                    }
                }
                for (k = l; k <= n; k++)
                {
                    CHECK((1 <= i && i <= m) && (1 <= k && k <= n));
                    a[i][k] *= scale;
                }
            }
        }
        CHECK(1 <= i && i <= n);
        anorm = MAX(anorm, (ABS(w[i]) + ABS(rv1[i])));
    }
    /* Accumulation of right-hand transformations */
    for (i = n; i >= 1; i--)
    {
        if (i < n)
        {
            if (g)
            {
                for (j = l; j <= n; j++)    /* Double division to avoid possible under?ow. */
                {
                    CHECK((1 <= i && i <= n) && (1 <= j && j <= n) && (1 <= l && l <= n));
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                }
                for (j = l; j <= n; j++)
                {
                    for (s = 0.0, k = l; k <= n; k++)
                    {
                        CHECK((1 <= i && i <= m) && (1 <= k && k <= n));
                        CHECK((1 <= k && k <= n) && (1 <= j && j <= n));
                        s += a[i][k] * v[k][j];
                    }
                    for (k = l; k <= n; k++)
                    {
                        CHECK((1 <= k && k <= n) && (1 <= j && j <= n) && (1 <= i && i <= n));
                        v[k][j] += s * v[k][i];
                    }
                }
            }
            for (j = l; j <= n; j++)
            {
                CHECK((1 <= i && i <= n) && (1 <= j && j <= n));
                v[i][j] = v[j][i] = 0.0;
            }
        }
        CHECK(1 <= i && i <= n);
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    /* Accumulation of left-hand transformations. */
    for (i = MIN(m, n); i >= 1; i--)
    {
        l = i + 1;
        CHECK(1 <= i && i <= n);
        g = w[i];
        for (j = l; j <= n; j++)
        {
            CHECK((1 <= i && i <= m) && (1 <= j && j <= n));
            a[i][j] = 0.0;
        }
        if (g)
        {
            g = 1.0 / g;
            if (i < n)
            {
                for (j = l; j <= n; j++)
                {
                    for (s = 0.0, k = l; k <= m; k++)
                    {
                        CHECK((1 <= k && k <= m) && (1 <= i && i <= n) && (1 <= j && j <= n));
                        s += a[k][i] * a[k][j];
                    }
                    CHECK((1 <= i && i <= m) && (1 <= i && i <= n));
                    f = (s / a[i][i]) * g;
                    for (k = i; k <= m; k++)
                    {
                        CHECK((1 <= k && k <= m) && (1 <= j && j <= n) && (1 <= i && i <= n));
                        a[k][j] += f * a[k][i];
                    }
                }
            }
            for (j = i; j <= m; j++)
            {
                CHECK((1 <= j && j <= m) && (1 <= i && i <= n));
                a[j][i] *= g;
            }
        }
        else
        {
            for (j = i; j <= m; j++)
            {
                CHECK((1 <= j && j <= m) && (1 <= i && i <= n));
                a[j][i] = 0.0;
            }
        }
        CHECK((1 <= i && i <= m) && (1 <= i && i <= n));
        ++a[i][i];
    }
    /*
     * Diagonalization of the bidiagonal form: Loop over singular values,
     * and over allowed iterations.
     */
    for (k = n; k >= 1; k--)
    {
        const int computeIterations = 30;
        for (its = 1; its <= computeIterations; its++)
        {
            flag = 1;

            /* Test for splitting. */
            for (l = k; l >= 1; l--)
            {
                nm = l - 1;  /* Note that rv1[1] is always zero.  */
                CHECK(1 <= l && l <= n);
                if ((double)(ABS(rv1[l]) + anorm) == anorm)
                {
                    flag = 0;
                    break;
                }
                CHECK(1 <= nm && nm <= n);
                if ((double)(ABS(w[nm]) + anorm) == anorm) break;
            }
            if (flag)
            {
                c = 0.0;   /* Cancellation of rv1[l], if l > 1. */
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    CHECK(1 <= i && i <= n);
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ((double)(ABS(f) + anorm) == anorm) break;
                    g = w[i];
                    h = COMPUTE_PATHAG(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 1; j <= m; j++)
                    {
                        CHECK((1 <= j && j <= m) && (1 <= nm && nm <= n) && (1 <= i && i <= n));
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y * c + z * s;
                        a[j][i] = z * c - y * s;
                    }
                }
            }
            CHECK(1 <= k && k <= n);
            z = w[k];

            /* Convergence. */
            if (l == k)
            {
                if (z < 0.0)  /* Singular value is made nonnegative. */
                {
                    CHECK(1 <= k && k <= n);
                    w[k] = -z;
                    for (j = 1; j <= n; j++)
                    {
                        CHECK((1 <= j && j <= n) && (1 <= k && k <= n));
                        v[j][k] = -v[j][k];
                    }
                }
                break;
            }
            /* if (its == computeIterations) nerror("no convergence in 'computeIterations' ComputeSVD iterations"); */
            CHECK(1 <= l && l <= n);
            x = w[l];   /* Shift from bottom 2-by-2 minor. */
            nm = k - 1;
            CHECK((1 <= nm && nm <= n) && (1 <= k && k <= n));
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = COMPUTE_PATHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0;   /* Next QR transformation: */
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                CHECK(1 <= i && i <= n);
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = COMPUTE_PATHAG(f, h);
                CHECK(1 <= j && j <= n);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 1; jj <= n; jj++)
                {
                    CHECK((1 <= jj && jj <= n) && (1 <= j && j <= n) && (1 <= i && i <= n));
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = COMPUTE_PATHAG(f, h);
                CHECK(1 <= j && j <= n);
                w[j] = z;   /* Rotation can be arbitrary if z = 0. */
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;

                for (jj = 1; jj <= m; jj++)
                {
                    CHECK((1 <= jj && jj <= m) && (1 <= i && i <= n) && (1 <= j && j <= n));
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }
            CHECK((1 <= l && l <= n) && (1 <= k && k <= n));
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    SVDFreeVector(rv1, 1, n);
}


/**
 * Function to solve system using the singular value decomposition (function
 * from Numerical Recipes in C).
 *
 * Solves A?X = B for a vector X, where A is specied by the arrays
 * u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by ComputeSVD. m and n are the
 * dimensions of a, and will be equal for square matrices. b[1..m] is the input
 * right-hand side. x[1..n] is the output solution vector.  No input quantities
 * are destroyed, so the routine may be called sequentially with di.erent b?s.
 */
void BackSubstituteSVD(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
    int jj, j, i;
    double s, *tmp;
    tmp = SVDAllocVector(1, n);
    for (j = 1; j <= n; j++)  /* Calculate UT B. */
    {
        s = 0.0;
        CHECK(1 <= j && j <= n);
        if (w[j])   /* Nonzero result only if wj is nonzero. */
        {
            for (i = 1; i <= m; i++)
            {
                CHECK((1 <= i && i <= m) && (1 <= j && j <= n));
                s += u[i][j] * b[i];
            }
            CHECK(1 <= j && j <= n);
            s /= w[j];   /* This is the divide by wj. */
        }
        tmp[j] = s;
    }
    for (j = 1; j <= n; j++)   /* Matrix multiply by V to get answer. */
    {
        s = 0.0;
        for (jj = 1; jj <= n; jj++)
        {
            CHECK((1 <= j && j <= n) && (1 <= jj && jj <= n));
            s += v[j][jj] * tmp[jj];
        }
        CHECK(1 <= j && j <= n);
        x[j] = s;
    }
    SVDFreeVector(tmp, 1, n);
}


