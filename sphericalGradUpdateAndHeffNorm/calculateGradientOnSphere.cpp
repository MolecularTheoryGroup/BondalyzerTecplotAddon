#include "calculateGradientOnSphere.h"

namespace tpcsm {

// Begin Functions from Numerical Recipes in C from Tecplot: could substitute favorite linear algebra library or routines here

#define ALLOC_ARRAY(N,type) (type *)malloc((N)*sizeof(type))
#define ALLOC_ITEM(type) (type *)malloc(sizeof(type))
#ifdef _DEBUG
/* NOTE: the pointer is set to 0xFFFF after the free for debug   */
/*       versions in the hopes of catching invalid pointer usage */
#define FREE_ARRAY(ptr)  do { free((void *)(ptr)); *((void **)&(ptr)) = (void *)0xFFFF; } while (0)
#define FREE_ITEM(ptr)   do { free((void *)(ptr)); *((void **)&(ptr)) = (void *)0xFFFF; } while (0)
#else
#define FREE_ARRAY(ptr)  free((void *)(ptr))
#define FREE_ITEM(ptr)   free((void *)(ptr))
#endif

#define VALID_REF(ptr) (ptr!=NULL)
#define SMALLFLOAT     (1.0e-38)
#define LARGEFLOAT     (1.0e+38)
#define SIGN(a,b)      ((b)<0 ? -ABS((a)): ABS((a)))
#define COMPUTE_PATHAG(a,b) ((at = ABS(a)) > (bt = ABS(b)) ? \
                            (ct = bt/at, at*sqrt(1.0+ct*ct)): (bt ? (ct = at/bt, bt*sqrt(1.0+ct*ct)): 0.0))
#define ABS(a)         ((a) >= 0 ? (a) : -(a) )
#define MAX(a,b)       ((a) > (b) ? (a) : (b) )
#define MIN(a,b)       ((a) < (b) ? (a) : (b) )



/**
* Allocate memory for a vector (function from Numerical Recipes in C).
*/
static double *nrcAllocVector(int nl, int nh)
{
    double *v;
    double *Result;

    REQUIRE(nl <= nh);

    v = ALLOC_ARRAY((unsigned)(nh - nl + 1), double);
    Result = v - nl;

    ENSURE(VALID_REF(v));
    ENSURE(VALID_REF(Result));
    return Result;
}


/**
* Free memory allocated for nrcAllocVector (function from Numerical Recipes in C).
*/
static void nrcFreeVector(double *v, int nl, int ASSERT_ONLY(nh))
{
    double *Vector = v + nl;

    REQUIRE(VALID_REF(v));
    REQUIRE(VALID_REF(Vector));
    REQUIRE(nl <= nh);

    FREE_ARRAY(Vector);
}

/**
* Allocate memory for a matrix (function from Numerical Recipes in C).
*/
static double **nrcAllocMatrix(int nrl, int nrh, int ncl, int nch)
/* Allocates a double matrix with range [nrl..nrh][ncl..nch]. */
{
    REQUIRE(nrl <= nrh);
    REQUIRE(ncl <= nch);

    /* Allocate pointers to rows. */
    double **m = ALLOC_ARRAY((unsigned)(nrh - nrl + 1), double*);
    ENSURE(VALID_REF(m));

    m -= nrl;

    double **result = m;
    ENSURE(VALID_REF(result));

    /* Allocate rows and set pointers to them. */
    for (Index_t i = nrl; i <= nrh; i++)
    {
        m[i] = ALLOC_ARRAY((unsigned)(nch - ncl + 1), double);
        ENSURE(VALID_REF(m[i]));
        ENSURE(VALID_REF(m[i] - ncl));
        m[i] -= ncl;
    }

    return result;
}


/**
* Free memory allocated for nrcAllocMatrix (function from Numerical Recipes in C).
*/
static void nrcFreeMatrix(double **m, int nrl, int nrh, int ncl, int /*nch*/)
{
    REQUIRE(VALID_REF(m));
    REQUIRE(VALID_REF(m + nrl));

    for (Index_t i = nrh; i >= nrl; i--)
    {
        double* matrixRow = m[i] + ncl;
        REQUIRE(VALID_REF(m[i]));
        REQUIRE(VALID_REF(matrixRow));
        FREE_ARRAY(matrixRow);
    }

    double** matrix = m + nrl;
    FREE_ARRAY(matrix);
}

/**
* Function to perform the singular value decomposition (function from
* Numerical Recipes in C).
*
* Given a matrix a[1..m][1..n], this routine computes its singular value
* decomposition, A = U?W?V T. The matrix U replaces a on output. The diagonal
* matrix of singular values W is output as a vector w[1..n]. The matrix V (not
* the transpose V T ) is output as v[1..n][1..n].
*/
static void nrcComputeSVD(double **a, int m, int n, double w[], double **v)
{
    int flag, i, its, j, jj, k;
    int l = 0; /* ... quiet compiler */
    int nm = 0; /* ... quiet compiler */
    double anorm, c, f, g, h, s, scale, x, y, z, *rv1;
    double at, bt, ct;
    rv1 = nrcAllocVector(1, n);
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
            for (k = l; k <= n; k++)
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
                            CHECK((1 <= i && i <= m) && (1 <= j && j <= m) && (1 <= k && k <= n));
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
        for (its = 1; its <= 30; its++)
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
            /* if (its == 30) nerror("no convergence in 30 ComputeSVD iterations"); */
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
    nrcFreeVector(rv1, 1, n);
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
static void nrcBackSubstituteSVD(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
    int jj, j, i;
    double s, *tmp;
    tmp = nrcAllocVector(1, n);
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
    nrcFreeVector(tmp, 1, n);
}

// END Functions from Numerical Recipes in C from Tecplot

// BEGIN Functions from Tecplot to use NRC code above

struct _MLSDeriv_s // MOVING LEAST SQUARES
{
    double ddx;
    double ddy;
    double ddz;
    double d2dxx;
    double d2dyy;
    double d2dxy;
    double d2dzz;
    double d2dxz;
    double d2dyz;
    double PrimaryNodeGradContribCoef;
    bool IsSingular;
    Index_t maxNumInfluencerPts;
    int       maxNumBasisFunctions;
    double    **a;  /* MxNumInfNodesOrElems x MxNumDerivTypes */
    double    **v;  /* MxNumDerivTypes x MxNumDerivTypes */
    double    *b;   /* MxNumInfNodesOrElems */
    double    *x;   /* MxNumDerivTypes */
    double    *w;   /* MxNumDerivTypes */
};

typedef struct _MLSDeriv_s  *MLSDeriv_pa;

/*
*  Allocate the coefficient matrices and vectors needed for the Moving Least
*  Square derivatives.
*/
static MLSDeriv_pa allocMLSMatrices(
    Index_t numInfluencerPts,
    Index_t numDerivTypes)
{
    REQUIRE(numInfluencerPts > 2);
    REQUIRE(numDerivTypes > 2);

    MLSDeriv_pa result = ALLOC_ITEM(struct _MLSDeriv_s);

    if (result != NULL)
    {
        std::memset(result, 0, sizeof(*result));
        result->maxNumBasisFunctions = numDerivTypes;
        result->maxNumInfluencerPts = numInfluencerPts;

        result->a = nrcAllocMatrix(1, result->maxNumInfluencerPts, 1, result->maxNumBasisFunctions);
        result->v = nrcAllocMatrix(1, result->maxNumBasisFunctions, 1, result->maxNumBasisFunctions);

        result->b = nrcAllocVector(1, result->maxNumInfluencerPts);
        result->x = nrcAllocVector(1, result->maxNumBasisFunctions);
        result->w = nrcAllocVector(1, result->maxNumBasisFunctions);
    }

    ENSURE(VALID_REF(result));
    ENSURE(VALID_REF(result->a));
    ENSURE(VALID_REF(result->v));
    ENSURE(VALID_REF(result->b));
    ENSURE(VALID_REF(result->x));
    ENSURE(VALID_REF(result->w));

    return result;
}

/*
*  Free the coefficient matrices and vectors needed for the Moving Least
*  Square derivatives.
*/
bool FreeMLSMatrices(MLSDeriv_pa *mlsDerivPtr)     /* holds deriv values */
{
    bool isOk = true;
    REQUIRE(VALID_REF(mlsDerivPtr));
    REQUIRE(*mlsDerivPtr == NULL || VALID_REF(*mlsDerivPtr));
    REQUIRE(IMPLICATION(*mlsDerivPtr != NULL, (*mlsDerivPtr)->a == NULL || VALID_REF((*mlsDerivPtr)->a)));
    REQUIRE(IMPLICATION(*mlsDerivPtr != NULL, (*mlsDerivPtr)->v == NULL || VALID_REF((*mlsDerivPtr)->v)));
    REQUIRE(IMPLICATION(*mlsDerivPtr != NULL, (*mlsDerivPtr)->b == NULL || VALID_REF((*mlsDerivPtr)->b)));
    REQUIRE(IMPLICATION(*mlsDerivPtr != NULL, (*mlsDerivPtr)->x == NULL || VALID_REF((*mlsDerivPtr)->x)));
    REQUIRE(IMPLICATION(*mlsDerivPtr != NULL, (*mlsDerivPtr)->w == NULL || VALID_REF((*mlsDerivPtr)->w)));

    if (*mlsDerivPtr != NULL)
    {
        if ((*mlsDerivPtr)->a != NULL)
        {
            nrcFreeMatrix((*mlsDerivPtr)->a, 1, (*mlsDerivPtr)->maxNumInfluencerPts,
                1, (*mlsDerivPtr)->maxNumBasisFunctions);
            (*mlsDerivPtr)->a = NULL;
        }
        if ((*mlsDerivPtr)->v != NULL)
        {
            nrcFreeMatrix((*mlsDerivPtr)->v, 1, (*mlsDerivPtr)->maxNumBasisFunctions,
                1, (*mlsDerivPtr)->maxNumBasisFunctions);
            (*mlsDerivPtr)->v = NULL;
        }

        if ((*mlsDerivPtr)->b != NULL)
        {
            nrcFreeVector((*mlsDerivPtr)->b, 1, (*mlsDerivPtr)->maxNumInfluencerPts);
            (*mlsDerivPtr)->b = NULL;
        }

        if ((*mlsDerivPtr)->x != NULL)
        {
            nrcFreeVector((*mlsDerivPtr)->x, 1, (*mlsDerivPtr)->maxNumBasisFunctions);
            (*mlsDerivPtr)->x = NULL;
        }

        if ((*mlsDerivPtr)->w != NULL)
        {
            nrcFreeVector((*mlsDerivPtr)->w, 1, (*mlsDerivPtr)->maxNumBasisFunctions);
            (*mlsDerivPtr)->w = NULL;
        }

        (*mlsDerivPtr)->maxNumBasisFunctions = 0;
        (*mlsDerivPtr)->maxNumInfluencerPts = 0;

        FREE_ITEM(*mlsDerivPtr);
        *mlsDerivPtr = NULL;
    }
    return isOk;
}

/*
*  Reallocate the coefficient matrices and vectors needed for the Moving Least
*  Square derivatives.
*/
bool reallocMLSMatrices(
    Index_t      numInfluencerPts,
    Index_t      numDerivTypes,
    MLSDeriv_pa  mlsDerivPtr) /* holds derivatives */
{
    bool isOk = true;
    REQUIRE(numInfluencerPts > 2);
    REQUIRE(numDerivTypes > 2);
    REQUIRE(VALID_REF(mlsDerivPtr));

    if (numDerivTypes > mlsDerivPtr->maxNumBasisFunctions)
    {
        nrcFreeMatrix(mlsDerivPtr->a, 1, mlsDerivPtr->maxNumInfluencerPts, 1, mlsDerivPtr->maxNumBasisFunctions);
        nrcFreeMatrix(mlsDerivPtr->v, 1, mlsDerivPtr->maxNumBasisFunctions, 1, mlsDerivPtr->maxNumBasisFunctions);

        nrcFreeVector(mlsDerivPtr->b, 1, mlsDerivPtr->maxNumInfluencerPts);
        nrcFreeVector(mlsDerivPtr->x, 1, mlsDerivPtr->maxNumBasisFunctions);
        nrcFreeVector(mlsDerivPtr->w, 1, mlsDerivPtr->maxNumBasisFunctions);

        mlsDerivPtr->maxNumBasisFunctions = numDerivTypes;
        if (mlsDerivPtr->maxNumInfluencerPts <= 0) mlsDerivPtr->maxNumInfluencerPts = 1;
        while (mlsDerivPtr->maxNumInfluencerPts < numInfluencerPts)
            mlsDerivPtr->maxNumInfluencerPts = 2 * mlsDerivPtr->maxNumInfluencerPts;

        mlsDerivPtr->a = nrcAllocMatrix(1, mlsDerivPtr->maxNumInfluencerPts, 1, mlsDerivPtr->maxNumBasisFunctions);
        mlsDerivPtr->v = nrcAllocMatrix(1, mlsDerivPtr->maxNumBasisFunctions, 1, mlsDerivPtr->maxNumBasisFunctions);
        mlsDerivPtr->b = nrcAllocVector(1, mlsDerivPtr->maxNumInfluencerPts);
        mlsDerivPtr->x = nrcAllocVector(1, mlsDerivPtr->maxNumBasisFunctions);
        mlsDerivPtr->w = nrcAllocVector(1, mlsDerivPtr->maxNumBasisFunctions);
    }

    if (numInfluencerPts > mlsDerivPtr->maxNumInfluencerPts)
    {
        nrcFreeMatrix(mlsDerivPtr->a, 1, mlsDerivPtr->maxNumInfluencerPts, 1, mlsDerivPtr->maxNumBasisFunctions);
        nrcFreeVector(mlsDerivPtr->b, 1, mlsDerivPtr->maxNumInfluencerPts);

        if (mlsDerivPtr->maxNumInfluencerPts <= 0) mlsDerivPtr->maxNumInfluencerPts = 1;
        while (mlsDerivPtr->maxNumInfluencerPts < numInfluencerPts)
            mlsDerivPtr->maxNumInfluencerPts = 2 * mlsDerivPtr->maxNumInfluencerPts;

        mlsDerivPtr->a = nrcAllocMatrix(1, mlsDerivPtr->maxNumInfluencerPts, 1, mlsDerivPtr->maxNumBasisFunctions);
        mlsDerivPtr->b = nrcAllocVector(1, mlsDerivPtr->maxNumInfluencerPts);
    }
    return isOk;
}



bool calculateGradientOnSphere(
    std::vector<Vec3> const&     nodalXYZUnnormalized, // in: nodal values of coordinate variables
    std::vector<double> const&   cellScalars, // in: cell-centered
    std::vector<TriNodes> const& surfaceTriangles, // in: cell-based
    Vec3 const&                  sphereCenter,// in: single XYZ
    double                       sphereRadius, // in: single scalar
    GradientCalcMethod_e         gradientCalcMethod, // in: number of pts to consider for gradient calculation
    std::vector<Vec3>&           gradients, // out: cell-centered dx,dy,dz
    char const*&                 statusMessage) // OUT: contains description of error or "Triangulation Successful" if successful
{
    bool result = false;

    CHECK(sphereRadius > SMALLFLOAT);
    double denormalizationFactor = 1.0 / sphereRadius;

    // check inputs
    size_t const numNodes = nodalXYZUnnormalized.size();
    size_t const numTriangles = surfaceTriangles.size();
    REQUIRE(cellScalars.size() == numTriangles);

    // check connectivity
    for (auto pos = 0; pos < numTriangles; ++pos)
    {
        REQUIRE(surfaceTriangles[pos].v1() < numNodes);
        REQUIRE(surfaceTriangles[pos].v2() < numNodes);
        REQUIRE(surfaceTriangles[pos].v3() < numNodes);
    }

    try
    {
        // normalize XYZs of triangulation
        std::vector<Vec3> normalizedNodalXYZs;
        normalizedNodalXYZs.reserve(numNodes);

        for (auto iter = nodalXYZUnnormalized.begin(); iter != nodalXYZUnnormalized.end(); ++iter)
            normalizedNodalXYZs.push_back((*iter - sphereCenter).normalize());

        Index_t const numBasisFunctions = 6; // 1, x, y, x^2, y^2, xy
        Index_t numInfluencerPts = 4;
        
        // for GradientCalcMethod_ThreeNodesAndCellCenter
        std::vector<double> nodalScalars;

        // for GradientCalcMethod_AllConnectedCellCenters
        typedef std::pair<Index_t, Index_t> NodeCellPair; // one cell number and one node number
        std::vector<NodeCellPair> nodeCellList;
        std::vector<Index_t> attachedCellList;

        if (gradientCalcMethod == GradientCalcMethod_ThreeNodesAndCellCenter)
        {
            // need scalar values at nodes 
            nodalScalars.resize(numNodes, 0.0);

            std::vector<Index_t> cellCountForNode;
            cellCountForNode.resize(numNodes, 0);

            for (Index_t cell = 0; cell < numTriangles; ++cell)
            {
                double const scalarVal = cellScalars[cell];

                Index_t const n1 = surfaceTriangles[cell].v1();
                Index_t const n2 = surfaceTriangles[cell].v2();
                Index_t const n3 = surfaceTriangles[cell].v3();

                nodalScalars[n1] += scalarVal;
                cellCountForNode[n1]++;

                nodalScalars[n2] += scalarVal;
                cellCountForNode[n2]++;

                nodalScalars[n3] += scalarVal;
                cellCountForNode[n3]++;
            }
            for (Index_t node = 0; node < numNodes; node++)
            {
                CHECK(cellCountForNode[node] > 0);
                nodalScalars[node] /= (double)cellCountForNode[node];
            }
        }
        else if ( gradientCalcMethod == GradientCalcMethod_AllConnectedCellCenters || 
                  gradientCalcMethod == GradientCalcMethod_FourCellCenters )
        {
            nodeCellList.reserve(3 * numTriangles);
            for (Index_t cell = 0; cell < numTriangles; ++cell)
            {
                Index_t const n1 = surfaceTriangles[cell].v1();
                Index_t const n2 = surfaceTriangles[cell].v2();
                Index_t const n3 = surfaceTriangles[cell].v3();
                nodeCellList.push_back(NodeCellPair(n1, cell));
                nodeCellList.push_back(NodeCellPair(n2, cell));
                nodeCellList.push_back(NodeCellPair(n3, cell));
            }
            std::sort(nodeCellList.begin(), nodeCellList.end());
            numInfluencerPts = 20; // so we don't keep reallocating, set to a decent sizes to start
            attachedCellList.reserve(40); 
        }
        else
            CHECK(false);

        // calculate the derivatives via MLS and stuff into gradients array
        std::vector<double> xvalues;
        std::vector<double> yvalues;
        std::vector<double> scalars;
        xvalues.reserve(40); // so we don't keep reallocating, set to a decent sizes to start
        yvalues.reserve(40);
        scalars.reserve(40);

        gradients.reserve(numTriangles);
        MLSDeriv_pa mlsDerivPtr = allocMLSMatrices(numInfluencerPts, numBasisFunctions);
        for (Index_t cell = 0; cell < numTriangles; ++cell)
        {
            Index_t const n1 = surfaceTriangles[cell].v1();
            Index_t const n2 = surfaceTriangles[cell].v2();
            Index_t const n3 = surfaceTriangles[cell].v3();
            Vec3 const xyzN1 = normalizedNodalXYZs[n1];
            Vec3 const xyzN2 = normalizedNodalXYZs[n2];
            Vec3 const xyzN3 = normalizedNodalXYZs[n3];
            Vec3 const xyzCC = (xyzN1 + xyzN2 + xyzN3).normalize();

            // we use xyzCC as the z'-axis in a new coordinate system
            static Vec3 xaxis(1.0, 0.0, 0.0);
            static Vec3 yaxis(0.0, 1.0, 0.0);
            static Vec3 zaxis(0.0, 0.0, 1.0);

            Vec3 const zprime = xyzCC;
            Vec3 xprime;
            double const threshold = 0.5774; // slightly bigger than sqrt(1/3)
            if (ABS(xyzCC.x()) < threshold)
                xprime = xyzCC.cross(xaxis);
            else if (ABS(xyzCC.y()) < threshold)
                xprime = xyzCC.cross(yaxis);
            else if (ABS(xyzCC.z()) < threshold)
                xprime = xyzCC.cross(zaxis);
            else
                CHECK(false);
            xprime.normalize();
            Vec3 const yprime = xprime.cross(zprime);
            yprime.normalize();
            CHECK(zprime.isNormalized());

            xvalues.clear();
            yvalues.clear();
            scalars.clear();

            xvalues.push_back(xprime.dot(xyzCC)); // rest of algorithm assumes this is the point where we are evaluating the derivative
            yvalues.push_back(yprime.dot(xyzCC));
            scalars.push_back(cellScalars[cell]);

            if (gradientCalcMethod == GradientCalcMethod_ThreeNodesAndCellCenter)
            {
                xvalues.push_back(xprime.dot(xyzN1));
                yvalues.push_back(yprime.dot(xyzN1));
                scalars.push_back(nodalScalars[n1]);

                xvalues.push_back(xprime.dot(xyzN2));
                yvalues.push_back(yprime.dot(xyzN2));
                scalars.push_back(nodalScalars[n2]);

                xvalues.push_back(xprime.dot(xyzN3));
                yvalues.push_back(yprime.dot(xyzN3));
                scalars.push_back(nodalScalars[n3]);
            }
            else if (gradientCalcMethod == GradientCalcMethod_AllConnectedCellCenters ||
                     gradientCalcMethod == GradientCalcMethod_FourCellCenters)
            {
                numInfluencerPts = 1;
                attachedCellList.clear();
                Index_t const cornerNodes[3] = { n1, n2, n3 };
                for (Index_t corner = 0; corner < 3; ++corner)
                {
                    Index_t cornerNode = cornerNodes[corner];
                    for (auto nodeCellListIter = std::lower_bound(nodeCellList.begin(), nodeCellList.end(), NodeCellPair(cornerNode, 0));
                         nodeCellListIter != nodeCellList.end() && nodeCellListIter->first == cornerNode;
                         ++nodeCellListIter)
                    {
                        Index_t const attachedCell = nodeCellListIter->second;
                        if (attachedCell != cell) // already included as first entry
                        {
                            if (gradientCalcMethod == GradientCalcMethod_FourCellCenters)
                            {
                                // push only if two node are the same, which is means an adjacent triangle
                                Index_t const a1 = surfaceTriangles[attachedCell].v1();
                                Index_t const a2 = surfaceTriangles[attachedCell].v2();
                                Index_t const a3 = surfaceTriangles[attachedCell].v3();
                                if (((a1 == n1 || a1 == n2 || a1 == n3) && (a2 == n1 || a2 == n2 || a2 == n3)) ||
                                    ((a2 == n1 || a2 == n2 || a2 == n3) && (a3 == n1 || a3 == n2 || a3 == n3)) ||
                                    ((a3 == n1 || a3 == n2 || a3 == n3) && (a1 == n1 || a1 == n2 || a1 == n3)))
                                {
                                    attachedCellList.push_back(attachedCell);
                                }
                            }
                            else
                            {
                                // use every triangle that includes one one of the base cell
                                attachedCellList.push_back(attachedCell);
                            }
                        }
                    }
                }

                std::sort(attachedCellList.begin(), attachedCellList.end());
                auto last = std::unique(attachedCellList.begin(), attachedCellList.end());
                attachedCellList.erase(last, attachedCellList.end());

                bool add = false;
                for (auto attachedCellListIter = attachedCellList.begin();
                    attachedCellListIter != attachedCellList.end();
                    ++attachedCellListIter)
                {
                    Index_t const attachedCell = *attachedCellListIter;

                    Index_t const acN1 = surfaceTriangles[attachedCell].v1(); // attached cell node 1
                    Index_t const acN2 = surfaceTriangles[attachedCell].v2();
                    Index_t const acN3 = surfaceTriangles[attachedCell].v3();
                    Vec3 const xyzACN1 = normalizedNodalXYZs[acN1]; // xyz of attached cell node 1
                    Vec3 const xyzACN2 = normalizedNodalXYZs[acN2];
                    Vec3 const xyzACN3 = normalizedNodalXYZs[acN3];
                    Vec3 const xyzACC = (xyzACN1 + xyzACN2 + xyzACN3).normalize(); // xyz of attached cell center

                    xvalues.push_back(xprime.dot(xyzACC));
                    yvalues.push_back(yprime.dot(xyzACC));
                    scalars.push_back(cellScalars[attachedCell]);
                    ++numInfluencerPts;
                }
            }
            else
                CHECK(false);

            CHECK(numInfluencerPts >= 3); // minimum to do a curve fit for first derivatives
            CHECK(xvalues.size() == numInfluencerPts);
            CHECK(yvalues.size() == numInfluencerPts);
            CHECK(scalars.size() == numInfluencerPts);

            if (mlsDerivPtr->maxNumBasisFunctions < numBasisFunctions ||
                mlsDerivPtr->maxNumInfluencerPts < numInfluencerPts)
            {
                if (!reallocMLSMatrices(numInfluencerPts, numBasisFunctions, mlsDerivPtr))
                    throw std::bad_alloc();
            }

            double const x0 = xvalues[0];
            double const y0 = yvalues[0];
            CHECK(ABS(x0) < 1e-15);
            CHECK(ABS(y0) < 1e-15);

            double** a = mlsDerivPtr->a;
            double*  b = mlsDerivPtr->b;
            double r2sum = 0.0;
            double r2min = LARGEFLOAT;
            for (Index_t ii = 0; ii < numInfluencerPts; ii++)
            {
                double const xn = xvalues[ii];
                double const yn = yvalues[ii];

                double const dx = xn - x0;
                double const dy = yn - y0;

                a[ii + 1][1] = 1;
                a[ii + 1][2] = dx;
                a[ii + 1][3] = dy;
                a[ii + 1][4] = 0.5 * dx * dx;
                a[ii + 1][5] = 0.5 * dy * dy;
                a[ii + 1][6] = dx * dy;

                b[ii + 1] = scalars[ii];

                double const r2 = dx * dx + dy * dy;
                r2sum += r2;
                if (ii != 0)
                    r2min = MIN(r2, r2min);
            }

            /*
            *  Scale the rows to by the weighting parameter 1/(r**2/r2ave + 0.1).
            *  Use the fact that dx and dy are stored in the a matrix.
            */
            double const r2ave = r2sum / (double)numInfluencerPts;
            for (Index_t ii = 0; ii < numInfluencerPts; ii++)
            {
                double const dx = a[ii + 1][2];
                double const dy = a[ii + 1][3];

                double const r2 = dx * dx + dy * dy;
                double const weight = r2ave / (r2 + 0.1 * r2ave);

                a[ii + 1][1] *= weight;
                a[ii + 1][2] *= weight;
                a[ii + 1][3] *= weight;
                a[ii + 1][4] *= weight;
                a[ii + 1][5] *= weight;
                a[ii + 1][6] *= weight;

                b[ii + 1] *= weight;
            }

            double* w = mlsDerivPtr->w;

            /* Compute the singular value decomposition. */
            nrcComputeSVD(mlsDerivPtr->a, numInfluencerPts, numBasisFunctions, mlsDerivPtr->w, mlsDerivPtr->v);

            /* If singular values are less than a certain level, set them to zero. */
            double wmax = 0.0;
            for (Index_t ii = 1; ii <= 6; ii++)
                if (wmax < w[ii])
                    wmax = w[ii];
            double const wmin = 1.0e-6 * wmax;

            mlsDerivPtr->IsSingular = false;
            for (Index_t ii = 1; ii <= 6; ii++)
            {
                if (w[ii] < wmin)
                {
                    w[ii] = 0.0;
                    mlsDerivPtr->IsSingular = true;
                }
            }

            /*
            * Solve the system of equations for the derivatives (coefficients of
            * the curve fit).
            */
            nrcBackSubstituteSVD(mlsDerivPtr->a, mlsDerivPtr->w, mlsDerivPtr->v,
                numInfluencerPts, numBasisFunctions,
                mlsDerivPtr->b, mlsDerivPtr->x);

            double* x = mlsDerivPtr->x;
            double const dPdu = x[2];
            double const dPdv = x[3];
            // dPdw is 0.0, no change in direction normal to sphere

            /*
            *
            mlsDerivPtr->ddx = x[2];
            mlsDerivPtr->ddy = x[3];
            *
            * don't care about these
            mlsDerivPtr->d2dxx = x[4];
            mlsDerivPtr->d2dyy = x[5];
            mlsDerivPtr->d2dxy = x[6];
            */

            // get back to XYZ space from x',y',z' space
            double const ddx = dPdu * xprime.x() + dPdv * yprime.x() /*+ 0.0 * zprime.x()*/;
            double const ddy = dPdu * xprime.y() + dPdv * yprime.y() /*+ 0.0 * zprime.y()*/;
            double const ddz = dPdu * xprime.z() + dPdv * yprime.z() /*+ 0.0 * zprime.z()*/;

            // denormalize back into the sphere we are calculating from
            gradients.push_back(Vec3(ddx*denormalizationFactor, ddy*denormalizationFactor, ddz*denormalizationFactor));
        }


        result = true;
        statusMessage = "Calculation successful";
    }
    catch (char const* message)
    {
        statusMessage = message;
        result = false;
    }
    return result;
}

}
