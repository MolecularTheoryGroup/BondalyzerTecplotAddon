#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "TASSERT.h"
#include "GUIDEFS.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "LINEARALGEBRA.h"
// #include "ODERUNGEKUTTA.h"  // TODO: Is this OK?
#include "CRITPOINTS.h"
#include "ENGINE.h"
#include <string.h>


#define MAXNEWTONITER 100

static double TetraTolerance = 1.0005;
// static double RangeEpsilon = 0.0;



/*
 * Compute the interpolation weights for a trilinear brick.
 * Input:
 *    r,s,t: Natural coordinates - ranges is from -1 to 1 and cell
 *           center is at 0.
 * Output:
 *    W[]:   Averaging weights for cell (duplicate nodes removed).
 *           At this stage, only the weights for the face are set.
 *    Returns: True if point (r,s,t) is in the brick.
 */
Boolean_t BrickTrilinearWeight(double       r,
                               double       s,
                               double       t,
                               double       W[8])
{
    Boolean_t CellFound = FALSE;
    REQUIRE(VALID_REF(W));

    if (r > -TetraTolerance && r < TetraTolerance &&
        s > -TetraTolerance && s < TetraTolerance &&
        t > -TetraTolerance && t < TetraTolerance)
    {
        double    OneMinusR, OnePlusR;
        double    OneMinusS, OnePlusS;
        double    OneMinusT, OnePlusT;

        CellFound = TRUE;
        OneMinusR = 1.0 - r;
        OnePlusR  = 1.0 + r;
        OneMinusS = 1.0 - s;
        OnePlusS  = 1.0 + s;
        OneMinusT = 1.0 - t;
        OnePlusT  = 1.0 + t;
        W[0] = 0.125 * OneMinusR * OneMinusS * OneMinusT;
        W[1] = 0.125 * OnePlusR  * OneMinusS * OneMinusT;
        W[2] = 0.125 * OnePlusR  * OnePlusS  * OneMinusT;
        W[3] = 0.125 * OneMinusR * OnePlusS  * OneMinusT;
        W[4] = 0.125 * OneMinusR * OneMinusS * OnePlusT;
        W[5] = 0.125 * OnePlusR  * OneMinusS * OnePlusT;
        W[6] = 0.125 * OnePlusR  * OnePlusS  * OnePlusT;
        W[7] = 0.125 * OneMinusR * OnePlusS  * OnePlusT;
        // TODO delete
        // W[0] = 0.125 * OneMinusR * OneMinusS * OneMinusT;
        // W[1] = 0.125 * OnePlusR  * OneMinusS * OneMinusT;
        // W[2] = 0.125 * OneMinusR * OnePlusS  * OneMinusT;
        // W[3] = 0.125 * OnePlusR  * OnePlusS  * OneMinusT;
        // W[4] = 0.125 * OneMinusR * OneMinusS * OnePlusT;
        // W[5] = 0.125 * OnePlusR  * OneMinusS * OnePlusT;
        // W[6] = 0.125 * OneMinusR * OnePlusS  * OnePlusT;
        // W[7] = 0.125 * OnePlusR  * OnePlusS  * OnePlusT;
    }

    ENSURE(VALID_BOOLEAN(CellFound));
    return (CellFound);
}


/*
 * Compute the derivatives wrt the natural coordinates of the
 * interpolation weights for a trilinear brick.
 * Input:
 *    r,s,t: Natural coordinates - ranges is from -1 to 1 and cell
 *           center is at 0.
 * Output:
 *    W[]:    Interpolation weights.
 *    dwdr[]: Derivative of weights wrt r (first natural coordinate).
 *    dwds[]: Derivative of weights wrt s (second natural coordinate).
 *    dwdt[]: Derivative of weights wrt t (third natural coordinate).
 *
 *    Returns: True if point (r,s,t) is in the brick.
 */
Boolean_t BrickTrilinearWeightDerivatives(double       r,
                                          double       s,
                                          double       t,
                                          double       W[8],
                                          double       dwdr[8],
                                          double       dwds[8],
                                          double       dwdt[8])
{
    Boolean_t CellFound = FALSE;
    REQUIRE(VALID_REF(W));
    REQUIRE(VALID_REF(dwdr));
    REQUIRE(VALID_REF(dwds));
    REQUIRE(VALID_REF(dwdt));

    if (r > -TetraTolerance && r < TetraTolerance &&
        s > -TetraTolerance && s < TetraTolerance &&
        t > -TetraTolerance && t < TetraTolerance) CellFound = TRUE;

    {
        double    OneMinusR, OnePlusR;
        double    OneMinusS, OnePlusS;
        double    OneMinusT, OnePlusT;

        OneMinusR = 1.0 - r;
        OnePlusR  = 1.0 + r;
        OneMinusS = 1.0 - s;
        OnePlusS  = 1.0 + s;
        OneMinusT = 1.0 - t;
        OnePlusT  = 1.0 + t;
        W[0] = 0.125 * OneMinusR * OneMinusS * OneMinusT;
        W[1] = 0.125 * OnePlusR  * OneMinusS * OneMinusT;
        W[2] = 0.125 * OnePlusR  * OnePlusS  * OneMinusT;
        W[3] = 0.125 * OneMinusR * OnePlusS  * OneMinusT;
        W[4] = 0.125 * OneMinusR * OneMinusS * OnePlusT;
        W[5] = 0.125 * OnePlusR  * OneMinusS * OnePlusT;
        W[6] = 0.125 * OnePlusR  * OnePlusS  * OnePlusT;
        W[7] = 0.125 * OneMinusR * OnePlusS  * OnePlusT;
        dwdr[0] = -0.125 * OneMinusS * OneMinusT;
        dwdr[1] =  0.125 * OneMinusS * OneMinusT;
        dwdr[2] =  0.125 * OnePlusS  * OneMinusT;
        dwdr[3] = -0.125 * OnePlusS  * OneMinusT;
        dwdr[4] = -0.125 * OneMinusS * OnePlusT;
        dwdr[5] =  0.125 * OneMinusS * OnePlusT;
        dwdr[6] =  0.125 * OnePlusS  * OnePlusT;
        dwdr[7] = -0.125 * OnePlusS  * OnePlusT;
        dwds[0] = -0.125 * OneMinusR * OneMinusT;
        dwds[1] = -0.125 * OnePlusR  * OneMinusT;
        dwds[2] =  0.125 * OnePlusR  * OneMinusT;
        dwds[3] =  0.125 * OneMinusR * OneMinusT;
        dwds[4] = -0.125 * OneMinusR * OnePlusT;
        dwds[5] = -0.125 * OnePlusR  * OnePlusT;
        dwds[6] =  0.125 * OnePlusR  * OnePlusT;
        dwds[7] =  0.125 * OneMinusR * OnePlusT;
        dwdt[0] = -0.125 * OneMinusR * OneMinusS;
        dwdt[1] = -0.125 * OnePlusR  * OneMinusS;
        dwdt[2] = -0.125 * OnePlusR  * OnePlusS;
        dwdt[3] = -0.125 * OneMinusR * OnePlusS;
        dwdt[4] =  0.125 * OneMinusR * OneMinusS;
        dwdt[5] =  0.125 * OnePlusR  * OneMinusS;
        dwdt[6] =  0.125 * OnePlusR  * OnePlusS;
        dwdt[7] =  0.125 * OneMinusR * OnePlusS;
    }

    ENSURE(VALID_BOOLEAN(CellFound));
    return (CellFound);
}



/*
 * Compute the trilinear natural coordinates (r,s,t) corresponding
 * to a point(X,Y,Z) in the cell.
 * Input:
 *    X,Y,Z: Physical coordinates of a point (perhaps) in the cell.
 * Output:
 *    r,s,t: Natural coordinates - ranges is from -1 to 1 and cell
 *           center is at 0.
 *    W[]:    Interpolation weights.
 *
 *    Returns: True if point (X,Y,Z) is in the trilinear brick.
 */
Boolean_t TrilinearBrickNaturalCoord(double       X,
                                     double       Y,
                                     double       Z,
                                     double       XCell[8],
                                     double       YCell[8],
                                     double       ZCell[8],
                                     double       XCMin,
                                     double       YCMin,
                                     double       ZCMin,
                                     double       XCMax,
                                     double       YCMax,
                                     double       ZCMax,
                                     double       *r,
                                     double       *s,
                                     double       *t,
                                     double       W[8])
{
    Boolean_t CellFound = FALSE;
    double dwdr[8];
    double dwds[8];
    double dwdt[8];
    REQUIRE(VALID_REF(W));
    REQUIRE(VALID_REF(XCell));
    REQUIRE(VALID_REF(YCell));
    REQUIRE(VALID_REF(ZCell));
    REQUIRE(VALID_REF(r));
    REQUIRE(VALID_REF(s));
    REQUIRE(VALID_REF(t));

    /*
     * Use a Newton iteration to solve the non-linear equations for
     * r,s,t in terms of X,Y,Z.
     */
    if (X <= XCMax || X >= XCMin || Y <= YCMax || Y >= YCMin || Z <= ZCMax || Z >= ZCMin)
    {
        Boolean_t IsOk = TRUE;
        Boolean_t IsInside;
        int       ic;
        int       iter = 0;
        double    Xm, Ym, Zm, DeltaX, DeltaY, DeltaZ, Toler, TolerSqr, DeltaSqr;

        /* If computing derivatives of weights, must pass through while loop once */
        Boolean_t FirstPass = TRUE;

        /* Matrix and temporary array for LU decomposition */
        LgIndex_t indx[3];
        double    TmpSpace[3];
        double    DeltaXYZ[3];
        Matrix_s  Jacobian;
        Jacobian = CreateMatrix(3, 3);

        Toler = 0.1 * (TetraTolerance - 1.0) * (XCMax - XCMin + YCMax - YCMin + ZCMax - ZCMin);
        TolerSqr = Toler * Toler;

        /* Set initial values of r,s,t - it out of range. */
        if (*r < -1.0 || *r > 1.0) *r = 0.0;
        if (*s < -1.0 || *s > 1.0) *s = 0.0;
        if (*t < -1.0 || *t > 1.0) *t = 0.0;


        /* Compute error. */
        IsInside = BrickTrilinearWeightDerivatives(*r, *s, *t, W, dwdr, dwds, dwdt);

        Xm = 0.0;
        Ym = 0.0;
        Zm = 0.0;
        for (ic = 0; ic < 8; ic++)
        {
            Xm = Xm + W[ic] * XCell[ic];
            Ym = Ym + W[ic] * YCell[ic];
            Zm = Zm + W[ic] * ZCell[ic];
        }
        DeltaX = X - Xm;
        DeltaY = Y - Ym;
        DeltaZ = Z - Zm;
        DeltaSqr = DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ;


        /* Iterate */
        while ((IsOk && DeltaSqr > TolerSqr && iter < MAXNEWTONITER) ||
               (FirstPass))
        {
            double dxdr = 0.0;
            double dxds = 0.0;
            double dxdt = 0.0;
            double dydr = 0.0;
            double dyds = 0.0;
            double dydt = 0.0;
            double dzdr = 0.0;
            double dzds = 0.0;
            double dzdt = 0.0;

            iter++;

            /* Compute Jacobian Matrix */
            for (ic = 0; ic < 8; ic++)
            {
                dxdr = dxdr + dwdr[ic] * XCell[ic];
                dxds = dxds + dwds[ic] * XCell[ic];
                dxdt = dxdt + dwdt[ic] * XCell[ic];
                dydr = dydr + dwdr[ic] * YCell[ic];
                dyds = dyds + dwds[ic] * YCell[ic];
                dydt = dydt + dwdt[ic] * YCell[ic];
                dzdr = dzdr + dwdr[ic] * ZCell[ic];
                dzds = dzds + dwds[ic] * ZCell[ic];
                dzdt = dzdt + dwdt[ic] * ZCell[ic];
            }

            Jacobian.val[0][0] = dxdr;
            Jacobian.val[0][1] = dxds;
            Jacobian.val[0][2] = dxdt;
            Jacobian.val[1][0] = dydr;
            Jacobian.val[1][1] = dyds;
            Jacobian.val[1][2] = dydt;
            Jacobian.val[2][0] = dzdr;
            Jacobian.val[2][1] = dzds;
            Jacobian.val[2][2] = dzdt;

            /* Solve linear system (Jacobian) (DeltaRST) = (DeltaXYZ) for DeltaRST */
            IsOk = cmludcmp(Jacobian, 3, indx, TmpSpace);

            DeltaXYZ[0] = DeltaX;
            DeltaXYZ[1] = DeltaY;
            DeltaXYZ[2] = DeltaZ;
            if (IsOk) IsOk = cmlubksb(Jacobian, 3, indx, DeltaXYZ);
            if (IsOk)
            {
                *r = *r + 0.5 * DeltaXYZ[0];
                *s = *s + 0.5 * DeltaXYZ[1];
                *t = *t + 0.5 * DeltaXYZ[2];
            }

            /* Recompute DeltaX,DeltaY,DeltaZ for new r,s,t */
            IsInside = BrickTrilinearWeightDerivatives(*r, *s, *t, W, dwdr, dwds, dwdt);
            Xm = 0.0;
            Ym = 0.0;
            Zm = 0.0;
            for (ic = 0; ic < 8; ic++)
            {
                Xm = Xm + W[ic] * XCell[ic];
                Ym = Ym + W[ic] * YCell[ic];
                Zm = Zm + W[ic] * ZCell[ic];
            }
            DeltaX = X - Xm;
            DeltaY = Y - Ym;
            DeltaZ = Z - Zm;
            DeltaSqr = DeltaX * DeltaX + DeltaY * DeltaY + DeltaZ * DeltaZ;

            /* Been through once, set FirstPass = FALSE */
            FirstPass = FALSE;
        }

        /* Is it inside? */
        if (*r < TetraTolerance && *r > -TetraTolerance &&
            *s < TetraTolerance && *s > -TetraTolerance &&
            *t < TetraTolerance && *t > -TetraTolerance && IsOk)
        {
            CellFound = TRUE;
        }

        /* Clean up */
        FreeMatrix(&Jacobian);
    }


    ENSURE(VALID_BOOLEAN(CellFound));
    return (CellFound);
}




/*
 * Use Jacobi transformations to compute the eigen system of a
 * symmetric 3x3 matrix. The input matrix is InMat[][]. The eigenvalues
 * are stored in array EgnVal[], and the right eigen vectors are stored
 * in array EgnVec[][]. The number of Jacobi rotations are stored in
 * NumRot.
 */
#define ROTATE(a,it,jt,kt,lt) g=a[it][jt];h=a[kt][lt];a[it][jt]=g-s*(h+g*Tau);\
a[kt][lt]=h+s*(g-h*Tau);

#define MAXITERS 50
Boolean_t Jacobi3by3(float **InMat,
                     float *EgnVal,
                     float **EgnVec,
                     int    *NumRot)
{
    Boolean_t IsOk = TRUE;
    Boolean_t IsDone = FALSE;
    int i, j, jj, Iter;
    float b[3], z[3], g;

    REQUIRE(VALID_REF(InMat));
    REQUIRE(VALID_REF(EgnVal));
    REQUIRE(VALID_REF(EgnVec));
    REQUIRE(VALID_REF(NumRot));

    /*
     * Initialize NumRot to zero, EigenVector matrix to identity, and
     * EigenValue array to diagonal of the input matrix.
     */
    *NumRot = 0;
    for (i = 0; i < 3 && IsOk; i++)
    {
        for (j = 0; j < 3; j++) EgnVec[i][j] = 0.0;
        EgnVec[i][i] = 1.0;
        EgnVal[i]    = InMat[i][i];
        b[i]         = EgnVal[i];
        z[i]         = 0.0;
    }

    /* Iteration - do up to MAXITERS Jacobi rotations */
    for (Iter = 0; Iter < MAXITERS && IsOk && !IsDone; Iter++)
    {
        float SumOffDiag, Threshold;

        /* Test to see if off-diagonal terms are zero */
        SumOffDiag = (float)(fabs(InMat[0][1]) + fabs(InMat[0][2]) + fabs(InMat[1][2]));
        if (SumOffDiag == 0.0) IsDone = TRUE;

        // Threshold = (float)(0.025 * SumOffDiag);  /* Equation 11.1.25 in Numerical Recipes */
        Threshold = (float)(0.00025 * SumOffDiag);  /* Equation 11.1.25 in Numerical Recipes */

        /* Rotate to zero each diagonal a[i][j] */
        for (i = 0; i < 2 && !IsDone; i++)
        {
            for (j = i + 1; j < 3; j++)
            {
                float Theta, t, c, s, Tau, h;
                /*
                 * After 2 iterations, set the off-diagonal terms
                 * to zero if it is 2-orders of magnitude smaller
                 * than the precision of either corresponding diagonal
                 */
                // float Err = (float)(100.0*fabs(InMat[i][j]));
                float Err = (float)(10000.0 * fabs(InMat[i][j]));
                if (Iter > 2
                    && (float)(fabs(EgnVal[i]) + Err) == (float)fabs(EgnVal[i])
                    && (float)(fabs(EgnVal[j]) + Err) == (float)fabs(EgnVal[j]))
                    InMat[i][j] = 0.0;
                else if (fabs(InMat[i][j]) > Threshold)
                {
                    float DiagDiff = EgnVal[j] - EgnVal[i];
                    if ((float)(fabs(DiagDiff) + Err) == (float)(fabs(DiagDiff)))
                        t = (InMat[i][j]) / DiagDiff;
                    else
                    {
                        Theta = (float)(0.5 * DiagDiff / (InMat[i][j]));  /* Equation 11.1.10 in NRC */
                        t = (float)(1.0 / (fabs(Theta) + sqrt(1.0 + Theta * Theta)));
                        if (Theta < 0.0) t = -t;
                    }

                    c = (float)(1.0 / sqrt(1 + t * t));
                    s = t * c;
                    Tau = (float)(s / (1.0 + c));
                    h = t * InMat[i][j];
                    z[i] -= h;
                    z[j] += h;
                    EgnVal[i] -= h;
                    EgnVal[j] += h;
                    InMat[i][j] = 0.0;
                    for (jj = 0; jj <= i - 1; jj++)
                    {
                        ROTATE(InMat, jj, i, jj, j)
                    }
                    for (jj = i + 1; jj <= j - 1; jj++)
                    {
                        ROTATE(InMat, i, jj, jj, j)
                    }
                    for (jj = j + 1; jj < 3; jj++)
                    {
                        ROTATE(InMat, i, jj, j, jj)
                    }
                    for (jj = 0; jj < 3; jj++)
                    {
                        ROTATE(EgnVec, jj, i, jj, j)
                        // g=EgnVec[jj][i];
                        // h=EgnVec[jj][j];
                        // EgnVec[jj][i]=g-s*(h+g*Tau);
                        // EgnVec[jj][j]=h+s*(g-h*Tau);
                    }
                    ++(*NumRot);
                }
            }
        }
        for (i = 0; i < 3 && !IsDone; i++)
        {
            b[i] += z[i];
            EgnVal[i] = b[i];
            z[i] = 0.0;
        }
    }
    if (!IsDone) IsOk = FALSE;

    ENSURE(VALID_BOOLEAN(IsOk));

    return (IsOk);
}


Boolean_t ComputeCurvature(EntIndex_t   ZoneNum,
                           LgIndex_t    i,
                           LgIndex_t    j,
                           LgIndex_t    k,
                           FieldData_pa CDVarFDPtr,
                           double      *d2dx2,
                           double      *d2dxdy,
                           double      *d2dxdz,
                           double      *d2dy2,
                           double      *d2dydz,
                           double      *d2dz2)
{
    Boolean_t  IsOk = TRUE;
    LgIndex_t  IMax, JMax, KMax;
    EntIndex_t NumZones, NumVars;
    LgIndex_t  IJMax;

    REQUIRE(TecUtilDataSetIsAvailable());
    REQUIRE(VALID_REF(CDVarFDPtr));


    /* Get num zones in dataset */
    if (IsOk)
    {
        char *DatasetTitle = NULL;
        TecUtilDataSetGetInfo(&DatasetTitle, &NumZones, &NumVars);
        TecUtilStringDealloc(&DatasetTitle);
    }

    REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
    REQUIRE(TecUtilZoneIsOrdered(ZoneNum));

    TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL, NULL);


    REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */
    REQUIRE(i > 0 && i <= IMax);
    REQUIRE(j > 0 && j <= JMax);
    REQUIRE(k > 0 && k <= KMax);

    IJMax = IMax * JMax;

    /* d2dx2 */
    if (IsOk)
    {
        LgIndex_t IndxI, IndxIp1, IndxIm1;
        LgIndex_t ii  = i;
        if (i == 1) ii  = 2;
        if (i == IMax) ii  = IMax - 1;

        IndxI = ii + (j - 1) * IMax + (k - 1) * IJMax;
        IndxIp1 = IndxI + 1;
        IndxIm1 = IndxI - 1;

        *d2dx2 =       TecUtilDataValueGetByRef(CDVarFDPtr, IndxIp1)
                       - 2.0 * TecUtilDataValueGetByRef(CDVarFDPtr, IndxI)
                       + TecUtilDataValueGetByRef(CDVarFDPtr, IndxIm1);
    }

    /* d2dy2 */
    if (IsOk)
    {
        LgIndex_t IndxJ, IndxJp1, IndxJm1;
        LgIndex_t jj  = j;
        if (j == 1) jj  = 2;
        if (j == JMax) jj  = JMax - 1;

        IndxJ = i + (jj - 1) * IMax + (k - 1) * IJMax;
        IndxJp1 = IndxJ + IMax;
        IndxJm1 = IndxJ - IMax;

        *d2dy2 =       TecUtilDataValueGetByRef(CDVarFDPtr, IndxJp1)
                       - 2.0 * TecUtilDataValueGetByRef(CDVarFDPtr, IndxJ)
                       + TecUtilDataValueGetByRef(CDVarFDPtr, IndxJm1);
    }

    /* d2dz2 */
    if (IsOk)
    {
        LgIndex_t IndxK, IndxKp1, IndxKm1;
        LgIndex_t kk  = k;
        if (k == 1) kk  = 2;
        if (k == KMax) kk  = KMax - 1;

        IndxK = i + (j - 1) * IMax + (kk - 1) * IJMax;
        IndxKp1 = IndxK + IJMax;
        IndxKm1 = IndxK - IJMax;

        *d2dz2 =       TecUtilDataValueGetByRef(CDVarFDPtr, IndxKp1)
                       - 2.0 * TecUtilDataValueGetByRef(CDVarFDPtr, IndxK)
                       + TecUtilDataValueGetByRef(CDVarFDPtr, IndxKm1);
    }

    /* d2dxdy, d2dxdz, and d2dydz */
    if (IsOk)
    {
        LgIndex_t ip1 = MIN(i + 1, IMax);
        LgIndex_t im1 = MAX(i - 1, 1);
        LgIndex_t jp1 = MIN(j + 1, JMax);
        LgIndex_t jm1 = MAX(j - 1, 1);
        LgIndex_t kp1 = MIN(k + 1, KMax);
        LgIndex_t km1 = MAX(k - 1, 1);

        *d2dxdy = (TecUtilDataValueGetByRef(CDVarFDPtr, ip1 + (jp1 - 1) * IMax + (k - 1) * IJMax)
                   - TecUtilDataValueGetByRef(CDVarFDPtr, ip1 + (jm1 - 1) * IMax + (k - 1) * IJMax)
                   - TecUtilDataValueGetByRef(CDVarFDPtr, im1 + (jp1 - 1) * IMax + (k - 1) * IJMax)
                   + TecUtilDataValueGetByRef(CDVarFDPtr, im1 + (jm1 - 1) * IMax + (k - 1) * IJMax))
                  / (float)((ip1 - im1) * (jp1 - jm1)) ;

        *d2dxdz = (TecUtilDataValueGetByRef(CDVarFDPtr, ip1 + (j - 1) * IMax + (kp1 - 1) * IJMax)
                   - TecUtilDataValueGetByRef(CDVarFDPtr, ip1 + (j - 1) * IMax + (km1 - 1) * IJMax)
                   - TecUtilDataValueGetByRef(CDVarFDPtr, im1 + (j - 1) * IMax + (kp1 - 1) * IJMax)
                   + TecUtilDataValueGetByRef(CDVarFDPtr, im1 + (j - 1) * IMax + (km1 - 1) * IJMax))
                  / (float)((ip1 - im1) * (kp1 - km1)) ;


        *d2dydz = (TecUtilDataValueGetByRef(CDVarFDPtr, i + (jp1 - 1) * IMax + (kp1 - 1) * IJMax)
                   - TecUtilDataValueGetByRef(CDVarFDPtr, i + (jp1 - 1) * IMax + (km1 - 1) * IJMax)
                   - TecUtilDataValueGetByRef(CDVarFDPtr, i + (jm1 - 1) * IMax + (kp1 - 1) * IJMax)
                   + TecUtilDataValueGetByRef(CDVarFDPtr, i + (jm1 - 1) * IMax + (km1 - 1) * IJMax))
                  / (float)((jp1 - jm1) * (kp1 - km1)) ;
    }

    return(IsOk);
}







Boolean_t CriticalPointInCell(EntIndex_t  ZoneNum,
                              EntIndex_t  UVarNum,
                              EntIndex_t  VVarNum,
                              EntIndex_t  WVarNum,
                              EntIndex_t  ChrgDensVarNum,
                              EntIndex_t  TypeVarNum,
                              LgIndex_t   IIndex,
                              LgIndex_t   JIndex,
                              LgIndex_t   KIndex,
                              double     *XCrtPt,
                              double     *YCrtPt,
                              double     *ZCrtPt,
                              double     *ChrgDensCrtPt,
                              EntIndex_t *TypeCrtPt,
                              float      *PrincDirX,
                              float      *PrincDirY,
                              float      *PrincDirZ)
{
    Boolean_t IsOk = TRUE;
    Boolean_t CellFound = FALSE;
    LgIndex_t IMax, JMax, KMax;
    double    UCell[8], VCell[8], WCell[8];  /* Velocities at cell corners */
    LgIndex_t Index[8];
    double    UCMin, UCMax, VCMin, VCMax, WCMin, WCMax;  /* Min and Max vel for cell */

    double    W[8];     /* Trilinear weight for critical point */
    double    r, s, t;  /* Cell natural coordinates of critical point */
    EntIndex_t NumZones, NumVars;
    FieldData_pa XVarFDPtr, YVarFDPtr, ZVarFDPtr;
    FieldData_pa UVarFDPtr = NULL;
    FieldData_pa VVarFDPtr = NULL;
    FieldData_pa WVarFDPtr = NULL;
    FieldData_pa CDVarFDPtr = NULL;
    FieldData_pa TypeVarFDPtr = NULL;

    int ic;

    REQUIRE(TecUtilDataSetIsAvailable());
    REQUIRE(VALID_REF(XCrtPt));
    REQUIRE(VALID_REF(YCrtPt));
    REQUIRE(VALID_REF(ZCrtPt));
    REQUIRE(VALID_REF(ChrgDensCrtPt));
    REQUIRE(VALID_REF(TypeCrtPt));
    REQUIRE(VALID_REF(PrincDirX));
    REQUIRE(VALID_REF(PrincDirY));
    REQUIRE(VALID_REF(PrincDirZ));


    /* Get num zones in dataset */
    if (IsOk)
    {
        char *DatasetTitle = NULL;
        TecUtilDataSetGetInfo(&DatasetTitle, &NumZones, &NumVars);
        TecUtilStringDealloc(&DatasetTitle);
    }

    REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
    REQUIRE(TecUtilZoneIsOrdered(ZoneNum));
    REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
    REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
    REQUIRE(WVarNum > 0 && WVarNum <= NumVars);
    REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
    REQUIRE(TypeVarNum > 0 && TypeVarNum <= NumVars);

    TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
                       NULL, NULL, NULL, NULL, NULL, NULL, NULL);


    REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */
    REQUIRE(IIndex > 0 && IIndex <= IMax);
    REQUIRE(JIndex > 0 && JIndex <= JMax);
    REQUIRE(KIndex > 0 && KIndex <= KMax);
    REQUIRE(VALID_REF(XVarFDPtr));
    REQUIRE(VALID_REF(YVarFDPtr));
    REQUIRE(VALID_REF(ZVarFDPtr));

    UVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, UVarNum);
    VVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, VVarNum);
    WVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, WVarNum);
    CDVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, ChrgDensVarNum);
    TypeVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, TypeVarNum);
    if (UVarFDPtr == NULL || VVarFDPtr == NULL || WVarFDPtr == NULL || CDVarFDPtr == NULL ||
        TypeVarFDPtr == NULL)
        IsOk = FALSE;




    /*
     * FE Brick and ZoneType_Ordered Data:
     *                                                            *
     *    7         6                                             *
     *    +---------+                                             *
     *   /|        /|                                             *
     * 4/ |      5/ |                                             *
     * +---------+  |                                             *
     * |  +------|--+                                             *
     * | /3      | /2                                             *
     * |/        |/                                               *
     * +---------+                                                *
     * 0         1                                                *
     */

    /* Set Velocities at corner of cell & compute vel min/max */
    if (IsOk)
    {

        Index[0] = IIndex + (JIndex - 1) * IMax + (KIndex - 1) * IMax * JMax;
        Index[1] = IIndex + 1 + (JIndex - 1) * IMax + (KIndex - 1) * IMax * JMax;
        Index[2] = IIndex + 1 + JIndex * IMax + (KIndex - 1) * IMax * JMax;
        Index[3] = IIndex + JIndex * IMax + (KIndex - 1) * IMax * JMax;
        Index[4] = IIndex + (JIndex - 1) * IMax + KIndex * IMax * JMax;
        Index[5] = IIndex + 1 + (JIndex - 1) * IMax + KIndex * IMax * JMax;
        Index[6] = IIndex + 1 + JIndex * IMax + KIndex * IMax * JMax;
        Index[7] = IIndex + JIndex * IMax + KIndex * IMax * JMax;

        for (ic = 0; ic < 8; ic++)
        {
            UCell[ic] = TecUtilDataValueGetByRef(UVarFDPtr, Index[ic]);
            VCell[ic] = TecUtilDataValueGetByRef(VVarFDPtr, Index[ic]);
            WCell[ic] = TecUtilDataValueGetByRef(WVarFDPtr, Index[ic]);
        }

        /* Compute the minima and maxima of gradients in the cell */
        UCMin = UCell[0];
        VCMin = VCell[0];
        WCMin = WCell[0];
        UCMax = UCell[0];
        VCMax = VCell[0];
        WCMax = WCell[0];
        for (ic = 1; ic < 8; ic++)
        {
            UCMin = MIN(UCMin, UCell[ic]);
            UCMax = MAX(UCMax, UCell[ic]);
            VCMin = MIN(VCMin, VCell[ic]);
            VCMax = MAX(VCMax, VCell[ic]);
            WCMin = MIN(WCMin, WCell[ic]);
            WCMax = MAX(WCMax, WCell[ic]);
        }

        if (UCMin <= 0.0 && UCMax > 0.0 &&
            VCMin <= 0.0 && VCMax > 0.0 &&
            WCMin <= 0.0 && WCMax > 0.0)
        {
            /* Atomic Nucleus */
            if (UCell[0] >  0.0  &&  VCell[0] >  0.0  &&  WCell[0] >  0.0  &&
                UCell[1] <= 0.0  &&  VCell[1] >  0.0  &&  WCell[1] >  0.0  &&
                UCell[2] <= 0.0  &&  VCell[2] <= 0.0  &&  WCell[2] >  0.0  &&
                UCell[3] >  0.0  &&  VCell[3] <= 0.0  &&  WCell[3] >  0.0  &&
                UCell[4] >  0.0  &&  VCell[4] >  0.0  &&  WCell[4] <= 0.0  &&
                UCell[5] <= 0.0  &&  VCell[5] >  0.0  &&  WCell[5] <= 0.0  &&
                UCell[6] <= 0.0  &&  VCell[6] <= 0.0  &&  WCell[6] <= 0.0  &&
                UCell[7] >  0.0  &&  VCell[7] <= 0.0  &&  WCell[7] <= 0.0)
            {
                /* Find corner with maximum rho */
                double RhoMax = TecUtilDataValueGetByRef(CDVarFDPtr, Index[0]);
                double rr[] = { -0.8,  0.8,  0.8, -0.8, -0.8,  0.8,  0.8, -0.8};
                double ss[] = { -0.8, -0.8,  0.8,  0.8, -0.8, -0.8,  0.8,  0.8};
                double tt[] = { -0.8, -0.8, -0.8, -0.8,  0.8,  0.8,  0.8,  0.8};
                int    icmax = 0;
                for (ic = 1; ic < 7; ic++)
                {
                    double RhoCrn = TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]);
                    if (RhoCrn > RhoMax)
                    {
                        RhoMax = RhoCrn;
                        icmax = ic;
                    }
                }
                r = rr[icmax];
                s = ss[icmax];
                t = tt[icmax];
            }
            else
            {
                r = 0.0;
                s = 0.0;
                t = 0.0;
            }

            CellFound =  TrilinearBrickNaturalCoord(0.0, 0.0, 0.0, UCell, VCell, WCell,
                                                    UCMin, VCMin, WCMin, UCMax, VCMax, WCMax,
                                                    &r, &s, &t, W);

            /* If Newton iteration failed, try some other initial conditions for r,s,t */
            if (!CellFound)
            {
                double rr[] = { -0.8,  0.8,  0.8, -0.8, -0.8,  0.8,  0.8, -0.8};
                double ss[] = { -0.8, -0.8,  0.8,  0.8, -0.8, -0.8,  0.8,  0.8};
                double tt[] = { -0.8, -0.8, -0.8, -0.8,  0.8,  0.8,  0.8,  0.8};
                int ictry = 0;
                while (!CellFound && ictry < 8)
                {
                    CellFound =  TrilinearBrickNaturalCoord(0.0, 0.0, 0.0, UCell, VCell, WCell,
                                                            UCMin, VCMin, WCMin, UCMax, VCMax, WCMax,
                                                            &r, &s, &t, W);
                    if (!CellFound)
                        ictry++;
                }
            }

            /* If Newton iteration failed, one last test for max */
            if (!CellFound)
            {
                LgIndex_t IMx, JMx, KMx;
                double    CellMaxChrgDens = TecUtilDataValueGetByRef(CDVarFDPtr, Index[0]);
                int       icmax = 0;
                LgIndex_t IndexMx;
                int       Type;

                /* Find the corner of the cell with the largest charge density */
                for (ic = 1; ic < 8; ic++)
                {
                    double ChrgDensCrn = TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]);
                    if (ChrgDensCrn > CellMaxChrgDens)
                    {
                        icmax = ic;
                        CellMaxChrgDens = ChrgDensCrn;
                    }
                }

                IMx = IIndex;
                if (icmax == 1 || icmax == 2 || icmax == 5 || icmax == 6) IMx++;

                JMx = JIndex;
                if (icmax == 2 || icmax == 3 || icmax == 6 || icmax == 7) JMx++;

                KMx = KIndex;
                if (icmax > 3) KMx++;

                IndexMx = IMx + (JMx - 1) * IMax + (KMx - 1) * IMax * JMax;
                Type = (int)TecUtilDataValueGetByRef(TypeVarFDPtr, IndexMx);

                /*
                 * Max if current point larger than the surrounding 26 points, and this
                 * max point hasn't already been found while testing another cell (Type>3).
                 */
                if (Type > -3 && IMx > 1 && IMx < IMax && JMx > 1 && JMx < JMax && KMx > 1 && KMx < KMax)
                {
                    Boolean_t MaxSoFar = TRUE;
                    Boolean_t MinSoFar = TRUE;
                    LgIndex_t ii, jj, kk;
                    LgIndex_t Indx;
                    double ChargeDensIni = TecUtilDataValueGetByRef(CDVarFDPtr, IndexMx);
                    double CDNeighbor;
                    for (kk = KMx - 1; MaxSoFar && kk <= KMx + 1; kk++)
                    {
                        for (jj = JMx - 1; MaxSoFar && jj <= JMx + 1; jj++)
                        {
                            for (ii = IMx - 1; MaxSoFar && ii <= IMx + 1; ii++)
                            {
                                if (!(ii == IMx && jj == JMx && kk == KMx))
                                {
                                    Indx = ii + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
                                    CDNeighbor = TecUtilDataValueGetByRef(CDVarFDPtr, Indx);
                                    if (CDNeighbor > ChargeDensIni) MaxSoFar = FALSE;
                                }
                            }
                        }
                    }
                    if (MaxSoFar)
                    {
                        CellFound = TRUE;
                        for (ic = 0; ic < 8; ic++) W[ic] = 0.0;
                        W[icmax] = 1.0;
                        TecUtilDataValueSetByRef(TypeVarFDPtr, Index[icmax], -3);
                    }
                }
            }


            /* If Newton iteration failed, one last test for min */
            if (!CellFound)
            {
                LgIndex_t IMn, JMn, KMn;
                double    CellMinChrgDens = TecUtilDataValueGetByRef(CDVarFDPtr, Index[0]);
                int       icmin = 0;
                LgIndex_t IndexMn;
                int       Type;

                /* Find the corner of the cell with the smallest charge density */
                for (ic = 1; ic < 8; ic++)
                {
                    double ChrgDensCrn = TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]);
                    if (ChrgDensCrn < CellMinChrgDens)
                    {
                        icmin = ic;
                        CellMinChrgDens = ChrgDensCrn;
                    }
                }

                IMn = IIndex;
                if (icmin == 1 || icmin == 2 || icmin == 5 || icmin == 6) IMn++;

                JMn = JIndex;
                if (icmin == 2 || icmin == 3 || icmin == 6 || icmin == 7) JMn++;

                KMn = KIndex;
                if (icmin > 3) KMn++;

                IndexMn = IMn + (JMn - 1) * IMax + (KMn - 1) * IMax * JMax;
                Type = (int)TecUtilDataValueGetByRef(TypeVarFDPtr, IndexMn);

                /*
                 * Min if current point smaller than the surrounding 26 points, and this
                 * min point hasn't already been found while testing another cell (Type>3).
                 */
                if (Type < 3 && IMn > 1 && IMn < IMax && JMn > 1 && JMn < JMax && KMn > 1 && KMn < KMax)
                {
                    Boolean_t MinSoFar = TRUE;
                    LgIndex_t ii, jj, kk;
                    LgIndex_t Indx;
                    double ChargeDensIni = TecUtilDataValueGetByRef(CDVarFDPtr, IndexMn);
                    double CDNeighbor;
                    for (kk = KMn - 1; MinSoFar && kk <= KMn + 1; kk++)
                    {
                        for (jj = JMn - 1; MinSoFar && jj <= JMn + 1; jj++)
                        {
                            for (ii = IMn - 1; MinSoFar && ii <= IMn + 1; ii++)
                            {
                                if (!(ii == IMn && jj == JMn && kk == KMn))
                                {
                                    Indx = ii + (jj - 1) * IMax + (kk - 1) * IMax * JMax;
                                    CDNeighbor = TecUtilDataValueGetByRef(CDVarFDPtr, Indx);
                                    if (CDNeighbor < ChargeDensIni) MinSoFar = FALSE;
                                }
                            }
                        }
                    }
                    if (MinSoFar)
                    {
                        CellFound = TRUE;
                        for (ic = 0; ic < 8; ic++) W[ic] = 0.0;
                        W[icmin] = 1.0;
                        TecUtilDataValueSetByRef(TypeVarFDPtr, Index[icmin], 3);
                    }
                }
            }




            /* Compute XYZ location of critical point */
            if (CellFound)
            {
                double D2RhoDx2  = 0.0;
                double D2RhoDxDy = 0.0;
                double D2RhoDxDz = 0.0;
                double D2RhoDy2  = 0.0;
                double D2RhoDyDz = 0.0;
                double D2RhoDz2  = 0.0;

                *XCrtPt = 0.0;
                *YCrtPt = 0.0;
                *ZCrtPt = 0.0;
                *ChrgDensCrtPt = 0.0;
                for (ic = 0; ic < 8; ic++)
                {
                    *XCrtPt += TecUtilDataValueGetByRef(XVarFDPtr, Index[ic]) * W[ic];
                    *YCrtPt += TecUtilDataValueGetByRef(YVarFDPtr, Index[ic]) * W[ic];
                    *ZCrtPt += TecUtilDataValueGetByRef(ZVarFDPtr, Index[ic]) * W[ic];
                    *ChrgDensCrtPt += TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]) * W[ic];
                }

                /* Compute second derivatives for curvature - classification.
                 * For now, assume rectangular grid with spacing 1 in each direction.
                 */
                for (ic = 0; ic < 8; ic++)
                {
                    double d2dx2, d2dxdy, d2dxdz, d2dy2, d2dydz, d2dz2;
                    LgIndex_t i, j, k;

                    i = IIndex;
                    if (ic == 1 || ic == 2 || ic == 5 || ic == 6) i++;

                    j = JIndex;
                    if (ic == 2 || ic == 3 || ic == 6 || ic == 7) j++;

                    k = KIndex;
                    if (ic > 3) k++;

                    IsOk = ComputeCurvature(ZoneNum, i, j, k, CDVarFDPtr,
                                            &d2dx2, &d2dxdy, &d2dxdz,
                                            &d2dy2, &d2dydz, &d2dz2);

                    D2RhoDx2  += W[ic] * d2dx2;
                    D2RhoDxDy += W[ic] * d2dxdy;
                    D2RhoDxDz += W[ic] * d2dxdz;
                    D2RhoDy2  += W[ic] * d2dy2;
                    D2RhoDyDz += W[ic] * d2dydz;
                    D2RhoDz2  += W[ic] * d2dz2;
                }

                /* Compute eigensystem of curvature tensor */
                if (IsOk)
                {
                    Boolean_t SortEgnV = TRUE;
                    float **InMat, **EgnVec, a[3][3], EgnVal[3], b[3][3];
                    int   NumRot = 0;
                    EntIndex_t CritPointType;

                    /* Set up matrices (pointers to pointers) */
                    if (IsOk)
                    {
                        InMat = ALLOC_ARRAY(3, float*, "InMat");
                        if (InMat == NULL)
                            IsOk = FALSE;
                        else
                        {
                            InMat[0] = a[0];
                            InMat[1] = a[1];
                            InMat[2] = a[2];
                        }
                    }

                    if (IsOk)
                        EgnVec = ALLOC_ARRAY(3, float*, "EgnVec");
                    if (EgnVec == NULL)
                        IsOk = FALSE;
                    else
                    {
                        EgnVec[0] = b[0];
                        EgnVec[1] = b[1];
                        EgnVec[2] = b[2];
                    }

                    a[0][0] = (float)D2RhoDx2;
                    a[0][1] = (float)D2RhoDxDy;
                    a[0][2] = (float)D2RhoDxDz;
                    a[1][0] = a[0][1];
                    a[1][1] = (float)D2RhoDy2;
                    a[1][2] = (float)D2RhoDyDz;
                    a[2][0] = a[0][2];
                    a[2][1] = a[1][2];
                    a[2][2] = (float)D2RhoDz2;

                    IsOk = Jacobi3by3(InMat, EgnVal, EgnVec, &NumRot);

                    /* Sort by eigenvalue value */
                    if (SortEgnV && IsOk)
                    {
                        if (EgnVal[1] < EgnVal[2])
                        {
                            int row;
                            float EgnVTmp = EgnVal[1];
                            EgnVal[1] = EgnVal[2];
                            EgnVal[2] = EgnVTmp;

                            for (row = 0; row < 3; row++)
                            {
                                EgnVTmp = EgnVec[row][1];
                                EgnVec[row][1] = EgnVec[row][2];
                                EgnVec[row][2] = EgnVTmp;
                            }
                        }
                        if (EgnVal[0] < EgnVal[1])
                        {
                            int row;
                            float EgnVTmp = EgnVal[0];
                            EgnVal[0] = EgnVal[1];
                            EgnVal[1] = EgnVTmp;

                            for (row = 0; row < 3; row++)
                            {
                                EgnVTmp = EgnVec[row][0];
                                EgnVec[row][0] = EgnVec[row][1];
                                EgnVec[row][1] = EgnVTmp;
                            }
                        }
                        if (EgnVal[1] < EgnVal[2])
                        {
                            int row;
                            float EgnVTmp = EgnVal[1];
                            EgnVal[1] = EgnVal[2];
                            EgnVal[2] = EgnVTmp;

                            for (row = 0; row < 3; row++)
                            {
                                EgnVTmp = EgnVec[row][1];
                                EgnVec[row][1] = EgnVec[row][2];
                                EgnVec[row][2] = EgnVTmp;
                            }
                        }
                    }

                    /* Determine critical point classifications */
                    if (EgnVal[2] < 0.0)
                    {
                        if (EgnVal[1] < 0.0)
                        {
                            if (EgnVal[0] < 0.0)   /* Atom */
                                CritPointType = -3;
                            else                   /* Bond */
                                CritPointType = -1;
                        }
                        else                       /* Ring */
                            CritPointType = 1;
                    }
                    else                           /* Cage */
                        CritPointType = 3;

                    *TypeCrtPt = CritPointType;

                    /* If saddle point, set Principle Direction vector. */
                    if (CritPointType == -1)
                    {
                        *PrincDirX = EgnVec[0][0];
                        *PrincDirY = EgnVec[1][0];
                        *PrincDirZ = EgnVec[2][0];
                    }
                    else if (CritPointType == 1)
                    {
                        *PrincDirX = EgnVec[0][2];
                        *PrincDirY = EgnVec[1][2];
                        *PrincDirZ = EgnVec[2][2];
                    }


                    /* If close to node - mark as max or min already found */
                    if (IsOk && CellFound)
                    {
                        double WMax = -1.0;
                        int icmax = -1;
                        for (ic = 0; ic < 8; ic++)
                        {
                            if (W[ic] > WMax)
                            {
                                WMax = W[ic];
                                icmax = ic;
                            }
                        }
                        if (CritPointType == -3 || CritPointType == 3 && WMax > 0.9)
                            TecUtilDataValueSetByRef(TypeVarFDPtr, Index[icmax], CritPointType);
                    }
                    if (InMat != NULL)
                        FREE_ARRAY(InMat, "InMat");
                    if (EgnVec != NULL)
                        FREE_ARRAY(EgnVec, "EgnVec");
                }
            }
        }
        else
            CellFound = FALSE;
    }

    return(CellFound);
}





/**
 * Determine if the CritPoints handle is sane.
 *
 * param CritPoints
 *     CritPoints structure in question.
 *
 * return
 *     TRUE if the CritPoints structure is valid, otherwise FALSE.
 */
Boolean_t CritPointsIsValid(CritPoints_pa CritPoints)
{
    Boolean_t IsValid = FALSE;

    IsValid = (VALID_REF(CritPoints) &&
               VALID_REF(CritPoints->X) && ArrListIsValid(CritPoints->X) &&
               VALID_REF(CritPoints->Y) && ArrListIsValid(CritPoints->Y) &&
               VALID_REF(CritPoints->Z) && ArrListIsValid(CritPoints->Z) &&
               VALID_REF(CritPoints->ChrgDens)  && ArrListIsValid(CritPoints->ChrgDens) &&
               VALID_REF(CritPoints->Type)      && ArrListIsValid(CritPoints->Type) &&
               VALID_REF(CritPoints->PrincDirX) && ArrListIsValid(CritPoints->PrincDirX) &&
               VALID_REF(CritPoints->PrincDirY) && ArrListIsValid(CritPoints->PrincDirY) &&
               VALID_REF(CritPoints->PrincDirZ) && ArrListIsValid(CritPoints->PrincDirZ));

    /* Require the same count for each array list in CritPoints structure. */
    if (IsValid)
    {
        LgIndex_t Count = ArrListGetCount(CritPoints->X);
        IsValid = (ArrListGetCount(CritPoints->Y) == Count);
        if (IsValid) IsValid = (ArrListGetCount(CritPoints->Z) == Count);
        if (IsValid) IsValid = (ArrListGetCount(CritPoints->ChrgDens) == Count);
        if (IsValid) IsValid = (ArrListGetCount(CritPoints->Type) == Count);
        if (IsValid) IsValid = (ArrListGetCount(CritPoints->PrincDirX) == Count);
        if (IsValid) IsValid = (ArrListGetCount(CritPoints->PrincDirY) == Count);
        if (IsValid) IsValid = (ArrListGetCount(CritPoints->PrincDirZ) == Count);
        if (IsValid) IsValid = (Count == (CritPoints->NumCrtPtsM3 + CritPoints->NumCrtPtsM1 +
                                              CritPoints->NumCrtPtsP1 + CritPoints->NumCrtPtsP3));
        if (IsValid) IsValid = (CritPoints->NumCrtPts == Count);
    }


    ENSURE(VALID_BOOLEAN(IsValid));
    return IsValid;
}




/**
 * Deallocates the CritPoints handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a CritPoints handle.
 */
void CritPointsDealloc(CritPoints_pa *CritPoints)
{
    REQUIRE(VALID_REF(CritPoints));
    REQUIRE(CritPointsIsValid(*CritPoints) || *CritPoints == NULL);

    if (*CritPoints != NULL)
    {
        /* release the ArrList's */
        if ((*CritPoints)->X != NULL) ArrListDealloc(&((*CritPoints)->X));
        if ((*CritPoints)->Y != NULL) ArrListDealloc(&((*CritPoints)->Y));
        if ((*CritPoints)->Z != NULL) ArrListDealloc(&((*CritPoints)->Z));
        if ((*CritPoints)->ChrgDens  != NULL) ArrListDealloc(&((*CritPoints)->ChrgDens));
        if ((*CritPoints)->Type      != NULL) ArrListDealloc(&((*CritPoints)->Type));
        if ((*CritPoints)->PrincDirX != NULL) ArrListDealloc(&((*CritPoints)->PrincDirX));
        if ((*CritPoints)->PrincDirY != NULL) ArrListDealloc(&((*CritPoints)->PrincDirY));
        if ((*CritPoints)->PrincDirZ != NULL) ArrListDealloc(&((*CritPoints)->PrincDirZ));

        /* Notify Tecplot of memory usage change */
        TecUtilMemoryChangeNotify(- (Int64_t)((*CritPoints)->MemUsageReported));
        (*CritPoints)->MemUsageReported = 0;

        /* release the list structure itself */
        FREE_ITEM(*CritPoints, "GradPath structure");
        *CritPoints = NULL;
    }

    ENSURE(*CritPoints == NULL);
}





/**
 * Empties the CritPoints structure of all Critical Points and resets the
 * principle directions and number of Critical Points to 0.
 *
 *
 * param CritPoints
 *     CritPoints to clear.
 */
void CritPointsClear(CritPoints_pa CritPoints)
{
    REQUIRE(CritPointsIsValid(CritPoints));

    CritPoints->NumCrtPts   = 0;
    CritPoints->NumCrtPtsM3 = 0;
    CritPoints->NumCrtPtsM1 = 0;
    CritPoints->NumCrtPtsP1 = 0;
    CritPoints->NumCrtPtsP3 = 0;

    ArrListClear(CritPoints->X);
    ArrListClear(CritPoints->Y);
    ArrListClear(CritPoints->Z);
    ArrListClear(CritPoints->ChrgDens);
    ArrListClear(CritPoints->Type);
    ArrListClear(CritPoints->PrincDirX);
    ArrListClear(CritPoints->PrincDirY);
    ArrListClear(CritPoints->PrincDirZ);

    ENSURE(CritPointsIsValid(CritPoints) && CritPointsGetCount(CritPoints) == 0);
}




/**
 * Critical point arrays are ordered by their critical point type: all atoms
 * (Type=-3) first, then Bonds (Type=-1), Rings (Type=1), and Cages (Type=3).
 * Get the beginning offset of Type into CritPoints arrays.
 *
 *
 * param CritPoints
 *     Critical Point structure containing the arrays of critical points.
 * param ItemOffset
 *     Type of critcal point.
 *
 * return
 *     Offset of beginning of critical point Type.
 */
LgIndex_t CritPointsGetBegOffset(const CritPoints_pa CritPoints,
                                 const char          Type)
{
    LgIndex_t BegOffset = 0;

    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(Type == -3 || Type == -1 || Type == 1 || Type == 3);

    /* Compute Offsets to the various critical point types */
    if (Type > -3) BegOffset += CritPoints->NumCrtPtsM3;
    if (Type > -1) BegOffset += CritPoints->NumCrtPtsM1;
    if (Type >  1) BegOffset += CritPoints->NumCrtPtsP1;

    ENSURE(0 <= BegOffset && BegOffset <= CritPointsGetCount(CritPoints));
    return BegOffset;
}





/**
 * Critical point arrays are ordered by their critical point type: all atoms
 * (Type=-3) first, then Bonds (Type=-1), Rings (Type=1), and Cages (Type=3).
 * Get the end offset (+1) of Type in CritPoints arrays.
 *
 *
 * param CritPoints
 *     Critical Point structure containing the arrays of critical points.
 * param ItemOffset
 *     Type of critcal point.
 *
 * return
 *     Offset of beginning of critical point Type.
 */
LgIndex_t CritPointsGetEndOffset(const CritPoints_pa CritPoints,
                                 const char          Type)
{
    LgIndex_t EndOffset = CritPoints->NumCrtPtsM3;

    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(Type == -3 || Type == -1 || Type == 1 || Type == 3);

    /* Compute Offsets to the various critical point types */
    if (Type > -3) EndOffset += CritPoints->NumCrtPtsM1;
    if (Type > -1) EndOffset += CritPoints->NumCrtPtsP1;
    if (Type >  1) EndOffset += CritPoints->NumCrtPtsP3;

    ENSURE(0 <= EndOffset && EndOffset <= CritPointsGetCount(CritPoints));
    return EndOffset;
}





/**
 * Get the Type of the critical point given it's offset into CritPoints.
 * Critical point arrays are ordered by their critical point type: all atoms
 * (Type=-3) first, then Bonds (Type=-1), Rings (Type=1), and Cages (Type=3).
 *
 *
 * param CritPoints
 *     Critical Point structure containing the arrays of critical points.
 * param PointOffset
 *     Offset of critcal point.
 *
 * return
 *     Critical point Type. (-3 = Atom,  -1 = Bond,  1 = Ring,  3 = Cage)
 */
char CritPointsGetType(const CritPoints_pa CritPoints,
                       const LgIndex_t     PointOffset)
{
    char Type = -3;

    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(0 <= PointOffset && PointOffset < CritPointsGetCount(CritPoints));

    /* Compute type based on Offsets */
    if (PointOffset >= CritPointsGetBegOffset(CritPoints, -1))
    {
        if (PointOffset >= CritPointsGetBegOffset(CritPoints, 1))
        {
            if (PointOffset >= CritPointsGetBegOffset(CritPoints, 3))
                Type = 3;
            else
                Type = 1;
        }
        else
            Type = -1;
    }

    ENSURE(Type == -3 || Type == -1 || Type == 1 || Type == 3);
    return Type;
}




/**
 * Removes a point from the CritPoints array. The members following the item
 * removed are shifted down accordingly to fill the vacated space.
 *
 *
 * param CritPoints
 *     Critical Point structure containing the point to remove.
 * param ItemOffset
 *     Offset to the point.
 *
 * return
 *     TRUE if successful, FALSE otherwise.
 */
Boolean_t CritPointsRemovePoint(CritPoints_pa CritPoints,
                                LgIndex_t     PointOffset)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t OffsetM1, OffsetP1, OffsetP3, NumCrtPts;

    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(0 <= PointOffset && PointOffset <= CritPointsGetCount(CritPoints) - 1);

    /* Compute Offsets to the various critical point types */
    NumCrtPts = CritPoints->NumCrtPts;
    OffsetM1  = CritPointsGetBegOffset(CritPoints, -1);
    OffsetP1  = CritPointsGetBegOffset(CritPoints,  1);
    OffsetP3  = CritPointsGetBegOffset(CritPoints,  3);

    /* Decrement the appropriate CrtPt count */
    (CritPoints->NumCrtPts)--;
    if (PointOffset <  OffsetM1)(CritPoints->NumCrtPtsM3)--;
    if (OffsetM1 <= PointOffset && PointOffset < OffsetP1)(CritPoints->NumCrtPtsM1)--;
    if (OffsetP1 <= PointOffset && PointOffset < OffsetP3)(CritPoints->NumCrtPtsP1)--;
    if (OffsetP3 <= PointOffset)(CritPoints->NumCrtPtsP3)--;

    /* Remove the items for the array lists */
    ArrListRemoveItem(CritPoints->X, PointOffset);
    ArrListRemoveItem(CritPoints->Y, PointOffset);
    ArrListRemoveItem(CritPoints->Z, PointOffset);
    ArrListRemoveItem(CritPoints->ChrgDens, PointOffset);
    ArrListRemoveItem(CritPoints->Type, PointOffset);
    ArrListRemoveItem(CritPoints->PrincDirX, PointOffset);
    ArrListRemoveItem(CritPoints->PrincDirY, PointOffset);
    ArrListRemoveItem(CritPoints->PrincDirZ, PointOffset);

    IsOk = CritPointsIsValid(CritPoints);

    ENSURE(CritPointsIsValid(CritPoints));
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}



/**
 * Allocates a CritPoints handle with a suitable default
 * capacity for the contained ArrList's.
 *
 *
 * return
 *     CritPoints handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
CritPoints_pa CritPointsAlloc()
{
    CritPoints_pa Result = NULL;

    Result = ALLOC_ITEM(CritPoints_s, "CritPoints structure");
    if (Result != NULL)
    {
        Result->MemUsageReported = 0;
        Result->NumCrtPts   = 0;
        Result->NumCrtPtsM3 = 0;
        Result->NumCrtPtsM1 = 0;
        Result->NumCrtPtsP1 = 0;
        Result->NumCrtPtsP3 = 0;
        Result->X           = ArrListAlloc(60, ArrListType_Double);
        Result->Y           = ArrListAlloc(60, ArrListType_Double);
        Result->Z           = ArrListAlloc(60, ArrListType_Double);
        Result->ChrgDens    = ArrListAlloc(60, ArrListType_Double);
        Result->Type        = ArrListAlloc(60, ArrListType_Char);
        Result->PrincDirX   = ArrListAlloc(60, ArrListType_Double);
        Result->PrincDirY   = ArrListAlloc(60, ArrListType_Double);
        Result->PrincDirZ   = ArrListAlloc(60, ArrListType_Double);

        /* If it failed to allocate any of the array lists, clean-up and exit. */
        if (Result->X == NULL || Result->Y == NULL || Result->Z == NULL ||
            Result->ChrgDens  == NULL || Result->Type == NULL ||
            Result->PrincDirX == NULL || Result->PrincDirY == NULL ||
            Result->PrincDirZ == NULL)
        {
            if (Result->X != NULL) ArrListDealloc(&(Result->X));
            if (Result->Y != NULL) ArrListDealloc(&(Result->Y));
            if (Result->Z != NULL) ArrListDealloc(&(Result->Z));
            if (Result->ChrgDens  != NULL) ArrListDealloc(&(Result->ChrgDens));
            if (Result->Type      != NULL) ArrListDealloc(&(Result->Type));
            if (Result->PrincDirX != NULL) ArrListDealloc(&(Result->PrincDirX));
            if (Result->PrincDirY != NULL) ArrListDealloc(&(Result->PrincDirY));
            if (Result->PrincDirZ != NULL) ArrListDealloc(&(Result->PrincDirZ));
            FREE_ITEM(Result, "CritPoints structure");
            Result = NULL;
        }
    }

    ENSURE(CritPointsIsValid(Result) || Result == NULL);
    return Result;
}




/**
 * Gets the number of points (nodes) currently in the CritPoints
 * (maintained by the CritPoints coordinate arrays).
 *
 * param
 *     CritPoints structure in question.
 *
 * return
 *     Number of points (nodes) in the CritPoints.
 */
LgIndex_t CritPointsGetCount(CritPoints_pa CritPoints)
{
    LgIndex_t Result = 0;

    REQUIRE(CritPointsIsValid(CritPoints));

    Result = ArrListGetCount(CritPoints->X);

    ENSURE(Result >= 0);
    return Result;
}






/**
 * Places critical point components at the specified offset. If the
 * offset is beyond the end of the array lists, they are sized
 * accordingly and the intervening array values between the last
 * node of the original state and the last node of the new state are
 * assigned 0.
 *
 * note
 *     If components already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param CritPoints
 *     CritPoints target in which to set the coordinates.
 * param PointOffset
 *     Offset of the point/node.
 * param X,Y,Z
 *     Coordinates to set at the specified offset.
 * param ChrgDens, Type
 *     Charge density and critical point type to set at the specified
 *     offset.
 * param PrincDirX, ...Y, ...Z
 *     Principle direction components to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t CritPointsSetPoint(CritPoints_pa CritPoints,
                             LgIndex_t     PointOffset,
                             double        X,
                             double        Y,
                             double        Z,
                             double        ChrgDens,
                             char          Type,
                             double        PrincDirX,
                             double        PrincDirY,
                             double        PrincDirZ)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(PointOffset >= 0);

    Item.Double = X;
    IsOk = ArrListSetItem(CritPoints->X, PointOffset, Item);

    if (IsOk)
    {
        Item.Double = Y;
        IsOk = ArrListSetItem(CritPoints->Y, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = Z;
        IsOk = ArrListSetItem(CritPoints->Z, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = ChrgDens;
        IsOk = ArrListSetItem(CritPoints->ChrgDens, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Char = Type;
        IsOk = ArrListSetItem(CritPoints->Type, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = PrincDirX;
        IsOk = ArrListSetItem(CritPoints->PrincDirX, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = PrincDirY;
        IsOk = ArrListSetItem(CritPoints->PrincDirY, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = PrincDirZ;
        IsOk = ArrListSetItem(CritPoints->PrincDirZ, PointOffset, Item);
    }

    // TODO: What if point is already set? Need to modify following logic.
    /* Update the NumCrtPts counters. */
    if (IsOk)
    {
        switch (Type)
        {
            case -3:
                (CritPoints->NumCrtPtsM3)++;
                break;
            case -1:
                (CritPoints->NumCrtPtsM1)++;
                break;
            case 1:
                (CritPoints->NumCrtPtsP1)++;
                break;
            case 3:
                (CritPoints->NumCrtPtsP3)++;
                break;
        }
        (CritPoints->NumCrtPts)++;
    }

    /* Require the same count for each coordinate array (checked by GradPathIsValid). */
    if (IsOk) IsOk = CritPointsIsValid(CritPoints);

    ENSURE(CritPointsIsValid(CritPoints));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Inserts critical point components at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 *
 * param CritPoints
 *     CritPoints target in which to set the coordinates.
 * param PointOffset
 *     Offset at which to insert the point/node.
 * param X,Y,Z
 *     Coordinates to set at the specified offset.
 * param ChrgDens, Type
 *     Charge density and critical point type to set at the specified
 *     offset.
 * param PrincDirX, ...Y, ...Z
 *     Principle direction components to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t CritPointsInsertPoint(CritPoints_pa CritPoints,
                                LgIndex_t     PointOffset,
                                double        X,
                                double        Y,
                                double        Z,
                                double        ChrgDens,
                                char          Type,
                                double        PrincDirX,
                                double        PrincDirY,
                                double        PrincDirZ)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;


    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(0 <= PointOffset && PointOffset <= CritPointsGetCount(CritPoints));

    Item.Double = X;
    IsOk = ArrListInsertItem(CritPoints->X, PointOffset, Item);

    if (IsOk)
    {
        Item.Double = Y;
        IsOk = ArrListInsertItem(CritPoints->Y, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = Z;
        IsOk = ArrListInsertItem(CritPoints->Z, PointOffset, Item);
    }
    if (IsOk)
    {
        Item.Double = ChrgDens;
        IsOk = ArrListInsertItem(CritPoints->ChrgDens, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Char = Type;
        IsOk = ArrListInsertItem(CritPoints->Type, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = PrincDirX;
        IsOk = ArrListInsertItem(CritPoints->PrincDirX, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = PrincDirY;
        IsOk = ArrListInsertItem(CritPoints->PrincDirY, PointOffset, Item);
    }

    if (IsOk)
    {
        Item.Double = PrincDirZ;
        IsOk = ArrListInsertItem(CritPoints->PrincDirZ, PointOffset, Item);
    }

    /* Update the NumCrtPts counters. */
    if (IsOk)
    {
        switch (Type)
        {
            case -3:
                (CritPoints->NumCrtPtsM3)++;
                break;
            case -1:
                (CritPoints->NumCrtPtsM1)++;
                break;
            case 1:
                (CritPoints->NumCrtPtsP1)++;
                break;
            case 3:
                (CritPoints->NumCrtPtsP3)++;
                break;
        }
        (CritPoints->NumCrtPts)++;
    }

    /* Require the same count for each coordinate array (checked by GradPathIsValid). */
    if (IsOk) IsOk = CritPointsIsValid(CritPoints);

    ENSURE(CritPointsIsValid(CritPoints));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/**
 * Appends the critical point components to CritPoints. The array lists
 * will be expanded to accommodate the additional items.
 *
 * param CritPoints
 *     CritPoints target to which the point/node is to be appended.
 * param X,Y,Z
 *     Coordinates of node to append to the CritPoints.
 * param ChrgDens, Type
 *     Charge density and critical point type to set at the specified
 *     offset.
 * param PrincDirX, ...Y, ...Z
 *     Principle direction components to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t CritPointsAppendPoint(CritPoints_pa  CritPoints,
                                double         X,
                                double         Y,
                                double         Z,
                                double         ChrgDens,
                                char           Type,
                                double         PrincDirX,
                                double         PrincDirY,
                                double         PrincDirZ)
{
    Boolean_t IsOk = TRUE;
    LgIndex_t Count;

    REQUIRE(CritPointsIsValid(CritPoints));

    Count = CritPointsGetCount(CritPoints);

    IsOk = CritPointsInsertPoint(CritPoints, Count, X, Y, Z, ChrgDens, Type,
                                 PrincDirX, PrincDirY, PrincDirZ);

    ENSURE(CritPointsIsValid(CritPoints));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}






/**
 * Gets the components of the CritPoints structure for the
 * critical point at the specified offset.
 *
 * param CritPoints
 *     CritPoints structure containing the desired item.
 * param PointOffset
 *     Offset to the coordinates in the CritPoints.
 * param *X, *Y, *Z
 *     Pointers to coordinates of the node
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */

Boolean_t CritPointsGetPoint(CritPoints_pa  CritPoints,
                             LgIndex_t      PointOffset,
                             double        *X,
                             double        *Y,
                             double        *Z,
                             double        *ChrgDens,
                             char          *Type,
                             double        *PrincDirX,
                             double        *PrincDirY,
                             double        *PrincDirZ)
{
    Boolean_t IsOk = TRUE;
    ArrListItem_u Item;

    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(PointOffset >= 0 && PointOffset < CritPointsGetCount(CritPoints));
    REQUIRE(VALID_REF(X));
    REQUIRE(VALID_REF(Y));
    REQUIRE(VALID_REF(Z));
    REQUIRE(VALID_REF(ChrgDens));
    REQUIRE(VALID_REF(Type));
    REQUIRE(VALID_REF(PrincDirX));
    REQUIRE(VALID_REF(PrincDirY));
    REQUIRE(VALID_REF(PrincDirZ));

    Item = ArrListGetItem(CritPoints->X, PointOffset);
    *X = Item.Double;

    Item = ArrListGetItem(CritPoints->Y, PointOffset);
    *Y = Item.Double;

    Item = ArrListGetItem(CritPoints->Z, PointOffset);
    *Z = Item.Double;

    Item = ArrListGetItem(CritPoints->ChrgDens, PointOffset);
    *ChrgDens = Item.Double;

    Item = ArrListGetItem(CritPoints->Type, PointOffset);
    *Type = Item.Char;

    Item = ArrListGetItem(CritPoints->PrincDirX, PointOffset);
    *PrincDirX = Item.Double;

    Item = ArrListGetItem(CritPoints->PrincDirY, PointOffset);
    *PrincDirY = Item.Double;

    Item = ArrListGetItem(CritPoints->PrincDirZ, PointOffset);
    *PrincDirZ = Item.Double;

    ENSURE(VALID_REF(X));
    ENSURE(VALID_REF(Y));
    ENSURE(VALID_REF(Z));
    ENSURE(VALID_BOOLEAN(IsOk));
    ENSURE(VALID_REF(ChrgDens));
    ENSURE(VALID_REF(Type));
    ENSURE(VALID_REF(PrincDirX));
    ENSURE(VALID_REF(PrincDirY));
    ENSURE(VALID_REF(PrincDirZ));

    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}





/*
 * Function to remove duplicate critical points from the lists
 *
 * TestForNode: If true, the one of the duplicate CP's must be located at a node.
 *              The CP at the node is the one removed from the list.
 * Tolerance:   Separation distance below which the nodes are considered duplicate.
 * Type:        Search for duplicates on within the critical point type:
 *              -3 = Atom, -1 = Bond, 1 = Ring, 3 = Cage
 */
Boolean_t RemDuplicateCritPoints(Boolean_t     TestForNode,
                                 double        Tolerance,
                                 char          Type,
                                 CritPoints_pa CritPoints)
{
    Boolean_t  IsOk = TRUE;
    LgIndex_t itst, ii;
    LgIndex_t NumCrtPts = CritPointsGetCount(CritPoints);

    REQUIRE(VALID_BOOLEAN(TestForNode));
    REQUIRE(Tolerance >= 0.0);
    REQUIRE(VALID_REF(CritPoints));
    REQUIRE(CritPointsIsValid(CritPoints));

    if (NumCrtPts > 1)
    {
        /* Compute begin and end critical point numbers */
        LgIndex_t BeginCrtPt = CritPointsGetBegOffset(CritPoints, Type);
        LgIndex_t EndCrtPt   = CritPointsGetEndOffset(CritPoints, Type);

        /* Loop through critical points, testing for duplicates */
        for (itst = BeginCrtPt; itst < CritPointsGetEndOffset(CritPoints, Type); itst++)
        {
            double XDistFromNode, YDistFromNode, ZDistFromNode;
            double DistFromNodeSq;
            double Xtst, Ytst, Ztst, dummy;
            char   cdummy;

            IsOk = CritPointsGetPoint(CritPoints, itst, &Xtst, &Ytst, &Ztst,
                                      &dummy, &cdummy, &dummy, &dummy, &dummy);

            XDistFromNode = ROUND2(Xtst) - Xtst;
            YDistFromNode = ROUND2(Ytst) - Ytst;
            ZDistFromNode = ROUND2(Ztst) - Ztst;
            DistFromNodeSq = XDistFromNode * XDistFromNode + YDistFromNode * YDistFromNode
                             + ZDistFromNode * ZDistFromNode;

            /* If TestForNode, remove only critical points that are at the node. */
            if (!TestForNode || DistFromNodeSq < 1.0e-7)
            {
                Boolean_t DupFound = FALSE;
                for (ii = BeginCrtPt; !DupFound && ii < CritPointsGetEndOffset(CritPoints, Type); ii++)
                {
                    if (ii != itst)
                    {
                        double Xii, Yii, Zii, DelX, DelY, DelZ;
                        IsOk = CritPointsGetPoint(CritPoints, ii, &Xii, &Yii, &Zii,
                                                  &dummy, &cdummy, &dummy, &dummy, &dummy);

                        DelX = Xii - Xtst;
                        DelY = Yii - Ytst;
                        DelZ = Zii - Ztst;

                        /* Points are duplicate if they are in the same quadrant of cell */
                        if (ABS(DelX) < Tolerance && ABS(DelY) < Tolerance && ABS(DelZ) < Tolerance)
                        {
                            // LgIndex_t idlt;
                            DupFound = TRUE;

                            /* Delete the itst node */
                            IsOk = CritPointsRemovePoint(CritPoints, itst);

                            /* Do not increment itst if the current one is deleted */
                            itst--;
                        }
                    }
                }
            }
        }
    }

    ENSURE(VALID_BOOLEAN(IsOk));
    return (IsOk);
}








Boolean_t ExtractCriticalPoints(const EntIndex_t    ZoneNum,
                                const EntIndex_t    UVarNum,
                                const EntIndex_t    VVarNum,
                                const EntIndex_t    WVarNum,
                                const EntIndex_t    ChrgDensVarNum,
                                const EntIndex_t    TypeVarNum,
                                CritPoints_pa       CritPoints)
{
    Boolean_t  IsOk = TRUE;
    Boolean_t  IsStop = FALSE;
    Boolean_t  ShowStatusBar = TRUE;
    EntIndex_t NumZones, NumVars;
    LgIndex_t  IMax, JMax, KMax;
    LgIndex_t  ii, jj, kk;

    REQUIRE(TecUtilDataSetIsAvailable());

    /* Get num zones in dataset */
    if (IsOk)
    {
        char *DatasetTitle = NULL;
        TecUtilDataSetGetInfo(&DatasetTitle, &NumZones, &NumVars);
        TecUtilStringDealloc(&DatasetTitle);
    }

    REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
    REQUIRE(TecUtilZoneIsOrdered(ZoneNum));
    REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
    REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
    REQUIRE(WVarNum > 0 && WVarNum <= NumVars);
    REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
    REQUIRE(CritPointsIsValid(CritPoints));


    TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */

    /* Set-up status bar */
    if (ShowStatusBar)
        TecUtilStatusStartPercentDone("Finding Critical Points", TRUE, TRUE);

    /* Set SourceZoneNum */
    CritPoints->SourceZoneNum = ZoneNum;

    /* Set Critical Point Type Var to zero for source ZoneNum */
    /*
    if (IsOk)
      {
        FieldData_pa TypeVarFDPtr = NULL;
        LgIndex_t    ii;
        LgIndex_t    NumNodex = IMax * JMax * KMax;

        TecUtilDataLoadBegin();
        TypeVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, TypeVarNum);

        for (ii=0; ii<NumNodex; ii++)
          {
            TecUtilDataValueSetByRef(TypeVarFDPtr, ii+1, 0);
          }
        TecUtilDataLoadEnd();
      }
      */

    /* Inform Tecplot that major data operation is beginning */
    TecUtilDataLoadBegin();

    /* Find the Critical Points */
    for (kk = 1; !IsStop && kk < KMax; kk++)
    {
        for (jj = 1; jj < JMax; jj++)
        {
            for (ii = 1; ii < IMax; ii++)
            {
                double XCrtPt, YCrtPt, ZCrtPt, ChrgDensCrtPt;
                EntIndex_t TypeCrtPt;
                float PrincDirX, PrincDirY, PrincDirZ;

                if (CriticalPointInCell(ZoneNum, UVarNum, VVarNum, WVarNum, ChrgDensVarNum, TypeVarNum,
                                        ii, jj, kk, &XCrtPt, &YCrtPt, &ZCrtPt,
                                        &ChrgDensCrtPt, &TypeCrtPt, &PrincDirX, &PrincDirY, &PrincDirZ))
                {

                    switch (TypeCrtPt)
                    {
                        case -3:
                        {
                            LgIndex_t PointOffset = CritPoints->NumCrtPtsM3;
                            IsOk =  CritPointsInsertPoint(CritPoints, PointOffset, XCrtPt, YCrtPt, ZCrtPt,
                                                          ChrgDensCrtPt, (char)TypeCrtPt, 0.0, 0.0, 0.0);
                        }
                        break;
                        case -1:
                        {
                            LgIndex_t PointOffset = CritPoints->NumCrtPtsM3 + CritPoints->NumCrtPtsM1;
                            IsOk =  CritPointsInsertPoint(CritPoints, PointOffset, XCrtPt, YCrtPt, ZCrtPt,
                                                          ChrgDensCrtPt, (char)TypeCrtPt, PrincDirX, PrincDirY, PrincDirZ);
                        }
                        break;
                        case 1:
                        {
                            LgIndex_t PointOffset = CritPoints->NumCrtPtsM3 + CritPoints->NumCrtPtsM1
                                                    + CritPoints->NumCrtPtsP1;
                            IsOk =  CritPointsInsertPoint(CritPoints, PointOffset, XCrtPt, YCrtPt, ZCrtPt,
                                                          ChrgDensCrtPt, (char)TypeCrtPt, PrincDirX, PrincDirY, PrincDirZ);
                        }
                        break;
                        case 3:
                        {
                            LgIndex_t PointOffset = CritPoints->NumCrtPtsM3 + CritPoints->NumCrtPtsM1
                                                    + CritPoints->NumCrtPtsP1 + CritPoints->NumCrtPtsP3;
                            IsOk =  CritPointsInsertPoint(CritPoints, PointOffset, XCrtPt, YCrtPt, ZCrtPt,
                                                          ChrgDensCrtPt, (char)TypeCrtPt, 0.0, 0.0, 0.0);
                        }
                        break;
                    }
                }
            }
        }

        /* Update status bar */
        if (ShowStatusBar)
        {
            char PercentDoneText[45];
            int  PercentDone;
            PercentDone = (int)((100 * kk) / KMax);
            sprintf_s(PercentDoneText, "Finding Critical Points: %d Percent Done", PercentDone);
            TecUtilStatusSetPercentDoneText(PercentDoneText);
            if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
            if (IsStop) IsOk = FALSE;
        }
    }


    /* Eliminate duplicate Cage critical points. */
    if (IsOk && CritPoints->NumCrtPtsP3 > 1)
        IsOk = RemDuplicateCritPoints(TRUE, 0.4, 3, CritPoints);
    if (IsOk && CritPoints->NumCrtPtsP3 > 1)
        IsOk = RemDuplicateCritPoints(FALSE, 1.0, 3, CritPoints);

    /* Eliminate duplicate Ring critical points. */
    if (IsOk && CritPoints->NumCrtPtsP1 > 1)
        IsOk = RemDuplicateCritPoints(FALSE, 1.0, 1, CritPoints);

    /* Eliminate duplicate Bond critical points. */
    if (IsOk && CritPoints->NumCrtPtsM1 > 1)
        IsOk = RemDuplicateCritPoints(FALSE, 1.0, -1, CritPoints);

    /* Eliminate duplicate Atom critical points. */
    if (IsOk && CritPoints->NumCrtPtsM3 > 1)
        IsOk = RemDuplicateCritPoints(TRUE, 0.4, -3, CritPoints);
    if (IsOk && CritPoints->NumCrtPtsM3 > 1)
        IsOk = RemDuplicateCritPoints(FALSE, 1.0, -3, CritPoints);

    if (ShowStatusBar) TecUtilStatusFinishPercentDone();

    /* Inform Tecplot that major data operation is ending */
    TecUtilDataLoadEnd();


    ENSURE(VALID_BOOLEAN(IsOk));
    return(IsOk);
}






/**
 * Create a Tecplot 1-D Scatter-point zone containing the critical
 * points components to CritPoints.
 *
 * param CritPoints
 *     CritPoints structure to be written to Tecplot Zone.
 *
 * param ChrgDensVarNum, TypeVarNum, UVarNum, VVarNum, WVarNum
 *     Variable numbers for the charge density, critical point type,
 *     and X-, Y-, and Z-Components of the charge density gradient
 *     vector.
 *
 * NOTE: For Bond and Ring critical points, this function stores the
 *       components of the principle direction in the UVarNum, VVarNum,
 *       and WVarNum variables.
 *
 * return
 *     Number of Tecplot zone if successful, otherwise 0.
 */
EntIndex_t CritPointsCreateTPZone(CritPoints_pa  CritPoints,
                                  EntIndex_t     ChrgDensVarNum,
                                  EntIndex_t     TypeVarNum,
                                  EntIndex_t     UVarNum,
                                  EntIndex_t     VVarNum,
                                  EntIndex_t     WVarNum)
{
    EntIndex_t CPZoneNum = TECUTILSETNOTMEMBER;
    Boolean_t  IsOk = TRUE;
    EntIndex_t NumZones = 0;
    EntIndex_t NumVars  = 0;
    LgIndex_t  NumCrtPts = CritPointsGetCount(CritPoints);

    REQUIRE(CritPointsIsValid(CritPoints));

    /* Get num zones and num vars in dataset */
    if (IsOk)
    {
        char *DatasetTitle = NULL;
        TecUtilDataSetGetInfo(&DatasetTitle, &NumZones, &NumVars);
        TecUtilStringDealloc(&DatasetTitle);
    }

    REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
    REQUIRE(TypeVarNum > 0 && TypeVarNum <= NumVars);
    REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
    REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
    REQUIRE(WVarNum > 0 && WVarNum <= NumVars);

    /* Inform Tecplot that major data operation is beginning */
    TecUtilDataLoadBegin();

    /* Create Critical Point Zone */
    if (IsOk && NumCrtPts > 0)
    {
        char       ZoneName[100];
        EntIndex_t SourceZoneNum = CritPoints->SourceZoneNum;

        sprintf_s(ZoneName, "Critical Points Zone %d", SourceZoneNum);

        // AveZoneNum = TUZoneGetNumByName(ZoneName);

        /* Set FieldDataType of all variables of CP zone */
        FieldDataType_e  *VarDataType = ALLOC_ARRAY(NumVars, FieldDataType_e, VarDataType);
        for (EntIndex_t iv = 0; iv < NumVars; iv++)
            VarDataType[iv] = FieldDataType_Double;

        /* Create zone if it doesn't already exist. */
        if (CPZoneNum == TECUTILSETNOTMEMBER)
        {
            ArgList_pa ArgList;
            TecUtilLockStart(AddOnID);
            ArgList = TecUtilArgListAlloc();
            TecUtilArgListAppendString(ArgList, SV_NAME, ZoneName);

            TecUtilArgListAppendInt(ArgList, SV_ZONETYPE,
                                    (ArbParam_t)ZoneType_Ordered);

            TecUtilArgListAppendInt(ArgList, SV_IMAX, NumCrtPts);
            TecUtilArgListAppendInt(ArgList, SV_JMAX, 1);
            TecUtilArgListAppendInt(ArgList, SV_KMAX, 1);

            TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)VarDataType);

            IsOk = TecUtilDataSetAddZoneX(ArgList);

            TecUtilArgListDealloc(&ArgList);
            TecUtilLockFinish(AddOnID);
            if (IsOk)
            {
                NumZones++;
                CPZoneNum = NumZones;
            }

            /* Inform Tecplot of new zone. */
            if (IsOk)
            {
                Set_pa ZSet = TecUtilSetAlloc(TRUE);
                TecUtilSetAddMember(ZSet, NumZones, TRUE);
                TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)ZSet);
                TecUtilSetDealloc(&ZSet);
            }
            if (VarDataType != NULL)
                FREE_ARRAY(VarDataType, "VarDataType");
        }

        /* Set Critical Point Locations into zone */
        if (IsOk)
        {
            FieldData_pa XVarFDPtr, YVarFDPtr, ZVarFDPtr;
            FieldData_pa UVarFDPtr, VVarFDPtr, WVarFDPtr;
            FieldData_pa CDVarFDPtr = NULL;
            FieldData_pa TypeVarFDPtr = NULL;
            LgIndex_t    ii, IJunk, JJunk, KJunk;

            TecUtilLockStart(AddOnID);

            TecUtilZoneGetInfo(CPZoneNum, &IJunk, &JJunk, &KJunk,  &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            CDVarFDPtr   = TecUtilDataValueGetWritableRef(CPZoneNum, ChrgDensVarNum);
            TypeVarFDPtr = TecUtilDataValueGetWritableRef(CPZoneNum, TypeVarNum);
            UVarFDPtr    = TecUtilDataValueGetWritableRef(CPZoneNum, UVarNum);
            VVarFDPtr    = TecUtilDataValueGetWritableRef(CPZoneNum, VVarNum);
            WVarFDPtr    = TecUtilDataValueGetWritableRef(CPZoneNum, WVarNum);

            for (ii = 0; IsOk && ii < NumCrtPts; ii++)
            {
                double XCrtPt, YCrtPt, ZCrtPt, ChrgDensCrtPt, PrincDirX, PrincDirY, PrincDirZ;
                char   TypeCrtPt;

                IsOk = CritPointsGetPoint(CritPoints, ii, &XCrtPt, &YCrtPt, &ZCrtPt,
                                          &ChrgDensCrtPt, &TypeCrtPt, &PrincDirX, &PrincDirY, &PrincDirZ);

                TecUtilDataValueSetByRef(XVarFDPtr, ii + 1, XCrtPt);
                TecUtilDataValueSetByRef(YVarFDPtr, ii + 1, YCrtPt);
                TecUtilDataValueSetByRef(ZVarFDPtr, ii + 1, ZCrtPt);
                TecUtilDataValueSetByRef(CDVarFDPtr, ii + 1, ChrgDensCrtPt);
                TecUtilDataValueSetByRef(TypeVarFDPtr, ii + 1, TypeCrtPt);
                // Store principle direction components in U, V, W vars
                TecUtilDataValueSetByRef(UVarFDPtr, ii + 1, PrincDirX);
                TecUtilDataValueSetByRef(VVarFDPtr, ii + 1, PrincDirY);
                TecUtilDataValueSetByRef(WVarFDPtr, ii + 1, PrincDirZ);
            }

            TecUtilLockFinish(AddOnID);
        }

        /* Set Zone aux data to completely transfer necessary information to Tecplot */
        if (IsOk)
        {
            AuxData_pa ZnAuxDataRef = TecUtilAuxDataZoneGetRef(CPZoneNum);
            if (ZnAuxDataRef != NULL)
            {
                char NumCrtPtStr[10], BaseZoneNumStr[10];

                IsOk = TecUtilAuxDataSetStrItem(ZnAuxDataRef, "CompChem.ZoneType", "CriticalPoints", TRUE);

                sprintf_s(NumCrtPtStr, "%d", CritPoints->NumCrtPtsM3);
                IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.NumCrtPtAtom", (ArbParam_t)NumCrtPtStr,
                                             AuxDataType_String, TRUE);

                sprintf_s(NumCrtPtStr, "%d", CritPoints->NumCrtPtsM1);
                if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.NumCrtPtBond",
                                                           (ArbParam_t)NumCrtPtStr, AuxDataType_String, TRUE);

                sprintf_s(NumCrtPtStr, "%d", CritPoints->NumCrtPtsP1);
                if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.NumCrtPtRing",
                                                           (ArbParam_t)NumCrtPtStr, AuxDataType_String, TRUE);

                sprintf_s(NumCrtPtStr, "%d", CritPoints->NumCrtPtsP3);
                if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.NumCrtPtCage",
                                                           (ArbParam_t)NumCrtPtStr, AuxDataType_String, TRUE);

                sprintf_s(BaseZoneNumStr, "%d", SourceZoneNum);
                if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BaseZoneNum",
                                                           (ArbParam_t)BaseZoneNumStr, AuxDataType_String, TRUE);
            }
            else IsOk = FALSE;
        }

    }

    /* Inform Tecplot that major data operation is ending */
    TecUtilDataLoadEnd();


    if (IsOk == FALSE) CPZoneNum = 0;


    ENSURE(VALID_BOOLEAN(IsOk));

    ENSURE(CPZoneNum >= 0 && CPZoneNum <= NumZones);
    return CPZoneNum;
}





/**
 * Create and populate a CritPoints structure using data from a
 * CompChem.ZoneType=CriticalPoints Tecplot zone.
 *
 *
 * param TPZoneNum
 *     Tecplot zone number of the CompChem.ZoneType=CriticalPoints zone.
 *
 * param ChrgDensVarNum, TypeVarNum, UVarNum, VVarNum, WVarNum
 *     Variable numbers for the charge density, critical point type,
 *     and X-, Y-, and Z-Components of the charge density gradient
 *     vector.
 *
 * NOTE: For Bond and Ring critical points, the components of the
 *       principle direction are stored in the UVarNum, VVarNum,
 *       and WVarNum variables.
 *
 * return
 *     CritPoints handle if successful, otherwise NULL.
 */
CritPoints_pa CritPointsGetFromTPZone(EntIndex_t TPZoneNum,
                                      EntIndex_t ChrgDensVarNum,
                                      EntIndex_t TypeVarNum,
                                      EntIndex_t UVarNum,
                                      EntIndex_t VVarNum,
                                      EntIndex_t WVarNum)
{
    CritPoints_pa Result = NULL;
    Boolean_t     IsOk   = TRUE;
    EntIndex_t    NumZones, NumVars;

    /* Get num zones and num vars in dataset */
    if (IsOk)
    {
        char *DatasetTitle = NULL;
        TecUtilDataSetGetInfo(&DatasetTitle, &NumZones, &NumVars);
        TecUtilStringDealloc(&DatasetTitle);
    }

    REQUIRE(TPZoneNum > 0 && TPZoneNum <= NumZones);
    REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
    REQUIRE(TypeVarNum > 0 && TypeVarNum <= NumVars);
    REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
    REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
    REQUIRE(WVarNum > 0 && WVarNum <= NumVars);

    /* Verify that TPZoneNum is a CriticalPoints zones */
    if (IsOk)
    {
        ArbParam_t Value;
        AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(TPZoneNum);
        if (AuxDataRef != NULL)
        {
            AuxDataType_e Type;
            Boolean_t     Retain;
            if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.ZoneType",
                                            &Value, &Type, &Retain))
            {
                if (Type == AuxDataType_String) // currently the only type supported
                {
                    char *ValueString = (char *)Value;
                    if (strcmp("CriticalPoints", ValueString) != 0) IsOk = FALSE;
                    TecUtilStringDealloc(&ValueString);
                }
                else
                    IsOk = FALSE;
            }
            else
                IsOk = FALSE;
        }
        else
            IsOk = FALSE;
        if (IsOk == FALSE)
            TecUtilDialogErrMsg("Incorrectly specified CriticalPoints zone.");
    }


    if (IsOk)
    {
        Result = CritPointsAlloc();
        if (Result != NULL)
        {

            // Extract the source zone number
            ArbParam_t Value;
            AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(TPZoneNum);
            if (AuxDataRef != NULL)
            {
                AuxDataType_e Type;
                Boolean_t     Retain;

                // Extract number of Atoms from Aux Data
                if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.BaseZoneNum",
                                                &Value, &Type, &Retain))
                {
                    if (Type == AuxDataType_String) // currently the only type supported
                    {
                        char      *ValueString = (char *)Value;
                        Result->SourceZoneNum = atoi(ValueString);
                        TecUtilStringDealloc(&ValueString);
                    }
                    else
                        IsOk = FALSE;
                }
                else
                    IsOk = FALSE;
            }
        }
        /*
              if (IsOk && Result != NULL)
                {

                  // Extract the aux data and set the cooresponding CritPoints structure elements
                  ArbParam_t Value;
                  AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(TPZoneNum );
                  if (AuxDataRef != NULL)
                    {
                      AuxDataType_e Type;
                      Boolean_t     Retain;

                      // Extract number of Atoms from Aux Data
                      if (TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.NumCrtPtAtom",
                                                   &Value, &Type, &Retain))
                        {
                          if (Type == AuxDataType_String) // currently the only type supported
                            {
                              char      *ValueString = (char *)Value;
                              Result->NumCrtPtsM3 = atoi(ValueString);
                              TecUtilStringDealloc(&ValueString);
                            }
                          else
                            IsOk = FALSE;
                        }
                      else
                        IsOk = FALSE;

                      // Extract number of Bonds from Aux Data
                      if (IsOk && TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.NumCrtPtBond",
                                                   &Value, &Type, &Retain))
                        {
                          if (Type == AuxDataType_String) // currently the only type supported
                            {
                              char      *ValueString = (char *)Value;
                              Result->NumCrtPtsM1 = atoi(ValueString);
                              TecUtilStringDealloc(&ValueString);
                            }
                          else
                            IsOk = FALSE;
                        }
                      else
                        IsOk = FALSE;

                      // Extract number of Rings from Aux Data
                      if (IsOk && TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.NumCrtPtRing",
                                                   &Value, &Type, &Retain))
                        {
                          if (Type == AuxDataType_String) // currently the only type supported
                            {
                              char      *ValueString = (char *)Value;
                              Result->NumCrtPtsP1 = atoi(ValueString);
                              TecUtilStringDealloc(&ValueString);
                            }
                          else
                            IsOk = FALSE;
                        }
                      else
                        IsOk = FALSE;

                      // Extract number of Cages from Aux Data
                      if (IsOk && TecUtilAuxDataGetItemByName(AuxDataRef, "CompChem.NumCrtPtCage",
                                                   &Value, &Type, &Retain))
                        {
                          if (Type == AuxDataType_String) // currently the only type supported
                            {
                              char      *ValueString = (char *)Value;
                              Result->NumCrtPtsP3 = atoi(ValueString);
                              TecUtilStringDealloc(&ValueString);
                            }
                          else
                            IsOk = FALSE;
                        }
                      else
                        IsOk = FALSE;


                      // Compute and set number of critical points.
                      if (IsOk) Result->NumCrtPts = Result->NumCrtPtsM3
                                                  + Result->NumCrtPtsM1
                                                  + Result->NumCrtPtsP1
                                                  + Result->NumCrtPtsP3;

                    }
                }
                */
    }


    /* Extract Critical Points from zone and set in CritPoints structure arrays */
    if (IsOk)
    {
        FieldData_pa XVarFDPtr, YVarFDPtr, ZVarFDPtr;
        FieldData_pa UVarFDPtr, VVarFDPtr, WVarFDPtr;
        FieldData_pa CDVarFDPtr = NULL;
        FieldData_pa TypeVarFDPtr = NULL;
        LgIndex_t    ii, IMax, JJunk, KJunk;

        TecUtilLockStart(AddOnID);

        TecUtilZoneGetInfo(TPZoneNum, &IMax, &JJunk, &KJunk,  &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
                           NULL, NULL, NULL, NULL, NULL, NULL, NULL);

        // CHECK(Result->NumCrtPts == IMax);

        CDVarFDPtr   = TecUtilDataValueGetReadableRef(TPZoneNum, ChrgDensVarNum);
        TypeVarFDPtr = TecUtilDataValueGetReadableRef(TPZoneNum, TypeVarNum);
        UVarFDPtr    = TecUtilDataValueGetReadableRef(TPZoneNum, UVarNum);
        VVarFDPtr    = TecUtilDataValueGetReadableRef(TPZoneNum, VVarNum);
        WVarFDPtr    = TecUtilDataValueGetReadableRef(TPZoneNum, WVarNum);

        for (ii = 0; IsOk && ii < IMax; ii++)
        {
            double XCrtPt, YCrtPt, ZCrtPt, ChrgDensCrtPt, PrincDirX, PrincDirY, PrincDirZ;
            char   TypeCrtPt;

            XCrtPt = TecUtilDataValueGetByRef(XVarFDPtr, ii + 1);
            YCrtPt = TecUtilDataValueGetByRef(YVarFDPtr, ii + 1);
            ZCrtPt = TecUtilDataValueGetByRef(ZVarFDPtr, ii + 1);
            ChrgDensCrtPt = TecUtilDataValueGetByRef(CDVarFDPtr, ii + 1);
            TypeCrtPt = (char)TecUtilDataValueGetByRef(TypeVarFDPtr, ii + 1);
            // Principle direction components stored in U, V, W vars
            PrincDirX = TecUtilDataValueGetByRef(UVarFDPtr, ii + 1);
            PrincDirY = TecUtilDataValueGetByRef(VVarFDPtr, ii + 1);
            PrincDirZ = TecUtilDataValueGetByRef(WVarFDPtr, ii + 1);

            IsOk = CritPointsAppendPoint(Result, XCrtPt, YCrtPt, ZCrtPt,
                                         ChrgDensCrtPt, TypeCrtPt, PrincDirX, PrincDirY, PrincDirZ);
        }

        TecUtilLockFinish(AddOnID);
    }


    /* Notify Tecplot of memory usage change */
    if (IsOk)
    {
        Result->MemUsageReported = (30 + 60 * Result->NumCrtPts) / 1024;
        if (Result->MemUsageReported > 0)
            TecUtilMemoryChangeNotify((Int64_t)(Result->MemUsageReported));
    }

    if (!IsOk)
    {
        CritPointsDealloc(&Result);
        Result = NULL;
    }


    ENSURE(CritPointsIsValid(Result) || Result == NULL);
    return Result;
}

