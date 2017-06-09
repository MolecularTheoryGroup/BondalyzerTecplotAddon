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
#include "LINEARALGEBRA.h"
#include "ELEMSHAPEFUNC.h"
#include <string.h>

#define MAXNEWTONITER 100

static double TetraTolerance = 1.0005;
// static double RangeEpsilon = 0.0;

/*
 * RECTANGULAR GRID WITH (X,Y,Z) ALIGNED WITH (I,J,K) ONLY. 
 * For given coordinates (X,Y,Z) compute the cell indices (i,j,k) 
 * and the natural coordinates (r,s,t) within the cell.
 *
 * Input:
 *    X,Y,Z: Physical coordinates of a point (perhaps) in the zone.
 * Input:
 *    [XYZ][BegEnd]Zone: Coordinate ranges of the rectangular zone.
 * Input:
 *    PeriodicBC: Can exceed XYZ limits of zone, repeat IJK/xyz sequence.
 * Output:
 *    i,j,k: Indices of cell containing the point.
 *    r,s,t: Natural coordinates - ranges is from -1 to 1 and cell
 *           center is at 0.
 *
 *    Returns: TRUE if point (X,Y,Z) is in the zone, FALSE otherwise.
 */
Boolean_t RectGridBrickXYZtoIJKRST(double          X,
								   double          Y,
								   double          Z,
								   LgIndex_t       IMax,
								   LgIndex_t       JMax,
								   LgIndex_t       KMax,
								   double          XBegZone,
								   double          XEndZone,
								   double          YBegZone,
								   double          YEndZone,
								   double          ZBegZone,
								   double          ZEndZone,
								   Boolean_t       PeriodicBC,
								   LgIndex_t      *i,
								   LgIndex_t      *j,
								   LgIndex_t      *k,
								   double         *r,
								   double         *s,
								   double         *t)
{
	Boolean_t IsInZone = TRUE;

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1);
	REQUIRE(XEndZone > XBegZone && YEndZone > YBegZone && ZEndZone > ZBegZone);
	REQUIRE(VALID_BOOLEAN(PeriodicBC));
	REQUIRE(VALID_REF(i));
	REQUIRE(VALID_REF(j));
	REQUIRE(VALID_REF(k));
	REQUIRE(VALID_REF(r));
	REQUIRE(VALID_REF(s));
	REQUIRE(VALID_REF(t));

	double DXZoneBuffer = 0.0;
	double DYZoneBuffer = 0.0;
	double DZZoneBuffer = 0.0;
	LgIndex_t ZoneCoordBuffer = 0;


	// Compute buffer, assuming XYZ are aligned with IJK respectively
	if (PeriodicBC)
	{
		DXZoneBuffer = (XEndZone - XBegZone) / (double)(MAX(IMax-1, 1));
		DYZoneBuffer = (YEndZone - YBegZone) / (double)(MAX(JMax-1, 1));
		DZZoneBuffer = (ZEndZone - ZBegZone) / (double)(MAX(KMax-1, 1));
		ZoneCoordBuffer = 1;
	}

	/* Is it in the solution domain */
	if (X < XBegZone - DXZoneBuffer || X > XEndZone + DXZoneBuffer ||
		Y < YBegZone - DYZoneBuffer || Y > YEndZone + DYZoneBuffer ||
		Z < ZBegZone - DZZoneBuffer || Z > ZEndZone + DZZoneBuffer ) IsInZone = FALSE;

	// Compute the cell and natural coordinates. Assume equally-spaced rectangular-
	// grid where X,Y,Z is aligned with I,J,K.
	if (IsInZone)
	{
		double ICoordPos = 1.0 + (double)(IMax - 1) * (X - XBegZone) / (XEndZone - XBegZone);
		double JCoordPos = 1.0 + (double)(JMax - 1) * (Y - YBegZone) / (YEndZone - YBegZone);
		double KCoordPos = 1.0 + (double)(KMax - 1) * (Z - ZBegZone) / (ZEndZone - ZBegZone);

		*i = MAX(MIN((LgIndex_t)ICoordPos, IMax - 1 + ZoneCoordBuffer), 1 - ZoneCoordBuffer);
		*j = MAX(MIN((LgIndex_t)JCoordPos, JMax - 1 + ZoneCoordBuffer), 1 - ZoneCoordBuffer);
		*k = MAX(MIN((LgIndex_t)KCoordPos, KMax - 1 + ZoneCoordBuffer), 1 - ZoneCoordBuffer);

		*r = 2.0 * (ICoordPos - (double)(*i)) - 1.0;
		*s = 2.0 * (JCoordPos - (double)(*j)) - 1.0;
		*t = 2.0 * (KCoordPos - (double)(*k)) - 1.0;
	}

	ENSURE(!IsInZone || (*i >= 1 - ZoneCoordBuffer && *i <= IMax + ZoneCoordBuffer));
	ENSURE(!IsInZone || (*j >= 1 - ZoneCoordBuffer && *j <= JMax + ZoneCoordBuffer));
	ENSURE(!IsInZone || (*k >= 1 - ZoneCoordBuffer && *k <= KMax + ZoneCoordBuffer));

	ENSURE(!IsInZone || (*r >= -1.0 && *r <= 1.0));
	ENSURE(!IsInZone || (*s >= -1.0 && *s <= 1.0));
	ENSURE(!IsInZone || (*t >= -1.0 && *t <= 1.0));

	ENSURE(VALID_BOOLEAN(IsInZone));
	return IsInZone;
}





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
	// TIM: For ChemBond tool, X,Y,Z can be replaced with U,V,W respectively,
	// and represent the values of the gradient.
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
	// TIM: Can be any X,Y,Z because CPs are CPs so long as at least 1 dimension has a CP in cell.
	// What type of CP can be determined later based on the eigenvalues at the point.
	if (X <= XCMax || X >= XCMin || Y <= YCMax || Y >= YCMin || Z <= ZCMax || Z >= ZCMin)
	{
		Boolean_t IsOk = TRUE;
		Boolean_t IsInside;
		Boolean_t WellOutside = FALSE;
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

		/* Set initial values of r,s,t - if out of range. */
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
		while ((IsOk && !WellOutside && DeltaSqr > TolerSqr && iter < MAXNEWTONITER) ||
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
				double dr = DeltaXYZ[0];
				double ds = DeltaXYZ[1];
				double dt = DeltaXYZ[2];

				// Compute an underrelaxation factor that will take you (almost because of 0.999 factor) to the boundary of the cell
				double UnderRelax = 1.0;
				/*
				if (*r + dr >  1.0) UnderRelax = MIN(MAX((1.0 - *r) / MAX(ABS(dr), 0.01), 0.1), 1.0);
				if (*r + dr < -1.0) UnderRelax = MIN(MAX((1.0 + *r) / MAX(ABS(dr), 0.01), 0.1), 1.0); 
				if (*s + ds >  1.0) UnderRelax = MIN(MAX((1.0 - *s) / MAX(ABS(ds), 0.01), 0.1), 1.0);
				if (*s + ds < -1.0) UnderRelax = MIN(MAX((1.0 + *s) / MAX(ABS(ds), 0.01), 0.1), 1.0); 
				if (*t + dt >  1.0) UnderRelax = MIN(MAX((1.0 - *t) / MAX(ABS(dt), 0.01), 0.1), 1.0);
				if (*t + dt < -1.0) UnderRelax = MIN(MAX((1.0 + *t) / MAX(ABS(dt), 0.01), 0.1), 1.0); 
				*/

				if (*r + dr >  1.0) UnderRelax = MIN(MAX((0.999 - *r) / MAX(ABS(dr), 0.01), 0.1), 1.0);
				if (*r + dr < -1.0) UnderRelax = MIN(MAX((0.999 + *r) / MAX(ABS(dr), 0.01), 0.1), 1.0); 
				if (*s + ds >  1.0) UnderRelax = MIN(MAX((0.999 - *s) / MAX(ABS(ds), 0.01), 0.1), 1.0);
				if (*s + ds < -1.0) UnderRelax = MIN(MAX((0.999 + *s) / MAX(ABS(ds), 0.01), 0.1), 1.0); 
				if (*t + dt >  1.0) UnderRelax = MIN(MAX((0.999 - *t) / MAX(ABS(dt), 0.01), 0.1), 1.0);
				if (*t + dt < -1.0) UnderRelax = MIN(MAX((0.999 + *t) / MAX(ABS(dt), 0.01), 0.1), 1.0); 

				*r = *r + UnderRelax * dr;
				*s = *s + UnderRelax * ds;
				*t = *t + UnderRelax * dt;
				/*
				*r = *r + 0.5 * DeltaXYZ[0];
				*s = *s + 0.5 * DeltaXYZ[1];
				*t = *t + 0.5 * DeltaXYZ[2];
				*/
			}
			else{
				CellFound = FALSE;
				return CellFound;
			}

			/* Is it well outside the cell? */
			if (*r < -200.0 || *r > 200.0 || *s < -200.0 || *s > 200.0 || *t < -200.0 || *t > 200.0)
			{
				WellOutside = TRUE;
			}

			if (!WellOutside)
			{
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
			}

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
 * Compute the interpolation weights for a bilinear quad.
 * Input:
 *    r,s: Natural coordinates - ranges is from -1 to 1 and cell
 *         center is at 0.
 * Output:
 *    W[]:   Averaging weights for cell (duplicate nodes removed).
 *           At this stage, only the weights for the face are set.
 *    Returns: True if point (r,s) is in the quad.
 */
Boolean_t QuadBilinearWeight(double       r,
							 double       s,
							 double       W[4])
{
	Boolean_t CellFound = FALSE;
	REQUIRE(VALID_REF(W));

	if (r > -TetraTolerance && r < TetraTolerance &&
		s > -TetraTolerance && s < TetraTolerance)
	{
		double    OneMinusR, OnePlusR;
		double    OneMinusS, OnePlusS;

		CellFound = TRUE;
		OneMinusR = 1.0 - r;
		OnePlusR  = 1.0 + r;
		OneMinusS = 1.0 - s;
		OnePlusS  = 1.0 + s;
		W[0] = 0.25 * OneMinusR * OneMinusS;
		W[1] = 0.25 * OnePlusR  * OneMinusS;
		W[2] = 0.25 * OnePlusR  * OnePlusS;
		W[3] = 0.25 * OneMinusR * OnePlusS;
	}

	ENSURE(VALID_BOOLEAN(CellFound));
	return (CellFound);
}


/*
 * Compute the derivatives wrt the natural coordinates of the
 * interpolation weights for a bilinear quad.
 * Input:
 *    r,s: Natural coordinates - ranges is from -1 to 1 and cell
 *         center is at 0.
 * Output:
 *    W[]:    Interpolation weights.
 *    dwdr[]: Derivative of weights wrt r (first natural coordinate).
 *    dwds[]: Derivative of weights wrt s (second natural coordinate).
 *
 *    Returns: True if point (r,s) is in the quad.
 */
Boolean_t QuadBilinearWeightDerivatives(double       r,
										double       s,
										double       W[4],
										double       dwdr[4],
										double       dwds[4])
{
	Boolean_t CellFound = FALSE;
	REQUIRE(VALID_REF(W));
	REQUIRE(VALID_REF(dwdr));
	REQUIRE(VALID_REF(dwds));

	if (r > -TetraTolerance && r < TetraTolerance &&
		s > -TetraTolerance && s < TetraTolerance) CellFound = TRUE;

	{
		double    OneMinusR, OnePlusR;
		double    OneMinusS, OnePlusS;

		OneMinusR = 1.0 - r;
		OnePlusR  = 1.0 + r;
		OneMinusS = 1.0 - s;
		OnePlusS  = 1.0 + s;
		W[0] = 0.25 * OneMinusR * OneMinusS;
		W[1] = 0.25 * OnePlusR  * OneMinusS;
		W[2] = 0.25 * OnePlusR  * OnePlusS;
		W[3] = 0.25 * OneMinusR * OnePlusS;
		dwdr[0] = -0.25 * OneMinusS;
		dwdr[1] =  0.25 * OneMinusS;
		dwdr[2] =  0.25 * OnePlusS;
		dwdr[3] = -0.25 * OnePlusS;
		dwds[0] = -0.25 * OneMinusR;
		dwds[1] = -0.25 * OnePlusR;
		dwds[2] =  0.25 * OnePlusR;
		dwds[3] =  0.25 * OneMinusR;
	}

	ENSURE(VALID_BOOLEAN(CellFound));
	return (CellFound);
}






/*
 * Compute the linear natural coordinates (r,s) corresponding to a vector
 * solution (U,V) in the 2D triangle cell. If it is a 3D surface, (X,Y)
 * should be single-values coordinates - perhaps a surface tangent
 * coordinate system.
 * Input:
 *    U,V: Vector solution at the desired point (perhaps) in the cell.
 *    UCell,VCell: Vector solution at the nodes.
 *    UCMin, UCMax, VCMin, VCMax: Coordinate range of cell to use in preliminary test.
 *    XCell,YCell: Physical surface-tangent coordinates of the nodes.
 * Output:
 *    X,Y: X,Y cordinates at which U,V is found.
 *    W[]: Interpolation weights.
 *
 *    Returns: True if point with vector (U,V) is in the linear triangle.
 */
Boolean_t LinearTriangleNaturalCoord(double       U,
									 double       V,
									 double       UCell[3],
									 double       VCell[3],
									 double       UCMin,
									 double       VCMin,
									 double       UCMax,
									 double       VCMax,
									 double       XCell[3],
									 double       YCell[3],
									 double       *X,
									 double       *Y,
									 double       W[3])
{
	Boolean_t CellFound = FALSE;
	double x, y;

	REQUIRE(VALID_REF(W));
	REQUIRE(VALID_REF(UCell));
	REQUIRE(VALID_REF(VCell));
	REQUIRE(VALID_REF(XCell));
	REQUIRE(VALID_REF(YCell));
	REQUIRE(VALID_REF(X));
	REQUIRE(VALID_REF(Y));

	/*
	 * Solve the following linear system of equations for x,y
	 *
	 *  U = L1(x,y) * UCell[0] + L2(x,y) * UCell[1] + L3(x,y) * UCell[2]
	 *  V = L1(x,y) * VCell[0] + L2(x,y) * VCell[1] + L3(x,y) * VCell[2]
	 *
	 * Where
	 *   L1(x,y) = (a1 + b1 * x + c1 * y) / (2 * Area)
	 *   L2(x,y) = (a2 + b2 * x + c2 * y) / (2 * Area)
	 *   L3(x,y) = (a3 + b3 * x + c3 * y) / (2 * Area)
	 *
	 * Resulting equations  M X = B
	 *    _                                                              _   _ _        _                                            _
	 *   | (b1 * U1 + b2 * U2 + b3 * U3)    (c1 * U1 + c2 * U2 + c3 * U3) | | x |      | U * 2 * Area - (a1 * U1 + a2 * U2 + a3 * U3) |
	 *   |                                                                | |   |  =   |                                              |
	 *   |_(b1 * V1 + b2 * V2 + b3 * V3)    (c1 * V1 + c2 * V2 + c3 * V3)_| |_y_|      |_V * 2 * Area - (a1 * V1 + a2 * V2 + a3 * V3)_|
	 *
	 */
	if (U <= UCMax && U >= UCMin && V <= VCMax && V >= VCMin)
	{
		Boolean_t IsOk = TRUE;

		double L1, L2, L3;

		double a1 = XCell[1] * YCell[2] - XCell[2] * YCell[1];
		double b1 = YCell[1] - YCell[2];
		double c1 = XCell[2] - XCell[1];

		double a2 = XCell[2] * YCell[0] - XCell[0] * YCell[2];
		double b2 = YCell[2] - YCell[0];
		double c2 = XCell[0] - XCell[2];

		double a3 = XCell[0] * YCell[1] - XCell[1] * YCell[0];
		double b3 = YCell[0] - YCell[1];
		double c3 = XCell[1] - XCell[0];

		double TwoDelta = XCell[1] * YCell[2] - XCell[2] * YCell[1]
						  + XCell[0] * (YCell[1] - YCell[2])
						  + YCell[0] * (XCell[2] - XCell[1]);

		// Compute the components of the system matrix M
		double M11 = b1 * UCell[0] + b2 * UCell[1] + b3 * UCell[2];
		double M12 = c1 * UCell[0] + c2 * UCell[1] + c3 * UCell[2];
		double M21 = b1 * VCell[0] + b2 * VCell[1] + b3 * VCell[2];
		double M22 = c1 * VCell[0] + c2 * VCell[1] + c3 * VCell[2];
		double Det = M11 * M22 - M21 * M12;

		// Compute the components of the RHS vector B
		double B1 = U * TwoDelta - (a1 * UCell[0] + a2 * UCell[1] + a3 * UCell[2]);
		double B2 = V * TwoDelta - (a1 * VCell[0] + a2 * VCell[1] + a3 * VCell[2]);

		// Solve system (indeterminant if Det==0.0)
		if (IsOk && ABS(Det) > SMALLFLOAT)
		{
			double RDet = 1.0 / Det;
			x = RDet * (M22 * B1 - M12 * B2);
			y = RDet * (-M21 * B1 + M11 * B2);

			L1 = (a1 + b1 * x + c1 * y) / TwoDelta;
			L2 = (a2 + b2 * x + c2 * y) / TwoDelta;
			L3 = 1.0 - L1 - L2;

			if (L1 >= 0.0 && L1 <= 1.0  &&  L2 >= 0.0 && L2 <= 1.0  &&
				L3 >= 0.0 && L3 <= 1.0)
			{
				CellFound = TRUE;

				*X = x;
				*Y = y;

				W[0] = L1;
				W[1] = L2;
				W[2] = L3;
			}
		}
	}

	ENSURE(VALID_BOOLEAN(CellFound));
	return (CellFound);
}





/*
 * Compute the bilinear natural coordinates (r,s) corresponding
 * to a point(X,Y) in the 2D cell. If it is a 3D surface, (X,Y) should
 * be single-values coordinates - perhaps a surface tangent coordinate
 * system.
 * Input:
 *    X,Y: Physical surface-tangent coordinates of a point (perhaps) in the cell.
 *    XCMin, XCMax, YCMin, YCMax: Coordinate range of cell to use in preliminary test.
 * Output:
 *    r,s: Natural coordinates - ranges is from -1 to 1 and cell
 *           center is at 0.
 *    W[]:    Interpolation weights.
 *
 *    Returns: True if point (X,Y) is in the bilinear brick.
 */
Boolean_t BilinearQuadNaturalCoord(double       X,
								   double       Y,
								   double       XCell[4],
								   double       YCell[4],
								   double       XCMin,
								   double       YCMin,
								   double       XCMax,
								   double       YCMax,
								   double       *r,
								   double       *s,
								   double       W[4])
{
	Boolean_t CellFound = FALSE;
	double dwdr[4];
	double dwds[4];
	REQUIRE(VALID_REF(W));
	REQUIRE(VALID_REF(XCell));
	REQUIRE(VALID_REF(YCell));
	REQUIRE(VALID_REF(r));
	REQUIRE(VALID_REF(s));

	/*
	 * Use a Newton iteration to solve the non-linear equations for
	 * r,s in terms of X,Y.
	 */
	if (X <= XCMax && X >= XCMin && Y <= YCMax && Y >= YCMin)
	{
		Boolean_t IsOk = TRUE;
		Boolean_t IsInside;
		int       ic;
		int       iter = 0;
		double    Xm, Ym, DeltaX, DeltaY, Toler, TolerSqr, DeltaSqr;

		/* If computing derivatives of weights, must pass through while loop once */
		Boolean_t FirstPass = TRUE;

		/* Matrix and temporary array for LU decomposition */
		// LgIndex_t indx[3];
		// double    TmpSpace[2];
		// double    DeltaXY[2];
		// Matrix_s  Jacobian;
		// Jacobian = CreateMatrix(2,2);

		Toler = 0.1 * (TetraTolerance - 1.0) * (XCMax - XCMin + YCMax - YCMin);
		TolerSqr = Toler * Toler;

		/* Set initial values of r,s,t - it out of range. */
		if (*r < -1.0 || *r > 1.0) *r = 0.0;
		if (*s < -1.0 || *s > 1.0) *s = 0.0;


		/* Compute error. */
		IsInside = QuadBilinearWeightDerivatives(*r, *s, W, dwdr, dwds);

		Xm = 0.0;
		Ym = 0.0;
		for (ic = 0; ic < 4; ic++)
		{
			Xm = Xm + W[ic] * XCell[ic];
			Ym = Ym + W[ic] * YCell[ic];
		}
		DeltaX = X - Xm;
		DeltaY = Y - Ym;
		DeltaSqr = DeltaX * DeltaX + DeltaY * DeltaY;


		/* Iterate */
		while ((IsOk && DeltaSqr > TolerSqr && iter < MAXNEWTONITER) ||
			   (FirstPass))
		{
			double dxdr = 0.0;
			double dxds = 0.0;
			double dydr = 0.0;
			double dyds = 0.0;
			double DeltaR, DeltaS, Det;

			iter++;

			/* Compute Jacobian Matrix */
			for (ic = 0; ic < 4; ic++)
			{
				dxdr = dxdr + dwdr[ic] * XCell[ic];
				dxds = dxds + dwds[ic] * XCell[ic];
				dydr = dydr + dwdr[ic] * YCell[ic];
				dyds = dyds + dwds[ic] * YCell[ic];
			}

			// Jacobian.val[0][0] = dxdr;
			// Jacobian.val[0][1] = dxds;
			// Jacobian.val[1][0] = dydr;
			// Jacobian.val[1][1] = dyds;

			/* Solve linear system (Jacobian) (DeltaRS) = (DeltaXY) for DeltaRS */
			// IsOk = cmludcmp(Jacobian, 2, indx, TmpSpace);

			// DeltaXY[0] = DeltaX;
			// DeltaXY[1] = DeltaY;
			// if (IsOk) IsOk = cmlubksb(Jacobian, 2, indx, DeltaXY);

			Det = dxdr * dyds - dxds * dydr;
			if (ABS(Det) > SMALLFLOAT)
			{
				double RDet = 1.0 / Det;
				DeltaR = RDet * (dyds * DeltaX - dxds * DeltaY);
				DeltaS = RDet * (-dydr * DeltaX + dxdr * DeltaY);
			}
			else
				IsOk = FALSE;

			// if (IsOk)
			//   {
			//     *r = *r + 0.5 * DeltaXY[0];
			//     *s = *s + 0.5 * DeltaXY[1];
			//   }
			if (IsOk)
			{
				*r = *r + 0.5 * DeltaR;
				*s = *s + 0.5 * DeltaS;
			}

			/* Recompute DeltaX,DeltaY for new r,s */
			IsInside = QuadBilinearWeightDerivatives(*r, *s, W, dwdr, dwds);
			Xm = 0.0;
			Ym = 0.0;
			for (ic = 0; ic < 4; ic++)
			{
				Xm = Xm + W[ic] * XCell[ic];
				Ym = Ym + W[ic] * YCell[ic];
			}
			DeltaX = X - Xm;
			DeltaY = Y - Ym;
			DeltaSqr = DeltaX * DeltaX + DeltaY * DeltaY;

			/* Been through once, set FirstPass = FALSE */
			FirstPass = FALSE;
		}

		/* Is it inside? */
		if (*r < TetraTolerance && *r > -TetraTolerance &&
			*s < TetraTolerance && *s > -TetraTolerance && IsOk)
		{
			CellFound = TRUE;
		}

		/* Clean up */
		// FreeMatrix(&Jacobian);
	}


	ENSURE(VALID_BOOLEAN(CellFound));
	return (CellFound);
}










/**
 * Compute the suface-tangent coordinates given deltas of X,Y,Z coordinates fron new origin.
 *
 * param PsiAlignment
 *     Coordinate direction (1==X, 2==Y, 3==Z) that most closely aligns with
 *     the psi surface-tangent coordinate direction. This, in combination with
 *     the normal vector, is sufficient to compute the transformation between
 *     surface-tangent and X,Y,Z coordinates.
 * param dX, dY, dZ
 *     Differnce of X, Y, Z coordinates form origin of psi,eta surface coords
 * param *dPsi, *dEta
 *     Pointers to surface tangent coordinates.
 *
 * return
 *     TRUE if operation is successful, otherwise FALSE.
 */
Boolean_t XYZtoPsiEta(int     PsiAlign,
					  XYZ_s   Normal,
					  double  dX,
					  double  dY,
					  double  dZ,
					  double *dPsi,
					  double *dEta)
{
	Boolean_t IsOk = TRUE;
	double t1x, t1y, t1z;
	double t2x, t2y, t2z;
	double rmagnitude;

	REQUIRE(PsiAlign > 0 && PsiAlign < 4);
	REQUIRE(VALID_REF(dPsi));
	REQUIRE(VALID_REF(dEta));

	// Compute the first tangent coordinate vector
	switch (PsiAlign)
	{
		case 1:  // normal is nearly aligned with z so t1 (psi) is nearly aligned with x
			// Cross product of n with i-coordinate vector gives t2 vector
			t2x =  0.0;
			t2y =  Normal.Z;
			t2z = -Normal.Y;
			rmagnitude = 1.0 / sqrt(t2y * t2y + t2z * t2z);
			break;
		case 2:  // Normal is nearly aligned with x, so t1 (psi) is nearly aligned with y
			// Cross product of n with j-coordinate vector gives t2 vector
			t2x = -Normal.Z;
			t2y =  0.0;
			t2z =  Normal.X;
			rmagnitude = 1.0 / sqrt(t2x * t2x + t2z * t2z);
			break;
		case 3:  // Normal is nearly aligned with y, so t1 (psi) is nearly aligned with z
			// Cross product of n with k-coordinate vector gives t2 vector
			t2x =  Normal.Y;
			t2y = -Normal.X;
			t2z =  0.0;
			rmagnitude = 1.0 / sqrt(t2x * t2x + t2y * t2y);
			break;

		default:
			IsOk = FALSE;
			break;
	}

	// Normalize t2
	t2x = t2x * rmagnitude;
	t2y = t2y * rmagnitude;
	t2z = t2z * rmagnitude;

	// Cross product of t2 with Normal gives t1
	t1x = t2y * Normal.Z - t2z * Normal.Y;
	t1y = t2z * Normal.X - t2x * Normal.Z;
	t1z = t2x * Normal.Y - t2y * Normal.X;

	// Compute the surface tanget coordinates
	*dPsi = dX * t1x + dY * t1y + dZ * t1z;
	*dEta = dX * t2x + dY * t2y + dZ * t2z;

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}









/**
 * Compute the deltas of X, Y, Z coordinates given the surface-tangent coordinates.
 *
 * param PsiAlignment
 *     Coordinate direction (1==X, 2==Y, 3==Z) that most closely aligns with
 *     the psi surface-tangent coordinate direction. This, in combination with
 *     the normal vector, is sufficient to compute the transformation between
 *     surface-tangent and X,Y,Z coordinates.
 * param dPsi, dEta
 *     Differnce of Psi, Eta surface-tangent coordinates
 * param *dX, *dY, *dZ
 *     Pointers to deltas of X, Y, Z coordinates.
 *
 * return
 *     TRUE if operation is successful, otherwise FALSE.
 */
Boolean_t PsiEtatoXYZ(int     PsiAlign,
					  XYZ_s   Normal,
					  double  dPsi,
					  double  dEta,
					  double *dX,
					  double *dY,
					  double *dZ)
{
	Boolean_t IsOk = TRUE;
	double t1x, t1y, t1z;
	double t2x, t2y, t2z;
	double rmagnitude;

	REQUIRE(PsiAlign > 0 && PsiAlign < 4);
	REQUIRE(VALID_REF(dX));
	REQUIRE(VALID_REF(dY));
	REQUIRE(VALID_REF(dZ));

	// Compute the first tangent coordinate vector
	switch (PsiAlign)
	{
		case 1:  // normal is nearly aligned with z so t1 (psi) is nearly aligned with x
			// Cross product of n with i-coordinate vector gives t2 vector
			t2x =  0.0;
			t2y =  Normal.Z;
			t2z = -Normal.Y;
			rmagnitude = 1.0 / sqrt(t2y * t2y + t2z * t2z);
			break;
		case 2:  // Normal is nearly aligned with x, so t1 (psi) is nearly aligned with y
			// Cross product of n with j-coordinate vector gives t2 vector
			t2x = -Normal.Z;
			t2y =  0.0;
			t2z =  Normal.X;
			rmagnitude = 1.0 / sqrt(t2x * t2x + t2z * t2z);
			break;
		case 3:  // Normal is nearly aligned with y, so t1 (psi) is nearly aligned with z
			// Cross product of n with k-coordinate vector gives t2 vector
			t2x =  Normal.Y;
			t2y = -Normal.X;
			t2z =  0.0;
			rmagnitude = 1.0 / sqrt(t2x * t2x + t2y * t2y);
			break;

		default:
			IsOk = FALSE;
			break;
	}

	// Normalize t2
	t2x = t2x * rmagnitude;
	t2y = t2y * rmagnitude;
	t2z = t2z * rmagnitude;

	// Cross product of t2 with Normal gives t1
	t1x = t2y * Normal.Z - t2z * Normal.Y;
	t1y = t2z * Normal.X - t2x * Normal.Z;
	t1z = t2x * Normal.Y - t2y * Normal.X;

	// Compute the surface tanget coordinates
	*dX = dPsi * t1x + dEta * t2x;
	*dY = dPsi * t1y + dEta * t2y;
	*dZ = dPsi * t1z + dEta * t2z;

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}







/**
 *
 * Compute PsiAlignment: the coordinate direction (1==X, 2==Y, 3==Z) that most closely
 * aligns with the psi surface-tangent coordinate direction. This, in combination with
 * the normal vector, is sufficient to compute the transformation between
 * surface-tangent and X,Y,Z coordinates.
 *
 * param IsQuad
 *     TRUE for FEQuad, FALSE for FETriangle
 * param XCell[4], YCell[4], ZCell[4]
 *     Arrays of coordinates of the cell corners.
 *
 * return (through pointers) XAve, YAve, ZAve
 *      Coordinates of the cell center.
 * return
 *     CellNormal: cell normal vector
 * return
 *     PsiAlign: 1==Psi->X, 2==Psi->Y, 3==Psi->Z).
 */
int CalcPsiAlign(Boolean_t    IsQuad,
				 double       XCell[4],
				 double       YCell[4],
				 double       ZCell[4],
				 double      *XAve,
				 double      *YAve,
				 double      *ZAve,
				 XYZ_pa       CellNormal)
{
	int PsiAlign = 0;
	double dxdpsi, dxdeta, dydpsi, dydeta, dzdpsi, dzdeta, RCellNormalTot;

	REQUIRE(VALID_BOOLEAN(IsQuad));
	REQUIRE(VALID_REF(XCell));
	REQUIRE(VALID_REF(YCell));
	REQUIRE(VALID_REF(ZCell));
	REQUIRE(VALID_REF(CellNormal));

	// Compute surface normal for cell
	if (IsQuad)
	{
		double x1 = XCell[0];
		double x2 = XCell[1];
		double x3 = XCell[2];
		double x4 = XCell[3];
		double y1 = YCell[0];
		double y2 = YCell[1];
		double y3 = YCell[2];
		double y4 = YCell[3];
		double z1 = ZCell[0];
		double z2 = ZCell[1];
		double z3 = ZCell[2];
		double z4 = ZCell[3];

		*XAve = 0.25 * (x1 + x2 + x3 + x4);
		*YAve = 0.25 * (y1 + y2 + y3 + y4);
		*ZAve = 0.25 * (z1 + z2 + z3 + z4);

		dxdpsi = 0.25 * (x2 - x1 + x3 - x2);
		dxdeta = 0.25 * (x3 - x2 + x4 - x1);
		dydpsi = 0.25 * (y2 - y1 + y3 - y2);
		dydeta = 0.25 * (y3 - y2 + y4 - y1);
		dzdpsi = 0.25 * (z2 - z1 + z3 - z2);
		dzdeta = 0.25 * (z3 - z2 + z4 - z1);
	}
	else   // collapsed to triangle
	{
		double OneThrd = 1.0 / 3.0;
		double x1 = XCell[0];
		double x2 = XCell[1];
		double x3 = XCell[2];
		double y1 = YCell[0];
		double y2 = YCell[1];
		double y3 = YCell[2];
		double z1 = ZCell[0];
		double z2 = ZCell[1];
		double z3 = ZCell[2];

		*XAve = OneThrd * (x1 + x2 + x3);
		*YAve = OneThrd * (y1 + y2 + y3);
		*ZAve = OneThrd * (z1 + z2 + z3);

		dxdpsi = 0.5 * (x2 - x1);
		dxdeta = 0.5 * (x3 - x1);
		dydpsi = 0.5 * (y2 - y1);
		dydeta = 0.5 * (y3 - y1);
		dzdpsi = 0.5 * (z2 - z1);
		dzdeta = 0.5 * (z3 - z1);
	}

	CellNormal->X = dydpsi * dzdeta - dydeta * dzdpsi;
	CellNormal->Y = dzdpsi * dxdeta - dzdeta * dxdpsi;
	CellNormal->Z = dxdpsi * dydeta - dxdeta * dydpsi;
	RCellNormalTot = 1.0 / sqrt(CellNormal->X * CellNormal->X + CellNormal->Y * CellNormal->Y +
								CellNormal->Z * CellNormal->Z);
	CellNormal->X *= RCellNormalTot;
	CellNormal->Y *= RCellNormalTot;
	CellNormal->Z *= RCellNormalTot;

	// If normal most closely aligns with Z, PsiAlignment is 1.
	if (ABS(CellNormal->Z) >= ABS(CellNormal->X) && ABS(CellNormal->Z) >= ABS(CellNormal->Y))
		PsiAlign = 1;
	// If normal most closely aligns with X, PsiAlignment is 2.
	else if (ABS(CellNormal->X) >= ABS(CellNormal->Y) && ABS(CellNormal->X) > ABS(CellNormal->Z))
		PsiAlign = 2;
	// If normal most closely aligns with Y, PsiAlignment is 3.
	else if (ABS(CellNormal->Y) > ABS(CellNormal->X) && ABS(CellNormal->Y) > ABS(CellNormal->Z))
		PsiAlign = 3;

	ENSURE(PsiAlign >= 0 && PsiAlign < 4);
	return PsiAlign;
}





