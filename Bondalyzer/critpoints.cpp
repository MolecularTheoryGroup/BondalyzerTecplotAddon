/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

//#ifdef _WIN32 
//	#define _COMPLEX_DEFINED 
//#endif

#include "TECADDON.h"

//#include "blaswrap.h"
//#include "f2c.h"
//#include "clapack.h"

#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "TASSERT.h"
#include "GUIDEFS.h"
#include "ENGINE.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "ARRLISTTOOLS.h"
#include "LINEARALGEBRA.h"
#include "ODERUNGEKUTTA.h"
#include "SURFELEMMAP.h"
#include "NORMALS.h"
#include "SURFACEFIT.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
#include "SVD.h"
// #include "SURFACEFIT.h"
#include "ELEMSHAPEFUNC.h"
#include "GRADPATH.h"
#include "GEOMTOOLS.h"
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>

#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"
#include "CSM_CALC_VARS.h"

#include <armadillo>
using namespace arma;


#define MAXNEWTONITER 100

static double TetraTolerance = 1.0005;
// static double RangeEpsilon = 0.0;




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





/**
 * Compute the linear index into the Tecplot arrays (1 based) from the I, J, K
 * indicies when there are periodic boundary conditions. In other words, if the
 * I, J, K, indicies exceed the limits of the zone, find the equivalent point inside
 * the zone.
 *
 * param I, J, K
 *     I, J, K indices for a 3-D ordered zone.
 * param IMax, JMax, KMax
 *     I, J, K index ranges for the 3D ordered zone.
 *
 * return
 *     Linear (1-based) index into the Tecplot field-data arrays.
 */
LgIndex_t PeriodicIndexFromIJK(LgIndex_t    I,
							   LgIndex_t    J,
							   LgIndex_t    K,
							   LgIndex_t IMax,
							   LgIndex_t JMax,
							   LgIndex_t KMax)
{
	LgIndex_t Result;
	LgIndex_t IIndex = I;
	LgIndex_t JIndex = J;
	LgIndex_t KIndex = K;

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1);
	REQUIRE(IIndex >= -1 && IIndex <= IMax + 3);
	REQUIRE(JIndex >= -1 && JIndex <= JMax + 3);
	REQUIRE(KIndex >= -1 && KIndex <= KMax + 3);

	if (IIndex < 1) IIndex = IIndex + IMax;
	if (JIndex < 1) JIndex = JIndex + JMax;
	if (KIndex < 1) KIndex = KIndex + KMax;
	if (IIndex > IMax) IIndex = IIndex - IMax;
	if (JIndex > JMax) JIndex = JIndex - JMax;
	if (KIndex > KMax) KIndex = KIndex - KMax;

	Result = IIndex + (JIndex - 1) * IMax + (KIndex - 1) * IMax * JMax;

	ENSURE(Result <= IMax * JMax * KMax + 1);
	return Result;
}








/**
 * Compute the linear index into the Tecplot arrays (1 based) from the I, J, K
 * indicies, both with and without periodic boundary conditions. In the periodic case, if the
 * I, J, K, indicies exceed the limits of the zone, find the equivalent point inside
 * the zone.
 *
 * param I, J, K
 *     I, J, K indices for a 3-D ordered zone.
 * param IMax, JMax, KMax
 *     I, J, K index ranges for the 3D ordered zone.
 * param PeriodicBC: True if periodic, FALSE otherwise
 *
 * return
 *     Linear (1-based) index into the Tecplot field-data arrays.
 */
LgIndex_t IndexFromIJK(LgIndex_t    I,
					   LgIndex_t    J,
					   LgIndex_t    K,
					   LgIndex_t IMax,
					   LgIndex_t JMax,
					   LgIndex_t KMax,
					   Boolean_t PeriodicBC)
{
	LgIndex_t Result;
	LgIndex_t IIndex = I;
	LgIndex_t JIndex = J;
	LgIndex_t KIndex = K;

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1);
	if (PeriodicBC)
	{
		REQUIRE(I >= -1 && I <= IMax + 3);
		REQUIRE(J >= -1 && J <= JMax + 3);
		REQUIRE(K >= -1 && K <= KMax + 3);
	}
	else
	{
		REQUIRE(I >= 1 && I <= IMax);
		REQUIRE(J >= 1 && J <= JMax);
		REQUIRE(K >= 1 && K <= KMax);
	}

	if (IIndex < 1) IIndex = IIndex + IMax;
	if (JIndex < 1) JIndex = JIndex + JMax;
	if (KIndex < 1) KIndex = KIndex + KMax;
	if (IIndex > IMax) IIndex = IIndex - IMax;
	if (JIndex > JMax) JIndex = JIndex - JMax;
	if (KIndex > KMax) KIndex = KIndex - KMax;

	Result = IIndex + (JIndex - 1) * IMax + (KIndex - 1) * IMax * JMax;

	ENSURE(Result <= IMax * JMax * KMax + 1);
	return Result;
}







Boolean_t ComputeCurvature(EntIndex_t   ZoneNum,
						   LgIndex_t    i,
						   LgIndex_t    j,
						   LgIndex_t    k,
						   FieldData_pa CDVarFDPtr,
						   Boolean_t    PeriodicBC,
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
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));

	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);


	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */
	if (PeriodicBC)
	{
		REQUIRE(i >= 0 && i <= IMax + 2);
		REQUIRE(j >= 0 && j <= JMax + 2);
		REQUIRE(k >= 0 && k <= KMax + 2);
	}
	else
	{
		REQUIRE(i > 0 && i <= IMax);
		REQUIRE(j > 0 && j <= JMax);
		REQUIRE(k > 0 && k <= KMax);
	}

	IJMax = IMax * JMax;

	/* d2dx2 */
	if (IsOk)
	{
		LgIndex_t IndxI, IndxIp1, IndxIm1;
		LgIndex_t ii  = i;

		if (PeriodicBC)
		{
			IndxI   = PeriodicIndexFromIJK(i    , j, k, IMax, JMax, KMax);
			IndxIp1 = PeriodicIndexFromIJK(i + 1, j, k, IMax, JMax, KMax);
			IndxIm1 = PeriodicIndexFromIJK(i - 1, j, k, IMax, JMax, KMax);
		}
		else
		{
			if (i == 1) ii  = 2;
			if (i == IMax) ii  = IMax - 1;

			IndxI = ii + (j - 1) * IMax + (k - 1) * IJMax;
			IndxIp1 = IndxI + 1;
			IndxIm1 = IndxI - 1;
		}

		*d2dx2 =       TecUtilDataValueGetByRef(CDVarFDPtr, IndxIp1)
					   - 2.0 * TecUtilDataValueGetByRef(CDVarFDPtr, IndxI)
					   + TecUtilDataValueGetByRef(CDVarFDPtr, IndxIm1);
	}

	/* d2dy2 */
	if (IsOk)
	{
		LgIndex_t IndxJ, IndxJp1, IndxJm1;
		LgIndex_t jj  = j;

		if (PeriodicBC)
		{
			IndxJ   = PeriodicIndexFromIJK(i, j    , k, IMax, JMax, KMax);
			IndxJp1 = PeriodicIndexFromIJK(i, j + 1, k, IMax, JMax, KMax);
			IndxJm1 = PeriodicIndexFromIJK(i, j - 1, k, IMax, JMax, KMax);
		}
		else
		{
			if (j == 1) jj  = 2;
			if (j == JMax) jj  = JMax - 1;

			IndxJ = i + (jj - 1) * IMax + (k - 1) * IJMax;
			IndxJp1 = IndxJ + IMax;
			IndxJm1 = IndxJ - IMax;
		}

		*d2dy2 =       TecUtilDataValueGetByRef(CDVarFDPtr, IndxJp1)
					   - 2.0 * TecUtilDataValueGetByRef(CDVarFDPtr, IndxJ)
					   + TecUtilDataValueGetByRef(CDVarFDPtr, IndxJm1);
	}

	/* d2dz2 */
	if (IsOk)
	{
		LgIndex_t IndxK, IndxKp1, IndxKm1;
		LgIndex_t kk  = k;

		if (PeriodicBC)
		{
			IndxK   = PeriodicIndexFromIJK(i, j, k    , IMax, JMax, KMax);
			IndxKp1 = PeriodicIndexFromIJK(i, j, k + 1, IMax, JMax, KMax);
			IndxKm1 = PeriodicIndexFromIJK(i, j, k - 1, IMax, JMax, KMax);
		}
		else
		{
			if (k == 1) kk  = 2;
			if (k == KMax) kk  = KMax - 1;

			IndxK = i + (j - 1) * IMax + (kk - 1) * IJMax;
			IndxKp1 = IndxK + IJMax;
			IndxKm1 = IndxK - IJMax;
		}

		*d2dz2 =       TecUtilDataValueGetByRef(CDVarFDPtr, IndxKp1)
					   - 2.0 * TecUtilDataValueGetByRef(CDVarFDPtr, IndxK)
					   + TecUtilDataValueGetByRef(CDVarFDPtr, IndxKm1);
	}

	/* d2dxdy, d2dxdz, and d2dydz */
	if (IsOk)
	{
		if (PeriodicBC)
		{
			LgIndex_t IndxIp1Jp1, IndxIm1Jp1, IndxIp1Jm1, IndxIm1Jm1;
			LgIndex_t IndxIp1Kp1, IndxIm1Kp1, IndxIp1Km1, IndxIm1Km1;
			LgIndex_t IndxJp1Kp1, IndxJm1Kp1, IndxJp1Km1, IndxJm1Km1;

			IndxIp1Jp1 = PeriodicIndexFromIJK(i + 1, j + 1, k, IMax, JMax, KMax);
			IndxIm1Jp1 = PeriodicIndexFromIJK(i - 1, j + 1, k, IMax, JMax, KMax);
			IndxIp1Jm1 = PeriodicIndexFromIJK(i + 1, j - 1, k, IMax, JMax, KMax);
			IndxIm1Jm1 = PeriodicIndexFromIJK(i - 1, j - 1, k, IMax, JMax, KMax);
			*d2dxdy = 0.25 * (TecUtilDataValueGetByRef(CDVarFDPtr, IndxIp1Jp1)
							  - TecUtilDataValueGetByRef(CDVarFDPtr, IndxIp1Jm1)
							  - TecUtilDataValueGetByRef(CDVarFDPtr, IndxIm1Jp1)
							  + TecUtilDataValueGetByRef(CDVarFDPtr, IndxIm1Jm1)) ;

			IndxIp1Kp1 = PeriodicIndexFromIJK(i + 1, j, k + 1, IMax, JMax, KMax);
			IndxIm1Kp1 = PeriodicIndexFromIJK(i - 1, j, k + 1, IMax, JMax, KMax);
			IndxIp1Km1 = PeriodicIndexFromIJK(i + 1, j, k - 1, IMax, JMax, KMax);
			IndxIm1Km1 = PeriodicIndexFromIJK(i - 1, j, k - 1, IMax, JMax, KMax);
			*d2dxdz = 0.25 * (TecUtilDataValueGetByRef(CDVarFDPtr, IndxIp1Kp1)
							  - TecUtilDataValueGetByRef(CDVarFDPtr, IndxIp1Km1)
							  - TecUtilDataValueGetByRef(CDVarFDPtr, IndxIm1Kp1)
							  + TecUtilDataValueGetByRef(CDVarFDPtr, IndxIm1Km1)) ;

			IndxJp1Kp1 = PeriodicIndexFromIJK(i, j + 1, k + 1, IMax, JMax, KMax);
			IndxJm1Kp1 = PeriodicIndexFromIJK(i, j - 1, k + 1, IMax, JMax, KMax);
			IndxJp1Km1 = PeriodicIndexFromIJK(i, j + 1, k - 1, IMax, JMax, KMax);
			IndxJm1Km1 = PeriodicIndexFromIJK(i, j - 1, k - 1, IMax, JMax, KMax);
			*d2dydz = 0.25 * (TecUtilDataValueGetByRef(CDVarFDPtr, IndxJp1Kp1)
							  - TecUtilDataValueGetByRef(CDVarFDPtr, IndxJp1Km1)
							  - TecUtilDataValueGetByRef(CDVarFDPtr, IndxJm1Kp1)
							  + TecUtilDataValueGetByRef(CDVarFDPtr, IndxJm1Km1)) ;
		}
		else
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
	}

	return(IsOk);
}






/**
 * For a given node of an FEQuad or FETriangle zone, compute the list of nodes
 * it connects to via edges. Also compute the alloted angle for each edge. The
 * alloted angle will be half of the corner angles for the two elements containing
 * this edge.
 *
 * param ZoneNum
 *     Number of target FE-surface zone number.
 * param NodeNum
 *     Number of the node of interest.
 * param SurfElemMap
 *     SurfElemMap for the target surface zone.
 * param EdgeOtherNodeList
 *     List of nodes connected to NodeNum via edges of the FE surface zone.
 *     Note: EdgeOtherNodeList is cleared of any existing contents first.
 * param EdgeAngleList
 *     ArrList containing the angles alloted to each edge
 *
 * return
 *     TRUE if operation is successful, otherwise FALSE.
 */
Boolean_t ComputeEdgeOtherNodeListFromSEM(EntIndex_t     ZoneNum,
										  LgIndex_t      NodeNum,
										  SurfElemMap_pa SurfElemMap,
										  ArrList_pa     EdgeOtherNodeList,
										  ArrList_pa     EdgeAngleList)
{
	Boolean_t  IsOk = TRUE;
	EntIndex_t NumZones, NumVars;
	LgIndex_t  NumNodes, NumElems, NumNodesPerElem;
	LgIndex_t  NumElemsForNode;
	NodeMap_pa   NodeMap = NULL;

	REQUIRE(SurfElemMapIsValid(SurfElemMap));
	REQUIRE(ArrListIsValid(EdgeOtherNodeList));
	REQUIRE(ArrListGetType(EdgeOtherNodeList) == ArrListType_Long);
	REQUIRE(ArrListIsValid(EdgeAngleList));
	REQUIRE(ArrListGetType(EdgeAngleList) == ArrListType_Double);
	REQUIRE(TecUtilDataSetIsAvailable());
	REQUIRE(NodeNum > 0);

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);

	TecUtilZoneGetInfo(ZoneNum, &NumNodes, &NumElems, &NumNodesPerElem, NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(NodeNum > 0 && NodeNum <= NumNodes && NumElems > 1);
	REQUIRE(TecUtilZoneIsFiniteElement(ZoneNum));

	NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
	CHECK(VALID_REF(NodeMap));

	// Clear EdgeOtherNodeList
	ArrListClear(EdgeOtherNodeList);

	NumElemsForNode = SurfElemMapGetElemCountForNode(SurfElemMap, NodeNum);

	// Find all edges radiating out from the node. Store them in EdgeOtherNodeList.
	if (IsOk)
	{
		LgIndex_t    ElemOffset;
		EntIndex_t   XVarNum = TecUtilVarGetNumByAssignment('X');
		EntIndex_t   YVarNum = TecUtilVarGetNumByAssignment('Y');
		EntIndex_t   ZVarNum = TecUtilVarGetNumByAssignment('Z');
		FieldData_pa XVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, XVarNum);
		FieldData_pa YVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, YVarNum);
		FieldData_pa ZVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, ZVarNum);

		for (ElemOffset = 0; ElemOffset < NumElemsForNode; ElemOffset++)
		{
			LgIndex_t Edge, CornerNodeBeg, CornerNodeEnd;
			LgIndex_t Elem = SurfElemMapGetElemNumForNodeOffset(SurfElemMap, NodeNum, ElemOffset);



			// Find the edges radiating out from the node. Store them by "other" node number
			// TEMP
			if (Elem < 1 || Elem > NumElems)
				TecUtilDialogMessageBox("Elem error in ComputeEdgeOtherNodeListFromSEM",
										MessageBox_Warning);

			CornerNodeEnd = TecUtilDataNodeGetByRef(NodeMap, Elem, 1);
			for (Edge = 1; Edge <= NumNodesPerElem; Edge++)
			{
				LgIndex_t NumEdgeOtherNodeList = ArrListGetCount(EdgeOtherNodeList);
				LgIndex_t CornerEnd = Edge + 1;
				LgIndex_t ItemOffset = -1;
				if (CornerEnd > NumNodesPerElem) CornerEnd -= NumNodesPerElem;

				CornerNodeBeg = CornerNodeEnd;
				CornerNodeEnd = TecUtilDataNodeGetByRef(NodeMap, Elem, CornerEnd);

				if (CornerNodeBeg == NodeNum && CornerNodeEnd != NodeNum)
					ItemOffset = ArrListAppendUniqueLongItem(EdgeOtherNodeList, CornerNodeEnd);

				if (CornerNodeEnd == NodeNum && CornerNodeBeg != NodeNum)
					ItemOffset = ArrListAppendUniqueLongItem(EdgeOtherNodeList, CornerNodeBeg);

				// Save the two element angles contribution to the edge. Use it for weights in SVD.
				if (ItemOffset >= 0)
				{
					ArrListItem_u Item;
					double PreviousAngle = 0.0;
					double XNode = TecUtilDataValueGetByRef(XVarFDPtr, NodeNum);
					double YNode = TecUtilDataValueGetByRef(YVarFDPtr, NodeNum);
					double ZNode = TecUtilDataValueGetByRef(ZVarFDPtr, NodeNum);
					double DXEdge1, DYEdge1, DZEdge1;
					double DXEdge2, DYEdge2, DZEdge2;
					double DotProd, MagEdge1, MagEdge2, Angle;

					if (CornerNodeBeg == NodeNum)
					{
						// Compute angle between edge and previous edge.
						LgIndex_t Corner2NodeEnd;
						LgIndex_t Corner2End = Edge - 1;
						if (Corner2End < 1) Corner2End += NumNodesPerElem;

						Corner2NodeEnd = TecUtilDataNodeGetByRef(NodeMap, Elem, Corner2End);

						DXEdge1 = TecUtilDataValueGetByRef(XVarFDPtr, CornerNodeEnd) - XNode;
						DYEdge1 = TecUtilDataValueGetByRef(YVarFDPtr, CornerNodeEnd) - YNode;
						DZEdge1 = TecUtilDataValueGetByRef(ZVarFDPtr, CornerNodeEnd) - ZNode;
						DXEdge2 = TecUtilDataValueGetByRef(XVarFDPtr, Corner2NodeEnd) - XNode;
						DYEdge2 = TecUtilDataValueGetByRef(YVarFDPtr, Corner2NodeEnd) - YNode;
						DZEdge2 = TecUtilDataValueGetByRef(ZVarFDPtr, Corner2NodeEnd) - ZNode;
					}
					else // CornerNodeEnd == NodeNum
					{
						// Compute angle between Edge and the following edge.
						LgIndex_t Corner2NodeEnd;
						LgIndex_t Corner2End = Edge - 1;
						if (Corner2End < 1) Corner2End += NumNodesPerElem;

						Corner2NodeEnd = TecUtilDataNodeGetByRef(NodeMap, Elem, Corner2End);

						DXEdge1 = TecUtilDataValueGetByRef(XVarFDPtr, CornerNodeBeg) - XNode;
						DYEdge1 = TecUtilDataValueGetByRef(YVarFDPtr, CornerNodeBeg) - YNode;
						DZEdge1 = TecUtilDataValueGetByRef(ZVarFDPtr, CornerNodeBeg) - ZNode;
						DXEdge2 = TecUtilDataValueGetByRef(XVarFDPtr, Corner2NodeEnd) - XNode;
						DYEdge2 = TecUtilDataValueGetByRef(YVarFDPtr, Corner2NodeEnd) - YNode;
						DZEdge2 = TecUtilDataValueGetByRef(ZVarFDPtr, Corner2NodeEnd) - ZNode;
					}
					DotProd = DXEdge1 * DXEdge2 + DYEdge1 * DYEdge2 + DZEdge1 * DZEdge2;
					MagEdge1 = sqrt(DXEdge1 * DXEdge1 + DYEdge1 * DYEdge1 + DZEdge1 * DZEdge1);
					MagEdge2 = sqrt(DXEdge2 * DXEdge2 + DYEdge2 * DYEdge2 + DZEdge2 * DZEdge2);
					Angle = 0.5 * acos(DotProd / (MagEdge1 * MagEdge2));

					// Sum to edge angle.
					if (ItemOffset < NumEdgeOtherNodeList)  // Second contribution
					{
						Item = ArrListGetItem(EdgeAngleList, ItemOffset);
						PreviousAngle = Item.Double;
					}
					Item.Double = Angle + PreviousAngle;
					IsOk = ArrListSetItem(EdgeAngleList, ItemOffset, Item);
				}
			}
		}
		CHECK((ArrListGetCount(EdgeOtherNodeList) == NumElemsForNode ||
			   ArrListGetCount(EdgeOtherNodeList) == NumElemsForNode + 1) &&
			  ArrListGetCount(EdgeAngleList) == ArrListGetCount(EdgeOtherNodeList));
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}








/**
 * Compute the curvature in of the scalar, Rho, for a specific (r,s) point in a
 * specific (ElemNum) element of an FEQuad surface.
 *
 * param ZoneNum
 *     Number of desired FEQuad surface zone.
 * param ElemNum
 *     Number of FEQuad element in ZoneNum. Repeated nodes may degenerate this
 *     element into a triangle.
 * param PsiAlign
 *     1 if Xp most nearly aligned with X, 2 if most nearly aligned with Y,
 *     3 if most nearly aligned with Z
 * param Normal
 *     XYZ_s normal vector to the surface for the element.
 * param UPCell, VPCell
 *     Xp and Yp Gradients of Rho at the nodes of the element.
 * param r, s
 *     Natural coordinates of the point in the cell.
 *
 *     4 ------------------------- 3
 *       |                       |
 *       |                       |
 *       |           ^s          |
 *       |           |           |             -1 <= r <= 1
 *       |           |           |             -1 <= s <= 1
 *       |           ---->r      |
 *       |                       |
 *       |                       |
 *       |                       |
 *       |                       |
 *     1 ------------------------- 2
 *
 * return (in pointers) D2RhoDx2, D2RhoDxDy, D2RhoDy2
 *     Curvature in the Xp, Yp surface-normal coordinate system
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t ComputeFEQuadCurvature(// EntIndex_t   ZoneNum,
	// LgIndex_t    ElemNum,
	// int          PsiAlign,
	// XYZ_s        Normal,
	double       XPCell[4],
	double       YPCell[4],
	double       UPCell[4],
	double       VPCell[4],
	double       r,
	double       s,
	double      *D2RhoDx2,
	double      *D2RhoDxDy,
	double      *D2RhoDy2)
{
	Boolean_t  IsOk = TRUE;
	// LgIndex_t  NumNodes, NumElems, NumNodesPerElem;
	// EntIndex_t NumZones, NumVars;

	double     OneMnsR, OnePlsR, OneMnsS, OnePlsS;
	double     DRDXp,  DSDXp,  DRDYp,  DSDYp;
	double     DUpDXp, DUpDYp, DVpDXp, DVpDYp;

	// REQUIRE(TecUtilDataSetIsAvailable());
	// REQUIRE(PsiAlign >= 1 && PsiAlign <= 3);
	// REQUIRE(Normal.X * Normal.X + Normal.Y * Normal.Y + Normal.Z * Normal.Z > 0.9);
	REQUIRE(VALID_REF(UPCell));
	REQUIRE(VALID_REF(VPCell));
	REQUIRE(r > -TetraTolerance && r < TetraTolerance);
	REQUIRE(s > -TetraTolerance && r < TetraTolerance);
	REQUIRE(VALID_REF(D2RhoDx2));
	REQUIRE(VALID_REF(D2RhoDxDy));
	REQUIRE(VALID_REF(D2RhoDy2));



	/* Get num zones in dataset */
	// if (IsOk)
	//     IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	// REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	// REQUIRE(TecUtilZoneGetType(ZoneNum) == ZoneType_FEQuad);

	// TecUtilZoneGetInfo(ZoneNum, &NumNodes, &NumElems, &NumNodesPerElem, NULL, NULL, NULL,
	//                    NULL, NULL, NULL, NULL, NULL, NULL, NULL);


	// REQUIRE(NumNodes > 1 && NumElems > 1 && NumNodesPerElem == 4); /* Must be FEQuad surface */
	// REQUIRE(ElemNum > 0 && ElemNum <= NumElems);

	// Basics for derivatives
	OneMnsR = 1 - r;
	OnePlsR = 1 + r;
	OneMnsS = 1 - s;
	OnePlsS = 1 + s;

	// Calculate coordinate transformation from r,s derivatives to xp,yp derivatives
	if (IsOk)
	{
		double DXpDR, DXpDS, DYpDR, DYpDS;
		double Denom, RDenom;

		// Reverse coord tranformation from chain rule
		DXpDR = 0.25 * (OneMnsS * (XPCell[1] - XPCell[0]) + OnePlsS * (XPCell[2] - XPCell[3]));
		DXpDS = 0.25 * (OneMnsR * (XPCell[3] - XPCell[0]) + OnePlsR * (XPCell[2] - XPCell[1]));
		DYpDR = 0.25 * (OneMnsS * (YPCell[1] - YPCell[0]) + OnePlsS * (YPCell[2] - YPCell[3]));
		DYpDS = 0.25 * (OneMnsR * (YPCell[3] - YPCell[0]) + OnePlsR * (YPCell[2] - YPCell[1]));

		// Invert reverse-coord-trans matrix
		Denom = DXpDR * DYpDS - DXpDS * DYpDR;
		if (ABS(Denom) > 1.0e-12)
		{
			RDenom = 1.0 / Denom;
			DRDXp =   RDenom * DYpDS;
			DSDXp = - RDenom * DYpDR;
			DRDYp = - RDenom * DXpDS;
			DSDYp =   RDenom * DXpDR;
		}
		else
			IsOk = FALSE;
	}

	CHECK(IsOk);

	// Compute derivatives of Up and Vp
	if (IsOk)
	{
		double DUpDR, DUpDS, DVpDR, DVpDS;
		DUpDR = 0.25 * (OneMnsS * (UPCell[1] - UPCell[0]) + OnePlsS * (UPCell[2] - UPCell[3]));
		DUpDS = 0.25 * (OneMnsR * (UPCell[3] - UPCell[0]) + OnePlsR * (UPCell[2] - UPCell[1]));
		DVpDR = 0.25 * (OneMnsS * (VPCell[1] - VPCell[0]) + OnePlsS * (VPCell[2] - VPCell[3]));
		DVpDS = 0.25 * (OneMnsR * (VPCell[3] - VPCell[0]) + OnePlsR * (VPCell[2] - VPCell[1]));

		DUpDXp = DUpDR * DRDXp + DUpDS * DSDXp;
		DUpDYp = DUpDR * DRDYp + DUpDS * DSDYp;
		DVpDXp = DVpDR * DRDXp + DVpDS * DSDXp;
		DVpDYp = DVpDR * DRDYp + DVpDS * DSDYp;
	}

	// Compute curvatures from derivatives of Up and Vp
	if (IsOk)
	{
		*D2RhoDx2  = DUpDXp;
		*D2RhoDxDy = 0.5 * (DUpDYp + DVpDXp);
		*D2RhoDy2  = DVpDYp;
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return(IsOk);
}











/**
 * Compute the curvature in of the scalar, Rho, for a specific (r,s) point in a
 * specific (ElemNum) element of an FEQuad surface.
 *
 * param XPCell, YPCell
 *     Xp and Yp location of nodes in a surface-tangent coordinate system.
 * param UPCell, VPCell
 *     Xp and Yp Gradients of Rho at the nodes of the element.
 *
 *
 *
 *                              /| 3
 *                            /  |
 *                          /    |
 *                        /      |
 *                      /        |
 *                    /          |
 *                  /            |
 *                /              |
 *              /                |
 *            /                  |
 *          /                    |
 *        /                      |
 *     1 ------------------------- 2
 *
 * return (in pointers) D2RhoDx2, D2RhoDxDy, D2RhoDy2
 *     Curvature in the Xp, Yp surface-tangent coordinate system
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t ComputeFETriangleCurvature(double       XPCell[3],
									 double       YPCell[3],
									 double       UPCell[3],
									 double       VPCell[3],
									 double      *D2RhoDx2,
									 double      *D2RhoDxDy,
									 double      *D2RhoDy2)
{
	Boolean_t  IsOk = TRUE;

	double     DUpDXp, DUpDYp, DVpDXp, DVpDYp;


	REQUIRE(VALID_REF(XPCell));
	REQUIRE(VALID_REF(YPCell));
	REQUIRE(VALID_REF(UPCell));
	REQUIRE(VALID_REF(VPCell));
	REQUIRE(VALID_REF(D2RhoDx2));
	REQUIRE(VALID_REF(D2RhoDxDy));
	REQUIRE(VALID_REF(D2RhoDy2));


	// Compute derivatives of Up and Vp
	if (IsOk)
	{
		double b1 = YPCell[1] - YPCell[2];
		double c1 = XPCell[2] - XPCell[1];

		double b2 = YPCell[2] - YPCell[0];
		double c2 = XPCell[0] - XPCell[2];

		double b3 = YPCell[0] - YPCell[1];
		double c3 = XPCell[1] - XPCell[0];

		double TwoDelta = XPCell[1] * YPCell[2] - XPCell[2] * YPCell[1]
						  + XPCell[0] * (YPCell[1] - YPCell[2])
						  + YPCell[0] * (XPCell[2] - XPCell[1]);

		double RTwoDelta = 1.0 / TwoDelta;

		DUpDXp = RTwoDelta * (b1 * UPCell[0] + b2 * UPCell[1] + b3 * UPCell[2]);
		DUpDYp = RTwoDelta * (c1 * UPCell[0] + c2 * UPCell[1] + c3 * UPCell[2]);

		DVpDXp = RTwoDelta * (b1 * VPCell[0] + b2 * VPCell[1] + b3 * VPCell[2]);
		DVpDYp = RTwoDelta * (c1 * VPCell[0] + c2 * VPCell[1] + c3 * VPCell[2]);
	}

	// Compute curvatures from derivatives of Up and Vp
	if (IsOk)
	{
		*D2RhoDx2  = DUpDXp;
		*D2RhoDxDy = 0.5 * (DUpDYp + DVpDXp);
		*D2RhoDy2  = DVpDYp;
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return(IsOk);
}










/**
 * Compute the critical point (position and type) using a quadratic surface fit of the
 * scalar, Rho, for a specific (ElemNum) element of an FEQuad surface.
 *
 * param ZoneNum
 *     Number of desired FEQuad surface zone.
 * param UVarNum, VVarNum, WVarNum
 *     Variable numbers fo the scalar-variable gradients. Not really used in this function
 * param ChrgDensVarNum
 *     Variable number of the scalar variable whose critical points we are looking for.
 * param TypeVarNum
 *     Variable number fo the critical point type (-2 = max, 0 = saddle, 2 = min)
 * param Normals
 *     Pointer to Normals_pa structure containing node-normals
 * param SurfElemMap
 *     Pointer to SurfElemMap_pa structure containing list of elements surrounding each node
 * param PsiAlignment
 *     Pointer to ArrList containing the PsiAlign for each node
 *     1 if Xp most nearly aligned with X, 2 if most nearly aligned with Y,
 *     3 if most nearly aligned with Z
 * param ElemNum
 *     Element Number
 * return (in pointers) XCrtPt, YCrtPt, ZCrtPt, ChrgDensCrtPt, TypeCrtPt
 *     Location, value of scalar, and type at the critical point.
 * return (in pointers) PrincDirX, PrincDirY, PrincDirZ
 *     Components of principle direction vector.
 *
 * return
 *     TRUE if critical point is in cell, FALSE otherwise.
 */
Boolean_t CriticalPointInQuadCellSF(EntIndex_t  ZoneNum,
									EntIndex_t  UVarNum,
									EntIndex_t  VVarNum,
									EntIndex_t  WVarNum,
									EntIndex_t  ChrgDensVarNum,
									EntIndex_t  TypeVarNum,
									Normals_pa  Normals,
									SurfElemMap_pa SurfElemMap,
									ArrList_pa  PsiAlignment,
									LgIndex_t   ElemNum,
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
	Boolean_t IsQuad = TRUE;
	int       NumNodesInElem = 4;
	LgIndex_t IMax, JMax, KMax;
	LgIndex_t Index[4];
	int       PsiAlign;
	XYZ_s     CellNormal;
	double XAve, YAve, ZAve;
	double RhoMinZone, RhoMaxZone;
	double RhoMinRegion =  LARGEFLOAT;
	double RhoMaxRegion = -LARGEFLOAT;

	EntIndex_t NumZones, NumVars;
	EntIndex_t XVarNum, YVarNum, ZVarNum;
	FieldData_pa XVarFDPtr = NULL;
	FieldData_pa YVarFDPtr = NULL;
	FieldData_pa ZVarFDPtr = NULL;
	FieldData_pa UVarFDPtr = NULL;
	FieldData_pa VVarFDPtr = NULL;
	FieldData_pa WVarFDPtr = NULL;
	FieldData_pa CDVarFDPtr = NULL;
	FieldData_pa TypeVarFDPtr = NULL;
	SurfaceFit_pa SurfaceFit = NULL;

	ArrList_pa NodesInSurfaceFit = NULL;

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


	// Get num zones in dataset
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	REQUIRE(TecUtilZoneGetType(ZoneNum) == ZoneType_FEQuad);
	REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
	REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
	REQUIRE(WVarNum > 0 && WVarNum <= NumVars);
	REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
	REQUIRE(TypeVarNum > 0 && TypeVarNum <= NumVars);

	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);


	REQUIRE(IMax > 1 && JMax > 1); // Must be 2D surface
	REQUIRE(ElemNum > 0 && ElemNum <= JMax);

	XVarNum   = TecUtilVarGetNumByAssignment('X');
	YVarNum   = TecUtilVarGetNumByAssignment('Y');
	ZVarNum   = TecUtilVarGetNumByAssignment('Z');

	XVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, XVarNum);
	YVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, YVarNum);
	ZVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, ZVarNum);
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

	// Find the min/max of the scalar varaible for the zone
	TecUtilDataValueGetMinMaxByRef(CDVarFDPtr, &RhoMinZone, &RhoMaxZone);


	/*
	 * FE Quad and ZoneType_Ordered Data:
	 *                                                            *
	 * 4         3                                                *
	 * +---------+                                                *
	 * |         |                                                *
	 * |         |                                                *
	 * |         |                                                *
	 * +---------+                                                *
	 * 0         1                                                *
	 */

	/* Set Velocities at corner of cell & compute vel min/max */
	if (IsOk)
	{
		ZoneType_e ZoneType = TecUtilZoneGetType(ZoneNum);

		if (ZoneType == ZoneType_FEQuad)
		{
			Index[0] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, 1);
			Index[1] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, 2);
			Index[2] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, 3);
			Index[3] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, 4);
		}
		/*
		else if(ZoneType == ZoneType_Ordered)
		  {
			Index[0] = IIndex + (JIndex - 1) * IMax;
			Index[1] = IIndex + 1 + (JIndex - 1) * IMax;
			Index[2] = IIndex + 1 + JIndex * IMax;
			Index[3] = IIndex + JIndex * IMax;
		  } */
		else
			IsOk = FALSE;
	}

	// Look for degenerate elements (triangles),
	// Adjust Index[] so first three items are Nodes of triangle
	NumNodesInElem = 4;
	if (Index[3] == Index[0] || Index[2] == Index[3])
	{
		IsQuad = FALSE;
		NumNodesInElem--;
	}
	if (Index[2] == Index[1])
	{
		IsQuad = FALSE;
		Index[2] = Index[3];
		NumNodesInElem--;
	}
	if (Index[1] == Index[0])
	{
		IsQuad = FALSE;
		Index[1] = Index[2];
		Index[2] = Index[3];
		NumNodesInElem--;
	}
	// Can only handle Quads or Triangles
	if (NumNodesInElem < 3) IsOk = FALSE;


	// Determine PsiAlignment for cell.
	if (IsOk)
	{
		double x1 = TecUtilDataValueGetByRef(XVarFDPtr, Index[0]);
		double x2 = TecUtilDataValueGetByRef(XVarFDPtr, Index[1]);
		double x3 = TecUtilDataValueGetByRef(XVarFDPtr, Index[2]);
		double x4 = TecUtilDataValueGetByRef(XVarFDPtr, Index[3]);
		double y1 = TecUtilDataValueGetByRef(YVarFDPtr, Index[0]);
		double y2 = TecUtilDataValueGetByRef(YVarFDPtr, Index[1]);
		double y3 = TecUtilDataValueGetByRef(YVarFDPtr, Index[2]);
		double y4 = TecUtilDataValueGetByRef(YVarFDPtr, Index[3]);
		double z1 = TecUtilDataValueGetByRef(ZVarFDPtr, Index[0]);
		double z2 = TecUtilDataValueGetByRef(ZVarFDPtr, Index[1]);
		double z3 = TecUtilDataValueGetByRef(ZVarFDPtr, Index[2]);
		double z4 = TecUtilDataValueGetByRef(ZVarFDPtr, Index[3]);

		double dxdpsi, dxdeta, dydpsi, dydeta, dzdpsi, dzdeta, RCellNormalTot;


		// Compute surface normal for cell
		if (IsQuad)
		{
			XAve = 0.25 * (x1 + x2 + x3 + x4);
			YAve = 0.25 * (y1 + y2 + y3 + y4);
			ZAve = 0.25 * (z1 + z2 + z3 + z4);

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

			XAve = OneThrd * (x1 + x2 + x3);
			YAve = OneThrd * (y1 + y2 + y3);
			ZAve = OneThrd * (z1 + z2 + z3);

			dxdpsi = 0.5 * (x2 - x1);
			dxdeta = 0.5 * (x3 - x1);
			dydpsi = 0.5 * (y2 - y1);
			dydeta = 0.5 * (y3 - y1);
			dzdpsi = 0.5 * (z2 - z1);
			dzdeta = 0.5 * (z3 - z1);
		}

		CellNormal.X = dydpsi * dzdeta - dydeta * dzdpsi;
		CellNormal.Y = dzdpsi * dxdeta - dzdeta * dxdpsi;
		CellNormal.Z = dxdpsi * dydeta - dxdeta * dydpsi;
		RCellNormalTot = 1.0 / sqrt(CellNormal.X * CellNormal.X + CellNormal.Y * CellNormal.Y +
									CellNormal.Z * CellNormal.Z);
		CellNormal.X *= RCellNormalTot;
		CellNormal.Y *= RCellNormalTot;
		CellNormal.Z *= RCellNormalTot;

		// If normal most closely aligns with Z, PsiAlignment is 1.
		if (ABS(CellNormal.Z) >= ABS(CellNormal.X) && ABS(CellNormal.Z) >= ABS(CellNormal.Y))
			PsiAlign = 1;
		// If normal most closely aligns with X, PsiAlignment is 2.
		else if (ABS(CellNormal.X) >= ABS(CellNormal.Y) && ABS(CellNormal.X) > ABS(CellNormal.Z))
			PsiAlign = 2;
		// If normal most closely aligns with Y, PsiAlignment is 3.
		else if (ABS(CellNormal.Y) > ABS(CellNormal.X) && ABS(CellNormal.Y) > ABS(CellNormal.Z))
			PsiAlign = 3;

		// Needed later in curvature calculation
		// DXCell[0] = x1 - XAve;
		// DXCell[1] = x2 - XAve;
		// DXCell[2] = x3 - XAve;
		// DXCell[3] = x4 - XAve;
		// DYCell[0] = y1 - YAve;
		// DYCell[1] = y2 - YAve;
		// DYCell[2] = y3 - YAve;
		// DYCell[3] = y4 - YAve;
		// DZCell[0] = z1 - ZAve;
		// DZCell[1] = z2 - ZAve;
		// DZCell[2] = z3 - ZAve;
		// DZCell[3] = z4 - ZAve;
	}

	// Allocate the SurfaceFit structure
	SurfaceFit = SurfaceFitAlloc();
	if (SurfaceFit == NULL) IsOk = FALSE;

	// Allocate the ArrList for the unique nodes in the surface fit
	NodesInSurfaceFit = ArrListAlloc(20, ArrListType_Long);
	if (NodesInSurfaceFit == NULL) IsOk = FALSE;

	// Find the elements surrounding the current element, unwrap (flatten) the local
	// surface segment, and add the nodes as points to the surface fit
	if (IsOk)
	{
		// First add nodes in element
		for (ic = 0; IsOk && ic < NumNodesInElem; ic++)
		{
			LgIndex_t    Node = Index[ic];
			double       DelX, DelY, DelZ, Psi, Eta, Rho;

			// First node are not already in the list, add it to the list
			// and the SurfaceFit
			ArrListAppendUniqueLongItem(NodesInSurfaceFit, Node);

			DelX = TecUtilDataValueGetByRef(XVarFDPtr, Node) - XAve;
			DelY = TecUtilDataValueGetByRef(YVarFDPtr, Node) - YAve;
			DelZ = TecUtilDataValueGetByRef(ZVarFDPtr, Node) - ZAve;

			Rho = TecUtilDataValueGetByRef(CDVarFDPtr, Node);

			// Find min/max of Rho for region used in surface fit
			RhoMinRegion = MIN(Rho, RhoMinRegion);
			RhoMaxRegion = MAX(Rho, RhoMaxRegion);

			// Using PsiAlign and Normal for the cell, compute surface-coord position of node
			IsOk = XYZtoPsiEta(PsiAlign, CellNormal, DelX, DelY, DelZ, &Psi, &Eta);

			if (IsOk) IsOk = SurfaceFitAppendPointAtEnd(SurfaceFit, Psi, Eta, Rho);
		}

		// Now add nodes in surrounding elements
		for (ic = 0; IsOk && ic < NumNodesInElem; ic++)
		{
			int ie;
			LgIndex_t    Node = Index[ic];
			LgIndex_t    NumElemForNode = SurfElemMapGetElemCountForNode(SurfElemMap, Node);
			for (ie = 0; ie < NumElemForNode; ie++)
			{
				int Corner;
				LgIndex_t ElemNum = SurfElemMapGetElemNumForNodeOffset(SurfElemMap, Node, ie);
				for (Corner = 1; Corner <= 4; Corner++)
				{
					LgIndex_t    ElemNode = (LgIndex_t)TecUtilDataNodeGetByZone(ZoneNum, ElemNum, Corner);
					double       DelX, DelY, DelZ, Psi, Eta, Rho;

					if (ArrListIsUniqueLongItem(NodesInSurfaceFit, ElemNode))
					{
						DelX = TecUtilDataValueGetByRef(XVarFDPtr, ElemNode) - XAve;
						DelY = TecUtilDataValueGetByRef(YVarFDPtr, ElemNode) - YAve;
						DelZ = TecUtilDataValueGetByRef(ZVarFDPtr, ElemNode) - ZAve;

						Rho = TecUtilDataValueGetByRef(CDVarFDPtr, ElemNode);

						// Find min/max of Rho for region used in surface fit
						RhoMinRegion = MIN(Rho, RhoMinRegion);
						RhoMaxRegion = MAX(Rho, RhoMaxRegion);

						ArrListAppendUniqueLongItem(NodesInSurfaceFit, ElemNode);

						// Using PsiAlign and Normal for the cell, compute surface-coord position of node
						IsOk = XYZtoPsiEta(PsiAlign, CellNormal, DelX, DelY, DelZ, &Psi, &Eta);

						if (IsOk) IsOk = SurfaceFitAppendPointAtEnd(SurfaceFit, Psi, Eta, Rho);
					}
				}
			}
		}
	}

	if (IsOk)
		IsOk = SurfaceFitCompute(SurfaceFit, SurfFitType_Quadratic);

	if (IsOk)
		IsOk = SurfaceFitCompCritPoint(SurfaceFit);

	if (IsOk)
	{
		Boolean_t IsInside = TRUE;
		double PsiCrtPt = SurfaceFit->c1_crt;
		double EtaCrtPt = SurfaceFit->c2_crt;

		double Psi1, Eta1, Psi2, Eta2, Psi3, Eta3, Psi4, Eta4;
		double Rho;

		IsOk = SurfaceFitGetPointFromOffset(SurfaceFit, 0, &Psi1, &Eta1, &Rho);
		IsOk = SurfaceFitGetPointFromOffset(SurfaceFit, 1, &Psi2, &Eta2, &Rho);

		// TEMP: Expand size of cell slightly
		Psi1 = 1.05 * Psi1;
		Eta1 = 1.05 * Eta1;
		Psi2 = 1.05 * Psi2;
		Eta2 = 1.05 * Eta2;

		// Is this inside the cell?
		// For each edge, Edge X Node->CrtPt should be positive (in Psi, Eta coord)
		if (IsInside)
		{
			if (((Psi2 - Psi1) *(EtaCrtPt - Eta1) - (Eta2 - Eta1) *(PsiCrtPt - Psi1)) <= 0.0) IsInside = FALSE;
		}
		if (IsInside)
		{
			IsOk = SurfaceFitGetPointFromOffset(SurfaceFit, 2, &Psi3, &Eta3, &Rho);

			// TEMP: Expand size of cell slightly
			Psi3 = 1.05 * Psi3;
			Eta3 = 1.05 * Eta3;

			if (((Psi3 - Psi2) *(EtaCrtPt - Eta2) - (Eta3 - Eta2) *(PsiCrtPt - Psi2)) <= 0.0) IsInside = FALSE;
		}
		if (IsInside && NumNodesInElem > 3)
		{
			IsOk = SurfaceFitGetPointFromOffset(SurfaceFit, 3, &Psi4, &Eta4, &Rho);

			// TEMP: Expand size of cell slightly
			Psi4 = 1.05 * Psi4;
			Eta4 = 1.05 * Eta4;

			if (((Psi4 - Psi3) *(EtaCrtPt - Eta3) - (Eta4 - Eta3) *(PsiCrtPt - Psi3)) <= 0.0) IsInside = FALSE;
			if (((Psi1 - Psi4) *(EtaCrtPt - Eta4) - (Eta1 - Eta4) *(PsiCrtPt - Psi4)) <= 0.0) IsInside = FALSE;
		}
		if (IsInside && NumNodesInElem == 3)
		{
			if (((Psi1 - Psi3) *(EtaCrtPt - Eta3) - (Eta1 - Eta3) *(PsiCrtPt - Psi3)) <= 0.0) IsInside = FALSE;
		}

		if (IsInside)
		{
			double dX, dY, dZ;

			CellFound = TRUE;

			IsOk = PsiEtatoXYZ(PsiAlign, CellNormal, PsiCrtPt, EtaCrtPt, &dX, &dY, &dZ);
			*XCrtPt = XAve + dX;
			*YCrtPt = YAve + dY;
			*ZCrtPt = ZAve + dZ;

			*TypeCrtPt = -2;
			if (SurfaceFit->lambda1 > 0.0 || SurfaceFit->lambda2 > 0.0)
			{
				if (SurfaceFit->lambda1 > 0.0 && SurfaceFit->lambda2 > 0.0)
					*TypeCrtPt = 2;
				else
				{
					double Angle, SinAngle, CosAngle;
					double XCompEV, YCompEV, ZCompEV, RTotEV;

					*TypeCrtPt = 0;

					Angle = SurfaceFit->RotationAngle;
					SinAngle = sin(Angle);
					CosAngle = cos(Angle);

					IsOk = PsiEtatoXYZ(PsiAlign, CellNormal, CosAngle, SinAngle, &XCompEV, &YCompEV, &ZCompEV);

					RTotEV = 1.0 / sqrt(XCompEV * XCompEV + YCompEV * YCompEV + ZCompEV * ZCompEV);
					*PrincDirX = (float)(XCompEV * RTotEV);
					*PrincDirY = (float)(YCompEV * RTotEV);
					*PrincDirZ = (float)(ZCompEV * RTotEV);
				}
			}
		}
	}

	// Max if scalar value at current node is larger than the values in all nodes of
	// ____ surrounding layers of cells.
	if (IsOk && !CellFound)
	{
		// TODO
		// for (nn=1; nn<=IMax; nn++)


	}

	// Min if scalar value at current node is smaller than the values at all nodes of
	// ____ surrounding layers of cells.
	if (IsOk && !CellFound)
	{
		// TODO:
	}




	return(CellFound);
}








Boolean_t CriticalPointInQuadCell(EntIndex_t  ZoneNum,
								  EntIndex_t  UVarNum,
								  EntIndex_t  VVarNum,
								  EntIndex_t  WVarNum,
								  EntIndex_t  GradMagVarNum,
								  EntIndex_t  TypeVarNum,
								  Normals_pa  Normals,
								  SurfElemMap_pa SurfElemMap,
								  ArrList_pa  PsiAlignment,
								  LgIndex_t   ElemNum,
								  double     *XCrtPt,
								  double     *YCrtPt,
								  double     *ZCrtPt,
								  double     *GradMagCrtPt,
								  EntIndex_t *TypeCrtPt,
								  float      *PrincDirX,
								  float      *PrincDirY,
								  float      *PrincDirZ)
{
	Boolean_t IsOk = TRUE;
	Boolean_t CellFound = FALSE;
	Boolean_t IsQuad = TRUE;
	int       NumNodesInElem = 4;
	LgIndex_t IMax, JMax, KMax;
	double    DXCell[4], DYCell[4], DZCell[4]; // Diff between node and cell center coordinates
	double    XPCell[4], YPCell[4];         // Surface-tangent coordinates of cell corners
	double    UCell[4], VCell[4], WCell[4]; // Velocities at cell corners
	double    UPCell[4], VPCell[4];         // Surface-tangent velocities at cell corners
	LgIndex_t Index[4];
	double    UPCMin, UPCMax, VPCMin, VPCMax;  // Min and Max vel for cell
	int       PsiAlign;
	XYZ_s     CellNormal;
	//	Next two are for data agitation to eliminate spurious CPs
	const double	AgitationFactor = 0.000001, SaddleFactor = 1.0;
	const int		AgitationNumber = 200;

	double    W[4];     // Trilinear weight for critical point
	double    r, s;  // Cell natural coordinates of critical point
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
	REQUIRE(VALID_REF(GradMagCrtPt));
	REQUIRE(VALID_REF(TypeCrtPt));
	REQUIRE(VALID_REF(PrincDirX));
	REQUIRE(VALID_REF(PrincDirY));
	REQUIRE(VALID_REF(PrincDirZ));

	// Get num zones in dataset
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	REQUIRE(TecUtilZoneGetType(ZoneNum) == ZoneType_FEQuad);
	REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
	REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
	REQUIRE(WVarNum > 0 && WVarNum <= NumVars);
	REQUIRE(GradMagVarNum > 0 && GradMagVarNum <= NumVars);
	REQUIRE(TypeVarNum > 0 && TypeVarNum <= NumVars);

	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr, // TODO: bad way to get FD pointers
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);


	REQUIRE(IMax > 1 && JMax > 1); // Must be 2D surface
	REQUIRE(ElemNum > 0 && ElemNum <= JMax);
	REQUIRE(VALID_REF(XVarFDPtr));
	REQUIRE(VALID_REF(YVarFDPtr));
	REQUIRE(VALID_REF(ZVarFDPtr));

	//	Initialize random time seed for data agitation
	srand((unsigned int)time(NULL));

	UVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, UVarNum);
	VVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, VVarNum);
	WVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, WVarNum);
	CDVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, GradMagVarNum);
	TypeVarFDPtr = TecUtilDataValueGetWritableRef(ZoneNum, TypeVarNum);
	if (UVarFDPtr == NULL || VVarFDPtr == NULL || WVarFDPtr == NULL || CDVarFDPtr == NULL ||
		TypeVarFDPtr == NULL)
		IsOk = FALSE;




	/*
	 * FE Quad and ZoneType_Ordered Data:
	 *                                                            *
	 * 4         3                                                *
	 * +---------+                                                *
	 * |         |                                                *
	 * |         |                                                *
	 * |         |                                                *
	 * +---------+                                                *
	 * 0         1                                                *
	 */

	/* Set Velocities at corner of cell & compute vel min/max */
	if (IsOk)
	{
		ZoneType_e ZoneType = TecUtilZoneGetType(ZoneNum);

		if (ZoneType == ZoneType_FEQuad)
		{
			Index[0] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, 1);
			Index[1] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, 2);
			Index[2] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, 3);
			Index[3] = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, 4);
		}
		/*
		else if(ZoneType == ZoneType_Ordered)
		  {
			Index[0] = IIndex + (JIndex - 1) * IMax;
			Index[1] = IIndex + 1 + (JIndex - 1) * IMax;
			Index[2] = IIndex + 1 + JIndex * IMax;
			Index[3] = IIndex + JIndex * IMax;
		  } */
		else
			IsOk = FALSE;
	}

	// Look for degenerate elements (triangles),
	// Adjust Index[] so first three items are Nodes of triangle
	NumNodesInElem = 4;
	if (Index[3] == Index[0] || Index[2] == Index[3])
	{
		IsQuad = FALSE;
		NumNodesInElem--;
	}
	if (Index[2] == Index[1])
	{
		IsQuad = FALSE;
		Index[2] = Index[3];
		NumNodesInElem--;
	}
	if (Index[1] == Index[0])
	{
		IsQuad = FALSE;
		Index[1] = Index[2];
		Index[2] = Index[3];
		NumNodesInElem--;
	}
	// Can only handle Quads or Triangles
	if (NumNodesInElem < 3) IsOk = FALSE;

	// Determine PsiAlignment for cell.
	if (IsOk)
	{
		int ic;
		double XCell[4], YCell[4], ZCell[4];

		// double x1 = TecUtilDataValueGetByRef(XVarFDPtr, Index[0]);
		// double x2 = TecUtilDataValueGetByRef(XVarFDPtr, Index[1]);
		// double x3 = TecUtilDataValueGetByRef(XVarFDPtr, Index[2]);
		// double x4 = TecUtilDataValueGetByRef(XVarFDPtr, Index[3]);
		// double y1 = TecUtilDataValueGetByRef(YVarFDPtr, Index[0]);
		// double y2 = TecUtilDataValueGetByRef(YVarFDPtr, Index[1]);
		// double y3 = TecUtilDataValueGetByRef(YVarFDPtr, Index[2]);
		// double y4 = TecUtilDataValueGetByRef(YVarFDPtr, Index[3]);
		// double z1 = TecUtilDataValueGetByRef(ZVarFDPtr, Index[0]);
		// double z2 = TecUtilDataValueGetByRef(ZVarFDPtr, Index[1]);
		// double z3 = TecUtilDataValueGetByRef(ZVarFDPtr, Index[2]);
		// double z4 = TecUtilDataValueGetByRef(ZVarFDPtr, Index[3]);

		double XAve, YAve, ZAve;

		// double dxdpsi, dxdeta, dydpsi, dydeta, dzdpsi, dzdeta, RCellNormalTot;

		// Collect the cell corner coordinates
		for (ic = 0; ic < 4; ic++)
		{
			LgIndex_t Node = Index[ic];
			XCell[ic] = TecUtilDataValueGetByRef(XVarFDPtr, Node);
			YCell[ic] = TecUtilDataValueGetByRef(YVarFDPtr, Node);
			ZCell[ic] = TecUtilDataValueGetByRef(ZVarFDPtr, Node);
			W[ic] = 0.0;
		}


		// Compute surface normal for cell
		/*
		if (IsQuad)
		  {
			XAve = 0.25 * (x1 + x2 + x3 + x4);
			YAve = 0.25 * (y1 + y2 + y3 + y4);
			ZAve = 0.25 * (z1 + z2 + z3 + z4);

			dxdpsi = 0.25 * ( x2 - x1 + x3 - x2);
			dxdeta = 0.25 * ( x3 - x2 + x4 - x1);
			dydpsi = 0.25 * ( y2 - y1 + y3 - y2);
			dydeta = 0.25 * ( y3 - y2 + y4 - y1);
			dzdpsi = 0.25 * ( z2 - z1 + z3 - z2);
			dzdeta = 0.25 * ( z3 - z2 + z4 - z1);
		  }
		else   // collapsed to triangle
		  {
			double OneThrd = 1.0/3.0;

			XAve = OneThrd * (x1 + x2 + x3);
			YAve = OneThrd * (y1 + y2 + y3);
			ZAve = OneThrd * (z1 + z2 + z3);

			dxdpsi = 0.5 * ( x2 - x1);
			dxdeta = 0.5 * ( x3 - x1);
			dydpsi = 0.5 * ( y2 - y1);
			dydeta = 0.5 * ( y3 - y1);
			dzdpsi = 0.5 * ( z2 - z1);
			dzdeta = 0.5 * ( z3 - z1);
		  }

		CellNormal.X = dydpsi * dzdeta - dydeta * dzdpsi;
		CellNormal.Y = dzdpsi * dxdeta - dzdeta * dxdpsi;
		CellNormal.Z = dxdpsi * dydeta - dxdeta * dydpsi;
		RCellNormalTot = 1.0 / sqrt(CellNormal.X * CellNormal.X + CellNormal.Y * CellNormal.Y +
									CellNormal.Z * CellNormal.Z);
		CellNormal.X *= RCellNormalTot;
		CellNormal.Y *= RCellNormalTot;
		CellNormal.Z *= RCellNormalTot;

		// If normal most closely aligns with Z, PsiAlignment is 1.
		if ( ABS(CellNormal.Z) >= ABS(CellNormal.X) && ABS(CellNormal.Z) >= ABS(CellNormal.Y) )
		  PsiAlign = 1;
		// If normal most closely aligns with X, PsiAlignment is 2.
		else if ( ABS(CellNormal.X) >= ABS(CellNormal.Y) && ABS(CellNormal.X) > ABS(CellNormal.Z) )
		  PsiAlign = 2;
		// If normal most closely aligns with Y, PsiAlignment is 3.
		else if ( ABS(CellNormal.Y) > ABS(CellNormal.X) && ABS(CellNormal.Y) > ABS(CellNormal.Z) )
		  PsiAlign = 3;
		  */

		PsiAlign = CalcPsiAlign(IsQuad, XCell, YCell, ZCell, &XAve, &YAve, &ZAve, &CellNormal);

		// Needed later in curvature calculation
		for (ic = 0; ic < 4; ic++)
		{
			DXCell[ic] = XCell[ic] - XAve;
			DYCell[ic] = YCell[ic] - YAve;
			DZCell[ic] = ZCell[ic] - ZAve;
		}
		/*
		DXCell[0] = x1 - XAve;
		DXCell[1] = x2 - XAve;
		DXCell[2] = x3 - XAve;
		DXCell[3] = x4 - XAve;
		DYCell[0] = y1 - YAve;
		DYCell[1] = y2 - YAve;
		DYCell[2] = y3 - YAve;
		DYCell[3] = y4 - YAve;
		DZCell[0] = z1 - ZAve;
		DZCell[1] = z2 - ZAve;
		DZCell[2] = z3 - ZAve;
		DZCell[3] = z4 - ZAve;
		*/
	}

	//	Check for CP. If found, agitate node values and see if CP remains. 
	//	If so, then use original values for CP
	
	for (int CheckNum = 0 ; IsOk && CheckNum < AgitationNumber ; CheckNum++)
	{

		// Compute the surface-tangent velocities (UPCell, VPCell) at the cell corners
		if (IsOk)
		{
			for (ic = 0; IsOk && ic < NumNodesInElem; ic++)
			{
				LgIndex_t     Node = Index[ic];
				XYZ_s         Normal = NormalsGetNormalForNode(Normals, Node);
				UCell[ic] = TecUtilDataValueGetByRef(UVarFDPtr, Node);
				VCell[ic] = TecUtilDataValueGetByRef(VVarFDPtr, Node);
				WCell[ic] = TecUtilDataValueGetByRef(WVarFDPtr, Node);

				if (CellFound && CheckNum >= 1 && CheckNum < AgitationNumber - 1)
				{
					double RandNum[3];
					for (int j = 0 ; j < 3 ; j++)
						RandNum[j] = (-1.0 + 0.0001 * (double)(rand() % 9999 + 1));
					if (*TypeCrtPt == 0)
					{
						UCell[ic] += RandNum[0] * AgitationFactor * SaddleFactor * *GradMagCrtPt;
						VCell[ic] += RandNum[1] * AgitationFactor * SaddleFactor * *GradMagCrtPt;
						WCell[ic] += RandNum[2] * AgitationFactor * SaddleFactor * *GradMagCrtPt;
					}
					else
					{
						UCell[ic] += RandNum[0] * AgitationFactor * *GradMagCrtPt;
						VCell[ic] += RandNum[1] * AgitationFactor * *GradMagCrtPt;
						WCell[ic] += RandNum[2] * AgitationFactor * *GradMagCrtPt;
					}
				}

				// Using PsiAlign for the cell and Normal for the Node, compute surface-tangent gradients
				IsOk = XYZtoPsiEta(PsiAlign, Normal, UCell[ic], VCell[ic], WCell[ic], &(UPCell[ic]), &(VPCell[ic]));
				IsOk = XYZtoPsiEta(PsiAlign, Normal, DXCell[ic], DYCell[ic], DZCell[ic], &(XPCell[ic]), &(YPCell[ic]));
			}
		}
		
		/* Compute the minima and maxima of surface-tangent gradients in the cell */
		if (IsOk)
		{
			UPCMin = UPCell[0];
			VPCMin = VPCell[0];
			UPCMax = UPCell[0];
			VPCMax = VPCell[0];
			for (ic = 1; ic < NumNodesInElem; ic++)
			{
				UPCMin = MIN(UPCMin, UPCell[ic]);
				UPCMax = MAX(UPCMax, UPCell[ic]);
				VPCMin = MIN(VPCMin, VPCell[ic]);
				VPCMax = MAX(VPCMax, VPCell[ic]);
			}
		}

		if (IsOk)
		{
			// if (UPCMin <= 0.0 && UPCMax > 0.0 &&
			//     VPCMin <= 0.0 && VPCMax > 0.0)
			if (1)
			{
				/* Maxima */
				// if (UPCell[0] >  0.0  &&  VPCell[0] >  0.0 &&
				//     UPCell[1] <= 0.0  &&  VPCell[1] >  0.0 &&          TODO: Only applies for othogonal grids.
				//     UPCell[2] <= 0.0  &&  VPCell[2] <= 0.0 &&                Need better culling criteria.
				//     UPCell[3] >  0.0  &&  VPCell[3] <= 0.0 )
				// if (IsOk)
				//   {
				//     /* Find corner with maximum of scalar variable */
				//     double RhoMax = TecUtilDataValueGetByRef(CDVarFDPtr, Index[0]);
				//     double rr[] = {-0.8,  0.8,  0.8, -0.8,};
				//     double ss[] = {-0.8, -0.8,  0.8,  0.8,};
				//     int    icmax = 0;
				//     for (ic=1; ic<NumNodesInElem; ic++)
				//      {
				//         double RhoCrn = TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]);
				//         if ( RhoCrn > RhoMax )
				//           {
				//             RhoMax = RhoCrn;
				//             icmax = ic;
				//           }
				//       }
				//     r = rr[icmax];
				//     s = ss[icmax];
				//   }
				// else
				{
					r = 0.0;
					s = 0.0;
				}

				if (NumNodesInElem == 4)   // Quad
				{
					CellFound =  BilinearQuadNaturalCoord(0.0, 0.0, UPCell, VPCell,
														  UPCMin, VPCMin, UPCMax, VPCMax,
														  &r, &s, W);
				}
				else if (NumNodesInElem == 3)  // Triangle
				{
					// CellFound = LinearTriangleNaturalCoord(0.0, 0.0, UPCell, VPCell,
					//                               UPCMin, VPCMin, UPCMax, VPCMax,
					//                               &r, &s, W);
					CellFound = LinearTriangleNaturalCoord(0.0, 0.0, UPCell, VPCell,
														   UPCMin, VPCMin, UPCMax, VPCMax, XPCell, YPCell,
														   &r, &s, W);
				}
			}
		}

		if (IsOk)
		{

			// TEMP
			// if (CellFound)
			//   {
			//     char Message[200];
			//     sprintf(Message, "After first test: r,s=%g,%g", r, s);
			//     TecUtilDialogMessageBox(Message, MessageBox_Warning);
			//   }

			/* If Newton iteration failed, try some other initial conditions for r,s */
			// TODO: This looks wrong - doesn't set r=rr[ictry], etc.
			if (!CellFound && NumNodesInElem == 4)
			{
				double rr[] = { -0.8,  0.8,  0.8, -0.8};
				double ss[] = { -0.8, -0.8,  0.8,  0.8};
				int ictry = 0;
				while (!CellFound && ictry < 4)
				{
					// if (NumNodesInElem == 4)   // Quad
					//   {
					CellFound =  BilinearQuadNaturalCoord(0.0, 0.0, UPCell, VPCell,
														  UPCMin, VPCMin, UPCMax, VPCMax,
														  &r, &s, W);
					//   }
					// else if (NumNodesInElem == 3)  // Triangle
					//   {
					//     CellFound = LinearTriangleNaturalCoord(0.0, 0.0, UPCell, VPCell,
					//                                   UPCMin, VPCMin, UPCMax, VPCMax, XPCell, YPCell,
					//                                   &r, &s, W);
					//   }
					if (!CellFound)
						ictry++;
				}
			}

			/* If Newton iteration failed, one last test for max */
			if (!CellFound)
			{
				LgIndex_t NodeMx;
				double    CellMaxChrgDens = TecUtilDataValueGetByRef(CDVarFDPtr, Index[0]);
				int       icmax = 0;
				int       Type;

				/* Find the corner of the cell with the largest charge density */
				for (ic = 1; ic < NumNodesInElem; ic++)
				{
					double ChrgDensCrn = TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]);
					if (ChrgDensCrn > CellMaxChrgDens)
					{
						icmax = ic;
						CellMaxChrgDens = ChrgDensCrn;
					}
				}


				NodeMx = Index[icmax];

				Type = (int)TecUtilDataValueGetByRef(TypeVarFDPtr, NodeMx);


				// >>>> TODO: the following is new and needs to be tested

				/*
				 * Max if current point larger than value at all nodes in surrounding elements and this
				 * max point hasn't already been found while testing another cell (Type>-2).
				 */

				#if(0)
					// if ( Type > -2 )
				{
					Boolean_t MaxSoFar = TRUE;
					Boolean_t MinSoFar = TRUE;
					LgIndex_t NodeListOffset;
					LgIndex_t ElemOffset;
					double ChargeDensIni = TecUtilDataValueGetByRef(CDVarFDPtr, NodeMx);
					double CDNeighbor;

					ArrList_pa NeighborNodes = ArrListAlloc(20, ArrListType_Long);
					LgIndex_t NumNeighborElems = SurfElemMapGetElemCountForNode(SurfElemMap, NodeMx);


					// For each node in element, make a list of nodes in surrounding element
					// First add the current node to the list, then delete it at the end
					ArrListAppendUniqueLongItem(NeighborNodes, NodeMx);
					for (ElemOffset = 0; ElemOffset < NumNeighborElems; ElemOffset++)
					{
						LgIndex_t ElemNum = SurfElemMapGetElemNumForNodeOffset(SurfElemMap, NodeMx, ElemOffset);

						for (ic = 1; ic < 4; ic++)
						{
							LgIndex_t Node = TecUtilDataNodeGetByZone(ZoneNum, ElemNum, ic);
							ArrListAppendUniqueLongItem(NeighborNodes, Node);
						}
					}

					// Now go through list and see if any node has a larger value than NodeMx
					for (NodeListOffset = 0; MaxSoFar && NodeListOffset < ArrListGetCount(NeighborNodes); NodeListOffset++)
					{
						LgIndex_t Node;
						ArrListItem_u Item = ArrListGetItem(NeighborNodes, NodeListOffset);
						Node = Item.Long;
						CDNeighbor = TecUtilDataValueGetByRef(CDVarFDPtr, Node);
						if (CDNeighbor > ChargeDensIni) MaxSoFar = FALSE;
					}

					if (MaxSoFar)
					{
						CellFound = TRUE;
						for (ic = 0; ic < 4; ic++) W[ic] = 0.0;
						W[icmax] = 1.0;
						TecUtilDataValueSetByRef(TypeVarFDPtr, NodeMx, -3);
						r = -1.0;
						if (icmax == 2 && icmax == 3) r = 1.0;
						s = -1.0;
						if (icmax >= 3) s = 1.0;
					}

					ArrListDealloc(&NeighborNodes);
				}
				#endif
			}
		}

		if (IsOk)
		{

			/* If Newton iteration failed, one last test for min */
			if (!CellFound)
			{
				LgIndex_t NodeMn;
				double    CellMaxChrgDens = TecUtilDataValueGetByRef(CDVarFDPtr, Index[0]);
				int       icmin = 0;
				int       Type;

				/* Find the corner of the cell with the smallest charge density */
				for (ic = 1; ic < NumNodesInElem; ic++)
				{
					double ChrgDensCrn = TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]);
					if (ChrgDensCrn > CellMaxChrgDens)
					{
						icmin = ic;
						CellMaxChrgDens = ChrgDensCrn;
					}
				}

				NodeMn = Index[icmin];

				Type = (int)TecUtilDataValueGetByRef(TypeVarFDPtr, NodeMn);

				/*
				 * Min if current point smaller than the value at all nodes in surrounding elements, and this
				 * min point hasn't already been found while testing another cell (Type<2).
				 */
				/* TODO
				if ( Type < 2 )
				  {
					Boolean_t MinSoFar = TRUE;
					LgIndex_t ii, jj;
					LgIndex_t Indx;
					double ChargeDensIni = TecUtilDataValueGetByRef(CDVarFDPtr, IndexMn);
					double CDNeighbor;
					for (jj=JMn-1; MinSoFar && jj<=JMn+1; jj++)
					  {
					  for (ii=IMn-1; MinSoFar && ii<=IMn+1; ii++)
						  {
							if ( !(ii == IMn && jj == JMn) )
							  {
								Indx = ii + (jj - 1) * IMax;
								CDNeighbor = TecUtilDataValueGetByRef(CDVarFDPtr, Indx);
								if (CDNeighbor < ChargeDensIni) MinSoFar = FALSE;
							 }
						  }
					  }
					if (MinSoFar)
					  {
						CellFound = TRUE;
						for (ic=0; ic<4; ic++) W[ic] = 0.0;
						W[icmin] = 1.0;
						TecUtilDataValueSetByRef(TypeVarFDPtr, Index[icmin], 3);
					  }
				  }
				  */
			}


			// TEMP
			// if (CellFound)
			//   {
			//     char Message[200];
			//     sprintf(Message, "After all test: r,s=%g,%g, W=%g,%g,%g,%g", r, s, W[0], W[1], W[2], W[3]);
			//     TecUtilDialogMessageBox(Message, MessageBox_Warning);
			//   }
		}

		/* Compute XYZ location and type of critical point */
		if (IsOk && CellFound)
		{
			double D2RhoDx2  = 0.0;
			double D2RhoDxDy = 0.0;
			double D2RhoDy2  = 0.0;
			
			*XCrtPt = 0.0;
			*YCrtPt = 0.0;
			*ZCrtPt = 0.0;
			*GradMagCrtPt = 0.0;
			for (ic = 0; ic < 4; ic++)
			{
				*XCrtPt += TecUtilDataValueGetByRef(XVarFDPtr, Index[ic]) * W[ic];
				*YCrtPt += TecUtilDataValueGetByRef(YVarFDPtr, Index[ic]) * W[ic];
				*ZCrtPt += TecUtilDataValueGetByRef(ZVarFDPtr, Index[ic]) * W[ic];
				*GradMagCrtPt += TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]) * W[ic];
			}
			
			/* Compute second derivatives for curvature - classification.
			 * Utilize a bilinear variation of UPCell = DRhoDXp and VPCell = DRhoDYp,
			 * where Xp and Yp are the local surface-tangent coordinates computed using
			 * PsiAlign and Normal.
			 */
			if (IsOk)
			{
				// Find surface-tangent coords of cell corners
				for (ic = 0; IsOk && ic < 4; ic++)
					IsOk = XYZtoPsiEta(PsiAlign, CellNormal, DXCell[ic], DYCell[ic], DZCell[ic], &(XPCell[ic]), &(YPCell[ic]));

				// Compute the curvature in surface-tangent coordinates
				// >>>>>> TODO: Needs testing
				if (NumNodesInElem == 4)   // Quad
				{
					IsOk = ComputeFEQuadCurvature(// ZoneNum, ElemNum, PsiAlign, CellNormal,
							   XPCell, YPCell, UPCell, VPCell, r, s,
							   &D2RhoDx2, &D2RhoDxDy, &D2RhoDy2);
				}
				else if (NumNodesInElem == 3)  // Triangle
				{
					IsOk = ComputeFETriangleCurvature(XPCell, YPCell, UPCell, VPCell,
													  &D2RhoDx2, &D2RhoDxDy, &D2RhoDy2);
				}

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
				a[1][0] = a[0][1];
				a[1][1] = (float)D2RhoDy2;

				// Temporary, because we are doing 2x2 with 3x3
				a[0][2] = 0.0;
				a[2][0] = 0.0;
				a[1][2] = 0.0;
				a[2][1] = 0.0;
				a[2][2] = 1.0;

				// TODO: Problems here! 3x3 versus 2x2
				IsOk = Jacobi3by3(InMat, EgnVal, EgnVec, &NumRot);

				/* Sort by eigenvalue value */
				if (SortEgnV && IsOk)
				{
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
				}

				/* Determine critical point classifications */
				if (EgnVal[1] < 0.0)
				{
					if (EgnVal[0] < 0.0)   /* Maxima */
						CritPointType = -2;
					else                   /* Saddle */
						CritPointType = 0;
				}
				else                           /* Minima */
					CritPointType = 2;

				*TypeCrtPt = CritPointType;

				/* If saddle point, set Principle Direction vector. */
				if (CritPointType == 0)
				{
					double XPrincDir, YPrincDir, ZPrincDir;
					double XPPrincDir = (double)EgnVec[0][0];
					double YPPrincDir = (double)EgnVec[1][0];

					IsOk = PsiEtatoXYZ(PsiAlign, CellNormal, XPPrincDir, YPPrincDir, &XPrincDir, &YPrincDir, &ZPrincDir);
					*PrincDirX = (float)XPrincDir;
					*PrincDirY = (float)YPrincDir;
					*PrincDirZ = (float)ZPrincDir;
				}


				/* If close to node - mark as max or min already found */
				if (IsOk && CellFound)
				{
					double WMax = -1.0;
					int icmax = -1;
					for (ic = 0; ic < 4; ic++)
					{
						if (W[ic] > WMax)
						{
							WMax = W[ic];
							icmax = ic;
						}
					}
					if ((CritPointType == -2 || CritPointType == 2) && WMax > 0.9)
						TecUtilDataValueSetByRef(TypeVarFDPtr, Index[icmax], CritPointType);
				}

				if (InMat!=NULL)
					FREE_ARRAY(InMat, "InMat");
				if (EgnVec!=NULL)
					FREE_ARRAY(EgnVec, "EgnVec");
			}
		}
		else
			CellFound = FALSE;

		if (!CellFound) break;
	}

	return(CellFound);
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
							  Boolean_t   PeriodicBC,
							  double     *XCrtPt,
							  double     *YCrtPt,
							  double     *ZCrtPt,
							  double     *ChrgDensCrtPt,
							  EntIndex_t *TypeCrtPt,
							  float      *PrincDirX,
							  float      *PrincDirY,
							  float      *PrincDirZ,
							  ZoneVarInfo_pa ZoneVarInfo)
{
	Boolean_t IsOk = TRUE;
	Boolean_t CellFound = FALSE;
	LgIndex_t IMax, JMax, KMax;
	LgIndex_t IBeg, IEnd, JBeg, JEnd, KBeg, KEnd;
	double    UCell[8], VCell[8], WCell[8];  /* Gradients at cell corners */
	double    XCell[8], YCell[8], ZCell[8];  /* Coordinates at cell corners */
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
	//	Next two are for data agitation to eliminate spurious CPs
	const double	AgitationFactor = 0.000005, AtomFactor = 0.0, SaddleFactor = 1.0;
	const int		AgitationNumber = 100;

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

	REQUIRE(VALID_BOOLEAN(PeriodicBC));

	//	Initialize random time seed for data agitation
	srand((unsigned int)time(NULL));

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

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

	if (PeriodicBC)
	{
		IBeg = 0;
		IEnd = IMax + 1;
		JBeg = 0;
		JEnd = JMax + 1;
		KBeg = 0;
		KEnd = KMax + 1;
	}
	else
	{
		IBeg = 1;
		IEnd = IMax;
		JBeg = 1;
		JEnd = JMax;
		KBeg = 1;
		KEnd = KMax;
	}
	REQUIRE(IIndex >= IBeg && IIndex <= IEnd);
	REQUIRE(JIndex >= JBeg && JIndex <= JEnd);
	REQUIRE(KIndex >= KBeg && KIndex <= KEnd);

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

	/* Set Gradients at corner of cell & compute grad min/max */
	if (IsOk)
	{

		if (PeriodicBC)
		{
			Index[0] = PeriodicIndexFromIJK(IIndex    , JIndex    , KIndex    , IMax, JMax, KMax);
			Index[1] = PeriodicIndexFromIJK(IIndex + 1, JIndex    , KIndex    , IMax, JMax, KMax);
			Index[2] = PeriodicIndexFromIJK(IIndex + 1, JIndex + 1, KIndex    , IMax, JMax, KMax);
			Index[3] = PeriodicIndexFromIJK(IIndex    , JIndex + 1, KIndex    , IMax, JMax, KMax);
			Index[4] = PeriodicIndexFromIJK(IIndex    , JIndex    , KIndex + 1, IMax, JMax, KMax);
			Index[5] = PeriodicIndexFromIJK(IIndex + 1, JIndex    , KIndex + 1, IMax, JMax, KMax);
			Index[6] = PeriodicIndexFromIJK(IIndex + 1, JIndex + 1, KIndex + 1, IMax, JMax, KMax);
			Index[7] = PeriodicIndexFromIJK(IIndex    , JIndex + 1, KIndex + 1, IMax, JMax, KMax);
		}
		else
		{
			Index[0] = IIndex + (JIndex - 1) * IMax + (KIndex - 1) * IMax * JMax;
			Index[1] = IIndex + 1 + (JIndex - 1) * IMax + (KIndex - 1) * IMax * JMax;
			Index[2] = IIndex + 1 + JIndex * IMax + (KIndex - 1) * IMax * JMax;
			Index[3] = IIndex + JIndex * IMax + (KIndex - 1) * IMax * JMax;
			Index[4] = IIndex + (JIndex - 1) * IMax + KIndex * IMax * JMax;
			Index[5] = IIndex + 1 + (JIndex - 1) * IMax + KIndex * IMax * JMax;
			Index[6] = IIndex + 1 + JIndex * IMax + KIndex * IMax * JMax;
			Index[7] = IIndex + JIndex * IMax + KIndex * IMax * JMax;
		}
	}

	for (int CheckNum = 0 ; IsOk && CheckNum < AgitationNumber ; CheckNum++)
	{		
		for (ic = 0; ic < 8; ic++)
		{
			UCell[ic] = TecUtilDataValueGetByRef(UVarFDPtr, Index[ic]);
			VCell[ic] = TecUtilDataValueGetByRef(VVarFDPtr, Index[ic]);
			WCell[ic] = TecUtilDataValueGetByRef(WVarFDPtr, Index[ic]);
			XCell[ic] = TecUtilDataValueGetByRef(XVarFDPtr, Index[ic]);
			YCell[ic] = TecUtilDataValueGetByRef(YVarFDPtr, Index[ic]);
			ZCell[ic] = TecUtilDataValueGetByRef(ZVarFDPtr, Index[ic]);

			if (CellFound && CheckNum >= 1 && CheckNum < AgitationNumber - 1)
			{
				double RandNum[3];
				for (int j = 0 ; j < 3 ; j++)
					RandNum[j] = (-1.0 + 0.0001 * (double)(rand() % 9999 + 1));
				if (*TypeCrtPt == 1 || *TypeCrtPt == -1)
				{
					UCell[ic] += RandNum[0] * AgitationFactor * SaddleFactor * *ChrgDensCrtPt;
					VCell[ic] += RandNum[1] * AgitationFactor * SaddleFactor * *ChrgDensCrtPt;
					WCell[ic] += RandNum[2] * AgitationFactor * SaddleFactor * *ChrgDensCrtPt;
				}
				else if (*TypeCrtPt == -3)
				{
					UCell[ic] += RandNum[0] * AgitationFactor * AtomFactor * *ChrgDensCrtPt;
					VCell[ic] += RandNum[1] * AgitationFactor * AtomFactor * *ChrgDensCrtPt;
					WCell[ic] += RandNum[2] * AgitationFactor * AtomFactor * *ChrgDensCrtPt;
				}
				else
				{
					UCell[ic] += RandNum[0] * AgitationFactor * *ChrgDensCrtPt;
					VCell[ic] += RandNum[1] * AgitationFactor * *ChrgDensCrtPt;
					WCell[ic] += RandNum[2] * AgitationFactor * *ChrgDensCrtPt;
				}
			}
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
				for (ic = 1; ic < 8; ic++)
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
			// TODO: This looks wrong - doesn't set r=rr[ictry], etc.
			if (!CellFound)
			{
				double rr[] = { -0.8,  0.8,  0.8, -0.8, -0.8,  0.8,  0.8, -0.8};
				double ss[] = { -0.8, -0.8,  0.8,  0.8, -0.8, -0.8,  0.8,  0.8};
				double tt[] = { -0.8, -0.8, -0.8, -0.8,  0.8,  0.8,  0.8,  0.8};
				int ictry = 0;
				while (!CellFound && ictry < 8)
				{
					//	TIM:	set r,s,t to values just inside the corners
					r = rr[ictry];
					s = ss[ictry];
					t = tt[ictry];
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

				if (PeriodicBC)
					IndexMx = PeriodicIndexFromIJK(IMx, JMx, KMx, IMax, JMax, KMax);
				else
					IndexMx = IMx + (JMx - 1) * IMax + (KMx - 1) * IMax * JMax;

				Type = (int)TecUtilDataValueGetByRef(TypeVarFDPtr, IndexMx);

				/*
				 * Max if current point larger than the surrounding 26 points, and this
				 * max point hasn't already been found while testing another cell (Type>3).
				 */
				if (Type > -3 && IMx > IBeg && IMx < IEnd && JMx > JBeg && JMx < JEnd && KMx > KBeg && KMx < KEnd)
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
									if (PeriodicBC)
										Indx = PeriodicIndexFromIJK(ii, jj, kk, IMax, JMax, KMax);
									else
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

				if (PeriodicBC)
					IndexMn = PeriodicIndexFromIJK(IMn, JMn, KMn, IMax, JMax, KMax);
				else
					IndexMn = IMn + (JMn - 1) * IMax + (KMn - 1) * IMax * JMax;

				Type = (int)TecUtilDataValueGetByRef(TypeVarFDPtr, IndexMn);

				/*
				 * Min if current point smaller than the surrounding 26 points, and this
				 * min point hasn't already been found while testing another cell (Type>3).
				 */
				if (Type < 3 && IMn > IBeg && IMn < IEnd && JMn > JBeg && JMn < JEnd && KMn > KBeg && KMn < KEnd)
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
									if (PeriodicBC)
										Indx = PeriodicIndexFromIJK(ii, jj, kk, IMax, JMax, KMax);
									else
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
			// TIM: The exact charge density was computed, as well as the weights to apply to x,y,z in order to
			// yield that density, so use those weights again to find the actual position of the CP
			// within the cell.
			if (CellFound)
			{
				double D2RhoDx2  = 0.0;
				double D2RhoDxDy = 0.0;
				double D2RhoDxDz = 0.0;
				double D2RhoDy2  = 0.0;
				double D2RhoDyDz = 0.0;
				double D2RhoDz2  = 0.0;

				double D2rhoDxDx = 0.0;
				double D2rhoDxDy = 0.0;
				double D2rhoDxDz = 0.0;
				double D2rhoDyDx = 0.0;
				double D2rhoDyDy = 0.0;
				double D2rhoDyDz = 0.0;
				double D2rhoDzDx = 0.0;
				double D2rhoDzDy = 0.0;
				double D2rhoDzDz = 0.0;

				double DrDx, DsDy, DtDz, DxDr, DyDs, DzDt;
				double dwdr[8], dwds[8], dwdt[8];

				*XCrtPt = 0.0;
				*YCrtPt = 0.0;
				*ZCrtPt = 0.0;
				*ChrgDensCrtPt = 0.0;
				for (ic = 0; ic < 8; ic++)
				{
					/*
					XNode = TecUtilDataValueGetByRef(XVarFDPtr, Index[ic]);
					if (XNode < (double)(IIndex))     XNode = XNode + (double)IMax;
					if (XNode > (double)(IIndex + 1)) XNode = XNode - (double)IMax;

					YNode = TecUtilDataValueGetByRef(YVarFDPtr, Index[ic]);
					if (YNode < (double)(JIndex))     YNode = YNode + (double)JMax;
					if (YNode > (double)(JIndex + 1)) YNode = YNode - (double)JMax;

					ZNode = TecUtilDataValueGetByRef(ZVarFDPtr, Index[ic]);
					if (ZNode < (double)(KIndex))     ZNode = ZNode + (double)KMax;
					if (ZNode > (double)(KIndex + 1)) ZNode = ZNode - (double)KMax;
					*/
					// *XCrtPt += TecUtilDataValueGetByRef(XVarFDPtr, Index[ic]) * W[ic];
					// *YCrtPt += TecUtilDataValueGetByRef(YVarFDPtr, Index[ic]) * W[ic];
					// *ZCrtPt += TecUtilDataValueGetByRef(ZVarFDPtr, Index[ic]) * W[ic];

					*XCrtPt += XCell[ic] * W[ic];
					*YCrtPt += YCell[ic] * W[ic];
					*ZCrtPt += ZCell[ic] * W[ic];
					// TIM: ASK does this assume a continuium in the charge density between corners and the point?
					// Yes, but the density changes linearly along the edges of a cell, so since it has different
					// slop along each edge, it is not linear as you pass through the cell.
					*ChrgDensCrtPt += TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]) * W[ic];
				}

				/* Compute second derivatives for curvature - classification.
				 * For now, assume rectangular grid with spacing 1 in each direction.
				 */
				// New - compute the derivatives of the gradients using the derivatives
				// of the shape functions (weights)
				// Grid spacing - assume rectangular with x aligned with r, y aligned with s, and z aligned with t
				DxDr = TecUtilDataValueGetByRef(XVarFDPtr, Index[1]) - TecUtilDataValueGetByRef(XVarFDPtr, Index[0]);
				if (ABS(DxDr) > SMALLFLOAT) DrDx = 1.0 / DxDr;
				DyDs = TecUtilDataValueGetByRef(YVarFDPtr, Index[2]) - TecUtilDataValueGetByRef(YVarFDPtr, Index[1]);
				if (ABS(DyDs) > SMALLFLOAT) DsDy = 1.0 / DyDs;
				DzDt = TecUtilDataValueGetByRef(ZVarFDPtr, Index[4]) - TecUtilDataValueGetByRef(ZVarFDPtr, Index[0]);
				if (ABS(DzDt) > SMALLFLOAT) DtDz = 1.0 / DzDt;

				// compute the derivatives of the natural coordinates
				IsOk = BrickTrilinearWeightDerivatives(r, s, t, W, dwdr, dwds, dwdt);

				double DRhoDx, DRhoDy, DRhoDz;

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

					// Old
					IsOk = ComputeCurvature(ZoneNum, i, j, k, CDVarFDPtr, PeriodicBC,
											&d2dx2, &d2dxdy, &d2dxdz,
											&d2dy2, &d2dydz, &d2dz2);

					D2RhoDx2  += W[ic] * d2dx2;
					D2RhoDxDy += W[ic] * d2dxdy;
					D2RhoDxDz += W[ic] * d2dxdz;
					D2RhoDy2  += W[ic] * d2dy2;
					D2RhoDyDz += W[ic] * d2dydz;
					D2RhoDz2  += W[ic] * d2dz2;

					// New - compute the derivatives of the gradients using the derivatives
					// of the shape functions (weights)
					DRhoDx = TecUtilDataValueGetByRef(UVarFDPtr, Index[ic]);
					DRhoDy = TecUtilDataValueGetByRef(VVarFDPtr, Index[ic]);
					DRhoDz = TecUtilDataValueGetByRef(WVarFDPtr, Index[ic]);
					//Assume rectangular grid with r aligned with x, s aligned with y, and t alighned with z
					D2rhoDxDx += dwdr[ic] * DRhoDx * DrDx;
					D2rhoDxDy += dwds[ic] * DRhoDx * DsDy;
					D2rhoDxDz += dwdt[ic] * DRhoDx * DtDz;
					D2rhoDyDx += dwdr[ic] * DRhoDy * DrDx;
					D2rhoDyDy += dwds[ic] * DRhoDy * DsDy;
					D2rhoDyDz += dwdt[ic] * DRhoDy * DtDz;
					D2rhoDzDx += dwdr[ic] * DRhoDz * DrDx;
					D2rhoDzDy += dwds[ic] * DRhoDz * DsDy;
					D2rhoDzDz += dwdt[ic] * DRhoDz * DtDz;

				}
				// Temp
				// D2RhoDx2  = D2rhoDxDx;
				// D2RhoDxDy = D2rhoDxDy;
				// D2RhoDxDz = D2rhoDxDz;
				// D2RhoDy2  = D2rhoDyDy;
				// D2RhoDyDz = D2rhoDyDz;
				// D2RhoDz2  = D2rhoDzDz;

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
					a[1][0] = (float)D2RhoDxDy;//a[0][1];
					a[1][1] = (float)D2RhoDy2;
					a[1][2] = (float)D2RhoDyDz;
					a[2][0] = (float)D2RhoDxDz;//a[0][2];
					a[2][1] = (float)D2RhoDyDz;//a[1][2];
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
					float theta;
					float phi;
					const float pi = 3.141592653589793;
					std::string Type;
					if (EgnVal[2] < 0.0)
					{
						if (EgnVal[1] < 0.0)
						{
							if (EgnVal[0] < 0.0){   /* Atom */
								CritPointType = -3;
								Type = "Atom";
								theta = 0;
								phi = 0;
							}
							else{                   /* Bond */
								CritPointType = -1;
								Type = "Bond";
								theta = atan(sqrt(abs(EgnVal[0]/EgnVal[1])));
								phi = atan(sqrt(abs(EgnVal[0]/EgnVal[2])));
								// LOG
								
							}
						}
						else{                       /* Ring */
							CritPointType = 1;
							Type = "Ring";
							theta = atan(sqrt(abs(EgnVal[2]/EgnVal[1])));
							phi = atan(sqrt(abs(EgnVal[2]/EgnVal[0])));
						}
					}
					else{                           /* Cage */
						CritPointType = 3;
						Type = "Cage";
						theta = atan(sqrt(abs(EgnVal[0]/EgnVal[1])));
						phi = atan(sqrt(abs(EgnVal[0]/EgnVal[2])));
					}

					if (CritPointType != -3 && CheckNum == AgitationNumber - 1 && ZoneVarInfo->LogBondInfo){ //bond has passed data agitation tests
						/*
						 *	Bondalyzer finds bogus curvature values, and rather than modify/break the code trying
						 *	to get legit values, I'll just calculate the curvature and eigensystem again with
						 *	my own functions
						 */
						vec3 Point(vector<double>({ *XCrtPt, *YCrtPt, *ZCrtPt }).data()), EigVals, Grad;
						VolExtentIndexWeights_s VolInfo;
						mat33 Hess, EigVecs;
						FieldDataPointer_c RhoPtr;
						vector<FieldDataPointer_c> VarReadPtrs(3);
						GPType_e CalcType = GPType_Classic;
						MultiRootParams_s Params;

						RhoPtr.GetReadPtr(ZoneNum, ChrgDensVarNum);
						VarReadPtrs[0].GetReadPtr(ZoneNum, UVarNum);
						VarReadPtrs[1].GetReadPtr(ZoneNum, VVarNum);
						VarReadPtrs[2].GetReadPtr(ZoneNum, WVarNum);

						VolInfo.AddOnID = AddOnID;
						VolInfo.IsPeriodic = PeriodicBC;
						vector<EntIndex_t> XYZVarNums(3, -1);
						TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
						if (XYZVarNums[0] <= 0 || XYZVarNums[1] <= 0 || XYZVarNums[2] <= 0){
							vector<vector<string> > TmpStrs = {
								{ "X", "Y", "Z" }, { "I", "J", "K" }
							};
							Boolean_t VarsFound = FALSE;
							for (int i = 0; i < 2 && !VarsFound; ++i){
								for (int j = 0; j < 3; ++j)
									XYZVarNums[j] = TecUtilVarGetNumByName(TmpStrs[i][j].c_str());
								VarsFound = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2] > 0);
							}
							if (!VarsFound){
								TecUtilDialogErrMsg("Failed to get x, y, z variable numbers.");
								TecUtilLockFinish(AddOnID);
								return FALSE;
							}
						}
						for (int i = 0; i < 3; ++i){
							TecUtilVarGetMinMax(XYZVarNums[i], &VolInfo.MinXYZ[i], &VolInfo.MaxXYZ[i]);
						}

						GetVolInfo(ZoneNum, XYZVarNums, PeriodicBC, VolInfo);

						Params.CalcType = GPType_Classic;
						Params.BasisVectors = &VolInfo.BasisNormalized;
						Params.HasHess = FALSE;
						Params.IsPeriodic = PeriodicBC;
						Params.RhoPtr = &RhoPtr;
						Params.GradPtrs = &VarReadPtrs;
						Params.HessPtrs = NULL;
						Params.VolInfo = new VolExtentIndexWeights_s;
						*Params.VolInfo = VolInfo;

						CalcGradForPoint(Point, VolInfo.DelXYZ, VolInfo, Params.BasisVectors, 0, PeriodicBC, Grad, RhoPtr, CalcType, reinterpret_cast<void *>(&Params));
						CalcHessFor3DPoint(Point, VolInfo.DelXYZ, VolInfo, PeriodicBC, Hess, VarReadPtrs, CalcType, reinterpret_cast<void *>(&Params));
						CalcEigenSystemForPoint(Point, EigVals, EigVecs, Params);

						switch (CritPointType){
							case 3:
								theta = atan(sqrt(abs(EigVals[0] / EigVals[1])));
								phi = atan(sqrt(abs(EigVals[0] / EigVals[2])));
								break;
							case 1:
								theta = atan(sqrt(abs(EigVals[2] / EigVals[1])));
								phi = atan(sqrt(abs(EigVals[2] / EigVals[0])));
								break;
							case -1:
								theta = atan(sqrt(abs(EigVals[0] / EigVals[1])));
								phi = atan(sqrt(abs(EigVals[0] / EigVals[2])));
								break;
						}

						std::ofstream* BondLog = ZoneVarInfo->BondLog;
						// Find theta and phi, the angles of the bond's directionality.
						// Eigenvalues are sorted largest to smallest, and correspond with order of eigenvectors.

						// I,J,K,rho,theta [rad],phi [rad],theta [degree],phi [degree],lambda 1,lambda 2,lambda 3,3x3 eigenmatrix [1][1],[1][2],[1][3],[2][1],[2][2],[2][3],[3][1],[3][2],[3][3]

						PrintTimeDate(*BondLog);
						*BondLog << ",," << Type << "," << *XCrtPt << "," << *YCrtPt << "," << *ZCrtPt 
							<< "," << *ChrgDensCrtPt << "," << theta << "," << phi
							<< "," << theta * 180.0/pi << "," << phi * 180.0/pi;
						for (int erow = 0; erow < 3; ++erow){
// 							*BondLog << "," << EgnVal[erow];
							*BondLog << "," << EigVals[erow];
						}
						for (int erow = 0 ; erow < 3 ; ++erow)
							for (int ecol = 0; ecol < 3; ++ecol){
// 								*BondLog << "," << EgnVec.at(erow, ecol);
								*BondLog << "," << EigVecs.at(erow, ecol);
							}
// 						*BondLog << "," << DRhoDx << "," << DRhoDy << "," << DRhoDz
// 							<< "," << sqrt(DRhoDx*DRhoDx + DRhoDy*DRhoDy + DRhoDz*DRhoDz);
						*BondLog << "," << Grad[0] << "," << Grad[1] << "," << Grad[2]
							<< "," << sqrt(Grad[0] * Grad[0] + Grad[1] * Grad[1] + Grad[2] * Grad[2]);
						for (int arow = 0 ; arow < 3 ; ++arow)
							for (int acol = 0 ; acol < 3 ; ++acol)
								if (acol <= arow){
// 									*BondLog << "," << a.at(arow, acol);
									*BondLog << "," << Hess.at(arow, acol);
								}
						*BondLog << "\n";
					}

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

		if (!CellFound) break;
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
		if (IsValid) IsValid = (Count == (CritPoints->NumCrtPtsM3   + CritPoints->NumCrtPtsM1 +
											  CritPoints->NumCrtPtsP1   + CritPoints->NumCrtPtsP3 +
											  CritPoints->NumFFCrtPtsP1 + CritPoints->NumFFCrtPtsP3));
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
		if ((*CritPoints)->MemUsageReported > 0)
		{
			TecUtilMemoryChangeNotify(- (Int64_t)((*CritPoints)->MemUsageReported));
			(*CritPoints)->MemUsageReported = 0;
		}

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
	CritPoints->NumFFCrtPtsP1 = 0;
	CritPoints->NumFFCrtPtsP3 = 0;

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
	REQUIRE((CritPoints->Dimensions == 3 && (Type == -3 || Type == -1 ||
											 Type ==  1 || Type ==  3 ||
											 Type == 11 || Type == 13)) ||
			(CritPoints->Dimensions == 2 &&  Type >= -3 && Type <= 3));

	/* Compute Offsets to the various critical point types */
	if (Type > -2) BegOffset += CritPoints->NumCrtPtsM3;
	if (Type >  0) BegOffset += CritPoints->NumCrtPtsM1;
	if (Type >  1) BegOffset += CritPoints->NumCrtPtsP1;
	/* For Farfield Critical Points */
	if (Type >  3) BegOffset += CritPoints->NumCrtPtsP3;
	if (Type > 11) BegOffset += CritPoints->NumFFCrtPtsP1;

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
	REQUIRE((CritPoints->Dimensions == 3 && (Type == -3 || Type == -1 ||
											 Type ==  1 || Type ==  3 ||
											 Type == 11 || Type == 13)) ||
			(CritPoints->Dimensions == 2 &&  Type >= -3 && Type <= 3));

	/* Compute Offsets to the various critical point types */
	if (Type > -2) EndOffset += CritPoints->NumCrtPtsM1;
	if (Type >  0) EndOffset += CritPoints->NumCrtPtsP1;
	if (Type >  1) EndOffset += CritPoints->NumCrtPtsP3;
	/* For Farfield Critical Points */
	if (Type >  3) EndOffset += CritPoints->NumFFCrtPtsP1;
	if (Type > 11) EndOffset += CritPoints->NumFFCrtPtsP3;

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
 *     Critical point Type.
 *           For Volume, Type is  Max     (Atom) = -3,
 *                                Saddle1 (Bond) = -1,
 *                                Saddle2 (Ring) =  1, and
 *                                Min     (Cage) =  3,
 *                                FFRing         = 11,
 *                                FFCage         = 13.
 *           For Surface, Type is Max    = -2,
 *                                Saddle =  0,
 *                                Min    =  2.
 */
char CritPointsGetType(const CritPoints_pa CritPoints,
					   const LgIndex_t     PointOffset)
{
	char Type = -3;

	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(0 <= PointOffset && PointOffset < CritPointsGetCount(CritPoints));

	if (CritPoints->Dimensions == 3)
	{
		/* Compute type based on Offsets */
		if (PointOffset >= CritPointsGetBegOffset(CritPoints, -1))
		{
			if (PointOffset >= CritPointsGetBegOffset(CritPoints, 1))
			{
				if (PointOffset >= CritPointsGetBegOffset(CritPoints, 3))
				{
					if (PointOffset >= CritPointsGetBegOffset(CritPoints, 11))
					{
						if (PointOffset >= CritPointsGetBegOffset(CritPoints, 13))
							Type = 13;
						else
							Type = 11;
					}
					else
						Type = 3;
				}
				else
					Type = 1;
			}
			else
				Type = -1;
		}
	}
	else if (CritPoints->Dimensions == 2)
	{
		/* Compute type based on Offsets */
		Type = -2;
		if (PointOffset >= CritPointsGetBegOffset(CritPoints, -1))
		{
			if (PointOffset >= CritPointsGetBegOffset(CritPoints, 1))
			{
				if (PointOffset >= CritPointsGetBegOffset(CritPoints, 3))
					Type = 2;
				else
					Type = 1;  // shouldn't happen
			}
			else
				Type = 0;
		}

	}

	ENSURE((CritPoints->Dimensions == 3 && (Type == -3 || Type == -1 ||
											Type ==  1 || Type ==  3 ||
											Type == 11 || Type == 13)) ||
		   (CritPoints->Dimensions == 2 &&  Type >= -3 && Type <= 3));
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
	LgIndex_t OffsetM1, OffsetP1, OffsetP3, OffsetFFP1, OffsetFFP3, NumCrtPts;

	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(0 <= PointOffset && PointOffset <= CritPointsGetCount(CritPoints) - 1);

	/* Compute Offsets to the various critical point types */
	NumCrtPts  = CritPoints->NumCrtPts;
	OffsetM1   = CritPointsGetBegOffset(CritPoints, -1);
	OffsetP1   = CritPointsGetBegOffset(CritPoints,  1);
	OffsetP3   = CritPointsGetBegOffset(CritPoints,  3);

	// So test for decrementing CrtPt count will work
	OffsetFFP1 = CritPointsGetEndOffset(CritPoints,  3) + 1;

	/* Far-Field CPs are 3D only */
	if (CritPoints->Dimensions == 3)
	{
		OffsetFFP1 = CritPointsGetBegOffset(CritPoints, 11);
		OffsetFFP3 = CritPointsGetBegOffset(CritPoints, 13);
	}

	/* Decrement the appropriate CrtPt count */
	(CritPoints->NumCrtPts)--;
	if (PointOffset <  OffsetM1)(CritPoints->NumCrtPtsM3)--;
	if (OffsetM1 <= PointOffset && PointOffset <   OffsetP1)(CritPoints->NumCrtPtsM1)--;
	if (OffsetP1 <= PointOffset && PointOffset <   OffsetP3)(CritPoints->NumCrtPtsP1)--;
	if (OffsetP3 <= PointOffset && PointOffset < OffsetFFP1)(CritPoints->NumCrtPtsP3)--;

	/* Far-Field CPs are 3D only */
	if (CritPoints->Dimensions == 3)
	{
		if (OffsetFFP1 <= PointOffset && PointOffset < OffsetFFP3)(CritPoints->NumFFCrtPtsP1)--;
		if (OffsetFFP3 <= PointOffset)(CritPoints->NumFFCrtPtsP3)--;
	}

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
		Result->NumCrtPts     = 0;
		Result->NumCrtPtsM3   = 0;
		Result->NumCrtPtsM1   = 0;
		Result->NumCrtPtsP1   = 0;
		Result->NumCrtPtsP3   = 0;
		Result->NumFFCrtPtsP1 = 0;
		Result->NumFFCrtPtsP3 = 0;
		Result->X             = ArrListAlloc(60, ArrListType_Double);
		Result->Y             = ArrListAlloc(60, ArrListType_Double);
		Result->Z             = ArrListAlloc(60, ArrListType_Double);
		Result->ChrgDens      = ArrListAlloc(60, ArrListType_Double);
		Result->Type          = ArrListAlloc(60, ArrListType_Char);
		Result->PrincDirX     = ArrListAlloc(60, ArrListType_Double);
		Result->PrincDirY     = ArrListAlloc(60, ArrListType_Double);
		Result->PrincDirZ     = ArrListAlloc(60, ArrListType_Double);

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
			case -3:   // 3D max
				(CritPoints->NumCrtPtsM3)++;
				break;
			case -2:   // 2D max
				(CritPoints->NumCrtPtsM3)++;
				break;
			case -1:   // 3D saddle, type 1
				(CritPoints->NumCrtPtsM1)++;
				break;
			case 0:    // 2D saddle
				(CritPoints->NumCrtPtsM1)++;
				break;
			case 1:    // 3D saddle, type 2
				(CritPoints->NumCrtPtsP1)++;
				break;
			case 2:    // 2D minimum
				(CritPoints->NumCrtPtsP3)++;
				break;
			case 3:    // 3D minimum
				(CritPoints->NumCrtPtsP3)++;
				break;
			case 11:   // 3D Far-Field Ring ("Saddle")
				(CritPoints->NumFFCrtPtsP1)++;
				break;
			case 13:   // 3D Far-Field Cage ("Min")
				(CritPoints->NumFFCrtPtsP1)++;
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
			case -3:   // 3D max
				(CritPoints->NumCrtPtsM3)++;
				break;
			case -2:   // 2D max
				(CritPoints->NumCrtPtsM3)++;
				break;
			case -1:   // 3D saddle, type 1
				(CritPoints->NumCrtPtsM1)++;
				break;
			case 0:    // 2D saddle
				(CritPoints->NumCrtPtsM1)++;
				break;
			case 1:    // 3D saddle, type 2
				(CritPoints->NumCrtPtsP1)++;
				break;
			case 2:    // 2D minimum
				(CritPoints->NumCrtPtsP3)++;
				break;
			case 3:    // 3D minimum
				(CritPoints->NumCrtPtsP3)++;
				break;
			case 11:   // 3D Far-Field Ring ("Saddle")
				(CritPoints->NumFFCrtPtsP1)++;
				break;
			case 13:   // 3D Far-Field Cage ("Max")
				(CritPoints->NumFFCrtPtsP3)++;
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
Boolean_t RemDuplicateCritPoints(const Boolean_t     TestForNode,
								 const double        Tolerance,
								 const double        GridSpacing,
								 const char          Type,
								 CritPoints_pa       CritPoints)
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

			XDistFromNode = GridSpacing * ROUND2(Xtst/GridSpacing) - Xtst;
			YDistFromNode = GridSpacing * ROUND2(Ytst/GridSpacing) - Ytst;
			ZDistFromNode = GridSpacing * ROUND2(Ztst/GridSpacing) - Ztst;
			DistFromNodeSq = XDistFromNode * XDistFromNode + YDistFromNode * YDistFromNode
							 + ZDistFromNode * ZDistFromNode;

			/* If TestForNode, remove only critical points that are at the node. */
			if (!TestForNode || DistFromNodeSq < 1.0e-7)
			{
				Boolean_t DupFound = FALSE;
				// TODO
				// for (ii = BeginCrtPt; !DupFound && ii < CritPointsGetEndOffset(CritPoints, Type); ii++)
				for (ii = 0; !DupFound && ii < CritPointsGetEndOffset(CritPoints, Type); ii++)
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
						if (ABS(DelX) < Tolerance * GridSpacing && ABS(DelY) < Tolerance * GridSpacing && ABS(DelZ) < Tolerance  * GridSpacing)
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









/**
 * Compute the tangential derivatives of a scalar variable on a FE-surface.
 *
 * param ZoneNum
 *     Number of target FE-surface zone number.
 * param VarNum
 *     Number of the scalar variable to be differentiated.
 * param SurfElemMap
 *     SurfElemMap for the target surface zone.
 * param Normals
 *     Normal vectors (nodal) for the target surface zone.
 * param PsiDeriv, EtaDeriv
 *     Partial derivatives (nodal) in the surface tangent psi, eta directions.
 * param PsiAlignment
 *     Coordinate direction (1==X, 2==Y, 3==Z) that most closely aligns with
 *     the psi surface-tangent coordinate direction. This, in combination with
 *     the normal vector, is sufficient to compute the transformation between
 *     surface-tangent and X,Y,Z coordinates.
 *
 * return
 *     TRUE if operation is successful, otherwise FALSE.
 */
Boolean_t FENodeDerivCalc(EntIndex_t     ZoneNum,
						  EntIndex_t     UVarNum,
						  EntIndex_t     VVarNum,
						  EntIndex_t     WVarNum,
						  EntIndex_t     VarNum,
						  SurfElemMap_pa SurfElemMap,
						  Normals_pa     Normals,
						  ArrList_pa     PsiDeriv,
						  ArrList_pa     EtaDeriv,
						  ArrList_pa     PsiAlignment)
{
	Boolean_t    IsOk = TRUE;
	EntIndex_t   NumZones, NumVars;
	LgIndex_t    NumNodes, NumElems, NumNodesPerElem, Node;
	NodeMap_pa   NodeMap = NULL;
	FieldData_pa XVarFDPtr = NULL;
	FieldData_pa YVarFDPtr = NULL;
	FieldData_pa ZVarFDPtr = NULL;
	FieldData_pa VarFDPtr  = NULL;
	FieldData_pa UVarFDPtr = NULL;
	FieldData_pa VVarFDPtr = NULL;
	FieldData_pa WVarFDPtr = NULL;

	ArrList_pa   EdgeOtherNodeList = NULL;
	ArrList_pa   EdgeAngleList  = NULL;
	ArrList_pa   EdgeEndPsiList = NULL;
	ArrList_pa   EdgeEndEtaList = NULL;


	REQUIRE(ZoneNum > 0);
	REQUIRE(VarNum > 0);
	REQUIRE(SurfElemMapIsValid(SurfElemMap));
	REQUIRE(NormalsIsValid(Normals));
	REQUIRE(ArrListIsValid(PsiDeriv));
	REQUIRE(ArrListIsValid(EtaDeriv));
	REQUIRE(ArrListIsValid(PsiAlignment));

	REQUIRE(TecUtilDataSetIsAvailable());

	// Test the SurfaceFit routines
	REQUIRE(SurfaceFitTest());

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	REQUIRE(VarNum  > 0 && VarNum  <= NumVars);

	TecUtilZoneGetInfo(ZoneNum, &NumNodes, &NumElems, &NumNodesPerElem, &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(NumNodes > 1 && NumElems > 1);
	REQUIRE(TecUtilZoneIsFiniteElement(ZoneNum));
	REQUIRE(VALID_REF(XVarFDPtr));
	REQUIRE(VALID_REF(YVarFDPtr));
	REQUIRE(VALID_REF(ZVarFDPtr));

	NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
	VarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, VarNum);

	UVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, UVarNum);
	VVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, VVarNum);
	WVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, WVarNum);
	CHECK(VALID_REF(UVarFDPtr) && VALID_REF(VVarFDPtr) && VALID_REF(WVarFDPtr));

	// Determine the number of nodes per element
	switch (TecUtilZoneGetType(ZoneNum))
	{
		case ZoneType_FETriangle:
			NumNodesPerElem = 3;
			break;
		case ZoneType_FEQuad:
			NumNodesPerElem = 4;
			break;
		default:
			IsOk = FALSE;
			break;
	}


	EdgeOtherNodeList = ArrListAlloc(10, ArrListType_Long);
	EdgeAngleList     = ArrListAlloc(10, ArrListType_Double);
	EdgeEndPsiList    = ArrListAlloc(10, ArrListType_Double);
	EdgeEndEtaList    = ArrListAlloc(10, ArrListType_Double);

	for (Node = 1; Node <= NumNodes; Node++)
	{
		ArrListItem_u Item;
		int           PsiAlign = 0;
		XYZ_s         Normal = NormalsGetNormalForNode(Normals, Node);
		LgIndex_t     NumElemsForNode;
		// LgIndex_t     NumElemsForNode, ElemOffset;

		// Determine PsiAlignment.
		// If normal most closely aligns with Z, PsiAlignment is 1.
		if (ABS(Normal.Z) >= ABS(Normal.X) && ABS(Normal.Z) >= ABS(Normal.Y))
			PsiAlign = 1;
		// If normal most closely aligns with X, PsiAlignment is 2.
		else if (ABS(Normal.X) >= ABS(Normal.Y) && ABS(Normal.X) > ABS(Normal.Z))
			PsiAlign = 2;
		// If normal most closely aligns with Y, PsiAlignment is 3.
		else if (ABS(Normal.Y) > ABS(Normal.X) && ABS(Normal.Y) > ABS(Normal.Z))
			PsiAlign = 3;

		Item.Short = PsiAlign;
		IsOk = ArrListSetItem(PsiAlignment, Node - 1, Item);


		NumElemsForNode = SurfElemMapGetElemCountForNode(SurfElemMap, Node);

		// Find all edges radiating out from the node. Store them in EdgeOtherNodeList.
		if (IsOk)
			IsOk = ComputeEdgeOtherNodeListFromSEM(ZoneNum, Node, SurfElemMap, EdgeOtherNodeList, EdgeAngleList);
		/*
		if (IsOk)
		  {
			for (ElemOffset = 0; ElemOffset < NumElemsForNode; ElemOffset++)
			  {
				LgIndex_t Edge, CornerNodeBeg, CornerNodeEnd;
				LgIndex_t Elem = SurfElemMapGetElemNumForNodeOffset(SurfElemMap, Node, ElemOffset);


				// Find the edges radiating out from the node. Store them by "other" node number
				CornerNodeEnd = TecUtilDataNodeGetByRef(NodeMap, Elem, 1);
				for (Edge = 1; Edge <= NumNodesPerElem; Edge++)
				  {
					LgIndex_t CornerEnd = Edge + 1;
					if (CornerEnd > NumNodesPerElem) CornerEnd -= NumNodesPerElem;

					CornerNodeBeg = CornerNodeEnd;
					CornerNodeEnd = TecUtilDataNodeGetByRef(NodeMap, Elem, CornerEnd);

					if (CornerNodeBeg == Node && CornerNodeEnd != Node)
					   ArrListAppendUniqueLongItem(EdgeOtherNodeList, CornerNodeEnd);

					if (CornerNodeEnd == Node && CornerNodeBeg != Node)
					   ArrListAppendUniqueLongItem(EdgeOtherNodeList, CornerNodeBeg);
				  }
			  }
			CHECK(ArrListGetCount(EdgeOtherNodeList) == NumElemsForNode ||
				  ArrListGetCount(EdgeOtherNodeList) == NumElemsForNode + 1);
		  } */



		// For each element containing the node, project the lines radiating from the node
		// onto the tangent plane
		if (IsOk)
		{
			LgIndex_t     EdgeOffset;
			LgIndex_t     EdgeOtherNode;
			double        x0, y0, z0, dx, dy, dz, Psi, Eta;
			XYZ_s         Normal;
			double        dXDotN, Length, Ratio;
			LgIndex_t     NumEdges = ArrListGetCount(EdgeOtherNodeList);
			ArrListItem_u Item;

			x0 = TecUtilDataValueGetByRef(XVarFDPtr, Node);
			y0 = TecUtilDataValueGetByRef(YVarFDPtr, Node);
			z0 = TecUtilDataValueGetByRef(ZVarFDPtr, Node);

			for (EdgeOffset = 0; EdgeOffset < NumEdges; EdgeOffset++)
			{
				Item = ArrListGetItem(EdgeOtherNodeList, EdgeOffset);
				EdgeOtherNode = Item.Long;

				// Compute DeltaX (X is coordinate vector) along edge line
				dx = TecUtilDataValueGetByRef(XVarFDPtr, EdgeOtherNode) - x0;
				dy = TecUtilDataValueGetByRef(YVarFDPtr, EdgeOtherNode) - y0;
				dz = TecUtilDataValueGetByRef(ZVarFDPtr, EdgeOtherNode) - z0;
				Length = sqrt(dx * dx + dy * dy + dz * dz);

				// Subtract (DeltaX -dot- Normal) to get projection of DeltaX along surface
				Normal = NormalsGetNormalForNode(Normals, Node);
				dXDotN = dx * Normal.X + dy * Normal.Y + dz * Normal.Z;

				dx = dx - dXDotN * Normal.X;
				dy = dy - dXDotN * Normal.Y;
				dz = dz - dXDotN * Normal.Z;

				Ratio = Length / sqrt(dx * dx + dy * dy + dz * dz);


				// Rescale the length of the radiating line to match non-tangent length
				dx = dx * Ratio;
				dy = dy * Ratio;
				dz = dz * Ratio;

				// Convert to surface-tangent coordinates (psi,eta)
				IsOk = XYZtoPsiEta(PsiAlign, Normal, dx, dy, dz, &Psi, &Eta);

				// Save the Psi,Eta coordinates for use in derivative comp
				if (IsOk)
				{
					Item.Double = Psi;
					IsOk = ArrListAppendItem(EdgeEndPsiList,  Item);
				}
				if (IsOk)
				{
					Item.Double = Eta;
					IsOk = ArrListAppendItem(EdgeEndEtaList,  Item);
				}
			}
		}


		// TODO(Finish this!) -- Function needs to be tested




		// Compute the derivatives of the scalar Var using linear least-square surface-fit
		if (IsOk)
		{
			// int           ne, i;
			int           ne;
			ArrListItem_u Item;
			// double      **a = NULL;
			// double      **v = NULL;
			// double       *b = NULL;
			// double       *w = NULL;
			// double       *x = NULL;
			// double        wmin, wmax;
			double        Psi, Eta, Rho;
			double        DRhoDPsi, DRhoDEta;
			LgIndex_t     NumEdges = ArrListGetCount(EdgeOtherNodeList);
			LgIndex_t     EdgeOtherNode;
			SurfaceFit_pa SurfaceFit = SurfaceFitAlloc();;

			// a = SVDAllocMatrix(1, NumEdges+1, 1, 3);
			// v = SVDAllocMatrix(1, 3, 1, 3);
			// b = SVDAllocVector(1, NumEdges+1);
			// w = SVDAllocVector(1, NumEdges+1);
			// x = SVDAllocVector(1, NumEdges+1);

			// Matrix contribution of Node
			// a[1][1] = 1.0;
			// a[1][2] = 0.0;
			// a[1][3] = 0.0;
			// b[1]    = TecUtilDataValueGetByRef(VarFDPtr, Node);
			Rho    = TecUtilDataValueGetByRef(VarFDPtr, Node);
			IsOk = SurfaceFitAppendPointAtEnd(SurfaceFit, 0.0, 0.0, Rho);
			for (ne = 1; IsOk && ne <= NumEdges; ne++)
			{
				Item = ArrListGetItem(EdgeEndPsiList, ne - 1);
				Psi  = Item.Double;
				Item = ArrListGetItem(EdgeEndEtaList, ne - 1);
				Eta  = Item.Double;
				// a[ne+1][1] = 1.0;
				// a[ne+1][2] = Psi;
				// a[ne+1][3] = Eta;

				Item = ArrListGetItem(EdgeOtherNodeList, ne - 1);
				EdgeOtherNode = Item.Long;
				// b[ne+1] = TecUtilDataValueGetByRef(VarFDPtr, EdgeOtherNode);
				Rho = TecUtilDataValueGetByRef(VarFDPtr, EdgeOtherNode);
				IsOk = SurfaceFitAppendPointAtEnd(SurfaceFit, Psi, Eta, Rho);
			}

			// Compute the planar surface fit coefficients
			if (IsOk) IsOk = SurfaceFitCompute(SurfaceFit, SurfFitType_Planar);

			// Compute the singular Value decomposition
			// ComputeSVD(a, NumEdges+1, 3, w, v);

			// If singular values are less than a certain level, set them to zero.
			// wmax = -LARGEFLOAT;
			// for (i = 1; i <= 3; i++) if (wmax < w[i]) wmax = w[i];
			// wmin = 1.0e-6 * wmax;
			// for (i = 1; i <= 3; i++)
			//   {
			//     if (w[i] < wmin)
			//       {
			//         w[i] = 0.0;
			//       }
			//   }

			// Solve the system of equations for the coefficients of the curve fit.
			// BackSubstituteSVD(a, w, v, NumEdges+1, 3, b, x);

			// Set the PsiDeriv and EtaDeriv list items for Node
			DRhoDPsi = SurfaceFitGetCoefFromOffset(SurfaceFit, 1);
			DRhoDEta = SurfaceFitGetCoefFromOffset(SurfaceFit, 2);
			// Item.Double = x[2];
			Item.Double = DRhoDPsi;
			IsOk = ArrListSetItem(PsiDeriv, Node - 1, Item);
			// Item.Double = x[3];
			Item.Double = DRhoDEta;
			if (IsOk) IsOk = ArrListSetItem(EtaDeriv, Node - 1, Item);


			// Set XDeriv, YDeriv, and ZDeriv variables
			if (IsOk)
			{
				double XDeriv, YDeriv, ZDeriv;
				// First transform from psi, eta to X, Y, Z
				// IsOk = PsiEtatoXYZ(PsiAlign, Normal, x[2], x[3], &XDeriv, &YDeriv, &ZDeriv);
				IsOk = PsiEtatoXYZ(PsiAlign, Normal, DRhoDPsi, DRhoDEta, &XDeriv, &YDeriv, &ZDeriv);

				// Set the nodal derivative variable for the zone
				TecUtilDataValueSetByRef(UVarFDPtr, Node, XDeriv);
				TecUtilDataValueSetByRef(VVarFDPtr, Node, YDeriv);
				TecUtilDataValueSetByRef(WVarFDPtr, Node, ZDeriv);
			}

			// Free temporary matrices
			// SVDFreeMatrix(a, 1, NumEdges+1, 1, 3);
			// SVDFreeMatrix(v, 1, 3, 1, 3);
			// SVDFreeVector(b, 1, NumEdges+1);
			// SVDFreeVector(w, 1, NumEdges+1);
			// SVDFreeVector(x, 1, NumEdges+1);

			// Clean up
			SurfaceFitDealloc(&SurfaceFit);
		}

		// Clean up
		ArrListClear(EdgeOtherNodeList);
		ArrListClear(EdgeAngleList);
		ArrListClear(EdgeEndPsiList);
		ArrListClear(EdgeEndEtaList);
	}

	ArrListDealloc(&EdgeOtherNodeList);
	ArrListDealloc(&EdgeAngleList);
	ArrListDealloc(&EdgeEndPsiList);
	ArrListDealloc(&EdgeEndEtaList);


	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}





// Boolean_t ExtractCriticalPoints(const EntIndex_t    ZoneNum,
//                                 const EntIndex_t    UVarNum,
//                                 const EntIndex_t    VVarNum,
//                                 const EntIndex_t    WVarNum,
//                                 const EntIndex_t    ChrgDensVarNum,
//                                 const EntIndex_t    TypeVarNum,
//                                 const Boolean_t     PeriodicBC,
//                                 CritPoints_pa       CritPoints)
Boolean_t ExtractCriticalPoints(const ZoneVarInfo_pa  VolZoneVarInfo,
								const ZoneVarInfo_pa  SurfZoneVarInfo,
								const Boolean_t       IsSurf,
								CritPoints_pa         CritPoints)
{
	Boolean_t  IsOk = TRUE;
	Boolean_t  IsStop = FALSE;
	Boolean_t  ShowStatusBar = TRUE;
	EntIndex_t NumZones, NumVars;
	LgIndex_t  IMax, JMax, KMax;
	LgIndex_t  ii, jj, kk;
	EntIndex_t ZoneNum, UVarNum, VVarNum, WVarNum, ChrgDensVarNum, TypeVarNum;
	Boolean_t  PeriodicBC;

	// LOG
	std::ofstream*	BondLog;
	Boolean_t	LogBondInfo = FALSE;

	REQUIRE(TecUtilDataSetIsAvailable());
	REQUIRE(VALID_BOOLEAN(IsSurf));

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	// Set the appropriate variables
	if (IsSurf)
	{
		REQUIRE(VALID_REF(VolZoneVarInfo));
		REQUIRE(VALID_REF(SurfZoneVarInfo));

		ZoneNum        = SurfZoneVarInfo->ZoneNum;
		UVarNum        = SurfZoneVarInfo->UVarNum;
		VVarNum        = SurfZoneVarInfo->VVarNum;
		WVarNum        = SurfZoneVarInfo->WVarNum;
		ChrgDensVarNum = SurfZoneVarInfo->ChrgDensVarNum;
		TypeVarNum     = SurfZoneVarInfo->TypeVarNum;
		PeriodicBC     = SurfZoneVarInfo->PeriodicBC;
	}
	else
	{
		REQUIRE(VALID_REF(VolZoneVarInfo));
		REQUIRE(SurfZoneVarInfo == NULL);

		ZoneNum        = VolZoneVarInfo->ZoneNum;
		UVarNum        = VolZoneVarInfo->UVarNum;
		VVarNum        = VolZoneVarInfo->VVarNum;
		WVarNum        = VolZoneVarInfo->WVarNum;
		ChrgDensVarNum = VolZoneVarInfo->ChrgDensVarNum;
		TypeVarNum     = VolZoneVarInfo->TypeVarNum;
		PeriodicBC     = VolZoneVarInfo->PeriodicBC;

		// LOG
		BondLog = VolZoneVarInfo->BondLog;
		LogBondInfo = VolZoneVarInfo->LogBondInfo;
	}


	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	// REQUIRE(TecUtilZoneIsOrdered(ZoneNum));
	REQUIRE(UVarNum == 0 || UVarNum > 0 && UVarNum <= NumVars);
	REQUIRE(VVarNum == 0 || VVarNum > 0 && VVarNum <= NumVars);
	REQUIRE(WVarNum == 0 || WVarNum > 0 && WVarNum <= NumVars);
	REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
	REQUIRE(CritPointsIsValid(CritPoints));


	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */


	/* Set-up status bar */
	if (ShowStatusBar && !IsSurf)
	{
		TecUtilStatusSuspend(FALSE);
		TecUtilStatusStartPercentDone("Finding Critical Points", TRUE, TRUE);
		TecUtilStatusSuspend(TRUE);
	}

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
	if (TecUtilZoneIsOrdered(ZoneNum))
	{
		LgIndex_t IBeg = 1;
		LgIndex_t IEnd = IMax;
		LgIndex_t JBeg = 1;
		LgIndex_t JEnd = JMax;
		LgIndex_t KBeg = 1;
		LgIndex_t KEnd = KMax;

		// Set it to 3D
		CritPoints->Dimensions = 3;

		if (PeriodicBC)
		{
			IBeg = 0;
			IEnd = IMax + 1; //IMax + 1;
			JBeg = 0;
			JEnd = JMax + 1; //JMax + 1;
			KBeg = 0;
			KEnd = KMax + 1; //KMax + 1;
		}

		// TIM: Looping over every cell in the system checking for critical point.
		// All the pointOffset stuff is to keep the array of CPs sorted.
		for (kk = KBeg; !IsStop && kk < KEnd; kk++)
		{
			for (jj = JBeg; jj < JEnd; jj++)
			{
				for (ii = IBeg; ii < IEnd; ii++)
				{
					double XCrtPt, YCrtPt, ZCrtPt, ChrgDensCrtPt;
					EntIndex_t TypeCrtPt;
					float PrincDirX, PrincDirY, PrincDirZ;

					if (CriticalPointInCell(ZoneNum, UVarNum, VVarNum, WVarNum, ChrgDensVarNum, TypeVarNum,
											ii, jj, kk, PeriodicBC, &XCrtPt, &YCrtPt, &ZCrtPt,
											&ChrgDensCrtPt, &TypeCrtPt, &PrincDirX, &PrincDirY, &PrincDirZ,VolZoneVarInfo))
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
			if (ShowStatusBar && !IsSurf)
			{
				char PercentDoneText[200];
				int  PercentDone;
				PercentDone = (int)((100 * kk) / KEnd);
				sprintf(PercentDoneText, "Finding Critical Points: %d Percent Done", PercentDone);
				TecUtilStatusSuspend(FALSE);
				TecUtilStatusSetPercentDoneText(PercentDoneText);
				if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
				TecUtilStatusSuspend(TRUE);
				if (IsStop) IsOk = FALSE;
			}
		}
	}
	//	This case is for when ExtractCriticalPoints is called for surface critical points
	else  // Zone is FE
	{
		if (TecUtilZoneGetType(ZoneNum) == ZoneType_FEQuad)
		{
			LgIndex_t ne;
			Boolean_t IsOk = TRUE;
			Normals_pa     Normals = NULL;
			SurfElemMap_pa SurfElemMap = NULL;

			// Array lists to store surface derivative data
			ArrList_pa     PsiDeriv = NULL;
			ArrList_pa     EtaDeriv = NULL;
			// Array list giving coordinate direction Psi is most closely aligned with:  1==X, 2==Y, 3==Z
			ArrList_pa     PsiAlignment = NULL;

			PsiDeriv     = ArrListAlloc(IMax + 1, ArrListType_Double);
			EtaDeriv     = ArrListAlloc(IMax + 1, ArrListType_Double);
			PsiAlignment = ArrListAlloc(IMax + 1, ArrListType_Short);

			// Set it to 2D
			CritPoints->Dimensions = 2;

			// Check that derivatives are to be computed (for now)
			// CHECK(UVarNum == 0 && VVarNum == 0 && WVarNum == 0);

			// Compute the surface element map
			if (IsOk)
			{
				// Compute the element map
				SurfElemMap = SurfElemMapAlloc(IMax);
				IsOk = SurfElemMapCompute(SurfElemMap, ZoneNum);
			}

			// Compute the surface normals
			if (IsOk)
			{
				Normals_pa NormalsOld = NormalsAlloc(IMax);
				double     MaxAngleDiff;

				Normals = NormalsAlloc(IMax);
				IsOk = NormalsComputeUsingSEM(NormalsOld, ZoneNum, SurfElemMap);
				IsOk = NormalsComputeUsingIsoVarGrad(Normals, ZoneNum, VolZoneVarInfo->UVarNum, VolZoneVarInfo->VVarNum, VolZoneVarInfo->WVarNum);
				MaxAngleDiff = NormalsCompare(Normals, NormalsOld);
				NormalsDealloc(&NormalsOld);
			}

			// Compute the derivatives at the nodes
			// TEMP
			/*
			if (IsOk)
			  {
				IsOk = FENodeDerivCalc(ZoneNum, UVarNum, VVarNum, WVarNum, ChrgDensVarNum, SurfElemMap, Normals, PsiDeriv, EtaDeriv, PsiAlignment);
			  }
			  */

			// Loop over elements, testing for the presence of critical points
			for (ne = 1; !IsStop && ne <= JMax; ne++)
			{
				double XCrtPt, YCrtPt, ZCrtPt, ChrgDensCrtPt;
				EntIndex_t TypeCrtPt;
				float PrincDirX, PrincDirY, PrincDirZ;

				// For each element, find the critical points
				if (CriticalPointInQuadCell(ZoneNum, UVarNum, VVarNum, WVarNum, ChrgDensVarNum, TypeVarNum,
											Normals, SurfElemMap, PsiAlignment, ne, &XCrtPt, &YCrtPt, &ZCrtPt,
											&ChrgDensCrtPt, &TypeCrtPt, &PrincDirX, &PrincDirY, &PrincDirZ))
					// if( CriticalPointInQuadCellSF(ZoneNum, UVarNum, VVarNum, WVarNum, ChrgDensVarNum, TypeVarNum,
					//                        Normals, SurfElemMap, PsiAlignment, ne, &XCrtPt, &YCrtPt, &ZCrtPt,
					//                        &ChrgDensCrtPt, &TypeCrtPt, &PrincDirX, &PrincDirY, &PrincDirZ) )
				{

					switch (TypeCrtPt)
					{
						case -2:
						{
							LgIndex_t PointOffset = CritPoints->NumCrtPtsM3;
							IsOk =  CritPointsInsertPoint(CritPoints, PointOffset, XCrtPt, YCrtPt, ZCrtPt,
														  ChrgDensCrtPt, (char)TypeCrtPt, 0.0, 0.0, 0.0);
						}
						break;
						case 0:
						{
							LgIndex_t PointOffset = CritPoints->NumCrtPtsM3 + CritPoints->NumCrtPtsM1;
							IsOk =  CritPointsInsertPoint(CritPoints, PointOffset, XCrtPt, YCrtPt, ZCrtPt,
														  ChrgDensCrtPt, (char)TypeCrtPt, PrincDirX, PrincDirY, PrincDirZ);
						}
						break;
						case 2:
						{
							LgIndex_t PointOffset = CritPoints->NumCrtPtsM3 + CritPoints->NumCrtPtsM1
													+ CritPoints->NumCrtPtsP1 + CritPoints->NumCrtPtsP3;
							IsOk =  CritPointsInsertPoint(CritPoints, PointOffset, XCrtPt, YCrtPt, ZCrtPt,
														  ChrgDensCrtPt, (char)TypeCrtPt, 0.0, 0.0, 0.0);
						}
						break;
					}
				}

				/* Update status bar */
				if (ShowStatusBar && !IsSurf)
				{
					char PercentDoneText[200];
					int  PercentDone;
					PercentDone = (int)((100 * ne) / JMax);
					sprintf(PercentDoneText, "Finding Critical Points: %d Percent Done", PercentDone);
					TecUtilStatusSuspend(FALSE);
					TecUtilStatusSetPercentDoneText(PercentDoneText);
					if (!TecUtilStatusCheckPercentDone(PercentDone)) IsStop = TRUE;
					TecUtilStatusSuspend(TRUE);
					if (IsStop) IsOk = FALSE;
				}
			}

			// Clean Up
			NormalsDealloc(&Normals);
			SurfElemMapDealloc(&SurfElemMap);
			ArrListDealloc(&PsiDeriv);
			ArrListDealloc(&EtaDeriv);
			ArrListDealloc(&PsiAlignment);
		}
		else  // FE but not Quad
			IsOk = FALSE;
	}

	//	Duplicate CP check isn't needed becuase data agitation is used to prevent them
#if 0

	// TIM: ASK The way duplicates are removed arbitrarily assumes one CP to be the correct one,
	// then removes the rest. "Duplicates" are probably so close to the arbitrarily assumed "correct" 
	// point that the combination of points is almost a single point anyways, so if this 
	// duplicate test catches a duplicate CP, the duplicate is no more incorrect than the CP that is kept.
	// Therefore this method is fine, since there's no way to know which CP is corret anyways.

	// Only eleminate duplicate CP tests for 3D ordered zones
	if (TecUtilZoneIsOrdered(ZoneNum))
	{
		EntIndex_t   XVarNum = TecUtilVarGetNumByAssignment('X');
		EntIndex_t   YVarNum = TecUtilVarGetNumByAssignment('Y');
		EntIndex_t   ZVarNum = TecUtilVarGetNumByAssignment('Z');
		FieldData_pa XFDPtr  = TecUtilDataValueGetReadableNativeRef(ZoneNum, XVarNum); 
		FieldData_pa YFDPtr  = TecUtilDataValueGetReadableNativeRef(ZoneNum, YVarNum); 
		FieldData_pa ZFDPtr  = TecUtilDataValueGetReadableNativeRef(ZoneNum, ZVarNum); 

		// Assume constant, non-unit, grid spacing
		double GridSpacing = /*2*/MAX( MAX( ABS(TecUtilDataValueGetByRef(XFDPtr, 2) - TecUtilDataValueGetByRef(XFDPtr, 1)), 
									   ABS(TecUtilDataValueGetByRef(YFDPtr, 2) - TecUtilDataValueGetByRef(YFDPtr, 1)) ),
									   ABS(TecUtilDataValueGetByRef(ZFDPtr, 2) - TecUtilDataValueGetByRef(ZFDPtr, 1))    );

		/* Eliminate duplicate Cage critical points. */
		if (IsOk && CritPoints->NumCrtPtsP3 > 1)
			IsOk = RemDuplicateCritPoints(TRUE, 0.4, GridSpacing, 3, CritPoints);
		if (IsOk && CritPoints->NumCrtPtsP3 > 1)
			IsOk = RemDuplicateCritPoints(FALSE, 1.0, GridSpacing, 3, CritPoints);

		/* Eliminate duplicate Ring critical points. */
		if (IsOk && CritPoints->NumCrtPtsP1 > 1)
			IsOk = RemDuplicateCritPoints(FALSE, 1.0, GridSpacing, 1, CritPoints);

		/* Eliminate duplicate Bond critical points. */
		if (IsOk && CritPoints->NumCrtPtsM1 > 1)
			IsOk = RemDuplicateCritPoints(FALSE, 1.0, GridSpacing, -1, CritPoints);

		/* Eliminate duplicate Atom critical points. */
		if (IsOk && CritPoints->NumCrtPtsM3 > 1)
			IsOk = RemDuplicateCritPoints(TRUE, 0.4, GridSpacing, -3, CritPoints);
		if (IsOk && CritPoints->NumCrtPtsM3 > 1)
			IsOk = RemDuplicateCritPoints(FALSE, 1.0, GridSpacing, -3, CritPoints);
	}
	else
	{
		// TODO: Is this OK?
		double GridSpacing = 1.0;

		// Eliminate duplicate Max critical points
		if (IsOk && CritPoints->NumCrtPtsM3 > 1)
			IsOk = RemDuplicateCritPoints(FALSE, 0.1, GridSpacing, -2, CritPoints);
		// Eliminate duplicate saddle critical points
		if (IsOk && CritPoints->NumCrtPtsM1 > 1)
			IsOk = RemDuplicateCritPoints(FALSE, 0.1, GridSpacing, 0, CritPoints);
		// Eliminate duplicate Min critical points
		if (IsOk && CritPoints->NumCrtPtsP3 > 1)
			IsOk = RemDuplicateCritPoints(FALSE, 0.1, GridSpacing, 2, CritPoints);
	}
#endif

	if (ShowStatusBar && !IsSurf) 
	{
		TecUtilStatusSuspend(FALSE);
		TecUtilStatusFinishPercentDone();
		TecUtilStatusSuspend(TRUE);
	}

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
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

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
		char              ZoneName[200];
		EntIndex_t        SourceZoneNum = CritPoints->SourceZoneNum;
		LgIndex_t         iv;

		sprintf(ZoneName, "Critical Points Zone %d", SourceZoneNum);

		// AveZoneNum = TUZoneGetNumByName(ZoneName);

		/* Set FieldDataType of all variables of CP zone */
		FieldDataType_e  *VarDataType = ALLOC_ARRAY(NumVars, FieldDataType_e, "VarDataType");
		for (iv = 0; iv < NumVars; iv++)
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
		}
		if ( VarDataType != NULL )
			FREE_ARRAY(VarDataType, "VarDataType");

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
				char NumCrtPtStr[200];
				char BaseZoneNumStr[200];

				IsOk = TecUtilAuxDataSetStrItem(ZnAuxDataRef, "CompChem.ZoneType", "CriticalPoints", TRUE);

				sprintf(NumCrtPtStr, "%d", CritPoints->NumCrtPtsM3);
				IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.NumCrtPtAtom", (ArbParam_t)NumCrtPtStr,
											 AuxDataType_String, TRUE);

				sprintf(NumCrtPtStr, "%d", CritPoints->NumCrtPtsM1);
				if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.NumCrtPtBond",
														   (ArbParam_t)NumCrtPtStr, AuxDataType_String, TRUE);

				sprintf(NumCrtPtStr, "%d", CritPoints->NumCrtPtsP1);
				if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.NumCrtPtRing",
														   (ArbParam_t)NumCrtPtStr, AuxDataType_String, TRUE);

				sprintf(NumCrtPtStr, "%d", CritPoints->NumCrtPtsP3);
				if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.NumCrtPtCage",
														   (ArbParam_t)NumCrtPtStr, AuxDataType_String, TRUE);

				sprintf(BaseZoneNumStr, "%d", SourceZoneNum);
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
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

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

		/*
		if (Result != NULL)
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








/*
 * Verify the unit test of LinearTriangleNaturalCoord for analytical distribution
 * and a simple triangle.
 *
 *                   (0,1/2)
 *                     /\ Node2
 *                   /    \
 *                 /        \
 *               /            \
 *             /                \
 *     Node0 /____________________\  Node1
 *     (-1/2,-1/2)                   (1/2,-1/2)
 *
 * return
 *     TRUE if successful, FALSE  if it fails.
 */
Boolean_t LinearTriangleNaturalCoordTest()
{
	Boolean_t  IsOk = TRUE;
	double     U, V, UCMin, UCMax, VCMin, VCMax;
	double     UCell[3], VCell[3];
	double     XCell[3], YCell[3];
	double     X, Y;
	double     Weight[3];
	double     SmallDiff = 1.0e-6;

	// Set up triangle
	XCell[0] = -0.5;
	YCell[0] = -0.5;

	XCell[1] =  0.5;
	YCell[1] = -0.5;

	XCell[2] =  0.0;
	YCell[2] =  0.5;

	// Velocities - assume  Rho = X**2 - Y**2;
	UCell[0] = -1.0;
	VCell[0] =  1.0;

	UCell[1] =  1.0;
	VCell[1] =  1.0;

	UCell[2] =  0.0;
	VCell[2] = -1.0;

	// Set min,max
	UCMin = -1.0;
	UCMax =  1.0;
	VCMin = -1.0;
	VCMax =  1.0;

	// Look for U,V=0,0
	U = 0.0;
	V = 0.0;


	IsOk = LinearTriangleNaturalCoord(U, V, UCell, VCell, UCMin, VCMin, UCMax, VCMax, XCell, YCell, &X, &Y, Weight);


	// Test the results: X, Y should be (0.0, 0.0)
	if (IsOk) IsOk = (ABS(X) < SmallDiff);
	if (IsOk) IsOk = (ABS(Y) < SmallDiff);

	if (IsOk)
	{
		double D2RhoDx2, D2RhoDxDy, D2RhoDy2;

		IsOk = ComputeFETriangleCurvature(XCell, YCell, UCell, VCell,
										  &D2RhoDx2, &D2RhoDxDy, &D2RhoDy2);

		if (IsOk) IsOk = (ABS(D2RhoDxDy) < SmallDiff);
		if (IsOk) IsOk = (D2RhoDx2 >  1.99999 && D2RhoDx2 <  2.00001);
		if (IsOk) IsOk = (D2RhoDy2 > -2.00001 && D2RhoDy2 < -1.99999);
	}


	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}








/*
 * Verify the unit test of BilinearQuadNaturalCoord for analytical distribution
 * and a simple quad.
 *
 *            (-1/2,1/2)______(1/2,1/2)
 *              Node3  /      \ Node2
 *                   /          \
 *                 /              \
 *               /                  \
 *             /                      \
 *     Node0 /__________________________\  Node1
 *     (-1,-1/2)                     (1,-1/2)
 *
 * return
 *     TRUE if successful, FALSE  if it fails.
 */
Boolean_t BilinearQuadNaturalCoordTest()
{
	Boolean_t  IsOk = TRUE;
	double     U, V, UCMin, UCMax, VCMin, VCMax;
	double     UCell[4], VCell[4];
	double     XCell[4], YCell[4];
	double     r, s;
	double     Weight[4];
	double     SmallDiff = 2.0e-4;

	// Set up quad
	XCell[0] = -1.0;
	YCell[0] = -0.5;

	XCell[1] =  1.0;
	YCell[1] = -0.5;

	XCell[2] =  0.5;
	YCell[2] =  0.5;

	XCell[3] = -0.5;
	YCell[3] =  0.5;

	// Velocities - assume  Rho = X**2 - Y**2;
	UCell[0] = -2.0;
	VCell[0] =  1.0;

	UCell[1] =  2.0;
	VCell[1] =  1.0;

	UCell[2] =  1.0;
	VCell[2] = -1.0;

	UCell[3] = -1.0;
	VCell[3] = -1.0;

	// Set min,max
	UCMin = -2.0;
	UCMax =  2.0;
	VCMin = -1.0;
	VCMax =  1.0;

	// Look for U,V=0,0
	U = 0.0;
	V = 0.0;

	r = 0.5;
	s = 0.5;
	IsOk = BilinearQuadNaturalCoord(U, V, UCell, VCell, UCMin, VCMin,
									UCMax, VCMax, &r, &s, Weight);

	// IsOk = LinearTriangleNaturalCoord(U, V, UCell, VCell, UCMin, VCMin,
	//                                   UCMax, VCMax, XCell, YCell, &X, &Y, Weight);


	// Test the results: X, Y should be (0.0, 0.0)
	if (IsOk) IsOk = (ABS(r) < SmallDiff);
	if (IsOk) IsOk = (ABS(s) < SmallDiff);

	if (IsOk)
	{
		double D2RhoDx2, D2RhoDxDy, D2RhoDy2;

		IsOk = ComputeFEQuadCurvature(XCell, YCell, UCell, VCell, r, s,
									  &D2RhoDx2, &D2RhoDxDy, &D2RhoDy2);

		if (IsOk) IsOk = (ABS(D2RhoDxDy) < SmallDiff);
		if (IsOk) IsOk = (D2RhoDx2 >  1.99999 && D2RhoDx2 <  2.00001);
		if (IsOk) IsOk = (D2RhoDy2 > -2.00001 && D2RhoDy2 < -1.99999);
	}


	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}


/*
 *  Find volume CP corresponding to IsoTopo surface CP. We know the volume
 *  atom, so we just need to search the volume VolCP-atom connector lines for
 *  one that passes within some tolerance of the SurfCP.
 *
 * param VolCPType
 *     Type of the volume CP being sought
 * param CritPoints
 *     Volume critical point structure
 * param SurfCritPoints
 *     IsoTopo surface critical points structure
 * param Atom
 *     Atom (volume) number (zero-based offset)
 * param SurfCPType
 *     Type of the surface CP being tested against
 * param SurfCPOffset
 *     SurfCP (IsoTopo surface) number (zero-based offset)
 * param Tolerance
 *     Maximum Bond-GradPath distance that is considered a "hit"
 *
 * return
 *     VolCP (volume, zero-based offset) number whose VolCP-Atom connector passes
 *     within Tolerance distance of the SurfCP.
 *     -1 if no VolCP-Atom line passed within Tolerance of the SurfCP (assumed FF)
 *     -2 if multiple VolCP-Atom lines passed within Tolerance of the SurfCP
 */
LgIndex_t VolCPFromAtomIsoTopoSurfCP(const char          VolCPType,
									 const CritPoints_pa CritPoints,
									 const CritPoints_pa SurfCritPoints,
									 const LgIndex_t     AtomNum,
									 const char          SurfCPType,
									 const LgIndex_t     SurfCPOffset,
									 const double        Tolerance)
{
	LgIndex_t Result = -1;

	Boolean_t IsOk = TRUE;
	LgIndex_t NumAtoms, NumVolCPs;
	LgIndex_t NumSurfCPs;
	LgIndex_t AtomCPNum, SurfCPNum;
	XYZ_s     SurfCPXYZ;

	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CritPointsIsValid(SurfCritPoints));
	REQUIRE(CritPoints->Dimensions == 3);
	REQUIRE(SurfCritPoints->Dimensions == 2);

	NumAtoms = CritPointsGetEndOffset(CritPoints, (char)(-3))
			   - CritPointsGetBegOffset(CritPoints, (char)(-3));
	NumVolCPs = CritPointsGetEndOffset(CritPoints, VolCPType)
				- CritPointsGetBegOffset(CritPoints, VolCPType);
	NumSurfCPs = CritPointsGetEndOffset(SurfCritPoints, SurfCPType)
				 - CritPointsGetBegOffset(SurfCritPoints, SurfCPType);

	REQUIRE(AtomNum >= 0 && AtomNum < NumAtoms);
	REQUIRE(SurfCPOffset >= 0 && SurfCPOffset < NumSurfCPs);
	REQUIRE(Tolerance > SMALLFLOAT);

	if (NumVolCPs > 0)
	{
		LgIndex_t  VolCP;
		double     XCP, YCP, ZCP, Rho, Px, Py, Pz;
		char       Type;
		EntIndex_t VolSourceZoneNum = CritPoints->SourceZoneNum;

		AtomCPNum = CritPointsGetBegOffset(CritPoints, (char)(-3)) + AtomNum;

		// Get the location of the IsoTopo surface Min CP
		SurfCPNum = CritPointsGetBegOffset(SurfCritPoints, SurfCPType) + SurfCPOffset;

		IsOk = CritPointsGetPoint(SurfCritPoints, SurfCPNum, &XCP, &YCP, &ZCP, &Rho,
								  &Type, &Px, &Py, &Pz);

		// Verify that this CP is a min
		if (IsOk && Type != SurfCPType) IsOk = FALSE;

		// Set the location in the SurfCPXYZ structure
		SurfCPXYZ.X = XCP;
		SurfCPXYZ.Y = YCP;
		SurfCPXYZ.Z = ZCP;

		// Cycle through Bonds, looking for a VolCPType-AtomNum GradPath that passes
		// within Tolerance distance of MinCPNum
		for (VolCP = 0; IsOk && VolCP < NumVolCPs; VolCP++)
		{
			EntIndex_t TPZoneOfGP;

			LgIndex_t  VolCPNum = CritPointsGetBegOffset(CritPoints, VolCPType) + VolCP;

			TPZoneOfGP = GradPathTPZoneFromBegEndCP(VolCPNum, AtomCPNum, VolSourceZoneNum);

			// Bond includes specified AtomNum if TPZone is found
			if (TPZoneOfGP > 0)
			{
				LgIndex_t I, J, K;
				TecUtilZoneGetIJK(TPZoneOfGP, &I, &J, &K);

				if (I > 1 && J == 1 && K == 1){
					double      DistPathToMinCP = -1.0;
					GradPath_pa GradPath = GradPathGetFromTPZone(TPZoneOfGP);


					DistPathToMinCP = GradPathGetDistToPt(GradPath, SurfCPXYZ);

					// Found!!
					if (DistPathToMinCP >= 0.0 && DistPathToMinCP < Tolerance)
					{
						// First Bond found
						if (Result == -1)
						{
							Result = VolCP;
						}
						else if (Result >= 0)  // Second VolCP found
						{
							Result = -2;
						}
						else  // more than 2 VolCPs found
						{
							Result = -2;
						}
					}

					GradPathDealloc(&GradPath);
				}
			}
		}
	}

	if (IsOk == FALSE) Result = -1;

	ENSURE(Result >= -2 && Result < NumVolCPs);
	return Result;
}

									 /*
 *  Find volume CP corresponding to IsoTopo surface CP. We know the volume
 *  atom, so we just need to search the volume VolCP-atom connector lines for
 *  one that passes within some tolerance of the SurfCP.
 *	This version of the function returns the volume CP of the closest volume
 *	GP, rather than returning a special value if more than one is found.
 *
 * param VolCPType
 *     Type of the volume CP being sought
 * param CritPoints
 *     Volume critical point structure
 * param SurfCritPoints
 *     IsoTopo surface critical points structure
 * param Atom
 *     Atom (volume) number (zero-based offset)
 * param SurfCPType
 *     Type of the surface CP being tested against
 * param SurfCPOffset
 *     SurfCP (IsoTopo surface) number (zero-based offset)
 * param Tolerance
 *     Maximum Bond-GradPath distance that is considered a "hit"
 *
 * return
 *     VolCP (volume, zero-based offset) number whose VolCP-Atom connector passes
 *     within Tolerance distance of the SurfCP.
 *     -1 if no VolCP-Atom line passed within Tolerance of the SurfCP (assumed FF)
 *     -2 if multiple VolCP-Atom lines passed within Tolerance of the SurfCP
 */
LgIndex_t VolCPFromAtomIsoTopoSurfCPClosest(const char          VolCPType,
											 const CritPoints_pa CritPoints,
											 const CritPoints_pa SurfCritPoints,
											 const LgIndex_t     AtomNum,
											 const char          SurfCPType,
											 const LgIndex_t     SurfCPOffset,
											 const double        Tolerance)
{
	LgIndex_t Result = -1;

	Boolean_t IsOk = TRUE;
	LgIndex_t NumAtoms, NumVolCPs;
	LgIndex_t NumSurfCPs;
	LgIndex_t AtomCPNum, SurfCPNum;
	XYZ_s     SurfCPXYZ;

	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CritPointsIsValid(SurfCritPoints));
	REQUIRE(CritPoints->Dimensions == 3);
	REQUIRE(SurfCritPoints->Dimensions == 2);

	NumAtoms = CritPointsGetEndOffset(CritPoints, (char)(-3))
			   - CritPointsGetBegOffset(CritPoints, (char)(-3));
	NumVolCPs = CritPointsGetEndOffset(CritPoints, VolCPType)
				- CritPointsGetBegOffset(CritPoints, VolCPType);
	NumSurfCPs = CritPointsGetEndOffset(SurfCritPoints, SurfCPType)
				 - CritPointsGetBegOffset(SurfCritPoints, SurfCPType);

	REQUIRE(AtomNum >= 0 && AtomNum < NumAtoms);
	REQUIRE(SurfCPOffset >= 0 && SurfCPOffset < NumSurfCPs);
	REQUIRE(Tolerance > SMALLFLOAT);

	if (NumVolCPs > 0)
	{
		LgIndex_t  VolCP = -1;
		double     XCP, YCP, ZCP, Rho, Px, Py, Pz;
		char       Type;
		EntIndex_t VolSourceZoneNum = CritPoints->SourceZoneNum;

		AtomCPNum = CritPointsGetBegOffset(CritPoints, (char)(-3)) + AtomNum;

		// Get the location of the IsoTopo surface CP
		SurfCPNum = CritPointsGetBegOffset(SurfCritPoints, SurfCPType) + SurfCPOffset;

		IsOk = CritPointsGetPoint(SurfCritPoints, SurfCPNum, &XCP, &YCP, &ZCP, &Rho,
								  &Type, &Px, &Py, &Pz);

		// Verify that this CP is of specified type
		if (IsOk && Type != SurfCPType) IsOk = FALSE;

		// Set the location in the SurfCPXYZ structure
		SurfCPXYZ.X = XCP;
		SurfCPXYZ.Y = YCP;
		SurfCPXYZ.Z = ZCP;

		// Cycle through Volume GPs, looking for a VolCPType-AtomNum GradPath that passes
		// within Tolerance distance of SurfCPNum
		double DistPathToMinCP2 = -1.0;

		for (VolCP = 0; IsOk && VolCP < NumVolCPs; VolCP++)
		{
			EntIndex_t TPZoneOfGP;

			LgIndex_t  VolCPNum = CritPointsGetBegOffset(CritPoints, VolCPType) + VolCP;

			TPZoneOfGP = GradPathTPZoneFromBegEndCP(VolCPNum, AtomCPNum, VolSourceZoneNum);

			// Volume GP includes specified AtomNum if TPZone is found
			if (TPZoneOfGP > 0)
			{
				double      DistPathToMinCP = -1;
				GradPath_pa GradPath = GradPathGetFromTPZone(TPZoneOfGP);


				DistPathToMinCP = GradPathGetDistToPt(GradPath, SurfCPXYZ);

				// Found!!
				if (DistPathToMinCP >= 0.0 && DistPathToMinCP < Tolerance)
				{
					if (DistPathToMinCP2 == -1)
					{
						DistPathToMinCP2 = DistPathToMinCP;
						Result = VolCP;
					}
					else
					{
						if (DistPathToMinCP < DistPathToMinCP2)
						{
							DistPathToMinCP2 = DistPathToMinCP;
							Result = VolCP;
						}
					}
				}

				GradPathDealloc(&GradPath);
			}
		}
	}

	if (IsOk == FALSE) Result = -1;


	ENSURE(Result >= -1 && Result < NumVolCPs);
	return Result;
}




/*
 * Verify the unit test CritPoints for a regular cube FEQuad zone (cube with ?)
 *
 * return
 *     TRUE if successful, FALSE  if it fails.
 */
Boolean_t CritPointsTest()
{
	Boolean_t      IsOk = TRUE;

	// Volume Zone



	// Surface Zone

	// First simple triangle
	IsOk = LinearTriangleNaturalCoordTest();

	// Simple quad
	IsOk = BilinearQuadNaturalCoordTest();

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}

/** TIM
 * Take a given surface grad path, if its beginning and/or end
 * point(s) (surf CPs) are sufficiently close to the point of 
 * isosurface-volume grad paths intersection then replace it with 
 * the intersection point. If the point is far from the intersection
 * then append that point to the surface grad path
 *
 * NOTE:	Currently assumes that a volume grad path exists for each
 *			surf CP at ends of supplied surf grad path.
 *
 * param SurfGradPath
 *		Surface grad path to have its ends adjusted.
 * param SurfCritPoints
 *		Isosurface critical points
 * param VolCritPoints
 *		Volume critical points
 * param AtomNum
 *		Index of atom CP isosurface is based on (do I need this, or
 *		can I take it from the SurfCritPoints->SourceZone property?)
 * param GPFoundTolorance
 *		Tolerance to determine if a volume grad path corresponds
 *		to a surface CP
 * param ReplaceTolorance
 *		Tolerance to determine whether the surface grad path's end
 *		point should be replaced with the intersection point or if
 *		the intersection point should be appended to the surf grad
 *		path
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t	SurfCritPtsMoveToVolGPIsoSurfIntersection(
				 const CritPoints_pa	SurfCritPoints,
				 const CritPoints_pa	VolCritPoints,
				 const LgIndex_t		AtomNum,
				 const double			GPFoundTolerance)
{
	Boolean_t	IsOk = TRUE;

	REQUIRE(CritPointsIsValid(SurfCritPoints));
	REQUIRE(CritPointsIsValid(VolCritPoints));
	REQUIRE(SurfCritPoints->Dimensions == 2);
	REQUIRE(VolCritPoints->Dimensions == 3);
	REQUIRE(AtomNum >= 0);
	REQUIRE(GPFoundTolerance > 0);

	const char	MinCP		= 2,
				MaxCP		= -2,
				SaddleCP	= 0;

	char		VolCPType,
				VolCPType2;

	const char	BondCP		= -1,
				RingCP		= 1,
				RingFFCP	= 11,
				CageCP		= 3,
				CageFFCP	= 13;

	XYZ_s		SurfCPXYZ;
	double		XCP, YCP, ZCP, Rho, Px, Py, Pz;
	char		SurfCPType = 1;

	//	Need atom offset
	LgIndex_t AtomCPNum = CritPointsGetBegOffset(VolCritPoints, (char)(-3)) + AtomNum;

	LgIndex_t NumOfSurfCPs = CritPointsGetCount(SurfCritPoints);
	IsOk = (NumOfSurfCPs > 0);

	//	Loop over all surf crit points 
	Boolean_t Skip = FALSE;
	for (int SurfCPOffset = 0 ; IsOk && !Skip && SurfCPOffset < NumOfSurfCPs ; SurfCPOffset++)
	{
		IsOk = CritPointsGetPoint(SurfCritPoints, SurfCPOffset, &XCP, &YCP, &ZCP, &Rho,
			&SurfCPType, &Px, &Py, &Pz);

		//if (SurfCPType == MaxCP || SurfCPType == MinCP) Skip = TRUE;
		
		if (!Skip)
		{
			if (IsOk)
			{
				SurfCPXYZ.X = XCP;
				SurfCPXYZ.Y = YCP;
				SurfCPXYZ.Z = ZCP;
			}

			if (IsOk)
			{
				//	Need surf CP number minus offset of its type
				LgIndex_t SurfCPNum = SurfCPOffset - CritPointsGetBegOffset(SurfCritPoints, SurfCPType);
				CHECK(SurfCPNum >= 0);

				//	Surface maxes and saddles could go to near or far field vol CPs
				//	so might need to check for both
				VolCPType	= 0;
				VolCPType2	= 0;

				if (SurfCPType == MinCP)
					VolCPType = BondCP;
				else if (SurfCPType == MaxCP)
				{
					VolCPType = CageCP;
					VolCPType2 = CageFFCP;
				}
				else if (SurfCPType == SaddleCP)
				{
					VolCPType = RingCP;
					VolCPType2 = RingFFCP;
				}
				else
					IsOk = FALSE;

				LgIndex_t VolCP;
				/*const double	Tolerance = 0.1;*/

				VolCP = VolCPFromAtomIsoTopoSurfCPClosest(VolCPType,
														VolCritPoints,
														SurfCritPoints,
														AtomNum,
														SurfCPType,
														SurfCPNum,
														GPFoundTolerance);

				Boolean_t	FFType = FALSE;
				if (VolCPType2 != 0 && VolCP < 0)
				{
					VolCP = VolCPFromAtomIsoTopoSurfCPClosest(VolCPType2,
															VolCritPoints,
															SurfCritPoints,
															AtomNum,
															SurfCPType,
															SurfCPNum,
															GPFoundTolerance);
					if (VolCP >= 0) FFType = TRUE;
				}
				if (VolCP < 0) IsOk = FALSE;

				if (IsOk)
				{
					//	Get the volume grad path, check the distance to the surf CP and replace
					LgIndex_t	VolCPNum;
					if (!FFType) VolCPNum = CritPointsGetBegOffset(VolCritPoints, VolCPType) + VolCP;
					else VolCPNum = CritPointsGetBegOffset(VolCritPoints, VolCPType2) + VolCP;

					/*EntIndex_t	VolGPZoneNum = GradPathTPZoneFromBegEndCP(VolCPNum, 
																		AtomCPNum, 
																		VolCritPoints->SourceZoneNum);*/
					GradPath_pa	VolGradPath = GradPathGetByBegEndCP(VolCPNum, AtomCPNum);
					ENSURE(GradPathIsValid(VolGradPath));
					IsOk = GradPathIsValid(VolGradPath);

					if (IsOk)
					{
						//double		DistPathToSurfCP = -1.0;

						//	Now need to find the intersection point of the volume grad
						//	path and the isosurface
						XYZ_s		LnBeg, LnEnd, Pt;
						double		Rho1, Rho2;
						LgIndex_t	NumOfVolGPPoints = GradPathGetCount(VolGradPath);
						if (NumOfVolGPPoints > 0)
						{
							double	RSquare1, RSquare2, Min1, Min2;
							Boolean_t IsFound = FALSE;
							for (LgIndex_t j = NumOfVolGPPoints - 1 ; j > 0 && !IsFound ; j--)
							{
								GradPathGetPoint(VolGradPath, j, 
													&LnBeg.X, 
													&LnBeg.Y, 
													&LnBeg.Z, 
													&Rho1);
								GradPathGetPoint(VolGradPath, j-1, 
													&LnEnd.X, 
													&LnEnd.Y, 
													&LnEnd.Z, 
													&Rho2);
								RSquare1 = DistanceSquaredXYZ(LnBeg, SurfCPXYZ);
								RSquare2 = DistanceSquaredXYZ(LnEnd, SurfCPXYZ);
								if (j == NumOfVolGPPoints)
								{
									Min1 = RSquare1;
									Min2 = RSquare2;
								}
								else
								{
									if (RSquare2 > Min2 && RSquare1 < Min1)
									{
										//	Interpolate to find intersection point.
										double DistanceRatio = RSquare1 / (RSquare1 + RSquare2);
										Pt.X = LnBeg.X + DistanceRatio * (LnEnd.X - LnBeg.X);
										Pt.Y = LnBeg.Y + DistanceRatio * (LnEnd.Y - LnBeg.Y);
										Pt.Z = LnBeg.Z + DistanceRatio * (LnEnd.Z - LnBeg.Z);

										//	Quadratic interpolation to find Rho at the point.
										double	Rho0, b1, b2;
										XYZ_s	Ln0;
										GradPathGetPoint(VolGradPath, j+1,
															&Ln0.X,
															&Ln0.Y,
															&Ln0.Z,
															&Rho0);
										b1 = (Rho1 - Rho0) / sqrt(DistanceSquaredXYZ(Ln0, LnBeg));
										b2 = ((Rho2 - Rho1) / sqrt(DistanceSquaredXYZ(LnBeg, LnEnd)) - b1) 
											/ sqrt(DistanceSquaredXYZ(LnEnd, Ln0));
										double TempDist = sqrt(DistanceSquaredXYZ(Ln0, Pt));
										Rho = Rho0 + b1 * TempDist + b2 * TempDist
												* sqrt(DistanceSquaredXYZ(LnBeg, Pt));

										double SurfX, SurfY, SurfZ;
										LgIndex_t SurfElem;

										SurfElem = GeomToolsClosestElemNum(SurfCritPoints->SourceZoneNum, 
																			Pt.X, Pt.Y, Pt.Z,
																			&SurfX, &SurfY, &SurfZ);
										if (SurfElem > 0)
										{
											Pt.X = SurfX;
											Pt.Y = SurfY;
											Pt.Z = SurfZ;
										}
										else IsOk = FALSE;

										IsFound = TRUE;
									}
									else
									{
										//	Still approaching surface, so reset Min variables
										//	and continue
										Min1 = RSquare1;
										Min2 = RSquare2;
									}
								}
							}
							if (!IsFound) IsOk = FALSE;
						}
						else IsOk = FALSE;

						//	Now replace the surf CP with the intersection point just found.

						if (IsOk)
						{
							IsOk = CritPointsRemovePoint(SurfCritPoints, SurfCPOffset);
							if (IsOk) IsOk = CritPointsInsertPoint(SurfCritPoints, SurfCPOffset, 
													Pt.X, Pt.Y, Pt.Z, Rho, SurfCPType, Px, Py, Pz);
						}
					}

					GradPathDealloc(&VolGradPath);

				}
			}
			//if (!IsOk) break;
		}
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}	//	Boolean_t	SurfCritPtsMoveToVolGPIsoSurfIntersection

Boolean_t CritPointGetTypeString(const LgIndex_t CritPointNum,
								 const char      CritPointType,
								 char           *CritPointTypeString)
{
	Boolean_t IsOk = TRUE;
	int Size = 0;

	REQUIRE(CritPointNum >= -1);
	// REQUIRE(CritPointNum == -1 ||
	//        (CritPointType == -3 || CritPointType == -1 || CritPointType == 1 || CritPointType == 3) );
	REQUIRE(VALID_REF(CritPointTypeString));

	if (CritPointNum < 0)
	{
		Size = sprintf(CritPointTypeString, "FarField");
		if (Size != 8) IsOk = FALSE;
	}
	else
	{
		switch (CritPointType)
		{
		case -3:
			Size = sprintf(CritPointTypeString, "Atom");
			break;
		case -2:
			Size = sprintf(CritPointTypeString, "Max");
			break;
		case -1:
			Size = sprintf(CritPointTypeString, "Bond");
			break;
		case 0:
			Size = sprintf(CritPointTypeString, "Saddle");
			break;
		case 1:
			Size = sprintf(CritPointTypeString, "Ring");
			break;
		case 2:
			Size = sprintf(CritPointTypeString, "Min");
			break;
		case 3:
			Size = sprintf(CritPointTypeString, "Cage");
			break;
		case 11:
			Size = sprintf(CritPointTypeString, "RingFF");
			break;
		case 13:
			Size = sprintf(CritPointTypeString, "CageFF");
			break;
		}
		if (!(Size == 4 || Size == 3 || Size == 6)) IsOk = FALSE;
	}
	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}

/*
 *  Find volume bond corresponding to IsoTopo surface Min. We know the volume
 *  atom, so we just need to search the volume bond-atom connector lines for
 *  one that passes within some tolerance of the Min CP.
 *
 * param CritPoints
 *     Volume critical point structure
 * param SurfCritPoints
 *     IsoTopo surface critical points structure
 * param Atom
 *     Atom (volume) number (zero-based offset)
 * param Min
 *     Min (IsoTopo surface) number (zero-based offset)
 * param Tolerance
 *     Maximum Bond-GradPath distance that is considered a "hit"
 *
 * return
 *     Bond (volume, zero-based offset) number that is within Tolerance distance
 *     of the MinCP.
 *     -1 if no Bond-Atom line passed within Tolerance of the MinCP (assumed FF)
 *     -2 if multiple Bond-Atom lines passed within Tolerance of the MinCP
 */
LgIndex_t BondCPFromAtomIsoTopoMin(CritPoints_pa CritPoints,
								   CritPoints_pa SurfCritPoints,
								   LgIndex_t     AtomNum,
								   LgIndex_t     SurfMinNum,
								   double        Tolerance)
{
	LgIndex_t BondNum = -1;

	Boolean_t IsOk = TRUE;
	LgIndex_t NumAtoms, NumBonds;
	LgIndex_t NumSurfMins;
	LgIndex_t AtomCPNum, BondCPNum, MinCPNum;
	XYZ_s     MinCPXYZ;

	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CritPointsIsValid(SurfCritPoints));
	REQUIRE(CritPoints->Dimensions == 3);
	REQUIRE(SurfCritPoints->Dimensions == 2);

	NumAtoms = CritPointsGetEndOffset(CritPoints, (char)(-3))
			   - CritPointsGetBegOffset(CritPoints, (char)(-3));
	NumBonds = CritPointsGetEndOffset(CritPoints, (char)(-1))
			   - CritPointsGetBegOffset(CritPoints, (char)(-1));
	NumSurfMins = CritPointsGetEndOffset(SurfCritPoints, (char)(2))
				  - CritPointsGetBegOffset(SurfCritPoints, (char)(2));

	REQUIRE(AtomNum >= 0 && AtomNum < NumAtoms);
	REQUIRE(SurfMinNum >= 0 && SurfMinNum < NumSurfMins);
	REQUIRE(Tolerance > SMALLFLOAT);

	if (NumBonds > 0)
	{
		LgIndex_t  Bond;
		double     XCP, YCP, ZCP, Rho, Px, Py, Pz;
		char       Type;
		EntIndex_t SourceZoneNum = CritPoints->SourceZoneNum;

		AtomCPNum = CritPointsGetBegOffset(CritPoints, (char)(-3)) + AtomNum;

		// Get the location of the IsoTopo surface Min CP
		MinCPNum = CritPointsGetBegOffset(SurfCritPoints, (char)(2)) + SurfMinNum;

		IsOk = CritPointsGetPoint(SurfCritPoints, MinCPNum, &XCP, &YCP, &ZCP, &Rho,
								  &Type, &Px, &Py, &Pz);

		// Verify that this CP is a min
		if (IsOk && Type != (char)(2)) IsOk = FALSE;

		// Set the location in the MinCPXYZ structure
		MinCPXYZ.X = XCP;
		MinCPXYZ.Y = YCP;
		MinCPXYZ.Z = ZCP;

		// Cycle through Bonds, looking for a Bond-AtomNum GradPath that passes
		// within Tolerance distance of MinCPNum
		for (Bond = 0; IsOk && Bond < NumBonds; Bond++)
		{
			EntIndex_t TPZoneOfGP;

			BondCPNum = CritPointsGetBegOffset(CritPoints, (char)(-1)) + Bond;

			TPZoneOfGP = GradPathTPZoneFromBegEndCP(BondCPNum, AtomCPNum, SourceZoneNum);

			// Bond includes specified AtomNum if TPZone is found
			if (TPZoneOfGP > 0)
			{
				double      DistPathToMinCP = -1.0;
				GradPath_pa GradPath = GradPathGetFromTPZone(TPZoneOfGP);


				DistPathToMinCP = GradPathGetDistToPt(GradPath, MinCPXYZ);

				// Found!!
				if (DistPathToMinCP > 0.0 && DistPathToMinCP < Tolerance)
				{
					// First Bond found
					if (BondNum == -1)
					{
						BondNum = Bond;
					}
					else if (BondNum >= 0)  // Second bond found
					{
						BondNum = -2;
					}
					else  // more than 2 bonds found
					{
						BondNum = -2;
					}
				}

				GradPathDealloc(&GradPath);
			}
		}
	}

	if (IsOk == FALSE) BondNum = -1;

	ENSURE(BondNum > -2 && BondNum < NumBonds);
	return BondNum;
}

/*
	Loops over critical points in a given CP zone
	finding the minimum distance between CPs of two specified
	type, or between any two CPs of different type.

	returns true if operation succedded, false if not
*/
Boolean_t	CritPointsMinDistance(CritPoints_pa	CritPoints,
								  const Boolean_t	TypesSpecified,
								  const char	Type1,
								  const char	Type2,
								  double		*MinDist)
{
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(VALID_BOOLEAN(TypesSpecified));
	REQUIRE(VALID_REF(MinDist));
	if (TypesSpecified) REQUIRE(Type1 > -4 && Type1 < 14 &&
								Type2 > -4 && Type2 < 14);

	LgIndex_t	NumCrtPts = CritPointsGetCount(CritPoints);
	CHECK(NumCrtPts > 0);

	double	DXYZ = 0, djunk;
	char	CheckType1, CheckType2;
	XYZ_s	Pt1, Pt2;
	*MinDist = -1.0;

	for (LgIndex_t i = 0 ; i < NumCrtPts-1 ; i++)
	{
		CritPointsGetPoint(CritPoints, i, &Pt1.X, &Pt1.Y, &Pt1.Z, &djunk, 
			&CheckType1, &djunk, &djunk, &djunk);
		if (CheckType1 == Type1 || !TypesSpecified)
		{
			for (LgIndex_t j = 0 ; j < NumCrtPts ; j++)
			{
				if (j != i)
				{
					CritPointsGetPoint(CritPoints, j, &Pt2.X, &Pt2.Y, &Pt2.Z, 
						&djunk, &CheckType2, &djunk, &djunk, &djunk);
					if ((CheckType2 != CheckType1 && !TypesSpecified) || CheckType2 == Type2)
					{
						if (*MinDist <= 0) *MinDist = DistanceSquaredXYZ(Pt1, Pt2);
						else *MinDist = MIN(*MinDist, DistanceSquaredXYZ(Pt1, Pt2));
					}
				}
			}
		}
	}
	Boolean_t IsOk = (*MinDist > 0.0);
	if (IsOk) *MinDist = sqrt(*MinDist);

	return IsOk;
}

/*
	A simple way to check for duplicate
*/
Boolean_t	CritPointsSurfHasNoDuplicates()
{
	Boolean_t	HasDuplicates = FALSE;

	ENSURE(VALID_BOOLEAN(HasDuplicates));
	return HasDuplicates;
}