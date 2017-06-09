/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#include <cmath>
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "ARRLISTTOOLS.h"
#include "ODERUNGEKUTTA.h"
#include "NORMALS.h"
#include "SURFELEMMAP.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
#include "ENGINE.h"
#include "GEOMTOOLS.h"
#include "ELEMSHAPEFUNC.h"
#include "GRADPATH.h"

/*
 * Callback from RK2UpdateSolution to get the gradient values
 * using trilinear interpolation of nodal gradients.
 *
 * Input:
 *   *ClientData: Pointer to structure containing Zone/Var numbers
 *   *Position:   Current XYZ location.
 * Output:
 *   *Gradient:   Vector of gradient values in XYZ directions.
 *
 */
Boolean_t ChrgDensGrad(const void       *ClientData,
					   const double     *Position,
					   double           *Gradient)
{
	Boolean_t  IsOk = TRUE;
	EntIndex_t NumZones, NumVars;
	LgIndex_t  IMax, JMax, KMax;
	LgIndex_t  IIndex, JIndex, KIndex;
	EntIndex_t ZoneNum, UVarNum, VVarNum, WVarNum;
	Boolean_t  PeriodicBC;
	double     BufferLayer = 0.0;   // used for periodic
	// double     MinCoordinate = 1.0;
	double r, s, t;
	double W[8];

	ZoneVarInfo_pa ZoneVarInfo = (ZoneVarInfo_pa)(ClientData);
	double     XBeg = ZoneVarInfo->XBeg;
	double     XEnd = ZoneVarInfo->XEnd;
	double     YBeg = ZoneVarInfo->YBeg;
	double     YEnd = ZoneVarInfo->YEnd; 
	double     ZBeg = ZoneVarInfo->ZBeg;
	double     ZEnd = ZoneVarInfo->ZEnd;
	double     DXZoneBuffer = 0.0;
	double     DYZoneBuffer = 0.0;
	double     DZZoneBuffer = 0.0;

	REQUIRE(VALID_REF(ClientData));
	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	ZoneNum    = ZoneVarInfo->ZoneNum;
	UVarNum    = ZoneVarInfo->UVarNum;
	VVarNum    = ZoneVarInfo->VVarNum;
	WVarNum    = ZoneVarInfo->WVarNum;
	PeriodicBC = ZoneVarInfo->PeriodicBC;

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);

	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));

	REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
	REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
	REQUIRE(WVarNum > 0 && WVarNum <= NumVars);

	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */

	if (PeriodicBC)
	{
		DXZoneBuffer = (XEnd - XBeg) / (double)(MAX(IMax - 1, 1));
		DYZoneBuffer = (YEnd - YBeg) / (double)(MAX(JMax - 1, 1));
		DZZoneBuffer = (ZEnd - ZBeg) / (double)(MAX(KMax - 1, 1));
	}


	/* Is it in the solution domain */
	if (Position[0] < XBeg - DXZoneBuffer || Position[0] > XEnd + DXZoneBuffer ||
		Position[1] < YBeg - DYZoneBuffer || Position[1] > YEnd + DYZoneBuffer ||
		Position[2] < ZBeg - DZZoneBuffer || Position[2] > ZEnd + DZZoneBuffer ) IsOk = FALSE;

	if (IsOk)
	{
		// Note: Assumes rectangular with I,J,K aligning with X,Y,Z
		/*
		double ICoordPos = 1.0 + (double)(IMax - 1) * (Position[0] - XBeg) / (XEnd - XBeg);
		double JCoordPos = 1.0 + (double)(JMax - 1) * (Position[1] - YBeg) / (YEnd - YBeg);
		double KCoordPos = 1.0 + (double)(KMax - 1) * (Position[2] - ZBeg) / (ZEnd - ZBeg);

		IIndex = MAX(MIN((LgIndex_t)ICoordPos, IMax-1),1);
		JIndex = MAX(MIN((LgIndex_t)JCoordPos, JMax-1),1);
		KIndex = MAX(MIN((LgIndex_t)KCoordPos, KMax-1),1);

		r = 2.0 * (ICoordPos - (double)IIndex) - 1.0;
		s = 2.0 * (JCoordPos - (double)JIndex) - 1.0;
		t = 2.0 * (KCoordPos - (double)KIndex) - 1.0;
		*/
		IsOk = RectGridBrickXYZtoIJKRST(Position[0], Position[1], Position[2], IMax, JMax, KMax,
										XBeg, XEnd, YBeg, YEnd, ZBeg, ZEnd, PeriodicBC, 
										&IIndex, &JIndex, &KIndex, &r, &s, &t);

		IsOk = BrickTrilinearWeight(r, s, t, W);

		if (IsOk)
		{
			int          ic;
			LgIndex_t    Index[8];
			FieldData_pa UVarFDPtr = NULL;
			FieldData_pa VVarFDPtr = NULL;
			FieldData_pa WVarFDPtr = NULL;

			// begin-temp
			// double xtmp = 0.0;
			// double ytmp = 0.0;
			// double ztmp = 0.0;
			// FieldData_pa XVarFDPtr = NULL;
			// FieldData_pa YVarFDPtr = NULL;
			// FieldData_pa ZVarFDPtr = NULL;
			// XVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, 1);
			// YVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, 2);
			// ZVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, 3);
			// end-temp

			UVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, UVarNum);
			VVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, VVarNum);
			WVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, WVarNum);

			/*
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
			  */
			Index[0] = IndexFromIJK(IIndex    , JIndex    , KIndex    , IMax, JMax, KMax, PeriodicBC);
			Index[1] = IndexFromIJK(IIndex + 1, JIndex    , KIndex    , IMax, JMax, KMax, PeriodicBC);
			Index[2] = IndexFromIJK(IIndex + 1, JIndex + 1, KIndex    , IMax, JMax, KMax, PeriodicBC);
			Index[3] = IndexFromIJK(IIndex    , JIndex + 1, KIndex    , IMax, JMax, KMax, PeriodicBC);
			Index[4] = IndexFromIJK(IIndex    , JIndex    , KIndex + 1, IMax, JMax, KMax, PeriodicBC);
			Index[5] = IndexFromIJK(IIndex + 1, JIndex    , KIndex + 1, IMax, JMax, KMax, PeriodicBC);
			Index[6] = IndexFromIJK(IIndex + 1, JIndex + 1, KIndex + 1, IMax, JMax, KMax, PeriodicBC);
			Index[7] = IndexFromIJK(IIndex    , JIndex + 1, KIndex + 1, IMax, JMax, KMax, PeriodicBC);

			// begin-temp
			// for (ic=0; ic<8; ic++)
			//   {
			//     xtmp += W[ic] * TecUtilDataValueGetByRef(XVarFDPtr, Index[ic]);
			//     ytmp += W[ic] * TecUtilDataValueGetByRef(YVarFDPtr, Index[ic]);
			//     ztmp += W[ic] * TecUtilDataValueGetByRef(ZVarFDPtr, Index[ic]);
			//   }
			// end-temp

			Gradient[0] = W[0] * TecUtilDataValueGetByRef(UVarFDPtr, Index[0]);
			Gradient[1] = W[0] * TecUtilDataValueGetByRef(VVarFDPtr, Index[0]);
			Gradient[2] = W[0] * TecUtilDataValueGetByRef(WVarFDPtr, Index[0]);
			for (ic = 1; ic < 8; ic++)
			{
				Gradient[0] += W[ic] * TecUtilDataValueGetByRef(UVarFDPtr, Index[ic]);
				Gradient[1] += W[ic] * TecUtilDataValueGetByRef(VVarFDPtr, Index[ic]);
				Gradient[2] += W[ic] * TecUtilDataValueGetByRef(WVarFDPtr, Index[ic]);
			}
		}
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}








/*
 * Function to determine the coordinates and coordinate ranges for an
 * FEBrick designed to encompase a surface FEQuad. The brick extends
 * inward and outward along the node normals from the FEQuad. The brick
 * is used to compute surface gradient paths by computing volume gradient
 * paths
 *
 * Input:
 *   Index[4]: Four corner nodes of an element
 *   Normals:  Structure containing the normal vectors at all nodes
 *             of a zone.
 *   XVarFDPtr, YVarFDPtr, ZVarFDPtr: Field data pointers for the
 *             nodes of a zone.
 * Output:
 *   XCell[8], YCell[8], ZCell[8]: Coordinate of the nodes of the FEBrick.
 *   XMin, XMax, YMin, YMax, ZMin, ZMax: Coordinate ranges for the FEBrick.
 *
 */
Boolean_t SetUpBrickForSurfQuad(LgIndex_t    Index[4],
								Normals_pa   Normals,
								FieldData_pa XVarFDPtr,
								FieldData_pa YVarFDPtr,
								FieldData_pa ZVarFDPtr,
								double       XCell[8],
								double       YCell[8],
								double       ZCell[8],
								double      *XMin,
								double      *XMax,
								double      *YMin,
								double      *YMax,
								double      *ZMin,
								double      *ZMax)
{
	Boolean_t IsOk = TRUE;
	Boolean_t IsQuad = TRUE;
	int NumNodesInElem;
	Boolean_t InCell = FALSE;

	REQUIRE(VALID_REF(Index));
	REQUIRE(NormalsIsValid(Normals));
	REQUIRE(VALID_REF(XVarFDPtr));
	REQUIRE(VALID_REF(YVarFDPtr));
	REQUIRE(VALID_REF(ZVarFDPtr));
	REQUIRE(VALID_REF(XCell));
	REQUIRE(VALID_REF(YCell));
	REQUIRE(VALID_REF(ZCell));
	REQUIRE(VALID_REF(XMin));
	REQUIRE(VALID_REF(XMax));
	REQUIRE(VALID_REF(YMin));
	REQUIRE(VALID_REF(YMax));
	REQUIRE(VALID_REF(ZMin));
	REQUIRE(VALID_REF(ZMax));

	// Create the corner coordinates and ranges for an FEBrick encompasing
	// a surface FEQuad
	if (IsOk)
	{
		// Look for degenerate elements (triangles),
		// Adjust Index[] so first three items are Nodes of triangle
		// TODO: actually using FEBrick for all so don't adjust

		NumNodesInElem = 4;
		/*
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
		*/
	}

	if (IsOk)
	{
		int ic;
		for (ic = 0; ic < 4; ic++)
		{
			int icp4 = ic + 4;
			LgIndex_t Node = Index[ic];
			// TEMP
			double DN = 0.5;
			XYZ_s Normal = NormalsGetNormalForNode(Normals, Node);
			XCell[ic] = TecUtilDataValueGetByRef(XVarFDPtr, Node) + DN * Normal.X;
			YCell[ic] = TecUtilDataValueGetByRef(YVarFDPtr, Node) + DN * Normal.Y;
			ZCell[ic] = TecUtilDataValueGetByRef(ZVarFDPtr, Node) + DN * Normal.Z;
			XCell[icp4] = TecUtilDataValueGetByRef(XVarFDPtr, Node) - DN * Normal.X;
			YCell[icp4] = TecUtilDataValueGetByRef(YVarFDPtr, Node) - DN * Normal.Y;
			ZCell[icp4] = TecUtilDataValueGetByRef(ZVarFDPtr, Node) - DN * Normal.Z;
		}


		// Compute the min/max coords
		*XMin = XCell[0];
		*YMin = YCell[0];
		*ZMin = ZCell[0];
		*XMax = *XMin;
		*YMax = *YMin;
		*ZMax = *ZMin;

		for (ic = 0; IsOk && ic < 8; ic++)
		{
			*XMin = MIN(*XMin, XCell[ic]);
			*YMin = MIN(*YMin, YCell[ic]);
			*ZMin = MIN(*ZMin, ZCell[ic]);
			*XMax = MAX(*XMax, XCell[ic]);
			*YMax = MAX(*YMax, YCell[ic]);
			*ZMax = MAX(*ZMax, ZCell[ic]);
		}
	}
	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}





/*
 * Callback from RK2UpdateSolution to get the gradient values
 * using trilinear interpolation of nodal gradients.
 *
 * Input:
 *   *ClientData: Pointer to structure containing Zone/Var numbers
 *   *Position:   Current XYZ location.
 * Output:
 *   *Gradient:   Vector of gradient values in XYZ directions.
 *
 */
Boolean_t SurfTotGradGrad(const void       *ClientData,
						  const double     *Position,
						  double           *Gradient)
{
	Boolean_t  IsOk = TRUE;
	EntIndex_t NumZones, NumVars;
	LgIndex_t  IMax, JMax, KMax;
	EntIndex_t ZoneNum, UVarNum, VVarNum, WVarNum, XVarNum, YVarNum, ZVarNum;
	ZoneType_e ZoneType;
	Normals_pa Normals = NULL;
	SurfElemMap_pa SurfElemMap = NULL;
	FieldData_pa XVarFDPtr = NULL;
	FieldData_pa YVarFDPtr = NULL;
	FieldData_pa ZVarFDPtr = NULL;
	FieldData_pa UVarFDPtr = NULL;
	FieldData_pa VVarFDPtr = NULL;
	FieldData_pa WVarFDPtr = NULL;
	Boolean_t  PeriodicBC;
	Boolean_t  ReverseDir = FALSE;
	LgIndex_t  LastElemUsed;
	double     BufferLayer = 0.0;   // used for periodic
	double     MinCoordinate = 1.0;
	double r = 0.0;
	double s = 0.0;
	double t = 0.0;
	double W[8];

	ZoneVarInfo_pa ZoneVarInfo = (ZoneVarInfo_pa)(ClientData);

	REQUIRE(VALID_REF(ClientData));
	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	ZoneNum    = ZoneVarInfo->ZoneNum;
	UVarNum    = ZoneVarInfo->UVarNum;
	VVarNum    = ZoneVarInfo->VVarNum;
	WVarNum    = ZoneVarInfo->WVarNum;
	PeriodicBC = ZoneVarInfo->PeriodicBC;
	LastElemUsed = ZoneVarInfo->LastElemUsed;
	if (NormalsIsValid(ZoneVarInfo->Normals))
		Normals = ZoneVarInfo->Normals;
	if (SurfElemMapIsValid(ZoneVarInfo->SurfElemMap))
		SurfElemMap = ZoneVarInfo->SurfElemMap;

	// Impacts the sign of the "redirect to face" term in gradient
	if (ZoneVarInfo->PathDir == StreamDir_Reverse) ReverseDir = TRUE;

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);

	// REQUIRE(TecUtilZoneIsOrdered(ZoneNum));

	REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
	REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
	REQUIRE(WVarNum > 0 && WVarNum <= NumVars);

	REQUIRE(NormalsIsValid(Normals));

	ZoneType = TecUtilZoneGetType(ZoneNum);

	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	if (TecUtilZoneIsOrdered(ZoneNum))
	{
		REQUIRE(IMax > 1 && JMax > 1 && KMax == 1); /* Must be 2D surface */
	}
	else  // Zone is FE - make sure it is a surface
	{
		REQUIRE(IMax > 2 && JMax > 0);
		REQUIRE(ZoneType == ZoneType_FETriangle || ZoneType == ZoneType_FEQuad);
	}


	// Find field data pointers
	UVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, UVarNum);
	VVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, VVarNum);
	WVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, WVarNum);

	XVarNum = TecUtilVarGetNumByAssignment('X');
	YVarNum = TecUtilVarGetNumByAssignment('Y');
	ZVarNum = TecUtilVarGetNumByAssignment('Z');

	if (XVarNum > 0 && XVarNum <= NumVars)
		XVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, XVarNum);
	else
		IsOk = FALSE;

	if (YVarNum > 0 && YVarNum <= NumVars)
		YVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, YVarNum);
	else
		IsOk = FALSE;

	if (YVarNum > 0 && YVarNum <= NumVars)
		ZVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, ZVarNum);
	else
		IsOk = FALSE;

	// IJ-Ordered surface
	if (TecUtilZoneIsOrdered(ZoneNum))
	{
		// Not yet implemented
		IsOk = FALSE;
		CHECK(0);
	}
	else  // FE-Triangle or FE-Quad surface
	{
		Boolean_t IsQuad = TRUE;
		LgIndex_t Index[4];
		int NumNodesInElem;
		double XCell[8],  YCell[8],  ZCell[8];
		double XMin, XMax, YMin, YMax, ZMin, ZMax;
		Boolean_t InCell = FALSE;


		// Find the closest cell and the weights
		// In LastElemUsed?
		if (IsOk)
		{
			ZoneType_e ZoneType = TecUtilZoneGetType(ZoneNum);

			if (ZoneType == ZoneType_FEQuad)
			{
				Index[0] = TecUtilDataNodeGetByZone(ZoneNum, LastElemUsed, 1);
				Index[1] = TecUtilDataNodeGetByZone(ZoneNum, LastElemUsed, 2);
				Index[2] = TecUtilDataNodeGetByZone(ZoneNum, LastElemUsed, 3);
				Index[3] = TecUtilDataNodeGetByZone(ZoneNum, LastElemUsed, 4);
			}
			else
				IsOk = FALSE;
		}
		if (IsOk)
		{
			// Look for degenerate elements (triangles),
			// Adjust Index[] so first three items are Nodes of triangle
			// TODO: actually using FEBrick for all so don't adjust

			NumNodesInElem = 4;
			/*
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
			*/
		}

		if (IsOk)
		{
			int ic;
			for (ic = 0; ic < 4; ic++)
			{
				int icp4 = ic + 4;
				LgIndex_t Node = Index[ic];
				// TEMP
				double DN = 0.5;
				XYZ_s Normal = NormalsGetNormalForNode(Normals, Node);
				XCell[ic] = TecUtilDataValueGetByRef(XVarFDPtr, Node) + DN * Normal.X;
				YCell[ic] = TecUtilDataValueGetByRef(YVarFDPtr, Node) + DN * Normal.Y;
				ZCell[ic] = TecUtilDataValueGetByRef(ZVarFDPtr, Node) + DN * Normal.Z;
				XCell[icp4] = TecUtilDataValueGetByRef(XVarFDPtr, Node) - DN * Normal.X;
				YCell[icp4] = TecUtilDataValueGetByRef(YVarFDPtr, Node) - DN * Normal.Y;
				ZCell[icp4] = TecUtilDataValueGetByRef(ZVarFDPtr, Node) - DN * Normal.Z;
			}


			// Compute the min/max coords
			XMin = XCell[0];
			YMin = YCell[0];
			ZMin = ZCell[0];
			XMax = XMin;
			YMax = YMin;
			ZMax = ZMin;

			for (ic = 0; IsOk && ic < 8; ic++)
			{
				XMin = MIN(XMin, XCell[ic]);
				YMin = MIN(YMin, YCell[ic]);
				ZMin = MIN(ZMin, ZCell[ic]);
				XMax = MAX(XMax, XCell[ic]);
				YMax = MAX(YMax, YCell[ic]);
				ZMax = MAX(ZMax, ZCell[ic]);
			}

			if (IsOk)
			{

				// if (IsQuad)
				//   {
				InCell = TrilinearBrickNaturalCoord(Position[0], Position[1], Position[2],
													XCell, YCell, ZCell, XMin, YMin, ZMin,
													XMax, YMax, ZMax, &r, &s, &t, W);
				// InCell = BilinearQuadNaturalCoord(XP, YP, XPCell, YPCell, XPMin,YPMin, XPMax,YPMax, &r, &s, W);
				//   }
				// else
				//   {
				// InCell = LinearTriangleNaturalCoord(XP, YP, XPCell, YPCell, XPMin,YPMin, XPMax, YPMax, XPCell, YPCell, &r, &s, W);
				//   }
			}
		}
		// If not InCell, look in first layer of neighboring cells
		if (IsOk && !InCell)
		{
			int ic;
			// ArrList of elements checked
			ArrList_pa ElemChecked = ArrListAlloc(20, ArrListType_Long);
			ArrListItem_u Item;

			// First elem in list is LastElemUsed
			Item.Long = LastElemUsed;
			IsOk = ArrListAppendItem(ElemChecked, Item);

			// For each node in LastElemUsed, cycle through the elements
			// in it's SurfElemMap checking elements that have not already
			// been checked.
			for (ic = 0; ic < 4; ic++)
			{
				LgIndex_t Node = Index[ic];
				LgIndex_t NumElems = SurfElemMapGetElemCountForNode(SurfElemMap, Node);
				LgIndex_t ie;
				for (ie = 0; !InCell && IsOk && ie < NumElems; ie++)
				{
					LgIndex_t Elem = SurfElemMapGetElemNumForNodeOffset(SurfElemMap, Node, ie);
					if (ArrListIsUniqueLongItem(ElemChecked, Elem))
					{
						// Try the element
						ZoneType_e ZoneType = TecUtilZoneGetType(ZoneNum);
						LgIndex_t IndexElem[4];
						ArrListAppendUniqueLongItem(ElemChecked, Elem);

						if (ZoneType == ZoneType_FEQuad)
						{
							IndexElem[0] = TecUtilDataNodeGetByZone(ZoneNum, Elem, 1);
							IndexElem[1] = TecUtilDataNodeGetByZone(ZoneNum, Elem, 2);
							IndexElem[2] = TecUtilDataNodeGetByZone(ZoneNum, Elem, 3);
							IndexElem[3] = TecUtilDataNodeGetByZone(ZoneNum, Elem, 4);
						}
						else
							IsOk = FALSE;

						if (IsOk)
						{
							IsOk = SetUpBrickForSurfQuad(IndexElem, Normals, XVarFDPtr,
														 YVarFDPtr, ZVarFDPtr, XCell, YCell,
														 ZCell, &XMin, &XMax, &YMin, &YMax, &ZMin, &ZMax);
						}

						if (IsOk)
						{
						   InCell = TrilinearBrickNaturalCoord(Position[0], Position[1], Position[2],
																XCell, YCell, ZCell, XMin, YMin, ZMin,
																XMax, YMax, ZMax, &r, &s, &t, W);

							if (InCell)
							{
								LastElemUsed = Elem;
								ZoneVarInfo->LastElemUsed = LastElemUsed;
								Index[0] = IndexElem[0];
								Index[1] = IndexElem[1];
								Index[2] = IndexElem[2];
								Index[3] = IndexElem[3];
							}
						}
					}
				}
			}

			// If still not InCell, look in second layer of neighboring cells
			if (IsOk && !InCell)
			{
				int ic, ec;
				// ArrList of elements checked
				ArrList_pa Elem2ndChecked = ArrListAlloc(20, ArrListType_Long);
				ArrListItem_u Item;

				LgIndex_t Num1rstLayerElem = ArrListGetCount(ElemChecked);

				// First elem in list is LastElemUsed
				Item.Long = LastElemUsed;
				IsOk = ArrListAppendItem(Elem2ndChecked, Item);

				// For each node in the elements of ElemChecked, cycle through the elements
				// in it's SurfElemMap checking elements that have not already
				// been checked.
				for (ec = 1; IsOk && !InCell && ec < Num1rstLayerElem; ec++)
				{
					LgIndex_t Index1rstElem[4];
					ArrListItem_u Item;
					LgIndex_t Elem1rstLyr;

					Item = ArrListGetItem(ElemChecked, ec);
					Elem1rstLyr = Item.Long;

					if (ZoneType == ZoneType_FEQuad)
					{
						Index1rstElem[0] = TecUtilDataNodeGetByZone(ZoneNum, Elem1rstLyr, 1);
						Index1rstElem[1] = TecUtilDataNodeGetByZone(ZoneNum, Elem1rstLyr, 2);
						Index1rstElem[2] = TecUtilDataNodeGetByZone(ZoneNum, Elem1rstLyr, 3);
						Index1rstElem[3] = TecUtilDataNodeGetByZone(ZoneNum, Elem1rstLyr, 4);
					}
					else
						IsOk = FALSE;

					// For each node in Elem1rstLyr (number of a first layer element), cycle through
					// the elements in it's SurfElemMap checking elements that have not already
					// been checked. Compare against both the first layer list (ElemChecked) and
					// the second layer list (Elem2ndChecked).
					for (ic = 0; !InCell && IsOk && ic < 4; ic++)
					{
						LgIndex_t Node = Index1rstElem[ic];
						LgIndex_t NumElems = SurfElemMapGetElemCountForNode(SurfElemMap, Node);
						LgIndex_t ie;
						for (ie = 0; !InCell && IsOk && ie < NumElems; ie++)
						{
							LgIndex_t Elem = SurfElemMapGetElemNumForNodeOffset(SurfElemMap, Node, ie);
							if (ArrListIsUniqueLongItem(ElemChecked, Elem) &&
								ArrListIsUniqueLongItem(Elem2ndChecked, Elem))
							{
								// Try the element
								ZoneType_e ZoneType = TecUtilZoneGetType(ZoneNum);
								LgIndex_t IndexElem[4];
								ArrListAppendUniqueLongItem(Elem2ndChecked, Elem);

								if (ZoneType == ZoneType_FEQuad)
								{
									IndexElem[0] = TecUtilDataNodeGetByZone(ZoneNum, Elem, 1);
									IndexElem[1] = TecUtilDataNodeGetByZone(ZoneNum, Elem, 2);
									IndexElem[2] = TecUtilDataNodeGetByZone(ZoneNum, Elem, 3);
									IndexElem[3] = TecUtilDataNodeGetByZone(ZoneNum, Elem, 4);
								}
								else
									IsOk = FALSE;

								if (IsOk)
								{
									IsOk = SetUpBrickForSurfQuad(IndexElem, Normals, XVarFDPtr,
																 YVarFDPtr, ZVarFDPtr, XCell, YCell,
																 ZCell, &XMin, &XMax, &YMin, &YMax, &ZMin, &ZMax);
								}

								if (IsOk)
								{
									InCell = TrilinearBrickNaturalCoord(Position[0], Position[1], Position[2],
																		XCell, YCell, ZCell, XMin, YMin, ZMin,
																		XMax, YMax, ZMax, &r, &s, &t, W);

									if (InCell)
									{
										LastElemUsed = Elem;
										ZoneVarInfo->LastElemUsed = LastElemUsed;
										Index[0] = IndexElem[0];
										Index[1] = IndexElem[1];
										Index[2] = IndexElem[2];
										Index[3] = IndexElem[3];
									}
								}
							}
						}
					}
				}
				ArrListDealloc(&Elem2ndChecked);
			}

			// TODO: if not InCell, do a global test to find the closest element

			ArrListDealloc(&ElemChecked);
		}

		// Compute the gradient
		if (IsOk)
		{
			if (InCell)
			{
				int ic;
				double     SignFeedback = 1.0;
				if (ReverseDir) SignFeedback = -1.0;

				Gradient[0] = 0.0;
				Gradient[1] = 0.0;
				Gradient[2] = 0.0;
				for (ic = 0; IsOk && ic < 4; ic++)
				{
					LgIndex_t  Node = Index[ic];
					double     Weight = W[ic] + W[ic+4];
					// double     WeightDiff = 2.0 * SignFeedback * (W[ic+4] - W[ic]);
					double     WeightDiff = 4.0 * SignFeedback * (W[ic+4] - W[ic]);
					double     UVal = TecUtilDataValueGetByRef(UVarFDPtr, Node);
					double     VVal = TecUtilDataValueGetByRef(VVarFDPtr, Node);
					double     WVal = TecUtilDataValueGetByRef(WVarFDPtr, Node);
					double     GradTotNode = ABS(UVal) + ABS(VVal) + ABS(WVal);  // Approximate
					XYZ_s Normal = NormalsGetNormalForNode(Normals, Node);

					Gradient[0] += Weight * UVal;
					Gradient[1] += Weight * VVal;
					Gradient[2] += Weight * WVal;

					// Adjust to push line toward surface (t=0 of FEBrick where W[ic] == W[ic+4]
					Gradient[0] += WeightDiff * GradTotNode * Normal.X;
					Gradient[1] += WeightDiff * GradTotNode * Normal.Y;
					Gradient[2] += WeightDiff * GradTotNode * Normal.Z;

				}

			}
			else
			{
				// TODO: Not in LastCellUsed - expand to next layer of cells
				IsOk = FALSE;
			}
		}





		// Project the vector onto the surface (tangent vector) - at least for FE-Triangle
	}

	// If it is periodic, set a buffer layer of 1.0
	// if (PeriodicBC)
	//   {
	//     BufferLayer = 1.0;
	//     MinCoordinate = 1.0 - BufferLayer;
	//   }

	/* Is it in the solution domain */
	// if ( Position[0] < MinCoordinate || Position[0] > ((double)IMax + 2*BufferLayer) ||
	//      Position[1] < MinCoordinate || Position[1] > ((double)JMax + 2*BufferLayer) ||
	//      Position[2] < MinCoordinate || Position[2] > ((double)KMax + 2*BufferLayer) ) IsOk = FALSE;

	// if (IsOk)
	//   {
	//     IIndex = MAX( (LgIndex_t)Position[0], (LgIndex_t)MinCoordinate);
	//     JIndex = MAX( (LgIndex_t)Position[1], (LgIndex_t)MinCoordinate);
	//     KIndex = MAX( (LgIndex_t)Position[2], (LgIndex_t)MinCoordinate);

	//     r = 2.0 * ( Position[0] - (double)IIndex ) - 1.0;
	//     s = 2.0 * ( Position[1] - (double)JIndex ) - 1.0;
	//     t = 2.0 * ( Position[2] - (double)KIndex ) - 1.0;

	//     IsOk = BrickTrilinearWeight(r, s, t, W);

	//     if (IsOk)
	//       {
	//         int          ic;
	//         LgIndex_t    Index[8];
	//         FieldData_pa UVarFDPtr = NULL;
	//         FieldData_pa VVarFDPtr = NULL;
	//         FieldData_pa WVarFDPtr = NULL;


	//         UVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, UVarNum);
	//         VVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, VVarNum);
	//         WVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, WVarNum);

	//         Index[0] = IndexFromIJK(IIndex    , JIndex    , KIndex    , IMax, JMax, KMax, PeriodicBC);
	//         Index[1] = IndexFromIJK(IIndex + 1, JIndex    , KIndex    , IMax, JMax, KMax, PeriodicBC);
	//         Index[2] = IndexFromIJK(IIndex + 1, JIndex + 1, KIndex    , IMax, JMax, KMax, PeriodicBC);
	//         Index[3] = IndexFromIJK(IIndex    , JIndex + 1, KIndex    , IMax, JMax, KMax, PeriodicBC);
	//         Index[4] = IndexFromIJK(IIndex    , JIndex    , KIndex + 1, IMax, JMax, KMax, PeriodicBC);
	//         Index[5] = IndexFromIJK(IIndex + 1, JIndex    , KIndex + 1, IMax, JMax, KMax, PeriodicBC);
	//         Index[6] = IndexFromIJK(IIndex + 1, JIndex + 1, KIndex + 1, IMax, JMax, KMax, PeriodicBC);
	//         Index[7] = IndexFromIJK(IIndex    , JIndex + 1, KIndex + 1, IMax, JMax, KMax, PeriodicBC);

	//         Gradient[0] = W[0] * TecUtilDataValueGetByRef(UVarFDPtr, Index[0]);
	//         Gradient[1] = W[0] * TecUtilDataValueGetByRef(VVarFDPtr, Index[0]);
	//        Gradient[2] = W[0] * TecUtilDataValueGetByRef(WVarFDPtr, Index[0]);
	//         for (ic=1; ic<8; ic++)
	//           {
	//             Gradient[0] += W[ic] * TecUtilDataValueGetByRef(UVarFDPtr, Index[ic]);
	//             Gradient[1] += W[ic] * TecUtilDataValueGetByRef(VVarFDPtr, Index[ic]);
	//             Gradient[2] += W[ic] * TecUtilDataValueGetByRef(WVarFDPtr, Index[ic]);
	//           }
	//       }
	//   }


	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}





/*
 * Compute the value of charge density at an X, Y, Z, position
 * using trilinear interpolation of nodal charge densities.
 *
 * Parameters:
 *   *ClientData: Pointer to structure containing Zone/Var numbers
 *   *Position:   Current XYZ location.
 *
 * Return:
 *   ChrgDens:    Scalar value of charge density at Position. Negative
 *                if it doesn't work.
 *
 */
double    ChrgDensTrilinearInterp(const ZoneVarInfo_pa  ZoneVarInfo,
								  const double         *Position)
{
	double     ChrgDens = -1.0;
	Boolean_t  IsOk = TRUE;
	EntIndex_t NumZones, NumVars;
	LgIndex_t  IMax, JMax, KMax;
	LgIndex_t  IIndex, JIndex, KIndex;
	EntIndex_t ZoneNum, ChrgDensVarNum;
	Boolean_t  PeriodicBC;
	double     BufferLayer = 0.0;   // used for periodic
	double     MinCoordinate = 1.0;
	double r, s, t;
	double W[8];

	REQUIRE(VALID_REF(ZoneVarInfo));
	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	ZoneNum = ZoneVarInfo->ZoneNum;
	ChrgDensVarNum = ZoneVarInfo->ChrgDensVarNum;
	PeriodicBC = ZoneVarInfo->PeriodicBC;

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	REQUIRE(VALID_BOOLEAN(PeriodicBC));

	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));

	REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);

	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */

	// If it is periodic, set a buffer layer of 1.0
	if (PeriodicBC)
	{
		BufferLayer = 1.0;
		MinCoordinate = 1.0 - BufferLayer;
	}

	/* Is it in the solution domain */
	/*
	if (Position[0] < MinCoordinate || Position[0] > ((double)IMax + 2*BufferLayer) ||
		Position[1] < MinCoordinate || Position[1] > ((double)JMax + 2*BufferLayer) ||
		Position[2] < MinCoordinate || Position[2] > ((double)KMax + 2*BufferLayer)) IsOk = FALSE;
		*/
	IsOk = RectGridBrickXYZtoIJKRST(Position[0], Position[1], Position[2], IMax, JMax, KMax,
									ZoneVarInfo->XBeg, ZoneVarInfo->XEnd, ZoneVarInfo->YBeg, 
									ZoneVarInfo->YEnd, ZoneVarInfo->ZBeg, ZoneVarInfo->ZEnd,
									ZoneVarInfo->PeriodicBC,
									&IIndex, &JIndex, &KIndex, &r, &s, &t);

	if (IsOk)
	{
		/*
		IIndex = MAX((LgIndex_t)Position[0], (LgIndex_t)MinCoordinate);
		JIndex = MAX((LgIndex_t)Position[1], (LgIndex_t)MinCoordinate);
		KIndex = MAX((LgIndex_t)Position[2], (LgIndex_t)MinCoordinate);

		r = 2.0 * (Position[0] - (double)IIndex) - 1.0;
		s = 2.0 * (Position[1] - (double)JIndex) - 1.0;
		t = 2.0 * (Position[2] - (double)KIndex) - 1.0;
		*/

		IsOk = BrickTrilinearWeight(r, s, t, W);

		if (IsOk)
		{
			int          ic;
			LgIndex_t    Index[8];
			FieldData_pa CDVarFDPtr = NULL;

			CDVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, ChrgDensVarNum);

			Index[0] = IndexFromIJK(IIndex    , JIndex    , KIndex    , IMax, JMax, KMax, PeriodicBC);
			Index[1] = IndexFromIJK(IIndex + 1, JIndex    , KIndex    , IMax, JMax, KMax, PeriodicBC);
			Index[2] = IndexFromIJK(IIndex + 1, JIndex + 1, KIndex    , IMax, JMax, KMax, PeriodicBC);
			Index[3] = IndexFromIJK(IIndex    , JIndex + 1, KIndex    , IMax, JMax, KMax, PeriodicBC);
			Index[4] = IndexFromIJK(IIndex    , JIndex    , KIndex + 1, IMax, JMax, KMax, PeriodicBC);
			Index[5] = IndexFromIJK(IIndex + 1, JIndex    , KIndex + 1, IMax, JMax, KMax, PeriodicBC);
			Index[6] = IndexFromIJK(IIndex + 1, JIndex + 1, KIndex + 1, IMax, JMax, KMax, PeriodicBC);
			Index[7] = IndexFromIJK(IIndex    , JIndex + 1, KIndex + 1, IMax, JMax, KMax, PeriodicBC);

			ChrgDens = W[0] * TecUtilDataValueGetByRef(CDVarFDPtr, Index[0]);
			for (ic = 1; ic < 8; ic++)
			{
				ChrgDens += W[ic] * TecUtilDataValueGetByRef(CDVarFDPtr, Index[ic]);
			}
		}
	}

	if (IsOk == FALSE) ChrgDens = -1.0;

	ENSURE(ChrgDens > 0.0 || ChrgDens == -1.0);
	return ChrgDens;
}





/*
 * Compute the maximum gradient values in the current cell.
 *
 * Input:
 *   *ClientData: Pointer to structure containing Zone/Var numbers
 *   *Position:   Current XYZ location.
 * Output:
 *   *MaxGradient: Scalar value of maximum gradient.
 *
 */
Boolean_t MaxChrgDensGrad(const void       *ClientData,
						  const double     *Position,
						  double           *MaxGradient,
						  double           *AveGradient)
{
	Boolean_t  IsOk = TRUE;
	EntIndex_t NumZones, NumVars;
	LgIndex_t  IMax, JMax, KMax;
	LgIndex_t  IIndex, JIndex, KIndex;
	EntIndex_t ZoneNum, UVarNum, VVarNum, WVarNum;
	LgIndex_t  LastElemUsed;
	Boolean_t  PeriodicBC;
	double     BufferLayer = 0.0;   // used for periodic
	double     MinCoordinate = 1.0;
	FieldData_pa UVarFDPtr = NULL;
	FieldData_pa VVarFDPtr = NULL;
	FieldData_pa WVarFDPtr = NULL;

	ZoneVarInfo_pa ZoneVarInfo = (ZoneVarInfo_pa)(ClientData);

	REQUIRE(VALID_REF(ClientData));
	REQUIRE(VALID_REF(Position));
	REQUIRE(VALID_REF(MaxGradient));
	REQUIRE(VALID_REF(AveGradient));
	/* Get num zones in dataset */
	if (IsOk)

		// Getting information about all zones so can reference the zone of 
		// interest next
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	ZoneNum = ZoneVarInfo->ZoneNum;
	UVarNum = ZoneVarInfo->UVarNum;
	VVarNum = ZoneVarInfo->VVarNum;
	WVarNum = ZoneVarInfo->WVarNum;
	PeriodicBC = ZoneVarInfo->PeriodicBC;
	LastElemUsed = MAX(ZoneVarInfo->LastElemUsed, 1);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);

	// REQUIRE(TecUtilZoneIsOrdered(ZoneNum));

	REQUIRE(UVarNum > 0 && UVarNum <= NumVars);
	REQUIRE(VVarNum > 0 && VVarNum <= NumVars);
	REQUIRE(WVarNum > 0 && WVarNum <= NumVars);
	REQUIRE(VALID_BOOLEAN(PeriodicBC));

	// Fetching dimensional information about zone of interest
	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */

	// With zone of interest verified as a 3D volume, fetch the rho information
	// from the zone
	UVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, UVarNum);
	VVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, VVarNum);
	WVarFDPtr = TecUtilDataValueGetReadableRef(ZoneNum, WVarNum);

	if (TecUtilZoneIsOrdered(ZoneNum))
	{
		double XBeg = ZoneVarInfo->XBeg;
		double XEnd = ZoneVarInfo->XEnd;
		double YBeg = ZoneVarInfo->YBeg;
		double YEnd = ZoneVarInfo->YEnd;
		double ZBeg = ZoneVarInfo->ZBeg;
		double ZEnd = ZoneVarInfo->ZEnd;
		double DXZoneBuffer = 0.0;
		double DYZoneBuffer = 0.0;
		double DZZoneBuffer = 0.0;
		double r, s, t;
 
		// If it is periodic, set a buffer layer of one cell
		if (PeriodicBC)
		{
			DXZoneBuffer = (XEnd - XBeg) / (double)(MAX(IMax - 1, 1));
			DYZoneBuffer = (YEnd - YBeg) / (double)(MAX(JMax - 1, 1));
			DZZoneBuffer = (ZEnd - ZBeg) / (double)(MAX(KMax - 1, 1));
		}


		/* Is it in the solution domain */
		/*
		if (Position[0] < MinCoordinate || Position[0] > ((double)IMax + 2*BufferLayer) ||
			Position[1] < MinCoordinate || Position[1] > ((double)JMax + 2*BufferLayer) ||
			Position[2] < MinCoordinate || Position[2] > ((double)KMax + 2*BufferLayer)) IsOk = FALSE;
			*/
		/*
		if ( Position[0] < XBeg || Position[0] > XEnd ||
			 Position[1] < YBeg || Position[1] > YEnd ||
			 Position[2] < ZBeg || Position[2] > ZEnd ) IsOk = FALSE;
			 */
		if (Position[0] < XBeg - DXZoneBuffer || Position[0] > XEnd + DXZoneBuffer ||
			Position[1] < YBeg - DYZoneBuffer || Position[1] > YEnd + DYZoneBuffer ||
			Position[2] < ZBeg - DZZoneBuffer || Position[2] > ZEnd + DZZoneBuffer ) IsOk = FALSE;

		if (IsOk)
		  IsOk =  RectGridBrickXYZtoIJKRST(Position[0], Position[1], Position[2], IMax, JMax, KMax,
										   XBeg, XEnd, YBeg, YEnd, ZBeg, ZEnd, PeriodicBC,
										   &IIndex, &JIndex, &KIndex, &r, &s, &t);


		if (IsOk)
		{
			/*
			IIndex = (LgIndex_t)Position[0];
			JIndex = (LgIndex_t)Position[1];
			KIndex = (LgIndex_t)Position[2];
			*/

			if (IsOk)
			{
				int          ic;
				LgIndex_t    Index[8];

				// Gets all the indices for the 8 vertices of the cell
				Index[0] = IndexFromIJK(IIndex    , JIndex    , KIndex    , IMax, JMax, KMax, PeriodicBC);
				Index[1] = IndexFromIJK(IIndex + 1, JIndex    , KIndex    , IMax, JMax, KMax, PeriodicBC);
				Index[2] = IndexFromIJK(IIndex + 1, JIndex + 1, KIndex    , IMax, JMax, KMax, PeriodicBC);
				Index[3] = IndexFromIJK(IIndex    , JIndex + 1, KIndex    , IMax, JMax, KMax, PeriodicBC);
				Index[4] = IndexFromIJK(IIndex    , JIndex    , KIndex + 1, IMax, JMax, KMax, PeriodicBC);
				Index[5] = IndexFromIJK(IIndex + 1, JIndex    , KIndex + 1, IMax, JMax, KMax, PeriodicBC);
				Index[6] = IndexFromIJK(IIndex + 1, JIndex + 1, KIndex + 1, IMax, JMax, KMax, PeriodicBC);
				Index[7] = IndexFromIJK(IIndex    , JIndex + 1, KIndex + 1, IMax, JMax, KMax, PeriodicBC);

				*MaxGradient = 0.0;
				*AveGradient = 0.0;
				LgIndex_t	IBig = -1;
				for (ic = 0; ic < 8; ic++)
				{
					// Finding sum of x,y,z gradient at each corner
					double CornerGrad = ABS(TecUtilDataValueGetByRef(UVarFDPtr, Index[ic])) +
										ABS(TecUtilDataValueGetByRef(VVarFDPtr, Index[ic])) +
										ABS(TecUtilDataValueGetByRef(WVarFDPtr, Index[ic]));
					if (*MaxGradient < CornerGrad)
					{
						*MaxGradient = CornerGrad;
						IBig = ic;
					}
					//*MaxGradient = MAX(*MaxGradient, CornerGrad);
					// Average gradient is max gradient /8
					//*AveGradient += 0.125 * CornerGrad;
				}
				*MaxGradient = sqrt(TecUtilDataValueGetByRef(UVarFDPtr, Index[IBig])
									* TecUtilDataValueGetByRef(UVarFDPtr, Index[IBig])
									+ TecUtilDataValueGetByRef(VVarFDPtr, Index[IBig])
									* TecUtilDataValueGetByRef(VVarFDPtr, Index[IBig])
									+ TecUtilDataValueGetByRef(WVarFDPtr, Index[IBig])
									* TecUtilDataValueGetByRef(WVarFDPtr, Index[IBig]));
				*AveGradient = *MaxGradient * 0.125;
			}
		}
	}
	else if (TecUtilZoneGetType(ZoneNum) == ZoneType_FEQuad)
	{
		// Now we're doing things according to the previous element used in grad path (I think?)
		LgIndex_t ic;
		NodeMap_pa NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);

		*MaxGradient = 0.0;
		*AveGradient = 0.0;

		// TEMP
		// if (LastElemUsed<1 || LastElemUsed > JMax)
		//   {
		//     char Message[200];
		//     sprintf_s(Message, "Elem=%d error in MaxChrgDensGrad", LastElemUsed);
		//     TecUtilDialogMessageBox(Message,  MessageBox_Warning);
		//   }

		LgIndex_t IBig = -1;
		for (ic = 0; ic < 4; ic++)
		{
			LgIndex_t Node = TecUtilDataNodeGetByRef(NodeMap, LastElemUsed, ic + 1);
			double CornerGrad = ABS(TecUtilDataValueGetByRef(UVarFDPtr, Node)) +
								ABS(TecUtilDataValueGetByRef(VVarFDPtr, Node)) +
								ABS(TecUtilDataValueGetByRef(WVarFDPtr, Node));
			if (*MaxGradient < CornerGrad)
			{
				*MaxGradient = CornerGrad;
				IBig = ic;
			}
			/**MaxGradient = MAX(*MaxGradient, CornerGrad);
			*AveGradient += 0.25 * CornerGrad;*/
		}
		LgIndex_t Node = TecUtilDataNodeGetByRef(NodeMap, LastElemUsed, IBig + 1);
		*MaxGradient = sqrt(TecUtilDataValueGetByRef(UVarFDPtr, Node)
							* TecUtilDataValueGetByRef(UVarFDPtr, Node)
							+ TecUtilDataValueGetByRef(VVarFDPtr, Node)
							* TecUtilDataValueGetByRef(VVarFDPtr, Node)
							+ TecUtilDataValueGetByRef(WVarFDPtr, Node)
							* TecUtilDataValueGetByRef(WVarFDPtr, Node));
		*AveGradient = *MaxGradient * 0.25;
	}
	else
	{
		// Zonetype wasn't one we can work with, so trigger break
		IsOk = FALSE;
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}




Boolean_t CoincidentWithCP(double     XPos,
						   double     YPos,
						   double     ZPos,
						   const CritPoints_pa CritPoints,
						   const double        CPTolerance,
						   LgIndex_t  ExceptionCrtPtNum,
						   char       MinTestCrtPtType,
						   char       MaxTestCrtPtType,
						   LgIndex_t *CritPointNum)
{
	Boolean_t IsOk = TRUE;
	Boolean_t IsCoincident = FALSE;
	LgIndex_t NumCritPoints = CritPointsGetCount(CritPoints);
	char      CrtPtType;
	LgIndex_t ii;

	REQUIRE(VALID_REF(CritPoints));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CPTolerance > 0.0);
	REQUIRE(ExceptionCrtPtNum >= -1 && ExceptionCrtPtNum < NumCritPoints);
	REQUIRE(VALID_REF(CritPointNum));


	for (ii = 0; !IsCoincident && ii < NumCritPoints; ii++)
	{
		CrtPtType = CritPointsGetType(CritPoints, ii);

		if (ii != ExceptionCrtPtNum && CrtPtType >= MinTestCrtPtType && CrtPtType <= MaxTestCrtPtType)
		{
			double XCrtPt, YCrtPt, ZCrtPt, dummy;
			double DelX, DelY, DelZ;
			char   cdummy;

			IsOk = CritPointsGetPoint(CritPoints, ii, &XCrtPt, &YCrtPt, &ZCrtPt,
									  &dummy, &cdummy, &dummy, &dummy, &dummy);

			if (IsOk)
			{
				DelX = XPos - XCrtPt;
				DelY = YPos - YCrtPt;
				DelZ = ZPos - ZCrtPt;

				/* Simple test to discard most cases. */
				if (ABS(DelX) < CPTolerance && ABS(DelY) < CPTolerance && ABS(DelZ) < CPTolerance)
				{
					if ((DelX * DelX + DelY * DelY + DelZ * DelZ) <= (CPTolerance * CPTolerance))
					{
						IsCoincident = TRUE;
						*CritPointNum = ii;
					}
				}
			}
		}
	}

	ENSURE(VALID_BOOLEAN(IsCoincident));
	ENSURE(!IsCoincident || (*CritPointNum >= 0 && *CritPointNum < NumCritPoints));
	return IsCoincident;
}


Boolean_t CoincidentWithSpecificCP(const XYZ_s	TestPt,
						   const CritPoints_pa CritPoints,
						   const double        CPTolerance,
						   char       TestCrtPtType,
						   LgIndex_t *CritPointNum)
{
	Boolean_t IsOk = TRUE;
	Boolean_t IsCoincident = FALSE;
	char      CrtPtType;
	LgIndex_t ii;

	REQUIRE(VALID_REF(CritPoints));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CPTolerance > 0.0);
	REQUIRE(VALID_REF(CritPointNum));

	LgIndex_t	BegOffset = CritPointsGetBegOffset(CritPoints, TestCrtPtType);
	LgIndex_t	EndOffset = CritPointsGetEndOffset(CritPoints, TestCrtPtType);

	for (ii = BegOffset; !IsCoincident && ii < EndOffset; ii++)
	{
		XYZ_s	CrtPtPt;
		double DelX, DelY, DelZ, djunk;
		char   cjunk;

		IsOk = CritPointsGetPoint(CritPoints, ii, &CrtPtPt.X, &CrtPtPt.Y, &CrtPtPt.Z,
			&djunk, &cjunk, &djunk, &djunk, &djunk);

		if (IsOk)
		{
			DelX = TestPt.X - CrtPtPt.X;
			DelY = TestPt.Y - CrtPtPt.Y;
			DelZ = TestPt.Z - CrtPtPt.Z;

			/* Simple test to discard most cases. */
			if (ABS(DelX) < CPTolerance && ABS(DelY) < CPTolerance && ABS(DelZ) < CPTolerance)
			{
				if (DistanceSquaredXYZ(TestPt, CrtPtPt) <= (CPTolerance * CPTolerance))
				{
					IsCoincident = TRUE;
					*CritPointNum = ii;
				}
			}
		}
	}

	ENSURE(VALID_BOOLEAN(IsCoincident));
	ENSURE(!IsCoincident || (*CritPointNum >= 0 && *CritPointNum < CritPointsGetCount(CritPoints)));
	return IsCoincident;
}





/**
 * Determine if the GradPath handle is sane.
 *
 * param GradPath
 *     GradPath structure in question.
 *
 * return
 *     TRUE if the GradPath structure is valid, otherwise FALSE.
 */
Boolean_t GradPathIsValid(GradPath_pa GradPath)
{
	Boolean_t IsValid = FALSE;

	IsValid = (VALID_REF(GradPath) &&
			   VALID_REF(GradPath->X) && ArrListIsValid(GradPath->X) &&
			   VALID_REF(GradPath->Y) && ArrListIsValid(GradPath->Y) &&
			   VALID_REF(GradPath->Z) && ArrListIsValid(GradPath->Z) &&
			   VALID_REF(GradPath->Rho) && ArrListIsValid(GradPath->Rho));

	/* Require the same count for each coordinate array. */
	if (IsValid) IsValid = (ArrListGetCount(GradPath->X) == ArrListGetCount(GradPath->Y));
	if (IsValid) IsValid = (ArrListGetCount(GradPath->X) == ArrListGetCount(GradPath->Z));
	if (IsValid) IsValid = (ArrListGetCount(GradPath->X) == ArrListGetCount(GradPath->Rho));

	ENSURE(VALID_BOOLEAN(IsValid));
	return IsValid;
}



/**
 * Deallocates the GradPath handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a GradPath handle.
 */
void GradPathDealloc(GradPath_pa *GradPath)
{
	REQUIRE(VALID_REF(GradPath));
	REQUIRE(GradPathIsValid(*GradPath) || *GradPath == NULL);

	if (*GradPath != NULL)
	{
		/* release the ArrList's */
		if ((*GradPath)->X != NULL) ArrListDealloc(&((*GradPath)->X));
		if ((*GradPath)->Y != NULL) ArrListDealloc(&((*GradPath)->Y));
		if ((*GradPath)->Z != NULL) ArrListDealloc(&((*GradPath)->Z));
		if ((*GradPath)->Rho != NULL) ArrListDealloc(&((*GradPath)->Rho));

		/* Notify Tecplot of memory usage change */
		if (ABS((*GradPath)->MemUsageReported) > 0)
		{
			TecUtilMemoryChangeNotify(- (Int64_t)((*GradPath)->MemUsageReported));
			(*GradPath)->MemUsageReported = 0;
		}

		/* release the list structure itself */
		FREE_ITEM(*GradPath, "GradPath structure");
		*GradPath = NULL;
	}

	ENSURE(*GradPath == NULL);
}




/**
 * Empties the GradPath of all points/nodes and resets the Critical Point
 * numbers to -1.
 *
 *
 * param GradPath
 *     GradPath to clear.
 */
void GradPathClear(GradPath_pa GradPath)
{
	REQUIRE(GradPathIsValid(GradPath));

	ArrListClear(GradPath->X);
	ArrListClear(GradPath->Y);
	ArrListClear(GradPath->Z);
	ArrListClear(GradPath->Rho);
	GradPath->BeginCrtPtNum = -1;
	GradPath->EndCrtPtNum   = -1;

	ENSURE(GradPathIsValid(GradPath) && GradPathGetCount(GradPath) == 0);
}


/**
 * Allocates a GradPath handle with a suitable default
 * capacity for the contained ArrList's.
 *
 *
 * return
 *     GradPath handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
GradPath_pa GradPathAlloc()
{
	GradPath_pa Result = NULL;

	Result = ALLOC_ITEM(GradPath_s, "GradPath structure");
	if (Result != NULL)
	{
		Result->MemUsageReported = 0;
		Result->BeginCrtPtNum = -1;
		Result->EndCrtPtNum   = -1;
		Result->X             = ArrListAlloc(60, ArrListType_Double);
		Result->Y             = ArrListAlloc(60, ArrListType_Double);
		Result->Z             = ArrListAlloc(60, ArrListType_Double);
		Result->Rho           = ArrListAlloc(60, ArrListType_Double);

		/* If it failed to allocate any of the array lists, clean-up and exit. */
		if (Result->X == NULL || Result->Y == NULL ||
			Result->Z == NULL)
		{
			if (Result->X != NULL) ArrListDealloc(&(Result->X));
			if (Result->Y != NULL) ArrListDealloc(&(Result->Y));
			if (Result->Z != NULL) ArrListDealloc(&(Result->Z));
			if (Result->Rho != NULL) ArrListDealloc(&(Result->Rho));
			FREE_ITEM(Result, "GradPath structure");
			Result = NULL;
		}
	}

	ENSURE(GradPathIsValid(Result) || Result == NULL);
	return Result;
}



/**
 * Gets the number of points (nodes) currently in the GradPath
 * (maintained by the GradPath coordinate arrays).
 *
 * param
 *     GradPath structure in question.
 *
 * return
 *     Number of points (nodes) in the GradPath.
 */
LgIndex_t GradPathGetCount(GradPath_pa GradPath)
{
	LgIndex_t Result = 0;

	REQUIRE(GradPathIsValid(GradPath));

	Result = ArrListGetCount(GradPath->X);

	ENSURE(Result >= 0);
	return Result;
}






/**
 * Places coordinate components at the specified offset. If the offset
 * is beyond the end of the coordinate array lists, they are sized
 * accordingly and the intervening coordinate values between the last
 * node of the original state and the last node of the new state are
 * assigned 0.
 *
 * note
 *     If coordinates already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param GradPath
 *     GradPath target in which to set the coordinates.
 * param PointOffset
 *     Offset of the point/node.
 * param X,Y,Z
 *     Coordinates to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t   GradPathSetPoint(GradPath_pa GradPath,
							 LgIndex_t   PointOffset,
							 double      X,
							 double      Y,
							 double      Z,
							 double      Rho)
{
	Boolean_t IsOk = TRUE;
	ArrListItem_u Item;


	REQUIRE(GradPathIsValid(GradPath));
	REQUIRE(PointOffset >= 0);
	// Confirming that respective arrays for x,y,z of gradpath are equal length
	REQUIRE(ArrListGetCount(GradPath->X) == ArrListGetCount(GradPath->Y));
	REQUIRE(ArrListGetCount(GradPath->X) == ArrListGetCount(GradPath->Z));

	// Item is a union type that has a property of each simple data type so 
	// it can be temporary storage of any type, but unlike a struct, a union can
	// only store ONE value at a time, where a struct can have a value assigned
	// to each data member at the same time
	Item.Double = X;
	IsOk = ArrListSetItem(GradPath->X, PointOffset, Item);

	if (IsOk)
	{
		Item.Double = Y;
		IsOk = ArrListSetItem(GradPath->Y, PointOffset, Item);
	}

	if (IsOk)
	{
		Item.Double = Z;
		IsOk = ArrListSetItem(GradPath->Z, PointOffset, Item);
	}

	if (IsOk)
	{
		Item.Double = Rho;
		IsOk = ArrListSetItem(GradPath->Rho, PointOffset, Item);
	}

	/* Require the same count for each coordinate array (checked by GradPathIsValid). */
	if (IsOk) IsOk = GradPathIsValid(GradPath);

	ENSURE(GradPathIsValid(GradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}







/**
 * Inserts coordinate components at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 * note
 *     If coordinates already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param GradPath
 *     GradPath target in which to set the coordinates.
 * param PointOffset
 *     Offset at which to insert the point/node.
 * param X,Y,Z
 *     Coordinates to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t   GradPathInsertPoint(GradPath_pa GradPath,
								LgIndex_t   PointOffset,
								double      X,
								double      Y,
								double      Z,
								double      Rho)
{
	Boolean_t IsOk = TRUE;
	ArrListItem_u Item;


	REQUIRE(GradPathIsValid(GradPath));
	REQUIRE(0 <= PointOffset && PointOffset <= GradPathGetCount(GradPath));

	Item.Double = X;
	IsOk = ArrListInsertItem(GradPath->X, PointOffset, Item);

	if (IsOk)
	{
		Item.Double = Y;
		IsOk = ArrListInsertItem(GradPath->Y, PointOffset, Item);
	}

	if (IsOk)
	{
		Item.Double = Z;
		IsOk = ArrListInsertItem(GradPath->Z, PointOffset, Item);
	}

	if (IsOk)
	{
		Item.Double = Rho;
		IsOk = ArrListInsertItem(GradPath->Rho, PointOffset, Item);
	}

	/* Require the same count for each coordinate array (checked by GradPathIsValid). */
	if (IsOk) IsOk = GradPathIsValid(GradPath);

	ENSURE(GradPathIsValid(GradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}


/**
 * Remove coordinate components at the specified offset. The members following the item
 * removed are shifted down accordingly to fill the vacated space.
 *
 * param GradPath
 *     GradPath target in which to set the coordinates.
 * param PointOffset
 *     Offset of the point/node.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t   GradPathRemovePoint(GradPath_pa GradPath,
								LgIndex_t   PointOffset)
{
	Boolean_t IsOk = TRUE;
	ArrListItem_u Item;

	REQUIRE(GradPathIsValid(GradPath));
	REQUIRE(PointOffset >= 0);

	Item = ArrListRemoveItem(GradPath->X, PointOffset);

	Item = ArrListRemoveItem(GradPath->Y, PointOffset);

	Item = ArrListRemoveItem(GradPath->Z, PointOffset);

	Item = ArrListRemoveItem(GradPath->Rho, PointOffset);


	ENSURE(GradPathIsValid(GradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}







/**
 * Inserts copies of the items from the source GradPath to the target
 * GradPath at the specified offset. The target lists will expand to
 * accommodate the additional items. The source GradPath remains
 * unchanged.
 *
 * param Target
 *     GradPath receiving the source items.
 * param PointOffset
 *     Offset at which to insert the source GradPath items.
 * param Source
 *     GradPath supplying the source items.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t   GradPathInsert(GradPath_pa Target,
						   LgIndex_t   PointOffset,
						   GradPath_pa Source)
{
	Boolean_t IsOk = TRUE;


	REQUIRE(GradPathIsValid(Target));
	REQUIRE(GradPathIsValid(Source));
	REQUIRE(0 <= PointOffset && PointOffset <= GradPathGetCount(Target));

	IsOk = ArrListInsert(Target->X, PointOffset, Source->X);

	if (IsOk)
	{
		IsOk = ArrListInsert(Target->Y, PointOffset, Source->Y);
	}

	if (IsOk)
	{
		IsOk = ArrListInsert(Target->Z, PointOffset, Source->Z);
	}

	if (IsOk)
	{
		IsOk = ArrListInsert(Target->Rho, PointOffset, Source->Rho);
	}

	/* Inform Tecplot of change in memory usage */
	if (IsOk)
	{
		/* Add estimated memory used by Source GradPath structure */
		Int64_t MemUsageChange = Source->MemUsageReported;
		if (MemUsageChange > 0)
		{
			Target->MemUsageReported += (LgIndex_t)MemUsageChange;
			TecUtilMemoryChangeNotify(MemUsageChange);
		}
	}

	ENSURE(GradPathIsValid(Target));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}








/**
 * Appends copies of the items from the source GradPath to the target
 * GradPath. The target lists will expand to accommodate the additional
 * items. The source GradPath remains unchanged.
 *
 * param Target
 *     GradPath receiving the source items.
 * param PointOffset
 *     Offset at which to insert the source GradPath items.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t   GradPathAppend(GradPath_pa Target,
						   GradPath_pa Source)
{
	Boolean_t IsOk = TRUE;

	REQUIRE(GradPathIsValid(Source));
	REQUIRE(GradPathIsValid(Target));

	IsOk = ArrListAppend(Target->X, Source->X);

	if (IsOk)
	{
		IsOk = ArrListAppend(Target->Y, Source->Y);
	}

	if (IsOk)
	{
		IsOk = ArrListAppend(Target->Z, Source->Z);
	}

	if (IsOk)
	{
		IsOk = ArrListAppend(Target->Rho, Source->Rho);
	}

	/* Inform Tecplot of change in memory usage */
	if (IsOk)
	{
		/* Add estimated memory used by Source GradPath structure */
		Int64_t MemUsageChange = Source->MemUsageReported;
		if (MemUsageChange > 0)
		{
			Target->MemUsageReported += (LgIndex_t)MemUsageChange;
			TecUtilMemoryChangeNotify(MemUsageChange);
		}
	}

	/* Change the EndCrtPtNum of Target */
	Target->EndCrtPtNum = Source->EndCrtPtNum;

	ENSURE(GradPathIsValid(Target));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}




/**
 * Appends the coordinates to the GradPath. The coordinate lists will be
 * expanded to accommodate the additional items.
 *
 * param GradPath
 *     GradPath target to which the point/node is to be appended.
 * param X,Y,Z
 *     Coordinates of node to append to the GradPath.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t GradPathAppendPoint(GradPath_pa    GradPath,
							  double         X,
							  double         Y,
							  double         Z,
							  double         Rho)
{
	Boolean_t IsOk = TRUE;
	LgIndex_t Count;

	REQUIRE(GradPathIsValid(GradPath));

	Count = GradPathGetCount(GradPath);

	IsOk = GradPathInsertPoint(GradPath, Count, X, Y, Z, Rho);

	ENSURE(GradPathIsValid(GradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}




/**
 * Gets the coordinates of the node at the specified offset in the GradPath.
 *
 * param GradPath
 *     GradPath structure containing the desired item.
 * param PointOffset
 *     Offset to the coordinates in the GradPath.
 * param *X, *Y, *Z
 *     Pointers to coordinates of the node
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */

Boolean_t   GradPathGetPoint(GradPath_pa GradPath,
							 LgIndex_t   PointOffset,
							 double     *X,
							 double     *Y,
							 double     *Z,
							 double     *Rho)
{
	Boolean_t IsOk = TRUE;
	ArrListItem_u Item;

	REQUIRE(GradPathIsValid(GradPath));
	//	TIM:	Here point offset is the number of points from the start of the grad path
	// along the grad path.
	REQUIRE(PointOffset >= 0 && PointOffset < GradPathGetCount(GradPath));
	REQUIRE(VALID_REF(X));
	REQUIRE(VALID_REF(Y));
	REQUIRE(VALID_REF(Z));
	REQUIRE(VALID_REF(Rho));

	//	TIM:	Here we're just fetching information about a point along the grad path
	// and assigning it to the x,y,z,rho passed by ref
	Item = ArrListGetItem(GradPath->X, PointOffset);
	*X = Item.Double;

	Item = ArrListGetItem(GradPath->Y, PointOffset);
	*Y = Item.Double;

	Item = ArrListGetItem(GradPath->Z, PointOffset);
	*Z = Item.Double;

	Item = ArrListGetItem(GradPath->Rho, PointOffset);
	*Rho = Item.Double;

	ENSURE(VALID_REF(X));
	ENSURE(VALID_REF(Y));
	ENSURE(VALID_REF(Z));
	ENSURE(VALID_REF(Rho));
	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}




/**
 * Get the charge density variable number from the Tecplot dataset
 * auxilliary data.
 *
 * param
 *     None.
 *
 * return
 *     ChrgDensVarNum if it works, 0 (zero) otherwise.
 */

EntIndex_t GradPathGetChrgDensVarNumFromTP()
{
	EntIndex_t ChrgDensVarNum = 0;
	EntIndex_t NumVars;
	Boolean_t  IsOk = TRUE;

	REQUIRE(TecUtilDataSetIsAvailable());

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, NULL, &NumVars);

	/* Query the DataSet aux data to find the ChrgDensVarNum */
	if (IsOk)
	{
		AuxData_pa AuxDataRef = TecUtilAuxDataDataSetGetRef();
		if (AuxDataRef != NULL)
		{
			char      *Value;
			Boolean_t  Retain;

			/* Get BegCrtPtNum */
			if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.ChrgDensVarNum",
											   &Value, &Retain))
			{
				if (Value != NULL)
				{
					/* Assign the parsed value to the structure */
					// if (sscanf_s(Value, "%d", &ChrgDensVarNum) != 1)
					//   IsOk = FALSE;
					ChrgDensVarNum = atoi(Value);

					/* release the allocated string copy */
					TecUtilStringDealloc(&Value);

					/* Test to verify that ChrgDensVarNum is within range */
					if (ChrgDensVarNum <= 0 || ChrgDensVarNum > NumVars)
						IsOk = FALSE;
				}
				else
					IsOk = FALSE;
			}
		}
	}

	if (IsOk == FALSE) ChrgDensVarNum = 0;

	ENSURE(ChrgDensVarNum >= 0);
	return ChrgDensVarNum;
}


/**
 * Find and return the TP Zone number for the GradPath with the
 * specified beginning and ending critical point numbers. Search
 * throught the Tepclot zones.
 *
 * param BegCPNum
 *     Number of the beginning critcial point number for the GradPath.
 *
 * param EndCPNum
 *     Number of the ending critcial point number for the GradPath.
 *
 * param BaseZone
 *     Number of the zone from which these topological structures were computed.
 *
 * return
 *     TP Zone number of GradPath, NULL otherwise.
 */
EntIndex_t     GradPathTPZoneFromBegEndCP(LgIndex_t     BegCrtPtNum,
										  LgIndex_t     EndCrtPtNum,
										  EntIndex_t    SourceZoneNum)
{
	EntIndex_t   TPZone = 0;
	Boolean_t    IsOk = TRUE;
	Boolean_t    IsFound = FALSE;
	EntIndex_t   NumZones, NumVars;
	LgIndex_t    BegCrtPtNumTry = -1;
	LgIndex_t    EndCrtPtNumTry = -1;
	LgIndex_t    BaseZoneNumTry = 0;

	REQUIRE(TecUtilDataSetIsAvailable());

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);


	/* Query the Auxiliary data for the zone to fill some of the structure */
	if (IsOk)
	{
		EntIndex_t ii;

		for (ii = 1; !IsFound && ii <= NumZones; ii++)
		{
			AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(ii);
			if (AuxDataRef != NULL)
			{
				char      *Value;
				Boolean_t  Retain;

				/* Get BaseZoneNum */
				if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BaseZoneNum",
												   &Value, &Retain))
				{
					if (Value != NULL)
					{
						/* Assign the parsed value to the structure */
						if (sscanf_s(Value, "%d", &BaseZoneNumTry) == 1)
						{
							IsOk = TRUE;
						}
						else
						{
							IsOk = FALSE;
						}
						// release the allocated string copy
						TecUtilStringDealloc(&Value);
					}
					else
					{
						IsOk = FALSE;
					}
				}


				if (BaseZoneNumTry == SourceZoneNum)
				{

					/* Get BegCrtPtNum */
					if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
													   &Value, &Retain))
					{
						if (Value != NULL)
						{
							/* Assign the parsed value to the structure */
							if (sscanf_s(Value, "%d", &BegCrtPtNumTry) == 1)
							{
								BegCrtPtNumTry--;
							}
							else
							{
								IsOk = FALSE;
							}
							// release the allocated string copy
							TecUtilStringDealloc(&Value);
						}
						else
						{
							IsOk = FALSE;
						}
					}

					/* Get EndCrtPtNum */
					if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
													   &Value, &Retain))
					{
						if (Value != NULL)
						{
							/* Assign the parsed value to the structure */
							if (sscanf_s(Value, "%d", &EndCrtPtNumTry) == 1)
							{
								EndCrtPtNumTry--;
							}
							else
							{
								IsOk = FALSE;
							}
							// release the allocated string copy
							TecUtilStringDealloc(&Value);
						}
						else
						{
							IsOk = FALSE;
						}
					}
				}
			}
			else  /* ZoneAuxDataRef is Null */
				IsOk = FALSE;

			if ((SourceZoneNum == BaseZoneNumTry) &&
				(BegCrtPtNum == BegCrtPtNumTry) &&
				(EndCrtPtNum == EndCrtPtNumTry))
			{
				TPZone = ii;
				IsFound = TRUE;
			}
		}
	}

	ENSURE(TPZone >= 0 && TPZone <= NumZones);
	return TPZone;
}






/**
 * Find and return the GradPath with the specified beginning and ending
 * critical point numbers. Search throught the Tepclot zones.
 *
 * param BegCPNum
 *     Number of the beginning critcial point number for the GradPath.
 *
 * param EndCPNum
 *     Number of the ending critcial point number for the GradPath.
 *
 * return
 *     Handle to GradPath if it works, NULL otherwise.
 */

GradPath_pa   GradPathGetByBegEndCP(LgIndex_t BegCrtPtNum,
									LgIndex_t EndCrtPtNum)
{
	GradPath_pa  Result = NULL;
	Boolean_t    IsFound = FALSE;
	Boolean_t    IsOk = TRUE;
	EntIndex_t   NumZones, NumVars, ZoneNum;
	EntIndex_t   ChrgDensVarNum = 0;
	LgIndex_t    IMax, JMax, KMax;
	LgIndex_t    BegCrtPtNumTry, EndCrtPtNumTry;
	FieldData_pa XVarFDPtr = NULL;
	FieldData_pa YVarFDPtr = NULL;
	FieldData_pa ZVarFDPtr = NULL;
	FieldData_pa CDVarFDPtr = NULL;
	ArrListItem_u Item;

	REQUIRE(TecUtilDataSetIsAvailable());

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	Result = GradPathAlloc();

	/* Query the Auxiliary data for the zone to fill some of the structure */
	if (IsOk)
	{
		EntIndex_t ii;

		for (ii = 1; !IsFound && ii <= NumZones; ii++)
		{
			AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(ii);
			if (AuxDataRef != NULL)
			{
				char      *Value;
				Boolean_t  Retain;

				/* Get BegCrtPtNum */
				if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
												   &Value, &Retain))
				{
					if (Value != NULL)
					{
						/* Assign the parsed value to the structure */
						if (sscanf_s(Value, "%d", &BegCrtPtNumTry) == 1)
						{
							BegCrtPtNumTry--;
						}
						else
						{
							IsOk = FALSE;
						}
						// release the allocated string copy
						TecUtilStringDealloc(&Value);
					}
					else
					{
						IsOk = FALSE;
					}
				}

				/* Get EndCrtPtNum */
				if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
												   &Value, &Retain))
				{
					if (Value != NULL)
					{
						/* Assign the parsed value to the structure */
						if (sscanf_s(Value, "%d", &EndCrtPtNumTry) == 1)
						{
							EndCrtPtNumTry--;
						}
						else
						{
							IsOk = FALSE;
						}
						// release the allocated string copy
						TecUtilStringDealloc(&Value);
					}
					else
					{
						IsOk = FALSE;
					}
				}
			}
			else  /* ZoneAuxDataRef is Null */
				IsOk = FALSE;

			if ((BegCrtPtNum == BegCrtPtNumTry) &&
				(EndCrtPtNum == EndCrtPtNumTry))
			{
				ZoneNum = ii;
				IsFound = TRUE;
			}
		}
	}

	/* Create the GradPath */
	if (IsFound == TRUE)
	{
		TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
						   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

		ChrgDensVarNum = GradPathGetChrgDensVarNumFromTP();
		CDVarFDPtr = TecUtilDataValueGetRef(ZoneNum, ChrgDensVarNum);

		/* Extract coordinate values and add to coordinate ArrLists */
		if (IsOk)
		{
			LgIndex_t ii;
			double X, Y, Z, Rho;

			for (ii = 1; IsOk && ii <= IMax; ii++)
			{
				X = TecUtilDataValueGetByRef(XVarFDPtr, ii);
				Y = TecUtilDataValueGetByRef(YVarFDPtr, ii);
				Z = TecUtilDataValueGetByRef(ZVarFDPtr, ii);
				Rho = TecUtilDataValueGetByRef(CDVarFDPtr, ii);

				Item.Double = X;
				IsOk = ArrListAppendItem(Result->X, Item);

				Item.Double = Y;
				IsOk = ArrListAppendItem(Result->Y, Item);

				Item.Double = Z;
				IsOk = ArrListAppendItem(Result->Z, Item);

				Item.Double = Rho;
				IsOk = ArrListAppendItem(Result->Rho, Item);
			}
		}

		/* Set the BeginCrtPtNum and EndCrtPtNum structure elements */
		Result->BeginCrtPtNum = BegCrtPtNum;
		Result->EndCrtPtNum   = EndCrtPtNum;

		if (IsOk == FALSE) GradPathDealloc(&Result);
	}
	else  // IsFound == FALSE
		GradPathDealloc(&Result);

	/* Inform Tecplot of memory usage */
	if (IsOk && IsFound)
	{
		/* Estimate of memory used by GradPath structured, 3 ArrList structures, and 3 double arrays */
		LgIndex_t MemUsageEstimate = (96 + 4 * 8 * (IMax + 1)) / 1024;
		if (MemUsageEstimate > 0 && Result->MemUsageReported < MemUsageEstimate)
		{
			Int64_t MemUsageChange = MemUsageEstimate - Result->MemUsageReported;
			Result->MemUsageReported = MemUsageEstimate;
			TecUtilMemoryChangeNotify(MemUsageChange);
		}
	}


	ENSURE(Result == NULL || GradPathIsValid(Result));
	return Result;
}







/**
 * Generate a GradPath from a GradPath that has been converted
 * into a Tecplot Zone.
 *
 * param TPZoneNum
 *     Tecplot Zone number from which to generate the GradPath.
 *
 * return
 *     Handle to GradPath if it works, NULL otherwise.
 */

GradPath_pa   GradPathGetFromTPZone(EntIndex_t ZoneNum)
{
	GradPath_pa  Result = NULL;
	Boolean_t    IsOk = TRUE;
	EntIndex_t   NumZones, NumVars;
	EntIndex_t   ChrgDensVarNum = 0;
	LgIndex_t    IMax, JMax, KMax;
	FieldData_pa XVarFDPtr = NULL;
	FieldData_pa YVarFDPtr = NULL;
	FieldData_pa ZVarFDPtr = NULL;
	FieldData_pa CDVarFDPtr = NULL;
	ArrListItem_u Item;

	REQUIRE(TecUtilDataSetIsAvailable());

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);
	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));

	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(IMax > 1 && JMax == 1 && KMax == 1); /* Must be I-ordered Line */

	REQUIRE(VALID_REF(XVarFDPtr));
	REQUIRE(VALID_REF(YVarFDPtr));
	REQUIRE(VALID_REF(ZVarFDPtr));

	Result = GradPathAlloc();

	/* Query the Auxiliary data for the zone to fill some of the structure */
	if (IsOk)
	{
		AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(ZoneNum);
		if (AuxDataRef != NULL)
		{
			char      *Value;
			Boolean_t  Retain;

			/* Get BegCrtPtNum */
			if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
											   &Value, &Retain))
			{
				if (Value != NULL)
				{
					/* Assign the parsed value to the structure */
					long BegCrtPtNum;
					if (sscanf_s(Value, "%d", &BegCrtPtNum) == 1)
					{
						Result->BeginCrtPtNum = (LgIndex_t)(BegCrtPtNum - 1);
					}
					else
					{
						IsOk = FALSE;
					}
					// release the allocated string copy
					TecUtilStringDealloc(&Value);
				}
				else
				{
					IsOk = FALSE;
				}
			}

			/* Get EndCrtPtNum */
			if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
											   &Value, &Retain))
			{
				if (Value != NULL)
				{
					/* Assign the parsed value to the structure */
					long EndCrtPtNum;
					if (sscanf_s(Value, "%d", &EndCrtPtNum) == 1)
					{
						Result->EndCrtPtNum = (LgIndex_t)(EndCrtPtNum - 1);
					}
					else
					{
						IsOk = FALSE;
					}
					// release the allocated string copy
					TecUtilStringDealloc(&Value);
				}
				else
				{
					IsOk = FALSE;
				}
			}
		}
		else  /* ZoneAuxDataRef is Null */
			IsOk = FALSE;
	}

	/* Get the ChrgDensVarNum from Tecplot DataSet Aux data */
	if (IsOk)
	{
		ChrgDensVarNum = GradPathGetChrgDensVarNumFromTP();
		if (ChrgDensVarNum == 0) IsOk = FALSE;
		if (IsOk) CDVarFDPtr = TecUtilDataValueGetRef(ZoneNum, ChrgDensVarNum);
	}

	/* Extract coordinate values and add to coordinate ArrLists */
	if (IsOk)
	{
		LgIndex_t ii;
		double X, Y, Z, Rho;

		for (ii = 1; IsOk && ii <= IMax; ii++)
		{
			X = TecUtilDataValueGetByRef(XVarFDPtr, ii);
			Y = TecUtilDataValueGetByRef(YVarFDPtr, ii);
			Z = TecUtilDataValueGetByRef(ZVarFDPtr, ii);
			Rho = TecUtilDataValueGetByRef(CDVarFDPtr, ii);

			Item.Double = X;
			IsOk = ArrListAppendItem(Result->X, Item);

			Item.Double = Y;
			IsOk = ArrListAppendItem(Result->Y, Item);

			Item.Double = Z;
			IsOk = ArrListAppendItem(Result->Z, Item);

			Item.Double = Rho;
			IsOk = ArrListAppendItem(Result->Rho, Item);
		}

		if (IsOk == FALSE) GradPathDealloc(&Result);
	}

	/* Inform Tecplot of memory usage */
	if (IsOk)
	{
		/* Estimate of memory used by GradPath structured, 3 ArrList structures, and 3 double arrays */
		LgIndex_t MemUsageEstimate = (96 + 4 * 8 * (IMax + 1)) / 1024;
		if (MemUsageEstimate > 0 && Result->MemUsageReported < MemUsageEstimate)
		{
			Int64_t MemUsageChange = MemUsageEstimate - Result->MemUsageReported;
			Result->MemUsageReported = MemUsageEstimate;
			TecUtilMemoryChangeNotify(MemUsageChange);
		}
	}


	ENSURE(GradPathIsValid(Result));
	return Result;
}





/**
 * Return the length of the GradPath start to finish.
 *
 * param GradPath
 *     GradPath structure containing the desired item.
 *
 * return
 *     Length if it works, 0.0 otherwise.
 */
double GradPathGetLength(GradPath_pa GradPath)
{
	Boolean_t IsOk = TRUE;
	double    Length = 0.0;
	LgIndex_t NumPoints = GradPathGetCount(GradPath);

	REQUIRE(GradPathIsValid(GradPath));

	if (NumPoints > 1)
	{
		LgIndex_t ii;
		double XIm1, YIm1, ZIm1, RhoIm1;

		IsOk = GradPathGetPoint(GradPath, 0, &XIm1, &YIm1, &ZIm1, &RhoIm1);

		for (ii = 1; IsOk && ii < NumPoints; ii++)
		{
			double X, Y, Z, Rho;
			IsOk = GradPathGetPoint(GradPath, ii, &X, &Y, &Z, &Rho);

			if (IsOk)
			{
				double DelX = X - XIm1;
				double DelY = Y - YIm1;
				double DelZ = Z - ZIm1;
				Length += sqrt(DelX * DelX + DelY * DelY + DelZ * DelZ);

				XIm1 = X;
				YIm1 = Y;
				ZIm1 = Z;
				RhoIm1 = Rho;
			}
		}
	}

	ENSURE(Length >= 0.0);
	return Length;
}







/**
 * Return the distance, at closest approach, between the GradPath
 * and a point in space.
 *
 * param GradPath
 *     GradPath structure containing the desired item.
 * param Point
 *     XYZ_s structure containing the coordinates of the point in space.
 *
 * return
 *     MinDistance if it works, -1.0 otherwise.
 */
double GradPathGetDistToPt(GradPath_pa GradPath,
						   XYZ_s       Point)
{
	Boolean_t IsOk = TRUE;
	double    MinDistance = LARGEFLOAT;

	REQUIRE(GradPathIsValid(GradPath));

	LgIndex_t NumPoints = GradPathGetCount(GradPath);


	if (NumPoints > 1)
	{
		LgIndex_t ii;
		double XIm1, YIm1, ZIm1, RhoIm1;
		double XPt, YPt, ZPt;
		double DXPtIm1, DYPtIm1, DZPtIm1;
		double DXPtI,   DYPtI,   DZPtI;
		double DistToPtFromIm1,    DistToPtFromI;
		double DistToPtFromIm1Sqr, DistToPtFromISqr;

		XPt = Point.X;
		YPt = Point.Y;
		ZPt = Point.Z;

		IsOk = GradPathGetPoint(GradPath, 0, &XIm1, &YIm1, &ZIm1, &RhoIm1);

		DXPtIm1 = XPt - XIm1;
		DYPtIm1 = YPt - YIm1;
		DZPtIm1 = ZPt - ZIm1;

		DistToPtFromIm1Sqr = DXPtIm1 * DXPtIm1 + DYPtIm1 * DYPtIm1 + DZPtIm1 * DZPtIm1;
		DistToPtFromIm1 = sqrt(DistToPtFromIm1Sqr);

		for (ii = 1; IsOk && ii < NumPoints; ii++)
		{
			double XI, YI, ZI, RhoI;
			IsOk = GradPathGetPoint(GradPath, ii, &XI, &YI, &ZI, &RhoI);

			if (IsOk)
			{
				double DXSeg = XI - XIm1;
				double DYSeg = YI - YIm1;
				double DZSeg = ZI - ZIm1;
				double LengthSeg = sqrt(DXSeg * DXSeg + DYSeg * DYSeg + DZSeg * DZSeg);

				double NormProjPtToLine = (DXPtIm1 * DXSeg + DYPtIm1 * DYSeg + DZPtIm1 * DZSeg) / LengthSeg;

				DXPtI = XPt - XI;
				DYPtI = YPt - YI;
				DZPtI = ZPt - ZI;

				DistToPtFromISqr = DXPtI * DXPtI + DYPtI * DYPtI + DZPtI * DZPtI;
				DistToPtFromI = sqrt(DistToPtFromISqr);

				// If the normal projection of the point to the line is outside the length
				// of the line, use the distance from the appropriate end point. If it is
				// within the length of the line, use the Pythagorean theorem.
				if (NormProjPtToLine <= 0.0)
				{
					MinDistance = MIN(MinDistance, DistToPtFromIm1);
				}
				else if (NormProjPtToLine >= LengthSeg)
				{
					MinDistance = MIN(MinDistance, DistToPtFromI);
				}
				else
				{
					double NormalDist = sqrt(DistToPtFromIm1Sqr - NormProjPtToLine * NormProjPtToLine);

					MinDistance = MIN(MinDistance, NormalDist);
				}

				// Incrementing I-Index. Reused previously computed I values for I-1
				XIm1 = XI;
				YIm1 = YI;
				ZIm1 = ZI;
				RhoIm1 = RhoI;

				DXPtIm1 = DXPtI;
				DYPtIm1 = DYPtI;
				DZPtIm1 = DZPtI;

				DistToPtFromIm1Sqr = DistToPtFromISqr;
				DistToPtFromIm1    = DistToPtFromI;
			}
		}
	}

	if (IsOk == FALSE)
		MinDistance = -1.0;

	ENSURE(MinDistance >= 0.0 || MinDistance == -1.0);
	return MinDistance;
}







/**
 * Modify the last point of a GradPath so that it is on the boundary
 * of the zone. Do this by interpolating/extrapolating the last line
 * segment.
 *
 * NOTE: Currently assumes IJK-block volume zone.
 *
 * param GradPath
 *     GradPath structure to be clipped.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t GradPathClipToSrcZoneDomain(ZoneVarInfo_pa ZoneVarInfo,
									  GradPath_pa    GradPath)
{
	Boolean_t  IsOk = TRUE;
	EntIndex_t NumZones, NumVars, ZoneNum;
	LgIndex_t  IMax, JMax, KMax;
	LgIndex_t  NumPoints;
	double XLm1 , YLm1 , ZLm1 , RhoLm1 ;
	double XL   , YL   , ZL   , RhoL   ;
	double XLNew, YLNew, ZLNew, RhoLNew;

	REQUIRE(ZoneVarInfo != NULL);

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);
	ZoneNum = ZoneVarInfo->ZoneNum;

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);

	REQUIRE(GradPathIsValid(GradPath));

	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	NumPoints = GradPathGetCount(GradPath);
	if (NumPoints < 2) IsOk = FALSE;

	if (IsOk)
		IsOk = GradPathGetPoint(GradPath, NumPoints - 2, &XLm1, &YLm1, &ZLm1, &RhoLm1);
	if (IsOk)
		IsOk = GradPathGetPoint(GradPath, NumPoints - 1, &XL  , &YL  , &ZL  , &RhoL);

	// Extrapolate/Interpolate the last segment to the closest intersection with a
	// boundary of the volume zone. NOTE: this assumes an IJK Block zone.
	if (IsOk)
	{
		double t = LARGEFLOAT;

		// t goes from zero at the Last-1 point, to one at the Last point. Find t
		// where it leaves the IJK block zone.
		if (XL > XLm1)
		{
			// t = MIN(t, ((double)IMax - XLm1) / MAX(XL - XLm1, 1.0e-6));
			t = MIN(t, (ZoneVarInfo->XEnd - XLm1) / MAX(XL - XLm1, 1.0e-6));
		}
		else
		{
			// t = MIN(t, (1.0 - XLm1) / MIN(XL - XLm1, -1.0e-6));
			t = MIN(t, (ZoneVarInfo->XBeg - XLm1) / MIN(XL - XLm1, -1.0e-6));
		}

		if (YL > YLm1)
		{
			// t = MIN(t, ((double)JMax - YLm1) / MAX(YL - YLm1, 1.0e-6));
			t = MIN(t, (ZoneVarInfo->YEnd - YLm1) / MAX(YL - YLm1, 1.0e-6));
		}
		else
		{
			// t = MIN(t, (1.0 - YLm1) / MIN(YL - YLm1, -1.0e-6));
			t = MIN(t, (ZoneVarInfo->YBeg - YLm1) / MIN(YL - YLm1, -1.0e-6));
		}

		if (ZL > ZLm1)
		{
			// t = MIN(t, ((double)KMax - ZLm1) / MAX(ZL - ZLm1, 1.0e-6));
			t = MIN(t, (ZoneVarInfo->ZEnd - ZLm1) / MAX(ZL - ZLm1, 1.0e-6));
		}
		else
		{
			// t = MIN(t, (1.0 - ZLm1) / MIN(ZL - ZLm1, -1.0e-6));
			t = MIN(t, (ZoneVarInfo->ZBeg - ZLm1) / MIN(ZL - ZLm1, -1.0e-6));
		}

		XLNew   = (1.0 - t) * XLm1   + t * XL;
		YLNew   = (1.0 - t) * YLm1   + t * YL;
		ZLNew   = (1.0 - t) * ZLm1   + t * ZL;
		RhoLNew = (1.0 - t) * RhoLm1 + t * RhoL;

		IsOk = GradPathSetPoint(GradPath, NumPoints - 1, XLNew, YLNew, ZLNew, RhoLNew);
	}

	ENSURE(GradPathIsValid(GradPath));
	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}





/**
 * Return the approximate maximum separation of the two GradPaths.
 * Separation is approximated as the normal distance between the
 * mid-points (in terms of number of points) of the two lines.
 *
 * param GradPath1
 *     GradPath structure containing the first of the two GradPaths
 *     to be compared.
 *
 * param GradPath2
 *     GradPath structure containing the second of the two GradPaths
 *     to be compared.
 *
 * return
 *     Length if it works, 0.0 otherwise.
 */
double GradPathGetSeparation(GradPath_pa GradPath1, GradPath_pa GradPath2)
{
	Boolean_t IsOk = TRUE;
	double    Separation = 0.0;
	LgIndex_t NumPoints1 = GradPathGetCount(GradPath1);
	LgIndex_t NumPoints2 = GradPathGetCount(GradPath2);

	REQUIRE(GradPathIsValid(GradPath1));
	REQUIRE(GradPathIsValid(GradPath2));

	if (NumPoints1 > 1 && NumPoints2 > 1)
	{
		double X1, Y1, Z1, X1Im1, Y1Im1, Z1Im1, Rho1;
		double X2, Y2, Z2, X2Im1, Y2Im1, Z2Im1, Rho2;
		double dx, dy, dz, tx, ty, tz, dtan;

		IsOk = GradPathGetPoint(GradPath1, NumPoints1 / 2, &X1, &Y1, &Z1, &Rho1);
		IsOk = GradPathGetPoint(GradPath2, NumPoints2 / 2, &X2, &Y2, &Z2, &Rho2);

		IsOk = GradPathGetPoint(GradPath1, NumPoints1 / 2 - 1, &X1Im1, &Y1Im1, &Z1Im1, &Rho1);
		IsOk = GradPathGetPoint(GradPath2, NumPoints2 / 2 - 1, &X2Im1, &Y2Im1, &Z2Im1, &Rho2);

		if (IsOk)
		{
			/* Compute GradPath tangent vector - averaged between the two GradPaths */
			tx = 0.5 * (X1 - X1Im1 + X2 - X2Im1);
			ty = 0.5 * (Y1 - Y1Im1 + Y2 - Y2Im1);
			tz = 0.5 * (Z1 - Z1Im1 + Z2 - Z2Im1);

			/* Compute the distance vector between the two points */
			dx = X2 - X1;
			dy = Y2 - Y1;
			dz = Z2 - Z1;

			/* Normal component of distance is computed from distance vector and tangent
			 * component of distance vector using Pythagorean theorem.
			 */
			dtan = dx * tx + dy * ty + dz * tz;
			Separation = sqrt(dx * dx + dy * dy + dz * dz - dtan * dtan);
		}
	}

	ENSURE(Separation >= 0.0);
	return Separation;
}








/**
 * Return the mean gradient along the GradPath = start to finish.
 *
 * param GradPath
 *     GradPath structure containing the desired item.
 *
 * return
 *     Mean gradient if it works, 0.0 otherwise.
 */
double GradPathGetMeanGrad(GradPath_pa GradPath)
{
	Boolean_t IsOk = TRUE;
	double    Length = 0.0;
	double    RhoChng = 0.0;
	double    Grad = 0.0;
	LgIndex_t NumPoints = GradPathGetCount(GradPath);

	REQUIRE(GradPathIsValid(GradPath));

	if (NumPoints > 1)
	{
		LgIndex_t ii;
		double XIm1, YIm1, ZIm1, RhoIm1;

		IsOk = GradPathGetPoint(GradPath, 0, &XIm1, &YIm1, &ZIm1, &RhoIm1);

		for (ii = 1; IsOk && ii < NumPoints; ii++)
		{
			double X, Y, Z, Rho;
			IsOk = GradPathGetPoint(GradPath, ii, &X, &Y, &Z, &Rho);

			if (IsOk && X > 1.0 && Y > 1.0 && Z > 1.0)
			{
				double DelX = X - XIm1;
				double DelY = Y - YIm1;
				double DelZ = Z - ZIm1;
				double DelRho = Rho - RhoIm1;
				double DelLength = sqrt(DelX * DelX + DelY * DelY + DelZ * DelZ);
				// double DelGrad = DelRho / ( MAX(DelLength, SMALLFLOAT) );

				/* This just gives the sum of the gradient of each piece. */
				Length += DelLength;
				RhoChng += DelRho;
				// Grad += DelGrad;

				XIm1 = X;
				YIm1 = Y;
				ZIm1 = Z;
				RhoIm1 = Rho;
			}
		}

		Grad = RhoChng / (MAX(Length, SMALLFLOAT));
	}

	return Grad;
}





/**
 * Resample the GradPath so that there are NPointsNew equally spaced
 * grid points along its length.
 *
 * param GradPath
 *     GradPath structure containing the desired item.
 * param NumPointsNew
 *     Pointers to coordinates of the node
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t GradPathResample(GradPath_pa *GradPath,
						   LgIndex_t    NumPointsNew,
						   double       segmentRatio)
{
	Boolean_t IsOk = TRUE;

	GradPath_pa GradPathNew  = NULL;
	LgIndex_t   NumPointsOld = GradPathGetCount(*GradPath);
	double      Length       = GradPathGetLength(*GradPath);

	REQUIRE(GradPathIsValid(*GradPath));
	REQUIRE(NumPointsNew > 1);

	if (NumPointsOld > 1)
	{
		// Count the total length of all segments. Note: there is one less
		// segment than there are points
		double curSegmentLength = 1.0;
		double totalLength = 0.0;
		for (int ii = 0; ii < NumPointsNew-1; ii++)
		{
			totalLength += curSegmentLength;
			curSegmentLength *= segmentRatio;
		}
		double delLength = Length / totalLength;

		LgIndex_t iold = 0;
		LgIndex_t inew;
		double    ArcLength = 0.0;
		double    XIm1, YIm1, ZIm1, RhoIm1;
		double    XI, YI, ZI, RhoI;
		double    ArcLengthIm1 = 0.0;
		double    ArcLengthI   = 0.0;

		GradPathNew = GradPathAlloc();

		/* The begin and end critical points are the same as before resampling. */
		GradPathNew->BeginCrtPtNum = (*GradPath)->BeginCrtPtNum;
		GradPathNew->EndCrtPtNum   = (*GradPath)->EndCrtPtNum;

		/* Set first point */
		IsOk = GradPathGetPoint(*GradPath, 0, &XI, &YI, &ZI, &RhoI);

		IsOk = GradPathAppendPoint(GradPathNew, XI, YI, ZI, RhoI);

		/* Initialize XIm1, YIm1, ZIm1 in case first ArcLength is <zero. */
		XIm1 = XI;
		YIm1 = YI;
		ZIm1 = ZI;
		RhoIm1 = RhoI;

		for (inew = 1; IsOk && inew < NumPointsNew - 1; inew++)
		{
			ArcLength += delLength;

			/*
			 * Increment index along old GradPath until find appropriate
			 * segment for interpolation.
			 */
			while (iold < NumPointsOld - 1 && ArcLengthI < ArcLength)
			{
				iold++;

				ArcLengthIm1 = ArcLengthI;
				XIm1 = XI;
				YIm1 = YI;
				ZIm1 = ZI;
				RhoIm1 = RhoI;

				IsOk = GradPathGetPoint(*GradPath, iold, &XI, &YI, &ZI, &RhoI);

				if (IsOk)
				{
					double DelX = XI - XIm1;
					double DelY = YI - YIm1;
					double DelZ = ZI - ZIm1;

					ArcLengthI += sqrt(DelX * DelX + DelY * DelY + DelZ * DelZ);
				}
			}

			/* Interpolate to get new coordinates and add to new GradPath */
			if (IsOk)
			{
				double Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
				double XNew = XIm1 + Ratio * (XI - XIm1);
				double YNew = YIm1 + Ratio * (YI - YIm1);
				double ZNew = ZIm1 + Ratio * (ZI - ZIm1);
				double RhoNew = RhoIm1 + Ratio * (RhoI - RhoIm1);
				IsOk = GradPathAppendPoint(GradPathNew, XNew, YNew, ZNew, RhoNew);
			}

			delLength *= segmentRatio;
		}

		/* Add Last Point (same as old last point) */
		IsOk = GradPathGetPoint(*GradPath, NumPointsOld - 1, &XI, &YI, &ZI, &RhoI);
		IsOk = GradPathAppendPoint(GradPathNew, XI, YI, ZI, RhoI);

		/* Replace GradPath */
		if (IsOk)
		{
			GradPathDealloc(GradPath);
			*GradPath = GradPathNew;
		}

		/* Inform Tecplot of memory usage */
		if (IsOk)
		{
			/* Estimate of memory used by GradPath structured, 3 ArrList structures, and 3 double arrays */
			LgIndex_t MemUsageEstimate = (96 + 4 * 8 * (NumPointsNew + 1)) / 1024;
			if (MemUsageEstimate > 0 && (*GradPath)->MemUsageReported < MemUsageEstimate)
			{
				Int64_t MemUsageChange = MemUsageEstimate - (*GradPath)->MemUsageReported;
				(*GradPath)->MemUsageReported = MemUsageEstimate;
				TecUtilMemoryChangeNotify(MemUsageChange);
			}
		}
	}

	ENSURE(GradPathGetCount(*GradPath) == NumPointsNew);
	ENSURE(GradPathIsValid(*GradPath));
	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}







/**
 * Reverse the order of points in GradPath.
 *
 * param GradPath
 *     GradPath structure containing the desired item.
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t GradPathReverse(GradPath_pa *GradPath)
{
	Boolean_t IsOk = TRUE;

	GradPath_pa GradPathNew  = NULL;
	LgIndex_t   NumPoints = GradPathGetCount(*GradPath);

	REQUIRE(VALID_REF(GradPath));
	REQUIRE(GradPathIsValid(*GradPath));

	if (NumPoints >= 1)
	{
		LgIndex_t iold = 0;
		LgIndex_t inew;
		double    XI, YI, ZI, RhoI;

		GradPathNew = GradPathAlloc();

		/* The begin and end critical points are also reversed. */
		GradPathNew->BeginCrtPtNum = (*GradPath)->EndCrtPtNum;
		GradPathNew->EndCrtPtNum   = (*GradPath)->BeginCrtPtNum;

		for (inew = 0; IsOk && inew < NumPoints; inew++)
		{
			iold = NumPoints - inew - 1;
			/*
			 * Increment index along old GradPath until find appropriate
			 * segment for interpolation.
			 */

			IsOk = GradPathGetPoint(*GradPath, iold, &XI, &YI, &ZI, &RhoI);

			/* Add to new GradPath */
			if (IsOk)
			{
				IsOk = GradPathAppendPoint(GradPathNew, XI, YI, ZI, RhoI);
			}
		}

		/* Replace GradPath */
		if (IsOk)
		{
			GradPathDealloc(GradPath);
			*GradPath = GradPathNew;
		}

		/* No change in memory usage, no need to inform Tecplot of memory usage */
	}

	ENSURE(GradPathGetCount(*GradPath) == NumPoints);
	ENSURE(GradPathIsValid(*GradPath));
	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}








/*
 * GradPath: Compute the gradient path (streamline). Start it at the given
 *  position and terminate it at a critical point or the boundary of the
 *  domain. Do this by solving the system of 3 ordinary differential equations:
 *
 *   d(XYZ)/dt = Gradient(XYZ)
 *
 */
Boolean_t GradPathAdd(const ZoneVarInfo_pa  ZoneVarInfo,
					  const LgIndex_t       MaxNumPathPoints,
					  const StreamDir_e     PathDir,
					  const CritPoints_pa   CritPoints,
					  const double          CPTolerance,
					  LgIndex_t      *NumPathPoints,
					  GradPath_pa     GradPath)
{
	Boolean_t     IsOk = TRUE;
	Boolean_t     IsFound = TRUE;
	Boolean_t     IsCritPoint = FALSE;
	Boolean_t     ReverseDir;
	double        DelTime;
	char          MinTestCrtPtType = (char)(-3);   // Will only test Gradpath proximity to CrtPts w/ Type >= this
	char          MaxTestCrtPtType = (char)(13);    // Will only test Gradpath proximity to CrtPts w/ Type <= this
	char          BegCrtPtType;


	double    b2  = 0.75;
	LgIndex_t ii = 0;

	float    GridSpacing = 1.0;

	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(VALID_REF(NumPathPoints));
	REQUIRE(*NumPathPoints >= 0);
	REQUIRE(MaxNumPathPoints > 0);
	REQUIRE(VALID_ENUM(PathDir, StreamDir_e));
	REQUIRE(VALID_REF(CritPoints));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CPTolerance > 0.0);
	REQUIRE(VALID_REF(GradPath));
	REQUIRE(GradPathIsValid(GradPath));
	REQUIRE((GradPath->BeginCrtPtNum >= 0 &&
			 GradPath->BeginCrtPtNum < CritPointsGetCount(CritPoints)) ||
			GradPath->BeginCrtPtNum == -1);

	/* Initialize EndCrtPtNum (-1 means doesn't end at Critical Point */
	GradPath->EndCrtPtNum = -1;

	/* Set ReverseDir based on PathDir */
	ReverseDir = FALSE;
	if (PathDir == StreamDir_Reverse)
		ReverseDir = TRUE;

	// Ideally you know the origin and destination CP types, so you can do the next line
	/* Determine range of critical point type to test against for GradPath termination */
	if (GradPath->BeginCrtPtNum > -1)
	{
		BegCrtPtType = CritPointsGetType(CritPoints, GradPath->BeginCrtPtNum);
		if (ReverseDir)
			MinTestCrtPtType = BegCrtPtType;
		else
			MaxTestCrtPtType = BegCrtPtType;
	}

	/* Compute initial time step based on local value of gradient. Want
	 * pathline to cross 1/4 of cell per step
	 */
	if (IsOk)
	{
		double Solution[3]; // Stores x,y,z solutions to diff eq
		double Rho;
		double MaxGradient, AveGradient;

		GridSpacing = (float)ZoneVarInfoGetDXFromZoneNum(ZoneVarInfo, ZoneVarInfo->ZoneNum);

		IsOk = GradPathGetPoint(GradPath, 0, &(Solution[0]), &(Solution[1]), &(Solution[2]), &Rho);

		// Here, if a new max charge density is found then base the initial time step on that,
		// otherwise base initial time step on the grid spacing and max value of gradient.
		// TODO: Make 1.0e-9 a defined constant.
		if (MaxChrgDensGrad((void *)ZoneVarInfo, Solution, &MaxGradient, &AveGradient))
			DelTime = 0.25 * GridSpacing / MAX(MaxGradient, 1.0e-9);
		// DelTime = 0.5/AveGradient;
		else
			DelTime = GridSpacing;
	}

	// This is where the grad path actually gets built.
	/* March streamline (gradient path) forward in time. */
	while (IsFound && *NumPathPoints < MaxNumPathPoints - 1 && !IsCritPoint)
	{
		double TmpSol[3], TmpFunc1[3], TmpFunc2[3], Solution[3], ChrgDens;
		Boolean_t TimeStepReduced = FALSE;

		// Get coordinate and density info about the new point.
		IsOk = GradPathGetPoint(GradPath, ii, &(Solution[0]), &(Solution[1]), &(Solution[2]), &ChrgDens);

		/* Update the position by one step */
		IsFound = RK2UpdateSolution(3, 0.75, &DelTime, ReverseDir, (void *)ZoneVarInfo,
			ChrgDensGrad, TmpSol, TmpFunc1, TmpFunc2, Solution,
			&TimeStepReduced);

		// This part only allows for a single time step reduction, when many might be necessary.
		// Need to loop over the above line until TimeStepReduced remains false.

		/* If time step is found to be too large, redo update with smaller time step. */
		if (TimeStepReduced)
		{
			IsOk = GradPathGetPoint(GradPath, ii, &(Solution[0]), &(Solution[1]), &(Solution[2]), &ChrgDens);

			IsFound = RK2UpdateSolution(3, 0.75, &DelTime, ReverseDir, (void *)ZoneVarInfo,
				ChrgDensGrad, TmpSol, TmpFunc1, TmpFunc2, Solution,
				&TimeStepReduced);
		}

		/* Update time step for next step. */
		if (IsOk)
		{
			double MaxGradient, AveGradient;

			GridSpacing = (float)ZoneVarInfoGetDXFromZoneNum(ZoneVarInfo, ZoneVarInfo->ZoneNum);

			// TODO: Make 1.0e-9 a defined constant.
			if (MaxChrgDensGrad((void *)ZoneVarInfo, Solution, &MaxGradient, &AveGradient))
				DelTime = 0.25 * GridSpacing / MAX(MaxGradient, 1.0e-9);
			// DelTime = 0.5/AveGradient;
			else
				DelTime = GridSpacing;

			// DelTime = 0.75*DelTime + 0.25*DelTimeNew;
			DelTime = DelTime;
		}

		/* Stop if !IsFound (left zone) or at critical point */
		if (IsFound)
		{
			LgIndex_t CritPointNum;
			double    ChrgDens;
			LgIndex_t NumInternalCPs = CritPoints->NumCrtPts - CritPoints->NumFFCrtPtsP1 - CritPoints->NumFFCrtPtsP3;

			ii++;
			*NumPathPoints = ii;

			ChrgDens = ChrgDensTrilinearInterp(ZoneVarInfo, Solution);

			IsOk = GradPathAppendPoint(GradPath, Solution[0], Solution[1], Solution[2], ChrgDens);

			IsCritPoint = CoincidentWithCP(Solution[0], Solution[1], Solution[2],
										   CritPoints, CPTolerance * GridSpacing, GradPath->BeginCrtPtNum,
										   MinTestCrtPtType, MaxTestCrtPtType, &CritPointNum);
			/* Set last pathline position to end critical point position */
			if (IsCritPoint)
				// if (IsCritPoint && CritPointNum < NumInternalCPs)
			{
				double XCrtPt, YCrtPt, ZCrtPt, dummy, SolTmp[3];
				char   cdummy;

				IsOk = CritPointsGetPoint(CritPoints, CritPointNum,
										  &XCrtPt, &YCrtPt, &ZCrtPt,
										  &dummy, &cdummy, &dummy, &dummy, &dummy);

				SolTmp[0] = XCrtPt;
				SolTmp[1] = YCrtPt;
				SolTmp[2] = ZCrtPt;
				// TODO: Temporary
				// ChrgDens = ChrgDensTrilinearInterp(ZoneVarInfo, SolTmp);

				IsOk = GradPathSetPoint(GradPath, ii, XCrtPt, YCrtPt, ZCrtPt, ChrgDens);

				GradPath->EndCrtPtNum = CritPointNum;
			}
		}

		/* If it has passed outside the zone, shorten the last line segment so
		 * that last point is on the boundary of the zone
		 */
		if (IsOk && (!IsFound || *NumPathPoints == (MaxNumPathPoints - 1)))
		{
			// IsOk = TRUE;
			IsOk = GradPathClipToSrcZoneDomain(ZoneVarInfo, GradPath);
		}

		/* Set first pathline position to beginning critical point position */
		if (IsOk && ii == 1 && GradPath->BeginCrtPtNum >= 0)
		{
			double XCrtPt, YCrtPt, ZCrtPt, dummy, SolTmp[3];
			char   cdummy;

			IsOk = CritPointsGetPoint(CritPoints, GradPath->BeginCrtPtNum,
									  &XCrtPt, &YCrtPt, &ZCrtPt, &dummy, &cdummy, &dummy, &dummy, &dummy);

			SolTmp[0] = XCrtPt;
			SolTmp[1] = YCrtPt;
			SolTmp[2] = ZCrtPt;
			ChrgDens = ChrgDensTrilinearInterp(ZoneVarInfo, SolTmp);

			IsOk = GradPathSetPoint(GradPath, 0, XCrtPt, YCrtPt, ZCrtPt, ChrgDens);
		}
	}

	/* Inform Tecplot of memory usage */
	if (IsOk)
	{
		/* Estimate of memory used by GradPath structured, 3 ArrList structures, and 3 double arrays */
		LgIndex_t MemUsageEstimate = (96 + 4 * 8 * (ii + 1)) / 1024;
		if (MemUsageEstimate > 0 && GradPath->MemUsageReported < MemUsageEstimate)
		{
			Int64_t MemUsageChange = MemUsageEstimate - GradPath->MemUsageReported;
			GradPath->MemUsageReported = MemUsageEstimate;
			TecUtilMemoryChangeNotify(MemUsageChange);
		}
	}

	*NumPathPoints = ii + 1;

	ENSURE(GradPathGetCount(GradPath) == *NumPathPoints);
	ENSURE(VALID_BOOLEAN(IsOk));
	return(IsOk);

} /* GradPathAdd() */








/*
 * GradPathAddMidway: Compute the gradient path (streamline). Start it at
 *  the given, integrate it both directions, and terminate it at a critical
 *  points or the boundary of the domain.
 *
 * param ZoneVarInfo
 *    Zone and Variable numbers needed by the function.
 * param MaxNumPathPoints
 *    Maximum number of path points. Integration is terminated it this is
 *    exceeded.
 * param CritPoints
 *    Structure describing the critical points.
 * param SeedPos
 *    XYZ_s position of the seed point.
 * param CPTolerance
 *    Distance from the critical point that is considered at the critical
 *    point. Needed since gradients are zero at the critical points.
 * param NumPathPoint
 *    Pointer to the output of the computed number of path points.
 * param GradPath
 *    Pointer to computed Gradient Path structure.
 *
 * return
 *    TRUE if it works, FALSE if it fails
 */
GradPath_pa GradPathAddMidway(const ZoneVarInfo_pa  ZoneVarInfo,
							  const LgIndex_t       MaxNumPathPoints,
							  const CritPoints_pa   CritPoints,
							  const XYZ_s           SeedPos,
							  const double          CPTolerance,
							  LgIndex_t            *NumPathPoints)
{
	Boolean_t     IsOk = TRUE;
	GradPath_pa GradPathForward  = GradPathAlloc();
	GradPath_pa GradPathBackward = GradPathAlloc();
	double XPos = SeedPos.X;
	double YPos = SeedPos.Y;
	double ZPos = SeedPos.Z;
	double junk = 0.0;

	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(VALID_REF(NumPathPoints));
	REQUIRE(*NumPathPoints >= 0);
	REQUIRE(MaxNumPathPoints > 0);
	REQUIRE(VALID_REF(CritPoints));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CPTolerance > 0.0);

	if (GradPathForward == NULL) IsOk = FALSE;
	if (GradPathBackward == NULL) IsOk = FALSE;

	// Clear GradPaths
	if (IsOk)
	{
		GradPathClear(GradPathForward);
		GradPathClear(GradPathBackward);
	}

	/* Seed first point */
	if (IsOk)
		IsOk = GradPathSetPoint(GradPathForward, 0, XPos, YPos, ZPos, junk);
	if (IsOk)
		IsOk = GradPathSetPoint(GradPathBackward, 0, XPos, YPos, ZPos, junk);

	/* Integrate forward gradient path line */
	if (IsOk)
	{
		GradPathForward->BeginCrtPtNum = -1;

		*NumPathPoints = 0;
		IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Forward,
						   CritPoints, 1.0, NumPathPoints, GradPathForward);
	}

	/* Integrate backward gradient path line */
	if (IsOk)
	{
		GradPathBackward->BeginCrtPtNum = -1;

		*NumPathPoints = 0;
		IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Reverse,
						   CritPoints, 1.0, NumPathPoints, GradPathBackward);
	}

	// Reverse the first GradPath and append the two
	if (IsOk)
	{
		IsOk = GradPathReverse(&GradPathForward);
	}
	if (IsOk)
	{
		LgIndex_t NumPoints = GradPathGetCount(GradPathForward);
		IsOk = GradPathRemovePoint(GradPathForward, NumPoints - 1);
	}
	if (IsOk)
		IsOk = GradPathAppend(GradPathForward, GradPathBackward);

	GradPathDealloc(&GradPathBackward);

	if (IsOk)
		*NumPathPoints = GradPathGetCount(GradPathForward);
	else
	{
		GradPathDealloc(&GradPathForward);
		GradPathForward = NULL;
	}

	ENSURE(GradPathForward == NULL || GradPathGetCount(GradPathForward) == *NumPathPoints);
	ENSURE(GradPathForward == NULL || GradPathIsValid(GradPathForward));
	return(GradPathForward);

} /* GradPathAddMidway() */











/*
 * GradPathAddSurf: Compute the gradient path (streamline). Start it at the given
 *  position and terminate it at a critical point or the boundary of the
 *  domain. Do this by solving the system of 3 ordinary differential equations:
 *
 *   d(XYZ)/dt = Gradient(XYZ)
 *
 *  and projecting back down to the surface on each step.
 *
 * param ZoneVarInfo
 *     ZoneVarInfo structure containing the surface zone number and gradient
 *     variable information.
 * param MaxNumPathPoints
 *     Terminate the GradPath after MaxNumGradPoints points are added.
 * param PathDir
 *     StreamDir_Forward (with the gradient), StreamDir_Reverse (against the gradient),
 *     or StreamDir_Both
 * param CritPoints
 *     Struture containing location, type, and orientation information about the
 *     critical points.
 * param CPToleranceOriginal
 *     Distance from a CP that is considered at the CP.
 *
 * returned through pointers
 *     NumPathPoints and GradPath structure.
 *
 * return
 *     TRUE if successful, FALSE otherwise.
 *
 */
Boolean_t GradPathAddSurf(const ZoneVarInfo_pa  ZoneVarInfo,
						  const LgIndex_t       MaxNumPathPoints,
						  const StreamDir_e     PathDir,
						  const CritPoints_pa   CritPoints,
						  const double          CPToleranceOriginal,
						  LgIndex_t            *NumPathPoints,
						  GradPath_pa           GradPath)
{
	Boolean_t		IsOk = TRUE;
	Boolean_t		IsFound = TRUE;
	Boolean_t		IsCritPoint = FALSE;
	Boolean_t		ReverseDir;
	double			DelTime;
	LgIndex_t		ElemLast = 1;
	
	//	Maximum ratio of grad path length to straight line distance of beg-end point
	
	XYZ_s			BeginXYZ, EndXYZ;

	char          MinTestCrtPtType = (char)(-2);   // Will only test Gradpath proximity to CrtPts w/ Type >= this
	char          MaxTestCrtPtType = (char)(2);    // Will only test Gradpath proximity to CrtPts w/ Type <= this
	char          BegCrtPtType;

	double    b2  = 0.75;
	LgIndex_t ii = 0;

	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(VALID_REF(NumPathPoints));
	REQUIRE(*NumPathPoints >= 0);
	REQUIRE(MaxNumPathPoints > 0);
	REQUIRE(VALID_ENUM(PathDir, StreamDir_e));
	REQUIRE(VALID_REF(CritPoints));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CPToleranceOriginal > 0.0);
	REQUIRE(VALID_REF(GradPath));
	REQUIRE(GradPathIsValid(GradPath));
	REQUIRE((GradPath->BeginCrtPtNum >= 0 &&
			 GradPath->BeginCrtPtNum < CritPointsGetCount(CritPoints)) ||
			GradPath->BeginCrtPtNum == -1);

	double			StepSizeFactor = 0.05 * CritPoints->MinCPDistance,
					StepSizeIncrement = 0.05 * CritPoints->MinCPDistance,
					StepSizeFactor2ndTry = 0.01 * CritPoints->MinCPDistance;
	int				StepSizeResizeMax = 10, StepSizeResizeMax2ndTry = 15;
	double			LengthRatioTolerance = 1.1;
	//if (CritPoints->NumCrtPts <= 4) LengthRatioTolerance = 3.141592653589793 / 2.0;
	/* Initialize EndCrtPtNum (-1 means doesn't end at Critical Point */
	GradPath->EndCrtPtNum = -1;

	/* Set ReverseDir based on PathDir */
	ReverseDir = FALSE;
	if (PathDir == StreamDir_Reverse)
		ReverseDir = TRUE;

  //  /* Determine range of critical point type to test against for GradPath termination */
  //  if (GradPath->BeginCrtPtNum > -1)
  //  {
  //      BegCrtPtType = CritPointsGetType(CritPoints, GradPath->BeginCrtPtNum);
  //      //	TIM:	all surf grad paths are seeded from saddles (0), so this is
		////			unnecessary flexibility.
		////			Changed so that only 2 and -2 are looked for, respectively
		//if (ReverseDir)
  //          MinTestCrtPtType = 2;
  //      else
  //          MaxTestCrtPtType = -2;
  //  }

	

	//	Outer loop to adjust stepsize in grad path to prevent
	//	erratic grad path behavior
	for (int StepSizeMult = 0 ; IsFound && !IsCritPoint && StepSizeMult < StepSizeResizeMax ; StepSizeMult++)
	{
		/* Compute initial time step based on local value of gradient. Want
		 * pathline to cross 1/4 of cell per step
		 */
		if (IsOk)
		{
			double Solution[3];
			double Rho;
			double MaxGradient, AveGradient;
			double XJunk, YJunk, ZJunk;

			IsOk = GradPathGetPoint(GradPath, 0, &(Solution[0]), &(Solution[1]), &(Solution[2]), &Rho);

			// Get nearest element number
			ElemLast = GeomToolsClosestElemNum(ZoneVarInfo->ZoneNum, Solution[0], Solution[1], Solution[2], &XJunk, &YJunk, &ZJunk);
			ZoneVarInfo->LastElemUsed = ElemLast;

			//// TEMP
			//// if (ElemLast == 0)
			////   {
			////     char Message[200];
			////     sprintf_s(Message, "ElemLast==0 at start of GradPathAddSurf");
			////     TecUtilDialogMessageBox(Message, MessageBox_Warning);
			////   }

			//// TODO: Make 1.0e-9 a defined constant.
			////	TIM: TODO: start the step size based on average length of 
			////				isosurface cell for surface.
			//if (MaxChrgDensGrad((void *)ZoneVarInfo, Solution, &MaxGradient, &AveGradient))
			//	DelTime = 0.25 / MAX(MaxGradient, 1.0e-9) * 
			//		(StepSizeFactor + (double)StepSizeMult * StepSizeIncrement);
			//// DelTime = 0.5/AveGradient;
			//else
			//	DelTime = 1.0;

			double GradGradVect[3], GradGradMag = 0.0;
			if (SurfTotGradGrad((void *)ZoneVarInfo, Solution, GradGradVect))
			{
				for (int j = 0 ; j < 3 ; j++)
					GradGradMag += GradGradVect[j] * GradGradVect[j];

				GradGradMag = sqrt(GradGradMag);

				DelTime = (StepSizeFactor + (double)StepSizeMult * StepSizeIncrement) / GradGradMag;
			}
			else DelTime = 1.0;
		}


		/* March streamline (gradient path) forward in time. */
		int			CPTolResizes = 0;
		const int	CPTolResizesMax = 5;
		double		CPTolerance = CPToleranceOriginal;
		ii = 0;
		while (IsFound && CPTolResizes < CPTolResizesMax - 1 && !IsCritPoint)
		{
			if (ii > 0 || *NumPathPoints > 0)
			{
				CPTolResizes++;
				CPTolerance *= (2.0 + CPTolResizes) / (1.0 + CPTolResizes);
				ii = 0;
				*NumPathPoints = ii;
				for (LgIndex_t jj = GradPathGetCount(GradPath) - 1 ; jj > 0 ; jj--)
					GradPathRemovePoint(GradPath, jj);
			}

			while (IsFound && *NumPathPoints < (MaxNumPathPoints / 10) - 1 && !IsCritPoint)
			{
				double TmpSol[3], TmpFunc1[3], TmpFunc2[3], Solution[3], ChrgDens;
				Boolean_t TimeStepReduced = FALSE;

				IsOk = GradPathGetPoint(GradPath, ii, &(Solution[0]), &(Solution[1]), &(Solution[2]), &ChrgDens);

				/* Update the position by one step */
				IsFound = RK2UpdateSolution(3, 0.75, &DelTime, ReverseDir, (void *)ZoneVarInfo,
					SurfTotGradGrad, TmpSol, TmpFunc1, TmpFunc2, Solution,
					&TimeStepReduced);

				/* If time step is found to be too large, redo update with smaller time step. */
				/*if (TimeStepReduced)
				{
					IsOk = GradPathGetPoint(GradPath, ii, &(Solution[0]), &(Solution[1]), &(Solution[2]), &ChrgDens);

					IsFound = RK2UpdateSolution(3, 0.75, &DelTime, ReverseDir, (void *)ZoneVarInfo,
						SurfTotGradGrad, TmpSol, TmpFunc1, TmpFunc2, Solution,
						&TimeStepReduced);
				}*/

				// Project back to surface
				// TODO
				
				XYZ_s SurfPt;

				if (IsOk)
				{
					// TODO IsOk = GradPathProjectPointToSurf(ZoneVarInfo, ii, GradPath, &ElemLast);
					IsOk = GradPathProjectPointToSurf(ZoneVarInfo, GradPath, ii, Solution, &SurfPt);

					if (IsOk)
					{
						Solution[0] = SurfPt.X;
						Solution[1] = SurfPt.Y;
						Solution[2] = SurfPt.Z;
					}
				}

				/* Update time step for next step. */
				if (IsOk)
				{
					//double MaxGradient, AveGradient;

					//// TODO: Make 1.0e-9 a defined constant.
					//if (MaxChrgDensGrad((void *)ZoneVarInfo, Solution, &MaxGradient, &AveGradient))
					//	DelTime = 0.25 / MAX(MaxGradient, 1.0e-9) * 
					//		(StepSizeFactor + (double)StepSizeMult * StepSizeIncrement);
					//// DelTime = 0.5/AveGradient;
					//else
					//	DelTime = 1.0;

					//// DelTime = 0.75*DelTime + 0.25*DelTimeNew;
					//DelTime = DelTime;

					double GradGradVect[3], GradGradMag = 0.0;
					if (SurfTotGradGrad((void *)ZoneVarInfo, Solution, GradGradVect))
					{
						for (int j = 0 ; j < 3 ; j++)
							GradGradMag += GradGradVect[j] * GradGradVect[j];

						GradGradMag = sqrt(GradGradMag);

						DelTime = (StepSizeFactor + (double)StepSizeMult * StepSizeIncrement) / GradGradMag;
					}
					else DelTime = 1.0;
				}

				/* Stop if !IsFound (left zone) or at critical point */
				if (IsFound)
				{
					LgIndex_t CritPointNum;
					double    ChrgDens = 0;

					ii++;
					*NumPathPoints = ii;

					//  ChrgDens = ChrgDensSurfInterp(ZoneVarInfo, Solution, ElemLast);

					IsOk = GradPathAppendPoint(GradPath, Solution[0], Solution[1], Solution[2], ChrgDens);

					IsCritPoint = CoincidentWithCP(Solution[0], Solution[1], Solution[2],
												   CritPoints, CPTolerance, GradPath->BeginCrtPtNum,
												   MinTestCrtPtType, MaxTestCrtPtType, &CritPointNum);
					
					/* Set last pathline position to end critical point position */
					if (IsCritPoint)
					{
						double XCrtPt, YCrtPt, ZCrtPt, dummy, SolTmp[3];
						char   cdummy;

						IsOk = CritPointsGetPoint(CritPoints, CritPointNum,
												  &XCrtPt, &YCrtPt, &ZCrtPt,
												  &dummy, &cdummy, &dummy, &dummy, &dummy);

						SolTmp[0] = XCrtPt;
						SolTmp[1] = YCrtPt;
						SolTmp[2] = ZCrtPt;

						EndXYZ.X = XCrtPt;
						EndXYZ.Y = YCrtPt;
						EndXYZ.Z = ZCrtPt;
						// TODO ChrgDens = ChrgDensSurfInterp(ZoneVarInfo, SolTmp, ElemLast);

						IsOk = GradPathSetPoint(GradPath, ii, XCrtPt, YCrtPt, ZCrtPt, ChrgDens);

						GradPath->EndCrtPtNum = CritPointNum;
					}
				}

				/* Set first pathline position to beginning critical point position */
				if (IsOk && ii == 1 && GradPath->BeginCrtPtNum >= 0)
				{
					double XCrtPt, YCrtPt, ZCrtPt, dummy, SolTmp[3];
					char   cdummy;

					IsOk = CritPointsGetPoint(CritPoints, GradPath->BeginCrtPtNum,
											  &XCrtPt, &YCrtPt, &ZCrtPt, &dummy, &cdummy, &dummy, &dummy, &dummy);

					SolTmp[0] = XCrtPt;
					SolTmp[1] = YCrtPt;
					SolTmp[2] = ZCrtPt;

					BeginXYZ.X = XCrtPt;
					BeginXYZ.Y = YCrtPt;
					BeginXYZ.Z = ZCrtPt;

					// TODO ChrgDens = ChrgDensSurfInterp(ZoneVarInfo, SolTmp, ElemLast);

					IsOk = GradPathSetPoint(GradPath, 0, XCrtPt, YCrtPt, ZCrtPt, ChrgDens);
				}

			}
		}
		if (IsFound && IsCritPoint)
		{
			if (GradPathGetLength(GradPath) / 
				SurfaceLength(ZoneVarInfo->Radius, sqrt(DistanceSquaredXYZ(EndXYZ, BeginXYZ))) 
					> LengthRatioTolerance)
			{
				for (LgIndex_t jj = GradPathGetCount(GradPath) - 1 ; jj > 0 ; jj--)
					GradPathRemovePoint(GradPath, jj);
				IsCritPoint = FALSE;
			}
		}
		if (!IsCritPoint && StepSizeMult == StepSizeResizeMax - 1)
		{
			StepSizeFactor = StepSizeFactor2ndTry;
			StepSizeResizeMax = StepSizeResizeMax2ndTry;
			StepSizeMult = -1;
		}
	}

	/* Inform Tecplot of memory usage */
	if (IsOk)
	{
		/* Estimate of memory used by GradPath structured, 3 ArrList structures, and 3 double arrays */
		LgIndex_t MemUsageEstimate = (96 + 4 * 8 * (ii + 1)) / 1024;
		if (MemUsageEstimate > 0 && GradPath->MemUsageReported < MemUsageEstimate)
		{
			Int64_t MemUsageChange = MemUsageEstimate - GradPath->MemUsageReported;
			GradPath->MemUsageReported = MemUsageEstimate;
			TecUtilMemoryChangeNotify(MemUsageChange);
		}
	}

	*NumPathPoints = ii + 1;

	ENSURE(GradPathGetCount(GradPath) == *NumPathPoints);
	ENSURE(VALID_BOOLEAN(IsOk));
	return(IsOk);

} /* GradPathAddSurf() */









/*
 * GradPathSurfAddMidway: Compute the gradient path (streamline). Start it at
 *  the given, integrate it both directions, and terminate it at a critical
 *  points or the boundary of the domain.
 *
 * param ZoneVarInfo
 *    Zone and Variable numbers needed by the function.
 * param MaxNumPathPoints
 *    Maximum number of path points. Integration is terminated it this is
 *    exceeded.
 * param CritPoints
 *    Structure describing the critical points.
 * param SeedPos
 *    XYZ_s position of the seed point.
 * param CPTolerance
 *    Distance from the critical point that is considered at the critical
 *    point. Needed since gradients are zero at the critical points.
 * param NumPathPoint
 *    Pointer to the output of the computed number of path points.
 * param GradPath
 *    Pointer to computed Gradient Path structure.
 *
 * return
 *    TRUE if it works, FALSE if it fails
 */
GradPath_pa GradPathSurfAddMidway(const ZoneVarInfo_pa  ZoneVarInfo,
								  const LgIndex_t       MaxNumPathPoints,
								  const CritPoints_pa   CritPoints,
								  const XYZ_s           SeedPos,
								  const double          CPTolerance,
								  LgIndex_t            *NumPathPoints)
{
	Boolean_t     IsOk = TRUE;
	GradPath_pa GradPathForward  = GradPathAlloc();
	GradPath_pa GradPathBackward = GradPathAlloc();
	double XPos = SeedPos.X;
	double YPos = SeedPos.Y;
	double ZPos = SeedPos.Z;
	double junk = 0.0;

	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(VALID_REF(NumPathPoints));
	REQUIRE(*NumPathPoints >= 0);
	REQUIRE(MaxNumPathPoints > 0);
	REQUIRE(VALID_REF(CritPoints));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CPTolerance > 0.0);

	if (GradPathForward == NULL) IsOk = FALSE;
	if (GradPathBackward == NULL) IsOk = FALSE;

	// Clear GradPaths
	if (IsOk)
	{
		GradPathClear(GradPathForward);
		GradPathClear(GradPathBackward);
	}

	/* Seed first point */
	if (IsOk)
		IsOk = GradPathSetPoint(GradPathForward, 0, XPos, YPos, ZPos, junk);
	if (IsOk)
		IsOk = GradPathSetPoint(GradPathBackward, 0, XPos, YPos, ZPos, junk);

	/* Integrate forward gradient path line */
	if (IsOk)
	{
		GradPathForward->BeginCrtPtNum = -1;

		ZoneVarInfo->PathDir = StreamDir_Forward;

		*NumPathPoints = 0;
		IsOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Forward,
							   CritPoints, 1.0, NumPathPoints, GradPathForward);
	}

	/* Integrate backward gradient path line */
	if (IsOk)
	{
		GradPathBackward->BeginCrtPtNum = -1;

		ZoneVarInfo->PathDir = StreamDir_Reverse;

		*NumPathPoints = 0;
		IsOk = GradPathAddSurf(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Reverse,
							   CritPoints, 1.0, NumPathPoints, GradPathBackward);
	}

	// Reverse the first GradPath and append the two
	if (IsOk)
	{
		IsOk = GradPathReverse(&GradPathForward);
	}
	if (IsOk)
	{
		LgIndex_t NumPoints = GradPathGetCount(GradPathForward);
		IsOk = GradPathRemovePoint(GradPathForward, NumPoints - 1);
	}
	if (IsOk)
		IsOk = GradPathAppend(GradPathForward, GradPathBackward);

	GradPathDealloc(&GradPathBackward);

	if (IsOk)
		*NumPathPoints = GradPathGetCount(GradPathForward);
	else
	{
		GradPathDealloc(&GradPathForward);
		GradPathForward = NULL;
	}

	ENSURE(GradPathForward == NULL || GradPathGetCount(GradPathForward) == *NumPathPoints);
	ENSURE(GradPathForward == NULL || GradPathIsValid(GradPathForward));
	return(GradPathForward);

} /* GradPathSurfAddMidway() */





/*
 * GradPathAdd2PtMinLen: Compute the minimum length gradient path
 *  (streamline) between two critical points. Start with a straight
 *  line between the points. Solve a system of three equations for the
 *  position of each node along the line, using a combination of a force
 *  pushing the line parallel to the velocity, and an elastic force
 *  pulling the line straight.
 *
 * parameter ZoneVarInfo:
 *   Zone and variable numbers.
 * parameter NumPathPoints:
 *   The number of path points desired.
 * parameter CritPoints
 *   Structure containing the information about critical points
 * parameter GradPath
 *   Gradient path structure to be populated by this function
 *
 *   d(XYZ)/dt = Gradient(XYZ)
 *
 * Returns
 *   True if it works, False if it fails.
 */

Boolean_t GradPathAdd2PtMinLen(const ZoneVarInfo_pa  ZoneVarInfo,
							   const LgIndex_t       NumPathPoints,
							   const CritPoints_pa   CritPoints,
							   GradPath_pa     GradPath)
{
	Boolean_t     IsOk = TRUE;
	// double        DelTime;

	double    b2  = 0.75;
	LgIndex_t ii = 0;

	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(NumPathPoints > 0);
	REQUIRE(VALID_REF(CritPoints));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(VALID_REF(GradPath));
	REQUIRE(GradPathIsValid(GradPath));
	REQUIRE(GradPath->BeginCrtPtNum >= 0 &&
			GradPath->BeginCrtPtNum < CritPointsGetCount(CritPoints));
	REQUIRE(GradPath->BeginCrtPtNum >= 0 &&
			GradPath->BeginCrtPtNum < CritPointsGetCount(CritPoints));

	/* Start with an equally spaced straight line between the beginning
	 * and ending critical points.
	 */
	if (IsOk)
	{
		double    XCrtBeg = 0.0;
		double    YCrtBeg = 0.0;
		double    ZCrtBeg = 0.0;
		double    XCrtEnd = 0.0;
		double    YCrtEnd = 0.0;
		double    ZCrtEnd = 0.0;
		double    RNumPts = 1.0 / ((double)(NumPathPoints - 1));
		double    DelX, DelY, DelZ, Length;
		double    junk;
		char      cjunk;
		LgIndex_t ii;

		IsOk = CritPointsGetPoint(CritPoints, GradPath->BeginCrtPtNum,
								  &XCrtBeg, &YCrtBeg, &ZCrtBeg, &junk,
								  &cjunk, &junk, &junk, &junk);
		IsOk = CritPointsGetPoint(CritPoints, GradPath->EndCrtPtNum,
								  &XCrtEnd, &YCrtEnd, &ZCrtEnd, &junk,
								  &cjunk, &junk, &junk, &junk);

		DelX = XCrtEnd - XCrtBeg;
		DelY = YCrtEnd - YCrtBeg;
		DelZ = ZCrtEnd - ZCrtBeg;
		Length = sqrt(DelX * DelX + DelY * DelY + DelZ * DelZ);

		for (ii = 0; IsOk && ii < NumPathPoints; ii++)
		{
			double X, Y, Z;
			double Rho = 0.0;
			X = XCrtBeg + ((double)ii) * DelX * RNumPts;
			Y = YCrtBeg + ((double)ii) * DelY * RNumPts;
			Z = ZCrtBeg + ((double)ii) * DelZ * RNumPts;
			IsOk = GradPathAppendPoint(GradPath, X, Y, Z, Rho);
		}
	}

	/*
	 * Iterate to move each point according to opposing forces - velocity
	 * and elastic stretching.
	 */
	if (IsOk)
	{
		int iter;
		for (iter = 0; IsOk && iter < 100; iter++)
		{
			double MaxVelContr, DXVel, DYVel, DZVel;
			double Residual = 0.0;
			double GradTotMax = 0.0;
			LgIndex_t ii;
			for (ii = 1; IsOk && ii < NumPathPoints - 1; ii++)
			{
				double X = 0.0;
				double Y = 0.0;
				double Z = 0.0;
				double Rho = 0.0;
				double XIm1 = 0.0;
				double YIm1 = 0.0;
				double ZIm1 = 0.0;
				double RhoIm1 = 0.0;
				double XIp1 = 0.0;
				double YIp1 = 0.0;
				double ZIp1 = 0.0;
				double RhoIp1 = 0.0;
				double XAve, YAve, ZAve, DXAve, DYAve, DZAve;
				double Position[3], GradIm[3], GradIp[3], GradNormIm[3], GradNormIp[3];
				double GradTot, GTRatio;
				double DelXYZ[3];

				IsOk = GradPathGetPoint(GradPath, ii, &X, &Y, &Z, &Rho);
				IsOk = GradPathGetPoint(GradPath, ii - 1, &XIm1, &YIm1, &ZIm1, &RhoIm1);
				IsOk = GradPathGetPoint(GradPath, ii + 1, &XIp1, &YIp1, &ZIp1, &RhoIp1);

				/* Elasticity - ii moves toward average of ii-1 and ii+1 */
				XAve = 0.5 * (XIm1 + XIp1);
				YAve = 0.5 * (YIm1 + YIp1);
				ZAve = 0.5 * (ZIm1 + ZIp1);
				DXAve = XAve - X;
				DYAve = YAve - Y;
				DZAve = ZAve - Z;

				/* Velocity - ii moves to be parallel to velocity */
				Position[0] = 0.5 * (X + XIm1);
				Position[1] = 0.5 * (Y + YIm1);
				Position[2] = 0.5 * (Z + ZIm1);
				IsOk = ChrgDensGrad(ZoneVarInfo, Position, GradIm);
				GradTot = sqrt(GradIm[0] * GradIm[0] + GradIm[1] * GradIm[1] + GradIm[2] * GradIm[2]);
				GradTotMax = MAX(GradTotMax, GradTot);
				GradIm[0] = GradIm[0] / GradTot;
				GradIm[1] = GradIm[1] / GradTot;
				GradIm[2] = GradIm[2] / GradTot;

				Position[0] = 0.5 * (X + XIp1);
				Position[1] = 0.5 * (Y + YIp1);
				Position[2] = 0.5 * (Z + ZIp1);
				IsOk = ChrgDensGrad(ZoneVarInfo, Position, GradIp);
				GradTot = sqrt(GradIp[0] * GradIp[0] + GradIp[1] * GradIp[1] + GradIp[2] * GradIp[2]);
				GradTotMax = MAX(GradTotMax, GradTot);
				GradIp[0] = GradIp[0] / GradTot;
				GradIp[1] = GradIp[1] / GradTot;
				GradIp[2] = GradIp[2] / GradTot;

				if (IsOk)
				{
					int jj;
					double GradTang = 0.0;
					double Length = 0.0;

					DelXYZ[0] = X - XIm1;
					DelXYZ[1] = Y - YIm1;
					DelXYZ[2] = Z - ZIm1;

					for (jj = 0; IsOk && jj < 3; jj++)
					{
						GradTang += DelXYZ[jj] * GradIm[jj];
						Length  += DelXYZ[jj] * DelXYZ[jj];
					}
					Length = sqrt(Length);
					GradTang = GradTang / Length;

					for (jj = 0; IsOk && jj < 3; jj++)
					{
						DelXYZ[jj] = DelXYZ[jj] / Length;
						GradNormIm[jj] = GradIm[jj] - GradTang * DelXYZ[jj];
						GradNormIm[jj] = GradNormIm[jj] * Length;
					}

					DelXYZ[0] = XIp1 - X;
					DelXYZ[1] = YIp1 - Y;
					DelXYZ[2] = ZIp1 - Z;
					GradTang = 0.0;
					Length = 0.0;

					for (jj = 0; IsOk && jj < 3; jj++)
					{
						GradTang += DelXYZ[jj] * GradIp[jj];
						Length  += DelXYZ[jj] * DelXYZ[jj];
					}
					Length = sqrt(Length);
					GradTang = GradTang / Length;
					for (jj = 0; IsOk && jj < 3; jj++)
					{
						DelXYZ[jj] = DelXYZ[jj] / Length;
						GradNormIp[jj] = GradIp[jj] - GradTang * DelXYZ[jj];
						GradNormIp[jj] = GradNormIp[jj] * Length;
					}
				}

				GTRatio = MIN((GradTot * 10.0 / GradTotMax), 1.0);
				DXVel = GTRatio * (GradNormIm[0] + GradNormIp[0]);
				DYVel = GTRatio * (GradNormIm[1] + GradNormIp[1]);
				DZVel = GTRatio * (GradNormIm[2] + GradNormIp[2]);
				// DXVel = GradNormIm[0];
				// DYVel = GradNormIm[1];
				// DZVel = GradNormIm[2];
				MaxVelContr = 0.1 * MAX((MAX(ABS(DelXYZ[0]), ABS(DelXYZ[1]))), ABS(DelXYZ[2]));

				// X += DXAve;
				// Y += DYAve;
				// Z += DZAve;
				X += 0.03 * DXVel + DXAve;
				Y += 0.03 * DYVel + DYAve;
				Z += 0.03 * DZVel + DZAve;
				IsOk = GradPathSetPoint(GradPath, ii, X, Y, Z, 0.0001);

				Residual += DXVel * DXVel + DYVel * DYVel + DZVel * DZVel +
							DXAve * DXAve + DYAve * DYAve + DZAve * DZAve;

			}
		}
	}


	ENSURE(GradPathGetCount(GradPath) == NumPathPoints);
	ENSURE(VALID_BOOLEAN(IsOk));
	return(IsOk);
}





/**
 * Test the GradPath module. Create a simple grad path and try some
 * utility funcitons.
 *
 * TODO: do the same thing for linear surface.
 *
 *
 * return
 *     TRUE if it worked. FALSE otherwise.
 */
Boolean_t GradPathTest()
{
	Boolean_t IsOk = TRUE;
	GradPath_pa GradPath = NULL;

	// Allocate the gradient path structure
	GradPath = GradPathAlloc();
	if (GradPath == NULL) IsOk = FALSE;

	CHECK(GradPathIsValid(GradPath));

	//  Simple Test: 4 point line
	//
	//                  ^ y = Coord2
	//                  |
	//                  |
	//                  |  sloping line
	//        (0,0.5,1) |___________ (1, 1,2)
	//                 /|
	//               /  |
	//             /    |
	//           /      |
	//(-1,0,0) /---------------------------> x = Coord1
	//        |         |z-coord out
	//        |         |
	//        |         |
	//        |         |
	//   (-1,-1,-1)
	//
	// Set first point
	if (IsOk) IsOk = GradPathSetPoint(GradPath, 0, -1.0, -1.0, -1.0, 0.0);

	// Append the third point (out of order)
	if (IsOk) IsOk = GradPathAppendPoint(GradPath, 0.0, 0.5, 1.0, 2.0);

	// Insert the second point
	if (IsOk) IsOk = GradPathInsertPoint(GradPath, 1, -1.0, 0.0, 0.0, 1.0);

	// Append an incorrect fourth point
	if (IsOk) IsOk = GradPathAppendPoint(GradPath, 4.0, 4.0, 4.0, -4.0);

	// Set the correct fourth point
	if (IsOk) IsOk = GradPathSetPoint(GradPath, 4, 1.0, 1.0, 2.0, 3.0);

	// Remove the incorrect fourth point
	if (IsOk) IsOk = GradPathRemovePoint(GradPath, 3);

	// Verify that the GradPath is correct
	if (IsOk)
	{
		double X, Y, Z, Rho;

		// First point should be (X,Y,Z,Rho) = (-1, -1, -1, 0)
		IsOk = GradPathGetPoint(GradPath, 0, &X, &Y, &Z, &Rho);
		if (IsOk &&
			(X < -1.000001 || X > -0.999999 || Y < -1.000001 || Y > -0.999999 ||
			 Z < -1.000001 || Z > -0.999999 || Rho < -0.000001 || Rho > 0.000001)) IsOk = FALSE;

		// Second point should be (X,Y,Z,Rho) = (-1, 0, 0, 1)
		if (IsOk) IsOk = GradPathGetPoint(GradPath, 1, &X, &Y, &Z, &Rho);
		if (IsOk &&
			(X < -1.000001 || X > -0.999999 || Y < -0.000001 || Y > 0.000001 ||
			 Z < -0.000001 || Z >  0.000001 || Rho < 0.999999 || Rho > 1.000001)) IsOk = FALSE;

		// Third point should be (X,Y,Z,Rho) = (0, 0.5, 1, 2)
		if (IsOk) IsOk = GradPathGetPoint(GradPath, 2, &X, &Y, &Z, &Rho);
		if (IsOk &&
			(X < -0.000001 || X >  0.000001 || Y <  0.499999 || Y > 0.500001 ||
			 Z <  0.999999 || Z >  1.000001 || Rho < 1.999999 || Rho > 2.000001)) IsOk = FALSE;

		// Fourth point should be (X,Y,Z,Rho) = (1, 1, 2, 3)
		if (IsOk) IsOk = GradPathGetPoint(GradPath, 3, &X, &Y, &Z, &Rho);
		if (IsOk &&
			(X <  0.999999 || X >  1.000001 || Y <  0.999999 || Y > 1.000001 ||
			 Z <  1.999999 || Z >  2.000001 || Rho < 2.999999 || Rho > 3.000001)) IsOk = FALSE;
	}

	// Test GradPathGetLength  (should be 3 + sqrt(2))
	if (IsOk)
	{
		double Length = GradPathGetLength(GradPath);
		double CorrectLength = 3.0 + sqrt(2.0);

		if (Length < 0.0) IsOk = FALSE;

		if (IsOk &&
			(Length < CorrectLength - 0.000001 || Length > CorrectLength + 0.000001)) IsOk = FALSE;
	}

	// Test GradPathGetDistToPt
	if (IsOk)
	{
		XYZ_s Point;
		double DistToPt = -1.0;

		Point.X = -1.0;
		Point.Y =  0.5;
		Point.Z =  1.0;

		DistToPt = GradPathGetDistToPt(GradPath, Point);

		if (DistToPt < 0.0) IsOk = FALSE;

		/// TODO: test that the distance is correct
	}

	GradPathDealloc(&GradPath);

	if (GradPath != NULL) IsOk = FALSE;


	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}


/**
 * Determine if the CircleGradPath handle is sane.
 *
 * param CircleGradPath
 *     CircleGradPath structure in question.
 *
 * return
 *     TRUE if the CircleGradPath structure is valid, otherwise FALSE.
 */
Boolean_t CircleGradPathIsValid(CircleGradPath_pa CircleGradPath)
{
	Boolean_t IsValid = FALSE;

	IsValid = (VALID_REF(CircleGradPath) &&
			   VALID_BOOLEAN(CircleGradPath->CompleteCircle) &&
			   VALID_REF(CircleGradPath->Angles) && ArrListIsValid(CircleGradPath->Angles) &&
			   VALID_REF(CircleGradPath->GradPaths) && ArrListIsValid(CircleGradPath->GradPaths) &&
			   VALID_REF(CircleGradPath->ConnectPathNums) && ArrListIsValid(CircleGradPath->ConnectPathNums));

	/* Require the same count for each coordinate array. */
	if (IsValid) IsValid = (ArrListGetCount(CircleGradPath->Angles) ==
								ArrListGetCount(CircleGradPath->GradPaths));

	ENSURE(VALID_BOOLEAN(IsValid));
	return IsValid;
}





/**
 * Gets the Index'th GradPath handle.
 *
 *
 * param CircleGradPath
 *     CircleGradPath to get the GradPath from.
 *
 * param Offset
 *     Offset into the GradPath pointer ArrList.
 *
 * return
 *     GradPath handle available, otherewise a handle to NULL.
 */
GradPath_pa CircleGradPathGetGP(CircleGradPath_pa CircleGradPath,
								LgIndex_t         Offset)
{
	GradPath_pa   Result = NULL;
	ArrListItem_u Item;
	LgIndex_t     NumGradPaths;

	REQUIRE(VALID_REF(CircleGradPath));
	REQUIRE(CircleGradPathIsValid(CircleGradPath) || CircleGradPath == NULL);
	// REQUIRE(Offset >= 0 && Offset < CircleGradPathGetCount(CircleGradPath) );
	REQUIRE(Offset >= 0);

	NumGradPaths = CircleGradPathGetCount(CircleGradPath);
	if (Offset >= NumGradPaths) Offset -= NumGradPaths;

	Item = ArrListGetItem(CircleGradPath->GradPaths, Offset);
	Result = (GradPath_pa)Item.VoidPtr;

	ENSURE(VALID_REF(Result));
	return Result;
}





/**
 * Gets the starting angle at for the Offset GridPath.
 *
 *
 * param CircleGradPath
 *     CircleGradPath from which to get the Angle.
 *
 * param Offset
 *     Offset into the Angle ArrList.
 *
 * return
 *     double Angle.
 */
double CircleGradPathGetAngle(CircleGradPath_pa CircleGradPath,
							  LgIndex_t         Offset)
{
	double  Result = 0.0;
	ArrListItem_u Item;
	LgIndex_t     NumGradPaths;

	REQUIRE(VALID_REF(CircleGradPath));
	REQUIRE(CircleGradPathIsValid(CircleGradPath) || CircleGradPath == NULL);
	REQUIRE(Offset >= 0 && Offset < CircleGradPathGetCount(CircleGradPath));
	REQUIRE(Offset >= 0);

	NumGradPaths = CircleGradPathGetCount(CircleGradPath);
	if (Offset >= NumGradPaths) Offset -= NumGradPaths;

	Item = ArrListGetItem(CircleGradPath->Angles, Offset);
	Result = Item.Double;

	return Result;
}






/**
 * Gets the connector path number the Offset into ConnectPathNums
 * array list.
 *
 *
 * param CircleGradPath
 *     CircleGradPath from which to get the Angle.
 *
 * param Offset
 *     Offset into the ConnectPathNums ArrList.
 *
 * return
 *     Connector path number.
 */
LgIndex_t CircleGradPathGetCPN(CircleGradPath_pa CircleGradPath,
							   LgIndex_t         Offset)
{
	LgIndex_t  Result = 0;
	ArrListItem_u Item;

	REQUIRE(VALID_REF(CircleGradPath));
	REQUIRE(CircleGradPathIsValid(CircleGradPath));
	REQUIRE(Offset >= 0 && Offset < CircleGradPathGetCPNCount(CircleGradPath));

	Item = ArrListGetItem(CircleGradPath->ConnectPathNums, Offset);
	Result = (LgIndex_t)Item.Long;

	return Result;
}






/**
 * Deallocates the CircleGradPath handle and set the handle to NULL.
 *
 * note
 *     Item destruction is the responsibility of the caller.
 *
 * param
 *     Reference to a CircleGradPath handle.
 */
void CircleGradPathDealloc(CircleGradPath_pa *CircleGradPath)
{

	REQUIRE(VALID_REF(CircleGradPath));
	REQUIRE(CircleGradPathIsValid(*CircleGradPath) || *CircleGradPath == NULL);

	if (*CircleGradPath != NULL)
	{
		LgIndex_t ii;
		LgIndex_t Count = ArrListGetCount((*CircleGradPath)->GradPaths);

		for (ii = 0; ii < Count; ii++)
		{
			GradPath_pa GradPath = CircleGradPathGetGP((*CircleGradPath), ii);
			GradPathDealloc(&GradPath);
		}

		/* release the ArrList's */
		if ((*CircleGradPath)->Angles != NULL) ArrListDealloc(&((*CircleGradPath)->Angles));
		if ((*CircleGradPath)->GradPaths != NULL) ArrListDealloc(&((*CircleGradPath)->GradPaths));
		if ((*CircleGradPath)->ConnectPathNums != NULL) ArrListDealloc(&((*CircleGradPath)->ConnectPathNums));

		/* release the list structure itself */
		FREE_ITEM(*CircleGradPath, "CircleGradPath structure");
		*CircleGradPath = NULL;
	}

	ENSURE(*CircleGradPath == NULL);
}




/**
 * Gets the number of Angles/GradPaths currently in the CircleGradPath
 * (maintained by the CircleGradPath array lists).
 *
 * param
 *     CircleGradPath structure in question.
 *
 * return
 *     Number of Angles/GradPaths in the CircleGradPath.
 */
LgIndex_t CircleGradPathGetCount(CircleGradPath_pa CircleGradPath)
{
	LgIndex_t Result = 0;

	REQUIRE(CircleGradPathIsValid(CircleGradPath));

	Result = ArrListGetCount(CircleGradPath->Angles);

	ENSURE(Result >= 0);
	return Result;
}





/**
 * Gets the number of connector GradPaths currently in the CircleGradPath
 * (maintained by a CircleGradPath ConnectPathNums array list).
 *
 * param
 *     CircleGradPath structure in question.
 *
 * return
 *     Number of ConnectGradNums in the CircleGradPath.
 */
LgIndex_t CircleGradPathGetCPNCount(CircleGradPath_pa CircleGradPath)
{
	LgIndex_t Result = 0;

	REQUIRE(CircleGradPathIsValid(CircleGradPath));

	Result = ArrListGetCount(CircleGradPath->ConnectPathNums);

	ENSURE(Result >= 0);
	return Result;
}





/**
 * Empties the CircleGradPath of all Angles and GradPaths and resets the
 * other variables.
 *
 *
 * param CircleGradPath
 *     CircleGradPath to clear.
 */
void CircleGradPathClear(CircleGradPath_pa CircleGradPath)
{
	LgIndex_t ii;
	LgIndex_t Count = ArrListGetCount(CircleGradPath->GradPaths);

	REQUIRE(CircleGradPathIsValid(CircleGradPath));

	for (ii = 0; ii < Count; ii++)
	{
		GradPath_pa GradPath = CircleGradPathGetGP(CircleGradPath, ii);
		GradPathDealloc(&GradPath);
	}

	ArrListClear(CircleGradPath->GradPaths);
	ArrListClear(CircleGradPath->Angles);
	ArrListClear(CircleGradPath->ConnectPathNums);

	CircleGradPath->Center[0]   = 0.0;
	CircleGradPath->Center[1]   = 0.0;
	CircleGradPath->Center[2]   = 0.0;
	CircleGradPath->PrincDir[0] = 0.0;
	CircleGradPath->PrincDir[1] = 0.0;
	CircleGradPath->PrincDir[2] = 0.0;
	CircleGradPath->Radius      = -1;
	CircleGradPath->CompleteCircle = TRUE;

	ENSURE(CircleGradPathIsValid(CircleGradPath) &&
		   CircleGradPathGetCount(CircleGradPath) == 0);
}




/**
 * Allocates a CircleGradPath handle with a suitable default
 * capacity for the contained ArrList's.
 *
 *
 * return
 *     CircleGradPath handle if sufficient memory was available,
 *     otherewise a handle to NULL.
 */
CircleGradPath_pa CircleGradPathAlloc()
{
	CircleGradPath_pa Result = NULL;

	Result = ALLOC_ITEM(CircleGradPath_s, "CircleGradPath structure");
	if (Result != NULL)
	{
		Result->Center[0]   = 0.0;
		Result->Center[1]   = 0.0;
		Result->Center[2]   = 0.0;
		Result->PrincDir[0] = 0.0;
		Result->PrincDir[1] = 0.0;
		Result->PrincDir[2] = 0.0;
		Result->Radius      = 0.0;
		Result->CompleteCircle = TRUE;
		Result->Angles      = ArrListAlloc(60, ArrListType_Double);
		Result->GradPaths   = ArrListAlloc(60, ArrListType_VoidPtr);
		Result->ConnectPathNums = ArrListAlloc(60, ArrListType_Long);

		/* If it failed to allocate any of the array lists, clean-up and exit. */
		if (Result->Angles == NULL || Result->GradPaths == NULL || Result->ConnectPathNums == NULL)
		{
			if (Result->Angles != NULL) ArrListDealloc(&(Result->Angles));
			if (Result->GradPaths != NULL) ArrListDealloc(&(Result->GradPaths));
			if (Result->ConnectPathNums != NULL) ArrListDealloc(&(Result->ConnectPathNums));
			FREE_ITEM(Result, "CircleGradPath structure");
			Result = NULL;
		}
	}

	ENSURE(CircleGradPathIsValid(Result) || Result == NULL);
	return Result;
}






/**
 * Places GradPath pointer and Angle at the specified offset. If the
 * offset is beyond the end of the array lists, they are sized
 * accordingly and the intervening array values between the last
 * node of the original state and the last node of the new state are
 * assigned 0.
 *
 * note
 *     If Angles/GradPaths already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param CircleGradPath
 *     CircleGradPath target in which to set the Angle/GradPath.
 * param Offset
 *     Offset into array lists.
 * param Angle
 *     Angle to set at the specified offset.
 * param GradPath
 *     GradPath handle to set at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t   CircleGradPathSetGP(CircleGradPath_pa CircleGradPath,
								LgIndex_t   Offset,
								double      Angle,
								GradPath_pa GradPath)
{
	Boolean_t IsOk = TRUE;
	ArrListItem_u Item;


	REQUIRE(CircleGradPathIsValid(CircleGradPath));
	REQUIRE(Offset >= 0);
	REQUIRE(ArrListGetCount(CircleGradPath->Angles) == ArrListGetCount(CircleGradPath->GradPaths));

	Item.Double = Angle;
	IsOk = ArrListSetItem(CircleGradPath->Angles, Offset, Item);

	if (IsOk)
	{
		Item.VoidPtr = (void *)GradPath;
		IsOk = ArrListSetItem(CircleGradPath->GradPaths, Offset, Item);
	}

	/* Require the same count for each coordinate array (checked by GradPathIsValid). */
	if (IsOk) IsOk = CircleGradPathIsValid(CircleGradPath);

	ENSURE(CircleGradPathIsValid(CircleGradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}






/**
 * Inserts Angle/GradPath at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 * note
 *     If Angle/GradPath already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param CircleGradPath
 *     CircleGradPath target in which to set the coordinates.
 * param Offset
 *     Offset at which to insert the Angle/GradPath.
 * param Angle
 *     Angle to insert at the specified offset.
 * param GradPath
 *     GradPath handle to insert at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t CircleGradPathInsertGP(CircleGradPath_pa CircleGradPath,
								 LgIndex_t         Offset,
								 double            Angle,
								 GradPath_pa       GradPath)
{
	Boolean_t IsOk = TRUE;
	ArrListItem_u Item;


	REQUIRE(CircleGradPathIsValid(CircleGradPath));
	REQUIRE(0 <= Offset && Offset <= CircleGradPathGetCount(CircleGradPath));

	Item.Double = Angle;
	IsOk = ArrListInsertItem(CircleGradPath->Angles, Offset, Item);

	if (IsOk)
	{
		Item.VoidPtr = (void *)GradPath;
		IsOk = ArrListInsertItem(CircleGradPath->GradPaths, Offset, Item);
	}


	/* Require the same count for each coordinate array (checked by CircleGradPathIsValid). */
	if (IsOk) IsOk = CircleGradPathIsValid(CircleGradPath);

	ENSURE(CircleGradPathIsValid(CircleGradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}






/**
 * Deletes Angle/GradPath at the specified offset.
 *
 *
 * param CircleGradPath
 *     CircleGradPath target in which to set the coordinates.
 * param Offset
 *     Offset at which to insert the Angle/GradPath.
 *
 * return
 *     TRUE if successful operation, otherwise FALSE.
 */

Boolean_t CircleGradPathRemoveGP(CircleGradPath_pa CircleGradPath,
								 LgIndex_t         Offset)
{
	Boolean_t IsOk = TRUE;
	ArrListItem_u Item;


	REQUIRE(CircleGradPathIsValid(CircleGradPath));
	REQUIRE(0 <= Offset && Offset <= CircleGradPathGetCount(CircleGradPath));

	Item = ArrListRemoveItem(CircleGradPath->Angles, Offset);

	if (IsOk)
	{
		GradPath_pa GradPath;
		Item = ArrListRemoveItem(CircleGradPath->GradPaths, Offset);
		GradPath = (GradPath_pa)Item.VoidPtr;
		GradPathDealloc(&GradPath);
	}


	/* Require the same count for each coordinate array (checked by CircleGradPathIsValid). */
	if (IsOk) IsOk = CircleGradPathIsValid(CircleGradPath);

	ENSURE(CircleGradPathIsValid(CircleGradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}




/**
 * Appends the Angle/GradPath to the array lists in CircleGradPath. The
 * array lists will be expanded to accommodate the additional items.
 *
 * param CircleGradPath
 *     CircleGradPath target to which the Angle/GradPath is to be appended.
 * param Angle
 *     Angle to append to the CircleGradPath.
 * param GradPath
 *     GradPath handle to append to the CircleGradPath.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t CircleGradPathAppendGP(CircleGradPath_pa  CircleGradPath,
								 double              Angle,
								 GradPath_pa         GradPath)
{
	Boolean_t IsOk = TRUE;
	LgIndex_t Count;

	REQUIRE(CircleGradPathIsValid(CircleGradPath));

	Count = CircleGradPathGetCount(CircleGradPath);

	IsOk = CircleGradPathInsertGP(CircleGradPath, Count, Angle, GradPath);

	ENSURE(CircleGradPathIsValid(CircleGradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}







/**
 * Inserts ConnectPathNums item at the specified offset. The arrays
 * will be expanded to accomodate the additional value.
 *
 * note
 *     If Angle/GradPath already exists at the specified location they are
 *     replaced. Therefore, item destruction is the responsibility of
 *     the caller.
 *
 * param CircleGradPath
 *     CircleGradPath target in which to set the coordinates.
 * param Offset
 *     Offset at which to insert the Angle/GradPath.
 * param ConnectPathNum
 *     Number of connector path to insert at the specified offset.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */

Boolean_t CircleGradPathInsertCPN(CircleGradPath_pa CircleGradPath,
								  LgIndex_t         Offset,
								  LgIndex_t         ConnectPathNum)
{
	Boolean_t IsOk = TRUE;
	ArrListItem_u Item;


	REQUIRE(CircleGradPathIsValid(CircleGradPath));
	REQUIRE(0 <= Offset && Offset <= CircleGradPathGetCPNCount(CircleGradPath));

	Item.Long = (long)ConnectPathNum;
	IsOk = ArrListInsertItem(CircleGradPath->ConnectPathNums, Offset, Item);

	ENSURE(CircleGradPathIsValid(CircleGradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}




/**
 * Appends the ConnectPathNums item to the array list in CircleGradPath. The
 * array lists will be expanded to accommodate the additional items.
 *
 * param CircleGradPath
 *     CircleGradPath target to which the ConnectPathNum is to be appended.
 * param ConnectPathNum
 *     ConnectPathNum to append to ConnectPathNums array list in the
 *     CircleGradPath structure.
 *
 * return
 *     TRUE if sufficient memory permitted the operation, otherwise FALSE.
 */
Boolean_t CircleGradPathAppendCPN(CircleGradPath_pa CircleGradPath,
								  LgIndex_t         ConnectPathNum)
{
	Boolean_t IsOk = TRUE;
	LgIndex_t Count;

	REQUIRE(CircleGradPathIsValid(CircleGradPath));

	Count = CircleGradPathGetCPNCount(CircleGradPath);

	IsOk = CircleGradPathInsertCPN(CircleGradPath, Count, ConnectPathNum);

	ENSURE(CircleGradPathIsValid(CircleGradPath));

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}







/**
 * Return the ending critical point number for the GradPath at Offset
 * in the array list.
 *
 * param CircleGradPath
 *     CircleGradPath target to extract the EndCrtPtNum from.
 * param Offset
 *     Offset into array list of GradPath from which to extract EndCrtPtNum.
 *
 * return
 *     EndCrtPtNum for the path.
 */
LgIndex_t CircleGradPathGetEndCPNum(CircleGradPath_pa  CircleGradPath,
									LgIndex_t          Offset)
{
	LgIndex_t EndCrtPtNum = -2;

	GradPath_pa GradPath = NULL;
	ArrListItem_u Item;

	REQUIRE(CircleGradPathIsValid(CircleGradPath));
	REQUIRE(0 <= Offset && Offset < CircleGradPathGetCount(CircleGradPath));

	Item = ArrListGetItem(CircleGradPath->GradPaths, Offset);
	GradPath = (GradPath_pa)Item.VoidPtr;

	EndCrtPtNum = GradPath->EndCrtPtNum;

	ENSURE(-1 <= EndCrtPtNum);
	return EndCrtPtNum;
}


/**
 * Search the GradPath zones to find a Bond critical point (CP) that
 * connect to two specified atom CPs and a specified Ring CP.
 *
 * param CritPoints
 *     CritPoints structure containing the critical points.
 * param Atom1CrtPtNum
 *     Number of first Atom critical point.
 * param Atom2CrtPtNum
 *     Number of second Atom critical point.
 * param RingCrtPtNum
 *     Number of Ring critical point.
 *
 * return
 *     Number of Bond CP if it works, -1 otherwise.
 */

LgIndex_t  GetBondFromAtomsAndRing(CritPoints_pa CritPoints,
								   LgIndex_t     Atom1CrtPtNum,
								   LgIndex_t     Atom2CrtPtNum,
								   LgIndex_t     RingCrtPtNum)
{
	LgIndex_t  BondCrtPtNum = -1;
	LgIndex_t  BondCPBegOffset, BondCPEndOffset, ii;
	EntIndex_t TPZoneBondRing = 0;
	EntIndex_t TPZoneBondAtom1 = 0;
	EntIndex_t TPZoneBondAtom2 = 0;
	EntIndex_t SourceZoneNum = CritPoints->SourceZoneNum;

	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(Atom1CrtPtNum == -1 || CritPointsGetType(CritPoints, Atom1CrtPtNum) == -3);
	REQUIRE(Atom2CrtPtNum == -1 || CritPointsGetType(CritPoints, Atom2CrtPtNum) == -3);
	REQUIRE(CritPointsGetType(CritPoints, RingCrtPtNum) == 1);


	BondCPBegOffset = CritPointsGetBegOffset(CritPoints, (char)(-1));
	BondCPEndOffset = CritPointsGetEndOffset(CritPoints, (char)(-1));

	/* Search through the bonds critical point for one that meets the criteria */
	for (ii = BondCPBegOffset; TPZoneBondRing == 0 && ii <= BondCPEndOffset; ii++)
	{
		TPZoneBondAtom1 = GradPathTPZoneFromBegEndCP(ii, Atom1CrtPtNum, SourceZoneNum);
		if (TPZoneBondAtom1 > 0)
		{
			TPZoneBondAtom2 = GradPathTPZoneFromBegEndCP(ii, Atom2CrtPtNum, SourceZoneNum);
			if (TPZoneBondAtom2 > 0)
			{
				TPZoneBondRing = GradPathTPZoneFromBegEndCP(ii, RingCrtPtNum, SourceZoneNum);
				if (TPZoneBondRing > 0) BondCrtPtNum = ii;
			}
		}
	}

	ENSURE(BondCrtPtNum >= -1 && BondCrtPtNum < CritPointsGetCount(CritPoints));
	return BondCrtPtNum;
}


/**
 * Resample the CircleGradPath so that there are NPointsNew equally
 * spaced grid points along each GradPath length.
 *
 * param CircleGradPath
 *     CircleGradPath structure containing the GradPaths to be
 *     resampled.
 * param NumPointsNew
 *     Pointers to coordinates of the node
 *
 * return
 *     TRUE if it works, FALSE otherwise.
 */
Boolean_t CircleGradPathResample(CircleGradPath_pa CircleGradPath,
								 LgIndex_t         NumPointsNew)
{
	Boolean_t IsOk = TRUE;

	LgIndex_t   icgp;
	LgIndex_t   NumGradPaths = CircleGradPathGetCount(CircleGradPath);

	REQUIRE(CircleGradPathIsValid(CircleGradPath));
	REQUIRE(NumPointsNew > 1);


	for (icgp = 0; IsOk && icgp < NumGradPaths; icgp++)
	{
		GradPath_pa GradPath = CircleGradPathGetGP(CircleGradPath, icgp);
		double      Angle = CircleGradPathGetAngle(CircleGradPath, icgp);

		IsOk = GradPathResample(&GradPath, NumPointsNew);

		IsOk = CircleGradPathSetGP(CircleGradPath, icgp, Angle, GradPath);
	}

	ENSURE(CircleGradPathIsValid(CircleGradPath));
	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}




/*
 * GradPathCircleAdd: Compute the set of gradient paths (streamlines) resulting
 *  from a seed circle of radius Radius from point CircleCenter, normal to the vector
 *  PrincipleDir. Terminate each gradient path at a critical point or the boundary of the
 *  domain.
 */
Boolean_t ComputeGradPathAtPtAngle(const ZoneVarInfo_pa  ZoneVarInfo,
								   const double          XCrtPt,
								   const double          YCrtPt,
								   const double          ZCrtPt,
								   const double          Angle,
								   const double          n2x,
								   const double          n2y,
								   const double          n2z,
								   const double          n3x,
								   const double          n3y,
								   const double          n3z,
								   const CritPoints_pa   CritPoints,
								   const StreamDir_e     PathDir,
								   const double          CPTolerance,
								   const double          CircleRadius,
								   LgIndex_t   *NumPathPoints,
								   GradPath_pa  GradPath)
{
	Boolean_t IsOk = TRUE;
	double xp, yp, XPos, YPos, ZPos;
	double junk = 0.0;

	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(VALID_ENUM(PathDir, StreamDir_e));
	REQUIRE(CPTolerance > 0.0);
	REQUIRE(VALID_REF(NumPathPoints));
	REQUIRE(GradPathIsValid(GradPath));

	*NumPathPoints = 0;

	// xp = CPTolerance * cos(Angle);
	// yp = CPTolerance * sin(Angle);
	xp = CircleRadius * cos(Angle);
	yp = CircleRadius * sin(Angle);

	XPos = XCrtPt + xp * n2x + yp * n3x;
	YPos = YCrtPt + xp * n2y + yp * n3y;
	ZPos = ZCrtPt + xp * n2z + yp * n3z;

	/* Seed first point */
	IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos, junk); // corret ChrgDens set later

	/* Integrate gradient path line */
	IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, PathDir,
					   CritPoints, CPTolerance, NumPathPoints, GradPath);

	ENSURE(*NumPathPoints > 0);
	ENSURE(GradPathIsValid(GradPath));
	ENSURE(VALID_BOOLEAN(IsOk));

	return IsOk;
}



/*
 * Compute the beginning and ending angles for a circle with a given n1, n2, n3
 * coordinate system, and a point at a given X, Y, Z
 */
Boolean_t CircleGradPathGetAngleRange(const ZoneVarInfo_pa ZoneVarInfo,
									  const double         n1x,
									  const double         n1y,
									  const double         n1z,
									  const double         n2x,
									  const double         n2y,
									  const double         n2z,
									  const double         n3x,
									  const double         n3y,
									  const double         n3z,
									  const double         XCrtPt,
									  const double         YCrtPt,
									  const double         ZCrtPt,
									  double              *AngleBeg,
									  double              *AngleEnd)
{
	Boolean_t  IsBoundary = TRUE;
	EntIndex_t NumZones, NumVars;
	LgIndex_t  IMax, JMax, KMax;
	double NormBoundaryX  = 0.0;
	double NormBoundaryY  = 0.0;
	double NormBoundaryZ  = 0.0;
	Boolean_t  PeriodicBC = ZoneVarInfo->PeriodicBC;
	EntIndex_t ZoneNum    = ZoneVarInfo->ZoneNum;
	double XBeg = ZoneVarInfo->XBeg;
	double XEnd = ZoneVarInfo->XEnd;
	double YBeg = ZoneVarInfo->YBeg;
	double YEnd = ZoneVarInfo->YEnd;
	double ZBeg = ZoneVarInfo->ZBeg;
	double ZEnd = ZoneVarInfo->ZEnd;
	double DXCell, DYCell, DZCell;

	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(VALID_BOOLEAN(ZoneVarInfo->PeriodicBC));

	/* Get num zones in dataset */
	if (IsBoundary)
		TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	REQUIRE(ZoneNum > 0 && ZoneNum <= NumZones);


	TecUtilZoneGetInfo(ZoneNum, &IMax, &JMax, &KMax,  NULL, NULL, NULL,
					   NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	REQUIRE(IMax > 1 && JMax > 1 && KMax > 1); /* Must be 3D volume */

	DXCell = (XEnd - XBeg) / (double)(IMax - 1);
	DYCell = (YEnd - YBeg) / (double)(JMax - 1);
	DZCell = (ZEnd - ZBeg) / (double)(KMax - 1);

	// I-Minus face of zone with periodic bc
	if (XCrtPt - XBeg < 1.5 * DXCell)
		NormBoundaryX = 1.0;

	// I-Plus face of zone with periodic bc
	else if (XCrtPt - XBeg > ((double)IMax - 0.5) * DXCell)
		NormBoundaryX = -1.0;

	// J-Minus face of zone with periodic bc
	else if (YCrtPt - YBeg < 1.5 * DYCell)
		NormBoundaryY = 1.0;

	// J-Plus face of zone with periodic bc
	else if (YCrtPt - YBeg > ((double)JMax - 0.5) * DYCell)
		NormBoundaryY = -1.0;

	// K-Minus face of zone with periodic bc
	else if (ZCrtPt - ZBeg < 1.5 * DZCell)
		NormBoundaryZ = 1.0;

	// K-Plus face of zone with periodic bc
	else if (ZCrtPt - ZBeg > ((double)KMax - 0.5) * DZCell)
		NormBoundaryZ = -1.0;

	else
		IsBoundary = FALSE;

	if (IsBoundary)
	{
		double NormBoundaryDotN1, NormProjectedBoundX, NormProjectedBoundY, NormProjectedBoundZ;
		double NormProjBoundTot;
		double MedianAngle;
		// Project NormBoundary vector onto plane of circle  Normalize(NormBoundary - (NormBoundary * n1) n1)
		// Gives NormProjectedBound
		NormBoundaryDotN1   = NormBoundaryX * n1x + NormBoundaryY * n1y + NormBoundaryZ * n1z;
		NormProjectedBoundX = NormBoundaryX - NormBoundaryDotN1 * n1x;
		NormProjectedBoundY = NormBoundaryY - NormBoundaryDotN1 * n1y;
		NormProjectedBoundZ = NormBoundaryZ - NormBoundaryDotN1 * n1z;

		NormProjBoundTot = sqrt(NormProjectedBoundX * NormProjectedBoundX + NormProjectedBoundY * NormProjectedBoundY +
								NormProjectedBoundZ * NormProjectedBoundZ);
		if (NormProjBoundTot > 0.0)
		{
			NormProjectedBoundX = NormProjectedBoundX / NormProjBoundTot;
			NormProjectedBoundY = NormProjectedBoundY / NormProjBoundTot;
			NormProjectedBoundZ = NormProjectedBoundZ / NormProjBoundTot;
		}

		// ArcCosine of NormProjectedBound * n2 gives median angle. Beg and Ending angles are 90 degrees either side
		MedianAngle = RAD_TO_DEG(acos(MAX(MIN(NormProjectedBoundX * n2x + NormProjectedBoundY * n2y + NormProjectedBoundZ * n2z,1.0),-1.0)));

		// If the n3 coordinate vector is outward, add 180 to the MedianAngle
		if (NormBoundaryX * n3x + NormBoundaryY * n3y + NormBoundaryZ * n3z < 0.0) MedianAngle += 180.0;

		// Comput the beginning and ending angles
		*AngleBeg = MedianAngle - 90.0;
		*AngleEnd = MedianAngle + 90.0;
	}



	ENSURE(VALID_BOOLEAN(IsBoundary));
	return IsBoundary;
}





/**
 * Compute the offset of the minimum length gradient path.
 *
 * param CircleGradPath
 *     CircleGradPath_pa structure containing the gradient paths of interest..
 * param BegOffset, EndOffset
 *     Range of offsets of GradPaths in CircleGradPath. If CircleGradPath is a complete
 *     circle, EndOffset may wrap around and be less than BegOffset.
 *
 * return
 *     Offset of the minimum length gradient path in the range, -1 if it failed.
 */
LgIndex_t CircleGradPathMinLengthInRange(const CircleGradPath_pa CircleGradPath,
										 const LgIndex_t         BegOffset,
										 const LgIndex_t         EndOffset)
{
	LgIndex_t GPWithMinLength = 0;
	LgIndex_t ii;

	LgIndex_t NumCircleGradPaths;
	Boolean_t CompleteCircle;
	LgIndex_t EndOffsetAdjusted = EndOffset;

	double    MinLength       = LARGEDOUBLE;


	REQUIRE(CircleGradPathIsValid(CircleGradPath));

	NumCircleGradPaths = CircleGradPathGetCount(CircleGradPath);
	CompleteCircle = CircleGradPath->CompleteCircle;

	REQUIRE(BegOffset >= 0 && BegOffset < NumCircleGradPaths);
	REQUIRE(EndOffset >= 0 && EndOffset <= NumCircleGradPaths + 1 + CompleteCircle * NumCircleGradPaths);

	// Handle regions that cross branch cut
	if (CompleteCircle && EndOffset < BegOffset)
		EndOffsetAdjusted += NumCircleGradPaths;

	// Loop over range looking for minimum lenght grad path
	for (ii = BegOffset; ii < EndOffsetAdjusted; ii++)
	{
		GradPath_pa GradPath = NULL;
		double      Length = 0.0;

		/* Handle regions that cross the branch cut */
		LgIndex_t Offset = ii;
		if (CompleteCircle && Offset >= NumCircleGradPaths) Offset -= NumCircleGradPaths;

		GradPath = CircleGradPathGetGP(CircleGradPath, Offset);
		Length   = GradPathGetLength(GradPath);
		if (Length < MinLength)
		{
			MinLength = Length;
			GPWithMinLength = Offset;
		}
	}


	ENSURE(GPWithMinLength >= 0 && GPWithMinLength < NumCircleGradPaths);
	return GPWithMinLength;
}







/*
 * GradPathCircleAdd: Compute the set of gradient paths (streamlines) resulting
 *  from a seed circle of radius Radius from point CircleCenter, normal to the vector
 *  PrincipleDir. Terminate each gradient path at a critical point or the boundary of the
 *  domain.
 */
Boolean_t CircleGradPathAdd(const ZoneVarInfo_pa  ZoneVarInfo,
							const LgIndex_t       BeginCrtPtNum,
							const CritPoints_pa   CritPoints,
							const StreamDir_e     PathDir,
							const double          CPTolerance,
							const double          CircleRadius,
							CircleGradPath_pa     CircleGradPath)
{
	Boolean_t IsOk = TRUE;
	double XCrtPt, YCrtPt, ZCrtPt, PrincDirX, PrincDirY, PrincDirZ, dummy;
	char   cdummy;
	double n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z;
	double AngleBeg, AngleEnd;
	// float  CircleRadius;
	Boolean_t PeriodicBC = ZoneVarInfo->PeriodicBC;
	Boolean_t CompleteCircle = FALSE;

	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(VALID_BOOLEAN(ZoneVarInfo->PeriodicBC));
	REQUIRE(BeginCrtPtNum >= 0 && BeginCrtPtNum < CritPointsGetCount(CritPoints));
	REQUIRE(VALID_REF(CritPoints));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(VALID_ENUM(PathDir, StreamDir_e));
	REQUIRE(CPTolerance > 0.0);

	/* Extract Critical point location and principle direction */
	IsOk = CritPointsGetPoint(CritPoints, BeginCrtPtNum,
							  &XCrtPt, &YCrtPt, &ZCrtPt, &dummy, &cdummy,
							  &PrincDirX, &PrincDirY, &PrincDirZ);


	/*
	 * Find seed points on a circle surrounding the critical point.
	 * Circle is on a plane passing through the center of the critical
	 * point and perpendicular to the principle direction.
	 */
	if (IsOk)
	{
		/* Compute the unit vectors of a local orthonormal coordinate
		 * system. The Principle direction is the first of the unit
		 * vectors.
		 */
		double n2mag;

		n1x = PrincDirX;
		n1y = PrincDirY;
		n1z = PrincDirZ;

		if (ABS(n1x) < 0.9)
		{
			/* n1x is the magnitude of the projection of the x-unit-normal
			 * vector on the n1 vector.
			 */
			n2x = 1.0 - n1x * n1x;
			n2y =     - n1x * n1y;
			n2z =     - n1x * n1z;
		}
		else
		{
			n2x =     - n1y * n1x;
			n2y = 1.0 - n1y * n1y;
			n2z =     - n1y * n1z;
		}
		n2mag = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
		n2x = n2x / n2mag;
		n2y = n2y / n2mag;
		n2z = n2z / n2mag;

		n3x = n1y * n2z - n1z * n2y;
		n3y = n1z * n2x - n1x * n2z;
		n3z = n1x * n2y - n1y * n2x;
	}


	/*
	 * Compute the beginning and ending angle for periodic bc
	 */
	if (!CircleGradPathGetAngleRange(ZoneVarInfo, n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, XCrtPt, YCrtPt, ZCrtPt,
									 &AngleBeg, &AngleEnd))   
	{
		CompleteCircle = TRUE;
		AngleBeg = 0.0;
		AngleEnd = 360.0;
	}
	CircleGradPath->CompleteCircle = CompleteCircle;

	/*
	 * Cycle through the circular seed points finding the GradPaths. Place
	 * them in the CircleGradPath structure.
	 */
	if (IsOk)
	{
		int ii;
		for (ii = 0; ii < numCircularSeeds; ii++)
		{
			double Angle;
			LgIndex_t      NumPathPoints = 0;
			GradPath_pa    GradPath = NULL;

			// TEMP
			if (CompleteCircle)
				Angle = ii * DEG_TO_RAD(360.0/numCircularSeeds);
			else
				Angle = DEG_TO_RAD(AngleBeg + ii * (AngleEnd - AngleBeg) / (numCircularSeeds-1.0));

			GradPath = GradPathAlloc();

			// xp = CPTolerance * cos(Angle);
			// yp = CPTolerance * sin(Angle);

			// XPos = XCrtPt + xp * n2x + yp * n3x;
			// YPos = YCrtPt + xp * n2y + yp * n3y;
			// ZPos = ZCrtPt + xp * n2z + yp * n3z;

			GradPathClear(GradPath);

			/* Seed first point */
			// IsOk = GradPathSetPoint(GradPath, 0, XPos, YPos, ZPos);
			GradPath->BeginCrtPtNum = BeginCrtPtNum;

			/* Integrate gradient path line */
			// IsOk = GradPathAdd(ZoneVarInfo, MAXGRADPATHPOINTS, StreamDir_Reverse,
			//   CritPoints, 1.0, &NumPathPoints, GradPath);
			IsOk = ComputeGradPathAtPtAngle(ZoneVarInfo, XCrtPt, YCrtPt, ZCrtPt, Angle,
				n2x, n2y, n2z, n3x, n3y, n3z, CritPoints, PathDir, CPTolerance, CircleRadius,
											&NumPathPoints, GradPath);

			/* Analyze EndCrtPtNum's. Between every cage-cage, or cage-boundary, combination
			 * should be a Ring. If Ring is missing, refine in this region.
			 */

			/* Add GradPath and Angle to CircleGradPath */
			IsOk = CircleGradPathAppendGP(CircleGradPath, Angle, GradPath);
		}
	}

	/*
	 * Refine where divergence of adjacent gradient paths is excessive
	 */
	if (IsOk)
	{
		int level;
		double InitialDivergenceAngle = DEG_TO_RAD(360.0 / numCircularSeeds);

		// Perform multiple passes of divergence refinement
		for (level=1; level<=20; level++)
		{
			int ii;
			LgIndex_t NumGradPaths = CircleGradPathGetCount(CircleGradPath);
			GradPath_pa GradPath1 = NULL;
			GradPath_pa GradPath2 = NULL;
			double SeedAngle1, SeedAngle2;

			// Find divergence angles between adjacent gradient paths
			GradPath1 = CircleGradPathGetGP(CircleGradPath,0);
			SeedAngle1 = CircleGradPathGetAngle(CircleGradPath, 0);

			if (CompleteCircle)
			{
				SeedAngle1 += DEG_TO_RAD(360.0);
			}


			for (ii = NumGradPaths-1; ii >= 0; ii--)
			{
				double Angle;
				double x1_11, y1_11, z1_11, rho1_11;
				double x1_10, y1_10, z1_10, rho1_10;
				double x2_11, y2_11, z2_11, rho2_11;
				double x2_10, y2_10, z2_10, rho2_10;
				LgIndex_t      NumPathPoints = 0;
				GradPath_pa    GradPathNew = NULL;
				int NumRefinements = 0;
				double SeedAngleNew;

				GradPath2 = CircleGradPathGetGP(CircleGradPath,ii);
				SeedAngle2 = CircleGradPathGetAngle(CircleGradPath, ii);

				// Pick 20th segment, or half way to end, to compare angles
				if (IsOk)
				{
					LgIndex_t NumPathPoints1 = GradPathGetCount(GradPath1);
					LgIndex_t NumPathPoints2 = GradPathGetCount(GradPath2);
					LgIndex_t jj = MIN(MIN( 21, NumPathPoints1/2),NumPathPoints2/2);
					LgIndex_t jjm1 = jj - 1;
			  
					IsOk           = GradPathGetPoint(GradPath1, jj  , &x1_11, &y1_11, &z1_11, &rho1_11);
					if (IsOk) IsOk = GradPathGetPoint(GradPath1, jjm1, &x1_10, &y1_10, &z1_10, &rho1_10);
					if (IsOk) IsOk = GradPathGetPoint(GradPath2, jj  , &x2_11, &y2_11, &z2_11, &rho2_11);
					if (IsOk) IsOk = GradPathGetPoint(GradPath2, jjm1, &x2_10, &y2_10, &z2_10, &rho2_10);
				}
				
				// Compare angles (use acos of dot product divided by segment lengths)
				if (IsOk)
				{
					double dx1 = x1_11 - x1_10;
					double dy1 = y1_11 - y1_10;
					double dz1 = z1_11 - z1_10;
					double dx2 = x2_11 - x2_10;
					double dy2 = y2_11 - y2_10;
					double dz2 = z2_11 - z2_10;
					double length1 = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
					double length2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);

					double dotproduct = dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
					double cosangle = dotproduct / MAX(length1 * length2, SMALLFLOAT);

					cosangle = MAX(MIN(cosangle, 1.0), -1.0);
					Angle = acos(cosangle);
				}
// GradPath_pa GradPathAddMidway(const ZoneVarInfo_pa  ZoneVarInfo,
//                               const LgIndex_t       MaxNumPathPoints,
//                               const CritPoints_pa   CritPoints,
//                               const XYZ_s           SeedPos,
//                               const double          CPTolerance,
//                               LgIndex_t      *NumPathPoints);
				// If angle exceeds limit (some-factor * 2 * pi / number-of-segments) add a
				// grad path at an angle halfway in between
				if (IsOk && Angle > 2.0 * InitialDivergenceAngle)
				{
					SeedAngleNew = 0.5*(SeedAngle1 + SeedAngle2);

					GradPathNew = GradPathAlloc();
					GradPathClear(GradPathNew);

					GradPathNew->BeginCrtPtNum = BeginCrtPtNum;

					IsOk = ComputeGradPathAtPtAngle(ZoneVarInfo, XCrtPt, YCrtPt, ZCrtPt, SeedAngleNew,
													n2x, n2y, n2z, n3x, n3y, n3z, CritPoints, PathDir, CPTolerance, CircleRadius,
													&NumPathPoints, GradPathNew);
					// Insert GradPath and Angle into CircleGradPath
					if (IsOk)
					{
						IsOk = CircleGradPathInsertGP(CircleGradPath, ii+1, SeedAngleNew, GradPathNew);
						NumRefinements++;
					}
				}


				// Reassign pointer in preparation for next segment
				GradPath1 = GradPath2;
				SeedAngle1 = SeedAngle2;

			}
		}
	}


	/*
	 * Search through the GradPaths to find the best connector for each
	 * critical point. Also label each GradPath based on which bond-bundle
	 * surface it is a member of.
	 */
	if (IsOk)
	{
		int ii;
		int NumRefinements = 0;

		ArrList_pa DividerPathNums = ArrListAlloc(60,  ArrListType_Long);

		char      BegCrtPtType     = CritPointsGetType(CritPoints, BeginCrtPtNum);

		LgIndex_t LastEndCrtPtNum = CircleGradPathGetEndCPNum(CircleGradPath, 0);
		char      LastEndCrtPtType = 0;   /* For CrtPtNum = -1 (none), set Type = 0 */

		LgIndex_t FullCircleAddOne = 1;

		if (!CompleteCircle) FullCircleAddOne = 0;

		if (LastEndCrtPtNum >= 0) LastEndCrtPtType = CritPointsGetType(CritPoints, LastEndCrtPtNum);

		ArrListClear(DividerPathNums);

		/* Find the GradPaths that mark transition from cage-to-ring, ring-to-cage,
		 * ring-to-none (outer boundary), cage-to-none.
		 */
		for (ii = 1; IsOk && ii < CircleGradPathGetCount(CircleGradPath) + FullCircleAddOne; ii++)
		{
			/* Paths are periodic - added one Offset past end of array which is equivalent
			 * to Offset zero.
			 */
			LgIndex_t NumPaths = CircleGradPathGetCount(CircleGradPath);
			LgIndex_t Offset = ii - (ii / NumPaths) * NumPaths; // ii when ii<NumPaths, 0 when ii=NumPaths
			LgIndex_t EndCrtPtNum = CircleGradPathGetEndCPNum(CircleGradPath, Offset);

			// TEMP
			GradPath_pa GradPathTmp = CircleGradPathGetGP(CircleGradPath, Offset);
			double GradTmp = 0.0;

			if (EndCrtPtNum != LastEndCrtPtNum)
			{
				char EndCrtPtType = 0;  /* For CrtPtNum = -1 (none), set Type = 0 */
				if (EndCrtPtNum >= 0) EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);

				/* If beginning point is Bond, and switching directly from cage-to-cage,
				 * or cage-to-none, subdivide amgle until Bond-Ring line is found.
				 */
				if (BegCrtPtType == -1)
				{
					// if ( ( BegCrtPtType == -1 && LastEndCrtPtType ==  3 ) ||
					//      ( BegCrtPtType ==  1 && LastEndCrtPtType == -3 ) ||
					//       LastEndCrtPtType == 0 )
					//  {
					//    if ( ( BegCrtPtType == -1 && EndCrtPtType ==  3 ) ||
					//         ( BegCrtPtType ==  1 && EndCrtPtType == -3 ) ||
					//         EndCrtPtType == 0)
					//      {
					if (LastEndCrtPtType == 3 || LastEndCrtPtType == 0)
					{
						if (EndCrtPtType ==  3 || EndCrtPtType == 0)
						{
							/* Add Critical Path between last and current existing paths */
							double      Angle;
							LgIndex_t   NumPathPoints = 0;
							GradPath_pa GradPath = NULL;

							GradPath = GradPathAlloc();
							GradPathClear(GradPath);

							/* Seed first point */
							GradPath->BeginCrtPtNum = BeginCrtPtNum;

							/* Compute the Angle - half way between previous angles */
							{
								double AngleOffset = CircleGradPathGetAngle(CircleGradPath, Offset);
								if (Offset == 0 && ii > 0) AngleOffset += DEG_TO_RAD(360);

								Angle = 0.5 * (CircleGradPathGetAngle(CircleGradPath, ii - 1) +
											   AngleOffset);
							}

							/* Compute the Path Line */
							IsOk = ComputeGradPathAtPtAngle(ZoneVarInfo, XCrtPt, YCrtPt, ZCrtPt, Angle,
															n2x, n2y, n2z, n3x, n3y, n3z, CritPoints, PathDir, CPTolerance, CircleRadius,
															&NumPathPoints, GradPath);

							/* Insert GradPath and Angle into CircleGradPath */
							IsOk = CircleGradPathInsertGP(CircleGradPath, ii, Angle, GradPath);

							/* Decrement the ii index to redo this GradPath */
							ii--;
						}
						else
						{
							ArrListItem_u Item;

							/* Reset the "Last..." values */
							LastEndCrtPtNum = EndCrtPtNum;
							LastEndCrtPtType = EndCrtPtType;

							/* Add path to list of divider paths */
							Item.Long = Offset;
							IsOk = ArrListAppendItem(DividerPathNums, Item);
						}
					}
					// else if( (BegCrtPtType == -1 && LastEndCrtPtType ==  1) ||
					//          (BegCrtPtType ==  1 && LastEndCrtPtType == -1) )
					//   {
					//     if ( (BegCrtPtType == -1 && EndCrtPtType ==  1) ||
					//          (BegCrtPtType ==  1 && EndCrtPtType == -1) )
					//       {
					else if (LastEndCrtPtType ==  1)
					{
						if (EndCrtPtType ==  1)
						{
							/* Add Critical Path between last and current existing paths */
						}
						else
						{
							ArrListItem_u Item;

							/* Reset the "Last..." values */
							LastEndCrtPtNum  = EndCrtPtNum;
							LastEndCrtPtType = EndCrtPtType;

							/* Add path to list of divider paths */
							Item.Long = Offset;
							IsOk = ArrListAppendItem(DividerPathNums, Item);
						}
					}
				}
				/* If begining point is Ring, and switching directly from Atom to Atom,
				 * or Atom to Far-Field, find existing Bond-Ring line and insert.
				 */
				else if (BegCrtPtType == 1)
				{
					if (LastEndCrtPtType == -3 || LastEndCrtPtType == 0)
					{
						if (EndCrtPtType ==  -3 || EndCrtPtType == 0)
						{
							/* Add existing Critical Path between last and current existing paths */
							double      Angle;
							LgIndex_t   NumPathPoints = 0;
							LgIndex_t   BondCrtPtNum;

							Boolean_t  OldMethod = TRUE; // FALSE; // TEMP

							GradPath_pa GradPath = NULL;

							GradPath = GradPathAlloc();
							GradPathClear(GradPath);

							/* Seed first point */
							GradPath->BeginCrtPtNum = BeginCrtPtNum;

							/*
							 * Find the existing Bond-Ring GradPath by:
							 *  - Find the Bond critical point that connects to
							 *    atoms LastEndCrtPtNum and EndCrtPtNum, and to
							 *    the current ring (BeginCrtPtNum).
							 *  - Extract the GradPath that connects this Bond CP
							 *    to the current Ring CP.
							 */
							//if (!OldMethod)  //TMP
							if (NumRefinements >= MAXCGPREFINEMENTS)
							{
								BondCrtPtNum = GetBondFromAtomsAndRing(CritPoints,
																	   LastEndCrtPtNum,
																	   EndCrtPtNum,
																	   BeginCrtPtNum);

								if (BondCrtPtNum > -1 && BondCrtPtNum < CritPointsGetEndOffset(CritPoints, -1))
								{
									GradPath = GradPathGetByBegEndCP(BondCrtPtNum, BeginCrtPtNum);

									/* Reverse the GradPath direction and compute the Angle */
									if (GradPath != NULL)
									{
										IsOk = GradPathReverse(&GradPath);
									}
									else
										IsOk = FALSE;
								}
								else
								{
									// Temporarily disable
									// IsOk = FALSE;
								}
							}  // TEMP: !OldMethod


							/* Compute the Angle - half way between previous angles */
							if (IsOk)
							{
								double AngleOffset = CircleGradPathGetAngle(CircleGradPath, Offset);
								if (Offset == 0 && ii > 0) AngleOffset += DEG_TO_RAD(360);

								Angle = 0.5 * (CircleGradPathGetAngle(CircleGradPath, ii - 1) +
											   AngleOffset);
							}

							/* Compute the Path Line */
							// if (OldMethod)
							if (IsOk && NumRefinements < MAXCGPREFINEMENTS)
							{
								IsOk = ComputeGradPathAtPtAngle(ZoneVarInfo, XCrtPt, YCrtPt, ZCrtPt, Angle,
																n2x, n2y, n2z, n3x, n3y, n3z, CritPoints, PathDir, CPTolerance, CircleRadius,
																&NumPathPoints, GradPath);
								NumRefinements++;
							}  // TEMP: OldMethod

							// if (IsOk) Temporarily modify to handle case where GradPathGetByBegEndCP fails
							if (IsOk && (NumRefinements < MAXCGPREFINEMENTS || BondCrtPtNum > -1))
							{
								/* Insert GradPath and Angle into CircleGradPath */
								IsOk = CircleGradPathInsertGP(CircleGradPath, ii, Angle, GradPath);

								/* Decrement the ii index to redo this GradPath */
								ii--;
							}
						}
						else
						{
							ArrListItem_u Item;

							/* Reset the "Last..." values */
							LastEndCrtPtNum = EndCrtPtNum;
							LastEndCrtPtType = EndCrtPtType;

							/* Reset the refinement level to zero */
							NumRefinements = 0;

							/* Add path to list of divider paths */
							Item.Long = Offset;
							IsOk = ArrListAppendItem(DividerPathNums, Item);
						}
					}
					else if (LastEndCrtPtType ==  -1)
					{
						if (EndCrtPtType ==  -1)
						{
							/* Add Critical Path between last and current existing paths */
						}
						else
						{
							ArrListItem_u Item;

							/* Reset the "Last..." values */
							LastEndCrtPtNum  = EndCrtPtNum;
							LastEndCrtPtType = EndCrtPtType;

							/* Add path to list of divider paths */
							Item.Long = Offset;
							IsOk = ArrListAppendItem(DividerPathNums, Item);
						}
					}
				}
			}

			// TEMP -
			// GradTmp = GradPathGetMeanGrad(GradPathTmp);
		}

		/*
		 * Refine, if necessary, the Bond-Cage or Ring-Atom connectors, trying to find
		 * the minimum length.
		 */
		// TODO - STILL NOT WORKING
		if (IsOk && 0)
		{
			LgIndex_t dpn;
			LgIndex_t NumDivPaths = ArrListGetCount(DividerPathNums);
			for (dpn = 0; dpn < NumDivPaths; dpn++)
			{
				ArrListItem_u Item;
				LgIndex_t ii, itmp;
				LgIndex_t BegDivPathNum, EndDivPathNum;
				LgIndex_t EndCrtPtNum;
				LgIndex_t NumCircleGradPaths = CircleGradPathGetCount(CircleGradPath);

				LgIndex_t dpnp1 = dpn + 1;
				if (dpn == NumDivPaths - 1) dpnp1 = 0;

				Item = ArrListGetItem(DividerPathNums, dpn);
				BegDivPathNum = Item.Long;

				Item = ArrListGetItem(DividerPathNums, dpnp1);
				EndDivPathNum = Item.Long;
				if (EndDivPathNum < BegDivPathNum)
					EndDivPathNum += NumCircleGradPaths;

				EndCrtPtNum = CircleGradPathGetEndCPNum(CircleGradPath, BegDivPathNum);

				/* If End critical point number == -1, it is going to outer boundary */
				if (EndCrtPtNum >= 0)
				{
					char      EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
					LgIndex_t NumPaths     = EndDivPathNum - BegDivPathNum;

					/* Add gradpaths if there aren't at least 4 lines to chose the min from */
					if ((EndCrtPtType == 3 || EndCrtPtType == -3)
						&& NumPaths < 4)
					{
						LgIndex_t EndPathOffset   = EndDivPathNum - (EndDivPathNum / NumCircleGradPaths) * NumCircleGradPaths;
						double    BegAngle        = CircleGradPathGetAngle(CircleGradPath, BegDivPathNum);
						double    EndAngle        = CircleGradPathGetAngle(CircleGradPath, EndPathOffset);
						double    Angle;
						double    MinLength       = LARGEDOUBLE;
						LgIndex_t GPWithMinLength = 0;
						//for (ii=EndDivPathNum; ii>BegDivPathNum; ii--)
						for (itmp = NumPaths; itmp > 0; itmp--)
						{
							// LgIndex_t   Offset = ii - (ii/NumCircleGradPaths) * NumCircleGradPaths; // ii when ii<NumPaths, 0 when ii=NumPaths
							// LgIndex_t   OffsetM1 = ii - 1 - ((ii - 1)/NumCircleGradPaths) * NumCircleGradPaths; // ii when ii<NumPaths, 0 when ii=NumPaths
							LgIndex_t   Offset,  OffsetM1;
							GradPath_pa GradPath = NULL;
							LgIndex_t   NumPathPoints = 0;
							double      Length = 0.0;

							// New >>>
							/* Recompute BegDivPathNum and NumCircleGradPaths  */
							Item = ArrListGetItem(DividerPathNums, dpn);
							BegDivPathNum = Item.Long;

							NumCircleGradPaths = CircleGradPathGetCount(CircleGradPath);

							// Item = ArrListGetItem( DividerPathNums, dpnp1 );
							// EndDivPathNum = Item.Long;
							// if (EndDivPathNum < BegDivPathNum)
							//   EndDivPathNum += NumCircleGradPaths;

							ii = BegDivPathNum + itmp;

							/* Offsets=ii when less than NumCircleGradPaths, ii+NumCircleGradPaths otherwise */
							Offset = ii - (ii / NumCircleGradPaths) * NumCircleGradPaths;
							OffsetM1 = ii - 1 - ((ii - 1) / NumCircleGradPaths) * NumCircleGradPaths;
							// <<<< New

							GradPath = GradPathAlloc();
							GradPathClear(GradPath);

							/* Seed first point */
							GradPath->BeginCrtPtNum = BeginCrtPtNum;

							/* Compute the Angle - half way between previous angles */
							{
								double AngleOffset   = CircleGradPathGetAngle(CircleGradPath, Offset);
								double AngleOffsetM1 = CircleGradPathGetAngle(CircleGradPath, OffsetM1);
								if (Offset == 0 && ii > 0) AngleOffset += DEG_TO_RAD(360);

								Angle = 0.5 * (AngleOffsetM1 + AngleOffset);
							}

							/* Compute the Path Line */
							IsOk = ComputeGradPathAtPtAngle(ZoneVarInfo, XCrtPt, YCrtPt, ZCrtPt, Angle,
															n2x, n2y, n2z, n3x, n3y, n3z, CritPoints, PathDir, CPTolerance, CircleRadius,
															&NumPathPoints, GradPath);

							/* Insert GradPath and Angle into CircleGradPath */
							IsOk = CircleGradPathInsertGP(CircleGradPath, Offset, Angle, GradPath);

							/* Adjust the divider path number that follow the added path */
							if (IsOk)
							{
								LgIndex_t dpnn;
								for (dpnn = 0; dpnn < NumDivPaths; dpnn++)
								{
									LgIndex_t DivPathNumNew;
									Item = ArrListGetItem(DividerPathNums, dpnn);
									DivPathNumNew = Item.Long;
									if (DivPathNumNew >= OffsetM1)
									{
										(Item.Long)++;
										IsOk = ArrListSetItem(DividerPathNums, dpnn, Item);
									}
								}
							}
						}
					}
				}
			}
		}

		/*
		 * Eliminate redundant bond-ring lines. Keep the shortest one.
		 */
		if (IsOk)
		{
			LgIndex_t dpn;
			LgIndex_t NumDivPaths = ArrListGetCount(DividerPathNums);
			for (dpn = 0; dpn < NumDivPaths; dpn++)
			{
				ArrListItem_u Item;
				LgIndex_t ii;
				LgIndex_t BegDivPathNum, EndDivPathNum;
				LgIndex_t EndCrtPtNum;
				LgIndex_t NumCircleGradPaths = CircleGradPathGetCount(CircleGradPath);

				LgIndex_t dpnp1 = dpn + 1;
				if (dpn == NumDivPaths - 1) dpnp1 = 0;

				Item = ArrListGetItem(DividerPathNums, dpn);
				BegDivPathNum = Item.Long;

				Item = ArrListGetItem(DividerPathNums, dpnp1);
				EndDivPathNum = Item.Long;
				if (EndDivPathNum < BegDivPathNum)
					EndDivPathNum += NumCircleGradPaths;

				EndCrtPtNum = CircleGradPathGetEndCPNum(CircleGradPath, BegDivPathNum);

				/* If End critical point number == -1, it is going to outer boundary */
				if (EndCrtPtNum >= 0)
				{
					char EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
					if ((EndCrtPtType == 1 || EndCrtPtType == -1)
						&& EndDivPathNum > BegDivPathNum + 1)
					{
						LgIndex_t NumGradPaths = EndDivPathNum - BegDivPathNum;

						double    MinLength       = LARGEDOUBLE;
						LgIndex_t GPWithMinLength = 0;
						for (ii = BegDivPathNum; ii < EndDivPathNum - 1; ii++)
						{
							GradPath_pa GradPath = NULL;
							double      Length = 0.0;

							/* Handle regions that cross the branch cut */
							LgIndex_t Offset = ii;
							if (Offset >= NumCircleGradPaths) Offset -= NumCircleGradPaths;

							GradPath = CircleGradPathGetGP(CircleGradPath, Offset);
							Length   = GradPathGetLength(GradPath);
							if (Length < MinLength)
							{
								MinLength = Length;
								GPWithMinLength = Offset;
							}
						}

						/* Delete all but GPWithMinLength */
						for (ii = 0; ii < NumGradPaths; ii++)
							//for (ii=EndDivPathNum-1; ii>=BegDivPathNum; ii--)
						{
							LgIndex_t Offset;

							NumCircleGradPaths = CircleGradPathGetCount(CircleGradPath);

							Item = ArrListGetItem(DividerPathNums, dpn);
							BegDivPathNum = Item.Long;

							Item = ArrListGetItem(DividerPathNums, dpnp1);
							EndDivPathNum = Item.Long;
							if (EndDivPathNum < BegDivPathNum)
								EndDivPathNum += NumCircleGradPaths;

							Offset = BegDivPathNum + ii;
							if (Offset >= NumCircleGradPaths) Offset -= NumCircleGradPaths;
							if (Offset != GPWithMinLength)
							{
								LgIndex_t dd;

								IsOk = CircleGradPathRemoveGP(CircleGradPath, Offset);

								/* Decrement the divider path numbers */
								for (dd = 0; dd < NumDivPaths; dd++)
								{
									LgIndex_t DivPathNum;
									Item = ArrListGetItem(DividerPathNums, dd);
									DivPathNum = Item.Long;
									if (DivPathNum > Offset)
									{
										Item.Long = DivPathNum - 1;
										IsOk = ArrListSetItem(DividerPathNums, dd, Item);
									}
								}
								/* Decrement CPWithMinLength, if appropriate.  */
								if (GPWithMinLength > Offset) GPWithMinLength--;

								/* Adjust range and index of loop */
								ii--;
								NumGradPaths--;
							}
						}
					}
				}
			}
		}

		/*
		 * Extract the Bond-Ring and Bond-Cage connectors and add the offset to the
		 * ConnectPathNums list.
		 */
		if (IsOk)
		{
			LgIndex_t NumCircleGradPaths = CircleGradPathGetCount(CircleGradPath);
			LgIndex_t dpn;
			LgIndex_t NumDivPaths = ArrListGetCount(DividerPathNums);
			LgIndex_t EndCrtPtNum;
			LgIndex_t BegDivPathNum, EndDivPathNum;
			ArrListItem_u Item;

			for (dpn = 0; dpn < NumDivPaths; dpn++)
			{
				// LgIndex_t ii;
				LgIndex_t dpnp1 = dpn + 1;
				if (dpn == NumDivPaths - 1) dpnp1 = 0;

				Item = ArrListGetItem(DividerPathNums, dpn);
				BegDivPathNum = Item.Long;

				if (CompleteCircle || dpnp1 > 0)
				{
					Item = ArrListGetItem(DividerPathNums, dpnp1);
					EndDivPathNum = Item.Long;
					if (EndDivPathNum < BegDivPathNum)
						EndDivPathNum += NumCircleGradPaths;
				}
				else  // (!CompleteCircle && dpnp1 == 0)  Last segment on partial circle
					// EndDivPathNum = NumCircleGradPaths + 2;
					EndDivPathNum = NumCircleGradPaths + 1;

				if (IsOk)
				{
					LgIndex_t Offset = BegDivPathNum;
					if (Offset >= NumCircleGradPaths) Offset -= NumCircleGradPaths;
					EndCrtPtNum = CircleGradPathGetEndCPNum(CircleGradPath, Offset);
				}

				/* If End critical point number == -1, it is going to outer boundary */
				if (EndCrtPtNum >= 0)
				{
					if (EndDivPathNum == BegDivPathNum + 1)
					{
						LgIndex_t Offset = BegDivPathNum;
						if (Offset >= NumCircleGradPaths) Offset -= NumCircleGradPaths;
						IsOk = CircleGradPathAppendCPN(CircleGradPath, Offset);
					}
					else
					{
						LgIndex_t GPWithMinLength = CircleGradPathMinLengthInRange(CircleGradPath, BegDivPathNum, EndDivPathNum - 1);
						IsOk = CircleGradPathAppendCPN(CircleGradPath, GPWithMinLength);
					}
					//   }
					// else
					//   {
					//     // TODO: Handle outer boundary conditions
					//     // IsOk = CircleGradPathAppendCPN(CircleGradPath, GPWithMinLength);
				}
			}

			// When not a complete circle, fill in missing piece before BegDivPathNum
			if (NumDivPaths > 0)
			{
				Item = ArrListGetItem(DividerPathNums, 0);
				BegDivPathNum = Item.Long;
				if (!CompleteCircle && BegDivPathNum > 0)
				{
					EndDivPathNum = BegDivPathNum;
					BegDivPathNum = 0;
					if (IsOk)
					{
						LgIndex_t Offset = BegDivPathNum;
						if (Offset >= NumCircleGradPaths) Offset -= NumCircleGradPaths;
						EndCrtPtNum = CircleGradPathGetEndCPNum(CircleGradPath, Offset);
					}

					/* If End critical point number == -1, it is going to outer boundary */
					if (EndCrtPtNum >= 0)
					{
						if (EndDivPathNum == BegDivPathNum + 1)
						{
							IsOk = CircleGradPathAppendCPN(CircleGradPath, BegDivPathNum);
						}
						else
						{
							LgIndex_t GPWithMinLength = CircleGradPathMinLengthInRange(CircleGradPath, BegDivPathNum, EndDivPathNum - 1);
							IsOk = CircleGradPathInsertCPN(CircleGradPath, 0, GPWithMinLength);
						}
					}
				}
			}
		}

		/* Find the Bond-Ring connectors and add the offset to the ConnectPathNums list.
		 */


		/* Dealloc temporary array list */
		if (DividerPathNums != NULL) ArrListDealloc(&DividerPathNums);
	}




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
Boolean_t	GradPathSurfConnectBegEndToVolGPs(GradPath_pa			SurfGradPath,
											 const CritPoints_pa	SurfCritPoints,
											 const CritPoints_pa	VolCritPoints,
											 const LgIndex_t		AtomNum)
{
	Boolean_t	IsOk = TRUE;

	REQUIRE(GradPathIsValid(SurfGradPath));
	REQUIRE(CritPointsIsValid(SurfCritPoints));
	REQUIRE(CritPointsIsValid(VolCritPoints));
	REQUIRE(SurfCritPoints->Dimensions == 2);
	REQUIRE(VolCritPoints->Dimensions == 3);
	REQUIRE(AtomNum >= 0);

	LgIndex_t	BegCPNum = SurfGradPath->BeginCrtPtNum;
	LgIndex_t	EndCPNum = SurfGradPath->EndCrtPtNum;

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
	
	//	Same actions performed on beginning and end, so loop over array
	LgIndex_t	SurfCPs[]	= {BegCPNum, EndCPNum};

	XYZ_s		SurfCPXYZ;
	double		XCP, YCP, ZCP, Rho, Px, Py, Pz;
	char		SurfCPType = 1;


	//	Need atom offset
	LgIndex_t AtomCPNum = CritPointsGetBegOffset(VolCritPoints, (char)(-3)) + AtomNum;

	
	//	Begin main loop (only goes twice)
	Boolean_t Skip = FALSE;
	for (int GradTip = 0 ; IsOk && !Skip && GradTip <= 1 ; GradTip++)
	{
		//	First, get the average spacing squared of the surf grad path,
		//	only using the middle points to determine it.
		LgIndex_t	NumOfSurfPts = GradPathGetCount(SurfGradPath);
		XYZ_s		Pt1, Pt2;
		double		AvgLen = 0, TotLen = 0, CheckLen = 0, junk;
		int			AvgNum = 0;
		
		for (LgIndex_t SurfGradPt = 1 ; IsOk && SurfGradPt < NumOfSurfPts - 1 ; SurfGradPt++)
		{
			GradPathGetPoint(SurfGradPath, SurfGradPt, &Pt1.X, &Pt1.Y, &Pt1.Z, &junk);
			GradPathGetPoint(SurfGradPath, SurfGradPt+1, &Pt2.X, &Pt2.Y, &Pt2.Z, &junk);
			TotLen += sqrt(DistanceSquaredXYZ(Pt1, Pt2));
			AvgNum++;
		}
		AvgLen = TotLen / (double)AvgNum;

		////	Now check the ends of the surf grad path to see if they're much longer than
		////	the average. If they are, then remove them.
		//GradPathGetPoint(SurfGradPath, NumOfSurfPts-1, &Pt1.X, &Pt1.Y, &Pt1.Z, &junk);
		//GradPathGetPoint(SurfGradPath, NumOfSurfPts-2, &Pt2.X, &Pt2.Y, &Pt2.Z, &junk);
		//CheckLenSqr = DistanceSquaredXYZ(Pt1, Pt2);
		//if (CheckLenSqr > AvgLenSqr * 2.25) GradPathRemovePoint(SurfGradPath, NumOfSurfPts-1);

		//GradPathGetPoint(SurfGradPath, (LgIndex_t)0, &Pt1.X, &Pt1.Y, &Pt1.Z, &junk);
		//GradPathGetPoint(SurfGradPath, (LgIndex_t)1, &Pt2.X, &Pt2.Y, &Pt2.Z, &junk);
		//CheckLenSqr = DistanceSquaredXYZ(Pt1, Pt2);
		//if (CheckLenSqr > AvgLenSqr * 2.25) GradPathRemovePoint(SurfGradPath, (LgIndex_t)0);

		
		IsOk = CritPointsGetPoint(SurfCritPoints, SurfCPs[GradTip], &XCP, &YCP, &ZCP, &Rho,
			&SurfCPType, &Px, &Py, &Pz);

		if (SurfCPType == MaxCP || SurfCPType == MinCP) Skip = TRUE;
		
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
				LgIndex_t SurfCPNum = SurfCPs[GradTip] - CritPointsGetBegOffset(SurfCritPoints, SurfCPType);
				CHECK(SurfCPNum >= 0);

				//	Surface maxes and saddles could go to near or far field vol CPs
				//	so might need to check for both
				VolCPType	= 0;
				VolCPType2	= 0;

				switch (SurfCPType)
				{
				case MinCP:
					VolCPType = BondCP;
					break;
				case MaxCP:
					VolCPType = CageCP;
					VolCPType2 = CageFFCP;
					break;
				case SaddleCP:
					VolCPType = RingCP;
					VolCPType2 = RingFFCP;
					break;
				default:
					IsOk = FALSE;
				}

				/*if (SurfCPType == MinCP)
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
					IsOk = FALSE;*/

				//	Get the volume grad path, check the distance to the surf CP and append/replace
				//	as necessary
				GradPath_pa	VolGradPath = GradPathVolGPFromAtomIsoTopoSurfCPClosest(VolCPType, VolCPType2, 
																					VolCritPoints, SurfCritPoints, 
																					AtomNum, SurfCPType, 
																					SurfCPNum, AvgLen*2.0);

				//	If no grad path is found, then assume that the surface saddle corresponds
				//	to a far field ring, in which case the atom-ringFF line was seeded from
				//	the surface CP anyways, so we don't need.
				IsOk = (GradPathIsValid(VolGradPath) && VolGradPath != NULL);

				if (IsOk)
				{
					double		DistPathToSurfCP = -1.0;
					//GradPath_pa	VolGradPath = GradPathGetFromTPZone(VolGPZoneNum);

					DistPathToSurfCP = GradPathGetDistToPt(VolGradPath, SurfCPXYZ);
					if (DistPathToSurfCP < 0) IsOk = FALSE;

					//	Now need to find the intersection point of the volume grad
					//	path and the isosurface
					XYZ_s		LnBeg, LnEnd, Pt;
					double		Rho1, Rho2;
					LgIndex_t	NumOfVolGPPoints = GradPathGetCount(VolGradPath);
					if (IsOk && NumOfVolGPPoints > 0)
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

// 									//	Check that no other surf CP of the same type is actually
// 									//	closer to the returned volume GP
// 									if (IsOk)
// 									{
// 										char CheckCPType = -4, cjunk;
// 										double	djunk;
// 										XYZ_s	CheckPt;
// 										LgIndex_t	SurfCPBegOffset = CritPointsGetBegOffset(SurfCritPoints, SurfCPType);
// 										LgIndex_t	SurfCPEndOffset = CritPointsGetEndOffset(SurfCritPoints, SurfCPType);
// 										Boolean_t	CloserPtFound = FALSE;
// 										for (LgIndex_t CurrentSurfCP = SurfCPBegOffset ; !CloserPtFound 
// 											&& CurrentSurfCP < SurfCPEndOffset ; CurrentSurfCP++)
// 										{
// 											if (CurrentSurfCP != SurfCPNum)
// 											{
// 												CritPointsGetPoint(SurfCritPoints, CurrentSurfCP, 
// 																		&CheckPt.X, &CheckPt.Y, &CheckPt.Z, &djunk, 
// 																		&cjunk, &djunk, &djunk, &djunk);
// 												if (GradPathGetDistToPt(VolGradPath, CheckPt) <= DistPathToSurfCP)
// 													CloserPtFound = TRUE;
// 											}
// 										}
// 										if (!CloserPtFound) IsFound = TRUE;
// 									}
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

					//	Now loop over each point in the surf grad path, starting 
					//	the original beginning point defined and check
					//	that it's between the intersection point and end point.
					//	Check by comparing each point's distance to the first-last
					//	distance. If two points pass in succession then assume the
					//	remaining will also pass (since closer to the end point than
					//	those checked before).

					if (IsOk && SurfCPType == SaddleCP)
					{
						double		TotalLen = 0.0, CheckLen = 0.0, junk;
						XYZ_s		CheckPt, EndPt;
						int			SuccessNum = 0;
						LgIndex_t	NumOfSurfPts = GradPathGetCount(SurfGradPath);

						GradPathGetPoint(SurfGradPath, NumOfSurfPts - 1, &EndPt.X, &EndPt.Y, &EndPt.Z, &junk);
						TotalLen = DistanceSquaredXYZ(Pt,EndPt);

						IsOk = (TotalLen > 0);

						for (LgIndex_t CheckPtOffset = 0 ; IsOk && SuccessNum < 1 
							&& CheckPtOffset < GradPathGetCount(SurfGradPath) ; CheckPtOffset++)
						{
							GradPathGetPoint(SurfGradPath, CheckPtOffset, &CheckPt.X, &CheckPt.Y, &CheckPt.Z, &junk);

							CheckLen = DistanceSquaredXYZ(CheckPt, EndPt);

							if (CheckLen < TotalLen) SuccessNum++;
							else
							{
								GradPathRemovePoint(SurfGradPath, CheckPtOffset);
								CheckPtOffset--;
							}
						}

						////	For debugging
						//LgIndex_t	NumOfSurfPts = GradPathGetCount(SurfGradPath);
						//XYZ_s		Pt1, Pt2;
						//double		MinLenSqr = 0, MaxLenSqr = 0, TotLenSqr = 0, CheckLenSqr = 0;
						//int			AvgNum = 0;

						//for (LgIndex_t SurfGradPt = 0 ; IsOk && SurfGradPt < NumOfSurfPts - 1 ; SurfGradPt++)
						//{
						//	GradPathGetPoint(SurfGradPath, SurfGradPt, &Pt1.X, &Pt1.Y, &Pt1.Z, &junk);
						//	GradPathGetPoint(SurfGradPath, SurfGradPt+1, &Pt2.X, &Pt2.Y, &Pt2.Z, &junk);
						//	CheckLenSqr = DistanceSquaredXYZ(Pt1, Pt2);
						//	MaxLenSqr = MAX(CheckLenSqr, MaxLenSqr);
						//	MinLenSqr = MIN(CheckLenSqr, MinLenSqr);
						//}

						//if (MaxLenSqr / MinLenSqr >= 2)
						//{
						//	IsOk = IsOk;
						//}

						NumOfSurfPts = GradPathGetCount(SurfGradPath);
						if (NumOfSurfPts > 0)
						{
							GradPathGetPoint(SurfGradPath, (LgIndex_t)0, &SurfCPXYZ.X, &SurfCPXYZ.Y, &SurfCPXYZ.Z, &junk);
							DistPathToSurfCP = GradPathGetDistToPt(VolGradPath, SurfCPXYZ);
							if (DistPathToSurfCP > AvgLen*1.1)
							{
								//	Surf CP is far from volume grad path-isosurface
								//	intersection, so append surfCP to surf gradpath
								if (DistPathToSurfCP <= AvgLen*2)
								{
									if (GradTip < 1)
									{
										IsOk = GradPathInsertPoint(SurfGradPath, 
													(LgIndex_t)0, Pt.X, Pt.Y, Pt.Z, Rho);
									}
									else
									{
										IsOk = GradPathAppendPoint(SurfGradPath, 
													Pt.X, Pt.Y, Pt.Z, Rho);
									}
								}
							}
							else
							{
								//	Surf CP is close to volume grad path-isosurface
								//	intersection, so replace tip of surf gradpath
								//	with surfcp
								if (GradTip < 1)
								{
									IsOk = GradPathRemovePoint(SurfGradPath, (LgIndex_t)1);
									if (IsOk) IsOk = GradPathInsertPoint(SurfGradPath, 
														(LgIndex_t)0, Pt.X, Pt.Y, Pt.Z, Rho);
								}
								else
								{
									IsOk = GradPathRemovePoint(SurfGradPath, NumOfSurfPts - 1);
									if (IsOk) IsOk = GradPathAppendPoint(SurfGradPath, 
															Pt.X, Pt.Y, Pt.Z, Rho);
								}
							}
						}
						else IsOk = FALSE;
					}
				}
				GradPathDealloc(&VolGradPath);
			}
			//if (!IsOk) break;
		}
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}	//	Boolean_t	GradPathSurfConnectBegEndToVolGPs



EntIndex_t CreateConnectorZone(EntIndex_t    BaseZoneNum,
							   EntIndex_t    ChrgDensVarNum,
							   EntIndex_t    TypeVarNum,
							   CritPoints_pa CritPoints,
							   GradPath_pa   GradPath)
{
	Boolean_t   IsOk = TRUE;
	EntIndex_t  ConnectorZoneNum = TECUTILSETNOTMEMBER;
	EntIndex_t  ConnectorType;
	EntIndex_t  NumZones = 0;
	EntIndex_t  NumVars  = 0;
	LgIndex_t   NumPathPoints = GradPathGetCount(GradPath);
	char        ZoneName[200];

	LgIndex_t  BegCrtPtNum, EndCrtPtNum;

	char       BegCrtPtType;
	char       EndCrtPtType = 0;

	LgIndex_t  BegCrtPtTypeOffset;
	LgIndex_t  EndCrtPtTypeOffset = 0;

	REQUIRE(TecUtilDataSetIsAvailable());

	/* Get num zones in dataset */
	if (IsOk)
		IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

	// REQUIRE(VALID_REF(ZoneName));
	// REQUIRE(strlen(ZoneName) > 0);
	REQUIRE(NumPathPoints > 1);
	REQUIRE(VALID_REF(GradPath));
	REQUIRE(GradPathIsValid(GradPath));
	// REQUIRE(GradPathGetCount(GradPath) == NumPathPoints);
	REQUIRE(CritPointsIsValid(CritPoints));

	/* Determine Zone Name */
	if (IsOk)
	{
		char       ZoneNameBeg[200];
		BegCrtPtNum  = GradPath->BeginCrtPtNum;
		EndCrtPtNum  = GradPath->EndCrtPtNum;

		BegCrtPtType = CritPointsGetType(CritPoints, BegCrtPtNum);

		BegCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, BegCrtPtType);

		if (EndCrtPtNum >= 0)
		{
			EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
			EndCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, EndCrtPtType);
		}

		/* Determine ConnectorType
		 * 3D:
		 *   4=Ring-Cage,
		 *   2=Bond-Cage,
		 *   1=Bond-Ring,
		 *  -1=Atom-Cage (or Atom-FFCage)
		 *  -2=Ring-Atom (or FFRing-Atom)
		 *  -4=Bond-Atom
		 * 2D surface:
		 *   2 = Saddle-Min
		 *   0 = Min-Max
		 *  -2 = Saddle-Max
		 */
		ConnectorType = EndCrtPtType + BegCrtPtType;
		if (ConnectorType == 0 && (EndCrtPtType == 1 || EndCrtPtType == -1)) ConnectorType =  1;
		if (ConnectorType == 0 && (EndCrtPtType == 3 || EndCrtPtType == -3)) ConnectorType = -1;
		if (ConnectorType == 10 && (EndCrtPtType == 13 || EndCrtPtType == -3)) ConnectorType = -1;
		if (ConnectorType == 8 && (EndCrtPtType == 11 || EndCrtPtType == -3)) ConnectorType = -2;

		/* Ring--Far-Field-Cage */
		if (BegCrtPtType == 1 && EndCrtPtNum == -1) ConnectorType = 4;

		/* Create Zone Name */
		/*
		if (BegCrtPtType == -1)
		  sprintf_s(ZoneNameBeg, "Bond%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
		else if (BegCrtPtType == 1)
		  sprintf_s(ZoneNameBeg, "Ring%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
		  */
		if (BegCrtPtNum >= 0)
		{
			switch (BegCrtPtType)
			{
				case -3:
					sprintf_s(ZoneNameBeg, "Atom%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
					break;
				case -2:
					sprintf_s(ZoneNameBeg, "Max%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
					break;
				case -1:
					sprintf_s(ZoneNameBeg, "Bond%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
					break;
				case 0:
					sprintf_s(ZoneNameBeg, "Saddle%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
					break;
				case 1:
					sprintf_s(ZoneNameBeg, "Ring%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
					break;
				case 2:
					sprintf_s(ZoneNameBeg, "Min%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
					break;
				case 3:
					sprintf_s(ZoneNameBeg, "Cage%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
					break;
				case 11:
					sprintf_s(ZoneNameBeg, "RingFF%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
					break;
				case 13:
					sprintf_s(ZoneNameBeg, "CageFF%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
					break;
			}
		}

		if (GradPath->EndCrtPtNum >= 0)
		{
			switch (EndCrtPtType)
			{
				case -3:
					sprintf_s(ZoneName, "%sAtom%d_Zone%d", ZoneNameBeg, EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
					break;
				case -2:
					sprintf_s(ZoneName, "%sMax%d_Zone%d", ZoneNameBeg, EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
					break;
				case -1:
					sprintf_s(ZoneName, "%sBond%d_Zone%d", ZoneNameBeg, EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
					break;
				case 0:
					sprintf_s(ZoneName, "%sSaddle%d_Zone%d", ZoneNameBeg, EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
					break;
				case 1:
					sprintf_s(ZoneName, "%sRing%d_Zone%d", ZoneNameBeg, EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
					break;
				case 2:
					sprintf_s(ZoneName, "%sMin%d_Zone%d", ZoneNameBeg, EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
					break;
				case 3:
					sprintf_s(ZoneName, "%sCage%d_Zone%d", ZoneNameBeg, EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
					break;
				case 11:
					sprintf_s(ZoneName, "%sRingFF%d_Zone%d", ZoneNameBeg, EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
					break;
				case 13:
					sprintf_s(ZoneName, "%sCageFF%d_Zone%d", ZoneNameBeg, EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
					break;
			}
		}
		else
			sprintf_s(ZoneName, "%sCage_FF_Zone_%d", ZoneNameBeg, BaseZoneNum);
	}

	if (IsOk)
	{
		TecUtilLockStart(AddOnID);

		// AveZoneNum = TUZoneGetNumByName(ZoneName);

		/* Set FieldDataType of all variables of CP zone */
		FieldDataType_e *VarDataType = ALLOC_ARRAY(NumVars, FieldDataType_e, "VarDataType");

		for (EntIndex_t iv = 0; iv < NumVars; iv++)
			VarDataType[iv] = FieldDataType_Double;
		VarDataType[TypeVarNum-1] = FieldDataType_Int16;

		/* Create zone if it doesn't already exist. */
		if (ConnectorZoneNum == TECUTILSETNOTMEMBER)
		{
			ArgList_pa ArgList;
			// TecUtilLockStart(AddOnID);
			ArgList = TecUtilArgListAlloc();
			TecUtilArgListAppendString(ArgList, SV_NAME, ZoneName);

			TecUtilArgListAppendInt(ArgList, SV_ZONETYPE,
									(ArbParam_t)ZoneType_Ordered);

			TecUtilArgListAppendInt(ArgList, SV_IMAX, NumPathPoints);
			TecUtilArgListAppendInt(ArgList, SV_JMAX, 1);
			TecUtilArgListAppendInt(ArgList, SV_KMAX, 1);

			TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)VarDataType);

			IsOk = TecUtilDataSetAddZoneX(ArgList);

			TecUtilArgListDealloc(&ArgList);
			// TecUtilLockFinish(AddOnID);
			if (IsOk)
			{
				NumZones++;
				ConnectorZoneNum = NumZones;
			}
		}
		if ( VarDataType != NULL )
			FREE_ARRAY(VarDataType, "VarDataType");

		/* Set Critical Point Connector Pathline Locations into zone */
		if (IsOk)
		{
			FieldData_pa XVarFDPtr, YVarFDPtr, ZVarFDPtr;
			FieldData_pa CDVarFDPtr = NULL;
			FieldData_pa TypeVarFDPtr = NULL;
			LgIndex_t    ii, IJunk, JJunk, KJunk;
			TecUtilZoneGetInfo(ConnectorZoneNum, &IJunk, &JJunk, &KJunk,  &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
							   NULL, NULL, NULL, NULL, NULL, NULL, NULL);
			CDVarFDPtr = TecUtilDataValueGetWritableRef(ConnectorZoneNum, ChrgDensVarNum);
			TypeVarFDPtr = TecUtilDataValueGetWritableRef(ConnectorZoneNum, TypeVarNum);

			for (ii = 0; IsOk && ii < NumPathPoints; ii++)
			{
				double X, Y, Z, Rho;

				IsOk = GradPathGetPoint(GradPath, ii, &X, &Y, &Z, &Rho);
				TecUtilDataValueSetByRef(XVarFDPtr, ii + 1, X);
				TecUtilDataValueSetByRef(YVarFDPtr, ii + 1, Y);
				TecUtilDataValueSetByRef(ZVarFDPtr, ii + 1, Z);
				TecUtilDataValueSetByRef(CDVarFDPtr, ii + 1, Rho);

				TecUtilDataValueSetByRef(TypeVarFDPtr, ii + 1, (double)ConnectorType);
			}
		}

		/* Set Zone aux data to completely transfer necessary information to Tecplot */
		if (IsOk)
		{
			AuxData_pa ZnAuxDataRef = TecUtilAuxDataZoneGetRef(ConnectorZoneNum);
			if (ZnAuxDataRef != NULL)
			{
				char CrtPtNumStr[200];
				char CrtPtTypeStr[200];
				char BaseZoneNumStr[200];

				sprintf_s(CrtPtNumStr, "%d", BegCrtPtNum + 1);
				IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BegCrtPtNum", (ArbParam_t)CrtPtNumStr,
											 AuxDataType_String, TRUE);

				IsOk = CritPointGetTypeString(BegCrtPtNum, BegCrtPtType, CrtPtTypeStr);
				if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BegCrtPtType",
														   (ArbParam_t)CrtPtTypeStr, AuxDataType_String, TRUE);

				sprintf_s(CrtPtNumStr, "%d", EndCrtPtNum + 1);
				if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.EndCrtPtNum",
														   (ArbParam_t)CrtPtNumStr, AuxDataType_String, TRUE);

				IsOk = CritPointGetTypeString(EndCrtPtNum, EndCrtPtType, CrtPtTypeStr);
				if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.EndCrtPtType",
														   (ArbParam_t)CrtPtTypeStr, AuxDataType_String, TRUE);

				sprintf_s(BaseZoneNumStr, "%d", BaseZoneNum);
				if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BaseZoneNum",
														   (ArbParam_t)BaseZoneNumStr, AuxDataType_String, TRUE);
			}
			else IsOk = FALSE;
		}


		/* Inform Tecplot of new zone. */
		if (IsOk)
		{
			Set_pa ZSet = TecUtilSetAlloc(TRUE);
			TecUtilSetAddMember(ZSet, NumZones, TRUE);
			TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)ZSet);
			TecUtilSetDealloc(&ZSet);
		}
		TecUtilLockFinish(AddOnID);
	}

	ENSURE(ConnectorZoneNum == TECUTILSETNOTMEMBER ||
		   (ConnectorZoneNum > 0 && ConnectorZoneNum <= NumZones));
	return (ConnectorZoneNum);
}

Boolean_t	GradPathProjectPointToSurf(const ZoneVarInfo_pa ZoneVarInfo,
									   const GradPath_pa	GradPath,
									   const LgIndex_t		CurrentGPPt,
									   double*				OriginalPtArray,
									   XYZ_s*				NewPt)
{
	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(VALID_REF(OriginalPtArray));
	REQUIRE(VALID_REF(NewPt));
	REQUIRE(VALID_REF(GradPath));
	REQUIRE(GradPathIsValid(GradPath));
	REQUIRE(CurrentGPPt >= 0);

	Boolean_t	IsOk = TRUE;
	
	LgIndex_t	ElemNum = -1;
	XYZ_s		OldPt;

	OldPt.X = OriginalPtArray[0];
	OldPt.Y = OriginalPtArray[1];
	OldPt.Z	= OriginalPtArray[2];
	
	ElemNum = GeomToolsClosestElemNum(ZoneVarInfo->ZoneNum, 
										OldPt.X, OldPt.Y, OldPt.Z, 
										&NewPt->X, &NewPt->Y, &NewPt->Z);

	if (ElemNum <= 0 && CurrentGPPt >= 1)
	{
		XYZ_s	GP1, GP2;
		double  PointDist = -1.0, PathDist = -1, djunk;

		IsOk = GradPathGetPoint(GradPath, CurrentGPPt, &GP1.X, &GP1.Y, &GP1.Z, &djunk);
		if (IsOk) IsOk = GradPathGetPoint(GradPath, CurrentGPPt-1, &GP2.X, &GP2.Y, &GP2.Z, &djunk);

		if (IsOk) IsOk = (DistanceSquaredXYZ(*NewPt, OldPt) <= DistanceSquaredXYZ(GP1, GP2));
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}

Boolean_t	GradPathProjectPointToSurf(const ZoneVarInfo_pa ZoneVarInfo,
									   double				CheckDistanceSqr,
									   XYZ_s				OldPt,
									   XYZ_s*				NewPt)
{
	REQUIRE(VALID_REF(ZoneVarInfo));
	REQUIRE(VALID_REF(NewPt));
	REQUIRE(CheckDistanceSqr > 0);

	Boolean_t	IsOk = TRUE;

	LgIndex_t	ElemNum = -1;

	ElemNum = GeomToolsClosestElemNum(ZoneVarInfo->ZoneNum, 
		OldPt.X, OldPt.Y, OldPt.Z, 
		&NewPt->X, &NewPt->Y, &NewPt->Z);

	if (ElemNum <= 0)
	{
		if (IsOk) IsOk = (DistanceSquaredXYZ(*NewPt, OldPt) <= CheckDistanceSqr);
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
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
GradPath_pa GradPathVolGPFromAtomIsoTopoSurfCPClosest(const char		VolCPType1,
													const char			VolCPType2,
													const CritPoints_pa CritPoints,
													const CritPoints_pa SurfCritPoints,
													const LgIndex_t     AtomNum,
													const char          SurfCPType,
													const LgIndex_t     SurfCPOffset,
													const double        Tolerance)
{
	GradPath_pa GradPath = GradPathAlloc();

	Boolean_t IsOk = TRUE;
	LgIndex_t NumAtoms, NumNonAtomCPs;
	LgIndex_t NumSurfCPs;
	LgIndex_t AtomCPNum, SurfCPNum;
	LgIndex_t NonAtomCPBegOffset, NonAtomCPEndOffset;
	XYZ_s     SurfCPXYZ;
	char	  VolCPType;

	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(CritPointsIsValid(SurfCritPoints));
	REQUIRE(CritPoints->Dimensions == 3);
	REQUIRE(SurfCritPoints->Dimensions == 2);

	NumAtoms = CritPointsGetEndOffset(CritPoints, (char)(-3))
			   - CritPointsGetBegOffset(CritPoints, (char)(-3));

	NonAtomCPBegOffset = CritPointsGetBegOffset(CritPoints, (char)(-1));
	NonAtomCPEndOffset = CritPointsGetEndOffset(CritPoints, (char)13);
	NumNonAtomCPs = NonAtomCPEndOffset - NonAtomCPBegOffset;

	NumSurfCPs = CritPointsGetEndOffset(SurfCritPoints, SurfCPType)
				 - CritPointsGetBegOffset(SurfCritPoints, SurfCPType);

	REQUIRE(AtomNum >= 0 && AtomNum < NumAtoms);
	REQUIRE(SurfCPOffset >= 0 && SurfCPOffset < NumSurfCPs);
	REQUIRE(Tolerance > SMALLFLOAT);

	

	if (NumNonAtomCPs > 0)
	{
		LgIndex_t  VolCPNum = -1;
		double     XCP, YCP, ZCP, Rho, Px, Py, Pz;
		char       Type;
		EntIndex_t VolSourceZoneNum = CritPoints->SourceZoneNum;

		AtomCPNum = CritPointsGetBegOffset(CritPoints, (char)(-3)) + AtomNum;

		// Get the location of the IsoTopo surface  CP
		SurfCPNum = CritPointsGetBegOffset(SurfCritPoints, SurfCPType) + SurfCPOffset;

		IsOk = CritPointsGetPoint(SurfCritPoints, SurfCPNum, &XCP, &YCP, &ZCP, &Rho,
								  &Type, &Px, &Py, &Pz);

		// Verify that this CP is the right type
		if (IsOk && Type != SurfCPType) IsOk = FALSE;

		// Set the location in the SurfCPXYZ structure
		SurfCPXYZ.X = XCP;
		SurfCPXYZ.Y = YCP;
		SurfCPXYZ.Z = ZCP;

		// Cycle through volume GPs, looking for a VolCPType-AtomNum GradPath that passes
		// within Tolerance distance of MinCPNum
		double DistPathToMinCP2 = -1.0;

		for (VolCPNum = NonAtomCPBegOffset; IsOk && VolCPNum < NonAtomCPEndOffset; VolCPNum++)
		{
			VolCPType = -4;
			VolCPType = CritPointsGetType(CritPoints, VolCPNum);
			if (VolCPType == VolCPType1 || VolCPType == VolCPType2)
			{
				EntIndex_t TPZoneOfGP;

				TPZoneOfGP = GradPathTPZoneFromBegEndCP(VolCPNum, AtomCPNum, VolSourceZoneNum);

				// Bond includes specified AtomNum if TPZone is found
				if (TPZoneOfGP > 0)
				{
					double      DistPathToMinCP = -1.0;
					GradPath = GradPathGetFromTPZone(TPZoneOfGP);


					DistPathToMinCP = GradPathGetDistToPt(GradPath, SurfCPXYZ);

					// Found!!
					if (DistPathToMinCP >= 0.0 && DistPathToMinCP < Tolerance)
					{
						if (DistPathToMinCP2 == -1.0)
						{
							DistPathToMinCP2 = DistPathToMinCP;
						}
						else
						{
							DistPathToMinCP2 = MIN(DistPathToMinCP2, DistPathToMinCP);
						}
					}
				}
			}
		}
		if (DistPathToMinCP2 < 0) IsOk = FALSE;
	}

	if (IsOk == FALSE) GradPathDealloc(&GradPath);

	if (IsOk) ENSURE(GradPathIsValid(GradPath));
	return GradPath;
}

/*
 * Create the concatenated GradPath consisting of a specified GradPath in a CircleGradPath,
 * and an unspecified GradPath that starts at the end of the specified GradPath and
 * ends with a specified EndCrtPtNum.
 * CGPWeight is the number of points to shift so that they'll be on the circle grad path
 * instead of the one found by its end CP number.
 *		0  to let the length ratio hold
 *		+1 to add 1 point to circle grad path (and subtract 1 from other)
 *		-1 to add 1 point to other grad path (and subtract one from circle grad path)
 */
GradPath_pa CircleGradPathConcatenatePaths(CircleGradPath_pa CircleGradPath,
									 CritPoints_pa     CritPoints,
									 LgIndex_t         CGPNum,
									 LgIndex_t         EndCPNum,
									 LgIndex_t		   CGPWeight,
									 LgIndex_t         IMax)
{
	Boolean_t   IsOk = TRUE;
	GradPath_pa ConnectPath;


	/* Create a GradPath that contains the double line (Bond-Ring, Ring-Cage) */
	if (IsOk)
	{
		GradPath_pa Path1       = NULL;
		LgIndex_t   NumPts1     = 0;
		double      Length1     = 0.0;
		LgIndex_t   EndCPNum1   = 0;

		LgIndex_t   EndCPNum2 = 0;
		GradPath_pa Path2     = NULL;
		LgIndex_t   NumPts2   = 0;
		double      Length2   = 0.0;

		LgIndex_t   NewNumPts1, NewNumPts2;

		ConnectPath = GradPathAlloc();  /* Remember to dealloc later */

		Path1 = CircleGradPathGetGP(CircleGradPath, CGPNum);

		NumPts1     = GradPathGetCount(Path1);
		Length1     = GradPathGetLength(Path1);
		EndCPNum1   = Path1->EndCrtPtNum;

		Path2     = GradPathGetByBegEndCP(EndCPNum1, EndCPNum);
		if (Path2 == NULL)
		{
			IsOk = FALSE;
		}
		else
		{
			NumPts2   = GradPathGetCount(Path2);
			Length2   = GradPathGetLength(Path2);

			NewNumPts1 = (LgIndex_t)(IMax * (Length1 / (Length1 + Length2)) + CGPWeight);
			NewNumPts2 = IMax + 1 - NewNumPts1;

			/* Resample both GradPaths so that they have nearly the same spacing
			 * and the number of unique points adds up to IMAX. Then append the
			 * GradPaths and add the combined line to the surface zone.
			 */
			IsOk = GradPathAppend(ConnectPath, Path1);
			ConnectPath->BeginCrtPtNum = Path1->BeginCrtPtNum;
			ConnectPath->EndCrtPtNum = Path1->EndCrtPtNum;

			IsOk = GradPathResample(&ConnectPath, NewNumPts1);

			IsOk = GradPathResample(&Path2, NewNumPts2);
			IsOk = GradPathRemovePoint(Path2, 0);  /* Remove the first (duplicate) point */
			IsOk = GradPathAppend(ConnectPath, Path2);
			ConnectPath->EndCrtPtNum = Path2->EndCrtPtNum;
		}

		GradPathDealloc(&Path2);
	}

	if (IsOk)
		IsOk = GradPathIsValid(ConnectPath);

	if (IsOk == FALSE)
	{
		GradPathDealloc(&ConnectPath);
		ConnectPath = NULL;
	}

	return ConnectPath;
}


/*
	sets the directional gradient gradient variables for streamtraces
	to be used in grad path generation
	Returns true if successful, false if otherwise.
*/

Boolean_t	GradPathSTSetVariables(const ZoneVarInfo_pa  VolZoneVarInfo)
{
	Boolean_t IsOk = TRUE;
	SetValueReturnCode_e Result;

	EntIndex_t	UVarNum = VolZoneVarInfo->DGradXNum;
	EntIndex_t	VVarNum = VolZoneVarInfo->DGradYNum;
	EntIndex_t	WVarNum = VolZoneVarInfo->DGradZNum;

	REQUIRE(UVarNum > 0 && UVarNum <= VolZoneVarInfo->NumVars);
	REQUIRE(VVarNum > 0 && VVarNum <= VolZoneVarInfo->NumVars);
	REQUIRE(VVarNum > 0 && VVarNum <= VolZoneVarInfo->NumVars);

	ArgList_pa argList = TecUtilArgListAlloc();
	
	//	Set the x,y,z variables of streamtraces
	TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALTHREEDVECTOR);
	TecUtilArgListAppendString(argList, SV_P2, SV_UVAR);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, UVarNum);
	Result = TecUtilStyleSetLowLevelX(argList);

	IsOk = (Result == SetValueReturnCode_Ok || Result == SetValueReturnCode_DuplicateValue
		|| Result == SetValue_Ok || Result == SetValue_DuplicateValue);

	if (IsOk)
	{
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALTHREEDVECTOR);
		TecUtilArgListAppendString(argList, SV_P2, SV_VVAR);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, VVarNum);
		Result = TecUtilStyleSetLowLevelX(argList);

		IsOk = (Result == SetValueReturnCode_Ok || Result == SetValueReturnCode_DuplicateValue
			|| Result == SetValue_Ok || Result == SetValue_DuplicateValue);
	}

	if (IsOk)
	{
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_GLOBALTHREEDVECTOR);
		TecUtilArgListAppendString(argList, SV_P2, SV_WVAR);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, WVarNum);
		Result = TecUtilStyleSetLowLevelX(argList);

		IsOk = (Result == SetValueReturnCode_Ok || Result == SetValueReturnCode_DuplicateValue
			|| Result == SetValue_Ok || Result == SetValue_DuplicateValue);
	}

	//	Enable streamtraces
	if (IsOk)
	{
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_STREAMTRACELAYERS);
		TecUtilArgListAppendString(argList, SV_P2, SV_SHOW);
		TecUtilArgListAppendArbParam(argList, SV_IVALUE, TRUE);
		Result = TecUtilStyleSetLowLevelX(argList);

		IsOk = (Result == SetValueReturnCode_Ok || Result == SetValueReturnCode_DuplicateValue
			|| Result == SetValue_Ok || Result == SetValue_DuplicateValue);
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}


/*
	Sets the max number of steps, step size and minimum step size properties
	of streamtraces.
	Returns true if successful, false if otherwise.
*/
Boolean_t	GradPathSTSetProperties(const LgIndex_t	MaxSteps,
									const double	StepSize,
									const double	MinStepSize)
{
	Boolean_t	IsOk = TRUE;
	SetValueReturnCode_e Result;

	REQUIRE(MaxSteps > 1);
	REQUIRE(StepSize > 0);
	REQUIRE(MinStepSize > 0);

	ArgList_pa	argList = TecUtilArgListAlloc();
	//	Set the max number of steps
	TecUtilArgListClear(argList);
	TecUtilArgListAppendString(argList, SV_P1, SV_STREAMATTRIBUTES);
	TecUtilArgListAppendString(argList, SV_P2, SV_MAXSTEPS);
	TecUtilArgListAppendArbParam(argList, SV_IVALUE, MaxSteps);
	Result = TecUtilStyleSetLowLevelX(argList);

	IsOk = (Result == SetValueReturnCode_Ok || Result == SetValueReturnCode_DuplicateValue
			|| Result == SetValue_Ok || Result == SetValue_DuplicateValue);

	//	Set the step size (number of cells per step)
	if (IsOk)
	{
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_STREAMATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_CELLFRACTION);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, StepSize);
		Result = TecUtilStyleSetLowLevelX(argList);

		IsOk = (Result == SetValueReturnCode_Ok || Result == SetValueReturnCode_DuplicateValue
			|| Result == SetValue_Ok || Result == SetValue_DuplicateValue);
	}

	//	Set the minimum stepsize
	if (IsOk)
	{
		TecUtilArgListClear(argList);
		TecUtilArgListAppendString(argList, SV_P1, SV_STREAMATTRIBUTES);
		TecUtilArgListAppendString(argList, SV_P2, SV_MINCELLFRACTION);
		TecUtilArgListAppendDouble(argList, SV_DVALUE, MinStepSize);
		Result = TecUtilStyleSetLowLevelX(argList);

		IsOk = (Result == SetValueReturnCode_Ok || Result == SetValueReturnCode_DuplicateValue
			|| Result == SetValue_Ok || Result == SetValue_DuplicateValue);
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}


/*
	Uses tecplot streamtraces to generate a gradpath. 
	Generates the streamtrace,
	extracts the streamtrace as a zone,
	checks the end of the streamtrace to confirm proper terminating CP type,
	adds each point of the streamtrace to the gradpath,
	checks that gradpath is not erratic (ratio of beg-end straight line distance
		and gradpath length is within tolerace),
	returns completed gradpath.
*/
Boolean_t	GradPathSTAddSurf(const ZoneVarInfo_pa  VolZoneVarInfo,
							  const StreamDir_e     PathDir,
							  const CritPoints_pa   CritPoints,
							  LgIndex_t            *NumPathPoints,
							  const double			CPTolerance,
							  GradPath_pa           GradPath)
{
	Boolean_t		IsOk = TRUE;
	Boolean_t		IsCritPoint = FALSE;
	
	char          TestCrtPtType;    // Will only test Gradpath proximity to these CrtPts
	char          BegCrtPtType;

	REQUIRE(VALID_REF(VolZoneVarInfo));
	REQUIRE(VALID_REF(NumPathPoints));
	REQUIRE(*NumPathPoints >= 0);
	REQUIRE(VALID_ENUM(PathDir, StreamDir_e));
	REQUIRE(VALID_REF(CritPoints));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(VALID_REF(GradPath));
	REQUIRE(GradPathIsValid(GradPath));
	REQUIRE((GradPath->BeginCrtPtNum >= 0 &&
			 GradPath->BeginCrtPtNum < CritPointsGetCount(CritPoints)) ||
			GradPath->BeginCrtPtNum == -1);
	REQUIRE(VolZoneVarInfo->XVarNum > 0 && VolZoneVarInfo->XVarNum <= VolZoneVarInfo->NumVars);
	REQUIRE(VolZoneVarInfo->YVarNum > 0 && VolZoneVarInfo->YVarNum <= VolZoneVarInfo->NumVars);
	REQUIRE(VolZoneVarInfo->ZVarNum > 0 && VolZoneVarInfo->ZVarNum <= VolZoneVarInfo->NumVars);
	REQUIRE(VolZoneVarInfo->ChrgDensVarNum > 0 && VolZoneVarInfo->ChrgDensVarNum <= VolZoneVarInfo->NumVars);

	if (PathDir == StreamDir_Forward)
		TestCrtPtType = (char)-2;
	else
		TestCrtPtType = (char)2;

	//	Maximum ratio of grad path length to straight line distance of beg-end point
	double			LengthRatioTolerance = 2.0;

	GradPath->EndCrtPtNum = -1;

	XYZ_s  SeedPt[3];
	if (IsOk)
	{
		double Rho;
		double XJunk, YJunk, ZJunk;

		IsOk = GradPathGetPoint(GradPath, 0, &SeedPt->X, &SeedPt->Y, &SeedPt->Z, &Rho);
		*NumPathPoints = 1;
	}

	//	Seed the streamtrace and extract it as a zone.
	//	The proper variables (gradgrad_x,y,z) have already been set for streamtrace generation

	//	First delete any existent streamtraces
	if (IsOk)
	{
		IsOk = TecUtilStreamtraceDeleteAll();
	}
	
	//	Add the new streamtrace on the surface and extract it
	if (IsOk)
	{
		ArgList_pa	StreamArgList = TecUtilArgListAlloc();
		
		TecUtilArgListAppendInt(StreamArgList, SV_STREAMTYPE, Streamtrace_SurfaceLine);
		TecUtilArgListAppendInt(StreamArgList, SV_DISTRIBUTIONREGION, DistributionRegion_Point);
		TecUtilArgListAppendInt(StreamArgList, SV_STREAMDIRECTION, PathDir);
		TecUtilArgListAppendDouble(StreamArgList, SV_X1, SeedPt->X);
		TecUtilArgListAppendDouble(StreamArgList, SV_Y1, SeedPt->Y);
		TecUtilArgListAppendDouble(StreamArgList, SV_Z1, SeedPt->Z);

		IsOk = TecUtilStreamtraceAddX(StreamArgList);
		
		if (IsOk)
		{
			TecUtilArgListDealloc(&StreamArgList);
			IsOk = TecUtilCreateStreamZones(FALSE);

			//	Delete the created streamtrace
			TecUtilStreamtraceDeleteAll();
		}
	}

	EntIndex_t		StreamZoneNum = -1;
	FieldData_pa	XPtr, YPtr, ZPtr, RhoPtr;
	LgIndex_t		NumOfStreamPoints = -1;

	//	Streamtrace has been made, so now extract as zone and get pointers to its
	//	x,y,z, and rho field data
	if (IsOk) 
	{
		StreamZoneNum	= TecUtilDataSetGetNumZones();

		ENSURE(StreamZoneNum > 0);

		XPtr = TecUtilDataValueGetReadableNativeRef(StreamZoneNum, VolZoneVarInfo->XVarNum);
		YPtr = TecUtilDataValueGetReadableNativeRef(StreamZoneNum, VolZoneVarInfo->YVarNum);
		ZPtr = TecUtilDataValueGetReadableNativeRef(StreamZoneNum, VolZoneVarInfo->ZVarNum);
		RhoPtr = TecUtilDataValueGetReadableNativeRef(StreamZoneNum, VolZoneVarInfo->ChrgDensVarNum);

		ENSURE(VALID_REF(XPtr));
		ENSURE(VALID_REF(YPtr));
		ENSURE(VALID_REF(ZPtr));
		ENSURE(VALID_REF(RhoPtr));

		NumOfStreamPoints = TecUtilDataValueGetCountByRef(XPtr);

		ENSURE(NumOfStreamPoints > 0);
		ENSURE(NumOfStreamPoints == TecUtilDataValueGetCountByRef(YPtr));
		ENSURE(NumOfStreamPoints == TecUtilDataValueGetCountByRef(ZPtr));
		ENSURE(NumOfStreamPoints == TecUtilDataValueGetCountByRef(RhoPtr));
	}

	XYZ_s		StreamPoint;
	double		StreamRho;
	LgIndex_t	CritPointNum = -1;

	//	First check that the last point of the streamtrace terminated at the right CP
	//	If not incident with CP, then function fails and is called again with a 
	//	different seed point.
	if (IsOk)
	{
		LgIndex_t LastPt;
		if (PathDir == StreamDir_Forward)
		{
			//	Here need to find the closest CP-coincident point, so sweep over
			//	streamtrace points to find the closest pass.
			LgIndex_t i, ClosestStreamPtNum = -1;
			XYZ_s	CrtPtPt, ClosestPt, StreamPtOld;
			double	rhojunk = -1.0, djunk = -1.0, MinDistSqr = -1.0, CheckDistSqr = -1.0;
			double	TolSqr = CritPoints->MinCPDistance * CritPoints->MinCPDistance * 0.5 * 0.5;
			char	cjunk;
			Boolean_t	ClosestFound = FALSE;

			for (i = NumOfStreamPoints ; !ClosestFound && i >= 1 ; i--)
			{
				StreamPoint.X = TecUtilDataValueGetByRef(XPtr, i);
				StreamPoint.Y = TecUtilDataValueGetByRef(YPtr, i);
				StreamPoint.Z = TecUtilDataValueGetByRef(ZPtr, i);

				if (StreamPoint.X != StreamPtOld.X || StreamPoint.Y != StreamPtOld.Y
					|| StreamPoint.Z != StreamPtOld.Z)
				{
					if (!IsCritPoint) IsCritPoint = CoincidentWithSpecificCP(StreamPoint, CritPoints, 
																		CPTolerance, TestCrtPtType, &CritPointNum);
					if (IsCritPoint)
					{
						// The if(djunk) makes the CP only get fetched once
						if (rhojunk < 0) IsOk = CritPointsGetPoint(CritPoints, CritPointNum, 
														&CrtPtPt.X, &CrtPtPt.Y, &CrtPtPt.Z, &rhojunk, 
														&cjunk, &djunk, &djunk, &djunk);
						CheckDistSqr = DistanceSquaredXYZ(CrtPtPt, StreamPoint);

						if (MinDistSqr < 0 || CheckDistSqr <= MinDistSqr)
						{
							MinDistSqr = CheckDistSqr;
							ClosestStreamPtNum = i;
						}
						else
							if (CheckDistSqr > TolSqr)
							{
								if (ClosestStreamPtNum > 0)
									ClosestFound = TRUE;
							}
					}
				}
				StreamPtOld.X = StreamPoint.X;
				StreamPtOld.Y = StreamPoint.Y;
				StreamPtOld.Z = StreamPoint.Z;
			}

			if (ClosestFound) 
				NumOfStreamPoints = ClosestStreamPtNum;
			else
				IsCritPoint = FALSE;
		}
		else
		{
			StreamPoint.X = TecUtilDataValueGetByRef(XPtr, 1);
			StreamPoint.Y = TecUtilDataValueGetByRef(YPtr, 1);
			StreamPoint.Z = TecUtilDataValueGetByRef(ZPtr, 1);

			IsCritPoint = CoincidentWithSpecificCP(StreamPoint, CritPoints, 
										CPTolerance, TestCrtPtType, &CritPointNum);
		}
	}

	//	Streamtrace terminates at a correct CP, so loop over each point, appending
	//	to the grad path. This allows all the existing gradpath functions to still
	//	be employed to find things like length easily. Also it allows for a gradpath
	//	to be returned to work well with the rest of the program.
	
	if(IsOk && IsCritPoint)
	{
		if (PathDir == StreamDir_Forward)
		{
			for (LgIndex_t i = 1 ; IsOk && i <= NumOfStreamPoints ; i++)
			{
				StreamPoint.X = TecUtilDataValueGetByRef(XPtr, i);
				StreamPoint.Y = TecUtilDataValueGetByRef(YPtr, i);
				StreamPoint.Z = TecUtilDataValueGetByRef(ZPtr, i);
				StreamRho = TecUtilDataValueGetByRef(RhoPtr, i);

				IsOk = GradPathAppendPoint(GradPath, StreamPoint.X, StreamPoint.Y, StreamPoint.Z, StreamRho);
			}
		}
		else
		{
			for (LgIndex_t i = NumOfStreamPoints ; IsOk && i >= 1 ; i--)
			{
				StreamPoint.X = TecUtilDataValueGetByRef(XPtr, i);
				StreamPoint.Y = TecUtilDataValueGetByRef(YPtr, i);
				StreamPoint.Z = TecUtilDataValueGetByRef(ZPtr, i);
				StreamRho = TecUtilDataValueGetByRef(RhoPtr, i);

				IsOk = GradPathAppendPoint(GradPath, StreamPoint.X, StreamPoint.Y, StreamPoint.Z, StreamRho);
			}
		}
	}

	XYZ_s	BegCPPoint, EndCPPoint;
	double	Rho, djunk;
	char	cjunk;

	//	Grad path has been generated, so replace the original seed point
	//	with the beginning CP and append the end CP
	if (IsOk && IsCritPoint)
	{
		IsOk = CritPointsGetPoint(CritPoints, CritPointNum, 
										&EndCPPoint.X, &EndCPPoint.Y, &EndCPPoint.Z, 
										&Rho, &cjunk, &djunk, &djunk, &djunk);
		if (IsOk) IsOk = GradPathAppendPoint(GradPath, EndCPPoint.X, EndCPPoint.Y, EndCPPoint.Z, Rho);

		if (IsOk) IsOk = CritPointsGetPoint(CritPoints, GradPath->BeginCrtPtNum, 
										&BegCPPoint.X, &BegCPPoint.Y, &BegCPPoint.Z, 
										&Rho, &cjunk, &djunk, &djunk, &djunk);
		if (IsOk) IsOk = GradPathSetPoint(GradPath, (LgIndex_t)0, BegCPPoint.X, 
											BegCPPoint.Y, BegCPPoint.Z, Rho);

		*NumPathPoints = NumOfStreamPoints + 2;
	}

	if (IsOk){
		//	Delete the streamtrace zone
		Set_pa	StreamZoneSet = TecUtilSetAlloc(TRUE);
		TecUtilSetAddMember(StreamZoneSet, (SetIndex_t)StreamZoneNum, TRUE);
		TecUtilDataSetDeleteZone(StreamZoneSet);

		ENSURE(GradPathGetCount(GradPath) == *NumPathPoints);
	}

	/* Inform Tecplot of memory usage */
	if (IsOk && IsCritPoint)
	{
		GradPath->EndCrtPtNum = CritPointNum;
		
		/* Estimate of memory used by GradPath structured, 3 ArrList structures, and 3 double arrays */
		LgIndex_t MemUsageEstimate = (96 + 4 * 8 * (*NumPathPoints)) / 1024;
		if (MemUsageEstimate > 0 && GradPath->MemUsageReported < MemUsageEstimate)
		{
			Int64_t MemUsageChange = MemUsageEstimate - GradPath->MemUsageReported;
			GradPath->MemUsageReported = MemUsageEstimate;
			TecUtilMemoryChangeNotify(MemUsageChange);
		}
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return(IsOk);

} /* GradPathSTAddSurf() */


Boolean_t	GradPathCheckLenRatio(const GradPath_pa	GradPath,
								  const CritPoints_pa CritPoints,
								  const double Radius)
{
	REQUIRE(GradPathIsValid(GradPath));
	REQUIRE(CritPointsIsValid(CritPoints));
	REQUIRE(Radius > 0);

	Boolean_t	GradPathIsGood = TRUE;
	XYZ_s		BegCPPoint, EndCPPoint;
	double		GradPathLen = -1.0, GradPathEndPtDist = -1.0, djunk;
	double		LengthRatioTolerance = 1.2;
	char		cjunk;

	if (CritPoints->NumCrtPts < 9) LengthRatioTolerance /= log(sqrt((double)(CritPoints->NumCrtPts - 1)));

	GradPathIsGood = CritPointsGetPoint(CritPoints, GradPath->EndCrtPtNum, 
		&EndCPPoint.X, &EndCPPoint.Y, &EndCPPoint.Z, 
		&djunk, &cjunk, &djunk, &djunk, &djunk);

	if (GradPathIsGood) GradPathIsGood = CritPointsGetPoint(CritPoints, GradPath->BeginCrtPtNum, 
		&BegCPPoint.X, &BegCPPoint.Y, &BegCPPoint.Z, 
		&djunk, &cjunk, &djunk, &djunk, &djunk);

	GradPathLen = GradPathGetLength(GradPath);
	GradPathEndPtDist = SurfaceLength(Radius, sqrt(DistanceSquaredXYZ(EndCPPoint, BegCPPoint)));

	GradPathIsGood = (GradPathLen > 0 && GradPathEndPtDist > 0);

	if (GradPathIsGood)
		GradPathIsGood = (GradPathLen / GradPathEndPtDist < LengthRatioTolerance);

	ENSURE(VALID_BOOLEAN(GradPathIsGood));
	return GradPathIsGood;
}