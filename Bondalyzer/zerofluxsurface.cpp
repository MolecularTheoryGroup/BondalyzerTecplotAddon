/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#include <vector>
#include <string.h>
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "ADDONVER.h"
#include "TASSERT.h"
#include "GUIDEFS.h"
#include "ENGINE.h"
#include "ADKUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "GLOBAL.h"
#include "SURFELEMMAP.h"
#include "NORMALS.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
#include "GRADPATH.h"
#include "ZEROFLUXSURFACE.h"

/*
 *  Find the zone number of the zero-flux surface connecting the specified
 *  Atom-Ring-Cage.
 *
 * param BaseZoneNum
 *     Tecplot zone number of the original volume zone upon which topology is based
 * param Atom, Ring, Cage
 *     Number of Atom, Ring, and Cage for surface
 *
 * return
 *     Zero-flux surface zone number connecting the specified atom, ring, and cage,
 *       0 if not found.
 */
EntIndex_t ZFSurfZoneNumFromAtomRingCage(const EntIndex_t BaseZoneNum,
                                         const LgIndex_t  AtomNum,
                                         const LgIndex_t  RingNum,
                                         const LgIndex_t  CageNum)
{
    EntIndex_t SurfZoneNum = 0;
    EntIndex_t NumZones, ZnNum;
    Boolean_t  IsOk = TRUE;

    REQUIRE(AtomNum > 0 && RingNum > 0 && CageNum > 0);

    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, NULL);

    REQUIRE(BaseZoneNum > 0 && BaseZoneNum <= NumZones);


    // Now search for surfaces that contain the Atom-Ring-Cage numbers corresponding to
    // AtomNum, RingNum, and CageNum.
    for (ZnNum = 1; IsOk && SurfZoneNum == 0 && ZnNum <= NumZones; ZnNum++)
    {
        Boolean_t IsSurface = FALSE;
        // Is ZnNum a surface zone?
        if (TecUtilZoneIsOrdered(ZnNum))
        {
            LgIndex_t IMax, JMax, KMax;
            // Use NULL for values we're not interested in
            TecUtilZoneGetInfo(ZnNum, &IMax, &JMax, &KMax, NULL, NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, NULL);
            if (IMax > 1 && JMax > 1 && KMax == 1) IsSurface = TRUE;
        }

        // Read Aux Data and compare to Bond-Ring-Cage or Atom-Ring-Cage combinations
        if (IsOk && IsSurface)
        {
            AuxData_pa AuxDataRef = TecUtilAuxDataZoneGetRef(ZnNum);
            if (AuxDataRef != NULL)
            {
                char      *Value;
                Boolean_t  Retain;
                LgIndex_t  RingIsFirst = FALSE;

                /* Check for correct BaseZoneNum */
                if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BaseZoneNum",
                                                   &Value, &Retain))
                {
                    if (Value != NULL)
                    {
                        /* Assign the parsed value to the structure */
                        long ZoneBaseZoneNum;
                        if (sscanf_s(Value, "%d", &ZoneBaseZoneNum) == 1)
                        {
                            if (ZoneBaseZoneNum != BaseZoneNum) IsSurface = FALSE;
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

                /* Check BegCrtPtType against "ATOM" */
                if (IsOk && IsSurface)
                {
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtType",
                                                       &Value, &Retain))
                    {
                        if (strcmp(Value, "Atom") != 0) IsSurface = FALSE;
                    }
                }


                /* Check BegCrtPtNum against AtomNum */
                if (IsOk && IsSurface)
                {
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.BegCrtPtNum",
                                                       &Value, &Retain))
                    {
                        if (Value != NULL)
                        {
                            /* Assign the parsed value to the structure */
                            long BegCrtPtNum;
                            if (sscanf_s(Value, "%d", &BegCrtPtNum) == 1)
                            {
                                if (AtomNum != BegCrtPtNum) IsSurface = FALSE;
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

                /* Check EndCrtPtType against "Ring", "RingFF", "Cage", or "CageFF" */
                if (IsOk && IsSurface)
                {
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtType",
                                                       &Value, &Retain))
                    {
                        if (strcmp(Value, "Ring") == 0 || strcmp(Value, "RingFF") == 0)
                            RingIsFirst = TRUE;
                        else if (strcmp(Value, "Cage") == 0 || strcmp(Value, "CageFF") == 0)
                            RingIsFirst = FALSE;
                        else
                            IsSurface = FALSE;
                    }
                }

                /* Check ThrdCrtPtType against "Ring", "RingFF", "Cage", or "CageFF" */
                if (IsOk && IsSurface)
                {
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.ThrdCrtPtType",
                                                       &Value, &Retain))
                    {
                        if (RingIsFirst)
                        {
                            if (strcmp(Value, "Cage") != 0 && strcmp(Value, "CageFF") != 0)
                                IsSurface = FALSE;
                        }
                        else
                        {
                            if (strcmp(Value, "Ring") != 0 && strcmp(Value, "RingFF") != 0)
                                IsSurface = FALSE;
                        }
                    }
                }


                /* Get EndCrtPtNum */
                if (IsOk && IsSurface)
                {
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.EndCrtPtNum",
                                                       &Value, &Retain))
                    {
                        if (Value != NULL)
                        {
                            /* Assign the parsed value to the structure */
                            long EndCrtPtNum;
                            if (sscanf_s(Value, "%d", &EndCrtPtNum) == 1)
                            {
                                if (RingIsFirst)
                                {
                                    if (EndCrtPtNum != RingNum) IsSurface = FALSE;
                                }
                                else
                                {
                                    if (EndCrtPtNum != CageNum) IsSurface = FALSE;
                                }
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


                /* Get ThrdCrtPtNum */
                if (IsOk && IsSurface)
                {
                    if (TecUtilAuxDataGetStrItemByName(AuxDataRef, "CompChem.ThrdCrtPtNum",
                                                       &Value, &Retain))
                    {
                        if (Value != NULL)
                        {
                            /* Assign the parsed value to the structure */
                            long ThrdCrtPtNum;
                            if (sscanf_s(Value, "%d", &ThrdCrtPtNum) == 1)
                            {
                                if (RingIsFirst)
                                {
                                    if (ThrdCrtPtNum != CageNum) IsSurface = FALSE;
                                }
                                else
                                {
                                    if (ThrdCrtPtNum != RingNum) IsSurface = FALSE;
                                }
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
        }

        if (IsOk && IsSurface == TRUE) SurfZoneNum = ZnNum;
    }

    ENSURE(SurfZoneNum >= 0 && SurfZoneNum < NumZones);
    return SurfZoneNum;
}

/*
 * Create a Bond-Ring-Cage or Ring-Bond-Atom connector surface
 *
 * Parameter
 *   BaseZoneNum
 *   TypeVarNum
 *   CritPoints
 *   CircleGradPath
 *   IMax
 *   GPNumBeg
 *   GPNumEnd
 *
 * Returns TRUE if it works, FALSE otherwise.
 */
EntIndex_t CreateBRSurfaceZone(EntIndex_t        BaseZoneNum,
                               EntIndex_t        ChrgDensVarNum,
                               EntIndex_t        TypeVarNum,
                               CritPoints_pa     CritPoints,
                               CircleGradPath_pa CircleGradPath,
                               LgIndex_t         IMax,
                               LgIndex_t         CGPNumBeg,
                               LgIndex_t         CGPNumEnd)
{
    Boolean_t   IsOk = TRUE;
    EntIndex_t  SurfaceZoneNum = TECUTILSETNOTMEMBER;
    EntIndex_t  SurfaceType;
    EntIndex_t  NumZones = 0;
    EntIndex_t  NumVars  = 0;
    LgIndex_t   NumGradPaths = CGPNumEnd - CGPNumBeg + 1;
    // GradPath_pa ConnectPath = NULL; /* Path containing two concatonated paths */
    GradPath_pa ConnectPathBeg = NULL; /* Path containing two concatonated paths */
    GradPath_pa ConnectPathEnd = NULL; /* Path containing two concatonated paths */
    char        ZoneName[200];

    // New
    LgIndex_t  FirstCrtPtNum,      LastCrtPtNum;
    char       FirstCrtPtType,     LastCrtPtType;
    // End New
    LgIndex_t  BegCrtPtNum,        EndCrtPtNum,        ThrdCrtPtNum;
    char       BegCrtPtType,       EndCrtPtType,       ThrdCrtPtType;
    LgIndex_t  BegCrtPtTypeOffset, EndCrtPtTypeOffset, ThrdCrtPtTypeOffset;
    LgIndex_t  npbeg, npend;
    LgIndex_t  CGPNumBegP1, CGPNumEndM1;
    LgIndex_t  CircleGradPathCount = CircleGradPathGetCount(CircleGradPath);

    FieldData_pa XVarFDPtr = NULL;
    FieldData_pa YVarFDPtr = NULL;
    FieldData_pa ZVarFDPtr = NULL;
    FieldData_pa CDVarFDPtr = NULL;
    FieldData_pa TypeVarFDPtr = NULL;

    Boolean_t IsFarField = FALSE;

    REQUIRE(TecUtilDataSetIsAvailable());

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    REQUIRE(BaseZoneNum > 0 && BaseZoneNum <= NumZones);
    REQUIRE(TypeVarNum > 0 && TypeVarNum <= NumVars);
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(CircleGradPathIsValid(CircleGradPath));
    REQUIRE(IMax > 1);
    REQUIRE(CGPNumBeg >= 0 && CGPNumBeg < CircleGradPathCount);
    REQUIRE(CGPNumEnd >= 0 && CGPNumEnd <= 2*CircleGradPathCount);

    /* Handle case where surface crosses branch cut */
    CGPNumBegP1 = CGPNumBeg + 1;
    //tmp if (CGPNumBegP1 >= CircleGradPathCount) CGPNumBegP1 -= CircleGradPathCount;

    CGPNumEndM1 = CGPNumEnd - 1;
    //tmp if (CGPNumEndM1 < 0) CGPNumEndM1 += CircleGradPathCount;

    //tmp if (NumGradPaths < 0) NumGradPaths += CircleGradPathCount;


    /* Find the three critical points in the surface */
    if (IsOk)
    {
        GradPath_pa GradPath1 = CircleGradPathGetGP(CircleGradPath, CGPNumBeg);
        GradPath_pa GradPath2 = CircleGradPathGetGP(CircleGradPath, CGPNumEnd);

        BegCrtPtNum = GradPath1->BeginCrtPtNum;
        CHECK(BegCrtPtNum == GradPath2->BeginCrtPtNum);
        BegCrtPtType = CritPointsGetType(CritPoints, BegCrtPtNum);
        BegCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, BegCrtPtType);

        // TODO: In Progress: Refactoring determination of surface corner points
        FirstCrtPtNum  = GradPath1->EndCrtPtNum;
        if (FirstCrtPtNum >= 0)
            FirstCrtPtType = CritPointsGetType(CritPoints, FirstCrtPtNum);

        LastCrtPtNum   = GradPath2->EndCrtPtNum;
        if (LastCrtPtNum >= 0)
            LastCrtPtType  = CritPointsGetType(CritPoints, LastCrtPtNum);

        // If either first or last goes to outer boundary, it is a FarField bundle
        if (GradPath1->EndCrtPtNum == -1 || GradPath2->EndCrtPtNum == -1)
            IsFarField = TRUE;

        if (NumGradPaths > 2)
        {
            GradPath_pa GradPathMid;
            LgIndex_t   CGPNumMid = (CGPNumBeg + CGPNumEnd) / 2;

            /* Handle group that crosses branch line */
            /* tmp
            if (CGPNumEnd < CGPNumBeg)
              {
                CGPNumMid = (CGPNumBeg + CGPNumEnd + CircleGradPathCount)/2;
                if (CGPNumMid >= CircleGradPathCount) CGPNumMid -= CircleGradPathCount;
              } */

            GradPathMid = CircleGradPathGetGP(CircleGradPath, CGPNumMid);

            // If mid path goes to outer boundary, it is a FarField bundle
            if (GradPathMid->EndCrtPtNum == -1)
                IsFarField = TRUE;

            if (!IsFarField)
            {
                EndCrtPtNum = GradPathMid->EndCrtPtNum;
                EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
                EndCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, EndCrtPtType);

                /* Determine the third critical point, the one in the middle of the dog-leg line */
                if (FirstCrtPtNum == EndCrtPtNum)
                {
                    ThrdCrtPtNum  = LastCrtPtNum;
                    ThrdCrtPtType = LastCrtPtType;
                    ThrdCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, ThrdCrtPtType);
                }
                else if (LastCrtPtNum == EndCrtPtNum)
                {
                    ThrdCrtPtNum  = FirstCrtPtNum;
                    ThrdCrtPtType = FirstCrtPtType;
                    ThrdCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, ThrdCrtPtType);
                }
                else
                    IsOk = FALSE;

            }
            else   // IsFarField
            {
                EndCrtPtNum  = -1;
                EndCrtPtType = 0;
                EndCrtPtTypeOffset = -1;
                ThrdCrtPtNum  = -1;
                ThrdCrtPtType = 0;
                ThrdCrtPtTypeOffset = -1;

                // Even though it is far-field, beginning or ending lines may go
                // to valid critical point
                if (FirstCrtPtNum >= 0)
                {
                    EndCrtPtNum  = FirstCrtPtNum;
                    EndCrtPtType = FirstCrtPtType;
                    EndCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, EndCrtPtType);
                    if (LastCrtPtNum >= 0 && LastCrtPtNum != FirstCrtPtNum)
                    {
                        ThrdCrtPtNum  = LastCrtPtNum;
                        ThrdCrtPtType = LastCrtPtType;
                        ThrdCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, ThrdCrtPtType);
                    }
                }
                else
                {
                    if (LastCrtPtNum >= 0)
                    {
                        EndCrtPtNum  = LastCrtPtNum;
                        EndCrtPtType = LastCrtPtType;
                        EndCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, EndCrtPtType);
                    }
                }
            }
        }
        else
        {
            if (FirstCrtPtType == -3 || FirstCrtPtType == 3)
                EndCrtPtNum = FirstCrtPtNum;
            else
                EndCrtPtNum = LastCrtPtNum;
            EndCrtPtType = CritPointsGetType(CritPoints, EndCrtPtNum);
            EndCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, EndCrtPtType);
        }
    }

    /* Determine Zone Name */
    if (IsOk)
    {
        char       ZoneNameBeg[200];

        /* Determine SurfaceType = sum of component CritPoint types
         *   4=Ring-FarField (not sum)
         *   3=Bond-Ring-Cage
         *   1=Ring-Atom-Cage
         *  -1=Bond-Atom-Cage
         *  -3=Ring-Bond-Atom
         *  -4=Bond-FarField (not sum)
         */
        SurfaceType = EndCrtPtType + BegCrtPtType + ThrdCrtPtType;

        if (IsFarField)
        {
            if (BegCrtPtType == -1)
                SurfaceType = -4;
            else
                SurfaceType = 4;
        }

        /* Ring--Far-Field-Cage */
        // TODO handle farfield
        // if (BegCrtPtType == 1 && EndCrtPtNum == -1) ConnectorType = 4;

        /* Create Zone Name */
        if (BegCrtPtType == -1)
            sprintf_s(ZoneNameBeg, "Bond%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);
        else if (BegCrtPtType == 1)
            sprintf_s(ZoneNameBeg, "Ring%d_", BegCrtPtNum - BegCrtPtTypeOffset + 1);

        if (!IsFarField)
        {
            switch (EndCrtPtType)
            {
                case -3:
                    sprintf_s(ZoneName, "%sBond%d_Atom%d_Zone%d", ZoneNameBeg,
                            ThrdCrtPtNum - ThrdCrtPtTypeOffset + 1,
                            EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 3:
                    sprintf_s(ZoneName, "%sRing%d_Cage%d_Zone%d", ZoneNameBeg,
                            ThrdCrtPtNum - ThrdCrtPtTypeOffset + 1,
                            EndCrtPtNum - EndCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
            }
        }
        else
            sprintf_s(ZoneName, "%s_FF_Zone_%d", ZoneNameBeg, BaseZoneNum);
    }

    /* Create the Tecplot Zone */
    TecUtilLockStart(AddOnID);

    if (IsOk)
    {
        // TecUtilLockStart(AddOnID);

        // AveZoneNum = TUZoneGetNumByName(ZoneName);

        /* Set FieldDataType of all variables of CP zone */
        FieldDataType_e *VarDataType = ALLOC_ARRAY(NumVars, FieldDataType_e, "VarDataType");

        for (EntIndex_t iv = 0; iv < NumVars; iv++)
            VarDataType[iv] = FieldDataType_Float;
        VarDataType[TypeVarNum-1] = FieldDataType_Int16;

        /* Create zone if it doesn't already exist. */
        if (SurfaceZoneNum == TECUTILSETNOTMEMBER)
        {
            ArgList_pa ArgList;
            // TecUtilLockStart(AddOnID);
            ArgList = TecUtilArgListAlloc();
            TecUtilArgListAppendString(ArgList, SV_NAME, ZoneName);

            TecUtilArgListAppendInt(ArgList, SV_ZONETYPE,
                                    (ArbParam_t)ZoneType_Ordered);

            TecUtilArgListAppendInt(ArgList, SV_IMAX, IMax);
            TecUtilArgListAppendInt(ArgList, SV_JMAX, NumGradPaths);
            TecUtilArgListAppendInt(ArgList, SV_KMAX, 1);

            TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)VarDataType);

            IsOk = TecUtilDataSetAddZoneX(ArgList);

            TecUtilArgListDealloc(&ArgList);
            // TecUtilLockFinish(AddOnID);
            if (IsOk)
            {
                NumZones++;
                SurfaceZoneNum = NumZones;
            }
        }
        if ( VarDataType != NULL )
            FREE_ARRAY(VarDataType, "VarDataType");
    }

    /* Determine beginning and ending GradPaths, before considering special dog-leg */
    npbeg = CGPNumBeg;
    npend = CGPNumEnd;

    /* Identify the special coupled pathlines (Bond-Ring-Cage or Ring-Bond-Atom),
     * at the beginning and/or last line, and compute the combined line(s).
     */
    // if (IsOk && EndCrtPtNum >= 0)
    if (IsOk)
    {
        /*
        GradPath_pa ConnectPathTst = CircleGradPathGetGP(CircleGradPath, CGPNumBeg);
        LgIndex_t   EndCPNum = ConnectPathTst->EndCrtPtNum;
        char        EndCPType = CritPointsGetType(CritPoints, EndCPNum);
        */

        /*
        if (EndCPType == 1 || EndCPType == -1)
          {
            npbeg = CGPNumBegP1;
            npend = CGPNumEnd;
          }
        else
          {
            npbeg = CGPNumBeg;
            npend = CGPNumEndM1;
          }
          */
        if (!IsFarField)
        {
            if (FirstCrtPtNum != EndCrtPtNum)
            {
                npbeg = CGPNumBegP1;
                ConnectPathBeg = CircleGradPathConcatenatePaths(CircleGradPath, CritPoints,
                                                          CGPNumBeg, EndCrtPtNum, 1, IMax);
				if (ConnectPathBeg == NULL) IsOk = FALSE;
            }
            if (LastCrtPtNum != EndCrtPtNum)
            {
                npend = CGPNumEndM1;
                ConnectPathEnd = CircleGradPathConcatenatePaths(CircleGradPath, CritPoints,
                                                          CGPNumEnd, EndCrtPtNum, 1, IMax);
				if (ConnectPathEnd == NULL) IsOk = FALSE;
            }
        }
		else
		{
			/* Create a GradPath that contains the double lines (Bond-Ring, Ring-Cage) */
			if (FirstCrtPtType == 1 && LastCrtPtType == 1)
			{
				ConnectPathBeg = CircleGradPathConcatenatePaths(CircleGradPath, 
									CritPoints, CGPNumBeg, -1, 2, IMax);
				if (ConnectPathBeg == NULL) IsOk = FALSE;

				if (IsOk) ConnectPathEnd = CircleGradPathConcatenatePaths(CircleGradPath, 
											CritPoints, CGPNumEnd, -1, 2, IMax);
				if (ConnectPathEnd == NULL) IsOk = FALSE;
			}
		}

        /* Create a GradPath that contains the double lines (Bond-Ring, Ring-Cage) */
      //  
      //  if (IsOk)
      //    {
      //      GradPath_pa Path1       = NULL,
						//ConnectPath = NULL;
      //      LgIndex_t   NumPts1     = 0;
      //      double      Length1     = 0.0;
      //      LgIndex_t   EndCPNum1   = 0;

      //      LgIndex_t   EndCPNum2 = 0;
      //      GradPath_pa Path2     = NULL;
      //      LgIndex_t   NumPts2   = 0;
      //      double      Length2   = 0.0;

      //      LgIndex_t   NewNumPts1, NewNumPts2;

      //      ConnectPath = GradPathAlloc();  // Remember to dealloc later

      //      if (npbeg != CGPNumBeg)
      //        Path1 = CircleGradPathGetGP(CircleGradPath, CGPNumBeg);
      //      else
      //        Path1 = CircleGradPathGetGP(CircleGradPath, CGPNumEnd);

      //      NumPts1     = GradPathGetCount(Path1);
      //      Length1     = GradPathGetLength(Path1);
      //      EndCPNum1   = Path1->EndCrtPtNum;

      //      if (npbeg != CGPNumBeg)
      //        EndCPNum2 = (CircleGradPathGetGP(CircleGradPath, CGPNumBegP1))->EndCrtPtNum;
      //      else
      //        EndCPNum2 = (CircleGradPathGetGP(CircleGradPath, CGPNumEndM1))->EndCrtPtNum;

      //      Path2     = GradPathGetByBegEndCP(EndCPNum1, EndCPNum2);
      //      NumPts2   = GradPathGetCount(Path2);
      //      Length2   = GradPathGetLength(Path2);

      //      NewNumPts1 = (LgIndex_t)(IMax * (Length1/(Length1 + Length2)) + 1);
      //      NewNumPts2 = IMax + 1 - NewNumPts1;

      //      // Resample both GradPaths so that they have nearly the same spacing
      //      // and the number of unique points adds up to IMAX. Then append the
      //      // GradPaths and add the combined line to the surface zone.
      //      //
      //      IsOk = GradPathAppend(ConnectPath, Path1);
      //      ConnectPath->BeginCrtPtNum = Path1->BeginCrtPtNum;
      //      ConnectPath->EndCrtPtNum = Path1->EndCrtPtNum;

      //      IsOk = GradPathResample(&ConnectPath, NewNumPts1);

      //      IsOk = GradPathResample(&Path2, NewNumPts2);
      //      IsOk = GradPathRemovePoint(Path2, 0);  // Remove the first (duplicate) point
      //      IsOk = GradPathAppend(ConnectPath, Path2);
      //      ConnectPath->EndCrtPtNum = Path2->EndCrtPtNum;

      //      //TEMP - add J=1 row to IJ-ordered surface zone
      //      // IsOk = AddSurfZoneLine(SurfZoneNum, 1, ConnectPath);

      //      GradPathDealloc(&Path2);
      //    }
      //    
    }

    /* Add GradPath lines as J=const lines of the surface zone */

    /* Get FieldData pointers */
    if (IsOk)
    {
        LgIndex_t IJunk, JJunk, KJunk;
        TecUtilZoneGetInfo(SurfaceZoneNum, &IJunk, &JJunk, &KJunk,  &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
                           NULL, NULL, NULL, NULL, NULL, NULL, NULL);
        CDVarFDPtr = TecUtilDataValueGetWritableRef(SurfaceZoneNum, ChrgDensVarNum);
        TypeVarFDPtr = TecUtilDataValueGetWritableRef(SurfaceZoneNum, TypeVarNum);
    }

    /* If first line is the special line, add it first */
    if (IsOk && (npbeg != CGPNumBeg || (IsFarField && FirstCrtPtType == 1 && LastCrtPtType == 1)))
    {
        LgIndex_t ii;
        for (ii = 0; IsOk && ii < IMax; ii++)
        {
            double X, Y, Z, Rho;

            IsOk = GradPathGetPoint(ConnectPathBeg, ii, &X, &Y, &Z, &Rho);

            TecUtilDataValueSetByRef(XVarFDPtr, ii + 1, X);
            TecUtilDataValueSetByRef(YVarFDPtr, ii + 1, Y);
            TecUtilDataValueSetByRef(ZVarFDPtr, ii + 1, Z);

            TecUtilDataValueSetByRef(CDVarFDPtr, ii + 1, Rho);

            TecUtilDataValueSetByRef(TypeVarFDPtr, ii + 1, (double)SurfaceType);
        }
    }

    /* Set Critical Point Connector Pathline Locations into zone */
    if (IsOk)
    {
        LgIndex_t ii, npt, NPEndTmp;

        /* When crossing the branch line, NPEnd may be less than NPBeg.
         * Adjust so that you still get the correct np path numbers
         */
        NPEndTmp = npend;
        if (NPEndTmp < npbeg) NPEndTmp += CircleGradPathCount;

        for (npt = npbeg; IsOk && npt <= NPEndTmp; npt++)
        {
			if ((IsFarField && FirstCrtPtType == 1 && LastCrtPtType == 1) && npt == CGPNumBeg)
				npt++;

			GradPath_pa GradPath;
            LgIndex_t np = npt;
            double Angle;
            LgIndex_t jj = npt - CGPNumBeg;

			if (jj < 0) jj += CircleGradPathCount;

            /* Crossing branch cut */
            if (np >= CircleGradPathCount) np -= CircleGradPathCount;

            Angle = CircleGradPathGetAngle(CircleGradPath, np);

            GradPath = CircleGradPathGetGP(CircleGradPath, np);

            GradPathResample(&GradPath, IMax, atomISegmentRatio);

            IsOk = CircleGradPathSetGP(CircleGradPath, np, Angle, GradPath);

            for (ii = 0; IsOk && ii < IMax; ii++)
            {
                double X, Y, Z, Rho;

                LgIndex_t FDOffset = ii + 1 + jj * IMax;

                IsOk = GradPathGetPoint(GradPath, ii, &X, &Y, &Z, &Rho);
                TecUtilDataValueSetByRef(XVarFDPtr, FDOffset, X);
                TecUtilDataValueSetByRef(YVarFDPtr, FDOffset, Y);
                TecUtilDataValueSetByRef(ZVarFDPtr, FDOffset, Z);

                TecUtilDataValueSetByRef(CDVarFDPtr, FDOffset, Rho);

                TecUtilDataValueSetByRef(TypeVarFDPtr, FDOffset, (double)SurfaceType);
            }
        }
    }

    /* If last line is the special line, add it last */
    if (IsOk && (npend != CGPNumEnd) || (IsFarField && FirstCrtPtType == 1 && LastCrtPtType == 1))
    {
        LgIndex_t ii;
        for (ii = 0; IsOk && ii < IMax; ii++)
        {
            double X, Y, Z, Rho;

            LgIndex_t FDOffset = ii + 1 + (NumGradPaths - 1) * IMax;

            IsOk = GradPathGetPoint(ConnectPathEnd, ii, &X, &Y, &Z, &Rho);

            TecUtilDataValueSetByRef(XVarFDPtr, FDOffset, X);
            TecUtilDataValueSetByRef(YVarFDPtr, FDOffset, Y);
            TecUtilDataValueSetByRef(ZVarFDPtr, FDOffset, Z);

            TecUtilDataValueSetByRef(CDVarFDPtr, FDOffset, Rho);

            TecUtilDataValueSetByRef(TypeVarFDPtr, FDOffset, (double)SurfaceType);
        }
    }

    // if (ConnectPath != NULL) GradPathDealloc(&ConnectPath);
    if (ConnectPathBeg != NULL) GradPathDealloc(&ConnectPathBeg);
    if (ConnectPathEnd != NULL) GradPathDealloc(&ConnectPathEnd);


    /* Set Zone aux data to completely transfer necessary information to Tecplot */
    if (IsOk)
    {
        AuxData_pa ZnAuxDataRef = TecUtilAuxDataZoneGetRef(SurfaceZoneNum);
        if (ZnAuxDataRef != NULL)
        {
            char CrtPtNumStr[200];
            char CrtPtTypeStr[200];
            char BaseZoneNumStr[200];

            sprintf_s(CrtPtNumStr, "%d", BegCrtPtNum + 1);
            IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BegCrtPtNum", (ArbParam_t)CrtPtNumStr,
                                         AuxDataType_String, TRUE);

            if (IsOk)
            {
                char CrtPtFFStr[200];
                if (!IsFarField)
                    sprintf_s(CrtPtFFStr, "FALSE");
                else
                    sprintf_s(CrtPtFFStr, "TRUE");

                IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.IsFarField",
                                             (ArbParam_t)CrtPtFFStr, AuxDataType_String, TRUE);
            }

            IsOk = CritPointGetTypeString(BegCrtPtNum, BegCrtPtType, CrtPtTypeStr);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BegCrtPtType",
                                                       (ArbParam_t)CrtPtTypeStr, AuxDataType_String, TRUE);

            sprintf_s(CrtPtNumStr, "%d", EndCrtPtNum + 1);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.EndCrtPtNum",
                                                       (ArbParam_t)CrtPtNumStr, AuxDataType_String, TRUE);

            IsOk = CritPointGetTypeString(EndCrtPtNum, EndCrtPtType, CrtPtTypeStr);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.EndCrtPtType",
                                                       (ArbParam_t)CrtPtTypeStr, AuxDataType_String, TRUE);

            sprintf_s(CrtPtNumStr, "%d", ThrdCrtPtNum + 1);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.ThrdCrtPtNum",
                                                       (ArbParam_t)CrtPtNumStr, AuxDataType_String, TRUE);

            IsOk = CritPointGetTypeString(ThrdCrtPtNum, ThrdCrtPtType, CrtPtTypeStr);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.ThrdCrtPtType",
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

    ENSURE(SurfaceZoneNum == TECUTILSETNOTMEMBER ||
           (SurfaceZoneNum > 0 && SurfaceZoneNum <= NumZones));
    return (SurfaceZoneNum);
}

/**
 * Create a three CP (example, Atom-Ring-Cage) zero-flux connector surface
 * zone from the three CPs that make its vertices and the grad path list
 * that defines the surface.
 *
 * param BaseZoneNum
 *     Zone number upon which topology (CritPoints) is based
 * param ChrgDensVarNum, TypeVarNum
 *     Variable numbers for electron charge density and CP type
 * param CritPoints
 *     Critical points structure containing CPs connected by surface
 * param BegCrtPtNum, FirstCrtPtNum, LastCrtPtNum
 *	   The offsets of the three CPs in the critpoints variable
 * param GradPathList
 *     ArrList of pointers to a series of GradPaths which will make
 *     up the I-lines of a 2D ordered zone.
 * param IMax
 *     NumPoints in GradPath's (after resampling). Also the I-dimension
 *     of the 2D ordered zone.
 *
 * return
 *     Surface Tecplot zone number if it works, zero otherwise.
 */
EntIndex_t CreateZFSurfZoneFrom3CPsGradpathList(EntIndex_t        ChrgDensVarNum,
												  EntIndex_t        TypeVarNum,
												  CritPoints_pa     CritPoints,
												  LgIndex_t			BegCrtPtNum,
												  LgIndex_t			FirstCrtPtNum,
												  LgIndex_t			LastCrtPtNum,
												  ArrList_pa        GradPathList,
												  LgIndex_t         IMax)
{
    Boolean_t   IsOk = TRUE;
    EntIndex_t  SurfaceZoneNum = TECUTILSETNOTMEMBER;
    EntIndex_t  SurfaceType;
    EntIndex_t  NumZones = 0;
    EntIndex_t  NumVars  = 0;
    LgIndex_t   IDim, JDim;
    char        ZoneName[200];
    char        BegCrtPtType,       FirstCrtPtType,       LastCrtPtType;
    LgIndex_t   BegCrtPtTypeOffset, FirstCrtPtTypeOffset, LastCrtPtTypeOffset;
    EntIndex_t  BaseZoneNum;

    FieldData_pa XVarFDPtr = NULL;
    FieldData_pa YVarFDPtr = NULL;
    FieldData_pa ZVarFDPtr = NULL;
    FieldData_pa CDVarFDPtr = NULL;
    FieldData_pa TypeVarFDPtr = NULL;

    TecUtilLockStart(AddOnID);

    REQUIRE(TecUtilDataSetIsAvailable());

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
    REQUIRE(TypeVarNum > 0 && TypeVarNum <= NumVars);
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(ArrListIsValid(GradPathList));
    REQUIRE(IMax > 1);

    BaseZoneNum = CritPoints->SourceZoneNum;
    REQUIRE(BaseZoneNum > 0 && BaseZoneNum <= NumZones);

    // Determine dimensions of the Tecplot surface zone and create it
    if (IsOk)
    {
        IDim = IMax;
        JDim = ArrListGetCount(GradPathList);
    }

    // Determine Zone Name
    if (IsOk)
    {
        char          ZoneNameBeg[200];
		char          ZoneNameMid[200];
		ArrListItem_u Item;
		GradPath_pa   GradPathFirst = NULL;
		GradPath_pa   GradPathLast = NULL;

		Item = ArrListGetItem(GradPathList, 0);
		GradPathFirst = (GradPath_pa)Item.VoidPtr;

		Item = ArrListGetItem(GradPathList, JDim - 1);
		GradPathLast = (GradPath_pa)Item.VoidPtr;

		
		CHECK(BegCrtPtNum == GradPathLast->BeginCrtPtNum
			&& LastCrtPtNum == GradPathLast->EndCrtPtNum);

        BegCrtPtType = CritPointsGetType(CritPoints, BegCrtPtNum);
        BegCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, BegCrtPtType);

        if (FirstCrtPtNum >= 0)
        {
            FirstCrtPtType = CritPointsGetType(CritPoints, FirstCrtPtNum);
            FirstCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, FirstCrtPtType);
        }

        if (LastCrtPtNum >= 0)
        {
            LastCrtPtType  = CritPointsGetType(CritPoints, LastCrtPtNum);
            LastCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, LastCrtPtType);
        }

        /* Determine SurfaceType = sum of component CritPoint types
         *   4=Ring-FarField (not sum)
         *   3=Bond-Ring-Cage
         *   1=Ring-Atom-Cage
         *  -1=Bond-Atom-Cage
         *  -3=Ring-Bond-Atom
         *  -4=Bond-FarField (not sum)
         */
        SurfaceType = BegCrtPtType + FirstCrtPtType + LastCrtPtType;

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

        if (FirstCrtPtNum >= 0)
        {
            switch (FirstCrtPtType)
            {
                case -3:
                    sprintf_s(ZoneNameMid, "%sAtom%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case -2:
                    sprintf_s(ZoneNameMid, "%sMax%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case -1:
                    sprintf_s(ZoneNameMid, "%sBond%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 0:
                    sprintf_s(ZoneNameMid, "%sSaddle%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 1:
                    sprintf_s(ZoneNameMid, "%sRing%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 2:
                    sprintf_s(ZoneNameMid, "%sMin%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 3:
                    sprintf_s(ZoneNameMid, "%sCage%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 11:
                    sprintf_s(ZoneNameMid, "%sRingFF%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 13:
                    sprintf_s(ZoneNameMid, "%sCageFF%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
            }
        }

        if (LastCrtPtNum >= 0)
        {
            switch (LastCrtPtType)
            {
                case -3:
                    sprintf_s(ZoneName, "%sAtom%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case -2:
                    sprintf_s(ZoneName, "%sMax%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case -1:
                    sprintf_s(ZoneName, "%sBond%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 0:
                    sprintf_s(ZoneName, "%sSaddle%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 1:
                    sprintf_s(ZoneName, "%sRing%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 2:
                    sprintf_s(ZoneName, "%sMin%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 3:
                    sprintf_s(ZoneName, "%sCage%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 11:
                    sprintf_s(ZoneName, "%sRingFF%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 13:
                    sprintf_s(ZoneName, "%sCageFF%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
            }
        }


        // Create Zone
        if (JDim > 1)
        {
            // TecUtilLockStart(AddOnID);

            // AveZoneNum = TUZoneGetNumByName(ZoneName);

            /* Set FieldDataType of all variables of CP zone */
            FieldDataType_e *VarDataType = ALLOC_ARRAY(NumVars, FieldDataType_e, "VarDataType");

            for (EntIndex_t iv = 0; iv < NumVars; iv++)
                VarDataType[iv] = FieldDataType_Float;
            VarDataType[TypeVarNum-1] = FieldDataType_Int16;

            /* Create zone if it doesn't already exist. */
            if (SurfaceZoneNum == TECUTILSETNOTMEMBER)
            {
                ArgList_pa ArgList;
                // TecUtilLockStart(AddOnID);
                ArgList = TecUtilArgListAlloc();
                TecUtilArgListAppendString(ArgList, SV_NAME, ZoneName);

                TecUtilArgListAppendInt(ArgList, SV_ZONETYPE,
                                        (ArbParam_t)ZoneType_Ordered);

                TecUtilArgListAppendInt(ArgList, SV_IMAX, IDim);
                TecUtilArgListAppendInt(ArgList, SV_JMAX, JDim);
                TecUtilArgListAppendInt(ArgList, SV_KMAX, 1);

                TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)VarDataType);

                IsOk = TecUtilDataSetAddZoneX(ArgList);

                TecUtilArgListDealloc(&ArgList);
                // TecUtilLockFinish(AddOnID);
                if (IsOk)
                {
                    NumZones++;
                    SurfaceZoneNum = NumZones;
                }
            }
            if ( VarDataType != NULL )
                FREE_ARRAY(VarDataType, "VarDataType");
        }

        /* Get FieldData pointers */
        if (IsOk)
        {
            LgIndex_t IJunk, JJunk, KJunk;
            TecUtilZoneGetInfo(SurfaceZoneNum, &IJunk, &JJunk, &KJunk,  &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            CDVarFDPtr = TecUtilDataValueGetWritableRef(SurfaceZoneNum, ChrgDensVarNum);
            TypeVarFDPtr = TecUtilDataValueGetWritableRef(SurfaceZoneNum, TypeVarNum);
        }

        /* Set the values for the zone */
        if (IsOk)
        {
            LgIndex_t ii, jj;
            for (jj = 0; IsOk && jj < JDim; jj++)
            {
                GradPath_pa GradPath = NULL;
                ArrListItem_u Item;

                Item = ArrListGetItem(GradPathList, jj);
                GradPath = (GradPath_pa)Item.VoidPtr;

                if (GradPathGetCount(GradPath) != IMax)
                    IsOk = GradPathResample(&GradPath, IMax);

                for (ii = 0; IsOk && ii < IMax; ii++)
                {
                    double X, Y, Z, Rho;

                    LgIndex_t FDOffset = ii + 1 + jj * IMax;

                    IsOk = GradPathGetPoint(GradPath, ii, &X, &Y, &Z, &Rho);
                    TecUtilDataValueSetByRef(XVarFDPtr, FDOffset, X);
                    TecUtilDataValueSetByRef(YVarFDPtr, FDOffset, Y);
                    TecUtilDataValueSetByRef(ZVarFDPtr, FDOffset, Z);

                    TecUtilDataValueSetByRef(CDVarFDPtr, FDOffset, Rho);

                    TecUtilDataValueSetByRef(TypeVarFDPtr, FDOffset, (double)SurfaceType);
                }
            }
        }
    }


    /* Set Zone aux data to completely transfer necessary information to Tecplot */
    if (IsOk)
    {
        AuxData_pa ZnAuxDataRef = TecUtilAuxDataZoneGetRef(SurfaceZoneNum);
        if (ZnAuxDataRef != NULL)
        {
            char CrtPtNumStr[200];
            char CrtPtTypeStr[200];
            char BaseZoneNumStr[200];

            sprintf_s(CrtPtNumStr, "%d", BegCrtPtNum + 1);
            IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BegCrtPtNum", (ArbParam_t)CrtPtNumStr,
                                         AuxDataType_String, TRUE);

            if (IsOk)
            {
                char CrtPtFFStr[200];
                Boolean_t IsFarField = FALSE;
                if (FirstCrtPtType > 10) IsFarField = TRUE;

                if (!IsFarField)
                    sprintf_s(CrtPtFFStr, "FALSE");
                else
                    sprintf_s(CrtPtFFStr, "TRUE");

                IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.IsFarField",
                                             (ArbParam_t)CrtPtFFStr, AuxDataType_String, TRUE);
            }

            IsOk = CritPointGetTypeString(BegCrtPtNum, BegCrtPtType, CrtPtTypeStr);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BegCrtPtType",
                                                       (ArbParam_t)CrtPtTypeStr, AuxDataType_String, TRUE);

            sprintf_s(CrtPtNumStr, "%d", FirstCrtPtNum + 1);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.EndCrtPtNum",
                                                       (ArbParam_t)CrtPtNumStr, AuxDataType_String, TRUE);

            IsOk = CritPointGetTypeString(FirstCrtPtNum, FirstCrtPtType, CrtPtTypeStr);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.EndCrtPtType",
                                                       (ArbParam_t)CrtPtTypeStr, AuxDataType_String, TRUE);

            sprintf_s(CrtPtNumStr, "%d", LastCrtPtNum + 1);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.ThrdCrtPtNum",
                                                       (ArbParam_t)CrtPtNumStr, AuxDataType_String, TRUE);

            IsOk = CritPointGetTypeString(LastCrtPtNum, LastCrtPtType, CrtPtTypeStr);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.ThrdCrtPtType",
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


    ENSURE(SurfaceZoneNum == TECUTILSETNOTMEMBER || SurfaceZoneNum > 0);
    return SurfaceZoneNum;
}

/**
 * Create a three CP (example, Atom-Ring-Cage) zero-flux connector surface
 * zone from an array list of GradPath pointers originating at Atom.
 *
 * param BaseZoneNum
 *     Zone number upon which topology (CritPoints) is based
 * param ChrgDensVarNum, TypeVarNum
 *     Variable numbers for electron charge density and CP type
 * param CritPoints
 *     Critical points structure containing CPs connected by surface
 * param GradPathList
 *     ArrList of pointers to a series of GradPaths which will make
 *     up the I-lines of a 2D ordered zone.
 * param IMax
 *     NumPoints in GradPath's (after resampling). Also the I-dimension
 *     of the 2D ordered zone.
 *
 * return
 *     Surface Tecplot zone number if it works, zero otherwise.
 */
EntIndex_t CreateZFSurfZoneFromGPList(EntIndex_t        ChrgDensVarNum,
                                      EntIndex_t        TypeVarNum,
                                      CritPoints_pa     CritPoints,
                                      ArrList_pa        GradPathList,
                                      LgIndex_t         IMax)
{
    Boolean_t   IsOk = TRUE;
    EntIndex_t  SurfaceZoneNum = TECUTILSETNOTMEMBER;
    EntIndex_t  SurfaceType;
    EntIndex_t  NumZones = 0;
    EntIndex_t  NumVars  = 0;
    LgIndex_t   IDim, JDim;
    char        ZoneName[200];
    LgIndex_t   BegCrtPtNum,        FirstCrtPtNum,        LastCrtPtNum;
    char        BegCrtPtType,       FirstCrtPtType,       LastCrtPtType;
    LgIndex_t   BegCrtPtTypeOffset, FirstCrtPtTypeOffset, LastCrtPtTypeOffset;
    EntIndex_t  BaseZoneNum;

    FieldData_pa XVarFDPtr = NULL;
    FieldData_pa YVarFDPtr = NULL;
    FieldData_pa ZVarFDPtr = NULL;
    FieldData_pa CDVarFDPtr = NULL;
    FieldData_pa TypeVarFDPtr = NULL;

    TecUtilLockStart(AddOnID);

    REQUIRE(TecUtilDataSetIsAvailable());

    /* Get num zones in dataset */
    if (IsOk)
        IsOk = TecUtilDataSetGetInfo(NULL, &NumZones, &NumVars);

    REQUIRE(ChrgDensVarNum > 0 && ChrgDensVarNum <= NumVars);
    REQUIRE(TypeVarNum > 0 && TypeVarNum <= NumVars);
    REQUIRE(CritPointsIsValid(CritPoints));
    REQUIRE(ArrListIsValid(GradPathList));
    REQUIRE(IMax > 1);

    BaseZoneNum = CritPoints->SourceZoneNum;
    REQUIRE(BaseZoneNum > 0 && BaseZoneNum <= NumZones);

    // Determine dimensions of the Tecplot surface zone and create it
    if (IsOk)
    {
        IDim = IMax;
        JDim = ArrListGetCount(GradPathList);
    }

    // Determine Zone Name
    if (IsOk)
    {
        char          ZoneNameBeg[200];
        char          ZoneNameMid[200];
        ArrListItem_u Item;
        GradPath_pa   GradPathFirst = NULL;
        GradPath_pa   GradPathLast = NULL;

        Item = ArrListGetItem(GradPathList, 0);
        GradPathFirst = (GradPath_pa)Item.VoidPtr;

        Item = ArrListGetItem(GradPathList, JDim - 1);
        GradPathLast = (GradPath_pa)Item.VoidPtr;


        BegCrtPtNum = GradPathFirst->BeginCrtPtNum;
        // CHECK(BegCrtPtNum == GradPathLast->BeginCrtPtNum);
        BegCrtPtType = CritPointsGetType(CritPoints, BegCrtPtNum);
        BegCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, BegCrtPtType);

        FirstCrtPtNum  = GradPathFirst->EndCrtPtNum;
        if (FirstCrtPtNum >= 0)
        {
            FirstCrtPtType = CritPointsGetType(CritPoints, FirstCrtPtNum);
            FirstCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, FirstCrtPtType);
        }

        LastCrtPtNum   = GradPathLast->EndCrtPtNum;
        if (LastCrtPtNum >= 0)
        {
            LastCrtPtType  = CritPointsGetType(CritPoints, LastCrtPtNum);
            LastCrtPtTypeOffset = CritPointsGetBegOffset(CritPoints, LastCrtPtType);
        }

        /* Determine SurfaceType = sum of component CritPoint types
         *   4=Ring-FarField (not sum)
         *   3=Bond-Ring-Cage
         *   1=Ring-Atom-Cage
         *  -1=Bond-Atom-Cage
         *  -3=Ring-Bond-Atom
         *  -4=Bond-FarField (not sum)
         */
        SurfaceType = BegCrtPtType + FirstCrtPtType + LastCrtPtType;

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

        if (FirstCrtPtNum >= 0)
        {
            switch (FirstCrtPtType)
            {
                case -3:
                    sprintf_s(ZoneNameMid, "%sAtom%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case -2:
                    sprintf_s(ZoneNameMid, "%sMax%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case -1:
                    sprintf_s(ZoneNameMid, "%sBond%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 0:
                    sprintf_s(ZoneNameMid, "%sSaddle%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 1:
                    sprintf_s(ZoneNameMid, "%sRing%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 2:
                    sprintf_s(ZoneNameMid, "%sMin%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 3:
                    sprintf_s(ZoneNameMid, "%sCage%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 11:
                    sprintf_s(ZoneNameMid, "%sRingFF%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
                case 13:
                    sprintf_s(ZoneNameMid, "%sCageFF%d_", ZoneNameBeg, FirstCrtPtNum - FirstCrtPtTypeOffset + 1);
                    break;
            }
        }

        if (LastCrtPtNum >= 0)
        {
            switch (LastCrtPtType)
            {
                case -3:
                    sprintf_s(ZoneName, "%sAtom%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case -2:
                    sprintf_s(ZoneName, "%sMax%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case -1:
                    sprintf_s(ZoneName, "%sBond%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 0:
                    sprintf_s(ZoneName, "%sSaddle%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 1:
                    sprintf_s(ZoneName, "%sRing%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 2:
                    sprintf_s(ZoneName, "%sMin%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 3:
                    sprintf_s(ZoneName, "%sCage%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 11:
                    sprintf_s(ZoneName, "%sRingFF%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
                case 13:
                    sprintf_s(ZoneName, "%sCageFF%d_Zone%d", ZoneNameMid, LastCrtPtNum - LastCrtPtTypeOffset + 1, BaseZoneNum);
                    break;
            }
        }


        // Create Zone
        if (JDim > 1)
        {
            // TecUtilLockStart(AddOnID);

            // AveZoneNum = TUZoneGetNumByName(ZoneName);

            /* Set FieldDataType of all variables of CP zone */
            FieldDataType_e *VarDataType = ALLOC_ARRAY(NumVars, FieldDataType_e, "VarDataType");

            for (EntIndex_t iv = 0; iv < NumVars; iv++)
                VarDataType[iv] = FieldDataType_Float;
            VarDataType[TypeVarNum-1] = FieldDataType_Int16;

            /* Create zone if it doesn't already exist. */
            if (SurfaceZoneNum == TECUTILSETNOTMEMBER)
            {
                ArgList_pa ArgList;
                // TecUtilLockStart(AddOnID);
                ArgList = TecUtilArgListAlloc();
                TecUtilArgListAppendString(ArgList, SV_NAME, ZoneName);

                TecUtilArgListAppendInt(ArgList, SV_ZONETYPE,
                                        (ArbParam_t)ZoneType_Ordered);

                TecUtilArgListAppendInt(ArgList, SV_IMAX, IDim);
                TecUtilArgListAppendInt(ArgList, SV_JMAX, JDim);
                TecUtilArgListAppendInt(ArgList, SV_KMAX, 1);

                TecUtilArgListAppendArray(ArgList, SV_VARDATATYPE, (void *)VarDataType);

                IsOk = TecUtilDataSetAddZoneX(ArgList);

                TecUtilArgListDealloc(&ArgList);
                // TecUtilLockFinish(AddOnID);
                if (IsOk)
                {
                    NumZones++;
                    SurfaceZoneNum = NumZones;
                }
            }
            if ( VarDataType != NULL )
                FREE_ARRAY(VarDataType, "VarDataType");
        }

        /* Get FieldData pointers */
        if (IsOk)
        {
            LgIndex_t IJunk, JJunk, KJunk;
            TecUtilZoneGetInfo(SurfaceZoneNum, &IJunk, &JJunk, &KJunk,  &XVarFDPtr, &YVarFDPtr, &ZVarFDPtr,
                               NULL, NULL, NULL, NULL, NULL, NULL, NULL);
            CDVarFDPtr = TecUtilDataValueGetWritableRef(SurfaceZoneNum, ChrgDensVarNum);
            TypeVarFDPtr = TecUtilDataValueGetWritableRef(SurfaceZoneNum, TypeVarNum);
        }

        /* Set the values for the zone */
        if (IsOk)
        {
            LgIndex_t ii, jj;
            for (jj = 0; IsOk && jj < JDim; jj++)
            {
                GradPath_pa GradPath = NULL;
                ArrListItem_u Item;

                Item = ArrListGetItem(GradPathList, jj);
                GradPath = (GradPath_pa)Item.VoidPtr;

                if (GradPathGetCount(GradPath) != IMax)
                    IsOk = GradPathResample(&GradPath, IMax);

                for (ii = 0; IsOk && ii < IMax; ii++)
                {
                    double X, Y, Z, Rho;

                    LgIndex_t FDOffset = ii + 1 + jj * IMax;

                    IsOk = GradPathGetPoint(GradPath, ii, &X, &Y, &Z, &Rho);
                    TecUtilDataValueSetByRef(XVarFDPtr, FDOffset, X);
                    TecUtilDataValueSetByRef(YVarFDPtr, FDOffset, Y);
                    TecUtilDataValueSetByRef(ZVarFDPtr, FDOffset, Z);

                    TecUtilDataValueSetByRef(CDVarFDPtr, FDOffset, Rho);

                    TecUtilDataValueSetByRef(TypeVarFDPtr, FDOffset, (double)SurfaceType);
                }
            }
        }
    }


    /* Set Zone aux data to completely transfer necessary information to Tecplot */
    if (IsOk)
    {
        AuxData_pa ZnAuxDataRef = TecUtilAuxDataZoneGetRef(SurfaceZoneNum);
        if (ZnAuxDataRef != NULL)
        {
            char CrtPtNumStr[200];
            char CrtPtTypeStr[200];
            char BaseZoneNumStr[200];

            sprintf_s(CrtPtNumStr, "%d", BegCrtPtNum + 1);
            IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BegCrtPtNum", (ArbParam_t)CrtPtNumStr,
                                         AuxDataType_String, TRUE);

            if (IsOk)
            {
                char CrtPtFFStr[200];
                Boolean_t IsFarField = FALSE;
                if (FirstCrtPtType > 10) IsFarField = TRUE;

                if (!IsFarField)
                    sprintf_s(CrtPtFFStr, "FALSE");
                else
                    sprintf_s(CrtPtFFStr, "TRUE");

                IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.IsFarField",
                                             (ArbParam_t)CrtPtFFStr, AuxDataType_String, TRUE);
            }

            IsOk = CritPointGetTypeString(BegCrtPtNum, BegCrtPtType, CrtPtTypeStr);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.BegCrtPtType",
                                                       (ArbParam_t)CrtPtTypeStr, AuxDataType_String, TRUE);

            sprintf_s(CrtPtNumStr, "%d", FirstCrtPtNum + 1);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.EndCrtPtNum",
                                                       (ArbParam_t)CrtPtNumStr, AuxDataType_String, TRUE);

            IsOk = CritPointGetTypeString(FirstCrtPtNum, FirstCrtPtType, CrtPtTypeStr);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.EndCrtPtType",
                                                       (ArbParam_t)CrtPtTypeStr, AuxDataType_String, TRUE);

            sprintf_s(CrtPtNumStr, "%d", LastCrtPtNum + 1);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.ThrdCrtPtNum",
                                                       (ArbParam_t)CrtPtNumStr, AuxDataType_String, TRUE);

            IsOk = CritPointGetTypeString(LastCrtPtNum, LastCrtPtType, CrtPtTypeStr);
            if (IsOk) IsOk = TecUtilAuxDataSetItem(ZnAuxDataRef, "CompChem.ThrdCrtPtType",
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


    ENSURE(SurfaceZoneNum == TECUTILSETNOTMEMBER || SurfaceZoneNum > 0);
    return SurfaceZoneNum;
}



