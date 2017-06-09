
#include "TECADDON.h"
#include "ADDGLBL.h"
#include "GUIDEFS.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

#include <string>
#include <vector>

#include "ZONEVARINFO.h"
#include "CSM_DATA_TYPES.h"

#include "STREAMTRACEOPERATIONS.h"

#include <armadillo>
using namespace arma;

using std::string;
using std::stringstream;
using std::vector;




Boolean_t StreamtraceResample(EntIndex_t STZoneNum, LgIndex_t NumSTPoints){

	TecUtilLockStart(AddOnID);

	Boolean_t IsOk = TecUtilDataSetAddZone("Resampled Streamtrace", NumSTPoints, 1, 1, ZoneType_Ordered, NULL);

	EntIndex_t TmpZoneNum = TecUtilDataSetGetNumZones();
	if (IsOk){
		Set_pa NewZoneSet = TecUtilSetAlloc(FALSE);
		IsOk = TecUtilSetAddMember(NewZoneSet, TmpZoneNum, FALSE);
		if (IsOk){
			TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)NewZoneSet);
		}
		TecUtilSetDealloc(&NewZoneSet);
	}
	EntIndex_t NumVars = TecUtilDataSetGetNumVars();
	EntIndex_t XYZVarNums[3] = { -1 };
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	IsOk = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2]);
	if (!IsOk){
		TecUtilDialogErrMsg("Couldn't get XYZ var nums");
		TecUtilLockFinish(AddOnID);
		return IsOk;
	}
	EntIndex_t VolZoneNum = ZoneNumByName(string("Full Volume"));
	vec3 MinVolPt;
	IsOk = (VolZoneNum > 0);
	if (IsOk){
		double DJunk;
		for (int i = 0; i < 3; ++i){
			TecUtilDataValueGetMinMaxByZoneVar(VolZoneNum, XYZVarNums[i], &MinVolPt[i], &DJunk);
		}
	}

	/*
	*	Get get and set functions for direct access to streamtrace
	*	zone and new zone's data.
	*	This looks overcomplicated, but it's to ensure that the XYZ
	*	data pointers are the first three in the vectors.
	*/
	vector<FieldData_pa> STGetPtrs(NumVars, NULL);
	vector<FieldValueGetFunction_pf> STGetFuncts(NumVars, NULL);
	vector<FieldData_pa> STSetPtrs(NumVars, NULL);
	vector<FieldValueSetFunction_pf> STSetFuncts(NumVars, NULL);
	if (IsOk){
		for (int i = 0; i < 3 && IsOk; ++i){
			STGetPtrs[i] = TecUtilDataValueGetReadableRef(STZoneNum, XYZVarNums[i]);
			STSetPtrs[i] = TecUtilDataValueGetWritableRef(TmpZoneNum, XYZVarNums[i]);
			IsOk = (VALID_REF(STGetPtrs[i]) && VALID_REF(STSetPtrs[i]));
			if (IsOk){
				STGetFuncts[i] = TecUtilDataValueRefGetGetFunc(STGetPtrs[i]);
				STSetFuncts[i] = TecUtilDataValueRefGetSetFunc(STSetPtrs[i]);
				IsOk = (VALID_REF(STGetFuncts[i]) && VALID_REF(STSetFuncts[i]));
			}
		}
		EntIndex_t VarNum = 3;
		for (int i = 1; i <= NumVars && IsOk; ++i){
			Boolean_t IsXYZVarNum = FALSE;
			for (int j = 0; j < 3 && !IsXYZVarNum; ++j){
				IsXYZVarNum = (XYZVarNums[j] == i);
			}
			if (!IsXYZVarNum){
				STGetPtrs[VarNum] = TecUtilDataValueGetReadableRef(STZoneNum, i);
				STSetPtrs[VarNum] = TecUtilDataValueGetWritableRef(TmpZoneNum, i);
				IsOk = (VALID_REF(STGetPtrs[VarNum]) && VALID_REF(STSetPtrs[VarNum]));
				if (IsOk){
					STGetFuncts[VarNum] = TecUtilDataValueRefGetGetFunc(STGetPtrs[VarNum]);
					STSetFuncts[VarNum] = TecUtilDataValueRefGetSetFunc(STSetPtrs[VarNum]);
					IsOk = (VALID_REF(STGetFuncts[VarNum]) && VALID_REF(STSetFuncts[VarNum]));
				}
				VarNum++;
			}
		}
		if (!IsOk){
			TecUtilDialogErrMsg("Failed to get streamtrace zone pointers for resampling");
			TecUtilLockFinish(AddOnID);
			return IsOk;
		}
	}

	/*
	*	Get length of streamtrace zone to be resampled and
	*	distance between each point of resampled zone.
	*	Here we also find the REAL max I value for the streamtrace,
	*	since it usually stops moving well before it reaches its
	*	max number of points.
	*/
	TecUtilDataLoadBegin();
	double STLength = 0.0, DelLength;
	LgIndex_t STIJKMax[3];
	if (IsOk){
		TecUtilZoneGetIJK(STZoneNum, &STIJKMax[0], &STIJKMax[1], &STIJKMax[2]);
		IsOk = (STIJKMax[0] > 1 && STIJKMax[1] == 1 && STIJKMax[2] == 1);
	}
	if (IsOk){
		vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;
		for (int i = 0; i < 3; ++i)
			PtIm1[i] = STGetFuncts[i](STGetPtrs[i], 0);
		Boolean_t ValuesStartedChanging = FALSE;
		for (LgIndex_t i = 1; i < STIJKMax[0]; ++i){
			for (int j = 0; j < 3; ++j)
				PtI[j] = STGetFuncts[j](STGetPtrs[j], i);
			STLength += Distance(PtI, PtIm1);
			if (!ValuesStartedChanging && sum(PtIm1 != PtI) == 3){
				ValuesStartedChanging = TRUE;
			}
			else if (sum(PtIm1 == PtI) == 3 && ValuesStartedChanging){
				STIJKMax[0] = i - 1;
				break;
			}
			PtIm1 = PtI;
		}
		DelLength = STLength / double(NumSTPoints);
	}
	TecUtilDataLoadEnd();

	/*
	*	Now step down streamtrace zone again, adding a point to the
	*	resampled zone each time a distance of DelLength is reached,
	*	interpolating values of all variables each time. First and
	*	last points will always be included.
	*/
	if (IsOk){
		vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;
		vector<double> ValsI(NumVars, MinVolPt[0] - 10);
		vector<double> ValsIm1(NumVars, MinVolPt[0] - 10);
		double ArcLength = 0.0,
			ArcLengthI = 0.0,
			ArcLengthIm1 = 0.0;

		//	Set first point and get values needed for first step in loop
		for (int i = 0; i < NumVars; ++i){
			ValsI[i] = STGetFuncts[i](STGetPtrs[i], 0);
			STSetFuncts[i](STSetPtrs[i], 0, ValsI[i]);
			if (i < 3){
				PtI[i] = ValsI[i];
			}
		}
		PtIm1 = PtI;
		ValsIm1 = ValsI;
		LgIndex_t OldI = 0;

		TecUtilDataLoadBegin();

		for (int NewI = 1; NewI < NumSTPoints; ++NewI){
			ArcLength += DelLength;

			/*
			*	Move along old path until a segment of DelLength has
			*	been traversed, then interpolate to get new values to
			*	add to new streamtrace zone
			*/
			while (OldI < STIJKMax[0] && ArcLengthI < ArcLength){
				OldI++;

				ArcLengthIm1 = ArcLengthI;
				PtIm1 = PtI;
				ValsIm1 = ValsI;

				for (int i = 0; i < NumVars; ++i){
					ValsI[i] = STGetFuncts[i](STGetPtrs[i], OldI);
					if (i < 3){
						PtI[i] = ValsI[i];
					}
				}

				ArcLengthI += Distance(PtI, PtIm1);
			}

			/*
			*	DelLength has been reached, so time to add a new point
			*	to the new streamtrace zone
			*/
			double Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
			for (int i = 0; i < NumVars; ++i){
				STSetFuncts[i](STSetPtrs[i], NewI,
					ValsIm1[i] + Ratio * (ValsI[i] - ValsIm1[i])
					);
			}
		}
		TecUtilDataLoadEnd();
		//	Set last point

		// 					for (int i = 0; i < NumVars; ++i){
		// 						ValsI[i] = STGetFuncts[i](STGetPtrs[i], OldI);
		// 						STSetFuncts[i](STSetPtrs[i], NumSTPoints,
		// 							ValsI[i]
		// 							);
		// 					}
	}
	Set_pa DeleteZoneSet = TecUtilSetAlloc(FALSE);
	IsOk = TecUtilSetAddMember(DeleteZoneSet, STZoneNum, FALSE);
	if (IsOk){
		IsOk = TecUtilDataSetDeleteZone(DeleteZoneSet);
// 		if (IsOk){
// 			TecUtilStateChanged(StateChange_ZonesDeleted, (ArbParam_t)DeleteZoneSet);
// 		}
		TecUtilSetDealloc(&DeleteZoneSet);
	}
	else{
		TecUtilDialogErrMsg("Failed to delete streamtrace zone");
		TecUtilLockFinish(AddOnID);
		return IsOk;
	}
	

	TecUtilLockFinish(AddOnID);

	return IsOk;
}

Boolean_t StreamtraceResampleToExistingZone(EntIndex_t STZoneNum, 
	EntIndex_t NewZoneNum, 
	LgIndex_t NumSTPoints, 
	EntIndex_t NumVars,
	EntIndex_t* XYZVarNums, 
	vec3 & MinVolPt)
{

	TecUtilLockStart(AddOnID);

	Boolean_t IsOk = TRUE;

	/*
	*	Get get and set functions for direct access to streamtrace
	*	zone and new zone's data.
	*	This looks overcomplicated, but it's to ensure that the XYZ
	*	data pointers are the first three in the vectors.
	*/
	vector<FieldData_pa> STGetPtrs(NumVars, NULL);
	vector<FieldValueGetFunction_pf> STGetFuncts(NumVars, NULL);
	vector<FieldData_pa> STSetPtrs(NumVars, NULL);
	vector<FieldValueSetFunction_pf> STSetFuncts(NumVars, NULL);
	if (IsOk){
		for (int i = 0; i < 3 && IsOk; ++i){
			STGetPtrs[i] = TecUtilDataValueGetReadableRef(STZoneNum, XYZVarNums[i]);
			STSetPtrs[i] = TecUtilDataValueGetWritableRef(NewZoneNum, XYZVarNums[i]);
			IsOk = (VALID_REF(STGetPtrs[i]) && VALID_REF(STSetPtrs[i]));
			if (IsOk){
				STGetFuncts[i] = TecUtilDataValueRefGetGetFunc(STGetPtrs[i]);
				STSetFuncts[i] = TecUtilDataValueRefGetSetFunc(STSetPtrs[i]);
				IsOk = (VALID_REF(STGetFuncts[i]) && VALID_REF(STSetFuncts[i]));
			}
		}
		EntIndex_t VarNum = 3;
		for (int i = 1; i <= NumVars && IsOk; ++i){
			Boolean_t IsXYZVarNum = FALSE;
			for (int j = 0; j < 3 && !IsXYZVarNum; ++j){
				IsXYZVarNum = (XYZVarNums[j] == i);
			}
			if (!IsXYZVarNum){
				STGetPtrs[VarNum] = TecUtilDataValueGetReadableRef(STZoneNum, i);
				STSetPtrs[VarNum] = TecUtilDataValueGetWritableRef(NewZoneNum, i);
				IsOk = (VALID_REF(STGetPtrs[VarNum]) && VALID_REF(STSetPtrs[VarNum]));
				if (IsOk){
					STGetFuncts[VarNum] = TecUtilDataValueRefGetGetFunc(STGetPtrs[VarNum]);
					STSetFuncts[VarNum] = TecUtilDataValueRefGetSetFunc(STSetPtrs[VarNum]);
					IsOk = (VALID_REF(STGetFuncts[VarNum]) && VALID_REF(STSetFuncts[VarNum]));
				}
				VarNum++;
			}
		}
		if (!IsOk){
			TecUtilDialogErrMsg("Failed to get streamtrace zone pointers for resampling");
			TecUtilLockFinish(AddOnID);
			return IsOk;
		}
	}

	/*
	*	Get length of streamtrace zone to be resampled and
	*	distance between each point of resampled zone.
	*	Here we also find the REAL max I value for the streamtrace,
	*	since it usually stops moving well before it reaches its
	*	max number of points.
	*/
	double STLength = 0.0, DelLength;
	LgIndex_t STIJKMax[3];
	if (IsOk){
		TecUtilZoneGetIJK(STZoneNum, &STIJKMax[0], &STIJKMax[1], &STIJKMax[2]);
		IsOk = (STIJKMax[0] > 1 && STIJKMax[1] == 1 && STIJKMax[2] == 1);
	}
	if (IsOk){
		vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;
		for (int i = 0; i < 3; ++i)
			PtIm1[i] = STGetFuncts[i](STGetPtrs[i], 0);
		Boolean_t ValuesStartedChanging = FALSE;
		for (LgIndex_t i = 1; i < STIJKMax[0]; ++i){
			for (int j = 0; j < 3; ++j)
				PtI[j] = STGetFuncts[j](STGetPtrs[j], i);
			STLength += Distance(PtI, PtIm1);
			if (!ValuesStartedChanging && sum(PtIm1 != PtI) == 3){
				ValuesStartedChanging = TRUE;
			}
			else if (sum(PtIm1 == PtI) == 3 && ValuesStartedChanging){
				STIJKMax[0] = i - 1;
				break;
			}
			PtIm1 = PtI;
		}
		DelLength = STLength / double(NumSTPoints);
	}

	/*
	*	Now step down streamtrace zone again, adding a point to the
	*	resampled zone each time a distance of DelLength is reached,
	*	interpolating values of all variables each time. First and
	*	last points will always be included.
	*/
	if (IsOk){
		vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;
		vector<double> ValsI(NumVars, MinVolPt[0] - 10);
		vector<double> ValsIm1(NumVars, MinVolPt[0] - 10);
		double ArcLength = 0.0,
			ArcLengthI = 0.0,
			ArcLengthIm1 = 0.0;

		//	Set first point and get values needed for first step in loop
		for (int i = 0; i < NumVars; ++i){
			ValsI[i] = STGetFuncts[i](STGetPtrs[i], 0);
			STSetFuncts[i](STSetPtrs[i], 0, ValsI[i]);
			if (i < 3){
				PtI[i] = ValsI[i];
			}
		}
		PtIm1 = PtI;
		ValsIm1 = ValsI;
		LgIndex_t OldI = 0;


		for (int NewI = 1; NewI < NumSTPoints; ++NewI){
			ArcLength += DelLength;

			/*
			*	Move along old path until a segment of DelLength has
			*	been traversed, then interpolate to get new values to
			*	add to new streamtrace zone
			*/
			while (OldI < STIJKMax[0] && ArcLengthI < ArcLength){
				OldI++;

				ArcLengthIm1 = ArcLengthI;
				PtIm1 = PtI;
				ValsIm1 = ValsI;

				for (int i = 0; i < NumVars; ++i){
					ValsI[i] = STGetFuncts[i](STGetPtrs[i], OldI);
					if (i < 3){
						PtI[i] = ValsI[i];
					}
				}

				ArcLengthI += Distance(PtI, PtIm1);
			}

			/*
			*	DelLength has been reached, so time to add a new point
			*	to the new streamtrace zone
			*/
			double Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
			for (int i = 0; i < NumVars; ++i){
				STSetFuncts[i](STSetPtrs[i], NewI,
					ValsIm1[i] + Ratio * (ValsI[i] - ValsIm1[i])
					);
			}
		}
		//	Set last point

		// 					for (int i = 0; i < NumVars; ++i){
		// 						ValsI[i] = STGetFuncts[i](STGetPtrs[i], OldI);
		// 						STSetFuncts[i](STSetPtrs[i], NumSTPoints,
		// 							ValsI[i]
		// 							);
		// 					}
	}

	TecUtilLockFinish(AddOnID);

	return IsOk;
}


Boolean_t StreamtraceConcatenateResample(EntIndex_t STZoneNum, EntIndex_t VolGPNum, LgIndex_t NumSTPoints){

	TecUtilLockStart(AddOnID);

	Boolean_t IsOk = TecUtilDataSetAddZone("Resampled Streamtrace", NumSTPoints, 1, 1, ZoneType_Ordered, NULL);
	

	EntIndex_t TmpZoneNum = TecUtilDataSetGetNumZones();
	if (IsOk){
		Set_pa NewZoneSet = TecUtilSetAlloc(FALSE);
		IsOk = TecUtilSetAddMember(NewZoneSet, TmpZoneNum, FALSE);
		if (IsOk){
			TecUtilStateChanged(StateChange_ZonesAdded, (ArbParam_t)NewZoneSet);
		}
		TecUtilSetDealloc(&NewZoneSet);
	}
	EntIndex_t NumVars = TecUtilDataSetGetNumVars();
	EntIndex_t XYZVarNums[3] = { -1 };
	TecUtilAxisGetVarAssignments(&XYZVarNums[0], &XYZVarNums[1], &XYZVarNums[2]);
	IsOk = (XYZVarNums[0] > 0 && XYZVarNums[1] > 0 && XYZVarNums[2]);
	if (!IsOk){
		TecUtilDialogErrMsg("Couldn't get XYZ var nums");
		TecUtilLockFinish(AddOnID);
		return IsOk;
	}
	EntIndex_t VolZoneNum = ZoneNumByName(string("Full Volume"));
	vec3 MinVolPt;
	IsOk = (VolZoneNum > 0);
	if (IsOk){
		double DJunk;
		for (int i = 0; i < 3; ++i){
			TecUtilDataValueGetMinMaxByZoneVar(VolZoneNum, XYZVarNums[i], &MinVolPt[i], &DJunk);
		}
	}

	/*
	*	Get get and set functions for direct access to streamtrace
	*	zone and new zone's data.
	*	This looks overcomplicated, but it's to ensure that the XYZ
	*	data pointers are the first three in the vectors.
	*/
	EntIndex_t GetSTZoneNums[2] = { STZoneNum, VolGPNum };
	vector<vector<FieldData_pa>> STGetPtrs(2, vector<FieldData_pa>(NumVars, NULL));
	vector<vector<FieldValueGetFunction_pf> > STGetFuncts(2, vector<FieldValueGetFunction_pf>(NumVars, NULL));
	vector<FieldData_pa> STSetPtrs(NumVars, NULL);
	vector<FieldValueSetFunction_pf> STSetFuncts(NumVars, NULL);
	if (IsOk){
		for (int i = 0; i < 3 && IsOk; ++i){
			for (int j = 0; j < 2 && IsOk; ++j){
				STGetPtrs[j][i] = TecUtilDataValueGetReadableRef(GetSTZoneNums[j], XYZVarNums[i]);
				IsOk = (VALID_REF(STGetPtrs[j][i]));
			}
			STSetPtrs[i] = TecUtilDataValueGetWritableRef(TmpZoneNum, XYZVarNums[i]);
			IsOk = (VALID_REF(STSetPtrs[i]));
			if (IsOk){
				for (int j = 0; j < 2 && IsOk; ++j){
					STGetFuncts[j][i] = TecUtilDataValueRefGetGetFunc(STGetPtrs[j][i]);
					IsOk = (VALID_REF(STGetFuncts[j][i]));
				}
				STSetFuncts[i] = TecUtilDataValueRefGetSetFunc(STSetPtrs[i]);
				IsOk = (VALID_REF(STSetFuncts[i]));
			}
		}
		EntIndex_t VarNum = 3;
		for (int i = 1; i <= NumVars && IsOk; ++i){
			Boolean_t IsXYZVarNum = FALSE;
			for (int j = 0; j < 3 && !IsXYZVarNum; ++j){
				IsXYZVarNum = (XYZVarNums[j] == i);
			}
			if (!IsXYZVarNum){
				for (int j = 0; j < 2 && IsOk; ++j){
					STGetPtrs[j][VarNum] = TecUtilDataValueGetReadableRef(GetSTZoneNums[j], i);
					IsOk = (VALID_REF(STGetPtrs[j][VarNum]));
				}
				STSetPtrs[VarNum] = TecUtilDataValueGetWritableRef(TmpZoneNum, i);
				IsOk = (VALID_REF(STSetPtrs[VarNum]));
				if (IsOk){
					for (int j = 0; j < 2 && IsOk; ++j){
						STGetFuncts[j][VarNum] = TecUtilDataValueRefGetGetFunc(STGetPtrs[j][VarNum]);
						IsOk = (VALID_REF(STGetFuncts[j][VarNum]));
					}
					STSetFuncts[VarNum] = TecUtilDataValueRefGetSetFunc(STSetPtrs[VarNum]);
					IsOk = (VALID_REF(STSetFuncts[VarNum]));
				}
				VarNum++;
			}
		}
		if (!IsOk){
			TecUtilDialogErrMsg("Failed to get streamtrace zone pointers for resampling");
			TecUtilLockFinish(AddOnID);
			return IsOk;
		}
	}

	/*
	*	Get length of streamtrace zones to be resampled and
	*	distance between each point of resampled zones.
	*	Here we also find the REAL max I value for the streamtraces,
	*	since it usually stops moving well before it reaches its
	*	max number of points.
	*/
	double DelLength;
	double STLength[2] = { 0.0 };
	LgIndex_t STIJKMax[2][3];

	TecUtilDataLoadBegin();

	for (int PassNum = 0; PassNum < 2; ++PassNum){
		if (IsOk){
			TecUtilZoneGetIJK(GetSTZoneNums[PassNum], &STIJKMax[PassNum][0], &STIJKMax[PassNum][1], &STIJKMax[PassNum][2]);
			IsOk = (STIJKMax[PassNum][0] > 1 && STIJKMax[PassNum][1] == 1 && STIJKMax[PassNum][2] == 1);
		}
		if (IsOk){
			vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;
			for (int i = 0; i < 3; ++i)
				PtIm1[i] = STGetFuncts[PassNum][i](STGetPtrs[PassNum][i], 0);
			Boolean_t ValuesStartedChanging = FALSE;
			for (LgIndex_t i = 1; i < STIJKMax[PassNum][0]; ++i){
				for (int j = 0; j < 3; ++j)
					PtI[j] = STGetFuncts[PassNum][j](STGetPtrs[PassNum][j], i);
				STLength[PassNum] += Distance(PtI, PtIm1);
				if (!ValuesStartedChanging && sum(PtIm1 != PtI) == 3){
					ValuesStartedChanging = TRUE;
				}
				else if (sum(PtIm1 == PtI) == 3 && ValuesStartedChanging){
					STIJKMax[PassNum][0] = i - 1;
					break;
				}
				PtIm1 = PtI;
			}
		}
	}
	TecUtilDataLoadEnd();
	/*
	*	Need to make sure that we start at the far end of one of the
	*	streamtraces so that when we concatenate the two, its one
	*	continuous line.
	*/


	Boolean_t StartFrom1[2] = { TRUE };
	int HalfNumSTPoints[2];
	if (IsOk){
		vec3 EndPoints[4];
		for (int i = 0; i < 2; ++i){
			for (int j = 0; j < 3; ++j){
				EndPoints[2 * i][j] = STGetFuncts[i][j](STGetPtrs[i][j], 0);
				EndPoints[1 + 2 * i][j] = STGetFuncts[i][j](STGetPtrs[i][j], STIJKMax[i][0] - 1);
			}
		}
		int GapNum = 1, MinGapNum = 1;
		double MinGap = 1e50;
		for (int i = 0; i < 2; ++i){
			for (int j = 2; j < 4; ++j){
				double TmpGap = DistSqr(EndPoints[i], EndPoints[j]);
				if (TmpGap < MinGap){
					MinGap = TmpGap;
					MinGapNum = GapNum;
				}
				GapNum++;
			}
		}
		double TotalLength = STLength[0] + STLength[1] + sqrt(MinGap);
		HalfNumSTPoints[0] = int(STLength[0] / TotalLength * (double)NumSTPoints);
		HalfNumSTPoints[1] = NumSTPoints - HalfNumSTPoints[0];
		DelLength = TotalLength / double(NumSTPoints);
		switch (MinGapNum){
			case 1:
				/*
				*	point 1 of the two streamtraces are next to eachother, so go in
				*	order last1 first1 first2 last2
				*/
				StartFrom1[0] = FALSE;
				StartFrom1[1] = TRUE;
				break;
			case 2:
				/*
				*	point 1 of the first streamtrace is next to last of second one, so go in
				*	order last1 first1 last2 first 2
				*/
				StartFrom1[0] = FALSE;
				StartFrom1[1] = FALSE;
				break;
			case 3:
				/*
				*	last point of first streamtrace next to first point of second, so go in
				*	order first1 last1 first2 last2
				*/
				StartFrom1[0] = TRUE;
				StartFrom1[1] = TRUE;
				break;
			default:
				/*
				*	last point of first streamtrace next to last point of second, so go in
				*	order first1 last1 last2 first2
				*/
				StartFrom1[0] = TRUE;
				StartFrom1[1] = FALSE;
				break;
		}
	}

	/*
	*	Now step down streamtrace zone again, adding a point to the
	*	resampled zone each time a distance of DelLength is reached,
	*	interpolating values of all variables each time. First and
	*	last points will always be included.
	*/
	if (IsOk){
		vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;
		vector<double> ValsI(NumVars, MinVolPt[0] - 10);
		vector<double> ValsIm1(NumVars, MinVolPt[0] - 10);
		double ArcLength = 0.0,
			ArcLengthI = 0.0,
			ArcLengthIm1 = 0.0;
		int NewI = 0;

		TecUtilDataLoadBegin();

		for (int PassNum = 0; PassNum < 2; ++PassNum){
			int OldI = 0;
			int OldStep = 1;
			if (!StartFrom1[PassNum]){
				OldI = STIJKMax[PassNum][0] - 1;
				OldStep = -1;
			}

			//	Set first point and get values needed for first step in loop
			for (int i = 0; i < NumVars; ++i){
				ValsI[i] = STGetFuncts[PassNum][i](STGetPtrs[PassNum][i], OldI);
				STSetFuncts[i](STSetPtrs[i], NewI, ValsI[i]);
				if (i < 3){
					PtI[i] = ValsI[i];
				}
			}
			NewI++;
			OldI += OldStep;
			int NewIEnd = HalfNumSTPoints[PassNum];
			if (PassNum > 0){
				NewIEnd += HalfNumSTPoints[0];
				ArcLengthI += Distance(PtI, PtIm1);
			}
			PtIm1 = PtI;
			ValsIm1 = ValsI;
			for (NewI = NewI; NewI < NewIEnd; ++NewI){
				ArcLength += DelLength;

				/*
				*	Move along old path until a segment of DelLength has
				*	been traversed, then interpolate to get new values to
				*	add to new streamtrace zone
				*/
				while (OldI < STIJKMax[PassNum][0] && OldI > 0 && ArcLengthI < ArcLength){
					ArcLengthIm1 = ArcLengthI;
					PtIm1 = PtI;
					ValsIm1 = ValsI;

					for (int i = 0; i < NumVars; ++i){
						ValsI[i] = STGetFuncts[PassNum][i](STGetPtrs[PassNum][i], OldI);
						if (i < 3){
							PtI[i] = ValsI[i];
						}
					}

					ArcLengthI += Distance(PtI, PtIm1);

					OldI += OldStep;
				}

				/*
				*	DelLength has been reached, so time to add a new point
				*	to the new streamtrace zone
				*/
				double Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
				for (int i = 0; i < NumVars; ++i){
					STSetFuncts[i](STSetPtrs[i], NewI,
						ValsIm1[i] + Ratio * (ValsI[i] - ValsIm1[i])
						);
				}
			}
			//	Set last point

			// 					for (int i = 0; i < NumVars; ++i){
			// 						ValsI[i] = STGetFuncts[i](STGetPtrs[i], OldI);
			// 						STSetFuncts[i](STSetPtrs[i], NumSTPoints,
			// 							ValsI[i]
			// 							);
			// 					}
		}

		TecUtilDataLoadEnd();

		Set_pa DeleteZoneSet = TecUtilSetAlloc(FALSE);
		IsOk = TecUtilSetAddMember(DeleteZoneSet, STZoneNum, FALSE);
		if (IsOk){
			IsOk = TecUtilDataSetDeleteZone(DeleteZoneSet);
// 			if (IsOk){
// 				TecUtilStateChanged(StateChange_ZonesDeleted, (ArbParam_t)DeleteZoneSet);
// 			}
			TecUtilSetDealloc(&DeleteZoneSet);
		}
		else{
			TecUtilDialogErrMsg("Failed to delete streamtrace zone");
			TecUtilLockFinish(AddOnID);
			return IsOk;
		}
	}

	TecUtilLockFinish(AddOnID);
	return IsOk;
}

Boolean_t StreamtraceConcatenateResampleToExistingZone(EntIndex_t STZoneNum, 
	EntIndex_t NewZoneNum, 
	EntIndex_t VolGPNum, 
	LgIndex_t NumSTPoints, 
	EntIndex_t NumVars, 
	EntIndex_t* XYZVarNums, 
	vec3 & MinVolPt)
{

	TecUtilLockStart(AddOnID);

	Boolean_t IsOk = TRUE;

	/*
	*	Get get and set functions for direct access to streamtrace
	*	zone and new zone's data.
	*	This looks overcomplicated, but it's to ensure that the XYZ
	*	data pointers are the first three in the vectors.
	*/
	EntIndex_t GetSTZoneNums[2] = { STZoneNum, VolGPNum };
	vector<vector<FieldData_pa>> STGetPtrs(2, vector<FieldData_pa>(NumVars, NULL));
	vector<vector<FieldValueGetFunction_pf> > STGetFuncts(2, vector<FieldValueGetFunction_pf>(NumVars, NULL));
	vector<FieldData_pa> STSetPtrs(NumVars, NULL);
	vector<FieldValueSetFunction_pf> STSetFuncts(NumVars, NULL);
	if (IsOk){
		for (int i = 0; i < 3 && IsOk; ++i){
			for (int j = 0; j < 2 && IsOk; ++j){
				STGetPtrs[j][i] = TecUtilDataValueGetReadableRef(GetSTZoneNums[j], XYZVarNums[i]);
				IsOk = (VALID_REF(STGetPtrs[j][i]));
			}
			STSetPtrs[i] = TecUtilDataValueGetWritableRef(NewZoneNum, XYZVarNums[i]);
			IsOk = (VALID_REF(STSetPtrs[i]));
			if (IsOk){
				for (int j = 0; j < 2 && IsOk; ++j){
					STGetFuncts[j][i] = TecUtilDataValueRefGetGetFunc(STGetPtrs[j][i]);
					IsOk = (VALID_REF(STGetFuncts[j][i]));
				}
				STSetFuncts[i] = TecUtilDataValueRefGetSetFunc(STSetPtrs[i]);
				IsOk = (VALID_REF(STSetFuncts[i]));
			}
		}
		EntIndex_t VarNum = 3;
		for (int i = 1; i <= NumVars && IsOk; ++i){
			Boolean_t IsXYZVarNum = FALSE;
			for (int j = 0; j < 3 && !IsXYZVarNum; ++j){
				IsXYZVarNum = (XYZVarNums[j] == i);
			}
			if (!IsXYZVarNum){
				for (int j = 0; j < 2 && IsOk; ++j){
					STGetPtrs[j][VarNum] = TecUtilDataValueGetReadableRef(GetSTZoneNums[j], i);
					IsOk = (VALID_REF(STGetPtrs[j][VarNum]));
				}
				STSetPtrs[VarNum] = TecUtilDataValueGetWritableRef(NewZoneNum, i);
				IsOk = (VALID_REF(STSetPtrs[VarNum]));
				if (IsOk){
					for (int j = 0; j < 2 && IsOk; ++j){
						STGetFuncts[j][VarNum] = TecUtilDataValueRefGetGetFunc(STGetPtrs[j][VarNum]);
						IsOk = (VALID_REF(STGetFuncts[j][VarNum]));
					}
					STSetFuncts[VarNum] = TecUtilDataValueRefGetSetFunc(STSetPtrs[VarNum]);
					IsOk = (VALID_REF(STSetFuncts[VarNum]));
				}
				VarNum++;
			}
		}
		if (!IsOk){
			TecUtilDialogErrMsg("Failed to get streamtrace zone pointers for resampling");
			TecUtilLockFinish(AddOnID);
			return IsOk;
		}
	}

	/*
	*	Get length of streamtrace zones to be resampled and
	*	distance between each point of resampled zones.
	*	Here we also find the REAL max I value for the streamtraces,
	*	since it usually stops moving well before it reaches its
	*	max number of points.
	*/
	double DelLength;
	double STLength[2] = { 0.0 };
	LgIndex_t STIJKMax[2][3];

	for (int PassNum = 0; PassNum < 2; ++PassNum){
		if (IsOk){
			TecUtilZoneGetIJK(GetSTZoneNums[PassNum], &STIJKMax[PassNum][0], &STIJKMax[PassNum][1], &STIJKMax[PassNum][2]);
			IsOk = (STIJKMax[PassNum][0] > 1 && STIJKMax[PassNum][1] == 1 && STIJKMax[PassNum][2] == 1);
		}
		if (IsOk){
			vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;
			for (int i = 0; i < 3; ++i)
				PtIm1[i] = STGetFuncts[PassNum][i](STGetPtrs[PassNum][i], 0);
			Boolean_t ValuesStartedChanging = FALSE;
			for (LgIndex_t i = 1; i < STIJKMax[PassNum][0]; ++i){
				for (int j = 0; j < 3; ++j)
					PtI[j] = STGetFuncts[PassNum][j](STGetPtrs[PassNum][j], i);
				STLength[PassNum] += Distance(PtI, PtIm1);
				if (!ValuesStartedChanging && sum(PtIm1 != PtI) == 3){
					ValuesStartedChanging = TRUE;
				}
				else if (sum(PtIm1 == PtI) == 3 && ValuesStartedChanging){
					STIJKMax[PassNum][0] = i - 1;
					break;
				}
				PtIm1 = PtI;
			}
		}
	}
	/*
	*	Need to make sure that we start at the far end of one of the
	*	streamtraces so that when we concatenate the two, its one
	*	continuous line.
	*/


	Boolean_t StartFrom1[2] = { TRUE };
	int HalfNumSTPoints[2];
	if (IsOk){
		vec3 EndPoints[4];
		for (int i = 0; i < 2; ++i){
			for (int j = 0; j < 3; ++j){
				EndPoints[2 * i][j] = STGetFuncts[i][j](STGetPtrs[i][j], 0);
				EndPoints[1 + 2 * i][j] = STGetFuncts[i][j](STGetPtrs[i][j], STIJKMax[i][0] - 1);
			}
		}
		int GapNum = 1, MinGapNum = 1;
		double MinGap = 1e50;
		for (int i = 0; i < 2; ++i){
			for (int j = 2; j < 4; ++j){
				double TmpGap = DistSqr(EndPoints[i], EndPoints[j]);
				if (TmpGap < MinGap){
					MinGap = TmpGap;
					MinGapNum = GapNum;
				}
				GapNum++;
			}
		}
		double TotalLength = STLength[0] + STLength[1] + sqrt(MinGap);
		HalfNumSTPoints[0] = int(STLength[0] / TotalLength * (double)NumSTPoints);
		HalfNumSTPoints[1] = NumSTPoints - HalfNumSTPoints[0];
		DelLength = TotalLength / double(NumSTPoints);
		switch (MinGapNum){
			case 1:
				/*
				*	point 1 of the two streamtraces are next to eachother, so go in
				*	order last1 first1 first2 last2
				*/
				StartFrom1[0] = FALSE;
				StartFrom1[1] = TRUE;
				break;
			case 2:
				/*
				*	point 1 of the first streamtrace is next to last of second one, so go in
				*	order last1 first1 last2 first 2
				*/
				StartFrom1[0] = FALSE;
				StartFrom1[1] = FALSE;
				break;
			case 3:
				/*
				*	last point of first streamtrace next to first point of second, so go in
				*	order first1 last1 first2 last2
				*/
				StartFrom1[0] = TRUE;
				StartFrom1[1] = TRUE;
				break;
			default:
				/*
				*	last point of first streamtrace next to last point of second, so go in
				*	order first1 last1 last2 first2
				*/
				StartFrom1[0] = TRUE;
				StartFrom1[1] = FALSE;
				break;
		}
	}

	/*
	*	Now step down streamtrace zone again, adding a point to the
	*	resampled zone each time a distance of DelLength is reached,
	*	interpolating values of all variables each time. First and
	*	last points will always be included.
	*/
	if (IsOk){
		vec3 PtI = MinVolPt - 10, PtIm1 = MinVolPt - 10;
		vector<double> ValsI(NumVars, MinVolPt[0] - 10);
		vector<double> ValsIm1(NumVars, MinVolPt[0] - 10);
		double ArcLength = 0.0,
			ArcLengthI = 0.0,
			ArcLengthIm1 = 0.0;
		int NewI = 0;

		for (int PassNum = 0; PassNum < 2; ++PassNum){
			int OldI = 0;
			int OldStep = 1;
			if (!StartFrom1[PassNum]){
				OldI = STIJKMax[PassNum][0] - 1;
				OldStep = -1;
			}

			//	Set first point and get values needed for first step in loop
			for (int i = 0; i < NumVars; ++i){
				ValsI[i] = STGetFuncts[PassNum][i](STGetPtrs[PassNum][i], OldI);
				STSetFuncts[i](STSetPtrs[i], NewI, ValsI[i]);
				if (i < 3){
					PtI[i] = ValsI[i];
				}
			}
			NewI++;
			OldI += OldStep;
			int NewIEnd = HalfNumSTPoints[PassNum];
			if (PassNum > 0){
				NewIEnd += HalfNumSTPoints[0];
				ArcLengthI += Distance(PtI, PtIm1);
			}
			PtIm1 = PtI;
			ValsIm1 = ValsI;
			for (NewI = NewI; NewI < NewIEnd; ++NewI){
				ArcLength += DelLength;

				/*
				*	Move along old path until a segment of DelLength has
				*	been traversed, then interpolate to get new values to
				*	add to new streamtrace zone
				*/
				while (OldI < STIJKMax[PassNum][0] && OldI > 0 && ArcLengthI < ArcLength){
					ArcLengthIm1 = ArcLengthI;
					PtIm1 = PtI;
					ValsIm1 = ValsI;

					for (int i = 0; i < NumVars; ++i){
						ValsI[i] = STGetFuncts[PassNum][i](STGetPtrs[PassNum][i], OldI);
						if (i < 3){
							PtI[i] = ValsI[i];
						}
					}

					ArcLengthI += Distance(PtI, PtIm1);

					OldI += OldStep;
				}

				/*
				*	DelLength has been reached, so time to add a new point
				*	to the new streamtrace zone
				*/
				double Ratio = (ArcLength - ArcLengthIm1) / (ArcLengthI - ArcLengthIm1);
				for (int i = 0; i < NumVars; ++i){
					STSetFuncts[i](STSetPtrs[i], NewI,
						ValsIm1[i] + Ratio * (ValsI[i] - ValsIm1[i])
						);
				}
			}
			//	Set last point

			// 					for (int i = 0; i < NumVars; ++i){
			// 						ValsI[i] = STGetFuncts[i](STGetPtrs[i], OldI);
			// 						STSetFuncts[i](STSetPtrs[i], NumSTPoints,
			// 							ValsI[i]
			// 							);
			// 					}
		}
	}

	TecUtilLockFinish(AddOnID);
	return IsOk;
}