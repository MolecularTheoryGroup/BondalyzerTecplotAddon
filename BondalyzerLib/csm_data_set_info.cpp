
#include <string>
#include <sstream>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "TECADDON.h"
#include "CSM_DATA_SET_INFO.h"


using std::stringstream;


size_t getTotalSystemMemory()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx(&status);
	return status.ullTotalPhys;
#else
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
#endif
}


const int SearchVectorForString(const vector<string> & Vec, const string & SearchString, const bool & UseVectorStringLength){
	Boolean_t IsFound = FALSE;

	for (int i = 0; i < Vec.size(); ++i)
		if (UseVectorStringLength && SearchString.compare(0, Vec[i].length(), Vec[i]) == 0)
			return i;
		else if (!UseVectorStringLength && SearchString == Vec[i])
			return i;

	return -1;
}

const vector<string> SplitString(const string &s, char delim) {
	stringstream ss(s);
	string item;
	vector<string> tokens;
	while (getline(ss, item, delim)) {
		tokens.push_back(item);
	}
	return tokens;
}

const Boolean_t StringIsInt(const string & s){
	Boolean_t IsInt = TRUE;
	for (int i = 0; i < s.length() && IsInt; ++i)
		IsInt = isdigit(s[i]);

	return IsInt;
}

const Boolean_t StringIsFloat(const string & s){
	for (const auto & i : s) if (!isdigit(i) && i != '.') return FALSE;

	return TRUE;
}



/*
	*	To guarantee valid aux data names (the name, not the contents).
	*	A valid name may only contain '_', '.', letters and numbers, and must
	*	start with a letter or '_'.
	*	This will replace ANY invalid characters with '_', and if the first character
	*	is not a letter or '_'. a '_' will be prepended.
	*	
	*	If a blank string is given, it will ouput 'Blank_Name'.
	*	
	*	ASCII numbers for checking validity of characters:
	*		sym		num
	*		
	*		.		46
	*		_		95
	*		0-9		48-57
	*		A-Z		65-90
	*		a-z		97-122
	*/

const string AuxDataMakeStringValidName(string Str){
	if (Str.length() <= 0) return "Blank_Name";

	if (!(Str[0] == 95 || (Str[0] >= 65 && Str[0] <= 90) || (Str[0] >= 97 && Str[0] <= 122)))
		Str = "_" + Str;

	for (char & Chr : Str){
		if (!(Chr == 95
			|| Chr == 46
			|| (Chr >= 48 && Chr <= 57)
			|| (Chr >= 65 && Chr <= 90)
			|| (Chr >= 97 && Chr <= 122)))
		{
			Chr = '_';
		}
	}

	return Str;
}

/*
	*	Begin aux data wrapper functions
	*/

const Boolean_t AuxDataZoneHasItem(const int & ZoneNum, const string & CheckString)
{
	Boolean_t ZoneOK = TRUE;

	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	ZoneOK = VALID_REF(TempAuxData);
	if (ZoneOK){
		char* TempCStr;
		AuxDataType_e ADTJunk;
		Boolean_t BJunk;
		ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, CheckString.c_str(), reinterpret_cast<ArbParam_t*>(&TempCStr), &ADTJunk, &BJunk);
		TecUtilStringDealloc(&TempCStr);
	}

	return ZoneOK;
}
const Boolean_t AuxDataZoneItemMatches(const int & ZoneNum, const string & AuxDataName, const string & AuxDataValue)
{
	Boolean_t ZoneOK = TRUE;

	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	ZoneOK = VALID_REF(TempAuxData);
	if (ZoneOK){
		char* TempCStr;
		AuxDataType_e ADTJunk;
		Boolean_t BJunk;
		ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataName.c_str(), reinterpret_cast<ArbParam_t*>(&TempCStr), &ADTJunk, &BJunk);
		if (ZoneOK){
			ZoneOK = (AuxDataValue.compare(TempCStr) == 0);
		}
		TecUtilStringDealloc(&TempCStr);
	}

	return ZoneOK;
}

const string AuxDataZoneGetItem(const int & ZoneNum, const string & AuxDataName)
{
	Boolean_t ZoneOK = TRUE;

	char* TempCStr = NULL;
	string TmpStr = "";

	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	ZoneOK = VALID_REF(TempAuxData);
	if (ZoneOK){
		AuxDataType_e ADTJunk;
		Boolean_t BJunk;
		ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), reinterpret_cast<ArbParam_t*>(&TempCStr), &ADTJunk, &BJunk);
		if (ZoneOK){
			TmpStr = TempCStr;
		}
		TecUtilStringDealloc(&TempCStr);
	}

	return TmpStr;
}

const Boolean_t AuxDataZoneGetItem(const int & ZoneNum, const string & AuxDataName, string & Value)
{
	Boolean_t ZoneOK = TRUE;

	char* TempCStr = NULL;

	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	ZoneOK = VALID_REF(TempAuxData);
	if (ZoneOK){
		AuxDataType_e ADTJunk;
		Boolean_t BJunk;
		ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), reinterpret_cast<ArbParam_t*>(&TempCStr), &ADTJunk, &BJunk);
		if (ZoneOK){
			Value = TempCStr;
		}
		TecUtilStringDealloc(&TempCStr);
	}

	return ZoneOK;
}


const Boolean_t AuxDataZoneSetItem(const int & ZoneNum, const string & AuxDataName, const string & AuxDataValue)
{
	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	Boolean_t IsOk = VALID_REF(TempAuxData);
	if (IsOk){
		IsOk = TecUtilAuxDataSetStrItem(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), AuxDataValue.c_str(), TRUE);
	}
	return IsOk;
}

const Boolean_t AuxDataZoneDeleteItemByName(const int & ZoneNum, const string & AuxDataName)
{
	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	Boolean_t IsOk = VALID_REF(TempAuxData);
	if (IsOk){
		IsOk = TecUtilAuxDataDeleteItemByName(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str());
	}
	return IsOk;
}

const Boolean_t AuxDataDataSetDeleteItemByName(const string & AuxDataName)
{
	AuxData_pa TempAuxData = TecUtilAuxDataDataSetGetRef();
	Boolean_t IsOk = VALID_REF(TempAuxData);
	if (IsOk){
		IsOk = TecUtilAuxDataDeleteItemByName(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str());
	}
	return IsOk;
}

const string AuxDataDataSetGetItem(const string & AuxDataName)
{
	Boolean_t ZoneOK = TRUE;

	char* TempCStr = NULL;
	string TmpStr;

	AuxData_pa TempAuxData = TecUtilAuxDataDataSetGetRef();
	ZoneOK = VALID_REF(TempAuxData);
	if (ZoneOK){
		AuxDataType_e ADTJunk;
		Boolean_t BJunk;
		ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), reinterpret_cast<ArbParam_t*>(&TempCStr), &ADTJunk, &BJunk);
		if (ZoneOK){
			TmpStr = TempCStr;
		}
		TecUtilStringDealloc(&TempCStr);
	}

	return TmpStr;
}

const Boolean_t AuxDataDataSetGetItem(const string & AuxDataName, string & Value)
{
	Boolean_t ZoneOK = TRUE;

	char* TempCStr = NULL;

	AuxData_pa TempAuxData = TecUtilAuxDataDataSetGetRef();
	ZoneOK = VALID_REF(TempAuxData);
	if (ZoneOK){
		AuxDataType_e ADTJunk;
		Boolean_t BJunk;
		ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), reinterpret_cast<ArbParam_t*>(&TempCStr), &ADTJunk, &BJunk);
		if (ZoneOK){
			Value = TempCStr;
		}
		TecUtilStringDealloc(&TempCStr);
	}

	return ZoneOK;
}


const Boolean_t AuxDataDataSetSetItem(const string & AuxDataName, const string & AuxDataValue)
{
	AuxData_pa TempAuxData = TecUtilAuxDataDataSetGetRef();
	Boolean_t IsOk = VALID_REF(TempAuxData);
	if (IsOk){
		IsOk = TecUtilAuxDataSetStrItem(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), AuxDataValue.c_str(), TRUE);
	}
	return IsOk;
}

/*
	*	Begin var/zone wrapper functions
	*/
SmInteger_t VarNumByNameList(const vector<string> & VarNameList, const bool PartialMatch)
{
	SmInteger_t VarNum = -1;

	EntIndex_t NumberOfVars = TecUtilDataSetGetNumVars();

	for (EntIndex_t CurrentVar = NumberOfVars; CurrentVar >= 1 && VarNum < 0; CurrentVar--){
		VarName_t CurrentVarName;
		if (TecUtilVarGetName(CurrentVar, &CurrentVarName)){
			for (const auto & s : VarNameList){
				if (!strncmp(CurrentVarName, s.c_str(), s.size())){
					VarNum = CurrentVar;
				}
			}
			TecUtilStringDealloc(&CurrentVarName);
		}
	}

	return VarNum;
}

SmInteger_t VarNumByName(const string & VarName, const bool PartialMatch)
{
	SmInteger_t VarNum = -1;

	EntIndex_t NumberOfVars = TecUtilDataSetGetNumVars();

	// 	for (EntIndex_t CurrentVar = NumberOfVars; CurrentVar >= 1 && VarNum < 0; CurrentVar--){
	for (EntIndex_t CurrentVar = 1; CurrentVar <= NumberOfVars && VarNum < 0; CurrentVar++){
		VarName_t CurrentVarName;
		if (TecUtilVarGetName(CurrentVar, &CurrentVarName)){
			if (!strncmp(CurrentVarName, VarName.c_str(), VarName.size())){
				VarNum = CurrentVar;
			}
			TecUtilStringDealloc(&CurrentVarName);
		}
	}

	return VarNum;
}

SmInteger_t ZoneNumByNameList(const vector<string> & ZoneNameList)
{
	SmInteger_t ZoneNum = -1;

	EntIndex_t NumberOfZones = TecUtilDataSetGetNumZones();

	for (EntIndex_t CurrentZone = 1; CurrentZone <= NumberOfZones && ZoneNum < 0; CurrentZone++){
		VarName_t CurrentZoneName;
		if (TecUtilZoneGetName(CurrentZone, &CurrentZoneName)){
			for (int CurrentNameCheck = 0; CurrentNameCheck < ZoneNameList.size() && ZoneNum < 0; CurrentNameCheck++){
				if (!strncmp(CurrentZoneName, ZoneNameList[CurrentNameCheck].c_str(), ZoneNameList[CurrentNameCheck].size())){
					ZoneNum = CurrentZone;
				}
			}
			TecUtilStringDealloc(&CurrentZoneName);
		}
	}

	return ZoneNum;
}

SmInteger_t ZoneNumByName(const string & ZoneName, const bool ActiveOnly, const bool PartialMatch)
{
	SmInteger_t ZoneNum = -1;

	EntIndex_t NumberOfZones = TecUtilDataSetGetNumZones();

	for (EntIndex_t CurrentZone = 1; CurrentZone <= NumberOfZones && ZoneNum < 0; CurrentZone++){
		ZoneName_t CurrentZoneName;
		if (TecUtilZoneIsActive(CurrentZone) || !ActiveOnly){
			if (TecUtilZoneGetName(CurrentZone, &CurrentZoneName)){
				if (!strncmp(CurrentZoneName, ZoneName.c_str(), ZoneName.size())){
					ZoneNum = CurrentZone;
				}
				TecUtilStringDealloc(&CurrentZoneName);
			}
		}
	}

	return ZoneNum;
}


const vec3 GetDelXYZ_Ordered3DZone(const vector<int> & XYZVarNums, const int & ZoneNum){
	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));
	REQUIRE(XYZVarNums.size() == 3);
	REQUIRE(ZoneNum >= 1 && ZoneNum <= TecUtilDataSetGetNumZones());
	for (int i = 0; i < 3; ++i)
		REQUIRE(XYZVarNums[i] >= 1 && XYZVarNums[i] <= TecUtilDataSetGetNumVars());

	vec3 DelXYZ;

	LgIndex_t IJK[3];
	TecUtilZoneGetIJK(ZoneNum, &IJK[0], &IJK[1], &IJK[2]);

	vec3 Pt1, Pt2;

	for (int i = 0; i < 3; ++i){
		Pt1[i] = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, 1, 1, IJK[0], IJK[1]));
		Pt2[i] = TecUtilDataValueGetByZoneVar(ZoneNum,
											XYZVarNums[i],
											IndexFromIJK(MIN(2, IJK[0]),
														MIN(2, IJK[1]),
														MIN(2, IJK[2]),
														IJK[0], IJK[1]));
	}

	for (int i = 0; i < 3; ++i)
		DelXYZ[i] = ABS(Pt2[i] - Pt1[i]);

// 		for (int i = 0; i < 3; ++i){
// 			if (IJK[i] > 1){
// 				double pt1, pt2;
// 				pt1 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1,1,1,IJK[0],IJK[1]));
// 				switch (i){
// 					case 0:
// 						pt2 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(2,1,1,IJK[0],IJK[1]));
// 						break;
// 					case 1:
// 						pt2 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, 2, 1, IJK[0], IJK[1]));
// 						break;
// 					default:
// 						pt2 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, 1, 2, IJK[0], IJK[1]));
// 						break;
// 				}
// 				DelXYZ[i] = ABS(pt2-pt1);
// 			}
// 			else
// 				DelXYZ[i] = 0.0;
// 		}

	return DelXYZ;
}


const vector<vec3> ZoneXYZVarGetMinMax_Ordered3DZone(const vector<int> & XYZVarNums, const int & ZoneNum){

	vector<vec3> MinMaxPts(2);

	ZoneXYZVarGetMinMax_Ordered3DZone(XYZVarNums, ZoneNum, MinMaxPts[0], MinMaxPts[1]);

	return MinMaxPts;
}


void ZoneXYZVarGetMinMax_Ordered3DZone(const vector<int> & XYZVarNums, const int & ZoneNum, vec3 & MinXYZ, vec3 & MaxXYZ){
	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));
	REQUIRE(XYZVarNums.size() == 3);
	REQUIRE(ZoneNum >= 1 && ZoneNum <= TecUtilDataSetGetNumZones());
	for (int i = 0; i < 3; ++i)
		REQUIRE(XYZVarNums[i] >= 1 && XYZVarNums[i] <= TecUtilDataSetGetNumVars());

	LgIndex_t IJK[2][3];
	IJK[0][0] = IJK[0][1] = IJK[0][2] = 1;
	TecUtilZoneGetIJK(ZoneNum, &IJK[1][0], &IJK[1][1], &IJK[1][2]);

	MinXYZ.fill(1e150);
	MaxXYZ.fill(-1e150);

	double  Value;

	for (int i = 0; i < 2; ++i){
		for (int j = 0; j < 2; ++j){
			for (int k = 0; k < 2; ++k){
				for (int dir = 0; dir < 3; ++dir){
					Value = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[dir], IndexFromIJK(IJK[i][0], IJK[j][1], IJK[k][2], IJK[1][0], IJK[1][1]));
					MinXYZ[dir] = MIN(MinXYZ[dir], Value);
					MaxXYZ[dir] = MAX(MaxXYZ[dir], Value);
				}
			}
		}
	}

	// 		for (int i = 0; i < 3; ++i){
	// 			MinMaxPts[0][i] = MIN(MinMaxPts[0][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, 1, 1, IJK[0], IJK[1])));
	// 			MinMaxPts[0][i] = MIN(MinMaxPts[0][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0], 1, 1, IJK[0], IJK[1])));
	// 			MinMaxPts[0][i] = MIN(MinMaxPts[0][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, IJK[1], 1, IJK[0], IJK[1])));
	// 			MinMaxPts[0][i] = MIN(MinMaxPts[0][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, 1, IJK[2], IJK[0], IJK[1])));
	// 			MinMaxPts[0][i] = MIN(MinMaxPts[0][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0], IJK[1], 1, IJK[0], IJK[1])));
	// 			MinMaxPts[0][i] = MIN(MinMaxPts[0][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0], 1, IJK[2], IJK[0], IJK[1])));
	// 			MinMaxPts[0][i] = MIN(MinMaxPts[0][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, IJK[1], IJK[2], IJK[0], IJK[1])));
	// 			MinMaxPts[0][i] = MIN(MinMaxPts[0][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0], IJK[1], IJK[2], IJK[0], IJK[1])));
	// 
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1,			1,		1,		IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0],	1,		1,		IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1,			IJK[1], 1,		IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1,			1,		IJK[2], IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0],	IJK[1], 1,		IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0],	1,		IJK[2],	IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1,			IJK[1], IJK[2],	IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0],	IJK[1], IJK[2], IJK[0], IJK[1])));
	// 
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1,			1,		1,		IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1,			1,		IJK[2], IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1,			IJK[1], 1,		IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1,			IJK[1], IJK[2],	IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0],	1,		1,		IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0],	1,		IJK[2],	IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0],	IJK[1], 1,		IJK[0], IJK[1])));
	// 			MinMaxPts[1][i] = MAX(MinMaxPts[1][i], TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[0],	IJK[1], IJK[2], IJK[0], IJK[1])));
	// 		}
	// 
	// 		for (int j = 0; j < 2; ++j){
	// 			for (int i = 0; i < 3; ++i){
	// 				if (IJK[i] > 1){
	// 					double pt1, pt2, pt3;
	// 					pt1 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, 1, 1, IJK[0], IJK[1]));
	// 					switch (i){
	// 						case 0:
	// 							pt2 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(2, 1, 1, IJK[0], IJK[1]));
	// 							pt3 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(IJK[i], 1, 1, IJK[0], IJK[1]));
	// 							break;
	// 						case 1:
	// 							pt2 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, 2, 1, IJK[0], IJK[1]));
	// 							pt3 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, IJK[i], 1, IJK[0], IJK[1]));
	// 							break;
	// 						default:
	// 							pt2 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, 1, 2, IJK[0], IJK[1]));
	// 							pt3 = TecUtilDataValueGetByZoneVar(ZoneNum, XYZVarNums[i], IndexFromIJK(1, 1, IJK[i], IJK[0], IJK[1]));
	// 							break;
	// 					}
	// 					DelXYZ[i] = ABS(pt2 - pt1);
	// 				}
	// 				else
	// 					DelXYZ[i] = 0.0;
	// 			}
	// 		}
}

/*
	*	Begin state change wrapper functions
	*/




/*
	*	Begin other functions
	*/


const Boolean_t SetPercent(unsigned int CurrentNum, unsigned int TotalNum, const string & ProgresssText, const AddOn_pa & AddOnID){
	unsigned int Percent = MAX(0, MIN(static_cast<int>(static_cast<double>(CurrentNum) / static_cast<double>(TotalNum)* 100.), 100));

	stringstream ss;
	ss << ProgresssText << "  (" << Percent << "% Complete)";

	TecUtilLockStart(AddOnID);
	Boolean_t IsOk = TecUtilDialogCheckPercentDone(Percent);
	TecUtilDialogSetPercentDoneText(ss.str().c_str());

	TecUtilLockFinish(AddOnID);

	return IsOk;
}


const LgIndex_t IndexFromIJK(const LgIndex_t &I,
	const LgIndex_t &J,
	const LgIndex_t &K,
	const LgIndex_t &MaxI,
	const LgIndex_t &MaxJ)
{
	return I + (J - 1) * MaxI + (K - 1) * MaxI * MaxJ;
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
const LgIndex_t IndexFromIJK(const LgIndex_t &    I,
	const LgIndex_t &    J,
	const LgIndex_t &    K,
	const LgIndex_t & IMax,
	const LgIndex_t & JMax,
	const LgIndex_t & KMax,
	const Boolean_t & PeriodicBC)
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

const vector<LgIndex_t> IJKFromIndex(const LgIndex_t & Index,
									const vector<LgIndex_t>& IJKMax){
	vector<LgIndex_t> IJK(3);
	Boolean_t IsOk = (Index <= IJKMax[0] * IJKMax[1] * IJKMax[2] + 1);
	if (IsOk){
		LgIndex_t IJ = IJKMax[0] * IJKMax[1];
		IJK[2] = Index / IJ + 1;
		IJK[1] = (Index % IJ) / IJKMax[0] + 1;
		IJK[0] = (Index % (IJKMax[0] - 1)) + 1;
	}

	return IJK;
}


const Boolean_t SaveVec3VecAsScatterZone(const vector<vec3> & VecVec, 
	const string & ZoneName,
	const ColorIndex_t & Color,
	const vector<EntIndex_t> & XYZVarNums){
	Boolean_t IsOk = VecVec.size() > 0;
	if (IsOk){

		TecUtilDataSetAddZone(ZoneName.c_str(), VecVec.size(), 1, 1, ZoneType_Ordered, NULL);
// 		for (int i = 0; i < VecVec.size(); ++i){
// 			for (int dir = 0; dir < 3; ++dir){
// 				TecUtilDataValueSetByZoneVar(TecUtilDataSetGetNumZones(), XYZVarNums[dir], i + 1, VecVec[i][dir]);
// 			}
// 		}
		int ZoneNum = TecUtilDataSetGetNumZones();

		vector<vector<double> > TmpValues(3, vector<double>(VecVec.size()));
		for (int i = 0; i < VecVec.size(); ++i){
			for (int j = 0; j < 3; ++j)
				TmpValues[j][i] = VecVec[i][j];
		}

		for (int i = 0; i < 3 && IsOk; ++i){
			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
			IsOk = VALID_REF(SetFDPtr);
			if (IsOk){
				TecUtilDataValueArraySetByRef(SetFDPtr, 1, VecVec.size(), reinterpret_cast<void*>(const_cast<double*>(TmpValues[i].data())));
			}
		}
		Set_pa TmpSet = TecUtilSetAlloc(FALSE);
		TecUtilSetAddMember(TmpSet, TecUtilDataSetGetNumZones(), FALSE);
		TecUtilZoneSetMesh(SV_SHOW, TmpSet, 0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, TmpSet, 0, FALSE);
		TecUtilZoneSetEdgeLayer(SV_SHOW, TmpSet, 0, FALSE);
		TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0, TRUE);
		TecUtilZoneSetScatter(SV_COLOR, TmpSet, 0, Color);
		TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, TmpSet, GeomShape_Sphere);
		TecUtilZoneSetScatter(SV_FRAMESIZE, TmpSet, 1, 0);

		TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);

		TecUtilSetDealloc(&TmpSet);
	}
	return IsOk;
}