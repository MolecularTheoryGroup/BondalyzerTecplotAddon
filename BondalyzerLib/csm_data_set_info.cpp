
#include <string>
#include <sstream>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <windows.h>
#else
#include <unistd.h>
#endif

// for profiling
#include <chrono>
#include <ctime>
#include <ratio>

#include "TECADDON.h"
#include "Set.h"

#include "CSM_DATA_SET_INFO.h"


//for profiling
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

using std::stringstream;

using tecplot::toolbox::Set;

int StatusNumValues;
int StatusMaxNumValues = 20;
vector<duration<double> > StatusMeanList(StatusMaxNumValues, duration<double>(0));

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

// template <class T>
// int const VectorGetElementNum(vector<T> const & SearchVec, T const & Item){
// 	int ElemNum;
// 
// 	vector<T>::const_iterator itr = std::find(SearchVec.cbegin(), SearchVec.cend(), Item);
// 
// 	if (itr == SearchVec.cend()) return -1;
// 
// 	ElemNum = itr - SearchVec.cbegin();
// 	return ElemNum;
// }

int SearchVectorForString(vector<string> const & Vec, string const & SearchString, bool UseVectorStringLength){
	Boolean_t IsFound = FALSE;

	for (int i = 0; i < Vec.size(); ++i)
		if (UseVectorStringLength && SearchString.compare(0, Vec[i].length(), Vec[i]) == 0)
			return i;
		else if (!UseVectorStringLength && SearchString == Vec[i])
			return i;

	return -1;
}

/*
 *	template split string function, except it doesn't know the type of the items to be split in the string
 *	so it can't know the type of the return vector. I hoped it would get that from the variable
 *	receiving the return of the function call, but I guess not, and even if that worked, you could break it
 *	by returning to an auto type variable.
 */
// template <class T>
// vector<T> SplitString(string const &s, string const & delim, bool RemoveAllBlanks, bool RemoveBlankAtEnd)
// {
// 	stringstream ss1(s), ss2;
// 	string tok;
// 	vector<T> out;
// 	if (RemoveAllBlanks) {
// 		while (std::getline(ss1, tok, delim)) {
// 			if (tok.length() > 0){
// 				ss2.str(tok);
// 				out.push_back();
// 				ss2 >> out.back();
// 			}
// 		}
// 	}
// 	else{
// 		while (std::getline(ss1, tok, delim)) {
// 			ss2.str(tok);
// 			out.push_back();
// 			ss2 >> out.back();
// 		}
// 	}
// 	if (RemoveBlankAtEnd && out.back() == T())
// 		out.pop_back();
// 
// 	return out;
// }

vector<string> SplitString(string const &s, string const & delim, bool RemoveAllBlanks, bool RemoveBlankAtEnd) {
  	vector<string> tokens;
  	auto start = 0U;
  	auto end = s.find(delim);
	if (RemoveAllBlanks) {
		RemoveBlankAtEnd = true;
  		while (end != std::string::npos) {
  			if (end - start > 0) {
  				tokens.push_back(s.substr(start, end - start));
  			}
  			start = end + delim.length();
  			end = s.find(delim, start);
  		}
  	}
  	else {
  		while (end != std::string::npos) {
  			tokens.push_back(s.substr(start, end - start));
  			start = end + delim.length();
  			end = s.find(delim, start);
  		}
  	}
  	if (!RemoveBlankAtEnd || s.length() - start > 0)
  		tokens.push_back(s.substr(start, s.length() - start));
  // 	if (tokens.size() == 0) tokens.push_back(s);
  	return tokens;
}

vector<int> SplitStringInt(string const &s, string const & delim) {
	vector<int> tokens;
	auto start = 0U;
	auto end = s.find(delim);
	while (end != std::string::npos) {
		if (end - start > 0) {
			tokens.push_back(stoi(s.substr(start, end - start)));
		}
		start = end + delim.length();
		end = s.find(delim, start);
	}

	if (s.length() - start > 0)
		tokens.push_back(stoi(s.substr(start, s.length() - start)));
	// 	if (tokens.size() == 0) tokens.push_back(s);
	return tokens;
}

vector<double> SplitStringDbl(string const &s, string const & delim) {
	vector<double> tokens;
	auto start = 0U;
	auto end = s.find(delim);
	while (end != std::string::npos) {
		if (end - start > 0) {
			tokens.push_back(stod(s.substr(start, end - start)));
		}
		start = end + delim.length();
		end = s.find(delim, start);
	}
	if (s.length() - start > 0)
		tokens.push_back(stod(s.substr(start, s.length() - start)));
	// 	if (tokens.size() == 0) tokens.push_back(s);
	return tokens;
}


Boolean_t StringIsInt(string const & s){
	Boolean_t IsInt = TRUE;
	for (int i = 0; i < s.length() && IsInt; ++i)
		IsInt = isdigit(s[i]);

	return IsInt;
}

Boolean_t StringIsFloat(string const & s){
	for (auto const & i : s) if (!isdigit(i) && i != '.') return FALSE;

	return TRUE;
}


string StringReplaceSubString(string const & InStr, string const & OldStr, string const & NewStr){
	string OutStr = InStr;

	size_t Index = 0;

	while (true){
		Index = OutStr.find(OldStr, Index);

		if (Index == string::npos) break;

		OutStr.replace(Index, OldStr.length(), NewStr);

		Index += NewStr.length();
	}

	return OutStr;
}

string StringRemoveSubString(string const & InString,
	string const & SubString)
{
	string NewString = InString;

	size_t Pos = NewString.find_first_of(SubString);
	while (Pos != string::npos) {
		NewString.erase(Pos, SubString.length());
		Pos = NewString.find_first_of(SubString);
	}

	return NewString;
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

string AuxDataMakeStringValidName(string Str){
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

void AuxDataCopy(int SourceNum, int DestNum, bool IsZone){
	AuxData_pa SourcePtr = nullptr, DestPtr = nullptr;

	if (IsZone){
		SourcePtr = TecUtilAuxDataZoneGetRef(SourceNum);
		DestPtr = TecUtilAuxDataZoneGetRef(DestNum);
	}
	else{
		SourcePtr = TecUtilAuxDataVarGetRef(SourceNum);
		DestPtr = TecUtilAuxDataVarGetRef(DestNum);
	}

	if (SourcePtr != nullptr && DestPtr != nullptr){
		int NumItems = TecUtilAuxDataGetNumItems(SourcePtr);
		for (int i = 1; i <= NumItems; ++i){
			char *name, *value;
			Boolean_t retain;
			AuxDataType_e type;
			TecUtilAuxDataGetItemByIndex(SourcePtr, i, &name, reinterpret_cast<ArbParam_t*>(&value), &type, &retain);
			TecUtilAuxDataSetItem(DestPtr, const_cast<const char*>(name), reinterpret_cast<ArbParam_t>(value), type, retain);

			TecUtilStringDealloc(&name);
			TecUtilStringDealloc(&value);
		}
	}
}

Boolean_t AuxDataZoneHasItem(int ZoneNum, string const & CheckString)
{
	Boolean_t ZoneOK = TRUE;

	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	ZoneOK = (TempAuxData != nullptr);
	if (ZoneOK){
		char* TempCStr;
		AuxDataType_e ADTJunk;
		Boolean_t BJunk;
		ZoneOK = TecUtilAuxDataGetItemByName(TempAuxData, CheckString.c_str(), reinterpret_cast<ArbParam_t*>(&TempCStr), &ADTJunk, &BJunk);
		TecUtilStringDealloc(&TempCStr);
	}

	return ZoneOK;
}
Boolean_t AuxDataZoneItemMatches(int ZoneNum, string const & AuxDataName, string const & AuxDataValue)
{
	Boolean_t ZoneOK = TRUE;

	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	ZoneOK = (TempAuxData != nullptr);
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

string AuxDataZoneGetItem(int ZoneNum, string const & AuxDataName)
{
	Boolean_t ZoneOK = TRUE;

	char* TempCStr = nullptr;
	string TmpStr = "-1";

	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	ZoneOK = (TempAuxData != nullptr);
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

Boolean_t AuxDataZoneGetItem(int ZoneNum, string const & AuxDataName, string & Value)
{
	Boolean_t ZoneOK = TRUE;

	char* TempCStr = nullptr;

	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	ZoneOK = (TempAuxData != nullptr);
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

vector<string> AuxDataZoneGetList(int ZoneNum, string const & AuxDataName, string & Value){
	return SplitString(AuxDataZoneGetItem(ZoneNum, AuxDataName), ",,");
}


Boolean_t AuxDataZoneSetItem(int ZoneNum, string const & AuxDataName, string const & AuxDataValue)
{
	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	Boolean_t IsOk = (TempAuxData != nullptr);
	if (IsOk){
		IsOk = TecUtilAuxDataSetStrItem(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), AuxDataValue.c_str(), TRUE);
	}
	return IsOk;
}

Boolean_t AuxDataZoneDeleteItemByName(int ZoneNum, string const & AuxDataName)
{
	AuxData_pa TempAuxData = TecUtilAuxDataZoneGetRef(ZoneNum);
	Boolean_t IsOk = (TempAuxData != nullptr);
	if (IsOk){
		IsOk = TecUtilAuxDataDeleteItemByName(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str());
	}
	return IsOk;
}

Boolean_t AuxDataDataSetDeleteItemByName(string const & AuxDataName)
{
	AuxData_pa TempAuxData = TecUtilAuxDataDataSetGetRef();
	Boolean_t IsOk = (TempAuxData != nullptr);
	if (IsOk){
		IsOk = TecUtilAuxDataDeleteItemByName(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str());
	}
	return IsOk;
}

string AuxDataDataSetGetItem(string const & AuxDataName)
{
	Boolean_t ZoneOK = TRUE;

	char* TempCStr = nullptr;
	string TmpStr;

	AuxData_pa TempAuxData = TecUtilAuxDataDataSetGetRef();
	ZoneOK = (TempAuxData != nullptr);
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

Boolean_t AuxDataDataSetGetItem(string const & AuxDataName, string & Value)
{
	Boolean_t ZoneOK = TRUE;

	char* TempCStr = nullptr;

	AuxData_pa TempAuxData = TecUtilAuxDataDataSetGetRef();
	ZoneOK = (TempAuxData != nullptr);
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


Boolean_t AuxDataDataSetSetItem(string const & AuxDataName, string const & AuxDataValue)
{
	AuxData_pa TempAuxData = TecUtilAuxDataDataSetGetRef();
	Boolean_t IsOk = (TempAuxData != nullptr);
	if (IsOk){
		IsOk = TecUtilAuxDataSetStrItem(TempAuxData, AuxDataMakeStringValidName(AuxDataName).c_str(), AuxDataValue.c_str(), TRUE);
	}
	return IsOk;
}

/*
	*	Begin var/zone wrapper functions
	*/
SmInteger_t VarNumByNameList(vector<string> const & VarNameList, bool const PartialMatch)
{
	SmInteger_t VarNum = -1;

	EntIndex_t NumberOfVars = TecUtilDataSetGetNumVars();

	for (EntIndex_t CurrentVar = NumberOfVars; CurrentVar >= 1 && VarNum < 0; CurrentVar--){
		VarName_t CurrentVarName;
		if (TecUtilVarGetName(CurrentVar, &CurrentVarName)){
			for (auto const & s : VarNameList){
				if (!strncmp(CurrentVarName, s.c_str(), s.size())){
					VarNum = CurrentVar;
				}
			}
			TecUtilStringDealloc(&CurrentVarName);
		}
	}

	return VarNum;
}

SmInteger_t VarNumByName(string const & VarName, bool const PartialMatch)
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

SmInteger_t ZoneNumByNameList(vector<string> const & ZoneNameList)
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

SmInteger_t ZoneNumByName(string const & ZoneName, bool const ActiveOnly, bool const PartialMatch)
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


vec3 GetDelXYZ_Ordered3DZone(vector<int> const & XYZVarNums, int ZoneNum){
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


vector<vec3> ZoneXYZVarGetMinMax_Ordered3DZone(vector<int> const & XYZVarNums, int ZoneNum){

	vector<vec3> MinMaxPts(2);

	ZoneXYZVarGetMinMax_Ordered3DZone(XYZVarNums, ZoneNum, MinMaxPts[0], MinMaxPts[1]);

	return MinMaxPts;
}


void ZoneXYZVarGetBasisVectors_Ordered3DZone(vector<int> const & XYZVarNums, int ZoneNum, mat33 & BasisVectors, vec3 & BVExtent){
	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));
	REQUIRE(XYZVarNums.size() == 3);
	REQUIRE(ZoneNum >= 1 && ZoneNum <= TecUtilDataSetGetNumZones());
	for (int i = 0; i < 3; ++i)
		REQUIRE(XYZVarNums[i] >= 1 && XYZVarNums[i] <= TecUtilDataSetGetNumVars());

	LgIndex_t IJK[3];
	TecUtilZoneGetIJK(ZoneNum, &IJK[0], &IJK[1], &IJK[2]);

	vec3 LowXYZ, HighXYZ;

	FieldVecPointer_c XYZPtr;

	TecUtilDataLoadBegin();

	if (XYZPtr.GetReadPtr(ZoneNum, XYZVarNums)){
		LowXYZ = XYZPtr[0];

		for (int dir = 0; dir < 3; ++dir){
			HighXYZ = XYZPtr[IndexFromIJK(
				dir == 0 ? IJK[0] : 1,
				dir == 1 ? IJK[1] : 1,
				dir == 2 ? IJK[2] : 1,
				IJK[0], IJK[1]) - 1];
			BasisVectors.col(dir) = HighXYZ - LowXYZ;
			BVExtent[dir] = norm(BasisVectors.col(dir));
		}
	}

	TecUtilDataLoadEnd();
}


void ZoneXYZVarGetMinMax_Ordered3DZone(vector<int> const & XYZVarNums, int ZoneNum, vec3 & MinXYZ, vec3 & MaxXYZ){
	REQUIRE(TecUtilZoneIsOrdered(ZoneNum));
	REQUIRE(XYZVarNums.size() == 3);
	REQUIRE(ZoneNum >= 1 && ZoneNum <= TecUtilDataSetGetNumZones());
	for (int i = 0; i < 3; ++i)
		REQUIRE(XYZVarNums[i] >= 1 && XYZVarNums[i] <= TecUtilDataSetGetNumVars());

	LgIndex_t IJK[2][3];
	IJK[0][0] = IJK[0][1] = IJK[0][2] = 1;
	TecUtilZoneGetIJK(ZoneNum, &IJK[1][0], &IJK[1][1], &IJK[1][2]);

	MinXYZ.fill(DBL_MAX);
	MaxXYZ.fill(DBL_MIN);

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

void StatusLaunch(string const & StatusStr, AddOn_pa const & AddOnID, Boolean_t ShowScale, Boolean_t ShowButton){
	TecUtilLockStart(AddOnID);
// 	TecUtilDrawGraphics(TRUE);
	TecUtilStatusSuspend(FALSE);

	StatusNumValues = 0;
	StatusMeanList.resize(StatusMaxNumValues, duration<double>(0));

	TecUtilStatusStartPercentDone(StatusStr.c_str(), ShowButton, ShowScale);
// 	TecUtilDialogLaunchPercentDone(StatusStr.c_str(), ShowScale);

	TecUtilStatusSuspend(TRUE);
// 	TecUtilDrawGraphics(FALSE);
	TecUtilLockFinish(AddOnID);
}

void StatusDrop(AddOn_pa const & AddOnID){
	TecUtilLockStart(AddOnID);
// 	TecUtilDrawGraphics(TRUE);
	TecUtilStatusSuspend(FALSE);

	TecUtilStatusFinishPercentDone();
// 	TecUtilDialogDropPercentDone();

	TecUtilLockFinish(AddOnID);
}

string PrintDuration(duration<double> dur){
	using namespace std::chrono;
	using day_t = duration<long, std::ratio<3600 * 24>>;
	
	auto d = duration_cast<day_t>(dur);
	auto h = duration_cast<hours>(dur -= d);
	auto m = duration_cast<minutes>(dur -= h);
	auto s = duration_cast<seconds>(dur -= m);
	auto ms = duration_cast<seconds>(dur -= s);

	string str = 
		(d.count() > 0 ? to_string(d.count()) + "-" : "")
		+ (h.count() > 0 ? (h.count() > 9 ? "" : "0") + to_string(h.count()) + ":" : "")
		+ (m.count() > 9 ? "" : "0") + to_string(m.count()) + ":"
		+ (s.count() > 9 ? "" : "0") + to_string(s.count());

	return str;
}

Boolean_t StatusUpdate(unsigned int CurrentNum,
	unsigned int TotalNum,
	string const & ProgresssText,
	AddOn_pa const & AddOnID,
	high_resolution_clock::time_point StartTime)
{
	unsigned int Percent = MAX(0, MIN(static_cast<int>((static_cast<double>(CurrentNum) + 0.5) / static_cast<double>(TotalNum)* 100.), 100));

	TecUtilLockStart(AddOnID);


// 	TecUtilDrawGraphics(TRUE);
	TecUtilStatusSuspend(FALSE);
	if (ProgresssText != string("")){
		stringstream ss;
		ss << ProgresssText << "  (" << Percent << "% Complete)";
		if (StartTime != high_resolution_clock::time_point() && Percent > 0 && Percent < 100)
		{
			if (CurrentNum > (int)(0.04 * (double)TotalNum)){
				high_resolution_clock::time_point CurrentTime = high_resolution_clock::now();
				duration<double> Elapsed = duration_cast<duration<double>>(CurrentTime - StartTime);
				// 			StatusMeanList[StatusNumValues % StatusMaxNumValues] = Elapsed / (StatusNumValues + 1.);
				// 			StatusNumValues++;
				// 			duration<double> Mean = duration<double>(0);
				// 			for (int i = 0; i < MIN(StatusNumValues, StatusMaxNumValues); ++i){
				// 				Mean += StatusMeanList[i];
				// 			}
				// 			Mean /= (double)MIN(StatusNumValues, StatusMaxNumValues);
				// 			duration<double> Remaining = Mean * (double)(TotalNum - CurrentNum);
				duration<double> Remaining = (Elapsed / (double)CurrentNum) * (double)(TotalNum - CurrentNum);
				ss << " (Ela. " << PrintDuration(Elapsed) << ". Rem. " << PrintDuration(Remaining) << ". Est. Tot. " << PrintDuration(Elapsed + Remaining) << ")";
			}
			else{
				ss << " Please wait";
			}
		}
		TecUtilStatusSetPercentDoneText(ss.str().c_str());
		// 		TecUtilDialogSetPercentDoneText(ss.str().c_str());
	}
	Boolean_t IsOk = TecUtilStatusCheckPercentDone(Percent);
// 	Boolean_t IsOk = TecUtilDialogCheckPercentDone(Percent);
	TecUtilStatusSuspend(TRUE);
// 	TecUtilDrawGraphics(FALSE);

	TecUtilLockFinish(AddOnID);

	return IsOk;
}


LgIndex_t IndexFromIJK(LgIndex_t I,
	LgIndex_t J,
	LgIndex_t K,
	LgIndex_t MaxI,
	LgIndex_t MaxJ)
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

vector<LgIndex_t> IJKFromIndex(LgIndex_t Index,
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


Boolean_t SaveVec3VecAsScatterZone(vector<vec3> const & VecVec, 
	string const & ZoneName,
	ColorIndex_t const & Color,
	vector<EntIndex_t> const & XYZVarNums){
	Boolean_t IsOk = VecVec.size() > 0;
	if (IsOk){

		TecUtilDataSetAddZone(ZoneName.c_str(), VecVec.size(), 1, 1, ZoneType_Ordered, nullptr);
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

Boolean_t SaveTetVec3VecAsFEZone(vector<vec3> const & Verts,
	string const & ZoneName,
	ColorIndex_t const & Color,
	vector<EntIndex_t> const & XYZVarNums){
	Boolean_t IsOk = Verts.size() == 4;
	if (IsOk){

		TecUtilDataSetAddZone(ZoneName.c_str(), Verts.size(), 1, 1, ZoneType_FETetra, nullptr);
		// 		for (int i = 0; i < VecVec.size(); ++i){
		// 			for (int dir = 0; dir < 3; ++dir){
		// 				TecUtilDataValueSetByZoneVar(TecUtilDataSetGetNumZones(), XYZVarNums[dir], i + 1, VecVec[i][dir]);
		// 			}
		// 		}
		int ZoneNum = TecUtilDataSetGetNumZones();

		vector<vector<double> > TmpValues(3, vector<double>(Verts.size()));
		for (int i = 0; i < Verts.size(); ++i){
			for (int j = 0; j < 3; ++j)
				TmpValues[j][i] = Verts[i][j];
		}

		for (int i = 0; i < 3 && IsOk; ++i){
			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
			IsOk = VALID_REF(SetFDPtr);
			if (IsOk){
				TecUtilDataValueArraySetByRef(SetFDPtr, 1, Verts.size(), reinterpret_cast<void*>(const_cast<double*>(TmpValues[i].data())));
			}
		}

		NodeMap_pa NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
		IsOk = VALID_REF(NodeMap);
		for (int i = 1; i <= 4; ++i){
			TecUtilDataNodeSetByRef(NodeMap, 1, i, i);
		}

		Set_pa TmpSet = TecUtilSetAlloc(FALSE);
		TecUtilSetAddMember(TmpSet, TecUtilDataSetGetNumZones(), FALSE);
		TecUtilZoneSetMesh(SV_SHOW, TmpSet, 0, FALSE);
		TecUtilZoneSetContour(SV_SHOW, TmpSet, 0, FALSE);
		TecUtilZoneSetEdgeLayer(SV_SHOW, TmpSet, 0, TRUE);
		TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0, FALSE);
		TecUtilZoneSetShade(SV_SHOW, TmpSet, 0, TRUE);
		TecUtilZoneSetShade(SV_COLOR, TmpSet, 0, Color);

		TecUtilZoneSetActive(TmpSet, AssignOp_PlusEquals);

		TecUtilSetDealloc(&TmpSet);
	}
	return IsOk;
}


void SetZoneStyle(vector<int> ZoneNums,
	CSMZoneStyle_e ZoneStyle,
	ColorIndex_t Color,
	double const Size)
{
	if (ZoneNums.size() == 0) ZoneNums.push_back(0);
	for (auto ZoneNum : ZoneNums) {
		if (ZoneNum <= 0 || ZoneNum > TecUtilDataSetGetNumZones()) {
			ZoneNum = TecUtilDataSetGetNumZones();
		}
		if (ZoneNum <= 0)
			return;

		int IMaxPathCutoff = 100;

		if (ZoneStyle == ZoneStyle_Invalid) {
			int I, J, K;
			TecUtilZoneGetIJK(ZoneNum, &I, &J, &K);
			if (TecUtilZoneIsOrdered(ZoneNum)) {
				if (J > 1)
					if (K > 1)
						ZoneStyle = ZoneStyle_Volume;
					else
						ZoneStyle = ZoneStyle_Surface;
				else if (I > IMaxPathCutoff)
					ZoneStyle = ZoneStyle_Path;
				else
					ZoneStyle = ZoneStyle_Points;
			}
			else if (TecUtilZoneIsFiniteElement(ZoneNum)) {
				if (K == 1)
					ZoneStyle = ZoneStyle_Points;
				else if (K == 2)
					ZoneStyle = ZoneStyle_Path;
				else if (K >= 3) {
					ZoneType_e ZoneType = TecUtilZoneGetType(ZoneNum);
					if (ZoneType == ZoneType_FETriangle
						|| ZoneType == ZoneType_FEQuad
						|| ZoneType == ZoneType_FEPolygon)
					{
						ZoneStyle = ZoneStyle_Surface;
					}
					else if (ZoneType == ZoneType_FETetra
						|| ZoneType == ZoneType_FEBrick
						|| ZoneType == ZoneType_FEPolyhedron)
					{
						ZoneStyle = ZoneStyle_Volume;
					}
				}
			}
		}

		ColorIndex_t ColorDefault = Black_C;

		Set ZoneSet(ZoneNum);

		if (ZoneStyle <= ZoneStyle_Points) {
			double SizeDefault = 0.5;
			TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, ZoneSet.getRef(), GeomShape_Sphere);
			TecUtilZoneSetScatter(SV_FRAMESIZE, ZoneSet.getRef(), (Size > 0 ? Size : SizeDefault), 0);
			TecUtilZoneSetScatter(SV_COLOR, ZoneSet.getRef(), 0, (Color != InvalidColor_C ? Color : ColorDefault));
		}
		else if (ZoneStyle <= ZoneStyle_Path) {
			double SizeDefault = 0.2;
			TecUtilZoneSetScatter(SV_SHOW, ZoneSet.getRef(), 0, FALSE);
			TecUtilZoneSetMesh(SV_SHOW, ZoneSet.getRef(), 0, TRUE);
			TecUtilZoneSetMesh(SV_COLOR, ZoneSet.getRef(), 0, (Color != InvalidColor_C ? Color : ColorDefault));
			TecUtilZoneSetMesh(SV_LINETHICKNESS, ZoneSet.getRef(), (Size > 0 ? Size : SizeDefault), 0);
		}
		else if (ZoneStyle <= ZoneStyle_Surface) {
			ColorDefault = Cyan_C;
			TecUtilZoneSetScatter(SV_SHOW, ZoneSet.getRef(), 0, FALSE);
			TecUtilZoneSetShade(SV_SHOW, ZoneSet.getRef(), 0, TRUE);
			TecUtilZoneSetShade(SV_COLOR, ZoneSet.getRef(), 0, (Color != InvalidColor_C ? Color : ColorDefault));
		}
		else if (ZoneStyle <= ZoneStyle_Volume) {
			double SizeDefault = 0.1;
			TecUtilZoneSetScatter(SV_SHOW, ZoneSet.getRef(), 0, FALSE);
			TecUtilZoneSetEdgeLayer(SV_SHOW, ZoneSet.getRef(), 0, TRUE);
			TecUtilZoneSetEdgeLayer(SV_LINETHICKNESS, ZoneSet.getRef(), (Size > 0 ? Size : SizeDefault), 0);
			TecUtilZoneSetEdgeLayer(SV_COLOR, ZoneSet.getRef(), 0, (Color != InvalidColor_C ? Color : ColorDefault));
		}
	}
}