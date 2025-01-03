
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>

#include "TECADDON.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_DATA_LOADER_FUNCTIONS.h"

#include "Set.h"
#include "StyleValue.h"

#include <armadillo>
using namespace arma;

using tecplot::toolbox::Set;


using std::vector;
using std::string;

std::string VarType_str[] = {
	"SCF",
	"FRAG",
	"ORTHO",
	"TRANS",
	"GEOMETRY",
	"NOTYPE"
};

std::string ab_str[] = {
	"A",
	"B",
	"NOAB"
};

string GetAtomStr(int AtomNum){
	if (AtomNum <= 0 || AtomNum > ElementSymbolList.size())
		return "";
	return ElementSymbolList[AtomNum - 1];
}

/*
* tw_datatype methods
*/
T41Var_s::T41Var_s(){
	VarType = NOTYPE;
	DTAB = NOAB;
}
T41Var_s::T41Var_s(const char* iName){
	NameStr = iName;
	VarType = NOTYPE;
	DTAB = NOAB;
}
T41Var_s::T41Var_s(const char* iT41Name, const char* iName, const char* iBaseName, Boolean_t UseForIsoSurfaceDisplay){
	this->Set(iT41Name, iName, iBaseName, UseForIsoSurfaceDisplay);
}
T41Var_s::~T41Var_s(){}

void T41Var_s::Set(const char* iT41Name, const char* iName, const char* iBaseName, Boolean_t UseForIsoSurfaceDisplay)
{
	T41NameStr = iT41Name;
	T41FullNameStr = iT41Name;
	NameStr = iName;
	BaseNameStr = iBaseName;
	VarType = NOTYPE;
	DTAB = NOAB;
	UseForIsosurface = UseForIsoSurfaceDisplay;
}

void T41Var_s::SetAB(const std::string &TempVarName){
	for (int i = 0; i < 2; ++i){
		std::string TempStr = T41NameStr;
		TempStr += "_";
		TempStr += ab_str[i];
		if (!TempVarName.compare(0, TempStr.length(), TempStr.c_str())){
			DTAB = ab_e(i);
			T41FullNameStr += "_";
			T41FullNameStr += ab_str[i];
			NameStr += " ";
			NameStr += ab_str[i];
			break;
		}
	}
}

void T41Var_s::SetAB(ab_e In){
	for (int i = 0; i < 2; ++i){
		if (In == ab_e(i)){
			DTAB = ab_e(i);
			T41FullNameStr += "_";
			T41FullNameStr += ab_str[i];
			NameStr += " ";
			NameStr += ab_str[i];
			break;
		}
	}
}

void T41Var_s::SetType(const std::string &TypeStr){
	std::vector<std::string> ChkStrList = { "SCF", "SumFrag", "Ortho", "Trans", "Geometry" };

	NameStr += " ";
	T41FullNameStr += "_";
	for (int i = 0; i < ChkStrList.size(); ++i)
		if (!TypeStr.compare(0, ChkStrList[i].length(), ChkStrList[i])){
		VarType = VarType_e(i);
		NameStr += VarType_str[i];
		T41FullNameStr += VarType_str[i];
		break;
		}
}

void T41Var_s::SetType(VarType_e In){
	NameStr += " ";
	T41FullNameStr += "_";
	for (int i = 0; i < 5; ++i)
		if (In == VarType_e(i)){
		VarType = VarType_e(i);
		NameStr += VarType_str[i];
		T41FullNameStr += VarType_str[i];
		break;
		}
}

Boolean_t T41Var_s::operator== (T41Var_s const & T2) const{
	return (T41FullNameStr == T2.T41FullNameStr);
}

Boolean_t T41Var_s::operator!= (T41Var_s const & T2) const{
	return !(*this == T2);
}

AtomGroup_s::AtomGroup_s(const std::string & AtomTypeName){
	Name = AtomTypeName;
	Count = 1;
}

void AtomGroup_s::AddPosition(vec3 & Position){
	for (int i = 0; i < 3; ++i)
		Positions[i].push_back(Position[i]);
}

/*
* Other functions
*/

Boolean_t IsInDataTypeList(const std::vector<T41Var_s> &SearchVec, T41Var_s const &ToFind){
	Boolean_t IsFound = FALSE;
	for (std::vector<T41Var_s>::const_iterator it = SearchVec.cbegin(); it != SearchVec.cend() && !IsFound; it++)
		if (*it == ToFind)
			IsFound = TRUE;
	return IsFound;
}

int IsVarLine(const std::string &CheckLine){
	std::vector<std::string> StrList = { "x values", "y values", "z values",
		"SCF", "SumFrag", "Ortho", "Trans", "Geometry" };
	for (int i = 0; i < StrList.size(); ++i)
		if (!CheckLine.compare(0, StrList[i].length(), StrList[i]))
			return i + 1;
	return -1;
}


Boolean_t VectorContainsString(const std::vector<std::string> &SearchVec, const std::string &ToFind){
	Boolean_t IsFound = FALSE;

	for (std::vector<std::string>::const_iterator it = SearchVec.cbegin(); it != SearchVec.cend() && !IsFound; it++)
		if (*it == ToFind) 
			return TRUE;

	return IsFound;
}

ColorIndex_t GetAtomColor(const std::vector<AtomColor_s> & AtomColorList, const std::string & AtomName){
	for (int i = 1; i < AtomColorList.size(); ++i)
		if (VectorContainsString(AtomColorList[i].Names, AtomName))
			return AtomColorList[i].Color;

	return AtomColorList[0].Color;
}

void PopulateAtomColorList(std::vector<AtomColor_s> &AtomColorList){
	/*
	*	The majority of species are orange, so the are ommited
	*	and when a species isn't found the default color will
	*	be orange.
	*/
	std::vector<std::string> TempStrList;
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom3_C));

	//	Light Gray
	TempStrList = { "H", "Na", "Pt", "Hg", "Nuclear CPs" };
	AtomColorList.push_back(AtomColor_s(TempStrList, White_C));

	//	Dark Gray
	TempStrList = { "Li", "Be", "Mg", "Al", "Ga", "Ti", "K", "Ag" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom1_C));

	//	Blue
	TempStrList = { "N" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Blue_C));

	//	Red
	TempStrList = { "O", "Bond CPs" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Red_C));

	//	Green (but darker than Green_C)
	TempStrList = { "F" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom12_C));

	//	Tangerine
	TempStrList = { "Si", "P" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom11_C));

	//	Yellow
	TempStrList = { "S" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Yellow_C));

	//	Pastel Green
	TempStrList = { "Cl" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom13_C));

	//	Black (for lack of a darker dark gray)
	TempStrList = { "C", "Ni", "Sn" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Black_C));

	//	Pastel yellow
	TempStrList = { "Cr" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom18_C));

	//	Teal
	TempStrList = { "Mn" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom45_C));

	//	Navy blue
	TempStrList = { "Co" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom38_C));

	//	Dark orange (brown really)
	TempStrList = { "Cu", "Ne" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom32_C));

	//	Pastel pink
	TempStrList = { "Zn" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom50_C));

	//	Forest green
	TempStrList = { "As" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom25_C));

	//	Crimson
	TempStrList = { "Se", "Br" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom9_C));

	//	Pink-red
	TempStrList = { "Ru" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom31_C));

	//	Gold
	TempStrList = { "Au", "Cd" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom19_C));

	//	Pastel purple
	TempStrList = { "I" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom15_C));

	//	Cyan
	TempStrList = { "Os", "Cage CPs" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Cyan_C));

	//	Green
	TempStrList = { "Sg", "Hs", "Ring CPs" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Green_C));

	//	Green (between forest and darker-than-regular)
	TempStrList = { "Bh" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom28_C));

	//	bubble-gum pink
	TempStrList = { "Nd" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom23_C));

	//	lipstick pink
	TempStrList = { "U" };
	AtomColorList.push_back(AtomColor_s(TempStrList, Custom39_C));
}

int PopulateTypeList(std::vector<T41Var_s> &TypeList){

	TypeList.push_back(T41Var_s("x values", "X", "", FALSE));

	TypeList.push_back(T41Var_s("y values", "Y", "", FALSE));

	TypeList.push_back(T41Var_s("z values", "Z", "", FALSE));

	TypeList.push_back(T41Var_s("Coulpot", "Coulomb Potential", "Coulomb Potential", FALSE));

	TypeList.push_back(T41Var_s("XCPot_fit", "Exchange-Correlation Potential", "Exchange-Correlation Potential", FALSE));

	TypeList.push_back(T41Var_s("DensityLap", "Laplacian", "Laplacian", FALSE));

	TypeList.push_back(T41Var_s("Kinetic Energy Density", "Kinetic Energy Density", "Kinetic Energy Density", FALSE));

	TypeList.push_back(T41Var_s("ELF", "Electron Localization Function", "Electron Localization Function", FALSE));

	/*
	*	The next three are only for getting information about the atomic count/types/positions/charges
	*/
	TypeList.push_back(T41Var_s("nnuc", "Atom Count", "Atom Count", FALSE));

	TypeList.push_back(T41Var_s("labels", "Atom Type Names", "Atom Type Names", FALSE));

	TypeList.push_back(T41Var_s("xyznuc", "Atom Positions", "Atom Positions", FALSE));

	TypeList.push_back(T41Var_s("qtch", "Atom Charges", "Atom Charges", FALSE));

	TypeList.push_back(T41Var_s("DensityGradX", "X Density Gradient", "Density Gradient", FALSE));

	TypeList.push_back(T41Var_s("DensityGradY", "Y Density Gradient", "Density Gradient", FALSE));

	TypeList.push_back(T41Var_s("DensityGradZ", "Z Density Gradient", "Density Gradient", FALSE));

	TypeList.push_back(T41Var_s("DensityHessXX", "XX Density Hessian", "Density Hessian", FALSE));

	TypeList.push_back(T41Var_s("DensityHessXY", "XY Density Hessian", "Density Hessian", FALSE));

	TypeList.push_back(T41Var_s("DensityHessXZ", "XZ Density Hessian", "Density Hessian", FALSE));

	TypeList.push_back(T41Var_s("DensityHessYY", "YY Density Hessian", "Density Hessian", FALSE));

	TypeList.push_back(T41Var_s("DensityHessYZ", "YZ Density Hessian", "Density Hessian", FALSE));

	TypeList.push_back(T41Var_s("DensityHessZZ", "ZZ Density Hessian", "Density Hessian", FALSE));

	TypeList.push_back(T41Var_s("Fitdensity", "Electron Fit Density", "Electron Fit Density", TRUE));

	TypeList.push_back(T41Var_s("Density", "Electron Density", "Electron Density", TRUE));

	return static_cast<int>(TypeList.size());
}

Boolean_t CreateAtomZonesFromAtomGroupList(vector<AtomGroup_s> const & AtomGroupList,
	vector<string> const & XYZVarNames,
	vector<FieldDataType_e> & VarDataTypes,
	int GroupIndex)
{
	/*
	*	Now create zones for each atom group and set their style settings accordingly
	*/
	Set_pa NewZones = TecUtilSetAlloc(TRUE);
	vector<int> NewZoneNums;
	vector<AtomColor_s> AtomColorList;
	PopulateAtomColorList(AtomColorList);

	Boolean_t IsOk = (AtomGroupList.size() > 0 && XYZVarNames.size() == 3);

	for (int GroupNum = 0; GroupNum < AtomGroupList.size() && IsOk; ++GroupNum){
		IsOk = TecUtilDataSetAddZone(StringMakeValidZoneName(NuclearPositionsZoneNameBase + AtomGroupList[GroupNum].Name).c_str(), AtomGroupList[GroupNum].Positions[0].size(), 1, 1, ZoneType_Ordered, VarDataTypes.size() > 0 ? VarDataTypes.data() : nullptr);

		EntIndex_t ZoneNum;

		if (IsOk) {
			ZoneNum = TecUtilDataSetGetNumZones();
			IsOk = AuxDataZoneSetItem(ZoneNum, DLZoneType, DLZoneTypeNuclearPositions);
		}
		if (IsOk) {
			IsOk = AuxDataZoneSetItem(ZoneNum, DLZoneAtomicSpecies, AtomGroupList[GroupNum].Name);
		}

		if (IsOk){
			FieldData_pa VarRef[3];
			for (int i = 0; i < 3; ++i)
				VarRef[i] = TecUtilDataValueGetWritableNativeRef(ZoneNum, TecUtilVarGetNumByName(XYZVarNames[i].c_str()));
			if (VALID_REF(VarRef[0]) && VALID_REF(VarRef[0]) && VALID_REF(VarRef[0])){
				for (int i = 0; i < 3; ++i){
					vector<ImportType_t> ImportBuffer(AtomGroupList[GroupNum].Positions[i].begin(), AtomGroupList[GroupNum].Positions[i].end());
					// 					TecUtilDataValueArraySetByRef(VarRef[i], 1, AtomGroupList[GroupNum].Count, AtomGroupList[GroupNum].Positions[i].data());
					TecUtilDataValueArraySetByRef(VarRef[i], 1, ImportBuffer.size(), ImportBuffer.data());
				}
			}

			NewZoneNums.push_back(ZoneNum);

			Set_pa TempSet = TecUtilSetAlloc(TRUE);
			TecUtilSetAddMember(TempSet, ZoneNum, TRUE);

			TecUtilZoneSetScatter(SV_SHOW, TempSet, 0.0, TRUE);
			TecUtilZoneSetScatter(SV_COLOR, TempSet, 0.0, GetAtomColor(AtomColorList, AtomGroupList[GroupNum].Name));
			TecUtilZoneSetScatter(SV_FRAMESIZE, TempSet, (AtomGroupList[GroupNum].Charges.size() > 0 ? MIN(2 * sqrt(AtomGroupList[GroupNum].Charges[0]), 8) : 4), 0);
			TecUtilZoneSetScatterSymbolShape(SV_GEOMSHAPE, TempSet, GeomShape_Sphere);

			TecUtilZoneSetMesh(SV_SHOW, TempSet, 0.0, FALSE);
			TecUtilZoneSetContour(SV_SHOW, TempSet, 0.0, FALSE);
			TecUtilZoneSetVector(SV_SHOW, TempSet, 0.0, FALSE);
			TecUtilZoneSetShade(SV_SHOW, TempSet, 0.0, FALSE);
			TecUtilZoneSetEdgeLayer(SV_SHOW, TempSet, 0.0, FALSE);

			TecUtilSetDealloc(&TempSet);
		}
	}

	if (!NewZoneNums.empty()){
		for (int zi : NewZoneNums){
			TecUtilSetAddMember(NewZones, zi, FALSE);
		}
		ArgList_pa ArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(ArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(ArgList, SV_P2, SV_GROUP);
		TecUtilArgListAppendSet(ArgList, SV_OBJECTSET, NewZones);
		TecUtilArgListAppendArbParam(ArgList, SV_IVALUE, (LgIndex_t)(GroupIndex));
		TecUtilStyleSetLowLevelX(ArgList);
		TecUtilArgListDealloc(&ArgList);

		TecUtilStateChanged(StateChange_ZonesAdded, reinterpret_cast<ArbParam_t>(NewZones));

		TecUtilZoneSetActive(NewZones, AssignOp_PlusEquals);

		TecUtilFrameSetPlotType(PlotType_Cartesian3D);
		DrawBondLinesBetweenNuclearPositions(NewZoneNums, { 1,2,3 }, 3.0, -1.0, 0.1, 0.4, false);


		tecplot::toolbox::StyleValue styleValue;
		styleValue.set((Boolean_t)1, SV_FIELDLAYERS, SV_SHOWSCATTER);
	}
	TecUtilSetDealloc(&NewZones);

	return IsOk;
}


void DrawBondLinesBetweenNuclearPositions(vector<int> const & ZoneNumList, vector<int> XYZVarNums, float BondCutoff, float WDVCutoff, float MinLineThickness, float MaxLineThickness, bool ReplaceOldBonds, bool HeteroBondsOnly){
	REQUIRE(!ZoneNumList.empty());
	for (int const & i : ZoneNumList){
		REQUIRE(i >= 1 && i <= TecUtilDataSetGetNumZones());
		if (!AuxDataZoneItemMatches(i, CSMAuxData.DL.ZoneType, CSMAuxData.DL.ZoneTypeNuclearPositions)){
			TecUtilDialogErrMsg("You can only select nuclear positions for this tool");
			return;
		}
	}
	REQUIRE(XYZVarNums.size() == 3);
	for (int const & i : XYZVarNums) {
		REQUIRE(i >= 1 && i <= TecUtilDataSetGetNumVars());
	}
	REQUIRE(BondCutoff > 0.0);

	float BondCutoffSqr = BondCutoff * BondCutoff;
	float WDVCutoffSqr = WDVCutoff * WDVCutoff;

	TecUtilDataLoadBegin();

	vector<FieldVecPointer_c> ZoneXYZPtrs(ZoneNumList.size());
	vector<string> ZoneNames(ZoneNumList.size());
	for (int zi = 0; zi < ZoneNumList.size(); ++zi){
		ZoneXYZPtrs[zi].InitializeReadPtr(ZoneNumList[zi], XYZVarNums);
		char * TmpChar;
		TecUtilZoneGetName(ZoneNumList[zi], &TmpChar);
		ZoneNames[zi] = TmpChar;
		ZoneNames[zi] = StringRemoveSubString(ZoneNames[zi], NuclearPositionsZoneNameBase);
		TecUtilStringDealloc(&TmpChar);
	}

	Set DeleteZones;

	// Get min distance between any two points
	float MinDistSqr = FLT_MAX, MaxDistSqr = FLT_MIN;
	for (int zi = 0; zi < ZoneXYZPtrs.size(); ++zi) {
		for (int pi = 0; pi < ZoneXYZPtrs[zi].Size(); ++pi) {
			vec3 p1 = ZoneXYZPtrs[zi][pi];
			for (int zj = 0; zj < ZoneXYZPtrs.size(); ++zj) {
				if ((HeteroBondsOnly && zj > zi) || (!HeteroBondsOnly && zj >= zi)) {
					for (int pj = (zj == zi ? pi + 1 : 0); pj < ZoneXYZPtrs[zj].Size(); ++pj) {
						if (zi != zj || pi != pj) {
							vec3 p2 = ZoneXYZPtrs[zj][pj];
							auto PtDistSqr = DistSqr(p1, p2);
							MinDistSqr = MIN(MinDistSqr, PtDistSqr);
							MaxDistSqr = MAX(MaxDistSqr, PtDistSqr);
						}
					}
				}
			}
		}
	}

	float MaxDist = BondCutoff, MinDist = sqrt(MinDistSqr);
	float DistRange = MaxDist - MinDist;
	float LineThicknessRange = MaxLineThickness - MinLineThickness;

	// Now make bond line zones for each pair
	for (int zi = 0; zi < ZoneXYZPtrs.size(); ++zi){
		for (int pi = 0; pi < ZoneXYZPtrs[zi].Size(); ++pi) {
			vec3 p1 = ZoneXYZPtrs[zi][pi];
			for (int zj = 0; zj < ZoneXYZPtrs.size(); ++zj) {
				if ((HeteroBondsOnly && zj > zi) || (!HeteroBondsOnly && zj >= zi)) {
					for (int pj = (zj == zi ? pi + 1 : 0); pj < ZoneXYZPtrs[zj].Size(); ++pj) {
						if (zi != zj || pi != pj) {
							vec3 p2 = ZoneXYZPtrs[zj][pj];
							auto PtDistSqr = DistSqr(p1, p2);

							if (PtDistSqr < BondCutoffSqr || (WDVCutoff > 0.0 && PtDistSqr < WDVCutoffSqr)) {
								string ZoneName = GUIBondZoneNameBase;
								int ZoneNum = 0;
								if (ZoneNames[zi][0] <= ZoneNames[zj][0]){
									ZoneName += ZoneNames[zi] + to_string(pi + 1) + "-" + ZoneNames[zj] + to_string(pj + 1);
								}
								else{
									ZoneName += ZoneNames[zj] + to_string(pj + 1) + "-" + ZoneNames[zi] + to_string(pi + 1);
								}
								int OldZoneNum = ZoneNumByName(ZoneName);
								if (OldZoneNum > 0) {
									if (ReplaceOldBonds) {
										DeleteZones.add(OldZoneNum);
									}
									else{
										continue;
									}
								}
								if (PtDistSqr < BondCutoffSqr) {
									// determine line thickness
									float PtDist = sqrt(PtDistSqr);
									float LineThickness = CLAMP(MinLineThickness + (MaxDist - PtDist) / DistRange * LineThicknessRange, MinLineThickness, MaxLineThickness);

									if (ZoneNames[zi] == "H" || ZoneNames[zj] == "H"){
										LineThickness = MinLineThickness;
									}

									// make zone
									ZoneNum = SaveVec3VecAsScatterZone({ p1, p2 }, ZoneName, Black_C, XYZVarNums);
									if (ZoneNum > 0){
										SetZoneStyle({ ZoneNum }, ZoneStyle_Path, Black_C, LineThickness);
									}
								}
								else{
									ZoneNum = SaveVec3VecAsScatterZone({ p1, p2 }, ZoneName, Black_C, XYZVarNums);
									if (ZoneNum > 0) {
										SetZoneStyle({ ZoneNum }, ZoneStyle_Path, Black_C, MinLineThickness);
										TecUtilZoneSetMesh(SV_LINEPATTERN, Set(ZoneNum).getRef(), 0.0, LinePattern_Dashed);
										TecUtilZoneSetMesh(SV_PATTERNLENGTH, Set(ZoneNum).getRef(), 0.8, 0);
									}
								}
								if (ZoneNum > 0){
									AuxDataZoneSetItem(ZoneNum, CSMAuxData.DL.ZoneType, CSMAuxData.DL.ZoneTypeGUIBond);
								}
							}
						}
					}
				}
			}
		}
	}

	if (!DeleteZones.isEmpty()){
		TecUtilDataSetDeleteZone(DeleteZones.getRef());
	}

	TecUtilDataLoadEnd();
}