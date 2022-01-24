#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <iomanip>

// #include "CSM_CRIT_POINTS.h"
#include "CSM_DATA_TYPES.h"


// for profiling
#include <chrono>
#include <ctime>
#include <ratio>

#include <armadillo>
using namespace arma;

//for profiling
using std::chrono::high_resolution_clock;

using std::string;
using std::vector;
using std::to_string;

double const DefaultRhoCutoff = 1e-5;

string const T21Prefix = "T21.";

string const NuclearPositionsZoneNameBase = "NUCLEI ";
string const GUIBondZoneNameBase = "GUIBOND ";

/*
*	These are the aux data tags and values (where applicable)
*	used by the Bondalyzer add-ons.
*/
struct CSMAuxData_s{
	/*
	*	Aux data tags and values for the DataLoader add-on
	*/
	struct DL_s{
		string const Prefix = "DataLoader.",
			ProgramName = Prefix + "ProgramName",
			LatticeConstants = Prefix + "LatticeConstants",
			Origin = Prefix + "Origin",
			ZoneType = Prefix + "ZoneType",
			ZoneTypeVolumeZone = "Volume Zone",
			ZoneTypeNuclearPositions = "Nuclear Positions",
			ZoneTypeGUIBond = "GUI Bond",
			ZoneAtomicSpecies = Prefix + "Atomic Species";
		vector<string> LatticeVecs;

		DL_s(){
			LatticeVecs = {
				Prefix + "LatticeVectorA",
				Prefix + "LatticeVectorB",
				Prefix + "LatticeVectorC"
			};
		}
	}; 
	DL_s const DL;
	/*
	*	Aux data tags and values for the GradientBundleAnalysis add-on
	*/
	struct GBA_s{
		string const Prefix = "GBA.",
			ElemNum = Prefix + "Elem",
			GPClosestPtNumToCP = Prefix + "ClosestPointNumToCP",
			SphereCPName = Prefix + "SphereCP",
			SphereCPNum = Prefix + "SphereCPNumber",
			SourceZoneNum = Prefix + "SourceZoneNumber",
			SourceNucleusName = Prefix + "SourceNucleusName",
			SphereSeedRadius = Prefix + "SphereSeedRadius",
			RadialApprxRadius = Prefix + "RadialApprxRadius",
			SphereSeedRadiusCPDistRatio = Prefix + "SphereSeedRadiusCPDistRatio",
			RadialApprxRadiusCPDistRatio = Prefix + "RadialApprxRadiusCPDistRatio",
			SphereConstrainedNodeNums = Prefix + "SphereConstrainedNodeNums",
			SphereConstrainedNodeIntersectCPTypes = Prefix + "SphereConstrainedNodeIntersectCPTypes",
			SphereConstrainedNodeIntersectCPNames = Prefix + "SphereConstrainedNodeIntersectCPNames",
			SphereConstrainedNodeIntersectCPTotalOffsetNames = Prefix + "SphereConstrainedNodeIntersectCPTotalOffsetNames",
			GPElemNums = Prefix + "SphereElemNums",
			GPNodeNum = Prefix + "SphereNodeNum",
			ZoneType = Prefix + "ZoneType",
			ZoneTypeSphereZone = "SphereMeshZone",
			ZoneTypeDGB = "DifferentialGradientBundle",
			ZoneTypeGradPath = "GradPath",
			ZoneTypeIBEdgeZone = "IBEdgeZone",
			ZoneTypeTopoCageWedge = "TopologicalCageWedge",
			ZoneTypeTruncIsosurface = "TruncatingIsosurface",
			ZoneTypeCondensedRepulsiveBasin = "CondensedRepulsiveBasin",
			ZoneTypeCondensedAttractiveBasin = "CondensedAttractiveBasin",
			ZoneTypeCondensedRepulsiveBasinWedge = "CondensedRepulsiveWedge",
			ZoneTypeCondensedAttractiveBasinWedge = "CondensedAttractiveWedge",
			CondensedBasinIsSpecialGradientBundle = Prefix + "IsSpecialGradientBundle",
			CondensedBasinName = Prefix + "CondensedBasinName",
			CondensedBasinInfo = Prefix + "CondensedBasinInfo",
			CondensedBasinSphereNodes = Prefix + "BasinNodes",
			CondensedBasinSphereElements = Prefix + "BasinElements",
			CondensedBasinCentralElementIndex = Prefix + "BasinCentralElementIndex",
			CondensedBasinDefiningVariable = Prefix + "CondensedBasinDefiningVariable",
			CondensedBasinBoundaryIntVals = Prefix + "CondensedBasinBoundaryIntVals",
			VolumeCPName = Prefix + "VolumeCP",
			IntPrecision = Prefix + "IntegrationPrecision",
			RadialSphereApproximation = Prefix + "RadialSphereApproximation",
			NumGBs = Prefix + "NumberOfGradientBundles",
			NumNodeGPs = Prefix + "NumberOfNodalGradientPaths",
			NumEdgeGPs = Prefix + "NumberOfEdgeGradientPaths",
			TotNumGPs = Prefix + "TotalNumberOfGradientPaths",
			NumWorkingGPs = Prefix + "NumberOfWorkingGradientPaths",
			PointsPerGP = Prefix + "PointsPerGradientPath",
			GPRhoCutoff = Prefix + "GradientPathRhoTruncationValue",
			SphereMinNumGPs = Prefix + "MinimumNumberOfGradientPathsSetting",
			SphereSubdivisionLevel = Prefix + "GBSubdivisionLevel",
			SphereSubdivisionTightness = Prefix + "GBSubdivisionTightness",
			SphereNumBondPathCoincidentGBs = Prefix + "NumBondPathCoincidentGBs",
			GBSurfaceGPMaxSpacing = Prefix + "GBSurfaceGPMaxSpacing",
			GPsPerGB = Prefix + "GradientPathsPerGradientBundle",
			IntWallTime = Prefix + "IntegrationWallTime",
			IntVarNames = Prefix + "IntegratedVariableNames",
			IntVarVals = Prefix + "IntegratedVariableValues",
			IntVarSphereInteriorVals = Prefix + "IntegratedVariableSphereInteriorValues",
			IntVarSphereInteriorMins = Prefix + "IntegratedVariableSphereInteriorMinValues",
			IntVarSphereInteriorMaxes = Prefix + "IntegratedVariableSphereInteriorMaxValues",
			IntVarSphereInteriorMeans = Prefix + "IntegratedVariableSphereInteriorMeanValues",
			IntVarSphereInteriorStdDevs = Prefix + "IntegratedVariableSphereInteriorValueStdDevs",
			IntVarSphereInteriorRanges = Prefix + "IntegratedVariableSphereInteriorValueRanges",
			AtomicBasinIntegrationVariables = Prefix + "AtomicBasinIntegrateVariables",
			AtomicBasinIntegrationValues = Prefix + "AtomicBasinIntegrationValues",
			AtomicBasinIntegrationTotalValues = Prefix + "AtomicBasinTotalIntegrationValues",
			TePartitionFunction = Prefix + "TePartitionFunction",
			VarType = Prefix + "VarType",
			VarTypeCondensedVariable = "CondensedVariable",
			SphereElemSymbol = Prefix + "SphereElementSymbol",
			SphereAtomicNumber = Prefix + "SphereAtomicNumber",
			SphereAtomicCoreECount = Prefix + "SphereAtomicCoreElectronCount",
			SphereOrigin = Prefix + "SphereOriginXYZ";
		vector<string> NodeNums;

		GBA_s(){
			NodeNums = {
				Prefix + "Node1",
				Prefix + "Node2",
				Prefix + "Node3"
			};
		}
	}; 
	GBA_s const GBA;
	/*
	*	Aux data tags and values for the Bondalyzer (CompChem) add-on
	*/
	struct CC_s{
		string const Prefix = "CompChem.",
		ZoneType = Prefix + "ZoneType",
		ZoneSubType = Prefix + "ZoneSubType",
		ZoneTypeCPs = "CriticalPoints",
		ZoneTypeCPsAll = "AllCPs",
		ZoneTypeGP = "GradientPath",
		ZoneTypeSZFS = "SpecialZeroFluxSurface",
		ZoneSubTypeIAS = "InteratomicSurface",
		ZoneSubTypeIASSegment = "InteratomicSurfaceSegment",
		ZoneSubTypeRS = "RingSurface",
		ZoneSubTypeRSSegment = "RingSurfaceSegment",
		ZoneSubTypeBBS = "BondBundleSurface",
		ZoneSubTypeRBS = "RingBundleSurface",
		ZoneSubTypeCBS = "CageBundleSurface",
		ZoneSubTypeBondPath = "BondPath",
		ZoneSubTypeBondPathSegment = "BondPathSegment",
		ZoneSubTypeBondCagePath = "CageBondPath",
		ZoneSubTypeRingLine = "RingLine",
		ZoneSubTypeRingLineSegment = "RingLineSegment",
		ZoneSubTypeCageNuclearPath = "CageNuclearPath",
		ZoneSubTypeRingNuclearPath = "RingNuclearPath",
		ZoneSubTypeRingBondPath = "RingBondPath",
		GPEndNumList = Prefix + "BegEndCrtPtNums",
		GPEndTypeList = Prefix + "BegEndCrtPtTypes",
		GPMiddleCPNum = Prefix + "MiddleCrtPtNum",
		GPMiddleCPNumForType = Prefix + "MiddleCrtPtNumForType",
		GPMiddleCPType = Prefix + "MiddleCrtPtType",
		ZFSCornerCPTypeList = Prefix + "CornerCrtPtTypes",
		ZFSCornerCPNumList = Prefix + "CornerCrtPtNums";
		vector<string> NumCPs,
			CPSubTypes,
			GPEndNumStrs,
			GPEndTypes,
			ZFSCornerCPNumStrs,
			ZFSCornerCPTypes,
			MolecularGraphGPSubTypes,
			OneSkeletonGPSubTypes;

		CC_s(){
			NumCPs = {
				Prefix + "NumCrtPtAtom",
				Prefix + "NumCrtPtBond",
				Prefix + "NumCrtPtRing",
				Prefix + "NumCrtPtCage",
				Prefix + "NumCrtPtRingFF",
				Prefix + "NumCrtPtCageFF"
			};
			CPSubTypes = {
				"NuclearCPs",
				"BondCPs",
				"RingCPs",
				"CageCPs",
				"RingFFCPs",
				"CageFFCPs"
			};
			// For gradient path zone
			GPEndNumStrs = {
				Prefix + "BegCrtPtNum",
				Prefix + "EndCrtPtNum"
			};
			GPEndTypes = {
				Prefix + "BegCrtPtType",
				Prefix + "EndCrtPtType"
			};
			MolecularGraphGPSubTypes = { ZoneSubTypeBondPathSegment };
			OneSkeletonGPSubTypes = {
				ZoneSubTypeBondPathSegment,
				ZoneSubTypeBondCagePath,
				ZoneSubTypeRingLineSegment,
				ZoneSubTypeCageNuclearPath,
				ZoneSubTypeRingNuclearPath,
				ZoneSubTypeRingBondPath
			};
			// for zero-flux surface zone
			for (int i = 0; i < 3; ++i){
				ZFSCornerCPNumStrs.push_back(Prefix + "CornerCrtPtNum" + to_string(i + 1));
				ZFSCornerCPTypes.push_back(Prefix + "CornerCrtPtType" + to_string(i + 1));
			}
		}
	}; 
	CC_s const CC;

	CSMAuxData_s(){}
}; 
CSMAuxData_s const CSMAuxData;

/*
 *	DataLoader aux data tags
 */

string const DLDataPrefix = "DataLoader.";
string const DLProgramName = DLDataPrefix + "ProgramName";
vector<string> const DLLatticeVecs = {
	DLDataPrefix + "LatticeVectorA",
	DLDataPrefix + "LatticeVectorB",
	DLDataPrefix + "LatticeVectorC"
};
string const DLLatticeConstants = DLDataPrefix + "LatticeConstants";
string const DLOrigin = DLDataPrefix + "Origin";
string const DLZoneType = DLDataPrefix + "ZoneType";
string const DLZoneTypeVolumeZone = "Volume Zone";
string const DLZoneTypeNuclearPositions = "Nuclear Positions";
string const DLZoneAtomicSpecies = DLDataPrefix + "Atomic Species";

/*
 *	Gradient bundle analysis aux data tags
 */



string const GBADataPrefix = "GBA.";
string const GBAElemNum = GBADataPrefix + "Elem";
vector<string> const GBANodeNums = {
	GBADataPrefix + "Node1",
	GBADataPrefix + "Node2",
	GBADataPrefix + "Node3"
};
string const GBAGPClosestPtNumToCP = GBADataPrefix + "ClosestPointNumToCP";
string const GBASphereCPName = GBADataPrefix + "SphereCP";
string const GBASphereCPNum = GBADataPrefix + "SphereCPNumber";
string const GBASphereConstrainedNodeNums = GBADataPrefix + "SphereConstrainedNodeNums";
string const GBASphereConstrainedNodeIntersectCPTypes = GBADataPrefix + "SphereConstrainedNodeIntersectCPTypes";
string const GBASphereConstrainedNodeIntersectCPNames = GBADataPrefix + "SphereConstrainedNodeIntersectCPNames";
string const GBAGPElemNums = GBADataPrefix + "SphereElemNums";
string const GBAGPNodeNum = GBADataPrefix + "SphereNodeNum";
string const GBAZoneType = GBADataPrefix + "ZoneType";
string const GBAZoneTypeSphereZone = "SphereMeshZone";
string const GBAZoneTypeFEVolumeZone = "FEVolumeZone";
string const GBAZoneTypeGradPath = "GradPath";
string const GBAZoneTypeIBEdgeZone = "IBEdgeZone";
string const GBAVolumeCPName = GBADataPrefix + "VolumeCP";

/*
*	Aux data names for bondalyzed files
*/
// For CP zone
string const CCDataNumCPs[4] = {
	"CompChem.NumCrtPtAtom",
	"CompChem.NumCrtPtBond",
	"CompChem.NumCrtPtRing",
	"CompChem.NumCrtPtCage"
};
// For gradiant path zone
string const CCDataGPEndNums[2] = {
	"CompChem.BegCrtPtNum",
	"CompChem.EndCrtPtNum"
};
string const CCDataGPEndTypes[2] = {
	"CompChem.BegCrtPtType",
	"CompChem.EndCrtPtType"
};

/*
 *	Strings for variable names
 */
struct CSMVarName_s{
	string const Dens = "Electron Density",
		DensGrad = "Density Gradient",
		DensHess = "Density Hessian",
		DensLap = "Density Laplacian",
		DensAvgCurv = DensLap + " (Avg. Curvature)",
		DeformDens = "Deformation Density",
		DeformDensFit = DeformDens + " Fit",
		DensFit = Dens + " Fit",
		DensExact = Dens + " Exact",
		Tau = "Kinetic Energy Density",
		CoulombPot = "Coulomb Potential",
		XCPot = "Exchange Correlation Potential",
		ELF = "Electron Localization Function",
		DensGradX = "X " + DensGrad,
		DensGradY = "Y " + DensGrad,
		DensGradZ = "Z " + DensGrad,
		DensGradMag = DensGrad + " Magnitude",
		DensHessXX = "XX " + DensHess,
		DensHessXY = "XY " + DensHess,
		DensHessXZ = "XZ " + DensHess,
		DensHessYY = "YY " + DensHess,
		DensHessYZ = "YZ " + DensHess,
		DensHessZZ = "ZZ " + DensHess,
		ElecStatPot = "Electrostatic Potential",
		ElecStatField = "Electrostatic Field Magnitude",
		ElecStatFieldGrad = ElecStatField + " Gradient (Avg. Curvature)",
		GaussCurvature = "Gaussian Curvature",
		DiaMagnShift = "Diamagnetic Shift (Shielding) (Avg. Curvature)",
		
		CritPointType = "CritPointType";
	vector<string> DensGradVec,
		DensHessTensor,
		EigVecs,
		EigVals,
		EigSys,
		GradDotEigVecs,
		EberlyFunc;

	CSMVarName_s(){
		DensGradVec = {
			DensGradX,
			DensGradY,
			DensGradZ
		};
		DensHessTensor = {
			DensHessXX,
			DensHessXY,
			DensHessXZ,
			DensHessYY,
			DensHessYZ,
			DensHessZZ
		};
		EigVecs = {
			"Eigenvector 1x",
			"Eigenvector 1y",
			"Eigenvector 1z",
			"Eigenvector 2x",
			"Eigenvector 2y",
			"Eigenvector 2z",
			"Eigenvector 3x",
			"Eigenvector 3y",
			"Eigenvector 3z"
		};
		EigVals = {
			"Eigenvalue 1",
			"Eigenvalue 2",
			"Eigenvalue 3"
		};
		EigSys = EigVecs;
		EigSys.insert(EigSys.end(), EigVals.cbegin(), EigVals.cend());
		GradDotEigVecs = {
			"Gradient dot Eigenvector 1",
			"Gradient dot Eigenvector 2",
			"Gradient dot Eigenvector 3"
		};
		EberlyFunc = {
			"Eberly 1-ridge Function",
			"Eberly 2-ridge Function"
		};
	}
};
CSMVarName_s const CSMVarName;


/*
*	Strings for variable names
*/
struct CSMZoneName_s{
	string const Delim = "_",
		FullVolume = "Full Volume Zone" + Delim,

		CriticalPoints = "Critical Points",

		GradientPath = "GP" + Delim,
		SpecialGradientPath = "SGP" + Delim,
		BondPath = "BP" + Delim,
		RingLine = "RL" + Delim,

		ZeroFluxSurface = "ZFS" + Delim,
		SpecialZeroFluxSurface = "SZFS" + Delim,
		InteratomicSurface = "IAS" + Delim,
		RingSurface = "RS" + Delim,
		
		FESurface = "FE Volume";

	vector<string> CPType = vector<string>({ "Nuclear", "Bond", "Ring", "Cage", "Ring FF", "Cage FF" });
	vector<double> CPScatterSize = { 2.5, 1.0, 1.0, 1.1, 0, 0 };
	vector<GeomShape_e> CPScatterSymbol = { GeomShape_Sphere, GeomShape_Sphere, GeomShape_Sphere, GeomShape_Octahedron, GeomShape_Invalid, GeomShape_Invalid };
	CSMZoneName_s(){
		for (string & i : CPType) i = CriticalPoints + Delim + i;
	}
};
CSMZoneName_s const CSMZoneName;


size_t getTotalSystemMemory();

template <class T>
int VectorGetElementNum(vector<T> const & SearchVec, T const & Item){
	int ElemNum;

	vector<T>::const_iterator itr = std::find(SearchVec.cbegin(), SearchVec.cend(), Item);

	if (itr == SearchVec.cend()) return -1;

	ElemNum = itr - SearchVec.cbegin();
	return ElemNum;
}
int SearchVectorForString(vector<string> const & Vec, string const & SearchString, bool UseVectorStringLength = true);

// template <class T>
// vector<T> SplitString(string const &s, string const & delim, bool RemoveAllBlanks = false, bool RemoveBlankAtEnd = false);
vector<string> SplitString(string const &s, string const & delim = " ", bool RemoveAllBlanks = false, bool RemoveBlankAtEnd = false);
vector<int> SplitStringInt(string const &s, string const & delim = ",");
vector<double> SplitStringDbl(string const &s, string const & delim = ",");
string StringJoin(vector<string> const & strVec, string const & delim = " ");

template <class T>
string VectorToString(vector<T> const & Items, string const & Delim = " "){
	//string OutStr;

	//for (int i = 0; i < Items.size(); ++i) {
	//	OutStr += to_string(Items[i]);
	//	if (i < Items.size() - 1) OutStr += Delim;
	//}

	//return OutStr;

	std::ostringstream oss;
	if (!Items.empty()) {
		oss.precision(16);
		std::copy(Items.begin(), Items.end() - 1, std::ostream_iterator<T>(oss, Delim.c_str()));
		oss << Items.back();
		return oss.str();
	}
	else return "";
}
string IntVectorToRangeString(vector<int> Items);
vector<int> RangeStringToIntVector(string const & Str);
Boolean_t StringIsInt(string const & s);
Boolean_t StringIsFloat(string const & s);
string StringReplaceSubString(string const & InStr, string const & OldStr, string const & NewStr);
string StringRemoveSubString(string const & InString, string const & SubString);
string DoubleToString(double const & Val, int Precision = 2);

void AuxDataCopy(int SourceNum, int DestNum, bool IsZone);

Boolean_t AuxDataZoneHasItem(int ZoneNum, string const & AuxDataName);
Boolean_t AuxDataZoneItemMatches(int ZoneNum, string const & AuxDataName, string const & AuxDataValue);
string AuxDataZoneGetItem(int ZoneNum, string const & AuxDataName);
Boolean_t AuxDataZoneGetItem(int ZoneNum, string const & AuxDataName, string & Value);
vector<string> AuxDataZoneGetList(int ZoneNum, string const & AuxDataName);
Boolean_t AuxDataZoneSetItem(int ZoneNum, string const & AuxDataName, string const & AuxDataValue);
template <class T>
Boolean_t AuxDataZoneSetList(int ZoneNum, string const & AuxDataName, vector<T> const & AuxDataValueList){
	AuxDataZoneSetItem(ZoneNum, AuxDataName, VectorToString(AuxDataValueList, ',,'));
}
Boolean_t AuxDataZoneDeleteItemByName(int ZoneNum, string const & AuxDataName);

vector<string> AuxDataDisassembleString(string const & Str, int StrLen = 32000);
string AuxDataAssembleString(vector<string> const & StrVec);

string AuxDataDataSetGetItem(string const & AuxDataName);
Boolean_t AuxDataDataSetGetItem(string const & AuxDataName, string & Value);
Boolean_t AuxDataDataSetSetItem(string const & AuxDataName, string const & AuxDataValue);
Boolean_t AuxDataDataSetDeleteItemByName(string const & AuxDataName);


SmInteger_t VarNumByNameList(vector<string> const & VarNameList, bool const PartialMatch = false);
SmInteger_t VarNumByName(string const & VarName, bool const PartialMatch = false);
SmInteger_t ZoneNumByNameList(vector<string> const & ZoneNameList);
SmInteger_t ZoneNumByName(string const & ZoneName, bool const ActiveOnly = false, bool const PartialMatch = false);

vec3 GetDelXYZ_Ordered3DZone(vector<int> const & XYZVarNums, int ZoneNum);
vector<vec3> ZoneXYZVarGetMinMax_Ordered3DZone(vector<int> const & XYZVarNums, int ZoneNum);
void ZoneXYZVarGetMinMax_Ordered3DZone(vector<int> const & XYZVarNums, int ZoneNum, vec3 & MinXYZ, vec3 & MaxXYZ);
void ZoneXYZVarGetBasisVectors_Ordered3DZone(vector<int> const & XYZVarNums, int ZoneNum, mat33 & BasisVectors, vec3 & BVExtent);

high_resolution_clock::time_point StatusLaunch(string const & StatusStr, AddOn_pa const & AddOnID, Boolean_t ShowScale = TRUE, Boolean_t ShowButton = TRUE, Boolean_t RelaunchStatus = FALSE);
void StatusDrop(AddOn_pa const & AddOnID);
Boolean_t StatusUpdate(unsigned int CurrentNum, 
	unsigned int TotalNum, 
	string const & ProgresssText, 
	AddOn_pa const & AddOnID, 
	high_resolution_clock::time_point StartTime = high_resolution_clock::time_point(),
	bool DisplayCount = true);

LgIndex_t IndexFromIJK(LgIndex_t I,
	LgIndex_t J,
	LgIndex_t K,
	LgIndex_t MaxI,
	LgIndex_t MaxJ);

LgIndex_t IndexFromIJK(LgIndex_t I,
	LgIndex_t    J,
	LgIndex_t    K,
	LgIndex_t IMax,
	LgIndex_t JMax,
	LgIndex_t KMax,
	Boolean_t PeriodicBC);

vector<LgIndex_t> IJKFromIndex(LgIndex_t Index,
	const vector<LgIndex_t>& IJKMax);

int SaveVec3VecAsScatterZone(vector<vec3> const & VecVec,
	string const & ZoneName = "IOrderedZone",
	ColorIndex_t const & Color = Black_C,
	vector<EntIndex_t> const & XYZVarNums = { 1,2,3 });

Boolean_t SaveTetVec3VecAsFEZone(vector<vec3> const & Verts,
	string const & ZoneName,
	ColorIndex_t const & Color,
	vector<EntIndex_t> const & XYZVarNums);

enum CSMZoneStyle_e{
	ZoneStyle_Invalid = -1,
	ZoneStyle_Points,
	ZoneStyle_Path,
	ZoneStyle_Surface,
	ZoneStyle_Volume
};
void SetZoneStyle(vector<int> ZoneNums = {},
	CSMZoneStyle_e ZoneStyle = ZoneStyle_Invalid,
	ColorIndex_t Color = InvalidColor_C,
	double const Size = -1);

int ZoneFinalSourceZoneNum(int ZoneNum, bool MaintainZoneType = false);

void SetZoneNum(int OldZoneNum = TecUtilDataSetGetNumZones(), int NewZoneNum = 1);