#pragma once
#ifndef CSMDATASETINFO_H_
#define CSMDATASETINFO_H_

#include <string>
#include <vector>

#include "CSM_DATA_TYPES.h"

#include <armadillo>
using namespace arma;


using std::string;
using std::vector;

const double DefaultRhoCutoff = 1e-3;

const string T21Prefix = "T21.";

const string CPTypeVarName = "CritPointType";

/*
*	These are the aux data tags and values (where applicable)
*	used by the Bondalyzer add-ons.
*/
struct CSMAuxData_s{
	/*
	*	Aux data tags and values for the DataLoader add-on
	*/
	struct DL_s{
		const string Prefix = "DataLoader.",
			ProgramName = Prefix + "ProgramName",
			LatticeConstants = Prefix + "LatticeConstants",
			Origin = Prefix + "Origin",
			ZoneType = Prefix + "ZoneType",
			ZoneTypeVolumeZone = "Volume Zone",
			ZoneTypeNuclearPositions = "Nuclear Positions",
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
	const DL_s DL;
	/*
	*	Aux data tags and values for the GradientBundleAnalysis add-on
	*/
	struct GBA_s{
		const string Prefix = "GBA.",
			ElemNum = Prefix + "Elem",
			GPClosestPtNumToCP = Prefix + "ClosestPointNumToCP",
			SphereCPName = Prefix + "SphereCP",
			SphereCPNum = Prefix + "SphereCPNumber",
			SourceZoneNum = Prefix + "SourceZoneNumber",
			SphereConstrainedNodeNums = Prefix + "SphereConstrainedNodeNums",
			SphereConstrainedNodeIntersectCPTypes = Prefix + "SphereConstrainedNodeIntersectCPTypes",
			SphereConstrainedNodeIntersectCPNames = Prefix + "SphereConstrainedNodeIntersectCPNames",
			GPElemNums = Prefix + "SphereElemNums",
			GPNodeNum = Prefix + "SphereNodeNum",
			ZoneType = Prefix + "ZoneType",
			ZoneTypeSphereZone = "SphereMeshZone",
			ZoneTypeFEVolumeZone = "FEVolumeZone",
			ZoneTypeGradPath = "GradPath",
			ZoneTypeIBEdgeZone = "IBEdgeZone",
			VolumeCPName = Prefix + "VolumeCP",
			IntPrecision = Prefix + "IntegrationPrecision",
			NumGBs = Prefix + "NumberOfGradientBundles",
			PointsPerGP = Prefix + "PointsPerGradientPath",
			GPsPerGB = Prefix + "GradientPathsPerGradientBundle",
			IntWallTime = Prefix + "IntegrationWallTime";
		vector<string> NodeNums;

		GBA_s(){
			NodeNums = {
				Prefix + "Node1",
				Prefix + "Node2",
				Prefix + "Node3"
			};
		}
	}; 
	const GBA_s GBA;
	/*
	*	Aux data tags and values for the Bondalyzer (CompChem) add-on
	*/
	struct CC_s{
		const string Prefix = "CompChem.",
		ZoneType = Prefix + "ZoneType",
		ZoneSubType = Prefix + "ZoneSubType",
		ZoneTypeCPs = Prefix + "CriticalPoints",
		ZoneTypeGP = Prefix + "GradientPath",
		ZoneSubTypeBondPath = Prefix + "BondPath",
		ZoneSubTypeRingLine = Prefix + "RingLine";
		vector<string> NumCPs,
			CPSubTypes,
			GPEndNums,
			GPEndTypes;

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
				Prefix + "NuclearCPs",
				Prefix + "BondCPs",
				Prefix + "RingCPs",
				Prefix + "CageCPs",
				Prefix + "RingFFCPs",
				Prefix + "CageFFCPs"
			};
			// For gradiant path zone
			GPEndNums = {
				Prefix + "BegCrtPtNum",
				Prefix + "EndCrtPtNum"
			};
			GPEndTypes = {
				Prefix + "BegCrtPtType",
				Prefix + "EndCrtPtType"
			};
		}
	}; 
	const CC_s CC;

	CSMAuxData_s(){}
}; 
const CSMAuxData_s CSMAuxData;

/*
 *	DataLoader aux data tags
 */

const string DLDataPrefix = "DataLoader.";
const string DLProgramName = DLDataPrefix + "ProgramName";
const vector<string> DLLatticeVecs = {
	DLDataPrefix + "LatticeVectorA",
	DLDataPrefix + "LatticeVectorB",
	DLDataPrefix + "LatticeVectorC"
};
const string DLLatticeConstants = DLDataPrefix + "LatticeConstants";
const string DLOrigin = DLDataPrefix + "Origin";
const string DLZoneType = DLDataPrefix + "ZoneType";
const string DLZoneTypeVolumeZone = "Volume Zone";
const string DLZoneTypeNuclearPositions = "Nuclear Positions";
const string DLZoneAtomicSpecies = DLDataPrefix + "Atomic Species";

/*
 *	Gradient bundle analysis aux data tags
 */



const string GBADataPrefix = "GBA.";
const string GBAElemNum = GBADataPrefix + "Elem";
const vector<string> GBANodeNums = {
	GBADataPrefix + "Node1",
	GBADataPrefix + "Node2",
	GBADataPrefix + "Node3"
};
const string GBAGPClosestPtNumToCP = GBADataPrefix + "ClosestPointNumToCP";
const string GBASphereCPName = GBADataPrefix + "SphereCP";
const string GBASphereCPNum = GBADataPrefix + "SphereCPNumber";
const string GBASphereConstrainedNodeNums = GBADataPrefix + "SphereConstrainedNodeNums";
const string GBASphereConstrainedNodeIntersectCPTypes = GBADataPrefix + "SphereConstrainedNodeIntersectCPTypes";
const string GBASphereConstrainedNodeIntersectCPNames = GBADataPrefix + "SphereConstrainedNodeIntersectCPNames";
const string GBAGPElemNums = GBADataPrefix + "SphereElemNums";
const string GBAGPNodeNum = GBADataPrefix + "SphereNodeNum";
const string GBAZoneType = GBADataPrefix + "ZoneType";
const string GBAZoneTypeSphereZone = "SphereMeshZone";
const string GBAZoneTypeFEVolumeZone = "FEVolumeZone";
const string GBAZoneTypeGradPath = "GradPath";
const string GBAZoneTypeIBEdgeZone = "IBEdgeZone";
const string GBAVolumeCPName = GBADataPrefix + "VolumeCP";

/*
*	Aux data names for bondalyzed files
*/
// For CP zone
const string CCDataNumCPs[4] = {
	"CompChem.NumCrtPtAtom",
	"CompChem.NumCrtPtBond",
	"CompChem.NumCrtPtRing",
	"CompChem.NumCrtPtCage"
};
// For gradiant path zone
const string CCDataGPEndNums[2] = {
	"CompChem.BegCrtPtNum",
	"CompChem.EndCrtPtNum"
};
const string CCDataGPEndTypes[2] = {
	"CompChem.BegCrtPtType",
	"CompChem.EndCrtPtType"
};

/*
 *	Strings for variable names
 */
struct CSMVarName_s{
	const string Dens = "Electron Density",
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
			"Eigenvector 1-1",
			"Eigenvector 1-2",
			"Eigenvector 1-3",
			"Eigenvector 2-1",
			"Eigenvector 2-2",
			"Eigenvector 2-3",
			"Eigenvector 3-1",
			"Eigenvector 3-2",
			"Eigenvector 3-3"
		};
		EigVals = {
			"Eigenvalue 1",
			"Eigenvalue 2",
			"Eigenvalue 3"
		};
		EigSys = EigVecs;
		EigSys.insert(EigSys.begin(), EigVals.cbegin(), EigVals.cend());
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
const CSMVarName_s CSMVarName;


size_t getTotalSystemMemory();

const int SearchVectorForString(const vector<string> & Vec, const string & SearchString, const bool & UseVectorStringLength = true);
const vector<string> SplitString(const string &s, char delim);
const Boolean_t StringIsInt(const string & s);
const Boolean_t StringIsFloat(const string & s);

const Boolean_t AuxDataZoneHasItem(const int & ZoneNum, const string & AuxDataName);
const Boolean_t AuxDataZoneItemMatches(const int & ZoneNum, const string & AuxDataName, const string & AuxDataValue);
const string AuxDataZoneGetItem(const int & ZoneNum, const string & AuxDataName);
const Boolean_t AuxDataZoneGetItem(const int & ZoneNum, const string & AuxDataName, string & Value);
const Boolean_t AuxDataZoneSetItem(const int & ZoneNum, const string & AuxDataName, const string & AuxDataValue);
const Boolean_t AuxDataZoneDeleteItemByName(const int & ZoneNum, const string & AuxDataName);

const string AuxDataDataSetGetItem(const string & AuxDataName);
const Boolean_t AuxDataDataSetGetItem(const string & AuxDataName, string & Value);
const Boolean_t AuxDataDataSetSetItem(const string & AuxDataName, const string & AuxDataValue);
const Boolean_t AuxDataDataSetDeleteItemByName(const string & AuxDataName);


SmInteger_t VarNumByNameList(const vector<string> & VarNameList, const bool PartialMatch = false);
SmInteger_t VarNumByName(const string & VarName, const bool PartialMatch = false);
SmInteger_t ZoneNumByNameList(const vector<string> & ZoneNameList);
SmInteger_t ZoneNumByName(const string & ZoneName, const bool ActiveOnly = false, const bool PartialMatch = false);

const vec3 GetDelXYZ_Ordered3DZone(const vector<int> & XYZVarNums, const int & ZoneNum);
const vector<vec3> ZoneXYZVarGetMinMax_Ordered3DZone(const vector<int> & XYZVarNums, const int & ZoneNum);
void ZoneXYZVarGetMinMax_Ordered3DZone(const vector<int> & XYZVarNums, const int & ZoneNum, vec3 & MinXYZ, vec3 & MaxXYZ);

const Boolean_t SetPercent(unsigned int CurrentNum, unsigned int TotalNum, const string & ProgresssText, const AddOn_pa & AddOnID);

const LgIndex_t IndexFromIJK(const LgIndex_t & I,
	const LgIndex_t & J,
	const LgIndex_t & K,
	const LgIndex_t & MaxI,
	const LgIndex_t & MaxJ);

const LgIndex_t IndexFromIJK(const LgIndex_t    & I,
	const LgIndex_t &    J,
	const LgIndex_t &    K,
	const LgIndex_t & IMax,
	const LgIndex_t & JMax,
	const LgIndex_t & KMax,
	const Boolean_t & PeriodicBC);

const vector<LgIndex_t> IJKFromIndex(const LgIndex_t & Index,
	const vector<LgIndex_t>& IJKMax);

const Boolean_t SaveVec3VecAsScatterZone(const vector<vec3> & VecVec,
	const string & ZoneName,
	const ColorIndex_t & Color,
	const vector<EntIndex_t> & XYZVarNums);

#endif