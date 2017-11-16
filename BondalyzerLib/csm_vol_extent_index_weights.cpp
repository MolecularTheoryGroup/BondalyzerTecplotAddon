
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include <omp.h>

//#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "TECADDON.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"

#include "CSM_VOL_EXTENT_INDEX_WEIGHTS.h"

#include <armadillo>
using namespace arma;

using std::vector;
using std::string;
using std::to_string;
using std::stringstream;


VolExtentIndexWeights_s & VolExtentIndexWeights_s::operator = (const VolExtentIndexWeights_s & rhs){
	if (this == &rhs)
		return *this;

	MaxIJK = rhs.MaxIJK;
	MaxXYZ = rhs.MaxXYZ;
	MinXYZ = rhs.MinXYZ;
	DelXYZ = rhs.DelXYZ;
	BasisVectors = rhs.BasisVectors;
	BasisNormalized = rhs.BasisNormalized;
	BasisInverse = rhs.BasisInverse;
	BasisExtent = rhs.BasisExtent;
	IsPeriodic = rhs.IsPeriodic;
	AddOnID = rhs.AddOnID;

	for (int i = 0; i < 8; ++i){
		Index[i] = rhs.Index[i];
		Weights[i] = rhs.Weights[i];
	}

	return *this;
}
const Boolean_t VolExtentIndexWeights_s::operator == (const VolExtentIndexWeights_s & rhs) const{
	Boolean_t AreEqual = (
		MaxIJK == rhs.MaxIJK &&
		sum(MaxXYZ == rhs.MaxXYZ) == 3 &&
		sum(MinXYZ == rhs.MinXYZ) == 3 &&
		sum(DelXYZ == rhs.DelXYZ ) == 3 &&
		sum(BasisExtent == rhs.BasisExtent) == 3 &&
		sum(sum(BasisVectors == rhs.BasisVectors)) == 9 &&
		sum(sum(BasisNormalized == rhs.BasisNormalized)) == 9 &&
		sum(sum(BasisInverse == rhs.BasisInverse)) == 9 &&
		IsPeriodic == rhs.IsPeriodic &&
		AddOnID == rhs.AddOnID
		);

	if (AreEqual){
		for (int i = 0; i < 8 && AreEqual; ++i){
			AreEqual = (
				Index[i] == rhs.Index[i] &&
				Weights[i] == rhs.Weights[i]
				);
		}
	}

	return AreEqual;
}

const Boolean_t GetVolInfo(const int & VolZoneNum,
	const vector<int> & XYZVarNums,
	const Boolean_t & IsPeriodic,
	VolExtentIndexWeights_s & VolInfo)
{
	TecUtilZoneGetIJK(VolZoneNum, &VolInfo.MaxIJK[0], &VolInfo.MaxIJK[1], &VolInfo.MaxIJK[2]);
	for (int i = 0; i < 3; ++i) REQUIRE(VolInfo.MaxIJK[i] >= 3); // if less than 3 points, can't do numerical gradients (if necessary)

	ZoneXYZVarGetMinMax_Ordered3DZone(XYZVarNums, VolZoneNum, VolInfo.MinXYZ, VolInfo.MaxXYZ);
	VolInfo.DelXYZ = GetDelXYZ_Ordered3DZone(XYZVarNums, VolZoneNum);
	ZoneXYZVarGetBasisVectors_Ordered3DZone(XYZVarNums, VolZoneNum, VolInfo.BasisVectors, VolInfo.BasisExtent);
	VolInfo.BasisInverse = mat33(VolInfo.BasisVectors.i());
	VolInfo.BasisNormalized = mat33(normalise(VolInfo.BasisVectors));
	VolInfo.IsPeriodic = IsPeriodic;

	return TRUE;
}

const Boolean_t SetIndexAndWeightsForPoint(vec3 Point, VolExtentIndexWeights_s & VolZoneInfo)
{

	Boolean_t IsOk = TRUE;
	/*
	*	Transform the point into the UVW coordinate system
	*/
	Point = VolZoneInfo.BasisInverse * (Point - VolZoneInfo.MinXYZ);

	if (VolZoneInfo.IsPeriodic){
		for (int i = 0; i < 3; ++i){
			if (Point[i] < 0.){
				Point[i] += 1.;
			}
			else if (Point[i] > 1.){
				Point[i] -= 1.;
			}
		}
	}
	else{
		/*
		*	Check that current position is in system bounds
		*/
		for (int i = 0; i < 3; ++i){
			if (Point[i] < 0. || Point[i] > 1.){
				Point[i] = MIN(1., MAX(Point[i], 0.));
			}
		}
	}

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
	/*
	*	Get cell and natural coordinates (the index values for the
	*	cell around the point, and the coordinates within the cell
	*	for the exact point for interpolation).
	*	Same as RectGridBrickXYZtoIJKRST in elemshapefunc.cpp
	*/
	double RST[3];
	if (IsOk){
		LgIndex_t IJK[3];
		for (int i = 0; i < 3; ++i){
// 			double TempCoord = 1.0 + static_cast<double>(VolZoneInfo.MaxIJK[i] - 1.0) * (Point[i] - VolZoneInfo.MinUVW[i]) / (VolZoneInfo.BasisExtent[i]);
			double TempCoord = 1.0 + static_cast<double>(VolZoneInfo.MaxIJK[i] - 1) * Point[i];
			IJK[i] = MAX(MIN(static_cast<LgIndex_t>(TempCoord), VolZoneInfo.MaxIJK[i] - 1), 1);
			RST[i] = 2.0 * (TempCoord - static_cast<double>(IJK[i])) - 1.0;

			if (IJK[i] < 1 || IJK[i] > VolZoneInfo.MaxIJK[i] || RST[i] < -1.0001 || RST[i] > 1.0001){
				IsOk = FALSE;
				break;
			}
		}
		VolZoneInfo.Index[0] = IndexFromIJK(IJK[0], IJK[1], IJK[2], VolZoneInfo.MaxIJK[0], VolZoneInfo.MaxIJK[1]) - 1;
		VolZoneInfo.Index[1] = IndexFromIJK(IJK[0] + 1, IJK[1], IJK[2], VolZoneInfo.MaxIJK[0], VolZoneInfo.MaxIJK[1]) - 1;
		VolZoneInfo.Index[2] = IndexFromIJK(IJK[0] + 1, IJK[1] + 1, IJK[2], VolZoneInfo.MaxIJK[0], VolZoneInfo.MaxIJK[1]) - 1;
		VolZoneInfo.Index[3] = IndexFromIJK(IJK[0], IJK[1] + 1, IJK[2], VolZoneInfo.MaxIJK[0], VolZoneInfo.MaxIJK[1]) - 1;
		VolZoneInfo.Index[4] = IndexFromIJK(IJK[0], IJK[1], IJK[2] + 1, VolZoneInfo.MaxIJK[0], VolZoneInfo.MaxIJK[1]) - 1;
		VolZoneInfo.Index[5] = IndexFromIJK(IJK[0] + 1, IJK[1], IJK[2] + 1, VolZoneInfo.MaxIJK[0], VolZoneInfo.MaxIJK[1]) - 1;
		VolZoneInfo.Index[6] = IndexFromIJK(IJK[0] + 1, IJK[1] + 1, IJK[2] + 1, VolZoneInfo.MaxIJK[0], VolZoneInfo.MaxIJK[1]) - 1;
		VolZoneInfo.Index[7] = IndexFromIJK(IJK[0], IJK[1] + 1, IJK[2] + 1, VolZoneInfo.MaxIJK[0], VolZoneInfo.MaxIJK[1]) - 1;
	}

	/*
	*	Get interpolation weights for interpolating
	*	values of X,Y,Z,Rho,Grad in cell.
	*	Same as BrickTrilinearWeight() in elemshapefunc.cpp
	*/
	if (IsOk){
		double TetraTol = 1.0005;
		double OneMinusRST[3];
		double OnePlusRST[3];
		for (int i = 0; i < 3; ++i){
			if (RST[i] > -TetraTol && RST[i] < TetraTol){
				OneMinusRST[i] = 1.0 - RST[i];
				OnePlusRST[i] = 1.0 + RST[i];
			}
			else{
				IsOk = FALSE;
				break;
			}
		}
		if (IsOk){
			VolZoneInfo.Weights[0] = 0.125 * OneMinusRST[0] * OneMinusRST[1] * OneMinusRST[2];
			VolZoneInfo.Weights[1] = 0.125 * OnePlusRST[0] * OneMinusRST[1] * OneMinusRST[2];
			VolZoneInfo.Weights[2] = 0.125 * OnePlusRST[0] * OnePlusRST[1] * OneMinusRST[2];
			VolZoneInfo.Weights[3] = 0.125 * OneMinusRST[0] * OnePlusRST[1] * OneMinusRST[2];
			VolZoneInfo.Weights[4] = 0.125 * OneMinusRST[0] * OneMinusRST[1] * OnePlusRST[2];
			VolZoneInfo.Weights[5] = 0.125 * OnePlusRST[0] * OneMinusRST[1] * OnePlusRST[2];
			VolZoneInfo.Weights[6] = 0.125 * OnePlusRST[0] * OnePlusRST[1] * OnePlusRST[2];
			VolZoneInfo.Weights[7] = 0.125 * OneMinusRST[0] * OnePlusRST[1] * OnePlusRST[2];
		}
	}


	return IsOk;
} //	const Boolean_t SetIndexAndWeightsForPoint() 

const vector<int> GetIJKForPoint(vec3 & Point, VolExtentIndexWeights_s & VolZoneInfo)
{
	Boolean_t IsOk = TRUE;

	vector<int> IJK(3);
	/*
	*	Transform the point into the UVW coordinate system
	*/
	Point = VolZoneInfo.BasisInverse * (Point - VolZoneInfo.MinXYZ);

	if (VolZoneInfo.IsPeriodic){
		for (int i = 0; i < 3; ++i){
			if (Point[i] < 0.){
				Point[i] += 1.;
			}
			else if (Point[i] > 1.){
				Point[i] -= 1.;
			}
		}
	}
	else{
		/*
		*	Check that current position is in system bounds
		*/
		for (int i = 0; i < 3; ++i){
			if (Point[i] < 0. || Point[i] > 1.){
				Point[i] = MIN(1., MAX(Point[i], 0.));
			}
		}
	}

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
	if (IsOk){
		for (int i = 0; i < 3; ++i){
			// 			double TempCoord = 1.0 + static_cast<double>(VolZoneInfo.MaxIJK[i] - 1.0) * (Point[i] - VolZoneInfo.MinUVW[i]) / (VolZoneInfo.BasisExtent[i]);
			double TempCoord = 1.0 + static_cast<double>(VolZoneInfo.MaxIJK[i] - 1) * Point[i];
			IJK[i] = MAX(MIN(static_cast<LgIndex_t>(TempCoord), VolZoneInfo.MaxIJK[i] - 1), 1);

			if (IJK[i] < 1 || IJK[i] > VolZoneInfo.MaxIJK[i]){
				IsOk = FALSE;
				break;
			}
		}
	}

	Point = VolZoneInfo.BasisVectors * Point + VolZoneInfo.MinXYZ;

	return IJK;

// 	Boolean_t IsOk = TRUE;
// 
// 	vector<int> IJK(3);
// 
// 	/*
// 	*	Check that current position is in system bounds
// 	*/
// 	for (int i = 0; i < 3; ++i){
// 		if (Point[i] < VolZoneInfo.MinXYZ[i] || Point[i] > VolZoneInfo.MaxXYZ[i]){
// 			Point[i] = MIN(VolZoneInfo.MaxXYZ[i], MAX(Point[i], VolZoneInfo.MinXYZ[i]));
// 		}
// 	}
// 
// 	if (IsOk){
// 		for (int i = 0; i < 3; ++i){
// 			double TempCoord = 1.0 + static_cast<double>(VolZoneInfo.MaxIJK[i] - 1.0) * (Point[i] - VolZoneInfo.MinXYZ[i]) / (VolZoneInfo.MaxXYZ[i] - VolZoneInfo.MinXYZ[i]);
// 			IJK[i] = MAX(MIN(static_cast<LgIndex_t>(TempCoord), VolZoneInfo.MaxIJK[i] - 1), 1);
// 
// 			if (IJK[i] < 1 || IJK[i] > VolZoneInfo.MaxIJK[i]){
// 				IsOk = FALSE;
// 				break;
// 			}
// 		}
// 	}
// 
// 	return IJK;
}

/*
*	Given the corner number, and the IJK values of the 0th corner,
*	increment the IJK values as necessary to get the IJK values of
*	the corner specified.
*/
void GetCellCornerIndices(const int & CornerNum, int & i, int & j, int & k){
	switch (CornerNum){
		case 0:
			break;
		case 1:
			i++;
			break;
		case 2:
			i++;
			j++;
			break;
		case 3:
			j++;
			break;
		case 4:
			k++;
			break;
		case 5:
			i++;
			k++;
			break;
		case 6:
			i++;
			j++;
			k++;
			break;
		case 7:
			j++;
			k++;
			break;
	}
}

const double ValByCurrentIndexAndWeightsFromRawPtr(const VolExtentIndexWeights_s & VolZoneInfo, const FieldDataPointer_c & FDPtr)
{
	double Value = 0.0;

	for (int i = 0; i < 8; ++i){
		//TecUtilDialogMessageBox(string("index " + to_string(VolZoneInfo.Index[i]) + ", weight " + to_string(VolZoneInfo.Weights[i]) + ", pointer good " + to_string(FDPtr.IsReady())).c_str(), MessageBoxType_Information);
		Value += VolZoneInfo.Weights[i] * static_cast<double>(FDPtr[VolZoneInfo.Index[i]]);
	}

	return Value;
}

const double ValAtPointByPtr(vec3 & Point, VolExtentIndexWeights_s & VolZoneInfo, const FieldDataPointer_c & FDPtr){
	if (SetIndexAndWeightsForPoint(Point, VolZoneInfo)){
		return ValByCurrentIndexAndWeightsFromRawPtr(VolZoneInfo, FDPtr);
	}
}