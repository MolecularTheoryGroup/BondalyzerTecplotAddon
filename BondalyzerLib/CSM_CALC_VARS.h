#pragma once
#ifndef CSMCALCVARS_H_
#define	CSMCALCVARS_H_

#include <vector>
#include <string>

#include "CSM_DATA_TYPES.h"

#include <armadillo>
using namespace arma;



using std::vector;
using std::string;

enum CalcVar_e{
	CalcGradientVectors = 0,
	CalcGradientMagnitude,
	CalcHessian,
	CalcEigenSystem,
	CalcLaplacian,
	CalcGaussianCurvature,
	CalcEigenVectorsDotGradient,
	CalcEberly1Ridge,
	CalcEberly2Ridge,
	CalcEigenRank,

	CalcInvalidVar = -1
};

static const vector<string> CalcVarTypeNames = {
	"Density gradient vector",
	"Density gradient magnitude",
	"Density Hessian",
	"Eigen system of density Hessian",
	"Laplacian of density",
	"Gaussian curvature of density",
	"Inner product of gradient vector with eigenvectors",
	"Eberly 1-ridge function",
	"Eberly 2-ridge function",
	"Rank of eigenvalues"
};

struct CalcVarsOptions_s{
	Boolean_t IsPeriodic = FALSE;
	Boolean_t CalcForAllZones = TRUE;
	EntIndex_t CalcZoneNum = -1;
	vector<CalcVar_e> CalcVarList;
	vector<Boolean_t> EberlyUseCutoff;
	vector<double> EberlyCutoff;
	EntIndex_t RhoVarNum = -1;
	vector<EntIndex_t> GradVarNums;
	vector<EntIndex_t> HessVarNums;
	Boolean_t HasGrad = FALSE;
	Boolean_t HasHess = FALSE;

	AddOn_pa AddOnID;

	CalcVarsOptions_s(){
		EberlyUseCutoff.resize(2, FALSE);
		EberlyCutoff.resize(2, 0.001);
		GradVarNums.resize(3, -1);
		HessVarNums.resize(6, -1);
	}
};

void CalcVars(CalcVarsOptions_s & Opt);

void CalcGradGradMagForDataset(Boolean_t IsPeriodic, const AddOn_pa & AddOnID);

const Boolean_t CalcGradForRegularVar(const vector<int> & IJKMax,
	const vec3 & DelXYZ,
	const Boolean_t & IsPeriodic,
	const FieldDataPointer_c & VarReadPtr,
	const vector<FieldDataPointer_c> & VarWritePtrs,
	const string & VarName,
	const AddOn_pa & AddOnID);

const Boolean_t CalcMagForRegularVectorVar(const vector<int> & IJKMax,
	const vector<FieldDataPointer_c> & VarReadPtrs,
	const FieldDataPointer_c & VarWritePtr,
	const string & VarName,
	const AddOn_pa & AddOnID);

void CalcGradForNode(const int & ii,
	const int & jj,
	const int & kk,
	const vec3 & DelXYZx2,
	const vec3 & DelXYZx12,
	vector<double> & Vals,
	vector<int> & DirInd,
	const int & StartDir,
	const vector<int> & IJKMax,
	const Boolean_t & IsPeriodic,
	const FieldDataPointer_c & VarReadPtr,
	vec3 & OutValues);

void CalcGradForPoint(const vec3 & Point,
	const vec3 & DelXYZ,
	VolExtentIndexWeights_s & VolInfo,
	const mat33 & DirVects,
	const int & StartDir,
	const Boolean_t & IsPeriodic,
	vec & OutValues,
	const FieldDataPointer_c & VarReadPtr,
	const GPType_e & CalcType,
	void * Params);

void CalcHessForNode(const int & ii,
	const int & jj,
	const int & kk,
	const vec3 & DelXYZx2,
	const vec3 & DelXYZx12,
	vector<double> & Vals,
	vector<int> & DirInd,
	const vector<int> & IJKMax,
	const Boolean_t & IsPeriodic,
	const FieldDataPointer_c & ScalarReadPtr,
	const vector<FieldDataPointer_c> & GradReadPtrs,
	mat33 & Hess);

void CalcHessForPoint(const vec3 & Point,
	const vec3 & DelXYZ,
	VolExtentIndexWeights_s & VolInfo,
	const mat33 & DirVects,
	const Boolean_t & IsPeriodic,
	mat & OutValues,
	const FieldDataPointer_c & VarReadPtr,
	const GPType_e & CalcType,
	void * Params);

void CalcHessFor3DPoint(const vec3 & Point,
	const vec3 & DelXYZ,
	VolExtentIndexWeights_s & VolInfo,
	const Boolean_t & IsPeriodic,
	mat33 & Hess,
	const vector<FieldDataPointer_c> & VarReadPtrs,
	const GPType_e & CalcType,
	void * Params);

void CalcHessForDataSet(Boolean_t IsPeriodic,
	const AddOn_pa & AddOnID);

const Boolean_t CalcEigenSystemForNode(const int & ii,
	const int & jj,
	const int & kk,
	vec3 & EigenValues,
	mat33 & EigenVectors,
	MultiRootParams_s & RootParams);

const Boolean_t CalcEigenSystemForPoint(vec3 & Point,
	vec & EigenValues,
	mat & EigenVectors,
	MultiRootParams_s & RootParams);

void CalcEigenSystemForDataSet(Boolean_t IsPeriodic,
	const AddOn_pa & AddOnID);

void CalcEigenvecDotGradForDataSet(Boolean_t IsPeriodic,
	const Boolean_t & NormalizeGrad,
	const AddOn_pa & AddOnID);

void CalcEigenvecDotGradForPoint(vec3 Point,
	vec3 & DotProducts,
	vec3 & EigenValues,
	mat33 & EigenVectors,
	const Boolean_t & NormalizeGrad,
	MultiRootParams_s & RootParams);

void CalcEigenvecDotGradForPoint(vec3 Point,
	vec3 & DotProducts,
	const Boolean_t & NormalizeGrad,
	MultiRootParams_s & RootParams);

const double Eberly1RidgeFunction(vec3 & Point,
	const double & RhoCutoff,
	const Boolean_t & NormalizeGrad,
	MultiRootParams_s & RootParams);
const double Eberly2RidgeFunction(vec3 & Point,
	const double & RhoCutoff,
	const Boolean_t & NormalizeGrad,
	MultiRootParams_s & RootParams);

const double NEBForceFunction(vec3 & Point,
	MultiRootParams_s & RootParams);

void CalcEberlyFunctions(Boolean_t IsPeriodic, const AddOn_pa & AddOnID, const double & RhoCutoff);

void MapAllVarsToAllZones(const AddOn_pa & AddOnID);

const vector<double> LogLevels(double Min, double Max);

void GaussianBlur(const Boolean_t & IsPeriodic,
	const AddOn_pa & AddOnID,
	const EntIndex_t & ZoneNum,
	const EntIndex_t & VarNum,
	const string & NewZoneName,
	const double & Sigma);

const vec3 Transform2dTo3d(const vec2 & TwoPt,
	const mat33 & BasisVectors,
	const vec3 & Origin);

#endif