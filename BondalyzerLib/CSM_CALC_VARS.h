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

static vector<string> const CalcVarTypeNames = {
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

void CalcGradGradMagForDataset(Boolean_t IsPeriodic, AddOn_pa const & AddOnID);

Boolean_t CalcGradForRegularVar(vector<int> const & IJKMax,
	vec3 const & DelXYZ,
	Boolean_t IsPeriodic,
	FieldDataPointer_c const & VarReadPtr,
	vector<FieldDataPointer_c> const & VarWritePtrs,
	string const & VarName,
	AddOn_pa const & AddOnID);

Boolean_t CalcMagForRegularVectorVar(vector<int> const & IJKMax,
	vector<FieldDataPointer_c> const & VarReadPtrs,
	FieldDataPointer_c const & VarWritePtr,
	string const & VarName,
	AddOn_pa const & AddOnID);

void CalcGradForNode(int ii,
	int jj,
	int kk,
	vec3 const & DelXYZx2,
	vec3 const & DelXYZx12,
	vector<double> & Vals,
	vector<int> & DirInd,
	int StartDir,
	vector<int> const & IJKMax,
	Boolean_t IsPeriodic,
	FieldDataPointer_c const & VarReadPtr,
	vec3 & OutValues);

void CalcGradForPoint(vec3 const & Point,
	vec3 const & DelXYZ,
	VolExtentIndexWeights_s & VolInfo,
	mat33 const & DirVects,
	int StartDir,
	Boolean_t IsPeriodic,
	vec & OutValues,
	FieldDataPointer_c const & VarReadPtr,
	GPType_e const & CalcType,
	void * Params);

void CalcHessForNode(int ii,
	int jj,
	int kk,
	vec3 const & DelXYZx2,
	vec3 const & DelXYZx12,
	vector<double> & Vals,
	vector<int> & DirInd,
	vector<int> const & IJKMax,
	Boolean_t IsPeriodic,
	FieldDataPointer_c const & ScalarReadPtr,
	vector<FieldDataPointer_c> const & GradReadPtrs,
	mat33 & Hess);

void CalcHessForPoint(vec3 const & Point,
	vec3 const & DelXYZ,
	VolExtentIndexWeights_s & VolInfo,
	mat33 const & DirVects,
	Boolean_t IsPeriodic,
	mat & OutValues,
	FieldDataPointer_c const & VarReadPtr,
	GPType_e const & CalcType,
	void * Params);

void CalcHessFor3DPoint(vec3 const & Point,
	vec3 const & DelXYZ,
	VolExtentIndexWeights_s & VolInfo,
	Boolean_t IsPeriodic,
	mat33 & Hess,
	vector<FieldDataPointer_c> const & VarReadPtrs,
	GPType_e const & CalcType,
	void * Params);

void CalcHessForDataSet(Boolean_t IsPeriodic,
	AddOn_pa const & AddOnID);

Boolean_t CalcEigenSystemForNode(int ii,
	int jj,
	int kk,
	vec3 & EigenValues,
	mat33 & EigenVectors,
	MultiRootParams_s & RootParams);

Boolean_t CalcEigenSystemForPoint(vec3 & Point,
	vec & EigenValues,
	mat & EigenVectors,
	MultiRootParams_s & RootParams);

void CalcEigenSystemForDataSet(Boolean_t IsPeriodic,
	AddOn_pa const & AddOnID);

void CalcEigenvecDotGradForDataSet(Boolean_t IsPeriodic,
	Boolean_t NormalizeGrad,
	AddOn_pa const & AddOnID);

void CalcEigenvecDotGradForPoint(vec3 Point,
	vec3 & DotProducts,
	vec3 & EigenValues,
	mat33 & EigenVectors,
	Boolean_t NormalizeGrad,
	MultiRootParams_s & RootParams);

void CalcEigenvecDotGradForPoint(vec3 Point,
	vec3 & DotProducts,
	Boolean_t NormalizeGrad,
	MultiRootParams_s & RootParams);

double Eberly1RidgeFunction(vec3 & Point,
	double const & RhoCutoff,
	Boolean_t NormalizeGrad,
	MultiRootParams_s & RootParams);
double Eberly2RidgeFunction(vec3 & Point,
	double const & RhoCutoff,
	Boolean_t NormalizeGrad,
	MultiRootParams_s & RootParams);

double NEBForceFunction(vec3 & Point,
	MultiRootParams_s & RootParams);

void CalcEberlyFunctions(Boolean_t IsPeriodic, AddOn_pa const & AddOnID, double const & RhoCutoff);

void MapAllVarsToAllZones(AddOn_pa const & AddOnID);

vector<double> LogLevels(double Min, double Max);

void GaussianBlur(Boolean_t IsPeriodic,
	AddOn_pa const & AddOnID,
	EntIndex_t const & ZoneNum,
	EntIndex_t const & VarNum,
	string const & NewZoneName,
	double const & Sigma);

vec3 Transform2dTo3d(vec2 const & TwoPt,
	mat33 const & BasisVectors,
	vec3 const & Origin);

#endif