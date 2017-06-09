#pragma once
#ifndef CSMFEVOLUME_H_
#define CSMFEVOLUME_H_

#include <vector>

#include "CSM_GRAD_PATH.h"


using std::vector;


class FEVolume_c
{
public:
	FEVolume_c();
	/*
	*	Constructor for making FEVolme_c from
	*	existing FE volume zone. Can be FEQuad or FETriangle.
	*	These are the FEVolume_c's that will integrate over them
	*	selves.
	*/
	FEVolume_c(const int & InZoneNum,
		const ZoneType_e & InZoneType,
		const vector<FieldDataPointer_c> & InXYZPtrs,
		const vec3 & InMaxXYZ,
		const vec3 & InMinXYZ,
		const vector<int> InMaxIJK,
		NodeMap_t* InConnectivityListPtr);
	/*
	*	This is a constructor for making a FEVolume_c
	*	from an existing zone. The FEVolume_c gets all the
	*	information it needs from Tec360.
	*/
	FEVolume_c(const int & ZoneNum,
		const int & VolZoneNum,
		const vector<int> & InXYZVarNums,
		const vector<int> & InIntVarNums);

	/*
	*	Make FEVolume from an existing zone.
	*	Can be IJ ordered zone or actual FE zone.
	*	If IJ ordered zone, need to convert to FE triangle zone representation
	*		with connectivity info and such.
	*	This assumes the structure of the surfaces made by Bondalyzer
	*/
	FEVolume_c(const int & ZoneNum,
		const vector<int> & InXYZVarNums);


	~FEVolume_c();

	/*
	*	Getters
	*/
	const int GetNumSides() const { return m_GPList.size(); }
	const int GetGPZoneNum(const int & i) const { REQUIRE(0 <= i && i < m_GPList.size()); return m_GPList[i].GetZoneNum(); }
	const Boolean_t IsMade() const { return m_FEVolumeMade; }
	const Boolean_t IntResultsReady() const { return m_IntegrationResultsReady; }
	const vector<double> GetIntResults() const;
	const int GetZoneNum() const { return m_ZoneNum; }
	const vector<vector<double> > GetTriSphereIntValsByElem(vector<double> * SphereTriangleAreas = NULL) const { return TriSphereIntValsByElem(SphereTriangleAreas); }

	/*
	*	Setters
	*/
	const Boolean_t Setup(const int & InZoneNum,
		const int & VolZoneNum,
		const vector<int> & InXYZVarNums,
		const vector<int> & InIntVarNums,
		const bool CopyData = false);
	const Boolean_t Setup(const int InZoneNum,
		const vector<int> & InXYZVarNums);

	void AddGP(const GradPath_c & GP){ m_GPList.push_back(GP); }

	const Boolean_t DoIntegration(const Boolean_t & IntegrateVolume, 
		const vector<vec> & stuW, 
		const vector<int> & SplitPtNums = vector<int>(),
		const vector<vec> & stuW2 = vector<vec>());
// 	const Boolean_t DoIntegration(const int & ResolutionScale, const Boolean_t & IntegrateVolume);
	const Boolean_t DoIntegrationNew(const int & ResolutionScale, const Boolean_t & IntegrateVolume);

	const Boolean_t GQIntegration(const int & NumGQPts, const vector<FieldDataPointer_c> & InIntFDPtrs, const Boolean_t & IntegrateVolume);

	const double IntVolume(const int & N, const vec3 & StartPoint) const;
	const vector<mat> GetIntegrationPointsWeights(const vec3 & StartPoint, const vector<vec> & stuW) const;
	const vector<mat> GetIntegrationPointsWeights(const vector<vec> & stuW) const;

	const Boolean_t Make(vector<GradPath_c*> GPs);

	const Boolean_t Refine();
	/*
	*	Two methods to make the 3- and 4-sided
	*	FE volumes from gradient paths.
	*/
	const Boolean_t Make(const GradPath_c & GP1,
		const GradPath_c & GP2,
		const GradPath_c & GP3);
	const Boolean_t Make(const GradPath_c & GP1,
		const GradPath_c & GP2,
		const GradPath_c & GP3,
		const GradPath_c & GP4);

	/*
	*	Other
	*/
	const int SaveAsTriFEZone(const vector<int> & XYZVarNums, string ZoneName = "");
	const int SaveAsTriFEZone(const string & ZoneName, 
		vector<FieldDataType_e> DataTypes,
		const vector<ValueLocation_e> & DataLocations,
		const vector<int> & XYZVarNums) const;
	const int SaveAsFEZone(
		vector<FieldDataType_e> DataTypes,
		const vector<int> & XYZVarNums,
		const int & RhoVarNum
		) const;

	friend class Domain_c;

private:

	const Boolean_t PointIsInterior(const vec3 & Point, const vector<vec3> & FarPoints) const;
	const int TriangleIntersect(const vec3 & T_P0,
		const vec3 & T_P1,
		const vec3 & T_P2,
		const vec3 & R_P1,
		const vec3 & R_P0) const;
	const vector<vector<double> > TriSphereIntValsByElem(vector<double> * SphereTriangleAreas = NULL) const;
	const Boolean_t CalcMaxNodeDistSqr();
	const Boolean_t DistSqrToSurfaceNodeWithinTolerance(const vec3 & CheckPt,
		double & NewDistSqrUnderTol,
		int & CloseNodeNum,
		const double & DistSqrTol = -1.0);
	const Boolean_t SubDivideIntegrateCellAtPoint(const vec3 & Point,
		const vector<vec3> & FarPoints,
		const vec3 & DelXYZ,
		const double & MinDistSqrToSurfaceNode,
		int MinDistNodeNum,
		const int & SubDivideLevel,
		const Boolean_t & IntegrateVolume);
	const vector<int> TriangleEdgeMidPointSubdivide(const int & TriNum);
	void RefineTriElems(const vector<int> & TriNumList);
	void TriPolyLines();
	void RemoveDupicateNodes();


	vector<GradPath_c> m_GPList;

	vector<vec3> m_XYZList;

	vector<double> m_RhoList;

	NodeMap_t* m_ConnectivityListPtr;
	vector<vector<LgIndex_t> > m_ConnectivityList;
	vector<double> m_MaxNeighborNodeDistSqr,
		m_MinNeighborNodeDistSqr;
	vector<vec3> m_RefinedXYZList;
	vector<vector<int> > m_ElemList;

	vector<FieldDataPointer_c> m_XYZPtrs;

	vector<FieldDataPointer_c> m_IntVarPtrs;

#ifdef _DEBUG
	vector<vec3> m_InteriorPts, m_CheckPts, m_InteriorSubPts, m_CheckSubPts;
#endif

	int m_ZoneNum;

	double m_MaxNodeDistanceSqr,
		m_MinNodeDistanceSqr;

	vector<double> m_IntValues;
	vector<int> m_IntVarNums;
	int m_NumIntVars;

	int m_NumGPs;
	int m_NumGPPts;
	int m_NumNodes;
	int m_NumElems;

	int m_NumNodesPerElem;

	vec3 m_ZoneMinXYZ;;
	vec3 m_ZoneMaxXYZ;

	VolExtentIndexWeights_s m_VolZoneInfo;

	Boolean_t m_FEVolumeMade;
	Boolean_t m_IntegrationResultsReady;
};


const vector<vec> GetWeightsPoints(const int & N);


class Domain_c
{
public:
	Domain_c(){}
	Domain_c(const vector<int> & V, FEVolume_c *Vol){ Setup(V, Vol); }
	~Domain_c(){ m_Vol = NULL; }
	void Setup(const vector<int> & V, FEVolume_c *Vol);

	const double Weight() const;
protected:
	
private:
	const double Split(const int & ei, vector<int> & t, Domain_c & D1, Domain_c & D2) const;
	const double WeightFunc(const vector<int> & t) const;


	vector<int> m_V;
	vector<vector<int> > m_E;
	FEVolume_c *m_Vol = NULL;
};



#endif