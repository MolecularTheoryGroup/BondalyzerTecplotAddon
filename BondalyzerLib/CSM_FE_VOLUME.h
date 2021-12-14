#pragma once

#include <vector>
#include <queue>

#include "CSM_GRAD_PATH.h"
#include "CSM_GEOMETRY.h"

using std::vector;

static int const DefaultProjectPointToSurfaceMaxBFSDepth = 10;
static int const DefaultProjectPointtoSurfaceNumConvergedBSFDepthLevels = 3;

class FESurface_c
{
public:
	FESurface_c();
	/*
	*	Constructor for making FEVolme_c from
	*	existing FE volume zone. Can be FEQuad or FETriangle.
	*	These are the FEVolume_c's that will integrate over them
	*	selves.
	*/
	FESurface_c(int InZoneNum,
		ZoneType_e const & InZoneType,
		vector<FieldDataPointer_c> const & InXYZPtrs,
		vec3 const & InMaxXYZ,
		vec3 const & InMinXYZ,
		vector<int> const InMaxIJK,
		NodeMap_t* InConnectivityListPtr);
	/*
	*	This is a constructor for making a FEVolume_c
	*	from an existing zone. The FEVolume_c gets all the
	*	information it needs from Tec360.
	*/
	FESurface_c(int ZoneNum,
		int VolZoneNum,
		vector<int> const & InXYZVarNums,
		vector<int> const & InIntVarNums,
		bool const CopyData = false);

	/*
	*	Make FEVolume from an existing zone.
	*	Can be IJ ordered zone or actual FE zone.
	*	If IJ ordered zone, need to convert to FE triangle zone representation
	*		with connectivity info and such.
	*	This assumes the structure of the surfaces made by Bondalyzer
	*/
	FESurface_c(int ZoneNum,
		vector<int> const & InXYZVarNums);

	/*
	 *	Make FEVolume from lists of nodes (points) and elements (lists of point indices).
	 */
	FESurface_c(vector<vec3> const & Nodes, vector<vector<int> > & Elements);


	~FESurface_c();

	/*
	 *	Operators
	 */

	FESurface_c & operator+=(FESurface_c const & rhs);
	FESurface_c operator+(FESurface_c const & rhs) const;

	/*
	*	Getters
	*/
	int GetNumSides() const { return m_GPList.size(); }
	int GetGPZoneNum(int i) const { REQUIRE(0 <= i && i < m_GPList.size()); return m_GPList[i].GetZoneNum(); }
	Boolean_t IsMade() const { return m_FEVolumeMade; }
	Boolean_t IntResultsReady() const { return m_IntegrationResultsReady; }
	vector<double> GetIntResults() const;
	int GetZoneNum() const { return m_ZoneNum; }
	vector<vector<double> > GetTriSphereIntValsByElem(vector<double> * SphereTriangleAreas = nullptr) const { return TriSphereIntValsByElem(SphereTriangleAreas); }
	vector<double> TriSphereElemSolidAngles(double * TotalAreaIn = nullptr, vec3 * Origin = nullptr) const;
	int GetNumElems() const {return m_ElemList.size(); }
	int GetNumNodes() const {return m_XYZList.size(); }

	vector<vec3> GetSphereIntersectionPath(vec3 const & SphereCenter, double const & SphereRadius);
	bool ProjectPointToSurface(vec3 const & OldPoint, vec3 & NewPoint, int & ProjectedElemIndex, bool & ProjectionIsInterior, int MaxBFSDepth = DefaultProjectPointToSurfaceMaxBFSDepth, bool StartWithBFS = true) const;

	/*
	*	Setters
	*/
	Boolean_t Setup(int InZoneNum,
		int VolZoneNum,
		vector<int> const & InXYZVarNums,
		vector<int> const & InIntVarNums,
		bool const CopyData = false);
	Boolean_t Setup(int const InZoneNum,
		vector<int> const & InXYZVarNums);

	void AddGP(GradPath_c const & GP){ m_GPList.push_back(GP); }

	Boolean_t DoIntegration(Boolean_t IntegrateVolume, 
		vector<vec> const & stuW, 
		vector<int> const & SplitPtNums = vector<int>(),
		vector<vec> const & stuW2 = vector<vec>());
// 	Boolean_t DoIntegration(int ResolutionScale, Boolean_t IntegrateVolume);
	Boolean_t DoIntegrationNew(int ResolutionScale, Boolean_t IntegrateVolume);
	void EdgeMidpointSubdivide(std::queue<int> & ElemsToDo, vector<vector<double> > & ElemVals, vec3 * CPPos = nullptr, double * SphereRadius = nullptr, vector<vector<int> > * OldElemsToNewElems = nullptr);
	void EdgeMidpointSubdivideParallel(vector<vector<double> > & ElemVals, vec3 * CPPos = nullptr, double * SphereRadius = nullptr, vector<vector<int> > * OldElemsToNewElems = nullptr);
	vector<vector<double> > CellCenteredToNodalVals(vector<vector<double> > & ElemVals, vec3 * SphereOrigin = nullptr);
	vector<vector<double> > NodalToCellCenteredVals(vector<vector<double> > & NodeVals, vec3 * SphereOrigin = nullptr);

	Boolean_t GQIntegration(int NumGQPts, vector<FieldDataPointer_c> const & InIntFDPtrs, Boolean_t IntegrateVolume);

	double IntVolume(int N, vec3 const & StartPoint) const;
	vector<mat> GetIntegrationPointsWeights(vec3 const & StartPoint, vector<vec> const & stuW) const;
	vector<mat> GetIntegrationPointsWeights(vector<vec> const & stuW) const;

	void GeneratePointElementDistanceCheckData();
	void GenerateElemMidpoints();
	void GetNodeConnectivityFromTecplot();
	void GenerateElemConnectivity(int numSharedNodes = 1);
	void GenerateNodeToElementList();
	vector<vector<int> > const * GetNodeToElementListPtr() const { return &m_NodeToElementList; }
	vector<vector<LgIndex_t> > const * GetNodeConnectivityListPtr() const { return &m_NodeConnectivityList; }
	vector<vector<LgIndex_t> > const * GetElemConnectivityListPtr() const { return &m_ElemConnectivityList; }
	vector<vector<int> > const * GetElemListPtr() const { return &m_ElemList; }
	vector<vec3> const * GetElemMidpointsPtr() const { return &m_ElemMidPoints; }
	vector<vec3> const * GetXYZListPtr() const { return &m_XYZList; }

// 	Boolean_t MakeGradientBundle(vector<GradPath_c*> GPs);
	Boolean_t MakeFromGPs(vector<GradPathBase_c const *> const & GPs, 
		bool ConnectBeginningAndEndGPs = false, 
		bool AddCapToSurface = false, 
		bool AddCapsBetweenGPs = true, 
		vector<vec3> const * CapPoints = nullptr, 
		int LeftPathDeviationPointInd = -1,
		int RightPathDeviationPointInd = -1,
		int CapMidPointInd = -1);
 	Boolean_t MakeFromGPs(vector<GradPath_c const *> const & GPs, 
		bool ConnectBeginningAndEndGPs = false, 
		bool AddCapToSurface = false, 
		bool AddCapsBetweenGPs = true, 
		vector<vec3> const * CapPoints = nullptr, 
		int LeftPathDeviationPointInd = -1,
		int RightPathDeviationPointInd = -1,
		int CapMidPointInd = -1)
 	{
 		vector<GradPathBase_c const *> TmpGPPtrs;
 		for (auto i : GPs)
 			TmpGPPtrs.push_back(reinterpret_cast<GradPathBase_c const *>(i));
 
 		return MakeFromGPs(TmpGPPtrs, ConnectBeginningAndEndGPs, AddCapToSurface, AddCapsBetweenGPs, CapPoints, LeftPathDeviationPointInd, RightPathDeviationPointInd, CapMidPointInd);
 	}

	Boolean_t MakeFromNodeElemList(vector<vec3> const & P, vector<vector<int> > const & T);

	Boolean_t Refine();

	int GetClosestNodeToPoint(vec3 const & Point, double * ClosestNodeDistance = nullptr) const;
	/*
	*	Two methods to make the 3- and 4-sided
	*	FE volumes from gradient paths.
	*/
// 	Boolean_t MakeGradientBundle(GradPath_c const & GP1,
// 		GradPath_c const & GP2,
// 		GradPath_c const & GP3);
// 	Boolean_t MakeGradientBundle(GradPath_c const & GP1,
// 		GradPath_c const & GP2,
// 		GradPath_c const & GP3,
// 		GradPath_c const & GP4);

	/*
	*	Other
	*/
	static int SetZoneStyle(int const ZoneNum,
		AssignOp_e const ZoneActive = AssignOp_PlusEquals,
		Boolean_t const ShowContour = TRUE,
		Boolean_t const ShowMesh = TRUE,
		Boolean_t const ShowScatter = FALSE,
		Boolean_t const ShowShade = FALSE);
	int SaveAsTriFEZone(vector<int> const & XYZVarNums, string ZoneName = "");
	int SaveAsTriFEZone(string const & ZoneName, 
		vector<FieldDataType_e> DataTypes,
		vector<ValueLocation_e> const & DataLocations,
		vector<int> const & XYZVarNums,
		int RhoVarNum = 4);
// 	int const SaveAsFEZone(
// 		vector<FieldDataType_e> DataTypes,
// 		vector<int> const & XYZVarNums,
// 		int RhoVarNum
// 		);

	friend class Domain_c;

private:

	double FESurface_c::PointDistanceToElementSquared(vec3 const & P, int e, vec3 & ClosestPoint)  const;
	Boolean_t PointIsInterior(vec3 const & Point, vector<vec3> const & FarPoints) const;
	int TriangleIntersect(vec3 const & T_P0,
		vec3 const & T_P1,
		vec3 const & T_P2,
		vec3 const & R_P1,
		vec3 const & R_P0) const;
	vector<vector<double> > TriSphereIntValsByElem(vector<double> * SphereTriangleAreas = nullptr) const;
	Boolean_t CalcMaxNodeDistSqr();
	Boolean_t DistSqrToSurfaceNodeWithinTolerance(vec3 const & CheckPt,
		double & NewDistSqrUnderTol,
		int & CloseNodeNum,
		double const & DistSqrTol = -1.0);
	Boolean_t SubDivideIntegrateCellAtPoint(vec3 const & Point,
		vector<vec3> const & FarPoints,
		vec3 const & DelXYZ,
		double const & MinDistSqrToSurfaceNode,
		int MinDistNodeNum,
		int SubDivideLevel,
		Boolean_t IntegrateVolume);
// 	vector<int> TriangleEdgeMidPointSubdivide(int TriNum);
	void RefineTriElems(vector<int> const & TriNumList);
// 	void TriPolyLines(bool const ConnectBeginningAndEndGPs = true);
	void RemoveDupicateNodes();


	vector<GradPath_c> m_GPList;

	vector<vec3> m_XYZList;

	// data for speeding up surface grad paths
	vector<vec3> m_ElemMidPoints, m_v21, m_v32, m_v13, m_normals,
		m_v21crossN, m_v32crossN, m_v13crossN;
	vector<double> m_oneOverMagSqrV21, m_oneOverMagSqrV32, m_oneOverMagSqrV13, m_oneOverMagSqrN;

	vector<vec3> m_edge0, m_edge1;
	vector<double> m_a, m_b, m_c, m_det, m_invDet, m_oneOverDenom, m_oneOverA, m_oneOverC;
	// end


	vector<double> m_RhoList;

	NodeMap_t* m_ConnectivityListPtr;
	vector<vector<LgIndex_t> > m_NodeConnectivityList,
		m_ElemConnectivityList;
	vector<double> m_MaxNeighborNodeDistSqr,
		m_MinNeighborNodeDistSqr;
	vector<vec3> m_RefinedXYZList;
	vector<vector<int> > m_ElemList,
		m_NodeToElementList;

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

// 	int m_NumGPs;
// 	int m_NumGPPts;
	int m_NumNodes;
	int m_NumElems;

	int m_NumNodesPerElem;

	vec3 m_ZoneMinXYZ;;
	vec3 m_ZoneMaxXYZ;

	VolExtentIndexWeights_s m_VolZoneInfo;

	Boolean_t m_FEVolumeMade;
	Boolean_t m_IntegrationResultsReady;
};


vector<vec> GetWeightsPoints(int N);


class Domain_c
{
public:
	Domain_c(){}
	Domain_c(vector<int> const & V, FESurface_c *Vol){ Setup(V, Vol); }
	~Domain_c(){ m_Vol = nullptr; }
	void Setup(vector<int> const & V, FESurface_c *Vol);

	double Weight() const;
protected:
	
private:
	double Split(int ei, vector<int> & t, Domain_c & D1, Domain_c & D2) const;
	double WeightFunc(vector<int> const & t) const;


	vector<int> m_V;
	vector<vector<int> > m_E;
	FESurface_c *m_Vol = nullptr;
};

void ResizeSphere(int ZoneNum, double const & SizeFactor, Boolean_t AbsoluteRadius, bool ScaleByVar = false, int ScaleVarNum = -1, double ScaleFactor = 1.0, bool LogScale = true);


void GetClosedIsoSurface(int IsoZoneNum, const std::vector<FieldDataPointer_c> & IsoReadPtrs, std::vector<int> & NodeNums, AddOn_pa AddOnID);