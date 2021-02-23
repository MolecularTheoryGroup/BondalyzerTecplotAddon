#include "TECADDON.h"
#if !defined (MSWIN)
#include <unistd.h>
#endif

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include <cmath>

#include <armadillo>

#include "ArgList.h"
#include "Set.h"

#include <omp.h>

#include "CSM_GRAD_PATH.h"
#include "CSM_DATA_TYPES.h"
#include "CSM_DATA_SET_INFO.h"
#include "CSM_FE_VOLUME.h"
#include "CSM_GEOMETRY.h"

#include "updateSphericalTriangulation.h"

// #define RECORD_INT_POINTS


using std::vector;
using std::string;
using std::to_string;
using std::stringstream;
using std::cout;
using std::endl;

using namespace arma;


/*
*	Public member functions
*/

/*
 *	Operators
 */

FESurface_c & FESurface_c::operator+=(FESurface_c const & rhs){
	int NodeCount = m_XYZList.size(),
		ElemCount = m_ElemList.size();

	// Copy nodes and elements from rhs
	m_XYZList.insert(m_XYZList.end(), rhs.m_XYZList.cbegin(), rhs.m_XYZList.cend());
	m_ElemList.insert(m_ElemList.end(), rhs.m_ElemList.cbegin(), rhs.m_ElemList.cend());

	// Fix new element numbers
	for (int e = ElemCount; e < m_ElemList.size(); ++e){
		for (int & i : m_ElemList[e]){
			i += NodeCount;
		}
	}

	m_NumElems = m_ElemList.size();
	m_NumNodes = m_XYZList.size();

	m_FEVolumeMade = (m_NumNodes > 0 && m_NumElems > 0);

	return *this;
}

FESurface_c FESurface_c::operator+(FESurface_c const & rhs) const{
	return FESurface_c(*this) += rhs;
}

/*
*	Constructor/destructor
*/

FESurface_c::FESurface_c()
{
// 	m_NumGPs = -1;
	m_FEVolumeMade = FALSE;
	m_IntegrationResultsReady = FALSE;
}

FESurface_c::FESurface_c(int InZoneNum,
	ZoneType_e const & InZoneType,
	vector<FieldDataPointer_c> const & InXYZPtrs,
	vec3 const & InMaxXYZ,
	vec3 const & InMinXYZ,
	vector<int> const InMaxIJK,
	NodeMap_t* InConnectivityListPtr)
{
	m_FEVolumeMade = TRUE;
	m_IntegrationResultsReady = FALSE;

	m_XYZPtrs = InXYZPtrs;

	m_FEVolumeMade = m_XYZPtrs.size() == 3;

	for (int i = 0; i < 3 && m_FEVolumeMade; ++i){
		m_FEVolumeMade = m_XYZPtrs[i].IsReady();
	}

	if (m_FEVolumeMade)
		m_FEVolumeMade = (InMaxIJK.size() == 3);

	if (m_FEVolumeMade)
		m_FEVolumeMade = (InConnectivityListPtr != nullptr);

	if (m_FEVolumeMade){
		TecUtilZoneGetIJK(InZoneNum, &m_NumNodes, &m_NumElems, &m_NumNodesPerElem);
		m_FEVolumeMade = ((InZoneType == ZoneType_FETriangle && m_NumNodesPerElem == 3)
			|| (InZoneType == ZoneType_FEQuad && m_NumNodesPerElem == 4));
		if (m_FEVolumeMade){
			m_VolZoneInfo.MaxIJK = InMaxIJK;
			m_VolZoneInfo.MaxXYZ = InMaxXYZ;
			m_VolZoneInfo.MinXYZ = InMinXYZ;
			m_ConnectivityListPtr = InConnectivityListPtr;
		}
	}
}

/*
*	Make FEVolume from lists of nodes (points) and elements (lists of point indices).
*/
FESurface_c::FESurface_c(vector<vec3> const & Nodes, vector<vector<int> > & Elements)
{
	m_FEVolumeMade = (Nodes.size() >= 3 && Elements.size() >= 1);
	if (m_FEVolumeMade){
		m_XYZList.clear();
		m_XYZList.insert(m_XYZList.begin(), Nodes.cbegin(), Nodes.cend());
		m_ElemList.insert(m_ElemList.begin(), Elements.cbegin(), Elements.cend());
		m_NumElems = m_ElemList.size();
		m_NumNodes = m_XYZList.size();
	}
}

FESurface_c::FESurface_c(int InZoneNum,
	int VolZoneNum,
	vector<int> const & InXYZVarNums,
	vector<int> const & InIntVarNums,
	bool const CopyData)
{
	Setup(InZoneNum, VolZoneNum, InXYZVarNums, InIntVarNums, CopyData);
}

FESurface_c::FESurface_c(int InZoneNum,
	vector<int> const & InXYZVarNums)
{
	Setup(InZoneNum, InXYZVarNums);
}


/*
	Make FEVolume from vector of pointers to GradPath_c's
*/
// Boolean_t FESurface_c::MakeGradientBundle(vector<GradPath_c*> GPs)
// {
// 	Boolean_t IsOk = TRUE;
// 	for (int i = 0; i < GPs.size() && IsOk; ++i){
// 		IsOk = GPs[i]->IsMade();
// 		if (IsOk && i > 0){
// 			IsOk = (GPs[i]->GetCount() == GPs[i-1]->GetCount());
// 		}
// 	}
// 	
// 	if (!IsOk){
// 		int MinCount = INT_MAX;
// 		Boolean_t AllMade = TRUE;
// 		for (auto *g : GPs){
// 			MinCount = MIN(MinCount, g->GetCount());
// 			AllMade = AllMade && g->IsMade();
// 		}
// 		if (AllMade){
// 			for (auto *g : GPs){
// 				g->Resample(MinCount);
// 			}
// 			IsOk = TRUE;
// 		}
// 	}
// 
// 	if (IsOk){
// 		m_NumGPs = GPs.size();
// 		m_NumGPPts = GPs[0]->GetCount();
// 		m_NumNodes = m_NumGPs * m_NumGPPts;
// 		m_NumElems = 2 * (m_NumGPs - 2) + 2 * m_NumGPs * (m_NumGPPts - 1);
// 
// 		vector<vector<vec3>::const_iterator> XYZIt(m_NumGPs);
// 		
// 		for (int i = 0; i < GPs.size(); ++i){
// 			XYZIt[i] = GPs[i]->m_XYZList.cbegin();
// 		}
// 
// 		m_XYZList.resize(m_NumGPs * m_NumGPPts);
// 
// 		for (int i = 0; i < m_NumGPs; ++i){
// 			for (int j = 0; j < m_NumGPPts; ++j){
// 				m_XYZList[i * m_NumGPPts + j] = *XYZIt[i];
// 				XYZIt[i]++;
// 			}
// 		}
// 
// 		/*
// 		 * Generate triangle list.
// 		 */
// 
// // 		vector<int> V(m_NumGPs);
// // 		for (int i = 0; i < m_NumGPs; ++i){
// // 			V[i] = i * m_NumGPPts;
// // 		}
// //  
// // 		Domain_c D(V, this);
// // 		D.Weight();
//  
// 		TriPolyLines();
//  
// // 		V = vector<int>();
// // 		V.reserve(m_NumGPs * 2);
// // 		int OldNumNodes = m_XYZList.size() - m_NumGPs;
// // 		for (int i = 0; i < m_NumGPs; ++i){
// // 			V.push_back((i + 1) * m_NumGPPts - 1);
// // 			V.push_back(OldNumNodes + i);
// // 		}
// // 		D.Setup(V, this);
// // 		D.Weight();
// 
// 		m_NumElems = m_ElemList.size();
// 
// 		m_FEVolumeMade = TRUE;
// 	}
// 
// 	return IsOk;
// }

/*
Make FEVolume from vector of pointers to GradPath_c's
*/
Boolean_t FESurface_c::MakeFromGPs(vector<GradPathBase_c const *> const & GPs, 
	bool ConnectBeginningAndEndGPs, 
	bool AddCap, 
	bool AddCapsBetweenGPs, 
	vector<vec3> const * CapPoints, 
	int LeftPathDeviationPointInd,
	int RightPathDeviationPointInd, 
	int CapMidPointInd)
{
	Boolean_t IsOk = TRUE;
	AddCapsBetweenGPs = false; // depricated; looks better the other way, and saving GBs as zones is purely for visualization.
// 	AddCap = true; // depricated; now we add "caps" along the entire GB for visual effect of being closed even when truncated.

	for (int i = 0; i < GPs.size() && IsOk; ++i){
	 	IsOk = GPs[i]->IsMade() && GPs[i]->GetCount() > 0;
	}



	
	if (IsOk){
		
		/*
		*	Get total number of points in all GPs
		*/
		m_NumNodes = 0;
		for (auto const & i : GPs)
			m_NumNodes += i->GetCount();

		/*
		 *	Save GP points to m_XYZList and save each GP's points'
		 *	indices to a 2d int vector
		 */
		m_XYZList.clear();
		m_XYZList.resize(m_NumNodes);
		m_RhoList.clear();
		vector<vector<int> > IndList;
		int ni = 0;
		for (auto * g : GPs){
			IndList.push_back(vector<int>(g->GetCount()));
			int gi = 0;
			for (auto const & p : g->m_XYZList){
				m_XYZList[ni] = p;
				IndList.back()[gi++] = ni++;
			}
			m_RhoList.insert(m_RhoList.end(), g->m_RhoList.cbegin(), g->m_RhoList.cend());
		}

		/*
		 * If connecting beginning and end paths, need to copy the first indList to the end of indList
		 */
		if (ConnectBeginningAndEndGPs) IndList.push_back(IndList[0]);
		/*
			*	Now need to generate the "cap" points that will be used to linearly
			*	fill in any excess space between the endpoints of neighboring GPs.
			*/

		/*
			* First, get the average spacing between GP points to come up with a decent
			*	spacing for the cap points
			*/
		double AvgCapSpacing = 0;
		if (AddCapsBetweenGPs) {
			int Denom = 0;
			for (auto * g : GPs) {
				for (int i = 0; i < g->GetCount() - 1; ++i) {
					AvgCapSpacing += Distance(g->XYZAt(i), g->XYZAt(i + 1));
					Denom++;
				}
			}
			AvgCapSpacing /= double(Denom);
		}

		/*
			*	Now check to see if each neighboring pair of paths needs a cap.
			*	If a cap is necessary (i.e. the end points of the paths are spaced
			*	farther apart than CapSpacing) then use the CapSpacing as a starting
			*	point to find the number of equidistant cap points, saving the points
			*	themselves to m_XYZList and the corresponding incides of the cap points
			*	to CapIndList.
			*/
		vector<vector<int> > CapIndList(IndList.size() - 1);
		// size - 1 because there are IndList.size() - 1 pairs

		if (AddCapsBetweenGPs) {
			for (int i = 0; i < IndList.size() - 1; ++i) {
				double CapDist = Distance(m_XYZList[IndList[i].back()], m_XYZList[IndList[i + 1].back()]);
				if (CapDist > AvgCapSpacing) {
					/*
						*	Cap points needed, so find NumCapPoints necessary
						*	and then make equidistant points.
						*/
					int NumCapPoints = CapDist / AvgCapSpacing;
					/*
						*	Always want an odd number of cap points to ensure that one
						*	is at the midpoint between the path end points.
						*/
					if (NumCapPoints > 0 && NumCapPoints % 2 == 1) NumCapPoints--;

					/*
						*	Now make the cap points, starting from the "left" path,
						*	saving the index of each new point to CapIndList.
						*/
					vec3 CapVec = (m_XYZList[IndList[i + 1].back()] - m_XYZList[IndList[i].back()]) / double(NumCapPoints);
					for (int j = 0; j < NumCapPoints; ++j) {
						m_XYZList.push_back(m_XYZList[IndList[i].back()] + CapVec * double(j));
						CapIndList[i].push_back(ni++);
					}
				}
			}
		}

		/*
		 * Now all the GP points, cap points, and their respective indices have been
		 * generated.
		 * Next call the path stitching function on each pair.
		 */
		m_ElemList.clear();
		if (CapPoints != nullptr){
			// The presence of a list of cap points means that this is a special case for
			// interatomic/ring surface generation where two surfaces going into the center CP
			// need to be stitched with a bond/ring path.
			// Assume that IndList is of length 2, so that CapIndList is of size 1.
			REQUIRE(IndList.size() == 2);
			for (auto const & v : *CapPoints){
				CapIndList.back().push_back(m_XYZList.size());
				m_XYZList.push_back(v);
			}
			AddCapsBetweenGPs = true;
		}

		for (int i = 0; i < IndList.size() - 1; ++i) {
			if (AddCapsBetweenGPs)
				StitchCapPaths(IndList[i], IndList[i + 1], CapIndList[i], m_XYZList, m_ElemList, LeftPathDeviationPointInd, RightPathDeviationPointInd, CapMidPointInd);
			else
				StitchPaths(IndList[i], IndList[i + 1], m_XYZList, m_ElemList);
		}



		if (AddCap) {
			int NumEdgePointsNormal = (GPs.size() - 3) / 3,
				NumEdgePointsSaddle = (GPs.size() - 4) / 4;
			int TotNumPointsNormal = 3 * NumEdgePointsNormal + 3,
				TotNumPointsSaddle = 4 * NumEdgePointsSaddle + 4;

// 			if (TotNumPointsSaddle > 4 || TotNumPointsNormal > 3) {
			if (GPs.size() > 4) {
				/*
				 *	GB edges have been connected into surfaces.
				 *	Now need to create surfaces for the initial and terminal 3d caps
				 *	to close off the GB.
				 *	The initial 2d cap is guaranteed to be a triangle since it coincides
				 *	with one of the triangular elements on the sphere.
				 *	The terminal 2d cap can be stitched by finding the pair of GP end
				 *	points of greatest separation and dividing all end points into
				 *	two paths and then stitching those.
				 */

				 /*
				 *	There are two kinds of GBs:
				 *	One for "normal" GBs and another for GBs that have a bond
				 *	path as an edge.
				 *
				 *	The normal GBs have a triangle of points at their base with
				 *	e edge points, so 3e+3 total points.
				 *
				 *	Bond path GBs are similar, but have multiple points at one of
				 *	the triangular vertices because of the e+2 paths on the far
				 *	side of the bond path (that is, along the interatomic surface).
				 *	So they have 4e+4 total points. The first 2e+2 and the last e+1 GPs of a bond path
				 *	GB (and the (2e+3)-th GP) describe the full set of points along
				 *	the triangular element, and the set [2e+3,3e+3] all originate
				 *	from the same point as the (3e+4)-th GP.
				 */


				 /*
				  *	If TotNumPoints == GPs.size() then it's a normal GB, if not it's
				  *	a bond path GB.
				  */
				  // 		vector<int> InitialCapInd;

// 				if (TotNumPointsSaddle == GPs.size()) {
// 					// 			for (int i = 0; i < 3; ++i) InitialCapInd.push_back(IndList[(NumEdgePointsNormal + 1) * i][0]);
// 					// 			InitialCapInd.back()++;
// 					m_ElemList.push_back({ IndList[0][0],
// 						IndList[NumEdgePointsSaddle + 1][0],
// 						IndList[(NumEdgePointsSaddle) * 2 + 2][0] });
// 				}
// 				else if (TotNumPointsNormal == GPs.size()) {
// 					// 			for (int i = 0; i < 3; ++i) InitialCapInd.push_back(IndList[(NumEdgePointsNormal + 1) * i][0]);
// 					m_ElemList.push_back({ IndList[0][0],
// 						IndList[NumEdgePointsNormal + 1][0],
// 						IndList[(NumEdgePointsNormal + 1) * 2][0] });
// 				}
// 				else
// 					TecUtilDialogMessageBox("Incorrect number of GPs for GB", MessageBoxType_Error);
				// 		m_ElemList.push_back(InitialCapInd);


				/*
				 *	Find the GP end points of maximum separation at each step down the GPs
				 *	which are assumed to have equispaced points, although the results wouldn't
				 *	be terrible even with nonuniform spacing.
				 *	We can neglect the presence of 1d cap points in CapIndList because
				 *	they are guaranteed to be straight lines, so by including the
				 *	GP endpoints we'll necessarily end up with triangular elements
				 *	in the resulting 2d cap that coincide with the existing 1d cap points.
				 */


// 				for (int gpi = 10; gpi < IndList[0].size(); ++gpi) {
					double MaxDistSqr = -1;
					int MaxDistsGPNums[2] = { -1, -1 };


					for (int i = 0; i < IndList.size() - (ConnectBeginningAndEndGPs ? 2 : 1); ++i) {

						for (int j = i + 1; j < GPs.size(); ++j) {

							double TmpDistSqr = DistSqr(GPs[i]->XYZAt(-1), GPs[j]->XYZAt(-1)); // old version that only capped end of GBs
// 							double TmpDistSqr = DistSqr(GPs[i]->XYZAt(MIN(gpi, GPs[i]->m_NumGPPoints - 1)), GPs[j]->XYZAt(MIN(gpi,GPs[j]->m_NumGPPoints-1)));
							if (TmpDistSqr > MaxDistSqr) {
								MaxDistSqr = TmpDistSqr;
								MaxDistsGPNums[0] = i;
								MaxDistsGPNums[1] = j;
							}
						}
					}


					if (MaxDistSqr > 0.0001) {
						vector<int> LVec, RVec;
						for (int i = MaxDistsGPNums[0]; i < MaxDistsGPNums[1]; ++i) {
							LVec.push_back(IndList[i].back());
// 							LVec.push_back(IndList[i][MIN(gpi, IndList[i].size() - 1)]);
						}
						for (int i = 0; i < GPs.size() - (MaxDistsGPNums[1] - MaxDistsGPNums[0]); ++i) {
							RVec.push_back(IndList[(MaxDistsGPNums[1] + i) % GPs.size()].back());
// 							RVec.push_back(IndList[(MaxDistsGPNums[1] + i) % GPs.size()][MIN(gpi, IndList[(MaxDistsGPNums[1] + i) % GPs.size()].size() - 1)]);
						}
						std::reverse(RVec.begin(), RVec.end());

						/*
						 *	Now stitch these two "paths" and add to the element list
						 */
						StitchPaths(LVec, RVec, m_XYZList, m_ElemList);
					}
// 				}
// 				else
// 					IsOk = FALSE;
			}
			else{
				// No edge GPs are present, so the cap connectivity is either the final points of the GPs,
				// if 3 total GPs, or splitting a quad into two triangles, if 4 total GPs.
				// 
				if (GPs.size() == 3){
					m_ElemList.push_back({ IndList[0][0], IndList[1][0], IndList[2][0] });
					m_ElemList.push_back({ IndList[0].back(), IndList[1].back(), IndList[2].back() });
				}
				else if (GPs.size() == 4) {
					// Check the two ways to split the quads initial and terminal caps and split
					// according to the shortest split edge.
					
					vector<vector<int> > TmpIndList = {
						{ IndList[0][0], IndList[1][0], IndList[2][0], IndList[3][0] },
						{ IndList[0].back(), IndList[1].back(), IndList[2].back(), IndList[3].back() } },
						Split02 = { {0,1,2},{0,2,3} },
						Split13 = { {0,1,3},{1,2,3} };

					for (auto const & inds : TmpIndList){
						if (DistSqr(m_XYZList[inds[0]], m_XYZList[inds[2]]) < DistSqr(m_XYZList[inds[1]], m_XYZList[inds[3]])){
							for (auto const & tri : Split02)
								m_ElemList.push_back({ inds[tri[0]], inds[tri[1]], inds[tri[2]] });
						}
						else {
							for (auto const & tri : Split13)
								m_ElemList.push_back({ inds[tri[0]], inds[tri[1]], inds[tri[2]] });
						}
					}
				}
				else{
					IsOk = FALSE;
				}
			}
		}
	}

	if (IsOk){
		m_NumNodes = m_XYZList.size();
		m_NumElems = m_ElemList.size();

		m_FEVolumeMade = TRUE;
	}

// 	RemoveDupicateNodes();

	return IsOk;


	/*
	 *	Old implementation that resampled paths to ensure they all had the same number of points,
	 *	which was stupid
	 */
// 	for (int i = 0; i < GPs.size() && IsOk; ++i){
// 		IsOk = GPs[i]->IsMade();
// 		if (IsOk && i > 0){
// 			IsOk = (GPs[i]->GetCount() == GPs[i - 1]->GetCount());
// 		}
// 	}
// 
// 	if (!IsOk){
// 		int MinCount = INT_MAX;
// 		Boolean_t AllMade = TRUE;
// 		for (auto *g : GPs){
// 			MinCount = MIN(MinCount, g->GetCount());
// 			AllMade = AllMade && g->IsMade();
// 		}
// 		if (AllMade){
// 			for (auto *g : GPs){
// 				g->Resample(MinCount);
// 			}
// 			IsOk = TRUE;
// 		}
// 	}
// 
// 	if (IsOk){
// 		m_NumGPs = GPs.size();
// 		m_NumGPPts = GPs[0]->GetCount();
// 		m_NumNodes = m_NumGPs * m_NumGPPts;
// 		m_NumElems = 2 * (m_NumGPs - 2) + 2 * m_NumGPs * (m_NumGPPts - 1);
// 
// 		vector<vector<vec3>::const_iterator> XYZIt(m_NumGPs);
// 
// 		m_XYZList.clear();
// 		for (GradPath_c* GP : GPs)
// 			m_XYZList.insert(m_XYZList.begin(), GP->m_XYZList.cbegin(), GP->m_XYZList.cend());
// 
// // 		for (int i = 0; i < GPs.size(); ++i){
// // 			XYZIt[i] = GPs[i]->m_XYZList.cbegin();
// // 		}
// // 
// // 		m_XYZList.resize(m_NumGPs * m_NumGPPts);
// // 
// // 		for (int i = 0; i < m_NumGPs; ++i){
// // 			for (int j = 0; j < m_NumGPPts; ++j){
// // 				m_XYZList[i * m_NumGPPts + j] = *XYZIt[i];
// // 				XYZIt[i]++;
// // 			}
// // 		}
// 
// 		/*
// 		* Generate triangle list.
// 		*/
// 
// 		TriPolyLines(ConnectBeginningAndEndGPs);
// 
// 		m_NumElems = m_ElemList.size();
// 
// 		m_FEVolumeMade = TRUE;
// 	}
// 
// 	return IsOk;
}

Boolean_t FESurface_c::MakeFromNodeElemList(vector<vec3> const & P, vector<vector<int> > const & T)
{
	Boolean_t IsOk = P.size() >= 3 && T.size() >= 1;

	if (IsOk){
		m_XYZList = P;
		m_ElemList = T;

		m_NumElems = m_ElemList.size();
		m_NumNodes = m_XYZList.size();
		m_FEVolumeMade = TRUE;
	}

	return IsOk;
}

/*
*	Two constructors for the 3- and 4-sided
*	FE volumes.
*/
// Boolean_t FESurface_c::MakeGradientBundle(GradPath_c const & GP1,
// 	GradPath_c const & GP2,
// 	GradPath_c const & GP3)
// {
// 	Boolean_t IsOk = (GP1.IsMade()
// 		&& GP2.IsMade()
// 		&& GP3.IsMade()
// 		&& GP1.GetCount() == GP2.GetCount()
// 		&& GP1.GetCount() == GP3.GetCount());
// 
// 
// 	if (IsOk){
// 		m_GPList = { GP1, GP2, GP3 };
// 
// 		m_NumGPs = 3;
// 		m_NumGPPts = GP1.GetCount();
// 		m_NumNodes = 3 * m_NumGPPts;
// 		m_NumElems = 3 * (m_NumGPPts - 1) + 2;
// 
// 		vector<vec3>::const_iterator XYZIt[3];
// 		vector<double>::const_iterator RhoIt[3];
// 
// 		int GPNum = 0;
// 		XYZIt[GPNum] = GP1.m_XYZList.cbegin();
// 
// 		RhoIt[GPNum] = GP1.m_RhoList.cbegin();
// 		GPNum++;
// 
// 
// 		XYZIt[GPNum] = GP2.m_XYZList.cbegin();
// 
// 		RhoIt[GPNum] = GP2.m_RhoList.cbegin();
// 		GPNum++;
// 
// 
// 		XYZIt[GPNum] = GP3.m_XYZList.cbegin();
// 
// 		RhoIt[GPNum] = GP3.m_RhoList.cbegin();
// 
// 		m_XYZList.resize(3 * m_NumGPPts);
// 
// 		for (int j = 0; j < m_NumGPPts; ++j){
// 			for (int i = 0; i < 3; ++i){
// 				m_XYZList[i + 3 * j] = *XYZIt[i];
// 				XYZIt[i]++;
// 			}
// 		}
// 
// 		m_RhoList.resize(3 * m_NumGPPts);
// 
// 		for (int j = 0; j < m_NumGPPts; ++j){
// 			for (int i = 0; i < 3; ++i){
// 				m_RhoList[i + 3 * j] = *RhoIt[i];
// 				RhoIt[i]++;
// 			}
// 		}
// 
// 		m_FEVolumeMade = TRUE;
// 	}
// 
// 	return IsOk;
// }
// Boolean_t FESurface_c::MakeGradientBundle(GradPath_c const & GP1,
// 	GradPath_c const & GP2,
// 	GradPath_c const & GP3,
// 	GradPath_c const & GP4)
// {
// 	Boolean_t IsOk = (GP1.IsMade()
// 		&& GP2.IsMade()
// 		&& GP3.IsMade()
// 		&& GP4.IsMade()
// 		&& GP1.GetCount() == GP2.GetCount()
// 		&& GP1.GetCount() == GP3.GetCount()
// 		&& GP1.GetCount() == GP4.GetCount());
// 
// 
// 	if (IsOk){
// 		m_GPList = { GP1, GP2, GP3, GP4 };
// 
// 		m_NumGPs = 4;
// 		m_NumGPPts = GP1.GetCount();
// 		m_NumNodes = 4 * m_NumGPPts;
// 		m_NumElems = 4 * (m_NumGPPts - 1) + 2;
// 
// 		vector<vec3>::const_iterator XYZIt[4];
// 		vector<double>::const_iterator RhoIt[4];
// 
// 		int GPNum = 0;
// 		XYZIt[GPNum] = GP1.m_XYZList.cbegin();
// 
// 		RhoIt[GPNum] = GP1.m_RhoList.cbegin();
// 		GPNum++;
// 
// 		XYZIt[GPNum] = GP2.m_XYZList.cbegin();
// 
// 		RhoIt[GPNum] = GP2.m_RhoList.cbegin();
// 		GPNum++;
// 
// 		XYZIt[GPNum] = GP3.m_XYZList.cbegin();
// 
// 		RhoIt[GPNum] = GP3.m_RhoList.cbegin();
// 		GPNum++;
// 
// 		XYZIt[GPNum] = GP4.m_XYZList.cbegin();
// 
// 		RhoIt[GPNum] = GP4.m_RhoList.cbegin();
// 
// 		m_XYZList.resize(4 * m_NumGPPts);
// 
// 		for (int j = 0; j < m_NumGPPts; ++j){
// 			for (int i = 0; i < 4; ++i){
// 				m_XYZList[i + 4 * j] = *XYZIt[i];
// 				XYZIt[i]++;
// 			}
// 		}
// 
// 		m_RhoList.resize(4 * m_NumGPPts);
// 
// 		for (int j = 0; j < m_NumGPPts; ++j){
// 			for (int i = 0; i < 4; ++i){
// 				m_RhoList[i + 3 * j] = *RhoIt[i];
// 				RhoIt[i]++;
// 			}
// 		}
// 
// 		m_FEVolumeMade = TRUE;
// 	}
// 
// 	return IsOk;
// }

FESurface_c::~FESurface_c()
{
}

/*
*	Getter methods
*/

vector<double> FESurface_c::GetIntResults() const
{
	if (m_IntegrationResultsReady){
		return m_IntValues;
	}
	else
		return vector<double>(1, 0);
}




double FESurface_c::PointDistanceToElementSquared(vec3 const & P, int elem, vec3 & ClosestPoint) const{
	/*
	 * from https://www.gamedev.net/forums/topic/552906-closest-point-on-triangle/
	 */
	vec3 v0 = m_XYZList[m_ElemList[elem][0]] - P;

	double d = dot(m_edge0[elem], v0),
		e = dot(m_edge1[elem], v0),
		s = m_b[elem] * e - m_c[elem] * d,
		t = m_b[elem] * d - m_a[elem] * e;

	if (s + t < m_det[elem])
	{
		if (s < 0.0)
		{
			if (t < 0.0)
			{
				if (d < 0.0)
				{
					s = CLAMP(-d * m_oneOverA[elem], 0.0, 1.0);
					t = 0.0;
				}
				else
				{
					s = 0.0;
					t = CLAMP(-e * m_oneOverC[elem], 0.0, 1.0);
				}
			}
			else
			{
				s = 0.0;
				t = CLAMP(-e * m_oneOverC[elem], 0.0, 1.0);
			}
		}
		else if (t < 0.0)
		{
			s = CLAMP(-d * m_oneOverA[elem], 0.0, 1.0);
			t = 0.0;
		}
		else
		{
			s *= m_invDet[elem];
			t *= m_invDet[elem];
		}
	}
	else
	{
		if (s < 0.0)
		{
			double tmp0 = m_b[elem] + d;
			double tmp1 = m_c[elem] + e;
			if (tmp1 > tmp0)
			{
				double numer = tmp1 - tmp0;
				s = CLAMP(numer * m_oneOverDenom[elem], 0.0, 1.0);
				t = 1 - s;
			}
			else
			{
				t = CLAMP(-e * m_oneOverC[elem], 0.0, 1.0);
				s = 0.0;
			}
		}
		else if (t < 0.0)
		{
			if (m_a[elem] + d > m_b[elem] + e)
			{
				double numer = m_c[elem] + e - m_b[elem] - d;
				s = CLAMP(numer * m_oneOverDenom[elem], 0.0, 1.0);
				t = 1 - s;
			}
			else
			{
				s = CLAMP(-e * m_oneOverC[elem], 0.0, 1.0);
				t = 0.0;
			}
		}
		else
		{
			double numer = m_c[elem] + e - m_b[elem] - d;
			s = CLAMP(numer * m_oneOverDenom[elem], 0.0, 1.0);
			t = 1.0 - s;
		}
	}

	ClosestPoint = m_XYZList[m_ElemList[elem][0]] + s * m_edge0[elem] + t * m_edge1[elem];
	return DistSqr(ClosestPoint, P);
}

bool FESurface_c::ProjectPointToSurface(vec3 const & OldPoint, vec3 & NewPoint, int & ProjectedElemIndex, bool & ProjectionIsInterior, int MaxBFSDepth, bool StartWithBFS) const{
	bool IsOk = m_XYZList.size() > 0
		&& m_ElemList.size() > 0;
	if (!IsOk) return IsOk;

	// find minimum distance element and try to project to that,
	// moving to a breadth first search starting from ProjectedElemIndex
	// to find an intersecting element if the closest element isn't a hit
	// or if ProjectedElemIndex >= 0.
	// 
	// 

	/*
	 *	quick point to triangle check for every element to get the minimum
	 *	from http://www.iquilezles.org/www/articles/triangledistance/triangledistance.htm
	 */


	 // 	if (ProjectedElemIndex < 0 || ProjectedElemIndex >= m_ElemList.size()){

// 	// find minimum distance element and try to project to that,
// 	// moving to a breadth first search starting from ProjectedElemIndex.

	double MinDistSqr = DBL_MAX,
		OldMinDistSqr = DBL_MAX,
		TmpDistSqr;
	int NumConvergedLevels = 0;
 	if (ProjectedElemIndex >= 0 && StartWithBFS){
 		// BFS with depth limit
 		std::queue<int> ToCheck;
 		vector<bool> Visited(m_ElemList.size(), false);
 		int CurDepth = 0;
 		ToCheck.push(ProjectedElemIndex);
 		ToCheck.push(-1); //Depth level marker
		Visited[ProjectedElemIndex] = true;
 		while (!ToCheck.empty()){
 			if (ToCheck.front() == -1){
 				// Moving down another depth
 				ToCheck.push(-1);
 				if (++CurDepth > MaxBFSDepth)
 					break;
				
				if (MinDistSqr < OldMinDistSqr) {
					NumConvergedLevels = 0;
					OldMinDistSqr = MinDistSqr;
				}
				else
					NumConvergedLevels++;

				if (NumConvergedLevels >= DefaultProjectPointtoSurfaceNumConvergedBSFDepthLevels)
					break;
 			}
 			else {
 				int e = ToCheck.front();
				TmpDistSqr = PointDistanceToElementSquared(OldPoint, e, NewPoint);
				// #pragma omp critical
				if (TmpDistSqr < MinDistSqr) {
					MinDistSqr = TmpDistSqr;
					ProjectedElemIndex = e;
				}
 				for (int n : m_ElemConnectivityList[e]) {
 					if (!Visited[n]) {
 						ToCheck.push(n);
						Visited[n] = true;
 					}
 				}
 			}
 			ToCheck.pop();
 		}
 	}

	if (MinDistSqr > 1e-12 || !StartWithBFS) {
		vector<std::pair<double, int> > MinDistSqrIndVector(omp_get_num_procs(), std::make_pair(MinDistSqr, ProjectedElemIndex));
		vector<vec3> JunkVecVector(omp_get_num_procs());
		int numElems = m_edge0.size();
#pragma omp parallel for
		for (int e = 0; e < numElems; ++e) {
			int numThreads = omp_get_num_threads();
			int ThreadNum = omp_get_thread_num();
			double ThreadDistSqr = PointDistanceToElementSquared(OldPoint, e, JunkVecVector[ThreadNum]);
			if (ThreadDistSqr < MinDistSqrIndVector[ThreadNum].first) {
				MinDistSqrIndVector[ThreadNum].first = ThreadDistSqr;
				MinDistSqrIndVector[ThreadNum].second = e;
			}
// 			TmpDistSqr = PointDistanceToElementSquared(OldPoint, e, NewPoint);
// 			if (TmpDistSqr < MinDistSqr) {
// 				MinDistSqr = TmpDistSqr;
// 				ProjectedElemIndex = e;
// 			}
		}

		auto MinDistSqrInd = std::min_element(MinDistSqrIndVector.begin(), MinDistSqrIndVector.end(), [](std::pair<double, int> a, std::pair<double, int> b) {return a.first < b.first; });
		ProjectedElemIndex = MinDistSqrInd->second;
	}

	if (ProjectedElemIndex >= 0) {
		PointDistanceToElementSquared(OldPoint, ProjectedElemIndex, NewPoint);
#if 0
		SaveVec3VecAsScatterZone({ OldPoint }, string("Element " + to_string(ProjectedElemIndex + 1) + " old point"), Red_C);
		TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.5, 0);
		TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
		SaveVec3VecAsScatterZone({ NewPoint }, string("Element " + to_string(ProjectedElemIndex + 1) + " new point"), Green_C);
		TecUtilZoneSetScatter(SV_FRAMESIZE, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.5, 0);
		TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);

		SaveVec3VecAsScatterZone({ OldPoint, NewPoint }, string("old to new point"), Blue_C);
		TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
		TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
		TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.2, 0);
		TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Blue_C);
		TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);

		SaveVec3VecAsScatterZone({ m_XYZList[m_ElemList[ProjectedElemIndex][0]], m_XYZList[m_ElemList[ProjectedElemIndex][1]] , m_XYZList[m_ElemList[ProjectedElemIndex][2]], m_XYZList[m_ElemList[ProjectedElemIndex][0]] }, string("Element " + to_string(ProjectedElemIndex + 1)));
		TecUtilZoneSetScatter(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, FALSE);
		TecUtilZoneSetMesh(SV_SHOW, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, TRUE);
		TecUtilZoneSetMesh(SV_LINETHICKNESS, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0.4, 0);
		TecUtilZoneSetMesh(SV_COLOR, tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), 0, Red_C);
		TecUtilZoneSetActive(tecplot::toolbox::Set(TecUtilDataSetGetNumZones()).getRef(), AssignOp_MinusEquals);
#endif
	}
	else
		TecUtilDialogErrMsg("Failed to find minimum triangle distance");
	// 	}

	IsOk = ProjectedElemIndex >= 0;
	if (!IsOk) return IsOk;

	return true;
}

int FESurface_c::SetZoneStyle(int const ZoneNum,
	AssignOp_e const ZoneActive,
	Boolean_t const ShowContour,
	Boolean_t const ShowMesh,
	Boolean_t const ShowScatter,
	Boolean_t const ShowShade)
{
	if (ZoneNum > 0 && ZoneNum <= TecUtilDataSetGetNumZones()){
		Set_pa ZoneSet = TecUtilSetAlloc(TRUE);
		TecUtilSetAddMember(ZoneSet, ZoneNum, TRUE);

		TecUtilZoneSetActive(ZoneSet, ZoneActive);
		TecUtilZoneSetMesh(SV_SHOW, ZoneSet, 0.0, ShowMesh);
		TecUtilZoneSetScatter(SV_SHOW, ZoneSet, 0.0, ShowScatter);
		TecUtilZoneSetContour(SV_SHOW, ZoneSet, 0.0, ShowContour);
		TecUtilZoneSetShade(SV_SHOW, ZoneSet, 0.0, ShowShade);
	}
	return ZoneNum;
}



/*
	Save an FEVolume_c that has it's elements and nodal structure
	already created.
	Only Saves the XYZ values.
*/
int FESurface_c::SaveAsTriFEZone(vector<int> const & XYZVarNums, 
	string ZoneName)
{
	Boolean_t IsOk = IsMade() && (m_XYZList.size() > 0 || m_RefinedXYZList.size() > 0);
	int ZoneNum = -1;

	if (ZoneName.length() == 0)
		ZoneName = "FE Volume";

	if (m_RefinedXYZList.size() <= 0) m_RefinedXYZList = m_XYZList;

	IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), m_RefinedXYZList.size(), m_ElemList.size(), 0, ZoneType_FETriangle, nullptr);

	Set_pa TmpSet = TecUtilSetAlloc(FALSE);

	if (IsOk){
		ZoneNum = TecUtilDataSetGetNumZones();
		m_ZoneNum = ZoneNum;
		TecUtilSetAddMember(TmpSet, ZoneNum, FALSE);
		TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
		TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);

		ArgList_pa CurrentArgList = TecUtilArgListAlloc();
		TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
		TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
		TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
		TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
		TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
		TecUtilStyleSetLowLevelX(CurrentArgList);
		TecUtilArgListDealloc(&CurrentArgList);

		TecUtilSetDealloc(&TmpSet);

		vector<vector<double> > TmpValues(3, vector<double>(m_RefinedXYZList.size()));
		for (int i = 0; i < m_RefinedXYZList.size(); ++i){
			for (int j = 0; j < 3; ++j)
				TmpValues[j][i] = m_RefinedXYZList[i][j];
		}

		for (int i = 0; i < 3 && IsOk; ++i){
			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
			IsOk = VALID_REF(SetFDPtr);
			if (IsOk){
				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_RefinedXYZList.size(), reinterpret_cast<void*>(const_cast<double*>(TmpValues[i].data())));
			}
		}

	}
	
	if (IsOk){
		NodeMap_pa NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
		IsOk = VALID_REF(NodeMap);
		for (int ei = 0; ei < m_ElemList.size(); ++ei){
			for (int i = 0; i < 3; ++i){
				TecUtilDataNodeSetByRef(NodeMap, ei + 1, i + 1, m_ElemList[ei][i] + 1);
			}
		}
	}

	TecUtilSetDealloc(&TmpSet);

	return ZoneNum;
}


/*
	Save a FEVolume_c that was made from multiple GPs
*/
int FESurface_c::SaveAsTriFEZone(string const & ZoneName, 
	vector<FieldDataType_e> DataTypes,
	vector<ValueLocation_e> const & DataLocations,
	vector<int> const & XYZVarNums,
	int RhoVarNum)
{
	Boolean_t IsOk = IsMade();
	int ZoneNum = -1;

	if (IsOk){
		if (DataTypes.size() > XYZVarNums[2]) for (int i = 0; i < 3; ++i)
			DataTypes[XYZVarNums[i] - 1] = FieldDataType_Double;

// 		IsOk = TecUtilDataSetAddZone(ZoneName.c_str(), m_XYZList.size(), m_ElemList.size(), 0, ZoneType_FETriangle, DataTypes.data());
		tecplot::toolbox::ArgList Args;
		Args.appendString(SV_NAME, ZoneName);
		Args.appendInt(SV_ZONETYPE, ZoneType_FETriangle);
		Args.appendInt(SV_IMAX, m_XYZList.size());
		Args.appendInt(SV_JMAX, m_ElemList.size());
		if (DataTypes.size() > 0)
			Args.appendArray(SV_VARDATATYPE, DataTypes.data());
		if (DataLocations.size() > 0)
			Args.appendArray(SV_VALUELOCATION, DataLocations.data());

		IsOk = TecUtilDataSetAddZoneX(Args.getRef());
	}

	tecplot::toolbox::Set TmpSet;

	if (IsOk){
		ZoneNum = TecUtilDataSetGetNumZones();
		m_ZoneNum = ZoneNum;
		TmpSet += ZoneNum;
		TecUtilZoneSetActive(TmpSet.getRef(), AssignOp_MinusEquals);
		TecUtilZoneSetScatter(SV_SHOW, TmpSet.getRef(), 0.0, FALSE);

		tecplot::toolbox::ArgList CurrentArgList;
		CurrentArgList.appendString(SV_P1, SV_FIELDMAP);
		CurrentArgList.appendString(SV_P2, SV_EFFECTS);
		CurrentArgList.appendString(SV_P3, SV_USETRANSLUCENCY);
		CurrentArgList.appendSet(SV_OBJECTSET, TmpSet);
		CurrentArgList.appendArbParam(SV_IVALUE, FALSE);
		TecUtilStyleSetLowLevelX(CurrentArgList.getRef());

		vector<vector<double> > TmpValues(3, vector<double>(m_XYZList.size()));
		for (int i = 0; i < m_XYZList.size(); ++i){
			for (int j = 0; j < 3; ++j)
				TmpValues[j][i] = m_XYZList[i][j];
		}

		for (int i = 0; i < 3 && IsOk; ++i){
			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
			IsOk = VALID_REF(SetFDPtr);
			if (IsOk){
				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_XYZList.size(), reinterpret_cast<void*>(const_cast<double*>(TmpValues[i].data())));
			}
		}
		if (IsOk && m_RhoList.size() == m_XYZList.size()){
			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, RhoVarNum);
			IsOk = VALID_REF(SetFDPtr);
			if (IsOk){
				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_NumNodes, reinterpret_cast<void*>(const_cast<double*>(m_RhoList.data())));
			}
		}

	}

	if (IsOk){
		/*
		*	Need to define the connectivity (nodal structure) of the FE volume.
		*	The order doesn't really matter here, just that the connectivity
		*	is correct between elements.
		*/
		NodeMap_pa NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
		IsOk = VALID_REF(NodeMap);
		if (IsOk){


			for (int ei = 0; ei < m_ElemList.size(); ++ei){
				for (int i = 0; i < 3; ++i){
					TecUtilDataNodeSetByRef(NodeMap, ei + 1, i + 1, m_ElemList[ei][i] + 1);
				}
			}
// 			int ei = 1;
// 			int ni = 1;
// 			vector<int> t;
// 			//	beginning cap
// 			for (int i = 1; i < m_NumGPs - 1; ++i){
// 				t = {ni, ni + i, ni + i + 1};
// 				for (int j = 0; j < 3; ++j){
// 					TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 				}
// 				ei++;
// 			}
// 			
// 			//	triangle's up sides of gradient bundle.
// 			//	Each while iteration is one more point moving
// 			//	down each GP in the GB.
// 			//	The ti for loop is the first m_NumGPs sides
// 			//	(2 triangular elements per side) of the GB
// 			//	and then after the loop the last side.
// 			while (ni < m_NumNodes - m_NumGPs && ei < m_NumElems){
// 				for (int ti = 0; ti < m_NumGPs - 1; ++ti){
// 					t = {ni, ni + m_NumGPs, ni + m_NumGPs + 1};
// 					for (int j = 0; j < 3; ++j){
// 						TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 					}
// 					ei++;
// 					t = {ni, ni + m_NumGPs + 1, ni + 1};
// 					for (int j = 0; j < 3; ++j){
// 						TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 					}
// 					ni++;
// 					ei++;
// 				}
// 				t = {ni, ni + m_NumGPs, ni + 1};
// 				for (int j = 0; j < 3; ++j){
// 					TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 				}
// 				ei++;
// 				t = {ni, ni + 1, ni - m_NumGPs + 1};
// 				for (int j = 0; j < 3; ++j){
// 					TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 				}
// 				ni++;
// 				ei++;
// 			}
// 			
// 			//	end cap
//   // 			ni -= m_NumGPs; 
// 			//	set ni back to the last node of the 
// 			//	"first" GP in the GB.
// 			for (int i = 1; i < m_NumGPs - 1; ++i){
// 				t = {ni, ni + i, ni + i + 1};
// 				for (int j = 0; j < 3; ++j){
// 					TecUtilDataNodeSetByRef(NodeMap, ei, j + 1, t[j]);
// 				}
// 				ei++;
// 			}
		}
	}

	return ZoneNum;
}


// int const FESurface_c::SaveAsFEZone(vector<FieldDataType_e> DataTypes,
// 	vector<int> const & XYZVarNums,
// 	int RhoVarNum)
// {
// 	Boolean_t IsOk = IsMade();
// 	int ZoneNum = -1;
// 
// 	if (IsOk){
// 		for (int i = 0; i < 3; ++i)
// 			DataTypes[XYZVarNums[i] - 1] = FieldDataType_Double;
// 
// 		DataTypes[RhoVarNum - 1] = FieldDataType_Double;
// 
// 		IsOk = TecUtilDataSetAddZone(CSMZoneName.FESurface.c_str(), m_XYZList.size(), m_NumElems, 0, ZoneType_FEQuad, DataTypes.data());
// 	}
// 
// 	Set_pa TmpSet = TecUtilSetAlloc(FALSE);
// 
// 	if (IsOk){
// 		ZoneNum = TecUtilDataSetGetNumZones();
// 		m_ZoneNum = ZoneNum;
// 		TecUtilSetAddMember(TmpSet, ZoneNum, FALSE);
// 		TecUtilZoneSetActive(TmpSet, AssignOp_MinusEquals);
// 		TecUtilZoneSetScatter(SV_SHOW, TmpSet, 0.0, FALSE);
// 
// 		ArgList_pa CurrentArgList = TecUtilArgListAlloc();
// 		TecUtilArgListAppendString(CurrentArgList, SV_P1, SV_FIELDMAP);
// 		TecUtilArgListAppendString(CurrentArgList, SV_P2, SV_EFFECTS);
// 		TecUtilArgListAppendString(CurrentArgList, SV_P3, SV_USETRANSLUCENCY);
// 		TecUtilArgListAppendSet(CurrentArgList, SV_OBJECTSET, TmpSet);
// 		TecUtilArgListAppendArbParam(CurrentArgList, SV_IVALUE, FALSE);
// 		TecUtilStyleSetLowLevelX(CurrentArgList);
// 		TecUtilArgListDealloc(&CurrentArgList);
// 
// 		TecUtilSetDealloc(&TmpSet);
// 
// 		vector<vector<double> > TmpValues(3, vector<double>(m_XYZList.size()));
// 		for (int i = 0; i < m_XYZList.size(); ++i){
// 			for (int j = 0; j < 3; ++j)
// 				TmpValues[j][i] = m_XYZList[i][j];
// 		}
// 
// 		for (int i = 0; i < 3 && IsOk; ++i){
// 			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, XYZVarNums[i]);
// 			IsOk = VALID_REF(SetFDPtr);
// 			if (IsOk){
// 				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_XYZList.size(), reinterpret_cast<void*>(const_cast<double*>(TmpValues[i].data())));
// 			}
// 		}
// 		if (IsOk){
// 			FieldData_pa SetFDPtr = TecUtilDataValueGetWritableNativeRef(ZoneNum, RhoVarNum);
// 			IsOk = VALID_REF(SetFDPtr);
// 			if (IsOk){
// 				TecUtilDataValueArraySetByRef(SetFDPtr, 1, m_XYZList.size(), reinterpret_cast<void*>(const_cast<double*>(m_RhoList.data())));
// 			}
// 		}
// 
// 	}
// 
// 	if (IsOk){
// 		/*
// 		*	Need to define the connectivity (nodal structure) of the FE volume.
// 		*	The order doesn't really matter here, just that the connectivity
// 		*	is correct between elements.
// 		*/
// 		NodeMap_pa NodeMap = TecUtilDataNodeGetWritableRef(ZoneNum);
// 		IsOk = VALID_REF(NodeMap);
// 		if (IsOk){
// 			if (m_NumGPs == 3){
// 				int ei = 1;
// 				//	Beginning triangular cap
// 				for (int i = 1; i <= 4; ++i){
// 					// 4th corner gets repeated
// 					TecUtilDataNodeSetByRef(NodeMap, ei, i, MIN(i, 3));
// 				}
// 				ei++;
// 				/*
// 				*	Now work down the FE volume, defining the elements starting
// 				*	at the sphere, wrapping around the FE volume as you work
// 				*	towards the FE volume's end.
// 				*/
// 				for (int i = 0; i < m_NumGPPts - 1 && ei < m_NumElems - 1; ++i){
// 					int ii = 1 + 3 * i;
// 					for (int j = 0; j < 2 && ei < m_NumElems - 1; ++j){
// 						/*
// 						*	First two elements of each iteration have
// 						*	the same relative nodal structure, offset by
// 						*	one for each element.
// 						*/
// 						int Vals[4] = { ii, ii + 1, ii + 4, ii + 3 };
// 						for (int k = 0; k < 4; ++k){
// 							TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
// 						}
// 						ii++;
// 						ei++;
// 					}
// 					/*
// 					*	The third element has a different nodal structure
// 					*/
// 					int Vals[4] = { ii, ii - 2, ii + 1, ii + 3 };
// 					for (int k = 0; k < 4; ++k){
// 						TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
// 					}
// 					ei++;
// 				}
// 				//	End triangular cap
// 				for (int i = 1; i <= 4; ++i){
// 					// 4th corner gets repeated
// 					TecUtilDataNodeSetByRef(NodeMap, ei, i, MIN(i + 3 * (m_NumGPPts - 1), m_NumNodes));
// 				}
// 			}
// 			else{
// 				/*
// 				*	Same thing for the 4-sided FE volumes, but the end cap
// 				*	is a quad and the "walls" of the volume have three with
// 				*	the same nodal structure instead of two.
// 				*/
// 				int ei = 1;
// 				for (int i = 1; i <= 4; ++i){
// 					TecUtilDataNodeSetByRef(NodeMap, ei, i, MAX(i % 4, 1));
// 				}
// 				ei++;
// 				for (int i = 0; i < m_NumGPPts - 1 && ei < m_NumElems - 1; ++i){
// 					int ii = 1 + 4 * i;
// 					for (int j = 0; j < 3 && ei < m_NumElems - 1; ++j){
// 						int Vals[4] = { ii, ii + 1, ii + 5, ii + 4 };
// 						for (int k = 0; k < 4; ++k){
// 							TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
// 						}
// 						ii++;
// 						ei++;
// 					}
// 					int Vals[4] = { ii, ii - 3, ii + 1, ii + 4 };
// 					for (int k = 0; k < 4; ++k){
// 						TecUtilDataNodeSetByRef(NodeMap, ei, k + 1, Vals[k]);
// 					}
// 					ei++;
// 				}
// 				//	End quadragonal cap
// 				for (int i = 1; i <= 4; ++i){
// 					// 4th corner not repeated
// 					TecUtilDataNodeSetByRef(NodeMap, ei, i, i + 4 * (m_NumGPPts - 1));
// 				}
// 			}
// 		}
// 	}
// 
// 	TecUtilSetDealloc(&TmpSet);
// 
// 	return ZoneNum;
// }



/*
*	Setter methods
*/

void FESurface_c::GetNodeConnectivityFromTecplot(){
	/*
	*	Get connectivity list for FE zone
	*/
	if (m_FEVolumeMade) {
		TecUtilDataNodeGetReadableRawPtr(m_ZoneNum, &m_ConnectivityListPtr);
		m_FEVolumeMade = (m_ConnectivityListPtr != nullptr);
	}
	if (m_FEVolumeMade) {
		NodeToElemMap_pa NodeToElemMap = TecUtilDataNodeToElemMapGetReadableRef(m_ZoneNum);
		m_NodeConnectivityList.resize(m_NumNodes);
		ElemFaceOffset_t FaceOffset = 0;
		LgIndex_t NumUniqueNodes, UniqueNodesSize = 0, *UniqueNodes = nullptr;
		for (LgIndex_t NodeNum = 0; NodeNum < m_NumNodes; ++NodeNum) {
			m_NodeConnectivityList[NodeNum].reserve(8);
			LgIndex_t NumElems = TecUtilDataNodeToElemMapGetNumElems(NodeToElemMap, NodeNum + 1);
			for (LgIndex_t ElemNum = 1; ElemNum <= NumElems; ++ElemNum) {
				LgIndex_t ElemIndex = TecUtilDataNodeToElemMapGetElem(NodeToElemMap, NodeNum + 1, ElemNum);
				TecUtilDataFECellGetUniqueNodes(m_ZoneNum, FaceOffset, ElemIndex, &NumUniqueNodes, &UniqueNodesSize, &UniqueNodes);
				for (LgIndex_t eNodeNum = 0; eNodeNum < UniqueNodesSize; ++eNodeNum) {
					if (UniqueNodes[eNodeNum] > 0 && NodeNum + 1 != UniqueNodes[eNodeNum] && std::find(m_NodeConnectivityList[NodeNum].begin(), m_NodeConnectivityList[NodeNum].end(), UniqueNodes[eNodeNum] - 1) == m_NodeConnectivityList[NodeNum].end()) {
						m_NodeConnectivityList[NodeNum].push_back(UniqueNodes[eNodeNum] - 1);
					}
				}
			}
		}
		TecUtilArrayDealloc((void**)(&UniqueNodes));
	}
}


void FESurface_c::GeneratePointElementDistanceCheckData() {
 	m_NumElems = m_ElemList.size();
 
 	if (m_edge0.size() != m_NumElems) {
 		m_edge0.clear(); m_edge0.resize(m_NumElems);
 		m_edge1.clear(); m_edge1.resize(m_NumElems);
 		m_a.clear(); m_a.resize(m_NumElems);
 		m_b.clear(); m_b.resize(m_NumElems);
 		m_c.clear(); m_c.resize(m_NumElems);
 		m_det.clear(); m_det.resize(m_NumElems);
 		m_oneOverA.clear(); m_oneOverA.resize(m_NumElems);
 		m_oneOverC.clear(); m_oneOverC.resize(m_NumElems);
 		m_invDet.clear(); m_invDet.resize(m_NumElems);
 		m_oneOverDenom.clear(); m_oneOverDenom.resize(m_NumElems);
 
 #pragma omp parallel for
 		for (int i = 0; i < m_NumElems; ++i) {
 			m_edge0[i] = m_XYZList[m_ElemList[i][1]] - m_XYZList[m_ElemList[i][0]];
 			m_edge1[i] = m_XYZList[m_ElemList[i][2]] - m_XYZList[m_ElemList[i][0]];
 			m_a[i] = dot2(m_edge0[i]);
 			m_b[i] = dot(m_edge0[i], m_edge1[i]);
 			m_c[i] = dot2(m_edge1[i]);
 			m_det[i] = m_a[i] * m_c[i] - m_b[i] * m_b[i];
 			m_oneOverA[i] = 1.0 / m_a[i];
 			m_oneOverC[i] = 1.0 / m_c[i];
 			m_invDet[i] = 1.0 / m_det[i];
 			m_oneOverDenom[i] = 1.0 / (m_a[i] - 2.0 * m_b[i] + m_c[i]);
 		}
 	}
	GenerateElemConnectivity();
}

void FESurface_c::GenerateElemMidpoints(){
	if (m_ElemMidPoints.size() != m_ElemList.size()) {
		m_ElemMidPoints.clear();
		m_ElemMidPoints.reserve(m_ElemList.size());
		for (auto const & e : m_ElemList) {
			if (e.size() > 1) {
				m_ElemMidPoints.push_back(m_XYZList[e[0]]);
				for (int i = 1; i < e.size(); ++i) m_ElemMidPoints.back() += m_XYZList[e[i]];
				m_ElemMidPoints.back() /= e.size();
			}
		}
	}
}

void FESurface_c::GenerateElemConnectivity() {
	if (m_ElemMidPoints.size() != m_ElemList.size()) {
		GetTriElementConnectivityList(&m_ElemList, m_ElemConnectivityList, 1);
	}
}

Boolean_t FESurface_c::Setup(int InZoneNum,
	int VolZoneNum,
	vector<int> const & InXYZVarNums,
	vector<int> const & InIntVarNums,
	bool const CopyData)
{
	m_FEVolumeMade = (InZoneNum > 0
		&& VolZoneNum > 0
		&& InZoneNum != VolZoneNum
		&& TecUtilZoneIsFiniteElement(InZoneNum)
		&& InXYZVarNums.size() == 3);

	if (m_FEVolumeMade){
		m_ZoneNum = InZoneNum;
		ZoneType_e ZoneType = TecUtilZoneGetType(InZoneNum);
		TecUtilZoneGetIJK(InZoneNum, &m_NumNodes, &m_NumElems, &m_NumNodesPerElem);
		m_FEVolumeMade = (m_NumNodes > 0 && m_NumElems > 0
			&& (ZoneType == ZoneType_FETriangle && m_NumNodesPerElem == 3)
			|| (ZoneType == ZoneType_FEQuad && m_NumNodesPerElem == 4));
	}

	if (m_FEVolumeMade){
		GetVolInfo(VolZoneNum, InXYZVarNums, FALSE, m_VolZoneInfo);
	}

	/*
	*	Get pointers to necessary data
	*/
	if (m_FEVolumeMade){
		m_IntVarNums = InIntVarNums;
		m_XYZPtrs.resize(3);
		m_NumIntVars = static_cast<int>(InIntVarNums.size());
		m_IntVarPtrs.resize(m_NumIntVars);

		for (int i = 0; i < 3 && m_FEVolumeMade; ++i)
			m_FEVolumeMade = m_XYZPtrs[i].InitializeReadPtr(InZoneNum, InXYZVarNums[i]);

		for (int i = 0; i < m_NumIntVars && m_FEVolumeMade; ++i)
			m_FEVolumeMade = m_IntVarPtrs[i].InitializeReadPtr(VolZoneNum, m_IntVarNums[i]);
	}


	if (m_FEVolumeMade && CopyData) {
		m_XYZList.resize(m_NumNodes);
		for (int i = 0; i < m_NumNodes; ++i){
			for (int j = 0; j < 3; ++j)
				m_XYZList[i][j] = m_XYZPtrs[j][i];
		}
		TecUtilDataNodeGetReadableRawPtr(m_ZoneNum, &m_ConnectivityListPtr);
		m_ElemList.resize(m_NumElems);
		for (int i = 0; i < m_NumElems; ++i){
			m_ElemList[i].reserve(m_NumNodesPerElem);
			for (int j = 0; j < m_NumNodesPerElem; ++j)
				m_ElemList[i].push_back(m_ConnectivityListPtr[i * m_NumNodesPerElem + j]);
		}
		if (m_NumNodesPerElem == 4){
			vector<vector<int> > NewElems;
			for (vector<int> const & e : m_ElemList){
				NewElems.push_back({ e[0], e[1], e[2] });
				NewElems.push_back({ e[0], e[2], e[3] });
			}
			m_ElemList = NewElems;
		}

// 		RemoveDupicateNodes();
	}

//  	 	GetNodeConnectivityFromTecplot();
	// 	GenerateElemConnectivity();

	m_IntegrationResultsReady = FALSE;

	return m_FEVolumeMade;
}

/*
 *	Make FEVolume from an existing zone.
 *	Can be IJ ordered zone or actual FE zone.
 *	If IJ ordered zone, need to convert to FE triangle zone representation
 *		with connectivity info and such. 
 *	This assumes the structure of the surfaces made by Bondalyzer
 */
Boolean_t FESurface_c::Setup(int const InZoneNum,
								vector<int> const & InXYZVarNums){
	REQUIRE(1 <= InZoneNum && InZoneNum <= TecUtilDataSetGetNumZones());
	REQUIRE(InXYZVarNums.size() == 3);
	for (int i : InXYZVarNums) 
		REQUIRE(1 <= i && i <= TecUtilDataSetGetNumVars());

	m_FEVolumeMade = FALSE;

// 	if (TecUtilZoneIsActive(InZoneNum)){
		vector<int> IJK(3);
		m_ZoneNum = InZoneNum;

		TecUtilZoneGetIJK(m_ZoneNum, &IJK[0], &IJK[1], &IJK[2]);

		if (TecUtilZoneIsOrdered(m_ZoneNum) && (1 < IJK[0] && 1 < IJK[1] && 1 == IJK[2])){
			m_XYZList.resize(IJK[0] * IJK[1]);
			m_ElemList.reserve(m_XYZList.size() * 2);

			m_XYZPtrs.resize(3);
			for (int i = 0; i < 3; ++i)
				m_XYZPtrs[i].InitializeReadPtr(m_ZoneNum, InXYZVarNums[i]);

			for (int j = 0; j < IJK[1]; ++j){
				for (int i = 0; i < IJK[0]; ++i){
					int Index = IJK[0] * j + i;
					for (int Dir = 0; Dir < 3; ++Dir)
						m_XYZList[Index][Dir] = m_XYZPtrs[Dir][Index];

					if (i < IJK[0] - 1 && j < IJK[1] - 1){
						vector<int> Inds;
						Inds.reserve(4);
						for (int jj = 0; jj < 2; ++jj){
							for (int ii = 0; ii < 2; ++ii){
								Inds.push_back(IJK[0] * (j + jj) + (i + ii));
							}
						}
						m_ElemList.push_back(vector<LgIndex_t>({ Inds[0], Inds[1], Inds[3] }));
						m_ElemList.push_back(vector<LgIndex_t>({ Inds[3], Inds[2], Inds[0] }));
					}
				}
			}
			m_NumNodes = m_XYZList.size();
			m_NumElems = m_ElemList.size();
			m_FEVolumeMade = TRUE;
			RemoveDupicateNodes();
			return TRUE;
		}
		else if (!TecUtilZoneIsOrdered(m_ZoneNum)){
			return Setup(InZoneNum, ZoneNumByName("Full Volume"), InXYZVarNums, vector<int>(), true);
		}
// 	}
	return FALSE;
}

/*
 *	Remove duplicate nodes (nodes with same position)
 *	and revise element list to reflect the removal of
 *	the nodes.
 */
void FESurface_c::RemoveDupicateNodes(){
	RemoveDupicateNodesFromMesh(m_XYZList, m_ElemList);
}


Boolean_t FESurface_c::Refine(){
	if (m_RefinedXYZList.size() == 0)
		m_RefinedXYZList = m_XYZList;

	m_RefinedXYZList.reserve(m_RefinedXYZList.size() * 2);
	m_ElemList.reserve(m_ElemList.size() * 4);

	for (int i = 0; i < m_NumElems; ++i)
		TriangleEdgeMidPointSubdivide(m_RefinedXYZList, m_ElemList, i, vector<int>());

	m_XYZList = m_RefinedXYZList;
	m_NumNodes = m_RefinedXYZList.size();
	m_NumElems = m_ElemList.size();

	RemoveDupicateNodes();

	m_RefinedXYZList = m_XYZList;
	m_NumNodes = m_RefinedXYZList.size();
	m_NumElems = m_ElemList.size();

	return TRUE;
}


Boolean_t FESurface_c::DoIntegration(Boolean_t IntegrateVolume,
	vector<vec> const & stuW,
	vector<int> const & SplitPtNums,
	vector<vec> const & stuW2)
{
	Boolean_t IsOk = m_FEVolumeMade && m_NumNodesPerElem == 4;

	if (IsOk){
		Boolean_t IsBond = m_GPList.size() > 3;
		Boolean_t IsCage;
		if (IsBond){

		}
		else{
			/*
			*	Check to see if gradient bundle ends at a cage point.
			*	First check area of triandle made by last three points
			*	(should be 0 if the come together at a cage point).
			*/
		}


		vector<mat> xyzW = GetIntegrationPointsWeights(stuW);

// 		string TempStr = "Negative weights at i = {", TempStr2 = TempStr;
// 		for (int i = 0; i < xyzW[0].n_rows; ++i){
// 			if (xyzW[1](i,0) < 0){
// 				TempStr += to_string(i) + ",";
// 			}
// 		}
// 
// 		if (TempStr != TempStr2){
// 			TecUtilDialogMessageBox(TempStr.c_str(), MessageBoxType_Error);
// 			TempStr = "";
// 		}

// 		vector<vector<double> > vxyzW(xyzW[0].n_rows, vector<double>(3));
// 		vector<double> vW(xyzW[1].n_rows);
// 		for (int i = 0; i < 3; ++i){
// 			for (int j = 0; j < xyzW[0].n_rows; ++j){
// 				vxyzW[j][i] = xyzW[0](j, i);
// 			}
// 		}
// 		for (int i = 0; i < vW.size(); ++i){
// 			vW[i] = xyzW[1](i, 0);
// 		}

		m_IntValues.resize(m_NumIntVars + int(IntegrateVolume), 0.0);

		for (int i = 0; i < xyzW[0].n_rows; ++i){
			vec3 TmpVec = xyzW[0].submat(i,0,i,2).t();
			SetIndexAndWeightsForPoint(TmpVec, m_VolZoneInfo);
			if (IntegrateVolume){
				m_IntValues[m_NumIntVars] += xyzW[1](i, 0);
			}
			for (int j = 0; j < m_NumIntVars; ++j){
				m_IntValues[j] += xyzW[1](i, 0) * ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[j]);
			}
		}

		m_IntegrationResultsReady = TRUE;
	}

	return IsOk;
}

// Boolean_t const FEVolume_c::DoIntegration(int ResolutionScale, Boolean_t IntegrateVolume)
// {
// 	Boolean_t IsOk = m_FEVolumeMade;
// 
// 	if (IsOk && !m_IntegrationResultsReady){
// 		m_ZoneMinXYZ = m_VolZoneInfo.MaxXYZ;
// 		m_ZoneMaxXYZ = m_VolZoneInfo.MinXYZ;
// 		/*
// 		*	Get FE zone's max and min XYZ values
// 		*/
// 		for (int i = 0; i < 3; ++i){
// 			for (int j = 0; j < m_NumNodes; ++j){
// 				if (m_XYZPtrs[i][j] < m_ZoneMinXYZ[i])
// 					m_ZoneMinXYZ[i] = m_XYZPtrs[i][j];
// 				else if (m_XYZPtrs[i][j] > m_ZoneMaxXYZ[i])
// 					m_ZoneMaxXYZ[i] = m_XYZPtrs[i][j];
// 			}
// 		}
// 
// 		IsOk = sum(m_ZoneMaxXYZ > m_ZoneMinXYZ) == 3;
// 	}
// 
// 	/*
// 	if (IsOk){
// 	/ *
// 	*	Loop over every point in the system checking if it's interior
// 	*	or not. If interior, add it to the integral value.
// 	* /
// 	vec3 PointSpacingV123, TmpPoint;
// 
// 	for (int i = 0; i < 3; ++i){
// 	PointSpacingV123[i] = (m_VolZoneInfo.MaxXYZ[i] - m_VolZoneInfo.MinXYZ[i]) / static_cast<double>(m_VolZoneInfo.MaxIJK[i] - 1);
// 	}
// 	m_IntValues.resize(m_NumIntVars, 0.0);
// 	int IJK[3];
// 	for (IJK[2] = 0; IJK[2] < m_VolZoneInfo.MaxIJK[2]; ++IJK[2]){
// 	for (IJK[1] = 0; IJK[1] < m_VolZoneInfo.MaxIJK[1]; ++IJK[1]){
// 	for (IJK[0] = 0; IJK[0] < m_VolZoneInfo.MaxIJK[0]; ++IJK[0]){
// 	for (int i = 0; i < 3; ++i){
// 	TmpPoint[i] = m_VolZoneInfo.MinXYZ[i] + static_cast<double>(IJK[i]) * PointSpacingV123[i];
// 	}
// 	if (PointIsInterior(TmpPoint)){
// 	SetIndexAndWeightsForPoint(TmpPoint, m_VolZoneInfo);
// 	for (int i = 0; i < m_NumIntVars; ++i){
// 	if (IntVarFDTypes[i] == FieldDataType_Float)
// 	m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, IntVarPtrsFlt[i]);
// 	else
// 	m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, IntVarPtrsDbl[i]);
// 	}
// 	}
// 	}
// 	}
// 	}
// 	}
// 	*/
// 
// 
// 	if (IsOk){
// 		/*
// 		*	Loop over each cell in the full system volume.
// 		*	Maintain list of bools to see which nodes are in/out
// 		*	of the FE volume.
// 		*	If whole cell is in FE volume, add whole cell to integral.
// 		*	If part of cell in FE volume, subsample cell and add pieces
// 		*	of it to integral if they're in FE volume.
// 		*/
// 		vec3 FarPoint = m_VolZoneInfo.MinXYZ - 1.;
// 		
// 		vec3 PointSpacingV123 = m_VolZoneInfo.PointSpacingV123, SubDelXYZ, TmpPoint, TmpSubPoint;
// 		double CellVolume = 1.0;
// 
// 		for (int i = 0; i < 3; ++i){
// // 			PointSpacingV123[i] = (m_VolZoneInfo.MaxXYZ[i] - m_VolZoneInfo.MinXYZ[i]) / static_cast<double>(m_VolZoneInfo.MaxIJK[i] - 1);
// 			CellVolume *= PointSpacingV123[i];
// 		}
// 
// 		SubDelXYZ = PointSpacingV123 / static_cast<double>(ResolutionScale);
// 		double SubCellFactor = 1.0 / static_cast<double>(ResolutionScale * ResolutionScale * ResolutionScale);
// 		double SubCellVolume = CellVolume * SubCellFactor;
// 
// 		if (IntegrateVolume)
// 			m_IntValues = vector<double>(m_NumIntVars + 1, 0.0);
// 		else
// 			m_IntValues = vector<double>(m_NumIntVars, 0.0);
// 
// 		vector<Boolean_t> VolNodeInFEVolume(m_VolZoneInfo.MaxIJK[0] * m_VolZoneInfo.MaxIJK[1] * m_VolZoneInfo.MaxIJK[2]);
// 		vector<Boolean_t> VolNodeChecked(m_VolZoneInfo.MaxIJK[0] * m_VolZoneInfo.MaxIJK[1] * m_VolZoneInfo.MaxIJK[2], FALSE);
// 		int IJK[3], TmpIJK[3];
// 		int Index[8];
// 		int NumInteriorNodes;
// 		for (IJK[2] = 1; IJK[2] < m_VolZoneInfo.MaxIJK[2]; ++IJK[2]){
// 			for (IJK[1] = 1; IJK[1] < m_VolZoneInfo.MaxIJK[1]; ++IJK[1]){
// 				for (IJK[0] = 1; IJK[0] < m_VolZoneInfo.MaxIJK[0]; ++IJK[0]){
// 					/*
// 					*	Get I-D index numbers for the nodes of the cell.
// 					*/
// 					Index[0] = IndexFromIJK(IJK[0], IJK[1], IJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[1] = IndexFromIJK(IJK[0] + 1, IJK[1], IJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[2] = IndexFromIJK(IJK[0] + 1, IJK[1] + 1, IJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[3] = IndexFromIJK(IJK[0], IJK[1] + 1, IJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[4] = IndexFromIJK(IJK[0], IJK[1], IJK[2] + 1, m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[5] = IndexFromIJK(IJK[0] + 1, IJK[1], IJK[2] + 1, m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[6] = IndexFromIJK(IJK[0] + 1, IJK[1] + 1, IJK[2] + 1, m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 					Index[7] = IndexFromIJK(IJK[0], IJK[1] + 1, IJK[2] + 1, m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
// 
// 					// 					if (IJK[0] == 40 && IJK[1] == 40 && IJK[2] == 40){
// 					// 						int a = 1;
// 					// 					}
// 
// 
// 					/*
// 					*	Check the number of cell nodes that are interior to the volume.
// 					*/
// 					NumInteriorNodes = 0;
// 
// 					for (int i = 0; i < 8; ++i){
// 						if (!VolNodeChecked[Index[i]]){
// 							/*
// 							*	Node has not been check as interior or not, so get the XYZ
// 							*	coordinates and check if interior.
// 							*/
// 							// Get IJK for corner of cell (subtract 1 since IJK were 1-based
// 							for (int j = 0; j < 3; ++j)
// 								TmpIJK[j] = IJK[j] - 1;
// 
// 							GetCellCornerIndices(i, TmpIJK[0], TmpIJK[1], TmpIJK[2]);
// 
// 							for (int j = 0; j < 3; ++j)
// 								TmpPoint[j] = m_VolZoneInfo.MinXYZ[j] + static_cast<double>(TmpIJK[j]) * PointSpacingV123[j];
// 
// 							//Boolean_t IsInside = VolNodeInFEVolume[Index[i]] = PointIsInterior(TmpPoint);
// 							VolNodeInFEVolume[Index[i]] = PointIsInterior(TmpPoint, FarPoint);
// 							VolNodeChecked[Index[i]] = TRUE;
// 						}
// 
// 						if (VolNodeInFEVolume[Index[i]]) NumInteriorNodes++;
// 					}
// 
// 					if (NumInteriorNodes == 8){
// 						/*
// 						*	Entire cell is in FE volume, so get midpoint of the cell and add that
// 						*	value to the integral.
// 						*/
// 
// 						/*
// 						*	Get the midpoint of the cell by taking the 0,0,0 corner and adding 0.5*dXYZ to it.
// 						*/
// 						for (int i = 0; i < 3; ++i)
// 							TmpPoint[i] = m_VolZoneInfo.MinXYZ[i] + (static_cast<double>(IJK[i]) - 0.5) * PointSpacingV123[i];
// 
// 						/*
// 						*	Add the value of each integrating variable for the cell midpoint to the
// 						*	integral value for that variable.
// 						*	This assumes that the value is already in the correct units for the system.
// 						*/
// 						if (SetIndexAndWeightsForPoint(TmpPoint, m_VolZoneInfo)){
// 							if (IntegrateVolume)
// 								m_IntValues[m_NumIntVars] += CellVolume;
// 							for (int i = 0; i < m_NumIntVars; ++i){
// 								m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * CellVolume;
// 							}
// 						}
// 						else
// 							TecUtilDialogErrMsg("Error: Indices out of bounds");
// 
// 					}
// 					else if (NumInteriorNodes > 0){
// 						/*
// 						*	Cell is partially contained in FE volume, so need to
// 						*	subsample the cell, adding values of subcells to integral
// 						*	as they're found.
// 						*
// 						*	Since we're testing points INSIDE the cell, need to start at 1/2 of
// 						*	the sub-resolution into the cell so that the subcells of density we add or don't add
// 						*	are entirely contained in the cell and span the cell.
// 						*/
// 						for (int i = 0; i < 3; ++i)
// 							TmpPoint[i] = m_VolZoneInfo.MinXYZ[i] + static_cast<double>(IJK[i] - 1) * PointSpacingV123[i] + (SubDelXYZ[i] / 2.0);
// 
// 						for (TmpIJK[0] = 0; TmpIJK[0] < ResolutionScale; ++TmpIJK[0]){
// 							for (TmpIJK[1] = 0; TmpIJK[1] < ResolutionScale; ++TmpIJK[1]){
// 								for (TmpIJK[2] = 0; TmpIJK[2] < ResolutionScale; ++TmpIJK[2]){
// 									for (int i = 0; i < 3; ++i)
// 										TmpSubPoint[i] = TmpPoint[i] + SubDelXYZ[i] * static_cast<double>(TmpIJK[i]);
// 
// 									if (PointIsInterior(TmpSubPoint, FarPoint)){
// 										if (SetIndexAndWeightsForPoint(TmpSubPoint, m_VolZoneInfo)){
// 											if (IntegrateVolume)
// 												m_IntValues[m_NumIntVars] += SubCellVolume;
// 											for (int i = 0; i < m_NumIntVars; ++i){
// 												m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * SubCellVolume;
// 											}
// 										}
// 										else
// 											TecUtilDialogErrMsg("Error: Indices out of bounds");
// 									}
// 								}
// 							}
// 						}
// 					}
// 				}
// 
// 			}
// 		}
// 
// 		m_IntegrationResultsReady = TRUE;
// 	}
// 
// 	return IsOk;
// }

Boolean_t FESurface_c::DoIntegrationNew(int ResolutionScale, Boolean_t IntegrateVolume)
{
	Boolean_t IsOk = m_FEVolumeMade;

	if (IsOk && !m_IntegrationResultsReady){
// 		m_ZoneMinXYZ = m_VolZoneInfo.MaxXYZ;
// 		m_ZoneMaxXYZ = m_VolZoneInfo.MinXYZ;
// 		/*
// 		*	Get FE zone's max and min XYZ values
// 		*/
// 		for (int i = 0; i < 3; ++i){
// 			for (int j = 0; j < m_NumNodes; ++j){
// 				if (m_XYZPtrs[i][j] < m_ZoneMinXYZ[i])
// 					m_ZoneMinXYZ[i] = m_XYZPtrs[i][j];
// 				else if (m_XYZPtrs[i][j] > m_ZoneMaxXYZ[i])
// 					m_ZoneMaxXYZ[i] = m_XYZPtrs[i][j];
// 			}
// 		}

		m_ZoneMinXYZ = ones<vec>(3);
		m_ZoneMaxXYZ = zeros<vec>(3);

		for (int n = 0; n < m_NumNodes; ++n){
			vec3 tmpPt;
			tmpPt << m_XYZPtrs[0][n] << m_XYZPtrs[1][n] << m_XYZPtrs[2][n];
			tmpPt = m_VolZoneInfo.BasisInverse * (tmpPt - m_VolZoneInfo.MinXYZ);
			for (int i = 0; i < 3; ++i){
				if (tmpPt[i] < m_ZoneMinXYZ[i])
					m_ZoneMinXYZ[i] = tmpPt[i];
				else if (tmpPt[i] > m_ZoneMaxXYZ[i])
					m_ZoneMaxXYZ[i] = tmpPt[i];
			}
		}

		m_ZoneMinXYZ = m_VolZoneInfo.BasisVectors * m_ZoneMinXYZ + m_VolZoneInfo.MinXYZ;
		m_ZoneMaxXYZ = m_VolZoneInfo.BasisVectors * m_ZoneMaxXYZ + m_VolZoneInfo.MinXYZ;

		IsOk = sum(m_ZoneMaxXYZ > m_ZoneMinXYZ) == 3;
	}

	vector<int> ZoneMinIJK, ZoneMaxIJK, ZoneNumIJK(3);
	int ZoneNumPts = 1;

	/*
	*	Get the IJK indices that correspond to the min and max XYZ for the zone
	*/
	int tmpZoneNum = m_XYZPtrs[0].ZoneNum();
	if (IsOk){
		ZoneMinIJK = GetIJKForPoint(m_ZoneMinXYZ, m_VolZoneInfo);
		ZoneMaxIJK = GetIJKForPoint(m_ZoneMaxXYZ, m_VolZoneInfo);
		for (int i = 0; i < 3; ++i){
			for (int j = 0; j < 2; ++j){
// 				if (ZoneMinIJK[i] > 1)
// 					ZoneMinIJK[i]--;
				if (ZoneMaxIJK[i] < m_VolZoneInfo.MaxIJK[i])
					ZoneMaxIJK[i]++;
			}
			ZoneNumIJK[i] = ZoneMaxIJK[i] - ZoneMinIJK[i] + 1;
			ZoneNumPts *= ZoneNumIJK[i];
		}
	}



	vector<Boolean_t> DoCheckNode(ZoneNumPts,FALSE), NodeIsInterior(ZoneNumPts, TRUE);
	vector<double> CheckNodeMinDistSqrToSurfaceNode(ZoneNumPts, 0.0);
	vector<int> CheckNodeMinDistNodeNum(ZoneNumPts, -1);

	m_IntValues = vector<double>(m_NumIntVars + static_cast<int>(IntegrateVolume), 0.0);

// 	if (IntegrateVolume)
// 		m_IntValues = vector<double>(m_NumIntVars + 1, 0.0);
// 	else
// 		m_IntValues = vector<double>(m_NumIntVars, 0.0);

	if (IsOk){
		CalcMaxNodeDistSqr();
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
		SaveAsTriFEZone(vector<int>({ 1, 2, 3 }));
// 		return IsOk;
#endif // _DEBUG

		vector<int> VolIJK(3), ZoneIJK(3);
		vec3 VolPt, NodePt, VolSubPt;
		vector<vec3> FarPoints(8), MinMaxXYZ = { m_VolZoneInfo.MinXYZ - 1, m_VolZoneInfo.MaxXYZ + 1 };
		for (int i = 0; i < 2; ++i){
			double x = MinMaxXYZ[i][0];
			for (int j = 0; j < 2; ++j){
				double y = MinMaxXYZ[j][1];
				for (int k = 0; k < 2; ++k){
					double z = MinMaxXYZ[k][2];
					FarPoints[4 * i + 2 * j + k] << x << y << z;
				}
			}
		}
		LgIndex_t VolIndex, ZoneIndex;
		ZoneIJK[2] = 1;
		for (VolIJK[2] = ZoneMinIJK[2]; VolIJK[2] <= ZoneMaxIJK[2]; ++VolIJK[2]){
			ZoneIJK[1] = 1;
			for (VolIJK[1] = ZoneMinIJK[1]; VolIJK[1] <= ZoneMaxIJK[1]; ++VolIJK[1]){
				ZoneIJK[0] = 1;
				for (VolIJK[0] = ZoneMinIJK[0]; VolIJK[0] <= ZoneMaxIJK[0]; ++VolIJK[0]){
// 					VolIndex = IndexFromIJK(VolIJK[0], VolIJK[1], VolIJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]);
					ZoneIndex = IndexFromIJK(ZoneIJK[0], ZoneIJK[1], ZoneIJK[2], ZoneNumIJK[0], ZoneNumIJK[1]) - 1;
					for (int ii = 0; ii < 3; ++ii){
						VolPt[ii] = m_VolZoneInfo.MinXYZ[ii] + static_cast<double>(VolIJK[ii]-1) * m_VolZoneInfo.PointSpacingV123[ii];
						NodeIsInterior[ZoneIndex] = NodeIsInterior[ZoneIndex] && (ZoneIJK[ii] > 1 && ZoneIJK[ii] < ZoneNumIJK[ii]);
					}
					if (ZoneIJK[0] == 11 && ZoneIJK[1] == 6 && ZoneIJK[2] == 6){
						int a = 1;
					}
					NodeIsInterior[ZoneIndex] = NodeIsInterior[ZoneIndex] && PointIsInterior(VolPt, FarPoints);
					if (!NodeIsInterior[ZoneIndex]){
						DoCheckNode[ZoneIndex] = DistSqrToSurfaceNodeWithinTolerance(VolPt, CheckNodeMinDistSqrToSurfaceNode[ZoneIndex], CheckNodeMinDistNodeNum[ZoneIndex]);
						/*for (int n = 0; n < m_NumNodes && !DoCheckNode[ZoneIndex]; ++n){
							for (int ii = 0; ii < 3; ++ii){
								NodePt[ii] = m_XYZPtrs[ii][n];
							}
							DoCheckNode[ZoneIndex] = (DistSqr(VolPt, NodePt) < m_MaxNeighborNodeDistSqr[n]);
							//DoCheckNode[ZoneIndex] = (DistSqr(VolPt, NodePt) < m_MaxNodeDistanceSqr);
						}*/
					}
					ZoneIJK[0]++;
				}
				ZoneIJK[1]++;
			}
			ZoneIJK[2]++;
		}


#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
		int NumCheckPoints = 0;
		for (auto const & i : DoCheckNode)
			NumCheckPoints += int(i);
		m_InteriorPts.reserve(NumCheckPoints);
		m_CheckPts.reserve(NumCheckPoints);
		m_InteriorSubPts.reserve(NumCheckPoints * (2 ^ ResolutionScale));
		m_CheckSubPts.reserve(NumCheckPoints * (2 ^ ResolutionScale));
#endif


		vec3 SubDelXYZ, TmpPoint, TmpSubPoint;

		vector<vec3> uvwVec(3);
		for (int i = 0; i < 3; ++i)
			uvwVec[i] = m_VolZoneInfo.BasisVectors.col(i) / (double)m_VolZoneInfo.MaxIJK[i];
		double CellVolume = ParallepipedVolume(uvwVec);

// 		double CellVolume = 1.0;

// 		for (int i = 0; i < 3; ++i){
// 			// 			PointSpacingV123[i] = (m_VolZoneInfo.MaxXYZ[i] - m_VolZoneInfo.MinXYZ[i]) / static_cast<double>(m_VolZoneInfo.MaxIJK[i] - 1);
// 			CellVolume *= m_VolZoneInfo.PointSpacingV123[i];
// 		}

		SubDelXYZ = m_VolZoneInfo.PointSpacingV123 / static_cast<double>(ResolutionScale);
		double SubDivideFactorUser = 1.0 / static_cast<double>(ResolutionScale * ResolutionScale * ResolutionScale);
		double SubCellVolumeUser = CellVolume * SubDivideFactorUser;
		double SubDivideFactor, SubCellVolume;
		vector<LgIndex_t> NeighborIJK(3);
		vector<Boolean_t> NodeIsInteriorCopy = NodeIsInterior;
		LgIndex_t CheckInt;
		ZoneIJK[2] = 1;

		for (VolIJK[2] = ZoneMinIJK[2]; VolIJK[2] <= ZoneMaxIJK[2]; ++VolIJK[2]){
		//for (VolIJK[2] = ZoneMinIJK[2]; VolIJK[2] <= ZoneMaxIJK[2] - (ZoneNumIJK[2] * 3 / 4); ++VolIJK[2]){
			ZoneIJK[1] = 1;
			for (VolIJK[1] = ZoneMinIJK[1]; VolIJK[1] <= ZoneMaxIJK[1]; ++VolIJK[1]){
				ZoneIJK[0] = 1;
				for (VolIJK[0] = ZoneMinIJK[0]; VolIJK[0] <= ZoneMaxIJK[0]; ++VolIJK[0]){
					ZoneIndex = IndexFromIJK(ZoneIJK[0], ZoneIJK[1], ZoneIJK[2], ZoneNumIJK[0], ZoneNumIJK[1]) - 1; 
					for (int ii = 0; ii < 3; ++ii){
						VolPt[ii] = m_VolZoneInfo.MinXYZ[ii] + static_cast<double>(VolIJK[ii] - 1) * m_VolZoneInfo.PointSpacingV123[ii];
					}
					if (ZoneIJK[0] == 11 && ZoneIJK[1] == 6 && ZoneIJK[2] == 6){
						int a = 1;
					}
					if (NodeIsInterior[ZoneIndex]){
						NeighborIJK = ZoneIJK;
						for (NeighborIJK[0] = ZoneIJK[0] - 1; NeighborIJK[0] < ZoneIJK[0] + 2; ++NeighborIJK[0]){
							for (NeighborIJK[1] = ZoneIJK[1] - 1; NeighborIJK[1] < ZoneIJK[1] + 2; ++NeighborIJK[1]){
								for (NeighborIJK[2] = ZoneIJK[2] - 1; NodeIsInterior[ZoneIndex] && NeighborIJK[2] < ZoneIJK[2] + 2; ++NeighborIJK[2]){
									NodeIsInterior[ZoneIndex] = NodeIsInteriorCopy[IndexFromIJK(NeighborIJK[0], NeighborIJK[1], NeighborIJK[2], ZoneNumIJK[0], ZoneNumIJK[1]) - 1];
								}
							}

						}
// 						for (int iDir = 0; iDir < 3; ++iDir){
// 							for (int d = -1; d < 2 && NodeIsInterior[ZoneIndex]; d += 2){
// 								CheckInt = NeighborIJK[iDir] + d;
// 								NodeIsInterior[ZoneIndex] = (CheckInt >= 1 && CheckInt <= ZoneNumIJK[iDir]);
// 								if (NodeIsInterior[ZoneIndex]){
// 									NeighborIJK[iDir] = CheckInt;
// 									NodeIsInterior[ZoneIndex] = NodeIsInteriorCopy[IndexFromIJK(NeighborIJK[0], NeighborIJK[1], NeighborIJK[2], ZoneNumIJK[0], ZoneNumIJK[1]) - 1];
// 									NeighborIJK[iDir] -= d;
// 								}
// 							}
// 						}
						DoCheckNode[ZoneIndex] = !NodeIsInterior[ZoneIndex];
						if (NodeIsInterior[ZoneIndex]){
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
							m_InteriorPts.push_back(VolPt);
#endif
							VolIndex = IndexFromIJK(VolIJK[0], VolIJK[1], VolIJK[2], m_VolZoneInfo.MaxIJK[0], m_VolZoneInfo.MaxIJK[1]) - 1;
							SubDivideFactor = 1.;
							for (int dir = 0; dir < 3; ++dir){
								if (VolIJK[dir] == 1 || VolIJK[dir] == m_VolZoneInfo.MaxIJK[dir])
									SubDivideFactor *= 0.5;
							}
							SubCellVolume = CellVolume * SubDivideFactor;
							if (IntegrateVolume)
								m_IntValues[m_NumIntVars] += SubCellVolume;
							for (int i = 0; i < m_NumIntVars; ++i){
								m_IntValues[i] += m_IntVarPtrs[i][VolIndex] * SubCellVolume;
							}
						}
						else{
							DistSqrToSurfaceNodeWithinTolerance(VolPt, CheckNodeMinDistSqrToSurfaceNode[ZoneIndex], CheckNodeMinDistNodeNum[ZoneIndex]);
						}
					}
					if (!NodeIsInterior[ZoneIndex] && DoCheckNode[ZoneIndex]){
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
						m_CheckPts.push_back(VolPt);
#endif
						SubDivideIntegrateCellAtPoint(VolPt, FarPoints, m_VolZoneInfo.PointSpacingV123 * 0.5, CheckNodeMinDistSqrToSurfaceNode[ZoneIndex], CheckNodeMinDistNodeNum[ZoneIndex], ResolutionScale - 1, IntegrateVolume);
// 						for (int ii = 0; ii < 3; ++ii){
// 							VolPt[ii] = m_VolZoneInfo.MinXYZ[ii] + static_cast<double>(VolIJK[ii]-1) * m_VolZoneInfo.PointSpacingV123[ii];
// 							VolPt[ii] = MAX(m_VolZoneInfo.MinXYZ[ii], MIN(m_VolZoneInfo.MaxXYZ[ii], VolPt[ii] - (m_VolZoneInfo.PointSpacingV123[ii] / 2.) + (SubDelXYZ[ii] / 2.)));
// 						}
// 						VolSubPt[0] = VolPt[0];
// 						for (int i = 0; i < ResolutionScale; ++i){
// 							VolSubPt[1] = VolPt[1];
// 							for (int j = 0; j < ResolutionScale; ++j){
// 								VolSubPt[2] = VolPt[2];
// 								for (int k = 0; k < ResolutionScale; ++k){
// #ifdef RECORD_INT_POINTS
// 									InteriorCheckPts.push_back(VolSubPt);
// #endif
// 									if (PointIsInterior(VolSubPt, FarPoint) && SetIndexAndWeightsForPoint(VolSubPt, m_VolZoneInfo)){
// #ifdef RECORD_INT_POINTS
// 										InteriorSubPts.push_back(VolSubPt);
// #endif
// 										if (IntegrateVolume)
// 											m_IntValues[m_NumIntVars] += SubCellVolumeUser;
// 										for (int i = 0; i < m_NumIntVars; ++i){
// 											m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * SubCellVolumeUser;
// 										}
// 									}
// 									VolSubPt[2] += SubDelXYZ[2];
// 								}
// 								VolSubPt[1] += SubDelXYZ[1];
// 							}
// 							VolSubPt[0] += SubDelXYZ[0];
// 						}
					}
					ZoneIJK[0]++;
				}
				ZoneIJK[1]++;
			}
			ZoneIJK[2]++;
		}
		m_IntegrationResultsReady = TRUE;
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
		vector<string> ZoneNames = { "Interior Points", "CheckPoints", "Interior Subpoints", "Sub Check Points" };
		for (auto & i : ZoneNames) i += " zone " + to_string(m_ZoneNum);
		vector<ColorIndex_t> Colors = { Green_C, Blue_C, Purple_C, Red_C };
		vector<vector<vec3> > Vecs = { m_InteriorPts, m_CheckPts, m_InteriorSubPts, m_CheckSubPts };
		vector<EntIndex_t> VarNums(3);
		for (int i = 0; i < 3; ++i)
			VarNums[i] = m_XYZPtrs[i].VarNum();
		for (int i = 0; i < Vecs.size(); ++i)
			SaveVec3VecAsScatterZone(Vecs[i], ZoneNames[i], Colors[i], VarNums);
#endif
	}

// 	m_ElemList = vector<vector<LgIndex_t> >();
// 	m_RefinedXYZList = vector<vec3>();
	
	return IsOk;
}



Boolean_t FESurface_c::CalcMaxNodeDistSqr(){
	Boolean_t IsOk = m_XYZPtrs.size() == 3;

	for (int i = 0; i < m_XYZPtrs.size() && IsOk; ++i){
		IsOk = m_XYZPtrs[i].IsReady();
	}
	
	m_XYZList.resize(m_NumNodes);
	for (int n = 0; n < m_NumNodes; ++n){
		for (int i = 0; i < 3; ++i){
			m_XYZList[n][i] = m_XYZPtrs[i][n];
		}
	}

// 	if (IsOk){
// 		m_MaxNodeDistanceSqr = 0.0;
// 		m_MinNodeDistanceSqr = 1e150;
// 		m_MaxNeighborNodeDistSqr.resize(m_NumNodes, 0.0);
// 		m_MinNeighborNodeDistSqr.resize(m_NumNodes, 1e100);
// 		double TempDist;
// 		for (int i = 0; i < m_NumNodes; ++i){
// 			for (auto const & j : m_ConnectivityList[i]){
// 				TempDist = DistSqr(m_XYZList[i], m_XYZList[j]);
// 				m_MaxNeighborNodeDistSqr[i] = MAX(m_MaxNeighborNodeDistSqr[i], TempDist);
// 				if (TempDist > 1e-5){
// 					m_MinNeighborNodeDistSqr[i] = MIN(m_MinNeighborNodeDistSqr[i], TempDist);
// 				}
// 			}
// 			m_MaxNodeDistanceSqr = MAX(m_MaxNodeDistanceSqr, m_MaxNeighborNodeDistSqr[i]);
// 			m_MinNodeDistanceSqr = MIN(m_MinNodeDistanceSqr, m_MinNeighborNodeDistSqr[i]);
// 		}
// 	}
// 	
	m_MinNodeDistanceSqr = sum(square(m_VolZoneInfo.PointSpacingV123));
	
	m_RefinedXYZList = m_XYZList;
	m_ElemList.resize(m_NumElems, vector<LgIndex_t>(3, -1));
	
	for (int i = 0; i < m_NumElems; ++i){
		for (int j = 0; j < 3; ++j){
			m_ElemList[i][j] = m_ConnectivityListPtr[i * m_NumNodesPerElem + j];
		}
	}

	vector<int> ElemNumList(m_NumElems);
	for (int i = 0; i < m_NumElems; ++i){
		ElemNumList[i] = i;
	}
	RefineTriElems(ElemNumList);

	m_MinNodeDistanceSqr *= 0.49; // 70% of distance (70%^2 of square distance)

	return IsOk;
}

// vector<int> FESurface_c::TriangleEdgeMidPointSubdivide(int TriNum){
// 	vector<LgIndex_t> NewNodes(3);
// 	vector<int> NewTris(4);
// 	NewTris[0] = TriNum;
// 	int ni = m_RefinedXYZList.size();
// 	int ti = m_ElemList.size();
// 	for (int i = 0; i < 3; ++i){
// 		m_RefinedXYZList.push_back(
// 			(
// 				m_RefinedXYZList[m_ElemList[TriNum][i]]
// 				+ m_RefinedXYZList[m_ElemList[TriNum][(i + 1) % 3]]
// 			) / 2.
// 		);
// 		NewNodes[i] = ni++;
// 		NewTris[i+1] = ti++;
// 	}
// 	m_ElemList.push_back(vector<LgIndex_t>({m_ElemList[TriNum][1], NewNodes[1], NewNodes[0]}));
// 	m_ElemList.push_back(vector<LgIndex_t>({m_ElemList[TriNum][2], NewNodes[2], NewNodes[1]}));
// 	m_ElemList.push_back(NewNodes);
// 	m_ElemList[TriNum] = vector<LgIndex_t>({m_ElemList[TriNum][0], NewNodes[0], NewNodes[2]});
// 	return NewTris;
// }

void FESurface_c::RefineTriElems(vector<int> const & TriNumList){
	for (auto const & t : TriNumList){
		double MaxNodeDistSqr = 0.0;
		for (int i = 0; i < 3; ++i){
			MaxNodeDistSqr = MAX(MaxNodeDistSqr, DistSqr(m_RefinedXYZList[m_ElemList[t][i]], m_RefinedXYZList[m_ElemList[t][(i+1) % 3]]));
		}
		if (MaxNodeDistSqr > m_MinNodeDistanceSqr){
			vector<int> refinedTriNums;
			TriangleEdgeMidPointSubdivide(m_RefinedXYZList, m_ElemList, t, refinedTriNums);
			RefineTriElems(refinedTriNums);
		}
	}
}


// void FESurface_c::TriPolyLines(bool const ConnectBeginningAndEndGPs)
// {
// 	vec3 NewNode, EdgeVec;
// 	for (int iGP = 0; iGP < m_NumGPs - (ConnectBeginningAndEndGPs ? 0 : 1); ++iGP){
// 		int lInd = iGP,
// 			rInd = (iGP + 1) % m_NumGPs;
// 		int li = lInd * m_NumGPPts,
// 			ri = rInd * m_NumGPPts;
// 		vector<int> lr = { li, ri }, 
// 			lrMid = { li,ri };
// 		vector<bool> MidFound(2, false);
// 		int NumEdgePts;
// 		double EdgePtSpacing, EdgeLen, MinNodeScore, TmpNodeScore, MinNodeNum, lLen, rLen, TmpLen;
// 		EdgeVec = m_XYZList[ri + m_NumGPPts - 1] - m_XYZList[li + m_NumGPPts - 1];
// 		EdgeLen = norm(EdgeVec);
// 		EdgePtSpacing = norm(m_XYZList[li] - m_XYZList[li + 1]);
// 		NumEdgePts = static_cast<int>(EdgeLen / EdgePtSpacing);
// 		bool HasFarEdge = (NumEdgePts > 2);
// 		if (HasFarEdge){
// 			EdgePtSpacing = EdgeLen / static_cast<double>(NumEdgePts - 1);
// 			EdgeVec = normalise(EdgeVec) * EdgePtSpacing;
// 			MinNodeScore = DBL_MAX;
// 			for (int i = 1; i < NumEdgePts - 1; ++i){
// 				NewNode = m_XYZList[li + m_NumGPPts - 1] + EdgeVec * static_cast<double>(i);
// 				lLen = DBL_MAX;
// 				rLen = DBL_MAX;
// 				for (int j = 0; j < m_NumGPPts; ++j){
// 					lLen = MIN(lLen, DistSqr(NewNode, m_XYZList[li + j]));
// 					rLen = MIN(rLen, DistSqr(NewNode, m_XYZList[ri + j]));
// 
// 				}
// 				TmpNodeScore = lLen + rLen;
// 				if (TmpNodeScore < MinNodeScore){
// 					MinNodeScore = TmpNodeScore;
// 					MinNodeNum = i;
// 				}
// 			}
// 			NewNode = m_XYZList[li + m_NumGPPts - 1] + EdgeVec * static_cast<double>(MinNodeNum);
// 		}
// 		else
// 			NewNode = (m_XYZList[li + m_NumGPPts - 1] + m_XYZList[ri + m_NumGPPts - 1]) / 2.0;
// 		if (HasFarEdge) m_XYZList.push_back(NewNode);
// 		int NewNodeNum = m_XYZList.size() - 1;
// 		while (lr[0] < (lInd + 1) * m_NumGPPts - 1 && lr[1] < (rInd + 1) * m_NumGPPts - 1){
// 			int FarPoint, MinSide, MinFarPoint;
// 			double MinLen = DBL_MAX, TmpLen;
// 			for (int i = 0; i < 2; ++i){
// 				for (int j = 0; j < 1 + int(HasFarEdge); ++j){
// 					if (j == 0)
// 						FarPoint = lr[(i + 1) % 2];
// 					else
// 						FarPoint = NewNodeNum;
// 					TmpLen = DistSqr(m_XYZList[lr[i] + 1], m_XYZList[FarPoint]);
// 					if (TmpLen < MinLen){
// 						MinSide = i;
// 						MinFarPoint = FarPoint;
// 						MinLen = TmpLen;
// 					}
// 				}
// 			}
// 			if (HasFarEdge && !MidFound[MinSide] && MinFarPoint == NewNodeNum){
// 				MidFound[MinSide] = true;
// 				lrMid[MinSide] = lr[MinSide];
// 			}
// 			m_ElemList.push_back({ lr[MinSide], lr[MinSide] + 1, MinFarPoint });
// 			lr[MinSide]++;
// 		}
// 		vector<int> lrInd = { lInd, rInd };
// 		for (int i = 0; i < 2; ++i){
// 			while (lr[i] < (lrInd[i] + 1) * m_NumGPPts - 1){
// 				if (DistSqr(m_XYZList[lr[i] + 1], m_XYZList[lr[(i + 1) % 2]]) < DistSqr(m_XYZList[lr[i] + 1], m_XYZList[NewNodeNum])){
// 					m_ElemList.push_back({ lr[i], lr[i] + 1, lr[(i + 1) % 2] });
// 				}
// 				else if (HasFarEdge){
// 					m_ElemList.push_back({ lr[i], lr[i] + 1, NewNodeNum });
// 					if (!MidFound[i]){
// 						MidFound[i] = true;
// 						lrMid[i] = lr[i];
// 					}
// 				}
// 				lr[i]++;
// 			}
// 		}
// 		if (MidFound[0] && MidFound[1]){
// 			m_ElemList.push_back({ lrMid[0], lrMid[1], NewNodeNum });
// 		}
// 	}
// }

Boolean_t FESurface_c::DistSqrToSurfaceNodeWithinTolerance(vec3 const & CheckPt,
																double & NewDistSqrUnderTol,
																int & CloseNodeNum,
																double const & DistSqrTol)
{
	Boolean_t IsFound = FALSE;
	// first check the old closest node to see if it's under tolerance
	if (CloseNodeNum >= 0){
		NewDistSqrUnderTol = DistSqr(CheckPt, m_RefinedXYZList[CloseNodeNum]);
		if (DistSqrTol > 0)
			IsFound = (NewDistSqrUnderTol <= DistSqrTol);
		else
			IsFound = (NewDistSqrUnderTol <= m_MinNodeDistanceSqr);
	}
	for (int CloseNodeNum = 0; CloseNodeNum < m_RefinedXYZList.size() && !IsFound; ++CloseNodeNum){
		NewDistSqrUnderTol = DistSqr(CheckPt, m_RefinedXYZList[CloseNodeNum]);
		if (DistSqrTol > 0)
			IsFound = (NewDistSqrUnderTol <= DistSqrTol);
		else
			IsFound = (NewDistSqrUnderTol <= m_MinNodeDistanceSqr);
		//DoCheckNode[ZoneIndex] = (DistSqr(VolPt, NodePt) < m_MaxNodeDistanceSqr);
	}
	return IsFound;
}

// Boolean_t FEVolume_c::DistSqrToSurfaceEdgeWithinTolerance(vec3 const & CheckPt,
// 	double & NewDistSqrUnderTol,
// 	double const & DistSqrTol)
// {
// 	Boolean_t IsFound = FALSE;
// 	vec3 NodePt;
// 
// 	for (int n = 0; n < m_NumNodes && !IsFound; ++n){
// 
// 	}
// 
// 	return IsFound;
// }

Boolean_t FESurface_c::SubDivideIntegrateCellAtPoint(vec3 const & Point,
														vector<vec3> const & FarPoints,
														vec3 const & DelXYZ,
														double const & MinDistSqrToSurfaceNode,
														int MinDistNodeNum,
														int SubDivideLevel,
														Boolean_t IntegrateVolume)
{
	Boolean_t IsOk = TRUE;
	vec3 TmpPoint;
	double TmpDistSqr;

	for (double xi = -1.0; xi < 1.1; xi += 2.0){
		TmpPoint[0] = Point[0] + DelXYZ[0] * xi * 0.5;
		for (double yi = -1.0; yi < 1.1; yi += 2.0){
			TmpPoint[1] = Point[1] + DelXYZ[1] * yi * 0.5;
			for (double zi = -1.0; zi < 1.1; zi += 2.0){
				TmpPoint[2] = Point[2] + DelXYZ[2] * zi * 0.5;
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
				m_CheckSubPts.push_back(TmpPoint);
				if (m_InteriorSubPts.size() == 3966){
					int a = 1;
				}
#endif // _DEBUG

				if (SubDivideLevel > 0){
					if (DistSqrToSurfaceNodeWithinTolerance(TmpPoint, TmpDistSqr, MinDistNodeNum, MinDistSqrToSurfaceNode)){
						IsOk = SubDivideIntegrateCellAtPoint(TmpPoint, FarPoints, DelXYZ * 0.5, TmpDistSqr, MinDistNodeNum, SubDivideLevel - 1, IntegrateVolume);
					}
					else if (PointIsInterior(TmpPoint, FarPoints) && SetIndexAndWeightsForPoint(TmpPoint, m_VolZoneInfo)){
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
						m_InteriorSubPts.push_back(TmpPoint);

#endif // _DEBUG

						double CellVolume = 1;
						for (int dir = 0; dir < 3; ++dir)
							CellVolume *= DelXYZ[dir];
						if (IntegrateVolume)
							m_IntValues[m_NumIntVars] += CellVolume;
						for (int i = 0; i < m_NumIntVars; ++i){
							m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * CellVolume;
						}
					}
				}
				else if (PointIsInterior(TmpPoint, FarPoints) && SetIndexAndWeightsForPoint(TmpPoint, m_VolZoneInfo)){
#if defined(RECORD_INT_POINTS) && defined(_DEBUG)
					m_InteriorSubPts.push_back(TmpPoint);
#endif // _DEBUG
					double CellVolume = 1;
					for (int dir = 0; dir < 3; ++dir)
						CellVolume *= DelXYZ[dir];
					if (IntegrateVolume)
						m_IntValues[m_NumIntVars] += CellVolume;
					for (int i = 0; i < m_NumIntVars; ++i){
						m_IntValues[i] += ValByCurrentIndexAndWeightsFromRawPtr(m_VolZoneInfo, m_IntVarPtrs[i]) * CellVolume;
					}
				}
			}
		}
	}

	return IsOk;
}


/*
	*	Begin implementation of Charles' integration code
	*/


vector<double> split(string const &s, char delim) {
	stringstream ss(s);
	string item;
	vector<double> tokens;
	while (getline(ss, item, delim)) {
		tokens.push_back(std::stod(item));
	}
	return tokens;
}

mat LoadFile(string const & Path, char Delim){
	mat Points(1, 3);
	std::ifstream File(Path.c_str());

	if (File.is_open()){
		string Line;
		while (!File.eof()){
			std::getline(File, Line);
			vector<double> tmp = split(Line, Delim);
			if (tmp.size() > 0)
				Points = join_cols(Points, mat(tmp).t());
		}
		File.close();
	}
	else{
		TecUtilDialogErrMsg("Failed to open file");
	}

	return Points.tail_rows(Points.n_rows - 1);
}



double FESurface_c::IntVolume(int N, vec3 const & StartPoint) const{
// 	vector<int> NumElemsList = { 80, 320, 1280, 5120 };

	mat Verts;
	Verts << 0 << 0 << 0 << endr
		<< 1 << 0 << 0 << endr
		<< 0 << 1 << 0 << endr
		<< 0 << 0 << 1;
	// 		mat Verts = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

	vector<vec> stuW = tetraquad(N, Verts);

	int NumW = stuW[3].n_elem;

	vector<mat> Points(3, mat(m_GPList[0].GetCount() + 1, 3));
	for (int i = 0; i < 3; ++i){
		Points[i].row(0) = StartPoint.t();
		for (int j = 0; j < m_GPList[i].GetCount(); ++j){
			Points[i].row(j + 1) = m_GPList[i][j].t();
		}
	}

	mat abc = cubTrans(Points[0], Points[1], Points[2]);
	mat xyz = cubTrans_func(abc, stuW[0], stuW[1], stuW[2]);
	vec WTemp(NumW);

	for (int i = 0; i < NumW; ++i){
		WTemp(i) = stuW[3](i) * cubJacobian(abc, stuW[0](i), stuW[1](i), stuW[2](i));
	}

	return sum(WTemp);

	return 0;
}

vector<mat> FESurface_c::GetIntegrationPointsWeights(vec3 const & StartPoint, vector<vec> const & stuW) const{

	if (stuW.size() == 4){

		int NumW = stuW[3].n_elem;

		vector<mat> Points(3, mat(m_GPList[0].GetCount() + 1, 3));
		for (int i = 0; i < 3; ++i){
			Points[i].row(0) = StartPoint.t();
			for (int j = 0; j < m_GPList[i].GetCount(); ++j){
				Points[i].row(j + 1) = m_GPList[i][j].t();
			}
		}

		mat abc = cubTrans(Points[0], Points[1], Points[2]);
		mat xyz = cubTrans_func(abc, stuW[0], stuW[1], stuW[2]);
		vec WTemp(NumW);

		for (int i = 0; i < NumW; ++i){
			WTemp(i) = stuW[3](i) * cubJacobian(abc, stuW[0](i), stuW[1](i), stuW[2](i));
		}
		return vector<mat>({ xyz, WTemp });
	}
	return vector<mat>();
}

vector<mat> FESurface_c::GetIntegrationPointsWeights(vector<vec> const & stuW) const{

	if (stuW.size() == 4){

		int NumW = stuW[3].n_elem;

		vector<mat> Points(3, mat(m_GPList[0].GetCount(), 3));
		for (int i = 0; i < 3; ++i){
			for (int j = 0; j < m_GPList[i].GetCount(); ++j){
				Points[i].row(j) = m_GPList[i][j].t();
			}
		}

		mat abc = cubTrans(Points[0], Points[1], Points[2]);
		mat xyz = cubTrans_func(abc, stuW[0], stuW[1], stuW[2]);
		vec WTemp(NumW);

		for (int i = 0; i < NumW; ++i){
			WTemp(i) = stuW[3](i) * cubJacobian(abc, stuW[0](i), stuW[1](i), stuW[2](i));
		}
		return vector<mat>({ xyz, WTemp });
	}
	return vector<mat>();
}




/*
*	Private member functions
*/

Boolean_t FESurface_c::PointIsInterior(vec3 const & Point, vector<vec3> const & FarPoints) const{
	/*
	*	Check that Point is in bounding cube of FEVolume
	*/
// 	Boolean_t IsInterior = (sum(Point >= m_ZoneMinXYZ) == 3 && sum(Point <= m_ZoneMaxXYZ) == 3);
	Boolean_t IsInterior;
	/*
	*	Now loop over each element, checking the triangles in that element
	*	for intersections between the triangle
	*	and a line segment between the Point and the FarPoints.
	*/
	int NumIntersections;
	bool FarPointGood;
	for (vec3 const & FarPoint : FarPoints){
		NumIntersections = 0;
		FarPointGood = true;
		vec3 T[3];
		for (int ElemNum = 0; ElemNum < m_NumElems && FarPointGood; ++ElemNum){
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					T[j][i] = m_XYZPtrs[i][m_ConnectivityListPtr[ElemNum * m_NumNodesPerElem + j]];
				}
			}
			switch (TriangleIntersect(T[0], T[1], T[2], FarPoint, Point)){
				case 1:
					NumIntersections++;
					break;
				case 2:
					FarPointGood = false;
					break;
			}
		}
		if (FarPointGood) break;
	}

	IsInterior = FarPointGood && (NumIntersections % 2 != 0);

	return IsInterior;
}
int FESurface_c::TriangleIntersect(vec3 const & T_P0,
	vec3 const & T_P1,
	vec3 const & T_P2,
	vec3 const & R_P1,
	vec3 const & R_P0) const
{
	if (sum(T_P0 == T_P1) == 3 || sum(T_P0 == T_P2) == 3 || sum(T_P1 == T_P2) == 3)
		return 0;

	vec3 u = T_P1 - T_P0;
	vec3 v = T_P2 - T_P0;
	vec3 n = cross(u, v);
	vec3 w0 = R_P0 - T_P0;
	vec3 D = R_P1 - R_P0;

	//TODO: check for degeneracy 

	//Check for ray-plane parallelism
	double a = abs(dot(n, D));
	if (abs(dot(n, D)) < 1e-10){
		return 0;
	}

	//Check if ray goes towards triangle plane
	double r = -dot(n, w0) / dot(n, D);
	if (r < 0)
		return 0;

	//Calc intersect or R and the plane
	vec3 P = R_P0 + (D * r);

	vec3 w = P - T_P0;
	double uv = dot(u, v);
	double wv = dot(w, v);
	double vv = dot(v, v);
	double wu = dot(w, u);
	double uu = dot(u, u);

	double Denom = ((uv * uv) - (uu * vv));

	double s = ((uv*wv) - (vv*wu)) / Denom;
	double t = ((uv*wu) - (uu*wv)) / Denom;

	//TODO : boundary cases
	if (s >= 0 && t >= 0 && s + t <= 1)
		return 1;
	else
		return 0;
}

vector<double> FESurface_c::TriSphereElemSolidAngles() const
{
	Boolean_t IsOk = TRUE;

	vector<double> SolidAngles;

	if (IsOk) {
		SolidAngles.resize(m_NumElems);
		double TotalArea = 0.0;
		/*
		*	Get total surface area of sphere, storing element areas
		*	as they're found.
		*/

// 		vec3 T[3];
// 		double A, B, C;

		for (int ElemNum = 0; ElemNum < m_NumElems; ++ElemNum) {
			/*
			*	Get the corners of the triangle element.
			*/
			
// 			for (int i = 0; i < 3; ++i) {
// 				T[i] = m_XYZList[m_ElemList[ElemNum][i]];
// 			}
// 
// 			/*
// 			*	Get the area of the element
// 			*/
// 			A = T[0][2] * (T[1][1] - T[2][1])
// 				+ T[1][2] * (T[2][1] - T[0][1])
// 				+ T[2][2] * (T[0][1] - T[1][1]);
// 
// 			B = T[0][0] * (T[1][2] - T[2][2])
// 				+ T[1][0] * (T[2][2] - T[0][2])
// 				+ T[2][0] * (T[0][2] - T[1][2]);
// 
// 			C = T[0][1] * (T[1][0] - T[2][0])
// 				+ T[1][1] * (T[2][0] - T[0][0])
// 				+ T[2][1] * (T[0][0] - T[1][0]);

// 			SolidAngles[ElemNum] = 0.5 * sqrt(A * A + B * B + C * C);
			SolidAngles[ElemNum] = TriArea(m_XYZList[m_ElemList[ElemNum][0]], m_XYZList[m_ElemList[ElemNum][1]], m_XYZList[m_ElemList[ElemNum][2]]);
			TotalArea += SolidAngles[ElemNum];
		}

#pragma omp parallel for
		for (int ElemNum = 0; ElemNum < m_NumElems; ++ElemNum) {
			SolidAngles[ElemNum] /= TotalArea;
		}
	}

	return SolidAngles;
}


vector<vector<double> > FESurface_c::TriSphereIntValsByElem(vector<double> * SphereTriangleAreas) const
{
	Boolean_t IsOk = m_IntegrationResultsReady && m_NumNodesPerElem == 3;

	vector<vector<double> > IntVals;

	if (IsOk){
		double TotalArea = 0.0;
		IntVals.resize(m_NumElems, vector<double>(m_IntValues.size()));
		/*
		*	Get total surface area of sphere, storing element areas
		*	as they're found.
		*/

		vector<double> AreaFactors(m_NumElems);
		vec3 T[3];
		double A, B, C;

		for (int ElemNum = 0; ElemNum < m_NumElems; ++ElemNum){
			/*
			*	Get the corners of the triangle element.
			*/
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					T[j][i] = m_XYZPtrs[i][m_ConnectivityListPtr[ElemNum * m_NumNodesPerElem + j]];
				}
			}

			/*
			*	Get the area of the element
			*/
			A = T[0][2]*(T[1][1] - T[2][1])
				+ T[1][2]*(T[2][1] - T[0][1])
				+ T[2][2]*(T[0][1] - T[1][1]);

			B = T[0][0]*(T[1][2] - T[2][2])
				+ T[1][0]*(T[2][2] - T[0][2])
				+ T[2][0]*(T[0][2] - T[1][2]);

			C = T[0][1]*(T[1][0] - T[2][0])
				+ T[1][1]*(T[2][0] - T[0][0])
				+ T[2][1]*(T[0][0] - T[1][0]);

			AreaFactors[ElemNum] = 0.5 * sqrt(A*A + B*B + C*C);
			TotalArea += AreaFactors[ElemNum];
		}

		if (SphereTriangleAreas != nullptr)
			*SphereTriangleAreas = AreaFactors;

		for (int ElemNum = 0; ElemNum < m_NumElems; ++ElemNum){
			AreaFactors[ElemNum] /= TotalArea;
			for (int VarNum = 0; VarNum < m_IntValues.size(); ++VarNum)
				IntVals[ElemNum][VarNum] = m_IntValues[VarNum] * AreaFactors[ElemNum];
		}
	}

	return IntVals;
}


vector<vec3> FESurface_c::GetSphereIntersectionPath(vec3 const & SphereCenter, double const & SphereRadius) {
	REQUIRE(SphereRadius > 0.0);
	REQUIRE(IsMade());
	REQUIRE(m_ElemList.size() > 0);

	/*
	 *	See if we'll be using the list of nodes or pointers to node locations
	 */
	bool UseXYZList = m_XYZList.size() > 0;
	FieldVecPointer_c XYZVecPtr;
	if (!UseXYZList){
		XYZVecPtr = FieldVecPointer_c(m_XYZPtrs);
		REQUIRE(XYZVecPtr.Size() == m_NumNodes);
	}

	/*
	 *	Get list of edges, stored as "n1,n2" where n1 < n2.
	 */
	std::map<Edge, std::set<int> > EdgeToElemMap;
	int elemNum = 0;
	for (int ti = 0; ti < m_ElemList.size(); ++ti){
		for (int ci = 0; ci < m_ElemList[ti].size(); ++ci){
			EdgeToElemMap[MakeEdge(m_ElemList[ti][ci], m_ElemList[ti][(ci + 1) % m_ElemList[ti].size()])].insert(ti);
		}
	}

	double RadSqr = SphereRadius * SphereRadius;
	bool InsideNodeFound = false, OutsideNodeFound = false;

	/*
	 *	Make sure the surface itself intersects the sphere by finding 
	 *	at least one node inside and outside the sphere.
	 */
	if (UseXYZList){
		for (auto const & node : m_XYZList){
			if (!InsideNodeFound && RadSqr > DistSqr(SphereCenter, node))
				InsideNodeFound = true;
			if (!OutsideNodeFound && RadSqr < DistSqr(SphereCenter, node))
				OutsideNodeFound = true;
			if (InsideNodeFound && OutsideNodeFound)
				break;
		}
	}
	else{
		for (int i = 0; i < m_NumNodes; ++i) {
			if (!InsideNodeFound && RadSqr > DistSqr(SphereCenter, XYZVecPtr[i]))
				InsideNodeFound = true;
			if (!OutsideNodeFound && RadSqr < DistSqr(SphereCenter, XYZVecPtr[i]))
				OutsideNodeFound = true;
			if (InsideNodeFound && OutsideNodeFound)
				break;
		}
	}

	/*
	 *	Now we can find the intersection points
	 */
	vector<vec3> IntPoints;
	if (InsideNodeFound && OutsideNodeFound) {
		// Remove duplicate nodes
// 		this->RemoveDupicateNodes();

		/*
		 *	There could be artifact edges that are huge compared to the rest,
		 *	so we'll check for that by comparing edge length to the edge
		 *	length standard deviation.
		 */
		vec EdgeLengths(EdgeToElemMap.size());
		double MaxLen = DBL_MIN;
		double MinLen = DBL_MAX;
		int ei = 0;
		for (auto const & edge : EdgeToElemMap){
			if (UseXYZList)
				EdgeLengths[ei] = Distance(m_XYZList[edge.first.first], m_XYZList[edge.first.second]);
			else
				EdgeLengths[ei] = Distance(XYZVecPtr[edge.first.first], XYZVecPtr[edge.first.second]);
			ei++;
		}

		double EdgeLenMean = mean(EdgeLengths),
			EdgeLenStdDev = stddev(EdgeLengths);

		EdgeLengths = abs(EdgeLengths - EdgeLenMean);

		/*
		 *	Now loop over the edges to check if they straddle the sphere radius.
		 *	When an intersecting edge is found, confirm that it's not an outlier in terms
		 *	of edge length (use stddevCheckFactor * stddev as the cutoff), then add it to the map of intersecting edges.
		 */
		std::map<Edge, std::set<int> > IntEdges;
		std::map<Edge, vec3> IntEdgesPoints;

		double StdDevCheckFactor = 100000;

		int edgeNum = 0;
		for (auto const & edge : EdgeToElemMap){
			if (EdgeLengths[edgeNum]< StdDevCheckFactor * EdgeLenStdDev) {
				std::pair<double, double> dist;
				if (UseXYZList) {
					dist.first = DistSqr(m_XYZList[edge.first.first], SphereCenter);
					dist.second = DistSqr(m_XYZList[edge.first.second], SphereCenter);
				}
				else {
					dist.first = DistSqr(XYZVecPtr[edge.first.first], SphereCenter);
					dist.second = DistSqr(XYZVecPtr[edge.first.second], SphereCenter);
				}
				if (abs(dist.second - dist.first) > 1e-16) {

					int closeNode = 0;
// 					if (dist.first > dist.second) {
// 						closeNode = 1;
// 						dist = std::make_pair(dist.second, dist.first);
// 					}

					if ((dist.first < RadSqr && dist.second > RadSqr)
						|| (dist.second < RadSqr && dist.first > RadSqr)) {
						vec3 v1, v2;
// 						if (UseXYZList) {
// 							if (closeNode = 0) {
								v1 = m_XYZList[edge.first.first];
								v2 = m_XYZList[edge.first.second];
// 							}
// 							else {
// 								v1 = m_XYZList[edge.first.second];
// 								v2 = m_XYZList[edge.first.first];
// 							}
// 						}
// 						else {
// 							if (closeNode = 0) {
// 								v1 = XYZVecPtr[edge.first.first];
// 								v2 = XYZVecPtr[edge.first.second];
// 							}
// 							else {
// 								v1 = XYZVecPtr[edge.first.second];
// 								v2 = XYZVecPtr[edge.first.first];
// 							}
// 						}
						IntEdges.insert(edge);
						dist.first = sqrt(dist.first);
						dist.second = sqrt(dist.second);
						IntEdgesPoints[edge.first] = (v1 + (v2 - v1) * (SphereRadius - dist.first) / (dist.second - dist.first));
					}
				}
			}

			edgeNum++;
		}

		// New method:
		// Assume the end point (first point) of the intersection path is the point farthest from
		// the midpoint of all the intersection points.
		// The second point is then the point closest to the first point.
		// The next point is then the closest point to the last point that hasn't
		// already been added to the path, until all points are added.
		vector<vec3> InitialIntPoints;
		InitialIntPoints.reserve(IntEdgesPoints.size());
		for (auto const & p : IntEdgesPoints) {
			InitialIntPoints.push_back(p.second);
		}
// 		RemoveDuplicatePointsFromVec3Vec(InitialIntPoints, EdgeLenMean * 0.1); // do duplicate search during path formation

		IntPoints.reserve(InitialIntPoints.size());
// 		vector<bool> PointUsed(InitialIntPoints.size(), false);
		vec3 MidPt = zeros(3);
		for (auto const & p : InitialIntPoints){
			MidPt += p;
		}
		MidPt /= (double)InitialIntPoints.size();

		double MaxDist = DBL_MIN;
		int MaxPtNum = -1;
		vector<vec3>::iterator PtIt = InitialIntPoints.end();
		for (auto it = InitialIntPoints.begin(); it != InitialIntPoints.end(); ++it) {
			double TmpDist = DistSqr(*it, MidPt);
			if (TmpDist > MaxDist) {
				MaxDist = TmpDist;
				PtIt = it;
			}
		}

		if (PtIt != InitialIntPoints.end()){
			auto CheckPoints = InitialIntPoints;

			IntPoints.push_back(*PtIt);

			InitialIntPoints.erase(PtIt);

			double MinDist;
			double CheckDist = 0.0005;
			CheckDist *= CheckDist;

			while (!InitialIntPoints.empty()){
				MinDist = DBL_MAX;
				for (auto it = InitialIntPoints.begin(); it != InitialIntPoints.end(); ++it){
					double TmpDist = DistSqr(*it, IntPoints.back());
					if (TmpDist < MinDist){
						MinDist = TmpDist;
						PtIt = it;
					}
				}
				if (MinDist > 0 && (InitialIntPoints.size() == 1 || MinDist > CheckDist)){
					IntPoints.push_back(*PtIt);
				}
				InitialIntPoints.erase(PtIt);
			}
		}

// 		for (int ni = 0; ni < InitialIntPoints.size(); ++ni) {
// 			double TmpDist = DistSqr(InitialIntPoints[ni], MidPt);
// 			if (TmpDist > MaxDist){
// 				MaxDist = TmpDist;
// 				MaxPtNum = ni;
// 			}
// 		}
// 
// 		if (MaxPtNum >= 0) {
// 			IntPoints.push_back(InitialIntPoints[MaxPtNum]);
// // 
//  			double CheckDist = EdgeLenMean * 0.05;
//  			CheckDist *= CheckDist;
// 
// // 			PointUsed[MaxPtNum] = true;
// 
// 			auto CheckPoints = InitialIntPoints;
// 
// 			while (IntPoints.size() < InitialIntPoints.size()) {
// 				double MinDist = DBL_MAX;
// 				int MinPtNum = -1;
// 				for (int ni = 0; ni < InitialIntPoints.size(); ++ni) {
// 					if (!PointUsed[ni]) {
// 						double TmpDist = DistSqr(IntPoints.back(), InitialIntPoints[ni]);
// 						if (TmpDist < MinDist) {
// 							MinDist = TmpDist;
// 							MinPtNum = ni;
// 						}
// 					}
// 				}
// 				if (MinPtNum >= 0 && (MinDist > CheckDist || InitialIntPoints.size())) {
// 					IntPoints.push_back(InitialIntPoints[MinPtNum]);
// 					PointUsed[MinPtNum] = true;
// 				}
// 			}
// 		}


		
		// Old method based on shared elements of edges.
// 		/*
// 		 *	Now we have all the valid intersecting edges, so we need to get the
// 		 *	intersection points themselves, and in the correct order.
// 		 *	We can get a correct order by identifying the endpoint of the 
// 		 *	intersection path.
// 		 *	For an edge that is interior to the intersection path, both of its
// 		 *	participating elements will be in the element set of other intersecting
// 		 *	edges. So if both elements can be found in the element sets of the other
// 		 *	intersecting edges, then it is an interior edge. Conversely, if only one
// 		 *	of an intersecting edge's elements can be found in the element sets of
// 		 *	the other intersecting edges, then it is an endpoint edge and can be used
// 		 *	as the first intersection point.
// 		 *	The remaining points can be deduced by shared element numbers.
// 		 */
// 		IntPoints.reserve(IntEdges.size());
// 		vector<bool> EdgeUsed(IntEdges.size(), false);
// 		edgeNum = 0;
// 		std::set<int> elemSet;
// 
// 		/*
// 		 *	First find an endpoint.
// 		 */
// 		for (auto const & edge : IntEdges){
// 			/*
// 			 * Search the other edges' elements for the elements of this
// 			 * edge. If the number of hits is 1 then it's an endpoint.
// 			 */
// 			int NumHits = 0;
// 			int otherEdgeNum = 0;
// 			for (auto const & elem : edge.second){
// 				for (auto const & otherEdge : IntEdges){
// 					if (edgeNum != otherEdgeNum++ && otherEdge.second.count(elem)){
// 						NumHits++;
// 						break;
// 					}
// 				}
// 				if (NumHits >= 2) break;
// 			}
// 
// 			if (NumHits == 1 || edge.second.size() == 1){
// 				IntPoints.push_back(IntEdgesPoints[edge.first]);
// 				elemSet = edge.second;
// 				EdgeUsed[edgeNum] = true;
// 				break;
// 			}
// 
// 			edgeNum++;
// 		}
// 
// 		/*
// 		 * Now loop over the remaining edges until all the intersection points are found.
// 		 */
// 		int iter = 0;
// 		while (IntPoints.size() < IntEdges.size() && iter < 1.5 * (double)IntEdges.size()){
// 			iter++;
// 
// 			edgeNum = 0;
// 			for (auto const & edge : IntEdges) {
// 				if (!EdgeUsed[edgeNum]) {
// 					for (auto const & elem : elemSet) {
// 						if (edge.second.count(elem)){
// 							IntPoints.push_back(IntEdgesPoints[edge.first]);
// 							EdgeUsed[edgeNum] = true;
// 							elemSet = edge.second;
// 							break;
// 						}
// 					}
// 				}
// 				edgeNum++;
// 			}
// 		}
// 
// 		if (iter > IntEdges.size()){
// 			TecUtilDialogErrMsg("Failed to generate ordered sphere-surface intersection points.");
// 		}

	}

	return IntPoints;
}

int FESurface_c::GetClosestNodeToPoint(vec3 const & Point, double * ClosestNodeDistance) const{
	int MinNodeInd = -1;
	if (this->IsMade()){
		double MinDistSqr = DBL_MAX,
			TmpDistSqr;

		for (int ni = 0; ni < this->m_XYZList.size(); ++ni){
			TmpDistSqr = DistSqr(Point, this->m_XYZList[ni]);
			if (TmpDistSqr < MinDistSqr){
				MinDistSqr = TmpDistSqr;
				MinNodeInd = ni;
			}
		}
		if (ClosestNodeDistance != nullptr)
			*ClosestNodeDistance = sqrt(MinDistSqr);
	}
	return MinNodeInd;
}



/*
 *	Begin Domain_c functions
 */


void Domain_c::Setup(vector<int> const & V, FESurface_c *Vol)
{
	m_V = V;
	m_Vol = Vol;
	m_E.resize(m_V.size(), vector<int>(2));
	for (int i = 0; i < m_V.size(); ++i){
		for (int j = 0; j < 2; ++j){
			m_E[i][j] = (i + j) % m_V.size();
		}
	}
}

double Domain_c::Weight() const
{
	if (m_V.size() <= 2)
		return 0.0;
	vector<double> MinTriAreas(m_E.size() - 2);
	vector<vector<int> > MinTris(m_E.size() - 2);
	vector<Domain_c> D1(m_E.size() - 2),
		D2(m_E.size() - 2);
	int tMin;
	double MinArea = DBL_MAX;
	for (int i = 0; i < MinTriAreas.size(); ++i){
		MinTriAreas[i] = Split(i, MinTris[i], D1[i], D2[i]);
		if (MinTriAreas[i] < MinArea){
			MinArea = MinTriAreas[i];
			tMin = i;
		}
	}
	m_Vol->m_ElemList.push_back(MinTris[tMin]);
	return MinTriAreas[tMin] + D1[tMin].Weight() + D2[tMin].Weight();
}

double Domain_c::Split(int ei, vector<int> & t, Domain_c & D1, Domain_c & D2) const
{
	vector<vector<vector<int> > > TriList(m_V.size() - 2 - ei, vector<vector<int>  >(2, vector<int>(3)));
	for (int i = 0; i < TriList.size(); ++i){
		for (int j = 0; j < 2; ++j){
			TriList[i][0][j] = m_V[m_E[ei][j]];
			TriList[i][1][j] = m_E[ei][j];
		}
		TriList[i][0][2] = m_V[(m_E[ei][1] + 1 + i) % m_V.size()];
		TriList[i][1][2] = (m_E[ei][1] + 1 + i) % m_V.size();
	}
	int tMin;
	double TempArea, MinArea = 1e150;
	for (int i = 0; i < TriList.size(); ++i){
		TempArea = WeightFunc(TriList[i][0]);
		if (TempArea < MinArea){
			MinArea = TempArea;
			tMin = i;
		}
	}
	t = TriList[tMin][0];
	vector<int> ti = TriList[tMin][1];
	if (abs(ti[0] - ti[2]) < abs(ti[1] - ti[2]) && ti[0] < ti[1]){
		vector<int> V(ti[0] - ti[2] + 1);
		for (int i = 0; i < V.size(); ++i){
			V[i] = m_V[ti[2] + i];
		}
		D1.Setup(V, m_Vol);
		V = vector<int>(m_V.size() - ti[1] + ti[2] + 1);
		for (int i = 0; i < V.size(); ++i){
			V[i] = m_V[(ti[1] + i) % m_V.size()];
		}
		D2.Setup(V, m_Vol);
	}
	else{
		vector<int> V(m_V.size() + ti[0] - ti[2] + 1);
		for (int i = 0; i < V.size(); ++i){
			V[i] = m_V[(ti[2] + i) % m_V.size()];
		}
		D1.Setup(V, m_Vol);
		V = vector<int>(ti[2] - ti[1] + 1);
		for (int i = 0; i < V.size(); ++i){
			V[i] = m_V[ti[1] + i];
		}
		D2.Setup(V, m_Vol);
	}

	return MinArea;
}

double Domain_c::WeightFunc(vector<int> const & t) const
{
	vec3 AB = m_Vol->m_XYZList[t[1]] - m_Vol->m_XYZList[t[0]],
		AC = m_Vol->m_XYZList[t[2]] - m_Vol->m_XYZList[t[0]];

	double Area = (pow(AB[1] * AC[2] - AB[2] * AC[1], 2)
		+ pow(AB[2] * AC[0] - AB[0] * AC[2], 2)
		+ pow(AB[0] * AC[1] - AB[1] * AC[0], 2));

	return Area;
}


