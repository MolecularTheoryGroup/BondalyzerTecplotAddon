#include "updateSphericalTriangulation.h"

namespace tpcsm {

namespace{
    class Arc
    {
    private:
        Vec3 m_p1;
        Vec3 m_p2;
    public:
        Arc(
            Vec3 const& p1,
            Vec3 const& p2)
            : m_p1(p1)
            , m_p2(p2)
        {
            static double const epsilon = 1.0e-7;
            REQUIRE(1.0 - epsilon < p1.getNorm() && p1.getNorm() < 1.0 + epsilon);
        }
        Vec3 const& p1() const { return m_p1; }
        void p1(Vec3 const& val) { m_p1 = val; }
        Vec3 const& p2() const { return m_p2; }
        void p2(Vec3 const& val) { m_p2 = val; }
    };
}

namespace{
    bool findArcIntersection(
        Arc const& e1,
        Arc const& e2,
        Vec3& intersection)
    {
        REQUIRE(e1.p1().isNormalized());
        REQUIRE(e1.p2().isNormalized());
        REQUIRE(e2.p1().isNormalized());
        REQUIRE(e2.p2().isNormalized());

        bool found = false;
        static double const epsilon = 1.0e-6;
        REQUIRE(e1.p1().dot(e1.p2()) > 0.5); // arc should be substantially less than entire circle
        REQUIRE(e2.p1().dot(e2.p2()) > 0.5);

        // avoid arcs too short to calculate intersection
        if (e1.p1().dot(e1.p2()) < 1.0 - epsilon &&
            e2.p1().dot(e2.p2()) < 1.0 - epsilon)
        {
            // find the normals for the planes of each segments (since these are great circles, a unique plane exists)
            Vec3 const e1PlaneNormal = e1.p1().cross(e1.p2()); // TODO: maybe precalculate this
            Vec3 const e2PlaneNormal = e2.p1().cross(e2.p2());

            // the intersection is in both planes, so perpendicular to both normals, so the cross product of those normals
            Vec3 planarIntersection = e1PlaneNormal.cross(e2PlaneNormal).normalize();

            // now the intersection might not be along a given segment
            double const aDotB = e1.p1().dot(e1.p2()); // TODO: maybe precalculate this?
            double aDotX = e1.p1().dot(planarIntersection);
            // we could end up with the intersection pointing the wrong direction
            if (aDotX < 0.0)
            {
                planarIntersection = -planarIntersection;
                aDotX = -aDotX;
            }
            double const bDotX = e1.p2().dot(planarIntersection);
            // bDotX could be <0 if the planar intersection is far away from segment

            // if aDotB is the lowest valued dot product, aDotB is the largest angle, and thus X must be between A&B
            // but add in an epsilon to account for numerical error
            if (aDotB < aDotX + epsilon && aDotB < bDotX + epsilon)
            {
                // now check other arc to be sure they are both on same side of sphere
                double const cDotD = e2.p1().dot(e2.p2()); // TODO: maybe precalculate this?
                double const cDotX = e2.p1().dot(planarIntersection);
                double const dDotX = e2.p2().dot(planarIntersection);

                if (cDotD < cDotX + epsilon && cDotD < dDotX + epsilon)
                {
                    intersection = planarIntersection;
                    found = true;
                }
            }
        }
        return found;
    }
}

namespace {
    bool coincident(Vec3 const& a, Vec3 const& b)
    {
        static double const epsilon = 1e-15;
        Vec3 const diff(a - b);
        double n2 = diff.getNormSquared(); // diff.max(-diff) < epsilon
        return n2 < epsilon;
    }
}

namespace {
    typedef Index3 TriNeighbors; // (triNeighborOppositeN1,triNeighborOppositeN2,triNeighborOppositeN3)
    Index_t changeNeighbor(
        TriNeighbors& triNeighbors,
        Index_t oldNeighbor,
        Index_t newNeighbor)
    {
        if (triNeighbors.v1() == oldNeighbor)
        {
            triNeighbors.setV1(newNeighbor);
            return true;
        }
        else if (triNeighbors.v2() == oldNeighbor)
        {
            triNeighbors.setV2(newNeighbor);
            return true;
        }
        else if (triNeighbors.v3() == oldNeighbor)
        {
            triNeighbors.setV3(newNeighbor);
            return true;
        }
        else
            return false;
    }

    void setEdgeTriNeighbor(
        Index_t         edgeNode1,
        Index_t         edgeNode2,
        TriNodes const& triNodes,
        Index_t         neighborTriPos,
        TriNeighbors&   triNeighbors)
    {
        if (edgeNode1 == triNodes.v1())
        {
            if (edgeNode2 == triNodes.v2())
                triNeighbors.setV3(neighborTriPos);
            else
            {
                CHECK(edgeNode2 == triNodes.v3());
                triNeighbors.setV2(neighborTriPos);
            }
        }
        else if (edgeNode1 == triNodes.v2())
        {
            if (edgeNode2 == triNodes.v1())
                triNeighbors.setV3(neighborTriPos);
            else
            {
                CHECK(edgeNode2 == triNodes.v3());
                triNeighbors.setV1(neighborTriPos);
            }
        }
        else
        {
            CHECK(edgeNode1 == triNodes.v3());
            if (edgeNode2 == triNodes.v1())
                triNeighbors.setV2(neighborTriPos);
            else
            {
                CHECK(edgeNode2 == triNodes.v2());
                triNeighbors.setV1(neighborTriPos);
            }
        }
    }
}

namespace {
    // remove coincident nodes and collapsed triangles from inTriangulation and return outTriangulation
    void simplifyTriangulation(
        std::vector<Vec3> const&     inNodeXYZs, // must be normalized
        std::vector<Index_t> const&  inNodeEdges,
        std::vector<TriNodes> const& inTriangleNodes,
        std::vector<Vec3>&           outNodeXYZs,
        std::vector<Index_t>&        outNodeEdges,
        std::vector<TriNodes>&       outTriangleNodes)
    {
        REQUIRE(inNodeXYZs.size() == inNodeEdges.size());
        // find coincident points
        static double const coincidentPtEpsilon = 2.0e-14; // empirically determined via a spherical zone generated within Tecplot (2e-14 works, 1e-14 failes)
        std::vector<Index_t> uniquePos;
        {
            Index_t const numTriangleNodes = Index_t(inNodeXYZs.size());
            uniquePos.reserve(numTriangleNodes);
            for (Index_t pos = 0; pos < numTriangleNodes; ++pos)
            {
                Vec3 const normalizedXYZ(inNodeXYZs[pos]);
                CHECK(normalizedXYZ.isNormalized());
                Index_t const nodeEdge = inNodeEdges[pos];

                // see if any coincident points are already in the list
                bool isCoincident = false;
                for (Index_t pos2 = 0; pos2 < pos; ++pos2) // this is O(n^2) TODO: Replace with faster algorithm
                {
                    Index_t const uniqueNodePos = uniquePos[pos2];
                    if (normalizedXYZ.dot(outNodeXYZs[uniqueNodePos]) > 1.0 - coincidentPtEpsilon)
                    {
                        isCoincident = true;
                        uniquePos.push_back(uniquePos[pos2]);
                        break;
                    }
                }
                if (!isCoincident)
                {
                    Index_t const uniqueNodePos = Index_t(outNodeXYZs.size());
                    uniquePos.push_back(uniqueNodePos);
                    outNodeXYZs.push_back(normalizedXYZ);
                    outNodeEdges.push_back(nodeEdge);
                }
            }
        }
        CHECK(uniquePos.size() == inNodeXYZs.size());

        // renumber nodes in triangles while we reorient them
        for (auto iter = inTriangleNodes.begin(); iter != inTriangleNodes.end(); ++iter)
        {
            Index_t const n1 = uniquePos[iter->v1()];
            Index_t const n2 = uniquePos[iter->v2()];
            Index_t const n3 = uniquePos[iter->v3()];

            if (n1 != n2 && n2 != n3 && n3 != n1) // skip degenerate triangles
            {
                Vec3 const p1 = outNodeXYZs[n1];
                Vec3 const p2 = outNodeXYZs[n2];
                Vec3 const p3 = outNodeXYZs[n3];

                // enforce consistent orientation
                static double const epsilon = 1e-15;
                double const tp = Vec3::tripleProduct(p1, p2, p3);
                TriNodes tri;
                if (tp < 0.0)
                {
                    CHECK(tp < -epsilon); // with coincident node removal, we should have weeded out degenerate triangle above
                    tri = TriNodes(n1, n3, n2);
                }
                else
                {
                    CHECK(tp > epsilon); // with coincident node removal, we should have weeded out degenerate triangle above
                    tri = TriNodes(n1, n2, n3);
                }
                Index_t const newTriPos = Index_t(outTriangleNodes.size());
                outTriangleNodes.push_back(tri);
            }
        }

        ENSURE(outNodeXYZs.size() == outNodeEdges.size());
    }
}

bool updateSphericalTriangulation(
    std::vector<Vec3> const&     startNodeXYZs, // IN: location of the initial triangulation's nodes
    std::vector<TriNodes> const& startTriangleNodes, // IN: connectivity of each triangle, defined in terms of points above
    Vec3 const&                  sphereCenter,
    double                       radius,
    std::vector<Vec3> const&     constraintNodes, // IN:
    std::vector<Edge> const&     constraintSegments, // IN: defined in terms of constraintNodes above
    bool                         includeConstraintSegmentsOnly, // IN: false means include intersections 
    std::vector<Vec3>&           nodeXYZs, // OUT: location of the resulting triangulation's nodes (includes all of those in initialNodeXYZs)
    std::vector<Index_t>&        nodeEdge, // OUT: for each node in nodeXYZ which constraintSegments it lines on or -1 if none
    std::vector<TriNodes>&       triangleNodes, // OUT: connectivity of each triangle, defined in terms of points above
    char const*&                 statusMessage) // OUT: contains description of error or "Triangulation Successful" if successful
{
    bool result = false;

    for (auto pos = 0; pos < startNodeXYZs.size(); ++pos)
    {
        REQUIRE(startTriangleNodes[pos].v1() < startNodeXYZs.size());
        REQUIRE(startTriangleNodes[pos].v2() < startNodeXYZs.size());
        REQUIRE(startTriangleNodes[pos].v3() < startNodeXYZs.size());
    }

    try
    {
        std::vector<Vec3> initialNodeXYZs(startNodeXYZs);
        std::vector<TriNodes> initialTriangleNodes(startTriangleNodes);

#if 0
        for (int pass = 0; pass < 1; pass++)
        {
            // refine existing mesh for performance testing
            std::vector<Vec3> refinedNodesXYZs;
            std::vector<TriNodes> refinedTriangleNodes;
            refinedNodesXYZs.reserve(initialNodeXYZs.size() * 2);
            refinedTriangleNodes.reserve(initialTriangleNodes.size() * 4);
            for (auto iter = initialTriangleNodes.begin(); iter != initialTriangleNodes.end(); ++iter)
            {
                Index_t const n1 = iter->v1();
                Index_t const n2 = iter->v2();
                Index_t const n3 = iter->v3();
                Vec3 const p1 = initialNodeXYZs[n1];
                Vec3 const p2 = initialNodeXYZs[n2];
                Vec3 const p3 = initialNodeXYZs[n3];
                Vec3 const p12 = (p1 + p2) / 2.0;
                Vec3 const p23 = (p2 + p3) / 2.0;
                Vec3 const p31 = (p3 + p1) / 2.0;
                Index_t const rnBase = Index_t(refinedNodesXYZs.size());
                refinedNodesXYZs.push_back(p1);
                refinedNodesXYZs.push_back(p2);
                refinedNodesXYZs.push_back(p3);
                refinedNodesXYZs.push_back(p12);
                refinedNodesXYZs.push_back(p23);
                refinedNodesXYZs.push_back(p31);
                refinedTriangleNodes.push_back(TriNodes(rnBase + 0, rnBase + 3, rnBase + 5));
                refinedTriangleNodes.push_back(TriNodes(rnBase + 3, rnBase + 1, rnBase + 4));
                refinedTriangleNodes.push_back(TriNodes(rnBase + 3, rnBase + 4, rnBase + 5));
                refinedTriangleNodes.push_back(TriNodes(rnBase + 5, rnBase + 4, rnBase + 2));
            }

            std::swap(initialNodeXYZs, refinedNodesXYZs);
            std::swap(initialTriangleNodes, refinedTriangleNodes);
        }
#endif

        // reserve space to avoid resizing as we go which could create O(n^2) performance
        nodeXYZs.reserve(initialNodeXYZs.size() + constraintNodes.size());
        triangleNodes.reserve(initialTriangleNodes.size() + 2 * constraintNodes.size() + constraintSegments.size());

        // normalize XYZs of initial triangulation
        std::vector<Vec3> normalizedXYZs;
        normalizedXYZs.reserve(nodeXYZs.size());
        std::vector<Index_t> initialNodeEdges;
        initialNodeEdges.reserve(nodeXYZs.size());
        for (auto iter = initialNodeXYZs.begin(); iter != initialNodeXYZs.end(); ++iter)
        {
            normalizedXYZs.push_back((*iter - sphereCenter).normalize());
            initialNodeEdges.push_back(-1); // part of initial triangulation and thus not on constraint edge
        }

        // simplify the triangulation
        simplifyTriangulation(normalizedXYZs, initialNodeEdges, initialTriangleNodes, nodeXYZs, nodeEdge, triangleNodes);

#if 0
        // collect edges in order to determine triangle neighbors
        typedef Index3 OneSidedEdge; // (edgeNode1, edgeNode2, tri), one-sided because it doesn't contain the other triangle
        std::vector<OneSidedEdge> oneSidedEdges;
        oneSidedEdges.reserve(triangleNodes.size() * 3);
        for (Index_t triPos = 0; triPos < triangleNodes.size(); ++triPos)
        {
            Index_t const n1 = triangleNodes[triPos].v1();
            Index_t const n2 = triangleNodes[triPos].v2();
            Index_t const n3 = triangleNodes[triPos].v3();

            // make sure edges go the same direction (low node to high node)
            if (n1 < n2)
                oneSidedEdges.push_back(OneSidedEdge(n1, n2, triPos));
            else
                oneSidedEdges.push_back(OneSidedEdge(n2, n1, triPos));

            if ( n2 < n3 )
                oneSidedEdges.push_back(OneSidedEdge(n2, n3, triPos));
            else
                oneSidedEdges.push_back(OneSidedEdge(n3, n2, triPos));

            if ( n3 < n1 )
                oneSidedEdges.push_back(OneSidedEdge(n3, n1, triPos));
            else
                oneSidedEdges.push_back(OneSidedEdge(n1, n3, triPos));
        }

        // sorting associates one-sided edges so we can find the adjacent triangles
        std::sort(oneSidedEdges.begin(), oneSidedEdges.end());

        std::vector<TriNeighbors> triNeighbors;
        static Index_t const BAD_TRIANGLE = Index_t(-1);
        triNeighbors.reserve(triangleNodes.capacity());
        triNeighbors.resize(triangleNodes.size(), TriNeighbors(BAD_TRIANGLE, BAD_TRIANGLE, BAD_TRIANGLE));

        Index_t numEdges = Index_t(oneSidedEdges.size());
        CHECK(numEdges % 2 == 0);
        for (Index_t edgePos = 0; edgePos < numEdges; edgePos+=2)
        {
            // oneSidedEdges contains two records for each edge, one for each triangle, sorted to be next to each other.
            CHECK(oneSidedEdges[edgePos].v1() == oneSidedEdges[edgePos + 1].v1());
            CHECK(oneSidedEdges[edgePos].v2() == oneSidedEdges[edgePos + 1].v2());

            Index_t const edge1Node1 = oneSidedEdges[edgePos].v1();
            Index_t const edge1Node2 = oneSidedEdges[edgePos].v2();
            Index_t const tri1Pos = oneSidedEdges[edgePos].v3();
            TriNodes tri1Nodes = triangleNodes[tri1Pos];

            Index_t const edge2Node1 = oneSidedEdges[edgePos + 1].v1();
            Index_t const edge2Node2 = oneSidedEdges[edgePos + 1].v2();
            Index_t const tri2Pos = oneSidedEdges[edgePos + 1].v3();
            TriNodes tri2Nodes = triangleNodes[tri2Pos];

            setEdgeTriNeighbor(edge1Node1, edge1Node2, tri1Nodes, tri2Pos, triNeighbors[tri1Pos]);
            setEdgeTriNeighbor(edge2Node1, edge2Node2, tri2Nodes, tri1Pos, triNeighbors[tri2Pos]);
        }

        // check triangle neighbors
        for (auto iter = triNeighbors.begin(); iter != triNeighbors.end(); ++iter)
        {
            CHECK(iter->v1() != BAD_TRIANGLE);
            CHECK(iter->v2() != BAD_TRIANGLE);
            CHECK(iter->v3() != BAD_TRIANGLE);
        }
#endif

        // add the new nodes splitting up the triangles as we go. This makes things simpler later.
        static double const TRIPLE_PRODUCT_EPSILON = 1.0e-14; // used to determine if point is inside of triangle (TODO:should this be the same as coincidentPtEpsilon?)

        // Because we are setting the edge flag for new points, we need to do this by looping over the constraints
        // and not the constraint nodes.

        // To make this much faster keep track of which points we already added instead of relying looping over the triangles and calcing triple products
        Index_t const numConstraintNodes = Index_t(constraintNodes.size());
        std::vector<bool> constraintNodeAdded;
        constraintNodeAdded.resize(numConstraintNodes, false);

        Index_t const numConstraintSegments = Index_t(constraintSegments.size());
        for (Index_t segPos = 0; segPos < numConstraintSegments; ++segPos)
        {
            for (int whichEnd = 0; whichEnd <= 1; ++whichEnd)
            {
                Index_t const constraintNodePos = (whichEnd == 0) ? constraintSegments[segPos].first : constraintSegments[segPos].second;
                //if (!constraintNodeAdded[constraintNodePos])
                {
                    constraintNodeAdded[constraintNodePos] = true;

                    Vec3 const pt = constraintNodes[constraintNodePos];

                    Index_t const newNode = Index_t(nodeXYZs.size());
                    Vec3 const newPt = (pt - sphereCenter).normalize();
                    nodeXYZs.push_back(newPt);
                    nodeEdge.push_back(segPos);

                    for (Index_t triPos = 0; triPos < triangleNodes.size(); ++triPos)
                    {
                        Index_t const n1 = triangleNodes[triPos].v1();
                        Index_t const n2 = triangleNodes[triPos].v2();
                        Index_t const n3 = triangleNodes[triPos].v3();
                        Vec3 const p1 = nodeXYZs[n1]; // newNodes are already normalized
                        Vec3 const p2 = nodeXYZs[n2];
                        Vec3 const p3 = nodeXYZs[n3];

                        if (Vec3::tripleProduct(p1, p2, p3) > TRIPLE_PRODUCT_EPSILON &&
                            Vec3::tripleProduct(p1, p2, newPt) > TRIPLE_PRODUCT_EPSILON &&
                            Vec3::tripleProduct(p1, newPt, p3) > TRIPLE_PRODUCT_EPSILON &&
                            Vec3::tripleProduct(newPt, p2, p3) > TRIPLE_PRODUCT_EPSILON)
                        {
                            // split existing triangle into three
                            Index_t const t1 = triPos;
                            Index_t const t2 = Index_t(triangleNodes.size());
                            Index_t const t3 = t2 + 1;

                            triangleNodes[triPos] = TriNodes(n1, n2, newNode);
                            triangleNodes.push_back(TriNodes(n2, n3, newNode));
                            triangleNodes.push_back(TriNodes(n3, n1, newNode));

#if 0 // neighbors were an attempted optimization cut due to development time
                            // fix up neighbors
                            TriNeighbors const existingTriNeighbors = triNeighbors[triPos];
                            triNeighbors[triPos] = TriNeighbors(t2, t3, existingTriNeighbors.v3());
                            CHECK(changeNeighbor(triNeighbors[existingTriNeighbors.v3()], triPos, t1)); // just check for neighbor triPos

                            triNeighbors.push_back(TriNeighbors(t3, triPos, existingTriNeighbors.v1()));
                            changeNeighbor(triNeighbors[existingTriNeighbors.v1()], triPos, t2);

                            triNeighbors.push_back(TriNeighbors(triPos, t2, existingTriNeighbors.v2()));
                            changeNeighbor(triNeighbors[existingTriNeighbors.v2()], triPos, t3);
#endif

                            break;
                        }
                    }
                }
            }
        }

#if 0 // check triangle neighbors again
        for (auto iter = triNeighbors.begin(); iter != triNeighbors.end(); ++iter)
        {
            CHECK(iter->v1() != BAD_TRIANGLE);
            CHECK(iter->v2() != BAD_TRIANGLE);
            CHECK(iter->v3() != BAD_TRIANGLE);
        }
#endif

        // includeConstraintSegments is a debugging tool to make sure the correct segments were added
        // this results in those segments as degenerate triangles in the new mesh, but no intersections
        if (includeConstraintSegmentsOnly)
        {
            Index_t offset = Index_t(initialNodeXYZs.size());
            for (auto iter = constraintSegments.begin(); iter != constraintSegments.end(); ++iter)
            {
                TriNodes tri(iter->first + offset, iter->second + offset, iter->second + offset);
                triangleNodes.push_back(tri);
            }
        }
        else
        {
            // find constraintSegment intersections, normalizing as we go
            Index_t const segPosSize = Index_t(constraintSegments.size());
            for (Index_t segPos = 0; segPos < segPosSize; ++segPos)
            {
                Index_t const first = constraintSegments[segPos].first;
                Index_t const second = constraintSegments[segPos].second;

                // normalize as we go
                Vec3 a = (constraintNodes[first] - sphereCenter).normalize();
                Vec3 b = (constraintNodes[second] - sphereCenter).normalize();
                Arc constraintArc(a, b);

                // For now we do this triangle by triangle which introduces duplicate points that we remove at the end.
                // TODO: get edge by edge version working
                // We can't use iterators to can't loop over newTriangles because we add triangles to the end 
                // and this might resize the vector and thus invalidate the iterator. So we use an index instead.
                for (size_t triPos = 0; triPos < triangleNodes.size(); ++triPos)
                {
                    Index_t const n1 = triangleNodes[triPos].v1();
                    Index_t const n2 = triangleNodes[triPos].v2();
                    Index_t const n3 = triangleNodes[triPos].v3();
                    Vec3 const p1 = nodeXYZs[n1]; // newNodes are already normalized
                    Vec3 const p2 = nodeXYZs[n2];
                    Vec3 const p3 = nodeXYZs[n3];

                    // weed out degenerate triangles: collapsed (repeated points) or "collinear" (on the same spherical arc)
                    double const tripleProduct = std::abs(Vec3::tripleProduct(p1, p2, p3));
                    if (tripleProduct < TRIPLE_PRODUCT_EPSILON)
                    {
                        // do we need special processing for these?
                    }
                    else
                    {
                        // out first question is, are the arc end points coincident with any edge
                        bool aIsP1 = coincident(a, p1);
                        bool bIsP1 = coincident(b, p1);
                        bool aIsP2 = coincident(a, p2);
                        bool bIsP2 = coincident(b, p2);
                        bool aIsP3 = coincident(a, p3);
                        bool bIsP3 = coincident(b, p3);

                        // to simplify the processing, reverse the direction of the arc if b is coincident
                        // so we can always assume a is coincident
                        if (!aIsP1 && !aIsP2 && !aIsP3)
                        {
                            CHECK(!aIsP1);
                            std::swap(a, b);
                            std::swap(aIsP1, bIsP1);
                            std::swap(aIsP2, bIsP2);
                            std::swap(aIsP3, bIsP3);
                            constraintArc = Arc(a, b);
                        }

                        Vec3 p12;
                        Vec3 p23;
                        Vec3 p31;
                        bool const i12 = findArcIntersection(constraintArc, Arc(p1, p2), p12);
                        bool const i23 = findArcIntersection(constraintArc, Arc(p2, p3), p23);
                        bool const i31 = findArcIntersection(constraintArc, Arc(p3, p1), p31);

                        if (aIsP1)
                        {
                            CHECK(!bIsP1); // a&b must be distinct
                            CHECK(!aIsP2 && !aIsP3); // triangle points must be distinct
                            //CHECK(i31 && i12); // intersection is at P1
                            if (i23) // is there an intersection on adjacent edge?
                            {
                                if (!bIsP2 && !bIsP3) // is B not one of the other points
                                {
                                    if (Vec3::tripleProduct(p1, p2, p23) > TRIPLE_PRODUCT_EPSILON &&
                                        Vec3::tripleProduct(p3, p1, p23) > TRIPLE_PRODUCT_EPSILON)
                                    {
                                        Index_t const n23 = Index_t(nodeXYZs.size());
                                        nodeXYZs.push_back(p23);
                                        nodeEdge.push_back(segPos);

                                        // record the two triangle split
                                        triangleNodes[triPos] = TriNodes(n1, n2, n23);
                                        triangleNodes.push_back(TriNodes(n3, n1, n23));

                                        // To record the edge as a collapsed triangle
                                        // triangleNodes.push_back(TriNodes(n1, n23, n23));
                                        // triNeighbors.push_back(TriNeighbors(BAD_TRIANGLE, BAD_TRIANGLE, BAD_TRIANGLE));
                                        continue;
                                    }
                                }
                            }
                        }
                        else if (aIsP2)
                        {
                            CHECK(!bIsP2); // a&b must be distinct
                            CHECK(!aIsP3 && !aIsP1); // triangle points must be distinct
                            //CHECK(i12 && i23); // intersection is at P2
                            if (i31) // is there an intersection on adjacent edge?
                            {
                                if (!bIsP3 && !bIsP1) // is B not one of the other points
                                {
                                    if (Vec3::tripleProduct(p2, p3, p31) > TRIPLE_PRODUCT_EPSILON &&
                                        Vec3::tripleProduct(p1, p2, p31) > TRIPLE_PRODUCT_EPSILON)
                                    {
                                        Index_t const n31 = Index_t(nodeXYZs.size());
                                        nodeXYZs.push_back(p31);
                                        nodeEdge.push_back(segPos);

                                        // record the two triangle split
                                        triangleNodes[triPos] = TriNodes(n2, n3, n31);
                                        triangleNodes.push_back(TriNodes(n1, n2, n31));

                                        // To record the edge as a collapsed triangle
                                        //triangleNodes.push_back(TriNodes(n2, n31, n31));
                                        //triNeighbors.push_back(TriNeighbors(BAD_TRIANGLE, BAD_TRIANGLE, BAD_TRIANGLE));

                                        continue;
                                    }
                                }
                            }
                        }
                        else if (aIsP3)
                        {
                            CHECK(!bIsP3); // a&b must be distinct
                            CHECK(!aIsP1 && !aIsP2); // triangle points must be distinct
                            //CHECK(i23 && i31); // intersection is at P3
                            if (i12) // is there an intersection on adjacent edge?
                            {
                                if (!bIsP1 && !bIsP2) // is B not one of the other points
                                {
                                    if (Vec3::tripleProduct(p3, p1, p12) > TRIPLE_PRODUCT_EPSILON &&
                                        Vec3::tripleProduct(p2, p3, p12) > TRIPLE_PRODUCT_EPSILON)
                                    {
                                        Index_t const n12 = Index_t(nodeXYZs.size());
                                        nodeXYZs.push_back(p12);
                                        nodeEdge.push_back(segPos);

                                        // record the two triangle split
                                        triangleNodes[triPos] = TriNodes(n3, n1, n12);
                                        triangleNodes.push_back(TriNodes(n2, n3, n12));

                                        // To record the edge as a collapsed triangle
                                        //triangleNodes.push_back(TriNodes(n3, n12, n12));
                                        //triNeighbors.push_back(TriNeighbors(BAD_TRIANGLE, BAD_TRIANGLE, BAD_TRIANGLE));

                                        continue;
                                    }
                                }
                            }
                        }
                        else if (i12 && i23 && i31)
                        {
                            //throw "bad triangles";
                        }
                        else if (i12 && i23)
                        {
                            if (!aIsP2 && !bIsP2)
                            {
                                if (Vec3::tripleProduct(p2, p23, p12) > TRIPLE_PRODUCT_EPSILON &&
                                    Vec3::tripleProduct(p3, p12, p23) > TRIPLE_PRODUCT_EPSILON &&
                                    Vec3::tripleProduct(p1, p12, p3) > TRIPLE_PRODUCT_EPSILON)
                                {
                                    Index_t const n12 = Index_t(nodeXYZs.size());
                                    nodeXYZs.push_back(p12);
                                    nodeEdge.push_back(segPos);

                                    Index_t const n23 = Index_t(nodeXYZs.size());
                                    nodeXYZs.push_back(p23);
                                    nodeEdge.push_back(segPos);

                                    triangleNodes[triPos] = TriNodes(n2, n23, n12);
                                    triangleNodes.push_back(TriNodes(n3, n12, n23)); // TODO: Choose between two potential splits based on a goodness of triangles formed
                                    triangleNodes.push_back(TriNodes(n1, n12, n3));

                                    // To record the edge as a collapsed triangle
                                    //triangleNodes.push_back(TriNodes(n12, n23, n23));
                                    //triNeighbors.push_back(TriNeighbors(BAD_TRIANGLE, BAD_TRIANGLE, BAD_TRIANGLE));
                                    continue;
                                }
                            }
                        }
                        else if (i23 && i31)
                        {
                            if (!aIsP3 && !bIsP3)
                            {
                                if (Vec3::tripleProduct(p3, p31, p23) > TRIPLE_PRODUCT_EPSILON &&
                                    Vec3::tripleProduct(p1, p23, p31) > TRIPLE_PRODUCT_EPSILON &&
                                    Vec3::tripleProduct(p2, p23, p1) > TRIPLE_PRODUCT_EPSILON)
                                {
                                    Index_t const n23 = Index_t(nodeXYZs.size());
                                    nodeXYZs.push_back(p23);
                                    nodeEdge.push_back(segPos);

                                    Index_t const n31 = Index_t(nodeXYZs.size());
                                    nodeXYZs.push_back(p31);
                                    nodeEdge.push_back(segPos);

                                    triangleNodes[triPos] = TriNodes(n3, n31, n23);
                                    triangleNodes.push_back(TriNodes(n1, n23, n31)); // TODO: Choose between two potential splits based on a goodness of triangles formed
                                    triangleNodes.push_back(TriNodes(n2, n23, n1));

                                    // To record the edge as a collapsed triangle
                                    //triangleNodes.push_back(TriNodes(n23, n31, n31));
                                    //triNeighbors.push_back(TriNeighbors(BAD_TRIANGLE, BAD_TRIANGLE, BAD_TRIANGLE));
                                    continue;
                                }
                            }
                        }
                        else if (i31 && i12)
                        {
                            if (!aIsP1 && !bIsP1)
                            {
                                if (Vec3::tripleProduct(p1, p12, p31) > TRIPLE_PRODUCT_EPSILON &&
                                    Vec3::tripleProduct(p2, p31, p12) > TRIPLE_PRODUCT_EPSILON &&
                                    Vec3::tripleProduct(p3, p31, p2) > TRIPLE_PRODUCT_EPSILON)
                                {
                                    Index_t const n31 = Index_t(nodeXYZs.size());
                                    nodeXYZs.push_back(p31);
                                    nodeEdge.push_back(segPos);

                                    Index_t const n12 = Index_t(nodeXYZs.size());
                                    nodeXYZs.push_back(p12);
                                    nodeEdge.push_back(segPos);

                                    triangleNodes[triPos] = TriNodes(n1, n12, n31);
                                    triangleNodes.push_back(TriNodes(n2, n31, n12)); // TODO: Choose between two potential splits based on a goodness of triangles formed
                                    triangleNodes.push_back(TriNodes(n3, n31, n2));

                                    // To record the edge as a collapsed triangle
                                    //triangleNodes.push_back(TriNodes(n31, n12, n12));
                                    //triNeighbors.push_back(TriNeighbors(BAD_TRIANGLE, BAD_TRIANGLE, BAD_TRIANGLE));
                                    continue;
                                }
                            }
                        }
                        else
                        {
                            // no intersections
                        }
                    }
                }
            }

            // re-simplify the triangulation
            std::vector<Vec3> simplifiedXYZs;
            std::vector<Index_t> simplifiedNodeEdge;
            std::vector<Index3> simplifiedTriangleNodes;
            simplifyTriangulation(nodeXYZs, nodeEdge, triangleNodes, simplifiedXYZs, simplifiedNodeEdge, simplifiedTriangleNodes);
            // put the simplified version back into the return values
            std::swap(nodeXYZs, simplifiedXYZs);
            std::swap(nodeEdge, simplifiedNodeEdge);
            std::swap(triangleNodes, simplifiedTriangleNodes);
        }

        // un-normalize the points. TODO: Should we just get the original values?
        for (auto iter = nodeXYZs.begin(); iter != nodeXYZs.end(); ++iter)
            *iter = (*iter * radius) + sphereCenter;

        result = true;
        statusMessage = "Triangulation update successful";
    }
    catch (char const* message)
    {
        statusMessage = message;
        result = false;
    }
    return result;
}

}
