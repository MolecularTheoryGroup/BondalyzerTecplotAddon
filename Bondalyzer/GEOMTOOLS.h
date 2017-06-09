/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef GEOMTOOLS_H_
#define GEOMTOOLS_H_

#define EPSILON_GT 0.000001

double DistanceSquaredXYZ(const XYZ_s A,
						  const XYZ_s B);
LgIndex_t GeomToolsClosestElemNum(EntIndex_t SurfZone,
                                  double     X,
                                  double     Y,
                                  double     Z,
                                  double    *XInter,
                                  double    *YInter,
                                  double    *ZInter);
Boolean_t GeomToolsLineSegTriangleIntersection(XYZ_s const & triCorner1,
                                               XYZ_s const & triCorner2,
                                               XYZ_s const & triCorner3,
                                               XYZ_s const & lineSegStart,
                                               XYZ_s const & lineSegEnd,
                                               XYZ_s *       naturalTUV,
                                               XYZ_s *       intersectionPt);
Boolean_t GeomToolsLineSegQuadIntersections(XYZ_s const & corner1,
                                            XYZ_s const & corner2,
                                            XYZ_s const & corner3,
                                            XYZ_s const & corner4,
                                            XYZ_s const & lineSegStart,
                                            XYZ_s const & lineSegEnd,
                                            LgIndex_t *   numIntersections,
                                            XYZ_s *       intersectionPts);
typedef void IntersectionCallback_pf(XYZ_s& intersectionPt,
                                     void*  clientData);
void GeomToolsLineSegSurfaceIntersection(EntIndex_t const        surfaceZoneNum, // 1-based
                                         XYZ_s const &           lineSegStartPt,
                                         XYZ_s const &           lineSegEndPt,
                                         IntersectionCallback_pf intersectionCallback,
                                         void *                  clientData);
Boolean_t GeomToolsTest();

double SurfaceLength(const double Radius, 
					 const double PointDistance);
#endif /* GEOMTOOLS_H_ */