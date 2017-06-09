/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef ELEMSHAPEFUNC_H_
#define ELEMSHAPEFUNC_H_
Boolean_t RectGridBrickXYZtoIJKRST(double          X,
                                   double          Y,
                                   double          Z,
                                   LgIndex_t       IMax,
                                   LgIndex_t       JMax,
                                   LgIndex_t       KMax,
                                   double          XBegZone,
                                   double          XEndZone,
                                   double          YBegZone,
                                   double          YEndZone,
                                   double          ZBegZone,
                                   double          ZEndZone,
								   Boolean_t       PeriodicBC,
                                   LgIndex_t      *i,
                                   LgIndex_t      *j,
                                   LgIndex_t      *k,
                                   double         *r,
                                   double         *s,
                                   double         *t);
Boolean_t BrickTrilinearWeight(double       r,
                               double       s,
                               double       t,
                               double       W[8]);
Boolean_t BrickTrilinearWeightDerivatives(double       r,
                                          double       s,
                                          double       t,
                                          double       W[8],
                                          double       dwdr[8],
                                          double       dwds[8],
                                          double       dwdt[8]);
Boolean_t TrilinearBrickNaturalCoord(double       X,
                                     double       Y,
                                     double       Z,
                                     double       XCell[8],
                                     double       YCell[8],
                                     double       ZCell[8],
                                     double       XCMin,
                                     double       YCMin,
                                     double       ZCMin,
                                     double       XCMax,
                                     double       YCMax,
                                     double       ZCMax,
                                     double       *r,
                                     double       *s,
                                     double       *t,
                                     double       W[8]);
Boolean_t QuadBilinearWeight(double       r,
                             double       s,
                             double       W[4]);
Boolean_t QuadBilinearWeightDerivatives(double       r,
                                        double       s,
                                        double       W[4],
                                        double       dwdr[4],
                                        double       dwds[4]);
Boolean_t LinearTriangleNaturalCoord(double       U,
                                     double       V,
                                     double       UCell[3],
                                     double       VCell[3],
                                     double       UCMin,
                                     double       VCMin,
                                     double       UCMax,
                                     double       VCMax,
                                     double       XCell[3],
                                     double       YCell[3],
                                     double       *X,
                                     double       *Y,
                                     double       W[3]);
Boolean_t BilinearQuadNaturalCoord(double       X,
                                   double       Y,
                                   double       XCell[4],
                                   double       YCell[4],
                                   double       XCMin,
                                   double       YCMin,
                                   double       XCMax,
                                   double       YCMax,
                                   double       *r,
                                   double       *s,
                                   double       W[4]);
Boolean_t XYZtoPsiEta(int     PsiAlign,
                      XYZ_s   Normal,
                      double  dX,
                      double  dY,
                      double  dZ,
                      double *dPsi,
                      double *dEta);
Boolean_t PsiEtatoXYZ(int     PsiAlign,
                      XYZ_s   Normal,
                      double  dPsi,
                      double  dEta,
                      double *dX,
                      double *dY,
                      double *dZ);
int CalcPsiAlign(Boolean_t    IsQuad,
                 double       XCell[4],
                 double       YCell[4],
                 double       ZCell[4],
                 double      *XAve,
                 double      *YAve,
                 double      *ZAve,
                 XYZ_pa       CellNormal);


#endif

