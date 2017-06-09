/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef SURFACEFIT_H_
#define SURFACEFIT_H_

typedef enum
{
    SurfFitType_NoFit,
    SurfFitType_Planar,
    SurfFitType_Quadratic,
    /* BEGINREMOVEFROMADDON */
    END_SurfFitType_e,
    /* ENDREMOVEFROMADDON */
    SurfFitType_Invalid = BadEnumValue
} SurfFitType_e;

/* private SurfaceFit structure */
typedef struct _SurfaceFit_s
{
    // Fit Type
    SurfFitType_e SurfFitType;
    // Temporary matrices,
    //double    **a;  // NumPoints x NumSurfFitParams
    //double    **v;  // NumSurfFitParams x NumSurfFitParams
    //double    *b;   // NumPoints
    //double    *x;   // NumSurfFitParams
    //double    *w;   // NumSurfFitParams
    // SurfaceFit Inputs: DepVar1, DepVar2, IndepVar. Index is point number.
    ArrList_pa Coord1;
    ArrList_pa Coord2;
    ArrList_pa Var;
    // SurfaceFit Outputs: Fit coefficients.
    ArrList_pa FitCoefs;
    // The following are for quadratic surface fit.
    // Critical Point
    double c1_crt;
    double c2_crt;
    // Eigenvalue of curvature matrix
    double lambda1;
    double lambda2;
    double RotationAngle;
} SurfaceFit_s;

typedef struct _SurfaceFit_s  *SurfaceFit_pa;

Boolean_t     SurfaceFitIsValid(SurfaceFit_pa SurfaceFit);
SurfaceFit_pa SurfaceFitAlloc();
void          SurfaceFitClear(SurfaceFit_pa SurfaceFit);
void          SurfaceFitDealloc(SurfaceFit_pa *SurfaceFit);
LgIndex_t     SurfaceFitGetPointCount(const SurfaceFit_pa SurfaceFit);
Boolean_t     SurfaceFitSetPointAtOffset(SurfaceFit_pa SurfaceFit,
                                         LgIndex_t  PointOffset,
                                         double     Coord1,
                                         double     Coord2,
                                         double     Var);
Boolean_t   SurfaceFitAppendPointAtEnd(SurfaceFit_pa SurfaceFit,
                                       double        Coord1,
                                       double        Coord2,
                                       double        Var);
Boolean_t SurfaceFitGetPointFromOffset(const SurfaceFit_pa  SurfaceFit,
                                       const LgIndex_t   PointOffset,
                                       double     *Coord1,
                                       double     *Coord2,
                                       double     *Var);
LgIndex_t     SurfaceFitGetCoefCount(const SurfaceFit_pa SurfaceFit);
Boolean_t     SurfaceFitSetCoefAtOffset(SurfaceFit_pa SurfaceFit,
                                        LgIndex_t  PointOffset,
                                        double     Coef);
Boolean_t   SurfaceFitAppendCoefAtEnd(SurfaceFit_pa SurfaceFit,
                                      double        Coef);
double  SurfaceFitGetCoefFromOffset(const SurfaceFit_pa  SurfaceFit,
                                    const LgIndex_t   PointOffset);
SurfFitType_e SurfaceFitGetType(const SurfaceFit_pa SurfaceFit);
//Boolean_t SurfaceFitGetCoefficients(SurfaceFit_pa  SurfaceFit,
//                                    double       **FitCoeffs);
Boolean_t SurfaceFitCompute(SurfaceFit_pa SurfaceFit,
                            SurfFitType_e SurfFitType);
double    SurfaceFitInterpolateVar(SurfaceFit_pa SurfaceFit,
                                   double        Coord1,
                                   double        Coord2);
Boolean_t SurfaceFitCompCritPoint(SurfaceFit_pa SurfaceFit);

Boolean_t SurfaceFitTest();



#endif /* SURFACEFIT_H_ */