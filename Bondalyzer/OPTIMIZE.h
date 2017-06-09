/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef OPTIMIZE_H_
#define OPTIMIZE_H_

typedef enum
{
    OptimizeMethod_NoOpt,
    OptimizeMethod_DownhillSimplex,
    OptimizeMethod_PowellsConjugateGrad,
    /* BEGINREMOVEFROMADDON */
    END_OptimizeMethod_e,
    /* ENDREMOVEFROMADDON */
    OptimizeMethod_Invalid = BadEnumValue
} OptimizeMethod_e;

/* private OptPoint structure */
typedef struct _OptPoint_s
{
    double    *Coords; // NDim Coordinates of a point in parameter space
    double     Value;  // Value of function at point
} OptPoint_s;

typedef struct _OptPoint_s  *OptPoint_pa;

/* private Optimize structure */
typedef struct _Optimize_s
{
    int         NDim;   // Dimension of parameter-space for optimization
    OptimizeMethod_e OptimizeMethod;
    double      Tolerance; // Convergence tolerance
    int         MaxIter;   // Maximum number of iterations allowed
    ArrList_pa  OptPoints; // NDim+1 points for DownhillSimplex
} Optimize_s;

typedef struct _Optimize_s  *Optimize_pa;


/* private OptTestInfo structure - Client data for test Function */
typedef struct _OptTestInfo_s
{
    int        NDims;   // NDim Coordinates of a point in parameter space
    double    *Origin;  // Origin of test function
} OptTestInfo_s;

typedef struct _OptTestInfo_s  *OptTestInfo_pa;


Boolean_t     OptPointIsValid(OptPoint_pa OptPoint);
OptPoint_pa   OptPointAlloc(LgIndex_t NDim);
void          OptPointDealloc(OptPoint_pa *OptPoint);
double*       OptPointGetCoordPointer(const OptPoint_pa OptPoint);

Boolean_t     OptimizeIsValid(Optimize_pa Optimize);
Optimize_pa   OptimizeAlloc(int NDim);
void          OptimizeClear(Optimize_pa Optimize);
void          OptimizeDealloc(Optimize_pa *Optimize);
LgIndex_t     OptimizeGetPointCount(const Optimize_pa Optimize);
LgIndex_t     OptimizeGetDimensions(const Optimize_pa Optimize);
Boolean_t     OptimizeSetOptPointAtOffset(Optimize_pa Optimize,
                                          LgIndex_t   PointOffset,
                                          OptPoint_pa OptPoint);
Boolean_t     OptimizeSetPointAtOffset(Optimize_pa Optimize,
                                       LgIndex_t   PointOffset,
                                       double     *Coords,
                                       double      Value);
double        OptimizeGetValueAtOffset(const Optimize_pa Optimize,
                                       const LgIndex_t   PointOffset);
Boolean_t     OptimizeSwapPointsAtOffsets(Optimize_pa Optimize,
                                          LgIndex_t   PointOffset1,
                                          LgIndex_t   PointOffset2);
OptPoint_pa OptimizeGetPointFromOffset(const Optimize_pa  Optimize,
                                       const LgIndex_t    PointOffset);
OptimizeMethod_e OptimizeGetMethod(const Optimize_pa Optimize);
OptPoint_pa OptimizeGetPointFromOffset(const Optimize_pa  Optimize,
                                       const LgIndex_t    PointOffset); Boolean_t OptimizeAveSimplexPoints(Optimize_pa Optimize);
LgIndex_t OptimizeCompute(Optimize_pa      Optimize,
                          OptimizeMethod_e OptimizeMethod,
                          void        *ClientData,
                          Boolean_t (*Function)(const void *, const double *, double *));
Boolean_t TestFunction(const void       *ClientData,
                       const double     *Coords,
                       double           *Value);


#endif /* OPTIMIZE_H_ */