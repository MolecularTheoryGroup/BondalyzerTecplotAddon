/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef ODERUNGEKUTTA_H_
#define ODERUNGEKUTTA_H_

Boolean_t RK2UpdateSolution(const EntIndex_t  NumSolVars,
							const double      b2,
							double      *DelTime,
							const Boolean_t   ReverseDir,
							const void       *ClientData,
							Boolean_t (*Function)(const void *, const double *, double *),
							double     *TmpSol,
							double     *TmpFunc1,
							double     *TmpFunc2,
							double     *Solution,
							Boolean_t  *TimeStepReduced);
Boolean_t RK45UpdateSolution(const EntIndex_t NumSolVars, 
							 double *DelTime, 
							 double *StepDistance, 
							 const double GridSpacing,
							 const double MinCPDist, 
							 const Boolean_t ReverseDir, 
							 const double StepSizeFactor, 
							 const int StepSizeMultiple, 
							 const void *ClientData, 
							 Boolean_t (*Function)(const void *, const double *, double *), 
							 double *Solution, 
							 const double Tolerance, 
							 Boolean_t *TimeStepReduced);
Boolean_t RK2UpdateSolutionSurf(const EntIndex_t NumSolVars, 
								const double b2, 
								double *DelTime, 
								double *StepDistance, 
								const Boolean_t ReverseDir, 
								const void *ClientData, 
								Boolean_t (*Function)(const void *, const double *, double *), 
								double *Solution,
								const double Tolerance, 
								Boolean_t *TimeStepReduced);
Boolean_t RK45UpdateSolutionSurf(const EntIndex_t NumSolVars,
								 double *DelTime, 
								 double *StepDistance, 
								 const Boolean_t ReverseDir, 
								 const void *ClientData, 
								 Boolean_t (*Function)(const void *, 
								 const double *, double *), 
								 double *Solution, 
								 const double Tolerance, 
								 Boolean_t *TimeStepReduced);
Boolean_t NormalizeVector(double *DataArray);
#endif /* ODERUNGEKUTTA_H_ */
