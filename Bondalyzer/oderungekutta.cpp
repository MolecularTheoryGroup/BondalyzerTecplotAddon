/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#include "TECADDON.h"
#include "GLOBAL.h"
#include "GEOMTOOLS.h"
#include "ODERUNGEKUTTA.h"
#include "ARRLIST.h"
#include "SURFELEMMAP.h"
#include "NORMALS.h"
#include "ZONEVARINFO.h"
#include "CRITPOINTS.h"
#include "GRADPATH.h"



/*
 * RK2UpdateSolution: Update the solution, by one time steps, of the ordinary
 * differential equation,
 *
 *   d(Solution)/dt = Function(Solution)
 *
 *  using a 2-step Runge Kutta method. The user can specify which method
 *  to use via the parameters b2.
 *
 *   NumSolVars   Number of variables in Solution array.
 *
 *   b2           Parameter describing the Runge-Kutta method.
 *                   1.0 is midpoint method.
 *                   0.75 is Atkinson method
 *
 *   DelTime      TimeStep
 *   ReverseDir   If TRUE, move in direction opposite the vector returned from Function
 *                if FALSE, move in direction of the vector returned from Function
 *   Function     Pointer to function in ODE.
 *   Solution     Solution vector in ODE
 *
 *   Boolean_t Function(const double *Solution, double *FunctionArray) returns
 *     TRUE if a value for the FunctionArray is found, and FALSE if it isn't.
 */
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
							Boolean_t  *TimeStepReduced)
{
	Boolean_t IsOk = TRUE;
	double    b1  = 1.0 - b2;
	double    c2  = 0.5 / b2;
	double    a21 = c2;
	double	  DelTimeNew;

	REQUIRE(NumSolVars > 0);
	REQUIRE(b2 > 0.0 && b2 <= 1.0);
	REQUIRE(ABS(b1 + b2 - 1.0) < 1.0e-5);  /* Consistent solution */
	REQUIRE(VALID_BOOLEAN(ReverseDir));
	REQUIRE(VALID_REF(ClientData));
	REQUIRE(VALID_REF(Function));
	REQUIRE(VALID_REF(TmpSol));
	REQUIRE(VALID_REF(TmpFunc1));
	REQUIRE(VALID_REF(TmpFunc2));
	REQUIRE(VALID_REF(Solution));
	REQUIRE(VALID_REF(TimeStepReduced));
	REQUIRE(VALID_REF(DelTime));

	/* Assume time step does not need to be reduced. */
	*TimeStepReduced = FALSE;

	/* Compute tmpSol from first stage of RK*/
	IsOk = Function(ClientData, Solution, TmpFunc1);
	if (IsOk)
	{
		int ii;
		if (ReverseDir)
		{
			for (ii = 0; ii < NumSolVars; ii++)
				TmpSol[ii] = Solution[ii] - a21 * *DelTime * TmpFunc1[ii];
		}
		else
		{
			for (ii = 0; ii < NumSolVars; ii++)
				TmpSol[ii] = Solution[ii] + a21 * *DelTime * TmpFunc1[ii];
		}
	}

	/* Evaluate function at TmpSol. If function increases dramatically, recommend reduced time step. */
	if (IsOk)
		IsOk = Function(ClientData, TmpSol, TmpFunc2);
	if (IsOk)
	{
		/* TODO - Generalize to other CFL */
		DelTimeNew = 0.25 / MAX(ABS(TmpFunc2[0]) + ABS(TmpFunc2[1]) + ABS(TmpFunc2[2]), 1.0e-9);
		if (DelTimeNew < 0.25 * *DelTime)
		{
			*TimeStepReduced = TRUE;
			*DelTime = DelTimeNew;
		}
	}

	/* Compute updated solution from second stage of RK */
	if (IsOk)
	{
		int ii;
		if (ReverseDir)
		{
			for (ii = 0; ii < NumSolVars; ii++)
				Solution[ii] = Solution[ii] - *DelTime * (b1 * TmpFunc1[ii]
			+ b2 * TmpFunc2[ii]);
		}
		else
		{
			for (ii = 0; ii < NumSolVars; ii++)
				Solution[ii] = Solution[ii] + *DelTime * (b1 * TmpFunc1[ii]
			+ b2 * TmpFunc2[ii]);
		}
	}

	ENSURE(VALID_BOOLEAN(IsOk));
	return(IsOk);

} /* cmludcmp() */


/*
 * RK45UpdateSolution: Update the solution, by one time steps, of the ordinary
 * differential equation,
 *
 *   d(Solution)/dt = Function(Solution)
 *
 *  using a 4-5 Cash-Karp Runge-Kutta Fehlberg method. 
 *
 *   NumSolVars   Number of variables in Solution array.
 *
 *   b2           Parameter describing the Runge-Kutta method.
 *                   1.0 is midpoint method.
 *                   0.75 is Atkinson method
 *
 *   DelTime      TimeStep
 *   ReverseDir   If TRUE, move in direction opposite the vector returned from Function
 *                if FALSE, move in direction of the vector returned from Function
 *   Function     Pointer to function in ODE.
 *   Solution     Solution vector in ODE
 *
 *   Boolean_t Function(const double *Solution, double *FunctionArray) returns
 *     TRUE if a value for the FunctionArray is found, and FALSE if it isn't.
 */
Boolean_t RK45UpdateSolution(const EntIndex_t	NumSolVars,
                            double				*DelTime,
							const Boolean_t		ReverseDir,
                            const void			*ClientData,
                            Boolean_t (*Function)(const void *, const double *, double *),
                            double				*Solution,
                            Boolean_t			*TimeStepReduced)
{
    Boolean_t IsOk = TRUE;

	double	DelTimeNew;

	REQUIRE(NumSolVars > 0);
	REQUIRE(VALID_REF(DelTime));
    REQUIRE(VALID_BOOLEAN(ReverseDir));
    REQUIRE(VALID_REF(ClientData));
    REQUIRE(VALID_REF(Function));
    REQUIRE(VALID_REF(Solution));
    REQUIRE(VALID_REF(TimeStepReduced));

	//	Constants for the 4th order estimate
	const double	aa[] = {37.0/378.0,
							0.0,
							250.0/621.0,
							125.0/594.0,
							0.0,
							512.0/1771.0};
	//	Constants for the 5th order estimate
	const double	bb[] = {2825.0/27648.0,
							0.0,
							18575.0/48384.0,
							13525.0/55296.0,
							277.0/14336.0,
							1.0/4.0};
	//	Constants for the x components of k_1-6
	//	(*k1x and k1y are 1, so not defined here.)
	const double	k2x = 1.0/5.0,
					k3x = 3.0/10.0,
					k4x = 3.0/5.0,
					k5x = 1.0,
					k6x = 7.0/8.0;
	//	Constants for the y components of k_1-6
	const double	k2y		=	1.0/5.0;
	const double	k3y[]	=	{3.0/40.0,
								9.0/40.0},
					k4y[]	=	{3.0/10.0,
								-9.0/10.0,
								6.0/5.0},
					k5y[]	=	{-11.0/54.0,
								5.0/2.0,
								-70.0/27.0,
								35.0/27.0},
					k6y[]	=	{1631.0/55296.0,
								175.0/512.0,
								575.0/13824.0,
								44275.0/110592.0,
								253.0/4096.0};

	//	Ideally the length of the arrays would be variable,
	//	but we're always dealing with 2 or 3 dimensions, so that
	//	flexibility isn't needed.
	double	TempSol[3],
			k1[3], k2[3], k3[3], 
			k4[3], k5[3], k6[3],
			Slope4[3], Slope5[3],
			Sol4[3], Sol5[3];
	//	Absolute error between 4th and 5th order approximations
	//double	SolError = 0.0;
	int		i;

	*TimeStepReduced = FALSE;
	
	//	Implement the 4-5 Cash-Karp Runge-Kutta Fehlberg method

	//	Compute k1
	//	If Function fails here then the grad path has left the system.
	IsOk = Function(ClientData, Solution, k1);

	//	Compute k2
	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - k2y*k1[i]*(*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + k2y*k1[i]*(*DelTime);

		IsOk = Function(ClientData, TempSol, k2);
	}

	//	Compute k3
	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - (k3y[0]*k1[i] + k3y[1]*k2[i]) * (*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + (k3y[0]*k1[i] + k3y[1]*k2[i]) * (*DelTime);

		IsOk = Function(ClientData, TempSol, k3);
	}

	//	Compute k4
	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - (k4y[0]*k1[i] + k4y[1]*k2[i] + k4y[2]*k3[i]) * (*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + (k4y[0]*k1[i] + k4y[1]*k2[i] + k4y[2]*k3[i]) * (*DelTime);

		IsOk = Function(ClientData, TempSol, k4);
	}

	//	Compute k5
	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - (k5y[0]*k1[i] + k5y[1]*k2[i] +
								k5y[2]*k3[i] + k5y[3]*k4[i]) * (*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + (k5y[0]*k1[i] + k5y[1]*k2[i] +
								k5y[2]*k3[i] + k5y[3]*k4[i]) * (*DelTime);

		IsOk = Function(ClientData, TempSol, k5);
	}

	//	Compute k6
	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - (k6y[0]*k1[i] + k6y[1]*k2[i] +
								k6y[2]*k3[i] + k6y[3]*k4[i] + k6y[4]*k5[i]) * (*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + (k6y[0]*k1[i] + k6y[1]*k2[i] +
								k6y[2]*k3[i] + k6y[3]*k4[i] + k6y[4]*k5[i]) * (*DelTime);

		IsOk = Function(ClientData, TempSol, k6);
	}

	//	Compute "slopes" of the 4th and 5th order approximations
	if (IsOk)
		for (i = 0; i < NumSolVars ; i++)
		{
			Slope4[i] = aa[0]*k1[i] + aa[2]*k3[i] + aa[3]*k4[i] + aa[5]*k6[i];
			Slope5[i] = bb[0]*k1[i] + bb[2]*k3[i] + bb[3]*k4[i] + bb[4]*k5[i] + bb[5]*k6[i];
		}

	//	Compute the 4th and 5th order approximations

	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
			{
				Sol4[i] = Solution[i] - Slope4[i] * (*DelTime);
				Sol5[i] = Solution[i] - Slope5[i] * (*DelTime);
			}
		else
			for (i = 0 ; i < NumSolVars ; i++)
			{
				Sol4[i] = Solution[i] + Slope4[i] * (*DelTime);
				Sol5[i] = Solution[i] + Slope5[i] * (*DelTime);
			}

		////	Compute absolute error between 4th and 5th order approximations
		//for (i = 0 ; i < NumSolVars ; i++)
		//	SolError = MAX(SolError, ABS(Sol4[i] - Sol5[i]));

		DelTimeNew = 0.25 / MAX(ABS(Slope5[0]) + ABS(Slope5[1]) + ABS(Slope5[2]), 1.0e-9);
		if (DelTimeNew < 0.25 * *DelTime)
		{
			*TimeStepReduced = TRUE;
			*DelTime = DelTimeNew;
		}
		
		//	Check if error is in tolerance and adjust
		//if (SolError <= Tolerance)
		if (!*TimeStepReduced)
		{
			/**DelTime = 0.9 * *DelTime * Tolerance / SolError;
			*TimeStepReduced = FALSE;*/
			
			////	Step was successful, so update solution and find 
			////	distance traversed
			//XYZ_s	InitPoint, EndPoint;

			//InitPoint.X = Solution[0];
			//InitPoint.Y = Solution[1];
			//InitPoint.Z = Solution[2];
			//EndPoint.X = Sol5[0];
			//EndPoint.Y = Sol5[1];
			//EndPoint.Z = Sol5[2];
			//*StepDistance = sqrt(DistanceSquaredXYZ(EndPoint, InitPoint));

			for (i = 0 ; i < NumSolVars ; i++)
				Solution[i] = Sol5[i];
		}
	}

	/*if (!IsOk)
	{
		*DelTime *= 0.1;
		*TimeStepReduced = TRUE;
		IsOk = TRUE;
	}*/

    ENSURE(VALID_BOOLEAN(IsOk));
    return(IsOk);

}	//	Boolean_t RK45UpdateSolution()


/*
 * RK2UpdateSolutionSurf: Update the solution, by one time steps, of the ordinary
 * differential equation,
 *
 *   d(Solution)/dt = Function(Solution)
 *
 *  using a 2-step Runge Kutta method. The user can specify which method
 *  to use via the parameters b2. 
 *	project back to the surface during each step, that way the gradients at
 *	at the surface points are used
 *
 *   NumSolVars   Number of variables in Solution array.
 *
 *   b2           Parameter describing the Runge-Kutta method.
 *                   1.0 is midpoint method.
 *                   0.75 is Atkinson method
 *
 *   DelTime      TimeStep
 *   ReverseDir   If TRUE, move in direction opposite the vector returned from Function
 *                if FALSE, move in direction of the vector returned from Function
 *   Function     Pointer to function in ODE.
 *   Solution     Solution vector in ODE
 *
 *   Boolean_t Function(const double *Solution, double *FunctionArray) returns
 *     TRUE if a value for the FunctionArray is found, and FALSE if it isn't.
 */
Boolean_t RK2UpdateSolutionSurf(const EntIndex_t  NumSolVars,
                            const double      b2,
							double			*DelTime,
							double				*StepDistance,
                            const Boolean_t   ReverseDir,
                            const void       *ClientData,
                            Boolean_t (*Function)(const void *, const double *, double *),
                            double     *Solution,
							const double	Tolerance,
                            Boolean_t  *TimeStepReduced)
{
    Boolean_t IsOk = TRUE;
	double    b1  = 1.0 - b2;
	double    c2  = 0.5 / b2;
	double    a21 = c2;
	double	  DelTimeNew, CheckDistanceSqr, SolError = 0.0;
	double		TmpFunc1[3], TmpFunc2[3], TmpSol1[3], TmpSol2[3];
	XYZ_s		Pt1, Pt2, SurfPt;
    int		i;

    REQUIRE(NumSolVars > 0);
    REQUIRE(b2 > 0.0 && b2 <= 1.0);
    REQUIRE(ABS(b1 + b2 - 1.0) < 1.0e-5);  /* Consistent solution */
    REQUIRE(VALID_BOOLEAN(ReverseDir));
    REQUIRE(VALID_REF(ClientData));
    REQUIRE(VALID_REF(Function));
    REQUIRE(VALID_REF(Solution));
    REQUIRE(VALID_REF(TimeStepReduced));

    /* Assume time step does not need to be reduced. */
    *TimeStepReduced = FALSE;

	Pt1.X = Solution[0];
	Pt1.Y = Solution[1];
	Pt1.Z = Solution[2];

    /* Compute tmpSol from first stage of RK*/
    IsOk = Function(ClientData, Solution, TmpFunc1);

	if (IsOk) IsOk = NormalizeVector(TmpFunc1);

    if (IsOk)
    {
        if (ReverseDir)
            for (i = 0; i < NumSolVars; i++)
                TmpSol1[i] = Solution[i] - a21 * *DelTime * TmpFunc1[i];
        else
            for (i = 0; i < NumSolVars; i++)
                TmpSol1[i] = Solution[i] + a21 * *DelTime * TmpFunc1[i];
    }

	/*Pt2.X = TmpSol1[0];
	Pt2.Y = TmpSol1[1];
	Pt2.Z = TmpSol1[2];

	CheckDistanceSqr = DistanceSquaredXYZ(Pt1,Pt2);
	CHECK(CheckDistanceSqr > 0);

	if (IsOk)
		IsOk = GradPathProjectPointToSurf((ZoneVarInfo_pa)ClientData, CheckDistanceSqr, Pt2, &SurfPt);
	if (IsOk)
	{
		TmpSol1[0] = SurfPt.X;
		TmpSol1[1] = SurfPt.Y;
		TmpSol1[2] = SurfPt.Z;
	}*/

    /* Evaluate function at TmpSol. If function increases dramatically, recommend reduced time step. */
    if (IsOk)
        IsOk = Function(ClientData, TmpSol1, TmpFunc2);
  //  if (IsOk)
  //  {
  //      /* TODO - Generalize to other CFL */
  //      DelTimeNew = 0.25 / MAX(ABS(TmpFunc2[0]) + ABS(TmpFunc2[1]) + ABS(TmpFunc2[2]), 1.0e-9);
  //      if (DelTimeNew < 0.25 * *DelTime)
		//{
		//	*DelTime *= 0.5;
		//	*TimeStepReduced = TRUE;
		//}
  //  }

    /* Compute updated solution from second stage of RK */
	if (IsOk) IsOk = NormalizeVector(TmpFunc2);

    if (IsOk/* && !*TimeStepReduced*/)
    {
        if (ReverseDir)
            for (i = 0; i < NumSolVars; i++)
                TmpSol2[i] = Solution[i] - *DelTime * (b1 * TmpFunc1[i]
                                                         + b2 * TmpFunc2[i]);
        else
            for (i = 0; i < NumSolVars; i++)
                TmpSol2[i] = Solution[i] + *DelTime * (b1 * TmpFunc1[i]
                                                         + b2 * TmpFunc2[i]);

		/*Pt2.X = TmpSol2[0];
		Pt2.Y = TmpSol2[1];
		Pt2.Z = TmpSol2[2];

		CheckDistanceSqr = DistanceSquaredXYZ(Pt1,Pt2);
		CHECK(CheckDistanceSqr > 0);

		if (IsOk)
			IsOk = GradPathProjectPointToSurf((ZoneVarInfo_pa)ClientData, CheckDistanceSqr, Pt2, &SurfPt);
		if (IsOk)
		{
			TmpSol2[0] = SurfPt.X;
			TmpSol2[1] = SurfPt.Y;
			TmpSol2[2] = SurfPt.Z;
		}*/
    }

	if ( 1 || !*TimeStepReduced)
	{
		for (i = 0 ; i < NumSolVars ; i++)
			Solution[i] = TmpSol2[i];
	}

	//for (i = 0 ; i < NumSolVars ; i++)
	//	SolError = MAX(SolError, ABS(TmpSol1[i] - TmpSol2[i]));

	//if (SolError <= Tolerance)
	//{
	//	*DelTime = 0.9 * *DelTime * Tolerance / SolError;
	//	*TimeStepReduced = FALSE;

	//	//	Step was successful, so update solution 

	//	*StepDistance = sqrt(DistanceSquaredXYZ(SurfPt, Pt1));

	//	for (i = 0 ; i < NumSolVars ; i++)
	//		Solution[i] = TmpSol2[i];
	//}
	//else
	//{
	//	*DelTime *= 0.5;
	//	*TimeStepReduced = TRUE;
	//}



    ENSURE(VALID_BOOLEAN(IsOk));
    return(IsOk);

}	//	Boolean_t RK2UpdateSolutionSurf()


/*
 * RK45UpdateSolution: Update the solution, by one time steps, of the ordinary
 * differential equation, 
 *
 *   d(Solution)/dt = Function(Solution)
 *
 *  using a 4-5 Cash-Karp Runge-Kutta Fehlberg method. 
 *	project back to the surface during each step, that way the gradients at
 *	at the surface points are used
 *
 *   NumSolVars   Number of variables in Solution array.
 *
 *   b2           Parameter describing the Runge-Kutta method.
 *                   1.0 is midpoint method.
 *                   0.75 is Atkinson method
 *
 *   DelTime      TimeStep
 *   ReverseDir   If TRUE, move in direction opposite the vector returned from Function
 *                if FALSE, move in direction of the vector returned from Function
 *   Function     Pointer to function in ODE.
 *   Solution     Solution vector in ODE
 *
 *   Boolean_t Function(const double *Solution, double *FunctionArray) returns
 *     TRUE if a value for the FunctionArray is found, and FALSE if it isn't.
 */
Boolean_t RK45UpdateSolutionSurf(const EntIndex_t	NumSolVars,
                            double				*DelTime,
							double				*StepDistance,
                            const Boolean_t		ReverseDir,
                            const void			*ClientData,
                            Boolean_t (*Function)(const void *, const double *, double *),
                            double				*Solution,
							const double		Tolerance,
                            Boolean_t			*TimeStepReduced)
{
	Boolean_t IsOk = TRUE;

	REQUIRE(NumSolVars > 0);
	REQUIRE(VALID_REF(DelTime));
	REQUIRE(VALID_BOOLEAN(ReverseDir));
	REQUIRE(VALID_REF(ClientData));
	REQUIRE(VALID_REF(Function));
	REQUIRE(VALID_REF(Solution));
	REQUIRE(VALID_REF(TimeStepReduced));

	//	Constants for the 4th order estimate
	const double	aa[] = {37.0/378.0,
							0.0,
							250.0/621.0,
							125.0/594.0,
							0.0,
							512.0/1771.0};
	//	Constants for the 5th order estimate
	const double	bb[] = {2825.0/27648.0,
							0.0,
							18575.0/48384.0,
							13525.0/55296.0,
							277.0/14336.0,
							1.0/4.0};
	//	Constants for the x components of k_1-6
	//	(*k1x and k1y are 1, so not defined here.)
	const double	k2x = 1.0/5.0,
					k3x = 3.0/10.0,
					k4x = 3.0/5.0,
					k5x = 1.0,
					k6x = 7.0/8.0;
	//	Constants for the y components of k_1-6
	const double	k2y		=	1.0/5.0;
	const double	k3y[]	=	{3.0/40.0,
								9.0/40.0},
					k4y[]	=	{3.0/10.0,
								-9.0/10.0,
								6.0/5.0},
					k5y[]	=	{-11.0/54.0,
								5.0/2.0,
								-70.0/27.0,
								35.0/27.0},
					k6y[]	=	{1631.0/55296.0,
								175.0/512.0,
								575.0/13824.0,
								44275.0/110592.0,
								253.0/4096.0};

	/* Assume time step does not need to be reduced. */
	*TimeStepReduced = FALSE;
	//	Ideally the length of the arrays would be variable,
	//	but we're always dealing with 2 or 3 dimensions, so that
	//	flexibility isn't needed.
	double	TempSol[3],
		k1[3], k2[3], k3[3], 
		k4[3], k5[3], k6[3],
		Slope4[3], Slope5[3],
		Sol4[3], Sol5[3];
	//	Absolute error between 4th and 5th order approximations
	double	SolError = 0.0, CheckDistanceSqr;
	int		i;
	XYZ_s	Pt1, Pt2, SurfPt;

	Pt1.X = Solution[0];
	Pt1.Y = Solution[1];
	Pt1.Z = Solution[2];

	//	Implement the 4-5 Cash-Karp Runge-Kutta Fehlberg method

	//	Compute k1
	IsOk = Function(ClientData, Solution, k1);

	//	Compute k2
	if (IsOk) IsOk = NormalizeVector(k1);

	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - k2y*k1[i]*(*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + k2y*k1[i]*(*DelTime);

		/*Pt2.X = TempSol[0];
		Pt2.Y = TempSol[1];
		Pt2.Z = TempSol[2];

		CheckDistanceSqr = DistanceSquaredXYZ(Pt1,Pt2);
		CHECK(CheckDistanceSqr > 0);

		if (IsOk)
			IsOk = GradPathProjectPointToSurf((ZoneVarInfo_pa)ClientData, CheckDistanceSqr, Pt2, &SurfPt);
		if (IsOk)
		{
			TempSol[0] = SurfPt.X;
			TempSol[1] = SurfPt.Y;
			TempSol[2] = SurfPt.Z;
		}*/

		IsOk = Function(ClientData, TempSol, k2);
	}

	//	Compute k3
	if (IsOk) IsOk = NormalizeVector(k2);

	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - (k3y[0]*k1[i] + k3y[1]*k2[i]) * (*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + (k3y[0]*k1[i] + k3y[1]*k2[i]) * (*DelTime);

		/*Pt2.X = TempSol[0];
		Pt2.Y = TempSol[1];
		Pt2.Z = TempSol[2];

		CheckDistanceSqr = DistanceSquaredXYZ(Pt1,Pt2);
		CHECK(CheckDistanceSqr > 0);

		if (IsOk)
			IsOk = GradPathProjectPointToSurf((ZoneVarInfo_pa)ClientData, CheckDistanceSqr, Pt2, &SurfPt);
		if (IsOk)
		{
			TempSol[0] = SurfPt.X;
			TempSol[1] = SurfPt.Y;
			TempSol[2] = SurfPt.Z;
		}*/

		IsOk = Function(ClientData, TempSol, k3);
	}

	//	Compute k4
	if (IsOk) IsOk = NormalizeVector(k3);

	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - (k4y[0]*k1[i] + k4y[1]*k2[i] + k4y[2]*k3[i]) * (*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + (k4y[0]*k1[i] + k4y[1]*k2[i] + k4y[2]*k3[i]) * (*DelTime);

		/*Pt2.X = TempSol[0];
		Pt2.Y = TempSol[1];
		Pt2.Z = TempSol[2];

		CheckDistanceSqr = DistanceSquaredXYZ(Pt1,Pt2);
		CHECK(CheckDistanceSqr > 0);

		if (IsOk)
			IsOk = GradPathProjectPointToSurf((ZoneVarInfo_pa)ClientData, CheckDistanceSqr, Pt2, &SurfPt);
		if (IsOk)
		{
			TempSol[0] = SurfPt.X;
			TempSol[1] = SurfPt.Y;
			TempSol[2] = SurfPt.Z;
		}*/

		IsOk = Function(ClientData, TempSol, k4);
	}

	//	Compute k5
	if (IsOk) IsOk = NormalizeVector(k4);

	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - (k5y[0]*k1[i] + k5y[1]*k2[i] +
				k5y[2]*k3[i] + k5y[3]*k4[i]) * (*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + (k5y[0]*k1[i] + k5y[1]*k2[i] +
				k5y[2]*k3[i] + k5y[3]*k4[i]) * (*DelTime);

		/*Pt2.X = TempSol[0];
		Pt2.Y = TempSol[1];
		Pt2.Z = TempSol[2];

		CheckDistanceSqr = DistanceSquaredXYZ(Pt1,Pt2);
		CHECK(CheckDistanceSqr > 0);

		if (IsOk)
			IsOk = GradPathProjectPointToSurf((ZoneVarInfo_pa)ClientData, CheckDistanceSqr, Pt2, &SurfPt);
		if (IsOk)
		{
			TempSol[0] = SurfPt.X;
			TempSol[1] = SurfPt.Y;
			TempSol[2] = SurfPt.Z;
		}*/

		IsOk = Function(ClientData, TempSol, k5);
	}

	//	Compute k6
	if (IsOk) IsOk = NormalizeVector(k5);

	if (IsOk)
	{
		if (ReverseDir)
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] - (k6y[0]*k1[i] + k6y[1]*k2[i] +
				k6y[2]*k3[i] + k6y[3]*k4[i] + k6y[4]*k5[i]) * (*DelTime);
		else
			for (i = 0 ; i < NumSolVars ; i++)
				TempSol[i] = Solution[i] + (k6y[0]*k1[i] + k6y[1]*k2[i] +
				k6y[2]*k3[i] + k6y[3]*k4[i] + k6y[4]*k5[i]) * (*DelTime);

		/*Pt2.X = TempSol[0];
		Pt2.Y = TempSol[1];
		Pt2.Z = TempSol[2];

		CheckDistanceSqr = DistanceSquaredXYZ(Pt1,Pt2);
		CHECK(CheckDistanceSqr > 0);

		if (IsOk)
			IsOk = GradPathProjectPointToSurf((ZoneVarInfo_pa)ClientData, CheckDistanceSqr, Pt2, &SurfPt);
		if (IsOk)
		{
			TempSol[0] = SurfPt.X;
			TempSol[1] = SurfPt.Y;
			TempSol[2] = SurfPt.Z;
		}*/

		IsOk = Function(ClientData, TempSol, k6);
	}

	//	Compute "slopes" of the 4th and 5th order approximations
	//if (IsOk) IsOk = NormalizeVector(k6);

	if (IsOk)
		for (i = 0; i < NumSolVars ; i++)
		{
			Slope4[i] = aa[0]*k1[i] + aa[2]*k3[i] + aa[3]*k4[i] + aa[5]*k6[i];
			Slope5[i] = bb[0]*k1[i] + bb[2]*k3[i] + bb[3]*k4[i] + bb[4]*k5[i] + bb[5]*k6[i];
		}

		//	Compute the 4th and 5th order approximations

		if (IsOk) IsOk = NormalizeVector(Slope4);
		if (IsOk) IsOk = NormalizeVector(Slope5);

		if (IsOk)
		{
			if (ReverseDir)
				for (i = 0 ; i < NumSolVars ; i++)
				{
					Sol4[i] = Solution[i] - Slope4[i] * (*DelTime);
					Sol5[i] = Solution[i] - Slope5[i] * (*DelTime);
				}
			else
				for (i = 0 ; i < NumSolVars ; i++)
				{
					Sol4[i] = Solution[i] + Slope4[i] * (*DelTime);
					Sol5[i] = Solution[i] + Slope5[i] * (*DelTime);
				}

				/*Pt2.X = Sol4[0];
				Pt2.Y = Sol4[1];
				Pt2.Z = Sol4[2];

				CheckDistanceSqr = DistanceSquaredXYZ(Pt1,Pt2);
				CHECK(CheckDistanceSqr > 0);

				if (IsOk)
					IsOk = GradPathProjectPointToSurf((ZoneVarInfo_pa)ClientData, CheckDistanceSqr, Pt2, &SurfPt);
				if (IsOk)
				{
					Sol4[0] = SurfPt.X;
					Sol4[1] = SurfPt.Y;
					Sol4[2] = SurfPt.Z;
				}

				Pt2.X = Sol5[0];
				Pt2.Y = Sol5[1];
				Pt2.Z = Sol5[2];

				CheckDistanceSqr = DistanceSquaredXYZ(Pt1,Pt2);
				CHECK(CheckDistanceSqr > 0);

				if (IsOk)
					IsOk = GradPathProjectPointToSurf((ZoneVarInfo_pa)ClientData, CheckDistanceSqr, Pt2, &SurfPt);
				if (IsOk)
				{
					Sol5[0] = SurfPt.X;
					Sol5[1] = SurfPt.Y;
					Sol5[2] = SurfPt.Z;
				}*/

				//	Compute absolute error between 4th and 5th order approximations
				for (i = 0 ; i < NumSolVars ; i++)
					SolError = MAX(SolError, ABS(Sol4[i] - Sol5[i]));

				//	Check if error is in tolerance and adjust
				if (SolError <= Tolerance)
				{
					*DelTime = 0.9 * *DelTime * Tolerance / SolError;
					*TimeStepReduced = FALSE;

					//	Step was successful, so update solution and find 
					//	distance traversed

					*StepDistance = sqrt(DistanceSquaredXYZ(SurfPt, Pt1));

					for (i = 0 ; i < NumSolVars ; i++)
						Solution[i] = Sol5[i];
				}
				else
				{
					*DelTime *= 0.5;
					*TimeStepReduced = TRUE;
				}
		}

		if (!IsOk)
		{
			*DelTime *= 0.1;
			*TimeStepReduced = TRUE;
			IsOk = TRUE;
		}

		ENSURE(VALID_BOOLEAN(IsOk));
		return(IsOk);

}	//	Boolean_t RK45UpdateSolutionSurf()

/*
	Normalizes an array of values

	Returns true if successful, false if failed
*/
Boolean_t	NormalizeVector(double *DataArray)
{
	Boolean_t IsOk = TRUE;
	REQUIRE(VALID_REF(DataArray));
	
	int i;
	double Normalizer = 0.0;

	for (i = 0 ; i < 3 ; i++)
		Normalizer += DataArray[i] * DataArray[i];

	Normalizer = sqrt(Normalizer);

	IsOk = (Normalizer > 0);

	if (IsOk)
		for (i = 0 ; i < 3 ; i++)
			DataArray[i] /= Normalizer;

	ENSURE(VALID_BOOLEAN(IsOk));
	return IsOk;
}