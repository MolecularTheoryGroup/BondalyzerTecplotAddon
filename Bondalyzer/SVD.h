/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/

#ifndef SVD_H_
#define SVD_H_

double *SVDAllocVector(int nl, int nh);
void SVDFreeVector(double *v, int nl, int nh);
double **SVDAllocMatrix(int nrl, int nrh, int ncl, int nch);
void SVDFreeMatrix(double **m, int nrl, int nrh, int ncl, int nch);
void ComputeSVD(double **a, int m, int n, double w[], double **v);
void BackSubstituteSVD(double **u, double w[], double **v, int m, int n, double b[], double x[]);


#endif /* SVD_H_ */