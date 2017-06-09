/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


#ifndef LINEARALGEBRA_H_
#define LINEARALGEBRA_H_

typedef struct _Matrix_s
{
    double **val;      /* appears to be 2D array of values */
    double  *raw_data; /* actual 1D array of dim1*dim2 values */
} Matrix_s;

void NullMatrix(Matrix_s *matrixptr);
void FreeMatrix(Matrix_s *matrixptr);
Matrix_s CreateMatrix(LgIndex_t dim1, LgIndex_t dim2);

Boolean_t cmludcmp(Matrix_s   matrix,
                   LgIndex_t  ndim,
                   LgIndex_t *indx,
                   double    *TmpSpace);

Boolean_t cmlubksb(Matrix_s   matrix,
                   LgIndex_t  ndim,
                   LgIndex_t *indx,
                   double    *bb);

#endif /* LINEARALGEBRA_H_ */
