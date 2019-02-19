/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNMatrix Dense module 
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_sunmatrix.h"


/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int          fails = 0;        /* counter for test failures  */
  sunindextype matrows, matcols; /* vector length              */
  N_Vector     x, y;             /* test vectors               */
  realtype     *xdata, *ydata;   /* pointers to vector data    */
  SUNMatrix    A, I;             /* test matrices              */
  realtype     *Adata, *Idata;   /* pointers to matrix data    */
  int          print_timing, square;
  sunindextype i, j, m, n;

  /* check input and set vector length */
  if (argc < 4){
    printf("ERROR: THREE (3) Input required: matrix rows, matrix cols, print timing \n");
    return(-1);
  }
  
  matrows = atol(argv[1]); 
  if (matrows <= 0) {
    printf("ERROR: number of rows must be a positive integer \n");
    return(-1); 
  }
  
  matcols = atol(argv[2]); 
  if (matcols <= 0) {
    printf("ERROR: number of cols must be a positive integer \n");
    return(-1); 
  }

  print_timing = atoi(argv[3]);
  SetTiming(print_timing);
  
  square = (matrows == matcols) ? 1 : 0;
  printf("\nDense matrix test: size %ld by %ld\n\n",
         (long int) matrows, (long int) matcols);

  /* Initialize vectors and matrices to NULL */
  x = NULL;
  y = NULL;
  A = NULL;
  I = NULL;
  
  /* Create vectors and matrices */
  x = N_VNew_Serial(matcols);
  y = N_VNew_Serial(matrows);
  A = SUNDenseMatrix(matrows, matcols);
  I = NULL;
  if (square)
    I = SUNDenseMatrix(matrows, matcols);
  
  /* Fill matrices and vectors */
  Adata = SUNDenseMatrix_Data(A);
  for(j=0; j < matcols; j++) {
    for(i=0; i < matrows; i++) {
      Adata[j*matrows + i] = (j+1)*(i+j);
    }
  }

  if (square) {
    Idata = SUNDenseMatrix_Data(I);
    for(i=0, j=0; i < matrows; i++, j++) {
      Idata[j*matrows + i] = ONE;
    }
  }

  xdata = N_VGetArrayPointer(x);
  for(i=0; i < matcols; i++) {
    xdata[i] = ONE / (i+1);
  }

  ydata = N_VGetArrayPointer(y);
  for(i=0; i < matrows; i++) {
    m = i;
    n = m + matcols - 1;
    ydata[i] = HALF*(n+1-m)*(n+m);
  }
    
  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(A, SUNMATRIX_DENSE, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  if (square) {
    fails += Test_SUNMatScaleAdd(A, I, 0);
    fails += Test_SUNMatScaleAddI(A, I, 0);
  }
  fails += Test_SUNMatMatvec(A, x, y, 0);
  fails += Test_SUNMatSpace(A, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
    printf("\nA =\n");
    SUNDenseMatrix_Print(A,stdout);
    if (square) {
      printf("\nI =\n");
      SUNDenseMatrix_Print(I,stdout);
    }
    printf("\nx =\n");
    N_VPrint_Serial(x);
    printf("\ny =\n");
    N_VPrint_Serial(y);
  } else {
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");
  }

  /* Free vectors and matrices */
  N_VDestroy_Serial(x);
  N_VDestroy_Serial(y);
  SUNMatDestroy(A);
  if (square)
    SUNMatDestroy(I);

  return(fails);
}

/* ----------------------------------------------------------------------
 * Check matrix
 * --------------------------------------------------------------------*/
int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  realtype *Adata, *Bdata;
  sunindextype Aldata, Bldata;
  sunindextype i;
  
  /* get data pointers */
  Adata = SUNDenseMatrix_Data(A);
  Bdata = SUNDenseMatrix_Data(B);

  /* get and check data lengths */
  Aldata = SUNDenseMatrix_LData(A);
  Bldata = SUNDenseMatrix_LData(B);

  if (Aldata != Bldata) {
    printf(">>> ERROR: check_matrix: Different data array lengths \n");
    return(1);
  }
  
  /* compare data */
  for(i=0; i < Aldata; i++){
    failure += FNEQ(Adata[i], Bdata[i], tol);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure = 0;
  realtype *Adata;
  sunindextype Aldata;
  sunindextype i;
  
  /* get data pointer */
  Adata = SUNDenseMatrix_Data(A);

  /* compare data */
  Aldata = SUNDenseMatrix_LData(A);
  for(i=0; i < Aldata; i++){
    failure += FNEQ(Adata[i], val, tol);
  }

  if (failure > ZERO) {
    printf("Check_matrix_entry failures:\n");
    for(i=0; i < Aldata; i++)
      if (FNEQ(Adata[i], val, tol) != 0)
        printf("  Adata[%ld] = %g != %g (err = %g)\n", (long int) i,
               Adata[i], val, SUNRabs(Adata[i]-val));
  }
  
  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_vector(N_Vector x, N_Vector y, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* get vector data */
  xdata = N_VGetArrayPointer(x);
  ydata = N_VGetArrayPointer(y);

  /* check data lengths */
  xldata = N_VGetLength_Serial(x);
  yldata = N_VGetLength_Serial(y);

  if (xldata != yldata) {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return(1);
  }

  /* check vector data */
  for(i=0; i < xldata; i++)
    failure += FNEQ(xdata[i], ydata[i], tol);

  if (failure > ZERO) {
    printf("Check_vector failures:\n");
    for(i=0; i < xldata; i++)
      if (FNEQ(xdata[i], ydata[i], tol) != 0)
        printf("  xdata[%ld] = %g != %g (err = %g)\n", (long int) i,
               xdata[i], ydata[i], SUNRabs(xdata[i]-ydata[i]));
  }
  
  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(SUNMatrix A)
{
  realtype *Adata = SUNDenseMatrix_Data(A);
  if (Adata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  if (SUNDenseMatrix_Rows(A) == SUNDenseMatrix_Columns(A))
    return SUNTRUE;
  else
    return SUNFALSE;
}
