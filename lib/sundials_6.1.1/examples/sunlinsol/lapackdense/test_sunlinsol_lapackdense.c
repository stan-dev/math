/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol LapackDense
 * module implementation.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_lapackdense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_sunlinsol.h"


/* ----------------------------------------------------------------------
 * SUNLinSol_LapackDense Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int             fails = 0;          /* counter for test failures  */
  sunindextype    cols, rows;         /* matrix columns, rows       */
  SUNLinearSolver LS;                 /* solver object              */
  SUNMatrix       A, B, I;            /* test matrices              */
  N_Vector        x, y, b;            /* test vectors               */
  int             print_timing;
  sunindextype    j, k;
  realtype        *colj, *xdata, *colIj;
  SUNContext      sunctx;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  /* check input and set matrix dimensions */
  if (argc < 3){
    printf("ERROR: TWO (2) Inputs required: matrix cols, print timing \n");
    return(-1);
  }

  cols = (sunindextype) atol(argv[1]);
  if (cols <= 0) {
    printf("ERROR: number of matrix columns must be a positive integer \n");
    return(-1);
  }

  rows = cols;

  print_timing = atoi(argv[2]);
  SetTiming(print_timing);

  printf("\nLapackDense linear solver test: size %ld\n\n",
         (long int) cols);

  /* Create matrices and vectors */
  A = SUNDenseMatrix(rows, cols, sunctx);
  B = SUNDenseMatrix(rows, cols, sunctx);
  I = SUNDenseMatrix(rows, cols, sunctx);
  x = N_VNew_Serial(cols, sunctx);
  y = N_VNew_Serial(cols, sunctx);
  b = N_VNew_Serial(cols, sunctx);

  /* Fill A matrix with uniform random data in [0,1/cols] */
  for (j=0; j<cols; j++) {
    colj = SUNDenseMatrix_Column(A, j);
    for (k=0; k<rows; k++)
      colj[k] = (realtype) rand() / (realtype) RAND_MAX / cols;
  }

  /* Create anti-identity matrix */
  j=cols-1;
  for (k=0; k<rows; k++) {
    colj = SUNDenseMatrix_Column(I,j);
    colj[k] = 1;
    j = j-1;
  }

  /* Add anti-identity to ensure the solver needs to do row-swapping */
  for (k=0; k<rows; k++){
    for(j=0; j<cols; j++){
      colj = SUNDenseMatrix_Column(A,j);
      colIj = SUNDenseMatrix_Column(I,j);
      colj[k]  = colj[k] + colIj[k];
   }
  }

  /* Fill x vector with uniform random data in [0,1] */
  xdata = N_VGetArrayPointer(x);
  for (j=0; j<cols; j++) {
    xdata[j] = (realtype) rand() / (realtype) RAND_MAX;
  }

  /* copy A and x into B and y to print in case of solver failure */
  SUNMatCopy(A, B);
  N_VScale(ONE, x, y);

  /* create right-hand side vector for linear solve */
  fails = SUNMatMatvec(A, x, b);
  if (fails) {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");
    return(1);
  }

  /* Create dense linear solver */
  LS = SUNLinSol_LapackDense(x, A, sunctx);

  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSetup(LS, A, 0);
  fails += Test_SUNLinSolSolve(LS, A, x, b, 100*UNIT_ROUNDOFF, SUNTRUE, 0);

  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
  fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_LAPACKDENSE, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
    printf("\nA (original) =\n");
    SUNDenseMatrix_Print(B,stdout);
    printf("\nA (factored) =\n");
    SUNDenseMatrix_Print(A,stdout);
    printf("\nx (original) =\n");
    N_VPrint_Serial(y);
    printf("\nx (computed) =\n");
    N_VPrint_Serial(x);
  } else {
    printf("SUCCESS: SUNLinSol module passed all tests \n \n");
  }

  /* Free solver, matrix and vectors */
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNMatDestroy(B);
  SUNMatDestroy(I);
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(b);

  SUNContext_Free(&sunctx);

  return(fails);
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_vector(N_Vector X, N_Vector Y, realtype tol)
{
  int failure = 0;
  sunindextype i, local_length;
  realtype *Xdata, *Ydata, maxerr;

  Xdata = N_VGetArrayPointer(X);
  Ydata = N_VGetArrayPointer(Y);
  local_length = N_VGetLength_Serial(X);

  /* check vector data */
  for(i=0; i < local_length; i++)
    failure += SUNRCompareTol(Xdata[i], Ydata[i], tol);

  if (failure > ZERO) {
    maxerr = ZERO;
    for(i=0; i < local_length; i++)
      maxerr = SUNMAX(SUNRabs(Xdata[i]-Ydata[i]), maxerr);
    printf("check err failure: maxerr = %g (tol = %g)\n",
	   maxerr, tol);
    return(1);
  }
  else
    return(0);
}

void sync_device()
{
}
