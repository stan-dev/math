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
 * This is the testing routine to check the SUNLinSol KLU module
 * implementation.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_sunlinsol.h"


/* ----------------------------------------------------------------------
 * SUNLinSol_KLU Linear Solver Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int             fails = 0;          /* counter for test failures  */
  sunindextype    N;                  /* matrix columns, rows       */
  SUNLinearSolver LS;                 /* linear solver object       */
  SUNMatrix       A, B;               /* test matrices              */
  N_Vector        x, y, b;            /* test vectors               */
  realtype        *matdata, *xdata;
  int             mattype, print_timing;
  sunindextype    i, j, k;
  sun_klu_symbolic *symbolic;
  sun_klu_numeric  *numeric;
  sun_klu_common   *common;
  SUNContext      sunctx;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  /* check input and set matrix dimensions */
  if (argc < 4){
    printf("ERROR: THREE (3) Inputs required: matrix size, matrix type (0/1), print timing \n");
    return(-1);
  }

  N = (sunindextype) atol(argv[1]);
  if (N <= 0) {
    printf("ERROR: matrix size must be a positive integer \n");
    return(-1);
  }

  mattype = atoi(argv[2]);
  if ((mattype != 0) && (mattype != 1)) {
    printf("ERROR: matrix type must be 0 or 1 \n");
    return(-1);
  }
  mattype = (mattype == 0) ? CSC_MAT : CSR_MAT;

  print_timing = atoi(argv[3]);
  SetTiming(print_timing);

  printf("\nKLU linear solver test: size %ld, type %i\n\n",
         (long int) N, mattype);

  /* Create matrices and vectors */
  B = SUNDenseMatrix(N, N, sunctx);
  x = N_VNew_Serial(N, sunctx);
  y = N_VNew_Serial(N, sunctx);
  b = N_VNew_Serial(N, sunctx);

  /* Fill matrix with uniform random data in [0,1/N] */
  for (k=0; k<5*N; k++) {
    i = rand() % N;
    j = rand() % N;
    matdata = SUNDenseMatrix_Column(B,j);
    matdata[i] = (realtype) rand() / (realtype) RAND_MAX / N;
  }

  /* Add identity to matrix */
  fails = SUNMatScaleAddI(ONE, B);
  if (fails) {
    printf("FAIL: SUNLinSol SUNMatScaleAddI failure\n");
    return(1);
  }

  /* Fill x vector with uniform random data in [0,1] */
  xdata = N_VGetArrayPointer(x);
  for (i=0; i<N; i++)
    xdata[i] = (realtype) rand() / (realtype) RAND_MAX;

  /* Create sparse matrix from dense, and destroy B */
  A = SUNSparseFromDenseMatrix(B, ZERO, mattype);
  SUNMatDestroy(B);

  /* copy x into y to print in case of solver failure */
  N_VScale(ONE, x, y);

  /* create right-hand side vector for linear solve */
  fails = SUNMatMatvec(A, x, b);
  if (fails) {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");
    return(1);
  }

  /* Create KLU linear solver */
  LS = SUNLinSol_KLU(x, A, sunctx);

  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSetup(LS, A, 0);
  fails += Test_SUNLinSolSolve(LS, A, x, b, 1000*UNIT_ROUNDOFF, SUNTRUE, 0);

  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
  fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_KLU, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);

  /* Test 'Get' routines */
  symbolic = SUNLinSol_KLUGetSymbolic(LS);
  if (symbolic->n != N) {
    printf("FAIL: SUNLinSol_KLUGetSymbolic failure\n");
    fails += 1;
  } else {
    printf("    PASSED test -- SUNLinSol_KLUGetSymbolic \n");
  }
  numeric = SUNLinSol_KLUGetNumeric(LS);
  if (numeric->n != N) {
    printf("FAIL: SUNLinSol_KLUGetNumeric failure\n");
    fails += 1;
  } else {
    printf("    PASSED test -- SUNLinSol_KLUGetNumeric \n");
  }
  common = SUNLinSol_KLUGetCommon(LS);
  if (common->singular_col != N) {
    printf("FAIL: SUNLinSol_KLUGetCommon failure\n");
    fails += 1;
  } else {
    printf("    PASSED test -- SUNLinSol_KLUGetCommon \n");
  }

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nx (original) =\n");
    N_VPrint_Serial(y);
    printf("\nb =\n");
    N_VPrint_Serial(b);
    printf("\nx (computed) =\n");
    N_VPrint_Serial(x);
  } else {
    printf("SUCCESS: SUNLinSol module passed all tests \n \n");
  }

  /* Free solver, matrix and vectors */
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
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
  sunindextype i, local_length, maxloc;
  realtype *Xdata, *Ydata, maxerr;

  Xdata = N_VGetArrayPointer(X);
  Ydata = N_VGetArrayPointer(Y);
  local_length = N_VGetLength_Serial(X);

  /* check vector data */
  for(i=0; i < local_length; i++)
    failure += SUNRCompareTol(Xdata[i], Ydata[i], tol);

  if (failure > ZERO) {
    maxerr = ZERO;
    maxloc = -1;
    for(i=0; i < local_length; i++) {
      if (SUNRabs(Xdata[i]-Ydata[i]) >  maxerr) {
        maxerr = SUNRabs(Xdata[i]-Ydata[i]);
        maxloc = i;
      }
    }
    printf("check err failure: maxerr = %g at loc %li (tol = %g)\n",
	   maxerr, (long int) maxloc, tol);
    return(1);
  }
  else
    return(0);
}

void sync_device()
{
}
