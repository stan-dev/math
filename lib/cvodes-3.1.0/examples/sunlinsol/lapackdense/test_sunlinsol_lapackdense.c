/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
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
 * SUNDenseLinearSolver Testing Routine
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

  /* check input and set matrix dimensions */
  if (argc < 3){
    printf("ERROR: TWO (2) Inputs required: matrix cols, print timing \n");
    return(-1);
  }

  cols = atol(argv[1]); 
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
  A = SUNDenseMatrix(rows, cols);
  B = SUNDenseMatrix(rows, cols);
  I = SUNDenseMatrix(rows, cols);
  x = N_VNew_Serial(cols);
  y = N_VNew_Serial(cols);
  b = N_VNew_Serial(cols);

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
  LS = SUNLapackDense(x, A);
  
  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSetup(LS, A, 0);
  fails += Test_SUNLinSolSolve(LS, A, x, b, RCONST(1.0e-15), 0);
 
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
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
  N_VDestroy(x);
  N_VDestroy(y);

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
    failure += FNEQ(Xdata[i], Ydata[i], tol);

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
