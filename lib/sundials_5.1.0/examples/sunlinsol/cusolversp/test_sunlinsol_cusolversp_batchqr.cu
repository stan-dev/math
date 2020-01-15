/* -----------------------------------------------------------------
 * Programmer(s): Cody J.Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol cuSolverSp
 * module  implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <nvector/nvector_cuda.h>
#include <sundials/sundials_math.h>
#include "test_sunlinsol.h"


/* ----------------------------------------------------------------------
 * SUNLinSol_KLU Linear Solver Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int             fails = 0;          /* counter for test failures     */
  sunindextype    N;                  /* matrix columns, rows          */
  int             block_size;         /* matrix block columns, rows    */
  int             block_nnz;          /* number of nonzeros in a block */
  int             block_nnz_max;      /* max nonzeros per block        */
  int             nblocks;            /* number of blocks              */
  SUNLinearSolver LS;                 /* linear solver object          */
  SUNMatrix       A, B;               /* test matrices                 */
  N_Vector        x, y, b;            /* test vectors                  */
  realtype        *matdata, *xdata;
  int             print_timing;
  sunindextype    i, j;

  /* check input and set matrix dimensions */
  if (argc < 4){
    printf("ERROR: THREE (3) Inputs required: matrix block size, number of blocks, print timing \n");
    return(-1);
  }

  block_size = atol(argv[1]);
  if (block_size <= 0) {
    printf("ERROR: matrix size must be a positive integer \n");
    return(-1);
  }

  block_nnz_max = block_size*block_size / 4;

  nblocks = atol(argv[2]);
  if (nblocks <= 0) {
    printf("ERROR: number of blocks must be a positive integer \n");
    return(-1);
  }

  /* calculate the size of the overall martrix */
  N = block_size * nblocks;

  print_timing = atoi(argv[3]);
  SetTiming(print_timing);

  printf("\ncuSolverSp linear solver test: size %ld, block size %ld, number of blocks %d\n\n",
    (long int) N, (long int) block_size, (long int) nblocks);

  /* Create matrices and vectors */
  B = SUNDenseMatrix(N, N);
  x = N_VNewManaged_Cuda(N);
  y = N_VNewManaged_Cuda(N);
  b = N_VNewManaged_Cuda(N);

  /* Zero the matrix */
  fails = SUNMatZero(B);

  /* Create sparsity pattern for a block. */
  sunindextype *cols = (sunindextype *) malloc(block_nnz_max*sizeof(sunindextype));
  sunindextype *rows = (sunindextype *) malloc(block_nnz_max*sizeof(sunindextype));
  for (i=0; i<block_nnz_max; i++) {
    cols[i] = rand() % block_size;
    rows[i] = rand() % block_size;
  }

  /* Fill matrix with uniform random data in [0,1/N] */
  for (i=0; i<nblocks; i++) {
    for (j=0; j<block_nnz_max; j++) {
      sunindextype col = cols[j] + block_size*i;
      sunindextype row = rows[j] + block_size*i;
      matdata = SUNDenseMatrix_Column(B,col);
      matdata[row] = (realtype) rand() / (realtype) RAND_MAX / N;
    }
  }

  /* Free temporary rows and cols variables */
  free(cols); free(rows);

  /* Add identity to B */
  fails = SUNMatScaleAddI(ONE, B);
  if (fails) {
    printf("FAIL: SUNMatScaleAddI failure\n");
    return(1);
  }

  /* Create sparse matrix from dense, and destroy B */
  A = SUNSparseFromDenseMatrix(B, ZERO, CSR_MAT);
  SUNMatDestroy(B);

  /* Calculate actual number of nonzeros per block */
  block_nnz = SUNSparseMatrix_NNZ(A) / nblocks;

  /* Fill x vector with uniform random data in [0,1] */
  xdata = N_VGetHostArrayPointer_Cuda(x);
  for (i=0; i<N; i++)
    xdata[i] = (realtype) rand() / (realtype) RAND_MAX;
  N_VCopyToDevice_Cuda(x);

  /* copy x into y to print in case of solver failure */
  N_VScale(ONE, x, y);

  /* create right-hand side vector for linear solve */
  fails = SUNMatMatvec(A, x, b);
  if (fails) {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");
    return(1);
  }

  /* Create cuSolverSp linear solver
   * The BatchedQR method allows you to solve many small subsystems in parallel.
   */
  LS = SUNLinSol_cuSolverSp_batchQR(x,           /* the overall system vector */
                                    A,           /* the overall system matrix */
                                    nblocks,     /* number of subsystems */
                                    block_size,  /* size of a subsystem  */
                                    block_nnz);  /* number of nonzeros in a subsystem */

  if (LS == NULL) {
    printf("FAIL: SUNLinSol_cuSolverSp_batchQR returned NULL\n");
    return(1);
  }

  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSetup(LS, A, 0);
  fails += Test_SUNLinSolSolve(LS, A, x, b, 1000*UNIT_ROUNDOFF, 0);

  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
  fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_CUSOLVERSP_BATCHQR, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nx (original) =\n");
    N_VPrint_Cuda(y);
    printf("\nb =\n");
    N_VPrint_Cuda(b);
    printf("\nx (computed) =\n");
    N_VPrint_Cuda(x);
  } else {
    printf("SUCCESS: SUNLinSol module passed all tests \n \n");
  }

  /* Free solver, matrix and vectors */
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(b);

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

  cudaDeviceSynchronize();

  Xdata = N_VGetHostArrayPointer_Cuda(X);
  Ydata = N_VGetHostArrayPointer_Cuda(Y);
  local_length = N_VGetLength(X);

  /* check vector data */
  for(i=0; i < local_length; i++)
    failure += FNEQ(Xdata[i], Ydata[i], tol);

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
