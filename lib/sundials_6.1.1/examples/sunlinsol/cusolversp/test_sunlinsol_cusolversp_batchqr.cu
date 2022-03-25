/* -----------------------------------------------------------------
 * Programmer(s): Cody J.Balos @ LLNL
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
 * This is the testing routine to check the SUNLinSol cuSolverSp
 * module  implementation.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunmatrix/sunmatrix_cusparse.h>
#include <nvector/nvector_cuda.h>
#include <nvector/nvector_serial.h>
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
  SUNMatrix       A, B, dA;           /* test matrices                 */
  N_Vector        x, b, d_x, d_xref, d_b;/* test vectors                  */
  realtype        *matdata, *xdata, *xrefdata;
  int             print_timing;
  sunindextype    i, j;
  SUNContext      sunctx;

  cusparseStatus_t cusp_status;
  cusolverStatus_t cusol_status;
  cusparseHandle_t cusp_handle;
  cusolverSpHandle_t cusol_handle;


  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

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

  printf("\ncuSolverSp linear solver test: size %ld, block size %ld, number of blocks %ld\n\n",
    (long int) N, (long int) block_size, (long int) nblocks);

  /* Initialize cuSPARSE */
  cusp_status = cusparseCreate(&cusp_handle);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS) {
    printf("ERROR: could not create cuSPARSE handle\n");
    return(-1);
  }

  /* Initialize cuSOLVER */
  cusol_status = cusolverSpCreate(&cusol_handle);
  if (cusol_status != CUSOLVER_STATUS_SUCCESS) {
    printf("ERROR: could not create cuSOLVER handle\n");
    return(-1);
  }

  /* Create matrices and vectors */
  B = SUNDenseMatrix(N, N, sunctx);
  d_x = N_VNew_Cuda(N, sunctx);
  d_xref = N_VNew_Cuda(N, sunctx);
  d_b = N_VNew_Cuda(N, sunctx);
  x = N_VMake_Serial(N, N_VGetHostArrayPointer_Cuda(d_x), sunctx);
  b = N_VMake_Serial(N, N_VGetHostArrayPointer_Cuda(d_b), sunctx);

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

  /* Create the device matrix */
  dA = SUNMatrix_cuSparse_NewBlockCSR(nblocks, block_size, block_size, block_nnz, cusp_handle, sunctx);
  if (dA == NULL) {
    printf("ERROR: could not create dA\n");
  }

  /* Copy data to device */
  fails = SUNMatrix_cuSparse_CopyToDevice(dA, SUNSparseMatrix_Data(A),
                                          SUNSparseMatrix_IndexPointers(A),
                                          SUNSparseMatrix_IndexValues(A));
  if (fails != 0) {
    printf("ERROR: could not copy A to the device\n");
    return(-1);
  }

  /* Fill x vector with uniform random data in [0,1] */
  xdata = N_VGetHostArrayPointer_Cuda(d_x);
  xrefdata = N_VGetHostArrayPointer_Cuda(d_xref);
  for (i=0; i<N; i++) {
    realtype tmp = (realtype) rand() / (realtype) RAND_MAX;
    xdata[i]    = tmp;
    xrefdata[i] = tmp;
  }
  N_VCopyToDevice_Cuda(d_x);
  N_VCopyToDevice_Cuda(d_xref);

  /* Synchronize before peforming dense operation on CPU */
  cudaDeviceSynchronize();

  /* create right-hand side vector for linear solve */
  fails = SUNMatMatvec(A, x, b);
  if (fails) {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");
    return(1);
  }
  N_VCopyToDevice_Cuda(d_b);

  /* Create cuSolverSp linear solver
   * The BatchedQR method allows you to solve many small subsystems in parallel.
   */
  LS = SUNLinSol_cuSolverSp_batchQR(d_x, dA, cusol_handle, sunctx);

  if (LS == NULL) {
    printf("FAIL: SUNLinSol_cuSolverSp_batchQR returned NULL\n");
    return(1);
  }

  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSetup(LS, dA, 0);
  fails += Test_SUNLinSolSolve(LS, dA, d_x, d_b, 1000*UNIT_ROUNDOFF, SUNTRUE, 0);

  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
  fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_CUSOLVERSP_BATCHQR, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);

    SUNMatrix_cuSparse_CopyFromDevice(dA, SUNSparseMatrix_Data(A), NULL, NULL);
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);

    N_VCopyFromDevice_Cuda(d_xref);
    printf("x (reference)\n");
    N_VPrint_Cuda(d_xref);

    N_VCopyFromDevice_Cuda(d_x); /* copy solution from device */
    printf("x (computed)\n");
    N_VPrint_Cuda(d_x);

    N_VCopyFromDevice_Cuda(d_b);
    printf("\nb = Ax (reference)\n");
    N_VPrint_Cuda(d_b);
  } else {
    printf("SUCCESS: SUNLinSol module passed all tests \n \n");
  }

  /* Free solver, matrix and vectors */
  SUNLinSolFree(LS);
  SUNMatDestroy(A); SUNMatDestroy(dA);
  N_VDestroy(x); N_VDestroy(d_x); N_VDestroy(d_xref);
  N_VDestroy(b); N_VDestroy(d_b);

  /* Destroy the cuSOLVER and cuSPARSE handles */
  cusparseDestroy(cusp_handle);
  cusolverSpDestroy(cusol_handle);

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

  cudaDeviceSynchronize();

  N_VCopyFromDevice_Cuda(X);
  N_VCopyFromDevice_Cuda(Y);

  Xdata = N_VGetHostArrayPointer_Cuda(X);
  Ydata = N_VGetHostArrayPointer_Cuda(Y);
  local_length = N_VGetLength(X);

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
  cudaDeviceSynchronize();
}
