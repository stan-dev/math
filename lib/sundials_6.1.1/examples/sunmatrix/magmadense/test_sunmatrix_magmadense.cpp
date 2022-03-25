/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This is the testing routine to check the SUNMatrix Dense module
 * implementation.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_magmadense.h>
#include <sundials/sundials_math.h>
#include "test_sunmatrix.h"

#if defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#define HIP_OR_CUDA(a,b) a
#elif defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#define HIP_OR_CUDA(a,b) b
#else
#define HIP_OR_CUDA(a,b) ((void)0);
#endif

#if defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#include <nvector/nvector_cuda.h>
#include <sunmemory/sunmemory_cuda.h>
#elif defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#include <nvector/nvector_hip.h>
#include <sunmemory/sunmemory_hip.h>
#endif


/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails = 0;        /* counter for test failures  */
  sunindextype matrows, matcols; /* matrix dimensions          */
  sunindextype nblocks;          /* number of blocks in matrix */
  N_Vector     x, y;             /* test vectors               */
  realtype     *xdata, *ydata;   /* pointers to vector data    */
  SUNMatrix    A, I;             /* test matrices              */
  realtype     *Adata, *Idata;   /* pointers to matrix data    */
  int          print_timing, square;
  sunindextype i, j, k, m, n;
  SUNContext   sunctx;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  SUNMemoryHelper memhelper = HIP_OR_CUDA( SUNMemoryHelper_Hip(sunctx);,
                                           SUNMemoryHelper_Cuda(sunctx); )

  /* check input and set vector length */
  if (argc < 5) {
    printf("ERROR: FOUR (4) Input required: matrix rows, matrix cols, number of matrix blocks, print timing \n");
    return(-1);
  }

  matrows = (sunindextype) atol(argv[1]);
  if (matrows <= 0) {
    printf("ERROR: number of rows must be a positive integer \n");
    return(-1);
  }

  matcols = (sunindextype) atol(argv[2]);
  if (matcols <= 0) {
    printf("ERROR: number of cols must be a positive integer \n");
    return(-1);
  }

  nblocks = (sunindextype) atol(argv[3]);
  if (nblocks <= 0) {
    printf("ERROR: number of blocks must be a positive integer \n");
    return(-1);
  }

  print_timing = atoi(argv[4]);
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
  x = HIP_OR_CUDA( N_VNew_Hip(matcols*nblocks, sunctx);,
                   N_VNew_Cuda(matcols*nblocks, sunctx); )
  y = HIP_OR_CUDA( N_VNew_Hip(matrows*nblocks, sunctx);,
                   N_VNew_Cuda(matrows*nblocks, sunctx); )
  A = SUNMatrix_MagmaDenseBlock(nblocks, matrows, matcols, SUNMEMTYPE_DEVICE, memhelper, NULL, sunctx);
  if (square)
    I = SUNMatClone(A);

  /* Allocate host data */
  Adata = (realtype*) malloc(sizeof(realtype)*SUNMatrix_MagmaDense_LData(A));
  if (square)
    Idata = (realtype*) malloc(sizeof(realtype)*SUNMatrix_MagmaDense_LData(I));

  /* Fill matrices and vectors */
  for(k=0; k < nblocks; k++) {
    for(j=0; j < matcols; j++) {
      for(i=0; i < matrows; i++) {
        Adata[k*matcols*matrows + j*matrows + i] = (j+1)*(i+j);
      }
    }
  }
  SUNMatrix_MagmaDense_CopyToDevice(A, Adata);

  if (square) {
    for(k=0; k < nblocks; k++) {
      for(j=0; j < matcols; j++) {
        for(i=0; i < matrows; i++) {
          Idata[k*matcols*matrows + j*matrows + i] = (j == i) ? ONE : ZERO;
        }
      }
    }
    SUNMatrix_MagmaDense_CopyToDevice(I, Idata);
  }

  xdata = N_VGetArrayPointer(x);
  for(k=0; k < nblocks; k++) {
    for(i=0; i < matcols; i++) {
      xdata[matcols*k+i] = ONE / (i+1);
    }
  }
  HIP_OR_CUDA(  N_VCopyToDevice_Hip(x);,
                N_VCopyToDevice_Cuda(x); )

  ydata = N_VGetArrayPointer(y);
  for(k=0; k < nblocks; k++) {
    for(i=0; i < matrows; i++) {
      m = i;
      n = m + matcols - 1;
      ydata[matrows*k+i] = HALF*(n+1-m)*(n+m);
    }
  }
  HIP_OR_CUDA(  N_VCopyToDevice_Hip(y);,
                N_VCopyToDevice_Cuda(y); )

  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(A, SUNMATRIX_MAGMADENSE, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  if (square) {
    fails += Test_SUNMatScaleAdd(A, I, 0);
    fails += Test_SUNMatScaleAddI(A, I, 0);
  }
  fails += Test_SUNMatMatvecSetup(A, 0);
  fails += Test_SUNMatMatvec(A, x, y, 0);
  fails += Test_SUNMatSpace(A, 0);

  /* Print result */
  if (fails)
    printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
  else
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");

  /* Free vectors and matrices */
  N_VDestroy(x);
  N_VDestroy(y);
  SUNMatDestroy(A);
  free(Adata);
  if (square) {
    SUNMatDestroy(I);
    free(Idata);
  }
  SUNContext_Free(&sunctx);

  return(fails);
}

/* ----------------------------------------------------------------------
 * Check matrix
 * --------------------------------------------------------------------*/
int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  sunindextype i = 0;
  sunindextype Aldata = SUNMatrix_MagmaDense_LData(A);
  sunindextype Bldata = SUNMatrix_MagmaDense_LData(B);
  realtype *Adata = (realtype*) malloc(sizeof(realtype)*Aldata);
  realtype *Bdata = (realtype*) malloc(sizeof(realtype)*Bldata);

  /* copy data to host */
  SUNMatrix_MagmaDense_CopyFromDevice(A, Adata);
  SUNMatrix_MagmaDense_CopyFromDevice(B, Bdata);

  /* check lengths */
  if (Aldata != Bldata) {
    printf(">>> ERROR: check_matrix: Different data array lengths \n");
    return(1);
  }

  /* compare data */
  for(i=0; i < Aldata; i++) {
    failure += SUNRCompareTol(Adata[i], Bdata[i], tol);
  }

  free(Adata);
  free(Bdata);

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure = 0;
  sunindextype i = 0;
  sunindextype Aldata = SUNMatrix_MagmaDense_LData(A);
  realtype *Adata = (realtype*) malloc(sizeof(realtype)*Aldata);

  /* copy data to host */
  SUNMatrix_MagmaDense_CopyFromDevice(A, Adata);

  /* compare data */
  for(i=0; i < Aldata; i++) {
    int check = SUNRCompareTol(Adata[i], val, tol);
    if (check) {
      printf("failed at %d\n", i);
      failure += check;
    }
  }

  free(Adata);

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_vector(N_Vector actual, N_Vector expected, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* copy vectors to host */
  HIP_OR_CUDA( N_VCopyFromDevice_Hip(actual);,
               N_VCopyFromDevice_Cuda(actual); )
  HIP_OR_CUDA( N_VCopyFromDevice_Hip(expected);,
               N_VCopyFromDevice_Cuda(expected); )

  /* get vector data */
  xdata = N_VGetArrayPointer(actual);
  ydata = N_VGetArrayPointer(expected);

  /* check data lengths */
  xldata = N_VGetLength(actual);
  yldata = N_VGetLength(expected);


  if (xldata != yldata) {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return(1);
  }

  /* check vector data */
  for(i=0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO) {
    printf("Check_vector failures:\n");
    for(i=0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        printf("  actual[%ld] = %g != %e (err = %g)\n", (long int) i,
               xdata[i], ydata[i], SUNRabs(xdata[i]-ydata[i]));
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(SUNMatrix A)
{
  realtype *Adata = SUNMatrix_MagmaDense_Data(A);
  if (Adata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  if (SUNMatrix_MagmaDense_Rows(A) == SUNMatrix_MagmaDense_Columns(A))
    return SUNTRUE;
  else
    return SUNFALSE;
}

void sync_device(SUNMatrix A)
{
  HIP_OR_CUDA( hipDeviceSynchronize();,
               cudaDeviceSynchronize(); )
}
