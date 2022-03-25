/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
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
 * This is the testing routine to check the SUNMatrix Sparse module
 * implementation.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_sunmatrix.h"


/* prototypes for custom tests */
int Test_SUNMatScaleAdd2(SUNMatrix A, SUNMatrix B, N_Vector x,
                         N_Vector y, N_Vector z);
int Test_SUNMatScaleAddI2(SUNMatrix A, N_Vector x, N_Vector y);
int Test_SUNSparseMatrixToCSC(SUNMatrix A);
int Test_SUNSparseMatrixToCSR(SUNMatrix A);

/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails=0;                    /* counter for test failures  */
  sunindextype matrows, matcols;           /* matrix dims                */
  int          mattype;                    /* matrix storage type        */
  N_Vector     x, y, z;                    /* test vectors               */
  realtype*    vecdata;                    /* pointers to vector data    */
  SUNMatrix    A, B, C, D, I;              /* test matrices              */
  realtype*    matdata;                    /* pointer to matrix data     */
  sunindextype i, j, k, kstart, kend, N, uband, lband;
  sunindextype *colptrs, *rowindices;
  sunindextype *rowptrs, *colindices;
  int          print_timing, square;
  SUNContext   sunctx;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  /* check input and set vector length */
  if (argc < 5){
    printf("ERROR: FOUR (4) Input required: matrix rows, matrix cols, matrix type (0/1), print timing \n");
    return(-1);
  }

  matrows = (sunindextype) atol(argv[1]);
  if (matrows < 1) {
    printf("ERROR: number of rows must be a positive integer\n");
    return(-1);
  }

  matcols = (sunindextype) atol(argv[2]);
  if (matcols < 1) {
    printf("ERROR: number of cols must be a positive integer\n");
    return(-1);
  }

  k = (sunindextype) atol(argv[3]);
  if ((k != 0) && (k != 1)) {
    printf("ERROR: matrix type must be 0 or 1\n");
    return(-1);
  }
  mattype = (k == 0) ? CSC_MAT : CSR_MAT;

  print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  square = (matrows == matcols) ? 1 : 0;
  printf("\nSparse matrix test: size %ld by %ld, type = %i\n\n",
         (long int) matrows, (long int) matcols, mattype);

  /* Initialize vectors and matrices to NULL */
  x = NULL;
  y = NULL;
  z = NULL;
  A = NULL;
  B = NULL;
  C = NULL;
  D = NULL;
  I = NULL;

  /* check creating sparse matrix from dense matrix */
  B = SUNDenseMatrix(5,6, sunctx);

  matdata = SUNDenseMatrix_Data(B);
  matdata[2]  = RCONST(1.0);    /* [ 0 2 0 0 7 0 ] */
  matdata[5]  = RCONST(2.0);    /* [ 0 0 4 0 8 0 ] */
  matdata[9]  = RCONST(3.0);    /* [ 1 0 0 0 0 0 ] */
  matdata[11] = RCONST(4.0);    /* [ 0 0 5 6 0 0 ] */
  matdata[13] = RCONST(5.0);    /* [ 0 3 0 0 0 9 ] */
  matdata[18] = RCONST(6.0);
  matdata[20] = RCONST(7.0);
  matdata[21] = RCONST(8.0);
  matdata[29] = RCONST(9.0);

  if (mattype == CSR_MAT) {

    /* Check CSR */
    C = SUNSparseMatrix(5, 6, 9, CSR_MAT, sunctx);
    rowptrs = SUNSparseMatrix_IndexPointers(C);
    colindices = SUNSparseMatrix_IndexValues(C);
    matdata = SUNSparseMatrix_Data(C);
    rowptrs[0] = 0;
    matdata[0] = RCONST(2.0);   colindices[0] = 1;
    matdata[1] = RCONST(7.0);   colindices[1] = 4;
    rowptrs[1] = 2;
    matdata[2] = RCONST(4.0);   colindices[2] = 2;
    matdata[3] = RCONST(8.0);   colindices[3] = 4;
    rowptrs[2] = 4;
    matdata[4] = RCONST(1.0);   colindices[4] = 0;
    rowptrs[3] = 5;
    matdata[5] = RCONST(5.0);   colindices[5] = 2;
    matdata[6] = RCONST(6.0);   colindices[6] = 3;
    rowptrs[4] = 7;
    matdata[7] = RCONST(3.0);   colindices[7] = 1;
    matdata[8] = RCONST(9.0);   colindices[8] = 5;
    rowptrs[5] = 9;

    A = SUNSparseFromDenseMatrix(B, ZERO, CSR_MAT);
    fails += check_matrix(A, C, 1e-15);

    if (fails) {
      printf("FAIL: SUNMatrix SparseFromDense CSR conversion failed\n");
      return(1);
    }

    SUNMatDestroy(A);
    SUNMatDestroy(C);

  } else {

    /* Check CSC */
    D = SUNSparseMatrix(5, 6, 9, CSC_MAT, sunctx);
    colptrs = SUNSparseMatrix_IndexPointers(D);
    rowindices = SUNSparseMatrix_IndexValues(D);
    matdata = SUNSparseMatrix_Data(D);
    colptrs[0] = 0;
    matdata[0] = RCONST(1.0);   rowindices[0] = 2;
    colptrs[1] = 1;
    matdata[1] = RCONST(2.0);   rowindices[1] = 0;
    matdata[2] = RCONST(3.0);   rowindices[2] = 4;
    colptrs[2] = 3;
    matdata[3] = RCONST(4.0);   rowindices[3] = 1;
    matdata[4] = RCONST(5.0);   rowindices[4] = 3;
    colptrs[3] = 5;
    matdata[5] = RCONST(6.0);   rowindices[5] = 3;
    colptrs[4] = 6;
    matdata[6] = RCONST(7.0);   rowindices[6] = 0;
    matdata[7] = RCONST(8.0);   rowindices[7] = 1;
    colptrs[5] = 8;
    matdata[8] = RCONST(9.0);   rowindices[8] = 4;
    colptrs[6] = 9;

    A = SUNSparseFromDenseMatrix(B, 1e-15, CSC_MAT);
    fails += check_matrix(A, D, 1e-15);

    if (fails) {
      printf("FAIL: SUNMatrix SparseFromDense CSC conversion failed\n");
      return(1);
    }

    SUNMatDestroy(A);
    SUNMatDestroy(D);

  }
  SUNMatDestroy(B);


  /* check creating sparse matrix from banded matrix */
  N = 7;
  uband = 1;
  lband = 2;                                   /* B(i,j) = j + (j-i) */
  B = SUNBandMatrix(N, uband, lband, sunctx);  /* B = [  0  2  0  0  0  0  0 ] */
  for (j=0; j<N; j++) {                        /*     [ -1  1  3  0  0  0  0 ] */
    matdata = SUNBandMatrix_Column(B, j);      /*     [ -2  0  2  4  0  0  0 ] */
    kstart = (j<uband) ? -j : -uband;          /*     [  0 -1  1  3  5  0  0 ] */
    kend = (j>N-1-lband) ? N-1-j: lband;       /*     [  0  0  0  2  4  6  0 ] */
    for (k=kstart; k<=kend; k++)               /*     [  0  0  0  1  3  5  7 ] */
      matdata[k] = j - k;                      /*     [  0  0  0  0  2  4  6 ] */
  }

  if (mattype == CSR_MAT) {

    /* CSR */
    C = SUNSparseMatrix(7, 7, 21, CSR_MAT, sunctx);
    rowptrs = SUNSparseMatrix_IndexPointers(C);
    colindices = SUNSparseMatrix_IndexValues(C);
    matdata = SUNSparseMatrix_Data(C);
    rowptrs[ 0] = 0;
    matdata[ 0] = RCONST(2.0);   colindices[ 0] = 1;
    rowptrs[ 1] = 1;
    matdata[ 1] = RCONST(-1.0);  colindices[ 1] = 0;
    matdata[ 2] = RCONST(1.0);   colindices[ 2] = 1;
    matdata[ 3] = RCONST(3.0);   colindices[ 3] = 2;
    rowptrs[ 2] = 4;
    matdata[ 4] = RCONST(-2.0);  colindices[ 4] = 0;
    matdata[ 5] = RCONST(2.0);   colindices[ 5] = 2;
    matdata[ 6] = RCONST(4.0);   colindices[ 6] = 3;
    rowptrs[ 3] = 7;
    matdata[ 7] = RCONST(-1.0);  colindices[ 7] = 1;
    matdata[ 8] = RCONST(1.0);   colindices[ 8] = 2;
    matdata[ 9] = RCONST(3.0);   colindices[ 9] = 3;
    matdata[10] = RCONST(5.0);   colindices[10] = 4;
    rowptrs[ 4] = 11;
    matdata[11] = RCONST(2.0);   colindices[11] = 3;
    matdata[12] = RCONST(4.0);   colindices[12] = 4;
    matdata[13] = RCONST(6.0);   colindices[13] = 5;
    rowptrs[ 5] = 14;
    matdata[14] = RCONST(1.0);   colindices[14] = 3;
    matdata[15] = RCONST(3.0);   colindices[15] = 4;
    matdata[16] = RCONST(5.0);   colindices[16] = 5;
    matdata[17] = RCONST(7.0);   colindices[17] = 6;
    rowptrs[ 6] = 18;
    matdata[18] = RCONST(2.0);   colindices[18] = 4;
    matdata[19] = RCONST(4.0);   colindices[19] = 5;
    matdata[20] = RCONST(6.0);   colindices[20] = 6;
    rowptrs[ 7] = 21;

    A = SUNSparseFromBandMatrix(B, ZERO, CSR_MAT);
    fails += check_matrix(A, C, 1e-15);

    if (fails) {
      printf("FAIL: SUNMatrix SparseFromBand CSR conversion failed\n");
      return(1);
    }

    SUNMatDestroy(A);
    SUNMatDestroy(C);

  } else {

    /* Check CSC */
    D = SUNSparseMatrix(7, 7, 21, CSC_MAT, sunctx);
    colptrs = SUNSparseMatrix_IndexPointers(D);
    rowindices = SUNSparseMatrix_IndexValues(D);
    matdata = SUNSparseMatrix_Data(D);
    colptrs[ 0] = 0;
    matdata[ 0] = RCONST(-1.0);  rowindices[ 0] = 1;
    matdata[ 1] = RCONST(-2.0);  rowindices[ 1] = 2;
    colptrs[ 1] = 2;
    matdata[ 2] = RCONST(2.0);   rowindices[ 2] = 0;
    matdata[ 3] = RCONST(1.0);   rowindices[ 3] = 1;
    matdata[ 4] = RCONST(-1.0);  rowindices[ 4] = 3;
    colptrs[ 2] = 5;
    matdata[ 5] = RCONST(3.0);   rowindices[ 5] = 1;
    matdata[ 6] = RCONST(2.0);   rowindices[ 6] = 2;
    matdata[ 7] = RCONST(1.0);   rowindices[ 7] = 3;
    colptrs[ 3] = 8;
    matdata[ 8] = RCONST(4.0);   rowindices[ 8] = 2;
    matdata[ 9] = RCONST(3.0);   rowindices[ 9] = 3;
    matdata[10] = RCONST(2.0);   rowindices[10] = 4;
    matdata[11] = RCONST(1.0);   rowindices[11] = 5;
    colptrs[ 4] = 12;
    matdata[12] = RCONST(5.0);   rowindices[12] = 3;
    matdata[13] = RCONST(4.0);   rowindices[13] = 4;
    matdata[14] = RCONST(3.0);   rowindices[14] = 5;
    matdata[15] = RCONST(2.0);   rowindices[15] = 6;
    colptrs[ 5] = 16;
    matdata[16] = RCONST(6.0);   rowindices[16] = 4;
    matdata[17] = RCONST(5.0);   rowindices[17] = 5;
    matdata[18] = RCONST(4.0);   rowindices[18] = 6;
    colptrs[ 6] = 19;
    matdata[19] = RCONST(7.0);   rowindices[19] = 5;
    matdata[20] = RCONST(6.0);   rowindices[20] = 6;
    colptrs[ 7] = 21;

    A = SUNSparseFromBandMatrix(B, 1e-15, CSC_MAT);
    fails += check_matrix(A, D, 1e-15);

    if (fails) {
      printf("FAIL: SUNMatrix SparseFromBand CSC conversion failed\n");
      return(1);
    }

    SUNMatDestroy(A);
    SUNMatDestroy(D);
  }

  SUNMatDestroy(B);


  /* Create/fill I matrix */
  I = NULL;
  if (square) {
    I = SUNSparseMatrix(matrows, matcols, matcols, mattype, sunctx);
    matdata    = SUNSparseMatrix_Data(I);
    colindices = SUNSparseMatrix_IndexValues(I);
    rowptrs    = SUNSparseMatrix_IndexPointers(I);
    for(i=0; i<matrows; i++) {
      matdata[i] = ONE;
      colindices[i] = i;
      rowptrs[i] = i;
    }
    rowptrs[matrows] = matrows;
  }

  /* Create/fill random dense matrices, create sparse from them */
  C = SUNDenseMatrix(matrows, matcols, sunctx);
  D = SUNDenseMatrix(matrows, matcols, sunctx);
  for (k=0; k<3*matrows; k++) {
    i = rand() % matrows;
    j = rand() % matcols;
    matdata = SUNDenseMatrix_Column(D,j);
    matdata[i] = (realtype) rand() / (realtype) RAND_MAX;
  }
  for (k=0; k<matrows; k++) {
    i = rand() % matrows;
    j = rand() % matcols;
    matdata = SUNDenseMatrix_Column(C,j);
    matdata[i] = (realtype) rand() / (realtype) RAND_MAX;
  }
  A = SUNSparseFromDenseMatrix(C, ZERO, mattype);
  B = SUNSparseFromDenseMatrix(D, ZERO, mattype);

  /* Create vectors and fill */
  x = N_VNew_Serial(matcols, sunctx);
  y = N_VNew_Serial(matrows, sunctx);
  z = N_VNew_Serial(matrows, sunctx);
  vecdata = N_VGetArrayPointer(x);
  for(i=0; i<matcols; i++)
    vecdata[i] = (realtype) rand() / (realtype) RAND_MAX;
  if (SUNMatMatvec(C, x, y) != 0) {
    printf("FAIL: SUNMatrix module Dense matvec failure \n \n");
    SUNMatDestroy(A);  SUNMatDestroy(B);
    SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(x);  N_VDestroy(y);  N_VDestroy(z);
    if (square)
      SUNMatDestroy(I);
    return(1);
  }
  if (SUNMatMatvec(D, x, z) != 0) {
    printf("FAIL: SUNMatrix module Dense matvec failure \n \n");
    SUNMatDestroy(A);  SUNMatDestroy(B);
    SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(x);  N_VDestroy(y);  N_VDestroy(z);
    if (square)
      SUNMatDestroy(I);
    return(1);
  }

  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(A, SUNMATRIX_SPARSE, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  fails += Test_SUNMatScaleAdd(A, I, 0);
  fails += Test_SUNMatScaleAdd2(A, B, x, y, z);
  if (square) {
    fails += Test_SUNMatScaleAddI(A, I, 0);
    fails += Test_SUNMatScaleAddI2(A, x, y);
  }
  fails += Test_SUNMatMatvec(A, x, y, 0);
  fails += Test_SUNMatSpace(A, 0);
  if (mattype == CSR_MAT) {
    fails += Test_SUNSparseMatrixToCSC(A);
  } else {
    fails += Test_SUNSparseMatrixToCSR(A);
  }

  /* Print result */
  if (fails) {
    printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nB =\n");
    SUNSparseMatrix_Print(B,stdout);
    if (square) {
      printf("\nI =\n");
      SUNSparseMatrix_Print(I,stdout);
    }
    printf("\nx =\n");
    N_VPrint_Serial(x);
    printf("\ny =\n");
    N_VPrint_Serial(y);
    printf("\nz =\n");
    N_VPrint_Serial(z);
  } else {
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");
  }

  /* Free vectors and matrices */
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(z);
  SUNMatDestroy(A);
  SUNMatDestroy(B);
  SUNMatDestroy(C);
  SUNMatDestroy(D);
  if (square)
    SUNMatDestroy(I);

  SUNContext_Free(&sunctx);

  return(fails);
}

/* ----------------------------------------------------------------------
 * Extra ScaleAdd tests for sparse matrices:
 *    A and B should have different sparsity patterns, and neither should
 *      contain sufficient storage to for their sum
 *    y should already equal A*x
 *    z should already equal B*x
 * --------------------------------------------------------------------*/
int Test_SUNMatScaleAdd2(SUNMatrix A, SUNMatrix B, N_Vector x,
                         N_Vector y, N_Vector z)
{
  int       failure;
  SUNMatrix C, D, E;
  N_Vector  u, v;
  realtype  tol=100*UNIT_ROUNDOFF;

  /* create clones for test */
  C = SUNMatClone(A);
  u = N_VClone(y);
  v = N_VClone(y);

  /* test 1: add A to B (output must be enlarged) */
  failure = SUNMatCopy(A, C);            /* C = A */
  if (failure) {
    printf(">>> FAILED test -- SUNMatCopy returned %d \n",
           failure);
    SUNMatDestroy(C);  N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  failure = SUNMatScaleAdd(ONE, C, B);   /* C = A+B */
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAdd returned %d \n",
           failure);
    SUNMatDestroy(C);  N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  failure = SUNMatMatvec(C, x, u);       /* u = Cx = Ax+Bx */
  if (failure) {
    printf(">>> FAILED test -- SUNMatMatvec returned %d \n",
           failure);
    SUNMatDestroy(C);  N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  N_VLinearSum(ONE,y,ONE,z,v);           /* v = y+z */
  failure = check_vector(u, v, tol);     /* u ?= v */
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAdd2 check 1 \n");
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nB =\n");
    SUNSparseMatrix_Print(B,stdout);
    printf("\nC =\n");
    SUNSparseMatrix_Print(C,stdout);
    printf("\nx =\n");
    N_VPrint_Serial(x);
    printf("\ny =\n");
    N_VPrint_Serial(y);
    printf("\nz =\n");
    N_VPrint_Serial(z);
    printf("\nu =\n");
    N_VPrint_Serial(u);
    printf("\nv =\n");
    N_VPrint_Serial(v);
    SUNMatDestroy(C);  N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  else {
    printf("    PASSED test -- SUNMatScaleAdd2 check 1 \n");
  }

  /* test 2: add A to a matrix with sufficient but misplaced storage */
  D = SUNMatClone(A);
  failure = SUNSparseMatrix_Reallocate(D, SM_NNZ_S(A)+SM_NNZ_S(B));
  failure = SUNMatCopy(A, D);            /* D = A */
  if (failure) {
    printf(">>> FAILED test -- SUNMatCopy returned %d \n",
           failure);
    SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  failure = SUNMatScaleAdd(ONE, D, B);   /* D = A+B */
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAdd returned %d \n",
           failure);
    SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  failure = SUNMatMatvec(D, x, u);       /* u = Cx = Ax+Bx */
  if (failure) {
    printf(">>> FAILED test -- SUNMatMatvec returned %d \n",
           failure);
    SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  N_VLinearSum(ONE,y,ONE,z,v);           /* v = y+z */
  failure = check_vector(u, v, tol);     /* u ?= v */
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAdd2 check 2 \n");
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nB =\n");
    SUNSparseMatrix_Print(B,stdout);
    printf("\nD =\n");
    SUNSparseMatrix_Print(D,stdout);
    printf("\nx =\n");
    N_VPrint_Serial(x);
    printf("\ny =\n");
    N_VPrint_Serial(y);
    printf("\nz =\n");
    N_VPrint_Serial(z);
    printf("\nu =\n");
    N_VPrint_Serial(u);
    printf("\nv =\n");
    N_VPrint_Serial(v);
    SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  else {
    printf("    PASSED test -- SUNMatScaleAdd2 check 2 \n");
  }


  /* test 3: add A to a matrix with the appropriate structure already in place */
  E = SUNMatClone(C);
  failure = SUNMatCopy(C, E);                /* E = A + B */
  if (failure) {
    printf(">>> FAILED test -- SUNMatCopy returned %d \n",
           failure);
    SUNMatDestroy(C);  SUNMatDestroy(D);  SUNMatDestroy(E);
    N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  failure = SUNMatScaleAdd(NEG_ONE, E, B);   /* E = -A */
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAdd returned %d \n",
           failure);
    SUNMatDestroy(C);  SUNMatDestroy(D);  SUNMatDestroy(E);
    N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  failure = SUNMatMatvec(E, x, u);           /* u = Ex = -Ax */
  if (failure) {
    printf(">>> FAILED test -- SUNMatMatvec returned %d \n",
           failure);
    SUNMatDestroy(C);  SUNMatDestroy(D);  SUNMatDestroy(E);
    N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  N_VLinearSum(NEG_ONE,y,ZERO,z,v);          /* v = -y */
  failure = check_vector(u, v, tol);         /* v ?= u */
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAdd2 check 3 \n");
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nB =\n");
    SUNSparseMatrix_Print(B,stdout);
    printf("\nC =\n");
    SUNSparseMatrix_Print(C,stdout);
    printf("\nE =\n");
    SUNSparseMatrix_Print(E,stdout);
    printf("\nx =\n");
    N_VPrint_Serial(x);
    printf("\ny =\n");
    N_VPrint_Serial(y);
    printf("\nu =\n");
    N_VPrint_Serial(u);
    printf("\nv =\n");
    N_VPrint_Serial(v);
    SUNMatDestroy(C);  SUNMatDestroy(D);  SUNMatDestroy(E);
    N_VDestroy(u);  N_VDestroy(v);  return(1);
  }
  else {
    printf("    PASSED test -- SUNMatScaleAdd2 check 3 \n");
  }

  SUNMatDestroy(C);
  SUNMatDestroy(D);
  SUNMatDestroy(E);
  N_VDestroy(u);
  N_VDestroy(v);
  return(0);
}

/* ----------------------------------------------------------------------
 * Extra ScaleAddI tests for sparse matrices:
 *    A should not contain values on the diagonal, nor should it contain
 *      sufficient storage to add those in
 *    y should already equal A*x
 * --------------------------------------------------------------------*/
int Test_SUNMatScaleAddI2(SUNMatrix A, N_Vector x, N_Vector y)
{
  int       failure;
  SUNMatrix B, C, D;
  N_Vector  w, z;
  realtype  tol=200*UNIT_ROUNDOFF;

  /* create clones for test */
  B = SUNMatClone(A);
  z = N_VClone(x);
  w = N_VClone(x);

  /* test 1: add I to a matrix with insufficient storage */
  failure = SUNMatCopy(A, B);
  if (failure) {
    printf(">>> FAILED test -- SUNMatCopy returned %d \n",
           failure);
    SUNMatDestroy(B);  N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  failure = SUNMatScaleAddI(NEG_ONE, B);   /* B = I-A */
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAddI returned %d \n",
           failure);
    SUNMatDestroy(B);  N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  failure = SUNMatMatvec(B, x, z);
  if (failure) {
    printf(">>> FAILED test -- SUNMatMatvec returned %d \n",
           failure);
    SUNMatDestroy(B);  N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  N_VLinearSum(ONE,x,NEG_ONE,y,w);  /* B x = (I - A) x = x - Ax = x - y = w */
  failure = check_vector(z, w, tol);
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAddI2 check 1 \n");
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nB =\n");
    SUNSparseMatrix_Print(B,stdout);
    printf("\nz =\n");
    N_VPrint_Serial(z);
    printf("\nw =\n");
    N_VPrint_Serial(w);
    SUNMatDestroy(B);  N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  else {
    printf("    PASSED test -- SUNMatScaleAddI2 check 1 \n");
  }

  /* test 2: add I to a matrix with sufficient but misplaced
     storage */
  C = SUNMatClone(A);
  failure = SUNSparseMatrix_Reallocate(C, SM_NNZ_S(A)+SM_ROWS_S(A));
  failure = SUNMatCopy(A, C);
  if (failure) {
    printf(">>> FAILED test -- SUNMatCopy returned %d \n",
           failure);
    SUNMatDestroy(B);  SUNMatDestroy(C);
    N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  failure = SUNMatScaleAddI(NEG_ONE, C);   /* C = I-A */
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAddI returned %d \n",
           failure);
    SUNMatDestroy(B);  SUNMatDestroy(C);
    N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  failure = SUNMatMatvec(C, x, z);
  if (failure) {
    printf(">>> FAILED test -- SUNMatMatvec returned %d \n",
           failure);
    SUNMatDestroy(B);  SUNMatDestroy(C);
    N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  N_VLinearSum(ONE,x,NEG_ONE,y,w);
  failure = check_vector(z, w, tol);
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAddI2 check 2 \n");
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nC =\n");
    SUNSparseMatrix_Print(C,stdout);
    printf("\nz =\n");
    N_VPrint_Serial(z);
    printf("\nw =\n");
    N_VPrint_Serial(w);
    SUNMatDestroy(B);  SUNMatDestroy(C);
    N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  else {
    printf("    PASSED test -- SUNMatScaleAddI2 check 2 \n");
  }


  /* test 3: add I to a matrix with appropriate structure already in place */
  D = SUNMatClone(C);
  failure = SUNMatCopy(C, D);
  if (failure) {
    printf(">>> FAILED test -- SUNMatCopy returned %d \n",
           failure);
    SUNMatDestroy(B);  SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  failure = SUNMatScaleAddI(NEG_ONE, D);   /* D = A */
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAddI returned %d \n",
           failure);
    SUNMatDestroy(B);  SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  failure = SUNMatMatvec(D, x, z);
  if (failure) {
    printf(">>> FAILED test -- SUNMatMatvec returned %d \n",
           failure);
    SUNMatDestroy(B);  SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  failure = check_vector(z, y, tol);
  if (failure) {
    printf(">>> FAILED test -- SUNMatScaleAddI2 check 3 \n");
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nD =\n");
    SUNSparseMatrix_Print(D,stdout);
    printf("\nz =\n");
    N_VPrint_Serial(z);
    printf("\ny =\n");
    N_VPrint_Serial(y);
    SUNMatDestroy(B);  SUNMatDestroy(C);  SUNMatDestroy(D);
    N_VDestroy(z);  N_VDestroy(w);  return(1);
  }
  else {
    printf("    PASSED test -- SUNMatScaleAddI2 check 3 \n");
  }

  SUNMatDestroy(B);
  SUNMatDestroy(C);
  SUNMatDestroy(D);
  N_VDestroy(z);
  N_VDestroy(w);
  return(0);
}

int Test_SUNSparseMatrixToCSR(SUNMatrix A)
{
  int       failure;
  SUNMatrix csc, csr;
  realtype  tol=200*UNIT_ROUNDOFF;

  failure = SUNSparseMatrix_ToCSR(A, &csr);

  if (failure) {
    printf(">>> FAILED test -- SUNSparseMatrix_ToCSR returned nonzero\n");
    SUNMatDestroy(csr);
    return(1);
  }

  /* check entries */
  if (SUNSparseMatrix_ToCSC(csr, &csc)) {
    printf(">>> FAILED test -- SUNSparseMatrix_ToCSC returned nonzero\n");
    SUNMatDestroy(csr);
    return(1);
  }

  if (check_matrix(A, csc, tol)) {
    printf(">>> FAILED test --  Test_SUNSparseMatrixToCSR check_matrix failed\n");
    SUNMatDestroy(csr);
    SUNMatDestroy(csc);
    return(1);
  }

  printf("    PASSED test -- SUNSparseMatrixToCSR\n");

  SUNMatDestroy(csr);
  SUNMatDestroy(csc);

  return(0);
}

int Test_SUNSparseMatrixToCSC(SUNMatrix A)
{
  int       failure;
  SUNMatrix csc=NULL, csr=NULL;
  realtype  tol=200*UNIT_ROUNDOFF;

  failure = SUNSparseMatrix_ToCSC(A, &csc);

  if (failure) {
    printf(">>> FAILED test -- SUNSparseMatrix_ToCSC returned nonzero\n");
    SUNMatDestroy(csr);
    return(1);
  }

  /* check entries */
  if (SUNSparseMatrix_ToCSR(csc, &csr)) {
    printf(">>> FAILED test -- SUNSparseMatrix_ToCSR returned nonzero\n");
    SUNMatDestroy(csc);
    return(1);
  }

  if (check_matrix(A, csr, tol)) {
    printf(">>> FAILED test --  Test_SUNSparseMatrixToCSC check_martrix failed\n");
    SUNMatDestroy(csr);
    SUNMatDestroy(csc);
    return(1);
  }

  printf("    PASSED test -- SUNSparseMatrixToCSC\n");

  SUNMatDestroy(csr);
  SUNMatDestroy(csc);

  return(0);
}


/* ----------------------------------------------------------------------
 * Check matrix
 * --------------------------------------------------------------------*/
int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  realtype *Adata, *Bdata;
  sunindextype *Aindexptrs, *Bindexptrs;
  sunindextype *Aindexvals, *Bindexvals;
  sunindextype i, ANP, BNP, Annz, Bnnz;

  /* get matrix pointers */
  Adata = SUNSparseMatrix_Data(A);
  Aindexptrs = SUNSparseMatrix_IndexPointers(A);
  Aindexvals = SUNSparseMatrix_IndexValues(A);
  ANP = SUNSparseMatrix_NP(A);
  Annz = Aindexptrs[ANP];

  Bdata = SUNSparseMatrix_Data(B);
  Bindexptrs = SUNSparseMatrix_IndexPointers(B);
  Bindexvals = SUNSparseMatrix_IndexValues(B);
  BNP = SUNSparseMatrix_NP(B);
  Bnnz = Bindexptrs[BNP];

  /* matrices must have same sparsetype, shape and actual data lengths */
  if (SUNMatGetID(A) != SUNMatGetID(B)) {
    printf(">>> ERROR: check_matrix: Different storage types (%d vs %d)\n",
           SUNMatGetID(A), SUNMatGetID(B));
    return(1);
  }
  if (SUNSparseMatrix_SparseType(A) != SUNSparseMatrix_SparseType(B)) {
    printf(">>> ERROR: check_matrix: Different storage types (%d vs %d)\n",
           SUNSparseMatrix_SparseType(A), SUNSparseMatrix_SparseType(B));
    return(1);
  }
  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Rows(B)) {
    printf(">>> ERROR: check_matrix: Different numbers of rows (%ld vs %ld)\n",
           (long int) SUNSparseMatrix_Rows(A), (long int) SUNSparseMatrix_Rows(B));
    return(1);
  }
  if (SUNSparseMatrix_Columns(A) != SUNSparseMatrix_Columns(B)) {
    printf(">>> ERROR: check_matrix: Different numbers of columns (%ld vs %ld)\n",
           (long int) SUNSparseMatrix_Columns(A),
           (long int) SUNSparseMatrix_Columns(B));
    return(1);
  }
  if (Annz != Bnnz) {
    printf(">>> ERROR: check_matrix: Different numbers of nonzeos (%ld vs %ld)\n",
           (long int) Annz, (long int) Bnnz);
    return(1);
  }

  /* compare sparsity patterns */
  for (i=0; i<ANP; i++)
    failure += (Aindexptrs[i] != Bindexptrs[i]);
  if (failure > ZERO) {
    printf(">>> ERROR: check_matrix: Different indexptrs \n");
    return(1);
  }
  for (i=0; i<Annz; i++)
    failure += (Aindexvals[i] != Bindexvals[i]);
  if (failure > ZERO) {
    printf(">>> ERROR: check_matrix: Different indexvals \n");
    return(1);
  }

  /* compare matrix values */
  for(i=0; i<Annz; i++)
    failure += SUNRCompareTol(Adata[i], Bdata[i], tol);
  if (failure > ZERO) {
    printf(">>> ERROR: check_matrix: Different entries \n");
    return(1);
  }

  return(0);
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure = 0;
  realtype *Adata;
  sunindextype *indexptrs;
  sunindextype i, NP;

  /* get data pointer */
  Adata = SUNSparseMatrix_Data(A);

  /* compare data */
  indexptrs = SUNSparseMatrix_IndexPointers(A);
  NP = SUNSparseMatrix_NP(A);
  for(i=0; i < indexptrs[NP]; i++){
    failure += SUNRCompareTol(Adata[i], val, tol);
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
  for(i=0; i < xldata; i++){
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(SUNMatrix A)
{
  realtype *Adata = SUNSparseMatrix_Data(A);
  if (Adata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  if (SUNSparseMatrix_Rows(A) == SUNSparseMatrix_Columns(A))
    return SUNTRUE;
  else
    return SUNFALSE;
}

void sync_device(SUNMatrix A)
{
  /* not running on GPU, just return */
  return;
}
