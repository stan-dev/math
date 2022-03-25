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
 * SUNMATRIX_CUSPARSE unit tests.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>
#include <stdio.h>

#include <cuda_runtime.h>

#include <nvector/nvector_cuda.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_cusparse.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "test_sunmatrix.h"
#include "dreadrb.h"

enum { IDENTITY, RANDOM, RBFILE };

/* Implementation specific test of SUNMatrix_cuSparse_SetKernelExecPolicy */
int Test_SetKernelExecPolicy(SUNMatrix A, int myid);

class ATestExecPolicy : public SUNCudaExecPolicy
{
public:
  ATestExecPolicy() : stream_(0) {}

  virtual size_t gridSize(size_t numWorkElements = 0, size_t blockDim = 0) const
  {
    return 1;
  }

  virtual size_t blockSize(size_t numWorkElements = 0, size_t gridDim = 0) const
  {
    return 1;
  }

  virtual const cudaStream_t* stream() const
  {
    return &stream_;
  }

  virtual SUNCudaExecPolicy* clone() const
  {
    return static_cast<SUNCudaExecPolicy*>(new ATestExecPolicy());
  }

private:
  const cudaStream_t stream_;
};

static SUNContext sunctx;

 /* ----------------------------------------------------------------------
  * Main SUNMatrix Testing Routine
  * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int          fails=0;                    /* counter for test failures  */
  sunindextype M, N;                       /* overall matrix dims        */
  sunindextype blkrows, blkcols;           /* block matrix dims          */
  int          nblocks;                    /* number of matrix blocks    */
  int          block_nnz_max;              /* max number of nnz in block */
  int          mattype;                    /* matrix storage type        */
  N_Vector     x, y, d_x, d_y;             /* test vectors               */
  realtype*    vecdata;                    /* pointers to vector data    */
  SUNMatrix    A, B, C, D, dA, dB, dI;     /* test matrices              */
  realtype*    matdata;                    /* pointer to matrix data     */
  int          print_timing, square;
  int          matrix_to_use;
  sunindextype i, j;
  FILE*        matrixfp;
  char*        filename;
  cusparseStatus_t cusp_status;
  cusparseHandle_t cusp_handle;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  /* initialize some input variables */
  blkrows = 0;
  blkcols = 0;
  nblocks = 0;
  square  = 0;

  /* check input */
  if (argc < 7) {
    printf("ERROR: SIX (6) inputs required: matrix (filename|random|identity), matrix rows, matrix cols, number of blocks, matrix type (CSR/BCSR), print timing (0/1)\n");
    return(-1);
  }

  /* determine what test matrix to use */
  if (!strcmp(argv[1], "random")) {
    matrix_to_use = RANDOM;
  } else if (!strcmp(argv[1], "identity")) {
    matrix_to_use = IDENTITY;
  } else {
    matrix_to_use = RBFILE;
    filename = argv[1];
  }

  /* if we are not reading from a file, verify that the dimension args are legal */
  if (matrix_to_use != RBFILE) {
    blkrows = (sunindextype) atol(argv[2]);
    if (blkrows <= 0) {
      printf("ERROR: number of rows must be a positive integer\n");
      return(-1);
    }

    blkcols = (sunindextype) atol(argv[3]);
    if (blkcols <= 0) {
      printf("ERROR: number of cols must be a positive integer\n");
      return(-1);
    }

    square = (blkrows == blkcols) ? 1 : 0;
  }

  nblocks = (sunindextype) atol(argv[4]);
  if (nblocks < 1) {
    printf("ERROR: number of blocks must be a positive integer\n");
    return(-1);
  }

  if (!strcmp(argv[5], "CSR")) {
    mattype = SUNMAT_CUSPARSE_CSR;
    if (nblocks != 1) {
        printf("ERROR: the CSR format only supports 1 block\n");
        return(-1);
    }
  } else if (!strcmp(argv[5], "BCSR")) {
    mattype = SUNMAT_CUSPARSE_BCSR;
    if (matrix_to_use == RBFILE) {
        printf("ERROR: cannot read BCSR format from a file\n");
    }
    if (!square) {
        printf("ERROR: the BCSR format only supports square block matrices\n");
        return(-1);
    }
  } else {
    printf("ERROR: matrix type must be CSR or BCSR\n");
    return(-1);
  }

  print_timing = atoi(argv[6]);
  SetTiming(print_timing);

  /* Initialize cuSPARSE */
  cusp_status = cusparseCreate(&cusp_handle);
  if (cusp_status != CUSPARSE_STATUS_SUCCESS) {
    printf("ERROR: could not create cuSPARSE handle\n");
    return(-1);
  }

  /* Initialize vectors and matrices to NULL */
  x  = NULL;
  y  = NULL;
  A  = NULL;
  B  = NULL;
  C  = NULL;
  D  = NULL;
  dA = NULL;
  dB = NULL;
  dI = NULL;

  if (matrix_to_use == RANDOM) {
    M = blkrows * nblocks;
    N = blkcols * nblocks;
    block_nnz_max = blkrows*blkcols / 2;

    /* Create sparsity pattern for a block. */
    sunindextype *cols = (sunindextype *) malloc(block_nnz_max*sizeof(sunindextype));
    sunindextype *rows = (sunindextype *) malloc(block_nnz_max*sizeof(sunindextype));
    for (i=0; i<block_nnz_max; i++) {
        cols[i] = rand() % blkcols;
        rows[i] = rand() % blkrows;
    }

    /* Fill matrix with uniform random data in [0,1/N] */
    D = SUNDenseMatrix(M, N, sunctx);
    for (i=0; i<nblocks; i++) {
        for (j=0; j<block_nnz_max; j++) {
          sunindextype col = cols[j] + blkcols*i;
          sunindextype row = rows[j] + blkrows*i;
          matdata = SUNDenseMatrix_Column(D,col);
          matdata[row] = (realtype) rand() / (realtype) RAND_MAX / N;
        }
    }
    if (SUNMatScaleAddI(RCONST(1.0), D)) {
      printf("ERROR: SUNMatScaleAddI failed for dense matrix D\n");
      return(-1);
    }

    /* Fill matrix with uniform random data in [0,1/N] */
    C = SUNDenseMatrix(M, N, sunctx);
    for (i=0; i<nblocks; i++) {
        for (j=0; j<block_nnz_max; j++) {
          sunindextype col = cols[j] + blkcols*i;
          sunindextype row = rows[j] + blkrows*i;
          matdata = SUNDenseMatrix_Column(C,col);
          matdata[row] = (realtype) rand() / (realtype) RAND_MAX / N;
        }
    }
    if (SUNMatScaleAddI(RCONST(1.0), C)) {
      printf("ERROR: SUNMatScaleAddI failed for dense matrix C\n");
      return(-1);
    }

    free(cols);
    free(rows);

    /* Create sparse matrices from dense */
    A = SUNSparseFromDenseMatrix(C, ZERO, CSR_MAT);
    if (A == NULL) {
      printf("ERROR: SUNSparseFromDenseMatrix returned NULL for A\n");
      return(-1);
    }
    B = SUNSparseFromDenseMatrix(D, ZERO, CSR_MAT);
    if (B == NULL) {
      printf("ERROR: SUNSparseFromDenseMatrix returned NULL B\n");
      return(-1);
    }
  } else if (matrix_to_use == IDENTITY) {
    M = blkrows * nblocks;
    N = blkcols * nblocks;

    D = SUNDenseMatrix(M, N, sunctx);
    SUNMatScaleAddI(RCONST(0.0), D);
    if (SUNMatScaleAddI(RCONST(0.0), D)) {
      printf("ERROR: SUNMatScaleAddI failed for dense matrix D\n");
      return(-1);
    }

    C = SUNDenseMatrix(M, N, sunctx);
    if (SUNMatScaleAddI(RCONST(0.0), C)) {
      printf("ERROR: SUNMatScaleAddI failed for dense matrix C\n");
      return(-1);
    }

    /* Create sparse matrices from dense */
    A = SUNSparseFromDenseMatrix(C, ZERO, CSR_MAT);
    if (A == NULL) {
      printf("ERROR: SUNSparseFromDenseMatrix returned NULL for A\n");
      return(-1);
    }
    B = SUNSparseFromDenseMatrix(D, ZERO, CSR_MAT);
    if (B == NULL) {
      printf("ERROR: SUNSparseFromDenseMatrix returned NULL B\n");
      return(-1);
    }
  } else {
    SUNMatrix cscA;

    matrixfp = fopen(filename, "r");
    dreadrb_dist(0, matrixfp, &cscA, sunctx);
    fclose(matrixfp);

    if (SUNSparseMatrix_ToCSR(cscA, &A)) {
      printf("ERROR: cannot convert matrix that was read to CSR\n");
      return(-1);
    }
    SUNMatDestroy(cscA);

    if (SUNMatScaleAddI(RCONST(1.0), A)) {
      printf("ERROR: SUNMatScaleAddI failed on matrix that read\n");
      return(-1);
    }

    blkrows = SUNSparseMatrix_Rows(A);
    blkcols = SUNSparseMatrix_Columns(A);
    square = (blkrows == blkcols) ? 1 : 0;
    nblocks = 1;
    M = blkrows * nblocks;
    N = blkcols * nblocks;

    B = SUNMatClone(A);
    if (B == NULL || (SUNMatCopy(A, B) != 0)) {
      printf("ERROR: failed to SUNMatClone and SUNMatCopy\n");
      return(-1);
    }
  }

  printf("cuSPARSE SUNMatrix test: size %ld by %ld, nblocks %ld, block size %ld by %ld, format = %i\n\n",
  (long int) M, (long int) N, (long int) nblocks, (long int) blkrows, (long int) blkcols, mattype);

  if (mattype == SUNMAT_CUSPARSE_CSR) {
    /* Create matrices that will be on the device */
    dA = SUNMatrix_cuSparse_NewCSR(SM_ROWS_S(A), SM_COLUMNS_S(A), SM_NNZ_S(A), cusp_handle, sunctx);
    if (dA == NULL) {
      printf("ERROR: SUNMatrix_cuSparse_NewCSR returned NULL for dA\n");
      return(-1);
    }
    dB = SUNMatrix_cuSparse_NewCSR(SM_ROWS_S(B), SM_COLUMNS_S(B), SM_NNZ_S(B), cusp_handle, sunctx);
    if (dB == NULL) {
      printf("ERROR: SUNMatrix_cuSparse_NewCSR returned NULL for dB\n");
      return(-1);
    }
  } else if (mattype == SUNMAT_CUSPARSE_BCSR) {
    sunindextype block_nnz;

    /* Calculate actual number of nonzeros per block */
    block_nnz = SUNSparseMatrix_NNZ(A) / nblocks;

    /* Create matrices that will be on the device */
    dA = SUNMatrix_cuSparse_NewBlockCSR(nblocks, blkrows, blkrows, block_nnz, cusp_handle, sunctx);
    if (dA == NULL) {
      printf("ERROR: SUNMatrix_cuSparse_NewCSR returned NULL for dA\n");
      return(-1);
    }
    dB = SUNMatrix_cuSparse_NewBlockCSR(nblocks, blkrows, blkrows, block_nnz, cusp_handle, sunctx);
    if (dB == NULL) {
      printf("ERROR: SUNMatrix_cuSparse_NewCSR returned NULL for dB\n");
      return(-1);
    }
  } else {
    printf("ERROR: unknown mattype\n");
    return(-1);
  }

  /* Copy data to device */
  fails = SUNMatrix_cuSparse_CopyToDevice(dA, SM_DATA_S(A), SM_INDEXPTRS_S(A), SM_INDEXVALS_S(A));
  if (fails != 0) {
    printf("ERROR: could not copy A to the device\n");
    return(-1);
  }
  fails = SUNMatrix_cuSparse_CopyToDevice(dB, SM_DATA_S(B), SM_INDEXPTRS_S(B), SM_INDEXVALS_S(B));
  if (fails != 0) {
    printf("ERROR: could not copy B to the device\n");
    return(-1);
  }

  /* Create/fill I matrix */
  dI = NULL;
  if (square) {
    dI = SUNMatClone_cuSparse(dA);
    if (dI == NULL) {
      printf("ERROR: SUNMatClone_cuSparse returned NULL\n");
      return(-1);
    }
    if (SUNMatCopy_cuSparse(dA, dI)) {
      printf("ERROR: SUNMatCopy_cuSparse failed\n");
      return(-1);
    }
    if (SUNMatScaleAddI_cuSparse(ZERO, dI)) {
      printf("ERROR: SUNMatScaleAddI_cuSparse failed\n");
      return(-1);
    }
  }

  /* Create vectors */
  d_x = N_VNew_Cuda(N, sunctx);
  d_y = N_VNew_Cuda(M, sunctx);
  if (d_x == NULL || d_y == NULL) {
    printf("ERROR: N_VNew_Cuda returned NULL\n");
    return(-1);
  }
  x = N_VMake_Serial(N, N_VGetHostArrayPointer_Cuda(d_x), sunctx);
  y = N_VMake_Serial(M, N_VGetHostArrayPointer_Cuda(d_y), sunctx);
  if (x == NULL || y == NULL) {
    printf("ERROR: N_VMake_Serial returned NULL\n");
    return(-1);
  }

  /* Zero the vectors on the host */
  N_VConst(ZERO, x);
  N_VConst(ZERO, y);

  /* Fill vector on the host */
  vecdata = N_VGetArrayPointer(x);
  for(i=0; i<N; i++)
    vecdata[i] = (realtype) rand() / (realtype) RAND_MAX;

  /* Compute reference y on the host */
  if (SUNMatMatvec(A, x, y)) {
    printf("FAIL: SUNSparseMatrix matvec failure \n \n");
    SUNMatDestroy(A);  SUNMatDestroy(B);
    SUNMatDestroy(C);  SUNMatDestroy(D);
    SUNMatDestroy(dA); SUNMatDestroy(dB);
    N_VDestroy(x);  N_VDestroy(y);
    N_VDestroy(d_x); N_VDestroy(d_y);
    if (square) {
        SUNMatDestroy(dI);
    }
    return(1);
  }

  /* Copy vectors to the device */
  N_VCopyToDevice_Cuda(d_x);
  N_VCopyToDevice_Cuda(d_y);

  printf("Setup complete\n");
  printf("Beginning tests\n\n");

  /* SUNMatrix Tests */
  fails += Test_SUNMatGetID(dA, SUNMATRIX_CUSPARSE, 0);
  fails += Test_SUNMatClone(dA, 0);
  fails += Test_SUNMatCopy(dA, 0);
  fails += Test_SUNMatZero(dA, 0);
  fails += Test_SUNMatScaleAdd(dA, dI, 0);
  if (square) fails += Test_SUNMatScaleAddI(dA, dI, 0);
  fails += Test_SUNMatMatvec(dA, d_x, d_y, 0);
  if (square) fails += Test_SetKernelExecPolicy(dI, 0);

  /* Print result */
  if (fails) {
    SUNMatrix_cuSparse_CopyFromDevice(dA, SM_DATA_S(A), NULL, NULL);
    SUNMatrix_cuSparse_CopyFromDevice(dB, SM_DATA_S(B), NULL, NULL);
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nB =\n");
    SUNSparseMatrix_Print(B,stdout);
    N_VCopyFromDevice_Cuda(d_x);
    N_VCopyFromDevice_Cuda(d_y);
    printf("\nx\n");
    N_VPrint_Cuda(d_x);
    printf("\ny = Ax (reference)\n");
    N_VPrint_Cuda(d_y);
  } else {
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");
  }

  printf("Beginning teardown\n");

  /* Free vectors and matrices */
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(d_x);
  N_VDestroy(d_y);
  SUNMatDestroy(A);
  SUNMatDestroy(B);
  SUNMatDestroy(C);
  SUNMatDestroy(D);
  SUNMatDestroy(dA);
  SUNMatDestroy(dB);
  if (square) {
    SUNMatDestroy(dI);
  }

  cusparseDestroy(cusp_handle);
  SUNContext_Free(&sunctx);

  printf("Teardown complete\n");

  return(fails);
 }

 /* ----------------------------------------------------------------------
  * Test the SUNMatrix_cuSparse_SetKernelExecPolicy function.
  * --------------------------------------------------------------------*/
int Test_SetKernelExecPolicy(SUNMatrix I, int myid)
{
  int print_all_ranks = 0;
  realtype  tol = 100*UNIT_ROUNDOFF;
  SUNMatrix B = SUNMatClone(I);

  /* check cloned matrix */
  if (B == NULL) {
    TEST_STATUS(">>> FAILED test -- SetKernelExecPolicy \n", myid);
    TEST_STATUS("    After SUNMatClone, B == NULL \n \n", myid);
    return(1);
  }

  /* copy data */
  if (SUNMatCopy(I, B)) {
    TEST_STATUS(">>> FAILED test -- SetKernelExecPolicy \n", myid);
    TEST_STATUS("    SUNMatCopy returned nonzero \n \n", myid);
    SUNMatDestroy(B);
    return(1);
  }

  /* set kernel exec policy */
  ATestExecPolicy exec_policy;
  SUNMatrix_cuSparse_SetKernelExecPolicy(B, &exec_policy);

  /* try out an operation */
  if (SUNMatScaleAddI(RCONST(-1.0), B)) {
    TEST_STATUS(">>> FAILED test -- SetKernelExecPolicy \n", myid);
    TEST_STATUS("    SUNMatScaleAddI returned nonzero \n \n", myid);
    SUNMatDestroy(B);
    return(1);
  }

  /* check matrix */
  if (check_matrix_entry(B, ZERO, tol)) {
    TEST_STATUS(">>> FAILED test -- SetKernelExecPolicy \n", myid);
    TEST_STATUS("    check_matrix_entry returned nonzero \n \n", myid);
    SUNMatDestroy(B);
    return(1);
  }

  TEST_STATUS("    PASSED test -- SetKernelExecPolicy \n", myid);

  SUNMatDestroy(B);

  return 0;
}

 /* ----------------------------------------------------------------------
  * Check matrix
  * --------------------------------------------------------------------*/
 int check_matrix(SUNMatrix dA, SUNMatrix dB, realtype tol)
 {
   int failure = 0;
   SUNMatrix A, B;
   realtype *Adata, *Bdata;
   sunindextype *Aindexptrs, *Bindexptrs;
   sunindextype *Aindexvals, *Bindexvals;
   sunindextype i, ANP, Annz, Bnnz;

   /* copy matrix data to host for the checks */
   A = SUNSparseMatrix(SUNMatrix_cuSparse_Rows(dA), SUNMatrix_cuSparse_Columns(dA),
                       SUNMatrix_cuSparse_NNZ(dA), CSR_MAT, sunctx);
   B = SUNSparseMatrix(SUNMatrix_cuSparse_Rows(dB), SUNMatrix_cuSparse_Columns(dB),
                       SUNMatrix_cuSparse_NNZ(dB), CSR_MAT, sunctx);

   failure = SUNMatrix_cuSparse_CopyFromDevice(dA, SM_DATA_S(A),
                                               SM_INDEXPTRS_S(A),
                                               SM_INDEXVALS_S(A));
   failure = SUNMatrix_cuSparse_CopyFromDevice(dB, SM_DATA_S(B),
                                               SM_INDEXPTRS_S(B),
                                               SM_INDEXVALS_S(B));
   cudaDeviceSynchronize();

   /* get matrix pointers */
   Adata = SUNSparseMatrix_Data(A);
   Aindexptrs = SUNSparseMatrix_IndexPointers(A);
   Aindexvals = SUNSparseMatrix_IndexValues(A);
   ANP = SUNSparseMatrix_NP(A);
   Annz = SUNSparseMatrix_NNZ(A);

   Bdata = SUNSparseMatrix_Data(B);
   Bindexptrs = SUNSparseMatrix_IndexPointers(B);
   Bindexvals = SUNSparseMatrix_IndexValues(B);
   Bnnz = SUNSparseMatrix_NNZ(B);

   /* matrices must have same sparsetype, shape and actual data lengths */
   if (SUNMatGetID(dA) != SUNMatGetID(dB)) {
     printf(">>> ERROR: check_matrix: Different storage types (%d vs %d)\n",
            SUNMatGetID(dA), SUNMatGetID(dB));
     SUNMatDestroy(dA); SUNMatDestroy(dB);
     return(1);
   }
   if (SUNMatrix_cuSparse_SparseType(A) != SUNMatrix_cuSparse_SparseType(B)) {
     printf(">>> ERROR: check_matrix: Different storage types (%d vs %d)\n",
            SUNMatrix_cuSparse_SparseType(A), SUNMatrix_cuSparse_SparseType(B));
     SUNMatDestroy(A); SUNMatDestroy(B);
     return(1);
   }
   if (SUNMatrix_cuSparse_Rows(dA) != SUNMatrix_cuSparse_Rows(dB)) {
     printf(">>> ERROR: check_matrix: Different numbers of rows (%ld vs %ld)\n",
            (long int) SUNMatrix_cuSparse_Rows(dA), (long int) SUNMatrix_cuSparse_Rows(dB));
     SUNMatDestroy(A); SUNMatDestroy(B);
     return(1);
   }
   if (SUNMatrix_cuSparse_Columns(dA) != SUNMatrix_cuSparse_Columns(dB)) {
     printf(">>> ERROR: check_matrix: Different numbers of columns (%ld vs %ld)\n",
            (long int) SUNMatrix_cuSparse_Columns(dA),
            (long int) SUNMatrix_cuSparse_Columns(dB));
     SUNMatDestroy(A); SUNMatDestroy(B);
     return(1);
   }
   if (Annz != Bnnz) {
     printf(">>> ERROR: check_matrix: Different numbers of nonzeros (%ld vs %ld)\n",
            (long int) Annz, (long int) Bnnz);
     SUNMatDestroy(A); SUNMatDestroy(B);
     return(1);
   }

   /* compare sparsity patterns */
   for (i=0; i<ANP; i++)
     failure += (Aindexptrs[i] != Bindexptrs[i]);
   if (failure > ZERO) {
     printf(">>> ERROR: check_matrix: Different indexptrs \n");
     SUNMatDestroy(A); SUNMatDestroy(B);
     return(1);
   }
   for (i=0; i<Annz; i++)
     failure += (Aindexvals[i] != Bindexvals[i]);
   if (failure > ZERO) {
     printf(">>> ERROR: check_matrix: Different indexvals \n");
     SUNMatDestroy(A); SUNMatDestroy(B);
     return(1);
   }

   /* compare matrix values */
   for(i=0; i<Annz; i++)
     failure += SUNRCompareTol(Adata[i], Bdata[i], tol);
   if (failure > ZERO) {
     printf(">>> ERROR: check_matrix: Different entries \n");
     SUNMatDestroy(A); SUNMatDestroy(B);
     return(1);
   }

   SUNMatDestroy(A); SUNMatDestroy(B);

   return(0);
 }

 int check_matrix_entry(SUNMatrix dA, realtype val, realtype tol)
 {
   int failure = 0;
   realtype *Adata;
   sunindextype i;

   /* copy matrix data to host for the checks */
   Adata = (realtype*) malloc(SUNMatrix_cuSparse_NNZ(dA)*sizeof(realtype));
   failure = SUNMatrix_cuSparse_CopyFromDevice(dA, Adata, NULL, NULL);
   cudaDeviceSynchronize();

   /* compare data */
   for(i=0; i < SUNMatrix_cuSparse_NNZ(dA); i++) {
     failure += SUNRCompareTol(Adata[i], val, tol);
   }

   free(Adata);

   if (failure > ZERO)
     return(1);
   else
     return(0);
 }

 int check_vector(N_Vector expected, N_Vector computed, realtype tol)
 {
   int failure = 0;
   realtype *xdata, *ydata;
   sunindextype xldata, yldata;
   sunindextype i;

   /* get vector data */
   xdata = N_VGetHostArrayPointer_Cuda(expected);
   ydata = N_VGetHostArrayPointer_Cuda(computed);

   /* copy data to host */
   N_VCopyFromDevice_Cuda(expected);
   N_VCopyFromDevice_Cuda(computed);
   cudaDeviceSynchronize();

   /* check data lengths */
   xldata = N_VGetLength_Cuda(expected);
   yldata = N_VGetLength_Cuda(computed);

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
   realtype *Adata = SUNMatrix_cuSparse_Data(A);
   if (Adata == NULL)
     return SUNFALSE;
   else
     return SUNTRUE;
 }

 booleantype is_square(SUNMatrix A)
 {
   if (SUNMatrix_cuSparse_Rows(A) == SUNMatrix_cuSparse_Columns(A))
     return SUNTRUE;
   else
     return SUNFALSE;
 }

void sync_device(SUNMatrix A)
{
  cudaDeviceSynchronize();
}
