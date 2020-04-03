/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This is the header file is for the cuSPARSE implementation of the
 * SUNMATRIX module.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_cuda.h>
#include <sunmatrix/sunmatrix_cusparse.h>

#include "sundials_cuda.h"
#include "sundials_debug.h"
#include "cusparse_kernels.cuh"

/* Use the namespace for the kernels */
using namespace sundials::device::sunmatrix_cusparse;

/* Constants */
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

#define MAX_THREAD_PER_BLOCK(val) ( (val > 16*CUDA_WARP_SIZE) ? (16*CUDA_WARP_SIZE) : (val) )

/* Private function prototypes */
static booleantype SMCompatible_cuSparse(SUNMatrix A, SUNMatrix B);
static SUNMatrix SUNMatrix_cuSparse_NewEmpty();

/* Macros for handling the different function names based on precision */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define cusparseXcsrmv cusparseDcsrmv
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define cusparseXcsrmv cusparseScsrmv
#endif

/* Content accessor macros */
#define SMCU_CONTENT_S(A)     ( (SUNMatrix_Content_cuSparse)(A->content) )
#define SMCU_ROWS_S(A)        ( SMCU_CONTENT_S(A)->M )
#define SMCU_COLUMNS_S(A)     ( SMCU_CONTENT_S(A)->N )
#define SMCU_NNZ_S(A)         ( SMCU_CONTENT_S(A)->NNZ )
#define SMCU_NBLOCKS_S(A)     ( SMCU_CONTENT_S(A)->nblocks )
#define SMCU_BLOCKROWS_S(A)   ( SMCU_CONTENT_S(A)->blockrows )
#define SMCU_BLOCKCOLS_S(A)   ( SMCU_CONTENT_S(A)->blockcols )
#define SMCU_BLOCKNNZ_S(A)    ( SMCU_CONTENT_S(A)->blocknnz )
#define SMCU_NP_S(A)          ( SMCU_CONTENT_S(A)->NP )
#define SMCU_SPARSETYPE_S(A)  ( SMCU_CONTENT_S(A)->sparse_type )
#define SMCU_OWNDATA_S(A)     ( SMCU_CONTENT_S(A)->own_data )
#define SMCU_DATA_S(A)        ( SMCU_CONTENT_S(A)->data )
#define SMCU_INDEXVALS_S(A)   ( SMCU_CONTENT_S(A)->colind )
#define SMCU_INDEXPTRS_S(A)   ( SMCU_CONTENT_S(A)->rowptrs )
#define SMCU_MATDESCR_S(A)    ( SMCU_CONTENT_S(A)->mat_descr )
#define SMCU_CUSPHANDLE_S(A)  ( SMCU_CONTENT_S(A)->cusp_handle )
#define SMCU_FIXEDPATTERN_S(A)( SMCU_CONTENT_S(A)->fixed_pattern )


/* ------------------------------------------------------------------
 * Constructors.
 * ------------------------------------------------------------------ */


SUNMatrix SUNMatrix_cuSparse_NewCSR(int M, int N, int NNZ, cusparseHandle_t cusp)
{
  /* return with NULL matrix on illegal input */
  if ( (M <= 0) || (N <= 0) || (NNZ < 0) )
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_NewCSR_cuSparse: illegal value(s) for M, N, or NNZ\n");
    return NULL;
  }

  SUNMatrix A = SUNMatrix_cuSparse_NewEmpty();
  if (A == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_NewCSR_cuSparse: SUNMatrix_cuSparse_NewEmpty returned NULL\n");
    return NULL;
  }

  /* Allocate device memory for the matrix */
  int *d_colind, *d_rowptr;
  realtype *d_values;

  d_colind = NULL;
  d_rowptr = NULL;
  d_values = NULL;

  cudaError_t cuerr;
  cuerr = cudaMalloc((void **) &d_colind, sizeof(*d_colind) * NNZ);
  if (!SUNDIALS_CUDA_VERIFY(cuerr))
  {
    SUNMatDestroy(A);
    return NULL;
  }
  cuerr = cudaMalloc((void **) &d_rowptr, sizeof(*d_rowptr) * (M+1));
  if (!SUNDIALS_CUDA_VERIFY(cuerr))
  {
    SUNMatDestroy(A);
    cudaFree(d_colind);
    return NULL;
  }
  cuerr = cudaMalloc((void **) &d_values, sizeof(*d_values) * NNZ);
  if (!SUNDIALS_CUDA_VERIFY(cuerr))
  {
    SUNMatDestroy(A);
    cudaFree(d_colind);
    cudaFree(d_rowptr);
    return NULL;
  }

  /* Choose sensible defaults */
  cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
  cusparseMatDescr_t mat_descr;
  cusparse_status = cusparseCreateMatDescr(&mat_descr);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNMatDestroy(A);
    return NULL;
  }

  cusparse_status = cusparseSetMatType(mat_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNMatDestroy(A);
    cusparseDestroyMatDescr(mat_descr);
    return NULL;
  }

  cusparse_status = cusparseSetMatIndexBase(mat_descr, CUSPARSE_INDEX_BASE_ZERO);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNMatDestroy(A);
    cusparseDestroyMatDescr(mat_descr);
    return NULL;
  }

  /* Fill the content */
  SMCU_CONTENT_S(A)->M             = M;
  SMCU_CONTENT_S(A)->N             = N;
  SMCU_CONTENT_S(A)->NNZ           = NNZ;
  SMCU_CONTENT_S(A)->nblocks       = 1;
  SMCU_CONTENT_S(A)->blockrows     = M;
  SMCU_CONTENT_S(A)->blockcols     = N;
  SMCU_CONTENT_S(A)->blocknnz      = NNZ;
  SMCU_CONTENT_S(A)->own_data      = SUNTRUE;
  SMCU_CONTENT_S(A)->sparse_type   = SUNMAT_CUSPARSE_CSR;
  SMCU_CONTENT_S(A)->colind        = d_colind;
  SMCU_CONTENT_S(A)->rowptrs       = d_rowptr;
  SMCU_CONTENT_S(A)->data          = d_values;
  SMCU_CONTENT_S(A)->mat_descr     = mat_descr;
  SMCU_CONTENT_S(A)->cusp_handle   = cusp;
  SMCU_CONTENT_S(A)->fixed_pattern = SUNFALSE;

  return A;
}


SUNMatrix SUNMatrix_cuSparse_MakeCSR(cusparseMatDescr_t mat_descr, int M, int N, int NNZ,
                                     int *rowptrs , int *colind , realtype *data,
                                     cusparseHandle_t cusp)
{
  /* return with NULL matrix on illegal input */
  if ( (M <= 0) || (N <= 0) || (NNZ < 0) )
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_MakeCSR_cuSparse: illegal value(s) for M, N, or NNZ\n");
    return NULL;
  }

  if ( (rowptrs == NULL) || (colind == NULL) || (data == NULL) )
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_MakeCSR_cuSparse: rowptrs, colind, or data is NULL\n");
    return NULL;
  }

  if (cusparseGetMatIndexBase(mat_descr) != CUSPARSE_INDEX_BASE_ZERO)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_MakeCSR_cuSparse: the cusparseMatDescr_t must have index base CUSPARSE_INDEX_BASE_ZERO\n");
    return NULL;
  }

  SUNMatrix A = SUNMatrix_cuSparse_NewEmpty();
  if (A == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_MakeCSR_cuSparse: SUNMatrix_cuSparse_NewEmpty returned NULL\n");
    return NULL;
  }

  /* Fill content */
  SMCU_CONTENT_S(A)->M             = M;
  SMCU_CONTENT_S(A)->N             = N;
  SMCU_CONTENT_S(A)->NNZ           = NNZ;
  SMCU_CONTENT_S(A)->nblocks       = 1;
  SMCU_CONTENT_S(A)->blockrows     = M;
  SMCU_CONTENT_S(A)->blockcols     = N;
  SMCU_CONTENT_S(A)->blocknnz      = NNZ;
  SMCU_CONTENT_S(A)->own_data      = SUNFALSE;
  SMCU_CONTENT_S(A)->sparse_type   = SUNMAT_CUSPARSE_CSR;
  SMCU_CONTENT_S(A)->colind        = colind;
  SMCU_CONTENT_S(A)->rowptrs       = rowptrs;
  SMCU_CONTENT_S(A)->data          = data;
  SMCU_CONTENT_S(A)->mat_descr     = mat_descr;
  SMCU_CONTENT_S(A)->cusp_handle   = cusp;
  SMCU_CONTENT_S(A)->fixed_pattern = SUNFALSE;

  return A;
}


SUNMatrix SUNMatrix_cuSparse_NewBlockCSR(int nblocks, int blockrows, int blockcols, int blocknnz, cusparseHandle_t cusp)
{
  int M, N, NNZ;

  /* Return with NULL matrix on illegal input */
  if (blockrows != blockcols)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_NewBlockCSR: matrix must be square for the BCSR format\n");
    return NULL;
  }

  M   = nblocks * blockrows;
  N   = M;
  NNZ = nblocks * blocknnz;

  /* Return with NULL matrix on illegal input */
  if ( (M <= 0) || (N <= 0) || (NNZ < 0) )
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_NewBlockCSR: illegal value(s) for M, N, or NNZ\n");
    return NULL;
  }

  /* Allocate the SUNMatrix object */
  SUNMatrix A = SUNMatrix_cuSparse_NewEmpty();
  if (A == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_NewBlockCSR: SUNMatrix_cuSparse_NewEmpty returned NULL\n");
    return NULL;
  }

  /* Allocate device memory for the matrix */
  int *d_colind, *d_rowptr;
  realtype *d_values;

  d_colind = NULL;
  d_rowptr = NULL;
  d_values = NULL;

  cudaError_t cuerr;
  cuerr = cudaMalloc((void **) &d_colind, sizeof(*d_colind) * blocknnz);
  if (!SUNDIALS_CUDA_VERIFY(cuerr))
  {
    SUNMatDestroy(A);
    return NULL;
  }
  cuerr = cudaMalloc((void **) &d_rowptr, sizeof(*d_rowptr) * (blockrows + 1));
  if (!SUNDIALS_CUDA_VERIFY(cuerr))
  {
    SUNMatDestroy(A);
    cudaFree(d_colind);
    return NULL;
  }
  cuerr = cudaMalloc((void **) &d_values, sizeof(*d_values) * blocknnz * nblocks);
  if (!SUNDIALS_CUDA_VERIFY(cuerr))
  {
    SUNMatDestroy(A);
    cudaFree(d_colind);
    cudaFree(d_rowptr);
    return NULL;
  }

  /* Choose sensible defaults */
  cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
  cusparseMatDescr_t mat_descr;
  cusparse_status = cusparseCreateMatDescr(&mat_descr);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNMatDestroy(A);
    return NULL;
  }

  cusparse_status = cusparseSetMatType(mat_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNMatDestroy(A);
    cusparseDestroyMatDescr(mat_descr);
    return NULL;
  }

  cusparse_status = cusparseSetMatIndexBase(mat_descr, CUSPARSE_INDEX_BASE_ZERO);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNMatDestroy(A);
    cusparseDestroyMatDescr(mat_descr);
    return NULL;
  }

  /* Fill the content */
  SMCU_CONTENT_S(A)->M             = M;
  SMCU_CONTENT_S(A)->N             = N;
  SMCU_CONTENT_S(A)->NNZ           = NNZ;
  SMCU_CONTENT_S(A)->nblocks       = nblocks;
  SMCU_CONTENT_S(A)->blockrows     = blockrows;
  SMCU_CONTENT_S(A)->blockcols     = blockrows;
  SMCU_CONTENT_S(A)->blocknnz      = blocknnz;
  SMCU_CONTENT_S(A)->own_data      = SUNTRUE;
  SMCU_CONTENT_S(A)->sparse_type   = SUNMAT_CUSPARSE_BCSR;
  SMCU_CONTENT_S(A)->colind        = d_colind;
  SMCU_CONTENT_S(A)->rowptrs       = d_rowptr;
  SMCU_CONTENT_S(A)->data          = d_values;
  SMCU_CONTENT_S(A)->mat_descr     = mat_descr;
  SMCU_CONTENT_S(A)->cusp_handle   = cusp;
  SMCU_CONTENT_S(A)->fixed_pattern = SUNFALSE;

  return A;
}

/* ------------------------------------------------------------------
 * Implementation specific routines.
 * ------------------------------------------------------------------ */

int SUNMatrix_cuSparse_SparseType(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_SPARSETYPE_S(A);
  else
    return SUNMAT_ILL_INPUT;
}

int SUNMatrix_cuSparse_Rows(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_ROWS_S(A);
  else
    return SUNMAT_ILL_INPUT;
}

int SUNMatrix_cuSparse_Columns(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_COLUMNS_S(A);
  else
    return SUNMAT_ILL_INPUT;
}

int SUNMatrix_cuSparse_NNZ(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_NNZ_S(A);
  else
    return SUNMAT_ILL_INPUT;
}

int* SUNMatrix_cuSparse_IndexPointers(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_INDEXPTRS_S(A);
  else
    return NULL;
}

int* SUNMatrix_cuSparse_IndexValues(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_INDEXVALS_S(A);
  else
    return NULL;
}

realtype* SUNMatrix_cuSparse_Data(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_DATA_S(A);
  else
    return NULL;
}

int SUNMatrix_cuSparse_NumBlocks(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_NBLOCKS_S(A);
  else
    return SUNMAT_ILL_INPUT;
}

int SUNMatrix_cuSparse_BlockRows(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_BLOCKROWS_S(A);
  else
    return SUNMAT_ILL_INPUT;
}

int SUNMatrix_cuSparse_BlockColumns(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_BLOCKCOLS_S(A);
  else
    return SUNMAT_ILL_INPUT;
}

int SUNMatrix_cuSparse_BlockNNZ(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_BLOCKNNZ_S(A);
  else
    return SUNMAT_ILL_INPUT;
}

realtype* SUNMatrix_cuSparse_BlockData(SUNMatrix A, int blockidx)
{
  realtype *matdata;
  int offset;

  if (SUNMatGetID(A) != SUNMATRIX_CUSPARSE)
    return NULL;

  if (blockidx >= SMCU_NBLOCKS_S(A))
    return NULL;

  matdata = SMCU_DATA_S(A);
  offset = SMCU_BLOCKNNZ_S(A)*blockidx;

  return (&matdata[offset]);
}

cusparseMatDescr_t SUNMatrix_cuSparse_MatDescr(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return SMCU_MATDESCR_S(A);
  else
    return NULL;
}

int SUNMatrix_cuSparse_SetFixedPattern(SUNMatrix A, booleantype yesno)
{
  if (SUNMatGetID(A) != SUNMATRIX_CUSPARSE)
    return SUNMAT_ILL_INPUT;

  SMCU_FIXEDPATTERN_S(A) = yesno;

  return SUNMAT_SUCCESS;
}


int SUNMatrix_cuSparse_CopyToDevice(SUNMatrix dA, realtype* h_data,
                                    int* h_idxptrs, int* h_idxvals)
{
  cudaError_t cuerr;
  cudaStream_t stream;
  cusparseStatus_t cusparse_status;
  int nidxvals, nidxptrs;

  if (SUNMatGetID(dA) != SUNMATRIX_CUSPARSE)
    return SUNMAT_ILL_INPUT;

  cusparse_status = cusparseGetStream(SMCU_CUSPHANDLE_S(dA), &stream);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status)) return SUNMAT_OPERATION_FAIL;

  if (h_data != NULL)
  {
    cuerr = cudaMemcpyAsync(SMCU_DATA_S(dA), h_data,
                            SMCU_NNZ_S(dA)*sizeof(realtype),
                            cudaMemcpyHostToDevice, stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;
  }

  switch(SMCU_SPARSETYPE_S(dA))
  {
    case SUNMAT_CUSPARSE_CSR:
      nidxptrs = SMCU_ROWS_S(dA)+1;
      nidxvals = SMCU_NNZ_S(dA);
      break;
    case SUNMAT_CUSPARSE_BCSR:
      nidxptrs = SMCU_BLOCKROWS_S(dA)+1;
      nidxvals = SMCU_BLOCKNNZ_S(dA);
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_CopyToDevice: unrecognized sparse type\n");
      return SUNMAT_ILL_INPUT;
  }

  if (h_idxptrs != NULL)
  {
    cuerr = cudaMemcpyAsync(SMCU_INDEXPTRS_S(dA), h_idxptrs,
                            nidxptrs*sizeof(int),
                            cudaMemcpyHostToDevice, stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;
  }

  if (h_idxvals != NULL)
  {
    cuerr = cudaMemcpyAsync(SMCU_INDEXVALS_S(dA), h_idxvals,
                            nidxvals*sizeof(int),
                            cudaMemcpyHostToDevice, stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;
  }

  return SUNMAT_SUCCESS;
}


int SUNMatrix_cuSparse_CopyFromDevice(SUNMatrix dA, realtype* h_data,
                                      int* h_idxptrs, int* h_idxvals)
{
  cudaError_t cuerr;
  cudaStream_t stream;
  cusparseStatus_t cusparse_status;
  int nidxvals, nidxptrs;

  if (SUNMatGetID(dA) != SUNMATRIX_CUSPARSE)
    return SUNMAT_ILL_INPUT;

  cusparse_status = cusparseGetStream(SMCU_CUSPHANDLE_S(dA), &stream);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status)) return SUNMAT_OPERATION_FAIL;

  if (h_data != NULL)
  {
    cuerr = cudaMemcpyAsync(h_data, SMCU_DATA_S(dA),
                            SMCU_NNZ_S(dA)*sizeof(realtype),
                            cudaMemcpyDeviceToHost, stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;
  }

  switch(SMCU_SPARSETYPE_S(dA))
  {
    case SUNMAT_CUSPARSE_CSR:
      nidxptrs = SMCU_ROWS_S(dA)+1;
      nidxvals = SMCU_NNZ_S(dA);
    case SUNMAT_CUSPARSE_BCSR:
      nidxptrs = SMCU_BLOCKROWS_S(dA)+1;
      nidxvals = SMCU_BLOCKNNZ_S(dA);
  }

  if (h_idxptrs != NULL)
  {
    cuerr = cudaMemcpyAsync(h_idxptrs, SMCU_INDEXPTRS_S(dA),
                            nidxptrs*sizeof(int),
                            cudaMemcpyDeviceToHost, stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;
  }

  if (h_idxvals != NULL)
  {
    cuerr = cudaMemcpyAsync(h_idxvals, SMCU_INDEXVALS_S(dA),
                            nidxvals*sizeof(int),
                            cudaMemcpyDeviceToHost, stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;
  }

  return SUNMAT_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */


SUNMatrix_ID SUNMatGetID_cuSparse(SUNMatrix A)
{
  return SUNMATRIX_CUSPARSE;
}

/* Returns a new matrix allocated to have the same structure as A,
   but it does not copy any nonzeros, column vals, or row pointers. */
SUNMatrix SUNMatClone_cuSparse(SUNMatrix A)
{
  SUNMatrix B;

  switch (SMCU_SPARSETYPE_S(A))
  {
    case SUNMAT_CUSPARSE_CSR:
      B = SUNMatrix_cuSparse_NewCSR(SMCU_ROWS_S(A), SMCU_COLUMNS_S(A), SMCU_NNZ_S(A),
                                    SMCU_CUSPHANDLE_S(A));
      break;
    case SUNMAT_CUSPARSE_BCSR:
      B = SUNMatrix_cuSparse_NewBlockCSR(SMCU_NBLOCKS_S(A), SMCU_BLOCKROWS_S(A), SMCU_BLOCKCOLS_S(A),
                                         SMCU_BLOCKNNZ_S(A), SMCU_CUSPHANDLE_S(A));
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMatClone_cuSparse: sparse type not recognized\n");
      B = NULL;
  }

  SMCU_FIXEDPATTERN_S(B) = SMCU_FIXEDPATTERN_S(A);

  return B;
}


/* Deallocates the SUNMatrix object and all data it owns */
void SUNMatDestroy_cuSparse(SUNMatrix A)
{
  if (A == NULL) return;

  /* free content */
  if (A->content != NULL)
  {
    if (SMCU_OWNDATA_S(A))
    {
      /* free data array */
      if (SMCU_DATA_S(A))
      {
        cudaFree(SMCU_DATA_S(A));
        SMCU_DATA_S(A) = NULL;
      }

      /* free index values array */
      if (SMCU_INDEXVALS_S(A))
      {
        cudaFree(SMCU_INDEXVALS_S(A));
        SMCU_INDEXVALS_S(A) = NULL;
      }

      /* free index pointers array */
      if (SMCU_INDEXPTRS_S(A))
      {
        cudaFree(SMCU_INDEXPTRS_S(A));
        SMCU_INDEXPTRS_S(A) = NULL;
      }

      /* free cusparseMatDescr_t */
      cusparseDestroyMatDescr(SMCU_MATDESCR_S(A));
    }

    /* free content struct */
    free(A->content);
    A->content = NULL;
  }

  /* free ops and matrix */
  if (A->ops) { free(A->ops); A->ops = NULL; }
  free(A); A = NULL;

  return;
}


/* Performs A_ij = 0 */
int SUNMatZero_cuSparse(SUNMatrix A)
{
  cudaError_t cuerr;
  cudaStream_t stream;

  cusparseGetStream(SMCU_CUSPHANDLE_S(A), &stream);

  /* set all data to zero */
  cuerr = cudaMemsetAsync(SMCU_DATA_S(A), 0, SMCU_NNZ_S(A)*sizeof(realtype), stream);
  if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;

  /* set all rowptrs to zero unless the sparsity pattern is fixed */
  if (!SMCU_FIXEDPATTERN_S(A))
  {
    cuerr = cudaMemsetAsync(SMCU_INDEXPTRS_S(A), 0,
                            (SMCU_BLOCKROWS_S(A)+1)*sizeof(int),
                            stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;

    /* set all colind to zero */
    cuerr = cudaMemsetAsync(SMCU_INDEXVALS_S(A), 0,
                            SMCU_BLOCKNNZ_S(A)*sizeof(int),
                            stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;
  }

  return SUNMAT_SUCCESS;
}


/* Copies the nonzeros, column vals, and row pointers into dst */
int SUNMatCopy_cuSparse(SUNMatrix src, SUNMatrix dst)
{
  cudaError_t cuerr;
  cudaStream_t stream;

  /* Verify that src and dst are compatible */
  if (!SMCompatible_cuSparse(src, dst))
    return SUNMAT_ILL_INPUT;

  cusparseGetStream(SMCU_CUSPHANDLE_S(src), &stream);

  /* Ensure that dst is allocated with at least as
     much memory as we have nonzeros in src */
  if (SMCU_NNZ_S(dst) < SMCU_NNZ_S(src))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatCopy_cuSparse: the destination matrix has less nonzeros than the source\n");
    return SUNMAT_ILL_INPUT;
  }

  /* Zero out dst so that copy works correctly */
  if (SUNMatZero_cuSparse(dst) != SUNMAT_SUCCESS)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatCopy_cuSparse: SUNMatZero_cuSparse failed\n");
    return SUNMAT_OPERATION_FAIL;
  }

  /* Copy the data over */
  cuerr = cudaMemcpyAsync(SMCU_DATA_S(dst), SMCU_DATA_S(src),
                          SMCU_NNZ_S(src)*sizeof(*SMCU_DATA_S(src)),
                          cudaMemcpyDeviceToDevice, stream);
  if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;

  /* Copy the row pointers over */
  cuerr = cudaMemcpyAsync(SMCU_INDEXPTRS_S(dst), SMCU_INDEXPTRS_S(src),
                          (SMCU_BLOCKROWS_S(src)+1)*sizeof(*SMCU_INDEXPTRS_S(src)),
                          cudaMemcpyDeviceToDevice, stream);
  if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;

  /* Copy the column indices over */
  cuerr = cudaMemcpyAsync(SMCU_INDEXVALS_S(dst), SMCU_INDEXVALS_S(src),
                          SMCU_BLOCKNNZ_S(src)*sizeof(*SMCU_INDEXVALS_S(src)),
                          cudaMemcpyDeviceToDevice, stream);
  if (!SUNDIALS_CUDA_VERIFY(cuerr)) return SUNMAT_OPERATION_FAIL;

  return SUNMAT_SUCCESS;
}


/* Performs A = cA + I. Requires the diagonal to be allocated already. */
int SUNMatScaleAddI_cuSparse(realtype c, SUNMatrix A)
{
  cudaStream_t stream;
  cusparseStatus_t cusparse_status;

  cusparse_status = cusparseGetStream(SMCU_CUSPHANDLE_S(A), &stream);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status)) return SUNMAT_OPERATION_FAIL;

  unsigned threadsPerBlock, gridSize;
  switch (SMCU_SPARSETYPE_S(A))
  {
    case SUNMAT_CUSPARSE_CSR:
      /* Choose the grid size to be the number of rows in the matrix,
        and then choose threadsPerBlock to be a multiple of the warp size
        that results in enough threads to have one per 2 columns. */
        threadsPerBlock = MAX_THREAD_PER_BLOCK(CUDA_WARP_SIZE*(SMCU_COLUMNS_S(A)/2 + CUDA_WARP_SIZE - 1)/CUDA_WARP_SIZE);
        gridSize = SMCU_ROWS_S(A);

      {
#ifdef SUNDIALS_CUDA_KERNEL_TIMING
        cudaEvent_t start, stop;
        float milliseconds = 0;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start);
#endif

        scaleAddIKernelCSR<realtype, int>
          <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_ROWS_S(A),
                                                     c,
                                                     SMCU_DATA_S(A),
                                                     SMCU_INDEXPTRS_S(A),
                                                     SMCU_INDEXVALS_S(A));

#ifdef SUNDIALS_CUDA_KERNEL_TIMING
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milliseconds, start, stop);
        fprintf(stdout, 
                "[performance] scaleAddIKernelCSR runtime (s): %22.15e\n",
                milliseconds/1000.0);
        /* scaleAddIKernelCSR reads 1 real, writes 1 real, reads 3 ints */
        fprintf(stdout,
                "[performance] scaleAddIKernelCSR effective bandwidth (GB/s): %f\n",
                (SMCU_NNZ_S(A)*(2*sizeof(realtype) + sizeof(int)) + 2*SMCU_ROWS_S(A)*sizeof(int))/milliseconds/1e6);
#endif
      }

      break;
    case SUNMAT_CUSPARSE_BCSR:
      /* Choose the grid size to be the number of blocks in the matrix,
         and then choose threadsPerBlock to be a multiple of the warp size
         that results in enough threads to have one per row of the block. */
      threadsPerBlock = MAX_THREAD_PER_BLOCK(CUDA_WARP_SIZE*(SMCU_BLOCKROWS_S(A) + CUDA_WARP_SIZE - 1)/CUDA_WARP_SIZE);
      gridSize = SMCU_NBLOCKS_S(A);

      {
#ifdef SUNDIALS_CUDA_KERNEL_TIMING
        cudaEvent_t start, stop;
        float milliseconds = 0;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start);
#endif

        scaleAddIKernelBCSR<realtype, int>
          <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_BLOCKROWS_S(A),
                                                     SMCU_NBLOCKS_S(A),
                                                     SMCU_BLOCKNNZ_S(A),
                                                     c,
                                                     SMCU_DATA_S(A),
                                                     SMCU_INDEXPTRS_S(A),
                                                     SMCU_INDEXVALS_S(A));

#ifdef SUNDIALS_CUDA_KERNEL_TIMING
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milliseconds, start, stop);
        fprintf(stdout, 
                "[performance] scaleAddIKernelBCSR runtime (s): %22.15e\n",
                milliseconds/1000.0);
        /* scaleAddIKernelBCSR reads 1 real, writes 1 real, reads 3 ints */
        fprintf(stdout,
                "[performance] scaleAddIKernelBCSR effective bandwidth (GB/s): %f\n",
                (SMCU_NNZ_S(A)*(2*sizeof(realtype) + sizeof(int)) + 2*SMCU_ROWS_S(A)*sizeof(int))/milliseconds/1e6);
#endif
      }
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMatScaleAddI_cuSparse: sparse type not recognized\n");
      return SUNMAT_ILL_INPUT;
  }

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return SUNMAT_OPERATION_FAIL;
#endif

  return SUNMAT_SUCCESS;
}


/* Performs A = cA + B */
int SUNMatScaleAdd_cuSparse(realtype c, SUNMatrix A, SUNMatrix B)
{
  cudaStream_t stream;
  cusparseStatus_t cusparse_status;

  if (!SMCompatible_cuSparse(A, B))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatScaleAdd_cuSparse: SUNMatScaleAdd_cuSparse failed\n");
    return SUNMAT_ILL_INPUT;
  }

  cusparse_status = cusparseGetStream(SMCU_CUSPHANDLE_S(A), &stream);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status)) return SUNMAT_OPERATION_FAIL;

  unsigned threadsPerBlock, gridSize;
  switch (SMCU_SPARSETYPE_S(A))
  {
    case SUNMAT_CUSPARSE_CSR:
      /* Choose the grid size to be the number of rows in the matrix,
        and then choose threadsPerBlock to be a multiple of the warp size
        that results in enough threads to have one per 2 columns. */
      threadsPerBlock = MAX_THREAD_PER_BLOCK(CUDA_WARP_SIZE*(SMCU_COLUMNS_S(A)/2 + CUDA_WARP_SIZE - 1)/CUDA_WARP_SIZE);
      gridSize = SMCU_ROWS_S(A);
     
      {
#ifdef SUNDIALS_CUDA_KERNEL_TIMING
        cudaEvent_t start, stop;
        float milliseconds = 0;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start);
#endif

        scaleAddKernelCSR<realtype, int>
          <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_NNZ_S(A),
                                                     c,
                                                     SMCU_DATA_S(A),
                                                     SMCU_DATA_S(B));

#ifdef SUNDIALS_CUDA_KERNEL_TIMING
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&milliseconds, start, stop);
        fprintf(stdout, 
                "[performance] scaleAddKernelCSR runtime (s): %22.15e\n",
                milliseconds/1000.0);
        /* scaleAddKernelCSR reads 2 realtype, and writes 1 realtype */
        fprintf(stdout,
                "[performance] scaleAddKernelCSR effective bandwidth (GB/s): %f\n",
                SMCU_NNZ_S(A)*sizeof(realtype)*3/milliseconds/1e6);
#endif
      }

      break;
    case SUNMAT_CUSPARSE_BCSR:
      /* Choose the grid size to be the number of blocks in the matrix,
         and then choose threadsPerBlock to be a multiple of the warp size
         that results in enough threads to have one per row of the block. */
      threadsPerBlock = MAX_THREAD_PER_BLOCK(CUDA_WARP_SIZE*(SMCU_BLOCKROWS_S(A) + CUDA_WARP_SIZE - 1)/CUDA_WARP_SIZE);
      gridSize = SMCU_NBLOCKS_S(A);

      {
#ifdef SUNDIALS_CUDA_KERNEL_TIMING
        cudaEvent_t start, stop;
        float milliseconds = 0;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start);
#endif

        scaleAddKernelCSR<realtype, int>
          <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_NNZ_S(A),
                                                     c,
                                                     SMCU_DATA_S(A),
                                                     SMCU_DATA_S(B));

#ifdef SUNDIALS_CUDA_KERNEL_TIMING
          cudaEventRecord(stop);
          cudaEventSynchronize(stop);
          cudaEventElapsedTime(&milliseconds, start, stop);
          fprintf(stdout, 
                  "[performance] scaleAddKernelCSR (BCSR format) runtime (s): %22.15e\n",
                  milliseconds/1000.0);
          /* scaleAddKernelCSR reads 2 realtype, and writes 1 realtype */
          fprintf(stdout,
                  "[performance] scaleAddKernelCSR (BCSR format) effective bandwidth (GB/s): %f\n",
                  SMCU_NNZ_S(A)*sizeof(realtype)*3/milliseconds/1e6);
#endif
      }

      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMatScaleAdd_cuSparse: sparse type not recognized\n");
      return SUNMAT_ILL_INPUT;
  }

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return SUNMAT_OPERATION_FAIL;
#endif

  return SUNMAT_SUCCESS;
}


/* Perform y = Ax */
int SUNMatMatvec_cuSparse(SUNMatrix A, N_Vector x, N_Vector y)
{
  /* Verify that the dimensions of A, x, and y agree */
  if ( (SMCU_COLUMNS_S(A) != N_VGetLength(x)) ||
       (SMCU_ROWS_S(A) != N_VGetLength(y)) )
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatMatvec_cuSparse: dimensions do not agree\n");
    return SUNMAT_ILL_INPUT;
  }

  realtype *d_xdata = N_VGetDeviceArrayPointer_Cuda(x);
  realtype *d_ydata = N_VGetDeviceArrayPointer_Cuda(y);

  if (SMCU_SPARSETYPE_S(A) == SUNMAT_CUSPARSE_CSR)
  {
    const realtype one = ONE;
    cusparseStatus_t cusparse_status;

    /* Zero result vector */
    N_VConst(ZERO, y);

    {
#ifdef SUNDIALS_CUDA_KERNEL_TIMING
      cudaEvent_t start, stop;
      float milliseconds = 0;
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      cudaEventRecord(start);
#endif

      cusparse_status = cusparseXcsrmv(SMCU_CUSPHANDLE_S(A),
                                       CUSPARSE_OPERATION_NON_TRANSPOSE,
                                       SMCU_ROWS_S(A),
                                       SMCU_COLUMNS_S(A),
                                       SMCU_NNZ_S(A),
                                       &one,
                                       SMCU_MATDESCR_S(A),
                                       SMCU_DATA_S(A),
                                       SMCU_INDEXPTRS_S(A),
                                       SMCU_INDEXVALS_S(A),
                                       d_xdata,
                                       &one,
                                       d_ydata);

#ifdef SUNDIALS_CUDA_KERNEL_TIMING
          cudaEventRecord(stop);
          cudaEventSynchronize(stop);
          cudaEventElapsedTime(&milliseconds, start, stop);
          fprintf(stdout, 
                  "[performance] cusparseXcsrmv untime (s): %22.15e\n",
                  milliseconds/1000.0);
#endif
    }

    if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status)) return SUNMAT_OPERATION_FAIL;
  }
  else if (SMCU_SPARSETYPE_S(A) == SUNMAT_CUSPARSE_BCSR)
  {
    cudaStream_t stream;
    cusparseStatus_t cusparse_status;
    unsigned gridSize, threadsPerBlock;

    cusparse_status = cusparseGetStream(SMCU_CUSPHANDLE_S(A), &stream);
    if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status)) return SUNMAT_OPERATION_FAIL;

    /* Choose the grid size to be the number of blocks in the matrix,
       and then choose threadsPerBlock to be a multiple of the warp size
       that results in enough threads to have one per row of the block. */
    threadsPerBlock = MAX_THREAD_PER_BLOCK(CUDA_WARP_SIZE*(SMCU_BLOCKROWS_S(A) + CUDA_WARP_SIZE - 1)/CUDA_WARP_SIZE);
    gridSize = SMCU_NBLOCKS_S(A);

    {
#ifdef SUNDIALS_CUDA_KERNEL_TIMING
        cudaEvent_t start, stop;
        float milliseconds = 0;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start);
#endif

      matvecBCSR<realtype, int>
        <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_BLOCKROWS_S(A),
                                                   SMCU_NBLOCKS_S(A),
                                                   SMCU_BLOCKNNZ_S(A),
                                                   SMCU_DATA_S(A),
                                                   SMCU_INDEXPTRS_S(A),
                                                   SMCU_INDEXVALS_S(A),
                                                   d_xdata,
                                                   d_ydata);

#ifdef SUNDIALS_CUDA_KERNEL_TIMING
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds, start, stop);
      fprintf(stdout, 
              "[performance] matvecBCSR runtime (s): %22.15e\n",
              milliseconds/1000.0);
      fprintf(stdout,
              "[performance] matvecBCSR effective bandwidth (GB/s): %f\n",
              (SMCU_NNZ_S(A)*(sizeof(realtype)*4 + sizeof(int)) + 2*SMCU_ROWS_S(A)*sizeof(int))/milliseconds/1e6);
#endif

    }

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
    cudaDeviceSynchronize();
    if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return SUNMAT_OPERATION_FAIL;
#endif
  }
  else
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatMatvec_cuSparse: sparse type not recognized\n");
    return SUNMAT_ILL_INPUT;
  }

  return SUNMAT_SUCCESS;
}


/*
 * =================================================================
 * private functions
 * =================================================================
 */


/* -----------------------------------------------------------------
 * Function to check compatibility of two sparse SUNMatrix objects
 */
static booleantype SMCompatible_cuSparse(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must be sparse */
  if ( (SUNMatGetID(A) != SUNMATRIX_CUSPARSE) ||
       (SUNMatGetID(B) != SUNMATRIX_CUSPARSE) )
    return SUNFALSE;

  /* both matrices must have the same shape and sparsity type */
  if (SMCU_ROWS_S(A) != SMCU_ROWS_S(B))
    return SUNFALSE;
  if (SMCU_COLUMNS_S(A) != SMCU_COLUMNS_S(B))
    return SUNFALSE;
  if (SMCU_SPARSETYPE_S(A) != SMCU_SPARSETYPE_S(B))
    return SUNFALSE;

  return SUNTRUE;
}

/* -----------------------------------------------------------------
 * Function to create empty SUNMatrix with ops attached and
 * the content structure allocated.
 */
SUNMatrix SUNMatrix_cuSparse_NewEmpty()
{
  /* Create an empty matrix object */
  SUNMatrix A = NULL;
  A = SUNMatNewEmpty();
  if (A == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_NewEmpty: SUNMatNewEmpty failed\n");
    return NULL;
  }

  /* Attach operations */
  A->ops->getid     = SUNMatGetID_cuSparse;
  A->ops->clone     = SUNMatClone_cuSparse;
  A->ops->destroy   = SUNMatDestroy_cuSparse;
  A->ops->zero      = SUNMatZero_cuSparse;
  A->ops->copy      = SUNMatCopy_cuSparse;
  A->ops->scaleadd  = SUNMatScaleAdd_cuSparse;
  A->ops->scaleaddi = SUNMatScaleAddI_cuSparse;
  A->ops->matvec    = SUNMatMatvec_cuSparse;

  /* Create content */
  SUNMatrix_Content_cuSparse content = NULL;
  content = (SUNMatrix_Content_cuSparse) malloc(sizeof *content);
  if (content == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_NewEmpty: failed to malloc content\n");
    SUNMatDestroy(A);
    return NULL;
  }

  /* Attach content */
  A->content = content;

  return A;
}
