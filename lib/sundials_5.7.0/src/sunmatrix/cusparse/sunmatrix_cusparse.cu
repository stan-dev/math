/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
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

#include <sunmemory/sunmemory_cuda.h>
#include <sunmatrix/sunmatrix_cusparse.h>

#include "sundials_cuda.h"
#include "sundials_debug.h"
#include "cusparse_kernels.cuh"


/* Use the namespace for the kernels */
using namespace sundials::sunmatrix_cusparse;

/* Constants */
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Private function prototypes */
static booleantype SMCompatible_cuSparse(SUNMatrix, SUNMatrix);
static SUNMatrix SUNMatrix_cuSparse_NewEmpty();
#if CUDART_VERSION >= 11000
static cusparseStatus_t CreateSpMatDescr(SUNMatrix, cusparseSpMatDescr_t*);
#endif

/* Macros for handling the different function names based on precision */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define cusparseXcsrmv cusparseDcsrmv
#define CUDA_R_XF CUDA_R_64F
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define cusparseXcsrmv cusparseScsrmv
#define CUDA_R_XF CUDA_R_32F
#endif

/* Content accessor macros */
#define SMCU_CONTENT(A)     ( (SUNMatrix_Content_cuSparse)(A->content) )
#define SMCU_ROWS(A)        ( SMCU_CONTENT(A)->M )
#define SMCU_COLUMNS(A)     ( SMCU_CONTENT(A)->N )
#define SMCU_NNZ(A)         ( SMCU_CONTENT(A)->NNZ )
#define SMCU_NBLOCKS(A)     ( SMCU_CONTENT(A)->nblocks )
#define SMCU_BLOCKROWS(A)   ( SMCU_CONTENT(A)->blockrows )
#define SMCU_BLOCKCOLS(A)   ( SMCU_CONTENT(A)->blockcols )
#define SMCU_BLOCKNNZ(A)    ( SMCU_CONTENT(A)->blocknnz )
#define SMCU_NP(A)          ( SMCU_CONTENT(A)->NP )
#define SMCU_SPARSETYPE(A)  ( SMCU_CONTENT(A)->sparse_type )
#define SMCU_OWNMATD(A)     ( SMCU_CONTENT(A)->own_matd )
#define SMCU_OWNEXEC(A)     ( SMCU_CONTENT(A)->own_exec )
#define SMCU_DATA(A)        ( SMCU_CONTENT(A)->data )
#define SMCU_DATAp(A)       ( (realtype*)SMCU_CONTENT(A)->data->ptr )
#define SMCU_INDEXVALS(A)   ( SMCU_CONTENT(A)->colind )
#define SMCU_INDEXPTRS(A)   ( SMCU_CONTENT(A)->rowptrs )
#define SMCU_INDEXVALSp(A)  ( (int*) SMCU_CONTENT(A)->colind->ptr )
#define SMCU_INDEXPTRSp(A)  ( (int*) SMCU_CONTENT(A)->rowptrs->ptr )
#define SMCU_MEMHELP(A)     ( SMCU_CONTENT(A)->mem_helper )
#define SMCU_MATDESCR(A)    ( SMCU_CONTENT(A)->mat_descr )
#define SMCU_CUSPHANDLE(A)  ( SMCU_CONTENT(A)->cusp_handle )
#define SMCU_FIXEDPATTERN(A)( SMCU_CONTENT(A)->fixed_pattern )
#define SMCU_EXECPOLICY(A)  ( SMCU_CONTENT(A)->exec_policy )


/* ------------------------------------------------------------------
 * Default execution policy definition.
 *
 * This policy tries to help us leverage the structure of the matrix.
 * It will choose block sizes which are a multiple of the warp size,
 * and it will choose a grid size to such that all work elements are
 * covered.
 * ------------------------------------------------------------------ */

class SUNCuSparseMatrixExecPolicy : public SUNCudaExecPolicy
{
public:
  SUNCuSparseMatrixExecPolicy(const cudaStream_t stream = 0)
    : stream_(stream)
  {}

  SUNCuSparseMatrixExecPolicy(const SUNCuSparseMatrixExecPolicy& ex)
    : stream_(ex.stream_)
  {}

  virtual size_t gridSize(size_t numWorkElements, size_t blockDim = 0) const
  {
    return(numWorkElements + blockDim - 1)/blockDim;
  }

  virtual size_t blockSize(size_t numWorkElements = 0, size_t gridDim = 0) const
  {
    return(max_block_size(CUDA_WARP_SIZE*(numWorkElements + CUDA_WARP_SIZE - 1)/CUDA_WARP_SIZE));
  }

  virtual const cudaStream_t* stream() const
  {
    return(&stream_);
  }

  virtual CudaExecPolicy* clone() const
  {
    return(static_cast<CudaExecPolicy*>(new SUNCuSparseMatrixExecPolicy(*this)));
  }

  static size_t max_block_size(int val)
  {
    return((val > MAX_CUDA_BLOCKSIZE) ? MAX_CUDA_BLOCKSIZE : val );
  }

private:
  const cudaStream_t stream_;
};

/* ------------------------------------------------------------------
 * Constructors.
 * ------------------------------------------------------------------ */

SUNMatrix SUNMatrix_cuSparse_NewCSR(int M, int N, int NNZ, cusparseHandle_t cusp)
{
  SUNMemory d_colind, d_rowptr, d_values;
  int alloc_fail = 0;

  /* return with NULL matrix on illegal input */
  if ( (M <= 0) || (N <= 0) || (NNZ < 0) )
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_NewCSR_cuSparse: illegal value(s) for M, N, or NNZ\n");
    return(NULL);
  }

  SUNMatrix A = SUNMatrix_cuSparse_NewEmpty();
  if (A == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_NewCSR_cuSparse: SUNMatrix_cuSparse_NewEmpty returned NULL\n");
    return(NULL);
  }

  SMCU_MEMHELP(A) = SUNMemoryHelper_Cuda();
  if (SMCU_MEMHELP(A) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_NewCSR_cuSparse: SUNMemoryHelper_Cuda returned NULL\n");
    SUNMatDestroy(A);
    return(NULL);
  }

  /* Allocate device memory for the matrix */
  alloc_fail += SUNMemoryHelper_Alloc(SMCU_MEMHELP(A), &d_colind,
                                      sizeof(int)*NNZ, SUNMEMTYPE_DEVICE);
  alloc_fail += SUNMemoryHelper_Alloc(SMCU_MEMHELP(A), &d_rowptr,
                                      sizeof(int)*(M+1), SUNMEMTYPE_DEVICE);
  alloc_fail += SUNMemoryHelper_Alloc(SMCU_MEMHELP(A), &d_values,
                                      sizeof(realtype)*NNZ, SUNMEMTYPE_DEVICE);
  if (alloc_fail)
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    SUNMatDestroy(A);
    return(NULL);
  }

  /* Choose sensible defaults */
  cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
  cusparseMatDescr_t mat_descr;
  cusparse_status = cusparseCreateMatDescr(&mat_descr);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    SUNMatDestroy(A);
    return(NULL);
  }

  cusparse_status = cusparseSetMatType(mat_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    cusparseDestroyMatDescr(mat_descr);
    SUNMatDestroy(A);
    return(NULL);
  }

  cusparse_status = cusparseSetMatIndexBase(mat_descr, CUSPARSE_INDEX_BASE_ZERO);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    cusparseDestroyMatDescr(mat_descr);
    SUNMatDestroy(A);
    return(NULL);
  }

  cudaStream_t stream;
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparseGetStream(cusp, &stream)))
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    cusparseDestroyMatDescr(mat_descr);
    SUNMatDestroy(A);
    return(NULL);
  }

  /* Fill the content */
  SMCU_CONTENT(A)->M              = M;
  SMCU_CONTENT(A)->N              = N;
  SMCU_CONTENT(A)->NNZ            = NNZ;
  SMCU_CONTENT(A)->nblocks        = 1;
  SMCU_CONTENT(A)->blockrows      = M;
  SMCU_CONTENT(A)->blockcols      = N;
  SMCU_CONTENT(A)->blocknnz       = NNZ;
  SMCU_CONTENT(A)->own_matd       = SUNTRUE;
  SMCU_CONTENT(A)->own_exec       = SUNTRUE;
  SMCU_CONTENT(A)->matvec_issetup = SUNFALSE;
  SMCU_CONTENT(A)->fixed_pattern  = SUNFALSE;
  SMCU_CONTENT(A)->sparse_type    = SUNMAT_CUSPARSE_CSR;
  SMCU_CONTENT(A)->colind         = d_colind;
  SMCU_CONTENT(A)->rowptrs        = d_rowptr;
  SMCU_CONTENT(A)->data           = d_values;
  SMCU_CONTENT(A)->mat_descr      = mat_descr;
  SMCU_CONTENT(A)->cusp_handle    = cusp;
  SMCU_CONTENT(A)->exec_policy    = new SUNCuSparseMatrixExecPolicy(stream);

#if CUDART_VERSION >= 11000
  cusparseSpMatDescr_t spmat_descr;
  if (!SUNDIALS_CUSPARSE_VERIFY(CreateSpMatDescr(A, &spmat_descr)))
  {
    SUNMatDestroy(A);
    return(NULL);
  }
  SMCU_CONTENT(A)->spmat_descr = spmat_descr;
  SMCU_CONTENT(A)->dBufferMem  = NULL;
  SMCU_CONTENT(A)->bufferSize  = 0;
  SMCU_CONTENT(A)->vecX        = NULL;
  SMCU_CONTENT(A)->vecY        = NULL;
#endif

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
    return(NULL);
  }

  if ( (rowptrs == NULL) || (colind == NULL) || (data == NULL) )
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_MakeCSR_cuSparse: rowptrs, colind, or data is NULL\n");
    return(NULL);
  }

  if (cusparseGetMatIndexBase(mat_descr) != CUSPARSE_INDEX_BASE_ZERO)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_MakeCSR_cuSparse: the cusparseMatDescr_t must have index base CUSPARSE_INDEX_BASE_ZERO\n");
    return(NULL);
  }

  SUNMatrix A = SUNMatrix_cuSparse_NewEmpty();
  if (A == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_MakeCSR_cuSparse: SUNMatrix_cuSparse_NewEmpty returned NULL\n");
    return(NULL);
  }

  SMCU_MEMHELP(A) = SUNMemoryHelper_Cuda();
  if (SMCU_MEMHELP(A) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_NewCSR_cuSparse: SUNMemoryHelper_Cuda returned NULL\n");
    SUNMatDestroy(A);
    return(NULL);
  }

  cudaStream_t stream;
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparseGetStream(cusp, &stream)))
  {
    SUNMatDestroy(A);
    return(NULL);
  }

  /* Fill content */
  SMCU_CONTENT(A)->M              = M;
  SMCU_CONTENT(A)->N              = N;
  SMCU_CONTENT(A)->NNZ            = NNZ;
  SMCU_CONTENT(A)->nblocks        = 1;
  SMCU_CONTENT(A)->blockrows      = M;
  SMCU_CONTENT(A)->blockcols      = N;
  SMCU_CONTENT(A)->blocknnz       = NNZ;
  SMCU_CONTENT(A)->own_matd       = SUNFALSE;
  SMCU_CONTENT(A)->own_exec       = SUNTRUE;
  SMCU_CONTENT(A)->matvec_issetup = SUNFALSE;
  SMCU_CONTENT(A)->fixed_pattern  = SUNFALSE;
  SMCU_CONTENT(A)->sparse_type    = SUNMAT_CUSPARSE_CSR;
  SMCU_CONTENT(A)->colind         = SUNMemoryHelper_Wrap(colind, SUNMEMTYPE_DEVICE);
  SMCU_CONTENT(A)->rowptrs        = SUNMemoryHelper_Wrap(rowptrs, SUNMEMTYPE_DEVICE);
  SMCU_CONTENT(A)->data           = SUNMemoryHelper_Wrap(data, SUNMEMTYPE_DEVICE);
  SMCU_CONTENT(A)->mat_descr      = mat_descr;
  SMCU_CONTENT(A)->cusp_handle    = cusp;

  SMCU_CONTENT(A)->exec_policy   = new SUNCuSparseMatrixExecPolicy(stream);

  if (SMCU_CONTENT(A)->colind == NULL ||
      SMCU_CONTENT(A)->rowptrs == NULL ||
      SMCU_CONTENT(A)->data == NULL)
  {
    SUNMatDestroy(A);
    return(NULL);
  }

#if CUDART_VERSION >= 11000
  cusparseSpMatDescr_t spmat_descr;
  if (!SUNDIALS_CUSPARSE_VERIFY(CreateSpMatDescr(A, &spmat_descr)))
  {
    SUNMatDestroy(A);
    return(NULL);
  }
  SMCU_CONTENT(A)->spmat_descr = spmat_descr;
  SMCU_CONTENT(A)->dBufferMem  = NULL;
  SMCU_CONTENT(A)->bufferSize  = 0;
  SMCU_CONTENT(A)->vecX        = NULL;
  SMCU_CONTENT(A)->vecY        = NULL;
#endif

  return(A);
}


SUNMatrix SUNMatrix_cuSparse_NewBlockCSR(int nblocks, int blockrows, int blockcols, int blocknnz, cusparseHandle_t cusp)
{
  SUNMemory d_colind, d_rowptr, d_values;
  int M, N, NNZ;
  int alloc_fail = 0;

  /* Return with NULL matrix on illegal input */
  if (blockrows != blockcols)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_NewBlockCSR: matrix must be square for the BCSR format\n");
    return(NULL);
  }

  M   = nblocks * blockrows;
  N   = M;
  NNZ = nblocks * blocknnz;

  /* Return with NULL matrix on illegal input */
  if ( (M <= 0) || (N <= 0) || (NNZ < 0) )
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_NewBlockCSR: illegal value(s) for M, N, or NNZ\n");
    return(NULL);
  }

  /* Allocate the SUNMatrix object */
  SUNMatrix A = SUNMatrix_cuSparse_NewEmpty();
  if (A == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_NewBlockCSR: SUNMatrix_cuSparse_NewEmpty returned NULL\n");
    return(NULL);
  }

  SMCU_MEMHELP(A) = SUNMemoryHelper_Cuda();
  if (SMCU_MEMHELP(A) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_NewCSR_cuSparse: SUNMemoryHelper_Cuda returned NULL\n");
    SUNMatDestroy(A);
    return(NULL);
  }

  /* Allocate device memory for the matrix */
  alloc_fail += SUNMemoryHelper_Alloc(SMCU_MEMHELP(A), &d_colind,
                                      sizeof(int)*blocknnz, SUNMEMTYPE_DEVICE);
  alloc_fail += SUNMemoryHelper_Alloc(SMCU_MEMHELP(A), &d_rowptr,
                                      sizeof(int)*(blockrows + 1),
                                      SUNMEMTYPE_DEVICE);
  alloc_fail += SUNMemoryHelper_Alloc(SMCU_MEMHELP(A), &d_values,
                                      sizeof(realtype)*blocknnz*nblocks,
                                      SUNMEMTYPE_DEVICE);
  if (alloc_fail)
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    SUNMatDestroy(A);
    return(NULL);
  }

  /* Choose sensible defaults */
  cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
  cusparseMatDescr_t mat_descr;
  cusparse_status = cusparseCreateMatDescr(&mat_descr);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    SUNMatDestroy(A);
    return(NULL);
  }

  cusparse_status = cusparseSetMatType(mat_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    cusparseDestroyMatDescr(mat_descr);
    SUNMatDestroy(A);
    return(NULL);
  }

  cusparse_status = cusparseSetMatIndexBase(mat_descr, CUSPARSE_INDEX_BASE_ZERO);
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparse_status))
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    cusparseDestroyMatDescr(mat_descr);
    SUNMatDestroy(A);
    return(NULL);
  }

  cudaStream_t stream;
  if (!SUNDIALS_CUSPARSE_VERIFY(cusparseGetStream(cusp, &stream)))
  {
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_colind);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_rowptr);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), d_values);
    cusparseDestroyMatDescr(mat_descr);
    SUNMatDestroy(A);
    return(NULL);
  }

  /* Fill the content */
  SMCU_CONTENT(A)->M              = M;
  SMCU_CONTENT(A)->N              = N;
  SMCU_CONTENT(A)->NNZ            = NNZ;
  SMCU_CONTENT(A)->nblocks        = nblocks;
  SMCU_CONTENT(A)->blockrows      = blockrows;
  SMCU_CONTENT(A)->blockcols      = blockrows;
  SMCU_CONTENT(A)->blocknnz       = blocknnz;
  SMCU_CONTENT(A)->own_matd       = SUNTRUE;
  SMCU_CONTENT(A)->own_exec       = SUNTRUE;
  SMCU_CONTENT(A)->matvec_issetup = SUNFALSE;
  SMCU_CONTENT(A)->cusp_handle    = cusp;
  SMCU_CONTENT(A)->fixed_pattern  = SUNFALSE;
  SMCU_CONTENT(A)->sparse_type    = SUNMAT_CUSPARSE_BCSR;
  SMCU_CONTENT(A)->colind         = d_colind;
  SMCU_CONTENT(A)->rowptrs        = d_rowptr;
  SMCU_CONTENT(A)->data           = d_values;
  SMCU_CONTENT(A)->mat_descr      = mat_descr;
  SMCU_CONTENT(A)->exec_policy    = new SUNCuSparseMatrixExecPolicy(stream);

#if CUDART_VERSION >= 11000
  cusparseSpMatDescr_t spmat_descr;
  if (!SUNDIALS_CUSPARSE_VERIFY(CreateSpMatDescr(A, &spmat_descr)))
  {
    SUNMatDestroy(A);
    return(NULL);
  }
  SMCU_CONTENT(A)->spmat_descr = spmat_descr;
  SMCU_CONTENT(A)->dBufferMem  = NULL;
  SMCU_CONTENT(A)->bufferSize  = 0;
  SMCU_CONTENT(A)->vecX        = NULL;
  SMCU_CONTENT(A)->vecY        = NULL;
#endif

  return(A);
}

/* ------------------------------------------------------------------
 * Implementation specific routines.
 * ------------------------------------------------------------------ */

int SUNMatrix_cuSparse_SparseType(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_SPARSETYPE(A));
  else
    return(SUNMAT_ILL_INPUT);
}

int SUNMatrix_cuSparse_Rows(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_ROWS(A));
  else
    return(SUNMAT_ILL_INPUT);
}

int SUNMatrix_cuSparse_Columns(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_COLUMNS(A));
  else
    return(SUNMAT_ILL_INPUT);
}

int SUNMatrix_cuSparse_NNZ(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_NNZ(A));
  else
    return(SUNMAT_ILL_INPUT);
}

int* SUNMatrix_cuSparse_IndexPointers(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_INDEXPTRSp(A));
  else
    return(NULL);
}

int* SUNMatrix_cuSparse_IndexValues(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_INDEXVALSp(A));
  else
    return(NULL);
}

realtype* SUNMatrix_cuSparse_Data(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_DATAp(A));
  else
    return(NULL);
}

int SUNMatrix_cuSparse_NumBlocks(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_NBLOCKS(A));
  else
    return(SUNMAT_ILL_INPUT);
}

int SUNMatrix_cuSparse_BlockRows(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_BLOCKROWS(A));
  else
    return(SUNMAT_ILL_INPUT);
}

int SUNMatrix_cuSparse_BlockColumns(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_BLOCKCOLS(A));
  else
    return(SUNMAT_ILL_INPUT);
}

int SUNMatrix_cuSparse_BlockNNZ(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_BLOCKNNZ(A));
  else
    return(SUNMAT_ILL_INPUT);
}

realtype* SUNMatrix_cuSparse_BlockData(SUNMatrix A, int blockidx)
{
  realtype *matdata;
  int offset;

  if (SUNMatGetID(A) != SUNMATRIX_CUSPARSE)
    return(NULL);

  if (blockidx >= SMCU_NBLOCKS(A))
    return(NULL);

  matdata = SMCU_DATAp(A);
  offset = SMCU_BLOCKNNZ(A)*blockidx;

  return(&matdata[offset]);
}

cusparseMatDescr_t SUNMatrix_cuSparse_MatDescr(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_CUSPARSE)
    return(SMCU_MATDESCR(A));
  else
    return(NULL);
}

int SUNMatrix_cuSparse_SetFixedPattern(SUNMatrix A, booleantype yesno)
{
  if (SUNMatGetID(A) != SUNMATRIX_CUSPARSE)
    return(SUNMAT_ILL_INPUT);

  SMCU_FIXEDPATTERN(A) = yesno;

  return(SUNMAT_SUCCESS);
}


int SUNMatrix_cuSparse_SetKernelExecPolicy(SUNMatrix A, SUNCudaExecPolicy* exec_policy)
{
  if (SUNMatGetID(A) != SUNMATRIX_CUSPARSE || exec_policy == NULL)
    return(SUNMAT_ILL_INPUT);

  if (SMCU_OWNEXEC(A)) delete SMCU_EXECPOLICY(A);
  SMCU_EXECPOLICY(A) = exec_policy;

  SMCU_OWNEXEC(A) = SUNFALSE;

  return(SUNMAT_SUCCESS);
}


int SUNMatrix_cuSparse_CopyToDevice(SUNMatrix dA, realtype* h_data,
                                    int* h_idxptrs, int* h_idxvals)
{
  int retval;
  SUNMemory _h_data, _h_idxptrs, _h_idxvals;
  const cudaStream_t* stream;
  int nidxvals, nidxptrs;

  if (SUNMatGetID(dA) != SUNMATRIX_CUSPARSE)
    return(SUNMAT_ILL_INPUT);

  stream  = SMCU_EXECPOLICY(dA)->stream();

  if (h_data != NULL)
  {
    _h_data = SUNMemoryHelper_Wrap(h_data, SUNMEMTYPE_HOST);
    retval  = SUNMemoryHelper_CopyAsync(SMCU_MEMHELP(dA),
                                        SMCU_DATA(dA),
                                        _h_data,
                                        SMCU_NNZ(dA)*sizeof(realtype),
                                        (void*) stream);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(dA), _h_data);
    if (retval != 0) return(SUNMAT_OPERATION_FAIL);
  }

  switch(SMCU_SPARSETYPE(dA))
  {
    case SUNMAT_CUSPARSE_CSR:
      nidxptrs = SMCU_ROWS(dA)+1;
      nidxvals = SMCU_NNZ(dA);
      break;
    case SUNMAT_CUSPARSE_BCSR:
      nidxptrs = SMCU_BLOCKROWS(dA)+1;
      nidxvals = SMCU_BLOCKNNZ(dA);
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_CopyToDevice: unrecognized sparse type\n");
      return(SUNMAT_ILL_INPUT);
  }

  if (h_idxptrs != NULL)
  {
    _h_idxptrs = SUNMemoryHelper_Wrap(h_idxptrs, SUNMEMTYPE_HOST);
    retval = SUNMemoryHelper_CopyAsync(SMCU_MEMHELP(dA),
                                       SMCU_INDEXPTRS(dA),
                                       _h_idxptrs,
                                       nidxptrs*sizeof(int),
                                       (void*) stream);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(dA), _h_idxptrs);
    if (retval != 0) return(SUNMAT_OPERATION_FAIL);
  }

  if (h_idxvals != NULL)
  {
    _h_idxvals = SUNMemoryHelper_Wrap(h_idxvals, SUNMEMTYPE_HOST);
    retval = SUNMemoryHelper_CopyAsync(SMCU_MEMHELP(dA),
                                       SMCU_INDEXVALS(dA),
                                       _h_idxvals,
                                       nidxvals*sizeof(int),
                                       (void*) stream);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(dA), _h_idxvals);
    if (retval != 0) return(SUNMAT_OPERATION_FAIL);
  }

  return(SUNMAT_SUCCESS);
}


int SUNMatrix_cuSparse_CopyFromDevice(SUNMatrix dA, realtype* h_data,
                                      int* h_idxptrs, int* h_idxvals)
{
  int retval;
  SUNMemory _h_data, _h_idxptrs, _h_idxvals;
  const cudaStream_t* stream;
  int nidxvals, nidxptrs;

  if (SUNMatGetID(dA) != SUNMATRIX_CUSPARSE)
    return(SUNMAT_ILL_INPUT);

  stream = SMCU_EXECPOLICY(dA)->stream();

  if (h_data != NULL)
  {
    _h_data = SUNMemoryHelper_Wrap(h_data, SUNMEMTYPE_HOST);
    retval  = SUNMemoryHelper_CopyAsync(SMCU_MEMHELP(dA),
                                        _h_data,
                                        SMCU_DATA(dA),
                                        SMCU_NNZ(dA)*sizeof(realtype),
                                        (void*) stream);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(dA), _h_data);
    if (retval != 0) return(SUNMAT_OPERATION_FAIL);
  }


  switch(SMCU_SPARSETYPE(dA))
  {
    case SUNMAT_CUSPARSE_CSR:
      nidxptrs = SMCU_ROWS(dA)+1;
      nidxvals = SMCU_NNZ(dA);
    case SUNMAT_CUSPARSE_BCSR:
      nidxptrs = SMCU_BLOCKROWS(dA)+1;
      nidxvals = SMCU_BLOCKNNZ(dA);
  }

  if (h_idxptrs != NULL)
  {
    _h_idxptrs = SUNMemoryHelper_Wrap(h_idxptrs, SUNMEMTYPE_HOST);
    retval = SUNMemoryHelper_CopyAsync(SMCU_MEMHELP(dA),
                                       _h_idxptrs,
                                       SMCU_INDEXPTRS(dA),
                                       nidxptrs*sizeof(int),
                                       (void*) stream);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(dA), _h_idxptrs);
    if (retval != 0) return(SUNMAT_OPERATION_FAIL);
  }

  if (h_idxvals != NULL)
  {
    _h_idxvals = SUNMemoryHelper_Wrap(h_idxvals, SUNMEMTYPE_HOST);
    retval = SUNMemoryHelper_CopyAsync(SMCU_MEMHELP(dA),
                                       _h_idxvals,
                                       SMCU_INDEXVALS(dA),
                                       nidxvals*sizeof(int),
                                       (void*) stream);
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(dA), _h_idxvals);
    if (retval != 0) return(SUNMAT_OPERATION_FAIL);
  }


  return(SUNMAT_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */


SUNMatrix_ID SUNMatGetID_cuSparse(SUNMatrix A)
{
  return(SUNMATRIX_CUSPARSE);
}

/* Returns a new matrix allocated to have the same structure as A,
   but it does not copy any nonzeros, column vals, or row pointers. */
SUNMatrix SUNMatClone_cuSparse(SUNMatrix A)
{
  SUNMatrix B;

  switch (SMCU_SPARSETYPE(A))
  {
    case SUNMAT_CUSPARSE_CSR:
      B = SUNMatrix_cuSparse_NewCSR(SMCU_ROWS(A), SMCU_COLUMNS(A), SMCU_NNZ(A),
                                    SMCU_CUSPHANDLE(A));
      break;
    case SUNMAT_CUSPARSE_BCSR:
      B = SUNMatrix_cuSparse_NewBlockCSR(SMCU_NBLOCKS(A), SMCU_BLOCKROWS(A), SMCU_BLOCKCOLS(A),
                                         SMCU_BLOCKNNZ(A), SMCU_CUSPHANDLE(A));
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMatClone_cuSparse: sparse type not recognized\n");
      B = NULL;
  }

  SMCU_FIXEDPATTERN(B) = SMCU_FIXEDPATTERN(A);
  delete SMCU_EXECPOLICY(B);
  SMCU_EXECPOLICY(B) = SMCU_EXECPOLICY(A)->clone();

  return(B);
}


/* Deallocates the SUNMatrix object and all data it owns */
void SUNMatDestroy_cuSparse(SUNMatrix A)
{
  if (A == NULL) return;

  /* free content */
  if (A->content != NULL)
  {
    if (SMCU_MEMHELP(A))
    {
      SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), SMCU_DATA(A));
      SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), SMCU_INDEXPTRS(A));
      SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), SMCU_INDEXVALS(A));
    }
    else
    {
      SUNDIALS_DEBUG_PRINT("WARNING in SUNMatDestroy_cuSparse: mem_helper was NULL when trying to dealloc data, this could result in a memory leak\n");
    }

    if (SMCU_OWNMATD(A))
    {
      /* free cusparseMatDescr_t */
      SUNDIALS_CUSPARSE_VERIFY( cusparseDestroyMatDescr(SMCU_MATDESCR(A)) );
    }

#if CUDART_VERSION >= 11000
    SUNDIALS_CUSPARSE_VERIFY( cusparseDestroyDnVec(SMCU_CONTENT(A)->vecX) );
    SUNDIALS_CUSPARSE_VERIFY( cusparseDestroyDnVec(SMCU_CONTENT(A)->vecY) );
    SUNDIALS_CUSPARSE_VERIFY( cusparseDestroySpMat(SMCU_CONTENT(A)->spmat_descr) );
    SUNMemoryHelper_Dealloc(SMCU_MEMHELP(A), SMCU_CONTENT(A)->dBufferMem);
#endif

    if (SMCU_EXECPOLICY(A) && SMCU_OWNEXEC(A))
    {
      delete SMCU_EXECPOLICY(A);
      SMCU_EXECPOLICY(A) = NULL;
    }

    SUNMemoryHelper_Destroy(SMCU_MEMHELP(A));

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

  stream = *SMCU_EXECPOLICY(A)->stream();

  /* set all data to zero */
  cuerr = cudaMemsetAsync(SMCU_DATAp(A), 0, SMCU_NNZ(A)*sizeof(realtype), stream);
  if (!SUNDIALS_CUDA_VERIFY(cuerr)) return(SUNMAT_OPERATION_FAIL);

  /* set all rowptrs to zero unless the sparsity pattern is fixed */
  if (!SMCU_FIXEDPATTERN(A))
  {
    cuerr = cudaMemsetAsync(SMCU_INDEXPTRSp(A), 0,
                            (SMCU_BLOCKROWS(A)+1)*sizeof(int),
                            stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return(SUNMAT_OPERATION_FAIL);

    /* set all colind to zero */
    cuerr = cudaMemsetAsync(SMCU_INDEXVALSp(A), 0,
                            SMCU_BLOCKNNZ(A)*sizeof(int),
                            stream);
    if (!SUNDIALS_CUDA_VERIFY(cuerr)) return(SUNMAT_OPERATION_FAIL);
  }

  return(SUNMAT_SUCCESS);
}


/* Copies the nonzeros, column vals, and row pointers into dst */
int SUNMatCopy_cuSparse(SUNMatrix src, SUNMatrix dst)
{
  int retval;
  const cudaStream_t* stream;

  /* Verify that src and dst are compatible */
  if (!SMCompatible_cuSparse(src, dst))
    return(SUNMAT_ILL_INPUT);

  stream = SMCU_EXECPOLICY(src)->stream();

  /* Ensure that dst is allocated with at least as
     much memory as we have nonzeros in src */
  if (SMCU_NNZ(dst) < SMCU_NNZ(src))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatCopy_cuSparse: the destination matrix has less nonzeros than the source\n");
    return(SUNMAT_ILL_INPUT);
  }

  /* Zero out dst so that copy works correctly */
  if (SUNMatZero_cuSparse(dst) != SUNMAT_SUCCESS)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatCopy_cuSparse: SUNMatZero_cuSparse failed\n");
    return(SUNMAT_OPERATION_FAIL);
  }

  /* Copy the data over */
  retval = SUNMemoryHelper_CopyAsync(SMCU_MEMHELP(src),
                                     SMCU_DATA(dst),
                                     SMCU_DATA(src),
                                     SMCU_NNZ(src)*sizeof(realtype),
                                     (void*) stream);
  if (retval) return(SUNMAT_OPERATION_FAIL);

  /* Copy the row pointers over */
  retval = SUNMemoryHelper_CopyAsync(SMCU_MEMHELP(src),
                                     SMCU_INDEXPTRS(dst),
                                     SMCU_INDEXPTRS(src),
                                     (SMCU_BLOCKROWS(src)+1)*sizeof(int),
                                     (void*) stream);
  if (retval) return(SUNMAT_OPERATION_FAIL);

  /* Copy the column indices over */
  retval = SUNMemoryHelper_CopyAsync(SMCU_MEMHELP(src),
                                     SMCU_INDEXVALS(dst),
                                     SMCU_INDEXVALS(src),
                                     SMCU_BLOCKNNZ(src)*sizeof(int),
                                     (void*) stream);
  if (retval) return(SUNMAT_OPERATION_FAIL);

  return(SUNMAT_SUCCESS);
}


/* Performs A = cA + I. Requires the diagonal to be allocated already. */
int SUNMatScaleAddI_cuSparse(realtype c, SUNMatrix A)
{
  unsigned threadsPerBlock, gridSize;
  cudaStream_t stream = *SMCU_EXECPOLICY(A)->stream();

  switch (SMCU_SPARSETYPE(A))
  {
    case SUNMAT_CUSPARSE_CSR:
      /* Choose the grid size to be the number of rows in the matrix,
        and then choose threadsPerBlock to be a multiple of the warp size
        that results in enough threads to have one per 2 columns. */
      threadsPerBlock = SMCU_EXECPOLICY(A)->blockSize(SMCU_COLUMNS(A)/2);
      gridSize = SMCU_EXECPOLICY(A)->gridSize(SMCU_ROWS(A)*SMCU_COLUMNS(A)/2, threadsPerBlock);
      scaleAddIKernelCSR<realtype, int>
        <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_ROWS(A),
                                                   c,
                                                   SMCU_DATAp(A),
                                                   SMCU_INDEXPTRSp(A),
                                                   SMCU_INDEXVALSp(A));
      break;
    case SUNMAT_CUSPARSE_BCSR:
      /* Choose the grid size to be the number of blocks in the matrix,
         and then choose threadsPerBlock to be a multiple of the warp size
         that results in enough threads to have one per row of the block. */
      threadsPerBlock = SMCU_EXECPOLICY(A)->blockSize(SMCU_BLOCKROWS(A));
      gridSize = SMCU_EXECPOLICY(A)->gridSize(SMCU_NBLOCKS(A)*SMCU_BLOCKROWS(A), threadsPerBlock);
      scaleAddIKernelBCSR<realtype, int>
        <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_BLOCKROWS(A),
                                                   SMCU_NBLOCKS(A),
                                                   SMCU_BLOCKNNZ(A),
                                                   c,
                                                   SMCU_DATAp(A),
                                                   SMCU_INDEXPTRSp(A),
                                                   SMCU_INDEXVALSp(A));
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMatScaleAddI_cuSparse: sparse type not recognized\n");
      return(SUNMAT_ILL_INPUT);
  }

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return(SUNMAT_OPERATION_FAIL);
#endif

  return(SUNMAT_SUCCESS);
}


/* Performs A = cA + B */
int SUNMatScaleAdd_cuSparse(realtype c, SUNMatrix A, SUNMatrix B)
{
  cudaStream_t stream;
  unsigned threadsPerBlock, gridSize;

  if (!SMCompatible_cuSparse(A, B))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatScaleAdd_cuSparse: SUNMatScaleAdd_cuSparse failed\n");
    return(SUNMAT_ILL_INPUT);
  }

  stream = *SMCU_EXECPOLICY(A)->stream();

  switch (SMCU_SPARSETYPE(A))
  {
    case SUNMAT_CUSPARSE_CSR:
      /* Choose the grid size to be the number of rows in the matrix,
        and then choose threadsPerBlock to be a multiple of the warp size
        that results in enough threads to have one per 2 columns. */
      threadsPerBlock = SMCU_EXECPOLICY(A)->blockSize(SMCU_COLUMNS(A)/2);
      gridSize = SMCU_EXECPOLICY(A)->gridSize(SMCU_ROWS(A)*SMCU_COLUMNS(A)/2, threadsPerBlock);
      scaleAddKernelCSR<realtype, int>
        <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_NNZ(A),
                                                   c,
                                                   SMCU_DATAp(A),
                                                   SMCU_DATAp(B));
      break;
    case SUNMAT_CUSPARSE_BCSR:
      /* Choose the grid size to be the number of blocks in the matrix,
         and then choose threadsPerBlock to be a multiple of the warp size
         that results in enough threads to have one per row of the block. */
      threadsPerBlock = SMCU_EXECPOLICY(A)->blockSize(SMCU_BLOCKROWS(A));
      gridSize = SMCU_EXECPOLICY(A)->gridSize(SMCU_NBLOCKS(A)*SMCU_BLOCKROWS(A), threadsPerBlock);
      scaleAddKernelCSR<realtype, int>
        <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_NNZ(A),
                                                   c,
                                                   SMCU_DATAp(A),
                                                   SMCU_DATAp(B));
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMatScaleAdd_cuSparse: sparse type not recognized\n");
      return(SUNMAT_ILL_INPUT);
  }

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return(SUNMAT_OPERATION_FAIL);
#endif

  return(SUNMAT_SUCCESS);
}

/* Setup buffers needed for Matvec */
int SUNMatMatvecSetup_cuSparse(SUNMatrix A)
{
#if CUDART_VERSION >= 11000
  realtype placeholder[1];
  const realtype one = ONE;

  /* Check if setup has already been done */
  if (!(SMCU_CONTENT(A)->matvec_issetup))
  {
    SUNDIALS_CUSPARSE_VERIFY( cusparseCreateDnVec(&SMCU_CONTENT(A)->vecX,
                                                  SMCU_COLUMNS(A),
                                                  placeholder, CUDA_R_XF) );
    SUNDIALS_CUSPARSE_VERIFY( cusparseCreateDnVec(&SMCU_CONTENT(A)->vecY,
                                                  SMCU_ROWS(A),
                                                  placeholder, CUDA_R_XF) );

    SUNDIALS_CUSPARSE_VERIFY(
      cusparseSpMV_bufferSize(SMCU_CUSPHANDLE(A),
                              CUSPARSE_OPERATION_NON_TRANSPOSE,
                              &one, SMCU_CONTENT(A)->spmat_descr,
                              SMCU_CONTENT(A)->vecX, &one, SMCU_CONTENT(A)->vecY,
                              CUDA_R_XF, CUSPARSE_MV_ALG_DEFAULT,
                              &SMCU_CONTENT(A)->bufferSize) );

    if ( SUNMemoryHelper_Alloc(SMCU_MEMHELP(A), &SMCU_CONTENT(A)->dBufferMem,
                               SMCU_CONTENT(A)->bufferSize, SUNMEMTYPE_DEVICE) )
      return(SUNMAT_OPERATION_FAIL);
  }
#endif
  SMCU_CONTENT(A)->matvec_issetup = SUNTRUE;
  return(SUNMAT_SUCCESS);
}

/* Perform y = Ax */
int SUNMatMatvec_cuSparse(SUNMatrix A, N_Vector x, N_Vector y)
{
  /* Verify that the dimensions of A, x, and y agree */
  if ( (SMCU_COLUMNS(A) != N_VGetLength(x)) ||
       (SMCU_ROWS(A) != N_VGetLength(y)) )
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatMatvec_cuSparse: dimensions do not agree\n");
    return(SUNMAT_ILL_INPUT);
  }

  realtype *d_xdata = N_VGetDeviceArrayPointer(x);
  realtype *d_ydata = N_VGetDeviceArrayPointer(y);

  if (SMCU_SPARSETYPE(A) == SUNMAT_CUSPARSE_CSR)
  {
    const realtype one = ONE;

    /* Zero result vector */
    N_VConst(ZERO, y);

#if CUDART_VERSION >= 11000
    {
      /* Setup matvec if it has not been done yet */
      if (!SMCU_CONTENT(A)->matvec_issetup && SUNMatMatvecSetup_cuSparse(A))
      {
        return(SUNMAT_OPERATION_FAIL);
      }

      SUNDIALS_CUSPARSE_VERIFY( cusparseDnVecSetValues(SMCU_CONTENT(A)->vecX,
                                                       d_xdata) );
      SUNDIALS_CUSPARSE_VERIFY( cusparseDnVecSetValues(SMCU_CONTENT(A)->vecY,
                                                       d_ydata) );

      SUNDIALS_CUSPARSE_VERIFY( cusparseSpMV(SMCU_CUSPHANDLE(A),
                                             CUSPARSE_OPERATION_NON_TRANSPOSE,
                                             &one, SMCU_CONTENT(A)->spmat_descr,
                                             SMCU_CONTENT(A)->vecX, &one,
                                             SMCU_CONTENT(A)->vecY, CUDA_R_XF,
                                             CUSPARSE_MV_ALG_DEFAULT,
                                             SMCU_CONTENT(A)->dBufferMem->ptr) );
    }
#else
    SUNDIALS_CUSPARSE_VERIFY(
      cusparseXcsrmv(SMCU_CUSPHANDLE(A), CUSPARSE_OPERATION_NON_TRANSPOSE,
                     SMCU_ROWS(A), SMCU_COLUMNS(A), SMCU_NNZ(A),
                     &one, SMCU_MATDESCR(A), SMCU_DATAp(A), SMCU_INDEXPTRSp(A),
                     SMCU_INDEXVALSp(A), d_xdata, &one, d_ydata) );
#endif
  }
  else if (SMCU_SPARSETYPE(A) == SUNMAT_CUSPARSE_BCSR)
  {
    cudaStream_t stream;
    unsigned gridSize, threadsPerBlock;

    stream = *SMCU_EXECPOLICY(A)->stream();

    /* Choose the grid size to be the number of blocks in the matrix,
       and then choose threadsPerBlock to be a multiple of the warp size
       that results in enough threads to have one per row of the block. */
    threadsPerBlock = SMCU_EXECPOLICY(A)->blockSize(SMCU_COLUMNS(A)/2);
    gridSize = SMCU_EXECPOLICY(A)->gridSize(SMCU_ROWS(A)*SMCU_COLUMNS(A)/2, threadsPerBlock);
    matvecBCSR<realtype, int>
      <<<gridSize, threadsPerBlock, 0, stream>>>(SMCU_BLOCKROWS(A),
                                                 SMCU_NBLOCKS(A),
                                                 SMCU_BLOCKNNZ(A),
                                                 SMCU_DATAp(A),
                                                 SMCU_INDEXPTRSp(A),
                                                 SMCU_INDEXVALSp(A),
                                                 d_xdata,
                                                 d_ydata);
  }
  else
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatMatvec_cuSparse: sparse type not recognized\n");
    return(SUNMAT_ILL_INPUT);
  }

#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
    cudaDeviceSynchronize();
    if (!SUNDIALS_CUDA_VERIFY(cudaGetLastError())) return(SUNMAT_OPERATION_FAIL);
#endif

  return(SUNMAT_SUCCESS);
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
    return(SUNFALSE);

  /* both matrices must have the same shape and sparsity type */
  if (SMCU_ROWS(A) != SMCU_ROWS(B))
    return(SUNFALSE);
  if (SMCU_COLUMNS(A) != SMCU_COLUMNS(B))
    return(SUNFALSE);
  if (SMCU_SPARSETYPE(A) != SMCU_SPARSETYPE(B))
    return(SUNFALSE);

  return(SUNTRUE);
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
    return(NULL);
  }

  /* Attach operations */
  A->ops->getid       = SUNMatGetID_cuSparse;
  A->ops->clone       = SUNMatClone_cuSparse;
  A->ops->destroy     = SUNMatDestroy_cuSparse;
  A->ops->zero        = SUNMatZero_cuSparse;
  A->ops->copy        = SUNMatCopy_cuSparse;
  A->ops->scaleadd    = SUNMatScaleAdd_cuSparse;
  A->ops->scaleaddi   = SUNMatScaleAddI_cuSparse;
  A->ops->matvecsetup = SUNMatMatvecSetup_cuSparse;
  A->ops->matvec      = SUNMatMatvec_cuSparse;

  /* Create content */
  SUNMatrix_Content_cuSparse content = NULL;
  content = (SUNMatrix_Content_cuSparse) malloc(sizeof *content);
  if (content == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMatrix_cuSparse_NewEmpty: failed to malloc content\n");
    SUNMatDestroy(A);
    return(NULL);
  }

  /* Attach content */
  A->content = content;
  content->mem_helper = NULL;

  return(A);
}

#if CUDART_VERSION >= 11000
cusparseStatus_t CreateSpMatDescr(SUNMatrix A, cusparseSpMatDescr_t *spmat_descr)
{
  /* CUDA 11 introduced the "Generic API" and removed the cusparseXcsrmv that
    works on the old cusparseMatDescr_t and raw data arrays. However,
    cuSolverSp stuff requires the cusparseMatDescr_t still. So, we have to
    create this cusparseSpMatDescr_t *and* the cusparseMatDescr_t. */
  return(cusparseCreateCsr(spmat_descr, SMCU_ROWS(A), SMCU_COLUMNS(A),
                           SMCU_NNZ(A), SMCU_INDEXPTRSp(A),
                           SMCU_INDEXVALSp(A), SMCU_DATAp(A),
                           CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                           CUSPARSE_INDEX_BASE_ZERO, CUDA_R_XF));
}
#endif
