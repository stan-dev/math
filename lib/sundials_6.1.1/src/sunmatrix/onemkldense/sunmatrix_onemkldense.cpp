/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * This is the implementation file for the dense implementation of the
 * SUNMATRIX class using the Intel oneAPI Math Kernel Library (oneMKL).
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <CL/sycl.hpp>
#include <oneapi/mkl/blas.hpp>

// SUNDIALS public headers
#include <sunmatrix/sunmatrix_onemkldense.h>
#include <sunmemory/sunmemory_sycl.h>

// SUNDIALS private headers
#include "sundials_debug.h"
#include "sundials_sycl.h"

// Check for a valid precision
#if defined(SUNDIALS_EXTENDED_PRECISION)
#error "oneMLK unsupported precision"
#endif

// Constants
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

// Content accessor macros
#define MAT_CONTENT(A)     ((SUNMatrixContent_OneMklDense) (A->content))
#define MAT_LAST_FLAG(A)   (MAT_CONTENT(A)->last_flag)
#define MAT_BLOCK_ROWS(A)  (MAT_CONTENT(A)->block_rows)
#define MAT_BLOCK_COLS(A)  (MAT_CONTENT(A)->block_cols)
#define MAT_ROWS(A)        (MAT_CONTENT(A)->rows)
#define MAT_COLS(A)        (MAT_CONTENT(A)->cols)
#define MAT_NBLOCKS(A)     (MAT_CONTENT(A)->num_blocks)
#define MAT_LDATA(A)       (MAT_CONTENT(A)->ldata)
#define MAT_DATA(A)        (MAT_CONTENT(A)->data)
#define MAT_DATAp(A)       ((realtype*) MAT_CONTENT(A)->data->ptr)
#define MAT_BLOCKS(A)      (MAT_CONTENT(A)->blocks)
#define MAT_BLOCKSp(A)     ((realtype**) MAT_CONTENT(A)->blocks->ptr)
#define MAT_EXECPOLICY(A)  (MAT_CONTENT(A)->exec_policy)
#define MAT_MEMTYPE(A)     (MAT_CONTENT(A)->mem_type)
#define MAT_MEMHELPER(A)   (MAT_CONTENT(A)->mem_helper)
#define MAT_QUEUE(A)       (MAT_CONTENT(A)->queue)

// Private function prototypes
static booleantype Compatible_AB(SUNMatrix A, SUNMatrix B);
static booleantype Compatible_Axy(SUNMatrix A, N_Vector x, N_Vector y);

// Kernel launch parameters
static int GetKernelParameters(SUNMatrix A, booleantype reduction,
                               size_t& nthreads_total,
                               size_t& nthreads_per_block);


/* --------------------------------------------------------------------------
 * Constructors
 * -------------------------------------------------------------------------- */


SUNMatrix SUNMatrix_OneMklDense(sunindextype M, sunindextype N,
                                SUNMemoryType mem_type,
                                SUNMemoryHelper mem_helper,
                                sycl::queue* queue, SUNContext sunctx)
{
  return SUNMatrix_OneMklDenseBlock(1, M, N, mem_type, mem_helper, queue,
                                    sunctx);
}


SUNMatrix SUNMatrix_OneMklDenseBlock(sunindextype num_blocks, sunindextype M,
                                     sunindextype N, SUNMemoryType mem_type,
                                     SUNMemoryHelper mem_helper,
                                     sycl::queue* queue, SUNContext sunctx)
{
  int retval;

  // Check inputs
  if ( (M <= 0) || (N <= 0) || (num_blocks <= 0) || (!mem_helper) ||
       ((mem_type != SUNMEMTYPE_UVM) && (mem_type != SUNMEMTYPE_DEVICE)))
  {
    SUNDIALS_DEBUG_ERROR("Illegal input\n");
    return NULL;
  }

  // Create an empty matrix object
  SUNMatrix A = SUNMatNewEmpty(sunctx);
  if (!A)
  {
    SUNDIALS_DEBUG_ERROR("SUNMatNewEmpty returned NULL\n");
    return NULL;
  }

  // Attach operations
  A->ops->getid     = SUNMatGetID_OneMklDense;
  A->ops->clone     = SUNMatClone_OneMklDense;
  A->ops->destroy   = SUNMatDestroy_OneMklDense;
  A->ops->zero      = SUNMatZero_OneMklDense;
  A->ops->copy      = SUNMatCopy_OneMklDense;
  A->ops->scaleadd  = SUNMatScaleAdd_OneMklDense;
  A->ops->scaleaddi = SUNMatScaleAddI_OneMklDense;
  A->ops->matvec    = SUNMatMatvec_OneMklDense;
  A->ops->space     = SUNMatSpace_OneMklDense;

  // Create content
  A->content = (SUNMatrixContent_OneMklDense) malloc(sizeof(_SUNMatrixContent_OneMklDense));
  if (!(A->content))
  {
    SUNDIALS_DEBUG_ERROR("Content allocation failed\n");
    SUNMatDestroy(A);
    return NULL;
  }

  // Fill content
  MAT_CONTENT(A)->block_rows  = M;
  MAT_CONTENT(A)->block_cols  = N;
  MAT_CONTENT(A)->rows        = num_blocks * M;
  MAT_CONTENT(A)->cols        = num_blocks * N;
  MAT_CONTENT(A)->num_blocks  = num_blocks;
  MAT_CONTENT(A)->ldata       = M * N * num_blocks;
  MAT_CONTENT(A)->data        = NULL;
  MAT_CONTENT(A)->blocks      = NULL;
  MAT_CONTENT(A)->mem_type    = mem_type;
  MAT_CONTENT(A)->mem_helper  = mem_helper;
  MAT_CONTENT(A)->exec_policy = new sundials::sycl::ThreadDirectExecPolicy(SYCL_BLOCKDIM(queue));
  MAT_CONTENT(A)->queue       = queue;

  // Allocate data
  retval = SUNMemoryHelper_Alloc(MAT_MEMHELPER(A), &(MAT_DATA(A)),
                                 sizeof(realtype) * MAT_LDATA(A), mem_type,
                                 queue);
  if (retval)
  {
    SUNDIALS_DEBUG_ERROR("SUNMemory allocation failed\n");
    SUNMatDestroy(A);
    return NULL;
  }

  if (MAT_NBLOCKS(A) > 1)
  {
    // Allocate array of pointers to block data
    retval = SUNMemoryHelper_Alloc(MAT_MEMHELPER(A), &(MAT_BLOCKS(A)),
                                   sizeof(realtype*) * MAT_NBLOCKS(A), mem_type,
                                   queue);
    if (retval)
    {
      SUNDIALS_DEBUG_ERROR("SUNMemory allocation failed\n");
      SUNMatDestroy(A);
      return NULL;
    }

    size_t nthreads_total, nthreads_per_block;

    if (GetKernelParameters(A, SUNFALSE, nthreads_total, nthreads_per_block))
    {
      SUNDIALS_DEBUG_ERROR("GetKernelParameters returned nonzero\n");
      SUNMatDestroy(A);
      return NULL;
    }

    realtype*  Adata   = MAT_DATAp(A);
    realtype** Ablocks = MAT_BLOCKSp(A);

    // Initialize array of pointers to block data
    SYCL_FOR(queue, nthreads_total, nthreads_per_block, item,
             GRID_STRIDE_XLOOP(item, i, num_blocks)
             {
               Ablocks[i] = Adata + i * M * N;
             });
  }

  return A;
}


/* --------------------------------------------------------------------------
 * Accessor functions
 * -------------------------------------------------------------------------- */


sunindextype SUNMatrix_OneMklDense_Rows(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_ONEMKLDENSE)
    return MAT_ROWS(A);
  else
    return SUNMAT_ILL_INPUT;
}


sunindextype SUNMatrix_OneMklDense_Columns(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_ONEMKLDENSE)
    return MAT_COLS(A);
  else
    return SUNMAT_ILL_INPUT;
}


sunindextype SUNMatrix_OneMklDense_NumBlocks(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_ONEMKLDENSE)
    return MAT_NBLOCKS(A);
  else
    return SUNMAT_ILL_INPUT;
}


sunindextype SUNMatrix_OneMklDense_BlockRows(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_ONEMKLDENSE)
    return MAT_BLOCK_ROWS(A);
  else
    return SUNMAT_ILL_INPUT;
}


sunindextype SUNMatrix_OneMklDense_BlockColumns(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_ONEMKLDENSE)
    return MAT_BLOCK_COLS(A);
  else
    return SUNMAT_ILL_INPUT;
}


sunindextype SUNMatrix_OneMklDense_LData(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_ONEMKLDENSE)
    return MAT_LDATA(A);
  else
    return SUNMAT_ILL_INPUT;
}


realtype* SUNMatrix_OneMklDense_Data(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_ONEMKLDENSE)
    return MAT_DATAp(A);
  else
    return NULL;
}


sunindextype SUNMatrix_OneMklDense_BlockLData(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_ONEMKLDENSE)
    return MAT_BLOCK_ROWS(A) * MAT_BLOCK_COLS(A);
  else
    return SUNMAT_ILL_INPUT;
}


realtype** SUNMatrix_OneMklDense_BlockData(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_ONEMKLDENSE)
    return MAT_BLOCKSp(A);
  else
    return NULL;
}


/* Functions that return pointers to the start of a block, column, or block
   column. These are defined as inline functions in sunmatrix_onemkldense.h, so
   we just mark them as extern here. */

extern realtype* SUNMatrix_OneMklDense_Block(SUNMatrix A, sunindextype k);

extern realtype* SUNMatrix_OneMklDense_Column(SUNMatrix A, sunindextype j);

extern realtype* SUNMatrix_OneMklDense_BlockColumn(SUNMatrix A, sunindextype k,
                                                   sunindextype j);


/* --------------------------------------------------------------------------
 * Utility functions
 * -------------------------------------------------------------------------- */


int SUNMatrix_OneMklDense_CopyToDevice(SUNMatrix A, realtype* h_data)
{
  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Illegal input\n");
    return SUNMAT_ILL_INPUT;
  }

  // Wrap the input pointer
  SUNMemory _h_data = SUNMemoryHelper_Wrap(h_data, SUNMEMTYPE_HOST);
  if (!_h_data)
  {
    SUNDIALS_DEBUG_ERROR("SUNMemory wrap failed\n");
    return SUNMAT_ILL_INPUT;
  }

  // Copy the data
  int copy_fail = SUNMemoryHelper_CopyAsync(MAT_MEMHELPER(A),
                                            MAT_DATA(A),
                                            _h_data,
                                            sizeof(realtype) * MAT_LDATA(A),
                                            MAT_QUEUE(A));

  // Sync with respect to host, but only this queue
  MAT_QUEUE(A)->wait_and_throw();

  int retval = SUNMemoryHelper_Dealloc(MAT_MEMHELPER(A), _h_data, MAT_QUEUE(A));
  if (retval)
  {
    SUNDIALS_DEBUG_ERROR("SUNMemory dealloc failed\n");
    return SUNMAT_MEM_FAIL;
  }

  return (copy_fail ? SUNMAT_MEM_FAIL : SUNMAT_SUCCESS);
}


int SUNMatrix_OneMklDense_CopyFromDevice(SUNMatrix A, realtype* h_data)
{
  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Illegal input\n");
    return SUNMAT_ILL_INPUT;
  }

  SUNMemory _h_data = SUNMemoryHelper_Wrap(h_data, SUNMEMTYPE_HOST);
  if (!_h_data)
  {
    SUNDIALS_DEBUG_ERROR("SUNMemory wrap failed\n");
    return SUNMAT_MEM_FAIL;
  }

  int copy_fail = SUNMemoryHelper_CopyAsync(MAT_MEMHELPER(A),
                                            _h_data,
                                            MAT_DATA(A),
                                            sizeof(realtype) * MAT_LDATA(A),
                                            MAT_QUEUE(A));

  // Sync with respect to host, but only this queue
  MAT_QUEUE(A)->wait_and_throw();

  int retval = SUNMemoryHelper_Dealloc(MAT_MEMHELPER(A), _h_data, MAT_QUEUE(A));
  if (retval)
  {
    SUNDIALS_DEBUG_ERROR("SUNMemory dealloc failed\n");
    return SUNMAT_MEM_FAIL;
  }

  return (copy_fail ? SUNMAT_MEM_FAIL : SUNMAT_SUCCESS);
}


/* --------------------------------------------------------------------------
 * Implementation of generic SUNMatrix operations.
 * -------------------------------------------------------------------------- */


SUNMatrix SUNMatClone_OneMklDense(SUNMatrix A)
{
  if (!A)
  {
    SUNDIALS_DEBUG_ERROR("Input matrix is NULL\n");
    return NULL;
  }

  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Invalid matrix ID\n");
    return NULL;
  }

  SUNMatrix B = SUNMatrix_OneMklDenseBlock(MAT_NBLOCKS(A),
                                           MAT_BLOCK_ROWS(A),
                                           MAT_BLOCK_COLS(A),
                                           MAT_DATA(A)->type,
                                           MAT_MEMHELPER(A),
                                           MAT_QUEUE(A),
                                           A->sunctx);

  if (!B)
  {
    SUNDIALS_DEBUG_ERROR("Output matrix is NULL\n");
    return NULL;
  }

  return B;
}


void SUNMatDestroy_OneMklDense(SUNMatrix A)
{
  if (!A)
  {
    SUNDIALS_DEBUG_ERROR("Input matrix is NULL\n");
    return;
  }

  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Invalid matrix ID\n");
    return;
  }

  // Free content
  if (A->content)
  {
    // Free data array(s)
    if (MAT_DATA(A)) SUNMemoryHelper_Dealloc(MAT_MEMHELPER(A), MAT_DATA(A),
                                             MAT_QUEUE(A));
    if (MAT_BLOCKS(A)) SUNMemoryHelper_Dealloc(MAT_MEMHELPER(A), MAT_BLOCKS(A),
                                               MAT_QUEUE(A));

    // Free content struct
    free(A->content);
    A->content = NULL;
  }

  // Free matrix
  SUNMatFreeEmpty(A);
  A = NULL;

  return;
}


int SUNMatZero_OneMklDense(SUNMatrix A)
{
  if (!A)
  {
    SUNDIALS_DEBUG_ERROR("Input matrix is NULL\n");
    return SUNMAT_ILL_INPUT;
  }

  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Invalid matrix ID\n");
    return SUNMAT_ILL_INPUT;
  }

  const sunindextype ldata  = MAT_LDATA(A);
  realtype           *Adata = MAT_DATAp(A);
  sycl::queue        *Q     = MAT_QUEUE(A);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(A, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_ERROR("GetKernelParameters returned nonzero\n");
    return SUNMAT_MEM_FAIL;
  }

  // Zero out matrix
  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, ldata)
           {
             Adata[i] = ZERO;
           });

  return SUNMAT_SUCCESS;
}


int SUNMatCopy_OneMklDense(SUNMatrix A, SUNMatrix B)
{
  if (!A || !B)
  {
    SUNDIALS_DEBUG_ERROR("An input matrix is NULL\n");
    return SUNMAT_ILL_INPUT;
  }

  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Invalid matrix ID\n");
    return SUNMAT_ILL_INPUT;
  }

  // Verify that A and B are compatible
  if (!Compatible_AB(A, B))
  {
    SUNDIALS_DEBUG_ERROR("Input matrices are incompatible\n");
    return SUNMAT_ILL_INPUT;
  }

  const sunindextype ldata  = MAT_LDATA(A);
  realtype           *Adata = MAT_DATAp(A);
  realtype           *Bdata = MAT_DATAp(B);
  sycl::queue        *Q     = MAT_QUEUE(A);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(A, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_ERROR("GetKernelParameters returned nonzero\n");
    return SUNMAT_MEM_FAIL;
  }

  // Copy A into B
  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, ldata)
           {
             Bdata[i] = Adata[i];
           });

  return SUNMAT_SUCCESS;
}


int SUNMatScaleAddI_OneMklDense(realtype c, SUNMatrix A)
{
  if (!A)
  {
    SUNDIALS_DEBUG_ERROR("Input matrix is NULL\n");
    return SUNMAT_ILL_INPUT;
  }

  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Invalid matrix ID\n");
    return SUNMAT_ILL_INPUT;
  }

  const size_t M     = static_cast<size_t>(MAT_BLOCK_ROWS(A));
  const size_t N     = static_cast<size_t>(MAT_BLOCK_COLS(A));
  const size_t B     = static_cast<size_t>(MAT_NBLOCKS(A));
  realtype*    Adata = MAT_DATAp(A);
  sycl::queue* Q     = MAT_QUEUE(A);

  // Compute A = c * A + I
  Q->submit([&](sycl::handler& h)
  {
    h.parallel_for(sycl::range{M, N, B}, [=](sycl::id<3> idx)
    {
      sunindextype i = idx[0];
      sunindextype j = idx[1];
      sunindextype k = idx[2];

      // Index into 1D data array
      sunindextype tid = k * M * N + j * M + i;

      if (i == j)
      {
        Adata[tid] = c * Adata[tid] + ONE;
      }
      else
      {
        Adata[tid] = c * Adata[tid];
      }
    });
  });

  return SUNMAT_SUCCESS;
}


int SUNMatScaleAdd_OneMklDense(realtype c, SUNMatrix A, SUNMatrix B)
{
  if (!A || !B)
  {
    SUNDIALS_DEBUG_ERROR("An input matrix is NULL\n");
    return SUNMAT_ILL_INPUT;
  }

  if (!Compatible_AB(A, B))
  {
    SUNDIALS_DEBUG_ERROR("Input matrices are incompatible\n");
    return SUNMAT_ILL_INPUT;
  }

  const sunindextype ldata  = MAT_LDATA(A);
  realtype           *Adata = MAT_DATAp(A);
  realtype           *Bdata = MAT_DATAp(B);
  sycl::queue        *Q     = MAT_QUEUE(A);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(A, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_ERROR("GetKernelParameters returned nonzero\n");
    return SUNMAT_MEM_FAIL;
  }

  // Compute A = c * A + B
  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, ldata)
           {
             Adata[i] = c * Adata[i] + Bdata[i];
           });

  return SUNMAT_SUCCESS;
}


int SUNMatMatvec_OneMklDense(SUNMatrix A, N_Vector x, N_Vector y)
{
  if (!A || !x || !y)
  {
    SUNDIALS_DEBUG_ERROR("Input matrix or vectors are NULL\n");
    return SUNMAT_ILL_INPUT;
  }

  if (!Compatible_Axy(A, x, y))
  {
    SUNDIALS_DEBUG_ERROR("Input matrix and vectors are incompatible\n");
    return SUNMAT_ILL_INPUT;
  }

  if (MAT_NBLOCKS(A) > 1)
  {
    sycl::queue* Q = MAT_QUEUE(A);
    sunindextype M = MAT_BLOCK_ROWS(A);
    sunindextype N = MAT_BLOCK_COLS(A);

    // TODO(DJG): Replace with batched function
    for (sunindextype i = 0; i < MAT_NBLOCKS(A); i++)
    {
      const realtype* Adata = MAT_DATAp(A) + i * M * N;
      const realtype* xdata = N_VGetDeviceArrayPointer(x) + i * N;
      realtype*       ydata = N_VGetDeviceArrayPointer(y) + i * M;

      // Copmute y = a * A * x + b * y
      oneapi::mkl::blas::gemv(*Q, oneapi::mkl::transpose::N, M, N, ONE, Adata,
                              M, xdata, 1, ZERO, ydata, 1);
    }
  }
  else
  {
    sycl::queue*    Q     = MAT_QUEUE(A);
    sunindextype    M     = MAT_ROWS(A);
    sunindextype    N     = MAT_COLS(A);
    const realtype* Adata = MAT_DATAp(A);
    const realtype* xdata = N_VGetDeviceArrayPointer(x);
    realtype*       ydata = N_VGetDeviceArrayPointer(y);

    // Copmute y = a * A * x + b * y
    oneapi::mkl::blas::gemv(*Q, oneapi::mkl::transpose::N, M, N, ONE, Adata, M,
                            xdata, 1, ZERO, ydata, 1);
  }

  return SUNMAT_SUCCESS;
}


int SUNMatSpace_OneMklDense(SUNMatrix A, long int *lenrw, long int *leniw)
{
  if (!A)
  {
    SUNDIALS_DEBUG_ERROR("Input matrix is NULL\n");
    return SUNMAT_ILL_INPUT;
  }

  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Invalid matrix ID\n");
    return SUNMAT_ILL_INPUT;
  }

  *lenrw = MAT_LDATA(A);
  *leniw = 4;

  return SUNMAT_SUCCESS;
}


/* --------------------------------------------------------------------------
 * Private functions
 * -------------------------------------------------------------------------- */


// Get the kernel launch parameters
static int GetKernelParameters(SUNMatrix A, booleantype reduction,
                               size_t& nthreads_total,
                               size_t& nthreads_per_block)
{
  if (!MAT_EXECPOLICY(A))
  {
    SUNDIALS_DEBUG_ERROR("The execution policy is NULL\n");
    return -1;
  }

  /* Get the number of threads per block and total number threads */
  nthreads_per_block = MAT_EXECPOLICY(A)->blockSize();
  nthreads_total     = nthreads_per_block *
                       MAT_EXECPOLICY(A)->gridSize(MAT_LDATA(A));

  if (nthreads_per_block == 0)
  {
    SUNDIALS_DEBUG_ERROR("The number of threads per block must be > 0\n");
    return -1;
  }

  if (nthreads_total == 0)
  {
    SUNDIALS_DEBUG_ERROR("The total number of threads must be > 0\n");
    return -1;
  }

  return 0;
}


static booleantype Compatible_AB(SUNMatrix A, SUNMatrix B)
{
  // Both matrices must have the SUNMATRIEX_MKLDENSE ID
  if (SUNMatGetID(A) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Illegal matrix ID\n");
    return SUNFALSE;
  }

  if (SUNMatGetID(B) != SUNMATRIX_ONEMKLDENSE)
  {
    SUNDIALS_DEBUG_ERROR("Illegal matrix ID\n");
    return SUNFALSE;
  }

  // Both matrices must have the same shape
  if (MAT_BLOCK_ROWS(A) != MAT_BLOCK_ROWS(B))
  {
    SUNDIALS_DEBUG_ERROR("Number of block rows do not match\n");
    return SUNFALSE;
  }

  if (MAT_BLOCK_COLS(A) != MAT_BLOCK_COLS(B))
  {
    SUNDIALS_DEBUG_ERROR("Number of block columns do not match\n");
    return SUNFALSE;
  }

  if (MAT_NBLOCKS(A) != MAT_NBLOCKS(B))
  {
    SUNDIALS_DEBUG_ERROR("Number of blocks do not match\n");
    return SUNFALSE;
  }

  return SUNTRUE;
}


static booleantype Compatible_Axy(SUNMatrix A, N_Vector x, N_Vector y)
{
  //  Vectors must implement N_VGetDeviceArrayPointer
  if (!(x->ops->nvgetdevicearraypointer) || !(y->ops->nvgetdevicearraypointer))
  {
    SUNDIALS_DEBUG_ERROR("Vectors do not have GetDeviceArrayPointer\n");
    return SUNFALSE;
  }

  // Inner dimensions must agree
  if (MAT_COLS(A) != N_VGetLength(x))
  {
    SUNDIALS_DEBUG_ERROR("Number of columns != input vectors length\n");
    return SUNFALSE;
  }

  // Outer dimensions must agree
  if (MAT_ROWS(A) != N_VGetLength(y))
  {
    SUNDIALS_DEBUG_ERROR("Number of rows != output vector length\n");
    return SUNFALSE;
  }

  return SUNTRUE;
}
