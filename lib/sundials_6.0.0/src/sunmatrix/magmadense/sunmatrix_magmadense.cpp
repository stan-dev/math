/* -----------------------------------------------------------------
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
 * This is the implementation file for the dense implementation of
 * the SUNMATRIX package based on MAGMA.
 * -----------------------------------------------------------------*/

#include <algorithm>
#include <sunmatrix/sunmatrix_magmadense.h>

#if defined(SUNDIALS_MAGMA_BACKENDS_CUDA)

#include "sundials_cuda.h"
#include "dense_cuda_kernels.cuh"
using namespace sundials::sunmatrix_gpudense::cuda;
#define SUNDIALS_HIP_OR_CUDA(a,b) b

#elif defined(SUNDIALS_MAGMA_BACKENDS_HIP)

#include "sundials_hip.h"
#include "dense_hip_kernels.hip.hpp"
using namespace sundials::sunmatrix_gpudense::hip;
#define SUNDIALS_HIP_OR_CUDA(a,b) a

#endif

/* Content accessor macro */
#define SMLD_CONTENT(A)  ( (SUNMatrixContent_MagmaDense) (A->content) )

/* Constants */
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Macros for magma operations based on precision */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define xgemv(q,...) magma_dgemv(__VA_ARGS__, q)
#define xgemv_batched(q,...) magmablas_dgemv_batched(__VA_ARGS__, q)
#define magma_xset_pointer(q,...) magma_dset_pointer(__VA_ARGS__, q)
#define xprint(q,...) magma_dprint_gpu(__VA_ARGS__, q)
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define xgemv(q,...) magma_sgemv(__VA_ARGS__, q)
#define xgemv_batched(q,...) magmablas_sgemv_batched(__VA_ARGS__, q)
#define magma_xset_pointer(q,...) magma_sset_pointer(__VA_ARGS__, q)
#define xprint(q,...) magma_sprint_gpu(__VA_ARGS__, q)
#else
#error unsupported precision
#endif

/* Private function prototypes */
static booleantype SMCompatible_MagmaDense(SUNMatrix A, SUNMatrix B);
static booleantype SMCompatible2_MagmaDense(SUNMatrix A, N_Vector x, N_Vector y);

/*
 * ----------------------------------------------------------------------------
 * Implementation specific routines
 * ----------------------------------------------------------------------------
 */

/*
 * Constructor functions
 */

SUNMatrix SUNMatrix_MagmaDense(sunindextype M, sunindextype N, SUNMemoryType memtype,
                               SUNMemoryHelper memhelper,
                               void* queue, SUNContext sunctx)
{
  return(SUNMatrix_MagmaDenseBlock(1, M, N, memtype, memhelper, queue, sunctx));
}

SUNMatrix SUNMatrix_MagmaDenseBlock(sunindextype nblocks, sunindextype M, sunindextype N,
                                    SUNMemoryType memtype, SUNMemoryHelper memhelper,
                                    void* queue, SUNContext sunctx)
{
  SUNMatrix Amat;
  SUNMatrixContent_MagmaDense A;
  int retval;

  /* Return with NULL matrix on illegal dimension input */
  if ( (M <= 0) || (N <= 0) || (nblocks <= 0))
    return(NULL);

  /* Check for valid memory type options */
  if ((memtype != SUNMEMTYPE_UVM) && (memtype != SUNMEMTYPE_DEVICE))
    return(NULL);

  /* Check for valid memory helper */
  if (memhelper == NULL)
    return(NULL);

  /* First thing we do is initialize magma */
  retval = magma_init();
  if (retval != MAGMA_SUCCESS) return(NULL);

  /* Create an empty matrix object */
  Amat = NULL;
  Amat = SUNMatNewEmpty(sunctx);
  if (Amat == NULL) return(NULL);

  /* Attach operations */
  Amat->ops->getid       = SUNMatGetID_MagmaDense;
  Amat->ops->clone       = SUNMatClone_MagmaDense;
  Amat->ops->destroy     = SUNMatDestroy_MagmaDense;
  Amat->ops->zero        = SUNMatZero_MagmaDense;
  Amat->ops->copy        = SUNMatCopy_MagmaDense;
  Amat->ops->scaleadd    = SUNMatScaleAdd_MagmaDense;
  Amat->ops->scaleaddi   = SUNMatScaleAddI_MagmaDense;
  Amat->ops->matvecsetup = SUNMatMatvecSetup_MagmaDense;
  Amat->ops->matvec      = SUNMatMatvec_MagmaDense;
  Amat->ops->space       = SUNMatSpace_MagmaDense;

  /* Create content */
  A = NULL;
  A = (SUNMatrixContent_MagmaDense) malloc(sizeof(*A));
  if (A == NULL) { SUNMatDestroy(Amat); return(NULL); }

  /* Attach content */
  Amat->content = A;

  /* Fill content */
  A->M       = M;
  A->N       = N;
  A->nblocks = nblocks;
  A->ldata   = M*N*nblocks;
  A->data    = NULL;
  A->blocks  = NULL;
  A->xblocks = NULL;
  A->yblocks = NULL;
  A->memhelp = memhelper;
  A->q       = NULL;

  magma_getdevice(&A->device_id);
  SUNDIALS_HIP_OR_CUDA(
    magma_queue_create_from_hip(A->device_id, (hipStream_t) queue, NULL, NULL, &A->q);,
    magma_queue_create_from_cuda(A->device_id, (cudaStream_t) queue, NULL, NULL, &A->q); )

  /* Allocate data */
  retval = SUNMemoryHelper_Alloc(A->memhelp, &A->data,
                                 sizeof(realtype) * A->ldata, memtype, nullptr);
  if (retval) { SUNMatDestroy(Amat); return(NULL); }

  if (A->nblocks > 1)
  {
    /* Allocate array of pointers to block data */
    retval = SUNMemoryHelper_Alloc(A->memhelp, &A->blocks,
                                   sizeof(realtype*) * A->nblocks, memtype,
                                   nullptr);
    if (retval) { SUNMatDestroy(Amat); return(NULL); }

    /* Initialize array of pointers to block data */
    magma_xset_pointer(A->q, (realtype**)A->blocks->ptr, (realtype*)A->data->ptr,
                       A->M, 0, 0, A->M*A->N, A->nblocks);
  }

  return(Amat);
}

/*
 * Accessor functions
 */

sunindextype SUNMatrix_MagmaDense_Rows(SUNMatrix Amat)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_MAGMADENSE)
    return(A->M * A->nblocks);
  else
    return(SUNMAT_ILL_INPUT);
}

sunindextype SUNMatrix_MagmaDense_Columns(SUNMatrix Amat)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_MAGMADENSE)
    return(A->N * A->nblocks);
  else
    return(SUNMAT_ILL_INPUT);
}

sunindextype SUNMatrix_MagmaDense_BlockRows(SUNMatrix Amat)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_MAGMADENSE)
    return(A->M);
  else
    return(SUNMAT_ILL_INPUT);
}

sunindextype SUNMatrix_MagmaDense_BlockColumns(SUNMatrix Amat)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_MAGMADENSE)
    return(A->N);
  else
    return(SUNMAT_ILL_INPUT);
}

sunindextype SUNMatrix_MagmaDense_NumBlocks(SUNMatrix Amat)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_MAGMADENSE)
    return(A->nblocks);
  else
    return(SUNMAT_ILL_INPUT);
}

sunindextype SUNMatrix_MagmaDense_LData(SUNMatrix Amat)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_MAGMADENSE)
    return(A->ldata);
  else
    return(SUNMAT_ILL_INPUT);
}

realtype* SUNMatrix_MagmaDense_Data(SUNMatrix Amat)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_MAGMADENSE)
    return((realtype*) A->data->ptr);
  else
    return(NULL);
}

realtype** SUNMatrix_MagmaDense_BlockData(SUNMatrix Amat)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (SUNMatGetID(Amat) == SUNMATRIX_MAGMADENSE)
    return((realtype**) A->blocks->ptr);
  else
    return(NULL);
}

extern realtype* SUNMatrix_MagmaDense_Block(SUNMatrix Amat, sunindextype k);

extern realtype* SUNMatrix_MagmaDense_Column(SUNMatrix Amat, sunindextype j);

extern realtype* SUNMatrix_MagmaDense_BlockColumn(SUNMatrix Amat, sunindextype k, sunindextype j);

/*
 * Utility functions
 */

void SUNMatrix_MagmaDense_Print(SUNMatrix Amat)
{
  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return;

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  for (sunindextype k = 0; k < A->nblocks; k++)
    xprint(A->q, A->M, A->N, SUNMatrix_MagmaDense_Block(Amat,k), A->M);
}

int SUNMatrix_MagmaDense_CopyToDevice(SUNMatrix Amat, realtype* h_data)
{
  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return(SUNMAT_ILL_INPUT);
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  int retval = 0;
  SUNMemory _h_data = SUNMemoryHelper_Wrap(h_data, SUNMEMTYPE_HOST);
  SUNDIALS_HIP_OR_CUDA( hipStream_t stream = magma_queue_get_hip_stream(A->q);,
                        cudaStream_t stream = magma_queue_get_cuda_stream(A->q); )

  retval = SUNMemoryHelper_CopyAsync(A->memhelp,
                                     A->data,
                                     _h_data,
                                     sizeof(realtype) * A->ldata,
                                     (void*) &stream);
  magma_queue_sync(A->q); /* sync with respect to host, but only this stream */

  SUNMemoryHelper_Dealloc(A->memhelp, _h_data, nullptr);
  return(retval == 0 ? SUNMAT_SUCCESS : SUNMAT_MEM_FAIL);
}

int SUNMatrix_MagmaDense_CopyFromDevice(SUNMatrix Amat, realtype* h_data)
{
  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return(SUNMAT_ILL_INPUT);
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  int retval = 0;
  SUNMemory _h_data = SUNMemoryHelper_Wrap(h_data, SUNMEMTYPE_HOST);
  SUNDIALS_HIP_OR_CUDA( hipStream_t stream = magma_queue_get_hip_stream(A->q);,
                        cudaStream_t stream = magma_queue_get_cuda_stream(A->q); )

  retval = SUNMemoryHelper_CopyAsync(A->memhelp,
                                     _h_data,
                                     A->data,
                                     sizeof(realtype) * A->ldata,
                                     (void*) &stream);
  magma_queue_sync(A->q); /* sync with respect to host, but only this stream */

  SUNMemoryHelper_Dealloc(A->memhelp, _h_data, nullptr);
  return(retval == 0 ? SUNMAT_SUCCESS : SUNMAT_MEM_FAIL);
}

/*
 * -----------------------------------------------------------------
 * Implementation of generic SUNMatrix operations.
 * -----------------------------------------------------------------
 */

SUNMatrix SUNMatClone_MagmaDense(SUNMatrix Amat)
{
  if (Amat == NULL) return(NULL);

  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return(NULL);

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);
  SUNMatrix B = NULL;
  SUNDIALS_HIP_OR_CUDA( hipStream_t stream = magma_queue_get_hip_stream(A->q);,
                        cudaStream_t stream = magma_queue_get_cuda_stream(A->q); )

  if (A->nblocks > 1)
    B = SUNMatrix_MagmaDenseBlock(A->nblocks, A->M, A->N, A->data->type,
                                  A->memhelp, stream, Amat->sunctx);
  else
    B = SUNMatrix_MagmaDense(A->M, A->N, A->data->type, A->memhelp, stream,
                             Amat->sunctx);

  return(B);
}

void SUNMatDestroy_MagmaDense(SUNMatrix Amat)
{
  if (Amat == NULL) return;

  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return;

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  /* Sync before destroying */
  magma_queue_sync(A->q);

  /* Free content */
  if (A)
  {
    /* Free data array(s) */
    if (A->data) SUNMemoryHelper_Dealloc(A->memhelp, A->data, nullptr);
    if (A->blocks) SUNMemoryHelper_Dealloc(A->memhelp, A->blocks, nullptr);
    if (A->xblocks) SUNMemoryHelper_Dealloc(A->memhelp, A->xblocks, nullptr);
    if (A->yblocks) SUNMemoryHelper_Dealloc(A->memhelp, A->yblocks, nullptr);
    magma_queue_destroy(A->q);
    /* Free content struct */
    free(A);
    Amat->content = NULL;
  }

  /* Free ops */
  if (Amat->ops)
  {
    free(Amat->ops);
    Amat->ops = NULL;
  }

  /* Free matrix */
  free(Amat);
  Amat = NULL;

  /* Call magma_finalize, but note that magma_finalize does
     nothing until it has been called the same number of times
     as magma_init */
  magma_finalize();

  return;
}

int SUNMatZero_MagmaDense(SUNMatrix Amat)
{
  if (Amat == NULL) return(SUNMAT_ILL_INPUT);

  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return(SUNMAT_ILL_INPUT);

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  /* Zero out matrix */
  SUNDIALS_LAUNCH_KERNEL(SUNDIALS_KERNEL_NAME(zeroKernel<realtype,sunindextype>),
    dim3(std::min<sunindextype>(A->nblocks,INT_MAX),1,1),
    SUNDIALS_HIP_OR_CUDA( dim3(1,16,16), dim3(1,16,32) ), /* We choose slightly larger thread blocks when using HIP since the warps are larger */
    0,
    SUNDIALS_HIP_OR_CUDA( magma_queue_get_hip_stream(A->q), magma_queue_get_cuda_stream(A->q) ),
    A->M,
    A->N,
    A->nblocks,
    (realtype*) A->data->ptr
  );

  return(SUNMAT_SUCCESS);
}

int SUNMatCopy_MagmaDense(SUNMatrix Amat, SUNMatrix Bmat)
{
  if ((Amat == NULL) || (Bmat == NULL)) return(SUNMAT_ILL_INPUT);

  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return(SUNMAT_ILL_INPUT);

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);
  SUNMatrixContent_MagmaDense B = SMLD_CONTENT(Bmat);

  /* Verify that A and B are compatible */
  if (!SMCompatible_MagmaDense(Amat, Bmat))
    return SUNMAT_ILL_INPUT;

  /* Copy A into B */
  SUNDIALS_LAUNCH_KERNEL(SUNDIALS_KERNEL_NAME(copyKernel<realtype,sunindextype>),
    dim3(std::min<sunindextype>(A->nblocks,INT_MAX),1,1),
    SUNDIALS_HIP_OR_CUDA( dim3(1,16,16), dim3(1,16,32) ),
    0,
    SUNDIALS_HIP_OR_CUDA( magma_queue_get_hip_stream(A->q), magma_queue_get_cuda_stream(A->q) ),
    A->M,
    A->N,
    A->nblocks,
    (const realtype*) A->data->ptr,
    (realtype*) B->data->ptr
  );

  return(SUNMAT_SUCCESS);
}

int SUNMatScaleAddI_MagmaDense(realtype c, SUNMatrix Amat)
{
  if (Amat == NULL) return(SUNMAT_ILL_INPUT);

  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return(SUNMAT_ILL_INPUT);

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);


  SUNDIALS_LAUNCH_KERNEL(SUNDIALS_KERNEL_NAME(scaleAddIKernel<realtype,sunindextype>),
    dim3(std::min<sunindextype>(A->nblocks,INT_MAX),1,1),
    SUNDIALS_HIP_OR_CUDA( dim3(1,16,16), dim3(1,16,32) ),
    0,
    SUNDIALS_HIP_OR_CUDA( magma_queue_get_hip_stream(A->q), magma_queue_get_cuda_stream(A->q) ),
    A->M,
    A->N,
    A->nblocks,
    c,
    (realtype*) A->data->ptr
  );

  return(SUNMAT_SUCCESS);
}

int SUNMatScaleAdd_MagmaDense(realtype c, SUNMatrix Amat, SUNMatrix Bmat)
{
  if ((Amat == NULL) || (Bmat == NULL)) return(SUNMAT_ILL_INPUT);

  if ((SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) ||
      !SMCompatible_MagmaDense(Amat, Bmat))
    return(SUNMAT_ILL_INPUT);

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);
  SUNMatrixContent_MagmaDense B = SMLD_CONTENT(Bmat);

  SUNDIALS_LAUNCH_KERNEL(SUNDIALS_KERNEL_NAME(scaleAddKernel<realtype,sunindextype>),
    dim3(std::min<sunindextype>(A->nblocks,INT_MAX),1,1),
    SUNDIALS_HIP_OR_CUDA( dim3(1,16,16), dim3(1,16,32) ),
    0,
    SUNDIALS_HIP_OR_CUDA( magma_queue_get_hip_stream(A->q), magma_queue_get_cuda_stream(A->q) ),
    A->M,
    A->N,
    A->nblocks,
    c,
    (realtype*) A->data->ptr,
    (const realtype*) B->data->ptr
  );

  return(SUNMAT_SUCCESS);
}

int SUNMatMatvecSetup_MagmaDense(SUNMatrix Amat)
{
  int retval = 0;

  if (Amat == NULL) return(SUNMAT_ILL_INPUT);

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (A->nblocks > 1)
  {
    /* Allocate array of pointers to blocks on device */
    if (A->xblocks == NULL)
      retval = SUNMemoryHelper_Alloc(A->memhelp, &A->xblocks,
                                     sizeof(realtype*) * A->nblocks,
                                     A->data->type, nullptr);
    if (retval) return(SUNMAT_MEM_FAIL);

    if (A->yblocks == NULL)
      retval = SUNMemoryHelper_Alloc(A->memhelp, &A->yblocks,
                                     sizeof(realtype*) * A->nblocks,
                                     A->data->type, nullptr);
    if (retval) return(SUNMAT_MEM_FAIL);
  }

  return(SUNMAT_SUCCESS);
}

int SUNMatMatvec_MagmaDense(SUNMatrix Amat, N_Vector x, N_Vector y)
{
  if ((Amat == NULL) || (x == NULL) || (y == NULL)) return(SUNMAT_ILL_INPUT);

  if ((SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) ||
      !SMCompatible2_MagmaDense(Amat, x, y))
    return(SUNMAT_ILL_INPUT);

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  if (A->nblocks > 1)
  {
    /* First, we need to create an array of pointers to the matrix and vector blocks */
    SUNDIALS_LAUNCH_KERNEL(SUNDIALS_KERNEL_NAME(getBlockPointers<realtype,sunindextype>),
      A->nblocks,
      256,
      0,
      SUNDIALS_HIP_OR_CUDA( magma_queue_get_hip_stream(A->q), magma_queue_get_cuda_stream(A->q) ),
      A->M,
      A->N,
      A->nblocks,
      (realtype*)A->data->ptr,
      (realtype**)A->blocks->ptr,
      (realtype*)N_VGetDeviceArrayPointer(x),
      (realtype**)A->xblocks->ptr,
      (realtype*)N_VGetDeviceArrayPointer(y),
      (realtype**)A->yblocks->ptr
    );

    /* Now we can use a batched gemv to do y = alpha*A*x + beta*y where A is block diagonal */
    xgemv_batched(
      A->q,         /* queue/stream to execute in */
      MagmaNoTrans, /* use A not A^T */
      A->M,         /* number of rows for a block */
      A->N,         /* number of cols for a block */
      ONE,          /* alpha */
      (realtype**)A->blocks->ptr,
      A->M,         /* leading dimension of A */
      (realtype**)A->xblocks->ptr,
      1,            /* increment (stride) of xblocks */
      ZERO,         /* beta */
      (realtype**)A->yblocks->ptr,
      1,            /* increment (stride) of yblocks */
      A->nblocks    /* number of blocks */
    );
  }
  else
  {
    /* Now we can use gemv to do y = alpha*A*x + beta*y */
    xgemv(
      A->q,         /* queue/stream to execute in */
      MagmaNoTrans, /* use A not A^T */
      A->M,         /* number of rows */
      A->N,         /* number of cols */
      ONE,          /* alpha */
      (const realtype*)A->data->ptr,
      A->M,         /* leading dimension of A */
      (const realtype*)N_VGetDeviceArrayPointer(x),
      1,            /* increment for x data */
      ZERO,         /* beta */
      (realtype*)N_VGetDeviceArrayPointer(y),
      1             /* increment for y data */
    );
  }

  return(SUNMAT_SUCCESS);
}

int SUNMatSpace_MagmaDense(SUNMatrix Amat, long int *lenrw, long int *leniw)
{
  if (Amat == NULL) return(SUNMAT_ILL_INPUT);

  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE) return(SUNMAT_ILL_INPUT);

  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  *lenrw = A->ldata;
  *leniw = 4;

  return(SUNMAT_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * Private functions
 * -----------------------------------------------------------------
 */

static booleantype SMCompatible_MagmaDense(SUNMatrix Amat, SUNMatrix Bmat)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);
  SUNMatrixContent_MagmaDense B = SMLD_CONTENT(Bmat);

  /* Both matrices must be SUNMATRIX_MAGMADENSE */
  if (SUNMatGetID(Amat) != SUNMATRIX_MAGMADENSE)
    return(SUNFALSE);
  if (SUNMatGetID(Bmat) != SUNMATRIX_MAGMADENSE)
    return(SUNFALSE);

  /* Both matrices must have the same shape */
  if (A->M != B->M)
    return(SUNFALSE);
  if (A->N != B->N)
    return(SUNFALSE);
  if (A->nblocks != B->nblocks)
    return(SUNFALSE);

  return(SUNTRUE);
}

static booleantype SMCompatible2_MagmaDense(SUNMatrix Amat, N_Vector x, N_Vector y)
{
  SUNMatrixContent_MagmaDense A = SMLD_CONTENT(Amat);

  /*  Vectors must implement N_VGetDeviceArrayPointer */
  if (x->ops->nvgetdevicearraypointer == NULL ||
      y->ops->nvgetdevicearraypointer == NULL)
    return(SUNFALSE);

  /* Inner dimensions must agree */
  if (A->N*A->nblocks != N_VGetLength(x))
    return(SUNFALSE);

  /* Outer dimensions must agree */
  if (A->M*A->nblocks != N_VGetLength(y))
    return(SUNFALSE);

  return(SUNTRUE);
}
