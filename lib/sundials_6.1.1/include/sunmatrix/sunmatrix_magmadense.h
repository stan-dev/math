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
 * This is the header file for the dense implementation of the
 * SUNMATRIX module, SUNMATRIX_MAGMADENSE.
 * -----------------------------------------------------------------
 */


#ifndef _SUNMATRIX_MAGMADENSE_H
#define _SUNMATRIX_MAGMADENSE_H

#include <stdio.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_memory.h>

#if defined(SUNDIALS_MAGMA_BACKENDS_CUDA)
#define HAVE_CUBLAS
#elif defined(SUNDIALS_MAGMA_BACKENDS_HIP)
#define HAVE_HIP
#endif
#include <magma_v2.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

struct _SUNMatrixContent_MagmaDense {
  int             last_flag; /* last error code returned by magma  */
  int             device_id; /* device ID used by magma            */
  sunindextype    M;         /* number of rows in block            */
  sunindextype    N;         /* number of columns in block         */
  sunindextype    nblocks;   /* number of blocks in matrix         */
  sunindextype    ldata;     /* length of data array               */
  SUNMemory       data;      /* matrix data; column-major          */
  SUNMemory       blocks;    /* device pointers to blocks of A     */
  SUNMemory       xblocks;   /* device pointers to blocks of x     */
  SUNMemory       yblocks;   /* device pointers to blocks of y     */
  SUNMemoryHelper memhelp;   /* memory helper                      */
  magma_queue_t   q;         /* operation queue (i.e. stream)      */
};

typedef struct _SUNMatrixContent_MagmaDense *SUNMatrixContent_MagmaDense;

/* ---------------------------------------
 * Implementation specific functions
 * ---------------------------------------*/

SUNDIALS_EXPORT SUNMatrix SUNMatrix_MagmaDense(sunindextype M, sunindextype N, SUNMemoryType memtype,
                                               SUNMemoryHelper memhelper, void* queue, SUNContext sunctx);
SUNDIALS_EXPORT SUNMatrix SUNMatrix_MagmaDenseBlock(sunindextype nblocks, sunindextype M, sunindextype N,
                                                    SUNMemoryType memtype, SUNMemoryHelper memhelper,
                                                    void* queue, SUNContext sunctx);
SUNDIALS_EXPORT void SUNMatrix_MagmaDense_Print(SUNMatrix A);
SUNDIALS_EXPORT realtype* SUNMatrix_MagmaDense_Data(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_LData(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_Columns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_BlockRows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_BlockColumns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNMatrix_MagmaDense_NumBlocks(SUNMatrix A);
SUNDIALS_EXPORT realtype** SUNMatrix_MagmaDense_BlockData(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatrix_MagmaDense_CopyToDevice(SUNMatrix A, realtype* h_data);
SUNDIALS_EXPORT int SUNMatrix_MagmaDense_CopyFromDevice(SUNMatrix A, realtype* h_data);

SUNDIALS_STATIC_INLINE
realtype* SUNMatrix_MagmaDense_Block(SUNMatrix Amat, sunindextype k)
{
  SUNMatrixContent_MagmaDense A = (SUNMatrixContent_MagmaDense) Amat->content;
  return( ((realtype*) A->data->ptr) + k*A->M*A->N );
}

SUNDIALS_STATIC_INLINE
realtype* SUNMatrix_MagmaDense_Column(SUNMatrix Amat, sunindextype j)
{
  SUNMatrixContent_MagmaDense A = (SUNMatrixContent_MagmaDense) Amat->content;
  return( ((realtype*) A->data->ptr) + j*A->M );
}

SUNDIALS_STATIC_INLINE
realtype* SUNMatrix_MagmaDense_BlockColumn(SUNMatrix Amat, sunindextype k, sunindextype j)
{
  SUNMatrixContent_MagmaDense A = (SUNMatrixContent_MagmaDense) Amat->content;
  return( ((realtype*) A->data->ptr) + k*A->M*A->N + j*A->M );
}


/* ---------------------------------------
 * SUNMatrix API functions
 * ---------------------------------------*/

SUNDIALS_STATIC_INLINE
SUNMatrix_ID SUNMatGetID_MagmaDense(SUNMatrix A) { return SUNMATRIX_MAGMADENSE; }

SUNDIALS_EXPORT SUNMatrix SUNMatClone_MagmaDense(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_MagmaDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_MagmaDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_MagmaDense(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd_MagmaDense(realtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_MagmaDense(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvecSetup_MagmaDense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec_MagmaDense(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace_MagmaDense(SUNMatrix A, long int *lenrw, long int *leniw);


#ifdef __cplusplus
}
#endif

#endif
