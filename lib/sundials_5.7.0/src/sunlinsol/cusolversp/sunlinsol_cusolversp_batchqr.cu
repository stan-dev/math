/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
 * Based on work by Donald Wilcox @ LBNL
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * Implementation file for cuSolverSp batched QR SUNLinearSolver interface.
 * ----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunmatrix/sunmatrix_cusparse.h>
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>

#include "sundials_cuda.h"
#include "sundials_debug.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/* macros for handling the different function names based on precision */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define _cusolverSpXcsrqrBufferInfoBatched cusolverSpDcsrqrBufferInfoBatched
#define _cusolverSpXcsrqrsvBatched cusolverSpDcsrqrsvBatched
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define _cusolverSpXcsrqrBufferInfoBatched cusolverSpScsrqrBufferInfoBatched
#define _cusolverSpXcsrqrsvBatched cusolverSpScsrqrsvBatched
#endif

/*
 * -----------------------------------------------------------------
 * cuSolverSp solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SUN_CUSP_CONTENT(S)        ( (SUNLinearSolverContent_cuSolverSp_batchQR)(S->content) )
#define SUN_CUSP_QRWORKSPACE(S)    ( SUN_CUSP_CONTENT(S)->workspace )
#define SUN_CUSP_FIRSTFACTORIZE(S) ( SUN_CUSP_CONTENT(S)->first_factorize )
#define SUN_CUSP_LASTFLAG(S)       ( SUN_CUSP_CONTENT(S)->last_flag )
#define SUN_CUSOL_HANDLE(S)        ( SUN_CUSP_CONTENT(S)->cusolver_handle )
#define SUN_CUSP_DESC(S)           ( SUN_CUSP_CONTENT(S)->desc )
#define SUN_CUSP_QRINFO(S)         ( SUN_CUSP_CONTENT(S)->info )
#define SUN_CUSP_INTERNAL_SIZE(S)  ( SUN_CUSP_CONTENT(S)->internal_size )
#define SUN_CUSP_WORK_SIZE(S)      ( SUN_CUSP_CONTENT(S)->workspace_size )

/*
 * ----------------------------------------------------------------------------
 *  Implementations of exported functions.
 * ----------------------------------------------------------------------------
 */

SUNLinearSolver SUNLinSol_cuSolverSp_batchQR(N_Vector y, SUNMatrix A, cusolverSpHandle_t cusol_handle)
{
  /* Check that required arguments are not NULL */
  if (y == NULL || A == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNLinSol_cuSolverSp_batchQR: y or A is null\n");
    return NULL;
  }

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_CUSPARSE || y->ops->nvgetdevicearraypointer == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNLinSol_cuSolverSp_batchQR: illegal type for y or A\n");
    return NULL;
  }

  /* Matrix and vector dimensions must agree */
  if (N_VGetLength(y) != SUNMatrix_cuSparse_Columns(A))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNLinSol_cuSolverSp_batchQR: matrix and vector dimensions don't agree\n");
    return NULL;
  }

  /* Create an empty linear solver */
  SUNLinearSolver S;

  S = NULL;
  S = SUNLinSolNewEmpty();
  if (S == NULL)
  {
    return NULL;
  }

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_cuSolverSp_batchQR;
  S->ops->getid      = SUNLinSolGetID_cuSolverSp_batchQR;
  S->ops->initialize = SUNLinSolInitialize_cuSolverSp_batchQR;
  S->ops->setup      = SUNLinSolSetup_cuSolverSp_batchQR;
  S->ops->solve      = SUNLinSolSolve_cuSolverSp_batchQR;
  S->ops->lastflag   = SUNLinSolLastFlag_cuSolverSp_batchQR;
  S->ops->free       = SUNLinSolFree_cuSolverSp_batchQR;

  /* Create content */
  SUNLinearSolverContent_cuSolverSp_batchQR content;

  content = NULL;
  content = (SUNLinearSolverContent_cuSolverSp_batchQR) malloc(sizeof(*content));
  if (content == NULL)
  {
    SUNLinSolFree(S);
    return NULL;
  }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag       = SUNLS_SUCCESS;
  content->first_factorize = SUNTRUE;
  content->internal_size   = 0;
  content->workspace_size  = 0;
  content->cusolver_handle = cusol_handle;
  content->info            = NULL;
  content->workspace       = NULL;
  content->desc            = NULL;

  return S;
}

/*
 * -----------------------------------------------------------------
 * Implementation of accessor and setter functions.
 * -----------------------------------------------------------------
 */

void SUNLinSol_cuSolverSp_batchQR_GetDescription(SUNLinearSolver S, const char** desc)
{
  *desc = SUN_CUSP_DESC(S);
}

void SUNLinSol_cuSolverSp_batchQR_SetDescription(SUNLinearSolver S, const char* desc)
{
  SUN_CUSP_DESC(S) = desc;
}

void SUNLinSol_cuSolverSp_batchQR_GetDeviceSpace(SUNLinearSolver S,
                                                 size_t* cuSolverInternal,
                                                 size_t* cuSolverWorkspace)
{
  /* size is in bytes */
  *cuSolverInternal  = SUN_CUSP_INTERNAL_SIZE(S); /* buffer for Q and R factors */
  *cuSolverWorkspace = SUN_CUSP_WORK_SIZE(S); /* numerical factorization buffer */
}

/*
 * -----------------------------------------------------------------
 * Implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_cuSolverSp_batchQR(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_cuSolverSp_batchQR(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_CUSOLVERSP_BATCHQR);
}

int SUNLinSolInitialize_cuSolverSp_batchQR(SUNLinearSolver S)
{
  SUN_CUSP_FIRSTFACTORIZE(S) = SUNTRUE;
  SUN_CUSP_LASTFLAG(S) = SUNLS_SUCCESS;
  return(SUN_CUSP_LASTFLAG(S));
}

int SUNLinSolSetup_cuSolverSp_batchQR(SUNLinearSolver S, SUNMatrix A)
{
  int blockrows, blockcols, blocknnz, nblock;
  int *d_rowptr, *d_colind;
  realtype *d_data;
  cusparseMatDescr_t mat_descr;
  cudaError_t cuerr;
  cusolverStatus_t status;

  if (SUN_CUSP_LASTFLAG(S) != SUNLS_SUCCESS)
    return SUN_CUSP_LASTFLAG(S);

  if (SUN_CUSP_FIRSTFACTORIZE(S))
  {

    /* Free old workspace and symbloic analysis */
    if (SUN_CUSP_QRWORKSPACE(S))
    {
      cudaFree(SUN_CUSP_QRWORKSPACE(S));
      cusolverSpDestroyCsrqrInfo(SUN_CUSP_QRINFO(S));
    }

    /* We must create a new csrqrinfo_t context every time we want to
       do a symbolic analysis. Trying to reuse it results in a
       CUSOLVER_STATUS_INVALID_VALUE error. */
    status = cusolverSpCreateCsrqrInfo(&SUN_CUSP_QRINFO(S));
    if (!SUNDIALS_CUSOLVER_VERIFY(status))
    {
      SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return SUN_CUSP_LASTFLAG(S);
    }

    nblock    = SUNMatrix_cuSparse_NumBlocks(A);
    blocknnz  = SUNMatrix_cuSparse_BlockNNZ(A);
    blockrows = SUNMatrix_cuSparse_BlockRows(A);
    blockcols = SUNMatrix_cuSparse_BlockColumns(A);
    d_data    = SUNMatrix_cuSparse_Data(A);
    d_rowptr  = SUNMatrix_cuSparse_IndexPointers(A);
    d_colind  = SUNMatrix_cuSparse_IndexValues(A);
    mat_descr = SUNMatrix_cuSparse_MatDescr(A);

    /* Perform symbolic analysis of sparsity structure */
    status = cusolverSpXcsrqrAnalysisBatched(SUN_CUSOL_HANDLE(S),
                                             blockrows,
                                             blockcols,
                                             blocknnz,
                                             mat_descr,
                                             d_rowptr,
                                             d_colind,
                                             SUN_CUSP_QRINFO(S));

    if (!SUNDIALS_CUSOLVER_VERIFY(status))
    {
      SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return SUN_CUSP_LASTFLAG(S);
    }

    /* Compute the workspace we will need */
    status = _cusolverSpXcsrqrBufferInfoBatched(SUN_CUSOL_HANDLE(S),
                                                blockrows,
                                                blockcols,
                                                blocknnz,
                                                mat_descr,
                                                d_data,
                                                d_rowptr,
                                                d_colind,
                                                nblock,
                                                SUN_CUSP_QRINFO(S),
                                                &SUN_CUSP_INTERNAL_SIZE(S),
                                                &SUN_CUSP_WORK_SIZE(S));

    if (!SUNDIALS_CUSOLVER_VERIFY(status))
    {
      SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return SUN_CUSP_LASTFLAG(S);
    }

    cuerr = cudaMalloc((void**) &SUN_CUSP_QRWORKSPACE(S), SUN_CUSP_WORK_SIZE(S));
    if (!SUNDIALS_CUDA_VERIFY(cuerr))
    {
      SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return SUN_CUSP_LASTFLAG(S);
    }

    SUN_CUSP_FIRSTFACTORIZE(S) = SUNFALSE;
  }

  SUN_CUSP_LASTFLAG(S) = SUNLS_SUCCESS;
  return(SUN_CUSP_LASTFLAG(S));
}

int SUNLinSolSolve_cuSolverSp_batchQR(SUNLinearSolver S, SUNMatrix A,
                                      N_Vector x, N_Vector b, realtype tol)
{
  cusolverStatus_t status;
  int blockrows, blockcols, blocknnz, nblock;
  int *d_rowptr, *d_colind;
  realtype *d_data;
  cusparseMatDescr_t mat_descr;

  if ((S == NULL) || (A == NULL) || (x == NULL) || (b == NULL))
    return SUNLS_MEM_NULL;

  SUN_CUSP_LASTFLAG(S) = SUNLS_SUCCESS;

  realtype* device_b = N_VGetDeviceArrayPointer(b);
  realtype* device_x = N_VGetDeviceArrayPointer(x);

  if (SUN_CUSP_LASTFLAG(S) != SUNLS_SUCCESS)
    return SUN_CUSP_LASTFLAG(S);

  /* solve the system */
  nblock    = SUNMatrix_cuSparse_NumBlocks(A);
  blocknnz  = SUNMatrix_cuSparse_BlockNNZ(A);
  blockrows = SUNMatrix_cuSparse_BlockRows(A);
  blockcols = SUNMatrix_cuSparse_BlockColumns(A);
  d_data    = SUNMatrix_cuSparse_Data(A);
  d_rowptr  = SUNMatrix_cuSparse_IndexPointers(A);
  d_colind  = SUNMatrix_cuSparse_IndexValues(A);
  mat_descr = SUNMatrix_cuSparse_MatDescr(A);

  status = _cusolverSpXcsrqrsvBatched(SUN_CUSOL_HANDLE(S),
                                      blockrows,
                                      blockcols,
                                      blocknnz,
                                      mat_descr,
                                      d_data,
                                      d_rowptr,
                                      d_colind,
                                      device_b,
                                      device_x,
                                      nblock,
                                      SUN_CUSP_QRINFO(S),
                                      SUN_CUSP_QRWORKSPACE(S));

  if (!SUNDIALS_CUSOLVER_VERIFY(status))
  {
    SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
    return SUN_CUSP_LASTFLAG(S);
  }

  return SUN_CUSP_LASTFLAG(S);
}

sunindextype SUNLinSolLastFlag_cuSolverSp_batchQR(SUNLinearSolver S)
{
  if (S == NULL) return -1;
  return SUN_CUSP_LASTFLAG(S);
}

int SUNLinSolFree_cuSolverSp_batchQR(SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == NULL) return SUNLS_SUCCESS;

  /* free stuff in the content structure */
  cusolverSpDestroyCsrqrInfo(SUN_CUSP_QRINFO(S));
  cudaFree(SUN_CUSP_QRWORKSPACE(S));

  /* free content structure */
  if (S->content) {
    free(S->content);
    S->content = NULL;
  }

  /* free ops structure */
  if (S->ops) {
    free(S->ops);
    S->ops = NULL;
  }

  /* free the actual SUNLinSol */
  free(S);
  S = NULL;

  return(SUNLS_SUCCESS);
}
