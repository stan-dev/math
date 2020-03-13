/* ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
 * Based on work by Donald Wilcox @ LBNL
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
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

#include <nvector/nvector_cuda.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>
#include <sundials/sundials_math.h>

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
#define SUN_CUSP_HANDLE(S)         ( SUN_CUSP_CONTENT(S)->cusolver_handle )
#define SUN_CUSP_SUBSYS_SIZE(S)    ( SUN_CUSP_CONTENT(S)->subsys_size )
#define SUN_CUSP_SUBSYS_NNZ(S)     ( SUN_CUSP_CONTENT(S)->subsys_nnz )
#define SUN_CUSP_MATDESC(S)        ( SUN_CUSP_CONTENT(S)->system_description )
#define SUN_CUSP_NUM_SUBSYS(S)     ( SUN_CUSP_CONTENT(S)->nsubsys )
#define SUN_CUSP_DROWPTR(S)        ( SUN_CUSP_CONTENT(S)->d_rowptr )
#define SUN_CUSP_DCOLIND(S)        ( SUN_CUSP_CONTENT(S)->d_colind )
#define SUN_CUSP_DVALUES(S)        ( SUN_CUSP_CONTENT(S)->d_values )
#define SUN_CUSP_DESC(S)           ( SUN_CUSP_CONTENT(S)->desc )
#define SUN_CUSP_QRINFO(S)         ( SUN_CUSP_CONTENT(S)->info )
#define SUN_CUSP_INTERNAL_SIZE(S)  ( SUN_CUSP_CONTENT(S)->internal_size )
#define SUN_CUSP_WORK_SIZE(S)      ( SUN_CUSP_CONTENT(S)->workspace_size )

/*
 * ----------------------------------------------------------------------------
 *  Implementations of exported functions.
 * ----------------------------------------------------------------------------
 */

SUNLinearSolver SUNLinSol_cuSolverSp_batchQR(N_Vector y, SUNMatrix A, int nsubsys,
                                             int subsys_size, int subsys_nnz)
{
  /* Check that required arguments are not NULL */
  if (y == NULL || A == NULL) return(NULL);

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE ||
      N_VGetVectorID(y) != SUNDIALS_NVEC_CUDA) return(NULL);

  /* Check that it is a CSR matrix */
  if (SUNSparseMatrix_SparseType(A) != CSR_MAT) return(NULL);

  /* Check that the vector is using managed memory */
  if (!N_VIsManagedMemory_Cuda(y)) return(NULL);

  /* Matrix must be square */
  if (SUNSparseMatrix_Columns(A) != SUNSparseMatrix_Rows(A)) return(NULL);

  /* Matrix and vector dimensions must agree */
  if (N_VGetLength(y) != SUNSparseMatrix_Columns(A)) return(NULL);

  /* All subsystems must be the same size */
  if (SUNSparseMatrix_Columns(A) != (subsys_size * nsubsys)) return(NULL);

  /* Number of nonzeros per subsys must be the same */
  if (SUNSparseMatrix_NNZ(A) != (subsys_nnz * nsubsys)) return(NULL);

  /* Allocate device memory for the matrix */
  int *d_colind, *d_rowptr;
  realtype *d_values;

  d_colind = NULL;
  d_rowptr = NULL;
  d_values = NULL;

  cudaError_t cuerr;
  cuerr = cudaMalloc((void **) &d_colind, sizeof(*d_colind) * subsys_nnz);
  if (cuerr != cudaSuccess) return(NULL);
  cuerr = cudaMalloc((void **) &d_rowptr, sizeof(*d_rowptr) * (subsys_size + 1));
  if (cuerr != cudaSuccess) { cudaFree(d_colind); return(NULL); }
  cuerr = cudaMalloc((void **) &d_values, sizeof(*d_values) * subsys_nnz * nsubsys);
  if (cuerr != cudaSuccess) { cudaFree(d_rowptr); cudaFree(d_colind); return(NULL); }

  /* Create an empty linear solver */
  SUNLinearSolver S;

  S = NULL;
  S = SUNLinSolNewEmpty();
  if (S == NULL) {
    cudaFree(d_rowptr); cudaFree(d_colind); cudaFree(d_values);
    return(NULL);
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
  content = (SUNLinearSolverContent_cuSolverSp_batchQR) malloc(sizeof *content);
  if (content == NULL) {
    cudaFree(d_rowptr); cudaFree(d_colind); cudaFree(d_values);
    SUNLinSolFree(S);
    return(NULL);
  }

  /* Attach content */
  S->content = content;

  /* Fill content */
  cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
  cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;

  cusolver_status = cusolverSpCreate(&content->cusolver_handle);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNLinSolFree(S);
    return(NULL);
  }

  cusparse_status = cusparseCreateMatDescr(&content->system_description);
  if (cusparse_status != CUSPARSE_STATUS_SUCCESS) {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNLinSolFree(S);
    return(NULL);
  }

  cusparse_status = cusparseSetMatType(content->system_description, CUSPARSE_MATRIX_TYPE_GENERAL);
  if (cusparse_status != CUSPARSE_STATUS_SUCCESS) {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNLinSolFree(S);
    return(NULL);
  }

  cusparse_status = cusparseSetMatIndexBase(content->system_description, CUSPARSE_INDEX_BASE_ZERO);
  if (cusparse_status != CUSPARSE_STATUS_SUCCESS) {
    cudaFree(d_rowptr); cudaFree(d_colind);
    cudaFree(d_values); SUNLinSolFree(S);
    return(NULL);
  }

  content->info        = NULL;
  content->workspace   = NULL;
  content->subsys_size = subsys_size;
  content->subsys_nnz  = subsys_nnz;
  content->nsubsys     = nsubsys;
  content->d_colind    = d_colind;
  content->d_rowptr    = d_rowptr;
  content->d_values    = d_values;
  content->desc        = NULL;

  return(S);
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
  cudaError_t cuerr;
  cusolverStatus_t status;

  /* copy matrix to the device */
  cuerr = cudaMemcpy(SUN_CUSP_DCOLIND(S), SUNSparseMatrix_IndexValues(A),
                     sizeof(int) * SUN_CUSP_SUBSYS_NNZ(S), cudaMemcpyHostToDevice);
  if (cuerr != cudaSuccess) SUN_CUSP_LASTFLAG(S) = SUNLS_MEM_FAIL;

  cuerr = cudaMemcpy(SUN_CUSP_DROWPTR(S), SUNSparseMatrix_IndexPointers(A),
                     sizeof(int) * (SUN_CUSP_SUBSYS_SIZE(S)+1), cudaMemcpyHostToDevice);
  if (cuerr != cudaSuccess) SUN_CUSP_LASTFLAG(S) = SUNLS_MEM_FAIL;

  cuerr = cudaMemcpy(SUN_CUSP_DVALUES(S), SUNSparseMatrix_Data(A),
                     sizeof(realtype) * SUN_CUSP_SUBSYS_NNZ(S) * SUN_CUSP_NUM_SUBSYS(S),
                     cudaMemcpyHostToDevice);
  if (cuerr != cudaSuccess) SUN_CUSP_LASTFLAG(S) = SUNLS_MEM_FAIL;

  if (SUN_CUSP_LASTFLAG(S) != SUNLS_SUCCESS)
    return(SUN_CUSP_LASTFLAG(S));

  if (SUN_CUSP_FIRSTFACTORIZE(S)) {

    /* Free old workspace and symbloic analysis */
    if (SUN_CUSP_QRWORKSPACE(S)) {
      cudaFree(SUN_CUSP_QRWORKSPACE(S));
      cusolverSpDestroyCsrqrInfo(SUN_CUSP_QRINFO(S));
    }

    /* We must create a new csrqrinfo_t context every time we want to
       do a symbolic analysis. Trying to reuse it results in a
       CUSOLVER_STATUS_INVALID_VALUE error. */
    status = cusolverSpCreateCsrqrInfo(&SUN_CUSP_QRINFO(S));
    if (status != CUSOLVER_STATUS_SUCCESS) {
      SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return(SUN_CUSP_LASTFLAG(S));
    }

    /* Perform symbolic analysis of sparsity structure */
    status = cusolverSpXcsrqrAnalysisBatched(SUN_CUSP_HANDLE(S),
                                             SUN_CUSP_SUBSYS_SIZE(S),
                                             SUN_CUSP_SUBSYS_SIZE(S),
                                             SUN_CUSP_SUBSYS_NNZ(S),
                                             SUN_CUSP_MATDESC(S),
                                             SUN_CUSP_DROWPTR(S),
                                             SUN_CUSP_DCOLIND(S),
                                             SUN_CUSP_QRINFO(S));

    if (status != CUSOLVER_STATUS_SUCCESS) {
      SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return(SUN_CUSP_LASTFLAG(S));
    }

    /* Compute the workspace we will need */
    status = _cusolverSpXcsrqrBufferInfoBatched(SUN_CUSP_HANDLE(S),
                                                SUN_CUSP_SUBSYS_SIZE(S),
                                                SUN_CUSP_SUBSYS_SIZE(S),
                                                SUN_CUSP_SUBSYS_NNZ(S),
                                                SUN_CUSP_MATDESC(S),
                                                SUN_CUSP_DVALUES(S),
                                                SUN_CUSP_DROWPTR(S),
                                                SUN_CUSP_DCOLIND(S),
                                                SUN_CUSP_NUM_SUBSYS(S),
                                                SUN_CUSP_QRINFO(S),
                                                &SUN_CUSP_INTERNAL_SIZE(S),
                                                &SUN_CUSP_WORK_SIZE(S));

    if (status != CUSOLVER_STATUS_SUCCESS) {
      SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return(SUN_CUSP_LASTFLAG(S));
    }

    cuerr = cudaMalloc((void**) &SUN_CUSP_QRWORKSPACE(S), SUN_CUSP_WORK_SIZE(S));
    if (cuerr != cudaSuccess) {
      SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return(SUN_CUSP_LASTFLAG(S));
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

  if ((S == NULL) || (A == NULL) || (x == NULL) || (b == NULL))
    return(SUNLS_MEM_NULL);

  SUN_CUSP_LASTFLAG(S) = SUNLS_SUCCESS;

  realtype* device_b = N_VGetDeviceArrayPointer_Cuda(b);
  realtype* device_x = N_VGetDeviceArrayPointer_Cuda(x);

  if (SUN_CUSP_LASTFLAG(S) != SUNLS_SUCCESS)
    return(SUN_CUSP_LASTFLAG(S));

  /* solve the system */
  status = _cusolverSpXcsrqrsvBatched(SUN_CUSP_HANDLE(S),
                                      SUN_CUSP_SUBSYS_SIZE(S),
                                      SUN_CUSP_SUBSYS_SIZE(S),
                                      SUN_CUSP_SUBSYS_NNZ(S),
                                      SUN_CUSP_MATDESC(S),
                                      SUN_CUSP_DVALUES(S),
                                      SUN_CUSP_DROWPTR(S),
                                      SUN_CUSP_DCOLIND(S),
                                      device_b,
                                      device_x,
                                      SUN_CUSP_NUM_SUBSYS(S),
                                      SUN_CUSP_QRINFO(S),
                                      SUN_CUSP_QRWORKSPACE(S));

  if (status != CUSOLVER_STATUS_SUCCESS) {
    SUN_CUSP_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
    return(SUN_CUSP_LASTFLAG(S));
  }

  return(SUN_CUSP_LASTFLAG(S));
}

sunindextype SUNLinSolLastFlag_cuSolverSp_batchQR(SUNLinearSolver S)
{
  if (S == NULL) return(-1);
  return(SUN_CUSP_LASTFLAG(S));
}

int SUNLinSolFree_cuSolverSp_batchQR(SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == NULL) return(SUNLS_SUCCESS);

  /* free stuff in the content structure */
  cusolverSpDestroy(SUN_CUSP_HANDLE(S));
  cusolverSpDestroyCsrqrInfo(SUN_CUSP_QRINFO(S));
  cudaFree(SUN_CUSP_DCOLIND(S));
  cudaFree(SUN_CUSP_DROWPTR(S));
  cudaFree(SUN_CUSP_DVALUES(S));
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
