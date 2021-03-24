/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on codes <solver>_lapack.c by: Radu Serban @ LLNL
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
 * This is the implementation file for the LAPACK dense
 * implementation of the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_lapackdense.h>
#include <sundials/sundials_math.h>

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * LapackDense solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define LAPACKDENSE_CONTENT(S) ( (SUNLinearSolverContent_LapackDense)(S->content) )
#define PIVOTS(S)              ( LAPACKDENSE_CONTENT(S)->pivots )
#define LASTFLAG(S)            ( LAPACKDENSE_CONTENT(S)->last_flag )

/*
 * -----------------------------------------------------------------
 * deprecated wrapper functions
 * -----------------------------------------------------------------
 */

SUNLinearSolver SUNLapackDense(N_Vector y, SUNMatrix A)
{ return(SUNLinSol_LapackDense(y, A)); }

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new LAPACK dense linear solver
 */

SUNLinearSolver SUNLinSol_LapackDense(N_Vector y, SUNMatrix A)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_LapackDense content;
  sunindextype MatrixRows;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE) return(NULL);

  if (SUNDenseMatrix_Rows(A) != SUNDenseMatrix_Columns(A)) return(NULL);

  if ( (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS) )
    return(NULL);

  MatrixRows = SUNDenseMatrix_Rows(A);
  if (MatrixRows != N_VGetLength(y)) return(NULL);

  /* Create linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty();
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_LapackDense;
  S->ops->getid      = SUNLinSolGetID_LapackDense;
  S->ops->initialize = SUNLinSolInitialize_LapackDense;
  S->ops->setup      = SUNLinSolSetup_LapackDense;
  S->ops->solve      = SUNLinSolSolve_LapackDense;
  S->ops->lastflag   = SUNLinSolLastFlag_LapackDense;
  S->ops->space      = SUNLinSolSpace_LapackDense;
  S->ops->free       = SUNLinSolFree_LapackDense;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_LapackDense) malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->N         = MatrixRows;
  content->last_flag = 0;
  content->pivots    = NULL;

  /* Allocate content */
  content->pivots = (sunindextype *) malloc(MatrixRows * sizeof(sunindextype));
  if (content->pivots == NULL) { SUNLinSolFree(S); return(NULL); }

  return(S);
}


/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_LapackDense(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}


SUNLinearSolver_ID SUNLinSolGetID_LapackDense(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_LAPACKDENSE);
}


int SUNLinSolInitialize_LapackDense(SUNLinearSolver S)
{
  /* all solver-specific memory has already been allocated */
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(SUNLS_SUCCESS);
}


int SUNLinSolSetup_LapackDense(SUNLinearSolver S, SUNMatrix A)
{
  sunindextype n, ier;

  /* check for valid inputs */
  if ( (A == NULL) || (S == NULL) )
    return(SUNLS_MEM_NULL);

  /* Ensure that A is a dense matrix */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE) {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(SUNLS_ILL_INPUT);
  }

  /* Call LAPACK to do LU factorization of A */
  n = SUNDenseMatrix_Rows(A);
  xgetrf_f77(&n, &n, SUNDenseMatrix_Data(A), &n, PIVOTS(S), &ier);
  LASTFLAG(S) = ier;
  if (ier > 0)
    return(SUNLS_LUFACT_FAIL);
  if (ier < 0)
    return(SUNLS_PACKAGE_FAIL_UNREC);
  return(SUNLS_SUCCESS);
}


int SUNLinSolSolve_LapackDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, realtype tol)
{
  sunindextype n, one, ier;
  realtype *xdata;

  if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) )
    return(SUNLS_MEM_NULL);

  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access x data array */
  xdata = N_VGetArrayPointer(x);
  if (xdata == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(SUNLS_MEM_FAIL);
  }

  /* Call LAPACK to solve the linear system */
  n = SUNDenseMatrix_Rows(A);
  one = 1;
  xgetrs_f77("N", &n, &one, SUNDenseMatrix_Data(A),
             &n, PIVOTS(S), xdata, &n, &ier);
  LASTFLAG(S) = ier;
  if (ier < 0)
    return(SUNLS_PACKAGE_FAIL_UNREC);

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(SUNLS_SUCCESS);
}


sunindextype SUNLinSolLastFlag_LapackDense(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return(-1);
  return(LASTFLAG(S));
}


int SUNLinSolSpace_LapackDense(SUNLinearSolver S,
                               long int *lenrwLS,
                               long int *leniwLS)
{
  *lenrwLS = 0;
  *leniwLS = 2 + LAPACKDENSE_CONTENT(S)->N;
  return(SUNLS_SUCCESS);
}

int SUNLinSolFree_LapackDense(SUNLinearSolver S)
{
  /* return if S is already free */
  if (S == NULL) return(SUNLS_SUCCESS);

  /* delete items from contents, then delete generic structure */
  if (S->content) {
    if (PIVOTS(S)) {
      free(PIVOTS(S));
      PIVOTS(S) = NULL;
    }
    free(S->content);
    S->content = NULL;
  }
  if (S->ops) {
    free(S->ops);
    S->ops = NULL;
  }
  free(S); S = NULL;
  return(SUNLS_SUCCESS);
}
