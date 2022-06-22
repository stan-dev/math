/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on codes <solver>_klu.c, written by Carol Woodward @ LLNL
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
 * This is the implementation file for the KLU implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_klu.h>
#include <sundials/sundials_math.h>

#define ZERO      RCONST(0.0)
#define ONE       RCONST(1.0)
#define TWO       RCONST(2.0)
#define TWOTHIRDS RCONST(0.666666666666666666666666666666667)

/*
 * -----------------------------------------------------------------
 * KLU solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define KLU_CONTENT(S)     ( (SUNLinearSolverContent_KLU)(S->content) )
#define LASTFLAG(S)        ( KLU_CONTENT(S)->last_flag )
#define FIRSTFACTORIZE(S)  ( KLU_CONTENT(S)->first_factorize )
#define SYMBOLIC(S)        ( KLU_CONTENT(S)->symbolic )
#define NUMERIC(S)         ( KLU_CONTENT(S)->numeric )
#define COMMON(S)          ( KLU_CONTENT(S)->common )
#define SOLVE(S)           ( KLU_CONTENT(S)->klu_solver )

/*
 * -----------------------------------------------------------------
 * typedef to handle pointer casts from sunindextype to KLU type
 * -----------------------------------------------------------------
 */

#if defined(SUNDIALS_INT64_T)
#define KLU_INDEXTYPE long int
#else
#define KLU_INDEXTYPE int
#endif

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new KLU linear solver
 */

SUNLinearSolver SUNLinSol_KLU(N_Vector y, SUNMatrix A, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_KLU content;
  int flag;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) return(NULL);

  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A)) return(NULL);

  if ( (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS) )
    return(NULL);

  if (SUNSparseMatrix_Rows(A) != N_VGetLength(y)) return(NULL);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_KLU;
  S->ops->getid      = SUNLinSolGetID_KLU;
  S->ops->initialize = SUNLinSolInitialize_KLU;
  S->ops->setup      = SUNLinSolSetup_KLU;
  S->ops->solve      = SUNLinSolSolve_KLU;
  S->ops->lastflag   = SUNLinSolLastFlag_KLU;
  S->ops->space      = SUNLinSolSpace_KLU;
  S->ops->free       = SUNLinSolFree_KLU;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_KLU) malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag       = 0;
  content->first_factorize = 1;
  content->symbolic        = NULL;
  content->numeric         = NULL;

#if defined(SUNDIALS_INT64_T)
  if (SUNSparseMatrix_SparseType(A) == CSC_MAT) {
    content->klu_solver = (KLUSolveFn) &klu_l_solve;
  } else {
    content->klu_solver = (KLUSolveFn) &klu_l_tsolve;
  }
#elif defined(SUNDIALS_INT32_T)
  if (SUNSparseMatrix_SparseType(A) == CSC_MAT) {
    content->klu_solver = &klu_solve;
  } else {
    content->klu_solver = &klu_tsolve;
  }
#else  /* incompatible sunindextype for KLU */
#error  Incompatible sunindextype for KLU
#endif

  flag = sun_klu_defaults(&(content->common));
  if (flag == 0) { SUNLinSolFree(S); return(NULL); }
  (content->common).ordering = SUNKLU_ORDERING_DEFAULT;

  return(S);
}


/* ----------------------------------------------------------------------------
 * Function to reinitialize a KLU linear solver
 */

int SUNLinSol_KLUReInit(SUNLinearSolver S, SUNMatrix A,
                        sunindextype nnz, int reinit_type)
{
  /* Check for non-NULL SUNLinearSolver */
  if ((S == NULL) || (A == NULL))
    return(SUNLS_MEM_NULL);

  /* Check for valid SUNMatrix */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE)
    return(SUNLS_ILL_INPUT);

  /* Check for valid reinit_type */
  if ((reinit_type != SUNKLU_REINIT_FULL) &&
      (reinit_type != SUNKLU_REINIT_PARTIAL))
    return(SUNLS_ILL_INPUT);

  /* Full re-initialization: reallocate matrix for updated storage */
  if (reinit_type == SUNKLU_REINIT_FULL)
    if (SUNSparseMatrix_Reallocate(A, nnz) != 0)
      return(SUNLS_MEM_FAIL);

  /* Free the prior factorazation and reset for first factorization */
  if( SYMBOLIC(S) != NULL)
    sun_klu_free_symbolic(&SYMBOLIC(S), &COMMON(S));
  if( NUMERIC(S) != NULL)
    sun_klu_free_numeric(&NUMERIC(S), &COMMON(S));
  FIRSTFACTORIZE(S) = 1;

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

/* ----------------------------------------------------------------------------
 * Function to set the ordering type for a KLU linear solver
 */

int SUNLinSol_KLUSetOrdering(SUNLinearSolver S, int ordering_choice)
{
  /* Check for legal ordering_choice */
  if ((ordering_choice < 0) || (ordering_choice > 2))
    return(SUNLS_ILL_INPUT);

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SUNLS_MEM_NULL);

  /* Set ordering_choice */
  COMMON(S).ordering = ordering_choice;

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


/*
 * -----------------------------------------------------------------
 * accessor functions
 * -----------------------------------------------------------------
 */

sun_klu_symbolic* SUNLinSol_KLUGetSymbolic(SUNLinearSolver S)
{
  return(SYMBOLIC(S));
}

sun_klu_numeric* SUNLinSol_KLUGetNumeric(SUNLinearSolver S)
{
  return(NUMERIC(S));
}

sun_klu_common* SUNLinSol_KLUGetCommon(SUNLinearSolver S)
{
  return(&(COMMON(S)));
}


/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_KLU(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}


SUNLinearSolver_ID SUNLinSolGetID_KLU(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_KLU);
}


int SUNLinSolInitialize_KLU(SUNLinearSolver S)
{
  /* Force factorization */
  FIRSTFACTORIZE(S) = 1;

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSetup_KLU(SUNLinearSolver S, SUNMatrix A)
{
  int retval;
  realtype uround_twothirds;

  uround_twothirds = SUNRpowerR(UNIT_ROUNDOFF,TWOTHIRDS);

  /* Ensure that A is a sparse matrix */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(LASTFLAG(S));
  }

  /* On first decomposition, get the symbolic factorization */
  if (FIRSTFACTORIZE(S)) {

    /* Perform symbolic analysis of sparsity structure */
    if (SYMBOLIC(S))
      sun_klu_free_symbolic(&SYMBOLIC(S), &COMMON(S));
    SYMBOLIC(S) = sun_klu_analyze(SUNSparseMatrix_NP(A),
                                  (KLU_INDEXTYPE*) SUNSparseMatrix_IndexPointers(A),
                                  (KLU_INDEXTYPE*) SUNSparseMatrix_IndexValues(A),
                                  &COMMON(S));
    if (SYMBOLIC(S) == NULL) {
      LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return(LASTFLAG(S));
    }

    /* ------------------------------------------------------------
       Compute the LU factorization of the matrix
       ------------------------------------------------------------*/
    if(NUMERIC(S))
      sun_klu_free_numeric(&NUMERIC(S), &COMMON(S));
    NUMERIC(S) = sun_klu_factor((KLU_INDEXTYPE*) SUNSparseMatrix_IndexPointers(A),
                                (KLU_INDEXTYPE*) SUNSparseMatrix_IndexValues(A),
                                SUNSparseMatrix_Data(A),
                                SYMBOLIC(S),
                                &COMMON(S));
    if (NUMERIC(S) == NULL) {
      LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return(LASTFLAG(S));
    }

    FIRSTFACTORIZE(S) = 0;

  } else {   /* not the first decomposition, so just refactor */

    retval = sun_klu_refactor((KLU_INDEXTYPE*) SUNSparseMatrix_IndexPointers(A),
                              (KLU_INDEXTYPE*) SUNSparseMatrix_IndexValues(A),
                              SUNSparseMatrix_Data(A),
                              SYMBOLIC(S),
                              NUMERIC(S),
                              &COMMON(S));
    if (retval == 0) {
      LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
      return(LASTFLAG(S));
    }

    /*-----------------------------------------------------------
      Check if a cheap estimate of the reciprocal of the condition
      number is getting too small.  If so, delete
      the prior numeric factorization and recompute it.
      -----------------------------------------------------------*/

    retval = sun_klu_rcond(SYMBOLIC(S), NUMERIC(S), &COMMON(S));
    if (retval == 0) {
      LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
      return(LASTFLAG(S));
    }

    if ( COMMON(S).rcond < uround_twothirds ) {

      /* Condition number may be getting large.
	 Compute more accurate estimate */
      retval = sun_klu_condest((KLU_INDEXTYPE*) SUNSparseMatrix_IndexPointers(A),
                               SUNSparseMatrix_Data(A),
                               SYMBOLIC(S),
                               NUMERIC(S),
                               &COMMON(S));
      if (retval == 0) {
	LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
        return(LASTFLAG(S));
      }

      if ( COMMON(S).condest > (ONE/uround_twothirds) ) {

	/* More accurate estimate also says condition number is
	   large, so recompute the numeric factorization */
	sun_klu_free_numeric(&NUMERIC(S), &COMMON(S));
	NUMERIC(S) = sun_klu_factor((KLU_INDEXTYPE*) SUNSparseMatrix_IndexPointers(A),
                                    (KLU_INDEXTYPE*) SUNSparseMatrix_IndexValues(A),
                                    SUNSparseMatrix_Data(A),
                                    SYMBOLIC(S),
                                    &COMMON(S));
	if (NUMERIC(S) == NULL) {
	  LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
          return(LASTFLAG(S));
	}
      }

    }
  }

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSolve_KLU(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                       N_Vector b, realtype tol)
{
  int flag;
  realtype *xdata;

  /* check for valid inputs */
  if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) )
    return(SUNLS_MEM_NULL);

  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access x data array */
  xdata = N_VGetArrayPointer(x);
  if (xdata == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }

  /* Call KLU to solve the linear system */
  flag = SOLVE(S)(SYMBOLIC(S), NUMERIC(S),
                  SUNSparseMatrix_NP(A), 1, xdata,
                  &COMMON(S));
  if (flag == 0) {
    LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
    return(LASTFLAG(S));
  }

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


sunindextype SUNLinSolLastFlag_KLU(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return(-1);
  return(LASTFLAG(S));
}


int SUNLinSolSpace_KLU(SUNLinearSolver S,
                       long int *lenrwLS,
                       long int *leniwLS)
{
  /* since the klu structures are opaque objects, we
     omit those from these results */
  *leniwLS = 2;
  *lenrwLS = 0;
  return(SUNLS_SUCCESS);
}

int SUNLinSolFree_KLU(SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == NULL) return(SUNLS_SUCCESS);

  /* delete items from the contents structure (if it exists) */
  if (S->content) {
    if (NUMERIC(S))
      sun_klu_free_numeric(&NUMERIC(S), &COMMON(S));
    if (SYMBOLIC(S))
      sun_klu_free_symbolic(&SYMBOLIC(S), &COMMON(S));
    free(S->content);
    S->content = NULL;
  }

  /* delete generic structures */
  if (S->ops) {
    free(S->ops);
    S->ops = NULL;
  }
  free(S); S = NULL;
  return(SUNLS_SUCCESS);
}
