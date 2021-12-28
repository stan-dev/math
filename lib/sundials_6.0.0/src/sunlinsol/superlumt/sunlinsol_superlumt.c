/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * -----------------------------------------------------------------
 * Based on codes <solver>_superlumt.c, written by
 * Carol S. Woodward @ LLNL
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
 * This is the implementation file for the SuperLUMT implementation
 * of the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_superlumt.h>
#include <sundials/sundials_math.h>

#define ZERO      RCONST(0.0)
#define ONE       RCONST(1.0)
#define TWO       RCONST(2.0)

/*
 * -----------------------------------------------------------------
 * SuperLUMT solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SLUMT_CONTENT(S)    ( (SUNLinearSolverContent_SuperLUMT)(S->content) )
#define LASTFLAG(S)         ( SLUMT_CONTENT(S)->last_flag )
#define FIRSTFACTORIZE(S)   ( SLUMT_CONTENT(S)->first_factorize )
#define SM_A(S)             ( SLUMT_CONTENT(S)->A )
#define SM_AC(S)            ( SLUMT_CONTENT(S)->AC )
#define SM_L(S)             ( SLUMT_CONTENT(S)->L )
#define SM_U(S)             ( SLUMT_CONTENT(S)->U )
#define SM_B(S)             ( SLUMT_CONTENT(S)->B )
#define GSTAT(S)            ( SLUMT_CONTENT(S)->Gstat )
#define PERMR(S)            ( SLUMT_CONTENT(S)->perm_r )
#define PERMC(S)            ( SLUMT_CONTENT(S)->perm_c )
#define SIZE(S)             ( SLUMT_CONTENT(S)->N )
#define NUMTHREADS(S)       ( SLUMT_CONTENT(S)->num_threads )
#define DIAGPIVOTTHRESH(S)  ( SLUMT_CONTENT(S)->diag_pivot_thresh )
#define ORDERING(S)         ( SLUMT_CONTENT(S)->ordering )
#define OPTIONS(S)          ( SLUMT_CONTENT(S)->options )

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SuperLUMT linear solver
 */

SUNLinearSolver SUNLinSol_SuperLUMT(N_Vector y, SUNMatrix A, int num_threads, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_SuperLUMT content;
  sunindextype MatrixRows;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) return(NULL);

  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A)) return(NULL);

  if ( (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS) )
    return(NULL);

  MatrixRows = SUNSparseMatrix_Rows(A);
  if (MatrixRows != N_VGetLength(y)) return(NULL);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_SuperLUMT;
  S->ops->getid      = SUNLinSolGetID_SuperLUMT;
  S->ops->initialize = SUNLinSolInitialize_SuperLUMT;
  S->ops->setup      = SUNLinSolSetup_SuperLUMT;
  S->ops->solve      = SUNLinSolSolve_SuperLUMT;
  S->ops->lastflag   = SUNLinSolLastFlag_SuperLUMT;
  S->ops->space      = SUNLinSolSpace_SuperLUMT;
  S->ops->free       = SUNLinSolFree_SuperLUMT;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SuperLUMT) malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->N                 = MatrixRows;
  content->last_flag         = 0;
  content->num_threads       = num_threads;
  content->diag_pivot_thresh = ONE;
  content->ordering          = SUNSLUMT_ORDERING_DEFAULT;
  content->perm_r            = NULL;
  content->perm_c            = NULL;
  content->Gstat             = NULL;
  content->A                 = NULL;
  content->AC                = NULL;
  content->L                 = NULL;
  content->U                 = NULL;
  content->B                 = NULL;
  content->options           = NULL;

  /* Allocate content */
  content->perm_r = (sunindextype *) malloc(MatrixRows*sizeof(sunindextype));
  if (content->perm_r == NULL) { SUNLinSolFree(S); return(NULL); }

  content->perm_c = (sunindextype *) malloc(MatrixRows*sizeof(sunindextype));
  if (content->perm_c == NULL) { SUNLinSolFree(S); return(NULL); }

  content->Gstat = (Gstat_t *) malloc(sizeof(Gstat_t));
  if (content->Gstat == NULL) { SUNLinSolFree(S); return(NULL); }

  content->A = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->A == NULL) { SUNLinSolFree(S); return(NULL); }
  content->A->Store = NULL;

  content->AC = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->AC == NULL) { SUNLinSolFree(S); return(NULL); }
  content->AC->Store = NULL;

  content->L = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->L == NULL) { SUNLinSolFree(S); return(NULL); }
  content->L->Store = NULL;

  content->U = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->U == NULL) { SUNLinSolFree(S); return(NULL); }
  content->U->Store = NULL;

  content->B = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->B == NULL) { SUNLinSolFree(S); return(NULL); }
  content->B->Store = NULL;
  xCreate_Dense_Matrix(content->B, MatrixRows, 1, NULL, MatrixRows, SLU_DN,
                       SLU_D, SLU_GE);

  content->options = (superlumt_options_t *) malloc(sizeof(superlumt_options_t));
  if (content->options == NULL) { SUNLinSolFree(S); return(NULL); }
  StatAlloc(MatrixRows, num_threads, sp_ienv(1), sp_ienv(2), content->Gstat);

  return(S);
}


/* ----------------------------------------------------------------------------
 * Function to set the ordering type for a SuperLUMT linear solver
 */

int SUNLinSol_SuperLUMTSetOrdering(SUNLinearSolver S, int ordering_choice)
{
  /* Check for legal ordering_choice */
  if ((ordering_choice < 0) || (ordering_choice > 3))
    return(SUNLS_ILL_INPUT);

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SUNLS_MEM_NULL);

  /* Set ordering_choice */
  ORDERING(S) = ordering_choice;

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SuperLUMT(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}


SUNLinearSolver_ID SUNLinSolGetID_SuperLUMT(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_SUPERLUMT);
}


int SUNLinSolInitialize_SuperLUMT(SUNLinearSolver S)
{
  /* force a first factorization */
  FIRSTFACTORIZE(S) = 1;

  /* Initialize statistics variables */
  StatInit(SIZE(S), NUMTHREADS(S), GSTAT(S));

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSetup_SuperLUMT(SUNLinearSolver S, SUNMatrix A)
{
  int_t retval;
  int panel_size, relax, lwork;
  double drop_tol;
  fact_t fact;
  trans_t trans;
  yes_no_t refact, usepr;
  void *work;

  /* Set option values for SuperLU_MT */
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  fact = EQUILIBRATE;
  trans = (SUNSparseMatrix_SparseType(A) == CSC_MAT) ? NOTRANS : TRANS;
  usepr = NO;
  drop_tol = ZERO;
  lwork = 0;
  work = NULL;

  /* free and reallocate sparse matrix */
  if (SM_A(S)->Store)
    SUPERLU_FREE(SM_A(S)->Store);
  xCreate_CompCol_Matrix(SM_A(S), SUNSparseMatrix_Rows(A),
			 SUNSparseMatrix_Columns(A),
                         SUNSparseMatrix_NNZ(A),
                         SUNSparseMatrix_Data(A),
                         (int_t*) SUNSparseMatrix_IndexValues(A),
                         (int_t*) SUNSparseMatrix_IndexPointers(A),
			 SLU_NC, SLU_D, SLU_GE);

  /* On first decomposition, set up reusable pieces */
  if (FIRSTFACTORIZE(S)) {

    /* Get column permutation vector perm_c[], according to ordering */
    get_perm_c(ORDERING(S), SM_A(S), (int_t *) PERMC(S));
    refact = NO;
    FIRSTFACTORIZE(S) = 0;

  } else {

    /* Re-initialize statistics variables */
    StatInit(SIZE(S), NUMTHREADS(S), GSTAT(S));
    Destroy_CompCol_Permuted(SM_AC(S));
    refact = YES;

  }

  /* Initialize the option structure using the user-input parameters.
     Subsequent calls will re-initialize options.  Apply perm_c to
     columns of original A to form AC */
  pxgstrf_init(NUMTHREADS(S), fact, trans, refact, panel_size, relax,
	       DIAGPIVOTTHRESH(S), usepr, drop_tol, (int_t *) PERMC(S), (int_t *) PERMR(S),
               work, lwork, SM_A(S), SM_AC(S), OPTIONS(S), GSTAT(S));

  /* Compute the LU factorization of A.
     The following routine will create num_threads threads. */
  pxgstrf(OPTIONS(S), SM_AC(S), (int_t *) PERMR(S), SM_L(S), SM_U(S),
          GSTAT(S), &retval);
  if (retval != 0) {
    LASTFLAG(S) = (retval < 0) ?
      SUNLS_PACKAGE_FAIL_UNREC : SUNLS_PACKAGE_FAIL_REC;
    return(LASTFLAG(S));
  }

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSolve_SuperLUMT(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                             N_Vector b, realtype tol)
{
  int_t retval;
  realtype *xdata;
  DNformat *Bstore;
  trans_t trans;

  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access x data array */
  xdata = N_VGetArrayPointer(x);
  if (xdata == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }

  Bstore = (DNformat *) (SM_B(S)->Store);
  Bstore->nzval = xdata;

  /* Call SuperLUMT to solve the linear system using L and U */
  trans = (SUNSparseMatrix_SparseType(A) == CSC_MAT) ? NOTRANS : TRANS;
  xgstrs(trans, SM_L(S), SM_U(S), (int_t *) PERMR(S), (int_t *) PERMC(S), SM_B(S), GSTAT(S), &retval);
  if (retval != 0) {
    LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
    return(LASTFLAG(S));
  }

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


sunindextype SUNLinSolLastFlag_SuperLUMT(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return(-1);
  return(LASTFLAG(S));
}


int SUNLinSolSpace_SuperLUMT(SUNLinearSolver S,
                             long int *lenrwLS,
                             long int *leniwLS)
{
  /* since the SuperLU_MT structures are opaque objects, we
     omit those from these results */
  *leniwLS = 5 + 2*SIZE(S);
  *lenrwLS = 1;
  return(SUNLS_SUCCESS);
}

int SUNLinSolFree_SuperLUMT(SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == NULL) return(SUNLS_SUCCESS);

  /* delete items from the contents structure (if it exists) */
  if (S->content) {
    if (OPTIONS(S) && SM_AC(S)) {
      pxgstrf_finalize(OPTIONS(S), SM_AC(S));
    }
    if (PERMR(S)) {
      free(PERMR(S));
      PERMR(S) = NULL;
    }
    if (PERMC(S)) {
      free(PERMC(S));
      PERMC(S) = NULL;
    }
    if (OPTIONS(S)) {
      free(OPTIONS(S));
      OPTIONS(S) = NULL;
    }
    if (SM_L(S)) {
      Destroy_SuperNode_SCP(SM_L(S));
      SM_L(S) = NULL;
    }
    if (SM_U(S)) {
      Destroy_CompCol_NCP(SM_U(S));
      SM_U(S) = NULL;
    }
    if (GSTAT(S)) {
      StatFree(GSTAT(S));
      free(GSTAT(S));
      GSTAT(S) = NULL;
    }
    if (SM_B(S)) {
      Destroy_SuperMatrix_Store(SM_B(S));
      SM_B(S) = NULL;
    }
    if (SM_A(S)) {
      if (SM_A(S)->Store) {
        SUPERLU_FREE(SM_A(S)->Store);
        SM_A(S)->Store = NULL;
      }
      free(SM_A(S));
      SM_A(S) = NULL;
    }
    if (SM_B(S)) {
      free(SM_B(S));
      SM_B(S) = NULL;
    }
    if (SM_AC(S)) {
      free(SM_AC(S));
      SM_AC(S) = NULL;
    }
    if (SM_L(S)) {
      free(SM_L(S));
      SM_L(S) = NULL;
    }
    if (SM_U(S)) {
      free(SM_U(S));
      SM_U(S) = NULL;
    }
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
