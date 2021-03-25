/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for the SuperLU-DIST implementation of the
 * SUNLINSOL package.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunmatrix/sunmatrix_slunrloc.h>
#include <sunlinsol/sunlinsol_superludist.h>
#include <sundials/sundials_math.h>
#include <superlu_ddefs.h>

#define ZERO      RCONST(0.0)
#define ONE       RCONST(1.0)
#define TWO       RCONST(2.0)

/*
 * -----------------------------------------------------------------
 * SuperLUDIST solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SLUDIST_CONTENT(S)    ( (SUNLinearSolverContent_SuperLUDIST)(S->content) )
#define SLU_FIRSTFACTORIZE(S) ( SLUDIST_CONTENT(S)->first_factorize )
#define SLU_BERR(S)           ( SLUDIST_CONTENT(S)->berr )
#define SLU_GRID(S)           ( SLUDIST_CONTENT(S)->grid )
#define SLU_LASTFLAG(S)       ( SLUDIST_CONTENT(S)->last_flag )
#define SLU_LU(S)             ( SLUDIST_CONTENT(S)->lu )
#define SLU_OPTIONS(S)        ( SLUDIST_CONTENT(S)->options )
#define SLU_SCALEPERM(S)      ( SLUDIST_CONTENT(S)->scaleperm )
#define SLU_SOLVESTRUCT(S)    ( SLUDIST_CONTENT(S)->solve )
#define SLU_STAT(S)           ( SLUDIST_CONTENT(S)->stat )
#define SLU_SIZE(S)           ( SLUDIST_CONTENT(S)->N )

/*
 * ----------------------------------------------------------------------------
 *  Implementations of exported functions.
 * ----------------------------------------------------------------------------
 */

SUNLinearSolver SUNLinSol_SuperLUDIST(N_Vector y, SUNMatrix A, gridinfo_t *grid,
                                      xLUstruct_t *lu, xScalePermstruct_t *scaleperm,
                                      xSOLVEstruct_t *solve, SuperLUStat_t *stat,
                                      superlu_dist_options_t *options)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_SuperLUDIST content;
  SuperMatrix *Asuper;
  NRformat_loc *Astore;

  /* Check that required arguments are not NULL */
  if (y == NULL || A == NULL || grid == NULL || lu == NULL ||
      scaleperm == NULL || solve == NULL || stat == NULL || options == NULL)
    return(NULL);

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SLUNRLOC) return(NULL);

  Asuper = SUNMatrix_SLUNRloc_SuperMatrix(A);
  Astore = (NRformat_loc *) Asuper->Store;

  /* Check that the matrix is square */
  if (Asuper->nrow != Asuper->ncol) return(NULL);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty();
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_SuperLUDIST;
  S->ops->getid      = SUNLinSolGetID_SuperLUDIST;
  S->ops->initialize = SUNLinSolInitialize_SuperLUDIST;
  S->ops->setup      = SUNLinSolSetup_SuperLUDIST;
  S->ops->solve      = SUNLinSolSolve_SuperLUDIST;
  S->ops->lastflag   = SUNLinSolLastFlag_SuperLUDIST;
  S->ops->space      = SUNLinSolSpace_SuperLUDIST;
  S->ops->free       = SUNLinSolFree_SuperLUDIST;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SuperLUDIST) malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->grid      = grid;
  content->lu        = lu;
  content->N         = Astore->m_loc;
  content->options   = options;
  content->scaleperm = scaleperm;
  content->solve     = solve;
  content->stat      = stat;
  content->berr      = ZERO;
  content->last_flag = 0;

  return(S);
}

/*
 * -----------------------------------------------------------------
 * Implementation of accessor functions.
 * -----------------------------------------------------------------
 */

realtype SUNLinSol_SuperLUDIST_GetBerr(SUNLinearSolver LS)
{
  return(SLU_BERR(LS));
}

gridinfo_t* SUNLinSol_SuperLUDIST_GetGridinfo(SUNLinearSolver LS)
{
  return(SLU_GRID(LS));
}

xLUstruct_t* SUNLinSol_SuperLUDIST_GetLUstruct(SUNLinearSolver LS)
{
  return(SLU_LU(LS));
}

superlu_dist_options_t* SUNLinSol_SuperLUDIST_GetSuperLUOptions(SUNLinearSolver LS)
{
  return(SLU_OPTIONS(LS));
}

xScalePermstruct_t* SUNLinSol_SuperLUDIST_GetScalePermstruct(SUNLinearSolver LS)
{
  return(SLU_SCALEPERM(LS));
}

xSOLVEstruct_t* SUNLinSol_SuperLUDIST_GetSOLVEstruct(SUNLinearSolver LS)
{
  return(SLU_SOLVESTRUCT(LS));
}

SuperLUStat_t* SUNLinSol_SuperLUDIST_GetSuperLUStat(SUNLinearSolver LS)
{
  return(SLU_STAT(LS));
}

/*
 * -----------------------------------------------------------------
 * Implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SuperLUDIST(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_SuperLUDIST(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_SUPERLUDIST);
}

int SUNLinSolInitialize_SuperLUDIST(SUNLinearSolver S)
{
  SLU_FIRSTFACTORIZE(S) = SUNTRUE;

  SLU_LASTFLAG(S) = SUNLS_SUCCESS;
  return(SLU_LASTFLAG(S));
}

int SUNLinSolSetup_SuperLUDIST(SUNLinearSolver S, SUNMatrix A)
{
  if (SLU_FIRSTFACTORIZE(S)) {
    /* Force factorization on next solve. */
    SLU_OPTIONS(S)->Fact = DOFACT;

    /* if the solve struct was already initialized, we need to
       finalize the last solve to avoid leaking memory */
    if (SLU_OPTIONS(S)->SolveInitialized == YES) {
      xDestroy_LU(SLU_SIZE(S), SLU_GRID(S), SLU_LU(S));
      dSolveFinalize(SLU_OPTIONS(S), SLU_SOLVESTRUCT(S));
    }
  } else {
    /* Force factorization on next solve, but reuse column permutation */
    SLU_OPTIONS(S)->Fact = SamePattern;
  }

  SLU_LASTFLAG(S) = SUNLS_SUCCESS;
  return(SLU_LASTFLAG(S));
}

int SUNLinSolSolve_SuperLUDIST(SUNLinearSolver S, SUNMatrix A,
                               N_Vector x, N_Vector b, realtype tol)
{
  int retval;
  realtype *xdata;
  SuperMatrix *Asuper;
  NRformat_loc *Astore;

  if ((S == NULL) || (A == NULL) || (x == NULL) || (b == NULL))
    return(SUNLS_MEM_NULL);

  Asuper = SUNMatrix_SLUNRloc_SuperMatrix(A);
  Astore = (NRformat_loc *) Asuper->Store;

  /* Copy b into x, and use xdata with SuperLU-DIST because
     SuperLU reuses b to store the solution vector. */
  N_VScale(ONE, b, x);
  xdata = N_VGetArrayPointer(x);
  if (xdata == NULL) {
    SLU_LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(SLU_LASTFLAG(S));
  }

  /* Call SuperLU-DIST to solve the linear system */
  pdgssvx(SLU_OPTIONS(S), Asuper, SLU_SCALEPERM(S), xdata, Astore->m_loc, 1,
          SLU_GRID(S), SLU_LU(S), SLU_SOLVESTRUCT(S), &SLU_BERR(S), SLU_STAT(S),
          &retval);

  if (retval != 0) {
    /* retval should never be < 0, and if it is > than ncol it means a memory
       allocation error occured. If 0 < retval <= ncol, then U is singular and
       the solve couldn't be completed. */
    if (retval < 0 || retval > Asuper->ncol)
      SLU_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
    else
      SLU_LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
    return(SLU_LASTFLAG(S));
  }

  /* Tell SuperLU-DIST that A is already factored so do not factorize */
  SLU_OPTIONS(S)->Fact = FACTORED;

  SLU_LASTFLAG(S) = SUNLS_SUCCESS;
  return(SLU_LASTFLAG(S));
}

sunindextype SUNLinSolLastFlag_SuperLUDIST(SUNLinearSolver S)
{
  if (S == NULL) return(-1);
  return(SLU_LASTFLAG(S));
}

int SUNLinSolFree_SuperLUDIST(SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == NULL) return(SUNLS_SUCCESS);

  /* Call SuperLU DIST destroy/finalize routines,
     but don't free the sturctures themselves - that is the user's job */
  xDestroy_LU(SLU_SIZE(S), SLU_GRID(S), SLU_LU(S));
  dSolveFinalize(SLU_OPTIONS(S), SLU_SOLVESTRUCT(S));

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

int SUNLinSolSpace_SuperLUDIST(SUNLinearSolver S, long int *leniwLS,
                               long int *lenrwLS)
{
  /* since the SuperLU_DIST structures are opaque objects, we
     omit those from these results */

  *leniwLS = 2; /* last_flag, N */
  *lenrwLS = 1; /* berr */

  return(SUNLS_SUCCESS);
}
