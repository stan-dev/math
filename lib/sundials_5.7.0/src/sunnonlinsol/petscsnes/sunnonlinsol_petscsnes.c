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
 * This is the implementation file for the SUNNonlinearSolver module
 * implementation that interfaces to the PETSc SNES nonlinear solvers.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <petscsnes.h>

#include <nvector/nvector_petsc.h>
#include <sunnonlinsol/sunnonlinsol_petscsnes.h>
#include <sundials/sundials_math.h>

#define SUNNLS_SNES_CONTENT(NLS) ( (SUNNonlinearSolverContent_PetscSNES)(NLS->content) )
#define SUNNLS_SNESOBJ(NLS)      ( SUNNLS_SNES_CONTENT(NLS)->snes )

/* private function which translates the SNESFunction form to the SUNNonlinSolSysFn form */
static PetscErrorCode PetscSysFn(SNES snes, Vec x, Vec f, void *ctx);

/*==============================================================================
  Constructor
  ============================================================================*/

/* create a SUNNonlinearSolver wrapper for the PETSc SNES context */
SUNNonlinearSolver SUNNonlinSol_PetscSNES(N_Vector y, SNES snes)
{
  int ierr;
  SUNNonlinearSolver NLS;
  SUNNonlinearSolverContent_PetscSNES content;

  /* check that the supplied SNES is non-NULL */
  if (snes == NULL || y == NULL) return NULL;

  /* check that the vector is the right type */
  if (N_VGetVectorID(y) != SUNDIALS_NVEC_PETSC) return NULL;

  /*
   * Create an empty nonlinear linear solver object
   */

  NLS = SUNNonlinSolNewEmpty();
  if (NLS == NULL) return NULL;

  /* Attach operations */
  NLS->ops->gettype         = SUNNonlinSolGetType_PetscSNES;
  NLS->ops->initialize      = SUNNonlinSolInitialize_PetscSNES;
  NLS->ops->solve           = SUNNonlinSolSolve_PetscSNES;
  NLS->ops->free            = SUNNonlinSolFree_PetscSNES;
  NLS->ops->setsysfn        = SUNNonlinSolSetSysFn_PetscSNES;
  NLS->ops->getnumiters     = SUNNonlinSolGetNumIters_PetscSNES;
  NLS->ops->getnumconvfails = SUNNonlinSolGetNumConvFails_PetscSNES;

  /*
   * Create content
   */

  content = NULL;
  content = (SUNNonlinearSolverContent_PetscSNES) malloc(sizeof *content);
  if (content == NULL) {
    SUNNonlinSolFree(NLS);
    return NULL;
  }

  /* Initialize all components of content to 0/NULL */
  memset(content, 0, sizeof(struct _SUNNonlinearSolverContent_PetscSNES));

  /* Attach content */
  NLS->content = content;

  /*
   * Fill content
   */

  content->snes = snes;

  /* Create all internal vectors */
  content->y = N_VCloneEmpty(y);
  if (content->y == NULL) {
    SUNNonlinSolFree(NLS);
    return NULL;
  }

  content->f = N_VCloneEmpty(y);
  if (content->f == NULL) {
    SUNNonlinSolFree(NLS);
    return NULL;
  }

  ierr = VecDuplicate(N_VGetVector_Petsc(y), &content->r);
  if (ierr != 0) {
    SUNNonlinSolFree(NLS);
    return NULL;
  }

  /* tell SNES about the sys function */
  ierr = SNESSetFunction(SUNNLS_SNESOBJ(NLS), SUNNLS_SNES_CONTENT(NLS)->r, PetscSysFn, NLS);
  if (ierr != 0) {
    SUNNonlinSolFree(NLS);
    return NULL;
  }

  return(NLS);
}


/*==============================================================================
  GetType, Initialize, Setup, Solve, and Free operations
  ============================================================================*/

/* get the type of SUNNonlinearSolver */
SUNNonlinearSolver_Type SUNNonlinSolGetType_PetscSNES(SUNNonlinearSolver NLS)
{
  return(SUNNONLINEARSOLVER_ROOTFIND);
}

/* performs any initialization needed */
int SUNNonlinSolInitialize_PetscSNES(SUNNonlinearSolver NLS)
{
  PetscErrorCode ptcerr;

  /* set the system function again to ensure it is the wrapper */
  ptcerr = SNESSetFunction(SUNNLS_SNESOBJ(NLS), SUNNLS_SNES_CONTENT(NLS)->r, PetscSysFn, NLS);
  if (ptcerr != 0) {
    SUNNLS_SNES_CONTENT(NLS)->petsc_last_err = ptcerr;
    return SUN_NLS_EXT_FAIL;
  }

  return SUN_NLS_SUCCESS;
}

/*------------------------------------------------------------------------------
  SUNNonlinSolSolve_PetscSNES: Performs the nonlinear solve F(y) = 0 or G(y) = y

  Successful solve return code:
    SUN_NLS_SUCCESS = 0

  Recoverable failure return codes (positive):
    SUN_NLS_CONV_RECVR
    *_RHSFUNC_RECVR (ODEs) or *_RES_RECVR (DAEs)

  Unrecoverable failure return codes (negative):
    *_MEM_NULL
    *_RHSFUNC_FAIL (ODEs) or *_RES_FAIL (DAEs)

  Note return values beginning with * are package specific values returned by
  the Sys function provided to the nonlinear solver.
  ----------------------------------------------------------------------------*/
int SUNNonlinSolSolve_PetscSNES(SUNNonlinearSolver NLS,
                                N_Vector y0, N_Vector y,
                                N_Vector w, realtype tol,
                                booleantype callLSetup, void* mem)
{
  /* local variables */
  PetscErrorCode ierr;
  SNESConvergedReason reason;
  int retval;

  /* check that the inputs are non-null */
  if ( (NLS == NULL) ||
       (y0  == NULL) ||
       (y   == NULL) ||
       (w   == NULL) )
    return SUN_NLS_MEM_NULL;

  /* store a pointer to the integrator memory so it can be
   * accessed in the system function */
  SUNNLS_SNES_CONTENT(NLS)->imem = mem;

  /* reset convergence failure count */
  SUNNLS_SNES_CONTENT(NLS)->nconvfails = 0;

  /* call petsc SNES solve */
  ierr = SNESSolve(SUNNLS_SNESOBJ(NLS), NULL, N_VGetVector_Petsc(y));

  /* check if the call to the system function failed */
  if (SUNNLS_SNES_CONTENT(NLS)->sysfn_last_err != SUN_NLS_SUCCESS)
    return SUNNLS_SNES_CONTENT(NLS)->sysfn_last_err;

  /* check if the SNESSolve had a failure elsewhere */
  if (ierr != 0) {
    SUNNLS_SNES_CONTENT(NLS)->petsc_last_err = ierr;
    return SUN_NLS_EXT_FAIL; /* ierr != 0 is not recoverable with PETSc */
  }

  /*
   * determine if and why the system converged/diverged so we can determine
   * if it is a recoverable failure or success
   */

  ierr = SNESGetConvergedReason(SUNNLS_SNESOBJ(NLS), &reason);
  if (ierr != 0) {
    SUNNLS_SNES_CONTENT(NLS)->petsc_last_err = ierr;
    return SUN_NLS_EXT_FAIL; /* ierr != 0 is not recoverable with PETSc */
  }

  if ( (reason == SNES_CONVERGED_ITERATING) ||
       (reason == SNES_CONVERGED_FNORM_ABS) ||
       (reason == SNES_CONVERGED_FNORM_RELATIVE) ||
       (reason == SNES_CONVERGED_SNORM_RELATIVE) ) {
    /* success */
    retval = SUN_NLS_SUCCESS;
  } else {
    /* recoverable failure */
    retval = SUN_NLS_CONV_RECVR;
    /* update convergence failure count */
    SUNNLS_SNES_CONTENT(NLS)->nconvfails++;
  }

  return retval;
}

/* free the SUNNonlinearSolver */
int SUNNonlinSolFree_PetscSNES(SUNNonlinearSolver NLS)
{
  /* return if NLS is already free */
  if (NLS == NULL)
    return SUN_NLS_SUCCESS;

  /* free items from contents, then the generic structure */
  if (NLS->content) {
    if (SUNNLS_SNES_CONTENT(NLS)->r) VecDestroy(&SUNNLS_SNES_CONTENT(NLS)->r);
    if (SUNNLS_SNES_CONTENT(NLS)->y) N_VDestroy_Petsc(SUNNLS_SNES_CONTENT(NLS)->y);
    if (SUNNLS_SNES_CONTENT(NLS)->f) N_VDestroy_Petsc(SUNNLS_SNES_CONTENT(NLS)->f);
    free(NLS->content);
    NLS->content = NULL;
  }

  /* free the ops structure */
  if (NLS->ops) {
    free(NLS->ops);
    NLS->ops = NULL;
  }

  /* free the nonlinear solver */
  free(NLS);

  return SUN_NLS_SUCCESS;
}


/*==============================================================================
  Set functions
  ============================================================================*/

/* set the system residual function */
int SUNNonlinSolSetSysFn_PetscSNES(SUNNonlinearSolver NLS, SUNNonlinSolSysFn SysFn)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  /* check that the nonlinear system function is non-null */
  if (SysFn == NULL)
    return(SUN_NLS_ILL_INPUT);

  SUNNLS_SNES_CONTENT(NLS)->Sys = SysFn;
  return SUN_NLS_SUCCESS;
}

/*==============================================================================
  Get functions
  ============================================================================*/


/* get the PETSc SNES context underneath the SUNNonlinearSolver object */
int SUNNonlinSolGetSNES_PetscSNES(SUNNonlinearSolver NLS, SNES *snes)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  /* return the SNES context */
  *snes = SUNNLS_SNESOBJ(NLS);
  return SUN_NLS_SUCCESS;
}

/* get the last error return by SNES */
int SUNNonlinSolGetPetscError_PetscSNES(SUNNonlinearSolver NLS, PetscErrorCode* err)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  /* return the last PETSc error code returned by SNES */
  *err = SUNNLS_SNES_CONTENT(NLS)->petsc_last_err;
  return SUN_NLS_SUCCESS;
}

/* get a pointer to the SUNDIALS integrator-provided system function F(y) */
int SUNNonlinSolGetSysFn_PetscSNES(SUNNonlinearSolver NLS, SUNNonlinSolSysFn* SysFn)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  /* return the nonlinear system defining function */
  *SysFn = SUNNLS_SNES_CONTENT(NLS)->Sys;
  return SUN_NLS_SUCCESS;
}

/* get the number of iterations performed in the last solve */
int SUNNonlinSolGetNumIters_PetscSNES(SUNNonlinearSolver NLS, long int* nni)
{
  int ierr;
  sunindextype niters;

  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  /* get iteration count */
  ierr = SNESGetIterationNumber(SUNNLS_SNESOBJ(NLS), &niters);
  if (ierr != 0) {
    SUNNLS_SNES_CONTENT(NLS)->petsc_last_err = ierr;
    return SUN_NLS_EXT_FAIL; /* ierr != 0 is not recoverable with PETSc */
  }

  *nni = (long int) niters;

  return SUN_NLS_SUCCESS;
}

/* get the total number of nonlinear convergence failures in the
   lifetime of this SUNNonlinearSolver object */
int SUNNonlinSolGetNumConvFails_PetscSNES(SUNNonlinearSolver NLS,
                                          long int* nconvfails)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  *nconvfails = SUNNLS_SNES_CONTENT(NLS)->nconvfails;

  return SUN_NLS_SUCCESS;
}

/*==============================================================================
  Private functions
  ============================================================================*/

static PetscErrorCode PetscSysFn(SNES snes, Vec x, Vec f, void *ctx)
{
  int retval;
  SUNNonlinearSolver NLS = (SUNNonlinearSolver) ctx;

  /* wrap up the petsc vectors in nvectors */
  N_VSetVector_Petsc(SUNNLS_SNES_CONTENT(NLS)->y, x);
  N_VSetVector_Petsc(SUNNLS_SNES_CONTENT(NLS)->f, f);

  /* now call the SUNDIALS format system function to evaluate F(y) */
  retval = SUNNLS_SNES_CONTENT(NLS)->Sys(SUNNLS_SNES_CONTENT(NLS)->y,
                                         SUNNLS_SNES_CONTENT(NLS)->f,
                                         SUNNLS_SNES_CONTENT(NLS)->imem);
  /* store the return value */
  SUNNLS_SNES_CONTENT(NLS)->sysfn_last_err = retval;

  /* if there was some sort of failure in the system function, return -1 to
   * indicate an error instead of retval so that we don't overlap with one of
   * the standard PETSc error codes */
  return (retval != 0) ? -1 : 0;
}
