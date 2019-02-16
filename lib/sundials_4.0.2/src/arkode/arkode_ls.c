/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Implementation file for ARKode's linear solver interface.
 *---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_ls_impl.h"
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

/* constants */
#define MIN_INC_MULT RCONST(1000.0)
#define MAX_DQITERS  3  /* max. # of attempts to recover in DQ J*v */
#define ZERO         RCONST(0.0)
#define PT25         RCONST(0.25)
#define ONE          RCONST(1.0)


/*===============================================================
  ARKLS utility routines (called by time-stepper modules)
  ===============================================================*/

/*---------------------------------------------------------------
  arkLSSetLinearSolver specifies the linear solver.
  ---------------------------------------------------------------*/
int arkLSSetLinearSolver(void *arkode_mem, SUNLinearSolver LS,
                         SUNMatrix A)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval, LSType;

  /* Return immediately if either arkode_mem or LS inputs are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKLS_MEM_NULL, "ARKLS",
                    "arkLSSetLinearSolver", MSG_LS_ARKMEM_NULL);
    return(ARKLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (LS == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetLinearSolver",
                    "LS must be non-NULL");
    return(ARKLS_ILL_INPUT);
  }

  /* Test if solver is compatible with LS interface */
  if ( (LS->ops->gettype == NULL) ||
       (LS->ops->initialize == NULL) ||
       (LS->ops->setup == NULL) ||
       (LS->ops->solve == NULL) ) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                   "arkLSSetLinearSolver",
                   "LS object is missing a required operation");
    return(ARKLS_ILL_INPUT);
  }

  /* Test if vector is compatible with LS interface */
  if ( (ark_mem->tempv1->ops->nvconst == NULL) ||
       (ark_mem->tempv1->ops->nvdotprod == NULL) ) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetLinearSolver", MSG_LS_BAD_NVECTOR);
    return(ARKLS_ILL_INPUT);
  }

  /* Retrieve the LS type */
  LSType = SUNLinSolGetType(LS);

  /* Check for compatible LS type, matrix and "atimes" support */
  if ((LSType == SUNLINEARSOLVER_ITERATIVE) && (LS->ops->setatimes == NULL)) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLSSetLinearSolver",
                    "Incompatible inputs: iterative LS must support ATimes routine");
    return(ARKLS_ILL_INPUT);
  }
  if ((LSType == SUNLINEARSOLVER_DIRECT) && (A == NULL)) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLSSetLinearSolver",
                    "Incompatible inputs: direct LS requires non-NULL matrix");
    return(ARKLS_ILL_INPUT);
  }
  if ((LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) && (A == NULL)) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLSSetLinearSolver",
                    "Incompatible inputs: matrix-iterative LS requires non-NULL matrix");
    return(ARKLS_ILL_INPUT);
  }


  /* Test whether time stepper module is supplied, with required routines */
  if ( (ark_mem->step_attachlinsol == NULL) ||
       (ark_mem->step_getlinmem == NULL) ||
       (ark_mem->step_getimplicitrhs == NULL) ||
       (ark_mem->step_getgammas == NULL) ) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetLinearSolver",
                    "Missing time step module or associated routines");
    return(ARKLS_ILL_INPUT);
  }

  /* Allocate memory for ARKLsMemRec */
  arkls_mem = NULL;
  arkls_mem = (ARKLsMem) malloc(sizeof(struct ARKLsMemRec));
  if (arkls_mem == NULL) {
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKLS",
                    "arkLSSetLinearSolver", MSG_LS_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }
  memset(arkls_mem, 0, sizeof(struct ARKLsMemRec));

  /* set SUNLinearSolver pointer */
  arkls_mem->LS = LS;

  /* Set defaults for Jacobian-related fields */
  if (A != NULL) {
    arkls_mem->jacDQ  = SUNTRUE;
    arkls_mem->jac    = arkLsDQJac;
    arkls_mem->J_data = ark_mem;
  } else {
    arkls_mem->jacDQ  = SUNFALSE;
    arkls_mem->jac    = NULL;
    arkls_mem->J_data = NULL;
  }
  arkls_mem->jtimesDQ = SUNTRUE;
  arkls_mem->jtsetup  = NULL;
  arkls_mem->jtimes   = arkLsDQJtimes;
  arkls_mem->Jt_data  = ark_mem;

  /* Set defaults for preconditioner-related fields */
  arkls_mem->pset   = NULL;
  arkls_mem->psolve = NULL;
  arkls_mem->pfree  = NULL;
  arkls_mem->P_data = ark_mem->user_data;

  /* Initialize counters */
  arkLsInitializeCounters(arkls_mem);

  /* Set default values for the rest of the LS parameters */
  arkls_mem->msbj      = ARKLS_MSBJ;
  arkls_mem->jbad      = SUNTRUE;
  arkls_mem->eplifac   = ARKLS_EPLIN;
  arkls_mem->last_flag = ARKLS_SUCCESS;

  /* If LS supports ATimes, attach ARKLs routine */
  if (LS->ops->setatimes) {
    retval = SUNLinSolSetATimes(LS, ark_mem, arkLsATimes);
    if (retval != SUNLS_SUCCESS) {
      arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKLS",
                      "arkLSSetLinearSolver",
                      "Error in calling SUNLinSolSetATimes");
      free(arkls_mem); arkls_mem = NULL;
      return(ARKLS_SUNLS_FAIL);
    }
  }

  /* If LS supports preconditioning, initialize pset/psol to NULL */
  if (LS->ops->setpreconditioner) {
    retval = SUNLinSolSetPreconditioner(LS, ark_mem, NULL, NULL);
    if (retval != SUNLS_SUCCESS) {
      arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKLS",
                      "arkLSSetLinearSolver",
                      "Error in calling SUNLinSolSetPreconditioner");
      free(arkls_mem); arkls_mem = NULL;
      return(ARKLS_SUNLS_FAIL);
    }
  }

  /* When using a non-NULL SUNMatrix object, store pointer to A and create saved_J */
  if (A != NULL) {
    arkls_mem->A = A;
    arkls_mem->savedJ = SUNMatClone(A);
    if (arkls_mem->savedJ == NULL) {
      arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKLS",
                      "arkLSSetLinearSolver", MSG_LS_MEM_FAIL);
      free(arkls_mem); arkls_mem = NULL;
      return(ARKLS_MEM_FAIL);
    }
  }

  /* Allocate memory for ytemp and x */
  arkls_mem->ytemp = N_VClone(ark_mem->tempv1);
  if (arkls_mem->ytemp == NULL) {
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKLS",
                    "arkLSSetLinearSolver", MSG_LS_MEM_FAIL);
    SUNMatDestroy(arkls_mem->savedJ);
    free(arkls_mem); arkls_mem = NULL;
    return(ARKLS_MEM_FAIL);
  }

  arkls_mem->x = N_VClone(ark_mem->tempv1);
  if (arkls_mem->x == NULL) {
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKLS",
                    "arkLSSetLinearSolver", MSG_LS_MEM_FAIL);
    N_VDestroy(arkls_mem->ytemp);
    SUNMatDestroy(arkls_mem->savedJ);
    free(arkls_mem); arkls_mem = NULL;
    return(ARKLS_MEM_FAIL);
  }

  /* For iterative LS, compute sqrtN from a dot product */
  if ( (LSType == SUNLINEARSOLVER_ITERATIVE) ||
       (LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) ) {
    N_VConst(ONE, arkls_mem->ytemp);
    arkls_mem->sqrtN = SUNRsqrt( N_VDotProd(arkls_mem->ytemp,
                                            arkls_mem->ytemp) );
  }

  /* Attach ARKLs interface to time stepper module */
  retval = ark_mem->step_attachlinsol(arkode_mem, arkLsInitialize,
                                      arkLsSetup, arkLsSolve,
                                      arkLsFree, 2, arkls_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKLS", "arkLSSetLinearSolver",
                    "Failed to attach to time stepper module");
    N_VDestroy(arkls_mem->x);
    N_VDestroy(arkls_mem->ytemp);
    SUNMatDestroy(arkls_mem->savedJ);
    free(arkls_mem); arkls_mem = NULL;
    return(retval);
  }

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSSetMassLinearSolver specifies the iterative mass-matrix
  linear solver and user-supplied routine to perform the
  mass-matrix-vector product.
  ---------------------------------------------------------------*/
int arkLSSetMassLinearSolver(void *arkode_mem, SUNLinearSolver LS,
                             SUNMatrix M, booleantype time_dep)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval, LSType;

  /* Return immediately if either arkode_mem or LS inputs are NULL */
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARKLS_MEM_NULL, "ARKLS",
                    "arkLSSetMassLinearSolver",
                    MSG_LS_ARKMEM_NULL);
    return(ARKLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (LS == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetMassLinearSolver",
                    "LS must be non-NULL");
    return(ARKLS_ILL_INPUT);
  }

  /* Test if solver is compatible with LS interface */
  if ( (LS->ops->gettype == NULL) ||
       (LS->ops->initialize == NULL) ||
       (LS->ops->setup == NULL) ||
       (LS->ops->solve == NULL) ) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                   "arkLSSetMassLinearSolver",
                   "LS object is missing a required operation");
    return(ARKLS_ILL_INPUT);
  }

  /* Test if vector is compatible with LS interface */
  if ( (ark_mem->tempv1->ops->nvconst == NULL) ||
       (ark_mem->tempv1->ops->nvdotprod == NULL) ){
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetMassLinearSolver", MSG_LS_BAD_NVECTOR);
    return(ARKLS_ILL_INPUT);
  }

  /* Retrieve the LS type */
  LSType = SUNLinSolGetType(LS);

  /* Check for compatible LS type, matrix and "atimes" support */
  if ((LSType == SUNLINEARSOLVER_ITERATIVE) && (LS->ops->setatimes == NULL)) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLSSetMassLinearSolver",
                    "Incompatible inputs: iterative LS must support ATimes routine");
    return(ARKLS_ILL_INPUT);
  }
  if ((LSType == SUNLINEARSOLVER_DIRECT) && (M == NULL)) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLSSetMassLinearSolver",
                    "Incompatible inputs: direct LS requires non-NULL matrix");
    return(ARKLS_ILL_INPUT);
  }
  if ((LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) && (M == NULL)) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLSSetMassLinearSolver",
                    "Incompatible inputs: matrix-iterative LS requires non-NULL matrix");
    return(ARKLS_ILL_INPUT);
  }

  /* Test whether time stepper module is supplied, with required routines */
  if ( (ark_mem->step_attachmasssol == NULL) ||
       (ark_mem->step_getmassmem == NULL) ) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetMassLinearSolver",
                    "Missing time step module or associated routines");
    return(ARKLS_ILL_INPUT);
  }

  /* Allocate memory for ARKLsMemRec */
  arkls_mem = NULL;
  arkls_mem = (ARKLsMassMem) malloc(sizeof(struct ARKLsMassMemRec));
  if (arkls_mem == NULL) {
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKLS",
                    "arkLSSetMassLinearSolver", MSG_LS_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }
  memset(arkls_mem, 0, sizeof(struct ARKLsMassMemRec));

  /* set SUNLinearSolver pointer; flag indicating time-dependence */
  arkls_mem->LS = LS;
  arkls_mem->time_dependent = time_dep;

  /* Set mass-matrix routines to NULL */
  arkls_mem->mass    = NULL;
  arkls_mem->mtsetup = NULL;
  arkls_mem->mtimes  = NULL;
  arkls_mem->mt_data = NULL;

  /* Set defaults for preconditioner-related fields */
  arkls_mem->pset   = NULL;
  arkls_mem->psolve = NULL;
  arkls_mem->pfree  = NULL;
  arkls_mem->P_data = ark_mem->user_data;

  /* Initialize counters */
  arkLsInitializeMassCounters(arkls_mem);

  /* Set default values for the rest of the LS parameters */
  arkls_mem->eplifac   = ARKLS_EPLIN;
  arkls_mem->last_flag = ARKLS_SUCCESS;

  /* If LS supports ATimes, attach ARKLs routine */
  if (LS->ops->setatimes) {
    retval = SUNLinSolSetATimes(LS, ark_mem, NULL);
    if (retval != SUNLS_SUCCESS) {
      arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKLS",
                      "arkLSSetMassLinearSolver",
                      "Error in calling SUNLinSolSetATimes");
      free(arkls_mem); arkls_mem = NULL;
      return(ARKLS_SUNLS_FAIL);
    }
  }

  /* If LS supports preconditioning, initialize pset/psol to NULL */
  if (LS->ops->setpreconditioner) {
    retval = SUNLinSolSetPreconditioner(LS, ark_mem, NULL, NULL);
    if (retval != SUNLS_SUCCESS) {
      arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKLS",
                      "arkLSSetMassLinearSolver",
                      "Error in calling SUNLinSolSetPreconditioner");
      free(arkls_mem); arkls_mem = NULL;
      return(ARKLS_SUNLS_FAIL);
    }
  }

  /* When using a non-NULL SUNMatrix object, store pointer to M and create M_lu */
  if (M != NULL) {
    arkls_mem->M = M;
    arkls_mem->M_lu = SUNMatClone(M);
    if (arkls_mem->M_lu == NULL) {
      arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKLS",
                      "arkLSSetMassLinearSolver", MSG_LS_MEM_FAIL);
      free(arkls_mem); arkls_mem = NULL;
      return(ARKLS_MEM_FAIL);
    }
  }

  /* Allocate memory for x */
  arkls_mem->x = N_VClone(ark_mem->tempv1);
  if (arkls_mem->x == NULL) {
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKLS",
                    "arkLSSetMassLinearSolver", MSG_LS_MEM_FAIL);
    SUNMatDestroy(arkls_mem->M_lu);
    free(arkls_mem); arkls_mem = NULL;
    return(ARKLS_MEM_FAIL);
  }

  /* For iterative LS, compute sqrtN from a dot product */
  if ( (LSType == SUNLINEARSOLVER_ITERATIVE) ||
       (LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) ) {
    N_VConst(ONE, arkls_mem->x);
    arkls_mem->sqrtN = SUNRsqrt( N_VDotProd(arkls_mem->x,
                                            arkls_mem->x) );
  }

  /* Attach ARKLs interface to time stepper module */
  retval = ark_mem->step_attachmasssol(arkode_mem, arkLsMassInitialize,
                                       arkLsMassSetup, arkLsMTimes,
                                       arkLsMassSolve, arkLsMassFree,
                                       2, arkls_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKLS", "arkLSSetMassLinearSolver",
                    "Failed to attach to time stepper module");
    N_VDestroy(arkls_mem->x);
    SUNMatDestroy(arkls_mem->M_lu);
    free(arkls_mem); arkls_mem = NULL;
    return(retval);
  }

  return(ARKLS_SUCCESS);
}


/*===============================================================
  Optional input/output (called by time-stepper modules)
  ===============================================================*/

/*---------------------------------------------------------------
  arkLSSetJacFn specifies the Jacobian function.
  ---------------------------------------------------------------*/
int arkLSSetJacFn(void *arkode_mem, ARKLsJacFn jac)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSSetJacFn",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* return with failure if jac cannot be used */
  if ((jac != NULL) && (arkls_mem->A == NULL)) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLSSetJacFn",
                    "Jacobian routine cannot be supplied for NULL SUNMatrix");
    return(ARKLS_ILL_INPUT);
  }

  /* set Jacobian routine pointer, and update relevant flags */
  if (jac != NULL) {
    arkls_mem->jacDQ  = SUNFALSE;
    arkls_mem->jac    = jac;
    arkls_mem->J_data = ark_mem->user_data;
  } else {
    arkls_mem->jacDQ  = SUNTRUE;
    arkls_mem->jac    = arkLsDQJac;
    arkls_mem->J_data = ark_mem;
  }

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSSetMassFn specifies the mass matrix function.
  ---------------------------------------------------------------*/
int arkLSSetMassFn(void *arkode_mem, ARKLsMassFn mass)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSSetMassFn",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* return with failure if mass cannot be used */
  if (mass == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLSSetMassFn",
                    "Mass-matrix routine must be non-NULL");
    return(ARKLS_ILL_INPUT);
  }
  if (arkls_mem->M == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLSSetMassFn",
                    "Mass-matrix routine cannot be supplied for NULL SUNMatrix");
    return(ARKLS_ILL_INPUT);
  }

  /* set mass matrix routine pointer and return */
  arkls_mem->mass = mass;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSSetEpsLin specifies the nonlinear -> linear tolerance
  scale factor.
  ---------------------------------------------------------------*/
int arkLSSetEpsLin(void *arkode_mem, realtype eplifac)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; store input and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSSetEpsLin",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  arkls_mem->eplifac = (eplifac <= ZERO) ? ARKLS_EPLIN : eplifac;

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSSetMaxStepsBetweenJac specifies the maximum number of
  time steps to wait before recomputing the Jacobian matrix
  and/or preconditioner.
  ---------------------------------------------------------------*/
int arkLSSetMaxStepsBetweenJac(void *arkode_mem, long int msbj)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; store input and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSSetMaxStepsBetweenJac",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  arkls_mem->msbj = (msbj <= ZERO) ? ARKLS_MSBJ : msbj;

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSSetPreconditioner specifies the user-supplied
  preconditioner setup and solve routines.
  ---------------------------------------------------------------*/
int arkLSSetPreconditioner(void *arkode_mem,
                           ARKLsPrecSetupFn psetup,
                           ARKLsPrecSolveFn psolve)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  PSetupFn  arkls_psetup;
  PSolveFn  arkls_psolve;
  int       retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSSetPreconditioner",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* issue error if LS object does not allow user-supplied preconditioning */
  if (arkls_mem->LS->ops->setpreconditioner == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetPreconditioner",
                    "SUNLinearSolver object does not support user-supplied preconditioning");
    return(ARKLS_ILL_INPUT);
  }

  /* store function pointers for user-supplied routines */
  arkls_mem->pset   = psetup;
  arkls_mem->psolve = psolve;

  /* notify linear solver to call ARKLs interface routines */
  arkls_psetup = (psetup == NULL) ? NULL : arkLsPSetup;
  arkls_psolve = (psolve == NULL) ? NULL : arkLsPSolve;
  retval = SUNLinSolSetPreconditioner(arkls_mem->LS, ark_mem,
                                      arkls_psetup, arkls_psolve);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKLS",
                    "arkLSSetPreconditioner",
                    "Error in calling SUNLinSolSetPreconditioner");
    return(ARKLS_SUNLS_FAIL);
  }

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSSetJacTimes specifies the user-supplied Jacobian-vector
  product setup and multiply routines.
  ---------------------------------------------------------------*/
int arkLSSetJacTimes(void *arkode_mem,
                     ARKLsJacTimesSetupFn jtsetup,
                     ARKLsJacTimesVecFn jtimes)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSSetJacTimes",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* issue error if LS object does not allow user-supplied ATimes */
  if (arkls_mem->LS->ops->setatimes == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetJacTimes",
                    "SUNLinearSolver object does not support user-supplied ATimes routine");
    return(ARKLS_ILL_INPUT);
  }

  /* store function pointers for user-supplied routines in ARKLs
     interface (NULL jtimes implies use of DQ default) */
  if (jtimes != NULL) {
    arkls_mem->jtimesDQ = SUNFALSE;
    arkls_mem->jtsetup  = jtsetup;
    arkls_mem->jtimes   = jtimes;
    arkls_mem->Jt_data  = ark_mem->user_data;
  } else {
    arkls_mem->jtimesDQ = SUNTRUE;
    arkls_mem->jtsetup  = NULL;
    arkls_mem->jtimes   = arkLsDQJtimes;
    arkls_mem->Jt_data  = ark_mem;
  }

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetWorkSpace returns the length of workspace allocated for
  the ARKLS linear solver interface.
  ---------------------------------------------------------------*/
int arkLSGetWorkSpace(void *arkode_mem, long int *lenrw,
                      long int *leniw)
{
  ARKodeMem    ark_mem;
  ARKLsMem     arkls_mem;
  sunindextype lrw1, liw1;
  long int     lrw, liw;
  int          retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetWorkSpace",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* start with fixed sizes plus vector/matrix pointers */
  *lenrw = 3;
  *leniw = 30;

  /* add NVector sizes */
  if (arkls_mem->x->ops->nvspace) {
    N_VSpace(arkls_mem->x, &lrw1, &liw1);
    *lenrw += 2*lrw1;
    *leniw += 2*liw1;
  }

  /* add SUNMatrix size (only account for the one owned by Ls interface) */
  if (arkls_mem->savedJ)
    if (arkls_mem->savedJ->ops->space) {
      retval = SUNMatSpace(arkls_mem->savedJ, &lrw, &liw);
      if (retval == 0) {
        *lenrw += lrw;
        *leniw += liw;
      }
    }

  /* add LS sizes */
  if (arkls_mem->LS->ops->space) {
    retval = SUNLinSolSpace(arkls_mem->LS, &lrw, &liw);
    if (retval == SUNLS_SUCCESS) {
      *lenrw += lrw;
      *leniw += liw;
    }
  }

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumJacEvals returns the number of Jacobian evaluations
  ---------------------------------------------------------------*/
int arkLSGetNumJacEvals(void *arkode_mem, long int *njevals)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; set output value and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetNumJacEvals",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *njevals = arkls_mem->nje;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumRhsEvals returns the number of calls to the ODE
  function needed for the DQ Jacobian approximation or J*v product
  approximation.
  ---------------------------------------------------------------*/
int arkLSGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; set output value and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetNumRhsEvals",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *nfevalsLS = arkls_mem->nfeDQ;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumPrecEvals returns the number of calls to the
  user- or ARKode-supplied preconditioner setup routine.
  ---------------------------------------------------------------*/
int arkLSGetNumPrecEvals(void *arkode_mem, long int *npevals)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; set output value and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetNumPrecEvals",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *npevals = arkls_mem->npe;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumPrecSolves returns the number of calls to the
  user- or ARKode-supplied preconditioner solve routine.
  ---------------------------------------------------------------*/
int arkLSGetNumPrecSolves(void *arkode_mem, long int *npsolves)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; set output value and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetNumPrecSolves",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *npsolves = arkls_mem->nps;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumLinIters returns the number of linear iterations
  (if accessible from the LS object).
  ---------------------------------------------------------------*/
int arkLSGetNumLinIters(void *arkode_mem, long int *nliters)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; set output value and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetNumLinIters",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *nliters = arkls_mem->nli;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumConvFails returns the number of linear solver
  convergence failures (as reported by the LS object).
  ---------------------------------------------------------------*/
int arkLSGetNumConvFails(void *arkode_mem, long int *nlcfails)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; set output value and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetNumConvFails",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *nlcfails = arkls_mem->ncfl;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumJTSetupEvals returns the number of calls to the
  user-supplied Jacobian-vector product setup routine.
  ---------------------------------------------------------------*/
int arkLSGetNumJTSetupEvals(void *arkode_mem, long int *njtsetups)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; set output value and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetNumJTSetupEvals",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *njtsetups = arkls_mem->njtsetup;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumJtimesEvals returns the number of calls to the
  Jacobian-vector product multiply routine.
  ---------------------------------------------------------------*/
int arkLSGetNumJtimesEvals(void *arkode_mem, long int *njvevals)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; set output value and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetNumJtimesEvals",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *njvevals = arkls_mem->njtimes;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetLastFlag returns the last flag set in a ARKLS
  function.
  ---------------------------------------------------------------*/
int arkLSGetLastFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  int       retval;

  /* access ARKLsMem structure; set output value and return */
  retval = arkLs_AccessLMem(arkode_mem, "arkLSGetLastFlag",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *flag = arkls_mem->last_flag;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetReturnFlagName translates from the integer error code
  returned by an ARKLs routine to the corresponding string
  equivalent for that flag
  ---------------------------------------------------------------*/
char *arkLSGetReturnFlagName(long int flag)
{
  char *name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case ARKLS_SUCCESS:
    sprintf(name,"ARKLS_SUCCESS");
    break;
  case ARKLS_MEM_NULL:
    sprintf(name,"ARKLS_MEM_NULL");
    break;
  case ARKLS_LMEM_NULL:
    sprintf(name,"ARKLS_LMEM_NULL");
    break;
  case ARKLS_ILL_INPUT:
    sprintf(name,"ARKLS_ILL_INPUT");
    break;
  case ARKLS_MEM_FAIL:
    sprintf(name,"ARKLS_MEM_FAIL");
    break;
  case ARKLS_MASSMEM_NULL:
    sprintf(name,"ARKLS_MASSMEM_NULL");
    break;
  case ARKLS_JACFUNC_UNRECVR:
    sprintf(name,"ARKLS_JACFUNC_UNRECVR");
    break;
  case ARKLS_JACFUNC_RECVR:
    sprintf(name,"ARKLS_JACFUNC_RECVR");
    break;
  case ARKLS_MASSFUNC_UNRECVR:
    sprintf(name,"ARKLS_MASSFUNC_UNRECVR");
    break;
  case ARKLS_MASSFUNC_RECVR:
    sprintf(name,"ARKLS_MASSFUNC_RECVR");
    break;
  case ARKLS_SUNMAT_FAIL:
    sprintf(name,"ARKLS_SUNMAT_FAIL");
    break;
  case ARKLS_SUNLS_FAIL:
    sprintf(name,"ARKLS_SUNLS_FAIL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*---------------------------------------------------------------
  arkLSSetMassEpsLin specifies the nonlinear -> linear tolerance
  scale factor for mass matrix linear systems.
  ---------------------------------------------------------------*/
int arkLSSetMassEpsLin(void *arkode_mem, realtype eplifac)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; store input and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSSetMassEpsLin",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  arkls_mem->eplifac = (eplifac <= ZERO) ? ARKLS_EPLIN : eplifac;

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSSetMassPreconditioner specifies the user-supplied
  preconditioner setup and solve routines.
  ---------------------------------------------------------------*/
int arkLSSetMassPreconditioner(void *arkode_mem,
                               ARKLsMassPrecSetupFn psetup,
                               ARKLsMassPrecSolveFn psolve)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  PSetupFn     arkls_mpsetup;
  PSolveFn     arkls_mpsolve;
  int          retval;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSSetMassPreconditioner",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* issue error if LS object does not allow user-supplied preconditioning */
  if (arkls_mem->LS->ops->setpreconditioner == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetMassPreconditioner",
                    "SUNLinearSolver object does not support user-supplied preconditioning");
    return(ARKLS_ILL_INPUT);
  }

  /* store function pointers for user-supplied routines in ARKLs interface */
  arkls_mem->pset   = psetup;
  arkls_mem->psolve = psolve;

  /* notify linear solver to call ARKLs interface routines */
  arkls_mpsetup = (psetup == NULL) ? NULL : arkLsMPSetup;
  arkls_mpsolve = (psolve == NULL) ? NULL : arkLsMPSolve;
  retval = SUNLinSolSetPreconditioner(arkls_mem->LS, ark_mem,
                                      arkls_mpsetup, arkls_mpsolve);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKLS",
                    "arkLSSetMassPreconditioner",
                    "Error in calling SUNLinSolSetPreconditioner");
    return(ARKLS_SUNLS_FAIL);
  }

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSSetMassTimes specifies the user-supplied mass
  matrix-vector product setup and multiply routines.
  ---------------------------------------------------------------*/
int arkLSSetMassTimes(void *arkode_mem,
                      ARKLsMassTimesSetupFn mtsetup,
                      ARKLsMassTimesVecFn mtimes,
                      void *mtimes_data)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSSetMassTimes",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* issue error if mtimes function is unusable */
  if (mtimes == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetMassTimes",
                    "non-NULL mtimes function must be supplied");
    return(ARKLS_ILL_INPUT);
  }

  /* issue error if LS object does not allow user-supplied ATimes */
  if (arkls_mem->LS->ops->setatimes == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLSSetMassTimes",
                    "SUNLinearSolver object does not support user-supplied ATimes routine");
    return(ARKLS_ILL_INPUT);
  }

  /* store pointers for user-supplied routines and data structure
     in ARKLs interface */
  arkls_mem->mtsetup = mtsetup;
  arkls_mem->mtimes  = mtimes;
  arkls_mem->mt_data = mtimes_data;

  /* notify linear solver to call ARKLs interface routine */
  retval = SUNLinSolSetATimes(arkls_mem->LS, ark_mem, arkLsMTimes);
  if (retval != SUNLS_SUCCESS) {
    arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKLS",
                    "arkLSSetMassTimes",
                    "Error in calling SUNLinSolSetATimes");
    return(ARKLS_SUNLS_FAIL);
  }

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetMassWorkSpace
  ---------------------------------------------------------------*/
int arkLSGetMassWorkSpace(void *arkode_mem, long int *lenrw,
                          long int *leniw)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  sunindextype lrw1, liw1;
  long int     lrw, liw;
  int          retval;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetMassWorkSpace",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* start with fixed sizes plus vector/matrix pointers */
  *lenrw = 2;
  *leniw = 23;

  /* add NVector sizes */
  if (ark_mem->tempv1->ops->nvspace) {
    N_VSpace(ark_mem->tempv1, &lrw1, &liw1);
    *lenrw += lrw1;
    *leniw += liw1;
  }

  /* add SUNMatrix size (only account for the one owned by Ls interface) */
  if (arkls_mem->M_lu)
    if (arkls_mem->M_lu->ops->space) {
      retval = SUNMatSpace(arkls_mem->M_lu, &lrw, &liw);
      if (retval == 0) {
        *lenrw += lrw;
        *leniw += liw;
      }
    }

  /* add LS sizes */
  if (arkls_mem->LS->ops->space) {
    retval = SUNLinSolSpace(arkls_mem->LS, &lrw, &liw);
    if (retval == SUNLS_SUCCESS) {
      *lenrw += lrw;
      *leniw += liw;
    }
  }

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumMassSetups returns the number of mass matrix
  solver 'setup' calls
  ---------------------------------------------------------------*/
int arkLSGetNumMassSetups(void *arkode_mem, long int *nmsetups)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; set output value and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetNumMassSetups",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *nmsetups = arkls_mem->nmsetups;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumMassMult returns the number of calls to the user-
  supplied or internal mass matrix-vector product multiply routine.
  ---------------------------------------------------------------*/
int arkLSGetNumMassMult(void *arkode_mem, long int *nmvevals)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; set output value and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetNumMassMult",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *nmvevals = arkls_mem->nmtimes;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumMassSolves returns the number of mass matrix
  solver 'solve' calls
  ---------------------------------------------------------------*/
int arkLSGetNumMassSolves(void *arkode_mem, long int *nmsolves)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; set output value and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetNumMassSolves",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *nmsolves = arkls_mem->nmsolves;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumMassPrecEvals returns the number of calls to the
  user- or ARKode-supplied preconditioner setup routine.
  ---------------------------------------------------------------*/
int arkLSGetNumMassPrecEvals(void *arkode_mem, long int *npevals)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; set output value and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetNumMassPrecEvals",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *npevals = arkls_mem->npe;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumMassPrecSolves returns the number of calls to the
  user- or ARKode-supplied preconditioner solve routine.
  ---------------------------------------------------------------*/
int arkLSGetNumMassPrecSolves(void *arkode_mem, long int *npsolves)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; set output value and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetNumMassPrecSolves",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *npsolves = arkls_mem->nps;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumMassIters returns the number of mass matrix solver
  linear iterations (if accessible from the LS object).
  ---------------------------------------------------------------*/
int arkLSGetNumMassIters(void *arkode_mem, long int *nmiters)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; set output value and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetNumMassIters",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *nmiters = arkls_mem->nli;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumMassConvFails returns the number of linear solver
  convergence failures (as reported by the LS object).
  ---------------------------------------------------------------*/
int arkLSGetNumMassConvFails(void *arkode_mem, long int *nmcfails)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; set output value and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetNumMassConvFails",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *nmcfails = arkls_mem->ncfl;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetNumMTSetups returns the number of calls to the
  user-supplied mass matrix-vector product setup routine.
  ---------------------------------------------------------------*/
int arkLSGetNumMTSetups(void *arkode_mem, long int *nmtsetups)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; set output value and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetNumMTSetups",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *nmtsetups = arkls_mem->nmtsetup;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLSGetLastMassFlag returns the last flag set in a ARKLS
  function.
  ---------------------------------------------------------------*/
int arkLSGetLastMassFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure; set output value and return */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLSGetLastMassFlag",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);
  *flag = arkls_mem->last_flag;
  return(ARKLS_SUCCESS);
}


/*===============================================================
  ARKLS Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkLsATimes:

  This routine generates the matrix-vector product z = Av, where
  A = M - gamma*J. The product M*v is obtained either by calling
  the mtimes routine or by just using v (if M=I).  The product
  J*v is obtained by calling the jtimes routine. It is then scaled
  by -gamma and added to M*v to obtain A*v. The return value is
  the same as the values returned by jtimes and mtimes --
  0 if successful, nonzero otherwise.
  ---------------------------------------------------------------*/
int arkLsATimes(void *arkode_mem, N_Vector v, N_Vector z)
{
  ARKodeMem   ark_mem;
  ARKLsMem    arkls_mem;
  void*       ark_step_massmem;
  int         retval;
  realtype    gamma, gamrat;
  booleantype dgamma_fail, *jcur;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLsATimes",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Access mass matrix solver (if it exists) */
  ark_step_massmem = NULL;
  if (ark_mem->step_getmassmem != NULL)
    ark_step_massmem = ark_mem->step_getmassmem(arkode_mem);

  /* get gamma values from time step module */
  retval = ark_mem->step_getgammas(arkode_mem, &gamma, &gamrat,
                                   &jcur, &dgamma_fail);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKLS", "arkLsATimes",
                    "An error occurred in ark_step_getgammas");
    return(retval);
  }

  /* call Jacobian-times-vector product routine
     (either user-supplied or internal DQ) */
  retval = arkls_mem->jtimes(v, z,
                             arkls_mem->tcur,
                             arkls_mem->ycur,
                             arkls_mem->fcur,
                             arkls_mem->Jt_data,
                             arkls_mem->ytemp);
  arkls_mem->njtimes++;
  if (retval != 0) return(retval);

  /* Compute mass matrix vector product and add to result */
  if (ark_step_massmem != NULL) {
    retval = arkLsMTimes(arkode_mem, v, arkls_mem->ytemp);
    if (retval != 0) return(retval);
    N_VLinearSum(ONE, arkls_mem->ytemp, -gamma, z, z);
  } else {
    N_VLinearSum(ONE, v, -gamma, z, z);
  }

  return(0);
}

/*---------------------------------------------------------------
  arkLsPSetup:

  This routine interfaces between the generic iterative linear
  solvers and the user's psetup routine.  It passes to psetup all
  required state information from arkode_mem.  Its return value
  is the same as that returned by psetup. Note that the generic
  iterative linear solvers guarantee that arkLsPSetup will only
  be called in the case that the user's psetup routine is non-NULL.
  ---------------------------------------------------------------*/
int arkLsPSetup(void *arkode_mem)
{
  ARKodeMem   ark_mem;
  ARKLsMem    arkls_mem;
  realtype    gamma, gamrat;
  booleantype dgamma_fail, *jcur;
  int         retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLsPSetup",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get gamma values from time step module */
  retval = ark_mem->step_getgammas(arkode_mem, &gamma, &gamrat,
                                   &jcur, &dgamma_fail);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKLS", "arkLsPSetup",
                    "An error occurred in ark_step_getgammas");
    return(retval);
  }

  /* Call user pset routine to update preconditioner and possibly
     reset jcur (pass !jbad as update suggestion) */
  retval = arkls_mem->pset(arkls_mem->tcur,
                           arkls_mem->ycur,
                           arkls_mem->fcur,
                           !(arkls_mem->jbad),
                           jcur, gamma,
                           arkls_mem->P_data);
  return(retval);
}

/*---------------------------------------------------------------
  arkLsPSolve:

  This routine interfaces between the generic SUNLinSolSolve
  routine and the user's psolve routine.  It passes to psolve all
  required state information from arkode_mem.  Its return value
  is the same as that returned by psolve. Note that the generic
  SUNLinSol solver guarantees that arkLsPSolve will not be
  called in the case in which preconditioning is not done. This
  is the only case in which the user's psolve routine is allowed
  to be NULL.
  ---------------------------------------------------------------*/
int arkLsPSolve(void *arkode_mem, N_Vector r, N_Vector z,
                realtype tol, int lr)
{
  ARKodeMem   ark_mem;
  ARKLsMem    arkls_mem;
  realtype    gamma, gamrat;
  booleantype dgamma_fail, *jcur;
  int         retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLsPSolve",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get gamma values from time step module */
  retval = ark_mem->step_getgammas(arkode_mem, &gamma, &gamrat,
                                   &jcur, &dgamma_fail);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKLS", "arkLsPSolve",
                    "An error occurred in ark_step_getgammas");
    return(retval);
  }

  /* call the user-supplied psolve routine, and accumulate count */
  retval = arkls_mem->psolve(arkls_mem->tcur,
                             arkls_mem->ycur,
                             arkls_mem->fcur, r, z,
                             gamma, tol, lr,
                             arkls_mem->P_data);
  arkls_mem->nps++;
  return(retval);
}

/*---------------------------------------------------------------
  arkLsMTimes:

  This routine generates the matrix-vector product z = Mv, where
  M is the system mass matrix, by calling the user-supplied mtimes
  routine. The return value is the same as the value returned
  by mtimes -- 0 if successful, nonzero otherwise.
  ---------------------------------------------------------------*/
int arkLsMTimes(void *arkode_mem, N_Vector v, N_Vector z)
{
  ARKodeMem       ark_mem;
  ARKLsMassMem    arkls_mem;
  int             retval;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLsMTimes",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* perform multiply by either calling the user-supplied routine
     (default), or asking the SUNMatrix to do the multiply */
  retval = -1;
  if (arkls_mem->mtimes) {

    /* call user-supplied mtimes routine and increment counter */
    retval = arkls_mem->mtimes(v, z, ark_mem->tcur,
                               arkls_mem->mt_data);

  } else if (arkls_mem->M) {

    if (arkls_mem->M->ops->matvec)
      retval = SUNMatMatvec(arkls_mem->M, v, z);

  }

  if (retval == 0) {
    arkls_mem->nmtimes++;
  } else {
    arkProcessError(ark_mem, retval, "ARKLS", "arkLsMTimes",
                    "Missing mass matrix-vector product routine");
  }
  return(retval);
}


/*---------------------------------------------------------------
  arkLsMPSetup:

  This routine interfaces between the generic linear solver and
  the user's mass matrix psetup routine.  It passes to psetup all
  required state information from arkode_mem.  Its return value
  is the same as that returned by psetup.  Note that the generic
  linear solvers guarantee that arkLsMPSetup will only be
  called if the user's psetup routine is non-NULL.
  ---------------------------------------------------------------*/
int arkLsMPSetup(void *arkode_mem)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLsMPSetup",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* only proceed if the mass matrix is time-independent or if
     pset has not been called previously */
  if (!arkls_mem->time_dependent && arkls_mem->npe)
    return(0);

  /* call user-supplied pset routine and increment counter */
  retval = arkls_mem->pset(ark_mem->tcur, arkls_mem->P_data);
  arkls_mem->npe++;
  return(retval);
}


/*---------------------------------------------------------------
  arkLsMPSolve:

  This routine interfaces between the generic LS routine and the
  user's mass matrix psolve routine.  It passes to psolve all
  required state information from arkode_mem.  Its return value is
  the same as that returned by psolve. Note that the generic
  solver guarantees that arkLsMPSolve will not be called in the
  case in which preconditioning is not done. This is the only case
  in which the user's psolve routine is allowed to be NULL.
  ---------------------------------------------------------------*/
int arkLsMPSolve(void *arkode_mem, N_Vector r, N_Vector z,
                 realtype tol, int lr)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLsMPSolve",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* call the user-supplied psolve routine, and accumulate count */
  retval = arkls_mem->psolve(ark_mem->tcur, r, z, tol,
                             lr, arkls_mem->P_data);
  arkls_mem->nps++;
  return(retval);
}


/*---------------------------------------------------------------
  arkLsDQJac:

  This routine is a wrapper for the Dense and Band
  implementations of the difference quotient Jacobian
  approximation routines.
  ---------------------------------------------------------------*/
int arkLsDQJac(realtype t, N_Vector y, N_Vector fy,
               SUNMatrix Jac, void *arkode_mem, N_Vector tmp1,
               N_Vector tmp2, N_Vector tmp3)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  ARKRhsFn  fi;
  int       retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLsDQJac",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* verify that Jac is non-NULL */
  if (Jac == NULL) {
    arkProcessError(ark_mem, ARKLS_LMEM_NULL, "ARKLS",
                    "arkLsDQJac", "SUNMatrix is NULL");
    return(ARKLS_LMEM_NULL);
  }

  /* Access implicit RHS function */
  fi = ark_mem->step_getimplicitrhs((void*) ark_mem);
  if (fi == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLsDQJac",
                    "Time step module is missing implicit RHS fcn");
    return(ARKLS_ILL_INPUT);
  }

  /* Verify that N_Vector supports required routines */
  if (ark_mem->tempv1->ops->nvcloneempty == NULL ||
      ark_mem->tempv1->ops->nvwrmsnorm == NULL ||
      ark_mem->tempv1->ops->nvlinearsum == NULL ||
      ark_mem->tempv1->ops->nvdestroy == NULL ||
      ark_mem->tempv1->ops->nvscale == NULL ||
      ark_mem->tempv1->ops->nvgetarraypointer == NULL ||
      ark_mem->tempv1->ops->nvsetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLsDQJac", MSG_LS_BAD_NVECTOR);
    return(ARKLS_ILL_INPUT);
  }

  /* Call the matrix-structure-specific DQ approximation routine */
  if (SUNMatGetID(Jac) == SUNMATRIX_DENSE) {
    retval = arkLsDenseDQJac(t, y, fy, Jac, ark_mem, arkls_mem,
                             fi, tmp1);
  } else if (SUNMatGetID(Jac) == SUNMATRIX_BAND) {
    retval = arkLsBandDQJac(t, y, fy, Jac, ark_mem, arkls_mem,
                            fi, tmp1, tmp2);
  } else {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLsDQJac",
                    "arkLsDQJac not implemented for this SUNMatrix type!");
    retval = ARKLS_ILL_INPUT;
  }
  return(retval);
}

/*---------------------------------------------------------------
  arkLsDenseDQJac:

  This routine generates a dense difference quotient approximation
  to the Jacobian of f(t,y). It assumes a dense SUNMatrix input
  (stored column-wise, and that elements within each column are
  contiguous). The address of the jth column of J is obtained via
  the function SUNDenseMatrix_Column() and this pointer is
  associated with an N_Vector using the
  N_VGetArrayPointer/N_VSetArrayPointer functions.  Finally, the
  actual computation of the jth column of the Jacobian is done
  with a call to N_VLinearSum.
  ---------------------------------------------------------------*/
int arkLsDenseDQJac(realtype t, N_Vector y, N_Vector fy,
                    SUNMatrix Jac, ARKodeMem ark_mem,
                    ARKLsMem arkls_mem, ARKRhsFn fi,
                    N_Vector tmp1)
{
  realtype     fnorm, minInc, inc, inc_inv, yjsaved, srur;
  realtype    *y_data, *ewt_data;
  N_Vector     ftemp, jthCol;
  sunindextype j, N;
  int          retval = 0;

  /* access matrix dimension */
  N = SUNDenseMatrix_Rows(Jac);

  /* Rename work vector for readibility */
  ftemp = tmp1;

  /* Create an empty vector for matrix column calculations */
  jthCol = N_VCloneEmpty(tmp1);

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(ark_mem->ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(ark_mem->uround);
  fnorm = N_VWrmsNorm(fy, ark_mem->rwt);
  minInc = (fnorm != ZERO) ?
    (MIN_INC_MULT * SUNRabs(ark_mem->h) * ark_mem->uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */
    N_VSetArrayPointer(SUNDenseMatrix_Column(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = SUNMAX(srur*SUNRabs(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = fi(t, y, ftemp, ark_mem->user_data);
    arkls_mem->nfeDQ++;
    if (retval != 0) break;

    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

  }

  /* Destroy jthCol vector */
  N_VSetArrayPointer(NULL, jthCol);  /* SHOULDN'T BE NEEDED */
  N_VDestroy(jthCol);

  return(retval);
}


/*---------------------------------------------------------------
  arkLsBandDQJac:

  This routine generates a banded difference quotient approximation
  to the Jacobian of f(t,y).  It assumes a band SUNMatrix input
  (stored column-wise, and that elements within each column are
  contiguous). This makes it possible to get the address
  of a column of J via the function SUNBandMatrix_Column() and to
  write a simple for loop to set each of the elements of a column
  in succession.
  ---------------------------------------------------------------*/
int arkLsBandDQJac(realtype t, N_Vector y, N_Vector fy,
                   SUNMatrix Jac, ARKodeMem ark_mem,
                   ARKLsMem arkls_mem, ARKRhsFn fi,
                   N_Vector tmp1, N_Vector tmp2)
{
  N_Vector     ftemp, ytemp;
  realtype     fnorm, minInc, inc, inc_inv, srur;
  realtype    *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  sunindextype group, i, j, width, ngroups, i1, i2;
  sunindextype N, mupper, mlower;
  int          retval = 0;

  /* access matrix dimensions */
  N = SUNBandMatrix_Columns(Jac);
  mupper = SUNBandMatrix_UpperBandwidth(Jac);
  mlower = SUNBandMatrix_LowerBandwidth(Jac);

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(ark_mem->ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(ark_mem->uround);
  fnorm = N_VWrmsNorm(fy, ark_mem->rwt);
  minInc = (fnorm != ZERO) ?
    (MIN_INC_MULT * SUNRabs(ark_mem->h) * ark_mem->uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = SUNMIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {

    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */
    retval = fi(ark_mem->tcur, ytemp, ftemp, ark_mem->user_data);
    arkls_mem->nfeDQ++;
    if (retval != 0) break;

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = SUNBandMatrix_Column(Jac, j);
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mupper);
      i2 = SUNMIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        SM_COLUMN_ELEMENT_B(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

  return(retval);
}


/*---------------------------------------------------------------
  arkLsDQJtimes:

  This routine generates a difference quotient approximation to
  the Jacobian-vector product fi_y(t,y) * v. The approximation is
  Jv = [fi(y + v*sig) - fi(y)]/sig, where sig = 1 / ||v||_WRMS,
  i.e. the WRMS norm of v*sig is 1.
  ---------------------------------------------------------------*/
int arkLsDQJtimes(N_Vector v, N_Vector Jv, realtype t,
                  N_Vector y, N_Vector fy, void *arkode_mem,
                  N_Vector work)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  ARKRhsFn  fi;
  realtype  sig, siginv;
  int       iter, retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLsDQJtimes",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Initialize perturbation to 1/||v|| */
  sig = ONE/N_VWrmsNorm(v, ark_mem->ewt);

  /* Access implicit RHS function */
  fi = ark_mem->step_getimplicitrhs(arkode_mem);
  if (fi == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLsDQJtimes",
                    "Time step module is missing implicit RHS fcn");
    return(ARKLS_ILL_INPUT);
  }

  for (iter=0; iter<MAX_DQITERS; iter++) {

    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* Set Jv = f(tn, y+sig*v) */
    retval = fi(t, work, Jv, ark_mem->user_data);
    arkls_mem->nfeDQ++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    /* If fi failed recoverably, shrink sig and retry */
    sig *= PT25;

  }

  /* If retval still isn't 0, return with a recoverable failure */
  if (retval > 0) return(+1);

  /* Replace Jv by (Jv - fy)/sig */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, fy, Jv);

  return(0);
}


/*---------------------------------------------------------------
  arkLsInitialize performs remaining initializations specific
  to the iterative linear solver interface (and solver itself)
  ---------------------------------------------------------------*/
int arkLsInitialize(void* arkode_mem)
{
  ARKodeMem    ark_mem;
  ARKLsMem     arkls_mem;
  ARKLsMassMem arkls_massmem;
  int          retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLsInitialize",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKLsMassMem (if applicable) */
  arkls_massmem = NULL;
  if (ark_mem->step_getmassmem != NULL)
    if (ark_mem->step_getmassmem(arkode_mem) != NULL) {
      retval = arkLs_AccessMassMem(arkode_mem, "arkLsInitialize",
                                   &ark_mem, &arkls_massmem);
      if (retval != ARK_SUCCESS)  return(retval);
    }


  /* Test for valid combinations of matrix & Jacobian routines: */
  if (arkls_mem->A == NULL) {

    /* If SUNMatrix A is NULL: ensure 'jac' function pointer is still NULL */
    arkls_mem->jacDQ  = SUNFALSE;
    arkls_mem->jac    = NULL;
    arkls_mem->J_data = NULL;

  } else if (arkls_mem->jacDQ) {

    /* If A is non-NULL, and 'jac' is not user-supplied:
       - if A is dense or band, ensure that our DQ approx. is used
       - otherwise => error */
    retval = 0;
    if (arkls_mem->A->ops->getid) {

      if ( (SUNMatGetID(arkls_mem->A) == SUNMATRIX_DENSE) ||
           (SUNMatGetID(arkls_mem->A) == SUNMATRIX_BAND) ) {
        arkls_mem->jac    = arkLsDQJac;
        arkls_mem->J_data = ark_mem;
      } else {
        retval++;
      }

    } else {
      retval++;
    }
    if (retval) {
      arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLsInitialize",
                      "No Jacobian constructor available for SUNMatrix type");
      arkls_mem->last_flag = ARKLS_ILL_INPUT;
      return(ARKLS_ILL_INPUT);
    }

  } else {

    /* If A is non-NULL, and 'jac' is user-supplied,
       reset J_data pointer (just in case) */
    arkls_mem->J_data = ark_mem->user_data;
  }


  /* Test for valid combination of system matrix and mass matrix (if applicable) */
  if (arkls_massmem) {

    /* A and M must both be NULL or non-NULL */
    if ( (arkls_mem->A==NULL) ^ (arkls_massmem->M==NULL) ) {
      arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLsInitialize",
                      "Cannot combine NULL and non-NULL System and mass matrices");
      arkls_mem->last_flag = ARKLS_ILL_INPUT;
      return(ARKLS_ILL_INPUT);
    }

    /* If A is non-NULL, A and M must have matching types (if accessible) */
    if (arkls_mem->A) {
      retval = 0;
      if ((arkls_mem->A->ops->getid==NULL) ^ (arkls_massmem->M->ops->getid==NULL))
        retval++;
      if (arkls_mem->A->ops->getid)
        if (SUNMatGetID(arkls_mem->A) != SUNMatGetID(arkls_massmem->M))
          retval++;
      if (retval) {
        arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLsInitialize",
                        "System and mass matrices have incompatible types");
        arkls_mem->last_flag = ARKLS_ILL_INPUT;
        return(ARKLS_ILL_INPUT);
      }
    }

    /* initialize mass matrix linear solver  */
    retval = arkLsMassInitialize(arkode_mem);
    if (retval != ARKLS_SUCCESS) {
      arkls_mem->last_flag = retval;
      return(retval);
    }
  }

  /* reset counters */
  arkLsInitializeCounters(arkls_mem);

  /* Set Jacobian-vector product fields, based on jtimesDQ */
  if (arkls_mem->jtimesDQ) {
    arkls_mem->jtsetup = NULL;
    arkls_mem->jtimes  = arkLsDQJtimes;
    arkls_mem->Jt_data = ark_mem;
  } else {
    arkls_mem->Jt_data = ark_mem->user_data;
  }

  /* if A is NULL and psetup is not present, then arkLsSetup does
     not need to be called, so set the lsetup function to NULL (if possible) */
  if ( (arkls_mem->A == NULL) &&
       (arkls_mem->pset == NULL) &&
       (ark_mem->step_disablelsetup != NULL) )
    ark_mem->step_disablelsetup(arkode_mem);

  /* Call LS initialize routine, and return result */
  arkls_mem->last_flag = SUNLinSolInitialize(arkls_mem->LS);
  return(arkls_mem->last_flag);
}


/*---------------------------------------------------------------
  arkLsSetup conditionally calls the LS 'setup' routine.

  When using a SUNMatrix object, this determines whether
  to update a Jacobian matrix (or use a stored version), based
  on heuristics regarding previous convergence issues, the number
  of time steps since it was last updated, etc.; it then creates
  the system matrix from this, the 'gamma' factor and the
  mass/identity matrix,
  A = M-gamma*J.

  This routine then calls the LS 'setup' routine with A.
  ---------------------------------------------------------------*/
int arkLsSetup(void* arkode_mem, int convfail, realtype tpred,
               N_Vector ypred, N_Vector fpred, booleantype *jcurPtr,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  ARKodeMem    ark_mem;
  ARKLsMem     arkls_mem;
  ARKLsMassMem arkls_massmem;
  void*        ark_step_massmem;
  realtype     gamma, gamrat;
  booleantype  dgamma_fail, *jcur;
  int          retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLsInitialize",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Set ARKLs time and N_Vector pointers to current time,
     solution and rhs */
  arkls_mem->tcur = tpred;
  arkls_mem->ycur = ypred;
  arkls_mem->fcur = fpred;

  /* get gamma values from time step module */
  arkls_mem->last_flag = ark_mem->step_getgammas(arkode_mem, &gamma, &gamrat,
                                                 &jcur, &dgamma_fail);
  if (arkls_mem->last_flag) {
    arkProcessError(ark_mem, retval, "ARKLS", "arkLsSetup",
                    "An error occurred in ark_step_getgammas");
    return(arkls_mem->last_flag);
  }

  /* Use nst, gamma/gammap, and convfail to set J/P eval. flag jok;
     Note: the "ARK_FAIL_BAD_J" test is asking whether the nonlinear
     solver converged due to a bad system Jacobian AND our gamma was
     fine, indicating that the J and/or P were invalid */
  arkls_mem->jbad = (ark_mem->nst == 0) ||
    (ark_mem->nst > arkls_mem->nstlj + arkls_mem->msbj) ||
    ((convfail == ARK_FAIL_BAD_J) && (!dgamma_fail)) ||
    (convfail == ARK_FAIL_OTHER);

  /* If using a NULL SUNMatrix, set jcur to jbad; otherwise update J as appropriate */
  if (arkls_mem->A == NULL) {

    *jcurPtr = arkls_mem->jbad;

  } else {

    /* If jbad = SUNFALSE, use saved copy of J */
    if (!arkls_mem->jbad) {

      *jcurPtr = SUNFALSE;
      retval = SUNMatCopy(arkls_mem->savedJ, arkls_mem->A);
      if (retval) {
        arkProcessError(ark_mem, ARKLS_SUNMAT_FAIL, "ARKLS",
                        "arkLsSetup",  MSG_LS_SUNMAT_FAILED);
        arkls_mem->last_flag = ARKLS_SUNMAT_FAIL;
        return(arkls_mem->last_flag);
      }

    /* If jbad = SUNTRUE, clear out J and call jac routine for new value */
    } else {

      arkls_mem->nje++;
      arkls_mem->nstlj = ark_mem->nst;
      *jcurPtr = SUNTRUE;
      retval = SUNMatZero(arkls_mem->A);
      if (retval) {
        arkProcessError(ark_mem, ARKLS_SUNMAT_FAIL, "ARKLS",
                        "arkLsSetup",  MSG_LS_SUNMAT_FAILED);
        arkls_mem->last_flag = ARKLS_SUNMAT_FAIL;
        return(arkls_mem->last_flag);
      }

      retval = arkls_mem->jac(tpred, ypred, fpred, arkls_mem->A,
                              arkls_mem->J_data, vtemp1, vtemp2, vtemp3);
      if (retval < 0) {
        arkProcessError(ark_mem, ARKLS_JACFUNC_UNRECVR, "ARKLS",
                        "arkLsSetup",  MSG_LS_JACFUNC_FAILED);
        arkls_mem->last_flag = ARKLS_JACFUNC_UNRECVR;
        return(-1);
      }
      if (retval > 0) {
        arkls_mem->last_flag = ARKLS_JACFUNC_RECVR;
        return(1);
      }

      retval = SUNMatCopy(arkls_mem->A, arkls_mem->savedJ);
      if (retval) {
        arkProcessError(ark_mem, ARKLS_SUNMAT_FAIL, "ARKLS",
                        "arkLsSetup",  MSG_LS_SUNMAT_FAILED);
        arkls_mem->last_flag = ARKLS_SUNMAT_FAIL;
        return(arkls_mem->last_flag);
      }

    }

    /* Scale and add mass matrix to get A = M-gamma*J*/
    ark_step_massmem = NULL;
    if (ark_mem->step_getmassmem)
      ark_step_massmem = ark_mem->step_getmassmem(arkode_mem);
    if (ark_step_massmem) {

      arkls_massmem = (ARKLsMassMem) ark_step_massmem;

      /* Setup mass matrix linear solver (including recomputation of mass matrix) */
      arkls_mem->last_flag = arkLsMassSetup(arkode_mem, vtemp1, vtemp2, vtemp3);
      if (retval) {
        arkProcessError(ark_mem, ARKLS_SUNMAT_FAIL, "ARKLS", "arkLsSetup",
                        "Error setting up mass-matrix linear solver");
        return(arkls_mem->last_flag);
      }

      /* Perform linear combination A = M-gamma*A */
      retval = SUNMatScaleAdd(-gamma, arkls_mem->A, arkls_massmem->M);

      /* or if M==I, set A = I-gamma*J*/
    } else {
      retval = SUNMatScaleAddI(-gamma, arkls_mem->A);
    }
    if (retval) {
      arkProcessError(ark_mem, ARKLS_SUNMAT_FAIL, "ARKLS",
                      "arkLsSetup",  MSG_LS_SUNMAT_FAILED);
      arkls_mem->last_flag = ARKLS_SUNMAT_FAIL;
      return(arkls_mem->last_flag);
    }

  }

  /* Call LS setup routine -- the LS may call arkLsPSetup, who will
     pass the heuristic suggestions above to the user code(s) */
  arkls_mem->last_flag = SUNLinSolSetup(arkls_mem->LS, arkls_mem->A);

  /* If the SUNMatrix was NULL, update heuristics flags */
  if (arkls_mem->A == NULL) {

    /* If user set jcur to SUNTRUE, increment npe and save nst value */
    if (*jcurPtr) {
      arkls_mem->npe++;
      arkls_mem->nstlj = ark_mem->nst;
    }

    /* Update jcurPtr flag if we suggested an update */
    if (arkls_mem->jbad) *jcurPtr = SUNTRUE;
  }

  return(arkls_mem->last_flag);
}

/*---------------------------------------------------------------
  arkLsSolve: interfaces between ARKode and the generic
  SUNLinearSolver object LS, by setting the appropriate tolerance
  and scaling vectors, calling the solver, and accumulating
  statistics from the solve for use/reporting by ARKode.

  When using a non-NULL SUNMatrix, this will additionally scale
  the solution appropriately when gamrat != 1.
  ---------------------------------------------------------------*/
int arkLsSolve(void* arkode_mem, N_Vector b, realtype tnow,
               N_Vector ynow, N_Vector fnow, realtype eRNrm, int mnewt)
{
  realtype    bnorm, resnorm;
  ARKodeMem   ark_mem;
  ARKLsMem    arkls_mem;
  realtype    gamma, gamrat, delta, deltar, ewt_mean;
  booleantype dgamma_fail, *jcur;
  int         nli_inc, nps_inc, retval, LSType;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "arkLsSolve",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Set scalar tcur and vectors ycur and fcur for use by the
     Atimes and Psolve interface routines */
  arkls_mem->tcur = tnow;
  arkls_mem->ycur = ynow;
  arkls_mem->fcur = fnow;

  /* Retrieve the LS type */
  LSType = SUNLinSolGetType(arkls_mem->LS);

  /* If the linear solver is iterative:
     test norm(b), if small, return x = 0 or x = b;
     set linear solver tolerance (in left/right scaled 2-norm) */
  if ( (LSType == SUNLINEARSOLVER_ITERATIVE) ||
       (LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) ) {
    deltar = arkls_mem->eplifac * eRNrm;
    bnorm = N_VWrmsNorm(b, ark_mem->rwt);
    if (bnorm <= deltar) {
      if (mnewt > 0) N_VConst(ZERO, b);
      arkls_mem->last_flag = ARKLS_SUCCESS;
      return(arkls_mem->last_flag);
    }
    delta = deltar * arkls_mem->sqrtN;
  } else {
    delta = ZERO;
  }

  /* Set initial guess x = 0 to LS */
  N_VConst(ZERO, arkls_mem->x);

  /* Set scaling vectors for LS to use (if applicable) */
  if (arkls_mem->LS->ops->setscalingvectors) {
    retval = SUNLinSolSetScalingVectors(arkls_mem->LS,
                                        ark_mem->ewt,
                                        ark_mem->rwt);
    if (retval != SUNLS_SUCCESS) {
      arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKLS", "arkLsSolve",
                      "Error in call to SUNLinSolSetScalingVectors");
      arkls_mem->last_flag = ARKLS_SUNLS_FAIL;
      return(arkls_mem->last_flag);
    }

  /* If solver is iterative and does not support scaling vectors, update the
     tolerance in an attempt to account for ewt/rwt vectors.  We make the
     following assumptions:
       1. rwt = ewt (i.e. the units of solution and residual are the same)
       2. ewt_i = ewt_mean, for i=0,...,n-1 (i.e. the solution units are identical)
       3. the linear solver uses a basic 2-norm to measure convergence
     Hence (using the notation from sunlinsol_spgmr.h, with S = diag(ewt)),
           || bbar - Abar xbar ||_2 < tol
       <=> || S b - S A x ||_2 < tol
       <=> || S (b - A x) ||_2 < tol
       <=> \sum_{i=0}^{n-1} (ewt_i (b - A x)_i)^2 < tol^2
       <=> ewt_mean^2 \sum_{i=0}^{n-1} (b - A x_i)^2 < tol^2
       <=> \sum_{i=0}^{n-1} (b - A x_i)^2 < tol^2 / ewt_mean^2
       <=> || b - A x ||_2 < tol / ewt_mean
     So we compute ewt_mean = ||ewt||_RMS = ||ewt||_2 / sqrt(n), and scale
     the desired tolerance accordingly. */
  } else if ( (LSType == SUNLINEARSOLVER_ITERATIVE) ||
              (LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) ) {

    ewt_mean = SUNRsqrt( N_VDotProd(ark_mem->ewt, ark_mem->ewt) ) / arkls_mem->sqrtN;
    delta /= ewt_mean;

  }

  /* Store previous nps value in nps_inc */
  nps_inc = arkls_mem->nps;

  /* If a user-provided jtsetup routine is supplied, call that here */
  if (arkls_mem->jtsetup) {
    arkls_mem->last_flag = arkls_mem->jtsetup(tnow, ynow, fnow,
                                              arkls_mem->Jt_data);
    arkls_mem->njtsetup++;
    if (arkls_mem->last_flag) {
      arkProcessError(ark_mem, retval, "ARKLS",
                      "arkLsSolve", MSG_LS_JTSETUP_FAILED);
      return(arkls_mem->last_flag);
    }
  }

  /* Call solver, and copy x to b */
  retval = SUNLinSolSolve(arkls_mem->LS, arkls_mem->A,
                          arkls_mem->x, b, delta);
  N_VScale(ONE, arkls_mem->x, b);

  /* If using a direct or matrix-iterative solver, scale the correction to
     account for change in gamma (this is only beneficial if M==I) */
  if ( (LSType == SUNLINEARSOLVER_DIRECT) ||
       (LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) ) {
    arkls_mem->last_flag = ark_mem->step_getgammas(arkode_mem, &gamma, &gamrat,
                                                   &jcur, &dgamma_fail);
    if (arkls_mem->last_flag != ARK_SUCCESS) {
      arkProcessError(ark_mem, arkls_mem->last_flag, "ARKLS", "arkLsSolve",
                      "An error occurred in ark_step_getgammas");
      return(arkls_mem->last_flag);
    }
    if (gamrat != ONE)  N_VScale(TWO/(ONE + gamrat), b, b);
  }

  /* Retrieve statistics from iterative linear solvers */
  resnorm = ZERO;
  nli_inc = 0;
  if ( (LSType == SUNLINEARSOLVER_ITERATIVE) ||
       (LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) ) {
    if (arkls_mem->LS->ops->resnorm)
      resnorm = SUNLinSolResNorm(arkls_mem->LS);
    if (arkls_mem->LS->ops->numiters)
      nli_inc = SUNLinSolNumIters(arkls_mem->LS);
  }

  /* Increment counters nli and ncfl */
  arkls_mem->nli += nli_inc;
  if (retval != SUNLS_SUCCESS) arkls_mem->ncfl++;

  /* Log solver statistics to diagnostics file (if requested) */
  if (ark_mem->report)
    fprintf(ark_mem->diagfp, "ARKLS  kry  %"RSYM"  %"RSYM"  %i  %i\n",
            bnorm, resnorm, nli_inc, (int) arkls_mem->nps - nps_inc);

  /* Interpret solver return value  */
  arkls_mem->last_flag = retval;

  switch(retval) {

  case SUNLS_SUCCESS:
    return(0);
    break;
  case SUNLS_RES_REDUCED:
    /* allow reduction but not solution on first nonlinear iteration,
       otherwise return with a recoverable failure */
    if (mnewt == 0) return(0);
    else            return(1);
    break;
  case SUNLS_CONV_FAIL:
  case SUNLS_ATIMES_FAIL_REC:
  case SUNLS_PSOLVE_FAIL_REC:
  case SUNLS_PACKAGE_FAIL_REC:
  case SUNLS_QRFACT_FAIL:
  case SUNLS_LUFACT_FAIL:
    return(1);
    break;
  case SUNLS_MEM_NULL:
  case SUNLS_ILL_INPUT:
  case SUNLS_MEM_FAIL:
  case SUNLS_GS_FAIL:
  case SUNLS_QRSOL_FAIL:
    return(-1);
    break;
  case SUNLS_PACKAGE_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_PACKAGE_FAIL_UNREC, "ARKLS",
                    "arkLsSolve",
                    "Failure in SUNLinSol external package");
    return(-1);
    break;
  case SUNLS_ATIMES_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_ATIMES_FAIL_UNREC, "ARKLS",
                    "arkLsSolve", MSG_LS_JTIMES_FAILED);
    return(-1);
    break;
  case SUNLS_PSOLVE_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_PSOLVE_FAIL_UNREC, "ARKLS",
                    "arkLsSolve", MSG_LS_PSOLVE_FAILED);
    return(-1);
    break;
  }

  return(0);
}


/*---------------------------------------------------------------
  arkLsFree frees memory associates with the ARKLs system
  solver interface.
  ---------------------------------------------------------------*/
int arkLsFree(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKLsMem  arkls_mem;
  void*     ark_step_lmem;

  /* Return immediately if ARKodeMem, ARKLsMem are NULL */
  if (arkode_mem == NULL)  return (ARKLS_SUCCESS);
  ark_mem = (ARKodeMem) arkode_mem;
  ark_step_lmem = ark_mem->step_getlinmem(arkode_mem);
  if (ark_step_lmem == NULL)  return(ARKLS_SUCCESS);
  arkls_mem = (ARKLsMem) ark_step_lmem;

  /* Free N_Vector memory */
  if (arkls_mem->ytemp) {
    N_VDestroy(arkls_mem->ytemp);
    arkls_mem->ytemp = NULL;
  }
  if (arkls_mem->x) {
    N_VDestroy(arkls_mem->x);
    arkls_mem->x = NULL;
  }

  /* Free savedJ memory */
  if (arkls_mem->savedJ) {
    SUNMatDestroy(arkls_mem->savedJ);
    arkls_mem->savedJ = NULL;
  }

  /* Nullify other N_Vector pointers */
  arkls_mem->ycur = NULL;
  arkls_mem->fcur = NULL;

  /* Nullify other SUNMatrix pointer */
  arkls_mem->A = NULL;

  /* Free preconditioner memory (if applicable) */
  if (arkls_mem->pfree)  arkls_mem->pfree(ark_mem);

  /* free ARKLs interface structure */
  free(arkls_mem);

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLsMassInitialize performs remaining initializations specific
  to the iterative linear solver interface (and solver itself)
  ---------------------------------------------------------------*/
int arkLsMassInitialize(void *arkode_mem)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          retval;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLsMassInitialize",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* reset counters */
  arkLsInitializeMassCounters(arkls_mem);

  /* perform checks for mass matrix constructor or mass matrix-vector product routine exist */
  if (arkls_mem->M == NULL) {
    if (arkls_mem->mtimes == NULL) {
      arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS", "arkLsMassInitialize",
                      "Missing user-provided mass matrix-vector product routine");
      arkls_mem->last_flag = ARKLS_ILL_INPUT;
      return(arkls_mem->last_flag);
    }
  } else {
    if (arkls_mem->mass == NULL) {
      arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                      "arkLsMassInitialize",
                      "Missing user-provided mass-matrix routine");
      arkls_mem->last_flag = ARKLS_ILL_INPUT;
      return(arkls_mem->last_flag);
    }
  }

  /* ensure that a mass matrix solver exists */
  if (arkls_mem->LS == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKLS",
                    "arkLsMassInitialize",
                    "Missing SUNLinearSolver object");
    arkls_mem->last_flag = ARKLS_ILL_INPUT;
    return(arkls_mem->last_flag);
  }

  /* if M is NULL and neither pset or mtsetup are present, then
     arkLsMassSetup does not need to be called, so set the
     msetup function to NULL */
  if ( (arkls_mem->M == NULL) &&
       (arkls_mem->pset == NULL) &&
       (arkls_mem->mtsetup == NULL) &&
       (ark_mem->step_disablemsetup != NULL) )
    ark_mem->step_disablemsetup(arkode_mem);

  /* Call LS initialize routine */
  arkls_mem->last_flag = SUNLinSolInitialize(arkls_mem->LS);
  return(arkls_mem->last_flag);
}


/*---------------------------------------------------------------
  arkLsMassSetup calls the LS 'setup' routine.
  ---------------------------------------------------------------*/
int arkLsMassSetup(void *arkode_mem, N_Vector vtemp1,
                   N_Vector vtemp2, N_Vector vtemp3)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  booleantype  call_mtsetup, call_lssetup;
  int          retval;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLsMassSetup",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Determine whether to call user-provided mtsetup routine */
  call_mtsetup = SUNFALSE;
  if ( (arkls_mem->mtsetup) &&
       (arkls_mem->time_dependent || (!arkls_mem->nmtsetup)) )
    call_mtsetup = SUNTRUE;

  /* call user-provided mtsetup routine if applicable */
  if (call_mtsetup) {
    arkls_mem->last_flag = arkls_mem->mtsetup(ark_mem->tcur,
                                              arkls_mem->mt_data);
    arkls_mem->nmtsetup++;
    if (arkls_mem->last_flag != 0) {
      arkProcessError(ark_mem, arkls_mem->last_flag, "ARKLS",
                      "arkLsMassSetup", MSG_LS_MTSETUP_FAILED);
      return(arkls_mem->last_flag);
    }
  }


  /* Perform user-facing setup based on whether this is matrix-free */
  if (arkls_mem->M == NULL) {

    /*** matrix-free -- only call LS setup if preconditioner setup exists ***/
    call_lssetup = (arkls_mem->pset != NULL);

  } else {

    /*** matrix-based ***/

    /* If mass matrix is not time dependent, and if it has been set up
       previously, just reuse existing M and M_lu */
    if (!arkls_mem->time_dependent && arkls_mem->nmsetups) {
      arkls_mem->last_flag = ARKLS_SUCCESS;
      return(arkls_mem->last_flag);
    }

    /* Update mass matrix */
    retval = SUNMatZero(arkls_mem->M);
    if (retval) {
      arkProcessError(ark_mem, ARKLS_SUNMAT_FAIL, "ARKLS",
                      "arkLsMassSetup",  MSG_LS_SUNMAT_FAILED);
      arkls_mem->last_flag = ARKLS_SUNMAT_FAIL;
      return(arkls_mem->last_flag);
    }

    retval = arkls_mem->mass(ark_mem->tcur, arkls_mem->M,
                             ark_mem->user_data,
                             vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      arkProcessError(ark_mem, ARKLS_MASSFUNC_UNRECVR, "ARKLS",
                      "arkLsMassSetup",  MSG_LS_MASSFUNC_FAILED);
      arkls_mem->last_flag = ARKLS_MASSFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arkls_mem->last_flag = ARKLS_MASSFUNC_RECVR;
      return(1);
    }

    /* Copy M into M_lu for factorization */
    retval = SUNMatCopy(arkls_mem->M, arkls_mem->M_lu);
    if (retval) {
      arkProcessError(ark_mem, ARKLS_SUNMAT_FAIL, "ARKLS",
                      "arkLsMassSetup",  MSG_LS_SUNMAT_FAILED);
      arkls_mem->last_flag = ARKLS_SUNMAT_FAIL;
      return(arkls_mem->last_flag);
    }

    /* signal call to LS setup routine */
    call_lssetup = SUNTRUE;

  }

  /* Call LS setup routine if applicable, and return */
  if (call_lssetup) {
    arkls_mem->last_flag = SUNLinSolSetup(arkls_mem->LS,
                                          arkls_mem->M_lu);
    arkls_mem->nmsetups++;
  }
  return(arkls_mem->last_flag);
}


/*---------------------------------------------------------------
  arkLsMassSolve: interfaces between ARKode and the generic
  SUNLinearSolver object LS, by setting the appropriate tolerance
  and scaling vectors, calling the solver, and accumulating
  statistics from the solve for use/reporting by ARKode.
  ---------------------------------------------------------------*/
int arkLsMassSolve(void *arkode_mem, N_Vector b, realtype nlscoef)
{
  realtype     resnorm, delta, rwt_mean;
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  int          nli_inc, nps_inc, retval, LSType;

  /* access ARKLsMassMem structure */
  retval = arkLs_AccessMassMem(arkode_mem, "arkLsMassSolve",
                               &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Retrieve the LS type */
  LSType = SUNLinSolGetType(arkls_mem->LS);

  /* Set input tolerance for iterative solvers */
  if ( (LSType == SUNLINEARSOLVER_ITERATIVE) ||
       (LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) ) {
    delta = arkls_mem->eplifac * nlscoef * arkls_mem->sqrtN;
  } else {
    delta = ZERO;
  }

  /* Set initial guess x = 0 for LS */
  N_VConst(ZERO, arkls_mem->x);

  /* Set scaling vectors for LS to use (if applicable) */
  if (arkls_mem->LS->ops->setscalingvectors) {
    retval = SUNLinSolSetScalingVectors(arkls_mem->LS,
                                        ark_mem->ewt,
                                        ark_mem->rwt);
    if (retval != SUNLS_SUCCESS) {
      arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKLS", "arkLsMassSolve",
                      "Error in call to SUNLinSolSetScalingVectors");
      arkls_mem->last_flag = ARKLS_SUNLS_FAIL;
      return(arkls_mem->last_flag);
    }

  /* If solver is iterative and does not support scaling vectors, update the
     tolerance in an attempt to account for rwt vector.  We make the
     following assumptions:
       1. rwt_i = rwt_mean, for i=0,...,n-1 (i.e. the solution units are identical)
       2. the linear solver uses a basic 2-norm to measure convergence
     Hence (using the notation from sunlinsol_spgmr.h, with S = diag(rwt)),
           || bbar - Abar xbar ||_2 < tol
       <=> || S b - S A x ||_2 < tol
       <=> || S (b - A x) ||_2 < tol
       <=> \sum_{i=0}^{n-1} (rwt_i (b - A x)_i)^2 < tol^2
       <=> rwt_mean^2 \sum_{i=0}^{n-1} (b - A x_i)^2 < tol^2
       <=> \sum_{i=0}^{n-1} (b - A x_i)^2 < tol^2 / rwt_mean^2
       <=> || b - A x ||_2 < tol / rwt_mean
     So we compute rwt_mean = ||rwt||_RMS = ||rwt||_2 / sqrt(n), and scale
     the desired tolerance accordingly. */
  } else if ( (LSType == SUNLINEARSOLVER_ITERATIVE) ||
              (LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) ) {

    rwt_mean = SUNRsqrt( N_VDotProd(ark_mem->rwt, ark_mem->rwt) ) / arkls_mem->sqrtN;
    delta /= rwt_mean;

  }

  /* Store previous nps value in nps_inc */
  nps_inc = arkls_mem->nps;

  /* Call solver, copy x to b, and increment mass solver counter */
  retval = SUNLinSolSolve(arkls_mem->LS, arkls_mem->M_lu,
                          arkls_mem->x, b, delta);
  N_VScale(ONE, arkls_mem->x, b);
  arkls_mem->nmsolves++;

  /* Retrieve statistics from iterative linear solvers */
  resnorm = ZERO;
  nli_inc = 0;
  if ( (LSType == SUNLINEARSOLVER_ITERATIVE) ||
       (LSType == SUNLINEARSOLVER_MATRIX_ITERATIVE) ) {
    if (arkls_mem->LS->ops->resnorm)
      resnorm = SUNLinSolResNorm(arkls_mem->LS);
    if (arkls_mem->LS->ops->numiters)
      nli_inc = SUNLinSolNumIters(arkls_mem->LS);
  }

  /* Increment counters nli and ncfl */
  arkls_mem->nli += nli_inc;
  if (retval != SUNLS_SUCCESS) arkls_mem->ncfl++;

  /* Log solver statistics to diagnostics file (if requested) */
  if (ark_mem->report)
    fprintf(ark_mem->diagfp, "ARKLS  mass  %"RSYM"  %i  %i\n",
            resnorm, nli_inc, (int) arkls_mem->nps - nps_inc);

  /* Interpret solver return value  */
  arkls_mem->last_flag = retval;

  switch(retval) {

  case SUNLS_SUCCESS:
    return(0);
    break;
  case SUNLS_RES_REDUCED:
  case SUNLS_CONV_FAIL:
  case SUNLS_ATIMES_FAIL_REC:
  case SUNLS_PSOLVE_FAIL_REC:
  case SUNLS_PACKAGE_FAIL_REC:
  case SUNLS_QRFACT_FAIL:
  case SUNLS_LUFACT_FAIL:
    return(1);
    break;
  case SUNLS_MEM_NULL:
  case SUNLS_ILL_INPUT:
  case SUNLS_MEM_FAIL:
  case SUNLS_GS_FAIL:
  case SUNLS_QRSOL_FAIL:
    return(-1);
    break;
  case SUNLS_PACKAGE_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_PACKAGE_FAIL_UNREC, "ARKLS",
                    "arkLsMassSolve",
                    "Failure in SUNLinSol external package");
    return(-1);
    break;
  case SUNLS_ATIMES_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_ATIMES_FAIL_UNREC, "ARKLS",
                    "arkLsMassSolve", MSG_LS_MTIMES_FAILED);
    return(-1);
    break;
  case SUNLS_PSOLVE_FAIL_UNREC:
    arkProcessError(ark_mem, SUNLS_PSOLVE_FAIL_UNREC, "ARKLS",
                    "arkLsMassSolve", MSG_LS_PSOLVE_FAILED);
    return(-1);
    break;
  }

  return(0);
}


/*---------------------------------------------------------------
  arkLsMassFree frees memory associates with the ARKLs mass
  matrix solver interface.
  ---------------------------------------------------------------*/
int arkLsMassFree(void *arkode_mem)
{
  ARKodeMem    ark_mem;
  ARKLsMassMem arkls_mem;
  void*        ark_step_massmem;

  /* Return immediately if ARKodeMem, ARKLsMassMem are NULL */
  if (arkode_mem == NULL)  return (ARKLS_SUCCESS);
  ark_mem = (ARKodeMem) arkode_mem;
  ark_step_massmem = ark_mem->step_getmassmem(arkode_mem);
  if (ark_step_massmem == NULL)  return(ARKLS_SUCCESS);
  arkls_mem = (ARKLsMassMem) ark_step_massmem;

  /* detach ARKLs interface routines from LS object (ignore return values) */
  if (arkls_mem->LS->ops->setatimes)
    SUNLinSolSetATimes(arkls_mem->LS, NULL, NULL);

  if (arkls_mem->LS->ops->setpreconditioner)
    SUNLinSolSetPreconditioner(arkls_mem->LS, NULL, NULL, NULL);

  /* Free N_Vector memory */
  if (arkls_mem->x) {
    N_VDestroy(arkls_mem->x);
    arkls_mem->x = NULL;
  }

  /* Free M_lu memory */
  if (arkls_mem->M_lu) {
    SUNMatDestroy(arkls_mem->M_lu);
    arkls_mem->M_lu = NULL;
  }

  /* Nullify other N_Vector pointers */
  arkls_mem->ycur = NULL;

  /* Nullify other SUNMatrix pointer */
  arkls_mem->M = NULL;

  /* Free preconditioner memory (if applicable) */
  if (arkls_mem->pfree)
    arkls_mem->pfree(ark_mem);

  /* free ARKLs interface structure */
  free(arkls_mem);

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  arkLsInitializeCounters and arkLsInitializeMassCounters:

  These routines reset all counters from an ARKLsMem or
  ARKLsMassMem structure.
  ---------------------------------------------------------------*/
int arkLsInitializeCounters(ARKLsMem arkls_mem)
{
  arkls_mem->nje      = 0;
  arkls_mem->nfeDQ    = 0;
  arkls_mem->nstlj    = 0;
  arkls_mem->npe      = 0;
  arkls_mem->nli      = 0;
  arkls_mem->nps      = 0;
  arkls_mem->ncfl     = 0;
  arkls_mem->njtsetup = 0;
  arkls_mem->njtimes  = 0;
  return(0);
}

int arkLsInitializeMassCounters(ARKLsMassMem arkls_mem)
{
  arkls_mem->nmsetups = 0;
  arkls_mem->nmsolves = 0;
  arkls_mem->nmtsetup = 0;
  arkls_mem->nmtimes  = 0;
  arkls_mem->npe      = 0;
  arkls_mem->nli      = 0;
  arkls_mem->nps      = 0;
  arkls_mem->ncfl     = 0;
  return(0);
}


/*---------------------------------------------------------------
  arkLs_AccessLMem and arkLs_AccessMassMem:

  Shortcut routines to unpack ark_mem, ls_mem and mass_mem
  structures from void* pointer.  If any is missing it returns
  ARKLS_MEM_NULL, ARKLS_LMEM_NULL or ARKLS_MASSMEM_NULL.
  ---------------------------------------------------------------*/
int arkLs_AccessLMem(void* arkode_mem, const char *fname,
                     ARKodeMem *ark_mem, ARKLsMem *arkls_mem)
{
  void* ark_step_lmem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARKLS_MEM_NULL, "ARKLS",
                    fname, MSG_LS_ARKMEM_NULL);
    return(ARKLS_MEM_NULL);
  }
  *ark_mem = (ARKodeMem) arkode_mem;
  ark_step_lmem = (*ark_mem)->step_getlinmem(arkode_mem);
  if (ark_step_lmem==NULL) {
    arkProcessError(*ark_mem, ARKLS_LMEM_NULL, "ARKLS",
                    fname, MSG_LS_LMEM_NULL);
    return(ARKLS_LMEM_NULL);
  }
  *arkls_mem = (ARKLsMem) ark_step_lmem;
  return(ARKLS_SUCCESS);
}

int arkLs_AccessMassMem(void* arkode_mem, const char *fname,
                        ARKodeMem *ark_mem, ARKLsMassMem *arkls_mem)
{
  void* ark_step_massmem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARKLS_MEM_NULL, "ARKLS",
                    fname, MSG_LS_ARKMEM_NULL);
    return(ARKLS_MEM_NULL);
  }
  *ark_mem = (ARKodeMem) arkode_mem;
  ark_step_massmem = (*ark_mem)->step_getmassmem(arkode_mem);
  if (ark_step_massmem==NULL) {
    arkProcessError(*ark_mem, ARKLS_MASSMEM_NULL, "ARKLS",
                    fname, MSG_LS_MASSMEM_NULL);
    return(ARKLS_MASSMEM_NULL);
  }
  *arkls_mem = (ARKLsMassMem) ark_step_massmem;
  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
