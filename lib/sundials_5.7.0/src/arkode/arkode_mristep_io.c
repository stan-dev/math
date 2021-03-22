/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
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
 * This is the implementation file for the optional input and output functions
 * for the ARKode MRIStep time stepper module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_mristep_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  DEPRECATED MRIStep optional I/O routines -- these are only
  here for backwards compatibility.
  ===============================================================*/

int MRIStepWriteButcher(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepWriteCoupling",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* return error flag/message to user */
  arkProcessError(ark_mem, ARK_WARNING, "ARKode::MRIStep",
                  "MRIStepWriteButcher",
                  "This routine is deprecated, and will be removed in a future release");
  return(ARK_WARNING);
}

int MRIStepGetCurrentButcherTables(void *arkode_mem, ARKodeButcherTable *B)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepWriteCoupling",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* return error flag/message to user */
  arkProcessError(ark_mem, ARK_WARNING, "ARKode::MRIStep",
                  "MRIStepGetCurrentButcherTables",
                  "This routine is deprecated, and will be removed in a future release");
  *B = NULL;
  return(ARK_WARNING);
}



/*===============================================================
  MRIStep Optional input functions (wrappers for generic ARKode
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int MRIStepSetDenseOrder(void *arkode_mem, int dord) {
  return(MRIStepSetInterpolantDegree(arkode_mem, dord)); }
int MRIStepSetInterpolantDegree(void *arkode_mem, int degree) {
  if (degree < 0) degree = ARK_INTERP_MAX_DEGREE;
  return(arkSetInterpolantDegree(arkode_mem, degree)); }
int MRIStepSetInterpolantType(void *arkode_mem, int itype) {
  return(arkSetInterpolantType(arkode_mem, itype)); }
int MRIStepSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                           void *eh_data) {
  return(arkSetErrHandlerFn(arkode_mem, ehfun, eh_data)); }
int MRIStepSetErrFile(void *arkode_mem, FILE *errfp) {
  return(arkSetErrFile(arkode_mem, errfp)); }
int MRIStepSetDiagnostics(void *arkode_mem, FILE *diagfp) {
  return(arkSetDiagnostics(arkode_mem, diagfp)); }
int MRIStepSetMaxNumSteps(void *arkode_mem, long int mxsteps) {
  return(arkSetMaxNumSteps(arkode_mem, mxsteps)); }
int MRIStepSetMaxHnilWarns(void *arkode_mem, int mxhnil) {
  return(arkSetMaxHnilWarns(arkode_mem, mxhnil)); }
int MRIStepSetStopTime(void *arkode_mem, realtype tstop) {
  return(arkSetStopTime(arkode_mem, tstop)); }
int MRIStepSetRootDirection(void *arkode_mem, int *rootdir) {
  return(arkSetRootDirection(arkode_mem, rootdir)); }
int MRIStepSetNoInactiveRootWarn(void *arkode_mem) {
  return(arkSetNoInactiveRootWarn(arkode_mem)); }
int MRIStepSetPostprocessStepFn(void *arkode_mem,
                                ARKPostProcessFn ProcessStep) {
  return(arkSetPostprocessStepFn(arkode_mem, ProcessStep)); }
int MRIStepSetPostprocessStageFn(void *arkode_mem,
                                 ARKPostProcessFn ProcessStage) {
  return(arkSetPostprocessStageFn(arkode_mem, ProcessStage)); }


/*---------------------------------------------------------------
  These wrappers for ARKLs module 'set' routines all are
  documented in arkode_mristep.h.
  ---------------------------------------------------------------*/
int MRIStepSetLinearSolver(void *arkode_mem, SUNLinearSolver LS,
                           SUNMatrix A) {
  return(arkLSSetLinearSolver(arkode_mem, LS, A)); }
int MRIStepSetJacFn(void *arkode_mem, ARKLsJacFn jac) {
  return(arkLSSetJacFn(arkode_mem, jac)); }
int MRIStepSetJacEvalFrequency(void *arkode_mem, long int msbj) {
  return(arkLSSetJacEvalFrequency(arkode_mem, msbj)); }
int MRIStepSetLinearSolutionScaling(void *arkode_mem, booleantype onoff) {
  return(arkLSSetLinearSolutionScaling(arkode_mem, onoff)); }
int MRIStepSetEpsLin(void *arkode_mem, realtype eplifac) {
  return(arkLSSetEpsLin(arkode_mem, eplifac)); }
int MRIStepSetLSNormFactor(void *arkode_mem, realtype nrmfac) {
  return(arkLSSetNormFactor(arkode_mem, nrmfac)); }
int MRIStepSetPreconditioner(void *arkode_mem, ARKLsPrecSetupFn psetup,
                             ARKLsPrecSolveFn psolve) {
  return(arkLSSetPreconditioner(arkode_mem, psetup, psolve)); }
int MRIStepSetJacTimes(void *arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                       ARKLsJacTimesVecFn jtimes) {
  return(arkLSSetJacTimes(arkode_mem, jtsetup, jtimes)); }
int MRIStepSetJacTimesRhsFn(void *arkode_mem, ARKRhsFn jtimesRhsFn) {
  return(arkLSSetJacTimesRhsFn(arkode_mem, jtimesRhsFn)); }
int MRIStepSetLinSysFn(void *arkode_mem, ARKLsLinSysFn linsys) {
  return(arkLSSetLinSysFn(arkode_mem, linsys)); }


/*===============================================================
  MRIStep Optional output functions (wrappers for generic ARKode
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int MRIStepGetNumSteps(void *arkode_mem, long int *nssteps) {
  return(arkGetNumSteps(arkode_mem, nssteps)); }
int MRIStepGetLastStep(void *arkode_mem, realtype *hlast) {
  return(arkGetLastStep(arkode_mem, hlast)); }
int MRIStepGetCurrentTime(void *arkode_mem, realtype *tcur) {
  return(arkGetCurrentTime(arkode_mem, tcur)); }
int MRIStepGetCurrentState(void *arkode_mem, N_Vector *state) {
  return(arkGetCurrentState(arkode_mem, state)); }
int MRIStepGetTolScaleFactor(void *arkode_mem, realtype *tolsfact) {
  return(arkGetTolScaleFactor(arkode_mem, tolsfact)); }
int MRIStepGetErrWeights(void *arkode_mem, N_Vector eweight) {
  return(arkGetErrWeights(arkode_mem, eweight)); }
int MRIStepGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw) {
  return(arkGetWorkSpace(arkode_mem, lenrw, leniw)); }
int MRIStepGetNumGEvals(void *arkode_mem, long int *ngevals) {
  return(arkGetNumGEvals(arkode_mem, ngevals)); }
int MRIStepGetRootInfo(void *arkode_mem, int *rootsfound) {
  return(arkGetRootInfo(arkode_mem, rootsfound)); }
char *MRIStepGetReturnFlagName(long int flag) {
  return(arkGetReturnFlagName(flag)); }

/*---------------------------------------------------------------
  These wrappers for ARKLs module 'get' routines all are
  documented in arkode_mristep.h.
  ---------------------------------------------------------------*/
int MRIStepGetLinWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS) {
  return(arkLSGetWorkSpace(arkode_mem, lenrwLS, leniwLS)); }
int MRIStepGetNumJacEvals(void *arkode_mem, long int *njevals) {
  return(arkLSGetNumJacEvals(arkode_mem, njevals)); }
int MRIStepGetNumPrecEvals(void *arkode_mem, long int *npevals) {
  return(arkLSGetNumPrecEvals(arkode_mem, npevals)); }
int MRIStepGetNumPrecSolves(void *arkode_mem, long int *npsolves) {
  return(arkLSGetNumPrecSolves(arkode_mem, npsolves)); }
int MRIStepGetNumLinIters(void *arkode_mem, long int *nliters) {
  return(arkLSGetNumLinIters(arkode_mem, nliters)); }
int MRIStepGetNumLinConvFails(void *arkode_mem, long int *nlcfails) {
  return(arkLSGetNumConvFails(arkode_mem, nlcfails)); }
int MRIStepGetNumJTSetupEvals(void *arkode_mem, long int *njtsetups) {
  return(arkLSGetNumJTSetupEvals(arkode_mem, njtsetups)); }
int MRIStepGetNumJtimesEvals(void *arkode_mem, long int *njvevals) {
  return(arkLSGetNumJtimesEvals(arkode_mem, njvevals)); }
int MRIStepGetNumLinRhsEvals(void *arkode_mem, long int *nfevalsLS) {
  return(arkLSGetNumRhsEvals(arkode_mem, nfevalsLS)); }
int MRIStepGetLastLinFlag(void *arkode_mem, long int *flag) {
  return(arkLSGetLastFlag(arkode_mem, flag)); }
char *MRIStepGetLinReturnFlagName(long int flag) {
  return(arkLSGetReturnFlagName(flag)); }



/*===============================================================
  MRIStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepSetUserData:

  Wrapper for generic arkSetUserData and arkLSSetUserData
  routines.
  ---------------------------------------------------------------*/
int MRIStepSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem        ark_mem;
  ARKodeMRIStepMem step_mem;
  int              retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetUserData",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* set user_data in ARKode mem */
  retval = arkSetUserData(arkode_mem, user_data);
  if (retval != ARK_SUCCESS) return(retval);

  /* set user data in ARKodeLS mem */
  if (step_mem->lmem != NULL) {
    retval = arkLSSetUserData(arkode_mem, user_data);
    if (retval != ARKLS_SUCCESS) return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetDefaults:

  Resets all MRIStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
int MRIStepSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set default values for integrator optional inputs */
  step_mem->q              = 3;              /* method order */
  step_mem->p              = 0;              /* embedding order */
  step_mem->predictor      = 0;              /* trivial predictor */
  step_mem->linear         = SUNFALSE;       /* nonlinear problem */
  step_mem->linear_timedep = SUNTRUE;        /* dfs/dy depends on t */
  step_mem->implicit       = SUNFALSE;       /* fs(t,y) is explicit */
  step_mem->maxcor         = MAXCOR;         /* max nonlinear iters/stage */
  step_mem->nlscoef        = NLSCOEF;        /* nonlinear tolerance coefficient */
  step_mem->crdown         = CRDOWN;         /* nonlinear convergence estimate coeff. */
  step_mem->rdiv           = RDIV;           /* nonlinear divergence tolerance */
  step_mem->dgmax          = DGMAX;          /* max gamma change before recomputing J or P */
  step_mem->msbp           = MSBP;           /* max steps between updates to J or P */
  step_mem->stages         = 0;              /* no stages */
  step_mem->istage         = 0;              /* current stage index */
  step_mem->MRIC           = NULL;           /* no slow->fast coupling */
  step_mem->NLS            = NULL;           /* no nonlinear solver object */
  step_mem->jcur           = SUNFALSE;
  step_mem->convfail       = ARK_NO_FAILURES;
  step_mem->stage_predict  = NULL;           /* no user-supplied stage predictor */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetLinear:

  Specifies that the implicit slow function, fs(t,y), is linear
  in y, and to tighten the linear solver tolerances while taking
  only one Newton iteration.  DO NOT USE IN COMBINATION WITH THE
  FIXED-POINT SOLVER.  Automatically tightens DeltaGammaMax
  to ensure that step size changes cause Jacobian recomputation.

  The argument should be 1 or 0, where 1 indicates that the
  Jacobian of fs with respect to y depends on time, and
  0 indicates that it is not time dependent.  Alternately, when
  using an iterative linear solver this flag denotes time
  dependence of the preconditioner.
  ---------------------------------------------------------------*/
int MRIStepSetLinear(void *arkode_mem, int timedepend)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetLinear",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set parameters */
  step_mem->linear = SUNTRUE;
  step_mem->linear_timedep = (timedepend == 1);
  step_mem->dgmax = RCONST(100.0)*UNIT_ROUNDOFF;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetNonlinear:

  Specifies that the implicit slow function, fs(t,y), is
  nonlinear in y.  Used to undo a previous call to
  MRIStepSetLinear.  Automatically loosens DeltaGammaMax back to
  default value.
  ---------------------------------------------------------------*/
int MRIStepSetNonlinear(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetNonlinear",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set parameters */
  step_mem->linear = SUNFALSE;
  step_mem->linear_timedep = SUNTRUE;
  step_mem->dgmax = DGMAX;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetCoupling:

  Specifies to use a customized coupling structure for the slow
  portion of the system.
  ---------------------------------------------------------------*/
int MRIStepSetCoupling(void *arkode_mem, MRIStepCoupling MRIC)
{
  int retval, is, stagetype;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  sunindextype Tlrw, Tliw;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetCoupling",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check for illegal inputs */
  if (MRIC == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetCoupling", MSG_ARK_NO_MEM);
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and coupling structure */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  MRIStepCoupling_Space(step_mem->MRIC, &Tliw, &Tlrw);
  MRIStepCoupling_Free(step_mem->MRIC);
  step_mem->MRIC = NULL;
  ark_mem->liw -= Tliw;
  ark_mem->lrw -= Tlrw;

  /* set the relevant parameters */
  step_mem->stages = MRIC->stages;
  step_mem->q = MRIC->q;
  step_mem->p = MRIC->p;

  /* copy the coupling structure in step memory */
  step_mem->MRIC = MRIStepCoupling_Copy(MRIC);
  if (step_mem->MRIC == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetCoupling", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  MRIStepCoupling_Space(step_mem->MRIC, &Tliw, &Tlrw);
  ark_mem->liw += Tliw;
  ark_mem->lrw += Tlrw;

  /* determine whether fs requires an implicit solver */
  step_mem->implicit = SUNFALSE;
  for (is=0; is<step_mem->stages; is++) {
    stagetype = mriStep_StageType(step_mem->MRIC, is);
    if ((stagetype == MRISTAGE_DIRK_FAST) || (stagetype == MRISTAGE_DIRK_NOFAST))
      step_mem->implicit = SUNTRUE;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetTable:

  Specifies to use a customized Butcher table for the slow
  portion of the system.
  ---------------------------------------------------------------*/
int MRIStepSetTable(void *arkode_mem, int q, ARKodeButcherTable B)
{
  int retval, is, stagetype;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  sunindextype Tlrw, Tliw;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetTable",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check for illegal inputs */
  if (B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetTable", MSG_ARK_NO_MEM);
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and coupling structure */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  MRIStepCoupling_Space(step_mem->MRIC, &Tliw, &Tlrw);
  MRIStepCoupling_Free(step_mem->MRIC);
  step_mem->MRIC = NULL;
  ark_mem->liw -= Tliw;
  ark_mem->lrw -= Tlrw;

  /* construct the MRI coupling table from B */
  step_mem->MRIC = MRIStepCoupling_MIStoMRI(B, q, 0);
  if (step_mem->MRIC == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                    "MRIStepSetTable",
                    "An error occurred in constructing coupling table.");
    return(ARK_MEM_FAIL);
  }

  /* set the relevant parameters */
  step_mem->stages = step_mem->MRIC->stages;
  step_mem->q = step_mem->MRIC->q;
  step_mem->p = step_mem->MRIC->p;

  /* determine whether fs requires an implicit solver */
  step_mem->implicit = SUNFALSE;
  for (is=0; is<step_mem->stages; is++) {
    stagetype = mriStep_StageType(step_mem->MRIC, is);
    if ((stagetype == MRISTAGE_DIRK_FAST) || (stagetype == MRISTAGE_DIRK_NOFAST))
      step_mem->implicit = SUNTRUE;
  }

  /* re-attach internal error weight function if necessary */
  if (step_mem->implicit) {
    if (!ark_mem->user_efun) {
      if (ark_mem->itol == ARK_SV && ark_mem->Vabstol != NULL)
        retval = arkSVtolerances(ark_mem, ark_mem->reltol, ark_mem->Vabstol);
      else
        retval = arkSStolerances(ark_mem, ark_mem->reltol, ark_mem->Sabstol);
      if (retval != ARK_SUCCESS) return(retval);
    }
  }

  MRIStepCoupling_Space(step_mem->MRIC, &Tliw, &Tlrw);
  ark_mem->liw += Tliw;
  ark_mem->lrw += Tlrw;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetTableNum:

  Specifies to use a pre-existing Butcher or MRI table for the
  slow portion of the problem, based on the integer flag
  indicating the ERK, DIRK or MRI table (see arkode_butcher_erk.h,
  arkode_butcher_dirk.h, and arkode_mristep.h for allowable values).
  ---------------------------------------------------------------*/
int MRIStepSetTableNum(void *arkode_mem, int itable)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  ARKodeButcherTable B;
  MRIStepCoupling C;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetTableNum",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get Butcher table based on argument */
  B = NULL;
  if (itable >= MIN_ERK_NUM && itable <= MAX_ERK_NUM) {
    B = ARKodeButcherTable_LoadERK(itable);
  } else if (itable >= MIN_DIRK_NUM && itable <= MAX_DIRK_NUM) {
    B = ARKodeButcherTable_LoadDIRK(itable);
  }

  /* set Butcher table into MRIStep -- FOR NOW ASSUME 2nd-ORDER COUPLING */
  if (B != NULL) {
    retval = MRIStepSetTable(arkode_mem, SUNMIN(B->q,2), B);
    ARKodeButcherTable_Free(B);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepSetTableNum",
                      "An error occurred in constructing coupling table.");
      return(ARK_MEM_FAIL);
    }
    return(ARK_SUCCESS);
  }

  /* if we've made it this far, then itable is either an MRI table,
     or is an illegal value */
  if (itable >= MIN_MRI_NUM && itable <= MAX_MRI_NUM) {
    C = NULL;
    C = MRIStepCoupling_LoadTable(itable);
    if (C == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepSetTableNum",
                      "An error occurred in constructing coupling table.");
      return(ARK_MEM_FAIL);
    }
    retval = MRIStepSetCoupling(arkode_mem, C);
    MRIStepCoupling_Free(C);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepSetTableNum",
                      "An error occurred in constructing coupling table.");
      return(ARK_MEM_FAIL);
    }
    return(ARK_SUCCESS);
  }

  /* if we've made it this far, then itable was an illegal value */
  arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                  "MRIStepSetTableNum",
                  "Illegal MRI table number");
  return(ARK_ILL_INPUT);
}

/*---------------------------------------------------------------
  MRIStepSetPreInnerFn:

  Sets the user-supplied function called BEFORE the inner evolve
  ---------------------------------------------------------------*/
int MRIStepSetPreInnerFn(void *arkode_mem, MRIStepPreInnerFn prefn)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set pre inner evolve function */
  step_mem->pre_inner_evolve = prefn;

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepSetPostInnerFn:

  Sets the user-supplied function called AFTER the inner evolve
  ---------------------------------------------------------------*/
int MRIStepSetPostInnerFn(void *arkode_mem, MRIStepPostInnerFn postfn)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set pre inner evolve function */
  step_mem->post_inner_evolve = postfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetFixedStep:

  Wrapper for generic arkSetFixedStep routine.  Additionally
  enforces current MRIStep constraint for fixed time-stepping.
  ---------------------------------------------------------------*/
int MRIStepSetFixedStep(void *arkode_mem, realtype hsfixed)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetFixedStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (hsfixed == ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepSetFixedStep",
                    "MIRStep does not support adaptive steps at this time.");
    return(ARK_ILL_INPUT);
  }

  /* call generic routine for remaining work */
  return(arkSetFixedStep(ark_mem, hsfixed));
}


/*---------------------------------------------------------------
  MRIStepSetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int MRIStepSetNonlinCRDown(void *arkode_mem, realtype crdown)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetNonlinCRDown",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (crdown <= ZERO) {
    step_mem->crdown = CRDOWN;
  } else {
    step_mem->crdown = crdown;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int MRIStepSetNonlinRDiv(void *arkode_mem, realtype rdiv)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetNonlinRDiv",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (rdiv <= ZERO) {
    step_mem->rdiv = RDIV;
  } else {
    step_mem->rdiv = rdiv;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply
  a reset to the default.
  ---------------------------------------------------------------*/
int MRIStepSetDeltaGammaMax(void *arkode_mem, realtype dgmax)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDeltaGammaMax",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (dgmax <= ZERO) {
    step_mem->dgmax = DGMAX;
  } else {
    step_mem->dgmax = dgmax;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetLSetupFrequency:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the frequency for calling lsetup;
  negative values imply recomputation of lsetup at each nonlinear
  solve; a zero value implies a reset to the default.
  ---------------------------------------------------------------*/
int MRIStepSetLSetupFrequency(void *arkode_mem, int msbp)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetLSetupFrequency",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (msbp == 0) {
    step_mem->msbp = MSBP;
  } else {
    step_mem->msbp = msbp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.
  Non-default choices are {1,2,3,4}, all others will use default
  (trivial) predictor.
  ---------------------------------------------------------------*/
int MRIStepSetPredictorMethod(void *arkode_mem, int pred_method)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetPredictorMethod",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set parameter */
  step_mem->predictor = pred_method;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int MRIStepSetMaxNonlinIters(void *arkode_mem, int maxcor)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetMaxNonlinIters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return error message if no NLS module is present */
  if (step_mem->NLS == NULL) {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::MRIStep",
                    "MRIStepSetMaxNonlinIters",
                    "No SUNNonlinearSolver object is present");
    return(ARK_ILL_INPUT);
  }

  /* argument <= 0 sets default, otherwise set input */
  if (maxcor <= 0) {
    step_mem->maxcor = MAXCOR;
  } else {
    step_mem->maxcor = maxcor;
  }

  /* send argument to NLS structure */
  retval = SUNNonlinSolSetMaxIters(step_mem->NLS, step_mem->maxcor);
  if (retval != SUN_NLS_SUCCESS) {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::MRIStep",
                    "MRIStepSetMaxNonlinIters",
                    "Error setting maxcor in SUNNonlinearSolver object");
    return(ARK_NLS_OP_ERR);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int MRIStepSetNonlinConvCoef(void *arkode_mem, realtype nlscoef)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetNonlinConvCoef",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) {
    step_mem->nlscoef = NLSCOEF;
  } else {
    step_mem->nlscoef = nlscoef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetStagePredictFn:  Specifies a user-provided step
  predictor function having type ARKStagePredictFn.  A
  NULL input function disables calls to this routine.
  ---------------------------------------------------------------*/
int MRIStepSetStagePredictFn(void *arkode_mem,
                             ARKStagePredictFn PredictStage)
{
  ARKodeMem        ark_mem;
  ARKodeMRIStepMem step_mem;
  int              retval;

  /* access ARKodeMRIStepMem structure and set function pointer */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetStagePredictFn",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* override predictor method 5 if non-NULL PredictStage is supplied */
  if ((step_mem->predictor == 5) && (PredictStage != NULL)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepSetStagePredictFn",
                    "User-supplied predictor is incompatible with predictor method 5");
    return(ARK_ILL_INPUT);
  }

  step_mem->stage_predict = PredictStage;
  return(ARK_SUCCESS);
}



/*===============================================================
  MRIStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepGetLastInnerStepFlag:

  Returns the last return value from the inner stepper.
  ---------------------------------------------------------------*/
int MRIStepGetLastInnerStepFlag(void *arkode_mem, int *flag)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetLastInnerStepFlag",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get the last return value from the inner stepper */
  *flag = step_mem->inner_retval;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetCurrentGamma: Returns the current value of gamma
  ---------------------------------------------------------------*/
int MRIStepGetCurrentGamma(void *arkode_mem, realtype *gamma)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  retval = mriStep_AccessStepMem(arkode_mem, NULL, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);
  *gamma = step_mem->gamma;
  return(retval);
}


/*---------------------------------------------------------------
  MRIStepGetNumRhsEvals:

  Returns the current number of calls to fs and ff
  ---------------------------------------------------------------*/
int MRIStepGetNumRhsEvals(void *arkode_mem, long int *nfs_evals)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetNumRhsEvals",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get number of fs evals from step_mem */
  *nfs_evals = step_mem->nfs;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int MRIStepGetNumLinSolvSetups(void *arkode_mem, long int *nlinsetups)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetNumLinSolvSetups",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get value from step_mem */
  *nlinsetups = step_mem->nsetups;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
int MRIStepGetNumNonlinSolvIters(void *arkode_mem, long int *nniters)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetNumNonlinSolvIters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  *nniters = step_mem->nls_iters;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int MRIStepGetNumNonlinSolvConvFails(void *arkode_mem, long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetNumNonlinSolvConvFails",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set output from step_mem */
  *nncfails = ark_mem->ncfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
int MRIStepGetNonlinSolvStats(void *arkode_mem, long int *nniters,
                              long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetNonlinSolvStats",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  *nniters  = step_mem->nls_iters;
  *nncfails = ark_mem->ncfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetCurrentCoupling:

  Sets pointer to the slow coupling structure currently in use.
  ---------------------------------------------------------------*/
int MRIStepGetCurrentCoupling(void *arkode_mem, MRIStepCoupling *MRIC)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetCurrentCoupling",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get coupling structure from step_mem */
  *MRIC = step_mem->MRIC;

  return(ARK_SUCCESS);
}


/*===============================================================
  MRIStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int MRIStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepWriteParameters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* output ARKode infrastructure parameters first */
  retval = arkWriteParameters(arkode_mem, fp);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepWriteParameters",
                    "Error writing ARKode infrastructure parameters");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepWriteCoupling:

  Outputs coupling structure to the provided file pointer.
  ---------------------------------------------------------------*/
int MRIStepWriteCoupling(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepWriteCoupling",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check that coupling structure is non-NULL (otherwise report error) */
  if (step_mem->MRIC == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepWriteCoupling", "Coupling structure is NULL");
    return(ARK_MEM_NULL);
  }

  /* write coupling structure to specified file */
  fprintf(fp, "\nMRIStep coupling structure:\n");
  MRIStepCoupling_Write(step_mem->MRIC, fp);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
