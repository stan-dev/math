/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for the optional input and
 * output functions for the ARKode ARKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here,
 * with slightly different names.  The code transition will be
 * minimal, but the documentation changes will be significant.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_arkstep_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  ARKStep Optional input functions (wrappers for generic ARKode
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int ARKStepSetDenseOrder(void *arkode_mem, int dord) {
  return(ARKStepSetInterpolantDegree(arkode_mem, dord)); }
int ARKStepSetInterpolantDegree(void *arkode_mem, int degree) {
  if (degree < 0) degree = ARK_INTERP_MAX_DEGREE;
  return(arkSetInterpolantDegree(arkode_mem, degree)); }
int ARKStepSetInterpolantType(void *arkode_mem, int itype) {
  return(arkSetInterpolantType(arkode_mem, itype)); }
int ARKStepSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                           void *eh_data) {
  return(arkSetErrHandlerFn(arkode_mem, ehfun, eh_data)); }
int ARKStepSetErrFile(void *arkode_mem, FILE *errfp) {
  return(arkSetErrFile(arkode_mem, errfp)); }
int ARKStepSetDiagnostics(void *arkode_mem, FILE *diagfp) {
  return(arkSetDiagnostics(arkode_mem, diagfp)); }
int ARKStepSetMaxNumSteps(void *arkode_mem, long int mxsteps) {
  return(arkSetMaxNumSteps(arkode_mem, mxsteps)); }
int ARKStepSetMaxHnilWarns(void *arkode_mem, int mxhnil) {
  return(arkSetMaxHnilWarns(arkode_mem, mxhnil)); }
int ARKStepSetInitStep(void *arkode_mem, realtype hin) {
  return(arkSetInitStep(arkode_mem, hin)); }
int ARKStepSetMinStep(void *arkode_mem, realtype hmin) {
  return(arkSetMinStep(arkode_mem, hmin)); }
int ARKStepSetMaxStep(void *arkode_mem, realtype hmax) {
  return(arkSetMaxStep(arkode_mem, hmax)); }
int ARKStepSetStopTime(void *arkode_mem, realtype tstop) {
  return(arkSetStopTime(arkode_mem, tstop)); }
int ARKStepSetRootDirection(void *arkode_mem, int *rootdir) {
  return(arkSetRootDirection(arkode_mem, rootdir)); }
int ARKStepSetNoInactiveRootWarn(void *arkode_mem) {
  return(arkSetNoInactiveRootWarn(arkode_mem)); }
int ARKStepSetConstraints(void *arkode_mem, N_Vector constraints) {
  return(arkSetConstraints(arkode_mem, constraints)); }
int ARKStepSetMaxNumConstrFails(void *arkode_mem, int maxfails) {
  return(arkSetMaxNumConstrFails(arkode_mem, maxfails)); }
int ARKStepSetPostprocessStepFn(void *arkode_mem,
                                ARKPostProcessFn ProcessStep) {
  return(arkSetPostprocessStepFn(arkode_mem, ProcessStep)); }
int ARKStepSetPostprocessStageFn(void *arkode_mem,
                                 ARKPostProcessFn ProcessStage) {
  return(arkSetPostprocessStageFn(arkode_mem, ProcessStage)); }
int ARKStepSetCFLFraction(void *arkode_mem, realtype cfl_frac) {
  return(arkSetCFLFraction(arkode_mem, cfl_frac)); }
int ARKStepSetSafetyFactor(void *arkode_mem, realtype safety) {
  return(arkSetSafetyFactor(arkode_mem, safety)); }
int ARKStepSetErrorBias(void *arkode_mem, realtype bias) {
  return(arkSetErrorBias(arkode_mem, bias)); }
int ARKStepSetMaxGrowth(void *arkode_mem, realtype mx_growth) {
  return(arkSetMaxGrowth(arkode_mem, mx_growth)); }
int ARKStepSetMinReduction(void *arkode_mem, realtype eta_min) {
  return(arkSetMinReduction(arkode_mem, eta_min)); }
int ARKStepSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub) {
  return(arkSetFixedStepBounds(arkode_mem, lb, ub)); }
int ARKStepSetAdaptivityMethod(void *arkode_mem, int imethod, int idefault,
                               int pq, realtype adapt_params[3]) {
  return(arkSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, adapt_params)); }
int ARKStepSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun, void *h_data) {
  return(arkSetAdaptivityFn(arkode_mem, hfun, h_data)); }
int ARKStepSetMaxFirstGrowth(void *arkode_mem, realtype etamx1) {
  return(arkSetMaxFirstGrowth(arkode_mem, etamx1)); }
int ARKStepSetMaxEFailGrowth(void *arkode_mem, realtype etamxf) {
  return(arkSetMaxEFailGrowth(arkode_mem, etamxf)); }
int ARKStepSetSmallNumEFails(void *arkode_mem, int small_nef) {
  return(arkSetSmallNumEFails(arkode_mem, small_nef)); }
int ARKStepSetMaxCFailGrowth(void *arkode_mem, realtype etacf) {
  return(arkSetMaxCFailGrowth(arkode_mem, etacf)); }
int ARKStepSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab, void *estab_data) {
  return(arkSetStabilityFn(arkode_mem, EStab, estab_data)); }
int ARKStepSetMaxErrTestFails(void *arkode_mem, int maxnef) {
  return(arkSetMaxErrTestFails(arkode_mem, maxnef)); }
int ARKStepSetMaxConvFails(void *arkode_mem, int maxncf) {
  return(arkSetMaxConvFails(arkode_mem, maxncf)); }
int ARKStepSetFixedStep(void *arkode_mem, realtype hfixed) {
  return(arkSetFixedStep(arkode_mem, hfixed)); }


/*---------------------------------------------------------------
  These wrappers for ARKLs module 'set' routines all are
  documented in arkode_arkstep.h.
  ---------------------------------------------------------------*/
int ARKStepSetLinearSolver(void *arkode_mem, SUNLinearSolver LS,
                           SUNMatrix A) {
  return(arkLSSetLinearSolver(arkode_mem, LS, A)); }
int ARKStepSetMassLinearSolver(void *arkode_mem, SUNLinearSolver LS,
                               SUNMatrix M, booleantype time_dep) {
  return(arkLSSetMassLinearSolver(arkode_mem, LS, M, time_dep)); }
int ARKStepSetJacFn(void *arkode_mem, ARKLsJacFn jac) {
  return(arkLSSetJacFn(arkode_mem, jac)); }
int ARKStepSetMassFn(void *arkode_mem, ARKLsMassFn mass) {
  return(arkLSSetMassFn(arkode_mem, mass)); }
int ARKStepSetJacEvalFrequency(void *arkode_mem, long int msbj) {
  return(arkLSSetJacEvalFrequency(arkode_mem, msbj)); }
int ARKStepSetLinearSolutionScaling(void *arkode_mem, booleantype onoff) {
  return(arkLSSetLinearSolutionScaling(arkode_mem, onoff)); }
int ARKStepSetEpsLin(void *arkode_mem, realtype eplifac) {
  return(arkLSSetEpsLin(arkode_mem, eplifac)); }
int ARKStepSetMassEpsLin(void *arkode_mem, realtype eplifac) {
  return(arkLSSetMassEpsLin(arkode_mem, eplifac)); }
int ARKStepSetLSNormFactor(void *arkode_mem, realtype nrmfac) {
  return(arkLSSetNormFactor(arkode_mem, nrmfac)); }
int ARKStepSetMassLSNormFactor(void *arkode_mem, realtype nrmfac) {
  return(arkLSSetMassNormFactor(arkode_mem, nrmfac)); }
int ARKStepSetPreconditioner(void *arkode_mem, ARKLsPrecSetupFn psetup,
                             ARKLsPrecSolveFn psolve) {
  return(arkLSSetPreconditioner(arkode_mem, psetup, psolve)); }
int ARKStepSetMassPreconditioner(void *arkode_mem, ARKLsMassPrecSetupFn psetup,
                                 ARKLsMassPrecSolveFn psolve) {
  return(arkLSSetMassPreconditioner(arkode_mem, psetup, psolve)); }
int ARKStepSetJacTimes(void *arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                       ARKLsJacTimesVecFn jtimes) {
  return(arkLSSetJacTimes(arkode_mem, jtsetup, jtimes)); }
int ARKStepSetJacTimesRhsFn(void *arkode_mem, ARKRhsFn jtimesRhsFn) {
  return(arkLSSetJacTimesRhsFn(arkode_mem, jtimesRhsFn)); }
int ARKStepSetMassTimes(void *arkode_mem, ARKLsMassTimesSetupFn msetup,
                        ARKLsMassTimesVecFn mtimes, void *mtimes_data) {
  return(arkLSSetMassTimes(arkode_mem, msetup, mtimes, mtimes_data)); }
int ARKStepSetLinSysFn(void *arkode_mem, ARKLsLinSysFn linsys) {
  return(arkLSSetLinSysFn(arkode_mem, linsys)); }

/* deprecated */
int ARKStepSetMaxStepsBetweenJac(void *arkode_mem, long int msbj) {
  return(arkLSSetJacEvalFrequency(arkode_mem, msbj)); }


/*===============================================================
  ARKStep Optional output functions (wrappers for generic ARKode
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int ARKStepGetNumStepAttempts(void *arkode_mem, long int *nstep_attempts) {
  return(arkGetNumStepAttempts(arkode_mem, nstep_attempts)); }
int ARKStepGetNumSteps(void *arkode_mem, long int *nsteps) {
  return(arkGetNumSteps(arkode_mem, nsteps)); }
int ARKStepGetActualInitStep(void *arkode_mem, realtype *hinused) {
  return(arkGetActualInitStep(arkode_mem, hinused)); }
int ARKStepGetLastStep(void *arkode_mem, realtype *hlast) {
  return(arkGetLastStep(arkode_mem, hlast)); }
int ARKStepGetCurrentStep(void *arkode_mem, realtype *hcur) {
  return(arkGetCurrentStep(arkode_mem, hcur)); }
int ARKStepGetCurrentTime(void *arkode_mem, realtype *tcur) {
  return(arkGetCurrentTime(arkode_mem, tcur)); }
int ARKStepGetCurrentState(void *arkode_mem, N_Vector *state) {
  return(arkGetCurrentState(arkode_mem, state)); }
int ARKStepGetTolScaleFactor(void *arkode_mem, realtype *tolsfact) {
  return(arkGetTolScaleFactor(arkode_mem, tolsfact)); }
int ARKStepGetErrWeights(void *arkode_mem, N_Vector eweight) {
  return(arkGetErrWeights(arkode_mem, eweight)); }
int ARKStepGetResWeights(void *arkode_mem, N_Vector rweight) {
  return(arkGetResWeights(arkode_mem, rweight)); }
int ARKStepGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw) {
  return(arkGetWorkSpace(arkode_mem, lenrw, leniw)); }
int ARKStepGetNumGEvals(void *arkode_mem, long int *ngevals) {
  return(arkGetNumGEvals(arkode_mem, ngevals)); }
int ARKStepGetRootInfo(void *arkode_mem, int *rootsfound) {
  return(arkGetRootInfo(arkode_mem, rootsfound)); }
int ARKStepGetStepStats(void *arkode_mem, long int *nsteps,
                        realtype *hinused, realtype *hlast,
                        realtype *hcur, realtype *tcur) {
  return(arkGetStepStats(arkode_mem, nsteps, hinused, hlast, hcur, tcur)); }
int ARKStepGetNumConstrFails(void *arkode_mem, long int *nconstrfails) {
  return(arkGetNumConstrFails(arkode_mem, nconstrfails)); }
int ARKStepGetNumExpSteps(void *arkode_mem, long int *nsteps) {
  return(arkGetNumExpSteps(arkode_mem, nsteps)); }
int ARKStepGetNumAccSteps(void *arkode_mem, long int *nsteps) {
  return(arkGetNumAccSteps(arkode_mem, nsteps)); }
int ARKStepGetNumErrTestFails(void *arkode_mem, long int *netfails) {
  return(arkGetNumErrTestFails(arkode_mem, netfails)); }
char *ARKStepGetReturnFlagName(long int flag) {
  return(arkGetReturnFlagName(flag)); }

/*---------------------------------------------------------------
  These wrappers for ARKLs module 'get' routines all are
  documented in arkode_arkstep.h.
  ---------------------------------------------------------------*/
int ARKStepGetLinWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS) {
  return(arkLSGetWorkSpace(arkode_mem, lenrwLS, leniwLS)); }
int ARKStepGetNumJacEvals(void *arkode_mem, long int *njevals) {
  return(arkLSGetNumJacEvals(arkode_mem, njevals)); }
int ARKStepGetNumPrecEvals(void *arkode_mem, long int *npevals) {
  return(arkLSGetNumPrecEvals(arkode_mem, npevals)); }
int ARKStepGetNumPrecSolves(void *arkode_mem, long int *npsolves) {
  return(arkLSGetNumPrecSolves(arkode_mem, npsolves)); }
int ARKStepGetNumLinIters(void *arkode_mem, long int *nliters) {
  return(arkLSGetNumLinIters(arkode_mem, nliters)); }
int ARKStepGetNumLinConvFails(void *arkode_mem, long int *nlcfails) {
  return(arkLSGetNumConvFails(arkode_mem, nlcfails)); }
int ARKStepGetNumJTSetupEvals(void *arkode_mem, long int *njtsetups) {
  return(arkLSGetNumJTSetupEvals(arkode_mem, njtsetups)); }
int ARKStepGetNumJtimesEvals(void *arkode_mem, long int *njvevals) {
  return(arkLSGetNumJtimesEvals(arkode_mem, njvevals)); }
int ARKStepGetNumLinRhsEvals(void *arkode_mem, long int *nfevalsLS) {
  return(arkLSGetNumRhsEvals(arkode_mem, nfevalsLS)); }
int ARKStepGetLastLinFlag(void *arkode_mem, long int *flag) {
  return(arkLSGetLastFlag(arkode_mem, flag)); }

int ARKStepGetMassWorkSpace(void *arkode_mem, long int *lenrwMLS, long int *leniwMLS) {
  return(arkLSGetMassWorkSpace(arkode_mem, lenrwMLS, leniwMLS)); }
int ARKStepGetNumMassSetups(void *arkode_mem, long int *nmsetups) {
  return(arkLSGetNumMassSetups(arkode_mem, nmsetups)); }
int ARKStepGetNumMassMultSetups(void *arkode_mem, long int *nmvsetups) {
  return(arkLSGetNumMassMatvecSetups(arkode_mem, nmvsetups)); }
int ARKStepGetNumMassMult(void *arkode_mem, long int *nmvevals) {
  return(arkLSGetNumMassMult(arkode_mem, nmvevals)); }
int ARKStepGetNumMassSolves(void *arkode_mem, long int *nmsolves) {
  return(arkLSGetNumMassSolves(arkode_mem, nmsolves)); }
int ARKStepGetNumMassPrecEvals(void *arkode_mem, long int *nmpevals) {
  return(arkLSGetNumMassPrecEvals(arkode_mem, nmpevals)); }
int ARKStepGetNumMassPrecSolves(void *arkode_mem, long int *nmpsolves) {
  return(arkLSGetNumMassPrecSolves(arkode_mem, nmpsolves)); }
int ARKStepGetNumMassIters(void *arkode_mem, long int *nmiters) {
  return(arkLSGetNumMassIters(arkode_mem, nmiters)); }
int ARKStepGetNumMassConvFails(void *arkode_mem, long int *nmcfails) {
  return(arkLSGetNumMassConvFails(arkode_mem, nmcfails)); }
int ARKStepGetNumMTSetups(void *arkode_mem, long int *nmtsetups) {
  return(arkLSGetNumMTSetups(arkode_mem, nmtsetups)); }
int ARKStepGetCurrentMassMatrix(void *arkode_mem, SUNMatrix *M) {
  return(arkLSGetCurrentMassMatrix(arkode_mem, M)); }
int ARKStepGetLastMassFlag(void *arkode_mem, long int *flag) {
  return(arkLSGetLastMassFlag(arkode_mem, flag)); }
char *ARKStepGetLinReturnFlagName(long int flag) {
  return(arkLSGetReturnFlagName(flag)); }



/*===============================================================
  ARKStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepSetUserData:

  Wrapper for generic arkSetUserData and arkLSSetUserData
  routines.
  ---------------------------------------------------------------*/
int ARKStepSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem        ark_mem;
  ARKodeARKStepMem step_mem;
  int              retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetUserData",
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

  /* set user data in ARKodeLSMass mem */
  if (step_mem->mass_mem != NULL) {
    retval = arkLSSetMassUserData(arkode_mem, user_data);
    if (retval != ARKLS_SUCCESS) return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetDefaults:

  Resets all ARKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.  Also leaves alone any data
  structures/options related to the ARKode infrastructure itself
  (e.g., root-finding and post-process step).
  ---------------------------------------------------------------*/
int ARKStepSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Set default ARKode infrastructure parameters */
  retval = arkSetDefaults(ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetDefaults",
                    "Error setting ARKode infrastructure defaults");
    return(retval);
  }

  /* Set default values for integrator optional inputs */
  step_mem->q                = Q_DEFAULT;      /* method order */
  step_mem->p                = 0;              /* embedding order */
  step_mem->predictor        = 0;              /* trivial predictor */
  step_mem->linear           = SUNFALSE;       /* nonlinear problem */
  step_mem->linear_timedep   = SUNTRUE;        /* dfi/dy depends on t */
  step_mem->explicit         = SUNTRUE;        /* fe(t,y) will be used */
  step_mem->implicit         = SUNTRUE;        /* fi(t,y) will be used */
  step_mem->maxcor           = MAXCOR;         /* max nonlinear iters/stage */
  step_mem->nlscoef          = NLSCOEF;        /* nonlinear tolerance coefficient */
  step_mem->crdown           = CRDOWN;         /* nonlinear convergence estimate coeff. */
  step_mem->rdiv             = RDIV;           /* nonlinear divergence tolerance */
  step_mem->dgmax            = DGMAX;          /* max step change before recomputing J or P */
  step_mem->msbp             = MSBP;           /* max steps between updates to J or P */
  step_mem->stages           = 0;              /* no stages */
  step_mem->istage           = 0;              /* current stage */
  step_mem->Be               = NULL;           /* no Butcher tables */
  step_mem->Bi               = NULL;
  step_mem->NLS              = NULL;           /* no nonlinear solver object */
  step_mem->jcur             = SUNFALSE;
  step_mem->convfail         = ARK_NO_FAILURES;
  step_mem->stage_predict    = NULL;           /* no user-supplied stage predictor */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetOptimalParams:

  Sets all adaptivity and solver parameters to our 'best guess'
  values, for a given ARKStep integration method (ERK, DIRK, ARK),
  a given method order, and a given nonlinear solver type.  Should
  only be called after the method order, solver, and integration
  method have been set, and only if time step adaptivity is
  enabled.
  ---------------------------------------------------------------*/
int ARKStepSetOptimalParams(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetOptimalParams",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (ark_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetOptimalParams",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = ark_mem->hadapt_mem;

  /* Choose values based on method, order */

  /*    explicit */
  if (step_mem->explicit && !step_mem->implicit) {
    hadapt_mem->imethod = ARK_ADAPT_PI;
    hadapt_mem->safety  = RCONST(0.99);
    hadapt_mem->bias    = RCONST(1.2);
    hadapt_mem->growth  = RCONST(25.0);
    hadapt_mem->k1      = RCONST(0.8);
    hadapt_mem->k2      = RCONST(0.31);
    hadapt_mem->etamxf  = RCONST(0.3);

  /*    implicit */
  } else if (step_mem->implicit && !step_mem->explicit) {
    switch (step_mem->q) {
    case 2:   /* just use standard defaults since better ones unknown */
      hadapt_mem->imethod   = ARK_ADAPT_PID;
      hadapt_mem->safety    = SAFETY;
      hadapt_mem->bias      = BIAS;
      hadapt_mem->growth    = GROWTH;
      hadapt_mem->etamxf    = ETAMXF;
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef     = RCONST(0.001);
      step_mem->maxcor      = 5;
      step_mem->crdown      = CRDOWN;
      step_mem->rdiv        = RDIV;
      step_mem->dgmax       = DGMAX;
      step_mem->msbp        = MSBP;
      break;
    case 3:
      hadapt_mem->imethod   = ARK_ADAPT_I;
      hadapt_mem->safety    = RCONST(0.957);
      hadapt_mem->bias      = RCONST(1.9);
      hadapt_mem->growth    = RCONST(17.6);
      hadapt_mem->etamxf    = RCONST(0.45);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef     = RCONST(0.22);
      step_mem->crdown      = RCONST(0.17);
      step_mem->rdiv        = RCONST(2.3);
      step_mem->dgmax       = RCONST(0.19);
      step_mem->msbp        = 60;
      break;
    case 4:
      hadapt_mem->imethod   = ARK_ADAPT_PID;
      hadapt_mem->safety    = RCONST(0.988);
      hadapt_mem->bias      = RCONST(1.2);
      hadapt_mem->growth    = RCONST(31.5);
      hadapt_mem->k1        = RCONST(0.535);
      hadapt_mem->k2        = RCONST(0.209);
      hadapt_mem->k3        = RCONST(0.148);
      hadapt_mem->etamxf    = RCONST(0.33);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef     = RCONST(0.24);
      step_mem->crdown      = RCONST(0.26);
      step_mem->rdiv        = RCONST(2.3);
      step_mem->dgmax       = RCONST(0.16);
      step_mem->msbp        = 31;
      break;
    case 5:
      hadapt_mem->imethod   = ARK_ADAPT_PID;
      hadapt_mem->safety    = RCONST(0.937);
      hadapt_mem->bias      = RCONST(3.3);
      hadapt_mem->growth    = RCONST(22.0);
      hadapt_mem->k1        = RCONST(0.56);
      hadapt_mem->k2        = RCONST(0.338);
      hadapt_mem->k3        = RCONST(0.14);
      hadapt_mem->etamxf    = RCONST(0.44);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef     = RCONST(0.25);
      step_mem->crdown      = RCONST(0.4);
      step_mem->rdiv        = RCONST(2.3);
      step_mem->dgmax       = RCONST(0.32);
      step_mem->msbp        = 31;
      break;
    }

  /*    imex */
  } else {
    switch (step_mem->q) {
    case 3:
      hadapt_mem->imethod   = ARK_ADAPT_PID;
      hadapt_mem->safety    = RCONST(0.965);
      hadapt_mem->bias      = RCONST(1.42);
      hadapt_mem->growth    = RCONST(28.7);
      hadapt_mem->k1        = RCONST(0.54);
      hadapt_mem->k2        = RCONST(0.36);
      hadapt_mem->k3        = RCONST(0.14);
      hadapt_mem->etamxf    = RCONST(0.46);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef     = RCONST(0.22);
      step_mem->crdown      = RCONST(0.17);
      step_mem->rdiv        = RCONST(2.3);
      step_mem->dgmax       = RCONST(0.19);
      step_mem->msbp        = 60;
      break;
    case 4:
      hadapt_mem->imethod   = ARK_ADAPT_PID;
      hadapt_mem->safety    = RCONST(0.97);
      hadapt_mem->bias      = RCONST(1.35);
      hadapt_mem->growth    = RCONST(25.0);
      hadapt_mem->k1        = RCONST(0.543);
      hadapt_mem->k2        = RCONST(0.297);
      hadapt_mem->k3        = RCONST(0.14);
      hadapt_mem->etamxf    = RCONST(0.47);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef     = RCONST(0.24);
      step_mem->crdown      = RCONST(0.26);
      step_mem->rdiv        = RCONST(2.3);
      step_mem->dgmax       = RCONST(0.16);
      step_mem->msbp        = 31;
      break;
    case 5:
      hadapt_mem->imethod   = ARK_ADAPT_PI;
      hadapt_mem->safety    = RCONST(0.993);
      hadapt_mem->bias      = RCONST(1.15);
      hadapt_mem->growth    = RCONST(28.5);
      hadapt_mem->k1        = RCONST(0.8);
      hadapt_mem->k2        = RCONST(0.35);
      hadapt_mem->etamxf    = RCONST(0.3);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef     = RCONST(0.25);
      step_mem->crdown      = RCONST(0.4);
      step_mem->rdiv        = RCONST(2.3);
      step_mem->dgmax       = RCONST(0.32);
      step_mem->msbp        = 31;
      break;
    }

  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetOrder:

  Specifies the method order

  ** Note in documentation that this should not be called along
  with ARKStepSetTable or ARKStepSetTableNum.  This routine
  is used to specify a desired method order using default Butcher
  tables, whereas any user-supplied table will have their own
  order associated with them.
  ---------------------------------------------------------------*/
int ARKStepSetOrder(void *arkode_mem, int ord)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetOrder",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) {
    step_mem->q = Q_DEFAULT;
  } else {
    step_mem->q = ord;
  }

  /* clear Butcher tables, since user is requesting a change in method
     or a reset to defaults.  Tables will be set in ARKInitialSetup. */
  step_mem->stages = 0;
  step_mem->istage = 0;
  step_mem->p = 0;

  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Be);
  step_mem->Be = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Bi);
  step_mem->Bi = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetLinear:

  Specifies that the implicit portion of the problem is linear,
  and to tighten the linear solver tolerances while taking only
  one Newton iteration.  DO NOT USE IN COMBINATION WITH THE
  FIXED-POINT SOLVER.  Automatically tightens DeltaGammaMax
  to ensure that step size changes cause Jacobian recomputation.

  The argument should be 1 or 0, where 1 indicates that the
  Jacobian of fi with respect to y depends on time, and
  0 indicates that it is not time dependent.  Alternately, when
  using an iterative linear solver this flag denotes time
  dependence of the preconditioner.
  ---------------------------------------------------------------*/
int ARKStepSetLinear(void *arkode_mem, int timedepend)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetLinear",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set parameters */
  step_mem->linear = SUNTRUE;
  step_mem->linear_timedep = (timedepend == 1);
  step_mem->dgmax = RCONST(100.0)*UNIT_ROUNDOFF;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetNonlinear:

  Specifies that the implicit portion of the problem is nonlinear.
  Used to undo a previous call to ARKStepSetLinear.  Automatically
  loosens DeltaGammaMax back to default value.
  ---------------------------------------------------------------*/
int ARKStepSetNonlinear(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetNonlinear",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set parameters */
  step_mem->linear = SUNFALSE;
  step_mem->linear_timedep = SUNTRUE;
  step_mem->dgmax = DGMAX;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetExplicit:

  Specifies that the implicit portion of the problem is disabled,
  and to use an explicit RK method.
  ---------------------------------------------------------------*/
int ARKStepSetExplicit(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetExplicit",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* ensure that fe is defined */
  if (step_mem->fe == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetExplicit", MSG_ARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  step_mem->explicit = SUNTRUE;
  step_mem->implicit = SUNFALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetImplicit:

  Specifies that the explicit portion of the problem is disabled,
  and to use an implicit RK method.
  ---------------------------------------------------------------*/
int ARKStepSetImplicit(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetImplicit",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* ensure that fi is defined */
  if (step_mem->fi == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetImplicit", MSG_ARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  step_mem->implicit = SUNTRUE;
  step_mem->explicit = SUNFALSE;

  /* re-attach internal error weight functions if necessary */
  if (!ark_mem->user_efun) {
    if (ark_mem->itol == ARK_SV && ark_mem->Vabstol != NULL)
      retval = arkSVtolerances(ark_mem, ark_mem->reltol, ark_mem->Vabstol);
    else
      retval = arkSStolerances(ark_mem, ark_mem->reltol, ark_mem->Sabstol);
    if (retval != ARK_SUCCESS) return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetImEx:

  Specifies that the specifies that problem has both implicit and
  explicit parts, and to use an ARK method (this is the default).
  ---------------------------------------------------------------*/
int ARKStepSetImEx(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetImEx",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* ensure that fe and fi are defined */
  if (step_mem->fe == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetImEx", MSG_ARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }
  if (step_mem->fi == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetImEx", MSG_ARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  step_mem->explicit = SUNTRUE;
  step_mem->implicit = SUNTRUE;

  /* re-attach internal error weight functions if necessary */
  if (!ark_mem->user_efun) {
    if (ark_mem->itol == ARK_SV && ark_mem->Vabstol != NULL)
      retval = arkSVtolerances(ark_mem, ark_mem->reltol, ark_mem->Vabstol);
    else
      retval = arkSStolerances(ark_mem, ark_mem->reltol, ark_mem->Sabstol);
    if (retval != ARK_SUCCESS) return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetTables:

  Specifies to use customized Butcher tables for the system.

  If Bi is NULL, then this sets the integrator in 'explicit' mode.

  If Be is NULL, then this sets the integrator in 'implicit' mode.

  Returns ARK_ILL_INPUT if both Butcher tables are not supplied.
  ---------------------------------------------------------------*/
int ARKStepSetTables(void *arkode_mem, int q, int p,
                     ARKodeButcherTable Bi, ARKodeButcherTable Be)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  sunindextype Blrw, Bliw;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetTables",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* check for illegal inputs */
  if ((Bi == NULL) && (Be == NULL)) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetTables",
                    "At least one complete table must be supplied");
    return(ARK_ILL_INPUT);
  }

  /* if both tables are set, check that they have the same number of stages */
  if ((Bi != NULL) && (Be != NULL)) {
    if (Bi->stages != Be->stages) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTables",
                      "Both tables must have the same number of stages");
      return(ARK_ILL_INPUT);
    }
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;

  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Be);
  step_mem->Be = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Bi);
  step_mem->Bi = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  /*
   * determine mode (implicit/explicit/ImEx), and perform appropriate actions
   */

  /* explicit */
  if (Bi == NULL) {

    /* set the relevant parameters (use table q and p) */
    step_mem->stages = Be->stages;
    step_mem->q = Be->q;
    step_mem->p = Be->p;

    /* copy the table in step memory */
    step_mem->Be = ARKodeButcherTable_Copy(Be);
    if (step_mem->Be == NULL) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTables", MSG_ARK_NO_MEM);
      return(ARK_MEM_NULL);
    }

    /* set method as purely explicit */
    retval = ARKStepSetExplicit(arkode_mem);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "ARKStepSetTables",
                      "Error in ARKStepSetExplicit");
      return(retval);
    }

  /* implicit */
  } else if (Be == NULL) {

    /* set the relevant parameters (use table q and p) */
    step_mem->stages = Bi->stages;
    step_mem->q = Bi->q;
    step_mem->p = Bi->p;

    /* copy the table in step memory */
    step_mem->Bi = ARKodeButcherTable_Copy(Bi);
    if (step_mem->Bi == NULL) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTables", MSG_ARK_NO_MEM);
      return(ARK_MEM_NULL);
    }

    /* set method as purely implicit */
    retval = ARKStepSetImplicit(arkode_mem);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "ARKStepSetTables",
                      "Error in ARKStepSetImplicit");
      return(ARK_ILL_INPUT);
    }

  /* ImEx */
  } else {

    /* set the relevant parameters (use input q and p) */
    step_mem->stages = Bi->stages;
    step_mem->q = q;
    step_mem->p = p;

    /* copy the explicit table into step memory */
    step_mem->Be = ARKodeButcherTable_Copy(Be);
    if (step_mem->Be == NULL) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTables", MSG_ARK_NO_MEM);
      return(ARK_MEM_NULL);
    }

    /* copy the implicit table into step memory */
    step_mem->Bi = ARKodeButcherTable_Copy(Bi);
    if (step_mem->Bi == NULL) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTables", MSG_ARK_NO_MEM);
      return(ARK_MEM_NULL);
    }

    /* set method as ImEx */
    retval = ARKStepSetImEx(arkode_mem);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "ARKStepSetTables",
                      "Error in ARKStepSetImEx");
      return(ARK_ILL_INPUT);
    }
  }

  /* note Butcher table space requirements */
  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetTableNum:

  Specifies to use pre-existing Butcher tables for the system,
  based on the integer flags passed to
  ARKodeButcherTable_LoadERK() and ARKodeButcherTable_LoadDIRK()
  within the files arkode_butcher_erk.c and arkode_butcher_dirk.c
  (automatically calls ARKStepSetImEx).

  If either argument is negative (illegal), then this disables the
  corresponding table (e.g. itable = -1  ->  explicit)
  ---------------------------------------------------------------*/
int ARKStepSetTableNum(void *arkode_mem, int itable, int etable)
{
  int flag, retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  sunindextype Blrw, Bliw;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetTableNum",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;

  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Be);
  step_mem->Be = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ARKodeButcherTable_Free(step_mem->Bi);
  step_mem->Bi = NULL;
  ark_mem->liw -= Bliw;
  ark_mem->lrw -= Blrw;

  /* determine mode (implicit/explicit/ImEx), and perform
     appropriate actions  */

  /*     illegal inputs */
  if ((itable < 0) && (etable < 0)) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetTableNum",
                    "At least one valid table number must be supplied");
    return(ARK_ILL_INPUT);


  /* explicit */
  } else if (itable < 0) {

    /* check that argument specifies an explicit table */
    if (etable<MIN_ERK_NUM || etable>MAX_ERK_NUM) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTableNum",
                      "Illegal ERK table number");
      return(ARK_ILL_INPUT);
    }

    /* fill in table based on argument */
    step_mem->Be = ARKodeButcherTable_LoadERK(etable);
    if (step_mem->Be == NULL) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTableNum",
                      "Error setting explicit table with that index");
      return(ARK_ILL_INPUT);
    }
    step_mem->stages = step_mem->Be->stages;
    step_mem->q = step_mem->Be->q;
    step_mem->p = step_mem->Be->p;

    /* set method as purely explicit */
    flag = ARKStepSetExplicit(arkode_mem);
    if (flag != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "ARKStepSetTableNum",
                      "Error in ARKStepSetExplicit");
      return(flag);
    }


  /* implicit */
  } else if (etable < 0) {

    /* check that argument specifies an implicit table */
    if (itable<MIN_DIRK_NUM || itable>MAX_DIRK_NUM) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTableNum",
                      "Illegal IRK table number");
      return(ARK_ILL_INPUT);
    }

    /* fill in table based on argument */
    step_mem->Bi = ARKodeButcherTable_LoadDIRK(itable);
    if (step_mem->Bi == NULL) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTableNum",
                      "Error setting table with that index");
      return(ARK_ILL_INPUT);
    }
    step_mem->stages = step_mem->Bi->stages;
    step_mem->q = step_mem->Bi->q;
    step_mem->p = step_mem->Bi->p;

    /* set method as purely implicit */
    flag = ARKStepSetImplicit(arkode_mem);
    if (flag != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "ARKStepSetTableNum",
                      "Error in ARKStepSetIxplicit");
      return(flag);
    }


  /* ImEx */
  } else {

    /* ensure that tables match */
    if ( !((etable == ARK324L2SA_ERK_4_2_3) && (itable == ARK324L2SA_DIRK_4_2_3)) &&
         !((etable == ARK436L2SA_ERK_6_3_4) && (itable == ARK436L2SA_DIRK_6_3_4)) &&
         !((etable == ARK437L2SA_ERK_7_3_4) && (itable == ARK437L2SA_DIRK_7_3_4)) &&
         !((etable == ARK548L2SA_ERK_8_4_5) && (itable == ARK548L2SA_DIRK_8_4_5)) &&
         !((etable == ARK548L2SAb_ERK_8_4_5) && (itable == ARK548L2SAb_DIRK_8_4_5)) ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "ARKStepSetTableNum",
                      "Incompatible Butcher tables for ARK method");
      return(ARK_ILL_INPUT);
    }

    /* fill in tables based on arguments */
    step_mem->Bi = ARKodeButcherTable_LoadDIRK(itable);
    step_mem->Be = ARKodeButcherTable_LoadERK(etable);
    if (step_mem->Bi == NULL) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTableNum",
                      "Illegal IRK table number");
      return(ARK_ILL_INPUT);
    }
    if (step_mem->Be == NULL) {
      arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                      "ARKStepSetTableNum",
                      "Illegal ERK table number");
      return(ARK_ILL_INPUT);
    }
    step_mem->stages = step_mem->Bi->stages;
    step_mem->q = step_mem->Bi->q;
    step_mem->p = step_mem->Bi->p;

    /* set method as ImEx */
    if (ARKStepSetImEx(arkode_mem) != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "ARKStepSetTableNum", MSG_ARK_MISSING_F);
      return(ARK_ILL_INPUT);
    }

  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetNonlinCRDown:

  Specifies the user-provided nonlinear convergence constant
  crdown.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetNonlinCRDown(void *arkode_mem, realtype crdown)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetNonlinCRDown",
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
  ARKStepSetNonlinRDiv:

  Specifies the user-provided nonlinear convergence constant
  rdiv.  Legal values are strictly positive; illegal values
  imply a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetNonlinRDiv(void *arkode_mem, realtype rdiv)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetNonlinRDiv",
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
  ARKStepSetDeltaGammaMax:

  Specifies the user-provided linear setup decision constant
  dgmax.  Legal values are strictly positive; illegal values imply
  a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetDeltaGammaMax(void *arkode_mem, realtype dgmax)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetDeltaGammaMax",
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
  ARKStepSetLSetupFrequency:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the frequency for calling lsetup;
  negative values imply recomputation of lsetup at each nonlinear
  solve; a zero value implies a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetLSetupFrequency(void *arkode_mem, int msbp)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetLSetupFrequency",
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

/* Deprecated */
int ARKStepSetMaxStepsBetweenLSet(void *arkode_mem, int msbp)
{
  return(ARKStepSetLSetupFrequency(arkode_mem, msbp));
}


/*---------------------------------------------------------------
  ARKStepSetPredictorMethod:

  Specifies the method to use for predicting implicit solutions.
  Non-default choices are {1,2,3,4}, all others will use default
  (trivial) predictor.
  ---------------------------------------------------------------*/
int ARKStepSetPredictorMethod(void *arkode_mem, int pred_method)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetPredictorMethod",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* return error if pred_method==5 and a non-NULL stage predictor function
     has been supplied */
  if ((pred_method == 5) && (step_mem->stage_predict != NULL)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep", "ARKStepSetPredictorMethod",
                    "predictor 5 cannot be combined with user-supplied stage predictor");
    return(ARK_ILL_INPUT);
  }

  /* set parameter */
  step_mem->predictor = pred_method;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxNonlinIters:

  Specifies the maximum number of nonlinear iterations during
  one solve.  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int ARKStepSetMaxNonlinIters(void *arkode_mem, int maxcor)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetMaxNonlinIters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return error message if no NLS module is present */
  if (step_mem->NLS == NULL) {
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::ARKStep",
                    "ARKStepSetMaxNonlinIters",
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
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::ARKStep",
                    "ARKStepSetMaxNonlinIters",
                    "Error setting maxcor in SUNNonlinearSolver object");
    return(ARK_NLS_OP_ERR);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetNonlinConvCoef:

  Specifies the coefficient in the nonlinear solver convergence
  test.  A non-positive input implies a reset to the default value.
  ---------------------------------------------------------------*/
int ARKStepSetNonlinConvCoef(void *arkode_mem, realtype nlscoef)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetNonlinConvCoef",
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
  ARKStepSetStagePredictFn:  Specifies a user-provided step
  predictor function having type ARKStagePredictFn.  A
  NULL input function disables calls to this routine.
  ---------------------------------------------------------------*/
int ARKStepSetStagePredictFn(void *arkode_mem,
                             ARKStagePredictFn PredictStage)
{
  ARKodeMem        ark_mem;
  ARKodeARKStepMem step_mem;
  int              retval;

  /* access ARKodeARKStepMem structure and set function pointer */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetStagePredictFn",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* override predictor method 5 if non-NULL PredictStage is supplied */
  if ((step_mem->predictor == 5) && (PredictStage != NULL)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetStagePredictFn",
                    "User-supplied predictor is incompatible with predictor method 5");
    return(ARK_ILL_INPUT);
  }

  step_mem->stage_predict = PredictStage;
  return(ARK_SUCCESS);
}


/*===============================================================
  ARKStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepGetCurrentGamma: Returns the current value of gamma
  ---------------------------------------------------------------*/
int ARKStepGetCurrentGamma(void *arkode_mem, realtype *gamma)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  retval = arkStep_AccessStepMem(arkode_mem, NULL, &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);
  *gamma = step_mem->gamma;
  return(retval);
}


/*---------------------------------------------------------------
  ARKStepGetNumRhsEvals:

  Returns the current number of calls to fe and fi
  ---------------------------------------------------------------*/
int ARKStepGetNumRhsEvals(void *arkode_mem, long int *fe_evals,
                          long int *fi_evals)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetNumRhsEvals",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get values from step_mem */
  *fe_evals = step_mem->nfe;
  *fi_evals = step_mem->nfi;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumLinSolvSetups:

  Returns the current number of calls to the lsetup routine
  ---------------------------------------------------------------*/
int ARKStepGetNumLinSolvSetups(void *arkode_mem, long int *nlinsetups)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetNumLinSolvSetups",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get value from step_mem */
  *nlinsetups = step_mem->nsetups;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetCurrentButcherTables:

  Sets pointers to the explicit and implicit Butcher tables
  currently in use.
  ---------------------------------------------------------------*/
int ARKStepGetCurrentButcherTables(void *arkode_mem,
                                   ARKodeButcherTable *Bi,
                                   ARKodeButcherTable *Be)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetCurrentButcherTables",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get tables from step_mem */
  *Bi = step_mem->Bi;
  *Be = step_mem->Be;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetEstLocalErrors: (updated to the correct vector, but
  need to verify that it is unchanged between filling the
  estimated error and the end of the time step)

  Returns an estimate of the local error
  ---------------------------------------------------------------*/
int ARKStepGetEstLocalErrors(void *arkode_mem, N_Vector ele)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetEstLocalErrors",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* copy vector to output */
  N_VScale(ONE, ark_mem->tempv1, ele);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int ARKStepGetTimestepperStats(void *arkode_mem, long int *expsteps,
                               long int *accsteps, long int *step_attempts,
                               long int *fe_evals, long int *fi_evals,
                               long int *nlinsetups, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetTimestepperStats",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set expsteps and accsteps from adaptivity structure */
  *expsteps = ark_mem->hadapt_mem->nst_exp;
  *accsteps = ark_mem->hadapt_mem->nst_acc;

  /* set remaining outputs */
  *step_attempts = ark_mem->nst_attempts;
  *fe_evals      = step_mem->nfe;
  *fi_evals      = step_mem->nfi;
  *nlinsetups    = step_mem->nsetups;
  *netfails      = ark_mem->netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumNonlinSolvIters:

  Returns the current number of nonlinear solver iterations
  ---------------------------------------------------------------*/
int ARKStepGetNumNonlinSolvIters(void *arkode_mem, long int *nniters)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetNumNonlinSolvIters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  *nniters = step_mem->nls_iters;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumNonlinSolvConvFails:

  Returns the current number of nonlinear solver convergence fails
  ---------------------------------------------------------------*/
int ARKStepGetNumNonlinSolvConvFails(void *arkode_mem, long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetNumNonlinSolvConvFails",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set output from step_mem */
  *nncfails = ark_mem->ncfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNonlinSolvStats:

  Returns nonlinear solver statistics
  ---------------------------------------------------------------*/
int ARKStepGetNonlinSolvStats(void *arkode_mem, long int *nniters,
                              long int *nncfails)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetNonlinSolvStats",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  *nniters  = step_mem->nls_iters;
  *nncfails = ark_mem->ncfn;

  return(ARK_SUCCESS);
}


/*===============================================================
  ARKStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int ARKStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int flag, retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepWriteParameters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* output ARKode infrastructure parameters first */
  flag = arkWriteParameters(ark_mem, fp);
  if (flag != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepWriteParameters",
                    "Error writing ARKode infrastructure parameters");
    return(flag);
  }

  /* print integrator parameters to file */
  fprintf(fp, "ARKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n",step_mem->q);
  if (step_mem->linear) {
    fprintf(fp, "  Linear implicit problem");
    if (step_mem->linear_timedep) {
      fprintf(fp, " (time-dependent Jacobian)\n");
    } else {
      fprintf(fp, " (time-independent Jacobian)\n");
    }
  }
  if (step_mem->explicit && step_mem->implicit) {
    fprintf(fp, "  ImEx integrator\n");
  } else if (step_mem->implicit) {
    fprintf(fp, "  Implicit integrator\n");
  } else {
    fprintf(fp, "  Explicit integrator\n");
  }

  if (step_mem->implicit) {
    fprintf(fp, "  Implicit predictor method = %i\n",step_mem->predictor);
    fprintf(fp, "  Implicit solver tolerance coefficient = %"RSYM"\n",step_mem->nlscoef);
    fprintf(fp, "  Maximum number of nonlinear corrections = %i\n",step_mem->maxcor);
    fprintf(fp, "  Nonlinear convergence rate constant = %"RSYM"\n",step_mem->crdown);
    fprintf(fp, "  Nonlinear divergence tolerance = %"RSYM"\n",step_mem->rdiv);
    fprintf(fp, "  Gamma factor LSetup tolerance = %"RSYM"\n",step_mem->dgmax);
    fprintf(fp, "  Number of steps between LSetup calls = %i\n",step_mem->msbp);
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepWriteButcher:

  Outputs Butcher tables to the provided file pointer.
  ---------------------------------------------------------------*/
int ARKStepWriteButcher(void *arkode_mem, FILE *fp)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepWriteButcher",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* check that Butcher table is non-NULL (otherwise report error) */
  if ((step_mem->Be == NULL) && (step_mem->Bi == NULL)) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepWriteButcher", "Butcher table memory is NULL");
    return(ARK_MEM_NULL);
  }

  /* print Butcher tables to file */
  fprintf(fp, "\nARKStep Butcher tables (stages = %i):\n", step_mem->stages);
  if (step_mem->explicit && (step_mem->Be != NULL)) {
    fprintf(fp, "  Explicit Butcher table:\n");
    ARKodeButcherTable_Write(step_mem->Be, fp);
  }
  fprintf(fp, "\n");
  if (step_mem->implicit && (step_mem->Bi != NULL)) {
    fprintf(fp, "  Implicit Butcher table:\n");
    ARKodeButcherTable_Write(step_mem->Bi, fp);
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
