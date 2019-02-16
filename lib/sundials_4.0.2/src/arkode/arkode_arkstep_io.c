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
  utility routines)
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepSetDenseOrder: Specifies the polynomial order for dense
  output.  Positive values are sent to the interpolation module;
  negative values imply to use the default.
  ---------------------------------------------------------------*/
int ARKStepSetDenseOrder(void *arkode_mem, int dord)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetDenseOrder", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetDenseOrder(ark_mem, dord));
}

/*---------------------------------------------------------------
  ARKStepSetErrHandlerFn: Specifies the error handler function
  ---------------------------------------------------------------*/
int ARKStepSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                           void *eh_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetErrHandlerFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetErrHandlerFn(ark_mem, ehfun, eh_data));
}

/*---------------------------------------------------------------
  ARKStepSetErrFile: Specifies the FILE pointer for output (NULL
  means no messages)
  ---------------------------------------------------------------*/
int ARKStepSetErrFile(void *arkode_mem, FILE *errfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetErrFile", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetErrFile(ark_mem, errfp));
}

/*---------------------------------------------------------------
  ARKStepSetUserData: Specifies the user data pointer for f
  ---------------------------------------------------------------*/
int ARKStepSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetUserData", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetUserData(ark_mem, user_data));
}

/*---------------------------------------------------------------
  ARKStepSetDiagnostics: Specifies to enable solver diagnostics,
  and specifies the FILE pointer for output (diagfp==NULL
  disables output)
  ---------------------------------------------------------------*/
int ARKStepSetDiagnostics(void *arkode_mem, FILE *diagfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetDiagnostics", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetDiagnostics(ark_mem, diagfp));
}

/*---------------------------------------------------------------
  ARKStepSetMaxNumSteps: Specifies the maximum number of
  integration steps
  ---------------------------------------------------------------*/
int ARKStepSetMaxNumSteps(void *arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxNumSteps(ark_mem, mxsteps));
}

/*---------------------------------------------------------------
  ARKStepSetMaxHnilWarns: Specifies the maximum number of warnings
  for small h
  ---------------------------------------------------------------*/
int ARKStepSetMaxHnilWarns(void *arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxHnilWarns", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxHnilWarns(ark_mem, mxhnil));
}

/*---------------------------------------------------------------
  ARKStepSetInitStep: Specifies the initial step size
  ---------------------------------------------------------------*/
int ARKStepSetInitStep(void *arkode_mem, realtype hin)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetInitStep(ark_mem, hin));
}

/*---------------------------------------------------------------
  ARKStepSetMinStep: Specifies the minimum step size
  ---------------------------------------------------------------*/
int ARKStepSetMinStep(void *arkode_mem, realtype hmin)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMinStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMinStep(ark_mem, hmin));
}

/*---------------------------------------------------------------
  ARKStepSetMaxStep: Specifies the maximum step size
  ---------------------------------------------------------------*/
int ARKStepSetMaxStep(void *arkode_mem, realtype hmax)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxStep(ark_mem, hmax));
}

/*---------------------------------------------------------------
  ARKStepSetStopTime: Specifies the time beyond which the
  integration is not to proceed.
  ---------------------------------------------------------------*/
int ARKStepSetStopTime(void *arkode_mem, realtype tstop)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetStopTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetStopTime(ark_mem, tstop));
}

/*---------------------------------------------------------------
  ARKStepSetFixedStep: Specifies to use a fixed time step size
  instead of performing any form of temporal adaptivity.  ARKStep
  will use this step size for all steps (unless tstop is set, in
  which case it may need to modify that last step approaching
  tstop.  If any solver failure occurs in the timestepping
  module, ARKStep will typically immediately return with an error
  message indicating that the selected step size cannot be used.

  Any nonzero argument will result in the use of that fixed step
  size; an argument of 0 will re-enable temporal adaptivity.
  ---------------------------------------------------------------*/
int ARKStepSetFixedStep(void *arkode_mem, realtype hfixed)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetFixedStep",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* allocate or free adaptivity memory as needed */
  if (hfixed != ZERO) {
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
    }
  } else if (step_mem->hadapt_mem == NULL) {
    step_mem->hadapt_mem = arkAdaptInit();
    if (step_mem->hadapt_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                      "ARKStepSetFixedStep",
                      "Allocation of Step Adaptivity Structure Failed");
      return(ARK_MEM_FAIL);
    }
  }

  return(arkSetFixedStep(ark_mem, hfixed));
}

/*---------------------------------------------------------------
  ARKStepSetRootDirection: Specifies the direction of zero-crossings
  to be monitored.  The default is to monitor both crossings.
  ---------------------------------------------------------------*/
int ARKStepSetRootDirection(void *arkode_mem, int *rootdir)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetRootDirection", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetRootDirection(ark_mem, rootdir));
}

/*---------------------------------------------------------------
  ARKStepSetNoInactiveRootWarn:  Disables issuing a warning if
  some root function appears to be identically zero at the
  beginning of the integration
  ---------------------------------------------------------------*/
int ARKStepSetNoInactiveRootWarn(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetNoInactiveRootWarn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetNoInactiveRootWarn(ark_mem));
}

/*---------------------------------------------------------------
  ARKStepSetPostprocessStepFn:  Specifies a user-provided step
  postprocessing function having type ARKPostProcessStepFn.  A
  NULL input function disables step postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.
  ---------------------------------------------------------------*/
int ARKStepSetPostprocessStepFn(void *arkode_mem,
                                ARKPostProcessStepFn ProcessStep)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetPostprocessStepFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetPostprocessStepFn(ark_mem, ProcessStep));
}

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
int ARKStepSetMaxStepsBetweenJac(void *arkode_mem, long int msbj) {
  return(arkLSSetMaxStepsBetweenJac(arkode_mem, msbj)); }
int ARKStepSetEpsLin(void *arkode_mem, realtype eplifac) {
  return(arkLSSetEpsLin(arkode_mem, eplifac)); }
int ARKStepSetMassEpsLin(void *arkode_mem, realtype eplifac) {
  return(arkLSSetMassEpsLin(arkode_mem, eplifac)); }
int ARKStepSetPreconditioner(void *arkode_mem, ARKLsPrecSetupFn psetup,
                             ARKLsPrecSolveFn psolve) {
  return(arkLSSetPreconditioner(arkode_mem, psetup, psolve)); }
int ARKStepSetMassPreconditioner(void *arkode_mem, ARKLsMassPrecSetupFn psetup,
                                 ARKLsMassPrecSolveFn psolve) {
  return(arkLSSetMassPreconditioner(arkode_mem, psetup, psolve)); }
int ARKStepSetJacTimes(void *arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                       ARKLsJacTimesVecFn jtimes) {
  return(arkLSSetJacTimes(arkode_mem, jtsetup, jtimes)); }
int ARKStepSetMassTimes(void *arkode_mem, ARKLsMassTimesSetupFn msetup,
                        ARKLsMassTimesVecFn mtimes, void *mtimes_data) {
  return(arkLSSetMassTimes(arkode_mem, msetup, mtimes, mtimes_data)); }



/*===============================================================
  ARKStep Optional output functions (wrappers for generic ARKode
  utility routines)
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepGetNumSteps:  Returns the current number of integration
  steps
  ---------------------------------------------------------------*/
int ARKStepGetNumSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetNumSteps(ark_mem, nsteps));
}

/*---------------------------------------------------------------
  ARKStepGetActualInitStep: Returns the step size used on the
  first step
  ---------------------------------------------------------------*/
int ARKStepGetActualInitStep(void *arkode_mem, realtype *hinused)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetActualInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetActualInitStep(ark_mem, hinused));
}

/*---------------------------------------------------------------
  ARKStepGetLastStep: Returns the step size used on the last
  successful step
  ---------------------------------------------------------------*/
int ARKStepGetLastStep(void *arkode_mem, realtype *hlast)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetLastStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetLastStep(ark_mem, hlast));
}

/*---------------------------------------------------------------
  ARKStepGetCurrentStep: Returns the step size to be attempted on
  the next step
  ---------------------------------------------------------------*/
int ARKStepGetCurrentStep(void *arkode_mem, realtype *hcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetCurrentStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetCurrentStep(ark_mem, hcur));
}

/*---------------------------------------------------------------
  ARKStepGetCurrentTime: Returns the current value of the
  independent variable
  ---------------------------------------------------------------*/
int ARKStepGetCurrentTime(void *arkode_mem, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetCurrentTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetCurrentTime(ark_mem, tcur));
}

/*---------------------------------------------------------------
  ARKStepGetTolScaleFactor: Returns a suggested factor for scaling
  tolerances
  ---------------------------------------------------------------*/
int ARKStepGetTolScaleFactor(void *arkode_mem, realtype *tolsfact)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetTolScaleFactor", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetTolScaleFactor(ark_mem, tolsfact));
}

/*---------------------------------------------------------------
  ARKStepGetErrWeights: This routine returns the current error
  weight vector.
  ---------------------------------------------------------------*/
int ARKStepGetErrWeights(void *arkode_mem, N_Vector eweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetErrWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetErrWeights(ark_mem, eweight));
}

/*---------------------------------------------------------------
  ARKStepGetResWeights: This routine returns the current residual
  weight vector.
  ---------------------------------------------------------------*/
int ARKStepGetResWeights(void *arkode_mem, N_Vector rweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetResWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetResWeights(ark_mem, rweight));
}

/*---------------------------------------------------------------
  ARKStepGetWorkSpace: Returns integrator work space requirements
  ---------------------------------------------------------------*/
int ARKStepGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetWorkSpace", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetWorkSpace(ark_mem, lenrw, leniw));
}

/*---------------------------------------------------------------
  ARKStepGetNumGEvals: Returns the current number of calls to g
  (for rootfinding)
  ---------------------------------------------------------------*/
int ARKStepGetNumGEvals(void *arkode_mem, long int *ngevals)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetNumGEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetNumGEvals(ark_mem, ngevals));
}

/*---------------------------------------------------------------
  ARKStepGetRootInfo: Returns pointer to array rootsfound showing
  roots found
  ---------------------------------------------------------------*/
int ARKStepGetRootInfo(void *arkode_mem, int *rootsfound)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetRootInfo", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetRootInfo(ark_mem, rootsfound));
}

/*---------------------------------------------------------------
  ARKStepGetStepStats: Returns step statistics
  ---------------------------------------------------------------*/
int ARKStepGetStepStats(void *arkode_mem, long int *nsteps,
                        realtype *hinused, realtype *hlast,
                        realtype *hcur, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetStepStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetStepStats(ark_mem, nsteps, hinused, hlast, hcur, tcur));
}

/*---------------------------------------------------------------
  ARKStepGetReturnFlagName: translates from return flags IDs to
  names
  ---------------------------------------------------------------*/
char *ARKStepGetReturnFlagName(long int flag)
{ return(arkGetReturnFlagName(flag)); }

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
int ARKStepGetLastMassFlag(void *arkode_mem, long int *flag) {
  return(arkLSGetLastMassFlag(arkode_mem, flag)); }

char *ARKStepGetLinReturnFlagName(long int flag) {
  return(arkLSGetReturnFlagName(flag)); }



/*===============================================================
  ARKStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepSetDefaults:

  Resets all ARKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.  Also leaves alone any data
  structures/options related to the ARKode infrastructure itself
  (e.g. root-finding).
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
  step_mem->hadapt_pq        = SUNFALSE;       /* use embedding order */
  step_mem->predictor        = 0;              /* trivial predictor */
  step_mem->linear           = SUNFALSE;       /* nonlinear problem */
  step_mem->linear_timedep   = SUNTRUE;        /* dfi/dy depends on t */
  step_mem->explicit         = SUNTRUE;        /* fe(t,y) will be used */
  step_mem->implicit         = SUNTRUE;        /* fi(t,y) will be used */
  if (step_mem->hadapt_mem != NULL) {
    step_mem->hadapt_mem->etamx1      = ETAMX1;     /* max change on first step */
    step_mem->hadapt_mem->etamxf      = ETAMXF;     /* max change on error-failed step */
    step_mem->hadapt_mem->small_nef   = SMALL_NEF;  /* num error fails before ETAMXF enforced */
    step_mem->hadapt_mem->etacf       = ETACF;      /* max change on convergence failure */
    step_mem->hadapt_mem->HAdapt      = NULL;       /* step adaptivity fn */
    step_mem->hadapt_mem->HAdapt_data = NULL;       /* step adaptivity data */
    step_mem->hadapt_mem->imethod     = 0;          /* PID controller */
    step_mem->hadapt_mem->cfl         = CFLFAC;     /* explicit stability factor */
    step_mem->hadapt_mem->safety      = SAFETY;     /* step adaptivity safety factor  */
    step_mem->hadapt_mem->bias        = BIAS;       /* step adaptivity error bias */
    step_mem->hadapt_mem->growth      = GROWTH;     /* step adaptivity growth factor */
    step_mem->hadapt_mem->lbound      = HFIXED_LB;  /* step adaptivity no-change lower bound */
    step_mem->hadapt_mem->ubound      = HFIXED_UB;  /* step adaptivity no-change upper bound */
    step_mem->hadapt_mem->k1          = AD0_K1;     /* step adaptivity parameter */
    step_mem->hadapt_mem->k2          = AD0_K2;     /* step adaptivity parameter */
    step_mem->hadapt_mem->k3          = AD0_K3;     /* step adaptivity parameter */
  }
  step_mem->maxcor           = MAXCOR;         /* max nonlinear iters/stage */
  step_mem->maxnef           = MAXNEF;         /* max error test fails */
  step_mem->maxncf           = MAXNCF;         /* max convergence fails */
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
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetOptimalParams",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* Choose values based on method, order */

  /*    explicit */
  if (step_mem->explicit && !step_mem->implicit) {
    hadapt_mem->imethod = 1;
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
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = SAFETY;
      hadapt_mem->bias      = BIAS;
      hadapt_mem->growth    = GROWTH;
      hadapt_mem->etamxf    = ETAMXF;
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef  = RCONST(0.001);
      step_mem->maxcor   = 5;
      step_mem->crdown   = CRDOWN;
      step_mem->rdiv     = RDIV;
      step_mem->dgmax    = DGMAX;
      step_mem->msbp     = MSBP;
      break;
    case 3:
      hadapt_mem->imethod   = 2;
      hadapt_mem->safety    = RCONST(0.957);
      hadapt_mem->bias      = RCONST(1.9);
      hadapt_mem->growth    = RCONST(17.6);
      hadapt_mem->etamxf    = RCONST(0.45);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef  = RCONST(0.22);
      step_mem->crdown   = RCONST(0.17);
      step_mem->rdiv     = RCONST(2.3);
      step_mem->dgmax    = RCONST(0.19);
      step_mem->msbp     = 60;
      break;
    case 4:
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = RCONST(0.988);
      hadapt_mem->bias      = RCONST(1.2);
      hadapt_mem->growth    = RCONST(31.5);
      hadapt_mem->k1        = RCONST(0.535);
      hadapt_mem->k2        = RCONST(0.209);
      hadapt_mem->k3        = RCONST(0.148);
      hadapt_mem->etamxf    = RCONST(0.33);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef  = RCONST(0.24);
      step_mem->crdown   = RCONST(0.26);
      step_mem->rdiv     = RCONST(2.3);
      step_mem->dgmax    = RCONST(0.16);
      step_mem->msbp     = 31;
      break;
    case 5:
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = RCONST(0.937);
      hadapt_mem->bias      = RCONST(3.3);
      hadapt_mem->growth    = RCONST(22.0);
      hadapt_mem->k1        = RCONST(0.56);
      hadapt_mem->k2        = RCONST(0.338);
      hadapt_mem->k3        = RCONST(0.14);
      hadapt_mem->etamxf    = RCONST(0.44);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef  = RCONST(0.25);
      step_mem->crdown   = RCONST(0.4);
      step_mem->rdiv     = RCONST(2.3);
      step_mem->dgmax    = RCONST(0.32);
      step_mem->msbp     = 31;
      break;
    }

    /*    imex */
  } else {
    switch (step_mem->q) {
    case 3:
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = RCONST(0.965);
      hadapt_mem->bias      = RCONST(1.42);
      hadapt_mem->growth    = RCONST(28.7);
      hadapt_mem->k1        = RCONST(0.54);
      hadapt_mem->k2        = RCONST(0.36);
      hadapt_mem->k3        = RCONST(0.14);
      hadapt_mem->etamxf    = RCONST(0.46);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef  = RCONST(0.22);
      step_mem->crdown   = RCONST(0.17);
      step_mem->rdiv     = RCONST(2.3);
      step_mem->dgmax    = RCONST(0.19);
      step_mem->msbp     = 60;
      break;
    case 4:
      hadapt_mem->imethod   = 0;
      hadapt_mem->safety    = RCONST(0.97);
      hadapt_mem->bias      = RCONST(1.35);
      hadapt_mem->growth    = RCONST(25.0);
      hadapt_mem->k1        = RCONST(0.543);
      hadapt_mem->k2        = RCONST(0.297);
      hadapt_mem->k3        = RCONST(0.14);
      hadapt_mem->etamxf    = RCONST(0.47);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef  = RCONST(0.24);
      step_mem->crdown   = RCONST(0.26);
      step_mem->rdiv     = RCONST(2.3);
      step_mem->dgmax    = RCONST(0.16);
      step_mem->msbp     = 31;
      break;
    case 5:
      hadapt_mem->imethod   = 1;
      hadapt_mem->safety    = RCONST(0.993);
      hadapt_mem->bias      = RCONST(1.15);
      hadapt_mem->growth    = RCONST(28.5);
      hadapt_mem->k1        = RCONST(0.8);
      hadapt_mem->k2        = RCONST(0.35);
      hadapt_mem->etamxf    = RCONST(0.3);
      hadapt_mem->small_nef = SMALL_NEF;
      hadapt_mem->etacf     = ETACF;
      step_mem->nlscoef  = RCONST(0.25);
      step_mem->crdown   = RCONST(0.4);
      step_mem->rdiv     = RCONST(2.3);
      step_mem->dgmax    = RCONST(0.32);
      step_mem->msbp     = 31;
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
  ARKodeButcherTable_Free(step_mem->Be);  step_mem->Be = NULL;
  ARKodeButcherTable_Free(step_mem->Bi);  step_mem->Bi = NULL;

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
  ARKodeButcherTable_Free(step_mem->Be);  step_mem->Be = NULL;
  ARKodeButcherTable_Free(step_mem->Bi);  step_mem->Bi = NULL;

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

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetTableNum",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  ARKodeButcherTable_Free(step_mem->Be);  step_mem->Be = NULL;
  ARKodeButcherTable_Free(step_mem->Bi);  step_mem->Bi = NULL;


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
         !((etable == ARK548L2SA_ERK_8_4_5) && (itable == ARK548L2SA_DIRK_8_4_5)) ) {
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
  ARKStepSetCFLFraction:

  Specifies the safety factor to use on the maximum explicitly-
  stable step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ARKStepSetCFLFraction(void *arkode_mem, realtype cfl_frac)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetCFLFraction",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetCFLFraction",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if (cfl_frac >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetCFLFraction", "Illegal CFL fraction");
    return(ARK_ILL_INPUT);
  }

  /* set positive-valued parameters, otherwise set default */
  if (cfl_frac <= ZERO) {
    hadapt_mem->cfl = CFLFAC;
  } else {
    hadapt_mem->cfl = cfl_frac;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetSafetyFactor:

  Specifies the safety factor to use on the error-based predicted
  time step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int ARKStepSetSafetyFactor(void *arkode_mem, realtype safety)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetSafetyFactor",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetSafetyFactoy",MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if (safety >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetSafetyFactor", "Illegal safety factor");
    return(ARK_ILL_INPUT);
  }

  /* set positive-valued parameters, otherwise set default */
  if (safety <= ZERO) {
    hadapt_mem->safety = SAFETY;
  } else {
    hadapt_mem->safety = safety;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetErrorBias:

  Specifies the error bias to use when performing adaptive-step
  error control.  Allowable values must be >= 1.0.  Any illegal
  value implies a reset to the default value.
  ---------------------------------------------------------------*/
int ARKStepSetErrorBias(void *arkode_mem, realtype bias)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetErrorBias",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetErrorBias", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (bias < 1.0) {
    hadapt_mem->bias = BIAS;
  } else {
    hadapt_mem->bias = bias;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxGrowth:

  Specifies the maximum step size growth factor to be allowed
  between successive integration steps.  Note: the first step uses
  a separate maximum growth factor.  Allowable values must be
  > 1.0.  Any illegal value implies a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetMaxGrowth(void *arkode_mem, realtype mx_growth)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetMaxGrowth",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* set allowed value, otherwise set default */
  if (mx_growth == ZERO) {
    hadapt_mem->growth = GROWTH;
  } else {
    hadapt_mem->growth = mx_growth;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetFixedStepBounds:

  Specifies the step size growth interval within which the step
  size will remain unchanged.  Allowable values must enclose the
  value 1.0.  Any illegal interval implies a reset to the default.
  ---------------------------------------------------------------*/
int ARKStepSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetFixedStepBounds",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetFixedStepBounds", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* set allowable interval, otherwise set defaults */
  if ((lb <= 1.0) && (ub >= 1.0)) {
    hadapt_mem->lbound = lb;
    hadapt_mem->ubound = ub;
  } else {
    hadapt_mem->lbound = HFIXED_LB;
    hadapt_mem->ubound = HFIXED_UB;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetAdaptivityMethod:

  Specifies the built-in time step adaptivity algorithm (and
  optionally, its associated parameters) to use.  All parameters
  will be checked for validity when used by the solver.
  ---------------------------------------------------------------*/
int ARKStepSetAdaptivityMethod(void *arkode_mem, int imethod,
                               int idefault, int pq,
                               realtype *adapt_params)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetAdaptivityMethod",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetAdaptivityMethod", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if ((imethod > 5) || (imethod < 0)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  /* set adaptivity method */
  hadapt_mem->imethod = imethod;

  /* set flag whether to use p or q */
  step_mem->hadapt_pq = (pq != 0);

  /* set method parameters */
  if (idefault == 1) {
    switch (hadapt_mem->imethod) {
    case (0):
      hadapt_mem->k1 = AD0_K1;
      hadapt_mem->k2 = AD0_K2;
      hadapt_mem->k3 = AD0_K3; break;
    case (1):
      hadapt_mem->k1 = AD1_K1;
      hadapt_mem->k2 = AD1_K2; break;
    case (2):
      hadapt_mem->k1 = AD2_K1; break;
    case (3):
      hadapt_mem->k1 = AD3_K1;
      hadapt_mem->k2 = AD3_K2; break;
    case (4):
      hadapt_mem->k1 = AD4_K1;
      hadapt_mem->k2 = AD4_K2; break;
    case (5):
      hadapt_mem->k1 = AD5_K1;
      hadapt_mem->k2 = AD5_K2;
      hadapt_mem->k3 = AD5_K3; break;
    }
  } else {
    hadapt_mem->k1 = adapt_params[0];
    hadapt_mem->k2 = adapt_params[1];
    hadapt_mem->k3 = adapt_params[2];
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetAdaptivityFn:

  Specifies the user-provided time step adaptivity function to use.
  ---------------------------------------------------------------*/
int ARKStepSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun,
                           void *h_data)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetAdaptivityFn",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetAdaptivityFn", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* NULL hfun sets default, otherwise set inputs */
  if (hfun == NULL) {
    hadapt_mem->HAdapt      = NULL;
    hadapt_mem->HAdapt_data = NULL;
    hadapt_mem->imethod     = 0;
  } else {
    hadapt_mem->HAdapt      = hfun;
    hadapt_mem->HAdapt_data = h_data;
    hadapt_mem->imethod     = -1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxFirstGrowth:

  Specifies the user-provided time step adaptivity constant
  etamx1.  Legal values are greater than 1.0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ARKStepSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetMaxFirstGrowth",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxFirstGrowth",MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* if argument legal set it, otherwise set default */
  if (etamx1 <= ONE) {
    hadapt_mem->etamx1 = ETAMX1;
  } else {
    hadapt_mem->etamx1 = etamx1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxEFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etamxf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ARKStepSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetMaxEFailGrowth",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxEFailGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* if argument legal set it, otherwise set default */
  if ((etamxf <= ZERO) || (etamxf > ONE)) {
    hadapt_mem->etamxf = ETAMXF;
  } else {
    hadapt_mem->etamxf = etamxf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetSmallNumEFails:

  Specifies the user-provided time step adaptivity constant
  small_nef.  Legal values are > 0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ARKStepSetSmallNumEFails(void *arkode_mem, int small_nef)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetSmallNumEFails",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetSmallNumEFails", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0) {
    hadapt_mem->small_nef = SMALL_NEF;
  } else {
    hadapt_mem->small_nef = small_nef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxCFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etacf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ARKStepSetMaxCFailGrowth(void *arkode_mem, realtype etacf)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetMaxCFailGrowth",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetMaxCFailGrowth", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* if argument legal set it, otherwise set default */
  if ((etacf <= ZERO) || (etacf > ONE)) {
    hadapt_mem->etacf = ETACF;
  } else {
    hadapt_mem->etacf = etacf;
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
  ARKStepSetMaxStepsBetweenLSet:

  Specifies the user-provided linear setup decision constant
  msbp.  Positive values give the number of time steps to wait
  before calling lsetup; negative values imply recomputation of
  lsetup at each nonlinear solve; a zero value implies a reset
  to the default.
  ---------------------------------------------------------------*/
int ARKStepSetMaxStepsBetweenLSet(void *arkode_mem, int msbp)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetMaxStepsBetweenLSet",
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

  /* set parameters */
  step_mem->predictor = pred_method;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetStabilityFn:

  Specifies the user-provided explicit time step stability
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int ARKStepSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab,
                          void *estab_data)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetStabilityFn",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* access ARKodeHAdaptMem structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSetStabilityFn", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* NULL argument sets default, otherwise set inputs */
  if (EStab == NULL) {
    hadapt_mem->expstab    = arkExpStab;
    hadapt_mem->estab_data = ark_mem;
  } else {
    hadapt_mem->expstab    = EStab;
    hadapt_mem->estab_data = estab_data;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ARKStepSetMaxErrTestFails(void *arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetMaxErrTestFails",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    step_mem->maxnef = MAXNEF;
  } else {
    step_mem->maxnef = maxnef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSetMaxConvFails:

  Specifies the maximum number of nonlinear convergence failures
  during one step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ARKStepSetMaxConvFails(void *arkode_mem, int maxncf)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetMaxConvFails",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* argument <= 0 sets default, otherwise set input */
  if (maxncf <= 0) {
    step_mem->maxncf = MAXNCF;
  } else {
    step_mem->maxncf = maxncf;
  }

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


/*===============================================================
  ARKStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepGetNumExpSteps:

  Returns the current number of stability-limited steps
  ---------------------------------------------------------------*/
int ARKStepGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetNumExpSteps",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if step adaptivity structure not allocated, just return 0 */
  if (step_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = step_mem->hadapt_mem->nst_exp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumAccSteps:

  Returns the current number of accuracy-limited steps
  ---------------------------------------------------------------*/
int ARKStepGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetNumAccSteps",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if step adaptivity structure not allocated, just return 0 */
  if (step_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = step_mem->hadapt_mem->nst_acc;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepGetNumStepAttempts:

  Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
int ARKStepGetNumStepAttempts(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetNumStepAttempts",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get value from step_mem */
  *nsteps = step_mem->nst_attempts;

  return(ARK_SUCCESS);
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
  ARKStepGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
int ARKStepGetNumErrTestFails(void *arkode_mem, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepGetNumErrTestFails",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* get value from step_mem */
  *netfails = step_mem->netf;

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

  /* if step adaptivity structure not allocated,
     just set expsteps and accsteps to 0 */
  if (step_mem->hadapt_mem == NULL) {
    *expsteps = 0;
    *accsteps = 0;
  } else {
    *expsteps = step_mem->hadapt_mem->nst_exp;
    *accsteps = step_mem->hadapt_mem->nst_acc;
  }

  /* set remaining outputs from step_mem */
  *step_attempts = step_mem->nst_attempts;
  *fe_evals      = step_mem->nfe;
  *fi_evals      = step_mem->nfi;
  *nlinsetups    = step_mem->nsetups;
  *netfails      = step_mem->netf;

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

  /* if a NLS object is present, set output from that; otherwise
     we took zero iterations */
  if (step_mem->NLS) {
    retval = SUNNonlinSolGetNumIters(step_mem->NLS, nniters);
    if (retval != SUN_NLS_SUCCESS) {
      arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::ARKStep",
                      "ARKStepGetNumNonlinSolvIters",
                      "Error retrieving nniters from SUNNonlinearSolver");
      return(ARK_NLS_OP_ERR);
    }
  } else {
    *nniters = 0;
  }

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
  *nncfails = step_mem->ncfn;

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

  /* set outputs from NLS module and step_mem structure (if present);
     otherwise there were zero iterations and no nonlinear failures */
  if (step_mem->NLS) {
    retval = SUNNonlinSolGetNumIters(step_mem->NLS, nniters);
    if (retval != SUN_NLS_SUCCESS) {
      arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode::ARKStep",
                      "ARKStepGetNonlinSolvStats",
                      "Error retrieving nniters from SUNNonlinearSolver");
      return(ARK_NLS_OP_ERR);
    }
    *nncfails = step_mem->ncfn;
  } else {
    *nniters = 0;
    *nncfails = 0;
  }

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
  if (step_mem->hadapt_mem != NULL) {
    fprintf(fp, "  Maximum step increase (first step) = %"RSYM"\n",
            step_mem->hadapt_mem->etamx1);
    fprintf(fp, "  Step reduction factor on multiple error fails = %"RSYM"\n",
            step_mem->hadapt_mem->etamxf);
    fprintf(fp, "  Minimum error fails before above factor is used = %i\n",
            step_mem->hadapt_mem->small_nef);
    fprintf(fp, "  Step reduction factor on nonlinear convergence failure = %"RSYM"\n",
            step_mem->hadapt_mem->etacf);
    if (step_mem->explicit)
      fprintf(fp, "  Explicit safety factor = %"RSYM"\n",
              step_mem->hadapt_mem->cfl);
    if (step_mem->hadapt_mem->HAdapt == NULL) {
      fprintf(fp, "  Time step adaptivity method %i\n", step_mem->hadapt_mem->imethod);
      fprintf(fp, "     Safety factor = %"RSYM"\n", step_mem->hadapt_mem->safety);
      fprintf(fp, "     Bias factor = %"RSYM"\n", step_mem->hadapt_mem->bias);
      fprintf(fp, "     Growth factor = %"RSYM"\n", step_mem->hadapt_mem->growth);
      fprintf(fp, "     Step growth lower bound = %"RSYM"\n", step_mem->hadapt_mem->lbound);
      fprintf(fp, "     Step growth upper bound = %"RSYM"\n", step_mem->hadapt_mem->ubound);
      fprintf(fp, "     k1 = %"RSYM"\n", step_mem->hadapt_mem->k1);
      fprintf(fp, "     k2 = %"RSYM"\n", step_mem->hadapt_mem->k2);
      fprintf(fp, "     k3 = %"RSYM"\n", step_mem->hadapt_mem->k3);
      if (step_mem->hadapt_mem->expstab == arkExpStab) {
        fprintf(fp, "  Default explicit stability function\n");
      } else {
        fprintf(fp, "  User provided explicit stability function\n");
      }
    } else {
      fprintf(fp, "  User provided time step adaptivity function\n");
    }
  }

  fprintf(fp, "  Maximum number of error test failures = %i\n",step_mem->maxnef);

  if (step_mem->implicit) {
    fprintf(fp, "  Maximum number of convergence test failures = %i\n",step_mem->maxncf);
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
