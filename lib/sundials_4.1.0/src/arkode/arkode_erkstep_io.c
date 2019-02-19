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
 * output functions for the ARKode ERKStep time stepper module.
 *
 * NOTE: many functions currently in arkode_io.c will move here,
 * with slightly different names.  The code transition will be
 * minimal, but the documentation changes will be significant.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_erkstep_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  ERKStep Optional input functions (wrappers for generic ARKode
  utility routines)
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepSetDenseOrder: Specifies the polynomial order for dense
  output.  Positive values are sent to the interpolation module;
  negative values imply to use the default.
  ---------------------------------------------------------------*/
int ERKStepSetDenseOrder(void *arkode_mem, int dord)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetDenseOrder", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetDenseOrder(ark_mem, dord));
}

/*---------------------------------------------------------------
  ERKStepSetErrHandlerFn: Specifies the error handler function
  ---------------------------------------------------------------*/
int ERKStepSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                           void *eh_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetErrHandlerFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetErrHandlerFn(ark_mem, ehfun, eh_data));
}

/*---------------------------------------------------------------
  ERKStepSetErrFile: Specifies the FILE pointer for output (NULL
  means no messages)
  ---------------------------------------------------------------*/
int ERKStepSetErrFile(void *arkode_mem, FILE *errfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetErrFile", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetErrFile(ark_mem, errfp));
}

/*---------------------------------------------------------------
  ERKStepSetUserData: Specifies the user data pointer for f
  ---------------------------------------------------------------*/
int ERKStepSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetUserData", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetUserData(ark_mem, user_data));
}

/*---------------------------------------------------------------
 ERKStepSetDiagnostics: Specifies to enable solver diagnostics,
 and specifies the FILE pointer for output (diagfp==NULL
 disables output)
---------------------------------------------------------------*/
int ERKStepSetDiagnostics(void *arkode_mem, FILE *diagfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetDiagnostics", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetDiagnostics(ark_mem, diagfp));
}

/*---------------------------------------------------------------
  ERKStepSetMaxNumSteps: Specifies the maximum number of
  integration steps
  ---------------------------------------------------------------*/
int ERKStepSetMaxNumSteps(void *arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetMaxNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxNumSteps(ark_mem, mxsteps));
}

/*---------------------------------------------------------------
 ERKStepSetMaxHnilWarns: Specifies the maximum number of warnings
 for small h
---------------------------------------------------------------*/
int ERKStepSetMaxHnilWarns(void *arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetMaxHnilWarns", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxHnilWarns(ark_mem, mxhnil));
}

/*---------------------------------------------------------------
  ERKStepSetInitStep: Specifies the initial step size
  ---------------------------------------------------------------*/
int ERKStepSetInitStep(void *arkode_mem, realtype hin)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetInitStep(ark_mem, hin));
}

/*---------------------------------------------------------------
  ERKStepSetMinStep: Specifies the minimum step size
  ---------------------------------------------------------------*/
int ERKStepSetMinStep(void *arkode_mem, realtype hmin)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetMinStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMinStep(ark_mem, hmin));
}

/*---------------------------------------------------------------
  ERKStepSetMaxStep: Specifies the maximum step size
  ---------------------------------------------------------------*/
int ERKStepSetMaxStep(void *arkode_mem, realtype hmax)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetMaxStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetMaxStep(ark_mem, hmax));
}

/*---------------------------------------------------------------
  ERKStepSetStopTime: Specifies the time beyond which the
  integration is not to proceed.
  ---------------------------------------------------------------*/
int ERKStepSetStopTime(void *arkode_mem, realtype tstop)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetStopTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetStopTime(ark_mem, tstop));
}

/*---------------------------------------------------------------
  ERKStepSetFixedStep: Specifies to use a fixed time step size
  instead of performing any form of temporal adaptivity.  ERKStep
  will use this step size for all steps (unless tstop is set, in
  which case it may need to modify that last step approaching
  tstop.  If any solver failure occurs in the timestepping
  module, ERKStep will typically immediately return with an error
  message indicating that the selected step size cannot be used.

  Any nonzero argument will result in the use of that fixed step
  size; an argument of 0 will re-enable temporal adaptivity.
  ---------------------------------------------------------------*/
int ERKStepSetFixedStep(void *arkode_mem, realtype hfixed)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetFixedStep",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* allocate or free adaptivity memory as needed */
  if (hfixed != ZERO) {
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
    }
  } else if (step_mem->hadapt_mem == NULL) {
    step_mem->hadapt_mem = arkAdaptInit();
    if (step_mem->hadapt_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ERKStep",
                      "ERKStepSetFixedStep",
                      "Allocation of Step Adaptivity Structure Failed");
      return(ARK_MEM_FAIL);
    }
  }

  return(arkSetFixedStep(ark_mem, hfixed));
}

/*---------------------------------------------------------------
  ERKStepSetRootDirection: Specifies the direction of zero-crossings
  to be monitored.  The default is to monitor both crossings.
  ---------------------------------------------------------------*/
int ERKStepSetRootDirection(void *arkode_mem, int *rootdir)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetRootDirection", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetRootDirection(ark_mem, rootdir));
}

/*---------------------------------------------------------------
  ERKStepSetNoInactiveRootWarn:  Disables issuing a warning if
  some root function appears to be identically zero at the
  beginning of the integration
  ---------------------------------------------------------------*/
int ERKStepSetNoInactiveRootWarn(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetNoInactiveRootWarn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetNoInactiveRootWarn(ark_mem));
}

/*---------------------------------------------------------------
  ERKStepSetPostprocessStepFn:  Specifies a user-provided step
  postprocessing function having type ARKPostProcessStepFn.  A
  NULL input function disables step postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.
  ---------------------------------------------------------------*/
int ERKStepSetPostprocessStepFn(void *arkode_mem,
                                ARKPostProcessStepFn ProcessStep)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetPostprocessStepFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetPostprocessStepFn(ark_mem, ProcessStep));
}



/*===============================================================
  ERKStep Optional output functions (wrappers for generic ARKode
  utility routines)
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepGetNumSteps:  Returns the current number of integration
  steps
  ---------------------------------------------------------------*/
int ERKStepGetNumSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetNumSteps(ark_mem, nsteps));
}

/*---------------------------------------------------------------
  ERKStepGetActualInitStep: Returns the step size used on the
  first step
  ---------------------------------------------------------------*/
int ERKStepGetActualInitStep(void *arkode_mem, realtype *hinused)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetActualInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetActualInitStep(ark_mem, hinused));
}

/*---------------------------------------------------------------
  ERKStepGetLastStep: Returns the step size used on the last
  successful step
  ---------------------------------------------------------------*/
int ERKStepGetLastStep(void *arkode_mem, realtype *hlast)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetLastStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetLastStep(ark_mem, hlast));
}

/*---------------------------------------------------------------
  ERKStepGetCurrentStep: Returns the step size to be attempted on
  the next step
  ---------------------------------------------------------------*/
int ERKStepGetCurrentStep(void *arkode_mem, realtype *hcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetCurrentStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetCurrentStep(ark_mem, hcur));
}

/*---------------------------------------------------------------
  ERKStepGetCurrentTime: Returns the current value of the
  independent variable
  ---------------------------------------------------------------*/
int ERKStepGetCurrentTime(void *arkode_mem, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetCurrentTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetCurrentTime(ark_mem, tcur));
}

/*---------------------------------------------------------------
  ERKStepGetTolScaleFactor: Returns a suggested factor for scaling
  tolerances
  ---------------------------------------------------------------*/
int ERKStepGetTolScaleFactor(void *arkode_mem, realtype *tolsfact)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetTolScaleFactor", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetTolScaleFactor(ark_mem, tolsfact));
}

/*---------------------------------------------------------------
  ERKStepGetErrWeights: This routine returns the current error
  weight vector.
  ---------------------------------------------------------------*/
int ERKStepGetErrWeights(void *arkode_mem, N_Vector eweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetErrWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetErrWeights(ark_mem, eweight));
}

/*---------------------------------------------------------------
  ERKStepGetWorkSpace: Returns integrator work space requirements
  ---------------------------------------------------------------*/
int ERKStepGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetWorkSpace", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetWorkSpace(ark_mem, lenrw, leniw));
}

/*---------------------------------------------------------------
  ERKStepGetNumGEvals: Returns the current number of calls to g
  (for rootfinding)
  ---------------------------------------------------------------*/
int ERKStepGetNumGEvals(void *arkode_mem, long int *ngevals)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetNumGEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetNumGEvals(ark_mem, ngevals));
}

/*---------------------------------------------------------------
  ERKStepGetRootInfo: Returns pointer to array rootsfound showing
  roots found
  ---------------------------------------------------------------*/
int ERKStepGetRootInfo(void *arkode_mem, int *rootsfound)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetRootInfo", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetRootInfo(ark_mem, rootsfound));
}

/*---------------------------------------------------------------
  ERKStepGetStepStats: Returns step statistics
  ---------------------------------------------------------------*/
int ERKStepGetStepStats(void *arkode_mem, long int *nsteps,
                        realtype *hinused, realtype *hlast,
                        realtype *hcur, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetStepStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetStepStats(ark_mem, nsteps, hinused, hlast, hcur, tcur));
}

/*---------------------------------------------------------------
  ERKStepGetReturnFlagName: translates from return flags IDs to
  names
  ---------------------------------------------------------------*/
char *ERKStepGetReturnFlagName(long int flag)
{ return(arkGetReturnFlagName(flag)); }



/*===============================================================
  ERKStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepSetDefaults:

  Resets all ERKStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
int ERKStepSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set default ARKode infrastructure parameters */
  retval = arkSetDefaults(arkode_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetDefaults",
                    "Error setting ARKode infrastructure defaults");
    return(retval);
  }

  /* Set default values for integrator optional inputs */
  step_mem->q         = Q_DEFAULT;      /* method order */
  step_mem->p         = 0;              /* embedding order */
  step_mem->hadapt_pq = SUNFALSE;       /* use embedding order */
  if (step_mem->hadapt_mem != NULL) {
    step_mem->hadapt_mem->etamx1      = ETAMX1;       /* max change on first step */
    step_mem->hadapt_mem->etamxf      = RCONST(0.3);  /* max change on error-failed step */
    step_mem->hadapt_mem->small_nef   = SMALL_NEF ;   /* num error fails before ETAMXF enforced */
    step_mem->hadapt_mem->etacf       = ETACF;        /* max change on convergence failure */
    step_mem->hadapt_mem->HAdapt      = NULL;         /* step adaptivity fn */
    step_mem->hadapt_mem->HAdapt_data = NULL;         /* step adaptivity data */
    step_mem->hadapt_mem->imethod     = 1;            /* PI controller */
    step_mem->hadapt_mem->cfl         = CFLFAC;       /* explicit stability factor */
    step_mem->hadapt_mem->safety      = RCONST(0.99); /* step adaptivity safety factor  */
    step_mem->hadapt_mem->bias        = RCONST(1.2);  /* step adaptivity error bias */
    step_mem->hadapt_mem->growth      = RCONST(25.0); /* step adaptivity growth factor */
    step_mem->hadapt_mem->lbound      = HFIXED_LB;    /* step adaptivity no-change lower bound */
    step_mem->hadapt_mem->ubound      = HFIXED_UB;    /* step adaptivity no-change upper bound */
    step_mem->hadapt_mem->k1          = RCONST(0.8);  /* step adaptivity parameter */
    step_mem->hadapt_mem->k2          = RCONST(0.31); /* step adaptivity parameter */
    step_mem->hadapt_mem->k3          = AD0_K3;       /* step adaptivity parameter */
  }
  step_mem->maxnef = MAXNEF;         /* max error test fails */
  step_mem->stages = 0;              /* no stages */
  step_mem->B      = NULL;           /* no Butcher table */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetOrder:

  Specifies the method order

  ** Note in documentation that this should not be called along
  with ERKStepSetTable or ERKStepSetTableNum.  This
  routine is used to specify a desired method order using
  default Butcher tables, whereas any user-supplied table will
  have their own order associated with them.
  ---------------------------------------------------------------*/
int ERKStepSetOrder(void *arkode_mem, int ord)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetOrder",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) {
    step_mem->q = Q_DEFAULT;
  } else {
    step_mem->q = ord;
  }

  /* clear Butcher tables, since user is requesting a change in method
     or a reset to defaults.  Tables will be set in ARKInitialSetup. */
  step_mem->stages = 0;
  step_mem->p = 0;
  ARKodeButcherTable_Free(step_mem->B);  step_mem->B = NULL;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetTable:

  Specifies to use a customized Butcher table for the explicit
  portion of the system.

  If d==NULL, then the method is automatically flagged as a
  fixed-step method; a user MUST also call either
  ERKStepSetFixedStep or ERKStepSetInitStep to set the desired
  time step size.
  ---------------------------------------------------------------*/
int ERKStepSetTable(void *arkode_mem, ARKodeButcherTable B)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetTable",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check for legal inputs */
  if (B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  ARKodeButcherTable_Free(step_mem->B); step_mem->B = NULL;

  /* set the relevant parameters */
  step_mem->stages = B->stages;
  step_mem->q = B->q;
  step_mem->p = B->p;

  /* copy the table into step memory */
  step_mem->B = ARKodeButcherTable_Copy(B);
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetTable", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetTableNum:

  Specifies to use a pre-existing Butcher table for the problem,
  based on the integer flag passed to ARKodeButcherTable_LoadERK()
  within the file arkode_butcher_erk.c.
  ---------------------------------------------------------------*/
int ERKStepSetTableNum(void *arkode_mem, int itable)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetTableNum",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check that argument specifies an explicit table */
  if (itable<MIN_ERK_NUM || itable>MAX_ERK_NUM) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetTableNum",
                    "Illegal ERK table number");
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  ARKodeButcherTable_Free(step_mem->B);  step_mem->B = NULL;

  /* fill in table based on argument */
  step_mem->B = ARKodeButcherTable_LoadERK(itable);
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetTableNum",
                    "Error setting table with that index");
    return(ARK_ILL_INPUT);
  }
  step_mem->stages = step_mem->B->stages;
  step_mem->q = step_mem->B->q;
  step_mem->p = step_mem->B->p;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepSetCFLFraction:

  Specifies the safety factor to use on the maximum explicitly-
  stable step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ERKStepSetCFLFraction(void *arkode_mem, realtype cfl_frac)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetCFLFraction",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetCFLFraction",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if (cfl_frac >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "ERKStepSetCFLFraction", "Illegal CFL fraction");
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
  ERKStepSetSafetyFactor:

  Specifies the safety factor to use on the error-based predicted
  time step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int ERKStepSetSafetyFactor(void *arkode_mem, realtype safety)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetSafetyFactor",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetSafetyFactoy",MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if (safety >= 1.0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "ERKStepSetSafetyFactor", "Illegal safety factor");
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
  ERKStepSetErrorBias:

  Specifies the error bias to use when performing adaptive-step
  error control.  Allowable values must be >= 1.0.  Any illegal
  value implies a reset to the default value.
  ---------------------------------------------------------------*/
int ERKStepSetErrorBias(void *arkode_mem, realtype bias)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetErrorBias",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetErrorBias", MSG_ARKADAPT_NO_MEM);
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
  ERKStepSetMaxGrowth:

  Specifies the maximum step size growth factor to be allowed
  between successive integration steps.  Note: the first step uses
  a separate maximum growth factor.  Allowable values must be
  > 1.0.  Any illegal value implies a reset to the default.
  ---------------------------------------------------------------*/
int ERKStepSetMaxGrowth(void *arkode_mem, realtype mx_growth)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetMaxGrowth",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetMaxGrowth", MSG_ARKADAPT_NO_MEM);
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
  ERKStepSetFixedStepBounds:

  Specifies the step size growth interval within which the step
  size will remain unchanged.  Allowable values must enclose the
  value 1.0.  Any illegal interval implies a reset to the default.
  ---------------------------------------------------------------*/
int ERKStepSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetFixedStepBounds",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetFixedStepBounds", MSG_ARKADAPT_NO_MEM);
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
  ERKStepSetAdaptivityMethod:

  Specifies the built-in time step adaptivity algorithm (and
  optionally, its associated parameters) to use.  All parameters
  will be checked for validity when used by the solver.
  ---------------------------------------------------------------*/
int ERKStepSetAdaptivityMethod(void *arkode_mem, int imethod,
                               int idefault, int pq,
                               realtype *adapt_params)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetAdaptivityMethod",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetAdaptivityMethod", MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* check for allowable parameters */
  if ((imethod > 5) || (imethod < 0)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "ERKStepSetAdaptivityMethod", "Illegal imethod");
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
  ERKStepSetAdaptivityFn:

  Specifies the user-provided time step adaptivity function to use.
  ---------------------------------------------------------------*/
int ERKStepSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun,
                           void *h_data)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetAdaptivityFn",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetAdaptivityFn", MSG_ARKADAPT_NO_MEM);
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
  ERKStepSetMaxFirstGrowth:

  Specifies the user-provided time step adaptivity constant
  etamx1.  Legal values are greater than 1.0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ERKStepSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetMaxFirstGrowth",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetMaxFirstGrowth",MSG_ARKADAPT_NO_MEM);
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
  ERKStepSetMaxEFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etamxf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ERKStepSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetMaxEFailGrowth",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetMaxEFailGrowth", MSG_ARKADAPT_NO_MEM);
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
  ERKStepSetSmallNumEFails:

  Specifies the user-provided time step adaptivity constant
  small_nef.  Legal values are > 0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int ERKStepSetSmallNumEFails(void *arkode_mem, int small_nef)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetSmallNumEFails",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetSmallNumEFails", MSG_ARKADAPT_NO_MEM);
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
  ERKStepSetStabilityFn:

  Specifies the user-provided explicit time step stability
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int ERKStepSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab,
                          void *estab_data)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  ARKodeHAdaptMem hadapt_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetStabilityFn",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* access structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSetStabilityFn", MSG_ARKADAPT_NO_MEM);
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
  ERKStepSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int ERKStepSetMaxErrTestFails(void *arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepSetMaxErrTestFails",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    step_mem->maxnef = MAXNEF;
  } else {
    step_mem->maxnef = maxnef;
  }

  return(ARK_SUCCESS);
}


/*===============================================================
  ERKStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepGetNumExpSteps:

  Returns the current number of stability-limited steps
  ---------------------------------------------------------------*/
int ERKStepGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetNumExpSteps",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* if step adaptivity structure not allocated, just return 0 */
  if (step_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = step_mem->hadapt_mem->nst_exp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetNumAccSteps:

  Returns the current number of accuracy-limited steps
  ---------------------------------------------------------------*/
int ERKStepGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetNumAccSteps",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* if step adaptivity structure not allocated, just return 0 */
  if (step_mem->hadapt_mem == NULL) {
    *nsteps = 0;
  } else {
    *nsteps = step_mem->hadapt_mem->nst_acc;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetNumStepAttempts:

  Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
int ERKStepGetNumStepAttempts(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetNumStepAttempts",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get value from step_mem */
  *nsteps = step_mem->nst_attempts;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetNumRhsEvals:

  Returns the current number of calls to fe and fi
  ---------------------------------------------------------------*/
int ERKStepGetNumRhsEvals(void *arkode_mem, long int *fevals)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetNumRhsEvals",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get values from step_mem */
  *fevals = step_mem->nfe;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
int ERKStepGetNumErrTestFails(void *arkode_mem, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetNumErrTestFails",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get value from step_mem */
  *netfails = step_mem->netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetCurrentButcherTable:

  Sets pointers to the Butcher table currently in use.
  ---------------------------------------------------------------*/
int ERKStepGetCurrentButcherTable(void *arkode_mem,
                                  ARKodeButcherTable *B)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetCurrentButcherTable",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get tables from step_mem */
  *B = step_mem->B;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetEstLocalErrors: (updated to the correct vector, but
  need to verify that it is unchanged between filling the
  estimated error and the end of the time step)

  Returns an estimate of the local error
  ---------------------------------------------------------------*/
int ERKStepGetEstLocalErrors(void *arkode_mem, N_Vector ele)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetEstLocalErrors",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* copy vector to output */
  N_VScale(ONE, ark_mem->tempv1, ele);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepGetTimestepperStats:

  Returns integrator statistics
  ---------------------------------------------------------------*/
int ERKStepGetTimestepperStats(void *arkode_mem, long int *expsteps,
                               long int *accsteps, long int *attempts,
                               long int *fevals, long int *netfails)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepGetTimestepperStats",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

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
  *attempts = step_mem->nst_attempts;
  *fevals   = step_mem->nfe;
  *netfails = step_mem->netf;

  return(ARK_SUCCESS);
}


/*===============================================================
  ERKStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int ERKStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepWriteParameters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* output ARKode infrastructure parameters first */
  retval = arkWriteParameters(arkode_mem, fp);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepWriteParameters",
                    "Error writing ARKode infrastructure parameters");
    return(retval);
  }

  /* print integrator parameters to file */
  fprintf(fp, "ERKStep time step module parameters:\n");
  fprintf(fp, "  Method order %i\n",step_mem->q);
  if (step_mem->hadapt_mem != NULL) {
    fprintf(fp, "  Maximum step increase (first step) = %"RSYM"\n",
            step_mem->hadapt_mem->etamx1);
    fprintf(fp, "  Step reduction factor on multiple error fails = %"RSYM"\n",
            step_mem->hadapt_mem->etamxf);
    fprintf(fp, "  Minimum error fails before above factor is used = %i\n",
            step_mem->hadapt_mem->small_nef);
    fprintf(fp, "  Step reduction factor on nonlinear convergence failure = %"RSYM"\n",
            step_mem->hadapt_mem->etacf);
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
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ERKStepWriteButcher:

  Outputs Butcher tables to the provided file pointer.
  ---------------------------------------------------------------*/
int ERKStepWriteButcher(void *arkode_mem, FILE *fp)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepWriteButcher",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check that Butcher table is non-NULL (otherwise report error) */
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepWriteButcher", "Butcher table memory is NULL");
    return(ARK_MEM_NULL);
  }

  /* print Butcher table to file */
  fprintf(fp, "\nERKStep Butcher table (stages = %i):\n", step_mem->stages);
  ARKodeButcherTable_Write(step_mem->B, fp);
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
