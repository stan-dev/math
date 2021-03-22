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
 * output functions for the ARKode infrastructure; these routines
 * should not be called directly by the user; instead they are
 * provided as utility routines for ARKode time-step modules
 * to use.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  ARKode optional input utility functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkSetDefaults:

  Resets all optional inputs to ARKode default values.  Does not
  change problem-defining function pointers fe and fi or
  user_data pointer.  Also leaves alone any data
  structures/options related to root-finding (those can be reset
  using ARKodeRootInit) or post-processing a step (ProcessStep).
  ---------------------------------------------------------------*/
int arkSetDefaults(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetDefaults", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Set default values for integrator optional inputs */
  ark_mem->fixedstep               = SUNFALSE;       /* default to use adaptive steps */
  ark_mem->reltol                  = RCONST(1.e-4);  /* relative tolerance */
  ark_mem->itol                    = ARK_SS;         /* scalar-scalar solution tolerances */
  ark_mem->ritol                   = ARK_SS;         /* scalar-scalar residual tolerances */
  ark_mem->Sabstol                 = RCONST(1.e-9);  /* solution absolute tolerance */
  ark_mem->atolmin0                = SUNFALSE;       /* min(abstol) > 0 */
  ark_mem->SRabstol                = RCONST(1.e-9);  /* residual absolute tolerance */
  ark_mem->Ratolmin0               = SUNFALSE;       /* min(Rabstol) > 0 */
  ark_mem->user_efun               = SUNFALSE;       /* no user-supplied ewt function */
  ark_mem->efun                    = arkEwtSetSS;    /* built-in scalar-scalar ewt function */
  ark_mem->e_data                  = ark_mem;        /* ewt function data */
  ark_mem->user_rfun               = SUNFALSE;       /* no user-supplied rwt function */
  ark_mem->rfun                    = arkRwtSet;      /* built-in rwt function */
  ark_mem->r_data                  = ark_mem;        /* rwt function data */
  ark_mem->ehfun                   = arkErrHandler;  /* default error handler fn */
  ark_mem->eh_data                 = ark_mem;        /* error handler data */
  ark_mem->errfp                   = stderr;         /* output stream for errors */
  ark_mem->mxstep                  = MXSTEP_DEFAULT; /* max number of steps */
  ark_mem->mxhnil                  = MXHNIL;         /* max warns of t+h==t */
  ark_mem->maxnef                  = MAXNEF;         /* max error test fails */
  ark_mem->maxncf                  = MAXNCF;         /* max convergence fails */
  ark_mem->maxconstrfails          = MAXCONSTRFAILS; /* max number of constraint fails */
  ark_mem->hin                     = ZERO;           /* determine initial step on-the-fly */
  ark_mem->hmin                    = ZERO;           /* no minimum step size */
  ark_mem->hmax_inv                = ZERO;           /* no maximum step size */
  ark_mem->tstopset                = SUNFALSE;       /* no stop time set */
  ark_mem->tstop                   = ZERO;           /* no fixed stop time */
  ark_mem->diagfp                  = NULL;           /* no solver diagnostics file */
  ark_mem->report                  = SUNFALSE;       /* don't report solver diagnostics */
  ark_mem->hadapt_mem->etamx1      = ETAMX1;         /* max change on first step */
  ark_mem->hadapt_mem->etamxf      = ETAMXF;         /* max change on error-failed step */
  ark_mem->hadapt_mem->etamin      = ETAMIN;         /* min bound on time step reduction */
  ark_mem->hadapt_mem->small_nef   = SMALL_NEF;      /* num error fails before ETAMXF enforced */
  ark_mem->hadapt_mem->etacf       = ETACF;          /* max change on convergence failure */
  ark_mem->hadapt_mem->HAdapt      = NULL;           /* step adaptivity fn */
  ark_mem->hadapt_mem->HAdapt_data = NULL;           /* step adaptivity data */
  ark_mem->hadapt_mem->imethod     = ARK_ADAPT_PID;  /* PID controller */
  ark_mem->hadapt_mem->cfl         = CFLFAC;         /* explicit stability factor */
  ark_mem->hadapt_mem->safety      = SAFETY;         /* step adaptivity safety factor  */
  ark_mem->hadapt_mem->bias        = BIAS;           /* step adaptivity error bias */
  ark_mem->hadapt_mem->growth      = GROWTH;         /* step adaptivity growth factor */
  ark_mem->hadapt_mem->lbound      = HFIXED_LB;      /* step adaptivity no-change lower bound */
  ark_mem->hadapt_mem->ubound      = HFIXED_UB;      /* step adaptivity no-change upper bound */
  ark_mem->hadapt_mem->k1          = AD0_K1;         /* step adaptivity parameter */
  ark_mem->hadapt_mem->k2          = AD0_K2;         /* step adaptivity parameter */
  ark_mem->hadapt_mem->k3          = AD0_K3;         /* step adaptivity parameter */
  ark_mem->hadapt_mem->pq          = SUNFALSE;       /* use embedding order */
  ark_mem->hadapt_mem->expstab     = arkExpStab;     /* internal explicit stability fn */
  ark_mem->hadapt_mem->estab_data  = NULL;           /* no explicit stability fn data */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetInterpolantType:

  Specifies use of the Lagrange or Hermite interpolation modules.
    itype == ARK_INTERP_HERMITE specifies the Hermite (nonstiff)
      interpolation module.
    itype == ARK_INTERP_LAGRANGE specifies the Lagrange (stiff)
      interpolation module.

  Return values:
     ARK_SUCCESS on success.
     ARK_MEM_NULL on NULL-valued arkode_mem input.
     ARK_MEM_FAIL if the interpolation module cannot be allocated.
     ARK_ILL_INPUT if the itype argument is not recognized.
  ---------------------------------------------------------------*/
int arkSetInterpolantType(void *arkode_mem, int itype)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetInterpolantType", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal itype input */
  if ((itype != ARK_INTERP_HERMITE) && (itype != ARK_INTERP_LAGRANGE)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSetInterpolantType",
                    "Illegal interpolation type input.");
    return(ARK_ILL_INPUT);
  }

  /* do not change type once the module has been initialized */
  if (ark_mem->initialized) {
    arkProcessError(ark_mem, ARK_INTERP_FAIL, "ARKode",
                    "arkSetInterpolantType",
                    "Type cannot be specified after module initialization.");
    return(ARK_ILL_INPUT);
  }

  /* delete any existing interpolation module */
  if (ark_mem->interp != NULL) {
    arkInterpFree(ark_mem, ark_mem->interp);
    ark_mem->interp = NULL;
  }

  /* create requested interpolation module, initially specifying
     the maximum possible interpolant degree. */
  if (itype == ARK_INTERP_HERMITE) {
    ark_mem->interp = arkInterpCreate_Hermite(arkode_mem, ARK_INTERP_MAX_DEGREE);
  } else {
    ark_mem->interp = arkInterpCreate_Lagrange(arkode_mem, ARK_INTERP_MAX_DEGREE);
  }
  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode", "arkSetInterpolantType",
                    "Unable to allocate interpolation structure");
    return(ARK_MEM_FAIL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetInterpolantDegree:

  Specifies the polynomial degree for the dense output
  interpolation module.

  Return values:
     ARK_SUCCESS on success.
     ARK_MEM_NULL on NULL-valued arkode_mem input or nonexistent
       interpolation module.
     ARK_INTERP_FAIL if the interpolation module is already
       initialized.
     ARK_ILL_INPUT if the degree is illegal.
  ---------------------------------------------------------------*/
int arkSetInterpolantDegree(void *arkode_mem, int degree)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetInterpolantDegree", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode",
                    "arkSetInterpolantDegree",
                    "Interpolation module is not yet allocated");
    return(ARK_MEM_NULL);
  }

  /* do not change degree once the module has been initialized */
  if (ark_mem->initialized) {
    arkProcessError(ark_mem, ARK_INTERP_FAIL, "ARKode",
                    "arkSetInterpolantType",
                    "Degree cannot be specified after module initialization.");
    return(ARK_ILL_INPUT);
  }

  /* pass 'degree' to interpolation module, returning its value */
  return(arkInterpSetDegree(ark_mem, ark_mem->interp, degree));
}


/*---------------------------------------------------------------
  arkSetErrHandlerFn:

  Specifies the error handler function
  ---------------------------------------------------------------*/
int arkSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                       void *eh_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetErrHandlerFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set user-provided values, or defaults, depending on argument */
  if (ehfun == NULL) {
    ark_mem->ehfun   = arkErrHandler;
    ark_mem->eh_data = ark_mem;
  } else {
    ark_mem->ehfun   = ehfun;
    ark_mem->eh_data = eh_data;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetErrFile:

  Specifies the FILE pointer for output (NULL means no messages)
  ---------------------------------------------------------------*/
int arkSetErrFile(void *arkode_mem, FILE *errfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetErrFile", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->errfp = errfp;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetUserData:

  Specifies the user data pointer for f
  ---------------------------------------------------------------*/
int arkSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetUserData", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->user_data = user_data;

  /* Set data for efun */
  if (ark_mem->user_efun)
    ark_mem->e_data = user_data;

  /* Set data for rfun */
  if (ark_mem->user_rfun)
    ark_mem->r_data = user_data;

  /* Set data for root finding */
  if (ark_mem->root_mem != NULL)
    ark_mem->root_mem->root_data = user_data;

  /* Set data for post-processing a step */
  if (ark_mem->ProcessStep != NULL)
    ark_mem->ps_data = user_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetDiagnostics:

  Specifies to enable solver diagnostics, and specifies the FILE
  pointer for output (diagfp==NULL disables output)
  ---------------------------------------------------------------*/
int arkSetDiagnostics(void *arkode_mem, FILE *diagfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetDiagnostics", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->diagfp = diagfp;
  if (diagfp != NULL) {
    ark_mem->report = SUNTRUE;
  } else {
    ark_mem->report = SUNFALSE;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxNumSteps:

  Specifies the maximum number of integration steps
  ---------------------------------------------------------------*/
int arkSetMaxNumSteps(void *arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMaxNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */
  if (mxsteps == 0)
    ark_mem->mxstep = MXSTEP_DEFAULT;
  else
    ark_mem->mxstep = mxsteps;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxHnilWarns:

  Specifies the maximum number of warnings for small h
  ---------------------------------------------------------------*/
int arkSetMaxHnilWarns(void *arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMaxHnilWarns", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing mxhnil=0 sets the default, otherwise use input. */
  if (mxhnil == 0) {
    ark_mem->mxhnil = 10;
  } else {
    ark_mem->mxhnil = mxhnil;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetInitStep:

  Specifies the initial step size
  ---------------------------------------------------------------*/
int arkSetInitStep(void *arkode_mem, realtype hin)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing hin=0 sets the default, otherwise use input. */
  if (hin == ZERO) {
    ark_mem->hin = ZERO;
  } else {
    ark_mem->hin = hin;
  }

  /* Clear previous initial step */
  ark_mem->h0u = ZERO;

  /* Clear error and step size history */
  ark_mem->hadapt_mem->ehist[0] = ONE;
  ark_mem->hadapt_mem->ehist[1] = ONE;
  ark_mem->hadapt_mem->hhist[0] = ZERO;
  ark_mem->hadapt_mem->hhist[1] = ZERO;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMinStep:

  Specifies the minimum step size
  ---------------------------------------------------------------*/
int arkSetMinStep(void *arkode_mem, realtype hmin)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMinStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing a value <= 0 sets hmin = 0 */
  if (hmin <= ZERO) {
    ark_mem->hmin = ZERO;
    return(ARK_SUCCESS);
  }

  /* check that hmin and hmax are agreeable */
  if (hmin * ark_mem->hmax_inv > ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSetMinStep", MSG_ARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->hmin = hmin;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxStep:

  Specifies the maximum step size
  ---------------------------------------------------------------*/
int arkSetMaxStep(void *arkode_mem, realtype hmax)
{
  realtype hmax_inv;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMaxStep", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing a value <= 0 sets hmax = infinity */
  if (hmax <= ZERO) {
    ark_mem->hmax_inv = ZERO;
    return(ARK_SUCCESS);
  }

  /* check that hmax and hmin are agreeable */
  hmax_inv = ONE/hmax;
  if (hmax_inv * ark_mem->hmin > ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSetMaxStep", MSG_ARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  /* set the value */
  ark_mem->hmax_inv = hmax_inv;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetStopTime:

  Specifies the time beyond which the integration is not to proceed.
  ---------------------------------------------------------------*/
int arkSetStopTime(void *arkode_mem, realtype tstop)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetStopTime", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* If ARKode was called at least once, test if tstop is legal
     (i.e. if it was not already passed).
     If arkSetStopTime is called before the first call to ARKode,
     tstop will be checked in ARKode. */
  if (ark_mem->nst > 0) {
    if ( (tstop - ark_mem->tcur) * ark_mem->h < ZERO ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                      "arkSetStopTime", MSG_ARK_BAD_TSTOP,
                      tstop, ark_mem->tcur);
      return(ARK_ILL_INPUT);
    }
  }

  ark_mem->tstop    = tstop;
  ark_mem->tstopset = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetFixedStep:

  Specifies to use a fixed time step size instead of performing
  any form of temporal adaptivity.  ARKode will use this step size
  for all steps (unless tstop is set, in which case it may need to
  modify that last step approaching tstop.  If any solver failure
  occurs in the timestepping module, ARKode will typically
  immediately return with an error message indicating that the
  selected step size cannot be used.

  Any nonzero argument will result in the use of that fixed step
  size; an argument of 0 will re-enable temporal adaptivity.
  ---------------------------------------------------------------*/
int arkSetFixedStep(void *arkode_mem, realtype hfixed)
{
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetFixedStep", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* re-attach internal error weight functions if necessary */
  if ((hfixed == ZERO) && (!ark_mem->user_efun)) {
    if (ark_mem->itol == ARK_SV && ark_mem->Vabstol != NULL)
      retval = arkSVtolerances(ark_mem, ark_mem->reltol, ark_mem->Vabstol);
    else
      retval = arkSStolerances(ark_mem, ark_mem->reltol, ark_mem->Sabstol);
    if (retval != ARK_SUCCESS) return(retval);
  }

  /* set ark_mem entry */
  if (hfixed != ZERO) {
    ark_mem->fixedstep = SUNTRUE;
    ark_mem->hin = hfixed;
  } else {
    ark_mem->fixedstep = SUNFALSE;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetRootDirection:

  Specifies the direction of zero-crossings to be monitored.
  The default is to monitor both crossings.
  ---------------------------------------------------------------*/
int arkSetRootDirection(void *arkode_mem, int *rootdir)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  int i;

  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetRootDirection", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode",
                    "arkSetRootDirection", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;

  if (ark_root_mem->nrtfn == 0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSetRootDirection", MSG_ARK_NO_ROOT);
    return(ARK_ILL_INPUT);
  }
  for(i=0; i<ark_root_mem->nrtfn; i++)
    ark_root_mem->rootdir[i] = rootdir[i];
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetNoInactiveRootWarn:

  Disables issuing a warning if some root function appears
  to be identically zero at the beginning of the integration
  ---------------------------------------------------------------*/
int arkSetNoInactiveRootWarn(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetNoInactiveRootWarn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode",
                    "arkSetNoInactiveRootWarn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;
  ark_root_mem->mxgnull = 0;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetPostprocessStepFn:

  Specifies a user-provided step postprocessing function having
  type ARKPostProcessFn.  A NULL input function disables step
  postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.
  ---------------------------------------------------------------*/
int arkSetPostprocessStepFn(void *arkode_mem,
                            ARKPostProcessFn ProcessStep)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetPostprocessStepFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* NULL argument sets default, otherwise set inputs */
  ark_mem->ProcessStep = ProcessStep;
  ark_mem->ps_data     = ark_mem->user_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetPostprocessStageFn:

  Specifies a user-provided stage postprocessing function having
  type ARKPostProcessFn.  A NULL input function disables
  stage postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.
  ---------------------------------------------------------------*/
int arkSetPostprocessStageFn(void *arkode_mem,
                             ARKPostProcessFn ProcessStage)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetPostprocessStageFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* NULL argument sets default, otherwise set inputs */
  ark_mem->ProcessStage = ProcessStage;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetConstraints:

  Activates or Deactivates inequality constraint checking.
  ---------------------------------------------------------------*/
int arkSetConstraints(void *arkode_mem, N_Vector constraints)
{
  realtype temptest;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetConstraints", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* If there are no constraints, destroy data structures */
  if (constraints == NULL) {
    arkFreeVec(ark_mem, &ark_mem->constraints);
    ark_mem->constraintsSet = SUNFALSE;
    return(ARK_SUCCESS);
  }

  /* Test if required vector ops. are defined */
  if (constraints->ops->nvdiv         == NULL ||
      constraints->ops->nvmaxnorm     == NULL ||
      constraints->ops->nvcompare     == NULL ||
      constraints->ops->nvconstrmask  == NULL ||
      constraints->ops->nvminquotient == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetConstraints", MSG_ARK_BAD_NVECTOR);
    return(ARK_ILL_INPUT);
  }

  /* Check the constraints vector */
  temptest = N_VMaxNorm(constraints);
  if ((temptest > RCONST(2.5)) || (temptest < HALF)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetConstraints", MSG_ARK_BAD_CONSTR);
    return(ARK_ILL_INPUT);
  }

  /* Allocate the internal constrains vector (if necessary) */
  if (!arkAllocVec(ark_mem, constraints, &ark_mem->constraints))
    return(ARK_MEM_FAIL);

  /* Load the constraints vector */
  N_VScale(ONE, constraints, ark_mem->constraints);
  ark_mem->constraintsSet = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxNumConstrFails:

  Set max number of allowed constraint failures in a step before
  returning an error
  ---------------------------------------------------------------*/
int arkSetMaxNumConstrFails(void *arkode_mem, int maxfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMaxNumConstrFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing maxfails = 0 sets the default, otherwise set to input */
  if (maxfails <= 0)
    ark_mem->maxconstrfails = MAXCONSTRFAILS;
  else
    ark_mem->maxconstrfails = maxfails;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetCFLFraction:

  Specifies the safety factor to use on the maximum explicitly-
  stable step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int arkSetCFLFraction(void *arkode_mem, realtype cfl_frac)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetCFLFraction",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* check for allowable parameters */
  if (cfl_frac >= ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSetCFLFraction", "Illegal CFL fraction");
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
  arkSetSafetyFactor:

  Specifies the safety factor to use on the error-based predicted
  time step size.  Allowable values must be within the open
  interval (0,1).  A non-positive input implies a reset to the
  default value.
  ---------------------------------------------------------------*/
int arkSetSafetyFactor(void *arkode_mem, realtype safety)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetSafetyFactor",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* check for allowable parameters */
  if (safety >= ONE) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSetSafetyFactor", "Illegal safety factor");
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
  arkSetErrorBias:

  Specifies the error bias to use when performing adaptive-step
  error control.  Allowable values must be >= 1.0.  Any illegal
  value implies a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetErrorBias(void *arkode_mem, realtype bias)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetErrorBias",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set allowed value, otherwise set default */
  if (bias < ONE) {
    hadapt_mem->bias = BIAS;
  } else {
    hadapt_mem->bias = bias;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxGrowth:

  Specifies the maximum step size growth factor to be allowed
  between successive integration steps.  Note: the first step uses
  a separate maximum growth factor.  Allowable values must be
  > 1.0.  Any illegal value implies a reset to the default.
  ---------------------------------------------------------------*/
int arkSetMaxGrowth(void *arkode_mem, realtype mx_growth)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMaxGrowth",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set allowed value, otherwise set default */
  if (mx_growth <= ONE) {
    hadapt_mem->growth = GROWTH;
  } else {
    hadapt_mem->growth = mx_growth;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMinReduction:

  Specifies the minimum possible step size reduction factor to be
  allowed between successive integration steps. Allowable values
  must be > 0.0 and < 1.0. Any illegal value implies a reset to
  the default.
  ---------------------------------------------------------------*/
int arkSetMinReduction(void *arkode_mem, realtype eta_min)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMinReduction",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* set allowed value, otherwise set default */
  if (eta_min >= ONE || eta_min <= ZERO) {
    hadapt_mem->etamin = ETAMIN;
  } else {
    hadapt_mem->etamin = eta_min;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetFixedStepBounds:

  Specifies the step size growth interval within which the step
  size will remain unchanged.  Allowable values must enclose the
  value 1.0.  Any illegal interval implies a reset to the default.
  ---------------------------------------------------------------*/
int arkSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetFixedStepBounds",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set allowable interval, otherwise set defaults */
  if ((lb <= ONE) && (ub >= ONE)) {
    hadapt_mem->lbound = lb;
    hadapt_mem->ubound = ub;
  } else {
    hadapt_mem->lbound = HFIXED_LB;
    hadapt_mem->ubound = HFIXED_UB;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetAdaptivityMethod:

  Specifies the built-in time step adaptivity algorithm (and
  optionally, its associated parameters) to use.  All parameters
  will be checked for validity when used by the solver.
  ---------------------------------------------------------------*/
int arkSetAdaptivityMethod(void *arkode_mem, int imethod, int idefault,
                           int pq, realtype adapt_params[3])
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetAdaptivityMethod",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* check for allowable parameters */
  if ((imethod > ARK_ADAPT_IMEX_GUS) || (imethod < ARK_ADAPT_PID)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  /* set adaptivity method */
  hadapt_mem->imethod = imethod;

  /* set flag whether to use p (embedding, 0) or q (method, 1) order */
  hadapt_mem->pq = (pq != 0);

  /* set method parameters */
  if (idefault == 1) {
    switch (hadapt_mem->imethod) {
    case (ARK_ADAPT_PID):
      hadapt_mem->k1 = AD0_K1;
      hadapt_mem->k2 = AD0_K2;
      hadapt_mem->k3 = AD0_K3; break;
    case (ARK_ADAPT_PI):
      hadapt_mem->k1 = AD1_K1;
      hadapt_mem->k2 = AD1_K2; break;
    case (ARK_ADAPT_I):
      hadapt_mem->k1 = AD2_K1; break;
    case (ARK_ADAPT_EXP_GUS):
      hadapt_mem->k1 = AD3_K1;
      hadapt_mem->k2 = AD3_K2; break;
    case (ARK_ADAPT_IMP_GUS):
      hadapt_mem->k1 = AD4_K1;
      hadapt_mem->k2 = AD4_K2; break;
    case (ARK_ADAPT_IMEX_GUS):
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
  arkSetAdaptivityFn:

  Specifies the user-provided time step adaptivity function to use.
  ---------------------------------------------------------------*/
int arkSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun, void *h_data)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetAdaptivityFn",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* NULL hfun sets default, otherwise set inputs */
  if (hfun == NULL) {
    hadapt_mem->HAdapt      = NULL;
    hadapt_mem->HAdapt_data = NULL;
    hadapt_mem->imethod     = ARK_ADAPT_PID;
  } else {
    hadapt_mem->HAdapt      = hfun;
    hadapt_mem->HAdapt_data = h_data;
    hadapt_mem->imethod     = ARK_ADAPT_CUSTOM;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxFirstGrowth:

  Specifies the user-provided time step adaptivity constant
  etamx1.  Legal values are greater than 1.0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMaxFirstGrowth",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (etamx1 <= ONE) {
    hadapt_mem->etamx1 = ETAMX1;
  } else {
    hadapt_mem->etamx1 = etamx1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxEFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etamxf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMaxEFailGrowth",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if ((etamxf <= ZERO) || (etamxf > ONE)) {
    hadapt_mem->etamxf = ETAMXF;
  } else {
    hadapt_mem->etamxf = etamxf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetSmallNumEFails:

  Specifies the user-provided time step adaptivity constant
  small_nef.  Legal values are > 0.  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetSmallNumEFails(void *arkode_mem, int small_nef)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetSmallNumEFails",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if (small_nef <= 0) {
    hadapt_mem->small_nef = SMALL_NEF;
  } else {
    hadapt_mem->small_nef = small_nef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxCFailGrowth:

  Specifies the user-provided time step adaptivity constant
  etacf. Legal values are in the interval (0,1].  Illegal values
  imply a reset to the default value.
  ---------------------------------------------------------------*/
int arkSetMaxCFailGrowth(void *arkode_mem, realtype etacf)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetMaxCFailGrowth",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if argument legal set it, otherwise set default */
  if ((etacf <= ZERO) || (etacf > ONE)) {
    hadapt_mem->etacf = ETACF;
  } else {
    hadapt_mem->etacf = etacf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetStabilityFn:

  Specifies the user-provided explicit time step stability
  function to use.  A NULL input function implies a reset to
  the default function (empty).
  ---------------------------------------------------------------*/
int arkSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab, void *estab_data)
{
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeMem ark_mem;
  retval = arkAccessHAdaptMem(arkode_mem, "arkSetStabilityFn",
                              &ark_mem, &hadapt_mem);
  if (retval != ARK_SUCCESS)  return(retval);

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
  arkSetMaxErrTestFails:

  Specifies the maximum number of error test failures during one
  step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int arkSetMaxErrTestFails(void *arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMaxErrTestFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    ark_mem->maxnef = MAXNEF;
  } else {
    ark_mem->maxnef = maxnef;
  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMaxConvFails:

  Specifies the maximum number of nonlinear convergence failures
  during one step try.  A non-positive input implies a reset to
  the default value.
  ---------------------------------------------------------------*/
int arkSetMaxConvFails(void *arkode_mem, int maxncf)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMaxConvFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* argument <= 0 sets default, otherwise set input */
  if (maxncf <= 0) {
    ark_mem->maxncf = MAXNCF;
  } else {
    ark_mem->maxncf = maxncf;
  }
  return(ARK_SUCCESS);
}



/*===============================================================
  ARKode optional output utility functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkGetNumStepAttempts:

   Returns the current number of steps attempted by the solver
  ---------------------------------------------------------------*/
int arkGetNumStepAttempts(void *arkode_mem, long int *nstep_attempts)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetNumStepAttempts", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nstep_attempts = ark_mem->nst_attempts;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumSteps:

  Returns the current number of integration steps
  ---------------------------------------------------------------*/
int arkGetNumSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->nst;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetActualInitStep:

  Returns the step size used on the first step
  ---------------------------------------------------------------*/
int arkGetActualInitStep(void *arkode_mem, realtype *hinused)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetActualInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *hinused = ark_mem->h0u;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetLastStep:

  Returns the step size used on the last successful step
  ---------------------------------------------------------------*/
int arkGetLastStep(void *arkode_mem, realtype *hlast)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetLastStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *hlast = ark_mem->hold;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetCurrentStep:

  Returns the step size to be attempted on the next step
  ---------------------------------------------------------------*/
int arkGetCurrentStep(void *arkode_mem, realtype *hcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetCurrentStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *hcur = ark_mem->next_h;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetCurrentState:

  Returns the current solution (before or after as step) or
  stage value (during step solve).
  ---------------------------------------------------------------*/
int arkGetCurrentState(void *arkode_mem, N_Vector *state)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetCurrentState", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *state = ark_mem->ycur;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetCurrentTime:

  Returns the current value of the independent variable
  ---------------------------------------------------------------*/
int arkGetCurrentTime(void *arkode_mem, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetCurrentTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *tcur = ark_mem->tcur;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetTolScaleFactor:

  Returns a suggested factor for scaling tolerances
  ---------------------------------------------------------------*/
int arkGetTolScaleFactor(void *arkode_mem, realtype *tolsfact)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetTolScaleFactor", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *tolsfact = ark_mem->tolsf;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetErrWeights:

  This routine returns the current error weight vector.
  ---------------------------------------------------------------*/
int arkGetErrWeights(void *arkode_mem, N_Vector eweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetErrWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->ewt, eweight);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetResWeights:

  This routine returns the current residual weight vector.
  ---------------------------------------------------------------*/
int arkGetResWeights(void *arkode_mem, N_Vector rweight)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetResWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->rwt, rweight);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetWorkSpace:

  Returns integrator work space requirements
  ---------------------------------------------------------------*/
int arkGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetWorkSpace", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *leniw = ark_mem->liw;
  *lenrw = ark_mem->lrw;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumGEvals:

  Returns the current number of calls to g (for rootfinding)
  ---------------------------------------------------------------*/
int arkGetNumGEvals(void *arkode_mem, long int *ngevals)
{
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetNumGEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode",
                    "arkGetNumGEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;
  *ngevals = ark_root_mem->nge;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetRootInfo:

  Returns pointer to array rootsfound showing roots found
  ---------------------------------------------------------------*/
int arkGetRootInfo(void *arkode_mem, int *rootsfound)
{
  int i;
  ARKodeMem ark_mem;
  ARKodeRootMem ark_root_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetRootInfo", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->root_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode",
                    "arkGetRootInfo", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_root_mem = (ARKodeRootMem) ark_mem->root_mem;
  for (i=0; i<ark_root_mem->nrtfn; i++)
    rootsfound[i] = ark_root_mem->iroots[i];
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetStepStats:

  Returns step statistics
  ---------------------------------------------------------------*/
int arkGetStepStats(void *arkode_mem, long int *nsteps,
                    realtype *hinused, realtype *hlast,
                    realtype *hcur, realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetStepStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps  = ark_mem->nst;
  *hinused = ark_mem->h0u;
  *hlast   = ark_mem->hold;
  *hcur    = ark_mem->next_h;
  *tcur    = ark_mem->tcur;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumConstrFails:

  Returns the current number of constraint fails
  ---------------------------------------------------------------*/
int arkGetNumConstrFails(void *arkode_mem, long int *nconstrfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetNumConstrFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nconstrfails = ark_mem->nconstrfails;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumExpSteps:

  Returns the current number of stability-limited steps
  ---------------------------------------------------------------*/
int arkGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetNumExpSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->hadapt_mem->nst_exp;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumAccSteps:

  Returns the current number of accuracy-limited steps
  ---------------------------------------------------------------*/
int arkGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetNumAccSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->hadapt_mem->nst_acc;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumErrTestFails:

  Returns the current number of error test failures
  ---------------------------------------------------------------*/
int arkGetNumErrTestFails(void *arkode_mem, long int *netfails)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetNumErrTestFails", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *netfails = ark_mem->netf;
  return(ARK_SUCCESS);
}

/*-----------------------------------------------------------------*/

char *arkGetReturnFlagName(long int flag)
{
  char *name;
  name = (char *)malloc(27*sizeof(char));

  switch(flag) {
  case ARK_SUCCESS:
    sprintf(name,"ARK_SUCCESS");
    break;
  case ARK_TSTOP_RETURN:
    sprintf(name,"ARK_TSTOP_RETURN");
    break;
  case ARK_ROOT_RETURN:
    sprintf(name,"ARK_ROOT_RETURN");
    break;
  case ARK_TOO_MUCH_WORK:
    sprintf(name,"ARK_TOO_MUCH_WORK");
    break;
  case ARK_TOO_MUCH_ACC:
    sprintf(name,"ARK_TOO_MUCH_ACC");
    break;
  case ARK_ERR_FAILURE:
    sprintf(name,"ARK_ERR_FAILURE");
    break;
  case ARK_CONV_FAILURE:
    sprintf(name,"ARK_CONV_FAILURE");
    break;
  case ARK_LINIT_FAIL:
    sprintf(name,"ARK_LINIT_FAIL");
    break;
  case ARK_LSETUP_FAIL:
    sprintf(name,"ARK_LSETUP_FAIL");
    break;
  case ARK_LSOLVE_FAIL:
    sprintf(name,"ARK_LSOLVE_FAIL");
    break;
  case ARK_RHSFUNC_FAIL:
    sprintf(name,"ARK_RHSFUNC_FAIL");
    break;
  case ARK_FIRST_RHSFUNC_ERR:
    sprintf(name,"ARK_FIRST_RHSFUNC_ERR");
    break;
  case ARK_REPTD_RHSFUNC_ERR:
    sprintf(name,"ARK_REPTD_RHSFUNC_ERR");
    break;
  case ARK_UNREC_RHSFUNC_ERR:
    sprintf(name,"ARK_UNREC_RHSFUNC_ERR");
    break;
  case ARK_RTFUNC_FAIL:
    sprintf(name,"ARK_RTFUNC_FAIL");
    break;
  case ARK_LFREE_FAIL:
    sprintf(name,"ARK_LFREE_FAIL");
    break;
  case ARK_MASSINIT_FAIL:
    sprintf(name,"ARK_MASSINIT_FAIL");
    break;
  case ARK_MASSSETUP_FAIL:
    sprintf(name,"ARK_MASSSETUP_FAIL");
    break;
  case ARK_MASSSOLVE_FAIL:
    sprintf(name,"ARK_MASSSOLVE_FAIL");
    break;
  case ARK_MASSFREE_FAIL:
    sprintf(name,"ARK_MASSFREE_FAIL");
    break;
  case ARK_MASSMULT_FAIL:
    sprintf(name,"ARK_MASSMULT_FAIL");
    break;
  case ARK_MEM_FAIL:
    sprintf(name,"ARK_MEM_FAIL");
    break;
  case ARK_MEM_NULL:
    sprintf(name,"ARK_MEM_NULL");
    break;
  case ARK_ILL_INPUT:
    sprintf(name,"ARK_ILL_INPUT");
    break;
  case ARK_NO_MALLOC:
    sprintf(name,"ARK_NO_MALLOC");
    break;
  case ARK_BAD_K:
    sprintf(name,"ARK_BAD_K");
    break;
  case ARK_BAD_T:
    sprintf(name,"ARK_BAD_T");
    break;
  case ARK_BAD_DKY:
    sprintf(name,"ARK_BAD_DKY");
    break;
  case ARK_TOO_CLOSE:
    sprintf(name,"ARK_TOO_CLOSE");
    break;
  case ARK_POSTPROCESS_STEP_FAIL:
    sprintf(name,"ARK_POSTPROCESS_STEP_FAIL");
    break;
  case ARK_POSTPROCESS_STAGE_FAIL:
    sprintf(name,"ARK_POSTPROCESS_STAGE_FAIL");
    break;
  case ARK_VECTOROP_ERR:
    sprintf(name,"ARK_VECTOROP_ERR");
    break;
  case ARK_NLS_INIT_FAIL:
    sprintf(name,"ARK_NLS_INIT_FAIL");
    break;
  case ARK_NLS_SETUP_FAIL:
    sprintf(name,"ARK_NLS_SETUP_FAIL");
    break;
  case ARK_NLS_OP_ERR:
    sprintf(name,"ARK_NLS_OP_ERR");
    break;
  case ARK_INNERSTEP_ATTACH_ERR:
    sprintf(name,"ARK_INNERSTEP_ATTACH_ERR");
    break;
  case ARK_INNERSTEP_FAIL:
    sprintf(name,"ARK_INNERSTEP_FAIL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}



/*===============================================================
  ARKode parameter output utility routine
  ===============================================================*/

/*---------------------------------------------------------------
  arkodeWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int arkWriteParameters(ARKodeMem ark_mem, FILE *fp)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkWriteParameters", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* print integrator parameters to file */
  fprintf(fp, "ARKode solver parameters:\n");
  if (ark_mem->hmin != ZERO)
    fprintf(fp, "  Minimum step size = %" RSYM"\n",ark_mem->hmin);
  if (ark_mem->hmax_inv != ZERO)
    fprintf(fp, "  Maximum step size = %" RSYM"\n",ONE/ark_mem->hmax_inv);
  if (ark_mem->fixedstep)
    fprintf(fp, "  Fixed time-stepping enabled\n");
  if (ark_mem->itol == ARK_WF) {
    fprintf(fp, "  User provided error weight function\n");
  } else {
    fprintf(fp, "  Solver relative tolerance = %" RSYM"\n", ark_mem->reltol);
    if (ark_mem->itol == ARK_SS) {
      fprintf(fp, "  Solver absolute tolerance = %" RSYM"\n", ark_mem->Sabstol);
    } else {
      fprintf(fp, "  Vector-valued solver absolute tolerance\n");
    }
  }
  if (!ark_mem->rwt_is_ewt) {
    if (ark_mem->ritol == ARK_WF) {
      fprintf(fp, "  User provided residual weight function\n");
    } else {
      if (ark_mem->ritol == ARK_SS) {
        fprintf(fp, "  Absolute residual tolerance = %" RSYM"\n", ark_mem->SRabstol);
      } else {
        fprintf(fp, "  Vector-valued residual absolute tolerance\n");
      }
    }
  }
  if (ark_mem->hin != ZERO)
    fprintf(fp, "  Initial step size = %" RSYM"\n",ark_mem->hin);
  fprintf(fp, "\n");
  fprintf(fp, "  Maximum step increase (first step) = %"RSYM"\n",
          ark_mem->hadapt_mem->etamx1);
  fprintf(fp, "  Step reduction factor on multiple error fails = %"RSYM"\n",
          ark_mem->hadapt_mem->etamxf);
  fprintf(fp, "  Minimum error fails before above factor is used = %i\n",
          ark_mem->hadapt_mem->small_nef);
  fprintf(fp, "  Step reduction factor on nonlinear convergence failure = %"RSYM"\n",
          ark_mem->hadapt_mem->etacf);
  fprintf(fp, "  Explicit safety factor = %"RSYM"\n",
          ark_mem->hadapt_mem->cfl);
  if (ark_mem->hadapt_mem->HAdapt == NULL) {
    fprintf(fp, "  Time step adaptivity method %i\n", ark_mem->hadapt_mem->imethod);
    fprintf(fp, "     Safety factor = %"RSYM"\n", ark_mem->hadapt_mem->safety);
    fprintf(fp, "     Bias factor = %"RSYM"\n", ark_mem->hadapt_mem->bias);
    fprintf(fp, "     Growth factor = %"RSYM"\n", ark_mem->hadapt_mem->growth);
    fprintf(fp, "     Step growth lower bound = %"RSYM"\n", ark_mem->hadapt_mem->lbound);
    fprintf(fp, "     Step growth upper bound = %"RSYM"\n", ark_mem->hadapt_mem->ubound);
    fprintf(fp, "     k1 = %"RSYM"\n", ark_mem->hadapt_mem->k1);
    fprintf(fp, "     k2 = %"RSYM"\n", ark_mem->hadapt_mem->k2);
    fprintf(fp, "     k3 = %"RSYM"\n", ark_mem->hadapt_mem->k3);
    if (ark_mem->hadapt_mem->expstab == arkExpStab) {
      fprintf(fp, "  Default explicit stability function\n");
    } else {
      fprintf(fp, "  User provided explicit stability function\n");
    }
  } else {
    fprintf(fp, "  User provided time step adaptivity function\n");
  }

  fprintf(fp, "  Maximum number of error test failures = %i\n",ark_mem->maxnef);
  fprintf(fp, "  Maximum number of convergence test failures = %i\n",ark_mem->maxncf);

  return(ARK_SUCCESS);
}


/*===============================================================
  ARKODE + XBraid interface utility functions
  ===============================================================*/


/*---------------------------------------------------------------
  arkSetForcePass:

  Ignore the value of kflag after the temporal error test and
  force the step to pass.
  ---------------------------------------------------------------*/
int arkSetForcePass(void *arkode_mem, booleantype force_pass)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetForcePass", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->force_pass = force_pass;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetLastKFlag:

  The last kflag value retured by the temporal error test.
  ---------------------------------------------------------------*/
int arkGetLastKFlag(void *arkode_mem, int *last_kflag)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetLastKFlag", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *last_kflag = ark_mem->last_kflag;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
