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
 * output functions for the ARKode infrastructure; these routines
 * should not be called directly by the user; instead they are
 * provided as utility routines for ARKode time-step modules
 * to use.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
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
  using ARKodeRootInit).
  ---------------------------------------------------------------*/
int arkSetDefaults(ARKodeMem ark_mem)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetDefaults", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* Set default values for integrator optional inputs */
  ark_mem->dense_q          = QDENSE_DEF;     /* dense output order */
  ark_mem->fixedstep        = SUNFALSE;       /* default to use adaptive steps */
  ark_mem->reltol           = 1.e-4;          /* relative tolerance */
  ark_mem->itol             = ARK_SS;         /* scalar-scalar solution tolerances */
  ark_mem->ritol            = ARK_SS;         /* scalar-scalar residual tolerances */
  ark_mem->Sabstol          = 1.e-9;          /* solution absolute tolerance */
  ark_mem->SRabstol         = 1.e-9;          /* residual absolute tolerance */
  ark_mem->user_efun        = SUNFALSE;       /* no user-supplied ewt function */
  ark_mem->efun             = arkEwtSet;      /* built-in ewt function */
  ark_mem->e_data           = NULL;           /* ewt function data */
  ark_mem->user_rfun        = SUNFALSE;       /* no user-supplied rwt function */
  ark_mem->rfun             = arkRwtSet;      /* built-in rwt function */
  ark_mem->r_data           = NULL;           /* rwt function data */
  ark_mem->ehfun            = arkErrHandler;  /* default error handler fn */
  ark_mem->eh_data          = ark_mem;        /* error handler data */
  ark_mem->errfp            = stderr;         /* output stream for errors */
  ark_mem->mxstep           = MXSTEP_DEFAULT; /* max number of steps */
  ark_mem->mxhnil           = MXHNIL;         /* max warns of t+h==t */
  ark_mem->hin              = ZERO;           /* determine initial step on-the-fly */
  ark_mem->hmin             = ZERO;           /* no minimum step size */
  ark_mem->hmax_inv         = ZERO;           /* no maximum step size */
  ark_mem->tstopset         = SUNFALSE;       /* no stop time set */
  ark_mem->tstop            = ZERO;           /* no fixed stop time */
  ark_mem->diagfp           = NULL;           /* no solver diagnostics file */
  ark_mem->report           = SUNFALSE;       /* don't report solver diagnostics */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetDenseOrder:

  Specifies the polynomial order for dense output.  Positive
  values are sent to the interpolation module; negative values
  imply to use the default.
  ---------------------------------------------------------------*/
int arkSetDenseOrder(ARKodeMem ark_mem, int dord)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetDenseOrder", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* set user-provided value, or default, depending on argument */
  if (dord < 0) {
    ark_mem->dense_q = QDENSE_DEF;
  } else {
    ark_mem->dense_q = dord;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetErrHandlerFn:

  Specifies the error handler function
  ---------------------------------------------------------------*/
int arkSetErrHandlerFn(ARKodeMem ark_mem, ARKErrHandlerFn ehfun,
                       void *eh_data)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetErrHandlerFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

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
int arkSetErrFile(ARKodeMem ark_mem, FILE *errfp)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetErrFile", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem->errfp = errfp;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetUserData:

  Specifies the user data pointer for f
  ---------------------------------------------------------------*/
int arkSetUserData(ARKodeMem ark_mem, void *user_data)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetUserData", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem->user_data = user_data;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetDiagnostics:

  Specifies to enable solver diagnostics, and specifies the FILE
  pointer for output (diagfp==NULL disables output)
  ---------------------------------------------------------------*/
int arkSetDiagnostics(ARKodeMem ark_mem, FILE *diagfp)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetDiagnostics", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

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
int arkSetMaxNumSteps(ARKodeMem ark_mem, long int mxsteps)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMaxNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

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
int arkSetMaxHnilWarns(ARKodeMem ark_mem, int mxhnil)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMaxHnilWarns", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

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
int arkSetInitStep(ARKodeMem ark_mem, realtype hin)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* Passing hin=0 sets the default, otherwise use input. */
  if (hin == ZERO) {
    ark_mem->hin = ZERO;
  } else {
    ark_mem->hin = hin;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSetMinStep:

  Specifies the minimum step size
  ---------------------------------------------------------------*/
int arkSetMinStep(ARKodeMem ark_mem, realtype hmin)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMinStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* Passing a value <= 0 sets hmax = infinity */
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
int arkSetMaxStep(ARKodeMem ark_mem, realtype hmax)
{
  realtype hmax_inv;
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetMaxStep", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }

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
int arkSetStopTime(ARKodeMem ark_mem, realtype tstop)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetStopTime", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
  }

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
int arkSetFixedStep(ARKodeMem ark_mem, realtype hfixed)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetFixedStep", MSG_ARK_NO_MEM);
    return (ARK_MEM_NULL);
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
int arkSetRootDirection(ARKodeMem ark_mem, int *rootdir)
{
  ARKodeRootMem ark_root_mem;
  int i;

  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetRootDirection", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
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
int arkSetNoInactiveRootWarn(ARKodeMem ark_mem)
{
  ARKodeRootMem ark_root_mem;
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetNoInactiveRootWarn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
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
  type ARKPostProcessStepFn.  A NULL input function disables step
  postprocessing.

  IF THE SUPPLIED FUNCTION MODIFIES ANY OF THE ACTIVE STATE DATA,
  THEN ALL THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND
  STABILITY ARE LOST.
  ---------------------------------------------------------------*/
int arkSetPostprocessStepFn(ARKodeMem ark_mem,
                            ARKPostProcessStepFn ProcessStep)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSetPostprocessStepFn", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* NULL argument sets default, otherwise set inputs */
  ark_mem->ProcessStep = ProcessStep;
  return(ARK_SUCCESS);
}


/*===============================================================
  ARKode optional output utility functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkGetNumSteps:

  Returns the current number of integration steps
  ---------------------------------------------------------------*/
int arkGetNumSteps(ARKodeMem ark_mem, long int *nsteps)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetNumSteps", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *nsteps = ark_mem->nst;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetActualInitStep:

  Returns the step size used on the first step
  ---------------------------------------------------------------*/
int arkGetActualInitStep(ARKodeMem ark_mem, realtype *hinused)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetActualInitStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *hinused = ark_mem->h0u;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetLastStep:

  Returns the step size used on the last successful step
  ---------------------------------------------------------------*/
int arkGetLastStep(ARKodeMem ark_mem, realtype *hlast)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetLastStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *hlast = ark_mem->hold;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetCurrentStep:

  Returns the step size to be attempted on the next step
  ---------------------------------------------------------------*/
int arkGetCurrentStep(ARKodeMem ark_mem, realtype *hcur)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetCurrentStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *hcur = ark_mem->next_h;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetCurrentTime:

  Returns the current value of the independent variable
  ---------------------------------------------------------------*/
int arkGetCurrentTime(ARKodeMem ark_mem, realtype *tcur)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetCurrentTime", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *tcur = ark_mem->tcur;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetTolScaleFactor:

  Returns a suggested factor for scaling tolerances
  ---------------------------------------------------------------*/
int arkGetTolScaleFactor(ARKodeMem ark_mem, realtype *tolsfact)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetTolScaleFactor", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *tolsfact = ark_mem->tolsf;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetErrWeights:

  This routine returns the current error weight vector.
  ---------------------------------------------------------------*/
int arkGetErrWeights(ARKodeMem ark_mem, N_Vector eweight)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetErrWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  N_VScale(ONE, ark_mem->ewt, eweight);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetResWeights:

  This routine returns the current residual weight vector.
  ---------------------------------------------------------------*/
int arkGetResWeights(ARKodeMem ark_mem, N_Vector rweight)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetResWeights", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  N_VScale(ONE, ark_mem->rwt, rweight);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetWorkSpace:

  Returns integrator work space requirements
  ---------------------------------------------------------------*/
int arkGetWorkSpace(ARKodeMem ark_mem, long int *lenrw, long int *leniw)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetWorkSpace", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *leniw = ark_mem->liw;
  *lenrw = ark_mem->lrw;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkGetNumGEvals:

  Returns the current number of calls to g (for rootfinding)
  ---------------------------------------------------------------*/
int arkGetNumGEvals(ARKodeMem ark_mem, long int *ngevals)
{
  ARKodeRootMem ark_root_mem;
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetNumGEvals", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
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
int arkGetRootInfo(ARKodeMem ark_mem, int *rootsfound)
{
  int i;
  ARKodeRootMem ark_root_mem;
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetRootInfo", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
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
int arkGetStepStats(ARKodeMem ark_mem, long int *nsteps,
                    realtype *hinused, realtype *hlast,
                    realtype *hcur, realtype *tcur)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkGetStepStats", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *nsteps  = ark_mem->nst;
  *hinused = ark_mem->h0u;
  *hlast   = ark_mem->hold;
  *hcur    = ark_mem->next_h;
  *tcur    = ark_mem->tcur;
  return(ARK_SUCCESS);
}


/*-----------------------------------------------------------------*/

char *arkGetReturnFlagName(long int flag)
{
  char *name;
  name = (char *)malloc(24*sizeof(char));

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
  case ARK_POSTPROCESS_FAIL:
    sprintf(name,"ARK_POSTPROCESS_FAIL");
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
  fprintf(fp, "  Dense output order %i\n",ark_mem->dense_q);
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

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
