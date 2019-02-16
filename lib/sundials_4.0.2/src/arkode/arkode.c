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
 * This is the implementation file for the main ARKode
 * infrastructure.  It is independent of the ARKode time step
 * module, nonlinear solver, linear solver and vector modules in
 * use.
 *--------------------------------------------------------------*/

/*===============================================================
  Import Header Files
  ===============================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#define NO_DEBUG_OUTPUT
/* #define DEBUG_OUTPUT */
#ifdef DEBUG_OUTPUT
#include <nvector/nvector_serial.h>
#endif

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif


/*===============================================================
  EXPORTED FUNCTIONS
  ===============================================================*/

/*---------------------------------------------------------------
  arkCreate:

  arkCreate creates an internal memory block for a problem to
  be solved by a time step module built on ARKode.  If successful,
  arkCreate returns a pointer to the problem memory. If an
  initialization error occurs, arkCreate prints an error message
  to standard err and returns NULL.
  ---------------------------------------------------------------*/
ARKodeMem arkCreate()
{
  int iret;
  ARKodeMem ark_mem;

  ark_mem = NULL;
  ark_mem = (ARKodeMem) malloc(sizeof(struct ARKodeMemRec));
  if (ark_mem == NULL) {
    arkProcessError(NULL, 0, "ARKode", "arkCreate",
                    MSG_ARK_ARKMEM_FAIL);
    return(NULL);
  }

  /* Zero out ark_mem */
  memset(ark_mem, 0, sizeof(struct ARKodeMemRec));

  /* Set uround */
  ark_mem->uround = UNIT_ROUNDOFF;

  /* Set default values for integrator optional inputs */
  iret = arkSetDefaults(ark_mem);
  if (iret != ARK_SUCCESS) {
    arkProcessError(NULL, 0, "ARKode", "arkCreate",
                    "Error setting default solver options");
    return(NULL);
  }

  /* Initialize time step module to NULL */
  ark_mem->step_attachlinsol = NULL;
  ark_mem->step_attachmasssol = NULL;
  ark_mem->step_disablelsetup = NULL;
  ark_mem->step_disablemsetup = NULL;
  ark_mem->step_getlinmem = NULL;
  ark_mem->step_getmassmem = NULL;
  ark_mem->step_getimplicitrhs = NULL;
  ark_mem->step_mmult = NULL;
  ark_mem->step_getgammas = NULL;
  ark_mem->step_init = NULL;
  ark_mem->step_fullrhs = NULL;
  ark_mem->step = NULL;
  ark_mem->step_mem = NULL;

  /* Initialize root finding variables */
  ark_mem->root_mem = NULL;

  /* Initialize diagnostics reporting variables */
  ark_mem->report  = SUNFALSE;
  ark_mem->diagfp  = NULL;

  /* Initialize lrw and liw */
  ark_mem->lrw = 18;
  ark_mem->liw = 39;  /* fcn/data ptr, int, long int, sunindextype, booleantype */

  /* No mallocs have been done yet */
  ark_mem->VabstolMallocDone  = SUNFALSE;
  ark_mem->VRabstolMallocDone = SUNFALSE;
  ark_mem->MallocDone         = SUNFALSE;

  /* No user-supplied step postprocessing function yet */
  ark_mem->ProcessStep = NULL;

  /* Return pointer to ARKode memory block */
  return(ark_mem);
}


/*---------------------------------------------------------------
  arkResize:

  arkResize re-initializes ARKode's memory for a problem with a
  changing vector size.  It is assumed that the problem dynamics
  before and after the vector resize will be comparable, so that
  all time-stepping heuristics prior to calling arkResize
  remain valid after the call.  If instead the dynamics should be
  re-calibrated, the ARKode memory structure should be deleted
  with a call to *StepFree, and re-created with a call to
  *StepCreate.

  To aid in the vector-resize operation, the user can supply a
  vector resize function, that will take as input an N_Vector with
  the previous size, and return as output a corresponding vector
  of the new size.  If this function (of type ARKVecResizeFn) is
  not supplied (i.e. is set to NULL), then all existing N_Vectors
  will be destroyed and re-cloned from the input vector.

  In the case that the dynamical time scale should be modified
  slightly from the previous time scale, an input "hscale" is
  allowed, that will re-scale the upcoming time step by the
  specified factor.  If a value <= 0 is specified, the default of
  1.0 will be used.

  Other arguments:
  ark_mem          Existing ARKode memory data structure.
  y0               The newly-sized solution vector, holding
                   the current dependent variable values.
  t0               The current value of the independent
                   variable.
  resize_data      User-supplied data structure that will be
                   passed to the supplied resize function.

  The return value is ARK_SUCCESS = 0 if no errors occurred, or
  a negative value otherwise.
  ---------------------------------------------------------------*/
int arkResize(ARKodeMem ark_mem, N_Vector y0, realtype hscale,
              realtype t0, ARKVecResizeFn resize, void *resize_data)
{
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int ier;

  /* Check ark_mem */
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkResize", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkResize", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkResize", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKode state */
  ark_mem->tcur = t0;
  ark_mem->tn   = t0;

  /* Update time-stepping parameters */
  /*   adjust upcoming step size depending on hscale */
  if (hscale < 0.0)  hscale = 1.0;
  if (hscale != 1.0) {

    /* Encode hscale into ark_mem structure */
    ark_mem->eta = hscale;
    ark_mem->hprime *= hscale;

    /* If next step would overtake tstop, adjust stepsize */
    if ( ark_mem->tstopset )
      if ( (ark_mem->tcur + ark_mem->hprime - ark_mem->tstop)*ark_mem->hprime > ZERO ) {
        ark_mem->hprime = (ark_mem->tstop-ark_mem->tcur) *
          (ONE-FOUR*ark_mem->uround);
        ark_mem->eta = ark_mem->hprime/ark_mem->h;
      }

  }

  /* Determing change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL)
    N_VSpace(y0, &lrw1, &liw1);
  lrw_diff = lrw1 - ark_mem->lrw1;
  liw_diff = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* Resize the ARKode vectors */
  /*     Vabstol */
  ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                     liw_diff, y0, &ark_mem->Vabstol);
  if (ier != ARK_SUCCESS)  return(ier);
  /*     VRabstol */
  ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                     liw_diff, y0, &ark_mem->VRabstol);
  if (ier != ARK_SUCCESS)  return(ier);
  /*     ewt */
  ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                     liw_diff, y0, &ark_mem->ewt);
  if (ier != ARK_SUCCESS)  return(ier);
  /*     rwt  */
  if (ark_mem->rwt_is_ewt) {      /* update pointer to ewt */
    ark_mem->rwt = ark_mem->ewt;
  } else {                            /* resize if distinct from ewt */
    ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &ark_mem->rwt);
    if (ier != ARK_SUCCESS)  return(ier);
  }
  /*     yn */
  ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                     liw_diff, y0, &ark_mem->yn);
  if (ier != ARK_SUCCESS)  return(ier);
  /*     tempv* */
  ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                     liw_diff, y0, &ark_mem->tempv1);
  if (ier != ARK_SUCCESS)  return(ier);
  ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                     liw_diff, y0, &ark_mem->tempv2);
  if (ier != ARK_SUCCESS)  return(ier);
  ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                     liw_diff, y0, &ark_mem->tempv3);
  if (ier != ARK_SUCCESS)  return(ier);
  ier = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                     liw_diff, y0, &ark_mem->tempv4);
  if (ier != ARK_SUCCESS)  return(ier);


  /* Resize interpolation structure memory */
  if (ark_mem->interp) {
    ier = arkInterpResize(ark_mem, ark_mem->interp, resize,
                          resize_data, lrw_diff, liw_diff, y0);
    if (ier != ARK_SUCCESS) {
      arkProcessError(ark_mem, ier, "ARKode", "arkResize",
                      "Interpolation module resize failure");
      return(ier);
    }
  }

  /* Copy y0 into ark_yn to set the current solution */
  N_VScale(ONE, y0, ark_mem->yn);

  /* Indicate that problem size is new */
  ark_mem->resized = SUNTRUE;
  ark_mem->firststage = SUNTRUE;

  /* Problem has been successfully re-sized */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkSStolerances, arkSVtolerances, arkWFtolerances:

  These functions specify the integration tolerances. One of them
  SHOULD be called before the first call to arkEvolve; otherwise
  default values of reltol=1e-4 and abstol=1e-9 will be used,
  which may be entirely incorrect for a specific problem.

  arkSStolerances specifies scalar relative and absolute
  tolerances.

  arkSVtolerances specifies scalar relative tolerance and a
  vector absolute tolerance (a potentially different absolute
  tolerance for each vector component).

  arkWFtolerances specifies a user-provides function (of type
  ARKEwtFn) which will be called to set the error weight vector.
  ---------------------------------------------------------------*/
int arkSStolerances(ARKodeMem ark_mem, realtype reltol, realtype abstol)
{
  /* Check inputs */
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSStolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkSStolerances", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }
  if (reltol < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSStolerances", MSG_ARK_BAD_RELTOL);
    return(ARK_ILL_INPUT);
  }
  if (abstol < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSStolerances", MSG_ARK_BAD_ABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  ark_mem->reltol  = reltol;
  ark_mem->Sabstol = abstol;
  ark_mem->itol    = ARK_SS;

  /* enforce use of arkEwtSet */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->efun      = arkEwtSet;
  ark_mem->e_data    = ark_mem;

  return(ARK_SUCCESS);
}


int arkSVtolerances(ARKodeMem ark_mem, realtype reltol, N_Vector abstol)
{
  /* Check inputs */
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkSVtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkSVtolerances", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }
  if (reltol < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSVtolerances", MSG_ARK_BAD_RELTOL);
    return(ARK_ILL_INPUT);
  }
  if (N_VMin(abstol) < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSVtolerances", MSG_ARK_BAD_ABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  if ( !(ark_mem->VabstolMallocDone) ) {
    ark_mem->Vabstol = N_VClone(ark_mem->ewt);
    ark_mem->lrw += ark_mem->lrw1;
    ark_mem->liw += ark_mem->liw1;
    ark_mem->VabstolMallocDone = SUNTRUE;
  }
  N_VScale(ONE, abstol, ark_mem->Vabstol);
  ark_mem->reltol = reltol;
  ark_mem->itol   = ARK_SV;

  /* enforce use of arkEwtSet */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->efun      = arkEwtSet;
  ark_mem->e_data    = ark_mem;

  return(ARK_SUCCESS);
}


int arkWFtolerances(ARKodeMem ark_mem, ARKEwtFn efun)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkWFtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkWFtolerances", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Copy tolerance data into memory */
  ark_mem->itol      = ARK_WF;
  ark_mem->user_efun = SUNTRUE;
  ark_mem->efun      = efun;
  ark_mem->e_data    = NULL; /* set to user_data in InitialSetup */

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkResStolerance, arkResVtolerance, arkResFtolerance:

  These functions specify the absolute residual tolerance.
  Specification of the absolute residual tolerance is only
  necessary for problems with non-identity mass matrices in which
  the units of the solution vector y dramatically differ from the
  units of the ODE right-hand side f(t,y).  If this occurs, one
  of these routines SHOULD be called before the first call to
  ARKode; otherwise the default value of rabstol=1e-9 will be
  used, which may be entirely incorrect for a specific problem.

  arkResStolerances specifies a scalar residual tolerance.

  arkResVtolerances specifies a vector residual tolerance
  (a potentially different absolute residual tolerance for
  each vector component).

  arkResFtolerances specifies a user-provides function (of
  type ARKRwtFn) which will be called to set the residual
  weight vector.
  ---------------------------------------------------------------*/
int arkResStolerance(ARKodeMem ark_mem, realtype rabstol)
{
  /* Check inputs */
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkResStolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkResStolerances", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }
  if (rabstol < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkResStolerances", MSG_ARK_BAD_RABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Allocate space for rwt if necessary */
  if (ark_mem->rwt_is_ewt) {
    ark_mem->rwt_is_ewt = SUNFALSE;
    ark_mem->rwt = N_VClone(ark_mem->ewt);
    ark_mem->lrw += ark_mem->lrw1;
    ark_mem->liw += ark_mem->liw1;
  }

  /* Copy tolerances into memory */
  ark_mem->SRabstol = rabstol;
  ark_mem->ritol    = ARK_SS;

  /* enforce use of arkRwtSet */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->rfun      = arkRwtSet;
  ark_mem->r_data    = ark_mem;

  return(ARK_SUCCESS);
}


int arkResVtolerance(ARKodeMem ark_mem, N_Vector rabstol)
{
  /* Check inputs */
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkResVtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkResVtolerances", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }
  if (N_VMin(rabstol) < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkResVtolerances", MSG_ARK_BAD_RABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Allocate space for rwt if necessary */
  if (ark_mem->rwt_is_ewt) {
    ark_mem->rwt_is_ewt = SUNFALSE;
    ark_mem->rwt = N_VClone(ark_mem->ewt);
    ark_mem->lrw += ark_mem->lrw1;
    ark_mem->liw += ark_mem->liw1;
  }

  /* Copy tolerances into memory */
  if ( !(ark_mem->VRabstolMallocDone) ) {
    ark_mem->VRabstol = N_VClone(ark_mem->rwt);
    ark_mem->lrw += ark_mem->lrw1;
    ark_mem->liw += ark_mem->liw1;
    ark_mem->VRabstolMallocDone = SUNTRUE;
  }
  N_VScale(ONE, rabstol, ark_mem->VRabstol);
  ark_mem->ritol = ARK_SV;


  /* enforce use of arkRwtSet */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->rfun      = arkRwtSet;
  ark_mem->r_data    = ark_mem;

  return(ARK_SUCCESS);
}


int arkResFtolerance(ARKodeMem ark_mem, ARKRwtFn rfun)
{
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkResFtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkResFtolerances", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Allocate space for rwt if necessary */
  if (ark_mem->rwt_is_ewt) {
    ark_mem->rwt_is_ewt = SUNFALSE;
    ark_mem->rwt = N_VClone(ark_mem->ewt);
    ark_mem->lrw += ark_mem->lrw1;
    ark_mem->liw += ark_mem->liw1;
  }

  /* Copy tolerance data into memory */
  ark_mem->ritol     = ARK_WF;
  ark_mem->user_rfun = SUNTRUE;
  ark_mem->rfun      = rfun;
  ark_mem->r_data    = NULL; /* set to user_data in InitialSetup */

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkEvolve:

  This routine is the main driver of ARKode-based integrators.

  It integrates over a time interval defined by the user, by
  calling the time step module to do internal time steps.

  The first time that arkEvolve is called for a successfully
  initialized problem, it computes a tentative initial step size.

  arkEvolve supports two modes as specified by itask: ARK_NORMAL and
  ARK_ONE_STEP.  In the ARK_NORMAL mode, the solver steps until
  it reaches or passes tout and then interpolates to obtain
  y(tout).  In the ARK_ONE_STEP mode, it takes one internal step
  and returns.  The behavior of both modes can be over-rided
  through user-specification of ark_tstop (through the
  *StepSetStopTime function), in which case if a solver step
  would pass tstop, the step is shortened so that it stops at
  exactly the specified stop time, and hence interpolation of
  y(tout) is not required.
  ---------------------------------------------------------------*/
int arkEvolve(ARKodeMem ark_mem, realtype tout, N_Vector yout,
              realtype *tret, int itask)
{
  long int nstloc;
  int retval, kflag, istate, ir, ier;
  int ewtsetOK;
  realtype troundoff, nrm;
  booleantype inactive_roots;


  /* Check and process inputs */

  /* Check if ark_mem exists */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode", "arkEvolve",
                    MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode", "arkEvolve",
                    MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check for yout != NULL */
  if ((ark_mem->ycur = yout) == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkEvolve",
                    MSG_ARK_YOUT_NULL);
    return(ARK_ILL_INPUT);
  }

  /* Check for tret != NULL */
  if (tret == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkEvolve",
                    MSG_ARK_TRET_NULL);
    return(ARK_ILL_INPUT);
  }

  /* Check for valid itask */
  if ( (itask != ARK_NORMAL) && (itask != ARK_ONE_STEP) ) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkEvolve",
                    MSG_ARK_BAD_ITASK);
    return(ARK_ILL_INPUT);
  }

  /* store copy of itask if using root-finding */
  if (ark_mem->root_mem != NULL) {
    if (itask == ARK_NORMAL) ark_mem->root_mem->toutc = tout;
    ark_mem->root_mem->taskc = itask;
  }


  /* perform first-step-specific initializations:
     - initialize tret values to initialization time
     - perform initial integrator setup  */
  if (ark_mem->nst == 0) {
    ark_mem->tretlast = *tret = ark_mem->tcur;
    ier = arkInitialSetup(ark_mem, tout);
    if (ier!= ARK_SUCCESS) return(ier);
  }


  /* perform first-step-after-resize initializations */
  if (ark_mem->nst > 0 && ark_mem->resized) {
    ier = arkPostResizeSetup(ark_mem);
    if (ier!= ARK_SUCCESS) return(ier);
  }


  /* perform stopping tests */
  if (ark_mem->nst > 0 && !ark_mem->resized)
    if (arkStopTests(ark_mem, tout, yout, tret, itask, &ier))
      return(ier);


  /*--------------------------------------------------
    Looping point for internal steps

    - update the ewt vector for the next step
    - check for errors (too many steps, too much
      accuracy requested, step size too small)
    - take a new step (via time stepper); stop on error
    - perform stop tests:
    - check for root in last step taken
    - check if tout was passed
    - check if close to tstop
    - check if in ONE_STEP mode (must return)
    --------------------------------------------------*/
  nstloc = 0;
  for(;;) {

    ark_mem->next_h = ark_mem->h;

    /* Reset and check ewt */
    if (ark_mem->nst > 0 && !ark_mem->resized) {
      ewtsetOK = ark_mem->efun(ark_mem->yn,
                               ark_mem->ewt,
                               ark_mem->e_data);
      if (ewtsetOK != 0) {
        if (ark_mem->itol == ARK_WF)
          arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkEvolve",
                          MSG_ARK_EWT_NOW_FAIL, ark_mem->tcur);
        else
          arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkEvolve",
                          MSG_ARK_EWT_NOW_BAD, ark_mem->tcur);

        istate = ARK_ILL_INPUT;
        ark_mem->tretlast = *tret = ark_mem->tcur;
        N_VScale(ONE, ark_mem->yn, yout);
        break;
      }
    }

    /* Reset and check rwt */
    if (!ark_mem->rwt_is_ewt) {
      if (ark_mem->nst > 0 && !ark_mem->resized) {
        ewtsetOK = ark_mem->rfun(ark_mem->yn,
                                 ark_mem->rwt,
                                 ark_mem->r_data);
        if (ewtsetOK != 0) {
          if (ark_mem->itol == ARK_WF)
            arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkEvolve",
                            MSG_ARK_RWT_NOW_FAIL, ark_mem->tcur);
          else
            arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkEvolve",
                            MSG_ARK_RWT_NOW_BAD, ark_mem->tcur);

          istate = ARK_ILL_INPUT;
          ark_mem->tretlast = *tret = ark_mem->tcur;
          N_VScale(ONE, ark_mem->yn, yout);
          break;
        }
      }
    }

    /* Check for too many steps */
    if ( (ark_mem->mxstep>0) && (nstloc >= ark_mem->mxstep) ) {
      arkProcessError(ark_mem, ARK_TOO_MUCH_WORK, "ARKode", "arkEvolve",
                      MSG_ARK_MAX_STEPS, ark_mem->tcur);
      istate = ARK_TOO_MUCH_WORK;
      ark_mem->tretlast = *tret = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, yout);
      break;
    }

    /* Check for too much accuracy requested */
    nrm = N_VWrmsNorm(ark_mem->yn, ark_mem->ewt);
    ark_mem->tolsf = ark_mem->uround * nrm;
    if (ark_mem->tolsf > ONE) {
      arkProcessError(ark_mem, ARK_TOO_MUCH_ACC, "ARKode", "arkEvolve",
                      MSG_ARK_TOO_MUCH_ACC, ark_mem->tcur);
      istate = ARK_TOO_MUCH_ACC;
      ark_mem->tretlast = *tret = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, yout);
      ark_mem->tolsf *= TWO;
      break;
    } else {
      ark_mem->tolsf = ONE;
    }

    /* Check for h below roundoff level in tn */
    if (ark_mem->tcur + ark_mem->h == ark_mem->tcur) {
      ark_mem->nhnil++;
      if (ark_mem->nhnil <= ark_mem->mxhnil)
        arkProcessError(ark_mem, ARK_WARNING, "ARKode", "arkEvolve",
                        MSG_ARK_HNIL, ark_mem->tcur, ark_mem->h);
      if (ark_mem->nhnil == ark_mem->mxhnil)
        arkProcessError(ark_mem, ARK_WARNING, "ARKode", "arkEvolve",
                        MSG_ARK_HNIL_DONE);
    }

    /* Update parameter for upcoming step size */
    if ((ark_mem->nst > 0) && (ark_mem->hprime != ark_mem->h)) {
      ark_mem->h = ark_mem->h * ark_mem->eta;
      ark_mem->next_h = ark_mem->h;
    }
    if (ark_mem->fixedstep) {
      ark_mem->h = ark_mem->hin;
      ark_mem->next_h = ark_mem->h;
    }

    /* Call time stepper module to take a step */
    kflag = ark_mem->step((void*) ark_mem);

    /* Process successful step, catch additional errors to send to arkHandleFailure */
    if (kflag == ARK_SUCCESS)
      kflag = arkCompleteStep(ark_mem);

    /* Process failed step cases, and exit loop */
    if (kflag != ARK_SUCCESS) {
      istate = arkHandleFailure(ark_mem, kflag);
      ark_mem->tretlast = *tret = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, yout);
      break;
    }

    nstloc++;

    /* Check for root in last step taken. */
    if (ark_mem->root_mem != NULL)
      if (ark_mem->root_mem->nrtfn > 0) {

        retval = arkRootCheck3((void*) ark_mem);
        if (retval == RTFOUND) {  /* A new root was found */
          ark_mem->root_mem->irfnd = 1;
          istate = ARK_ROOT_RETURN;
          ark_mem->tretlast = *tret = ark_mem->root_mem->tlo;
          break;
        } else if (retval == ARK_RTFUNC_FAIL) { /* g failed */
          arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKode", "arkEvolve",
                          MSG_ARK_RTFUNC_FAILED, ark_mem->root_mem->tlo);
          istate = ARK_RTFUNC_FAIL;
          break;
        }

        /* If we are at the end of the first step and we still have
           some event functions that are inactive, issue a warning
           as this may indicate a user error in the implementation
           of the root function. */
        if (ark_mem->nst==1) {
          inactive_roots = SUNFALSE;
          for (ir=0; ir<ark_mem->root_mem->nrtfn; ir++) {
            if (!ark_mem->root_mem->gactive[ir]) {
              inactive_roots = SUNTRUE;
              break;
            }
          }
          if ((ark_mem->root_mem->mxgnull > 0) && inactive_roots) {
            arkProcessError(ark_mem, ARK_WARNING, "ARKode", "arkEvolve",
                            MSG_ARK_INACTIVE_ROOTS);
          }
        }
      }

    /* In NORMAL mode, check if tout reached */
    if ( (itask == ARK_NORMAL) &&
         (ark_mem->tcur-tout)*ark_mem->h >= ZERO ) {
      istate = ARK_SUCCESS;
      ark_mem->tretlast = *tret = tout;
      (void) arkGetDky(ark_mem, tout, 0, yout);
      ark_mem->next_h = ark_mem->hprime;
      break;
    }

    /* Check if tn is at tstop or near tstop */
    if ( ark_mem->tstopset ) {
      troundoff = FUZZ_FACTOR*ark_mem->uround *
        (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->h));
      if ( SUNRabs(ark_mem->tcur - ark_mem->tstop) <= troundoff) {
        (void) arkGetDky(ark_mem, ark_mem->tstop, 0, yout);
        ark_mem->tretlast = *tret = ark_mem->tstop;
        ark_mem->tstopset = SUNFALSE;
        istate = ARK_TSTOP_RETURN;
        break;
      }
      if ( (ark_mem->tcur + ark_mem->hprime - ark_mem->tstop)*ark_mem->h > ZERO ) {
        ark_mem->hprime = (ark_mem->tstop - ark_mem->tcur) *
          (ONE-FOUR*ark_mem->uround);
        ark_mem->eta = ark_mem->hprime/ark_mem->h;
      }
    }

    /* In ONE_STEP mode, copy y and exit loop */
    if (itask == ARK_ONE_STEP) {
      istate = ARK_SUCCESS;
      ark_mem->tretlast = *tret = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, yout);
      ark_mem->next_h = ark_mem->hprime;
      break;
    }

  } /* end looping for internal steps */

  return(istate);
}


/*---------------------------------------------------------------
  arkGetDky:

  This routine computes the k-th derivative of the interpolating
  polynomial at the time t and stores the result in the vector
  dky. This routine internally calls arkInterpEvaluate to perform the
  interpolation.  We have the restriction that 0 <= k <= 3.  This
  routine uses an interpolating polynomial of degree
  max(ark_dense_q, k), i.e. it will form a polynomial of the
  degree requested by the user through ark_dense_q, unless
  higher-order derivatives are requested.

  This function is called by arkEvolve with k=0 and t=tout to perform
  interpolation of outputs, but may also be called indirectly by the
  user via time step module *StepGetDky calls.  Note: in all cases
  it will be called after ark_tcur has been updated to correspond
  with the end time of the last successful step.
  ---------------------------------------------------------------*/
int arkGetDky(ARKodeMem ark_mem, realtype t, int k, N_Vector dky)
{
  realtype s, tfuzz, tp, tn1;
  int retval;

  /* Check all inputs for legality */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode", "arkGetDky",
                    MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  if (dky == NULL) {
    arkProcessError(ark_mem, ARK_BAD_DKY, "ARKode", "arkGetDky",
                    MSG_ARK_NULL_DKY);
    return(ARK_BAD_DKY);
  }
  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode", "arkGetDky",
                    "Missing interpolation structure");
    return(ARK_MEM_NULL);
  }


  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * ark_mem->uround *
    (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->hold));
  if (ark_mem->hold < ZERO) tfuzz = -tfuzz;
  tp = ark_mem->tcur - ark_mem->hold - tfuzz;
  tn1 = ark_mem->tcur + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
    arkProcessError(ark_mem, ARK_BAD_T, "ARKode", "arkGetDky",
                    MSG_ARK_BAD_T, t, ark_mem->tcur-ark_mem->hold,
                    ark_mem->tcur);
    return(ARK_BAD_T);
  }

  /* call arkInterpEvaluate to evaluate result */
  s = (t - ark_mem->tcur) / ark_mem->h;
  retval = arkInterpEvaluate(ark_mem, ark_mem->interp, s,
                             k, ark_mem->dense_q, dky);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode", "arkGetDky",
                    "Error calling arkInterpEvaluate");
    return(retval);
  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkFree:

  This routine frees the ARKode infrastructure memory.
  ---------------------------------------------------------------*/
void arkFree(void **arkode_mem)
{
  ARKodeMem ark_mem;

  if (*arkode_mem == NULL) return;

  ark_mem = (ARKodeMem) (*arkode_mem);

  arkFreeVectors(ark_mem);
  if (ark_mem->interp != NULL)
    arkInterpFree(&(ark_mem->interp));

  if (ark_mem->root_mem != NULL)
    (void) arkRootFree(*arkode_mem);

  free(*arkode_mem);
  *arkode_mem = NULL;
}



/*===============================================================
  Internal functions that may be replaced by the user
  ===============================================================*/

/*---------------------------------------------------------------
  arkEwtSet

  This routine is responsible for setting the error weight vector ewt,
  according to tol_type, as follows:

  (1) ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol), i=0,...,neq-1
      if tol_type = ARK_SS
  (2) ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol[i]), i=0,...,neq-1
      if tol_type = ARK_SV

  arkEwtSet returns 0 if ewt is successfully set as above to a
  positive vector and -1 otherwise. In the latter case, ewt is
  considered undefined.

  All the real work is done in the routines arkEwtSetSS, arkEwtSetSV.
  ---------------------------------------------------------------*/
int arkEwtSet(N_Vector ycur, N_Vector weight, void *data)
{
  ARKodeMem ark_mem;
  int flag = 0;

  /* data points to ark_mem here */
  ark_mem = (ARKodeMem) data;

  switch(ark_mem->itol) {
  case ARK_SS:
    flag = arkEwtSetSS(ark_mem, ycur, weight);
    break;
  case ARK_SV:
    flag = arkEwtSetSV(ark_mem, ycur, weight);
    break;
  }

  return(flag);
}


/*---------------------------------------------------------------
  arkRwtSet

  This routine is responsible for setting the residual weight
  vector rwt, according to tol_type, as follows:

  (1) rwt[i] = 1 / (reltol * SUNRabs(M*ycur[i]) + rabstol), i=0,...,neq-1
      if tol_type = ARK_SS
  (2) rwt[i] = 1 / (reltol * SUNRabs(M*ycur[i]) + rabstol[i]), i=0,...,neq-1
      if tol_type = ARK_SV
  (3) unset if tol_type is any other value (occurs rwt=ewt)

  arkRwtSet returns 0 if rwt is successfully set as above to a
  positive vector and -1 otherwise. In the latter case, rwt is
  considered undefined.

  All the real work is done in the routines arkRwtSetSS, arkRwtSetSV.
  ---------------------------------------------------------------*/
int arkRwtSet(N_Vector y, N_Vector weight, void *data)
{
  ARKodeMem ark_mem;
  N_Vector My;
  int flag = 0;

  /* data points to ark_mem here */
  ark_mem = (ARKodeMem) data;

  /* return if rwt is just ewt */
  if (ark_mem->rwt_is_ewt)  return(0);

  /* put M*y into ark_tempv1 */
  My = ark_mem->tempv1;
  if (ark_mem->step_mmult != NULL) {
    flag = ark_mem->step_mmult((void *) ark_mem, y, My);
    if (flag != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
  } else {  /* this condition should not apply, but just in case */
    N_VScale(ONE, y, My);
  }

  /* call appropriate routine to fill rwt */
  switch(ark_mem->ritol) {
  case ARK_SS:
    flag = arkRwtSetSS(ark_mem, My, weight);
    break;
  case ARK_SV:
    flag = arkRwtSetSV(ark_mem, My, weight);
    break;
  }

  return(flag);
}


/*---------------------------------------------------------------
  arkErrHandler is the default error handling function.
  It sends the error message to the stream pointed to by ark_errfp
  ---------------------------------------------------------------*/
void arkErrHandler(int error_code, const char *module,
                   const char *function, char *msg, void *data)
{
  ARKodeMem ark_mem;
  char err_type[10];

  /* data points to ark_mem here */
  ark_mem = (ARKodeMem) data;

  if (error_code == ARK_WARNING)
    sprintf(err_type,"WARNING");
  else
    sprintf(err_type,"ERROR");

#ifndef NO_FPRINTF_OUTPUT
  if (ark_mem->errfp!=NULL) {
    fprintf(ark_mem->errfp,"\n[%s %s]  %s\n",module,err_type,function);
    fprintf(ark_mem->errfp,"  %s\n\n",msg);
  }
#endif

  return;
}



/*===============================================================
  Private Helper Functions
  ===============================================================*/

/*---------------------------------------------------------------
  arkInit:

  arkInit allocates and initializes memory for a problem. All
  inputs are checked for errors. If any error occurs during
  initialization, it is reported to the file whose file pointer
  is errfp and an error flag is returned. Otherwise, it returns
  ARK_SUCCESS.  This routine should be called by an ARKode
  timestepper module (not by the user).  This routine must be
  called prior to calling arkEvolve to evolve the problem.
  ---------------------------------------------------------------*/
int arkInit(ARKodeMem ark_mem, realtype t0, N_Vector y0)
{
  booleantype stepperOK, nvectorOK, allocOK;
  sunindextype lrw1, liw1;

  /* Check for legal input parameters */
  if (y0==NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Test if all required time stepper operations are implemented */
  stepperOK = arkCheckTimestepper(ark_mem);
  if (!stepperOK) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInit",
                    "Time stepper module is missing required functionality");
    return(ARK_ILL_INPUT);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = arkCheckNvector(y0);
  if (!nvectorOK) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInit", MSG_ARK_BAD_NVECTOR);
    return(ARK_ILL_INPUT);
  }

  /* Set space requirements for one N_Vector */
  if (y0->ops->nvspace != NULL) {
    N_VSpace(y0, &lrw1, &liw1);
  } else {
    lrw1 = 0;
    liw1 = 0;
  }
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;


  /* Allocate the solver vectors (using y0 as a template) */
  allocOK = arkAllocVectors(ark_mem, y0);
  if (!allocOK) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode",
                    "arkInit", MSG_ARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  /* Initialize the interpolation structure to NULL */
  ark_mem->interp = NULL;

  /* All error checking is complete at this point */

  /* Copy the input parameters into ARKode state */
  ark_mem->tcur = t0;
  ark_mem->tn   = t0;

  /* Set step parameters */
  ark_mem->hold     = ZERO;
  ark_mem->tolsf    = ONE;
  ark_mem->hmin     = ZERO;       /* no minimum step size */
  ark_mem->hmax_inv = ZERO;       /* no maximum step size */

  /* Initialize yn */
  N_VScale(ONE, y0, ark_mem->yn);

  /* Initialize all the counters */
  ark_mem->nst   = 0;
  ark_mem->nhnil = 0;

  /* Initialize other integrator optional outputs */
  ark_mem->h0u    = ZERO;
  ark_mem->next_h = ZERO;

  /* Initially, rwt should point to ewt */
  ark_mem->rwt_is_ewt = SUNTRUE;

  /* Indicate that problem size is new */
  ark_mem->resized    = SUNTRUE;
  ark_mem->firststage = SUNTRUE;

  /* Problem has been successfully initialized */
  ark_mem->MallocDone = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkReInit:

  arkReInit re-initializes ARKode's memory for a problem,
  assuming it has already been allocated in a prior arkInit
  call.  All problem specification inputs are checked for errors.
  If any error occurs during initialization, it is reported to
  the file whose file pointer is errfp.  This routine should only
  be called after arkInit, and only when the problem dynamics
  or desired solvers have changed dramatically, so that the
  problem integration should resume as if started from scratch.

  The return value is ARK_SUCCESS = 0 if no errors occurred, or
  a negative value otherwise.
  ---------------------------------------------------------------*/
int arkReInit(ARKodeMem ark_mem, realtype t0, N_Vector y0)
{
  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkReInit", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Copy the input parameters into ARKode state */
  ark_mem->tcur = t0;
  ark_mem->tn   = t0;

  /* Set step parameters */
  ark_mem->hold     = ZERO;
  ark_mem->tolsf    = ONE;
  ark_mem->hmin     = ZERO;       /* no minimum step size */
  ark_mem->hmax_inv = ZERO;       /* no maximum step size */

  /* Do not reset the linear solver addresses to NULL.  This means
     that if the user does not re-set these manually, we'll re-use
     the linear solver routines that were set during arkInit. */

  /* Initialize yn */
  N_VScale(ONE, y0, ark_mem->yn);

  /* Initialize all the counters */
  ark_mem->nst   = 0;
  ark_mem->nhnil = 0;

  /* Indicate that problem size is new */
  ark_mem->resized    = SUNTRUE;
  ark_mem->firststage = SUNTRUE;

  /* Initialize other integrator optional outputs */
  ark_mem->h0u    = ZERO;
  ark_mem->next_h = ZERO;

  /* Problem has been successfully re-initialized */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkPrintMem:

  This routine outputs the ark_mem structure to a specified file
  pointer.
  ---------------------------------------------------------------*/
void arkPrintMem(ARKodeMem ark_mem, FILE *outfile)
{
  /* output general values */
  fprintf(outfile, "ark_itol = %i\n", ark_mem->itol);
  fprintf(outfile, "ark_ritol = %i\n", ark_mem->ritol);
  fprintf(outfile, "ark_dense_q = %i\n", ark_mem->dense_q);
  fprintf(outfile, "ark_mxhnil = %i\n", ark_mem->mxhnil);
  fprintf(outfile, "ark_mxstep = %li\n", ark_mem->mxstep);
  fprintf(outfile, "ark_lrw1 = %li\n", (long int) ark_mem->lrw1);
  fprintf(outfile, "ark_liw1 = %li\n", (long int) ark_mem->liw1);
  fprintf(outfile, "ark_lrw = %li\n", (long int) ark_mem->lrw);
  fprintf(outfile, "ark_liw = %li\n", (long int) ark_mem->liw);
  fprintf(outfile, "ark_user_efun = %i\n", ark_mem->user_efun);
  fprintf(outfile, "ark_tstopset = %i\n", ark_mem->tstopset);
  fprintf(outfile, "ark_tstop = %" RSYM"\n", ark_mem->tstop);
  fprintf(outfile, "ark_report = %i\n", ark_mem->report);
  fprintf(outfile, "ark_VabstolMallocDone = %i\n", ark_mem->VabstolMallocDone);
  fprintf(outfile, "ark_MallocDone = %i\n", ark_mem->MallocDone);
  fprintf(outfile, "ark_resized = %i\n", ark_mem->resized);
  fprintf(outfile, "ark_firststage = %i\n", ark_mem->firststage);
  fprintf(outfile, "ark_uround = %" RSYM"\n", ark_mem->uround);
  fprintf(outfile, "ark_reltol = %" RSYM"\n", ark_mem->reltol);
  fprintf(outfile, "ark_Sabstol = %" RSYM"\n", ark_mem->Sabstol);
  fprintf(outfile, "ark_fixedstep = %i\n", ark_mem->fixedstep);
  fprintf(outfile, "ark_tolsf = %" RSYM"\n", ark_mem->tolsf);

  /* output counters */
  fprintf(outfile, "ark_nhnil = %i\n", ark_mem->nhnil);
  fprintf(outfile, "ark_nst = %li\n", ark_mem->nst);

  /* output time-stepping values */
  fprintf(outfile, "ark_hin = %" RSYM"\n", ark_mem->hin);
  fprintf(outfile, "ark_h = %" RSYM"\n", ark_mem->h);
  fprintf(outfile, "ark_hprime = %" RSYM"\n", ark_mem->hprime);
  fprintf(outfile, "ark_next_h = %" RSYM"\n", ark_mem->next_h);
  fprintf(outfile, "ark_eta = %" RSYM"\n", ark_mem->eta);
  fprintf(outfile, "ark_tcur = %" RSYM"\n", ark_mem->tcur);
  fprintf(outfile, "ark_tretlast = %" RSYM"\n", ark_mem->tretlast);
  fprintf(outfile, "ark_hmin = %" RSYM"\n", ark_mem->hmin);
  fprintf(outfile, "ark_hmax_inv = %" RSYM"\n", ark_mem->hmax_inv);
  fprintf(outfile, "ark_h0u = %" RSYM"\n", ark_mem->h0u);
  fprintf(outfile, "ark_tn = %" RSYM"\n", ark_mem->tn);
  fprintf(outfile, "ark_hold = %" RSYM"\n", ark_mem->hold);

  /* output root-finding quantities */
  if (ark_mem->root_mem != NULL)
    (void) arkPrintRootMem((void*) ark_mem, outfile);

  /* output interpolation quantities */
  if (ark_mem->interp != NULL)
    arkPrintInterpMem(ark_mem->interp, outfile);

#ifdef DEBUG_OUTPUT
  /* output vector quantities */
  if (ark_mem->Vabstol != NULL) {
    fprintf(outfile, "ark_Vapbsol:\n");
    N_VPrint_Serial(ark_mem->Vabstol);
  }
  if (ark_mem->ewt != NULL) {
    fprintf(outfile, "ark_ewt:\n");
    N_VPrint_Serial(ark_mem->ewt);
  }
  if (!ark_mem->rwt_is_ewt && ark_mem->rwt != NULL) {
    fprintf(outfile, "ark_rwt:\n");
    N_VPrint_Serial(ark_mem->rwt);
  }
  if (ark_mem->ycur != NULL) {
    fprintf(outfile, "ark_ycur:\n");
    N_VPrint_Serial(ark_mem->ycur);
  }
  if (ark_mem->yn != NULL) {
    fprintf(outfile, "ark_yn:\n");
    N_VPrint_Serial(ark_mem->yn);
  }
  if (ark_mem->tempv1 != NULL) {
    fprintf(outfile, "ark_tempv1:\n");
    N_VPrint_Serial(ark_mem->tempv1);
  }
  if (ark_mem->tempv2 != NULL) {
    fprintf(outfile, "ark_tempv2:\n");
    N_VPrint_Serial(ark_mem->tempv2);
  }
  if (ark_mem->tempv3 != NULL) {
    fprintf(outfile, "ark_tempv3:\n");
    N_VPrint_Serial(ark_mem->tempv3);
  }
  if (ark_mem->tempv4 != NULL) {
    fprintf(outfile, "ark_tempv4:\n");
    N_VPrint_Serial(ark_mem->tempv4);
  }
#endif

}


/*---------------------------------------------------------------
  arkCheckTimestepper:

  This routine checks if all required time stepper function
  pointers have been supplied.  If any of them is missing it
  returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype arkCheckTimestepper(ARKodeMem ark_mem)
{
  if ( (ark_mem->step_init == NULL) ||
       (ark_mem->step      == NULL) ||
       (ark_mem->step_mem  == NULL) )
    return(SUNFALSE);
  if ( (ark_mem->interp != NULL) &&
       (ark_mem->step_fullrhs == NULL) )
    return(SUNFALSE);
  return(SUNTRUE);
}


/*---------------------------------------------------------------
  arkCheckNvector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype arkCheckNvector(N_Vector tmpl)  /* to be updated?? */
{
  if ((tmpl->ops->nvclone     == NULL) ||
      (tmpl->ops->nvdestroy   == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) ||
      (tmpl->ops->nvconst     == NULL) ||
      (tmpl->ops->nvdiv       == NULL) ||
      (tmpl->ops->nvscale     == NULL) ||
      (tmpl->ops->nvabs       == NULL) ||
      (tmpl->ops->nvinv       == NULL) ||
      (tmpl->ops->nvaddconst  == NULL) ||
      (tmpl->ops->nvmaxnorm   == NULL) ||
      (tmpl->ops->nvwrmsnorm  == NULL) ||
      (tmpl->ops->nvmin       == NULL))
    return(SUNFALSE);
  else
    return(SUNTRUE);
}


/*---------------------------------------------------------------
  arkAllocVec:

  This routine allocates a single vector based on a template
  vector.  If the target vector already exists it is left alone;
  otherwise it is allocated by cloning the input vector. If the
  allocation is successful (or if the target vector already
  exists) then this returns SUNTRUE.  This routine also updates
  the optional outputs lrw and liw, which are (respectively) the
  lengths of the overall ARKode real and integer work spaces.
  ---------------------------------------------------------------*/
booleantype arkAllocVec(ARKodeMem ark_mem,
                        N_Vector tmpl,
                        N_Vector *v)
{
  if (*v == NULL) {
    *v = N_VClone(tmpl);
    if (*v == NULL) {
      arkFreeVectors(ark_mem);
      return(SUNFALSE);
    } else {
      ark_mem->lrw += ark_mem->lrw1;
      ark_mem->liw += ark_mem->liw1;
    }
  }
  return (SUNTRUE);
}


/*---------------------------------------------------------------
  arkFreeVec:

  This routine frees a single vector.  If the target vector is
  already NULL it is left alone; otherwise it is freed and the
  optional outputs lrw and liw are updated accordingly.
  ---------------------------------------------------------------*/
void arkFreeVec(ARKodeMem ark_mem, N_Vector *v)
{
  if (*v != NULL) {
    N_VDestroy(*v);
    *v = NULL;
    ark_mem->lrw -= ark_mem->lrw1;
    ark_mem->liw -= ark_mem->liw1;
  }
}


/*---------------------------------------------------------------
  arkResizeVec:

  This routine resizes a single vector based on a template
  vector.  If the ARKVecResizeFn function is non-NULL, then it
  calls that routine to perform the single-vector resize;
  otherwise it deallocates and reallocates the target vector based
  on the template vector.  If the resize is successful then this
  returns SUNTRUE.  This routine also updates the optional outputs
  lrw and liw, which are (respectively) the lengths of the overall
  ARKode real and integer work spaces.
  ---------------------------------------------------------------*/
int arkResizeVec(ARKodeMem ark_mem, ARKVecResizeFn resize,
                 void *resize_data, sunindextype lrw_diff,
                 sunindextype liw_diff, N_Vector tmpl, N_Vector *v)
{
  if (*v != NULL) {
    if (resize == NULL) {
      N_VDestroy(*v);
      *v = N_VClone(tmpl);
    } else {
      if (resize(*v, tmpl, resize_data)) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                        "arkResizeVec", MSG_ARK_RESIZE_FAIL);
        return(ARK_ILL_INPUT);
      }
    }
    ark_mem->lrw += lrw_diff;
    ark_mem->liw += liw_diff;
  }
  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  arkAllocVectors:

  This routine allocates the ARKode vectors ewt, yn, tempv* and
  ftemp.  If any of these vectors already exist, they are left
  alone.  Otherwise, it will allocate each vector by cloning the
  input vector. If all memory allocations are successful,
  arkAllocVectors returns SUNTRUE. Otherwise all vector memory
  is freed and arkAllocVectors returns SUNFALSE.  This routine
  also updates the optional outputs lrw and liw, which are
  (respectively) the lengths of the real and integer work spaces.
  ---------------------------------------------------------------*/
booleantype arkAllocVectors(ARKodeMem ark_mem, N_Vector tmpl)
{
  /* Allocate ewt if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->ewt))
    return(SUNFALSE);

  /* Set rwt to point at ewt */
  if (ark_mem->rwt_is_ewt)
    ark_mem->rwt = ark_mem->ewt;

  /* Allocate yn if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->yn))
    return(SUNFALSE);

  /* Allocate tempv1 if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->tempv1))
    return(SUNFALSE);

  /* Allocate tempv2 if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->tempv2))
    return(SUNFALSE);

  /* Allocate tempv3 if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->tempv3))
    return(SUNFALSE);

  /* Allocate tempv4 if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->tempv4))
    return(SUNFALSE);

  return(SUNTRUE);
}


/*---------------------------------------------------------------
  arkFreeVectors

  This routine frees the ARKode vectors allocated in both
  arkAllocVectors and arkAllocRKVectors.
  ---------------------------------------------------------------*/
void arkFreeVectors(ARKodeMem ark_mem)
{
  arkFreeVec(ark_mem, &ark_mem->ewt);
  if (!ark_mem->rwt_is_ewt)
    arkFreeVec(ark_mem, &ark_mem->rwt);
  arkFreeVec(ark_mem, &ark_mem->tempv1);
  arkFreeVec(ark_mem, &ark_mem->tempv2);
  arkFreeVec(ark_mem, &ark_mem->tempv3);
  arkFreeVec(ark_mem, &ark_mem->tempv4);
  arkFreeVec(ark_mem, &ark_mem->yn);
  arkFreeVec(ark_mem, &ark_mem->Vabstol);
}


/*---------------------------------------------------------------
  arkInitialSetup

  This routine performs all necessary items to prepare ARKode for
  the first internal step, including:
  - checks for valid initial step input or estimates first step
  - input consistency checks
  - checks the linear solver module (if applicable)
  - initializes linear solver (if applicable)
  ---------------------------------------------------------------*/
int arkInitialSetup(ARKodeMem ark_mem, realtype tout)
{
  int retval, hflag, istate, ier;
  realtype tout_hin, rh;

  /* Temporarily set ark_h */
  ark_mem->h = SUNRabs(tout - ark_mem->tcur);
  if (ark_mem->h == ZERO)  ark_mem->h = ONE;

  /* Set up the time stepper module */
  if (ark_mem->step_init == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInitialSetup", "Time stepper module is missing");
    return(ARK_ILL_INPUT);
  }
  retval = ark_mem->step_init(ark_mem, 0);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode", "arkInitialSetup",
                    "Error in initialization of time stepper module");
    return(retval);
  }

  /* Check that user has supplied an initial step size if fixedstep mode is on */
  if ( (ark_mem->fixedstep) && (ark_mem->hin == ZERO) ) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInitialSetup",
                    "Fixed step mode enabled, but no step size set");
    return(ARK_ILL_INPUT);
  }

  /* Set data for efun (if left unspecified) */
  if (ark_mem->user_efun)
    ark_mem->e_data = ark_mem->user_data;
  else
    ark_mem->e_data = ark_mem;

  /* Load initial error weights */
  ier = ark_mem->efun(ark_mem->yn,
                      ark_mem->ewt,
                      ark_mem->e_data);
  if (ier != 0) {
    if (ark_mem->itol == ARK_WF)
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                      "arkInitialSetup", MSG_ARK_EWT_FAIL);
    else
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                      "arkInitialSetup", MSG_ARK_BAD_EWT);
    return(ARK_ILL_INPUT);
  }

  /* Set data for rfun (if left unspecified) */
  if (ark_mem->user_rfun)
    ark_mem->r_data = ark_mem->user_data;
  else
    ark_mem->r_data = ark_mem;

  /* Load initial residual weights */
  if (ark_mem->rwt_is_ewt) {      /* update pointer to ewt */
    ark_mem->rwt = ark_mem->ewt;
  } else {
    ier = ark_mem->rfun(ark_mem->yn,
                        ark_mem->rwt,
                        ark_mem->r_data);
    if (ier != 0) {
      if (ark_mem->itol == ARK_WF)
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                        "arkInitialSetup", MSG_ARK_RWT_FAIL);
      else
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                        "arkInitialSetup", MSG_ARK_BAD_RWT);
      return(ARK_ILL_INPUT);
    }
  }

  /* Allocate interpolation memory (if unallocated, and if needed) */
  if (ark_mem->interp == NULL) {
    ark_mem->interp = arkInterpCreate(ark_mem);
    if (ark_mem->interp == NULL)
      return(ARK_MEM_FAIL);
  }

  /* Fill initial interpolation data (if needed) */
  if (ark_mem->interp != NULL) {
    ier = arkInterpInit(ark_mem, ark_mem->interp, ark_mem->tcur);
    if (ier != 0)  return(ier);
  }

  /* Test input tstop for legality. */
  if ( ark_mem->tstopset ) {
    if ( (ark_mem->tstop - ark_mem->tcur)*(tout - ark_mem->tcur) <= ZERO ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInitialSetup",
                      MSG_ARK_BAD_TSTOP, ark_mem->tstop, ark_mem->tcur);
      return(ARK_ILL_INPUT);
    }
  }

  /* Check input h for validity */
  ark_mem->h = ark_mem->hin;
  if ( (ark_mem->h != ZERO) &&
       ((tout-ark_mem->tcur)*ark_mem->h < ZERO) ) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInitialSetup",
                    MSG_ARK_BAD_H0);
    return(ARK_ILL_INPUT);
  }
  if ((ark_mem->hin == ZERO) && (ark_mem->fixedstep)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInitialSetup",
                    "nonzero step size must be supplied when using fixed-step mode");
    return(ARK_ILL_INPUT);
  }

  /* Estimate initial h if not set */
  if (ark_mem->h == ZERO) {
    /* Again, temporarily set ark_h for estimating an optimal value */
    ark_mem->h = SUNRabs(tout - ark_mem->tcur);
    if (ark_mem->h == ZERO)  ark_mem->h = ONE;
    /* Estimate the first step size */
    tout_hin = tout;
    if ( ark_mem->tstopset &&
         (tout-ark_mem->tcur)*(tout-ark_mem->tstop) > ZERO )
      tout_hin = ark_mem->tstop;
    hflag = arkHin(ark_mem, tout_hin);
    if (hflag != ARK_SUCCESS) {
      istate = arkHandleFailure(ark_mem, hflag);
      return(istate);
    }
  }

  /* Enforce step size bounds */
  rh = SUNRabs(ark_mem->h)*ark_mem->hmax_inv;
  if (rh > ONE) ark_mem->h /= rh;
  if (SUNRabs(ark_mem->h) < ark_mem->hmin)
    ark_mem->h *= ark_mem->hmin/SUNRabs(ark_mem->h);
  /* Check for approach to tstop */
  if (ark_mem->tstopset) {
    if ( (ark_mem->tcur + ark_mem->h - ark_mem->tstop)*ark_mem->h > ZERO ) {
      ark_mem->h = (ark_mem->tstop - ark_mem->tcur)*(ONE-FOUR*ark_mem->uround);
    }
  }

  /* Set initial time step factors */
  ark_mem->h0u    = ark_mem->h;
  ark_mem->hprime = ark_mem->h;

  /* Check for zeros of root function g at and near t0. */
  if (ark_mem->root_mem != NULL)
    if (ark_mem->root_mem->nrtfn > 0) {
      retval = arkRootCheck1((void*) ark_mem);

      if (retval == ARK_RTFUNC_FAIL) {
        arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKode", "arkRootCheck1",
                        MSG_ARK_RTFUNC_FAILED, ark_mem->tcur);
        return(ARK_RTFUNC_FAIL);
      }
    }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkPostResizeSetup

  This routine performs all necessary items to prepare ARKode for
  the first internal step after a resize() call, including:
  - re-initialize the linear solver
  - re-initialize the interpolation structure
  - check for approach to tstop
  - check for root near t0
  ---------------------------------------------------------------*/
int arkPostResizeSetup(ARKodeMem ark_mem)
{
  int retval, ier;

  /* Load updated error weights */
  ier = ark_mem->efun(ark_mem->yn,
                      ark_mem->ewt,
                      ark_mem->e_data);
  if (ier != 0) {
    if (ark_mem->itol == ARK_WF)
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                      "arkPostResizeSetup", MSG_ARK_EWT_FAIL);
    else
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                      "arkPostResizeSetup", MSG_ARK_BAD_EWT);
    return(ARK_ILL_INPUT);
  }

  /* Load updated residual weights */
  if (!ark_mem->rwt_is_ewt) {
    ier = ark_mem->rfun(ark_mem->yn,
                        ark_mem->rwt,
                        ark_mem->r_data);
    if (ier != 0) {
      if (ark_mem->itol == ARK_WF)
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                        "arkPostResizeSetup", MSG_ARK_RWT_FAIL);
      else
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                        "arkPostResizeSetup", MSG_ARK_BAD_RWT);
      return(ARK_ILL_INPUT);
    }
  }

  /* Fill initial interpolation data (if needed) */
  if (ark_mem->interp != NULL) {
    ier = arkInterpInit(ark_mem, ark_mem->interp, ark_mem->tcur);
    if (ier != 0)  return(ier);
  }

  /* Check for legal tstop (correct direction of integration) */
  if (ark_mem->tstopset) {
    if ( (ark_mem->tstop - ark_mem->tcur)*ark_mem->h < ZERO ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkPostResizeSetup",
                      MSG_ARK_BAD_TSTOP, ark_mem->tstop, ark_mem->tcur);
      return(ARK_ILL_INPUT);
    }
  }

  /* re-initialize the time stepper module */
  if (ark_mem->step_init == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkPostResizeSetup", "Time stepper module is missing");
    return(ARK_ILL_INPUT);
  }
  retval = ark_mem->step_init(ark_mem, 1);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode", "arkPostResizeSetup",
                    "Error in re-initialization of time stepper module");
    return(retval);
  }

  /* Check for zeros of root function g at and near t0. */
  if (ark_mem->root_mem != NULL)
    if (ark_mem->root_mem->nrtfn > 0) {
      retval = arkRootCheck1((void*) ark_mem);

      if (retval == ARK_RTFUNC_FAIL) {
        arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKode", "arkRootCheck1",
                        MSG_ARK_RTFUNC_FAILED, ark_mem->tcur);
        return(ARK_RTFUNC_FAIL);
      }
    }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStopTests

  This routine performs relevant stopping tests:
  - check for root in last step
  - check if we passed tstop
  - check if we passed tout (NORMAL mode)
  - check if current tn was returned (ONE_STEP mode)
  - check if we are close to tstop
  (adjust step size if needed)
  ---------------------------------------------------------------*/
int arkStopTests(ARKodeMem ark_mem, realtype tout, N_Vector yout,
                 realtype *tret, int itask, int *ier)
{
  int irfndp, retval;
  realtype troundoff;

  /* Estimate an infinitesimal time interval to be used as
     a roundoff for time quantities (based on current time
     and step size) */
  troundoff = FUZZ_FACTOR*ark_mem->uround *
    (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->h));

  /* First, check for a root in the last step taken, other than the
     last root found, if any.  If itask = ARK_ONE_STEP and y(tn) was not
     returned because of an intervening root, return y(tn) now.     */
  if (ark_mem->root_mem != NULL)
    if (ark_mem->root_mem->nrtfn > 0) {

      irfndp = ark_mem->root_mem->irfnd;

      retval = arkRootCheck2((void*) ark_mem);

      if (retval == CLOSERT) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkStopTests",
                        MSG_ARK_CLOSE_ROOTS, ark_mem->root_mem->tlo);
        *ier = ARK_ILL_INPUT;
        return(1);
      } else if (retval == ARK_RTFUNC_FAIL) {
        arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKode", "arkStopTests",
                        MSG_ARK_RTFUNC_FAILED, ark_mem->root_mem->tlo);
        *ier = ARK_RTFUNC_FAIL;
        return(1);
      } else if (retval == RTFOUND) {
        ark_mem->tretlast = *tret = ark_mem->root_mem->tlo;
        *ier = ARK_ROOT_RETURN;
        return(1);
      }

      /* If tn is distinct from tretlast (within roundoff),
         check remaining interval for roots */
      if ( SUNRabs(ark_mem->tcur - ark_mem->tretlast) > troundoff ) {

        retval = arkRootCheck3((void*) ark_mem);

        if (retval == ARK_SUCCESS) {     /* no root found */
          ark_mem->root_mem->irfnd = 0;
          if ((irfndp == 1) && (itask == ARK_ONE_STEP)) {
            ark_mem->tretlast = *tret = ark_mem->tcur;
            N_VScale(ONE, ark_mem->yn, yout);
            *ier = ARK_SUCCESS;
            return(1);
          }
        } else if (retval == RTFOUND) {  /* a new root was found */
          ark_mem->root_mem->irfnd = 1;
          ark_mem->tretlast = *tret = ark_mem->root_mem->tlo;
          *ier = ARK_ROOT_RETURN;
          return(1);
        } else if (retval == ARK_RTFUNC_FAIL) {  /* g failed */
          arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKode", "arkStopTests",
                          MSG_ARK_RTFUNC_FAILED, ark_mem->root_mem->tlo);
          *ier = ARK_RTFUNC_FAIL;
          return(1);
        }
      }

    } /* end of root stop check */

  /* In ARK_NORMAL mode, test if tout was reached */
  if ( (itask == ARK_NORMAL) &&
       ((ark_mem->tcur-tout)*ark_mem->h >= ZERO) ) {
    ark_mem->tretlast = *tret = tout;
    *ier = arkGetDky(ark_mem, tout, 0, yout);
    if (*ier != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                      "arkStopTests", MSG_ARK_BAD_TOUT, tout);
      *ier = ARK_ILL_INPUT;
      return(1);
    }
    *ier = ARK_SUCCESS;
    return(1);
  }

  /* In ARK_ONE_STEP mode, test if tn was returned */
  if ( itask == ARK_ONE_STEP &&
       SUNRabs(ark_mem->tcur - ark_mem->tretlast) > troundoff ) {
    ark_mem->tretlast = *tret = ark_mem->tcur;
    N_VScale(ONE, ark_mem->yn, yout);
    *ier = ARK_SUCCESS;
    return(1);
  }

  /* Test for tn at tstop or near tstop */
  if ( ark_mem->tstopset ) {

    if ( SUNRabs(ark_mem->tcur - ark_mem->tstop) <= troundoff) {
      *ier = arkGetDky(ark_mem, ark_mem->tstop, 0, yout);
      if (*ier != ARK_SUCCESS) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkStopTests",
                        MSG_ARK_BAD_TSTOP, ark_mem->tstop, ark_mem->tcur);
        *ier = ARK_ILL_INPUT;
        return(1);
      }
      ark_mem->tretlast = *tret = ark_mem->tstop;
      ark_mem->tstopset = SUNFALSE;
      *ier = ARK_TSTOP_RETURN;
      return(1);
    }

    /* If next step would overtake tstop, adjust stepsize */
    if ( (ark_mem->tcur + ark_mem->hprime - ark_mem->tstop)*ark_mem->h > ZERO ) {
      ark_mem->hprime = (ark_mem->tstop - ark_mem->tcur)*(ONE-FOUR*ark_mem->uround);
      ark_mem->eta = ark_mem->hprime/ark_mem->h;
    }
  }

  return(0);
}


/*---------------------------------------------------------------
  arkHin

  This routine computes a tentative initial step size h0.
  If tout is too close to tn (= t0), then arkHin returns
  ARK_TOO_CLOSE and h remains uninitialized. Note that here tout
  is either the value passed to arkEvolve at the first call or the
  value of tstop (if tstop is enabled and it is closer to t0=tn
  than tout). If the RHS function fails unrecoverably, arkHin
  returns ARK_RHSFUNC_FAIL. If the RHS function fails recoverably
  too many times and recovery is not possible, arkHin returns
  ARK_REPTD_RHSFUNC_ERR. Otherwise, arkHin sets h to the chosen
  value h0 and returns ARK_SUCCESS.

  The algorithm used seeks to find h0 as a solution of
  (WRMS norm of (h0^2 ydd / 2)) = 1,
  where ydd = estimated second derivative of y.  Although this
  choice is based on an error expansion of the Backward Euler
  method, and hence results in an overly-conservative time step
  for our higher-order ARK methods, it does find an order-of-
  magnitude estimate of the initial time scale of the solution.
  Since this method is only used on the first time step, the
  additional caution will not overly hinder solver efficiency.

  We start with an initial estimate equal to the geometric mean
  of the lower and upper bounds on the step size.

  Loop up to H0_ITERS times to find h0.
  Stop if new and previous values differ by a factor < 2.
  Stop if hnew/hg > 2 after one iteration, as this probably
  means that the ydd value is bad because of cancellation error.

  For each new proposed hg, we allow H0_ITERS attempts to
  resolve a possible recoverable failure from f() by reducing
  the proposed stepsize by a factor of 0.2. If a legal stepsize
  still cannot be found, fall back on a previous value if
  possible, or else return ARK_REPTD_RHSFUNC_ERR.

  Finally, we apply a bias (0.5) and verify that h0 is within
  bounds.
  ---------------------------------------------------------------*/
int arkHin(ARKodeMem ark_mem, realtype tout)
{
  int retval, sign, count1, count2;
  realtype tdiff, tdist, tround, hlb, hub;
  realtype hg, hgs, hs, hnew, hrat, h0, yddnrm;
  booleantype hgOK;

  /* If tout is too close to tn, give up */
  if ((tdiff = tout-ark_mem->tcur) == ZERO) return(ARK_TOO_CLOSE);

  sign = (tdiff > ZERO) ? 1 : -1;
  tdist = SUNRabs(tdiff);
  tround = ark_mem->uround * SUNMAX(SUNRabs(ark_mem->tcur), SUNRabs(tout));

  if (tdist < TWO*tround) return(ARK_TOO_CLOSE);

  /* Set lower and upper bounds on h0, and take geometric mean
     as first trial value.
     Exit with this value if the bounds cross each other. */
  hlb = H0_LBFACTOR * tround;
  hub = arkUpperBoundH0(ark_mem, tdist);

  hg  = SUNRsqrt(hlb*hub);

  if (hub < hlb) {
    if (sign == -1) ark_mem->h = -hg;
    else            ark_mem->h =  hg;
    return(ARK_SUCCESS);
  }

  /* Outer loop */
  hs = hg;     /* safeguard against 'uninitialized variable' warning */
  for(count1 = 1; count1 <= H0_ITERS; count1++) {

    /* Attempts to estimate ydd */
    hgOK = SUNFALSE;

    for (count2 = 1; count2 <= H0_ITERS; count2++) {
      hgs = hg*sign;
      retval = arkYddNorm(ark_mem, hgs, &yddnrm);
      /* If f() failed unrecoverably, give up */
      if (retval < 0) return(ARK_RHSFUNC_FAIL);
      /* If successful, we can use ydd */
      if (retval == ARK_SUCCESS) {hgOK = SUNTRUE; break;}
      /* f() failed recoverably; cut step size and test it again */
      hg *= POINT2;
    }

    /* If f() failed recoverably H0_ITERS times */
    if (!hgOK) {
      /* Exit if this is the first or second pass. No recovery possible */
      if (count1 <= 2) return(ARK_REPTD_RHSFUNC_ERR);
      /* We have a fall-back option. The value hs is a previous hnew which
         passed through f(). Use it and break */
      hnew = hs;
      break;
    }

    /* The proposed step size is feasible. Save it. */
    hs = hg;

    /* Propose new step size */
    hnew = (yddnrm*hub*hub > TWO) ? SUNRsqrt(TWO/yddnrm) : SUNRsqrt(hg*hub);
    
    /* If last pass, stop now with hnew */
    if (count1 == H0_ITERS) break;
    
    hrat = hnew/hg;

    /* Accept hnew if it does not differ from hg by more than a factor of 2 */
    if ((hrat > HALF) && (hrat < TWO)) break;

    /* After one pass, if ydd seems to be bad, use fall-back value. */
    if ((count1 > 1) && (hrat > TWO)) {
      hnew = hg;
      break;
    }

    /* Send this value back through f() */
    hg = hnew;
  }

  /* Apply bounds, bias factor, and attach sign */
  h0 = H0_BIAS*hnew;
  if (h0 < hlb) h0 = hlb;
  if (h0 > hub) h0 = hub;
  if (sign == -1) h0 = -h0;
  ark_mem->h = h0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkUpperBoundH0

  This routine sets an upper bound on abs(h0) based on
  tdist = tn - t0 and the values of y[i]/y'[i].
  ---------------------------------------------------------------*/
realtype arkUpperBoundH0(ARKodeMem ark_mem, realtype tdist)
{
  realtype hub_inv, hub;
  N_Vector temp1, temp2;

  /* Bound based on |y0|/|y0'| -- allow at most an increase of
   * H0_UBFACTOR in y0 (based on a forward Euler step). The weight
   * factor is used as a safeguard against zero components in y0. */
  temp1 = ark_mem->tempv1;
  temp2 = ark_mem->tempv2;

  N_VAbs(ark_mem->yn, temp2);
  ark_mem->efun(ark_mem->yn, temp1, ark_mem->e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(H0_UBFACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(ark_mem->interp->fnew, temp2);

  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);

  /* bound based on tdist -- allow at most a step of magnitude
   * H0_UBFACTOR * tdist */
  hub = H0_UBFACTOR*tdist;

  /* Use the smaller of the two */
  if (hub*hub_inv > ONE) hub = ONE/hub_inv;

  return(hub);
}


/*---------------------------------------------------------------
  arkYddNorm

  This routine computes an estimate of the second derivative of y
  using a difference quotient, and returns its WRMS norm.
  ---------------------------------------------------------------*/
int arkYddNorm(ARKodeMem ark_mem, realtype hg, realtype *yddnrm)
{
  int retval;

  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode", "arkYddNorm",
                    "Missing interpolation structure");
    return(ARK_MEM_NULL);
  }

  /* increment y with a multiple of f */
  N_VLinearSum(hg, ark_mem->interp->fnew, ONE,
               ark_mem->yn, ark_mem->ycur);

  /* compute y', via the ODE RHS routine */
  retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tcur+hg,
                                 ark_mem->ycur,
                                 ark_mem->tempv1, 2);
  if (retval != 0) return(ARK_RHSFUNC_FAIL);

  /* difference new f and original f to estimate y'' */
  N_VLinearSum(ONE/hg, ark_mem->tempv1, -ONE/hg,
	       ark_mem->interp->fnew, ark_mem->tempv1);

  /* reset ycur to equal yn (unnecessary?) */
  N_VScale(ONE, ark_mem->yn, ark_mem->ycur);
  
  /* compute norm of y'' */
  *yddnrm = N_VWrmsNorm(ark_mem->tempv1, ark_mem->ewt);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkCompleteStep

  This routine performs various update operations when the step
  solution is complete.  It is assumed that the timestepper
  module has stored the time-evolved solution in ark_mem->ycur,
  and the step that gave rise to this solution in ark_mem->h.
  We update the current time (tn), the current solution (yn),
  increment the overall step counter nst, record the values hold
  and tnew, reset the resized flag, allow for user-provided
  postprocessing, and update the interpolation structure.
  ---------------------------------------------------------------*/
int arkCompleteStep(ARKodeMem ark_mem)
{
  int retval;

  /* Set current time to the end of the step (in case the last
     stage time does not coincide with the step solution time).
     If tstop is enabled, it is possible for tn + h to be past
     tstop by roundoff, and in that case, we reset tn (after
     incrementing by h) to tstop. */
  ark_mem->tcur = ark_mem->tn + ark_mem->h;
  if (ark_mem->tstopset) {
    if ((ark_mem->tcur - ark_mem->tstop)*ark_mem->h > ZERO)
      ark_mem->tcur = ark_mem->tstop;
  }

  /* apply user-supplied step postprocessing function (if supplied) */
  if (ark_mem->ProcessStep != NULL) {
    retval = ark_mem->ProcessStep(ark_mem->tcur,
                                  ark_mem->ycur,
                                  ark_mem->user_data);
    if (retval != 0) return(ARK_POSTPROCESS_FAIL);
  }

  /* update interpolation structure */
  if (ark_mem->interp != NULL) {
    retval = arkInterpUpdate(ark_mem, ark_mem->interp,
                             ark_mem->tcur,
                             (ark_mem->ProcessStep != NULL));
    if (retval != ARK_SUCCESS)  return(retval);
  }

  /* update yn to current solution */
  N_VScale(ONE, ark_mem->ycur, ark_mem->yn);

  /* update scalar quantities */
  ark_mem->nst++;
  ark_mem->hold = ark_mem->h;
  ark_mem->tn   = ark_mem->tcur;

  /* turn off flag regarding resized problem */
  ark_mem->resized = SUNFALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkHandleFailure

  This routine prints error messages for all cases of failure by
  arkHin and ark_step. It returns to ARKode the value that ARKode
  is to return to the user.
  ---------------------------------------------------------------*/
int arkHandleFailure(ARKodeMem ark_mem, int flag)
{

  /* Depending on flag, print error message and return error flag */
  switch (flag) {
  case ARK_ERR_FAILURE:
    arkProcessError(ark_mem, ARK_ERR_FAILURE, "ARKode", "ARKode",
                    MSG_ARK_ERR_FAILS, ark_mem->tcur, ark_mem->h);
    break;
  case ARK_CONV_FAILURE:
    arkProcessError(ark_mem, ARK_CONV_FAILURE, "ARKode", "ARKode",
                    MSG_ARK_CONV_FAILS, ark_mem->tcur, ark_mem->h);
    break;
  case ARK_LSETUP_FAIL:
    arkProcessError(ark_mem, ARK_LSETUP_FAIL, "ARKode", "ARKode",
                    MSG_ARK_SETUP_FAILED, ark_mem->tcur);
    break;
  case ARK_LSOLVE_FAIL:
    arkProcessError(ark_mem, ARK_LSOLVE_FAIL, "ARKode", "ARKode",
                    MSG_ARK_SOLVE_FAILED, ark_mem->tcur);
    break;
  case ARK_RHSFUNC_FAIL:
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode", "ARKode",
                    MSG_ARK_RHSFUNC_FAILED, ark_mem->tcur);
    break;
  case ARK_UNREC_RHSFUNC_ERR:
    arkProcessError(ark_mem, ARK_UNREC_RHSFUNC_ERR, "ARKode", "ARKode",
                    MSG_ARK_RHSFUNC_UNREC, ark_mem->tcur);
    break;
  case ARK_REPTD_RHSFUNC_ERR:
    arkProcessError(ark_mem, ARK_REPTD_RHSFUNC_ERR, "ARKode", "ARKode",
                    MSG_ARK_RHSFUNC_REPTD, ark_mem->tcur);
    break;
  case ARK_RTFUNC_FAIL:
    arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKode", "ARKode",
                    MSG_ARK_RTFUNC_FAILED, ark_mem->tcur);
    break;
  case ARK_TOO_CLOSE:
    arkProcessError(ark_mem, ARK_TOO_CLOSE, "ARKode", "ARKode",
                    MSG_ARK_TOO_CLOSE);
    break;
  case ARK_MASSSOLVE_FAIL:
    arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKode", "ARKode",
                    MSG_ARK_MASSSOLVE_FAIL);
    break;
  case ARK_NLS_SETUP_FAIL:
    arkProcessError(ark_mem, ARK_NLS_SETUP_FAIL, "ARKode", "ARKode",
                    "At t = %Lg the nonlinear solver setup failed unrecoverably",
                    (long double) ark_mem->tcur);
    break;
  case ARK_VECTOROP_ERR:
    arkProcessError(ark_mem, ARK_VECTOROP_ERR, "ARKode", "ARKode",
                    MSG_ARK_VECTOROP_ERR, ark_mem->tcur);
    break;
  case ARK_INNERSTEP_FAIL:
    arkProcessError(ark_mem, ARK_INNERSTEP_FAIL, "ARKode", "ARKode",
                    MSG_ARK_INNERSTEP_FAILED, ark_mem->tcur);
    break;
  default:
    /* This return should never happen */
    arkProcessError(ark_mem, ARK_UNRECOGNIZED_ERROR, "ARKode", "ARKode",
                    "ARKode encountered an unrecognized error. Please report this to the Sundials developers at sundials-users@llnl.gov");
    return(ARK_UNRECOGNIZED_ERROR);
  }

  return(flag);
}


/*---------------------------------------------------------------
  arkEwtSetSS

  This routine sets ewt as decribed above in the case tol_type = ARK_SS.
  It tests for non-positive components before inverting. arkEwtSetSS
  returns 0 if ewt is successfully set to a positive vector
  and -1 otherwise. In the latter case, ewt is considered undefined.
  ---------------------------------------------------------------*/
int arkEwtSetSS(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, ark_mem->tempv1);
  N_VScale(ark_mem->reltol, ark_mem->tempv1, ark_mem->tempv1);
  N_VAddConst(ark_mem->tempv1, ark_mem->Sabstol, ark_mem->tempv1);
  if (N_VMin(ark_mem->tempv1) <= ZERO) return(-1);
  N_VInv(ark_mem->tempv1, weight);
  return(0);
}


/*---------------------------------------------------------------
  arkEwtSetSV

  This routine sets ewt as decribed above in the case tol_type = ARK_SV.
  It tests for non-positive components before inverting. arkEwtSetSV
  returns 0 if ewt is successfully set to a positive vector
  and -1 otherwise. In the latter case, ewt is considered undefined.
  ---------------------------------------------------------------*/
int arkEwtSetSV(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, ark_mem->tempv1);
  N_VLinearSum(ark_mem->reltol, ark_mem->tempv1, ONE,
               ark_mem->Vabstol, ark_mem->tempv1);
  if (N_VMin(ark_mem->tempv1) <= ZERO) return(-1);
  N_VInv(ark_mem->tempv1, weight);
  return(0);
}


/*---------------------------------------------------------------
  arkRwtSetSS

  This routine sets rwt as decribed above in the case tol_type = ARK_SS.
  It tests for non-positive components before inverting. arkRwtSetSS
  returns 0 if rwt is successfully set to a positive vector
  and -1 otherwise. In the latter case, rwt is considered undefined.
  ---------------------------------------------------------------*/
int arkRwtSetSS(ARKodeMem ark_mem, N_Vector My, N_Vector weight)
{
  N_VAbs(My, ark_mem->tempv1);
  N_VScale(ark_mem->reltol, ark_mem->tempv1, ark_mem->tempv1);
  N_VAddConst(ark_mem->tempv1, ark_mem->SRabstol, ark_mem->tempv1);
  if (N_VMin(ark_mem->tempv1) <= ZERO) return(-1);
  N_VInv(ark_mem->tempv1, weight);
  return(0);
}


/*---------------------------------------------------------------
  arkRwtSetSV

  This routine sets rwt as decribed above in the case tol_type = ARK_SV.
  It tests for non-positive components before inverting. arkRwtSetSV
  returns 0 if rwt is successfully set to a positive vector
  and -1 otherwise. In the latter case, rwt is considered undefined.
  ---------------------------------------------------------------*/
int arkRwtSetSV(ARKodeMem ark_mem, N_Vector My, N_Vector weight)
{
  N_VAbs(My, ark_mem->tempv1);
  N_VLinearSum(ark_mem->reltol, ark_mem->tempv1, ONE,
               ark_mem->VRabstol, ark_mem->tempv1);
  if (N_VMin(ark_mem->tempv1) <= ZERO) return(-1);
  N_VInv(ark_mem->tempv1, weight);
  return(0);
}


/*---------------------------------------------------------------
  arkExpStab is the default explicit stability estimation function
  ---------------------------------------------------------------*/
int arkExpStab(N_Vector y, realtype t, realtype *hstab, void *data)
{
  /* explicit stability not used by default,
     set to zero to disable */
  *hstab = RCONST(0.0);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkPredict_MaximumOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKode interpolation module.  This uses the
  highest-degree interpolant supported by the module (stored
  as dense_q in the ark_mem structure).
  ---------------------------------------------------------------*/
int arkPredict_MaximumOrder(ARKodeMem ark_mem, realtype tau, N_Vector yguess)
{

  /* verify that ark_mem and interpolation structure are provided */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkPredict_MaximumOrder",
                    "ARKodeMem structure is NULL");
    return(ARK_MEM_NULL);
  }
  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode",
                    "arkPredict_MaximumOrder",
                    "ARKodeInterpMem structure is NULL");
    return(ARK_MEM_NULL);
  }

  /* call the interpolation module to do the work */
  return(arkInterpEvaluate(ark_mem, ark_mem->interp, tau,
                           0, ark_mem->dense_q, yguess));
}


/*---------------------------------------------------------------
  arkPredict_VariableOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKode interpolation module.  The degree of the
  interpolant is based on the level of extrapolation outside the
  preceding time step.
  ---------------------------------------------------------------*/
int arkPredict_VariableOrder(ARKodeMem ark_mem, realtype tau, N_Vector yguess)
{
  int ord;
  realtype tau_tol  = 0.5;
  realtype tau_tol2 = 0.75;

  /* verify that ark_mem and interpolation structure are provided */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkPredict_VariableOrder",
                    "ARKodeMem structure is NULL");
    return(ARK_MEM_NULL);
  }
  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode",
                    "arkPredict_VariableOrder",
                    "ARKodeInterpMem structure is NULL");
    return(ARK_MEM_NULL);
  }

  /* set the polynomial order based on tau input */
  if (tau <= tau_tol) {
    ord = 3;
  } else if (tau <= tau_tol2) {
    ord = 2;
  } else {
    ord = 1;
  }

  /* call the interpolation module to do the work */
  return(arkInterpEvaluate(ark_mem, ark_mem->interp, tau,
                           0, ord, yguess));
}


/*---------------------------------------------------------------
  arkPredict_CutoffOrder

  This routine predicts the nonlinear implicit stage solution
  using the ARKode interpolation module.  If the level of
  extrapolation is small enough, it uses the maximum degree
  polynomial available (stored in ark_mem->dense_q); otherwise
  it uses a linear polynomial.
  ---------------------------------------------------------------*/
int arkPredict_CutoffOrder(ARKodeMem ark_mem, realtype tau, N_Vector yguess)
{
  int ord;
  realtype tau_tol = 0.5;

  /* verify that ark_mem and interpolation structure are provided */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkPredict_CutoffOrder",
                    "ARKodeMem structure is NULL");
    return(ARK_MEM_NULL);
  }
  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode",
                    "arkPredict_CutoffOrder",
                    "ARKodeInterpMem structure is NULL");
    return(ARK_MEM_NULL);
  }

  /* set the polynomial order based on tau input */
  if (tau <= tau_tol) {
    ord = ark_mem->dense_q;
  } else {
    ord = 1;
  }

  /* call the interpolation module to do the work */
  return(arkInterpEvaluate(ark_mem, ark_mem->interp, tau,
                           0, ord, yguess));
}


/*---------------------------------------------------------------
  arkPredict_Bootstrap

  This routine predicts the nonlinear implicit stage solution
  using a quadratic Hermite interpolating polynomial, based on
  the data {y_n, f(t_n,y_n), f(t_n+hj,z_j)}.

  Note: we assume that ftemp = f(t_n+hj,z_j) can be computed via 
     N_VLinearCombination(nvec, cvals, Xvecs, ftemp),
  i.e. the inputs cvals[0:nvec-1] and Xvecs[0:nvec-1] may be 
  combined to form f(t_n+hj,z_j).
  ---------------------------------------------------------------*/
int arkPredict_Bootstrap(ARKodeMem ark_mem, realtype hj,
                         realtype tau, int nvec, realtype *cvals,
                         N_Vector *Xvecs, N_Vector yguess)
{
  realtype a0, a1, a2;
  int i;

  /* verify that ark_mem and interpolation structure are provided */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkPredict_Bootstrap",
                    "ARKodeMem structure is NULL");
    return(ARK_MEM_NULL);
  }
  if (ark_mem->interp == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode",
                    "arkPredict_Bootstrap",
                    "ARKodeInterpMem structure is NULL");
    return(ARK_MEM_NULL);
  }

  /* set coefficients for Hermite interpolant */
  a0 = ONE;
  a2 = tau*tau/TWO/hj;
  a1 = tau - a2;

  /* set arrays for fused vector operation; shift inputs for 
     f(t_n+hj,z_j) to end of queue */
  for (i=0; i<nvec; i++) {
    cvals[2+i] = a2*cvals[i];
    Xvecs[2+i] = Xvecs[i];
  }
  cvals[0] = a0;
  Xvecs[0] = ark_mem->yn;
  cvals[1] = a1;
  Xvecs[1] = ark_mem->interp->fnew;

  /* call fused vector operation to compute prediction */
  return(N_VLinearCombination(nvec+2, cvals, Xvecs, yguess));
}


/*---------------------------------------------------------------
  arkProcessError is a high level error handling function
  - if ark_mem==NULL it prints the error message to stderr
  - otherwise, it sets-up and calls the error handling function
    pointed to by ark_ehfun
  ---------------------------------------------------------------*/
void arkProcessError(ARKodeMem ark_mem, int error_code,
                     const char *module, const char *fname,
                     const char *msgfmt, ...)
{
  va_list ap;
  char msg[256];

  /* Initialize the argument pointer variable
     (msgfmt is the last required argument to arkProcessError) */
  va_start(ap, msgfmt);

  /* Compose the message */
  vsprintf(msg, msgfmt, ap);

  if (ark_mem == NULL) {    /* We write to stderr */

#ifndef NO_FPRINTF_OUTPUT
    fprintf(stderr, "\n[%s ERROR]  %s\n  ", module, fname);
    fprintf(stderr, "%s\n\n", msg);
#endif

  } else {                 /* We can call ehfun */
    ark_mem->ehfun(error_code, module, fname, msg,
                   ark_mem->eh_data);
  }

  /* Finalize argument processing */
  va_end(ap);

  return;
}


/*===============================================================
  EOF
  ===============================================================*/
