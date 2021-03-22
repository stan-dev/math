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
#include "arkode_interp_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

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

  /* Initialize time step module to NULL */
  ark_mem->step_attachlinsol   = NULL;
  ark_mem->step_attachmasssol  = NULL;
  ark_mem->step_disablelsetup  = NULL;
  ark_mem->step_disablemsetup  = NULL;
  ark_mem->step_getlinmem      = NULL;
  ark_mem->step_getmassmem     = NULL;
  ark_mem->step_getimplicitrhs = NULL;
  ark_mem->step_mmult          = NULL;
  ark_mem->step_getgammas      = NULL;
  ark_mem->step_init           = NULL;
  ark_mem->step_fullrhs        = NULL;
  ark_mem->step                = NULL;
  ark_mem->step_mem            = NULL;

  /* Initialize root finding variables */
  ark_mem->root_mem = NULL;

  /* Initialize inequality constraints variables */
  ark_mem->constraintsSet = SUNFALSE;
  ark_mem->constraints    = NULL;

  /* Initialize diagnostics reporting variables */
  ark_mem->report  = SUNFALSE;
  ark_mem->diagfp  = NULL;

  /* Initialize lrw and liw */
  ark_mem->lrw = 18;
  ark_mem->liw = 39;  /* fcn/data ptr, int, long int, sunindextype, booleantype */

  /* No mallocs have been done yet */
  ark_mem->VabstolMallocDone     = SUNFALSE;
  ark_mem->VRabstolMallocDone    = SUNFALSE;
  ark_mem->MallocDone            = SUNFALSE;

  /* No user-supplied step postprocessing function yet */
  ark_mem->ProcessStep  = NULL;
  ark_mem->ps_data      = NULL;

  /* No user-supplied stage postprocessing function yet */
  ark_mem->ProcessStage = NULL;

  /* No user_data pointer yet */
  ark_mem->user_data = NULL;

  /* Allocate step adaptivity structure and note storage */
  ark_mem->hadapt_mem = arkAdaptInit();
  if (ark_mem->hadapt_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_FAIL, "ARKode", "arkCreate",
                    "Allocation of step adaptivity structure failed");
    return(NULL);
  }
  ark_mem->lrw += ARK_ADAPT_LRW;
  ark_mem->liw += ARK_ADAPT_LIW;

  /* Initialize the interpolation structure to NULL */
  ark_mem->interp = NULL;

  /* Initially, rwt should point to ewt */
  ark_mem->rwt_is_ewt = SUNTRUE;

  /* Indicate that evaluation of the full RHS is not required after each step,
     this flag is updated to SUNTRUE by the interpolation module initialization
     function and/or the stepper initialization function in arkInitialSetup */
  ark_mem->call_fullrhs = SUNFALSE;

  /* Indicate that the problem needs to be initialized */
  ark_mem->initsetup   = SUNTRUE;
  ark_mem->init_type   = FIRST_INIT;
  ark_mem->firststage  = SUNTRUE;
  ark_mem->initialized = SUNFALSE;

  /* Initial step size has not been determined yet */
  ark_mem->h   = ZERO;
  ark_mem->h0u = ZERO;

  /* Set default values for integrator optional inputs */
  iret = arkSetDefaults(ark_mem);
  if (iret != ARK_SUCCESS) {
    arkProcessError(NULL, 0, "ARKode", "arkCreate",
                    "Error setting default solver options");
    return(NULL);
  }

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
  booleantype resizeOK;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int retval;

  /* Check ark_mem */
  if (ark_mem == NULL) {
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
  if (hscale < ZERO)  hscale = ONE;
  if (hscale != ONE) {

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

  /* Resize the solver vectors (using y0 as a template) */
  resizeOK = arkResizeVectors(ark_mem, resize, resize_data,
                              lrw_diff, liw_diff, y0);
  if (!resizeOK) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode",
                    "arkResize", "Unable to resize vector");
    return(ARK_MEM_FAIL);
  }

  /* Resize the interpolation structure memory */
  if (ark_mem->interp != NULL) {
    retval = arkInterpResize(ark_mem, ark_mem->interp, resize,
                             resize_data, lrw_diff, liw_diff, y0);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, retval, "ARKode", "arkResize",
                      "Interpolation module resize failure");
      return(retval);
    }
  }

  /* Copy y0 into ark_yn to set the current solution */
  N_VScale(ONE, y0, ark_mem->yn);

  /* Disable constraints */
  ark_mem->constraintsSet = SUNFALSE;

  /* Indicate that problem needs to be initialized */
  ark_mem->initsetup  = SUNTRUE;
  ark_mem->init_type  = RESIZE_INIT;
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

  /* Set flag indicating whether abstol == 0 */
  ark_mem->atolmin0 = (abstol == ZERO);

  /* Copy tolerances into memory */
  ark_mem->reltol  = reltol;
  ark_mem->Sabstol = abstol;
  ark_mem->itol    = ARK_SS;

  /* enforce use of arkEwtSetSS */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->efun      = arkEwtSetSS;
  ark_mem->e_data    = ark_mem;

  return(ARK_SUCCESS);
}


int arkSVtolerances(ARKodeMem ark_mem, realtype reltol, N_Vector abstol)
{
  /* local variables */
  realtype abstolmin;

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
  if (abstol == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSVtolerances", MSG_ARK_NULL_ABSTOL);
    return(ARK_ILL_INPUT);
  }
  if (abstol->ops->nvmin == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkSVtolerances",
                    "Missing N_VMin routine from N_Vector");
    return(ARK_ILL_INPUT);
  }
  abstolmin = N_VMin(abstol);
  if (abstolmin < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkSVtolerances", MSG_ARK_BAD_ABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Set flag indicating whether min(abstol) == 0 */
  ark_mem->atolmin0 = (abstolmin == ZERO);

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

  /* enforce use of arkEwtSetSV */
  ark_mem->user_efun = SUNFALSE;
  ark_mem->efun      = arkEwtSetSV;
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
  ark_mem->e_data    = ark_mem->user_data;

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

  /* Set flag indicating whether rabstol == 0 */
  ark_mem->Ratolmin0 = (rabstol == ZERO);

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
  /* local variables */
  realtype rabstolmin;

  /* Check inputs */
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkResVtolerance", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkResVtolerance", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }
  if (rabstol == NULL) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode",
                    "arkResVtolerance", MSG_ARK_NULL_RABSTOL);
    return(ARK_NO_MALLOC);
  }
  if (rabstol->ops->nvmin == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkResVtolerance",
                   "Missing N_VMin routine from N_Vector");
    return(ARK_ILL_INPUT);
  }
  rabstolmin = N_VMin(rabstol);
  if (rabstolmin < ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkResVtolerance", MSG_ARK_BAD_RABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Set flag indicating whether min(abstol) == 0 */
  ark_mem->Ratolmin0 = (rabstolmin == ZERO);

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
  ark_mem->r_data    = ark_mem->user_data;

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
  int retval, kflag, istate, ir;
  int ewtsetOK;
  realtype troundoff, nrm;
  booleantype inactive_roots;
  realtype dsm;
  int nflag, ncf, nef, constrfails;

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
  if (ark_mem->initsetup) {
    ark_mem->tretlast = *tret = ark_mem->tcur;
    retval = arkInitialSetup(ark_mem, tout);
    if (retval!= ARK_SUCCESS) return(retval);
  }

  /* perform stopping tests */
  if (!ark_mem->initsetup)
    if (arkStopTests(ark_mem, tout, yout, tret, itask, &retval))
      return(retval);


  /*--------------------------------------------------
    Looping point for successful internal steps

    - update the ewt/rwt vectors for upcoming step
    - check for errors (too many steps, too much
      accuracy requested, step size too small)
    - loop over attempts at a new step:
      * try to take step (via time stepper module),
        handle solver convergence or other failures
      * perform constraint-handling (if selected)
      * check temporal error
      * if all of the above pass, complete step by
        updating current time, solution, error &
        stepsize history arrays.
    - perform stop tests:
      * check for root in last step taken
      * check if tout was passed
      * check if close to tstop
      * check if in ONE_STEP mode (must return)
    --------------------------------------------------*/
  nstloc = 0;
  for(;;) {

    ark_mem->next_h = ark_mem->h;

    /* Reset and check ewt and rwt */
    if (!ark_mem->initsetup) {

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

      if (!ark_mem->rwt_is_ewt) {
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
    if ( (ark_mem->mxstep > 0) && (nstloc >= ark_mem->mxstep) ) {
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
    if (ark_mem->hprime != ark_mem->h) {
      ark_mem->h = ark_mem->h * ark_mem->eta;
      ark_mem->next_h = ark_mem->h;
    }
    if (ark_mem->fixedstep) {
      ark_mem->h = ark_mem->hin;
      ark_mem->next_h = ark_mem->h;

      /* patch for 'fixedstep' + 'tstop' use case:
         limit fixed step size if step would overtake tstop */
      if ( ark_mem->tstopset ) {
        if ( (ark_mem->tcur + ark_mem->h - ark_mem->tstop)*ark_mem->h > ZERO ) {
          ark_mem->h = (ark_mem->tstop - ark_mem->tcur) *
            (ONE-FOUR*ark_mem->uround);
        }
      }
    }

    /* Looping point for step attempts */
    dsm = ZERO;
    ncf = nef = constrfails = ark_mem->last_kflag = 0;
    nflag = FIRST_CALL;
    for(;;) {

      /* increment attempt counter */
      ark_mem->nst_attempts++;

      /* Call time stepper module to attempt a step:
            0 => step completed successfully
           >0 => step encountered recoverable failure; reduce step if possible
           <0 => step encountered unrecoverable failure */
      kflag = ark_mem->step((void*) ark_mem, &dsm, &nflag);
      if (kflag < 0)  break;

      /* handle solver convergence failures */
      kflag = arkCheckConvergence(ark_mem, &nflag, &ncf);
      if (kflag < 0)  break;

      /* perform constraint-handling (if selected, and if solver check passed) */
      if (ark_mem->constraintsSet && (kflag == ARK_SUCCESS)) {
        kflag = arkCheckConstraints(ark_mem, &constrfails, &nflag);
        if (kflag < 0)  break;
      }

      /* when fixed time-stepping is enabled, 'success' == successful stage solves
         (checked in previous block), so just enforce no step size change */
      if (ark_mem->fixedstep) {
        ark_mem->eta = ONE;
        break;
      }

      /* check temporal error (if checks above passed) */
      if (kflag == ARK_SUCCESS) {
        kflag = arkCheckTemporalError(ark_mem, &nflag, &nef, dsm);
        if (kflag < 0)  break;
      }

      /* if we've made it here then no nonrecoverable failures occurred; someone above
         has recommended an 'eta' value for the next step -- enforce bounds on that value
         and set upcoming step size */
      ark_mem->eta = SUNMIN(ark_mem->eta, ark_mem->hadapt_mem->etamax);
      ark_mem->eta = SUNMAX(ark_mem->eta, ark_mem->hmin / SUNRabs(ark_mem->h));
      ark_mem->eta /= SUNMAX(ONE, SUNRabs(ark_mem->h) * ark_mem->hmax_inv*ark_mem->eta);

      /* if ignoring temporal error test result (XBraid) force step to pass */
      if (ark_mem->force_pass) {
        ark_mem->last_kflag = kflag;
        kflag = ARK_SUCCESS;
        break;
      }

      /* break attempt loop on successful step */
      if (kflag == ARK_SUCCESS)  break;

      /* unsuccessful step, if |h| = hmin, return ARK_ERR_FAILURE */
      if (SUNRabs(ark_mem->h) <= ark_mem->hmin*ONEPSM)  return(ARK_ERR_FAILURE);

      /* update h, hprime and next_h for next iteration */
      ark_mem->h *= ark_mem->eta;
      ark_mem->next_h = ark_mem->hprime = ark_mem->h;

    }

    /* If step attempt loop succeeded, complete step (update current time, solution,
       error stepsize history arrays; call user-supplied step postprocessing function)
       (added stuff from arkStep_PrepareNextStep -- revisit) */
    if (kflag == ARK_SUCCESS)  kflag = arkCompleteStep(ark_mem, dsm);

    /* If step attempt loop failed, process flag and return to user */
    if (kflag != ARK_SUCCESS) {
      istate = arkHandleFailure(ark_mem, kflag);
      ark_mem->tretlast = *tret = ark_mem->tcur;
      N_VScale(ONE, ark_mem->yn, yout);
      break;
    }

    nstloc++;

    /* Check for root in last step taken. */
    if (ark_mem->root_mem != NULL) {
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
      /* limit upcoming step if it will overcome tstop */
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
  dky. This routine internally calls arkInterpEvaluate to perform
  the interpolation.  We have the restriction that 0 <= k <= 3.
  This routine uses an interpolating polynomial of degree
  max(deg, k), i.e. it will form a polynomial of the degree
  available by the interpolation module and/or requested by
  the user through deg, unless higher-order derivatives are
  requested.

  This function is called by arkEvolve with k=0 and t=tout to
  perform interpolation of outputs, but may also be called
  indirectly by the user via time step module *StepGetDky calls.
  Note: in all cases it will be called after ark_tcur has been
  updated to correspond with the end time of the last successful
  step.
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
                             k, ARK_INTERP_MAX_DEGREE, dky);
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

  /* free vector storage */
  arkFreeVectors(ark_mem);

  /* free the time step adaptivity module */
  if (ark_mem->hadapt_mem != NULL) {
    free(ark_mem->hadapt_mem);
    ark_mem->hadapt_mem = NULL;
  }

  /* free the interpolation module */
  if (ark_mem->interp != NULL) {
    arkInterpFree(ark_mem, ark_mem->interp);
    ark_mem->interp = NULL;
  }

  /* free the root-finding module */
  if (ark_mem->root_mem != NULL) {
    (void) arkRootFree(*arkode_mem);
    ark_mem->root_mem = NULL;
  }

  free(*arkode_mem);
  *arkode_mem = NULL;
}



/*===============================================================
  Internal functions that may be replaced by the user
  ===============================================================*/

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
  called prior to calling arkEvolve to evolve the problem. The
  initialization type indicates if the values of internal counters
  should be reinitialized (FIRST_INIT) or retained (RESET_INIT).
  ---------------------------------------------------------------*/
int arkInit(ARKodeMem ark_mem, realtype t0, N_Vector y0,
            int init_type)
{
  booleantype stepperOK, nvectorOK, allocOK;
  sunindextype lrw1, liw1;

  /* Check ark_mem */
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Check if reset was called before the first Evolve call */
  if (init_type == RESET_INIT && !(ark_mem->initialized))
    init_type = FIRST_INIT;

  /* Check if allocations have been done i.e., is this first init call */
  if (ark_mem->MallocDone == SUNFALSE) {

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

    /* Create default Hermite interpolation module */
    ark_mem->interp = arkInterpCreate_Hermite(ark_mem, ARK_INTERP_MAX_DEGREE);
    if (ark_mem->interp == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode", "arkInit",
                      "Unable to allocate interpolation module");
      return(ARK_MEM_FAIL);
    }

    /* All allocations are complete */
    ark_mem->MallocDone = SUNTRUE;

  }

  /* All allocation and error checking is complete at this point */

  /* Copy the input parameters into ARKode state */
  ark_mem->tcur = t0;
  ark_mem->tn   = t0;

  /* Initialize yn */
  N_VScale(ONE, y0, ark_mem->yn);

  /* Initializations on (re-)initialization call, skip on reset */
  if (init_type == FIRST_INIT) {

    /* Counters */
    ark_mem->nst_attempts = 0;
    ark_mem->nst          = 0;
    ark_mem->nhnil        = 0;
    ark_mem->ncfn         = 0;
    ark_mem->netf         = 0;
    ark_mem->nconstrfails = 0;

    /* Initial, old, and next step sizes */
    ark_mem->h0u    = ZERO;
    ark_mem->hold   = ZERO;
    ark_mem->next_h = ZERO;

    /* Tolerance scale factor */
    ark_mem->tolsf = ONE;

    /* Adaptivity counters */
    ark_mem->hadapt_mem->nst_acc = 0;
    ark_mem->hadapt_mem->nst_exp = 0;

    /* Error and step size history */
    ark_mem->hadapt_mem->ehist[0] = ONE;
    ark_mem->hadapt_mem->ehist[1] = ONE;
    ark_mem->hadapt_mem->hhist[0] = ZERO;
    ark_mem->hadapt_mem->hhist[1] = ZERO;

    /* Indicate that evaluation of the full RHS is not required after each step,
       this flag is updated to SUNTRUE by the interpolation module initialization
       function and/or the stepper initialization function in arkInitialSetup */
    ark_mem->call_fullrhs = SUNFALSE;

    /* Indicate that initialization has not been done before */
    ark_mem->initialized = SUNFALSE;
  }

  /* Indicate initialization is needed */
  ark_mem->initsetup  = SUNTRUE;
  ark_mem->init_type  = init_type;
  ark_mem->firststage = SUNTRUE;

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
  fprintf(outfile, "itol = %i\n", ark_mem->itol);
  fprintf(outfile, "ritol = %i\n", ark_mem->ritol);
  fprintf(outfile, "mxhnil = %i\n", ark_mem->mxhnil);
  fprintf(outfile, "mxstep = %li\n", ark_mem->mxstep);
  fprintf(outfile, "lrw1 = %li\n", (long int) ark_mem->lrw1);
  fprintf(outfile, "liw1 = %li\n", (long int) ark_mem->liw1);
  fprintf(outfile, "lrw = %li\n", (long int) ark_mem->lrw);
  fprintf(outfile, "liw = %li\n", (long int) ark_mem->liw);
  fprintf(outfile, "user_efun = %i\n", ark_mem->user_efun);
  fprintf(outfile, "tstopset = %i\n", ark_mem->tstopset);
  fprintf(outfile, "tstop = %" RSYM"\n", ark_mem->tstop);
  fprintf(outfile, "report = %i\n", ark_mem->report);
  fprintf(outfile, "VabstolMallocDone = %i\n", ark_mem->VabstolMallocDone);
  fprintf(outfile, "MallocDone = %i\n", ark_mem->MallocDone);
  fprintf(outfile, "initsetup = %i\n", ark_mem->initsetup);
  fprintf(outfile, "init_type = %i\n", ark_mem->init_type);
  fprintf(outfile, "firststage = %i\n", ark_mem->firststage);
  fprintf(outfile, "uround = %" RSYM"\n", ark_mem->uround);
  fprintf(outfile, "reltol = %" RSYM"\n", ark_mem->reltol);
  fprintf(outfile, "Sabstol = %" RSYM"\n", ark_mem->Sabstol);
  fprintf(outfile, "fixedstep = %i\n", ark_mem->fixedstep);
  fprintf(outfile, "tolsf = %" RSYM"\n", ark_mem->tolsf);
  fprintf(outfile, "call_fullrhs = %i\n", ark_mem->call_fullrhs);

  /* output counters */
  fprintf(outfile, "nhnil = %i\n", ark_mem->nhnil);
  fprintf(outfile, "nst_attempts = %li\n", ark_mem->nst_attempts);
  fprintf(outfile, "nst = %li\n", ark_mem->nst);
  fprintf(outfile, "ncfn = %li\n", ark_mem->ncfn);
  fprintf(outfile, "netf = %li\n", ark_mem->netf);

  /* output time-stepping values */
  fprintf(outfile, "hin = %" RSYM"\n", ark_mem->hin);
  fprintf(outfile, "h = %" RSYM"\n", ark_mem->h);
  fprintf(outfile, "hprime = %" RSYM"\n", ark_mem->hprime);
  fprintf(outfile, "next_h = %" RSYM"\n", ark_mem->next_h);
  fprintf(outfile, "eta = %" RSYM"\n", ark_mem->eta);
  fprintf(outfile, "tcur = %" RSYM"\n", ark_mem->tcur);
  fprintf(outfile, "tretlast = %" RSYM"\n", ark_mem->tretlast);
  fprintf(outfile, "hmin = %" RSYM"\n", ark_mem->hmin);
  fprintf(outfile, "hmax_inv = %" RSYM"\n", ark_mem->hmax_inv);
  fprintf(outfile, "h0u = %" RSYM"\n", ark_mem->h0u);
  fprintf(outfile, "tn = %" RSYM"\n", ark_mem->tn);
  fprintf(outfile, "hold = %" RSYM"\n", ark_mem->hold);
  fprintf(outfile, "maxnef = %i\n", ark_mem->maxnef);
  fprintf(outfile, "maxncf = %i\n", ark_mem->maxncf);

  /* output time-stepping adaptivity structure */
  fprintf(outfile, "timestep adaptivity structure:\n");
  arkPrintAdaptMem(ark_mem->hadapt_mem, outfile);

  /* output inequality constraints quantities */
  fprintf(outfile, "constraintsSet = %i\n", ark_mem->constraintsSet);
  fprintf(outfile, "maxconstrfails = %i\n", ark_mem->maxconstrfails);

  /* output root-finding quantities */
  if (ark_mem->root_mem != NULL)
    (void) arkPrintRootMem((void*) ark_mem, outfile);

  /* output interpolation quantities */
  arkInterpPrintMem(ark_mem->interp, outfile);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */
  fprintf(outfile, "Vapbsol:\n");
  N_VPrintFile(ark_mem->Vabstol, outfile);
  fprintf(outfile, "ewt:\n");
  N_VPrintFile(ark_mem->ewt, outfile);
  if (!ark_mem->rwt_is_ewt) {
    fprintf(outfile, "rwt:\n");
    N_VPrintFile(ark_mem->rwt, outfile);
  }
  fprintf(outfile, "ycur:\n");
  N_VPrintFile(ark_mem->ycur, outfile);
  fprintf(outfile, "yn:\n");
  N_VPrintFile(ark_mem->yn, outfile);
  fprintf(outfile, "fn:\n");
  N_VPrintFile(ark_mem->fn, outfile);
  fprintf(outfile, "tempv1:\n");
  N_VPrintFile(ark_mem->tempv1, outfile);
  fprintf(outfile, "tempv2:\n");
  N_VPrintFile(ark_mem->tempv2, outfile);
  fprintf(outfile, "tempv3:\n");
  N_VPrintFile(ark_mem->tempv3, outfile);
  fprintf(outfile, "tempv4:\n");
  N_VPrintFile(ark_mem->tempv4, outfile);
  fprintf(outfile, "constraints:\n");
  N_VPrintFile(ark_mem->constraints, outfile);
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
  if ( (ark_mem->step_init    == NULL) ||
       (ark_mem->step         == NULL) ||
       (ark_mem->step_mem     == NULL) ||
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
      (tmpl->ops->nvwrmsnorm  == NULL))
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
  vector. If the ARKVecResizeFn function is non-NULL, then it
  calls that routine to perform the single-vector resize;
  otherwise it deallocates and reallocates the target vector based
  on the template vector. This routine also updates the optional
  outputs lrw and liw, which are (respectively) the lengths of the
  overall ARKode real and integer work spaces.

  If the resize is successful then this returns SUNTRUE,
  otherwise it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype arkResizeVec(ARKodeMem ark_mem, ARKVecResizeFn resize,
                         void *resize_data, sunindextype lrw_diff,
                         sunindextype liw_diff, N_Vector tmpl, N_Vector *v)
{
  if (*v != NULL) {
    if (resize == NULL) {
      N_VDestroy(*v);
      *v = NULL;
      *v = N_VClone(tmpl);
      if (*v == NULL) {
        arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode",
                        "arkResizeVec", "Unable to clone vector");
        return(SUNFALSE);
      }
    } else {
      if (resize(*v, tmpl, resize_data)) {
        arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode",
                        "arkResizeVec", MSG_ARK_RESIZE_FAIL);
        return(SUNFALSE);
      }
    }
    ark_mem->lrw += lrw_diff;
    ark_mem->liw += liw_diff;
  }
  return(SUNTRUE);
}


/*---------------------------------------------------------------
  arkAllocVectors:

  This routine allocates the ARKode vectors ewt, yn, tempv* and
  ftemp. If any of these vectors already exist, they are left
  alone. Otherwise, it will allocate each vector by cloning the
  input vector. This routine also updates the optional outputs
  lrw and liw, which are (respectively) the lengths of the real
  and integer work spaces.

  If all memory allocations are successful, arkAllocVectors
  returns SUNTRUE, otherwise it returns SUNFALSE.
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

  /* Allocate fn if needed */
  if (!arkAllocVec(ark_mem, tmpl, &ark_mem->fn))
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
  arkResizeVectors:

  This routine resizes all ARKode vectors if they exist,
  otherwise they are left alone. If a resize function is provided
  it is called to resize the vectors otherwise the vector is
  freed and a new vector is created by cloning in input vector.
  This routine also updates the optional outputs lrw and liw,
  which are (respectively) the lengths of the real and integer
  work spaces.

  If all memory allocations are successful, arkResizeVectors
  returns SUNTRUE, otherwise it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype arkResizeVectors(ARKodeMem ark_mem, ARKVecResizeFn resize,
                             void *resize_data, sunindextype lrw_diff,
                             sunindextype liw_diff, N_Vector tmpl)
{
  /* Vabstol */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->Vabstol))
    return(SUNFALSE);

  /* VRabstol */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->VRabstol))
    return(SUNFALSE);

  /* ewt */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->ewt))
    return(SUNFALSE);

  /* rwt  */
  if (ark_mem->rwt_is_ewt) {      /* update pointer to ewt */
    ark_mem->rwt = ark_mem->ewt;
  } else {                        /* resize if distinct from ewt */
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                      liw_diff, tmpl, &ark_mem->rwt))
      return(SUNFALSE);
  }

  /* yn */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->yn))
    return(SUNFALSE);

  /* fn */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->fn))
    return(SUNFALSE);

  /* tempv* */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->tempv1))
    return(SUNFALSE);

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->tempv2))
    return(SUNFALSE);

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->tempv3))
    return(SUNFALSE);

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->tempv4))
    return(SUNFALSE);

  /* constraints */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, tmpl, &ark_mem->constraints))
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
  arkFreeVec(ark_mem, &ark_mem->fn);
  arkFreeVec(ark_mem, &ark_mem->Vabstol);
  arkFreeVec(ark_mem, &ark_mem->constraints);
}


/*---------------------------------------------------------------
  arkInitialSetup

  This routine performs all necessary items to prepare ARKode for
  the first internal step after initialization, reinitialization,
  a reset() call, or a resize() call, including:
  - input consistency checks
  - (re)initializes the stepper
  - computes error and residual weights
  - (re)initialize the interpolation structure
  - checks for valid initial step input or estimates first step
  - checks for approach to tstop
  - checks for root near t0
  ---------------------------------------------------------------*/
int arkInitialSetup(ARKodeMem ark_mem, realtype tout)
{
  int retval, hflag, istate;
  realtype tout_hin, rh, htmp;
  booleantype conOK;

  /* Set up the time stepper module */
  if (ark_mem->step_init == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInitialSetup", "Time stepper module is missing");
    return(ARK_ILL_INPUT);
  }
  retval = ark_mem->step_init(ark_mem, ark_mem->init_type);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode", "arkInitialSetup",
                    "Error in initialization of time stepper module");
    return(retval);
  }

  /* Check that user has supplied an initial step size if fixedstep mode is on */
  if ( (ark_mem->fixedstep) && (ark_mem->hin == ZERO) ) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInitialSetup",
                    "Fixed step mode enabled, but no step size set");
    return(ARK_ILL_INPUT);
  }

  /* If using a built-in routine for error/residual weights with abstol==0,
     ensure that N_VMin is available */
  if ((!ark_mem->user_efun) && (ark_mem->atolmin0) && (!ark_mem->yn->ops->nvmin)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInitialSetup",
                    "N_VMin unimplemented (required by error-weight function)");
    return(ARK_ILL_INPUT);
  }
  if ( (!ark_mem->user_rfun) && (!ark_mem->rwt_is_ewt) &&
       (ark_mem->Ratolmin0) && (!ark_mem->yn->ops->nvmin) ) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInitialSetup",
                    "N_VMin unimplemented (required by residual-weight function)");
    return(ARK_ILL_INPUT);
  }

  /* Test input tstop for legality (correct direction of integration) */
  if ( ark_mem->tstopset ) {
    htmp = (ark_mem->h == ZERO) ? tout - ark_mem->tcur : ark_mem->h;
    if ( (ark_mem->tstop - ark_mem->tcur) * htmp <= ZERO ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInitialSetup",
                      MSG_ARK_BAD_TSTOP, ark_mem->tstop, ark_mem->tcur);
      return(ARK_ILL_INPUT);
    }
  }

  /* Check to see if y0 satisfies constraints */
  if (ark_mem->constraintsSet) {
    conOK = N_VConstrMask(ark_mem->constraints, ark_mem->yn, ark_mem->tempv1);
    if (!conOK) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInitialSetup", MSG_ARK_Y0_FAIL_CONSTR);
      return(ARK_ILL_INPUT);
    }
  }

  /* Load initial error weights */
  retval = ark_mem->efun(ark_mem->yn, ark_mem->ewt, ark_mem->e_data);
  if (retval != 0) {
    if (ark_mem->itol == ARK_WF)
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                      "arkInitialSetup", MSG_ARK_EWT_FAIL);
    else
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                      "arkInitialSetup", MSG_ARK_BAD_EWT);
    return(ARK_ILL_INPUT);
  }

  /* Load initial residual weights */
  if (ark_mem->rwt_is_ewt) {      /* update pointer to ewt */
    ark_mem->rwt = ark_mem->ewt;
  } else {
    retval = ark_mem->rfun(ark_mem->yn, ark_mem->rwt, ark_mem->r_data);
    if (retval != 0) {
      if (ark_mem->itol == ARK_WF)
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                        "arkInitialSetup", MSG_ARK_RWT_FAIL);
      else
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                        "arkInitialSetup", MSG_ARK_BAD_RWT);
      return(ARK_ILL_INPUT);
    }
  }

  /* If necessary, temporarily set h as it is used to compute the tolerance in a
     potential mass matrix solve when computing the full rhs */
  if (ark_mem->h == ZERO) ark_mem->h = ONE;

  /* Call fullrhs (used in estimating initial step, explicit steppers, Hermite
     interpolation module, and possibly (but not always) arkRootCheck1) */
  retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tcur,
                                 ark_mem->yn, ark_mem->fn, 0);
  if (retval != 0) return(ARK_RHSFUNC_FAIL);

  /* Fill initial interpolation data (if needed) */
  if (ark_mem->interp != NULL) {
    retval = arkInterpInit(ark_mem, ark_mem->interp, ark_mem->tcur);
    if (retval != 0)  return(retval);
  }

  /* initialization complete */
  ark_mem->initialized = SUNTRUE;

  /* Set initial step size */
  if (ark_mem->h0u == ZERO) {

    /* Check input h for validity */
    ark_mem->h = ark_mem->hin;
    if ( (ark_mem->h != ZERO) &&
         ((tout-ark_mem->tcur)*ark_mem->h < ZERO) ) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkInitialSetup",
                      MSG_ARK_BAD_H0);
      return(ARK_ILL_INPUT);
    }

    /* Estimate initial h if not set */
    if (ark_mem->h == ZERO) {
      /* Again, temporarily set h for estimating an optimal value */
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
      /* Use first step growth factor for estimated h */
      ark_mem->hadapt_mem->etamax = ark_mem->hadapt_mem->etamx1;
    } else if (ark_mem->nst == 0) {
      /* Use first step growth factor for user defined h */
      ark_mem->hadapt_mem->etamax = ark_mem->hadapt_mem->etamx1;
    } else {
      /* Use standard growth factor (e.g., for reset) */
      ark_mem->hadapt_mem->etamax = ark_mem->hadapt_mem->growth;
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
    ark_mem->eta    = ONE;
    ark_mem->hprime = ark_mem->h;
  }

  /* Check for zeros of root function g at and near t0. */
  if (ark_mem->root_mem != NULL) {
    if (ark_mem->root_mem->nrtfn > 0) {
      retval = arkRootCheck1((void*) ark_mem);

      if (retval == ARK_RTFUNC_FAIL) {
        arkProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKode", "arkRootCheck1",
                        MSG_ARK_RTFUNC_FAILED, ark_mem->tcur);
        return(ARK_RTFUNC_FAIL);
      }
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

      /* Shortcut to roots found in previous step */
      irfndp = ark_mem->root_mem->irfnd;

      /* If the full rhs was not computed in the last call to arkCompleteStep
         and roots were found in the previous step, then compute the full rhs
         for possible use in arkRootCheck2 (not always necessary) */
      if (!(ark_mem->call_fullrhs) && irfndp != 0) {
        retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tcur,
                                       ark_mem->yn, ark_mem->fn, 1);
        if (retval != 0) {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode", "arkStopTests",
                          MSG_ARK_RHSFUNC_FAILED);
          *ier = ARK_RHSFUNC_FAIL;
          return(1);
        }
      }

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
      hg *= RCONST(0.2);
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

  N_VAbs(ark_mem->fn, temp2);

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
  N_VLinearSum(hg, ark_mem->fn, ONE, ark_mem->yn, ark_mem->ycur);

  /* compute y', via the ODE RHS routine */
  retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tcur+hg,
                                 ark_mem->ycur,
                                 ark_mem->tempv1, 2);
  if (retval != 0) return(ARK_RHSFUNC_FAIL);

  /* difference new f and original f to estimate y'' */
  N_VLinearSum(ONE/hg, ark_mem->tempv1, -ONE/hg, ark_mem->fn, ark_mem->tempv1);

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
  and tnew, allow for user-provided postprocessing, and update
  the interpolation structure.
  ---------------------------------------------------------------*/
int arkCompleteStep(ARKodeMem ark_mem, realtype dsm)
{
  int retval, mode;

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
                                  ark_mem->ps_data);
    if (retval != 0) return(ARK_POSTPROCESS_STEP_FAIL);
  }

  /* update interpolation structure */
  if (ark_mem->interp != NULL) {
    retval = arkInterpUpdate(ark_mem, ark_mem->interp, ark_mem->tcur);
    if (retval != ARK_SUCCESS)  return(retval);
  }

  /* call fullrhs if needed */
  if (ark_mem->call_fullrhs) {
    mode = (ark_mem->ProcessStep != NULL) ? 0 : 1;
    retval = ark_mem->step_fullrhs(ark_mem, ark_mem->tcur,
                                   ark_mem->ycur,
                                   ark_mem->fn, mode);
    if (retval != 0) return(ARK_RHSFUNC_FAIL);
  }

  /* update yn to current solution */
  N_VScale(ONE, ark_mem->ycur, ark_mem->yn);

  /* Update step size and error history arrays */
  ark_mem->hadapt_mem->ehist[1] = ark_mem->hadapt_mem->ehist[0];
  ark_mem->hadapt_mem->ehist[0] = dsm*ark_mem->hadapt_mem->bias;
  ark_mem->hadapt_mem->hhist[1] = ark_mem->hadapt_mem->hhist[0];
  ark_mem->hadapt_mem->hhist[0] = ark_mem->h;

  /* update scalar quantities */
  ark_mem->nst++;
  ark_mem->hold   = ark_mem->h;
  ark_mem->tn     = ark_mem->tcur;
  ark_mem->hprime = ark_mem->h * ark_mem->eta;

  /* Reset growth factor for subsequent time step */
  ark_mem->hadapt_mem->etamax = ark_mem->hadapt_mem->growth;

  /* Turn off flag indicating initial step and first stage */
  ark_mem->initsetup  = SUNFALSE;
  ark_mem->firststage = SUNFALSE;

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
  case ARK_CONSTR_FAIL:
    arkProcessError(ark_mem, ARK_CONSTR_FAIL, "ARKode", "ARKode",
                    MSG_ARK_FAILED_CONSTR, ark_mem->tcur);
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
  case ARK_NLS_OP_ERR:
    arkProcessError(ark_mem, ARK_NLS_OP_ERR, "ARKode", "ARKode",
                    MSG_ARK_NLS_FAIL, ark_mem->tcur);
    break;
  case ARK_USER_PREDICT_FAIL:
    arkProcessError(ark_mem, ARK_USER_PREDICT_FAIL, "ARKode", "ARKode",
                    MSG_ARK_USER_PREDICT_FAIL, ark_mem->tcur);
    break;
  case ARK_POSTPROCESS_STEP_FAIL:
    arkProcessError(ark_mem, ARK_POSTPROCESS_STEP_FAIL, "ARKode", "ARKode",
                    MSG_ARK_POSTPROCESS_STEP_FAIL, ark_mem->tcur);
    break;
  case ARK_POSTPROCESS_STAGE_FAIL:
    arkProcessError(ark_mem, ARK_POSTPROCESS_STAGE_FAIL, "ARKode", "ARKode",
                    MSG_ARK_POSTPROCESS_STAGE_FAIL, ark_mem->tcur);
    break;
  case ARK_INTERP_FAIL:
    arkProcessError(ark_mem, ARK_INTERP_FAIL, "ARKode", "ARKode",
                    "At t = %Lg the interpolation module failed unrecoverably",
                    (long double) ark_mem->tcur);
    break;
  case ARK_INVALID_TABLE:
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode", "ARKode",
                    "ARKode was provided an invalid method table");
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

  This routine is responsible for setting the error weight vector
  ewt as follows:

  ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol), i=0,...,neq-1

  When the absolute tolerance is zero, it tests for non-positive
  components before inverting. arkEwtSetSS returns 0 if ewt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, ewt is considered undefined.
  ---------------------------------------------------------------*/
int arkEwtSetSS(N_Vector ycur, N_Vector weight, void* arkode_mem)
{
  ARKodeMem ark_mem = (ARKodeMem) arkode_mem;
  N_VAbs(ycur, ark_mem->tempv1);
  N_VScale(ark_mem->reltol, ark_mem->tempv1, ark_mem->tempv1);
  N_VAddConst(ark_mem->tempv1, ark_mem->Sabstol, ark_mem->tempv1);
  if (ark_mem->atolmin0) {
    if (N_VMin(ark_mem->tempv1) <= ZERO) return(-1);
  }
  N_VInv(ark_mem->tempv1, weight);
  return(0);
}


/*---------------------------------------------------------------
  arkEwtSetSV

  This routine is responsible for setting the error weight vector
  ewt as follows:

  ewt[i] = 1 / (reltol * SUNRabs(ycur[i]) + abstol[i]), i=0,...,neq-1

  When any absolute tolerance is zero, it tests for non-positive
  components before inverting. arkEwtSetSV returns 0 if ewt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, ewt is considered undefined.
  ---------------------------------------------------------------*/
int arkEwtSetSV(N_Vector ycur, N_Vector weight, void* arkode_mem)
{
  ARKodeMem ark_mem = (ARKodeMem) arkode_mem;
  N_VAbs(ycur, ark_mem->tempv1);
  N_VLinearSum(ark_mem->reltol, ark_mem->tempv1, ONE,
               ark_mem->Vabstol, ark_mem->tempv1);
  if (ark_mem->atolmin0) {
    if (N_VMin(ark_mem->tempv1) <= ZERO) return(-1);
  }
  N_VInv(ark_mem->tempv1, weight);
  return(0);
}


/*---------------------------------------------------------------
  arkEwtSetSmallReal

  This routine is responsible for setting the error weight vector
  ewt as follows:

  ewt[i] = SMALL_REAL

  This is routine is only used with explicit time stepping with
  a fixed step size to avoid a potential too much error return
  to the user.
  ---------------------------------------------------------------*/
int arkEwtSetSmallReal(N_Vector ycur, N_Vector weight, void* arkode_mem)
{
  N_VConst(SMALL_REAL, weight);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkRwtSetSS

  This routine sets rwt as decribed above in the case tol_type = ARK_SS.
  When the absolute tolerance is zero, it tests for non-positive
  components before inverting. arkRwtSetSS returns 0 if rwt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, rwt is considered undefined.
  ---------------------------------------------------------------*/
int arkRwtSetSS(ARKodeMem ark_mem, N_Vector My, N_Vector weight)
{
  N_VAbs(My, ark_mem->tempv1);
  N_VScale(ark_mem->reltol, ark_mem->tempv1, ark_mem->tempv1);
  N_VAddConst(ark_mem->tempv1, ark_mem->SRabstol, ark_mem->tempv1);
  if (ark_mem->Ratolmin0) {
    if (N_VMin(ark_mem->tempv1) <= ZERO) return(-1);
  }
  N_VInv(ark_mem->tempv1, weight);
  return(0);
}


/*---------------------------------------------------------------
  arkRwtSetSV

  This routine sets rwt as decribed above in the case tol_type = ARK_SV.
  When any absolute tolerance is zero, it tests for non-positive
  components before inverting. arkRwtSetSV returns 0 if rwt is
  successfully set to a positive vector and -1 otherwise. In the
  latter case, rwt is considered undefined.
  ---------------------------------------------------------------*/
int arkRwtSetSV(ARKodeMem ark_mem, N_Vector My, N_Vector weight)
{
  N_VAbs(My, ark_mem->tempv1);
  N_VLinearSum(ark_mem->reltol, ark_mem->tempv1, ONE,
               ark_mem->VRabstol, ark_mem->tempv1);
  if (ark_mem->Ratolmin0) {
    if (N_VMin(ark_mem->tempv1) <= ZERO) return(-1);
  }
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
  in the interpolation module).
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
                           0, ARK_INTERP_MAX_DEGREE, yguess));
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
  polynomial available (stored in the interpolation module
  structure); otherwise it uses a linear polynomial.
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
    ord = ARK_INTERP_MAX_DEGREE;
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
  int i, retval;

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
  Xvecs[1] = ark_mem->fn;

  /* call fused vector operation to compute prediction */
  retval = N_VLinearCombination(nvec+2, cvals, Xvecs, yguess);
  if (retval != 0)  return(ARK_VECTOROP_ERR);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkCheckConvergence

  This routine checks the return flag from the time-stepper's
  "step" routine for algebraic solver convergence issues.

  Returns ARK_SUCCESS (0) if successful, PREDICT_AGAIN (>0)
  on a recoverable convergence failure, or a relevant
  nonrecoverable failure flag (<0).
  --------------------------------------------------------------*/
int arkCheckConvergence(ARKodeMem ark_mem, int *nflagPtr, int *ncfPtr)
{
  ARKodeHAdaptMem hadapt_mem;

  if (*nflagPtr == ARK_SUCCESS) return(ARK_SUCCESS);

  /* The nonlinear soln. failed; increment ncfn */
  ark_mem->ncfn++;

  /* If fixed time stepping, then return with convergence failure */
  if (ark_mem->fixedstep) return(ARK_CONV_FAILURE);

  /* Otherwise, access adaptivity structure */
  if (ark_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode", "arkCheckConvergence",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = ark_mem->hadapt_mem;

  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (*nflagPtr < 0) {
    if (*nflagPtr == ARK_LSETUP_FAIL)       return(ARK_LSETUP_FAIL);
    else if (*nflagPtr == ARK_LSOLVE_FAIL)  return(ARK_LSOLVE_FAIL);
    else if (*nflagPtr == ARK_RHSFUNC_FAIL) return(ARK_RHSFUNC_FAIL);
    else                                    return(ARK_NLS_OP_ERR);
  }

  /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
  (*ncfPtr)++;
  hadapt_mem->etamax = ONE;

  /* If we had maxncf failures, or if |h| = hmin,
     return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */
  if ((*ncfPtr == ark_mem->maxncf) ||
      (SUNRabs(ark_mem->h) <= ark_mem->hmin*ONEPSM)) {
    if (*nflagPtr == CONV_FAIL)     return(ARK_CONV_FAILURE);
    if (*nflagPtr == RHSFUNC_RECVR) return(ARK_REPTD_RHSFUNC_ERR);
  }

  /* Reduce step size due to convergence failure */
  ark_mem->eta = hadapt_mem->etacf;

  /* Signal for Jacobian/preconditioner setup */
  *nflagPtr = PREV_CONV_FAIL;

  /* Return to reattempt the step */
  return(PREDICT_AGAIN);
}


/*---------------------------------------------------------------
  arkCheckConstraints

  This routine determines if the constraints of the problem
  are satisfied by the proposed step

  Returns ARK_SUCCESS if successful, otherwise CONSTR_RECVR
  --------------------------------------------------------------*/
int arkCheckConstraints(ARKodeMem ark_mem, int *constrfails, int *nflag)
{
  booleantype constraintsPassed;
  N_Vector mm  = ark_mem->tempv4;
  N_Vector tmp = ark_mem->tempv1;

  /* Check constraints and get mask vector mm for where constraints failed */
  constraintsPassed = N_VConstrMask(ark_mem->constraints, ark_mem->ycur, mm);
  if (constraintsPassed) return(ARK_SUCCESS);

  /* Constraints not met */

  /* Update total fails and fails in current step */
  ark_mem->nconstrfails++;
  (*constrfails)++;

  /* Return with error if reached max fails in a step */
  if (*constrfails == ark_mem->maxconstrfails) return(ARK_CONSTR_FAIL);

  /* Return with error if using fixed step sizes */
  if (ark_mem->fixedstep) return(ARK_CONSTR_FAIL);

  /* Return with error if |h| == hmin */
  if (SUNRabs(ark_mem->h) <= ark_mem->hmin*ONEPSM) return(ARK_CONSTR_FAIL);

  /* Reduce h by computing eta = h'/h */
  N_VLinearSum(ONE, ark_mem->yn, -ONE, ark_mem->ycur, tmp);
  N_VProd(mm, tmp, tmp);
  ark_mem->eta = RCONST(0.9)*N_VMinQuotient(ark_mem->yn, tmp);
  ark_mem->eta = SUNMAX(ark_mem->eta, TENTH);

  /* Signal for Jacobian/preconditioner setup */
  *nflag = PREV_CONV_FAIL;

  /* Return to reattempt the step */
  return(CONSTR_RECVR);
}


/*---------------------------------------------------------------
  arkCheckTemporalError

  This routine performs the local error test for the method.
  The weighted local error norm dsm is passed in.  This value is
  used to predict the next step to attempt based on dsm.
  The test dsm <= 1 is made, and if this fails then additional
  checks are performed based on the number of successive error
  test failures.

  Returns ARK_SUCCESS if the test passes.

  If the test fails:
    - if maxnef error test failures have occurred or if
      SUNRabs(h) = hmin, we return ARK_ERR_FAILURE.
    - otherwise: set *nflagPtr to PREV_ERR_FAIL, and
      return TRY_AGAIN.
  --------------------------------------------------------------*/
int arkCheckTemporalError(ARKodeMem ark_mem, int *nflagPtr, int *nefPtr, realtype dsm)
{
  int retval;
  realtype ttmp;
  long int nsttmp;
  ARKodeHAdaptMem hadapt_mem;

  /* Access hadapt_mem structure */
  if (ark_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode", "arkCheckTemporalError",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = ark_mem->hadapt_mem;

  /* consider change of step size for next step attempt (may be
     larger/smaller than current step, depending on dsm) */
  ttmp = (dsm <= ONE) ? ark_mem->tn + ark_mem->h : ark_mem->tn;
  nsttmp = (dsm <= ONE) ? ark_mem->nst+1 : ark_mem->nst;
  retval = arkAdapt((void*) ark_mem, hadapt_mem, ark_mem->ycur, ttmp,
                    ark_mem->h, dsm*ark_mem->hadapt_mem->bias, nsttmp);
  if (retval != ARK_SUCCESS)  return(ARK_ERR_FAILURE);

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */
  if (dsm <= ONE) return(ARK_SUCCESS);

  /* Test failed; increment counters, set nflag */
  (*nefPtr)++;
  ark_mem->netf++;
  *nflagPtr = PREV_ERR_FAIL;

  /* At maxnef failures, return ARK_ERR_FAILURE */
  if (*nefPtr == ark_mem->maxnef)  return(ARK_ERR_FAILURE);

  /* Set etamax=1 to prevent step size increase at end of this step */
  hadapt_mem->etamax = ONE;

  /* Enforce failure bounds on eta, update h, and return for retry of step */
  if (*nefPtr >= hadapt_mem->small_nef)
    ark_mem->eta = SUNMIN(ark_mem->eta, hadapt_mem->etamxf);
  return(TRY_AGAIN);
}


/*---------------------------------------------------------------
  arkAccessHAdaptMem:

  Shortcut routine to unpack ark_mem and hadapt_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int arkAccessHAdaptMem(void* arkode_mem, const char *fname,
                       ARKodeMem *ark_mem, ARKodeHAdaptMem *hadapt_mem)
{

  /* access ARKodeMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    fname, MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem) arkode_mem;
  if ((*ark_mem)->hadapt_mem==NULL) {
    arkProcessError(*ark_mem, ARK_MEM_NULL, "ARKode",
                    fname, MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *hadapt_mem = (ARKodeHAdaptMem) (*ark_mem)->hadapt_mem;
  return(ARK_SUCCESS);
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
