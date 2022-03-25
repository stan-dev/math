/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 *                Rujeko Chinomona @ SMU
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for ARKODE's MRI time stepper module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"
#include "arkode_interp_impl.h"
#include <sundials/sundials_math.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>


/*===============================================================
  SHORTCUTS
  ===============================================================*/

#define ARK_PROFILER ark_mem->sunctx->profiler

/*===============================================================
  MRIStep Exported functions -- Required
  ===============================================================*/

/*---------------------------------------------------------------
  Create MRIStep integrator memory struct
  ---------------------------------------------------------------*/
void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, realtype t0, N_Vector y0,
                    MRIStepInnerStepper stepper, SUNContext sunctx)
{
  ARKodeMem          ark_mem;         /* outer ARKode memory   */
  ARKodeMRIStepMem   step_mem;        /* outer stepper memory  */
  SUNNonlinearSolver NLS;             /* default nonlin solver */
  booleantype        nvectorOK;
  int                retval;

  /* Check that at least one of fse, fsi is supplied and is to be used*/
  if (fse == NULL && fsi == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_NULL_F);
    return(NULL);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_NULL_Y0);
    return(NULL);
  }

  /* Check that stepper is supplied */
  if (stepper == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepCreate",
                    "The inner stepper memory is NULL");
    return(NULL);
  }

  /* Check that context is supplied */
  if (!sunctx) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_NULL_SUNCTX);
    return(NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = mriStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_BAD_NVECTOR);
    return(NULL);
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate(sunctx);
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_NO_MEM);
    return(NULL);
  }

  /* Allocate ARKodeMRIStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeMRIStepMem) malloc(sizeof(struct ARKodeMRIStepMemRec));
  if (step_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                    "MRIStepCreate", MSG_ARK_ARKMEM_FAIL);
    MRIStepFree((void**) &ark_mem);  return(NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeMRIStepMemRec));

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_attachlinsol   = mriStep_AttachLinsol;
  ark_mem->step_disablelsetup  = mriStep_DisableLSetup;
  ark_mem->step_getlinmem      = mriStep_GetLmem;
  ark_mem->step_getimplicitrhs = mriStep_GetImplicitRHS;
  ark_mem->step_getgammas      = mriStep_GetGammas;
  ark_mem->step_init           = mriStep_Init;
  ark_mem->step_fullrhs        = mriStep_FullRHS;
  ark_mem->step                = mriStep_TakeStep;
  ark_mem->step_mem            = (void*) step_mem;

  /* Set default values for MRIStep optional inputs */
  retval = MRIStepSetDefaults((void *) ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::MRIStep",
                    "MRIStepCreate",
                    "Error setting default solver options");
    MRIStepFree((void**) &ark_mem);  return(NULL);
  }

  /* Allocate the general MRI stepper vectors using y0 as a template */
  /* NOTE: Fse, Fsi, inner_forcing, cvals, Xvecs, sdata, zpred and zcor will
     be allocated later on (based on the MRI method) */

  /* Copy the slow RHS functions into stepper memory */
  step_mem->fse = fse;
  step_mem->fsi = fsi;

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit_rhs = (fse == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit_rhs = (fsi == NULL) ? SUNFALSE : SUNTRUE;

  /* Update the ARKode workspace requirements */
  ark_mem->liw += 42;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->lrw += 10;

  /* Create a default Newton NLS object (just in case; will be deleted if
     the user attaches a nonlinear solver) */
  step_mem->NLS    = NULL;
  step_mem->ownNLS = SUNFALSE;

  if (step_mem->implicit_rhs) {
    NLS = SUNNonlinSol_Newton(y0, ark_mem->sunctx);
    if (!NLS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepCreate", "Error creating default Newton solver");
      MRIStepFree((void**) &ark_mem);  return(NULL);
    }
    retval = MRIStepSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepCreate", "Error attaching default Newton solver");
      MRIStepFree((void**) &ark_mem);  return(NULL);
    }
    step_mem->ownNLS = SUNTRUE;
  }

  /* Set the linear solver addresses to NULL (we check != NULL later) */
  step_mem->linit  = NULL;
  step_mem->lsetup = NULL;
  step_mem->lsolve = NULL;
  step_mem->lfree  = NULL;
  step_mem->lmem   = NULL;

  /* Initialize all the counters */
  step_mem->nfse      = 0;
  step_mem->nfsi      = 0;
  step_mem->nsetups   = 0;
  step_mem->nstlp     = 0;
  step_mem->nls_iters = 0;

  /* Initialize fused op work space */
  step_mem->cvals        = NULL;
  step_mem->Xvecs        = NULL;

  /* Initialize pre and post inner evolve functions */
  step_mem->pre_inner_evolve  = NULL;
  step_mem->post_inner_evolve = NULL;

  /* Initialize main ARKode infrastructure (allocates vectors) */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::MRIStep", "MRIStepCreate",
                    "Unable to initialize main ARKode infrastructure");
    MRIStepFree((void**) &ark_mem);  return(NULL);
  }

  /* Attach the inner stepper memory */
  step_mem->stepper = stepper;

  /* Check for required stepper functions */
  retval = mriStepInnerStepper_HasRequiredOps(step_mem->stepper);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepCreate",
                    "A required inner stepper function is NULL");
    MRIStepFree((void**) &ark_mem);
    return(NULL);
  }

  /* return ARKode memory */
  return((void*) ark_mem);
}


/*---------------------------------------------------------------
  MRIStepResize:

  This routine resizes the memory within the MRIStep module.
  It first resizes the main ARKode infrastructure memory, and
  then resizes its own data.
  ---------------------------------------------------------------*/
int MRIStepResize(void *arkode_mem, N_Vector y0, realtype t0,
                  ARKVecResizeFn resize, void *resize_data)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  SUNNonlinearSolver NLS;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepResize",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Determing change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL)
    N_VSpace(y0, &lrw1, &liw1);
  lrw_diff = lrw1 - ark_mem->lrw1;
  liw_diff = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* resize ARKode infrastructure memory (use hscale = 1.0) */
  retval = arkResize(ark_mem, y0, RCONST(1.0), t0, resize, resize_data);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::MRIStep", "MRIStepResize",
                    "Unable to resize main ARKode infrastructure");
    return(retval);
  }

  /* Resize Fse */
  if (step_mem->Fse) {
    if (!arkResizeVecArray(resize, resize_data,
                           step_mem->nstages_stored, y0, &(step_mem->Fse),
                           lrw_diff, &(ark_mem->lrw),
                           liw_diff, &(ark_mem->liw))) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepResize", "Unable to resize vector");
      return(ARK_MEM_FAIL);
    }
  }

  /* Resize Fsi */
  if (step_mem->Fsi) {
    if (!arkResizeVecArray(resize, resize_data,
                           step_mem->nstages_stored, y0, &(step_mem->Fsi),
                           lrw_diff, &(ark_mem->lrw),
                           liw_diff, &(ark_mem->liw))) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepResize", "Unable to resize vector");
      return(ARK_MEM_FAIL);
    }
  }

  /* Resize the nonlinear solver interface vectors (if applicable) */
  if (step_mem->sdata != NULL)
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                      liw_diff, y0, &step_mem->sdata)) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepResize", "Unable to resize vector");
      return(ARK_MEM_FAIL);
    }
  if (step_mem->zpred != NULL)
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                      liw_diff, y0, &step_mem->zpred)) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepResize", "Unable to resize vector");
      return(ARK_MEM_FAIL);
    }
  if (step_mem->zcor != NULL)
    if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                      liw_diff, y0, &step_mem->zcor)) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepResize", "Unable to resize vector");
      return(ARK_MEM_FAIL);
    }

  /* If a NLS object was previously used, destroy and recreate default Newton
     NLS object (can be replaced by user-defined object if desired) */
  if ((step_mem->NLS != NULL) && (step_mem->ownNLS)) {

    /* destroy existing NLS object */
    retval = SUNNonlinSolFree(step_mem->NLS);
    if (retval != ARK_SUCCESS)  return(retval);
    step_mem->NLS = NULL;
    step_mem->ownNLS = SUNFALSE;

    /* create new Newton NLS object */
    NLS = SUNNonlinSol_Newton(y0, ark_mem->sunctx);
    if (NLS == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepResize", "Error creating default Newton solver");
      return(ARK_MEM_FAIL);
    }

    /* attach new Newton NLS object to MRIStep */
    retval = MRIStepSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepResize", "Error attaching default Newton solver");
      return(ARK_MEM_FAIL);
    }
    step_mem->ownNLS = SUNTRUE;

  }

  /* Resize the inner stepper vectors */
  retval = mriStepInnerStepper_Resize(step_mem->stepper, resize, resize_data,
                                      lrw_diff, liw_diff, y0);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                    "MRIStepResize", "Unable to resize vector");
    return(ARK_MEM_FAIL);
  }

  /* reset nonlinear solver counters */
  if (step_mem->NLS != NULL)  step_mem->nsetups = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepReInit:

  This routine re-initializes the MRIStep module to solve a new
  problem of the same size as was previously solved (all counter
  values are set to 0).

  NOTE: the inner stepper needs to be reinitialized before
  calling this function.
  ---------------------------------------------------------------*/
int MRIStepReInit(void* arkode_mem, ARKRhsFn fse, ARKRhsFn fsi, realtype t0,
                  N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  SUNNonlinearSolver NLS;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepReInit",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode::MRIStep",
                    "MRIStepReInit", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check that at least one of fse, fsi is supplied and is to be used */
  if (fse == NULL && fsi == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepReInit", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit_rhs = (fse == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit_rhs = (fsi == NULL) ? SUNFALSE : SUNTRUE;

  /* Create a default Newton NLS object (just in case; will be deleted if
     the user attaches a nonlinear solver) */
  if (step_mem->implicit_rhs && !(step_mem->NLS)) {
    NLS = SUNNonlinSol_Newton(y0, ark_mem->sunctx);
    if (!NLS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepReInit", "Error creating default Newton solver");
      MRIStepFree((void**) &ark_mem); return(ARK_MEM_FAIL);
    }
    retval = MRIStepSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::MRIStep",
                      "MRIStepReInit", "Error attaching default Newton solver");
      MRIStepFree((void**) &ark_mem);  return(ARK_MEM_FAIL);
    }
    step_mem->ownNLS = SUNTRUE;
  }

  /* ReInitialize main ARKode infrastructure */
  retval = arkInit(arkode_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::MRIStep", "MRIStepReInit",
                    "Unable to reinitialize main ARKode infrastructure");
    return(retval);
  }

  /* Copy the input parameters into ARKode state */
  step_mem->fse = fse;
  step_mem->fsi = fsi;

  /* Initialize all the counters */
  step_mem->nfse      = 0;
  step_mem->nfsi      = 0;
  step_mem->nsetups   = 0;
  step_mem->nstlp     = 0;
  step_mem->nls_iters = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepReset:

  This routine resets the MRIStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
int MRIStepReset(void* arkode_mem, realtype tR, N_Vector yR)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepReset",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Initialize main ARKode infrastructure */
  retval = arkInit(ark_mem, tR, yR, RESET_INIT);

  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::MRIStep", "MRIStepReset",
                    "Unable to initialize main ARKode infrastructure");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSStolerances, MRIStepSVtolerances, MRIStepWFtolerances:

  These routines set integration tolerances (wrappers for general
  ARKode utility routines)
  ---------------------------------------------------------------*/
int MRIStepSStolerances(void *arkode_mem, realtype reltol, realtype abstol)
{
  /* unpack ark_mem, call arkSStolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSStolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSStolerances(ark_mem, reltol, abstol));
}

int MRIStepSVtolerances(void *arkode_mem, realtype reltol, N_Vector abstol)
{
  /* unpack ark_mem, call arkSVtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSVtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSVtolerances(ark_mem, reltol, abstol));
}

int MRIStepWFtolerances(void *arkode_mem, ARKEwtFn efun)
{
  /* unpack ark_mem, call arkWFtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepWFtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkWFtolerances(ark_mem, efun));
}


/*---------------------------------------------------------------
  MRIStepRootInit:

  Initialize (attach) a rootfinding problem to the stepper
  (wrappers for general ARKode utility routine)
  ---------------------------------------------------------------*/
int MRIStepRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)
{
  /* unpack ark_mem, call arkRootInit, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepRootInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkRootInit(ark_mem, nrtfn, g));
}


/*---------------------------------------------------------------
  MRIStepEvolve:

  This is the main time-integration driver (wrappers for general
  ARKode utility routine)
  ---------------------------------------------------------------*/
int MRIStepEvolve(void *arkode_mem, realtype tout, N_Vector yout,
                  realtype *tret, int itask)
{
  /* unpack ark_mem, call arkEvolve, and return */
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepEvolve", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkEvolve(ark_mem, tout, yout, tret, itask);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return(retval);
}


/*---------------------------------------------------------------
  MRIStepGetDky:

  This returns interpolated output of the solution or its
  derivatives over the most-recently-computed step (wrapper for
  generic ARKode utility routine)
  ---------------------------------------------------------------*/
int MRIStepGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  /* unpack ark_mem, call arkGetDky, and return */
  int retval;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepGetDky", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  SUNDIALS_MARK_FUNCTION_BEGIN(ARK_PROFILER);
  retval = arkGetDky(ark_mem, t, k, dky);
  SUNDIALS_MARK_FUNCTION_END(ARK_PROFILER);
  return(retval);
}

/*---------------------------------------------------------------
  MRIStepComputeState:

  Computes y based on the current prediction and given correction.
  ---------------------------------------------------------------*/
int MRIStepComputeState(void *arkode_mem, N_Vector zcor, N_Vector z)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepComputeState",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, z);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepFree frees all MRIStep memory, and then calls an ARKode
  utility routine to free the ARKode infrastructure memory.
  ---------------------------------------------------------------*/
void MRIStepFree(void **arkode_mem)
{
  sunindextype Cliw, Clrw;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL)  return;

  /* conditional frees on non-NULL MRIStep module */
  ark_mem = (ARKodeMem) (*arkode_mem);
  if (ark_mem->step_mem != NULL) {

    step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

    /* free the coupling structure and derived quantities */
    if (step_mem->MRIC != NULL) {
      MRIStepCoupling_Space(step_mem->MRIC, &Cliw, &Clrw);
      MRIStepCoupling_Free(step_mem->MRIC);
      step_mem->MRIC = NULL;
      ark_mem->liw -= Cliw;
      ark_mem->lrw -= Clrw;
      if (step_mem->stagetypes) {
        free(step_mem->stagetypes);
        step_mem->stagetypes = NULL;
        ark_mem->liw -= step_mem->stages;
      }
      if (step_mem->stage_map) {
        free(step_mem->stage_map);
        step_mem->stage_map = NULL;
        ark_mem->liw -= step_mem->stages;
      }
      if (step_mem->Ae_row) {
        free(step_mem->Ae_row);
        step_mem->Ae_row = NULL;
        ark_mem->lrw -= step_mem->stages;
      }
      if (step_mem->Ai_row) {
        free(step_mem->Ai_row);
        step_mem->Ai_row = NULL;
        ark_mem->lrw -= step_mem->stages;
      }
    }

    /* free the nonlinear solver memory (if applicable) */
    if ((step_mem->NLS != NULL) && (step_mem->ownNLS)) {
      SUNNonlinSolFree(step_mem->NLS);
      step_mem->ownNLS = SUNFALSE;
    }
    step_mem->NLS = NULL;

    /* free the linear solver memory */
    if (step_mem->lfree != NULL) {
      step_mem->lfree((void *) ark_mem);
      step_mem->lmem = NULL;
    }

    /* free the sdata, zpred and zcor vectors */
    if (step_mem->sdata != NULL) {
      arkFreeVec(ark_mem, &step_mem->sdata);
      step_mem->sdata = NULL;
    }
    if (step_mem->zpred != NULL) {
      arkFreeVec(ark_mem, &step_mem->zpred);
      step_mem->zpred = NULL;
    }
    if (step_mem->zcor != NULL) {
      arkFreeVec(ark_mem, &step_mem->zcor);
      step_mem->zcor = NULL;
    }

    /* free the RHS vectors */
    if (step_mem->Fse) {
      arkFreeVecArray(step_mem->nstages_stored, &(step_mem->Fse),
                      ark_mem->lrw1, &(ark_mem->lrw),
                      ark_mem->liw1, &(ark_mem->liw));
    }

    if (step_mem->Fsi) {
      arkFreeVecArray(step_mem->nstages_stored, &(step_mem->Fsi),
                      ark_mem->lrw1, &(ark_mem->lrw),
                      ark_mem->liw1, &(ark_mem->liw));
    }

    /* free the reusable arrays for fused vector interface */
    if (step_mem->cvals != NULL) {
      free(step_mem->cvals);
      step_mem->cvals = NULL;
      ark_mem->lrw -= (step_mem->nfusedopvecs);
    }
    if (step_mem->Xvecs != NULL) {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= (step_mem->nfusedopvecs);
    }
    step_mem->nfusedopvecs = 0;

    /* free the time stepper module itself */
    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;
  }

  /* free memory for overall ARKode infrastructure */
  arkFree(arkode_mem);
}


/*---------------------------------------------------------------
  MRIStepPrintMem:

  This routine outputs the memory from the MRIStep structure and
  the main ARKode infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
void MRIStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int i, retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepPrintMem",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return;

  /* if outfile==NULL, set it to stdout */
  if (outfile == NULL)  outfile = stdout;

  /* output data from main ARKode infrastructure */
  fprintf(outfile,"MRIStep Slow Stepper Mem:\n");
  arkPrintMem(ark_mem, outfile);

  /* output integer quantities */
  fprintf(outfile,"MRIStep: q = %i\n", step_mem->q);
  fprintf(outfile,"MRIStep: p = %i\n", step_mem->p);
  fprintf(outfile,"MRIStep: istage = %i\n", step_mem->istage);
  fprintf(outfile,"MRIStep: stages = %i\n", step_mem->stages);
  fprintf(outfile,"MRIStep: maxcor = %i\n", step_mem->maxcor);
  fprintf(outfile,"MRIStep: msbp = %i\n", step_mem->msbp);
  fprintf(outfile,"MRIStep: predictor = %i\n", step_mem->predictor);
  fprintf(outfile,"MRIStep: convfail = %i\n", step_mem->convfail);
  fprintf(outfile,"MRIStep: stagetypes =");
  for (i=0; i<step_mem->stages; i++)
    fprintf(outfile," %i",step_mem->stagetypes[i]);
  fprintf(outfile,"\n");

  /* output long integer quantities */
  fprintf(outfile,"MRIStep: nfse = %li\n", step_mem->nfse);
  fprintf(outfile,"MRIStep: nfsi = %li\n", step_mem->nfsi);
  fprintf(outfile,"MRIStep: nsetups = %li\n", step_mem->nsetups);
  fprintf(outfile,"MRIStep: nstlp = %li\n", step_mem->nstlp);
  fprintf(outfile,"MRIStep: nls_iters = %li\n", step_mem->nls_iters);

  /* output boolean quantities */
  fprintf(outfile,"MRIStep: user_linear = %i\n", step_mem->linear);
  fprintf(outfile,"MRIStep: user_linear_timedep = %i\n", step_mem->linear_timedep);
  fprintf(outfile,"MRIStep: user_explicit = %i\n", step_mem->explicit_rhs);
  fprintf(outfile,"MRIStep: user_implicit = %i\n", step_mem->implicit_rhs);
  fprintf(outfile,"MRIStep: jcur = %i\n", step_mem->jcur);
  fprintf(outfile,"MRIStep: ownNLS = %i\n", step_mem->ownNLS);

  /* output realtype quantities */
  fprintf(outfile,"MRIStep: Coupling structure:\n");
  MRIStepCoupling_Write(step_mem->MRIC, outfile);

  fprintf(outfile,"MRIStep: gamma = %"RSYM"\n", step_mem->gamma);
  fprintf(outfile,"MRIStep: gammap = %"RSYM"\n", step_mem->gammap);
  fprintf(outfile,"MRIStep: gamrat = %"RSYM"\n", step_mem->gamrat);
  fprintf(outfile,"MRIStep: crate = %"RSYM"\n", step_mem->crate);
  fprintf(outfile,"MRIStep: delp = %"RSYM"\n", step_mem->delp);
  fprintf(outfile,"MRIStep: eRNrm = %"RSYM"\n", step_mem->eRNrm);
  fprintf(outfile,"MRIStep: nlscoef = %"RSYM"\n", step_mem->nlscoef);
  fprintf(outfile,"MRIStep: crdown = %"RSYM"\n", step_mem->crdown);
  fprintf(outfile,"MRIStep: rdiv = %"RSYM"\n", step_mem->rdiv);
  fprintf(outfile,"MRIStep: dgmax = %"RSYM"\n", step_mem->dgmax);
  fprintf(outfile,"MRIStep: Ae_row =");
  for (i=0; i<step_mem->nstages_stored; i++)
    fprintf(outfile," %"RSYM,step_mem->Ae_row[i]);
  fprintf(outfile,"\n");
  fprintf(outfile,"MRIStep: Ai_row =");
  for (i=0; i<step_mem->nstages_stored; i++)
    fprintf(outfile," %"RSYM,step_mem->Ai_row[i]);
  fprintf(outfile,"\n");

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */
  fprintf(outfile, "MRIStep: sdata:\n");
  N_VPrintFile(step_mem->sdata, outfile);
  fprintf(outfile, "MRIStep: zpred:\n");
  N_VPrintFile(step_mem->zpred, outfile);
  fprintf(outfile, "MRIStep: zcor:\n");
  N_VPrintFile(step_mem->zcor, outfile);
  if (step_mem->Fse)
    for (i=0; i<step_mem->nstages_stored; i++) {
      fprintf(outfile,"MRIStep: Fse[%i]:\n", i);
      N_VPrintFile(step_mem->Fse[i], outfile);
    }
  if (step_mem->Fsi)
    for (i=0; i<step_mem->nstages_stored; i++) {
      fprintf(outfile,"MRIStep: Fsi[%i]:\n", i);
      N_VPrintFile(step_mem->Fsi[i], outfile);
    }
#endif

  /* print the inner stepper memory */
  mriStepInnerStepper_PrintMem(step_mem->stepper, outfile);

  return;
}



/*===============================================================
  MRIStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Interface routines supplied to ARKode
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  mriStep_AttachLinsol:

  This routine attaches the various set of system linear solver
  interface routines, data structure, and solver type to the
  MRIStep module.
  ---------------------------------------------------------------*/
int mriStep_AttachLinsol(void* arkode_mem, ARKLinsolInitFn linit,
                         ARKLinsolSetupFn lsetup,
                         ARKLinsolSolveFn lsolve,
                         ARKLinsolFreeFn lfree,
                         SUNLinearSolver_Type lsolve_type,
                         void *lmem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_AttachLinsol",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* free any existing system solver */
  if (step_mem->lfree != NULL)  step_mem->lfree(arkode_mem);

  /* Attach the provided routines, data structure and solve type */
  step_mem->linit  = linit;
  step_mem->lsetup = lsetup;
  step_mem->lsolve = lsolve;
  step_mem->lfree  = lfree;
  step_mem->lmem   = lmem;

  /* Reset all linear solver counters */
  step_mem->nsetups = 0;
  step_mem->nstlp   = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_DisableLSetup:

  This routine NULLifies the lsetup function pointer in the
  MRIStep module.
  ---------------------------------------------------------------*/
void mriStep_DisableLSetup(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_DisableLSetup",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return;

  /* nullify the lsetup function pointer */
  step_mem->lsetup = NULL;
}


/*---------------------------------------------------------------
  mriStep_GetLmem:

  This routine returns the system linear solver interface memory
  structure, lmem.
  ---------------------------------------------------------------*/
void* mriStep_GetLmem(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure, and return lmem */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_GetLmem",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(NULL);
  return(step_mem->lmem);
}


/*---------------------------------------------------------------
  mriStep_GetImplicitRHS:

  This routine returns the implicit RHS function pointer, fi.
  ---------------------------------------------------------------*/
ARKRhsFn mriStep_GetImplicitRHS(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure, and return fi */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_GetImplicitRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(NULL);
  if (step_mem->implicit_rhs) {
    return(step_mem->fsi);
  } else {
    return(NULL);
  }
}


/*---------------------------------------------------------------
  mriStep_GetGammas:

  This routine fills the current value of gamma, and states
  whether the gamma ratio fails the dgmax criteria.
  ---------------------------------------------------------------*/
int mriStep_GetGammas(void* arkode_mem, realtype *gamma,
                      realtype *gamrat, booleantype **jcur,
                      booleantype *dgamma_fail)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_GetGammas",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set outputs */
  *gamma  = step_mem->gamma;
  *gamrat = step_mem->gamrat;
  *jcur = &step_mem->jcur;
  *dgamma_fail = (SUNRabs(*gamrat - ONE) >= step_mem->dgmax);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  With initialization types FIRST_INIT this routine:
  - sets/checks the ARK Butcher tables to be used
  - allocates any memory that depends on the number of ARK
    stages, method order, or solver options
  - sets the call_fullrhs flag

  With other initialization types, this routine does nothing.
  ---------------------------------------------------------------*/
int mriStep_Init(void* arkode_mem, int init_type)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval, j;
  booleantype reset_efun;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_Init",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* immediately return if reset */
  if (init_type == RESET_INIT) return(ARK_SUCCESS);

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT) {

    /* enforce use of arkEwtSmallReal if using a fixed step size for
       an explicit method and an internal error weight function */
    reset_efun = SUNTRUE;
    if ( step_mem->implicit_rhs )  reset_efun = SUNFALSE;
    if ( ark_mem->user_efun )  reset_efun = SUNFALSE;
    if (reset_efun) {
      ark_mem->user_efun = SUNFALSE;
      ark_mem->efun      = arkEwtSetSmallReal;
      ark_mem->e_data    = ark_mem;
    }

    /* assume fixed outer step size */
    if (!ark_mem->fixedstep) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep", "mriStep_Init",
                      "Adaptive outer time stepping is not currently supported");
      return(ARK_ILL_INPUT);
    }

    /* Create coupling structure (if not already set) */
    retval = mriStep_SetCoupling(ark_mem);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep", "mriStep_Init",
                      "Could not create coupling table");
      return(ARK_ILL_INPUT);
    }

    /* Check that coupling structure is OK */
    retval = mriStep_CheckCoupling(ark_mem);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                      "mriStep_Init", "Error in coupling table");
      return(ARK_ILL_INPUT);
    }

    /* Retrieve/store method and embedding orders now that tables are finalized */
    step_mem->stages = step_mem->MRIC->stages;
    step_mem->q = step_mem->MRIC->q;
    step_mem->p = step_mem->MRIC->p;

    /* allocate/fill derived quantities from MRIC structure */

    /* stage map */
    if (step_mem->stage_map) {
      free(step_mem->stage_map);
      step_mem->stage_map = NULL;
      ark_mem->liw -= step_mem->stages;
    }
    step_mem->stage_map = (int *) calloc(step_mem->stages, sizeof(int));
    ark_mem->liw += step_mem->stages;

    retval = mriStepCoupling_GetStageMap(step_mem->MRIC,
                                         step_mem->stage_map,
                                         &(step_mem->nstages_stored));
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                      "mriStep_Init", "Error in coupling table");
      return(ARK_ILL_INPUT);
    }

    /* stage types */
    if (step_mem->stagetypes) {
      free(step_mem->stagetypes);
      step_mem->stagetypes = NULL;
      ark_mem->liw -= step_mem->stages;
    }
    step_mem->stagetypes = (int *) calloc(step_mem->stages, sizeof(int));
    ark_mem->liw += step_mem->stages;
    for (j=0; j<step_mem->stages; j++)
      step_mem->stagetypes[j] = mriStepCoupling_GetStageType(step_mem->MRIC, j);

    /* explicit RK coefficient row */
    if (step_mem->Ae_row) {
      free(step_mem->Ae_row);
      step_mem->Ae_row = NULL;
      ark_mem->lrw -= step_mem->stages;
    }
    step_mem->Ae_row = (realtype *) calloc(step_mem->stages,
                                           sizeof(realtype));
    ark_mem->lrw += step_mem->stages;

    /* implicit RK coefficient row */
    if (step_mem->Ai_row) {
      free(step_mem->Ai_row);
      step_mem->Ai_row = NULL;
      ark_mem->lrw -= step_mem->stages;
    }
    step_mem->Ai_row = (realtype *) calloc(step_mem->stages,
                                           sizeof(realtype));
    ark_mem->lrw += step_mem->stages;

    /* Allocate MRI RHS vector memory, update storage requirements */
    /*   Allocate Fse[0] ... Fse[nstages_stored - 1] if needed */
    if (step_mem->explicit_rhs) {
      if (!arkAllocVecArray(step_mem->nstages_stored,
                            ark_mem->ewt, &(step_mem->Fse),
                            ark_mem->lrw1, &(ark_mem->lrw),
                            ark_mem->liw1, &(ark_mem->liw)))
        return(ARK_MEM_FAIL);
    }

    /*   Allocate Fsi[0] ... Fsi[nstages_stored - 1] if needed */
    if (step_mem->implicit_rhs) {
      if (!arkAllocVecArray(step_mem->nstages_stored,
                            ark_mem->ewt, &(step_mem->Fsi),
                            ark_mem->lrw1, &(ark_mem->lrw),
                            ark_mem->liw1, &(ark_mem->liw)))
        return(ARK_MEM_FAIL);
    }

    /* if any slow stage is implicit, allocate sdata, zpred, zcor vectors;
       if all stages explicit, free default NLS object, and detach all
       linear solver routines.  Note: step_mem->implicit_rhs will only equal
       SUNTRUE if an implicit table has been user-provided. */
    if (step_mem->implicit_rhs) {
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->sdata)))
        return(ARK_MEM_FAIL);
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->zpred)))
        return(ARK_MEM_FAIL);
      if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->zcor)))
        return(ARK_MEM_FAIL);
    } else {
      if ((step_mem->NLS != NULL) && (step_mem->ownNLS)) {
        SUNNonlinSolFree(step_mem->NLS);
        step_mem->NLS = NULL;
        step_mem->ownNLS = SUNFALSE;
      }
      step_mem->linit  = NULL;
      step_mem->lsetup = NULL;
      step_mem->lsolve = NULL;
      step_mem->lfree  = NULL;
      step_mem->lmem   = NULL;
    }

    /* Allocate reusable arrays for fused vector interface */
    step_mem->nfusedopvecs = 2*step_mem->stages + 2;
    if (step_mem->cvals == NULL) {
      step_mem->cvals = (realtype *) calloc(step_mem->nfusedopvecs, sizeof(realtype));
      if (step_mem->cvals == NULL)  return(ARK_MEM_FAIL);
      ark_mem->lrw += (step_mem->nfusedopvecs);
    }
    if (step_mem->Xvecs == NULL) {
      step_mem->Xvecs = (N_Vector *) calloc(step_mem->nfusedopvecs, sizeof(N_Vector));
      if (step_mem->Xvecs == NULL)  return(ARK_MEM_FAIL);
      ark_mem->liw += (step_mem->nfusedopvecs);   /* pointers */
    }

    /* Allocate inner stepper data */
    retval = mriStepInnerStepper_AllocVecs(step_mem->stepper,
                                           step_mem->MRIC->nmat, ark_mem->ewt);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep", "mriStep_Init",
                      "Error allocating inner stepper memory");
      return(ARK_MEM_FAIL);
    }

    /* Limit interpolant degree based on method order (use negative
       argument to specify update instead of overwrite) */
    if (ark_mem->interp != NULL) {
      retval = arkInterpSetDegree(ark_mem, ark_mem->interp, -(step_mem->q-1));
      if (retval != ARK_SUCCESS) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep", "mriStep_Init",
                        "Unable to update interpolation polynomial degree");
        return(ARK_ILL_INPUT);
      }
    }

  }

  /* Call linit (if it exists) */
  if (step_mem->linit) {
    retval = step_mem->linit(ark_mem);
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_LINIT_FAIL, "ARKode::MRIStep", "mriStep_Init",
                      MSG_ARK_LINIT_FAIL);
      return(ARK_LINIT_FAIL);
    }
  }

  /* Initialize the nonlinear solver object (if it exists) */
  if (step_mem->NLS) {
    retval = mriStep_NlsInit(ark_mem);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_NLS_INIT_FAIL, "ARKode::MRIStep", "mriStep_Init",
                      "Unable to initialize SUNNonlinearSolver object");
      return(ARK_NLS_INIT_FAIL);
    }
  }

  /* Signal to shared arkode module that fullrhs is required after each step */
  ark_mem->call_fullrhs = SUNTRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS functions,
  f(t,y) = fse(t,y) + fsi(t,y)  + ff(t,y).

  This will be called in one of three 'modes':
    ARK_FULLRHS_START -> called at the beginning of a simulation
                         or after post processing at step
    ARK_FULLRHS_END   -> called at the end of a successful step
    ARK_FULLRHS_OTHER -> called elsewhere (e.g. for dense output)

  If it is called in ARK_FULLRHS_START mode, we store the vectors
  f(t,y) in F[0] for possible reuse in the first stage of the
  subsequent time step.

  If it is called in ARK_FULLRHS_END mode, we reevauate f(t,y). At
  this time no checks are made to see if the method coefficient
  support copying vectors F[stages] to fill f instead of calling f().

  ARK_FULLRHS_OTHER mode is only called for dense output in-between
  steps, so we strive to store the intermediate parts so that they
  do not interfere with the other two modes.
  ---------------------------------------------------------------*/
int mriStep_FullRHS(void* arkode_mem, realtype t, N_Vector y, N_Vector f,
                    int mode)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_FullRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* perform RHS functions contingent on 'mode' argument */
  switch(mode) {

  /* ARK_FULLRHS_START: called at the beginning of a simulation
     Store the vector fs(t,y) in F[0] for possible reuse
     in the first stage of the subsequent time step */
  case ARK_FULLRHS_START:

    /* call fse if the problem has an explicit component */
    if (step_mem->explicit_rhs) {
      retval = step_mem->fse(t, y, step_mem->Fse[0], ark_mem->user_data);
      step_mem->nfse++;
      if (retval != 0) {
         arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                         "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
         return(ARK_RHSFUNC_FAIL);
      }
    }

    /* call fsi if the problem has an implicit component */
    if (step_mem->implicit_rhs) {
      retval = step_mem->fsi(t, y, step_mem->Fsi[0], ark_mem->user_data);
      step_mem->nfsi++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                        "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* call ff (force new RHS computation) */
    retval = mriStepInnerStepper_FullRhs(step_mem->stepper, t, y, f,
                                         ARK_FULLRHS_OTHER);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                      "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* combine RHS vectors into output */
    if (step_mem->explicit_rhs && step_mem->implicit_rhs) { /* ImEx */
      N_VLinearSum(ONE, step_mem->Fse[0], ONE, f, f);
      N_VLinearSum(ONE, step_mem->Fsi[0], ONE, f, f);
    } else {
      if (step_mem->implicit_rhs) {         /* implicit */
        N_VLinearSum(ONE, step_mem->Fsi[0], ONE, f, f);
      } else {                          /* explicit */
        N_VLinearSum(ONE, step_mem->Fse[0], ONE, f, f);
      }
    }

    break;


  /* ARK_FULLRHS_END: called at the end of a successful step
     This always recomputes the full RHS (i.e., this is the
     same as case 0). */
  case ARK_FULLRHS_END:

    /* call fse if the problem has an explicit component */
    if (step_mem->explicit_rhs) {
      retval = step_mem->fse(t, y, step_mem->Fse[0], ark_mem->user_data);
      step_mem->nfse++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                        "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* call fsi if the problem has an implicit component */
    if (step_mem->implicit_rhs) {
      retval = step_mem->fsi(t, y, step_mem->Fsi[0], ark_mem->user_data);
      step_mem->nfsi++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                        "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* call ff (force new RHS computation) */
    retval = mriStepInnerStepper_FullRhs(step_mem->stepper, t, y, f,
                                         ARK_FULLRHS_OTHER);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                      "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* combine RHS vectors into output */
    if (step_mem->explicit_rhs && step_mem->implicit_rhs) { /* ImEx */
      N_VLinearSum(ONE, step_mem->Fse[0], ONE, f, f);
      N_VLinearSum(ONE, step_mem->Fsi[0], ONE, f, f);
    } else {
      if (step_mem->implicit_rhs) {         /* implicit */
        N_VLinearSum(ONE, step_mem->Fsi[0], ONE, f, f);
      } else {                          /* explicit */
        N_VLinearSum(ONE, step_mem->Fse[0], ONE, f, f);
      }
    }
    break;

  /*  ARK_FULLRHS_OTHER: called for dense output in-between steps
      store the intermediate calculations in such a way as to not
      interfere with the other two modes */
  case ARK_FULLRHS_OTHER:

    /* call fse if the problem has an explicit component (store in ark_tempv2) */
    if (step_mem->explicit_rhs) {
      retval = step_mem->fse(t, y, ark_mem->tempv2, ark_mem->user_data);
      step_mem->nfse++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                        "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }

    /* call fsi if the problem has an implicit component (store in sdata) */
    if (step_mem->implicit_rhs) {
      retval = step_mem->fsi(t, y, step_mem->sdata, ark_mem->user_data);
      step_mem->nfsi++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                        "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
    }


    /* call ff (force new RHS computation) */
    retval = mriStepInnerStepper_FullRhs(step_mem->stepper, t, y, f,
                                         ARK_FULLRHS_OTHER);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                      "mriStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* combine RHS vectors into output */
    if (step_mem->explicit_rhs && step_mem->implicit_rhs) { /* ImEx */
      N_VLinearSum(ONE, ark_mem->tempv2, ONE, f, f);
      N_VLinearSum(ONE, step_mem->sdata, ONE, f, f);
    } else {                   /* implicit */
      if (step_mem->implicit_rhs) {
        N_VLinearSum(ONE, step_mem->sdata, ONE, f, f);
      } else {                                           /* explicit */
        N_VLinearSum(ONE, ark_mem->tempv2, ONE, f, f);
      }
    }

    break;

  default:
    /* return with RHS failure if unknown mode is passed */
    arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::MRIStep",
                    "mriStep_FullRHS", "Unknown full RHS mode");
    return(ARK_RHSFUNC_FAIL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_TakeStep:

  This routine serves the primary purpose of the MRIStep module:
  it performs a single MRI step (with embedding, if possible).

  The output variable dsmPtr should contain estimate of the
  weighted local error if an embedding is present; otherwise it
  should be 0.

  The input/output variable nflagPtr is used to gauge convergence
  of any algebraic solvers within the step.  At the start of a new
  time step, this will initially have the value FIRST_CALL.  On
  return from this function, nflagPtr should have a value:
            0 => algebraic solve completed successfully
           >0 => solve did not converge at this step size
                 (but may with a smaller stepsize)
           <0 => solve encountered an unrecoverable failure

  The return value from this routine is:
            0 => step completed successfully
           >0 => step encountered recoverable failure;
                 reduce step and retry (if possible)
           <0 => step encountered unrecoverable failure
  ---------------------------------------------------------------*/
int mriStep_TakeStep(void* arkode_mem, realtype *dsmPtr, int *nflagPtr)
{
  ARKodeMem ark_mem;           /* outer ARKode memory        */
  ARKodeMRIStepMem step_mem;   /* outer stepper memory       */
  int is;                      /* current stage index        */
  int retval;                  /* reusable return flag       */

  /* initialize algebraic solver convergence flag to success;
     error estimate to zero */
  *nflagPtr = ARK_SUCCESS;
  *dsmPtr = ZERO;

  /* access the MRIStep mem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_TakeStep",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep step %li,  stage 0,  h = %"RSYM",  t_n = %"RSYM"\n",
         ark_mem->nst, ark_mem->h, ark_mem->tcur);
#endif

#ifdef SUNDIALS_DEBUG_PRINTVEC
    printf("    MRIStep slow stage 0 solution:\n");
    N_VPrint(ark_mem->ycur);
    if (step_mem->explicit_rhs)
    {
      printf("    MRIStep slow stage RHS Fe[0]:\n");
      N_VPrint(step_mem->Fse[0]);
    }
    if (step_mem->implicit_rhs)
    {
      printf("    MRIStep slow stage RHS Fi[0]:\n");
      N_VPrint(step_mem->Fsi[0]);
    }
#endif

  /* call nonlinear solver setup if it exists */
  if (step_mem->NLS)
    if ((step_mem->NLS)->ops->setup) {
      N_VConst(ZERO, ark_mem->tempv3);   /* set guess to 0 for predictor-corrector form */
      retval = SUNNonlinSolSetup(step_mem->NLS, ark_mem->tempv3, ark_mem);
      if (retval < 0) return(ARK_NLS_SETUP_FAIL);
      if (retval > 0) return(ARK_NLS_SETUP_RECVR);
    }

  /* The first stage is the previous time-step solution, so its RHS
     is the [already-computed] slow RHS from the start of the step */

  /* Loop over remaining stages */
  for (is = 1; is < step_mem->stages; is++) {

    /* Set current stage time  */
    ark_mem->tcur = ark_mem->tn + step_mem->MRIC->c[is]*ark_mem->h;

#ifdef SUNDIALS_DEBUG
    printf("    ------------------------------------"
           "----------------------------------------\n");
    printf("    MRIStep step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n",
           ark_mem->nst, is, ark_mem->h, ark_mem->tcur);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->report)
      fprintf(ark_mem->diagfp, "MRIStep  step  %li  %"RSYM"  %i  %"RSYM"\n",
              ark_mem->nst, ark_mem->h, is, ark_mem->tcur);

    /* Determine current stage type, and call corresponding routine; the
       vector ark_mem->ycur stores the previous stage solution on input, and
       should store the result of this stage solution on output. */
    switch (step_mem->stagetypes[is]) {
    case(MRISTAGE_ERK_FAST):
      retval = mriStep_StageERKFast(ark_mem, step_mem, is);
      break;
    case(MRISTAGE_ERK_NOFAST):
      retval = mriStep_StageERKNoFast(ark_mem, step_mem, is);
      break;
    case(MRISTAGE_DIRK_NOFAST):
      retval = mriStep_StageDIRKNoFast(ark_mem, step_mem, is, nflagPtr);
      break;
    case(MRISTAGE_DIRK_FAST):
      retval = mriStep_StageDIRKFast(ark_mem, step_mem, is, nflagPtr);
      break;
    }
    if (retval != ARK_SUCCESS)  return(retval);

#ifdef SUNDIALS_DEBUG_PRINTVEC
    printf("    MRIStep slow stage %i solution:\n",is);
    N_VPrint(ark_mem->ycur);
#endif

    /* apply user-supplied stage postprocessing function (if supplied) */
    if (ark_mem->ProcessStage != NULL) {
      retval = ark_mem->ProcessStage(ark_mem->tcur,
                                     ark_mem->ycur,
                                     ark_mem->user_data);
      if (retval != 0) return(ARK_POSTPROCESS_STAGE_FAIL);
    }

    /* conditionally reset the inner integrator with the modified stage solution */
    if ( (step_mem->stagetypes[is] != MRISTAGE_ERK_FAST) ||
         (ark_mem->ProcessStage != NULL) ) {
      retval = mriStepInnerStepper_Reset(step_mem->stepper,
                                         ark_mem->tcur, ark_mem->ycur);
      if (retval != ARK_SUCCESS)  return(ARK_INNERSTEP_FAIL);
    }

    /* Compute updated slow RHS except at last stage which is the new solution.
     * The new solution RHS evaluation happens in arkCompleteStep */
    if (is < step_mem->stages - 1 && step_mem->stage_map[is] > -1) {
      /* store explicit slow rhs */
      if (step_mem->explicit_rhs) {
        retval = step_mem->fse(ark_mem->tcur, ark_mem->ycur,
                               step_mem->Fse[step_mem->stage_map[is]],
                               ark_mem->user_data);
        step_mem->nfse++;
        if (retval < 0)  return(ARK_RHSFUNC_FAIL);
        if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);

#ifdef SUNDIALS_DEBUG_PRINTVEC
        printf("    MRIStep slow stage RHS Fe[%i]:\n",is);
        N_VPrint(step_mem->Fse[step_mem->stage_map[is]]);
#endif
      }

      /* store implicit slow rhs  */
      if (step_mem->implicit_rhs) {
        retval = step_mem->fsi(ark_mem->tcur, ark_mem->ycur,
                               step_mem->Fsi[step_mem->stage_map[is]],
                               ark_mem->user_data);
        step_mem->nfsi++;
        if (retval < 0)  return(ARK_RHSFUNC_FAIL);
        if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);

#ifdef SUNDIALS_DEBUG_PRINTVEC
        printf("    MRIStep slow stage RHS Fi[%i]:\n",is);
        N_VPrint(step_mem->Fsi[step_mem->stage_map[is]]);
#endif
      }
    } /* compute slow RHS */
  } /* loop over stages */

#ifdef SUNDIALS_DEBUG_PRINTVEC
    printf("    MRIStep updated solution:\n");
    N_VPrint(ark_mem->ycur);
#endif

  /* Solver diagnostics reporting */
  if (ark_mem->report)
    fprintf(ark_mem->diagfp, "MRIStep  etest  %li  %"RSYM"  %"RSYM"\n",
            ark_mem->nst, ark_mem->h, *dsmPtr);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  mriStep_AccessStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int mriStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeMRIStepMem *step_mem)
{

  /* access ARKodeMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    fname, MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem) arkode_mem;
  if ((*ark_mem)->step_mem==NULL) {
    arkProcessError(*ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    fname, MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *step_mem = (ARKodeMRIStepMem) (*ark_mem)->step_mem;
  return(ARK_SUCCESS);
}



/*---------------------------------------------------------------
  mriStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype mriStep_CheckNVector(N_Vector tmpl)
{
  if ( (tmpl->ops->nvclone     == NULL) ||
       (tmpl->ops->nvdestroy   == NULL) ||
       (tmpl->ops->nvlinearsum == NULL) ||
       (tmpl->ops->nvconst     == NULL) ||
       (tmpl->ops->nvscale     == NULL) ||
       (tmpl->ops->nvwrmsnorm  == NULL) )
    return(SUNFALSE);
  return(SUNTRUE);
}


/*---------------------------------------------------------------
  mriStep_SetCoupling

  This routine determines the MRI method to use, based on the
  desired accuracy.
  ---------------------------------------------------------------*/
int mriStep_SetCoupling(ARKodeMem ark_mem)
{
  ARKodeMRIStepMem step_mem;
  sunindextype Cliw, Clrw;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "mriStep_SetCoupling", MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* if coupling has already been specified, just return */
  if (step_mem->MRIC != NULL) return(ARK_SUCCESS);

  /* select method based on order and type */

  /**** ImEx methods ****/
  if (step_mem->implicit_rhs && step_mem->explicit_rhs) {

    switch (step_mem->q) {
    case 3:
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_IMEX_SD_3);
      break;
    case 4:
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_IMEX_SD_4);
      break;
    default:
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                      "mriStep_SetCoupling",
                      "No MRI method at requested order, using q=3.");
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_IMEX_SD_3);
      break;
    }

  /**** implicit methods ****/
  } else if (step_mem->implicit_rhs) {

    switch (step_mem->q) {
    case 2:
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_IMPL_SD_3);
      break;
    case 3:
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_IMPL_SD_3);
      break;
    case 4:
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_IMPL_SD_4);
      break;
    default:
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                      "mriStep_SetCoupling",
                      "No MRI method at requested order, using q=3.");
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_IMPL_SD_3);
      break;
    }

  /**** explicit methods ****/
  } else {

    switch (step_mem->q) {
    case 3:
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_EXPL_3);
      break;
    case 4:
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_EXPL_4);
      break;
    default:
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                      "mriStep_SetCoupling",
                      "No MRI method at requested order, using q=3.");
      step_mem->MRIC = MRIStepCoupling_LoadTable(MRISTEP_DEFAULT_EXPL_3);
      break;
    }

  }

  if (step_mem->MRIC == NULL) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                    "mriStep_SetCoupling",
                    "An error occurred in constructing coupling table.");
    return(ARK_INVALID_TABLE);
  }

  /* note coupling structure space requirements */
  MRIStepCoupling_Space(step_mem->MRIC, &Cliw, &Clrw);
  ark_mem->liw += Cliw;
  ark_mem->lrw += Clrw;

  /* set [redundant] stored values for stage numbers and
     method/embedding orders */
  step_mem->stages = step_mem->MRIC->stages;
  step_mem->q = step_mem->MRIC->q;
  step_mem->p = step_mem->MRIC->p;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_CheckCoupling

  This routine runs through the MRI coupling structure to ensure
  that it meets all necessary requirements, including:
    sorted abscissae, with c[0] = 0 and c[end] = 1
    lower-triangular (i.e., ERK or DIRK)
    all DIRK stages are solve-decoupled [temporarily]
    method order q > 0 (all)
    stages > 0 (all)

  Returns ARK_SUCCESS if it passes, ARK_INVALID_TABLE otherwise.
  ---------------------------------------------------------------*/
int mriStep_CheckCoupling(ARKodeMem ark_mem)
{
  int i, j, k;
  booleantype okay;
  ARKodeMRIStepMem step_mem;
  realtype Gabs, Wabs;
  const realtype tol = RCONST(100.0)*UNIT_ROUNDOFF;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "mriStep_CheckCoupling", MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* check that stages > 0 */
  if (step_mem->MRIC->stages < 1) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                    "mriStep_CheckCoupling", "stages < 1!");
    return(ARK_INVALID_TABLE);
  }

  /* check that method order q > 0 */
  if (step_mem->MRIC->q < 1) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                    "mriStep_CheckCoupling", "method order < 1");
    return(ARK_INVALID_TABLE);
  }

  /* check that embedding order p > 0 (if adaptive) */
  if ((step_mem->MRIC->p < 1) && (!ark_mem->fixedstep)) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                    "mriStep_CheckCoupling", "embedding order < 1");
    return(ARK_INVALID_TABLE);
  }

  /* Check that the matrices are defined appropriately */
  if (step_mem->implicit_rhs && step_mem->explicit_rhs) {
    /* ImEx */
    if (!(step_mem->MRIC->W) || !(step_mem->MRIC->G)) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                      "mriStep_CheckCoupling",
                      "Invalid coupling table for an IMEX problem!");
      return(ARK_ILL_INPUT);
    }
  } else if (step_mem->explicit_rhs) {
    /* Explicit */
    if (!(step_mem->MRIC->W) || step_mem->MRIC->G) {
      arkProcessError(ark_mem,ARK_ILL_INPUT, "ARKode::MRIStep",
                      "mriStep_CheckCoupling",
                      "Invalid coupling table for an explicit problem!");
      return(ARK_ILL_INPUT);
    }
  } else {
    /* Implicit */
    if (step_mem->MRIC->W || !(step_mem->MRIC->G)) {
      arkProcessError(ark_mem,ARK_ILL_INPUT, "ARKode::MRIStep",
                      "mriStep_CheckCoupling",
                      "Invalid coupling table fro an implicit problem!");
      return(ARK_ILL_INPUT);
    }
  }

  /* Check that W tables are strictly lower triangular */
  if (step_mem->MRIC->W) {
    Wabs = RCONST(0.0);
    for (k=0; k<step_mem->MRIC->nmat; k++)
      for (i=0; i<step_mem->MRIC->stages; i++)
        for (j=i; j<step_mem->MRIC->stages; j++)
          Wabs += SUNRabs(step_mem->MRIC->W[k][i][j]);
    if (Wabs > tol) {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                      "mriStep_CheckCoupling",
                      "Coupling can be up to ERK (at most)!");
      return(ARK_INVALID_TABLE);
    }
  }

  /* Check that G tables are lower triangular */
  if (step_mem->MRIC->G) {
    Gabs = RCONST(0.0);
    for (k=0; k<step_mem->MRIC->nmat; k++)
      for (i=0; i<step_mem->MRIC->stages; i++)
        for (j=i+1; j<step_mem->MRIC->stages; j++)
          Gabs += SUNRabs(step_mem->MRIC->G[k][i][j]);
    if (Gabs > tol) {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                      "mriStep_CheckCoupling",
                      "Coupling can be up to DIRK (at most)!");
      return(ARK_INVALID_TABLE);
    }
  }

  /* Check that no stage has MRISTAGE_DIRK_FAST type (for now) */
  okay = SUNTRUE;
  for (i=0; i<step_mem->MRIC->stages; i++)
    if (mriStepCoupling_GetStageType(step_mem->MRIC, i) == MRISTAGE_DIRK_FAST)
      okay = SUNFALSE;
  if (!okay) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                    "mriStep_CheckCoupling",
                    "solve-coupled DIRK stages not currently supported");
    return(ARK_INVALID_TABLE);
  }

  /* check that stage times are sorted */
  okay = SUNTRUE;
  for (i=1; i<step_mem->MRIC->stages; i++) {
    if ((step_mem->MRIC->c[i] - step_mem->MRIC->c[i-1]) < -tol)
      okay = SUNFALSE;
  }
  if (!okay) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                    "mriStep_CheckCoupling",
                    "Stage times must be sorted.");
    return(ARK_INVALID_TABLE);
  }

  /* check that the first stage is just the old step solution */
  Gabs = SUNRabs(step_mem->MRIC->c[0]);
  for (k=0; k<step_mem->MRIC->nmat; k++)
    for (j=0; j<step_mem->MRIC->stages; j++) {
      if (step_mem->MRIC->W)
        Gabs += SUNRabs(step_mem->MRIC->W[k][0][j]);
      if (step_mem->MRIC->G)
        Gabs += SUNRabs(step_mem->MRIC->G[k][0][j]);
    }
  if (Gabs > tol) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                    "mriStep_CheckCoupling",
                    "First stage must equal old solution.");
    return(ARK_INVALID_TABLE);
  }

  /* check that the last stage is at the final time */
  if (SUNRabs(ONE - step_mem->MRIC->c[step_mem->MRIC->stages-1]) > tol) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                    "mriStep_CheckCoupling",
                    "Final stage time must be equal 1.");
    return(ARK_INVALID_TABLE);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_StageERKFast

  This routine performs a single MRI stage with explicit slow
  time scale and fast time scale that requires evolution.
  ---------------------------------------------------------------*/
int mriStep_StageERKFast(ARKodeMem ark_mem,
                         ARKodeMRIStepMem step_mem, int is)
{
  realtype cdiff;    /* stage time increment */
  realtype t0;       /* start time for stage */
  int retval;        /* reusable return flag */

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep ERK fast stage\n");
#endif

  /* Set initial time for fast evolution */
  t0 = ark_mem->tn + step_mem->MRIC->c[is-1]*ark_mem->h;

  /* compute the inner forcing */
  cdiff = step_mem->MRIC->c[is] - step_mem->MRIC->c[is-1];
  retval = mriStep_ComputeInnerForcing(step_mem, is, cdiff);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set inner forcing time normalization constants */
  step_mem->stepper->tshift = t0;
  step_mem->stepper->tscale = cdiff * ark_mem->h;

  /* pre inner evolve function (if supplied) */
  if (step_mem->pre_inner_evolve) {
    retval = step_mem->pre_inner_evolve(t0, step_mem->stepper->forcing,
                                        step_mem->stepper->nforcing,
                                        ark_mem->user_data);
    if (retval != 0) return(ARK_OUTERTOINNER_FAIL);
  }

  /* advance inner method in time */
  retval = mriStepInnerStepper_Evolve(step_mem->stepper, t0, ark_mem->tcur,
                                      ark_mem->ycur);
  if (retval < 0) return(ARK_INNERSTEP_FAIL);

  /* post inner evolve function (if supplied) */
  if (step_mem->post_inner_evolve) {
    retval = step_mem->post_inner_evolve(ark_mem->tcur, ark_mem->ycur,
                                         ark_mem->user_data);
    if (retval != 0) return(ARK_INNERTOOUTER_FAIL);
  }

  /* return with success */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_StageERKNoFast

  This routine performs a single MRI stage with explicit slow
  time scale only (no fast time scale evolution).
  ---------------------------------------------------------------*/
int mriStep_StageERKNoFast(ARKodeMem ark_mem,
                           ARKodeMRIStepMem step_mem, int is)
{
  int retval, j, nvec;

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep ERK stage\n");
#endif

  /* determine effective ERK coefficients (store in cvals) */
  retval = mriStep_RKCoeffs(step_mem->MRIC, is, step_mem->stage_map,
                            step_mem->Ae_row, step_mem->Ai_row);
  if (retval != ARK_SUCCESS) { return(retval); }

  /* call fused vector operation to perform ERK update */
  step_mem->cvals[0] = ONE;
  step_mem->Xvecs[0] = ark_mem->ycur;
  nvec = 1;
  for (j = 0; j < is; j++) {
    if (step_mem->explicit_rhs && step_mem->stage_map[j] > -1) {
      step_mem->cvals[nvec] = ark_mem->h *
        step_mem->Ae_row[step_mem->stage_map[j]];
      step_mem->Xvecs[nvec] = step_mem->Fse[step_mem->stage_map[j]];
      nvec += 1;
    }
    if (step_mem->implicit_rhs && step_mem->stage_map[j] > -1) {
      step_mem->cvals[nvec] = ark_mem->h *
        step_mem->Ai_row[step_mem->stage_map[j]];
      step_mem->Xvecs[nvec] = step_mem->Fsi[step_mem->stage_map[j]];
      nvec += 1;
    }
  }
  /* Is there a case where we have an explicit update with Fsi? */

  retval = N_VLinearCombination(nvec, step_mem->cvals,
                                step_mem->Xvecs, ark_mem->ycur);
  if (retval != 0) return(ARK_VECTOROP_ERR);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_StageDIRKFast

  This routine performs a single stage of a "solve coupled"
  MRI method, i.e. a stage that is DIRK on the slow time scale
  and involves evolution of the fast time scale, in a
  fully-coupled fashion.
  ---------------------------------------------------------------*/
int mriStep_StageDIRKFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                          int is, int *nflagPtr)
{
#ifdef SUNDIALS_DEBUG
  printf("    MRIStep DIRK fast stage\n");
#endif

  /* this is not currently implemented */
  arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::MRIStep",
                  "mriStep_StageDIRKFast",
                  "This routine is not yet implemented.");
  return(ARK_INVALID_TABLE);
}


/*---------------------------------------------------------------
  mriStep_StageDIRKNoFast

  This routine performs a single MRI stage with implicit slow
  time scale only (no fast time scale evolution).
  ---------------------------------------------------------------*/
int mriStep_StageDIRKNoFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                            int is, int *nflagPtr)
{
  int retval;

#ifdef SUNDIALS_DEBUG
  printf("    MRIStep DIRK stage\n");
#endif

  /* store current stage index */
  step_mem->istage = is;

  /* Call predictor for current stage solution (result placed in zpred) */
  retval = mriStep_Predict(ark_mem, is, step_mem->zpred);
  if (retval != ARK_SUCCESS)  return (retval);

  /* If a user-supplied predictor routine is provided, call that here
     Note that mriStep_Predict is *still* called, so this user-supplied
     routine can just 'clean up' the built-in prediction, if desired. */
  if (step_mem->stage_predict) {
    retval = step_mem->stage_predict(ark_mem->tcur, step_mem->zpred,
                                     ark_mem->user_data);
    if (retval < 0)  return(ARK_USER_PREDICT_FAIL);
    if (retval > 0)  return(TRY_AGAIN);
  }

#ifdef SUNDIALS_DEBUG_PRINTVEC
  printf("    MRIStep predictor:\n");
  N_VPrint(step_mem->zpred);
#endif

  /* determine effective DIRK coefficients (store in cvals) */
  retval = mriStep_RKCoeffs(step_mem->MRIC, is, step_mem->stage_map,
                            step_mem->Ae_row, step_mem->Ai_row);
  if (retval != ARK_SUCCESS) { return(retval); }

  /* Set up data for evaluation of DIRK stage residual (data stored in sdata) */
  retval = mriStep_StageSetup(ark_mem);
  if (retval != ARK_SUCCESS)  return (retval);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  printf("    MRIStep rhs data:\n");
  N_VPrint(step_mem->sdata);
#endif

  /* perform implicit solve (result is stored in ark_mem->ycur); return
     with positive value on anything but success */
  *nflagPtr = mriStep_Nls(ark_mem, *nflagPtr);
  if (*nflagPtr != ARK_SUCCESS)  return(TRY_AGAIN);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_ComputeInnerForcing

  Constructs the 'coefficient' vectors for the forcing polynomial
  for a 'fast' outer MRI stage i:

  p_i(theta) = sum_{k=0}^{n-1} forcing[k] * theta^k

  where theta = (t - t0) / (cdiff*h) is the mapped 'time' for
  each 'fast' MRIStep evolution, with:
  * t0 -- the start of this outer MRIStep stage
  * cdiff*h, the temporal width of this MRIStep stage
  * n -- shorthand for MRIC->nmat

  explicit and solve-decoupled implicit or IMEX MRI-based methods
  define this forcing polynomial for each outer stage i > 0:

  p_i(theta) = w_i,0(theta) * fse_0 + ... + w_i,{i-1}(theta) * fse_{i-1}
             + g_i,0(theta) * fsi_0 + ... + g_i,{i-1}(theta) * fsi_{i-1}

  where

  w_i,j(theta) = w_0,i,j + w_1,i,j * theta + ... + w_n,i,j * theta^{n-1},
  w_k,i,j = 1/cdiff * MRIC->W[k][i][j]

  and

  g_i,j(theta) = g_0,i,j + g_1,i,j * theta + ... + g_n,i,j * theta^{n-1},
  g_k,i,j = 1/cdiff * MRIC->G[k][i][j]

  Converting to the appropriate form, we have

  p_i(theta) = ( w_0,i,0 * fse_0 + ... + w_0,i,{i-1} * fse_{i-1} +
                 g_0,i,0 * fsi_0 + ... + g_0,i,{i-1} * fsi_{i-1} ) * theta^0
             + ( w_1,i,0 * fse_0 + ... + w_1,i,{i-1} * fse_{i-1} +
                 g_1,i,0 * fsi_0 + ... + g_1,i,{i-1} * fsi_{i-1} ) * theta^1
                                    .
                                    .
                                    .
             + ( w_n,i,0 * fse_0 + ... + w_n,i,{i-1} * fse_{i-1} +
                 g_n,i,0 * fsi_0 + ... + g_n,i,{i-1} * fsi_{i-1} ) * theta^{n-1}

  Thus we define the forcing vectors for k = 0,...,nmat - 1

  forcing[k] = w_k,i,0 * fse_0 + ... + w_k,i,{i-1} * fse_{i-1}
             + g_k,i,0 * fsi_0 + ... + g_k,i,{i-1} * fsi_{i-1}

             = 1 / cdiff *
               ( W[k][i][0] * fse_0 + ... + W[k][i][i-1] * fse_{i-1} +
               ( G[k][i][0] * fsi_0 + ... + G[k][i][i-1] * fsi_{i-1} )

  This routine additionally returns a success/failure flag:
     ARK_SUCCESS -- successful evaluation
  ---------------------------------------------------------------*/

int mriStep_ComputeInnerForcing(ARKodeMRIStepMem step_mem,
                                int stage, realtype cdiff)
{
  realtype  rcdiff;
  int       j, k, nmat, nstore, retval;
  realtype* cvals;
  N_Vector* Xvecs;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* compute inner forcing vectors (assumes cdiff != 0) */
  nstore = 0;
  for (j = 0; j < stage; j++) {
    if (step_mem->explicit_rhs && step_mem->stage_map[j] > -1) {
      Xvecs[nstore] = step_mem->Fse[step_mem->stage_map[j]];
      nstore += 1;
    }
    if (step_mem->implicit_rhs && step_mem->stage_map[j] > -1) {
      Xvecs[nstore] = step_mem->Fsi[step_mem->stage_map[j]];
      nstore += 1;
    }
  }

  nmat   = step_mem->MRIC->nmat;
  rcdiff = ONE / cdiff;

  for (k = 0; k < nmat; k++) {
    nstore = 0;
    for (j = 0; j < stage; j++) {
      if (step_mem->stage_map[j] > -1) {
        if (step_mem->explicit_rhs && step_mem->implicit_rhs) {
          /* ImEx */
          cvals[nstore] = rcdiff * step_mem->MRIC->W[k][stage][j];
          nstore += 1;
          cvals[nstore] = rcdiff * step_mem->MRIC->G[k][stage][j];
          nstore += 1;
        } else if (step_mem->explicit_rhs) {
          /* explicit only */
          cvals[nstore] = rcdiff * step_mem->MRIC->W[k][stage][j];
          nstore += 1;
        } else {
          /* implicit only */
          cvals[nstore] = rcdiff * step_mem->MRIC->G[k][stage][j];
          nstore += 1;
        }
      }
    }

    retval = N_VLinearCombination(nstore, cvals, Xvecs,
                                  step_mem->stepper->forcing[k]);
    if (retval != 0) return(ARK_VECTOROP_ERR);
  }

#ifdef SUNDIALS_DEBUG_PRINTVEC
  for (k = 0; k < nmat; k++) {
    printf("    MRIStep forcing[%i]:\n", k);
    N_VPrint(step_mem->stepper->forcing[k]);
  }
#endif

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Compute/return the 'effective' RK coefficients for a 'nofast'
  stage.  It is assumed that the array 'A' has already been
  allocated to have length MRIC->stages.
  ---------------------------------------------------------------*/

int mriStep_RKCoeffs(MRIStepCoupling MRIC, int is, int *stage_map,
                     realtype *Ae_row, realtype *Ai_row)
{
  int j, k;
  realtype kconst;

  if (is < 1 || is >= MRIC->stages || !stage_map || !Ae_row || !Ai_row)
    return ARK_INVALID_TABLE;

  /* initialize RK coefficient array */
  for (j = 0; j < MRIC->stages; j++) {
    Ae_row[j] = ZERO;
    Ai_row[j] = ZERO;
  }

  /* compute RK coefficients */
  for (k = 0; k < MRIC->nmat; k++) {
    kconst = ONE/(k+ONE);
    if (MRIC->W) {
      for (j = 0; j < is; j++)
        if (stage_map[j] > -1)
          Ae_row[stage_map[j]] += (MRIC->W[k][is][j] * kconst);
    }
    if (MRIC->G) {
      for (j = 0; j <= is; j++)
        if (stage_map[j] > -1)
          Ai_row[stage_map[j]] += (MRIC->G[k][is][j] * kconst);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_Predict

  This routine computes the prediction for a specific internal
  stage solution, storing the result in yguess.  The
  prediction is done using the interpolation structure in
  extrapolation mode, hence stages "far" from the previous time
  interval are predicted using lower order polynomials than the
  "nearby" stages.
  ---------------------------------------------------------------*/
int mriStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess)
{
  int i, retval, jstage, nvec;
  realtype tau;
  realtype h;
  ARKodeMRIStepMem step_mem;
  realtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "mriStep_Predict", MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* verify that interpolation structure is provided */
  if ((ark_mem->interp == NULL) && (step_mem->predictor > 0)) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "mriStep_Predict",
                    "Interpolation structure is NULL");
    return(ARK_MEM_NULL);
  }

  /* local shortcuts for use with fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* if the first step (or if resized), use initial condition as guess */
  if (ark_mem->initsetup) {
    N_VScale(ONE, ark_mem->yn, yguess);
    return(ARK_SUCCESS);
  }

  /* set evaluation time tau as relative shift from previous successful time */
  tau = step_mem->MRIC->c[istage]*ark_mem->h/ark_mem->hold;

  /* use requested predictor formula */
  switch (step_mem->predictor) {

  case 1:

    /***** Interpolatory Predictor 1 -- all to max order *****/
    retval = arkPredict_MaximumOrder(ark_mem, tau, yguess);
    if (retval != ARK_ILL_INPUT)  return(retval);
    break;

  case 2:

    /***** Interpolatory Predictor 2 -- decrease order w/ increasing level of extrapolation *****/
    retval = arkPredict_VariableOrder(ark_mem, tau, yguess);
    if (retval != ARK_ILL_INPUT)  return(retval);
    break;

  case 3:

    /***** Cutoff predictor: max order interpolatory output for stages "close"
           to previous step, first-order predictor for subsequent stages *****/
    retval = arkPredict_CutoffOrder(ark_mem, tau, yguess);
    if (retval != ARK_ILL_INPUT)  return(retval);
    break;

  case 4:

    /***** Bootstrap predictor: if any previous stage in step has nonzero c_i,
           construct a quadratic Hermite interpolant for prediction; otherwise
           use the trivial predictor.  The actual calculations are performed in
           arkPredict_Bootstrap, but here we need to determine the appropriate
           stage, c_j, to use. *****/

    /* determine if any previous stages in step meet criteria */
    jstage = -1;
    for (i=0; i<istage; i++)
      jstage = (step_mem->MRIC->c[i] != ZERO) ? i : jstage;

    /* if using the trivial predictor, break */
    if (jstage == -1)  break;

    /* find the "optimal" previous stage to use */
    for (i=0; i<istage; i++)
      if ( (step_mem->MRIC->c[i] > step_mem->MRIC->c[jstage]) &&
           (step_mem->MRIC->c[i] != ZERO) && step_mem->stage_map[i] > -1)
        jstage = i;

    /* set stage time, stage RHS and interpolation values */
    h = ark_mem->h * step_mem->MRIC->c[jstage];
    tau = ark_mem->h * step_mem->MRIC->c[istage];
    nvec = 0;
    if (step_mem->implicit_rhs) {    /* Implicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fsi[step_mem->stage_map[jstage]];
      nvec += 1;
    }
    if (step_mem->explicit_rhs) {    /* Explicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fse[step_mem->stage_map[jstage]];
      nvec += 1;
    }

    /* call predictor routine */
    retval = arkPredict_Bootstrap(ark_mem, h, tau, nvec, cvals, Xvecs, yguess);
    if (retval != ARK_ILL_INPUT)  return(retval);
    break;

  }

  /* if we made it here, use the trivial predictor (previous step solution) */
  N_VScale(ONE, ark_mem->yn, yguess);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_StageSetup

  This routine sets up the stage data for computing the
  solve-decoupled MRI stage residual, along with the step- and
  method-related factors gamma, gammap and gamrat.

  At the ith stage, we compute the residual vector for
  z=z_i=zp+zc:
    r = z - z_{i-1} - h*sum_{j=0}^{i} A(i,j)*F(z_j)
    r = (zp + zc) - z_{i-1} - h*sum_{j=0}^{i} A(i,j)*F(z_j)
    r = (zc - gamma*F(z)) - data,
  where data = (z_{i-1} - zp + h*sum_{j=0}^{i-1} A(i,j)*F(z_j))
  corresponds to existing information.  This routine computes
  this 'data' vector and stores in step_mem->sdata.

  Note: on input, this row A(i,:) is already stored in rkcoeffs.
  ---------------------------------------------------------------*/
int mriStep_StageSetup(ARKodeMem ark_mem)
{
  /* local data */
  ARKodeMRIStepMem step_mem;
  int retval, i, j, nvec;
  realtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "mriStep_StageSetup", MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* Set shortcut to current stage index */
  i = step_mem->istage;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* Update gamma (if the method contains an implicit component) */
  step_mem->gamma = ark_mem->h * step_mem->Ai_row[step_mem->stage_map[i]];

  if (ark_mem->firststage)
    step_mem->gammap = step_mem->gamma;
  step_mem->gamrat = (ark_mem->firststage) ?
    ONE : step_mem->gamma / step_mem->gammap;  /* protect x/x != 1.0 */

  /* set cvals and Xvecs for setting stage data */
  cvals[0] = ONE;
  Xvecs[0] = ark_mem->ycur;
  cvals[1] = -ONE;
  Xvecs[1] = step_mem->zpred;
  nvec     = 2;

  for (j = 0; j < i; j++) {
    if (step_mem->explicit_rhs && step_mem->stage_map[j] > -1) {
      cvals[nvec] = ark_mem->h * step_mem->Ae_row[step_mem->stage_map[j]];
      Xvecs[nvec] = step_mem->Fse[step_mem->stage_map[j]];
      nvec += 1;
    }
    if (step_mem->implicit_rhs && step_mem->stage_map[j] > -1) {
      cvals[nvec] = ark_mem->h * step_mem->Ai_row[step_mem->stage_map[j]];
      Xvecs[nvec] = step_mem->Fsi[step_mem->stage_map[j]];
      nvec += 1;
    }
  }

  /* call fused vector operation to do the work */
  retval = N_VLinearCombination(nvec, cvals, Xvecs, step_mem->sdata);
  if (retval != 0) return(ARK_VECTOROP_ERR);

  /* return with success */
  return (ARK_SUCCESS);
}


/*---------------------------------------------------------------
  User-callable functions for a custom inner integrator
  ---------------------------------------------------------------*/


int MRIStepInnerStepper_Create(SUNContext sunctx, MRIStepInnerStepper *stepper)
{
  if (!sunctx) return ARK_ILL_INPUT;

  *stepper = NULL;
  *stepper = (MRIStepInnerStepper) malloc(sizeof(**stepper));
  if (*stepper == NULL) {
    arkProcessError(NULL, ARK_MEM_FAIL, "ARKode::MRIStep",
                    "MRIStepInnerStepper_Create",
                    MSG_ARK_ARKMEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  memset(*stepper, 0, sizeof(**stepper));

  (*stepper)->ops =
    (MRIStepInnerStepper_Ops) malloc(sizeof(*((*stepper)->ops)));
  if ((*stepper)->ops == NULL) {
    arkProcessError(NULL, ARK_MEM_FAIL, "ARKode::MRIStep",
                    "MRIStepInnerStepper_Create",
                    MSG_ARK_ARKMEM_FAIL);
    free(*stepper);
    return(ARK_MEM_FAIL);
  }
  memset((*stepper)->ops, 0, sizeof(*((*stepper)->ops)));

  /* initialize stepper data */
  (*stepper)->last_flag = ARK_SUCCESS;
  (*stepper)->sunctx    = sunctx;

  return(ARK_SUCCESS);
}


int MRIStepInnerStepper_Free(MRIStepInnerStepper *stepper)
{
  if (*stepper == NULL) return ARK_SUCCESS;

  /* free the inner forcing and fused op workspace vector */
  mriStepInnerStepper_FreeVecs(*stepper);

  /* free inner stepper mem */
  free(*stepper);
  *stepper = NULL;

  return(ARK_SUCCESS);
}


int MRIStepInnerStepper_SetContent(MRIStepInnerStepper stepper,
                                   void *content)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_SetContent",
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }
  stepper->content = content;

  return ARK_SUCCESS;
}


int MRIStepInnerStepper_GetContent(MRIStepInnerStepper stepper,
                                   void **content)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_GetContent",
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }
  *content = stepper->content;

  return ARK_SUCCESS;
}


int MRIStepInnerStepper_SetEvolveFn(MRIStepInnerStepper stepper,
                                    MRIStepInnerEvolveFn fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_SetEvolveFn",
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_SetEvolveFn",
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->evolve = fn;

  return ARK_SUCCESS;
}


int MRIStepInnerStepper_SetFullRhsFn(MRIStepInnerStepper stepper,
                                     MRIStepInnerFullRhsFn fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_SetFullRhsFn",
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_SetFullRhsFn",
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->fullrhs = fn;

  return ARK_SUCCESS;
}


int MRIStepInnerStepper_SetResetFn(MRIStepInnerStepper stepper,
                                   MRIStepInnerResetFn fn)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_SetResetFn",
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  if (stepper->ops == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_SetResetFn",
                    "Inner stepper operations structure is NULL");
    return ARK_ILL_INPUT;
  }

  stepper->ops->reset = fn;

  return ARK_SUCCESS;
}


int MRIStepInnerStepper_AddForcing(MRIStepInnerStepper stepper,
                                   realtype t, N_Vector f)
{
  realtype tau, taui;
  int i;

  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_AddForcing",
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  /* always append the constant forcing term */
  stepper->vals[0] = ONE;
  stepper->vecs[0] = f;

  /* compute normalized time tau and initialize tau^i */
  tau  = (t - stepper->tshift) / (stepper->tscale);
  taui = ONE;

  for (i = 0; i < stepper->nforcing; i++) {
    stepper->vals[i+1] = taui;
    stepper->vecs[i+1] = stepper->forcing[i];
    taui *= tau;
  }

  N_VLinearCombination(stepper->nforcing + 1,
                       stepper->vals,
                       stepper->vecs,
                       f);

  return ARK_SUCCESS;
}


int MRIStepInnerStepper_GetForcingData(MRIStepInnerStepper stepper,
                                       realtype *tshift, realtype *tscale,
                                       N_Vector **forcing, int *nforcing)
{
  if (stepper == NULL)
  {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepInnerStepper_GetForcingData",
                    "Inner stepper memory is NULL");
    return ARK_ILL_INPUT;
  }

  *tshift   = stepper->tshift;
  *tscale   = stepper->tscale;
  *forcing  = stepper->forcing;
  *nforcing = stepper->nforcing;

  return ARK_SUCCESS;
}


/*---------------------------------------------------------------
  Internal inner integrator functions
  ---------------------------------------------------------------*/


/* Check for required operations */
int mriStepInnerStepper_HasRequiredOps(MRIStepInnerStepper stepper)
{
  if (stepper == NULL) return ARK_ILL_INPUT;
  if (stepper->ops == NULL) return ARK_ILL_INPUT;

  if (stepper->ops->evolve && stepper->ops->fullrhs)
    return ARK_SUCCESS;
  else
    return ARK_ILL_INPUT;
}


/* Evolve the inner (fast) ODE */
int mriStepInnerStepper_Evolve(MRIStepInnerStepper stepper,
                               realtype t0, realtype tout, N_Vector y)
{
  if (stepper == NULL) return ARK_ILL_INPUT;
  if (stepper->ops == NULL) return ARK_ILL_INPUT;
  if (stepper->ops->evolve == NULL) return ARK_ILL_INPUT;

  stepper->last_flag = stepper->ops->evolve(stepper, t0, tout, y);
  return stepper->last_flag;
}


/* Compute the full RHS for inner (fast) time scale TODO(DJG): This function can
   be made optional when fullrhs is not called unconditionally by the ARKODE
   infrastructure e.g., in arkInitialSetup, arkYddNorm, and arkCompleteStep. */
int mriStepInnerStepper_FullRhs(MRIStepInnerStepper stepper,
                                realtype t, N_Vector y, N_Vector f,
                                int mode)
{
  if (stepper == NULL) return ARK_ILL_INPUT;
  if (stepper->ops == NULL) return ARK_ILL_INPUT;
  if (stepper->ops->fullrhs == NULL) return ARK_ILL_INPUT;

  stepper->last_flag = stepper->ops->fullrhs(stepper, t, y, f, mode);
  return stepper->last_flag;
}


/* Reset the inner (fast) stepper state */
int mriStepInnerStepper_Reset(MRIStepInnerStepper stepper,
                              realtype tR, N_Vector yR)
{
  if (stepper == NULL) return ARK_ILL_INPUT;
  if (stepper->ops == NULL) return ARK_ILL_INPUT;

  if (stepper->ops->reset) {
    stepper->last_flag = stepper->ops->reset(stepper, tR, yR);
    return stepper->last_flag;
  } else {
    /* assume stepper uses input state and does not need to be reset */
    return ARK_SUCCESS;
  }
}


/* Allocate MRI forcing and fused op workspace vectors if necessary */
int mriStepInnerStepper_AllocVecs(MRIStepInnerStepper stepper, int count,
                                  N_Vector tmpl)
{
  sunindextype lrw1, liw1;

  if (stepper == NULL) return ARK_ILL_INPUT;

  /* Set space requirements for one N_Vector */
  if (tmpl->ops->nvspace) {
    N_VSpace(tmpl, &lrw1, &liw1);
  } else {
    lrw1 = 0;
    liw1 = 0;
  }
  stepper->lrw1 = lrw1;
  stepper->liw1 = liw1;

  /* Set the number of forcing vectors and allocate vectors */
  stepper->nforcing = count;

  if (!arkAllocVecArray(count,
                        tmpl, &(stepper->forcing),
                        stepper->lrw1, &(stepper->lrw),
                        stepper->liw1, &(stepper->liw))) {
    mriStepInnerStepper_FreeVecs(stepper);
    return(ARK_MEM_FAIL);
  }

  /* Allocate fused operation workspace arrays */
  if (stepper->vecs == NULL) {
    stepper->vecs = (N_Vector *) calloc(count + 1, sizeof(N_Vector));
    if (stepper->vecs == NULL) {
      mriStepInnerStepper_FreeVecs(stepper);
      return(ARK_MEM_FAIL);
    }
  }

  if (stepper->vals == NULL) {
    stepper->vals = (realtype *) calloc(count + 1, sizeof(realtype));
    if (stepper->vals == NULL) {
      mriStepInnerStepper_FreeVecs(stepper);
      return(ARK_MEM_FAIL);
    }
  }

  return(ARK_SUCCESS);
}


/* Resize MRI forcing and fused op workspace vectors if necessary */
int mriStepInnerStepper_Resize(MRIStepInnerStepper stepper,
                               ARKVecResizeFn resize, void* resize_data,
                               sunindextype lrw_diff, sunindextype liw_diff,
                               N_Vector tmpl)
{
  int retval;

  if (stepper == NULL) return ARK_ILL_INPUT;

  retval = arkResizeVecArray(resize, resize_data,
                             stepper->nforcing, tmpl, &(stepper->forcing),
                             lrw_diff, &(stepper->lrw),
                             liw_diff, &(stepper->liw));
  if (retval != ARK_SUCCESS) return(ARK_MEM_FAIL);

  return(ARK_SUCCESS);
}


/* Free MRI forcing and fused op workspace vectors if necessary */
int mriStepInnerStepper_FreeVecs(MRIStepInnerStepper stepper)
{
  if (stepper == NULL) return ARK_ILL_INPUT;

  arkFreeVecArray(stepper->nforcing, &(stepper->forcing),
                  stepper->lrw1, &(stepper->lrw),
                  stepper->liw1, &(stepper->liw));

  if (stepper->vecs != NULL) {
    free(stepper->vecs);
    stepper->vecs = NULL;
  }

  if (stepper->vals != NULL) {
    free(stepper->vals);
    stepper->vals = NULL;
  }

  return(ARK_SUCCESS);
}


/* Print forcing vectors to output file */
void mriStepInnerStepper_PrintMem(MRIStepInnerStepper stepper,
                                  FILE* outfile)
{
#ifdef SUNDIALS_DEBUG_PRINTVEC
  int i;
#endif
  if (stepper == NULL) return;

  /* output data from the inner stepper */
  fprintf(outfile,"MRIStepInnerStepper Mem:\n");
  fprintf(outfile,"MRIStepInnerStepper: inner_nforcing = %i\n",
          stepper->nforcing);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  if (stepper->forcing != NULL) {
    for (i = 0; i < stepper->nforcing; i++) {
      fprintf(outfile,"MRIStep: inner_forcing[%i]:\n", i);
      N_VPrintFile(stepper->forcing[i], outfile);
    }
  }
#endif

  return;
}

/*===============================================================
  EOF
  ===============================================================*/
