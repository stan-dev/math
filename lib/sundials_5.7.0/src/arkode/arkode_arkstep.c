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
 * This is the implementation file for ARKode's ARK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_arkstep_impl.h"
#include "arkode_interp_impl.h"
#include <sundials/sundials_math.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

/* constants */
#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)

#define FIXED_LIN_TOL


/*===============================================================
  ARKStep Exported functions -- Required
  ===============================================================*/

void* ARKStepCreate(ARKRhsFn fe, ARKRhsFn fi, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  SUNNonlinearSolver NLS;
  booleantype nvectorOK;
  int retval;

  /* Check that at least one of fe, fi is supplied and is to be used */
  if (fe == NULL && fi == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_NULL_F);
    return(NULL);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_NULL_Y0);
    return(NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = arkStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_BAD_NVECTOR);
    return(NULL);
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate();
  if (ark_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_NO_MEM);
    return(NULL);
  }

  /* Allocate ARKodeARKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeARKStepMem) malloc(sizeof(struct ARKodeARKStepMemRec));
  if (step_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_ARKMEM_FAIL);
    return(NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeARKStepMemRec));

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_attachlinsol   = arkStep_AttachLinsol;
  ark_mem->step_attachmasssol  = arkStep_AttachMasssol;
  ark_mem->step_disablelsetup  = arkStep_DisableLSetup;
  ark_mem->step_disablemsetup  = arkStep_DisableMSetup;
  ark_mem->step_getlinmem      = arkStep_GetLmem;
  ark_mem->step_getmassmem     = arkStep_GetMassMem;
  ark_mem->step_getimplicitrhs = arkStep_GetImplicitRHS;
  ark_mem->step_mmult          = NULL;
  ark_mem->step_getgammas      = arkStep_GetGammas;
  ark_mem->step_init           = arkStep_Init;
  ark_mem->step_fullrhs        = arkStep_FullRHS;
  ark_mem->step                = arkStep_TakeStep_Z;
  ark_mem->step_mem            = (void*) step_mem;

  /* Set default values for ARKStep optional inputs */
  retval = ARKStepSetDefaults((void *)ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ARKStep",
                    "ARKStepCreate",
                    "Error setting default solver options");
    ARKStepFree((void**) &ark_mem);  return(NULL);
  }

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Allocate the general ARK stepper vectors using y0 as a template */
  /* NOTE: Fe, Fi, cvals and Xvecs will be allocated later on
     (based on the number of ARK stages) */

  /* Clone the input vector to create sdata, zpred and zcor */
  if (!arkAllocVec(ark_mem, y0, &(step_mem->sdata))) {
    ARKStepFree((void**) &ark_mem);  return(NULL); }
  if (!arkAllocVec(ark_mem, y0, &(step_mem->zpred))) {
    ARKStepFree((void**) &ark_mem);  return(NULL); }
  if (!arkAllocVec(ark_mem, y0, &(step_mem->zcor))) {
    ARKStepFree((void**) &ark_mem);  return(NULL); }

  /* Copy the input parameters into ARKode state */
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Update the ARKode workspace requirements */
  ark_mem->liw += 41;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->lrw += 10;

  /* If an implicit component is to be solved, create default Newton NLS object */
  step_mem->ownNLS = SUNFALSE;
  if (step_mem->implicit)  {
    NLS = NULL;
    NLS = SUNNonlinSol_Newton(y0);
    if (NLS == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                      "ARKStepCreate", "Error creating default Newton solver");
      ARKStepFree((void**) &ark_mem);  return(NULL);
    }
    retval = ARKStepSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                      "ARKStepCreate", "Error attaching default Newton solver");
      ARKStepFree((void**) &ark_mem);  return(NULL);
    }
    step_mem->ownNLS = SUNTRUE;
  }

  /* Set the linear solver addresses to NULL (we check != NULL later) */
  step_mem->linit       = NULL;
  step_mem->lsetup      = NULL;
  step_mem->lsolve      = NULL;
  step_mem->lfree       = NULL;
  step_mem->lmem        = NULL;
  step_mem->lsolve_type = -1;

  /* Set the mass matrix solver addresses to NULL */
  step_mem->minit       = NULL;
  step_mem->msetup      = NULL;
  step_mem->mmult       = NULL;
  step_mem->msolve      = NULL;
  step_mem->mfree       = NULL;
  step_mem->mass_mem    = NULL;
  step_mem->mass_type   = MASS_IDENTITY;
  step_mem->msolve_type = -1;

  /* Initialize initial error norm  */
  step_mem->eRNrm = ONE;

  /* Initialize all the counters */
  step_mem->nfe       = 0;
  step_mem->nfi       = 0;
  step_mem->nsetups   = 0;
  step_mem->nstlp     = 0;
  step_mem->nls_iters = 0;

  /* Initialize fused op work space */
  step_mem->cvals        = NULL;
  step_mem->Xvecs        = NULL;
  step_mem->nfusedopvecs = 0;

  /* Initialize external polynomial forcing data */
  step_mem->expforcing = SUNFALSE;
  step_mem->impforcing = SUNFALSE;
  step_mem->forcing    = NULL;
  step_mem->nforcing   = 0;

  /* Initialize main ARKode infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ARKStep", "ARKStepCreate",
                    "Unable to initialize main ARKode infrastructure");
    ARKStepFree((void**) &ark_mem);  return(NULL);
  }

  return((void *)ark_mem);
}


/*---------------------------------------------------------------
  ARKStepResize:

  This routine resizes the memory within the ARKStep module.
  It first resizes the main ARKode infrastructure memory, and
  then resizes its own data.
  ---------------------------------------------------------------*/
int ARKStepResize(void *arkode_mem, N_Vector y0, realtype hscale,
                  realtype t0, ARKVecResizeFn resize, void *resize_data)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  SUNNonlinearSolver NLS;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int i, retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepResize",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Determing change in vector sizes */
  lrw1 = liw1 = 0;
  if (y0->ops->nvspace != NULL)
    N_VSpace(y0, &lrw1, &liw1);
  lrw_diff = lrw1 - ark_mem->lrw1;
  liw_diff = liw1 - ark_mem->liw1;
  ark_mem->lrw1 = lrw1;
  ark_mem->liw1 = liw1;

  /* resize ARKode infrastructure memory */
  retval = arkResize(ark_mem, y0, hscale, t0, resize, resize_data);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ARKStep", "ARKStepResize",
                    "Unable to resize main ARKode infrastructure");
    return(retval);
  }

  /* Resize the sdata, zpred and zcor vectors */
  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, y0, &step_mem->sdata)) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "ARKStepResize",
                    "Unable to resize vector");
    return(ARK_MEM_FAIL);
  }

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, y0, &step_mem->zpred)) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "ARKStepResize",
                    "Unable to resize vector");
    return(ARK_MEM_FAIL);
  }

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, y0, &step_mem->zcor)) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "ARKStepResize",
                    "Unable to resize vector");
    return(ARK_MEM_FAIL);
  }

  /* Resize the ARKStep vectors */
  /*     Fe */
  if (step_mem->Fe != NULL) {
    for (i=0; i<step_mem->stages; i++) {
      if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                        liw_diff, y0, &step_mem->Fe[i])) {
        arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "ARKStepResize",
                        "Unable to resize vector");
        return(ARK_MEM_FAIL);
      }
    }
  }
  /*     Fi */
  if (step_mem->Fi != NULL) {
    for (i=0; i<step_mem->stages; i++) {
      if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                        liw_diff, y0, &step_mem->Fi[i])) {
        arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "ARKStepResize",
                        "Unable to resize vector");
        return(ARK_MEM_FAIL);
      }
    }
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
    NLS = NULL;
    NLS = SUNNonlinSol_Newton(y0);
    if (NLS == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                      "ARKStepResize", "Error creating default Newton solver");
      return(ARK_MEM_FAIL);
    }

    /* attach new Newton NLS object to ARKStep */
    retval = ARKStepSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                      "ARKStepResize", "Error attaching default Newton solver");
      return(ARK_MEM_FAIL);
    }
    step_mem->ownNLS = SUNTRUE;

  }

  /* reset nonlinear solver counters */
  if (step_mem->NLS != NULL)  step_mem->nsetups = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepReInit:

  This routine re-initializes the ARKStep module to solve a new
  problem of the same size as was previously solved. This routine
  should also be called when the problem dynamics or desired solvers
  have changed dramatically, so that the problem integration should
  resume as if started from scratch.

  Note all internal counters are set to 0 on re-initialization.
  ---------------------------------------------------------------*/
int ARKStepReInit(void* arkode_mem, ARKRhsFn fe,
                  ARKRhsFn fi, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepReInit",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Check if ark_mem was allocated */
  if (ark_mem->MallocDone == SUNFALSE) {
    arkProcessError(ark_mem, ARK_NO_MALLOC, "ARKode::ARKStep",
                    "ARKStepReInit", MSG_ARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check that at least one of fe, fi is supplied and is to be used */
  if (fe == NULL && fi == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepReInit", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Check that y0 is supplied */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Copy the input parameters into ARKode state */
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Initialize initial error norm  */
  step_mem->eRNrm = ONE;

  /* Initialize main ARKode infrastructure */
  retval = arkInit(ark_mem, t0, y0, FIRST_INIT);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ARKStep", "ARKStepReInit",
                    "Unable to reinitialize main ARKode infrastructure");
    return(retval);
  }

  /* Initialize all the counters */
  step_mem->nfe     = 0;
  step_mem->nfi     = 0;
  step_mem->nsetups = 0;
  step_mem->nstlp   = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepReset:

  This routine resets the ARKStep module state to solve the same
  problem from the given time with the input state (all counter
  values are retained).
  ---------------------------------------------------------------*/
int ARKStepReset(void* arkode_mem, realtype tR, N_Vector yR)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepReset",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Initialize main ARKode infrastructure */
  retval = arkInit(ark_mem, tR, yR, RESET_INIT);

  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ARKStep", "ARKStepReset",
                    "Unable to initialize main ARKode infrastructure");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepSStolerances, ARKStepSVtolerances, ARKStepWFtolerances,
  ARKStepResStolerance, ARKStepResVtolerance, ARKStepResFtolerance:

  These routines set integration tolerances (wrappers for general
  ARKode utility routines)
  ---------------------------------------------------------------*/
int ARKStepSStolerances(void *arkode_mem, realtype reltol, realtype abstol)
{
  /* unpack ark_mem, call arkSStolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSStolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSStolerances(ark_mem, reltol, abstol));
}

int ARKStepSVtolerances(void *arkode_mem, realtype reltol, N_Vector abstol)
{
  /* unpack ark_mem, call arkSVtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepSVtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSVtolerances(ark_mem, reltol, abstol));
}

int ARKStepWFtolerances(void *arkode_mem, ARKEwtFn efun)
{
  /* unpack ark_mem, call arkWFtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepWFtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkWFtolerances(ark_mem, efun));
}

int ARKStepResStolerance(void *arkode_mem, realtype rabstol)
{
  /* unpack ark_mem, call arkResStolerance, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepResStolerance", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkResStolerance(ark_mem, rabstol));
}

int ARKStepResVtolerance(void *arkode_mem, N_Vector rabstol)
{
  /* unpack ark_mem, call arkResVtolerance, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepResVtolerance", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkResVtolerance(ark_mem, rabstol));
}

int ARKStepResFtolerance(void *arkode_mem, ARKRwtFn rfun)
{
  /* unpack ark_mem, call arkResFtolerance, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepResFtolerance", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkResFtolerance(ark_mem, rfun));
}


/*---------------------------------------------------------------
  ARKStepRootInit:

  Initialize (attach) a rootfinding problem to the stepper
  (wrappers for general ARKode utility routine)
  ---------------------------------------------------------------*/
int ARKStepRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)
{
  /* unpack ark_mem, call arkRootInit, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepRootInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkRootInit(ark_mem, nrtfn, g));
}


/*---------------------------------------------------------------
  ARKStepEvolve:

  This is the main time-integration driver (wrappers for general
  ARKode utility routine)
  ---------------------------------------------------------------*/
int ARKStepEvolve(void *arkode_mem, realtype tout, N_Vector yout,
                  realtype *tret, int itask)
{
  /* unpack ark_mem, call arkEvolve, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepEvolve", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkEvolve(ark_mem, tout, yout, tret, itask));
}


/*---------------------------------------------------------------
  ARKStepGetDky:

  This returns interpolated output of the solution or its
  derivatives over the most-recently-computed step (wrapper for
  generic ARKode utility routine)
  ---------------------------------------------------------------*/
int ARKStepGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  /* unpack ark_mem, call arkGetDky, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "ARKStepGetDky", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetDky(ark_mem, t, k, dky));
}


/*---------------------------------------------------------------
  ARKStepComputeState:

  Computes y based on the current prediction and given correction.
  ---------------------------------------------------------------*/
int ARKStepComputeState(void *arkode_mem, N_Vector zcor, N_Vector z)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepComputeState",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, z);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepFree frees all ARKStep memory, and then calls an ARKode
  utility routine to free the ARKode infrastructure memory.
  ---------------------------------------------------------------*/
void ARKStepFree(void **arkode_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL)  return;

  /* conditional frees on non-NULL ARKStep module */
  ark_mem = (ARKodeMem) (*arkode_mem);
  if (ark_mem->step_mem != NULL) {

    step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

    /* free the Butcher tables */
    if (step_mem->Be != NULL) {
      ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
      ARKodeButcherTable_Free(step_mem->Be);
      step_mem->Be = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
    }
    if (step_mem->Bi != NULL) {
      ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
      ARKodeButcherTable_Free(step_mem->Bi);
      step_mem->Bi = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
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

    /* free the mass matrix solver memory */
    if (step_mem->mfree != NULL) {
      step_mem->mfree((void *) ark_mem);
      step_mem->mass_mem = NULL;
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
    if (step_mem->Fe != NULL) {
      for(j=0; j<step_mem->stages; j++)
        arkFreeVec(ark_mem, &step_mem->Fe[j]);
      free(step_mem->Fe);
      step_mem->Fe = NULL;
      ark_mem->liw -= step_mem->stages;
    }
    if (step_mem->Fi != NULL) {
      for(j=0; j<step_mem->stages; j++)
        arkFreeVec(ark_mem, &step_mem->Fi[j]);
      free(step_mem->Fi);
      step_mem->Fi = NULL;
      ark_mem->liw -= step_mem->stages;
    }

    /* free the reusable arrays for fused vector interface */
    if (step_mem->cvals != NULL) {
      free(step_mem->cvals);
      step_mem->cvals = NULL;
      ark_mem->lrw -= step_mem->nfusedopvecs;
    }
    if (step_mem->Xvecs != NULL) {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= step_mem->nfusedopvecs;
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
  ARKStepPrintMem:

  This routine outputs the memory from the ARKStep structure and
  the main ARKode infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
void ARKStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

#ifdef SUNDIALS_DEBUG_PRINTVEC
  int i;
#endif

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepPrintMem",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return;

  /* if outfile==NULL, set it to stdout */
  if (outfile == NULL)  outfile = stdout;

  /* output data from main ARKode infrastructure */
  arkPrintMem(ark_mem, outfile);

  /* output integer quantities */
  fprintf(outfile,"ARKStep: q = %i\n", step_mem->q);
  fprintf(outfile,"ARKStep: p = %i\n", step_mem->p);
  fprintf(outfile,"ARKStep: istage = %i\n", step_mem->istage);
  fprintf(outfile,"ARKStep: stages = %i\n", step_mem->stages);
  fprintf(outfile,"ARKStep: maxcor = %i\n", step_mem->maxcor);
  fprintf(outfile,"ARKStep: msbp = %i\n", step_mem->msbp);
  fprintf(outfile,"ARKStep: predictor = %i\n", step_mem->predictor);
  fprintf(outfile,"ARKStep: lsolve_type = %i\n", step_mem->lsolve_type);
  fprintf(outfile,"ARKStep: msolve_type = %i\n", step_mem->msolve_type);
  fprintf(outfile,"ARKStep: convfail = %i\n", step_mem->convfail);

  /* output long integer quantities */
  fprintf(outfile,"ARKStep: nfe = %li\n", step_mem->nfe);
  fprintf(outfile,"ARKStep: nfi = %li\n", step_mem->nfi);
  fprintf(outfile,"ARKStep: nsetups = %li\n", step_mem->nsetups);
  fprintf(outfile,"ARKStep: nstlp = %li\n", step_mem->nstlp);

  /* output boolean quantities */
  fprintf(outfile,"ARKStep: user_linear = %i\n", step_mem->linear);
  fprintf(outfile,"ARKStep: user_linear_timedep = %i\n", step_mem->linear_timedep);
  fprintf(outfile,"ARKStep: user_explicit = %i\n", step_mem->explicit);
  fprintf(outfile,"ARKStep: user_implicit = %i\n", step_mem->implicit);
  fprintf(outfile,"ARKStep: jcur = %i\n", step_mem->jcur);

  /* output realtype quantities */
  if (step_mem->Be != NULL) {
    fprintf(outfile,"ARKStep: explicit Butcher table:\n");
    ARKodeButcherTable_Write(step_mem->Be, outfile);
  }
  if (step_mem->Bi != NULL) {
    fprintf(outfile,"ARKStep: implicit Butcher table:\n");
    ARKodeButcherTable_Write(step_mem->Bi, outfile);
  }
  fprintf(outfile,"ARKStep: gamma = %"RSYM"\n", step_mem->gamma);
  fprintf(outfile,"ARKStep: gammap = %"RSYM"\n", step_mem->gammap);
  fprintf(outfile,"ARKStep: gamrat = %"RSYM"\n", step_mem->gamrat);
  fprintf(outfile,"ARKStep: crate = %"RSYM"\n", step_mem->crate);
  fprintf(outfile,"ARKStep: eRNrm = %"RSYM"\n", step_mem->eRNrm);
  fprintf(outfile,"ARKStep: nlscoef = %"RSYM"\n", step_mem->nlscoef);
  fprintf(outfile,"ARKStep: crdown = %"RSYM"\n", step_mem->crdown);
  fprintf(outfile,"ARKStep: rdiv = %"RSYM"\n", step_mem->rdiv);
  fprintf(outfile,"ARKStep: dgmax = %"RSYM"\n", step_mem->dgmax);

#ifdef SUNDIALS_DEBUG_PRINTVEC
  /* output vector quantities */
  fprintf(outfile, "ARKStep: sdata:\n");
  N_VPrintFile(step_mem->sdata, outfile);
  fprintf(outfile, "ARKStep: zpred:\n");
  N_VPrintFile(step_mem->zpred, outfile);
  fprintf(outfile, "ARKStep: zcor:\n");
  N_VPrintFile(step_mem->zcor, outfile);
  if (step_mem->Fe != NULL)
    for (i=0; i<step_mem->stages; i++) {
      fprintf(outfile,"ARKStep: Fe[%i]:\n", i);
      N_VPrintFile(step_mem->Fe[i], outfile);
    }
  if (step_mem->Fi != NULL)
    for (i=0; i<step_mem->stages; i++) {
      fprintf(outfile,"ARKStep: Fi[%i]:\n", i);
      N_VPrintFile(step_mem->Fi[i], outfile);
    }
#endif
}


/*===============================================================
  ARKStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Interface routines supplied to ARKode
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  arkStep_AttachLinsol:

  This routine attaches the various set of system linear solver
  interface routines, data structure, and solver type to the
  ARKStep module.
  ---------------------------------------------------------------*/
int arkStep_AttachLinsol(void* arkode_mem, ARKLinsolInitFn linit,
                         ARKLinsolSetupFn lsetup,
                         ARKLinsolSolveFn lsolve,
                         ARKLinsolFreeFn lfree,
                         SUNLinearSolver_Type lsolve_type,
                         void *lmem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_AttachLinsol",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* free any existing system solver */
  if (step_mem->lfree != NULL)  step_mem->lfree(arkode_mem);

  /* Attach the provided routines, data structure and solve type */
  step_mem->linit       = linit;
  step_mem->lsetup      = lsetup;
  step_mem->lsolve      = lsolve;
  step_mem->lfree       = lfree;
  step_mem->lmem        = lmem;
  step_mem->lsolve_type = lsolve_type;

  /* Reset all linear solver counters */
  step_mem->nsetups = 0;
  step_mem->nstlp   = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_AttachMasssol:

  This routine attaches the set of mass matrix linear solver
  interface routines, data structure, and solver type to the
  ARKStep module.
  ---------------------------------------------------------------*/
int arkStep_AttachMasssol(void* arkode_mem,
                          ARKMassInitFn minit,
                          ARKMassSetupFn msetup,
                          ARKMassMultFn mmult,
                          ARKMassSolveFn msolve,
                          ARKMassFreeFn mfree,
                          booleantype time_dep,
                          SUNLinearSolver_Type msolve_type,
                          void *mass_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_AttachMasssol",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* free any existing mass matrix solver */
  if (step_mem->mfree != NULL)  step_mem->mfree(arkode_mem);

  /* Attach the provided routines, data structure and solve type */
  step_mem->minit       = minit;
  step_mem->msetup      = msetup;
  step_mem->mmult       = mmult;
  step_mem->msolve      = msolve;
  step_mem->mfree       = mfree;
  step_mem->mass_mem    = mass_mem;
  step_mem->mass_type   = (time_dep) ? MASS_TIMEDEP : MASS_FIXED;
  step_mem->msolve_type = msolve_type;

  /* Attach mmult function pointer to ark_mem as well */
  ark_mem->step_mmult = mmult;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_DisableLSetup:

  This routine NULLifies the lsetup function pointer in the
  ARKStep module.
  ---------------------------------------------------------------*/
void arkStep_DisableLSetup(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL)  return;
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL)  return;
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* nullify the lsetup function pointer */
  step_mem->lsetup = NULL;
}


/*---------------------------------------------------------------
  arkStep_DisableMSetup:

  This routine NULLifies the msetup function pointer in the
  ARKStep module.
  ---------------------------------------------------------------*/
void arkStep_DisableMSetup(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (arkode_mem==NULL)  return;
  ark_mem = (ARKodeMem) arkode_mem;
  if (ark_mem->step_mem==NULL)  return;
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* nullify the msetup function pointer */
  step_mem->msetup = NULL;
}


/*---------------------------------------------------------------
  arkStep_GetLmem:

  This routine returns the system linear solver interface memory
  structure, lmem.
  ---------------------------------------------------------------*/
void* arkStep_GetLmem(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure, and return lmem */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_GetLmem",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(NULL);
  return(step_mem->lmem);
}


/*---------------------------------------------------------------
  arkStep_GetMassMem:

  This routine returns the mass matrix solver interface memory
  structure, mass_mem.
  ---------------------------------------------------------------*/
void* arkStep_GetMassMem(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure, and return mass_mem */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_GetMassMem",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(NULL);
  return(step_mem->mass_mem);
}


/*---------------------------------------------------------------
  arkStep_GetImplicitRHS:

  This routine returns the implicit RHS function pointer, fi.
  ---------------------------------------------------------------*/
ARKRhsFn arkStep_GetImplicitRHS(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure, and return fi */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_GetImplicitRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(NULL);
  return(step_mem->fi);
}


/*---------------------------------------------------------------
  arkStep_GetGammas:

  This routine fills the current value of gamma, and states
  whether the gamma ratio fails the dgmax criteria.
  ---------------------------------------------------------------*/
int arkStep_GetGammas(void* arkode_mem, realtype *gamma,
                      realtype *gamrat, booleantype **jcur,
                      booleantype *dgamma_fail)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_GetGammas",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* set outputs */
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;
  *gamma  = step_mem->gamma;
  *gamrat = step_mem->gamrat;
  *jcur = &step_mem->jcur;
  *dgamma_fail = (SUNRabs(*gamrat - ONE) >= step_mem->dgmax);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_Init:

  This routine is called just prior to performing internal time
  steps (after all user "set" routines have been called) from
  within arkInitialSetup.

  For all initialization types, this routine sets the relevant
  TakeStep routine based on the current problem configuration.

  With initialization type FIRST_INIT this routine:
  - sets/checks the ARK Butcher tables to be used
  - allocates any memory that depends on the number of ARK stages,
    method order, or solver options
  - checks for consistency between the system and mass matrix
    linear solvers (if applicable)
  - initializes and sets up the system and mass matrix linear
    solvers (if applicable)
  - initializes and sets up the nonlinear solver (if applicable)
  - allocates the interpolation data structure (if needed based
    on ARKStep solver options)
  - updates the call_fullrhs flag if necessary

  With initialization type FIRST_INIT or RESIZE_INIT, this routine:
  - sets the relevant TakeStep routine based on the current
    problem configuration
  - checks for consistency between the system and mass matrix
    linear solvers (if applicable)
  - initializes and sets up the system and mass matrix linear
    solvers (if applicable)
  - initializes and sets up the nonlinear solver (if applicable)

  With initialization type RESET_INIT, this routine does nothing.
  ---------------------------------------------------------------*/
int arkStep_Init(void* arkode_mem, int init_type)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int j, retval;
  booleantype reset_efun;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_Init",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* immediately return if reset */
  if (init_type == RESET_INIT) return(ARK_SUCCESS);

  /* initializations/checks for (re-)initialization call */
  if (init_type == FIRST_INIT) {

    /* enforce use of arkEwtSmallReal if using a fixed step size for
       an explicit method, an internal error weight function, and not
       using an iterative mass matrix solver with rwt=ewt */
    reset_efun = SUNTRUE;
    if ( step_mem->implicit )   reset_efun = SUNFALSE;
    if ( !ark_mem->fixedstep )  reset_efun = SUNFALSE;
    if ( ark_mem->user_efun )   reset_efun = SUNFALSE;
    if ( ark_mem->rwt_is_ewt && (step_mem->msolve_type == SUNLINEARSOLVER_ITERATIVE) )
      reset_efun = SUNFALSE;
    if ( ark_mem->rwt_is_ewt && (step_mem->msolve_type == SUNLINEARSOLVER_MATRIX_ITERATIVE) )
      reset_efun = SUNFALSE;
    if (reset_efun) {
      ark_mem->user_efun = SUNFALSE;
      ark_mem->efun      = arkEwtSetSmallReal;
      ark_mem->e_data    = ark_mem;
    }

    /* Create Butcher tables (if not already set) */
    retval = arkStep_SetButcherTables(ark_mem);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep", "arkStep_Init",
                      "Could not create Butcher table(s)");
      return(ARK_ILL_INPUT);
    }

    /* Check that Butcher tables are OK */
    retval = arkStep_CheckButcherTables(ark_mem);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_Init", "Error in Butcher table(s)");
      return(ARK_ILL_INPUT);
    }

    /* Retrieve/store method and embedding orders now that tables are finalized */
    if (step_mem->Bi != NULL) {
      step_mem->q = ark_mem->hadapt_mem->q = step_mem->Bi->q;
      step_mem->p = ark_mem->hadapt_mem->p = step_mem->Bi->p;
    } else {
      step_mem->q = ark_mem->hadapt_mem->q = step_mem->Be->q;
      step_mem->p = ark_mem->hadapt_mem->p = step_mem->Be->p;
    }

    /* Ensure that if adaptivity is enabled, then method includes embedding coefficients */
    if (!ark_mem->fixedstep && (step_mem->p == 0)) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep", "arkStep_Init",
                      "Adaptive timestepping cannot be performed without embedding coefficients");
      return(ARK_ILL_INPUT);
    }

    /* Allocate ARK RHS vector memory, update storage requirements */
    /*   Allocate Fe[0] ... Fe[stages-1] if needed */
    if (step_mem->explicit) {
      if (step_mem->Fe == NULL)
        step_mem->Fe = (N_Vector *) calloc(step_mem->stages, sizeof(N_Vector));
      for (j=0; j<step_mem->stages; j++) {
        if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->Fe[j])))
          return(ARK_MEM_FAIL);
      }
      ark_mem->liw += step_mem->stages;  /* pointers */
    }

    /*   Allocate Fi[0] ... Fi[stages-1] if needed */
    if (step_mem->implicit) {
      if (step_mem->Fi == NULL)
        step_mem->Fi = (N_Vector *) calloc(step_mem->stages, sizeof(N_Vector));
      for (j=0; j<step_mem->stages; j++) {
        if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->Fi[j])))
          return(ARK_MEM_FAIL);
      }
      ark_mem->liw += step_mem->stages;  /* pointers */
    }

    /* Allocate reusable arrays for fused vector operations */
    step_mem->nfusedopvecs = 2 * step_mem->stages + 2 + step_mem->nforcing;
    if (step_mem->cvals == NULL) {
      step_mem->cvals = (realtype *) calloc(step_mem->nfusedopvecs,
                                            sizeof(realtype));
      if (step_mem->cvals == NULL)  return(ARK_MEM_FAIL);
      ark_mem->lrw += step_mem->nfusedopvecs;
    }
    if (step_mem->Xvecs == NULL) {
      step_mem->Xvecs = (N_Vector *) calloc(step_mem->nfusedopvecs,
                                            sizeof(N_Vector));
      if (step_mem->Xvecs == NULL)  return(ARK_MEM_FAIL);
      ark_mem->liw += step_mem->nfusedopvecs;   /* pointers */
    }

    /* Limit interpolant degree based on method order (use negative
       argument to specify update instead of overwrite) */
    if (ark_mem->interp != NULL) {
      retval = arkInterpSetDegree(ark_mem, ark_mem->interp, -(step_mem->q-1));
      if (retval != ARK_SUCCESS) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep", "arkStep_Init",
                        "Unable to update interpolation polynomial degree");
        return(ARK_ILL_INPUT);
      }
    }

    /* If configured with either predictor 4 or 5 and a non-identity mass
       matrix, reset to trivial predictor */
    if (step_mem->mass_type != MASS_IDENTITY)
      if ((step_mem->predictor == 4) || (step_mem->predictor == 5))
        step_mem->predictor = 0;

    /* If the bootstrap predictor is enabled, signal to shared arkode module that
       fullrhs is required after each step */
    if (step_mem->predictor == 4)  ark_mem->call_fullrhs = SUNTRUE;
  }

  /* set appropriate TakeStep routine based on problem configuration */
  /*    (only one choice for now) */
  ark_mem->step = arkStep_TakeStep_Z;

  /* Check for consistency between mass system and system linear system modules
     (e.g., if lsolve is direct, msolve needs to match) */
  if ((step_mem->mass_type != MASS_IDENTITY) && step_mem->lmem) {
    if (step_mem->lsolve_type != step_mem->msolve_type) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep", "arkStep_Init",
                      "Incompatible linear and mass matrix solvers");
      return(ARK_ILL_INPUT);
    }
  }

  /* Perform mass matrix solver initialization and setup (if applicable) */
  if (step_mem->mass_type != MASS_IDENTITY) {

    /* Call minit (if it exists) */
    if (step_mem->minit != NULL) {
      retval = step_mem->minit((void *) ark_mem);
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_MASSINIT_FAIL, "ARKode::ARKStep",
                        "arkStep_Init", MSG_ARK_MASSINIT_FAIL);
        return(ARK_MASSINIT_FAIL);
      }
    }

    /* Call msetup (if it exists) */
    if (step_mem->msetup != NULL) {
      retval = step_mem->msetup((void *) ark_mem, ark_mem->tcur,
                                ark_mem->tempv1, ark_mem->tempv2,
                                ark_mem->tempv3);
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_MASSSETUP_FAIL, "ARKode::ARKStep",
                        "arkStep_Init", MSG_ARK_MASSSETUP_FAIL);
        return(ARK_MASSSETUP_FAIL);
      }
    }
  }

  /* Call linit (if it exists) */
  if (step_mem->linit) {
    retval = step_mem->linit(ark_mem);
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_LINIT_FAIL, "ARKode::ARKStep",
                      "arkStep_Init", MSG_ARK_LINIT_FAIL);
      return(ARK_LINIT_FAIL);
    }
  }

  /* Initialize the nonlinear solver object (if it exists) */
  if (step_mem->NLS) {
    retval = arkStep_NlsInit(ark_mem);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_NLS_INIT_FAIL, "ARKode::ARKStep", "arkStep_Init",
                      "Unable to initialize SUNNonlinearSolver object");
      return(ARK_NLS_INIT_FAIL);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_FullRHS:

  Rewriting the problem
    My' = fe(t,y) + fi(t,y)
  in the form
    y' = M^{-1}*[ fe(t,y) + fi(t,y) ],
  this routine computes the full right-hand side vector,
    f = M^{-1}*[ fe(t,y) + fi(t,y) ]

  This will be called in one of three 'modes':
    0 -> called at the beginning of a simulation
    1 -> called at the end of a successful step
    2 -> called elsewhere (e.g. for dense output)

  If it is called in mode 0, we store the vectors fe(t,y) and
  fi(t,y) in Fe[0] and Fi[0] for possible reuse in the first
  stage of the subsequent time step.

  If it is called in mode 1 and the ARK method coefficients
  support it, we may just copy vectors Fe[stages] and Fi[stages]
  to fill f instead of calling fe() and fi().

  Mode 2 is only called for dense output in-between steps, or
  when estimating the initial time step size, so we strive to
  store the intermediate parts so that they do not interfere
  with the other two modes.
  ---------------------------------------------------------------*/
int arkStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int nvec, retval;
  booleantype recomputeRHS;
  realtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_FullRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* local shortcuts for use with fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* setup mass-matrix if required (use output f as a temporary) */
  if ((step_mem->mass_type == MASS_TIMEDEP) && (step_mem->msetup != NULL)) {
    retval = step_mem->msetup((void *) ark_mem, t, f,
                              ark_mem->tempv2, ark_mem->tempv3);
    if (retval != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
  }

  /* perform RHS functions contingent on 'mode' argument */
  switch(mode) {

  /* Mode 0: called at the beginning of a simulation
     Store the vectors fe(t,y) and fi(t,y) in Fe[0] and Fi[0] for
     possible reuse in the first stage of the subsequent time step */
  case 0:

    /* call fe if the problem has an explicit component */
    if (step_mem->explicit) {
      retval = step_mem->fe(t, y, step_mem->Fe[0], ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
      /* apply external polynomial forcing */
      if (step_mem->expforcing) {
        cvals[0] = ONE;
        Xvecs[0] = step_mem->Fe[0];
        nvec     = 1;
        arkStep_ApplyForcing(step_mem, t, ONE, &nvec);
        N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fe[0]);
      }
    }

    /* call fi if the problem has an implicit component */
    if (step_mem->implicit) {
      retval = step_mem->fi(t, y, step_mem->Fi[0], ark_mem->user_data);
      step_mem->nfi++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
      /* apply external polynomial forcing */
      if (step_mem->impforcing) {
        cvals[0] = ONE;
        Xvecs[0] = step_mem->Fi[0];
        nvec     = 1;
        arkStep_ApplyForcing(step_mem, t, ONE, &nvec);
        N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fi[0]);
      }
    }

    /* combine RHS vector(s) into output */
    if (step_mem->explicit && step_mem->implicit) { /* ImEx */
      N_VLinearSum(ONE, step_mem->Fi[0], ONE, step_mem->Fe[0], f);
    } else if (step_mem->implicit) {                   /* implicit */
      N_VScale(ONE, step_mem->Fi[0], f);
    } else {                                           /* explicit */
      N_VScale(ONE, step_mem->Fe[0], f);
    }

    break;


  /* Mode 1: called at the end of a successful step
     If the ARK method coefficients support it, we just copy the last stage RHS vectors
     to fill f instead of calling fe() and fi().
     Copy the results to Fe[0] and Fi[0] if the ARK coefficients support it. */
  case 1:

    /* determine if explicit/implicit RHS functions need to be recomputed */
    recomputeRHS = SUNFALSE;
    if ( step_mem->explicit && (SUNRabs(step_mem->Be->c[step_mem->stages-1]-ONE)>TINY) )
      recomputeRHS = SUNTRUE;
    if ( step_mem->implicit && (SUNRabs(step_mem->Bi->c[step_mem->stages-1]-ONE)>TINY) )
      recomputeRHS = SUNTRUE;

    /* base RHS calls on recomputeRHS argument */
    if (recomputeRHS) {

      /* call fe if the problem has an explicit component */
      if (step_mem->explicit) {
        retval = step_mem->fe(t, y, step_mem->Fe[0], ark_mem->user_data);
        step_mem->nfe++;
        if (retval != 0) {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                          "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
          return(ARK_RHSFUNC_FAIL);
        }
        /* apply external polynomial forcing */
        if (step_mem->expforcing) {
          cvals[0] = ONE;
          Xvecs[0] = step_mem->Fe[0];
          nvec     = 1;
          arkStep_ApplyForcing(step_mem, t, ONE, &nvec);
          N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fe[0]);
        }
      }

      /* call fi if the problem has an implicit component */
      if (step_mem->implicit) {
        retval = step_mem->fi(t, y, step_mem->Fi[0], ark_mem->user_data);
        step_mem->nfi++;
        if (retval != 0) {
          arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                          "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
          return(ARK_RHSFUNC_FAIL);
        }
        /* apply external polynomial forcing */
        if (step_mem->impforcing) {
          cvals[0] = ONE;
          Xvecs[0] = step_mem->Fi[0];
          nvec     = 1;
          arkStep_ApplyForcing(step_mem, t, ONE, &nvec);
          N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fi[0]);
        }
      }
    } else {
      if (step_mem->explicit)
        N_VScale(ONE, step_mem->Fe[step_mem->stages-1], step_mem->Fe[0]);
      if (step_mem->implicit)
        N_VScale(ONE, step_mem->Fi[step_mem->stages-1], step_mem->Fi[0]);
    }

    /* combine RHS vector(s) into output */
    if (step_mem->explicit && step_mem->implicit) { /* ImEx */
      N_VLinearSum(ONE, step_mem->Fi[0], ONE, step_mem->Fe[0], f);
    } else if (step_mem->implicit) {                   /* implicit */
      N_VScale(ONE, step_mem->Fi[0], f);
    } else {                                           /* explicit */
      N_VScale(ONE, step_mem->Fe[0], f);
    }

    break;

  /*  Mode 2: called for dense output in-between steps or for estimation
      of the initial time step size, store the intermediate calculations
      in such a way as to not interfere with the other two modes */
  default:

    /* call fe if the problem has an explicit component (store in ark_tempv2) */
    if (step_mem->explicit) {
      retval = step_mem->fe(t, y, ark_mem->tempv2, ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
      /* apply external polynomial forcing */
      if (step_mem->expforcing) {
        cvals[0] = ONE;
        Xvecs[0] = ark_mem->tempv2;
        nvec     = 1;
        arkStep_ApplyForcing(step_mem, t, ONE, &nvec);
        N_VLinearCombination(nvec, cvals, Xvecs, ark_mem->tempv2);
      }
    }

    /* call fi if the problem has an implicit component (store in sdata) */
    if (step_mem->implicit) {
      retval = step_mem->fi(t, y, step_mem->sdata, ark_mem->user_data);
      step_mem->nfi++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ARKStep",
                        "arkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }
      /* apply external polynomial forcing */
      if (step_mem->impforcing) {
        cvals[0] = ONE;
        Xvecs[0] = step_mem->sdata;
        nvec     = 1;
        arkStep_ApplyForcing(step_mem, t, ONE, &nvec);
        N_VLinearCombination(nvec, cvals, Xvecs, step_mem->sdata);
      }
    }

    /* combine RHS vector(s) into output */
    if (step_mem->explicit && step_mem->implicit) { /* ImEx */
      N_VLinearSum(ONE, step_mem->sdata, ONE, ark_mem->tempv2, f);
    } else if (step_mem->implicit) {                   /* implicit */
      N_VScale(ONE, step_mem->sdata, f);
    } else {                                           /* explicit */
      N_VScale(ONE, ark_mem->tempv2, f);
    }

    break;
  }

  /* if M != I, then update f = M^{-1}*f */
  if (step_mem->mass_type != MASS_IDENTITY) {
    retval = step_mem->msolve((void *) ark_mem, f,
                              step_mem->nlscoef/ark_mem->h);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKode::ARKStep",
                      "arkStep_FullRHS", "Mass matrix solver failure");
      return(ARK_MASSSOLVE_FAIL);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_TakeStep_Z:

  This routine serves the primary purpose of the ARKStep module:
  it performs a single ARK step (with embedding, if possible).
  This version solves for each ARK stage vector, z_i.

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
int arkStep_TakeStep_Z(void* arkode_mem, realtype *dsmPtr, int *nflagPtr)
{
  int retval, is, nvec;
  booleantype implicit_stage;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  N_Vector zcor0;
  realtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_TakeStep_Z",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* local shortcuts for use with fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* if problem will involve no algebraic solvers, initialize nflagPtr to success */
  if ((!step_mem->implicit) && (step_mem->mass_type == MASS_IDENTITY))
    *nflagPtr = ARK_SUCCESS;

  /* call nonlinear solver setup if it exists */
  if (step_mem->NLS)
    if ((step_mem->NLS)->ops->setup) {
      zcor0 = ark_mem->tempv3;
      N_VConst(ZERO, zcor0);    /* set guess to all 0 (since using predictor-corrector form) */
      retval = SUNNonlinSolSetup(step_mem->NLS, zcor0, ark_mem);
      if (retval < 0) return(ARK_NLS_SETUP_FAIL);
      if (retval > 0) return(ARK_NLS_SETUP_RECVR);
    }

  /* loop over internal stages to the step */
  for (is=0; is<step_mem->stages; is++) {

    /* store current stage index */
    step_mem->istage = is;

    /* set current stage time(s) */
    if (step_mem->implicit)
      ark_mem->tcur = ark_mem->tn + step_mem->Bi->c[is]*ark_mem->h;
    else
      ark_mem->tcur = ark_mem->tn + step_mem->Be->c[is]*ark_mem->h;

#ifdef SUNDIALS_DEBUG
    printf("step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n",
           ark_mem->nst, is, ark_mem->h, ark_mem->tcur);
#endif

    /* setup time-dependent mass matrix */
    if ((step_mem->mass_type == MASS_TIMEDEP) && (step_mem->msetup != NULL)) {
      retval = step_mem->msetup((void *) ark_mem, ark_mem->tcur,
                                ark_mem->tempv1, ark_mem->tempv2,
                                ark_mem->tempv3);
      if (retval != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
    }

    /* determine whether implicit solve is required */
    implicit_stage = SUNFALSE;
    if (step_mem->implicit)
      if (SUNRabs(step_mem->Bi->A[is][is]) > TINY)
        implicit_stage = SUNTRUE;

    /* if implicit, call built-in and user-supplied predictors
       (results placed in zpred) */
    if (implicit_stage) {

      retval = arkStep_Predict(ark_mem, is, step_mem->zpred);
      if (retval != ARK_SUCCESS)  return (retval);

      /* if a user-supplied predictor routine is provided, call that here.
         Note that arkStep_Predict is *still* called, so this user-supplied
         routine can just 'clean up' the built-in prediction, if desired. */
      if (step_mem->stage_predict) {
        retval = step_mem->stage_predict(ark_mem->tcur, step_mem->zpred,
                                         ark_mem->user_data);
        if (retval < 0)  return(ARK_USER_PREDICT_FAIL);
        if (retval > 0)  return(TRY_AGAIN);
      }

    }

#ifdef SUNDIALS_DEBUG_PRINTVEC
    printf("predictor:\n");
    N_VPrint(step_mem->zpred);
#endif

    /* set up explicit data for evaluation of ARK stage (store in sdata) */
    retval = arkStep_StageSetup(ark_mem, implicit_stage);
    if (retval != ARK_SUCCESS)  return (retval);

#ifdef SUNDIALS_DEBUG_PRINTVEC
    printf("rhs data:\n");
    N_VPrint(step_mem->sdata);
#endif

    /* solver diagnostics reporting */
    if (ark_mem->report)
      fprintf(ark_mem->diagfp, "ARKStep  step  %li  %"RSYM"  %i  %"RSYM"\n",
              ark_mem->nst, ark_mem->h, is, ark_mem->tcur);

    /* perform implicit solve if required */
    if (implicit_stage) {

      /* implicit solve result is stored in ark_mem->ycur;
         return with positive value on anything but success */
      *nflagPtr = arkStep_Nls(ark_mem, *nflagPtr);
      if (*nflagPtr != ARK_SUCCESS)  return(TRY_AGAIN);

#ifdef SUNDIALS_DEBUG_PRINTVEC
      printf("nonlinear solution:\n");
      N_VPrint(ark_mem->ycur);
#endif

    /* otherwise no implicit solve is needed */
    } else {

      /* if M is fixed, solve with it to compute update (place back in sdata) */
      if (step_mem->mass_type == MASS_FIXED) {

        /* perform solve; return with positive value on anything but success */
        *nflagPtr = step_mem->msolve((void *) ark_mem, step_mem->sdata,
                                     step_mem->nlscoef);
        if (*nflagPtr != ARK_SUCCESS)  return(TRY_AGAIN);

      }

      /* set y to be yn + sdata (either computed in arkStep_StageSetup,
         or updated in prev. block) */
      N_VLinearSum(ONE, ark_mem->yn, ONE, step_mem->sdata, ark_mem->ycur);

#ifdef SUNDIALS_DEBUG_PRINTVEC
      printf("explicit solution:\n");
      N_VPrint(ark_mem->ycur);
#endif

    }

    /* apply user-supplied stage postprocessing function (if supplied) */
    /* NOTE: with internally inconsistent IMEX methods (c_i^E != c_i^I) the value
       of tcur corresponds to the stage time from the implicit table (c_i^I). */
    if (ark_mem->ProcessStage != NULL) {
      retval = ark_mem->ProcessStage(ark_mem->tcur,
                                     ark_mem->ycur,
                                     ark_mem->user_data);
      if (retval != 0) return(ARK_POSTPROCESS_STAGE_FAIL);
    }

    /* successful stage solve */
    /*    store implicit RHS (value in Fi[is] is from preceding nonlinear iteration) */
    if (step_mem->implicit) {
      retval = step_mem->fi(ark_mem->tcur, ark_mem->ycur,
                            step_mem->Fi[is], ark_mem->user_data);
      step_mem->nfi++;
      if (retval < 0)  return(ARK_RHSFUNC_FAIL);
      if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      /* apply external polynomial forcing */
      if (step_mem->impforcing) {
        cvals[0] = ONE;
        Xvecs[0] = step_mem->Fi[is];
        nvec     = 1;
        arkStep_ApplyForcing(step_mem, ark_mem->tcur, ONE, &nvec);
        N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fi[is]);
      }
    }

    /*    store explicit RHS */
    if (step_mem->explicit) {
        retval = step_mem->fe(ark_mem->tn + step_mem->Be->c[is]*ark_mem->h,
                              ark_mem->ycur, step_mem->Fe[is], ark_mem->user_data);
        step_mem->nfe++;
        if (retval < 0)  return(ARK_RHSFUNC_FAIL);
        if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
        /* apply external polynomial forcing */
        if (step_mem->expforcing) {
          cvals[0] = ONE;
          Xvecs[0] = step_mem->Fe[is];
          nvec     = 1;
          arkStep_ApplyForcing(step_mem, ark_mem->tn+step_mem->Be->c[is]*ark_mem->h,
                               ONE, &nvec);
          N_VLinearCombination(nvec, cvals, Xvecs, step_mem->Fe[is]);
        }
    }

    /* if using a time-dependent mass matrix, update Fe[is] and/or Fi[is] with M(t)^{-1} */
    if (step_mem->mass_type == MASS_TIMEDEP) {
      if (step_mem->implicit) {
        *nflagPtr = step_mem->msolve((void *) ark_mem, step_mem->Fi[is], step_mem->nlscoef);
        if (*nflagPtr != ARK_SUCCESS)  return(TRY_AGAIN);
      }
      if (step_mem->explicit) {
        *nflagPtr = step_mem->msolve((void *) ark_mem, step_mem->Fe[is], step_mem->nlscoef);
        if (*nflagPtr != ARK_SUCCESS)  return(TRY_AGAIN);
      }
    }

  } /* loop over stages */

  /* compute time-evolved solution (in ark_ycur), error estimate (in dsm).
     This can fail recoverably due to nonconvergence of the mass matrix solve,
     so handle that appropriately. */
  if (step_mem->mass_type == MASS_FIXED) {
    retval = arkStep_ComputeSolutions_MassFixed(ark_mem, dsmPtr);
  } else {
    retval = arkStep_ComputeSolutions(ark_mem, dsmPtr);
  }
  if (retval < 0)  return(retval);
  if (retval > 0) {
    *nflagPtr = retval;
    return(TRY_AGAIN);
  }

#ifdef SUNDIALS_DEBUG
  printf("error estimate = %"RSYM"\n", *dsmPtr);
#endif

  /* solver diagnostics reporting */
  if (ark_mem->report)
    fprintf(ark_mem->diagfp, "ARKStep  etest  %li  %"RSYM"  %"RSYM"\n",
            ark_mem->nst, ark_mem->h, *dsmPtr);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/


/*---------------------------------------------------------------
  arkStep_AccessStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int arkStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeARKStepMem *step_mem)
{

  /* access ARKodeMem structure */
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    fname, MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *ark_mem = (ARKodeMem) arkode_mem;
  if ((*ark_mem)->step_mem==NULL) {
    arkProcessError(*ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    fname, MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *step_mem = (ARKodeARKStepMem) (*ark_mem)->step_mem;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype arkStep_CheckNVector(N_Vector tmpl)
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
  arkStep_SetButcherTables

  This routine determines the ERK/DIRK/ARK method to use, based
  on the desired accuracy and information on whether the problem
  is explicit, implicit or imex.
  ---------------------------------------------------------------*/
int arkStep_SetButcherTables(ARKodeMem ark_mem)
{
  int etable, itable;
  ARKodeARKStepMem step_mem;
  sunindextype Blrw, Bliw;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_SetButcherTables", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* if tables have already been specified, just return */
  if ( (step_mem->Be != NULL) || (step_mem->Bi != NULL) )
    return(ARK_SUCCESS);

  /* initialize table numbers to illegal values */
  etable = itable = -1;

  /**** ImEx methods ****/
  if (step_mem->explicit && step_mem->implicit) {

    switch (step_mem->q) {

    case(2):
    case(3):
      etable = DEFAULT_ARK_ETABLE_3;
      itable = DEFAULT_ARK_ITABLE_3;
      break;
    case(4):
      etable = DEFAULT_ARK_ETABLE_4;
      itable = DEFAULT_ARK_ITABLE_4;
      break;
    case(5):
      etable = DEFAULT_ARK_ETABLE_5;
      itable = DEFAULT_ARK_ITABLE_5;
      break;
    default:    /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_SetButcherTables",
                      "No ImEx method at requested order, using q=5.");
      etable = DEFAULT_ARK_ETABLE_5;
      itable = DEFAULT_ARK_ITABLE_5;
      break;
    }

  /**** implicit methods ****/
  } else if (step_mem->implicit) {

    switch (step_mem->q) {
    case(2):
      itable = DEFAULT_DIRK_2;
      break;
    case(3):
      itable = DEFAULT_DIRK_3;
      break;
    case(4):
      itable = DEFAULT_DIRK_4;
      break;
    case(5):
      itable = DEFAULT_DIRK_5;
      break;
    default:    /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_SetButcherTables",
                      "No implicit method at requested order, using q=5.");
      itable = DEFAULT_DIRK_5;
      break;
    }

  /**** explicit methods ****/
  } else {

    switch (step_mem->q) {
    case(2):
      etable = DEFAULT_ERK_2;
      break;
    case(3):
      etable = DEFAULT_ERK_3;
      break;
    case(4):
      etable = DEFAULT_ERK_4;
      break;
    case(5):
      etable = DEFAULT_ERK_5;
      break;
    case(6):
      etable = DEFAULT_ERK_6;
      break;
    case(7):
    case(8):
      etable = DEFAULT_ERK_8;
      break;
    default:    /* no available method, set default */
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_SetButcherTables",
                      "No explicit method at requested order, using q=6.");
      etable = DEFAULT_ERK_6;
      break;
    }

  }

  if (etable > -1)
    step_mem->Be = ARKodeButcherTable_LoadERK(etable);
  if (itable > -1)
    step_mem->Bi = ARKodeButcherTable_LoadDIRK(itable);

  /* note Butcher table space requirements */
  ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  /* set [redundant] ARK stored values for stage numbers and method orders */
  if (step_mem->Be != NULL) {
    step_mem->stages = step_mem->Be->stages;
    step_mem->q = step_mem->Be->q;
    step_mem->p = step_mem->Be->p;
  }
  if (step_mem->Bi != NULL) {
    step_mem->stages = step_mem->Bi->stages;
    step_mem->q = step_mem->Bi->q;
    step_mem->p = step_mem->Bi->p;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_CheckButcherTables

  This routine runs through the explicit and/or implicit Butcher
  tables to ensure that they meet all necessary requirements,
  including:
    strictly lower-triangular (ERK)
    lower-triangular with some nonzeros on diagonal (IRK)
    method order q > 0 (all)
    embedding order q > 0 (all -- if adaptive time-stepping enabled)
    stages > 0 (all)

  Returns ARK_SUCCESS if tables pass, ARK_INVALID_TABLE otherwise.
  ---------------------------------------------------------------*/
int arkStep_CheckButcherTables(ARKodeMem ark_mem)
{
  int i, j;
  booleantype okay;
  ARKodeARKStepMem step_mem;
  realtype tol = RCONST(1.0e-12);

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* check that the expected tables are set */
  if (step_mem->explicit && step_mem->Be == NULL) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "explicit table is NULL!");
    return(ARK_INVALID_TABLE);
  }

  if (step_mem->implicit && step_mem->Bi == NULL) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "implicit table is NULL!");
    return(ARK_INVALID_TABLE);
  }

  /* check that stages > 0 */
  if (step_mem->stages < 1) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "stages < 1!");
    return(ARK_INVALID_TABLE);
  }

  /* check that method order q > 0 */
  if (step_mem->q < 1) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "method order < 1!");
    return(ARK_INVALID_TABLE);
  }

  /* check that embedding order p > 0 */
  if ((step_mem->p < 1) && (!ark_mem->fixedstep)) {
    arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "embedding order < 1!");
    return(ARK_INVALID_TABLE);
  }

  /* check that embedding exists */
  if ((step_mem->p > 0) && (!ark_mem->fixedstep)) {
    if (step_mem->implicit) {
      if (step_mem->Bi->d == NULL) {
        arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                        "arkStep_CheckButcherTables",
                        "no implicit embedding!");
        return(ARK_INVALID_TABLE);
      }
    }
    if (step_mem->explicit) {
      if (step_mem->Be->d == NULL) {
        arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                        "arkStep_CheckButcherTables",
                        "no explicit embedding!");
        return(ARK_INVALID_TABLE);
      }
    }
  }

  /* check that ERK table is strictly lower triangular */
  if (step_mem->explicit) {
    okay = SUNTRUE;
    for (i=0; i<step_mem->stages; i++)
      for (j=i; j<step_mem->stages; j++)
        if (SUNRabs(step_mem->Be->A[i][j]) > tol)
          okay = SUNFALSE;
    if (!okay) {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                      "arkStep_CheckButcherTables",
                      "Ae Butcher table is implicit!");
      return(ARK_INVALID_TABLE);
    }
  }

  /* check that IRK table is implicit and lower triangular */
  if (step_mem->implicit) {
    okay = SUNFALSE;
    for (i=0; i<step_mem->stages; i++)
      if (SUNRabs(step_mem->Bi->A[i][i]) > tol)
        okay = SUNTRUE;
    if (!okay) {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                      "arkStep_CheckButcherTables",
                      "Ai Butcher table is explicit!");
      return(ARK_INVALID_TABLE);
    }

    okay = SUNTRUE;
    for (i=0; i<step_mem->stages; i++)
      for (j=i+1; j<step_mem->stages; j++)
        if (SUNRabs(step_mem->Bi->A[i][j]) > tol)
          okay = SUNFALSE;
    if (!okay) {
      arkProcessError(ark_mem, ARK_INVALID_TABLE, "ARKode::ARKStep",
                      "arkStep_CheckButcherTables",
                      "Ai Butcher table has entries above diagonal!");
      return(ARK_INVALID_TABLE);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_Predict

  This routine computes the prediction for a specific internal
  stage solution, storing the result in yguess.  The
  prediction is done using the interpolation structure in
  extrapolation mode, hence stages "far" from the previous time
  interval are predicted using lower order polynomials than the
  "nearby" stages.
  ---------------------------------------------------------------*/
int arkStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess)
{
  int i, retval, jstage, nvec;
  realtype tau;
  realtype h;
  ARKodeARKStepMem step_mem;
  realtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Predict", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* verify that interpolation structure is provided */
  if ((ark_mem->interp == NULL) && (step_mem->predictor > 0) && (step_mem->predictor < 4)) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Predict",
                    "Interpolation structure is NULL");
    return(ARK_MEM_NULL);
  }

  /* local shortcuts for use with fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* if the first step, use initial condition as guess */
  if (ark_mem->initsetup) {
    N_VScale(ONE, ark_mem->yn, yguess);
    return(ARK_SUCCESS);
  }

  /* set evaluation time tau as relative shift from previous successful time */
  tau = step_mem->Bi->c[istage]*ark_mem->h/ark_mem->hold;

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
      jstage = (step_mem->Bi->c[i] != ZERO) ? i : jstage;

    /* if using the trivial predictor, break */
    if (jstage == -1)  break;

    /* find the "optimal" previous stage to use */
    for (i=0; i<istage; i++)
      if ( (step_mem->Bi->c[i] > step_mem->Bi->c[jstage]) &&
           (step_mem->Bi->c[i] != ZERO) )
        jstage = i;

    /* set stage time, stage RHS and interpolation values */
    h = ark_mem->h * step_mem->Bi->c[jstage];
    tau = ark_mem->h * step_mem->Bi->c[istage];
    nvec = 0;
    if (step_mem->implicit) {    /* Implicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fi[jstage];
      nvec += 1;
    }
    if (step_mem->explicit) {    /* Explicit piece */
      cvals[nvec] = ONE;
      Xvecs[nvec] = step_mem->Fe[jstage];
      nvec += 1;
    }

    /* call predictor routine */
    retval = arkPredict_Bootstrap(ark_mem, h, tau, nvec, cvals, Xvecs, yguess);
    if (retval != ARK_ILL_INPUT)  return(retval);
    break;

  case 5:

    /***** Minimal correction predictor: use all previous stage
           information in this step *****/

    /* set arrays for fused vector operation */
    nvec = 0;
    if (step_mem->explicit) {       /* Explicit pieces */
      for (jstage=0; jstage<istage; jstage++) {
        cvals[nvec] = ark_mem->h * step_mem->Be->A[istage][jstage];
        Xvecs[nvec] = step_mem->Fe[jstage];
        nvec += 1;
      }
    }
    if (step_mem->implicit) {      /* Implicit pieces */
      for (jstage=0; jstage<istage; jstage++) {
        cvals[nvec] = ark_mem->h * step_mem->Bi->A[istage][jstage];
        Xvecs[nvec] = step_mem->Fi[jstage];
        nvec += 1;
      }
    }
    cvals[nvec] = ONE;
    Xvecs[nvec] = ark_mem->yn;
    nvec += 1;

    /* compute predictor */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yguess);
    if (retval != 0) return(ARK_VECTOROP_ERR);
    return(ARK_SUCCESS);
    break;

  }

  /* if we made it here, use the trivial predictor (previous step solution) */
  N_VScale(ONE, ark_mem->yn, yguess);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_StageSetup

  This routine sets up the stage data for computing the RK
  residual, along with the step- and method-related factors
  gamma, gammap and gamrat.

  The internal behavior of this setup depends on two factors:
  (a) whether the stage is explicit or implicit, and
  (b) the form of mass matrix (I, M, M(t)).

  In each of the following:
  * yn is the previous time step solution,
  * r corresponds to the nonlinear residual,
  * z=zp+zc corresponds to the updated stage solution, where
    zp is the predictor (zp=yn for explicit stages), and zc
    is the corrector,
  * i is the current stage index, and
  * g = h*Ai(i,i) is the implicit coefficient.

  Explicit, I:
      z = yn + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
     <=>
      zc = h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes zc, stored in step_mem->sdata.

  Implicit, I:
      z = yn - h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
             - h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
     <=>
      r = zc - g*Fi(i) - s
      s = yn - zp + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes s, stored in step_mem->sdata.

  Explicit, M:
      M*z = M*yn + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
     <=>
      M*zc = s
      s = h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes s, stored in step_mem->sdata.

  Implicit, M:
      M*z = M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                 + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
     <=>
      r = M*zc - g*Fi(i) - s
      s = M*(yn - zp) + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes s, stored in step_mem->sdata.

  Explicit, M(t):
      z = yn + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j)+Ai(i,j)*Fi(j))
     <=>
      zc = h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j)+Ai(i,j)*Fi(j))
    This routine computes zc, stored in step_mem->sdata.

  Implicit, M(t):
      M(t)*z = M(t)*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                       + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
     <=>
      r = M(t)*(zc - s) - g*Fi(i)
      s = yn - zp + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
    This routine computes s, stored in step_mem->sdata.


  Thus _internal_ to this routine, we have 3 modes:

  Explicit (any):
    sdata = h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
  Implicit, M:
    sdata = M*(yn - zp) + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))
  Implicit, I or M(t):
    sdata = yn - zp + h*sum_{j=0}^{i-1} (Ae(i,j)*Fe(j) + Ai(i,j)*Fi(j))

  ---------------------------------------------------------------*/
int arkStep_StageSetup(ARKodeMem ark_mem, booleantype implicit)
{
  /* local data */
  ARKodeARKStepMem step_mem;
  int retval, i, j, nvec;
  realtype* cvals;
  N_Vector* Xvecs;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_StageSetup", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* Set shortcut to current stage index */
  i = step_mem->istage;

  /* If this is the first stage, and explicit, just set sdata=0 and return */
  if (!implicit && (i==0)) {
    N_VConst(ZERO, step_mem->sdata);
    return (ARK_SUCCESS);
  }

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* Update gamma if stage is implicit */
  if (implicit) {
    step_mem->gamma = ark_mem->h * step_mem->Bi->A[i][i];
    if (ark_mem->firststage)
      step_mem->gammap = step_mem->gamma;
    step_mem->gamrat = (ark_mem->firststage) ?
      ONE : step_mem->gamma / step_mem->gammap;  /* protect x/x != 1.0 */
  }

  /* If predictor==5, then sdata=0 (plus any implicit forcing).
     Set sdata appropriately and return */
  if (implicit && (step_mem->predictor == 5)) {

    /* apply external polynomial forcing (updates nvec, cvals, Xvecs) */
    if (step_mem->impforcing) {
      nvec = 0;
      arkStep_ApplyForcing(step_mem, ark_mem->tcur, step_mem->gamma, &nvec);
      retval = N_VLinearCombination(nvec, cvals, Xvecs, step_mem->sdata);
      if (retval != 0) return(ARK_VECTOROP_ERR);
    } else {
      N_VConst(ZERO, step_mem->sdata);
    }
    return (ARK_SUCCESS);

  }

  /* If implicit, initialize sdata to yn - zpred (here: zpred = zp), and set
     first entries for eventual N_VLinearCombination call */
  nvec = 0;
  if (implicit) {
    N_VLinearSum(ONE, ark_mem->yn, -ONE, step_mem->zpred, step_mem->sdata);
    cvals[0] = ONE;
    Xvecs[0] = step_mem->sdata;
    nvec = 1;
  }

  /* If implicit with fixed M!=I, update sdata with M*sdata */
  if (implicit && (step_mem->mass_type == MASS_FIXED)) {
    N_VScale(ONE, step_mem->sdata, ark_mem->tempv1);
    retval = step_mem->mmult((void *) ark_mem, ark_mem->tempv1, step_mem->sdata);
    if (retval != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
  }

  /* Update sdata with prior stage information */
  if (step_mem->explicit) {   /* Explicit pieces */
    for (j=0; j<i; j++) {
      cvals[nvec] = ark_mem->h * step_mem->Be->A[i][j];
      Xvecs[nvec] = step_mem->Fe[j];
      nvec += 1;
    }
  }
  if (step_mem->implicit) {   /* Implicit pieces */
    for (j=0; j<i; j++) {
      cvals[nvec] = ark_mem->h * step_mem->Bi->A[i][j];
      Xvecs[nvec] = step_mem->Fi[j];
      nvec += 1;
    }
  }

  /* apply external polynomial forcing (updates nvec, cvals, Xvecs) */
  if (step_mem->impforcing) {
    arkStep_ApplyForcing(step_mem, ark_mem->tcur,
                         ark_mem->h * step_mem->Bi->A[i][i], &nvec);
  }

  /* call fused vector operation to do the work */
  retval = N_VLinearCombination(nvec, cvals, Xvecs, step_mem->sdata);
  if (retval != 0) return(ARK_VECTOROP_ERR);

  /* return with success */
  return (ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_ComputeSolutions

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This routine
  also computes the error estimate ||y-ytilde||_WRMS, where ytilde
  is the embedded solution, and the norm weights come from
  ark_ewt.  This norm value is returned.  The vector form of this
  estimated error (y-ytilde) is stored in ark_mem->tempv1, in case
  the calling routine wishes to examine the error locations.

  This version assumes either an identity or time-dependent mass
  matrix (identical steps).
  ---------------------------------------------------------------*/
int arkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsmPtr)
{
  /* local data */
  int retval, j, nvec;
  N_Vector y, yerr;
  realtype* cvals;
  N_Vector* Xvecs;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_ComputeSolutions", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set N_Vector shortcuts, and shortcut to time at end of step */
  y    = ark_mem->ycur;
  yerr = ark_mem->tempv1;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* initialize output */
  *dsmPtr = ZERO;

  /* Compute time step solution */
  /*   set arrays for fused vector operation */
  cvals[0] = ONE;
  Xvecs[0] = ark_mem->yn;
  nvec = 1;
  for (j=0; j<step_mem->stages; j++) {
    if (step_mem->explicit) {      /* Explicit pieces */
      cvals[nvec] = ark_mem->h * step_mem->Be->b[j];
      Xvecs[nvec] = step_mem->Fe[j];
      nvec += 1;
    }
    if (step_mem->implicit) {      /* Implicit pieces */
      cvals[nvec] = ark_mem->h * step_mem->Bi->b[j];
      Xvecs[nvec] = step_mem->Fi[j];
      nvec += 1;
    }
  }

  /*   call fused vector operation to do the work */
  retval = N_VLinearCombination(nvec, cvals, Xvecs, y);
  if (retval != 0) return(ARK_VECTOROP_ERR);

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep) {

    /* set arrays for fused vector operation */
    nvec = 0;
    for (j=0; j<step_mem->stages; j++) {
      if (step_mem->explicit) {        /* Explicit pieces */
        cvals[nvec] = ark_mem->h * (step_mem->Be->b[j] - step_mem->Be->d[j]);
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
      if (step_mem->implicit) {        /* Implicit pieces */
        cvals[nvec] = ark_mem->h * (step_mem->Bi->b[j] - step_mem->Bi->d[j]);
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }
    }

    /* call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yerr);
    if (retval != 0) return(ARK_VECTOROP_ERR);

    /* fill error norm */
    *dsmPtr = N_VWrmsNorm(yerr, ark_mem->ewt);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_ComputeSolutions_MassFixed

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This routine
  also computes the error estimate ||y-ytilde||_WRMS, where ytilde
  is the embedded solution, and the norm weights come from
  ark_ewt.  This norm value is returned.  The vector form of this
  estimated error (y-ytilde) is stored in ark_mem->tempv1, in case
  the calling routine wishes to examine the error locations.

  This version assumes a fixed mass matrix.
  ---------------------------------------------------------------*/
int arkStep_ComputeSolutions_MassFixed(ARKodeMem ark_mem, realtype *dsmPtr)
{
  /* local data */
  int retval, j, nvec;
  N_Vector y, yerr;
  realtype* cvals;
  N_Vector* Xvecs;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_ComputeSolutions_MassFixed", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set N_Vector shortcuts, and shortcut to time at end of step */
  y    = ark_mem->ycur;
  yerr = ark_mem->tempv1;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* initialize output */
  *dsmPtr = ZERO;

  /* compute y RHS (store in y) */
  /*   set arrays for fused vector operation */
  nvec = 0;
  for (j=0; j<step_mem->stages; j++) {
    if (step_mem->explicit) {      /* Explicit pieces */
      cvals[nvec] = ark_mem->h * step_mem->Be->b[j];
      Xvecs[nvec] = step_mem->Fe[j];
      nvec += 1;
    }
    if (step_mem->implicit) {      /* Implicit pieces */
      cvals[nvec] = ark_mem->h * step_mem->Bi->b[j];
      Xvecs[nvec] = step_mem->Fi[j];
      nvec += 1;
    }
  }

  /*   call fused vector operation to compute RHS */
  retval = N_VLinearCombination(nvec, cvals, Xvecs, y);
  if (retval != 0) return(ARK_VECTOROP_ERR);

  /* solve for y update (stored in y) */
  retval = step_mem->msolve((void *) ark_mem, y, step_mem->nlscoef);
  if (retval < 0) {
    *dsmPtr = RCONST(2.0);   /* indicate too much error, step with smaller step */
    N_VScale(ONE, ark_mem->yn, y);      /* place old solution into y */
    return(CONV_FAIL);
  }

  /* compute y = yn + update */
  N_VLinearSum(ONE, ark_mem->yn, ONE, y, y);


  /* compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep) {

    /* compute yerr RHS vector */
    /*   set arrays for fused vector operation */
    nvec = 0;
    for (j=0; j<step_mem->stages; j++) {
      if (step_mem->explicit) {        /* Explicit pieces */
        cvals[nvec] = ark_mem->h * (step_mem->Be->b[j] - step_mem->Be->d[j]);
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
      if (step_mem->implicit) {        /* Implicit pieces */
        cvals[nvec] = ark_mem->h * (step_mem->Bi->b[j] - step_mem->Bi->d[j]);
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }
    }

    /*   call fused vector operation to compute yerr RHS */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yerr);
    if (retval != 0) return(ARK_VECTOROP_ERR);

    /* solve for yerr */
    retval = step_mem->msolve((void *) ark_mem, yerr, step_mem->nlscoef);
    if (retval < 0) {
      *dsmPtr = RCONST(2.0);  /* next attempt will reduce step by 'etacf';
                                 insert dsmPtr placeholder here */
      return(CONV_FAIL);
    }
    /* fill error norm */
    *dsmPtr = N_VWrmsNorm(yerr, ark_mem->ewt);
  }

  return(ARK_SUCCESS);
}


/*------------------------------------------------------------------------------
  arkStep_ApplyForcing

  Determines the linear combination coefficients and vectors to apply forcing
  at a given value of the independent variable (t).  This occurs through
  appending coefficients and N_Vector pointers to the underlying cvals and Xvecs
  arrays in the step_mem structure.  The dereferenced input *nvec should indicate
  the next available entry in the cvals/Xvecs arrays.  The input 's' is a
  scaling factor that should be applied to each of these coefficients.
  ----------------------------------------------------------------------------*/

void arkStep_ApplyForcing(ARKodeARKStepMem step_mem, realtype t,
                          realtype s, int *nvec)
{
  realtype tau, taui;
  int i;

  /* always append the constant forcing term */
  step_mem->cvals[*nvec] = s;
  step_mem->Xvecs[*nvec] = step_mem->forcing[0];
  (*nvec) += 1;

  /* compute normalized time tau and initialize tau^i */
  tau  = (t - step_mem->tshift) / (step_mem->tscale);
  taui = tau;
  for (i=1; i<step_mem->nforcing; i++) {
    step_mem->cvals[*nvec] = s*taui;
    step_mem->Xvecs[*nvec] = step_mem->forcing[i];
    taui *= tau;
    (*nvec) += 1;
  }
}

/*------------------------------------------------------------------------------
  arkStep_SetInnerForcing

  Sets an array of coefficient vectors for a time-dependent external polynomial
  forcing term in the ODE RHS i.e., y' = fe(t,y) + fi(t,y) + p(t). This
  function is primarily intended for use with multirate integration methods
  (e.g., MRIStep) where ARKStep is used to solve a modified ODE at a fast time
  scale. The polynomial is of the form

  p(t) = sum_{i = 0}^{nvecs - 1} forcing[i] * ((t - tshift) / (tscale))^i

  where tshift and tscale are used to normalize the time t (e.g., with MRIGARK
  methods).
  ----------------------------------------------------------------------------*/

int arkStep_SetInnerForcing(void* arkode_mem, realtype tshift, realtype tscale,
                            N_Vector* forcing, int nvecs)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_SetInnerForcing",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  if (nvecs > 0) {

    /* enable forcing */
    if (step_mem->explicit) {
      step_mem->expforcing = SUNTRUE;
      step_mem->impforcing = SUNFALSE;
    } else {
      step_mem->expforcing = SUNFALSE;
      step_mem->impforcing = SUNTRUE;
    }
    step_mem->tshift   = tshift;
    step_mem->tscale   = tscale;
    step_mem->forcing  = forcing;
    step_mem->nforcing = nvecs;

    /* If cvals and Xvecs are not allocated then arkStep_Init has not been
       called and the number of stages has not been set yet. These arrays will
       be allocated in arkStep_Init and take into account the value of nforcing.
       On subsequent calls will check if enough space has allocated in case
       nforcing has increased since the original allocation. */
    if (step_mem->cvals != NULL && step_mem->Xvecs != NULL) {

      /* check if there are enough reusable arrays for fused operations */
      if ((step_mem->nfusedopvecs - nvecs) < (2 * step_mem->stages + 2)) {

        /* free current work space */
        if (step_mem->cvals != NULL) {
          free(step_mem->cvals);
          ark_mem->lrw -= step_mem->nfusedopvecs;
        }
        if (step_mem->Xvecs != NULL) {
          free(step_mem->Xvecs);
          ark_mem->liw -= step_mem->nfusedopvecs;
        }

        /* allocate reusable arrays for fused vector operations */
        step_mem->nfusedopvecs = 2 * step_mem->stages + 2 + nvecs;

        step_mem->cvals = NULL;
        step_mem->cvals = (realtype *) calloc(step_mem->nfusedopvecs,
                                              sizeof(realtype));
        if (step_mem->cvals == NULL) return(ARK_MEM_FAIL);
        ark_mem->lrw += step_mem->nfusedopvecs;

        step_mem->Xvecs = NULL;
        step_mem->Xvecs = (N_Vector *) calloc(step_mem->nfusedopvecs,
                                              sizeof(N_Vector));
        if (step_mem->Xvecs == NULL) return(ARK_MEM_FAIL);
        ark_mem->liw += step_mem->nfusedopvecs;
      }
    }

  } else {

    /* disable forcing */
    step_mem->expforcing = SUNFALSE;
    step_mem->impforcing = SUNFALSE;
    step_mem->tshift     = ZERO;
    step_mem->tscale     = ONE;
    step_mem->forcing    = NULL;
    step_mem->nforcing   = 0;

  }

  return(0);
}

/*===============================================================
  EOF
  ===============================================================*/
