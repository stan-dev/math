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
 * This is the implementation file for ARKode's ARK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_arkstep_impl.h"
#include <sundials/sundials_math.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

#define NO_DEBUG_OUTPUT
/* #define DEBUG_OUTPUT */
#ifdef DEBUG_OUTPUT
#include <nvector/nvector_serial.h>
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
  ark_mem->step                = arkStep_TakeStep;
  ark_mem->step_mem            = (void*) step_mem;

  /* Set default values for ARKStep optional inputs */
  retval = ARKStepSetDefaults((void *)ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ARKStep",
                    "ARKStepCreate",
                    "Error setting default solver options");
    return(NULL);
  }

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Allocate the general ARK stepper vectors using y0 as a template */
  /* NOTE: Fe, Fi, cvals and Xvecs will be allocated later on
     (based on the number of ARK stages) */

  /* Clone the input vector to create sdata, zpred and zcor */
  if (!arkAllocVec(ark_mem, y0, &(step_mem->sdata)))
    return(NULL);
  if (!arkAllocVec(ark_mem, y0, &(step_mem->zpred)))
    return(NULL);
  if (!arkAllocVec(ark_mem, y0, &(step_mem->zcor)))
    return(NULL);

  /* Copy the input parameters into ARKode state */
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Update the ARKode workspace requirements */
  ark_mem->liw += 41;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->lrw += 10;

  /* Allocate step adaptivity structure, set default values, note storage */
  step_mem->hadapt_mem = arkAdaptInit();
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "ARKStepCreate",
                    "Allocation of step adaptivity structure failed");
    return(NULL);
  }
  ark_mem->lrw += ARK_ADAPT_LRW;
  ark_mem->liw += ARK_ADAPT_LIW;

  /* If an implicit component is to be solved, create default Newton NLS object */
  step_mem->ownNLS = SUNFALSE;
  if (step_mem->implicit)  {
    NLS = NULL;
    NLS = SUNNonlinSol_Newton(y0);
    if (NLS == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                      "ARKStepCreate", "Error creating default Newton solver");
      return(NULL);
    }
    retval = ARKStepSetNonlinearSolver(ark_mem, NLS);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep",
                      "ARKStepCreate", "Error attaching default Newton solver");
      return(NULL);
    }
    step_mem->ownNLS = SUNTRUE;
  }

  /* Set the linear solver addresses to NULL (we check != NULL later) */
  step_mem->linit  = NULL;
  step_mem->lsetup = NULL;
  step_mem->lsolve = NULL;
  step_mem->lfree  = NULL;
  step_mem->lmem   = NULL;
  step_mem->lsolve_type = -1;

  /* Set the mass matrix solver addresses to NULL */
  step_mem->minit    = NULL;
  step_mem->msetup   = NULL;
  step_mem->mmult    = NULL;
  step_mem->msolve   = NULL;
  step_mem->mfree    = NULL;
  step_mem->mass_mem = NULL;
  step_mem->msetuptime = -RCONST(99999999999.0);
  step_mem->msolve_type = -1;

  /* Initialize initial error norm  */
  step_mem->eRNrm = 1.0;

  /* Initialize all the counters */
  step_mem->nst_attempts = 0;
  step_mem->nfe          = 0;
  step_mem->nfi          = 0;
  step_mem->ncfn         = 0;
  step_mem->netf         = 0;
  step_mem->nsetups      = 0;
  step_mem->nstlp        = 0;

  /* Initialize main ARKode infrastructure */
  retval = arkInit(ark_mem, t0, y0);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ARKStep", "ARKStepCreate",
                    "Unable to initialize main ARKode infrastructure");
    return(NULL);
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
  if (step_mem->sdata != NULL) {
    retval = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &step_mem->sdata);
    if (retval != ARK_SUCCESS)  return(retval);
  }
  if (step_mem->zpred != NULL) {
    retval = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &step_mem->zpred);
    if (retval != ARK_SUCCESS)  return(retval);
  }
  if (step_mem->zcor != NULL) {
    retval = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &step_mem->zcor);
    if (retval != ARK_SUCCESS)  return(retval);
  }

  /* Resize the ARKStep vectors */
  /*     Fe */
  if (step_mem->Fe != NULL) {
    for (i=0; i<step_mem->stages; i++) {
      retval = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                         liw_diff, y0, &step_mem->Fe[i]);
      if (retval != ARK_SUCCESS)  return(retval);
    }
  }
  /*     Fi */
  if (step_mem->Fi != NULL) {
    for (i=0; i<step_mem->stages; i++) {
      retval = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                         liw_diff, y0, &step_mem->Fi[i]);
      if (retval != ARK_SUCCESS)  return(retval);
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
  if (step_mem->NLS != NULL) {
    step_mem->ncfn    = 0;
    step_mem->nsetups = 0;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  ARKStepReInit:

  This routine re-initializes the ARKStep module to solve a new
  problem of the same size as was previously solved.
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

  /* Check that at least one of fe, fi is supplied and is to be used */
  if (fe == NULL && fi == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepCreate", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* ReInitialize main ARKode infrastructure */
  retval = arkReInit(ark_mem, t0, y0);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ARKStep", "ARKStepReInit",
                    "Unable to initialize main ARKode infrastructure");
    return(retval);
  }

  /* Set implicit/explicit problem based on function pointers */
  step_mem->explicit = (fe == NULL) ? SUNFALSE : SUNTRUE;
  step_mem->implicit = (fi == NULL) ? SUNFALSE : SUNTRUE;

  /* Copy the input parameters into ARKode state */
  step_mem->fe = fe;
  step_mem->fi = fi;

  /* Destroy/Reinitialize time step adaptivity structure (if present) */
  if (step_mem->hadapt_mem != NULL) {
    free(step_mem->hadapt_mem);
    step_mem->hadapt_mem = arkAdaptInit();
    if (step_mem->hadapt_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "ARKStepReInit",
                      "Allocation of Step Adaptivity Structure Failed");
      return(ARK_MEM_FAIL);
    }
  }
  /* Initialize initial error norm  */
  step_mem->eRNrm = 1.0;

  /* Initialize all the counters */
  step_mem->nst_attempts = 0;
  step_mem->nfe          = 0;
  step_mem->nfi          = 0;
  step_mem->ncfn         = 0;
  step_mem->netf         = 0;
  step_mem->nsetups      = 0;
  step_mem->nstlp        = 0;

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

    /* free the time step adaptivity module */
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
      ark_mem->lrw -= ARK_ADAPT_LRW;
      ark_mem->liw -= ARK_ADAPT_LIW;
    }

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
      ark_mem->lrw -= (2*step_mem->stages + 1);
    }
    if (step_mem->Xvecs != NULL) {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= (2*step_mem->stages + 1);
    }

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
  fprintf(outfile,"ARKStep: mnewt = %i\n", step_mem->mnewt);
  fprintf(outfile,"ARKStep: maxcor = %i\n", step_mem->maxcor);
  fprintf(outfile,"ARKStep: maxnef = %i\n", step_mem->maxnef);
  fprintf(outfile,"ARKStep: maxncf = %i\n", step_mem->maxncf);
  fprintf(outfile,"ARKStep: msbp = %i\n", step_mem->msbp);
  fprintf(outfile,"ARKStep: predictor = %i\n", step_mem->predictor);
  fprintf(outfile,"ARKStep: lsolve_type = %i\n", step_mem->lsolve_type);
  fprintf(outfile,"ARKStep: msolve_type = %i\n", step_mem->msolve_type);
  fprintf(outfile,"ARKStep: convfail = %i\n", step_mem->convfail);

  /* output long integer quantities */
  fprintf(outfile,"ARKStep: nst_attempts = %li\n", step_mem->nst_attempts);
  fprintf(outfile,"ARKStep: nfe = %li\n", step_mem->nfe);
  fprintf(outfile,"ARKStep: nfi = %li\n", step_mem->nfi);
  fprintf(outfile,"ARKStep: ncfn = %li\n", step_mem->ncfn);
  fprintf(outfile,"ARKStep: netf = %li\n", step_mem->netf);
  fprintf(outfile,"ARKStep: nsetups = %li\n", step_mem->nsetups);
  fprintf(outfile,"ARKStep: nstlp = %li\n", step_mem->nstlp);

  /* output boolean quantities */
  fprintf(outfile,"ARKStep: user_linear = %i\n", step_mem->linear);
  fprintf(outfile,"ARKStep: user_linear_timedep = %i\n", step_mem->linear_timedep);
  fprintf(outfile,"ARKStep: user_explicit = %i\n", step_mem->explicit);
  fprintf(outfile,"ARKStep: user_implicit = %i\n", step_mem->implicit);
  fprintf(outfile,"ARKStep: hadapt_pq = %i\n", step_mem->hadapt_pq);
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
  if (step_mem->hadapt_mem != NULL) {
    fprintf(outfile,"ARKStep: timestep adaptivity structure:\n");
    arkPrintAdaptMem(step_mem->hadapt_mem, outfile);
  }
  fprintf(outfile,"ARKStep: crdown = %"RSYM"\n", step_mem->crdown);
  fprintf(outfile,"ARKStep: rdiv = %"RSYM"\n", step_mem->rdiv);
  fprintf(outfile,"ARKStep: dgmax = %"RSYM"\n", step_mem->dgmax);

#ifdef DEBUG_OUTPUT
  /* output vector quantities */
  if (step_mem->sdata != NULL) {
    fprintf(outfile, "ARKStep: sdata:\n");
    N_VPrint_Serial(step_mem->sdata);
  }
  if (step_mem->zpred != NULL) {
    fprintf(outfile, "ARKStep: zpred:\n");
    N_VPrint_Serial(step_mem->zpred);
  }
  if (step_mem->zcor != NULL) {
    fprintf(outfile, "ARKStep: zcor:\n");
    N_VPrint_Serial(step_mem->zcor);
  }
  if (step_mem->Fe != NULL)
    for (i=0; i<step_mem->stages; i++) {
      fprintf(outfile,"ARKStep: Fe[%i]:\n", i);
      N_VPrint_Serial(step_mem->Fe[i]);
    }
  if (step_mem->Fi != NULL)
    for (i=0; i<step_mem->stages; i++) {
      fprintf(outfile,"ARKStep: Fi[%i]:\n", i);
      N_VPrint_Serial(step_mem->Fi[i]);
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
                         int lsolve_type, void *lmem)
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
int arkStep_AttachMasssol(void* arkode_mem, ARKMassInitFn minit,
                          ARKMassSetupFn msetup,
                          ARKMassMultFn mmult,
                          ARKMassSolveFn msolve,
                          ARKMassFreeFn mfree,
                          int msolve_type, void *mass_mem)
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
  within arkInitialSetup (init_type == 0) or arkPostResizeSetup
  (init_type == 1).

  With init_type == 0, this routine:
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

  With init_type == 1, this routine:
  - checks for consistency between the system and mass matrix
    linear solvers (if applicable)
  - initializes and sets up the system and mass matrix linear
    solvers (if applicable)
  - initializes and sets up the nonlinear solver (if applicable)
  ---------------------------------------------------------------*/
int arkStep_Init(void* arkode_mem, int init_type)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int j, retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_Init",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* perform initializations specific to init_type 0 */
  if (init_type == 0) {

    /* destroy adaptivity structure if fixed-stepping is requested */
    if (ark_mem->fixedstep)
      if (step_mem->hadapt_mem != NULL) {
        free(step_mem->hadapt_mem);
        step_mem->hadapt_mem = NULL;
      }

    /* Set first step growth factor */
    if (step_mem->hadapt_mem != NULL)
      step_mem->hadapt_mem->etamax = step_mem->hadapt_mem->etamx1;

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

    /* note Butcher table space requirements */
    ARKodeButcherTable_Space(step_mem->Be, &Bliw, &Blrw);
    ark_mem->liw += Bliw;
    ark_mem->lrw += Blrw;
    ARKodeButcherTable_Space(step_mem->Bi, &Bliw, &Blrw);
    ark_mem->liw += Bliw;
    ark_mem->lrw += Blrw;

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

    /* Allocate reusable arrays for fused vector interface */
    j = (2*step_mem->stages+1 > 4) ? 2*step_mem->stages+1 : 4;
    if (step_mem->cvals == NULL) {
      step_mem->cvals = (realtype *) calloc(j, sizeof(realtype));
      if (step_mem->cvals == NULL)  return(ARK_MEM_FAIL);
      ark_mem->lrw += j;
    }
    if (step_mem->Xvecs == NULL) {
      step_mem->Xvecs = (N_Vector *) calloc(j, sizeof(N_Vector));
      if (step_mem->Xvecs == NULL)  return(ARK_MEM_FAIL);
      ark_mem->liw += j;   /* pointers */
    }

    /* Allocate interpolation memory (if unallocated, and needed) */
    if ((ark_mem->interp == NULL) && (step_mem->predictor > 0)) {
      ark_mem->interp = arkInterpCreate(ark_mem);
      if (ark_mem->interp == NULL) {
        arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ARKStep", "arkStep_Init",
                        "Unable to allocate interpolation structure");
        return(ARK_MEM_FAIL);
      }
    }

  }

  /* Check for consistency between linear system modules
       (e.g., if lsolve is direct, msolve needs to match) */
  if (step_mem->mass_mem != NULL) {  /* M != I */
    if (step_mem->lsolve_type != step_mem->msolve_type) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep", "arkStep_Init",
                      "Incompatible linear and mass matrix solvers");
      return(ARK_ILL_INPUT);
    }
  }

  /* Perform mass matrix solver initialization and setup (if applicable) */
  if (step_mem->mass_mem != NULL) {

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
      retval = step_mem->msetup((void *) ark_mem, ark_mem->tempv1,
                             ark_mem->tempv2, ark_mem->tempv3);
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_MASSSETUP_FAIL, "ARKode::ARKStep",
                        "arkStep_Init", MSG_ARK_MASSSETUP_FAIL);
        return(ARK_MASSSETUP_FAIL);
        step_mem->msetuptime = ark_mem->tcur;
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
  int retval;
  booleantype recomputeRHS;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_FullRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if the problem involves a non-identity mass matrix and setup is
     required, do so here (use output f as a temporary) */
  if ( (step_mem->mass_mem != NULL) && (step_mem->msetup != NULL) )
    if (SUNRabs(step_mem->msetuptime - t) > FUZZ_FACTOR*ark_mem->uround) {
      retval = step_mem->msetup((void *) ark_mem, f, ark_mem->tempv2,
                                ark_mem->tempv3);
      if (retval != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
      step_mem->msetuptime = t;
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
  if (step_mem->mass_mem != NULL) {
    retval = step_mem->msolve((void *) ark_mem, f, step_mem->nlscoef/ark_mem->h);
    if (retval != ARK_SUCCESS) {
      arkProcessError(ark_mem, ARK_MASSSOLVE_FAIL, "ARKode::ARKStep",
                      "arkStep_FullRHS", "Mass matrix solver failure");
      return(ARK_MASSSOLVE_FAIL);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_TakeStep:

  This routine serves the primary purpose of the ARKStep module:
  it performs a single successful ARK step (with embedding, if
  possible).  Multiple attempts may be taken in this process --
  once a step completes with successful (non)linear solves at
  each stage and passes the error estimate, the routine returns
  successfully.  If it cannot do so, it returns with an
  appropriate error flag.
  ---------------------------------------------------------------*/
int arkStep_TakeStep(void* arkode_mem)
{
  realtype dsm;
  int retval, ncf, nef, is, nflag, kflag, eflag;
  booleantype implicit_stage;
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  N_Vector zcor0;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_TakeStep",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  ncf = nef = 0;
  nflag = FIRST_CALL;
  eflag = ARK_SUCCESS;
  kflag = SOLVE_SUCCESS;

  /* Looping point for attempts to take a step */
  for(;;) {

    /* increment attempt counter */
    step_mem->nst_attempts++;

    /* call nonlinear solver setup if it exists */
    if (step_mem->NLS)
      if ((step_mem->NLS)->ops->setup) {
        zcor0 = ark_mem->tempv3;
        N_VConst(ZERO, zcor0);    /* set guess to all 0 (since ARKode uses predictor-corrector form) */
        retval = SUNNonlinSolSetup(step_mem->NLS, zcor0, ark_mem);
        if (retval < 0) return(ARK_NLS_SETUP_FAIL);
        if (retval > 0) return(ARK_NLS_SETUP_RECVR);
      }

    /* Loop over internal stages to the step */
    for (is=0; is<step_mem->stages; is++) {

      /* store current stage index */
      step_mem->istage = is;

      /* Set current stage time(s) */
      if (step_mem->implicit)
        ark_mem->tcur = ark_mem->tn + step_mem->Bi->c[is]*ark_mem->h;
      else
        ark_mem->tcur = ark_mem->tn + step_mem->Be->c[is]*ark_mem->h;

#ifdef DEBUG_OUTPUT
 printf("step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n",
         ark_mem->nst, is, ark_mem->h, ark_mem->tcur);
#endif

      /* determine whether implicit solve is required */
      implicit_stage = SUNFALSE;
      if (step_mem->implicit)
        if (SUNRabs(step_mem->Bi->A[is][is]) > TINY)
          implicit_stage = SUNTRUE;

      /* Call predictor for current stage solution (result placed in zpred) */
      if (implicit_stage) {
        eflag = arkStep_Predict(ark_mem, is, step_mem->zpred);
        if (eflag != ARK_SUCCESS)  return (eflag);
      } else {
        N_VScale(ONE, ark_mem->yn, step_mem->zpred);
      }

#ifdef DEBUG_OUTPUT
 printf("predictor:\n");
 N_VPrint_Serial(step_mem->zpred);
#endif

      /* Set up data for evaluation of ARK stage residual (data stored in sdata) */
      eflag = arkStep_StageSetup(ark_mem);
      if (eflag != ARK_SUCCESS)  return (eflag);

#ifdef DEBUG_OUTPUT
 printf("rhs data:\n");
 N_VPrint_Serial(step_mem->sdata);
#endif

      /* Solver diagnostics reporting */
      if (ark_mem->report)
        fprintf(ark_mem->diagfp, "ARKStep  step  %li  %"RSYM"  %i  %"RSYM"\n",
                ark_mem->nst, ark_mem->h, is, ark_mem->tcur);

      /* perform implicit solve if required */
      if (implicit_stage) {

        /* perform implicit solve (result is stored in ark_mem->ycur) */
        nflag = arkStep_Nls(ark_mem, nflag);

#ifdef DEBUG_OUTPUT
 printf("nonlinear solution:\n");
 N_VPrint_Serial(ark_mem->ycur);
#endif

        /* check for convergence (on failure, h will have been modified) */
        kflag = arkStep_HandleNFlag(ark_mem, &nflag, &ncf);

        /* If fixed time-stepping is used, then anything other than a
           successful solve must result in an error */
        if (ark_mem->fixedstep && (kflag != SOLVE_SUCCESS))
          return(kflag);

        /* If h reduced and step needs to be retried, break loop */
        if (kflag == PREDICT_AGAIN) break;

        /* Return if nonlinear solve failed and recovery not possible. */
        if (kflag != SOLVE_SUCCESS) return(kflag);

      /* otherwise no implicit solve is needed */
      } else {

        /* if M!=I, solve with M to compute update (place back in sdata) */
        if (step_mem->mass_mem != NULL) {

          /* perform mass matrix solve */
          nflag = step_mem->msolve((void *) ark_mem, step_mem->sdata,
                                   step_mem->nlscoef);

          /* check for convergence (on failure, h will have been modified) */
          kflag = arkStep_HandleNFlag(ark_mem, &nflag, &ncf);

          /* If fixed time-stepping is used, then anything other than a
             successful solve must result in an error */
          if (ark_mem->fixedstep && (kflag != SOLVE_SUCCESS))
            return(kflag);

          /* If h reduced and step needs to be retried, break loop */
          if (kflag == PREDICT_AGAIN) break;

          /* Return if solve failed and recovery not possible. */
          if (kflag != SOLVE_SUCCESS) return(kflag);

        /* if M==I, set y to be zpred + RHS data computed in arkStep_StageSetup */
        } else {
          N_VLinearSum(ONE, step_mem->sdata, ONE,
                       step_mem->zpred, ark_mem->ycur);
        }

#ifdef DEBUG_OUTPUT
 printf("explicit solution:\n");
 N_VPrint_Serial(ark_mem->ycur);
#endif

      }

      /* successful stage solve */
      /*    store implicit RHS (value in Fi[is] is from preceding nonlinear iteration) */
      if (step_mem->implicit) {
        retval = step_mem->fi(ark_mem->tcur, ark_mem->ycur,
                              step_mem->Fi[is], ark_mem->user_data);
        step_mem->nfi++;
        if (retval < 0)  return(ARK_RHSFUNC_FAIL);
        if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      }

      /*    store explicit RHS if necessary
            (already computed at first stage of purely explicit runs) */
      if (step_mem->explicit) {
          retval = step_mem->fe(ark_mem->tn + step_mem->Be->c[is]*ark_mem->h,
                                ark_mem->ycur, step_mem->Fe[is], ark_mem->user_data);
          step_mem->nfe++;
          if (retval < 0)  return(ARK_RHSFUNC_FAIL);
          if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
      }

    } /* loop over stages */

    /* if h has changed due to convergence failure and a new
       prediction is needed, continue to next attempt at step
       (cannot occur if fixed time stepping is enabled) */
    if (kflag == PREDICT_AGAIN)  continue;

    /* compute time-evolved solution (in ark_ycur), error estimate (in dsm) */
    retval = arkStep_ComputeSolutions(ark_mem, &dsm);
    if (retval < 0)  return(retval);

#ifdef DEBUG_OUTPUT
 printf("error estimate = %"RSYM"\n", dsm);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->report)
      fprintf(ark_mem->diagfp, "ARKStep  etest  %li  %"RSYM"  %"RSYM"\n",
              ark_mem->nst, ark_mem->h, dsm);

    /* Perform time accuracy error test (if failure, updates h for next try) */
    if (!ark_mem->fixedstep)
      eflag = arkStep_DoErrorTest(ark_mem, &nflag, &nef, dsm);

#ifdef DEBUG_OUTPUT
 printf("error test flag = %i\n", eflag);
#endif

    /* Restart step attempt (recompute all stages) if error test fails recoverably */
    if (eflag == TRY_AGAIN)  continue;

    /* Return if error test failed and recovery not possible. */
    if (eflag != ARK_SUCCESS)  return(eflag);

    /* Error test passed (eflag=ARK_SUCCESS), break from loop */
    break;

  } /* loop over step attempts */


  /* The step has completed successfully, clean up and
     consider change of step size */
  retval = arkStep_PrepareNextStep(ark_mem, dsm);
  if (retval != ARK_SUCCESS)  return(retval);

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
      itable = DEFAULT_DIRK_5;;
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

  Returns ARK_SUCCESS if tables pass, ARK_ILL_INPUT otherwise.
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

  /* check that stages > 0 */
  if (step_mem->stages < 1) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "stages < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that method order q > 0 */
  if (step_mem->q < 1) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "method order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that embedding order p > 0 */
  if ((step_mem->p < 1) && (!ark_mem->fixedstep)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_CheckButcherTables",
                    "embedding order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that embedding exists */
  if ((step_mem->p > 0) && (!ark_mem->fixedstep)) {
    if (step_mem->implicit) {
      if (step_mem->Bi->d == NULL) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                        "arkStep_CheckButcherTables",
                        "no implicit embedding!");
        return(ARK_ILL_INPUT);
      }
    }
    if (step_mem->explicit) {
      if (step_mem->Be->d == NULL) {
        arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                        "arkStep_CheckButcherTables",
                        "no explicit embedding!");
        return(ARK_ILL_INPUT);
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
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_CheckButcherTables",
                      "Ae Butcher table is implicit!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that IRK table is implicit and lower triangular */
  if (step_mem->implicit) {
    okay = SUNFALSE;
    for (i=0; i<step_mem->stages; i++)
      if (SUNRabs(step_mem->Bi->A[i][i]) > tol)
        okay = SUNTRUE;
    if (!okay) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_CheckButcherTables",
                      "Ai Butcher table is explicit!");
      return(ARK_ILL_INPUT);
    }

    okay = SUNTRUE;
    for (i=0; i<step_mem->stages; i++)
      for (j=i+1; j<step_mem->stages; j++)
        if (SUNRabs(step_mem->Bi->A[i][j]) > tol)
          okay = SUNFALSE;
    if (!okay) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                      "arkStep_CheckButcherTables",
                      "Ai Butcher table has entries above diagonal!");
      return(ARK_ILL_INPUT);
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_Predict

  This routine computes the prediction for a specific internal
  stage solution, storing the result in ypred.  The
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
  if ((ark_mem->interp == NULL) && (step_mem->predictor > 0)) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Predict",
                    "Interpolation structure is NULL");
    return(ARK_MEM_NULL);
  }

  /* local shortcuts to fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* if the first step (or if resized), use initial condition as guess */
  if (ark_mem->nst == 0 || ark_mem->resized) {
    N_VScale(ONE, ark_mem->yn, yguess);
    return(ARK_SUCCESS);
  }

  /* set evaluation time tau relative shift from previous successful time */
  tau = step_mem->Bi->c[istage]*ark_mem->h/ark_mem->hold;

  /* use requested predictor formula */
  switch (step_mem->predictor) {

  case 1:

    /***** Interpolatory Predictor 1 -- all to max order *****/
    retval = arkPredict_MaximumOrder(ark_mem, tau, yguess);
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 2:

    /***** Interpolatory Predictor 2 -- decrease order w/ increasing level of extrapolation *****/
    retval = arkPredict_VariableOrder(ark_mem, tau, yguess);
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 3:

    /***** Cutoff predictor: max order interpolatory output for stages "close"
           to previous step, first-order predictor for subsequent stages *****/
    retval = arkPredict_CutoffOrder(ark_mem, tau, yguess);
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 4:

    /***** Bootstrap predictor: if any previous stage in step has nonzero c_i,
           construct a quadratic Hermite interpolant for prediction; otherwise
           use the trivial predictor.  The actual calculations are performed in
           arkPredict_Bootstrap, but here we need to determine the appropriate
           stage, c_j, to use. *****/

    /* this approach will not work (for now) when using a non-identity mass matrix */
    if (step_mem->mass_mem) break;

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
    if (retval == ARK_SUCCESS)  return(ARK_SUCCESS);
    break;

  case 5:

    /***** Minimal correction predictor: use all previous stage
           information in this step *****/

    /* this approach will not work (for now) when using a non-identity mass matrix */
    if (step_mem->mass_mem != NULL)  {
      N_VScale(ONE, ark_mem->yn, yguess);
      break;
    }

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

  At the ith stage, we compute the residual vector:
    r = -M*z + M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                    + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
    r = -M*(zp + zc) + M*yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                            + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
    r = (-M*zc + gamma*Fi(zi)) + (M*(yn - zp) + data)
  where z = zp + zc.  In the above form of the residual,
  the first group corresponds to the current solution
  correction, and the second group corresponds to existing data.
  This routine computes this existing data, (M*(yn - zp) + data)
  and stores in step_mem->sdata.
  ---------------------------------------------------------------*/
int arkStep_StageSetup(ARKodeMem ark_mem)
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

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* If predictor==5, then sdata=0, otherwise set sdata appropriately */
  if ( (step_mem->predictor == 5) && (step_mem->mass_mem == NULL) ) {

    N_VConst(ZERO, step_mem->sdata);

  } else {

    /* Initialize sdata to ycur - zpred (here: ycur = yn and zpred = zp) */
    N_VLinearSum(ONE, ark_mem->yn, -ONE, step_mem->zpred,
                 step_mem->sdata);

    /* If M!=I, replace sdata with M*sdata, so that sdata = M*(yn-zpred) */
    if (step_mem->mass_mem != NULL) {
      N_VScale(ONE, step_mem->sdata, ark_mem->tempv1);
      retval = step_mem->mmult((void *) ark_mem, ark_mem->tempv1, step_mem->sdata);
      if (retval != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
    }

    /* Update rhs with prior stage information */
    /*   set arrays for fused vector operation */
    cvals[0] = ONE;
    Xvecs[0] = step_mem->sdata;
    nvec = 1;
    if (step_mem->explicit)    /* Explicit pieces */
      for (j=0; j<i; j++) {
        cvals[nvec] = ark_mem->h * step_mem->Be->A[i][j];
        Xvecs[nvec] = step_mem->Fe[j];
        nvec += 1;
      }
    if (step_mem->implicit)    /* Implicit pieces */
      for (j=0; j<i; j++) {
        cvals[nvec] = ark_mem->h * step_mem->Bi->A[i][j];
        Xvecs[nvec] = step_mem->Fi[j];
        nvec += 1;
      }

    /*   call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, step_mem->sdata);
    if (retval != 0) return(ARK_VECTOROP_ERR);

  }

  /* Update gamma (if the method contains an implicit component) */
  if (step_mem->implicit) {
    step_mem->gamma = ark_mem->h * step_mem->Bi->A[i][i];
    if (ark_mem->firststage)
      step_mem->gammap = step_mem->gamma;
    step_mem->gamrat = (ark_mem->firststage) ?
      ONE : step_mem->gamma / step_mem->gammap;  /* protect x/x != 1.0 */
  }

  /* return with success */
  return (ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_HandleNFlag

  This routine takes action on the return value nflag = *nflagPtr
  returned by arkNls, as follows:

  If arkStep_Nls succeeded in solving the nonlinear system, then
  arkHandleNFlag returns the constant SOLVE_SUCCESS, which tells
  arkStep it is safe to continue with other stage solves, or to
  perform the error test.

  If the nonlinear system was not solved successfully, then ncfn and
  ncf = *ncfPtr are incremented.

  If the solution of the nonlinear system failed due to an
  unrecoverable failure by setup, we return the value ARK_LSETUP_FAIL.

  If it failed due to an unrecoverable failure in solve, then we return
  the value ARK_LSOLVE_FAIL.

  If it failed due to an unrecoverable failure in rhs, then we return
  the value ARK_RHSFUNC_FAIL.

  Otherwise, a recoverable failure occurred when solving the
  nonlinear system (arkNls returned nflag == CONV_FAIL or RHSFUNC_RECVR).
  In this case, if using fixed time step sizes, or if ncf is now equal
  to maxncf, or if |h| = hmin, then we return the value ARK_CONV_FAILURE
  (if nflag=CONV_FAIL) or ARK_REPTD_RHSFUNC_ERR (if nflag=RHSFUNC_RECVR).
  If not, we set *nflagPtr = PREV_CONV_FAIL and return the value
  PREDICT_AGAIN, telling arkStep to reattempt the step.
  ---------------------------------------------------------------*/
int arkStep_HandleNFlag(ARKodeMem ark_mem, int *nflagPtr, int *ncfPtr)
{
  int nflag;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_HandleNFlag", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  nflag = *nflagPtr;

  if (nflag == ARK_SUCCESS) return(SOLVE_SUCCESS);

  /* The nonlinear soln. failed; increment ncfn */
  step_mem->ncfn++;

  /* If fixed time stepping, then return with convergence failure */
  if (ark_mem->fixedstep)    return(ARK_CONV_FAILURE);

  /* Otherwise, access adaptivity structure */
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep", "arkStep_HandleNFlag",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (nflag == ARK_LSETUP_FAIL)  return(ARK_LSETUP_FAIL);
  if (nflag == ARK_LSOLVE_FAIL)  return(ARK_LSOLVE_FAIL);
  if (nflag == ARK_RHSFUNC_FAIL) return(ARK_RHSFUNC_FAIL);

  /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
  (*ncfPtr)++;
  hadapt_mem->etamax = ONE;

  /* If we had maxncf failures, or if |h| = hmin,
     return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */
  if ((*ncfPtr == step_mem->maxncf) ||
      (SUNRabs(ark_mem->h) <= ark_mem->hmin*ONEPSM)) {
    if (nflag == CONV_FAIL)     return(ARK_CONV_FAILURE);
    if (nflag == RHSFUNC_RECVR) return(ARK_REPTD_RHSFUNC_ERR);
  }

  /* Reduce step size; return to reattempt the step */
  ark_mem->eta = SUNMAX(hadapt_mem->etacf,
                        ark_mem->hmin / SUNRabs(ark_mem->h));
  ark_mem->h *= ark_mem->eta;
  ark_mem->next_h = ark_mem->h;
  *nflagPtr = PREV_CONV_FAIL;

  return(PREDICT_AGAIN);
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
  ---------------------------------------------------------------*/
int arkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsm)
{
  /* local data */
  realtype tend;
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
  tend = ark_mem->tn + ark_mem->h;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* initialize output */
  *dsm = ZERO;

  /* Compute updated solution and error estimate based on whether
     a non-identity mass matrix is present */
  if (step_mem->mass_mem != NULL) {   /* M != I */

    /* setup mass matrix */
    if (step_mem->msetup != NULL)
      if (SUNRabs(step_mem->msetuptime - tend) > FUZZ_FACTOR*ark_mem->uround) {
        retval = step_mem->msetup((void *) ark_mem, ark_mem->tempv1,
                               ark_mem->tempv2, ark_mem->tempv3);
        if (retval != ARK_SUCCESS)  return(ARK_MASSSETUP_FAIL);
        step_mem->msetuptime = tend;
      }

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
      *dsm = 2.0;         /* indicate too much error, step with smaller step */
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
        *dsm = 2.0;         /* indicate too much error, step with smaller step */
        return(CONV_FAIL);
      }
      /* fill error norm */
      *dsm = N_VWrmsNorm(yerr, ark_mem->ewt);
    }

  } else {                          /* M == I */

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
      *dsm = N_VWrmsNorm(yerr, ark_mem->ewt);
    }

  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_DoErrorTest

  This routine performs the local error test for the ARK method.
  The weighted local error norm dsm is passed in, and
  the test dsm ?<= 1 is made.

  If the test passes, arkDoErrorTest returns ARK_SUCCESS.

  If the test fails, we revert to the last successful solution
  time, and:
    - if maxnef error test failures have occurred or if
      SUNRabs(h) = hmin, we return ARK_ERR_FAILURE.
    - otherwise: update time step factor eta based on local error
      estimate and reduce h.  Then set *nflagPtr to PREV_ERR_FAIL,
      and return TRY_AGAIN.
  ---------------------------------------------------------------*/
int arkStep_DoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
                        int *nefPtr, realtype dsm)
{
  realtype ehist2, hhist2;
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_DoErrorTest", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep", "arkDoErrorTest",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */
  if (dsm <= ONE) return(ARK_SUCCESS);

  /* Test failed; increment counters, set nflag */
  (*nefPtr)++;
  step_mem->netf++;
  *nflagPtr = PREV_ERR_FAIL;

  /* At |h| = hmin or maxnef failures, return ARK_ERR_FAILURE */
  if ((SUNRabs(ark_mem->h) <= ark_mem->hmin*ONEPSM) ||
      (*nefPtr == step_mem->maxnef))
    return(ARK_ERR_FAILURE);

  /* Set etamax=1 to prevent step size increase at end of this step */
  hadapt_mem->etamax = ONE;

  /* Temporarily update error history array for recomputation of h */
  ehist2 = hadapt_mem->ehist[2];
  hadapt_mem->ehist[2] = hadapt_mem->ehist[1];
  hadapt_mem->ehist[1] = hadapt_mem->ehist[0];
  hadapt_mem->ehist[0] = dsm*hadapt_mem->bias;

  /* Temporarily update step history array for recomputation of h */
  hhist2 = hadapt_mem->hhist[2];
  hadapt_mem->hhist[2] = hadapt_mem->hhist[1];
  hadapt_mem->hhist[1] = hadapt_mem->hhist[0];
  hadapt_mem->hhist[0] = ark_mem->h;

  /* Compute accuracy-based time step estimate (updated ark_eta) */
  retval = arkAdapt((void*) ark_mem, step_mem->hadapt_mem, ark_mem->ycur,
                    ark_mem->tcur, ark_mem->h, step_mem->q, step_mem->p,
                    step_mem->hadapt_pq, ark_mem->nst);
  if (retval != ARK_SUCCESS)  return(ARK_ERR_FAILURE);

  /* Revert error history array */
  hadapt_mem->ehist[0] = hadapt_mem->ehist[1];
  hadapt_mem->ehist[1] = hadapt_mem->ehist[2];
  hadapt_mem->ehist[2] = ehist2;

  /* Revert step history array */
  hadapt_mem->hhist[0] = hadapt_mem->hhist[1];
  hadapt_mem->hhist[1] = hadapt_mem->hhist[2];
  hadapt_mem->hhist[2] = hhist2;

  /* Enforce failure bounds on eta, update h, and return for retry of step */
  if (*nefPtr >= hadapt_mem->small_nef)
    ark_mem->eta = SUNMIN(ark_mem->eta, hadapt_mem->etamxf);
  ark_mem->h *= ark_mem->eta;
  ark_mem->next_h = ark_mem->h;
  return(TRY_AGAIN);
}


/*---------------------------------------------------------------
  arkStep_PrepareNextStep

  This routine handles ARK-specific updates following a successful
  step: copying the ARK result to the current solution vector,
  updating the error/step history arrays, and setting the
  prospective step size, hprime, for the next step.  Along with
  hprime, it sets the ratio eta=hprime/h.  It also updates other
  state variables related to a change of step size.
  ---------------------------------------------------------------*/
int arkStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm)
{
  int retval;
  ARKodeARKStepMem step_mem;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_PrepareNextStep", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* Update step size and error history arrays */
  if (step_mem->hadapt_mem != NULL) {
    step_mem->hadapt_mem->ehist[2] = step_mem->hadapt_mem->ehist[1];
    step_mem->hadapt_mem->ehist[1] = step_mem->hadapt_mem->ehist[0];
    step_mem->hadapt_mem->ehist[0] = dsm*step_mem->hadapt_mem->bias;
    step_mem->hadapt_mem->hhist[2] = step_mem->hadapt_mem->hhist[1];
    step_mem->hadapt_mem->hhist[1] = step_mem->hadapt_mem->hhist[0];
    step_mem->hadapt_mem->hhist[0] = ark_mem->h;
  }

  /* If fixed time-stepping requested, defer
     step size changes until next step */
  if (ark_mem->fixedstep){
    ark_mem->hprime = ark_mem->h;
    ark_mem->eta = ONE;
    return(ARK_SUCCESS);
  }

  /* If etamax = 1, defer step size changes until next step,
     and reset etamax */
  if (step_mem->hadapt_mem != NULL)
    if (step_mem->hadapt_mem->etamax == ONE) {
      ark_mem->hprime = ark_mem->h;
      ark_mem->eta = ONE;
      step_mem->hadapt_mem->etamax = step_mem->hadapt_mem->growth;
      return(ARK_SUCCESS);
    }

  /* Adjust ark_eta in arkAdapt */
  if (step_mem->hadapt_mem != NULL) {
    retval = arkAdapt((void*) ark_mem, step_mem->hadapt_mem,
                      ark_mem->ycur, ark_mem->tn + ark_mem->h,
                      ark_mem->h, step_mem->q, step_mem->p,
                      step_mem->hadapt_pq, ark_mem->nst+1);
    if (retval != ARK_SUCCESS)  return(ARK_ERR_FAILURE);
  }

  /* Set hprime value for next step size */
  ark_mem->hprime = ark_mem->h * ark_mem->eta;

  /* Reset growth factor for subsequent time step */
  if (step_mem->hadapt_mem != NULL)
    step_mem->hadapt_mem->etamax = step_mem->hadapt_mem->growth;

  return(ARK_SUCCESS);
}


/*===============================================================
  EOF
  ===============================================================*/
