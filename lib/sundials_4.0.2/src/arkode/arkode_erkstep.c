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
 * This is the implementation file for ARKode's ERK time stepper
 * module.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_erkstep_impl.h"
#include <sundials/sundials_math.h>

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



/*===============================================================
  ERKStep Exported functions -- Required
  ===============================================================*/

void* ERKStepCreate(ARKRhsFn f, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  booleantype nvectorOK;
  int retval;

  /* Check that f is supplied */
  if (f == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "ERKStepCreate", MSG_ARK_NULL_F);
    return(NULL);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "ERKStepCreate", MSG_ARK_NULL_Y0);
    return(NULL);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = erkStep_CheckNVector(y0);
  if (!nvectorOK) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "ERKStepCreate", MSG_ARK_BAD_NVECTOR);
    return(NULL);
  }

  /* Create ark_mem structure and set default values */
  ark_mem = arkCreate();
  if (ark_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepCreate", MSG_ARK_NO_MEM);
    return(NULL);
  }

  /* Allocate ARKodeERKStepMem structure, and initialize to zero */
  step_mem = NULL;
  step_mem = (ARKodeERKStepMem) malloc(sizeof(struct ARKodeERKStepMemRec));
  if (step_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ERKStep",
                    "ERKStepCreate", MSG_ARK_ARKMEM_FAIL);
    return(NULL);
  }
  memset(step_mem, 0, sizeof(struct ARKodeERKStepMemRec));

  /* Attach step_mem structure and function pointers to ark_mem */
  ark_mem->step_init    = erkStep_Init;
  ark_mem->step_fullrhs = erkStep_FullRHS;
  ark_mem->step         = erkStep_TakeStep;
  ark_mem->step_mem     = (void*) step_mem;

  /* Set default values for ERKStep optional inputs */
  retval = ERKStepSetDefaults((void *) ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ERKStep",
                    "ERKStepCreate",
                    "Error setting default solver options");
    return(NULL);
  }

  /* Allocate the general ERK stepper vectors using y0 as a template */
  /* NOTE: F, cvals and Xvecs will be allocated later on
     (based on the number of ERK stages) */

  /* Copy the input parameters into ARKode state */
  step_mem->f = f;

  /* Update the ARKode workspace requirements -- UPDATE */
  ark_mem->liw += 41;  /* fcn/data ptr, int, long int, sunindextype, booleantype */
  ark_mem->lrw += 10;

  /* Allocate step adaptivity structure, set default values, note storage */
  step_mem->hadapt_mem = arkAdaptInit();
  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ERKStep", "ERKStepCreate",
                    "Allocation of step adaptivity structure failed");
    return(NULL);
  }
  ark_mem->lrw += ARK_ADAPT_LRW;
  ark_mem->liw += ARK_ADAPT_LIW;

  /* Initialize all the counters */
  step_mem->nst_attempts = 0;
  step_mem->nfe          = 0;
  step_mem->netf         = 0;

  /* Initialize main ARKode infrastructure */
  retval = arkInit(ark_mem, t0, y0);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ERKStep", "ERKStepCreate",
                    "Unable to initialize main ARKode infrastructure");
    return(NULL);
  }

  return((void *)ark_mem);
}


/*---------------------------------------------------------------
  ERKStepResize:

  This routine resizes the memory within the ERKStep module.
  It first resizes the main ARKode infrastructure memory, and
  then resizes its own data.
  ---------------------------------------------------------------*/
int ERKStepResize(void *arkode_mem, N_Vector y0, realtype hscale,
                  realtype t0, ARKVecResizeFn resize, void *resize_data)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  sunindextype lrw1, liw1, lrw_diff, liw_diff;
  int i, retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepReSize",
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

  /* resize ARKode infrastructure memory */
  retval = arkResize(ark_mem, y0, hscale, t0, resize, resize_data);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ERKStep", "ERKStepResize",
                    "Unable to resize main ARKode infrastructure");
    return(retval);
  }

  /* Resize the RHS vectors */
  for (i=0; i<step_mem->stages; i++) {
    retval = arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                       liw_diff, y0, &step_mem->F[i]);
    if (retval != ARK_SUCCESS)  return(retval);
  }

  return(ARK_SUCCESS);
}


int ERKStepReInit(void* arkode_mem, ARKRhsFn f, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepReInit",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Check that f is supplied */
  if (f == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "ERKStepReInit", MSG_ARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "ERKStepReInit", MSG_ARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* ReInitialize main ARKode infrastructure */
  retval = arkReInit(arkode_mem, t0, y0);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, retval, "ARKode::ERKStep", "ERKStepReInit",
                    "Unable to initialize main ARKode infrastructure");
    return(retval);
  }

  /* Copy the input parameters into ARKode state */
  step_mem->f = f;

  /* Destroy/Reinitialize time step adaptivity structure (if present) */
  if (step_mem->hadapt_mem != NULL) {
    free(step_mem->hadapt_mem);
    step_mem->hadapt_mem = arkAdaptInit();
    if (step_mem->hadapt_mem == NULL) {
      arkProcessError(ark_mem, ARK_MEM_FAIL, "ARKode::ERKStep", "ERKStepReInit",
                      "Allocation of Step Adaptivity Structure Failed");
      return(ARK_MEM_FAIL);
    }
  }

  /* Initialize all the counters */
  step_mem->nst_attempts = 0;
  step_mem->nfe          = 0;
  step_mem->netf         = 0;

  return(ARK_SUCCESS);
}


int ERKStepSStolerances(void *arkode_mem, realtype reltol, realtype abstol)
{
  /* unpack ark_mem, call arkSStolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSStolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSStolerances(ark_mem, reltol, abstol));
}


int ERKStepSVtolerances(void *arkode_mem, realtype reltol, N_Vector abstol)
{
  /* unpack ark_mem, call arkSVtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepSVtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSVtolerances(ark_mem, reltol, abstol));
}


int ERKStepWFtolerances(void *arkode_mem, ARKEwtFn efun)
{
  /* unpack ark_mem, call arkWFtolerances, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepWFtolerances", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkWFtolerances(ark_mem, efun));
}


int ERKStepRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)
{
  /* unpack ark_mem, call arkRootInit, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepRootInit", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkRootInit(ark_mem, nrtfn, g));
}


int ERKStepEvolve(void *arkode_mem, realtype tout, N_Vector yout,
                  realtype *tret, int itask)
{
  /* unpack ark_mem, call arkEvolve, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepEvolve", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkEvolve(ark_mem, tout, yout, tret, itask));
}


int ERKStepGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  /* unpack ark_mem, call arkGetDky, and return */
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::ERKStep",
                    "ERKStepGetDky", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkGetDky(ark_mem, t, k, dky));
}


/*---------------------------------------------------------------
  ERKStepFree frees all ERKStep memory, and then calls an ARKode
  utility routine to free the ARKode infrastructure memory.
  ---------------------------------------------------------------*/
void ERKStepFree(void **arkode_mem)
{
  int j;
  sunindextype Bliw, Blrw;
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;

  /* nothing to do if arkode_mem is already NULL */
  if (*arkode_mem == NULL)  return;

  /* conditional frees on non-NULL ERKStep module */
  ark_mem = (ARKodeMem) (*arkode_mem);
  if (ark_mem->step_mem != NULL) {

    step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

    /* free the time step adaptivity module */
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
      ark_mem->lrw -= ARK_ADAPT_LRW;
      ark_mem->liw -= ARK_ADAPT_LIW;
    }

    /* free the Butcher table */
    if (step_mem->B != NULL) {
      ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
      ARKodeButcherTable_Free(step_mem->B);
      step_mem->B = NULL;
      ark_mem->liw -= Bliw;
      ark_mem->lrw -= Blrw;
    }

    /* free the RHS vectors */
    if (step_mem->F != NULL) {
      for(j=0; j<step_mem->stages; j++)
        arkFreeVec(ark_mem, &step_mem->F[j]);
      free(step_mem->F);
      step_mem->F = NULL;
      ark_mem->liw -= step_mem->stages;
    }

    /* free the reusable arrays for fused vector interface */
    if (step_mem->cvals != NULL) {
      free(step_mem->cvals);
      step_mem->cvals = NULL;
      ark_mem->lrw -= (step_mem->stages + 1);
    }
    if (step_mem->Xvecs != NULL) {
      free(step_mem->Xvecs);
      step_mem->Xvecs = NULL;
      ark_mem->liw -= (step_mem->stages + 1);
    }

    /* free the time stepper module itself */
    free(ark_mem->step_mem);
    ark_mem->step_mem = NULL;

  }

  /* free memory for overall ARKode infrastructure */
  arkFree(arkode_mem);
}


/*---------------------------------------------------------------
  ERKStepPrintMem:

  This routine outputs the memory from the ERKStep structure and
  the main ARKode infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
void ERKStepPrintMem(void* arkode_mem, FILE* outfile)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int retval;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "ERKStepPrintMem",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return;

  /* output data from main ARKode infrastructure */
  arkPrintMem(ark_mem, outfile);

  /* output integer quantities */
  fprintf(outfile,"ERKStep: q = %i\n", step_mem->q);
  fprintf(outfile,"ERKStep: p = %i\n", step_mem->p);
  fprintf(outfile,"ERKStep: stages = %i\n", step_mem->stages);
  fprintf(outfile,"ERKStep: maxnef = %i\n", step_mem->maxnef);

  /* output long integer quantities */
  fprintf(outfile,"ERKStep: nst_attempts = %li\n", step_mem->nst_attempts);
  fprintf(outfile,"ERKStep: nfe = %li\n", step_mem->nfe);
  fprintf(outfile,"ERKStep: netf = %li\n", step_mem->netf);

  /* output boolean quantities */
  fprintf(outfile,"ERKStep: hadapt_pq = %i\n", step_mem->hadapt_pq);

  /* output realtype quantities */
  fprintf(outfile,"ERKStep: Butcher table:\n");
  ARKodeButcherTable_Write(step_mem->B, outfile);
  if (step_mem->hadapt_mem != NULL) {
    fprintf(outfile,"ERKStep: timestep adaptivity structure:\n");
    arkPrintAdaptMem(step_mem->hadapt_mem, outfile);
  }

#ifdef DEBUG_OUTPUT
  /* output vector quantities */
  for (i=0; i<step_mem->stages; i++) {
    fprintf(outfile,"ERKStep: F[%i]:\n", i);
    N_VPrint_Serial(step_mem->F[i]);
  }
#endif
}



/*===============================================================
  ERKStep Private functions
  ===============================================================*/

/*---------------------------------------------------------------
  Interface routines supplied to ARKode
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  erkStep_Init:

  This routine is called just prior to performing internal time 
  steps (after all user "set" routines have been called) from 
  within arkInitialSetup (init_type == 0) or arkPostResizeSetup
  (init_type == 1).

  With init_type == 0, this routine:
  - sets/checks the ARK Butcher tables to be used
  - allocates any memory that depends on the number of ARK
    stages, method order, or solver options

  With init_type == 1, this routine does nothing.
  ---------------------------------------------------------------*/
int erkStep_Init(void* arkode_mem, int init_type)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  sunindextype Blrw, Bliw;
  int retval, j;

  /* immediately return if init_type == 1 */
  if (init_type == 1)  return(ARK_SUCCESS);
  
  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "erkStep_Init",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* destroy adaptivity structure if fixed-stepping is requested */
  if (ark_mem->fixedstep)
    if (step_mem->hadapt_mem != NULL) {
      free(step_mem->hadapt_mem);
      step_mem->hadapt_mem = NULL;
    }

  /* Set first step growth factor */
  if (step_mem->hadapt_mem != NULL)
    step_mem->hadapt_mem->etamax = step_mem->hadapt_mem->etamx1;

  /* Create Butcher table (if not already set) */
  retval = erkStep_SetButcherTable(ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep", "erkStep_Init",
                    "Could not create Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* Check that Butcher table are OK */
  retval = erkStep_CheckButcherTable(ark_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "erkStep_Init", "Error in Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* note Butcher table space requirements */
  ARKodeButcherTable_Space(step_mem->B, &Bliw, &Blrw);
  ark_mem->liw += Bliw;
  ark_mem->lrw += Blrw;

  /* Allocate ARK RHS vector memory, update storage requirements */
  /*   Allocate F[0] ... F[stages-1] if needed */
  if (step_mem->F == NULL)
    step_mem->F = (N_Vector *) calloc(step_mem->stages, sizeof(N_Vector));
  for (j=0; j<step_mem->stages; j++) {
    if (!arkAllocVec(ark_mem, ark_mem->ewt, &(step_mem->F[j])))
      return(ARK_MEM_FAIL);
  }
  ark_mem->liw += step_mem->stages;  /* pointers */

  /* Allocate reusable arrays for fused vector interface */
  if (step_mem->cvals == NULL) {
    step_mem->cvals = (realtype *) calloc(step_mem->stages+1, sizeof(realtype));
    if (step_mem->cvals == NULL)  return(ARK_MEM_FAIL);
    ark_mem->lrw += (step_mem->stages + 1);
  }
  if (step_mem->Xvecs == NULL) {
    step_mem->Xvecs = (N_Vector *) calloc(step_mem->stages+1, sizeof(N_Vector));
    if (step_mem->Xvecs == NULL)  return(ARK_MEM_FAIL);
    ark_mem->liw += (step_mem->stages + 1);   /* pointers */
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_FullRHS:

  This is just a wrapper to call the user-supplied RHS function,
  f(t,y).

  This will be called in one of three 'modes':
    0 -> called at the beginning of a simulation
    1 -> called at the end of a successful step
    2 -> called elsewhere (e.g. for dense output)

  If it is called in mode 0, we store the vectors f(t,y) in F[0]
  for possible reuse in the first stage of the subsequent time step.

  If it is called in mode 1 and the method coefficients
  support it, we may just copy vectors F[stages] to fill f instead
  of calling f().

  Mode 2 is only called for dense output in-between steps, so we
  strive to store the intermediate parts so that they do not
  interfere with the other two modes.
  ---------------------------------------------------------------*/
int erkStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode)
{
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;
  int i, s, retval;
  booleantype recomputeRHS;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "erkStep_FullRHS",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* perform RHS functions contingent on 'mode' argument */
  switch(mode) {

  /* Mode 0: called at the beginning of a simulation
     Store the vectors f(t,y) in F[0] for possible reuse
     in the first stage of the subsequent time step */
  case 0:

    /* call f */
    retval = step_mem->f(t, y, step_mem->F[0], ark_mem->user_data);
    step_mem->nfe++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ERKStep",
                      "erkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    /* copy RHS vector into output */
    N_VScale(ONE, step_mem->F[0], f);

    break;


  /* Mode 1: called at the end of a successful step
     If the method coefficients support it, we just copy the last stage RHS vectors
     to fill f instead of calling f(t,y).
     Copy the results to F[0] if the coefficients support it. */
  case 1:

    /* determine if explicit/implicit RHS functions need to be recomputed */
    recomputeRHS = SUNFALSE;
    s = step_mem->B->stages;
    for (i=0; i<s; i++)
      if (SUNRabs(step_mem->B->b[i] - step_mem->B->A[s-1][i])>TINY)
        recomputeRHS = SUNTRUE;

    /* base RHS calls on recomputeRHS argument */
    if (recomputeRHS) {

      /* call f */
      retval = step_mem->f(t, y, step_mem->F[0], ark_mem->user_data);
      step_mem->nfe++;
      if (retval != 0) {
        arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ERKStep",
                        "erkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
        return(ARK_RHSFUNC_FAIL);
      }

    } else {
      N_VScale(ONE, step_mem->F[step_mem->stages-1], step_mem->F[0]);
    }

    /* copy RHS vector into output */
    N_VScale(ONE, step_mem->F[0], f);

    break;

  /*  Mode 2: called for dense output in-between steps
      store the intermediate calculations in such a way as to not
      interfere with the other two modes */
  default:

    /* call f */
    retval = step_mem->f(t, y, f, ark_mem->user_data);
    step_mem->nfe++;
    if (retval != 0) {
      arkProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKode::ERKStep",
                      "erkStep_FullRHS", MSG_ARK_RHSFUNC_FAILED, t);
      return(ARK_RHSFUNC_FAIL);
    }

    break;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_TakeStep:

  This routine serves the primary purpose of the ERKStep module:
  it performs a single successful embedded ERK step (if possible).
  Multiple attempts may be taken in this process -- once a step
  passes the error estimate, the routine returns successfully.
  If it cannot do so, it returns with an appropriate error flag.
  ---------------------------------------------------------------*/
int erkStep_TakeStep(void* arkode_mem)
{
  realtype dsm;
  int retval, nef, is, eflag, js, nvec;
  realtype* cvals;
  N_Vector* Xvecs;
  ARKodeMem ark_mem;
  ARKodeERKStepMem step_mem;

  /* access ARKodeERKStepMem structure */
  retval = erkStep_AccessStepMem(arkode_mem, "erkStep_TakeStep",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  nef = 0;
  eflag = ARK_SUCCESS;

  /* Looping point for attempts to take a step */
  for(;;) {

    /* increment attempt counter */
    step_mem->nst_attempts++;

#ifdef DEBUG_OUTPUT
 printf("stage 0 RHS:\n");
 N_VPrint_Serial(step_mem->F[0]);
#endif

    /* Loop over internal stages to the step; since the method is explicit
       the first stage RHS is just the full RHS from the start of the step */
    for (is=1; is<step_mem->stages; is++) {

      /* Set current stage time(s) */
      ark_mem->tcur = ark_mem->tn + step_mem->B->c[is]*ark_mem->h;

#ifdef DEBUG_OUTPUT
 printf("step %li,  stage %i,  h = %"RSYM",  t_n = %"RSYM"\n",
         ark_mem->nst, is, ark_mem->h, ark_mem->tcur);
#endif

      /* Solver diagnostics reporting */
      if (ark_mem->report)
        fprintf(ark_mem->diagfp, "ERKStep  step  %li  %"RSYM"  %i  %"RSYM"\n",
                ark_mem->nst, ark_mem->h, is, ark_mem->tcur);

      /* Set ycur to current stage solution */
      nvec = 0;
      for (js=0; js<is; js++) {
        cvals[nvec] = ark_mem->h * step_mem->B->A[is][js];
        Xvecs[nvec] = step_mem->F[js];
        nvec += 1;
      }
      cvals[nvec] = ONE;
      Xvecs[nvec] = ark_mem->yn;
      nvec += 1;

      /*   call fused vector operation to do the work */
      retval = N_VLinearCombination(nvec, cvals, Xvecs, ark_mem->ycur);
      if (retval != 0) return(ARK_VECTOROP_ERR);

      /* compute updated RHS */
      retval = step_mem->f(ark_mem->tcur, ark_mem->ycur,
                              step_mem->F[is], ark_mem->user_data);
      step_mem->nfe++;
      if (retval < 0)  return(ARK_RHSFUNC_FAIL);
      if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);

#ifdef DEBUG_OUTPUT
 printf("RHS:\n");
 N_VPrint_Serial(step_mem->F[is]);
#endif

    } /* loop over stages */

    /* compute time-evolved solution (in ark_ycur), error estimate (in dsm) */
    retval = erkStep_ComputeSolutions(ark_mem, &dsm);
    if (retval < 0)  return(retval);    /* msetup failure */

#ifdef DEBUG_OUTPUT
 printf("error estimate = %"RSYM"\n", dsm);
 printf("updated solution:\n");
 N_VPrint_Serial(ark_mem->ycur);
#endif

    /* Solver diagnostics reporting */
    if (ark_mem->report)
      fprintf(ark_mem->diagfp, "ERKStep  etest  %li  %"RSYM"  %"RSYM"\n",
              ark_mem->nst, ark_mem->h, dsm);

    /* Perform time accuracy error test (if failure, updates h for next try) */
    if (!ark_mem->fixedstep)
      eflag = erkStep_DoErrorTest(ark_mem, &nef, dsm);

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
  retval = erkStep_PrepareNextStep(ark_mem, dsm);
  if (retval != ARK_SUCCESS)  return(retval);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Internal utility routines
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  erkStep_AccessStepMem:

  Shortcut routine to unpack ark_mem and step_mem structures from
  void* pointer.  If either is missing it returns ARK_MEM_NULL.
  ---------------------------------------------------------------*/
int erkStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeERKStepMem *step_mem)
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
                    fname, MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  *step_mem = (ARKodeERKStepMem) (*ark_mem)->step_mem;
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_CheckNVector:

  This routine checks if all required vector operations are
  present.  If any of them is missing it returns SUNFALSE.
  ---------------------------------------------------------------*/
booleantype erkStep_CheckNVector(N_Vector tmpl)
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
  erkStep_SetButcherTable

  This routine determines the ERK method to use, based on the
  desired accuracy.
  ---------------------------------------------------------------*/
int erkStep_SetButcherTable(ARKodeMem ark_mem)
{
  int etable;
  ARKodeERKStepMem step_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "erkStep_SetButcherTable", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* if table has already been specified, just return */
  if (step_mem->B != NULL)
    return(ARK_SUCCESS);

  /* initialize table number to illegal values */
  etable = -1;

  /* select method based on order */
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
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "erkStep_SetButcherTable",
                    "No explicit method at requested order, using q=6.");
    etable = DEFAULT_ERK_6;
    break;
  }

  if (etable > -1)
    step_mem->B = ARKodeButcherTable_LoadERK(etable);

  /* set [redundant] stored values for stage numbers and method orders */
  if (step_mem->B != NULL) {
    step_mem->stages = step_mem->B->stages;
    step_mem->q = step_mem->B->q;
    step_mem->p = step_mem->B->p;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_CheckButcherTable

  This routine runs through the explicit Butcher table to ensure
  that it meets all necessary requirements, including:
    strictly lower-triangular (ERK)
    method order q > 0 (all)
    embedding order q > 0 (all -- if adaptive time-stepping enabled)
    stages > 0 (all)

  Returns ARK_SUCCESS if tables pass, ARK_ILL_INPUT otherwise.
  ---------------------------------------------------------------*/
int erkStep_CheckButcherTable(ARKodeMem ark_mem)
{
  int i, j;
  booleantype okay;
  ARKodeERKStepMem step_mem;
  realtype tol = RCONST(1.0e-12);

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "erkStep_CheckButcherTable", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* check that stages > 0 */
  if (step_mem->stages < 1) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "erkStep_CheckButcherTable",
                    "stages < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that method order q > 0 */
  if (step_mem->q < 1) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "erkStep_CheckButcherTable",
                    "method order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that embedding order p > 0 */
  if ((step_mem->p < 1) && (!ark_mem->fixedstep)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "erkStep_CheckButcherTable",
                    "embedding order < 1!");
    return(ARK_ILL_INPUT);
  }

  /* check that embedding exists */
  if ((step_mem->p > 0) && (!ark_mem->fixedstep)) {
    if (step_mem->B->d == NULL) {
      arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                      "erkStep_CheckButcherTable",
                      "no embedding!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that ERK table is strictly lower triangular */
  okay = SUNTRUE;
  for (i=0; i<step_mem->stages; i++)
    for (j=i; j<step_mem->stages; j++)
      if (SUNRabs(step_mem->B->A[i][j]) > tol)
        okay = SUNFALSE;
  if (!okay) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ERKStep",
                    "erkStep_CheckButcherTable",
                    "Ae Butcher table is implicit!");
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_ComputeSolutions

  This routine calculates the final RK solution using the existing
  data.  This solution is placed directly in ark_ycur.  This routine
  also computes the error estimate ||y-ytilde||_WRMS, where ytilde
  is the embedded solution, and the norm weights come from
  ark_ewt.  This norm value is returned.  The vector form of this
  estimated error (y-ytilde) is stored in ark_tempv1, in case the
  calling routine wishes to examine the error locations.

  Note: at this point in the step, the vector ark_tempv1 may be
  used as a temporary vector.
  ---------------------------------------------------------------*/
int erkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsm)
{
  /* local data */
  int retval, j, nvec;
  N_Vector y, yerr;
  realtype* cvals;
  N_Vector* Xvecs;
  ARKodeERKStepMem step_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "erkStep_ComputeSolutions", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  /* set N_Vector shortcuts */
  y    = ark_mem->ycur;
  yerr = ark_mem->tempv1;

  /* local shortcuts for fused vector operations */
  cvals = step_mem->cvals;
  Xvecs = step_mem->Xvecs;

  /* initialize output */
  *dsm = ZERO;


  /* Compute time step solution */
  /*   set arrays for fused vector operation */
  nvec = 0;
  for (j=0; j<step_mem->stages; j++) {
    cvals[nvec] = ark_mem->h * step_mem->B->b[j];
    Xvecs[nvec] = step_mem->F[j];
    nvec += 1;
  }
  cvals[nvec] = ONE;
  Xvecs[nvec] = ark_mem->yn;
  nvec += 1;

  /*   call fused vector operation to do the work */
  retval = N_VLinearCombination(nvec, cvals, Xvecs, y);
  if (retval != 0) return(ARK_VECTOROP_ERR);

  /* Compute yerr (if step adaptivity enabled) */
  if (!ark_mem->fixedstep) {

    /* set arrays for fused vector operation */
    nvec = 0;
    for (j=0; j<step_mem->stages; j++) {
      cvals[nvec] = ark_mem->h * (step_mem->B->b[j] - step_mem->B->d[j]);
      Xvecs[nvec] = step_mem->F[j];
      nvec += 1;
    }

    /* call fused vector operation to do the work */
    retval = N_VLinearCombination(nvec, cvals, Xvecs, yerr);
    if (retval != 0) return(ARK_VECTOROP_ERR);

    /* fill error norm */
    *dsm = N_VWrmsNorm(yerr, ark_mem->ewt);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  erkStep_DoErrorTest

  This routine performs the local error test for the ARK method.
  The weighted local error norm dsm is passed in, and
  the test dsm ?<= 1 is made.

  If the test passes, arkDoErrorTest returns ARK_SUCCESS.

  If the test fails, we revert to the last successful solution
  time, and:
    - if maxnef error test failures have occurred or if
      SUNRabs(h) = hmin, we return ARK_ERR_FAILURE.
    - otherwise: update time step factor eta based on local error
      estimate and reduce h.
  ---------------------------------------------------------------*/
int erkStep_DoErrorTest(ARKodeMem ark_mem, int *nefPtr, realtype dsm)
{
  realtype ehist2, hhist2;
  int retval;
  ARKodeHAdaptMem hadapt_mem;
  ARKodeERKStepMem step_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "erkStep_DoErrorTest", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

  if (step_mem->hadapt_mem == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep", "arkDoErrorTest",
                    MSG_ARKADAPT_NO_MEM);
    return(ARK_MEM_NULL);
  }
  hadapt_mem = step_mem->hadapt_mem;

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */
  if (dsm <= ONE) return(ARK_SUCCESS);

  /* Test failed; increment counters */
  (*nefPtr)++;
  step_mem->netf++;

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
  erkStep_PrepareNextStep

  This routine handles ARK-specific updates following a successful
  step: copying the ARK result to the current solution vector,
  updating the error/step history arrays, and setting the
  prospective step size, hprime, for the next step.  Along with
  hprime, it sets the ratio eta=hprime/h.  It also updates other
  state variables related to a change of step size.
  ---------------------------------------------------------------*/
int erkStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm)
{
  int retval;
  ARKodeERKStepMem step_mem;

  /* access ARKodeERKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ERKStep",
                    "erkStep_PrepareNextStep", MSG_ERKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeERKStepMem) ark_mem->step_mem;

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
