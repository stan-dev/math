/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for the optional input and output functions
 * for the ARKode MRIStep time stepper module.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_mristep_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM "Lg"
#else
#define RSYM "g"
#endif


/*===============================================================
  MRIStep Optional input functions (wrappers for generic ARKode
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int MRIStepSetDenseOrder(void *arkode_mem, int dord) {
  return(MRIStepSetInterpolantDegree(arkode_mem, dord)); }
int MRIStepSetInterpolantDegree(void *arkode_mem, int degree) {
  if (degree < 0) degree = ARK_INTERP_MAX_DEGREE;
  return(arkSetInterpolantDegree(arkode_mem, degree)); }
int MRIStepSetInterpolantType(void *arkode_mem, int itype) {
  return(arkSetInterpolantType(arkode_mem, itype)); }
int MRIStepSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun,
                           void *eh_data) {
  return(arkSetErrHandlerFn(arkode_mem, ehfun, eh_data)); }
int MRIStepSetErrFile(void *arkode_mem, FILE *errfp) {
  return(arkSetErrFile(arkode_mem, errfp)); }
int MRIStepSetDiagnostics(void *arkode_mem, FILE *diagfp) {
  return(arkSetDiagnostics(arkode_mem, diagfp)); }
int MRIStepSetMaxNumSteps(void *arkode_mem, long int mxsteps) {
  return(arkSetMaxNumSteps(arkode_mem, mxsteps)); }
int MRIStepSetMaxHnilWarns(void *arkode_mem, int mxhnil) {
  return(arkSetMaxHnilWarns(arkode_mem, mxhnil)); }
int MRIStepSetStopTime(void *arkode_mem, realtype tstop) {
  return(arkSetStopTime(arkode_mem, tstop)); }
int MRIStepSetRootDirection(void *arkode_mem, int *rootdir) {
  return(arkSetRootDirection(arkode_mem, rootdir)); }
int MRIStepSetNoInactiveRootWarn(void *arkode_mem) {
  return(arkSetNoInactiveRootWarn(arkode_mem)); }
int MRIStepSetPostprocessStepFn(void *arkode_mem,
                                ARKPostProcessFn ProcessStep) {
  return(arkSetPostprocessStepFn(arkode_mem, ProcessStep)); }
int MRIStepSetPostprocessStageFn(void *arkode_mem,
                                 ARKPostProcessFn ProcessStage) {
  return(arkSetPostprocessStageFn(arkode_mem, ProcessStage)); }


/*===============================================================
  MRIStep Optional input functions (customized wrappers for
  generic ARKode utility routines).  All are documented in
  arkode_io.c and arkode_ls.c.
  ===============================================================*/

int MRIStepSetFixedStep(void *arkode_mem, realtype hsfixed)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetFixedStep", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (hsfixed == ZERO) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepSetFixedStep",
                    "MIRStep does not support adaptive steps at this time.");
    return(ARK_ILL_INPUT);
  }

  /* call generic routine for remaining work */
  return(arkSetFixedStep(ark_mem, hsfixed));
}


/*===============================================================
  MRIStep Optional output functions (wrappers for generic ARKode
  utility routines).  All are documented in arkode_io.c.
  ===============================================================*/
int MRIStepGetNumSteps(void *arkode_mem, long int *nssteps) {
  return(arkGetNumSteps(arkode_mem, nssteps)); }
int MRIStepGetLastStep(void *arkode_mem, realtype *hlast) {
  return(arkGetLastStep(arkode_mem, hlast)); }
int MRIStepGetCurrentTime(void *arkode_mem, realtype *tcur) {
  return(arkGetCurrentTime(arkode_mem, tcur)); }
int MRIStepGetCurrentState(void *arkode_mem, N_Vector *ycur) {
  return(arkGetCurrentState(arkode_mem, ycur)); }
int MRIStepGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw) {
  return(arkGetWorkSpace(arkode_mem, lenrw, leniw)); }
int MRIStepGetNumGEvals(void *arkode_mem, long int *ngevals) {
  return(arkGetNumGEvals(arkode_mem, ngevals)); }
int MRIStepGetRootInfo(void *arkode_mem, int *rootsfound) {
  return(arkGetRootInfo(arkode_mem, rootsfound)); }
char *MRIStepGetReturnFlagName(long int flag) {
  return(arkGetReturnFlagName(flag)); }


/*===============================================================
  MRIStep optional input functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepSetDefaults:

  Resets all MRIStep optional inputs to their default values.
  Does not change problem-defining function pointers or
  user_data pointer.
  ---------------------------------------------------------------*/
int MRIStepSetDefaults(void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set default values for integrator optional inputs */
  step_mem->q         = 3;              /* method order */
  step_mem->p         = 0;              /* embedding order */
  step_mem->stages    = 0;              /* no stages */
  step_mem->B         = NULL;           /* no Butcher table */

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetUserData: Specifies the user data pointer
  ---------------------------------------------------------------*/
int MRIStepSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetUserData", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  return(arkSetUserData(ark_mem, user_data));
}


/*---------------------------------------------------------------
  MRIStepSetTable:

  Specifies to use a customized Butcher table for the explicit
  portion of the system.
  ---------------------------------------------------------------*/
int MRIStepSetTable(void *arkode_mem, int q, ARKodeButcherTable B)
{
  int retval;
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetTable",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check for illegal inputs */
  if (B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetTables", MSG_ARK_NO_MEM);
    return(ARK_ILL_INPUT);
  }

  /* clear any existing parameters and Butcher tables */
  step_mem->stages = 0;
  step_mem->q = 0;
  step_mem->p = 0;
  ARKodeButcherTable_Free(step_mem->B);
  step_mem->B = NULL;

  /* set the relevant parameters */
  step_mem->stages = B->stages;
  step_mem->q = B->q;
  step_mem->p = 0; /* assume fixed stepping */

  /* copy the table in step memory */
  step_mem->B = ARKodeButcherTable_Copy(B);
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetTables", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepSetTableNum:

  Specifies to use pre-existing Butcher tables for the problem,
  based on the integer flags passed to ARKodeButcherTable_LoadERK()
  within the file arkode_butcher_erk.c.
  ---------------------------------------------------------------*/
int MRIStepSetTableNum(void *arkode_mem, int itable)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetTableNum",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check that argument specifies an explicit table (assume explicit) */
  if (itable < MIN_ERK_NUM || itable > MAX_ERK_NUM ) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetTableNum",
                    "Illegal MRI table number");
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
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepSetTableNum",
                    "Error setting table with that index");
    return(ARK_ILL_INPUT);
  }
  step_mem->stages = step_mem->B->stages;
  step_mem->q = step_mem->B->q;
  step_mem->p = step_mem->B->p;

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepSetPreInnerFn:

  Sets the user-supplied function called BEFORE the inner evolve
  ---------------------------------------------------------------*/
int MRIStepSetPreInnerFn(void *arkode_mem, MRIStepPreInnerFn prefn)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set pre inner evolve function */
  step_mem->pre_inner_evolve = prefn;

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
  MRIStepSetPostInnerFn:

  Sets the user-supplied function called AFTER the inner evolve
  ---------------------------------------------------------------*/
int MRIStepSetPostInnerFn(void *arkode_mem, MRIStepPostInnerFn postfn)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetDefaults",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* Set pre inner evolve function */
  step_mem->post_inner_evolve = postfn;

  return(ARK_SUCCESS);
}


/*===============================================================
  MRIStep optional output functions -- stepper-specific
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepGetLastInnerStepFlag:

  Returns the last return value from the inner stepper.
  ---------------------------------------------------------------*/
int MRIStepGetLastInnerStepFlag(void *arkode_mem, int *flag)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetLastInnerStepFlag",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get the last return value from the inner stepper */
  *flag = step_mem->inner_retval;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetNumRhsEvals:

  Returns the current number of calls to fs and ff
  ---------------------------------------------------------------*/
int MRIStepGetNumRhsEvals(void *arkode_mem, long int *nfs_evals)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetNumRhsEvals",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get number of fs evals from step_mem */
  *nfs_evals = step_mem->nfs;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetCurrentButcherTables:

  Sets pointers to the slow and fast Butcher tables currently in
  use.
  ---------------------------------------------------------------*/
int MRIStepGetCurrentButcherTables(void *arkode_mem, ARKodeButcherTable *B)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetCurrentButcherTable",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* get table from step_mem */
  *B = step_mem->B;

  return(ARK_SUCCESS);
}


/*===============================================================
  MRIStep parameter output
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepWriteParameters:

  Outputs all solver parameters to the provided file pointer.
  ---------------------------------------------------------------*/
int MRIStepWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepWriteParameters",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* output ARKode infrastructure parameters first */
  retval = arkWriteParameters(arkode_mem, fp);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepWriteParameters",
                    "Error writing ARKode infrastructure parameters");
    return(retval);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepWriteButcher:

  Outputs Butcher tables to the provided file pointer.
  ---------------------------------------------------------------*/
int MRIStepWriteButcher(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepWriteButcher",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS) return(retval);

  /* check that Butcher table is non-NULL (otherwise report error) */
  if (step_mem->B == NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "MRIStepWriteButcher", "Butcher table memory is NULL");
    return(ARK_MEM_NULL);
  }

  /* wrie outer Butcher table */
  fprintf(fp, "\nMRIStep Butcher tables:\n");
  if (step_mem->B != NULL) {
    fprintf(fp, "  Slow Butcher table (stages = %i):\n", step_mem->stages);
    ARKodeButcherTable_Write(step_mem->B, fp);
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
