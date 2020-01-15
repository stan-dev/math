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
 * Implementation header file for ARKode's MRI time stepper module.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_MRISTEP_IMPL_H
#define _ARKODE_MRISTEP_IMPL_H

#include "arkode/arkode_mristep.h"
#include "arkode/arkode_arkstep.h"

#include "arkode_impl.h"
#include "arkode_arkstep_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  MRI function types
  ===============================================================*/

typedef int (*MRIStepSetInnerForcingFn)(void* inner_arkode_mem, realtype tshift,
                                        realtype tscale, N_Vector *f,
                                        int nvecs);

typedef int (*MRIStepInnerEvolveFn)(void* inner_arkode_mem, realtype t0,
                                    N_Vector y0, realtype tout);


/*===============================================================
  MRI time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  The type ARKodeMRIStepMem is type pointer to struct
  ARKodeMRIStepMemRec. This structure contains fields to
  perform a MRI time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeMRIStepMemRec {

  /* MRI problem specification */
  ARKRhsFn fs;   /* y' = fs(t,y) + ff(t,y)  */

  /* Outer RK method storage and parameters */
  N_Vector *F;            /* slow RHS at each stage */
  int q;                  /* method order           */
  int p;                  /* embedding order        */
  int stages;             /* number of stages       */
  ARKodeButcherTable B;   /* MRI Butcher table      */

  /* Inner stepper data */
  void           *inner_mem;         /* inner stepper memory            */
  N_Vector       *inner_forcing;     /* RHS forcing vectors             */
  int             inner_num_forcing; /* number of RHS forcing vectors   */
  int             inner_retval;      /* last inner stepper return value */
  MRISTEP_ID      inner_stepper_id;  /* inner stepper identifier        */

  /* Inner-stepper-supplied functions */
  MRIStepSetInnerForcingFn inner_setforcing; /* set inner forcing data  */
  MRIStepInnerEvolveFn     inner_evolve;     /* inner evolve function   */
  ARKTimestepFullRHSFn     inner_fullrhs;    /* inner full RHS function */

  /* User-supplied pre and post inner evolve functions */
  MRIStepPreInnerFn  pre_inner_evolve;
  MRIStepPostInnerFn post_inner_evolve;

  /* Counters */
  long int nfs;  /* num fs calls */

  /* Reusable arrays for fused vector operations */
  realtype* cvals;
  N_Vector* Xvecs;

} *ARKodeMRIStepMem;


/*===============================================================
  MRI time step module private function prototypes
  ===============================================================*/

/* Create MRIStep memory structure */
void* mriStep_Create(ARKRhsFn fs, realtype t0, N_Vector y0);

/* Interface routines supplied to ARKode */
int mriStep_Init(void* arkode_mem, int init_type);
int mriStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode);
int mriStep_TakeStep(void* arkode_mem, realtype *dsmPtr, int *nflagPtr);

/* Internal utility routines */
int mriStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeMRIStepMem *step_mem);
booleantype mriStep_CheckNVector(N_Vector tmpl);
int mriStep_SetButcherTable(ARKodeMem ark_mem);
int mriStep_CheckButcherTable(ARKodeMem ark_mem);

/* Attach ARKStep inner stepper */
int mriStep_AttachARK(void* arkode_mem, void* inner_arkode_mem);

/* Compute forcing for inner stepper */
int mriStep_ComputeInnerForcing(ARKodeMRIStepMem step_mem, int stage,
                                realtype cdiff);

/* Evolve ARKStep inner stepper */
int mriStep_EvolveInnerARK(void* inner_arkode_mem, realtype t0,
                           N_Vector y0, realtype tout);

/*===============================================================
  Reusable MRIStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_MRISTEP_NO_MEM    "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
