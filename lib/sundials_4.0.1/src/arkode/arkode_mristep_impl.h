/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------------------
 * Implementation header file for ARKode's MRI time stepper module.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_MRISTEP_IMPL_H
#define _ARKODE_MRISTEP_IMPL_H

#include <arkode/arkode_mristep.h>
#include "arkode_impl.h"
#include "arkode/arkode_arkstep.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  MRI time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeMRIStepMemRec, ARKodeMRIStepMem
  ---------------------------------------------------------------
  The type ARKodeMRIStepMem is type pointer to struct
  ARKodeMRIStepMemRec. This structure contains fields to
  perform a MRI time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeMRIStepMemRec {

  /* MRI problem specification */
  ARKRhsFn fs;            /* y' = fs(t,y) + ff(t,y)     */
  ARKRhsFn ff;

  /* Outer RK method storage and parameters */
  N_Vector *F;            /* slow RHS at each stage */
  int q;                  /* method order           */
  int p;                  /* embedding order        */
  int stages;             /* number of stages       */
  ARKodeButcherTable B;   /* MRI Butcher table      */

  /* Inner stepper data */
  void     *inner_arkode_mem;  /* inner stepper memory            */
  N_Vector forcing;            /* RHS forcing vector              */
  realtype hf;                 /* inner step size                 */
  int      inner_retval;       /* last inner stepper return value */

  /* Counters */
  long int nfs;  /* num fe calls */
  long int nff;  /* num fe calls */

  /* Reusable arrays for fused vector operations */
  realtype* cvals;
  N_Vector* Xvecs;

} *ARKodeMRIStepMem;



/*===============================================================
  MRI time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKode */
int mriStep_Init(void* arkode_mem, int init_type);
int mriStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode);
int mriStep_TakeStep(void* arkode_mem);

/* Internal utility routines */
int mriStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeMRIStepMem *step_mem);
booleantype mriStep_CheckNVector(N_Vector tmpl);
int mriStep_SetButcherTable(ARKodeMem ark_mem);
int mriStep_CheckButcherTable(ARKodeMem ark_mem);

int mriStep_ComputeErrorEst(ARKodeMem ark_mem, realtype *dsm);
int mriStep_DoErrorTest(ARKodeMem ark_mem, int *nefPtr,
                        realtype dsm);
int mriStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm);

/* Internal inner stepper routines */
int mriStep_InnerRhsFn(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/*===============================================================
  Reusable MRIStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_MRISTEP_NO_MEM    "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
