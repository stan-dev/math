/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Implementation header file for ARKode's ERK time stepper
 * module.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_ERKSTEP_IMPL_H
#define _ARKODE_ERKSTEP_IMPL_H

#include <arkode/arkode_erkstep.h>
#include "arkode_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ERK time step module constants -- move many items here from
  arkode_impl.h
  ===============================================================*/



/*===============================================================
  ERK time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeERKStepMemRec, ARKodeERKStepMem
  ---------------------------------------------------------------
  The type ARKodeERKStepMem is type pointer to struct
  ARKodeERKStepMemRec.  This structure contains fields to
  perform an explicit Runge-Kutta time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeERKStepMemRec {

  /* ERK problem specification */
  ARKRhsFn f;             /* y' = f(t,y)                */

  /* ARK method storage and parameters */
  N_Vector *F;            /* explicit RHS at each stage */
  int q;                  /* method order               */
  int p;                  /* embedding order            */
  int stages;             /* number of stages           */
  ARKodeButcherTable B;   /* ERK Butcher table          */

  /* Time step adaptivity data */
  ARKodeHAdaptMem hadapt_mem;  /* time step adaptivity structure   */
  booleantype     hadapt_pq;   /* choice of using p (0) vs q (1)   */
  int             maxnef;      /* max error test fails in one step */

  /* Counters */
  long int nst_attempts;  /* num attempted steps                */
  long int nfe;           /* num fe calls                       */
  long int netf;          /* num error test failures            */

  /* Reusable arrays for fused vector operations */
  realtype* cvals;
  N_Vector* Xvecs;

} *ARKodeERKStepMem;


/*===============================================================
  ERK time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKode */
int erkStep_Init(void* arkode_mem, int init_type);
int erkStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode);
int erkStep_TakeStep(void* arkode_mem);

/* Internal utility routines */
int erkStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeERKStepMem *step_mem);
booleantype erkStep_CheckNVector(N_Vector tmpl);
int erkStep_SetButcherTable(ARKodeMem ark_mem);
int erkStep_CheckButcherTable(ARKodeMem ark_mem);

int erkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsm);
int erkStep_DoErrorTest(ARKodeMem ark_mem, int *nefPtr,
                        realtype dsm);
int erkStep_PrepareNextStep(ARKodeMem ark_mem, realtype dsm);

/*===============================================================
  Reusable ERKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_ERKSTEP_NO_MEM    "Time step module memory is NULL."

#ifdef __cplusplus
}
#endif

#endif
