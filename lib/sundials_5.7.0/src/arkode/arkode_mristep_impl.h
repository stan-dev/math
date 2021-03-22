/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
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
#include "arkode_ls_impl.h"
#include "arkode_arkstep_impl.h"   /* need to improve API for user-supplied
                                      inner stepper, so that direct access to
                                      arkStep_SetInnerForcing and arkStep_FullRHS
                                      are no longer needed, then remove this line */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* Stage type identifiers */
#define MRISTAGE_ERK_FAST    0
#define MRISTAGE_ERK_NOFAST  1
#define MRISTAGE_DIRK_NOFAST 2
#define MRISTAGE_DIRK_FAST   3

/* Implicit solver constants (duplicate from arkode_arkstep_impl.h) */
#define MAXCOR    3              /* max number of nonlinear iterations */
#define CRDOWN    RCONST(0.3)    /* constant to estimate the convergence
                                    rate for the nonlinear equation */
#define DGMAX     RCONST(0.2)    /* if |gamma/gammap-1| > DGMAX then call lsetup */
#define RDIV      RCONST(2.3)    /* declare divergence if ratio del/delp > RDIV */
#define MSBP      20             /* max no. of steps between lsetup calls */
#define NLSCOEF   RCONST(0.1)


/*===============================================================
  MRI function types
  ===============================================================*/

typedef int (*MRIStepSetInnerForcingFn)(void* inner_arkode_mem, realtype tshift,
                                        realtype tscale, N_Vector *f,
                                        int nvecs);

typedef int (*MRIStepInnerEvolveFn)(void* inner_arkode_mem, realtype t0,
                                    N_Vector y0, realtype tout);

typedef int (*MRIStepInnerResetFn)(void* inner_arkode_mem, realtype tR,
                                   N_Vector yR);

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
  ARKRhsFn     fs;               /* y' = fs(t,y) + ff(t,y)         */
  booleantype  linear;           /* SUNTRUE if fs is linear        */
  booleantype  linear_timedep;   /* SUNTRUE if dfs/dy depends on t */
  booleantype  implicit;         /* SUNTRUE if fs is implicit      */

  /* Outer RK method storage and parameters */
  N_Vector *F;                   /* slow RHS at each stage */
  int q;                         /* method order           */
  int p;                         /* embedding order        */
  int stages;                    /* number of stages       */
  MRIStepCoupling MRIC;          /* slow->fast coupling    */
  int *stagetypes;               /* type flags for stages  */
  realtype *rkcoeffs;            /* equivalent RK coeffs   */

  /* Algebraic solver data and parameters */
  N_Vector           sdata;      /* old stage data in residual               */
  N_Vector           zpred;      /* predicted stage solution                 */
  N_Vector           zcor;       /* stage correction                         */
  int                istage;     /* current stage index                      */
  SUNNonlinearSolver NLS;        /* generic SUNNonlinearSolver object        */
  booleantype        ownNLS;     /* flag indicating ownership of NLS         */
  realtype           gamma;      /* gamma = h * A(i,i)                       */
  realtype           gammap;     /* gamma at the last setup call             */
  realtype           gamrat;     /* gamma / gammap                           */
  realtype           dgmax;      /* call lsetup if |gamma/gammap-1| >= dgmax */
  int                predictor;  /* implicit prediction method to use        */
  realtype           crdown;     /* nonlinear conv rate estimation constant  */
  realtype           rdiv;       /* nonlin divergence if del/delp > rdiv     */
  realtype           crate;      /* estimated nonlin convergence rate        */
  realtype           delp;       /* norm of previous nonlinear solver update */
  realtype           eRNrm;      /* estimated residual norm, used in nonlin
                                    and linear solver convergence tests      */
  realtype           nlscoef;    /* coefficient in nonlin. convergence test  */

  int                msbp;       /* positive => max # steps between lsetup
                                    negative => call at each Newton iter     */
  long int           nstlp;      /* step number of last setup call           */

  int                maxcor;     /* max num iterations for solving the
                                    nonlinear equation                       */
  int                convfail;   /* NLS fail flag (for interface routines)   */
  booleantype        jcur;       /* is Jacobian info for lin solver current? */
  ARKStagePredictFn  stage_predict;  /* User-supplied stage predictor        */

  /* Linear Solver Data */
  ARKLinsolInitFn    linit;
  ARKLinsolSetupFn   lsetup;
  ARKLinsolSolveFn   lsolve;
  ARKLinsolFreeFn    lfree;
  void              *lmem;

  /* Inner stepper data */
  void              *inner_mem;         /* inner stepper memory            */
  N_Vector          *inner_forcing;     /* RHS forcing vectors             */
  int                inner_num_forcing; /* number of RHS forcing vectors   */
  int                inner_retval;      /* last inner stepper return value */
  MRISTEP_ID         inner_stepper_id;  /* inner stepper identifier        */

  /* Inner-stepper-supplied functions */
  MRIStepSetInnerForcingFn inner_setforcing; /* set inner forcing data  */
  MRIStepInnerEvolveFn     inner_evolve;     /* inner evolve function   */
  ARKTimestepFullRHSFn     inner_fullrhs;    /* inner full RHS function */
  MRIStepInnerResetFn      inner_reset;      /* inner reset function    */

  /* User-supplied pre and post inner evolve functions */
  MRIStepPreInnerFn  pre_inner_evolve;
  MRIStepPostInnerFn post_inner_evolve;

  /* Counters */
  long int nfs;       /* num fs calls                  */
  long int nsetups;   /* num linear solver setup calls */
  long int nls_iters; /* num nonlinear solver iters    */

  /* Reusable arrays for fused vector operations */
  realtype* cvals;
  N_Vector* Xvecs;

} *ARKodeMRIStepMem;


/*===============================================================
  MRI time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKode */
int mriStep_AttachLinsol(void* arkode_mem, ARKLinsolInitFn linit,
                         ARKLinsolSetupFn lsetup,
                         ARKLinsolSolveFn lsolve,
                         ARKLinsolFreeFn lfree,
                         SUNLinearSolver_Type lsolve_type,
                         void *lmem);
void mriStep_DisableLSetup(void* arkode_mem);
int mriStep_Init(void* arkode_mem, int init_type);
void* mriStep_GetLmem(void* arkode_mem);
ARKRhsFn mriStep_GetImplicitRHS(void* arkode_mem);
int mriStep_GetGammas(void* arkode_mem, realtype *gamma,
                      realtype *gamrat, booleantype **jcur,
                      booleantype *dgamma_fail);
int mriStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode);
int mriStep_TakeStep(void* arkode_mem, realtype *dsmPtr, int *nflagPtr);

/* Internal utility routines */
int mriStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeMRIStepMem *step_mem);
booleantype mriStep_CheckNVector(N_Vector tmpl);
int mriStep_SetCoupling(ARKodeMem ark_mem);
int mriStep_CheckCoupling(ARKodeMem ark_mem);
int mriStep_StageERKFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                         int is);
int mriStep_StageERKNoFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                           int is);
int mriStep_StageDIRKFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                          int is, int *nflagPtr);
int mriStep_StageDIRKNoFast(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                            int is, int *nflagPtr);
int mriStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess);
int mriStep_StageSetup(ARKodeMem ark_mem);
int mriStep_NlsInit(ARKodeMem ark_mem);
int mriStep_Nls(ARKodeMem ark_mem, int nflag);

/* private functions passed to nonlinear solver */
int mriStep_NlsResidual(N_Vector yy, N_Vector res, void* arkode_mem);
int mriStep_NlsFPFunction(N_Vector yy, N_Vector res, void* arkode_mem);
int mriStep_NlsLSetup(booleantype jbad, booleantype* jcur, void* arkode_mem);
int mriStep_NlsLSolve(N_Vector delta, void* arkode_mem);
int mriStep_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                        realtype tol, N_Vector ewt, void* arkode_mem);


/* Attach ARKStep inner stepper */
int mriStep_AttachARK(void* arkode_mem, void* inner_arkode_mem);

/* Compute forcing for inner stepper */
int mriStep_ComputeInnerForcing(ARKodeMRIStepMem step_mem, int stage,
                                realtype cdiff);

/* Evolve ARKStep inner stepper */
int mriStep_EvolveInnerARK(void* inner_arkode_mem, realtype t0,
                           N_Vector y0, realtype tout);

/* Return stage type (implicit/explicit + fast/nofast) */
int mriStep_StageType(MRIStepCoupling MRIC, int is);

/* Return effective RK coefficients (nofast stage) */
int mriStep_RKCoeffs(MRIStepCoupling MRIC, int is, realtype *Arow);

/*===============================================================
  Reusable MRIStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_MRISTEP_NO_MEM    "Time step module memory is NULL."
#define MSG_NLS_INIT_FAIL     "The nonlinear solver's init routine failed."

#ifdef __cplusplus
}
#endif

#endif
