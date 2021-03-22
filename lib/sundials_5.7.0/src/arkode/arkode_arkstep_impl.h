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
 * Implementation header file for ARKode's ARK time stepper
 * module.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_ARKSTEP_IMPL_H
#define _ARKODE_ARKSTEP_IMPL_H

#include <arkode/arkode_arkstep.h>
#include "arkode_impl.h"
#include "arkode_ls_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARK time step module constants
  ===============================================================*/

#define MAXCOR    3              /* max number of nonlinear iterations */
#define CRDOWN    RCONST(0.3)    /* constant to estimate the convergence
                                    rate for the nonlinear equation */
#define DGMAX     RCONST(0.2)    /* if |gamma/gammap-1| > DGMAX then call lsetup */
#define RDIV      RCONST(2.3)    /* declare divergence if ratio del/delp > RDIV */
#define MSBP      20             /* max no. of steps between lsetup calls */

/* Default solver tolerance factor */
/* #define NLSCOEF   RCONST(0.003) */  /* Hairer & Wanner constant */
/* #define NLSCOEF   RCONST(0.2)   */  /* CVODE constant */
#define NLSCOEF   RCONST(0.1)

/* Mass matrix types */
#define MASS_IDENTITY 0
#define MASS_FIXED    1
#define MASS_TIMEDEP  2


/*===============================================================
  ARK time step module data structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeARKStepMemRec, ARKodeARKStepMem
  ---------------------------------------------------------------
  The type ARKodeARKStepMem is type pointer to struct
  ARKodeARKStepMemRec.  This structure contains fields to
  perform an additive Runge-Kutta time step.
  ---------------------------------------------------------------*/
typedef struct ARKodeARKStepMemRec {

  /* ARK problem specification */
  ARKRhsFn     fe;              /* My' = fe(t,y) + fi(t,y)        */
  ARKRhsFn     fi;
  booleantype  linear;          /* SUNTRUE if fi is linear        */
  booleantype  linear_timedep;  /* SUNTRUE if dfi/dy depends on t */
  booleantype  explicit;        /* SUNTRUE if fe is enabled       */
  booleantype  implicit;        /* SUNTRUE if fi is enabled       */

  /* ARK method storage and parameters */
  N_Vector *Fe;           /* explicit RHS at each stage */
  N_Vector *Fi;           /* implicit RHS at each stage */
  N_Vector  sdata;        /* old stage data in residual */
  N_Vector  zpred;        /* predicted stage solution   */
  N_Vector  zcor;         /* stage correction           */
  int q;                  /* method order               */
  int p;                  /* embedding order            */
  int istage;             /* current stage              */
  int stages;             /* number of stages           */
  ARKodeButcherTable Be;  /* ERK Butcher table          */
  ARKodeButcherTable Bi;  /* IRK Butcher table          */

  /* User-supplied stage predictor routine */
  ARKStagePredictFn stage_predict;

  /* (Non)Linear solver parameters & data */
  SUNNonlinearSolver NLS;   /* generic SUNNonlinearSolver object     */
  booleantype     ownNLS;   /* flag indicating ownership of NLS      */
  realtype gamma;        /* gamma = h * A(i,i)                       */
  realtype gammap;       /* gamma at the last setup call             */
  realtype gamrat;       /* gamma / gammap                           */
  realtype dgmax;        /* call lsetup if |gamma/gammap-1| >= dgmax */

  int      predictor;    /* implicit prediction method to use        */
  realtype crdown;       /* nonlinear conv rate estimation constant  */
  realtype rdiv;         /* nonlin divergence if del/delp > rdiv     */
  realtype crate;        /* estimated nonlin convergence rate        */
  realtype delp;         /* norm of previous nonlinear solver update */
  realtype eRNrm;        /* estimated residual norm, used in nonlin
                            and linear solver convergence tests      */
  realtype nlscoef;      /* coefficient in nonlin. convergence test  */

  int      msbp;         /* positive => max # steps between lsetup
                            negative => call at each Newton iter     */
  long int nstlp;        /* step number of last setup call           */

  int      maxcor;       /* max num iterations for solving the
                            nonlinear equation                       */

  int      convfail;     /* NLS fail flag (for interface routines)   */
  booleantype jcur;      /* is Jacobian info for lin solver current? */

  /* Linear Solver Data */
  ARKLinsolInitFn      linit;
  ARKLinsolSetupFn     lsetup;
  ARKLinsolSolveFn     lsolve;
  ARKLinsolFreeFn      lfree;
  void                *lmem;
  SUNLinearSolver_Type lsolve_type;

  /* Mass matrix solver data */
  ARKMassInitFn        minit;
  ARKMassSetupFn       msetup;
  ARKMassMultFn        mmult;
  ARKMassSolveFn       msolve;
  ARKMassFreeFn        mfree;
  void*                mass_mem;
  int                  mass_type;  /* 0=identity, 1=fixed, 2=time-dep */
  SUNLinearSolver_Type msolve_type;

  /* Counters */
  long int nfe;       /* num fe calls               */
  long int nfi;       /* num fi calls               */
  long int nsetups;   /* num setup calls            */
  long int nls_iters; /* num nonlinear solver iters */

  /* Reusable arrays for fused vector operations */
  realtype *cvals;         /* scalar array for fused ops       */
  N_Vector *Xvecs;         /* array of vectors for fused ops   */
  int       nfusedopvecs;  /* length of cvals and Xvecs arrays */

  /* Data for using ARKStep with external polynomial forcing */
  booleantype expforcing;  /* add forcing to explicit RHS */
  booleantype impforcing;  /* add forcing to implicit RHS */
  realtype    tshift;      /* time normalization shift    */
  realtype    tscale;      /* time normalization scaling  */
  N_Vector*   forcing;     /* array of forcing vectors    */
  int         nforcing;    /* number of forcing vectors   */

} *ARKodeARKStepMem;


/*===============================================================
  ARK time step module private function prototypes
  ===============================================================*/

/* Interface routines supplied to ARKode */
int arkStep_AttachLinsol(void* arkode_mem, ARKLinsolInitFn linit,
                         ARKLinsolSetupFn lsetup,
                         ARKLinsolSolveFn lsolve,
                         ARKLinsolFreeFn lfree,
                         SUNLinearSolver_Type lsolve_type,
                         void *lmem);
int arkStep_AttachMasssol(void* arkode_mem,
                          ARKMassInitFn minit,
                          ARKMassSetupFn msetup,
                          ARKMassMultFn mmult,
                          ARKMassSolveFn msolve,
                          ARKMassFreeFn lfree,
                          booleantype time_dep,
                          SUNLinearSolver_Type msolve_type,
                          void *mass_mem);
void arkStep_DisableLSetup(void* arkode_mem);
void arkStep_DisableMSetup(void* arkode_mem);
int arkStep_Init(void* arkode_mem, int init_type);
void* arkStep_GetLmem(void* arkode_mem);
void* arkStep_GetMassMem(void* arkode_mem);
ARKRhsFn arkStep_GetImplicitRHS(void* arkode_mem);
int arkStep_GetGammas(void* arkode_mem, realtype *gamma,
                      realtype *gamrat, booleantype **jcur,
                      booleantype *dgamma_fail);
int arkStep_FullRHS(void* arkode_mem, realtype t,
                    N_Vector y, N_Vector f, int mode);
int arkStep_TakeStep_Z(void* arkode_mem, realtype *dsmPtr, int *nflagPtr);

/* Internal utility routines */
int arkStep_AccessStepMem(void* arkode_mem, const char *fname,
                          ARKodeMem *ark_mem, ARKodeARKStepMem *step_mem);
booleantype arkStep_CheckNVector(N_Vector tmpl);
int arkStep_SetButcherTables(ARKodeMem ark_mem);
int arkStep_CheckButcherTables(ARKodeMem ark_mem);
int arkStep_Predict(ARKodeMem ark_mem, int istage, N_Vector yguess);
int arkStep_StageSetup(ARKodeMem ark_mem, booleantype implicit);
int arkStep_NlsInit(ARKodeMem ark_mem);
int arkStep_Nls(ARKodeMem ark_mem, int nflag);
int arkStep_ComputeSolutions(ARKodeMem ark_mem, realtype *dsm);
int arkStep_ComputeSolutions_MassFixed(ARKodeMem ark_mem, realtype *dsm);
void arkStep_ApplyForcing(ARKodeARKStepMem step_mem, realtype t,
                          realtype s, int *nvec);

/* private functions passed to nonlinear solver */
int arkStep_NlsResidual_MassIdent(N_Vector zcor, N_Vector r, void* arkode_mem);
int arkStep_NlsResidual_MassFixed(N_Vector zcor, N_Vector r, void* arkode_mem);
int arkStep_NlsResidual_MassTDep(N_Vector zcor, N_Vector r, void* arkode_mem);
int arkStep_NlsFPFunction_MassIdent(N_Vector zcor, N_Vector g, void* arkode_mem);
int arkStep_NlsFPFunction_MassFixed(N_Vector zcor, N_Vector g, void* arkode_mem);
int arkStep_NlsFPFunction_MassTDep(N_Vector zcor, N_Vector g, void* arkode_mem);
int arkStep_NlsLSetup(booleantype jbad, booleantype* jcur, void* arkode_mem);
int arkStep_NlsLSolve(N_Vector delta, void* arkode_mem);
int arkStep_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                        realtype tol, N_Vector ewt, void* arkode_mem);

/* private functions used by MRIStep */
int arkStep_SetInnerForcing(void* arkode_mem, realtype tshift, realtype tscale,
                            N_Vector *f, int nvecs);

/*===============================================================
  Reusable ARKStep Error Messages
  ===============================================================*/

/* Initialization and I/O error messages */
#define MSG_ARKSTEP_NO_MEM    "Time step module memory is NULL."
#define MSG_NLS_INIT_FAIL     "The nonlinear solver's init routine failed."

/* Other error messages */
#define MSG_ARK_MISSING_FE     "Cannot specify that method is explicit without providing a function pointer to fe(t,y)."
#define MSG_ARK_MISSING_FI     "Cannot specify that method is implicit without providing a function pointer to fi(t,y)."
#define MSG_ARK_MISSING_F      "Cannot specify that method is ImEx without providing function pointers to fi(t,y) and fe(t,y)."

#ifdef __cplusplus
}
#endif

#endif
