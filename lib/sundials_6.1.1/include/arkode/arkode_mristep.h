/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the ARKode MRIStep module.
 * -----------------------------------------------------------------*/

#ifndef _MRISTEP_H
#define _MRISTEP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <arkode/arkode.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_butcher_dirk.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * MRIStep Constants
 * ----------------- */

/* MRIStep method types */
typedef enum {
  MRISTEP_EXPLICIT,
  MRISTEP_IMPLICIT,
  MRISTEP_IMEX
} MRISTEP_METHOD_TYPE;


/* MRI coupling table table accessor IDs:
     ERK:    0 -  99
     DIRK: 100 - 199
     MRI:  200 - 299          */

/* DEPRECATED MIS_KW3: use ARKODE_MIS_KW3 */
#define MIS_KW3             200
/* DEPRECATED MRI_GARK_ERK33a: use ARKODE_MRI_GARK_ERK33a */
#define MRI_GARK_ERK33a     201
/* DEPRECATED MRI_GARK_ERK45a: use ARKODE_MRI_GARK_ERK45a */
#define MRI_GARK_ERK45a     202
/* DEPRECATED MRI_GARK_IRK21a: use ARKODE_MRI_GARK_IRK21a */
#define MRI_GARK_IRK21a     203
/* DEPRECATED MRI_GARK_ESDIRK34a: use ARKODE_MRI_GARK_ESDIRK34a */
#define MRI_GARK_ESDIRK34a  204
/* DEPRECATED MRI_GARK_ESDIRK46a: use ARKODE_MRI_GARK_ESDIRK46a */
#define MRI_GARK_ESDIRK46a  205
/* DEPRECATED IMEX_MRI_GARK3a: use ARKODE_IMEX_MRI_GARK3a */
#define IMEX_MRI_GARK3a     206
/* DEPRECATED IMEX_MRI_GARK3b: use ARKODE_IMEX_MRI_GARK3b */
#define IMEX_MRI_GARK3b     207
/* DEPRECATED IMEX_MRI_GARK4: use ARKODE_IMEX_MRI_GARK4 */
#define IMEX_MRI_GARK4      208

/* Utility #defines to ensure valid input IDs for MRI tables */

/* DEPRECATED MIN_MRI_NUM: use ARKODE_MIN_MRI_NUM */
#define MIN_MRI_NUM         200
/* DEPRECATED MAX_MRI_NUM: use ARKODE_MAX_MRI_NUM */
#define MAX_MRI_NUM         208

typedef enum {
  ARKODE_MRI_NONE = -1, /* ensure enum is signed int */
  ARKODE_MIN_MRI_NUM = 200,
  ARKODE_MIS_KW3 = ARKODE_MIN_MRI_NUM,
  ARKODE_MRI_GARK_ERK33a,
  ARKODE_MRI_GARK_ERK45a,
  ARKODE_MRI_GARK_IRK21a,
  ARKODE_MRI_GARK_ESDIRK34a,
  ARKODE_MRI_GARK_ESDIRK46a,
  ARKODE_IMEX_MRI_GARK3a,
  ARKODE_IMEX_MRI_GARK3b,
  ARKODE_IMEX_MRI_GARK4,
  ARKODE_MAX_MRI_NUM = ARKODE_IMEX_MRI_GARK4
} ARKODE_MRITableID;


/* Default MRI coupling tables for each order */

static const int MRISTEP_DEFAULT_3         = ARKODE_MIS_KW3;
static const int MRISTEP_DEFAULT_EXPL_3    = ARKODE_MIS_KW3;
static const int MRISTEP_DEFAULT_EXPL_4    = ARKODE_MRI_GARK_ERK45a;
static const int MRISTEP_DEFAULT_IMPL_SD_2 = ARKODE_MRI_GARK_IRK21a;
static const int MRISTEP_DEFAULT_IMPL_SD_3 = ARKODE_MRI_GARK_ESDIRK34a;
static const int MRISTEP_DEFAULT_IMPL_SD_4 = ARKODE_MRI_GARK_ESDIRK46a;
static const int MRISTEP_DEFAULT_IMEX_SD_3 = ARKODE_IMEX_MRI_GARK3b;
static const int MRISTEP_DEFAULT_IMEX_SD_4 = ARKODE_IMEX_MRI_GARK4;

/* DEPRECATED DEFAULT_MRI_TABLE_3: use MRISTEP_DEFAULT_3 */
#define DEFAULT_MRI_TABLE_3          MRISTEP_DEFAULT_3     /* backwards-compatibility */
/* DEPRECATED DEFAULT_EXPL_MRI_TABLE_3: use MRISTEP_DEFAULT_EXPL_3 */
#define DEFAULT_EXPL_MRI_TABLE_3     MRISTEP_DEFAULT_EXPL_3
/* DEPRECATED DEFAULT_EXPL_MRI_TABLE_4: use MRISTEP_DEFAULT_EXPL_4 */
#define DEFAULT_EXPL_MRI_TABLE_4     MRISTEP_DEFAULT_EXPL_4
/* DEPRECATED DEFAULT_IMPL_SD_TABLE_2: use MRISTEP_DEFAULT_IMPL_SD_2 */
#define DEFAULT_IMPL_SD_MRI_TABLE_2  MRISTEP_DEFAULT_IMPL_SD_2
/* DEPRECATED DEFAULT_IMPL_SD_TABLE_3: use MRISTEP_DEFAULT_IMPL_SD_3 */
#define DEFAULT_IMPL_SD_MRI_TABLE_3  MRISTEP_DEFAULT_IMPL_SD_3
/* DEPRECATED DEFAULT_IMPL_SD_TABLE_4: use MRISTEP_DEFAULT_IMPL_SD_4 */
#define DEFAULT_IMPL_SD_MRI_TABLE_4  MRISTEP_DEFAULT_IMPL_SD_4
/* DEPRECATED DEFAULT_IMEX_SD_TABLE_3: use MRISTEP_DEFAULT_IMEX_SD_3 */
#define DEFAULT_IMEX_SD_MRI_TABLE_3  MRISTEP_DEFAULT_IMEX_SD_3
/* DEPRECATED DEFAULT_IMEX_SD_TABLE_4: use MRISTEP_DEFAULT_IMEX_SD_4 */
#define DEFAULT_IMEX_SD_MRI_TABLE_4  MRISTEP_DEFAULT_IMEX_SD_4

/* ------------------------------------
 * MRIStep Inner Stepper Function Types
 * ------------------------------------ */

typedef int (*MRIStepInnerEvolveFn)(MRIStepInnerStepper stepper,
                                    realtype t0, realtype tout, N_Vector y);

typedef int (*MRIStepInnerFullRhsFn)(MRIStepInnerStepper stepper,
                                     realtype t, N_Vector y, N_Vector f,
                                     int mode);

typedef int (*MRIStepInnerResetFn)(MRIStepInnerStepper stepper,
                                   realtype tR, N_Vector yR);

/*---------------------------------------------------------------
  MRI coupling data structure and associated utility routines
  ---------------------------------------------------------------*/
struct MRIStepCouplingMem
{
  int nmat;        /* number of MRI coupling matrices                   */
  int stages;      /* size of coupling matrices (stages * stages)       */
  int q;           /* method order of accuracy                          */
  int p;           /* embedding order of accuracy                       */
  realtype *c;     /* stage abscissae                                   */
  realtype ***W;   /* explicit coupling matrices [nmat][stages][stages] */
  realtype ***G;   /* implicit coupling matrices [nmat][stages][stages] */
};

typedef _SUNDIALS_STRUCT_ MRIStepCouplingMem *MRIStepCoupling;

/* Accessor routine to load built-in MRI table */
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_LoadTable(ARKODE_MRITableID imethod);

/* Utility routines to allocate/free/output coupling table structures */
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages,
                                                      MRISTEP_METHOD_TYPE type);
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_Create(int nmat,
                                                       int stages,
                                                       int q,
                                                       int p,
                                                       realtype *W,
                                                       realtype *G,
                                                       realtype *c);
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B,
                                                         int q, int p);
SUNDIALS_EXPORT MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling MRIC);
SUNDIALS_EXPORT void MRIStepCoupling_Space(MRIStepCoupling MRIC,
                                           sunindextype *liw,
                                           sunindextype *lrw);
SUNDIALS_EXPORT void MRIStepCoupling_Free(MRIStepCoupling MRIC);
SUNDIALS_EXPORT void MRIStepCoupling_Write(MRIStepCoupling MRIC,
                                           FILE *outfile);


/* ------------------------------
 * User-Supplied Function Types
 * ------------------------------ */

typedef int (*MRIStepPreInnerFn)(realtype t, N_Vector *f, int nvecs,
                                 void *user_data);

typedef int (*MRIStepPostInnerFn)(realtype t, N_Vector y, void *user_data);

/* -------------------
 * Exported Functions
 * ------------------- */

/* Create, Resize, and Reinitialization functions */
SUNDIALS_EXPORT void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, realtype t0,
                                    N_Vector y0, MRIStepInnerStepper stepper,
                                    SUNContext sunctx);

SUNDIALS_EXPORT int MRIStepResize(void *arkode_mem, N_Vector ynew,
                                  realtype t0, ARKVecResizeFn resize,
                                  void *resize_data);

SUNDIALS_EXPORT int MRIStepReInit(void* arkode_mem, ARKRhsFn fse, ARKRhsFn fsi,
                                  realtype t0, N_Vector y0);

SUNDIALS_EXPORT int MRIStepReset(void* arkode_mem, realtype tR, N_Vector yR);

/* Tolerance input functions */
SUNDIALS_EXPORT int MRIStepSStolerances(void *arkode_mem,
                                        realtype reltol,
                                        realtype abstol);
SUNDIALS_EXPORT int MRIStepSVtolerances(void *arkode_mem,
                                        realtype reltol,
                                        N_Vector abstol);
SUNDIALS_EXPORT int MRIStepWFtolerances(void *arkode_mem,
                                        ARKEwtFn efun);

/* Linear solver set function */
SUNDIALS_EXPORT int MRIStepSetLinearSolver(void *arkode_mem,
                                           SUNLinearSolver LS,
                                           SUNMatrix A);

/* Rootfinding initialization */
SUNDIALS_EXPORT int MRIStepRootInit(void *arkode_mem, int nrtfn,
                                    ARKRootFn g);

/* Optional input functions -- must be called AFTER MRIStepCreate */
SUNDIALS_EXPORT int MRIStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int MRIStepSetInterpolantType(void *arkode_mem, int itype);
SUNDIALS_EXPORT int MRIStepSetInterpolantDegree(void *arkode_mem, int degree);
SUNDIALS_EXPORT int MRIStepSetDenseOrder(void *arkode_mem, int dord);
SUNDIALS_EXPORT int MRIStepSetNonlinearSolver(void *arkode_mem,
                                              SUNNonlinearSolver NLS);
SUNDIALS_EXPORT int MRIStepSetNlsRhsFn(void *arkode_mem, ARKRhsFn nls_fs);
SUNDIALS_EXPORT int MRIStepSetLinear(void *arkode_mem, int timedepend);
SUNDIALS_EXPORT int MRIStepSetNonlinear(void *arkode_mem);
SUNDIALS_EXPORT int MRIStepSetCoupling(void *arkode_mem,
                                       MRIStepCoupling MRIC);
SUNDIALS_EXPORT int MRIStepSetMaxNumSteps(void *arkode_mem,
                                          long int mxsteps);
SUNDIALS_EXPORT int MRIStepSetNonlinCRDown(void *arkode_mem,
                                           realtype crdown);
SUNDIALS_EXPORT int MRIStepSetNonlinRDiv(void *arkode_mem,
                                         realtype rdiv);
SUNDIALS_EXPORT int MRIStepSetDeltaGammaMax(void *arkode_mem,
                                            realtype dgmax);
SUNDIALS_EXPORT int MRIStepSetLSetupFrequency(void *arkode_mem,
                                              int msbp);
SUNDIALS_EXPORT int MRIStepSetPredictorMethod(void *arkode_mem,
                                              int method);
SUNDIALS_EXPORT int MRIStepSetMaxNonlinIters(void *arkode_mem,
                                             int maxcor);
SUNDIALS_EXPORT int MRIStepSetNonlinConvCoef(void *arkode_mem,
                                             realtype nlscoef);
SUNDIALS_EXPORT int MRIStepSetMaxHnilWarns(void *arkode_mem,
                                           int mxhnil);
SUNDIALS_EXPORT int MRIStepSetStopTime(void *arkode_mem,
                                       realtype tstop);
SUNDIALS_EXPORT int MRIStepSetFixedStep(void *arkode_mem,
                                        realtype hsfixed);
SUNDIALS_EXPORT int MRIStepSetRootDirection(void *arkode_mem,
                                            int *rootdir);
SUNDIALS_EXPORT int MRIStepSetNoInactiveRootWarn(void *arkode_mem);
SUNDIALS_EXPORT int MRIStepSetErrHandlerFn(void *arkode_mem,
                                           ARKErrHandlerFn ehfun,
                                           void *eh_data);
SUNDIALS_EXPORT int MRIStepSetErrFile(void *arkode_mem,
                                      FILE *errfp);
SUNDIALS_EXPORT int MRIStepSetUserData(void *arkode_mem,
                                       void *user_data);
SUNDIALS_EXPORT int MRIStepSetDiagnostics(void *arkode_mem,
                                          FILE *diagfp);
SUNDIALS_EXPORT int MRIStepSetPostprocessStepFn(void *arkode_mem,
                                                ARKPostProcessFn ProcessStep);
SUNDIALS_EXPORT int MRIStepSetPostprocessStageFn(void *arkode_mem,
                                                 ARKPostProcessFn ProcessStage);
SUNDIALS_EXPORT int MRIStepSetPreInnerFn(void *arkode_mem,
                                         MRIStepPreInnerFn prefn);
SUNDIALS_EXPORT int MRIStepSetPostInnerFn(void *arkode_mem,
                                          MRIStepPostInnerFn postfn);
SUNDIALS_EXPORT int MRIStepSetStagePredictFn(void *arkode_mem,
                                             ARKStagePredictFn PredictStage);

/* Linear solver interface optional input functions -- must be called
   AFTER MRIStepSetLinearSolver */
SUNDIALS_EXPORT int MRIStepSetJacFn(void *arkode_mem, ARKLsJacFn jac);
SUNDIALS_EXPORT int MRIStepSetJacEvalFrequency(void *arkode_mem,
                                               long int msbj);
SUNDIALS_EXPORT int MRIStepSetLinearSolutionScaling(void *arkode_mem,
                                                    booleantype onoff);
SUNDIALS_EXPORT int MRIStepSetEpsLin(void *arkode_mem, realtype eplifac);
SUNDIALS_EXPORT int MRIStepSetLSNormFactor(void *arkode_mem,
                                           realtype nrmfac);
SUNDIALS_EXPORT int MRIStepSetPreconditioner(void *arkode_mem,
                                             ARKLsPrecSetupFn psetup,
                                             ARKLsPrecSolveFn psolve);
SUNDIALS_EXPORT int MRIStepSetJacTimes(void *arkode_mem,
                                       ARKLsJacTimesSetupFn jtsetup,
                                       ARKLsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int MRIStepSetJacTimesRhsFn(void *arkode_mem,
                                            ARKRhsFn jtimesRhsFn);
SUNDIALS_EXPORT int MRIStepSetLinSysFn(void *arkode_mem, ARKLsLinSysFn linsys);

/* Integrate the ODE over an interval in t */
SUNDIALS_EXPORT int MRIStepEvolve(void *arkode_mem, realtype tout,
                                  N_Vector yout, realtype *tret,
                                  int itask);

/* Computes the kth derivative of the y function at time t */
SUNDIALS_EXPORT int MRIStepGetDky(void *arkode_mem, realtype t,
                                  int k, N_Vector dky);

/* Utility function to update/compute y based on zcor */
SUNDIALS_EXPORT int MRIStepComputeState(void *arkode_mem, N_Vector zcor,
                                        N_Vector z);

/* Optional output functions */
SUNDIALS_EXPORT int MRIStepGetNumRhsEvals(void *arkode_mem,
                                          long int *nfse_evals,
                                          long int *nfsi_evals);
SUNDIALS_EXPORT int MRIStepGetNumLinSolvSetups(void *arkode_mem,
                                               long int *nlinsetups);
SUNDIALS_EXPORT int MRIStepGetCurrentCoupling(void *arkode_mem,
                                              MRIStepCoupling *MRIC);
SUNDIALS_EXPORT int MRIStepGetWorkSpace(void *arkode_mem,
                                        long int *lenrw,
                                        long int *leniw);
SUNDIALS_EXPORT int MRIStepGetNumSteps(void *arkode_mem,
                                       long int *nssteps);
SUNDIALS_EXPORT int MRIStepGetLastStep(void *arkode_mem,
                                       realtype *hlast);
SUNDIALS_EXPORT int MRIStepGetCurrentTime(void *arkode_mem,
                                          realtype *tcur);
SUNDIALS_EXPORT int MRIStepGetCurrentState(void *arkode_mem,
                                           N_Vector *state);
SUNDIALS_EXPORT int MRIStepGetCurrentGamma(void *arkode_mem,
                                           realtype *gamma);
SUNDIALS_EXPORT int MRIStepGetTolScaleFactor(void *arkode_mem,
                                             realtype *tolsfac);
SUNDIALS_EXPORT int MRIStepGetErrWeights(void *arkode_mem,
                                         N_Vector eweight);
SUNDIALS_EXPORT int MRIStepGetNumGEvals(void *arkode_mem,
                                        long int *ngevals);
SUNDIALS_EXPORT int MRIStepGetRootInfo(void *arkode_mem,
                                       int *rootsfound);
SUNDIALS_EXPORT int MRIStepGetLastInnerStepFlag(void *arkode_mem, int *flag);

SUNDIALS_EXPORT char *MRIStepGetReturnFlagName(long int flag);

SUNDIALS_EXPORT int MRIStepWriteParameters(void *arkode_mem, FILE *fp);

SUNDIALS_EXPORT int MRIStepWriteCoupling(void *arkode_mem, FILE *fp);

/* Nonlinear solver optional output functions */
SUNDIALS_EXPORT int MRIStepGetNonlinearSystemData(void *arkode_mem,
                                                  realtype *tcur,
                                                  N_Vector *zpred,
                                                  N_Vector *z,
                                                  N_Vector *F,
                                                  realtype *gamma,
                                                  N_Vector *sdata,
                                                  void     **user_data);
SUNDIALS_EXPORT int MRIStepGetNumNonlinSolvIters(void *arkode_mem,
                                                 long int *nniters);
SUNDIALS_EXPORT int MRIStepGetNumNonlinSolvConvFails(void *arkode_mem,
                                                     long int *nncfails);
SUNDIALS_EXPORT int MRIStepGetNonlinSolvStats(void *arkode_mem,
                                              long int *nniters,
                                              long int *nncfails);

/* Linear solver optional output functions */
SUNDIALS_EXPORT int MRIStepGetLinWorkSpace(void *arkode_mem,
                                           long int *lenrwLS,
                                           long int *leniwLS);
SUNDIALS_EXPORT int MRIStepGetNumJacEvals(void *arkode_mem,
                                          long int *njevals);
SUNDIALS_EXPORT int MRIStepGetNumPrecEvals(void *arkode_mem,
                                           long int *npevals);
SUNDIALS_EXPORT int MRIStepGetNumPrecSolves(void *arkode_mem,
                                            long int *npsolves);
SUNDIALS_EXPORT int MRIStepGetNumLinIters(void *arkode_mem,
                                          long int *nliters);
SUNDIALS_EXPORT int MRIStepGetNumLinConvFails(void *arkode_mem,
                                              long int *nlcfails);
SUNDIALS_EXPORT int MRIStepGetNumJTSetupEvals(void *arkode_mem,
                                              long int *njtsetups);
SUNDIALS_EXPORT int MRIStepGetNumJtimesEvals(void *arkode_mem,
                                             long int *njvevals);
SUNDIALS_EXPORT int MRIStepGetNumLinRhsEvals(void *arkode_mem,
                                             long int *nfevalsLS);
SUNDIALS_EXPORT int MRIStepGetLastLinFlag(void *arkode_mem,
                                          long int *flag);

SUNDIALS_EXPORT char *MRIStepGetLinReturnFlagName(long int flag);


/* Free function */
SUNDIALS_EXPORT void MRIStepFree(void **arkode_mem);

/* Output the MRIStep memory structure (useful when debugging) */
SUNDIALS_EXPORT void MRIStepPrintMem(void* arkode_mem, FILE* outfile);

/* Custom inner stepper functions */
SUNDIALS_EXPORT int MRIStepInnerStepper_Create(SUNContext sunctx,
                                               MRIStepInnerStepper *stepper);

SUNDIALS_EXPORT int MRIStepInnerStepper_Free(MRIStepInnerStepper *stepper);

SUNDIALS_EXPORT int MRIStepInnerStepper_SetContent(MRIStepInnerStepper stepper,
                                                   void *content);

SUNDIALS_EXPORT int MRIStepInnerStepper_GetContent(MRIStepInnerStepper stepper,
                                                   void **content);

SUNDIALS_EXPORT int MRIStepInnerStepper_SetEvolveFn(MRIStepInnerStepper stepper,
                                                    MRIStepInnerEvolveFn fn);

SUNDIALS_EXPORT int MRIStepInnerStepper_SetFullRhsFn(MRIStepInnerStepper stepper,
                                                     MRIStepInnerFullRhsFn fn);

SUNDIALS_EXPORT int MRIStepInnerStepper_SetResetFn(MRIStepInnerStepper stepper,
                                                   MRIStepInnerResetFn fn);

SUNDIALS_EXPORT int MRIStepInnerStepper_AddForcing(MRIStepInnerStepper stepper,
                                                   realtype t, N_Vector f);

SUNDIALS_EXPORT int MRIStepInnerStepper_GetForcingData(MRIStepInnerStepper stepper,
                                                       realtype *tshift,
                                                       realtype *tscale,
                                                       N_Vector **forcing,
                                                       int *nforcing);

#ifdef __cplusplus
}
#endif

#endif
