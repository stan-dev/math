/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the ARKode ERKStep module.
 * -----------------------------------------------------------------*/

#ifndef _ERKSTEP_H
#define _ERKSTEP_H

#include <sundials/sundials_nvector.h>
#include <arkode/arkode.h>
#include <arkode/arkode_butcher_erk.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * ERKStep Constants
 * ----------------- */

/* Default Butcher tables for each order */

#define DEFAULT_ERK_2           HEUN_EULER_2_1_2
#define DEFAULT_ERK_3           BOGACKI_SHAMPINE_4_2_3
#define DEFAULT_ERK_4           ZONNEVELD_5_3_4
#define DEFAULT_ERK_5           CASH_KARP_6_4_5
#define DEFAULT_ERK_6           VERNER_8_5_6
#define DEFAULT_ERK_8           FEHLBERG_13_7_8


/* -------------------
 * Exported Functions
 * ------------------- */

/* Create, Resize, and Reinitialization functions */
SUNDIALS_EXPORT void* ERKStepCreate(ARKRhsFn f, realtype t0,
                                    N_Vector y0);

SUNDIALS_EXPORT int ERKStepResize(void *arkode_mem, N_Vector ynew,
                                  realtype hscale, realtype t0,
                                  ARKVecResizeFn resize,
                                  void *resize_data);

SUNDIALS_EXPORT int ERKStepReInit(void* arkode_mem, ARKRhsFn f,
                                  realtype t0, N_Vector y0);

/* Tolerance input functions */
SUNDIALS_EXPORT int ERKStepSStolerances(void *arkode_mem,
                                        realtype reltol,
                                        realtype abstol);
SUNDIALS_EXPORT int ERKStepSVtolerances(void *arkode_mem,
                                        realtype reltol,
                                        N_Vector abstol);
SUNDIALS_EXPORT int ERKStepWFtolerances(void *arkode_mem,
                                        ARKEwtFn efun);

/* Rootfinding initialization */
SUNDIALS_EXPORT int ERKStepRootInit(void *arkode_mem, int nrtfn,
                                    ARKRootFn g);

/* Optional input functions -- must be called AFTER ERKStepCreate */
SUNDIALS_EXPORT int ERKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int ERKStepSetOrder(void *arkode_mem, int maxord);
SUNDIALS_EXPORT int ERKStepSetDenseOrder(void *arkode_mem, int dord);
SUNDIALS_EXPORT int ERKStepSetTable(void *arkode_mem,
                                    ARKodeButcherTable B);
SUNDIALS_EXPORT int ERKStepSetTableNum(void *arkode_mem, int itable);
SUNDIALS_EXPORT int ERKStepSetCFLFraction(void *arkode_mem,
                                          realtype cfl_frac);
SUNDIALS_EXPORT int ERKStepSetSafetyFactor(void *arkode_mem,
                                           realtype safety);
SUNDIALS_EXPORT int ERKStepSetErrorBias(void *arkode_mem,
                                        realtype bias);
SUNDIALS_EXPORT int ERKStepSetMaxGrowth(void *arkode_mem,
                                        realtype mx_growth);
SUNDIALS_EXPORT int ERKStepSetFixedStepBounds(void *arkode_mem,
                                              realtype lb, realtype ub);
SUNDIALS_EXPORT int ERKStepSetAdaptivityMethod(void *arkode_mem,
                                               int imethod,
                                               int idefault, int pq,
                                               realtype *adapt_params);
SUNDIALS_EXPORT int ERKStepSetAdaptivityFn(void *arkode_mem,
                                           ARKAdaptFn hfun,
                                           void *h_data);
SUNDIALS_EXPORT int ERKStepSetMaxFirstGrowth(void *arkode_mem,
                                             realtype etamx1);
SUNDIALS_EXPORT int ERKStepSetMaxEFailGrowth(void *arkode_mem,
                                             realtype etamxf);
SUNDIALS_EXPORT int ERKStepSetSmallNumEFails(void *arkode_mem,
                                             int small_nef);
SUNDIALS_EXPORT int ERKStepSetStabilityFn(void *arkode_mem,
                                          ARKExpStabFn EStab,
                                          void *estab_data);
SUNDIALS_EXPORT int ERKStepSetMaxErrTestFails(void *arkode_mem,
                                              int maxnef);
SUNDIALS_EXPORT int ERKStepSetMaxNumSteps(void *arkode_mem,
                                          long int mxsteps);
SUNDIALS_EXPORT int ERKStepSetMaxHnilWarns(void *arkode_mem,
                                           int mxhnil);
SUNDIALS_EXPORT int ERKStepSetInitStep(void *arkode_mem,
                                       realtype hin);
SUNDIALS_EXPORT int ERKStepSetMinStep(void *arkode_mem,
                                      realtype hmin);
SUNDIALS_EXPORT int ERKStepSetMaxStep(void *arkode_mem,
                                      realtype hmax);
SUNDIALS_EXPORT int ERKStepSetStopTime(void *arkode_mem,
                                       realtype tstop);
SUNDIALS_EXPORT int ERKStepSetFixedStep(void *arkode_mem,
                                        realtype hfixed);

SUNDIALS_EXPORT int ERKStepSetRootDirection(void *arkode_mem,
                                            int *rootdir);
SUNDIALS_EXPORT int ERKStepSetNoInactiveRootWarn(void *arkode_mem);

SUNDIALS_EXPORT int ERKStepSetErrHandlerFn(void *arkode_mem,
                                           ARKErrHandlerFn ehfun,
                                           void *eh_data);
SUNDIALS_EXPORT int ERKStepSetErrFile(void *arkode_mem,
                                      FILE *errfp);
SUNDIALS_EXPORT int ERKStepSetUserData(void *arkode_mem,
                                       void *user_data);
SUNDIALS_EXPORT int ERKStepSetDiagnostics(void *arkode_mem,
                                          FILE *diagfp);

SUNDIALS_EXPORT int ERKStepSetPostprocessStepFn(void *arkode_mem,
                                                ARKPostProcessStepFn ProcessStep);


/* Integrate the ODE over an interval in t */
SUNDIALS_EXPORT int ERKStepEvolve(void *arkode_mem, realtype tout,
                                  N_Vector yout, realtype *tret,
                                  int itask);

/* Computes the kth derivative of the y function at time t */
SUNDIALS_EXPORT int ERKStepGetDky(void *arkode_mem, realtype t,
                                  int k, N_Vector dky);

/* Optional output functions */
SUNDIALS_EXPORT int ERKStepGetNumExpSteps(void *arkode_mem,
                                          long int *expsteps);
SUNDIALS_EXPORT int ERKStepGetNumAccSteps(void *arkode_mem,
                                          long int *accsteps);
SUNDIALS_EXPORT int ERKStepGetNumStepAttempts(void *arkode_mem,
                                              long int *step_attempts);
SUNDIALS_EXPORT int ERKStepGetNumRhsEvals(void *arkode_mem,
                                          long int *nfevals);
SUNDIALS_EXPORT int ERKStepGetNumErrTestFails(void *arkode_mem,
                                              long int *netfails);
SUNDIALS_EXPORT int ERKStepGetCurrentButcherTable(void *arkode_mem,
                                                  ARKodeButcherTable *B);
SUNDIALS_EXPORT int ERKStepGetEstLocalErrors(void *arkode_mem,
                                             N_Vector ele);
SUNDIALS_EXPORT int ERKStepGetWorkSpace(void *arkode_mem,
                                        long int *lenrw,
                                        long int *leniw);
SUNDIALS_EXPORT int ERKStepGetNumSteps(void *arkode_mem,
                                       long int *nsteps);
SUNDIALS_EXPORT int ERKStepGetActualInitStep(void *arkode_mem,
                                             realtype *hinused);
SUNDIALS_EXPORT int ERKStepGetLastStep(void *arkode_mem,
                                       realtype *hlast);
SUNDIALS_EXPORT int ERKStepGetCurrentStep(void *arkode_mem,
                                          realtype *hcur);
SUNDIALS_EXPORT int ERKStepGetCurrentTime(void *arkode_mem,
                                          realtype *tcur);
SUNDIALS_EXPORT int ERKStepGetTolScaleFactor(void *arkode_mem,
                                             realtype *tolsfac);
SUNDIALS_EXPORT int ERKStepGetErrWeights(void *arkode_mem,
                                         N_Vector eweight);
SUNDIALS_EXPORT int ERKStepGetNumGEvals(void *arkode_mem,
                                        long int *ngevals);
SUNDIALS_EXPORT int ERKStepGetRootInfo(void *arkode_mem,
                                       int *rootsfound);
SUNDIALS_EXPORT char *ERKStepGetReturnFlagName(long int flag);

SUNDIALS_EXPORT int ERKStepWriteParameters(void *arkode_mem, FILE *fp);

SUNDIALS_EXPORT int ERKStepWriteButcher(void *arkode_mem, FILE *fp);


/* Grouped optional output functions */
SUNDIALS_EXPORT int ERKStepGetTimestepperStats(void *arkode_mem,
                                               long int *expsteps,
                                               long int *accsteps,
                                               long int *step_attempts,
                                               long int *nfevals,
                                               long int *netfails);
SUNDIALS_EXPORT int ERKStepGetStepStats(void *arkode_mem,
                                        long int *nsteps,
                                        realtype *hinused,
                                        realtype *hlast,
                                        realtype *hcur,
                                        realtype *tcur);


/* Free function */
SUNDIALS_EXPORT void ERKStepFree(void **arkode_mem);

/* Output the ERKStep memory structure (useful when debugging) */
SUNDIALS_EXPORT void ERKStepPrintMem(void* arkode_mem, FILE* outfile);


#ifdef __cplusplus
}
#endif

#endif
