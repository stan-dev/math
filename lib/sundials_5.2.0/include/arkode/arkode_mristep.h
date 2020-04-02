/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
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
#include <arkode/arkode.h>
#include <arkode/arkode_butcher_erk.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * MRIStep Constants
 * ----------------- */

/* Inner stepper module identifiers */
typedef enum {
  MRISTEP_ARKSTEP
} MRISTEP_ID;

/* Default slow (outer) Butcher tables for each order */
#define DEFAULT_MRI_TABLE_3  KNOTH_WOLKE_3_3

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
SUNDIALS_EXPORT void* MRIStepCreate(ARKRhsFn fs, realtype t0, N_Vector y0,
                                    MRISTEP_ID inner_step_id,
                                    void* inner_step_mem);

SUNDIALS_EXPORT int MRIStepResize(void *arkode_mem, N_Vector ynew,
                                  realtype t0, ARKVecResizeFn resize,
                                  void *resize_data);

SUNDIALS_EXPORT int MRIStepReInit(void* arkode_mem, ARKRhsFn fs, realtype t0,
                                  N_Vector y0);

/* Rootfinding initialization */
SUNDIALS_EXPORT int MRIStepRootInit(void *arkode_mem, int nrtfn,
                                    ARKRootFn g);

/* Optional input functions -- must be called AFTER MRIStepCreate */
SUNDIALS_EXPORT int MRIStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int MRIStepSetInterpolantType(void *arkode_mem, int itype);
SUNDIALS_EXPORT int MRIStepSetInterpolantDegree(void *arkode_mem, int degree);
SUNDIALS_EXPORT int MRIStepSetDenseOrder(void *arkode_mem, int dord);
SUNDIALS_EXPORT int MRIStepSetTable(void *arkode_mem, int q,
                                    ARKodeButcherTable B);
SUNDIALS_EXPORT int MRIStepSetTableNum(void *arkode_mem,
                                       int itable);
SUNDIALS_EXPORT int MRIStepSetMaxNumSteps(void *arkode_mem,
                                          long int mxsteps);
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


/* Integrate the ODE over an interval in t */
SUNDIALS_EXPORT int MRIStepEvolve(void *arkode_mem, realtype tout,
                                  N_Vector yout, realtype *tret,
                                  int itask);

/* Computes the kth derivative of the y function at time t */
SUNDIALS_EXPORT int MRIStepGetDky(void *arkode_mem, realtype t,
                                  int k, N_Vector dky);

/* Optional output functions */
SUNDIALS_EXPORT int MRIStepGetNumRhsEvals(void *arkode_mem,
                                          long int *nfs_evals);
SUNDIALS_EXPORT int MRIStepGetCurrentButcherTables(void *arkode_mem,
                                                   ARKodeButcherTable *B);
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
                                           N_Vector *ycur);
SUNDIALS_EXPORT int MRIStepGetNumGEvals(void *arkode_mem,
                                        long int *ngevals);
SUNDIALS_EXPORT int MRIStepGetRootInfo(void *arkode_mem,
                                       int *rootsfound);
SUNDIALS_EXPORT int MRIStepGetLastInnerStepFlag(void *arkode_mem, int *flag);

SUNDIALS_EXPORT char *MRIStepGetReturnFlagName(long int flag);

SUNDIALS_EXPORT int MRIStepWriteParameters(void *arkode_mem, FILE *fp);

SUNDIALS_EXPORT int MRIStepWriteButcher(void *arkode_mem, FILE *fp);

/* Free function */
SUNDIALS_EXPORT void MRIStepFree(void **arkode_mem);

/* Output the MRIStep memory structure (useful when debugging) */
SUNDIALS_EXPORT void MRIStepPrintMem(void* arkode_mem, FILE* outfile);


#ifdef __cplusplus
}
#endif

#endif
