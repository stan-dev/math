/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on code sundials_spfgmr.h by: Daniel R. Reynolds and
 *    Hilari C. Tiedeman @ SMU
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
 * This is the header file for the SPFGMR implementation of the
 * SUNLINSOL module, SUNLINSOL_SPFGMR.  The SPFGMR algorithm is based
 * on the Scaled Preconditioned FGMRES (Flexible Generalized Minimal
 * Residual) method [Y. Saad, SIAM J. Sci. Comput., 1993].
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SPFGMR_H
#define _SUNLINSOL_SPFGMR_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SPFGMR solver parameters */
#define SUNSPFGMR_MAXL_DEFAULT    5
#define SUNSPFGMR_MAXRS_DEFAULT   0
#define SUNSPFGMR_GSTYPE_DEFAULT  SUN_MODIFIED_GS

/* -----------------------------------------
 * SPFGMR Implementation of SUNLinearSolver
 * ----------------------------------------- */

struct _SUNLinearSolverContent_SPFGMR {
  int maxl;
  int pretype;
  int gstype;
  int max_restarts;
  booleantype zeroguess;
  int numiters;
  realtype resnorm;
  int last_flag;

  SUNATimesFn ATimes;
  void* ATData;
  SUNPSetupFn Psetup;
  SUNPSolveFn Psolve;
  void* PData;

  N_Vector s1;
  N_Vector s2;
  N_Vector *V;
  N_Vector *Z;
  realtype **Hes;
  realtype *givens;
  N_Vector xcor;
  realtype *yg;
  N_Vector vtemp;

  realtype *cv;
  N_Vector *Xv;

  int print_level;
  FILE* info_file;
};

typedef struct _SUNLinearSolverContent_SPFGMR *SUNLinearSolverContent_SPFGMR;

/* ----------------------------------------
 * Exported Functions for SUNLINSOL_SPFGMR
 * ---------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SPFGMR(N_Vector y,
                                                 int pretype,
                                                 int maxl,
                                                 SUNContext sunctx);
SUNDIALS_EXPORT int SUNLinSol_SPFGMRSetPrecType(SUNLinearSolver S,
                                                int pretype);
SUNDIALS_EXPORT int SUNLinSol_SPFGMRSetGSType(SUNLinearSolver S,
                                              int gstype);
SUNDIALS_EXPORT int SUNLinSol_SPFGMRSetMaxRestarts(SUNLinearSolver S,
                                                   int maxrs);
SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SPFGMR(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_SPFGMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_SPFGMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_SPFGMR(SUNLinearSolver S, void* A_data,
                                              SUNATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_SPFGMR(SUNLinearSolver S,
                                                      void* P_data,
                                                      SUNPSetupFn Pset,
                                                      SUNPSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_SPFGMR(SUNLinearSolver S,
                                                      N_Vector s1,
                                                      N_Vector s2);
SUNDIALS_EXPORT int SUNLinSolSetZeroGuess_SPFGMR(SUNLinearSolver S,
                                                 booleantype onoff);
SUNDIALS_EXPORT int SUNLinSolSetup_SPFGMR(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_SPFGMR(SUNLinearSolver S, SUNMatrix A,
                                          N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_SPFGMR(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_SPFGMR(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_SPFGMR(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_SPFGMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_SPFGMR(SUNLinearSolver S,
                                          long int *lenrwLS,
                                          long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_SPFGMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetInfoFile_SPFGMR(SUNLinearSolver LS,
                                                FILE* info_file);
SUNDIALS_EXPORT int SUNLinSolSetPrintLevel_SPFGMR(SUNLinearSolver LS,
                                                  int print_level);


#ifdef __cplusplus
}
#endif

#endif
