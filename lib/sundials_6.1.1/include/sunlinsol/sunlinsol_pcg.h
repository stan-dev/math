/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ SMU
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
 * This is the header file for the PCG implementation of the
 * SUNLINSOL module, SUNLINSOL_PCG.  The PCG algorithm is based
 * on the Preconditioned Conjugate Gradient.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_PCG_H
#define _SUNLINSOL_PCG_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default PCG solver parameters */
#define SUNPCG_MAXL_DEFAULT    5

/* --------------------------------------
 * PCG Implementation of SUNLinearSolver
 * -------------------------------------- */

struct _SUNLinearSolverContent_PCG {
  int maxl;
  int pretype;
  booleantype zeroguess;
  int numiters;
  realtype resnorm;
  int last_flag;

  SUNATimesFn ATimes;
  void* ATData;
  SUNPSetupFn Psetup;
  SUNPSolveFn Psolve;
  void* PData;

  N_Vector s;
  N_Vector r;
  N_Vector p;
  N_Vector z;
  N_Vector Ap;

  int print_level;
  FILE* info_file;
};

typedef struct _SUNLinearSolverContent_PCG *SUNLinearSolverContent_PCG;


/* -------------------------------------
 * Exported Functions for SUNLINSOL_PCG
 * ------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_PCG(N_Vector y,
                                              int pretype,
                                              int maxl,
                                              SUNContext sunctx);
SUNDIALS_EXPORT int SUNLinSol_PCGSetPrecType(SUNLinearSolver S,
                                             int pretype);
SUNDIALS_EXPORT int SUNLinSol_PCGSetMaxl(SUNLinearSolver S,
                                         int maxl);

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_PCG(SUNLinearSolver S, void* A_data,
                                           SUNATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_PCG(SUNLinearSolver S,
                                                   void* P_data,
                                                   SUNPSetupFn Pset,
                                                   SUNPSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_PCG(SUNLinearSolver S,
                                                   N_Vector s,
                                                   N_Vector nul);
SUNDIALS_EXPORT int SUNLinSolSetZeroGuess_PCG(SUNLinearSolver S,
                                              booleantype onoff);
SUNDIALS_EXPORT int SUNLinSolSetup_PCG(SUNLinearSolver S, SUNMatrix nul);
SUNDIALS_EXPORT int SUNLinSolSolve_PCG(SUNLinearSolver S, SUNMatrix nul,
                                       N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_PCG(SUNLinearSolver S,
                                       long int *lenrwLS,
                                       long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetInfoFile_PCG(SUNLinearSolver LS,
                                             FILE* info_file);
SUNDIALS_EXPORT int SUNLinSolSetPrintLevel_PCG(SUNLinearSolver LS,
                                               int print_level);

#ifdef __cplusplus
}
#endif

#endif
