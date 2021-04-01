/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on codes sundials_superlumt_impl.h and <solver>_superlumt.h
 *     written by Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the SuperLUMT implementation of the
 * SUNLINSOL module, SUNLINSOL_SUPERLUMT.
 *
 * Note:
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SLUMT_H
#define _SUNLINSOL_SLUMT_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>

/* Assume SuperLU_MT library was built with compatible index type */
#if defined(SUNDIALS_INT64_T)
#define _LONGINT
#endif

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SuperLU_MT solver parameters */
#define SUNSLUMT_ORDERING_DEFAULT  3     /* COLAMD */

/* Interfaces to match 'realtype' with the correct SuperLUMT functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#ifndef _SLUMT_H
#define _SLUMT_H
#include "slu_mt_ddefs.h"
#endif
#define xgstrs                  dgstrs
#define pxgstrf                 pdgstrf
#define pxgstrf_init            pdgstrf_init
#define xCreate_Dense_Matrix    dCreate_Dense_Matrix
#define xCreate_CompCol_Matrix  dCreate_CompCol_Matrix
#elif defined(SUNDIALS_SINGLE_PRECISION)
#ifndef _SLUMT_H
#define _SLUMT_H
#include "slu_mt_sdefs.h"
#endif
#define xgstrs                  sgstrs
#define pxgstrf                 psgstrf
#define pxgstrf_init            psgstrf_init
#define xCreate_Dense_Matrix    sCreate_Dense_Matrix
#define xCreate_CompCol_Matrix  sCreate_CompCol_Matrix
#else  /* incompatible sunindextype for SuperLUMT */
#error  Incompatible realtype for SuperLUMT
#endif


/* --------------------------------------------
 * SuperLUMT Implementation of SUNLinearSolver
 * -------------------------------------------- */

struct _SUNLinearSolverContent_SuperLUMT {
  int          last_flag;
  int          first_factorize;
  SuperMatrix  *A, *AC, *L, *U, *B;
  Gstat_t      *Gstat;
  sunindextype *perm_r, *perm_c;
  sunindextype N;
  int          num_threads;
  realtype     diag_pivot_thresh;
  int          ordering;
  superlumt_options_t *options;
};

typedef struct _SUNLinearSolverContent_SuperLUMT *SUNLinearSolverContent_SuperLUMT;


/* -------------------------------------------
 * Exported Functions for SUNLINSOL_SUPERLUMT
 * ------------------------------------------- */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SuperLUMT(N_Vector y,
                                                    SUNMatrix A,
                                                    int num_threads);
SUNDIALS_EXPORT int SUNLinSol_SuperLUMTSetOrdering(SUNLinearSolver S,
                                                   int ordering_choice);

/* deprecated */
SUNDIALS_EXPORT SUNLinearSolver SUNSuperLUMT(N_Vector y, SUNMatrix A,
                                             int num_threads);
/* deprecated */
SUNDIALS_EXPORT int SUNSuperLUMTSetOrdering(SUNLinearSolver S,
                                            int ordering_choice);

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SuperLUMT(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_SuperLUMT(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_SuperLUMT(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_SuperLUMT(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_SuperLUMT(SUNLinearSolver S, SUNMatrix A,
                                       N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_SuperLUMT(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_SuperLUMT(SUNLinearSolver S,
                                             long int *lenrwLS,
                                             long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_SuperLUMT(SUNLinearSolver S);


#ifdef __cplusplus
}
#endif

#endif
