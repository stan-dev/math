/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * This is the header file for the SUNLINEARSOLVER class implementation using
 * the Intel oneAPI Math Kernel Library (oneMKL).
 * ---------------------------------------------------------------------------*/

#ifndef _SUNLINSOL_ONEMKLDENSE_H
#define _SUNLINSOL_ONEMKLDENSE_H

#include <CL/sycl.hpp>

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

struct _SUNLinearSolverContent_OneMklDense
{
  int             last_flag;      /* last error code returned */
  sunindextype    rows;           /* number of rows in A      */
  SUNMemory       pivots;         /* pivots array             */
  sunindextype    f_scratch_size; /* num scratchpad elements  */
  SUNMemory       f_scratchpad;   /* scratchpad memory        */
  sunindextype    s_scratch_size; /* num scratchpad elements  */
  SUNMemory       s_scratchpad;   /* scratchpad memory        */
  SUNMemoryType   mem_type;       /* memory type              */
  SUNMemoryHelper mem_helper;     /* memory helper            */
  ::sycl::queue*  queue;          /* operation queue          */
};

typedef struct _SUNLinearSolverContent_OneMklDense *SUNLinearSolverContent_OneMklDense;

/* ---------------------------------------------------------------------------
 * Implementation specific functions
 * ---------------------------------------------------------------------------*/

SUNDIALS_EXPORT
SUNLinearSolver SUNLinSol_OneMklDense(N_Vector y, SUNMatrix A, SUNContext sunctx);

SUNDIALS_STATIC_INLINE
SUNLinearSolver_Type SUNLinSolGetType_OneMklDense(SUNLinearSolver S) { return SUNLINEARSOLVER_DIRECT; };

SUNDIALS_STATIC_INLINE
SUNLinearSolver_ID SUNLinSolGetID_OneMklDense(SUNLinearSolver S) { return SUNLINEARSOLVER_ONEMKLDENSE; };

SUNDIALS_EXPORT
int SUNLinSolInitialize_OneMklDense(SUNLinearSolver S);

SUNDIALS_EXPORT
int SUNLinSolSetup_OneMklDense(SUNLinearSolver S, SUNMatrix A);

SUNDIALS_EXPORT
int SUNLinSolSolve_OneMklDense(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                               N_Vector b, realtype tol);

SUNDIALS_EXPORT
sunindextype SUNLinSolLastFlag_OneMklDense(SUNLinearSolver S);

SUNDIALS_EXPORT
int SUNLinSolSpace_OneMklDense(SUNLinearSolver S, long int *lenrwLS,
                               long int *leniwLS);

SUNDIALS_EXPORT
int SUNLinSolFree_OneMklDense(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
