/* ----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * This is the header file for ARKode's linear solver interface.
 * ----------------------------------------------------------------*/

#ifndef _ARKLS_H
#define _ARKLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*=================================================================
  ARKLS Constants
  =================================================================*/

#define ARKLS_SUCCESS           0
#define ARKLS_MEM_NULL         -1
#define ARKLS_LMEM_NULL        -2
#define ARKLS_ILL_INPUT        -3
#define ARKLS_MEM_FAIL         -4
#define ARKLS_PMEM_NULL        -5
#define ARKLS_MASSMEM_NULL     -6
#define ARKLS_JACFUNC_UNRECVR  -7
#define ARKLS_JACFUNC_RECVR    -8
#define ARKLS_MASSFUNC_UNRECVR -9
#define ARKLS_MASSFUNC_RECVR   -10
#define ARKLS_SUNMAT_FAIL      -11
#define ARKLS_SUNLS_FAIL       -12


/*=================================================================
  ARKLS user-supplied function prototypes
  =================================================================*/

typedef int (*ARKLsJacFn)(realtype t, N_Vector y, N_Vector fy,
                          SUNMatrix Jac, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

typedef int (*ARKLsMassFn)(realtype t, SUNMatrix M, void *user_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

typedef int (*ARKLsPrecSetupFn)(realtype t, N_Vector y,
                                N_Vector fy, booleantype jok,
                                booleantype *jcurPtr,
                                realtype gamma, void *user_data);

typedef int (*ARKLsPrecSolveFn)(realtype t, N_Vector y,
                                N_Vector fy, N_Vector r,
                                N_Vector z, realtype gamma,
                                realtype delta, int lr,
                                void *user_data);

typedef int (*ARKLsJacTimesSetupFn)(realtype t, N_Vector y,
                                    N_Vector fy, void *user_data);

typedef int (*ARKLsJacTimesVecFn)(N_Vector v, N_Vector Jv,
                                  realtype t, N_Vector y,
                                  N_Vector fy, void *user_data,
                                  N_Vector tmp);

typedef int (*ARKLsMassTimesSetupFn)(realtype t, void *mtimes_data);


typedef int (*ARKLsMassTimesVecFn)(N_Vector v, N_Vector Mv,
                                   realtype t, void *mtimes_data);

typedef int (*ARKLsMassPrecSetupFn)(realtype t, void *user_data);

typedef int (*ARKLsMassPrecSolveFn)(realtype t, N_Vector r,
                                    N_Vector z, realtype delta,
                                    int lr, void *user_data);


#ifdef __cplusplus
}
#endif

#endif
