/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * Header file for the deprecated direct linear solver interface in
 * CVODES; these routines now just wrap the updated CVODE generic
 * linear solver interface in cvodes_ls.h.
 * -----------------------------------------------------------------*/

#ifndef _CVSDLS_H
#define _CVSDLS_H

#include <cvodes/cvodes_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*=================================================================
  Function Types (typedefs for equivalent types in cvodes_ls.h)
  =================================================================*/

typedef CVLsJacFn CVDlsJacFn;
typedef CVLsJacFnB CVDlsJacFnB;
typedef CVLsJacFnBS CVDlsJacFnBS;

/*====================================================================
  Exported Functions (wrappers for equivalent routines in cvodes_ls.h)
  ====================================================================*/

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetLinearSolver instead")
int CVDlsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS,
                         SUNMatrix A);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetJacFn instead")
int CVDlsSetJacFn(void *cvode_mem, CVDlsJacFn jac);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetLinWorkSpace instead")
int CVDlsGetWorkSpace(void *cvode_mem, long int *lenrwLS,
                      long int *leniwLS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetNumJacEvals instead")
int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetNumLinRhsEvals instead")
int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetLastLinFlag instead")
int CVDlsGetLastFlag(void *cvode_mem, long int *flag);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetLinReturnFlagName instead")
char *CVDlsGetReturnFlagName(long int flag);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetLinearSolverB instead")
int CVDlsSetLinearSolverB(void *cvode_mem, int which,
                          SUNLinearSolver LS, SUNMatrix A);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetJacFnB instead")
int CVDlsSetJacFnB(void *cvode_mem, int which, CVDlsJacFnB jacB);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetJacFnBS instead")
int CVDlsSetJacFnBS(void *cvode_mem, int which, CVDlsJacFnBS jacBS);


#ifdef __cplusplus
}
#endif

#endif
