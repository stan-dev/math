/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * Header file for the deprecated Scaled, Preconditioned Iterative
 * Linear Solver interface in CVODE; these routines now just wrap
 * the updated CVODE generic linear solver interface in cvode_ls.h.
 * -----------------------------------------------------------------*/

#ifndef _CVSPILS_H
#define _CVSPILS_H

#include <cvode/cvode_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  Function Types (typedefs for equivalent types in cvode_ls.h)
  ===============================================================*/

typedef CVLsPrecSetupFn CVSpilsPrecSetupFn;
typedef CVLsPrecSolveFn CVSpilsPrecSolveFn;
typedef CVLsJacTimesSetupFn CVSpilsJacTimesSetupFn;
typedef CVLsJacTimesVecFn CVSpilsJacTimesVecFn;

/*====================================================================
  Exported Functions (wrappers for equivalent routines in cvode_ls.h)
  ====================================================================*/

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetLinearSolver instead")
int CVSpilsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetEpsLin instead")
int CVSpilsSetEpsLin(void *cvode_mem, realtype eplifac);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetPreconditioner instead")
int CVSpilsSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn pset,
                             CVSpilsPrecSolveFn psolve);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetJacTimes instead")
int CVSpilsSetJacTimes(void *cvode_mem, CVSpilsJacTimesSetupFn jtsetup,
                       CVSpilsJacTimesVecFn jtimes);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetLinWorkSpace instead")
int CVSpilsGetWorkSpace(void *cvode_mem, long int *lenrwLS,
                        long int *leniwLS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetNumPrecEvals instead")
int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetNumPrecSolves instead")
int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetNumLinIters instead")
int CVSpilsGetNumLinIters(void *cvode_mem, long int *nliters);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetNumConvFails instead")
int CVSpilsGetNumConvFails(void *cvode_mem, long int *nlcfails);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetNumJTSetupEvals instead")
int CVSpilsGetNumJTSetupEvals(void *cvode_mem, long int *njtsetups);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetNumJtimesEvals instead")
int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetNumLinRhsEvals instead")
int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetLastLinFlag instead")
int CVSpilsGetLastFlag(void *cvode_mem, long int *flag);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeGetLinReturnFlagName instead")
char *CVSpilsGetReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif
