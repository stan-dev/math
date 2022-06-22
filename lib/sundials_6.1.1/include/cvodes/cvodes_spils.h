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
 * Linear Solver interface in CVODES; these routines now just wrap
 * the updated CVODES generic linear solver interface in cvodes_ls.h.
 * -----------------------------------------------------------------*/

#ifndef _CVSSPILS_H
#define _CVSSPILS_H

#include <cvodes/cvodes_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  Function Types (typedefs for equivalent types in cvodes_ls.h)
  ===============================================================*/

typedef CVLsPrecSetupFn CVSpilsPrecSetupFn;
typedef CVLsPrecSolveFn CVSpilsPrecSolveFn;
typedef CVLsJacTimesSetupFn CVSpilsJacTimesSetupFn;
typedef CVLsJacTimesVecFn CVSpilsJacTimesVecFn;
typedef CVLsPrecSetupFnB CVSpilsPrecSetupFnB;
typedef CVLsPrecSetupFnBS CVSpilsPrecSetupFnBS;
typedef CVLsPrecSolveFnB CVSpilsPrecSolveFnB;
typedef CVLsPrecSolveFnBS CVSpilsPrecSolveFnBS;
typedef CVLsJacTimesSetupFnB CVSpilsJacTimesSetupFnB;
typedef CVLsJacTimesSetupFnBS CVSpilsJacTimesSetupFnBS;
typedef CVLsJacTimesVecFnB CVSpilsJacTimesVecFnB;
typedef CVLsJacTimesVecFnBS CVSpilsJacTimesVecFnBS;

/*====================================================================
  Exported Functions (wrappers for equivalent routines in cvodes_ls.h)
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

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetLinearSolverB instead")
int CVSpilsSetLinearSolverB(void *cvode_mem, int which,
                            SUNLinearSolver LS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetEpsLinB instead")
int CVSpilsSetEpsLinB(void *cvode_mem, int which, realtype eplifacB);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetPreconditionerB instead")
int CVSpilsSetPreconditionerB(void *cvode_mem, int which,
                              CVSpilsPrecSetupFnB psetB,
                              CVSpilsPrecSolveFnB psolveB);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetPreconditionerBS instead")
int CVSpilsSetPreconditionerBS(void *cvode_mem, int which,
                               CVSpilsPrecSetupFnBS psetBS,
                               CVSpilsPrecSolveFnBS psolveBS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetJacTimesB instead")
int CVSpilsSetJacTimesB(void *cvode_mem, int which,
                        CVSpilsJacTimesSetupFnB jtsetupB,
                        CVSpilsJacTimesVecFnB jtimesB);

SUNDIALS_DEPRECATED_EXPORT_MSG("use CVodeSetJacTimesBS instead")
int CVSpilsSetJacTimesBS(void *cvode_mem, int which,
                         CVSpilsJacTimesSetupFnBS jtsetupBS,
                         CVSpilsJacTimesVecFnBS jtimesBS);


#ifdef __cplusplus
}
#endif

#endif
