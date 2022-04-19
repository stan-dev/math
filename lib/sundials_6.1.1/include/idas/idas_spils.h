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
 * Linear Solver interface in IDAS; these routines now just wrap
 * the updated IDA generic linear solver interface in idas_ls.h.
 * -----------------------------------------------------------------*/

#ifndef _IDASSPILS_H
#define _IDASSPILS_H

#include <idas/idas_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  Function Types (typedefs for equivalent types in idas_ls.h)
  ===============================================================*/

typedef IDALsPrecSetupFn IDASpilsPrecSetupFn;
typedef IDALsPrecSolveFn IDASpilsPrecSolveFn;
typedef IDALsJacTimesSetupFn IDASpilsJacTimesSetupFn;
typedef IDALsJacTimesVecFn IDASpilsJacTimesVecFn;
typedef IDALsPrecSetupFnB IDASpilsPrecSetupFnB;
typedef IDALsPrecSetupFnBS IDASpilsPrecSetupFnBS;
typedef IDALsPrecSolveFnB IDASpilsPrecSolveFnB;
typedef IDALsPrecSolveFnBS IDASpilsPrecSolveFnBS;
typedef IDALsJacTimesSetupFnB IDASpilsJacTimesSetupFnB;
typedef IDALsJacTimesSetupFnBS IDASpilsJacTimesSetupFnBS;
typedef IDALsJacTimesVecFnB IDASpilsJacTimesVecFnB;
typedef IDALsJacTimesVecFnBS IDASpilsJacTimesVecFnBS;

/*====================================================================
  Exported Functions (wrappers for equivalent routines in idas_ls.h)
  ====================================================================*/

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetLinearSolver instead")
int IDASpilsSetLinearSolver(void *ida_mem, SUNLinearSolver LS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetPreconditioner instead")
int IDASpilsSetPreconditioner(void *ida_mem, IDASpilsPrecSetupFn pset,
                              IDASpilsPrecSolveFn psolve);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetJacTimes instead")
int IDASpilsSetJacTimes(void *ida_mem, IDASpilsJacTimesSetupFn jtsetup,
                        IDASpilsJacTimesVecFn jtimes);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetEpsLin instead")
int IDASpilsSetEpsLin(void *ida_mem, realtype eplifac);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetIncrementFactor instead")
int IDASpilsSetIncrementFactor(void *ida_mem, realtype dqincfac);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetLinWorkSpace instead")
int IDASpilsGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetNumPrecEvals instead")
int IDASpilsGetNumPrecEvals(void *ida_mem, long int *npevals);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetNumPrecSolves instead")
int IDASpilsGetNumPrecSolves(void *ida_mem, long int *npsolves);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetNumLinIters instead")
int IDASpilsGetNumLinIters(void *ida_mem, long int *nliters);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetNumLinConvFails instead")
int IDASpilsGetNumConvFails(void *ida_mem, long int *nlcfails);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetNumJTSetupEvals instead")
int IDASpilsGetNumJTSetupEvals(void *ida_mem, long int *njtsetups);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetNumJtimesEvals instead")
int IDASpilsGetNumJtimesEvals(void *ida_mem, long int *njvevals);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetNumLinResEvals instead")
int IDASpilsGetNumResEvals(void *ida_mem, long int *nrevalsLS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetLastLinFlag instead")
int IDASpilsGetLastFlag(void *ida_mem, long int *flag);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDAGetLinReturnFlagName instead")
char *IDASpilsGetReturnFlagName(long int flag);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetLinearSolverB instead")
int IDASpilsSetLinearSolverB(void *ida_mem, int which, SUNLinearSolver LS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetEpsLinB instead")
int IDASpilsSetEpsLinB(void *ida_mem, int which, realtype eplifacB);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetIncrementFactorB instead")
int IDASpilsSetIncrementFactorB(void *ida_mem, int which, realtype dqincfacB);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetPreconditionerB instead")
int IDASpilsSetPreconditionerB(void *ida_mem, int which,
                               IDASpilsPrecSetupFnB psetB,
                               IDASpilsPrecSolveFnB psolveB);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetPreconditionerBS instead")
int IDASpilsSetPreconditionerBS(void *ida_mem, int which,
                                IDASpilsPrecSetupFnBS psetBS,
                                IDASpilsPrecSolveFnBS psolveBS);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetJacTimesB instead")
int IDASpilsSetJacTimesB(void *ida_mem, int which,
                         IDASpilsJacTimesSetupFnB jtsetupB,
                         IDASpilsJacTimesVecFnB jtimesB);

SUNDIALS_DEPRECATED_EXPORT_MSG("use IDASetJacTimesBS instead")
int IDASpilsSetJacTimesBS(void *ida_mem, int which,
                          IDASpilsJacTimesSetupFnBS jtsetupBS,
                          IDASpilsJacTimesVecFnBS jtimesBS);


#ifdef __cplusplus
}
#endif

#endif
