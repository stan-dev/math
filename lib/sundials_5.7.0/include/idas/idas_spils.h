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

SUNDIALS_EXPORT int IDASpilsSetLinearSolver(void *ida_mem, SUNLinearSolver LS);

SUNDIALS_EXPORT int IDASpilsSetPreconditioner(void *ida_mem, IDASpilsPrecSetupFn pset,
                                              IDASpilsPrecSolveFn psolve);

SUNDIALS_EXPORT int IDASpilsSetJacTimes(void *ida_mem, IDASpilsJacTimesSetupFn jtsetup,
                                        IDASpilsJacTimesVecFn jtimes);

SUNDIALS_EXPORT int IDASpilsSetEpsLin(void *ida_mem, realtype eplifac);

SUNDIALS_EXPORT int IDASpilsSetIncrementFactor(void *ida_mem, realtype dqincfac);

SUNDIALS_EXPORT int IDASpilsGetWorkSpace(void *ida_mem, long int *lenrwLS, long int *leniwLS);

SUNDIALS_EXPORT int IDASpilsGetNumPrecEvals(void *ida_mem, long int *npevals);

SUNDIALS_EXPORT int IDASpilsGetNumPrecSolves(void *ida_mem, long int *npsolves);

SUNDIALS_EXPORT int IDASpilsGetNumLinIters(void *ida_mem, long int *nliters);

SUNDIALS_EXPORT int IDASpilsGetNumConvFails(void *ida_mem, long int *nlcfails);

SUNDIALS_EXPORT int IDASpilsGetNumJTSetupEvals(void *ida_mem, long int *njtsetups);

SUNDIALS_EXPORT int IDASpilsGetNumJtimesEvals(void *ida_mem, long int *njvevals);

SUNDIALS_EXPORT int IDASpilsGetNumResEvals(void *ida_mem, long int *nrevalsLS);

SUNDIALS_EXPORT int IDASpilsGetLastFlag(void *ida_mem, long int *flag);

SUNDIALS_EXPORT char *IDASpilsGetReturnFlagName(long int flag);

SUNDIALS_EXPORT int IDASpilsSetLinearSolverB(void *ida_mem, int which,
                                             SUNLinearSolver LS);

SUNDIALS_EXPORT int IDASpilsSetEpsLinB(void *ida_mem, int which, realtype eplifacB);

SUNDIALS_EXPORT int IDASpilsSetIncrementFactorB(void *ida_mem, int which,
                                                realtype dqincfacB);

SUNDIALS_EXPORT int IDASpilsSetPreconditionerB(void *ida_mem, int which,
                                               IDASpilsPrecSetupFnB psetB,
                                               IDASpilsPrecSolveFnB psolveB);

SUNDIALS_EXPORT int IDASpilsSetPreconditionerBS(void *ida_mem, int which,
                                                IDASpilsPrecSetupFnBS psetBS,
                                                IDASpilsPrecSolveFnBS psolveBS);

SUNDIALS_EXPORT int IDASpilsSetJacTimesB(void *ida_mem, int which,
                                         IDASpilsJacTimesSetupFnB jtsetupB,
                                         IDASpilsJacTimesVecFnB jtimesB);

SUNDIALS_EXPORT int IDASpilsSetJacTimesBS(void *ida_mem, int which,
                                          IDASpilsJacTimesSetupFnBS jtsetupBS,
                                          IDASpilsJacTimesVecFnBS jtimesBS);


#ifdef __cplusplus
}
#endif

#endif
