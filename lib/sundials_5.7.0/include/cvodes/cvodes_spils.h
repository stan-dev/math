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

SUNDIALS_EXPORT int CVSpilsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS);

SUNDIALS_EXPORT int CVSpilsSetEpsLin(void *cvode_mem, realtype eplifac);

SUNDIALS_EXPORT int CVSpilsSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn pset,
                             CVSpilsPrecSolveFn psolve);

SUNDIALS_EXPORT int CVSpilsSetJacTimes(void *cvode_mem, CVSpilsJacTimesSetupFn jtsetup,
                       CVSpilsJacTimesVecFn jtimes);

SUNDIALS_EXPORT int CVSpilsGetWorkSpace(void *cvode_mem, long int *lenrwLS,
                        long int *leniwLS);

SUNDIALS_EXPORT int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals);

SUNDIALS_EXPORT int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves);

SUNDIALS_EXPORT int CVSpilsGetNumLinIters(void *cvode_mem, long int *nliters);

SUNDIALS_EXPORT int CVSpilsGetNumConvFails(void *cvode_mem, long int *nlcfails);

SUNDIALS_EXPORT int CVSpilsGetNumJTSetupEvals(void *cvode_mem, long int *njtsetups);

SUNDIALS_EXPORT int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals);

SUNDIALS_EXPORT int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS);

SUNDIALS_EXPORT int CVSpilsGetLastFlag(void *cvode_mem, long int *flag);

SUNDIALS_EXPORT char *CVSpilsGetReturnFlagName(long int flag);

SUNDIALS_EXPORT int CVSpilsSetLinearSolverB(void *cvode_mem, int which,
                            SUNLinearSolver LS);

SUNDIALS_EXPORT int CVSpilsSetEpsLinB(void *cvode_mem, int which, realtype eplifacB);

SUNDIALS_EXPORT int CVSpilsSetPreconditionerB(void *cvode_mem, int which,
                              CVSpilsPrecSetupFnB psetB,
                              CVSpilsPrecSolveFnB psolveB);

SUNDIALS_EXPORT int CVSpilsSetPreconditionerBS(void *cvode_mem, int which,
                               CVSpilsPrecSetupFnBS psetBS,
                               CVSpilsPrecSolveFnBS psolveBS);

SUNDIALS_EXPORT int CVSpilsSetJacTimesB(void *cvode_mem, int which,
                        CVSpilsJacTimesSetupFnB jtsetupB,
                        CVSpilsJacTimesVecFnB jtimesB);

SUNDIALS_EXPORT int CVSpilsSetJacTimesBS(void *cvode_mem, int which,
                         CVSpilsJacTimesSetupFnBS jtsetupBS,
                         CVSpilsJacTimesVecFnBS jtimesBS);


#ifdef __cplusplus
}
#endif

#endif
