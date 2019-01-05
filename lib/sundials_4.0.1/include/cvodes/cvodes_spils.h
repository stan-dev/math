/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
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
  CVSSPILS user-supplied function prototypes (typedefs for 
  equivalent types in cvodes_ls.h)
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

/*=================================================================
  CVSSPILS Exported functions (wrappers for equivalent routines in 
  cvodes_ls.h)
  =================================================================*/

int CVSpilsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS)
{ return(CVodeSetLinearSolver(cvode_mem, LS, NULL)); }

int CVSpilsSetEpsLin(void *cvode_mem, realtype eplifac)
{ return(CVodeSetEpsLin(cvode_mem, eplifac)); }

int CVSpilsSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn pset, 
                             CVSpilsPrecSolveFn psolve)
{ return(CVodeSetPreconditioner(cvode_mem, pset, psolve)); }

int CVSpilsSetJacTimes(void *cvode_mem, CVSpilsJacTimesSetupFn jtsetup,
                       CVSpilsJacTimesVecFn jtimes)
{ return(CVodeSetJacTimes(cvode_mem, jtsetup, jtimes)); }

int CVSpilsGetWorkSpace(void *cvode_mem, long int *lenrwLS,
                        long int *leniwLS)
{ return(CVodeGetLinWorkSpace(cvode_mem, lenrwLS, leniwLS)); }
  
int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals)
{ return(CVodeGetNumPrecEvals(cvode_mem, npevals)); }
  
int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves)
{ return(CVodeGetNumPrecSolves(cvode_mem, npsolves)); }
  
int CVSpilsGetNumLinIters(void *cvode_mem, long int *nliters)
{ return(CVodeGetNumLinIters(cvode_mem, nliters)); }
  
int CVSpilsGetNumConvFails(void *cvode_mem, long int *nlcfails)
{ return(CVodeGetNumLinConvFails(cvode_mem, nlcfails)); }
  
int CVSpilsGetNumJTSetupEvals(void *cvode_mem, long int *njtsetups)
{ return(CVodeGetNumJTSetupEvals(cvode_mem, njtsetups)); }
  
int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals)
{ return(CVodeGetNumJtimesEvals(cvode_mem, njvevals)); }
  
int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{ return(CVodeGetNumLinRhsEvals(cvode_mem, nfevalsLS)); }
  
int CVSpilsGetLastFlag(void *cvode_mem, long int *flag)
{ return(CVodeGetLastLinFlag(cvode_mem, flag)); }
  
char *CVSpilsGetReturnFlagName(long int flag)
{ return(CVodeGetLinReturnFlagName(flag)); }
  
int CVSpilsSetLinearSolverB(void *cvode_mem, int which,
                            SUNLinearSolver LS)
{ return(CVodeSetLinearSolverB(cvode_mem, which, LS, NULL)); }
  
int CVSpilsSetEpsLinB(void *cvode_mem, int which, realtype eplifacB)
{ return(CVodeSetEpsLinB(cvode_mem, which, eplifacB)); }
  
int CVSpilsSetPreconditionerB(void *cvode_mem, int which, 
                              CVSpilsPrecSetupFnB psetB,
                              CVSpilsPrecSolveFnB psolveB)
{ return(CVodeSetPreconditionerB(cvode_mem, which, psetB, psolveB)); }
  
int CVSpilsSetPreconditionerBS(void *cvode_mem, int which, 
                               CVSpilsPrecSetupFnBS psetBS,
                               CVSpilsPrecSolveFnBS psolveBS)
{ return(CVodeSetPreconditionerBS(cvode_mem, which, psetBS, psolveBS)); }
  
int CVSpilsSetJacTimesB(void *cvode_mem, int which, 
                        CVSpilsJacTimesSetupFnB jtsetupB,
                        CVSpilsJacTimesVecFnB jtimesB)
{ return(CVodeSetJacTimesB(cvode_mem, which, jtsetupB, jtimesB)); }
  
int CVSpilsSetJacTimesBS(void *cvode_mem, int which, 
                         CVSpilsJacTimesSetupFnBS jtsetupBS,
                         CVSpilsJacTimesVecFnBS jtimesBS)
{ return(CVodeSetJacTimesBS(cvode_mem, which, jtsetupBS, jtimesBS)); }
  

#ifdef __cplusplus
}
#endif

#endif
