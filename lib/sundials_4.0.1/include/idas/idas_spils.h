/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *-----------------------------------------------------------------
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
 *-----------------------------------------------------------------
 * Header file for the deprecated Scaled and Preconditioned 
 * Iterative Linear Solver interface in IDAS; these routines now just
 * wrap the updated IDA generic linear solver interface in idas_ls.h.
 *-----------------------------------------------------------------*/

#ifndef _IDASSPILS_H
#define _IDASSPILS_H

#include <idas/idas_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  Function Types (typedefs for equivalent types in idas_ls.h)
  =================================================================*/
  
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


/*=================================================================
  Exported Functions (wrappers for equivalent routines in idas_ls.h)
  =================================================================*/
  
int IDASpilsSetLinearSolver(void *ida_mem, SUNLinearSolver LS)
{ return(IDASetLinearSolver(ida_mem, LS, NULL)); }
  
int IDASpilsSetPreconditioner(void *ida_mem, IDASpilsPrecSetupFn pset, 
                              IDASpilsPrecSolveFn psolve)
{ return(IDASetPreconditioner(ida_mem, pset,  psolve)); }
  
int IDASpilsSetJacTimes(void *ida_mem, IDASpilsJacTimesSetupFn jtsetup,
                        IDASpilsJacTimesVecFn jtimes)
{ return(IDASetJacTimes(ida_mem, jtsetup, jtimes)); }
  
int IDASpilsSetEpsLin(void *ida_mem, realtype eplifac)
{ return(IDASetEpsLin(ida_mem, eplifac)); }
  
int IDASpilsSetIncrementFactor(void *ida_mem, realtype dqincfac)
{ return(IDASetIncrementFactor(ida_mem, dqincfac)); }
  
int IDASpilsGetWorkSpace(void *ida_mem, long int *lenrwLS,
                         long int *leniwLS)
{ return(IDAGetLinWorkSpace(ida_mem, lenrwLS, leniwLS)); }
  
int IDASpilsGetNumPrecEvals(void *ida_mem, long int *npevals)
{ return(IDAGetNumPrecEvals(ida_mem, npevals)); }
  
int IDASpilsGetNumPrecSolves(void *ida_mem, long int *npsolves)
{ return(IDAGetNumPrecSolves(ida_mem, npsolves)); }
  
int IDASpilsGetNumLinIters(void *ida_mem, long int *nliters)
{ return(IDAGetNumLinIters(ida_mem, nliters)); }
  
int IDASpilsGetNumConvFails(void *ida_mem, long int *nlcfails)
{ return(IDAGetNumLinConvFails(ida_mem, nlcfails)); }
  
int IDASpilsGetNumJTSetupEvals(void *ida_mem, long int *njtsetups)
{ return(IDAGetNumJTSetupEvals(ida_mem, njtsetups)); }
  
int IDASpilsGetNumJtimesEvals(void *ida_mem, long int *njvevals)
{ return(IDAGetNumJtimesEvals(ida_mem, njvevals)); }
  
int IDASpilsGetNumResEvals(void *ida_mem, long int *nrevalsLS)
{ return(IDAGetNumLinResEvals(ida_mem, nrevalsLS)); }
  
int IDASpilsGetLastFlag(void *ida_mem, long int *flag)
{ return(IDAGetLastLinFlag(ida_mem, flag)); }
  
char *IDASpilsGetReturnFlagName(long int flag)
{ return(IDAGetLinReturnFlagName(flag)); }
  
int IDASpilsSetLinearSolverB(void *ida_mem, int which,
                             SUNLinearSolver LS)
{ return(IDASetLinearSolverB(ida_mem, which, LS, NULL)); }
  
int IDASpilsSetEpsLinB(void *ida_mem, int which, realtype eplifacB)
{ return(IDASetEpsLinB(ida_mem, which, eplifacB)); }
  
int IDASpilsSetIncrementFactorB(void *ida_mem, int which,
                                realtype dqincfacB)
{ return(IDASetIncrementFactorB(ida_mem, which, dqincfacB)); }
  
int IDASpilsSetPreconditionerB(void *ida_mem, int which,
                               IDASpilsPrecSetupFnB psetB,
                               IDASpilsPrecSolveFnB psolveB)
{ return(IDASetPreconditionerB(ida_mem, which, psetB, psolveB)); }
  
int IDASpilsSetPreconditionerBS(void *ida_mem, int which,
                                IDASpilsPrecSetupFnBS psetBS,
                                IDASpilsPrecSolveFnBS psolveBS)
{ return(IDASetPreconditionerBS(ida_mem, which, psetBS, psolveBS)); }
  
int IDASpilsSetJacTimesB(void *ida_mem, int which,
                         IDASpilsJacTimesSetupFnB jtsetupB,
                         IDASpilsJacTimesVecFnB jtimesB)
{ return(IDASetJacTimesB(ida_mem, which, jtsetupB, jtimesB)); }
  
int IDASpilsSetJacTimesBS(void *ida_mem, int which,
                          IDASpilsJacTimesSetupFnBS jtsetupBS,
                          IDASpilsJacTimesVecFnBS jtimesBS)
{ return(IDASetJacTimesBS(ida_mem, which, jtsetupBS, jtimesBS)); }
  

 
#ifdef __cplusplus
}
#endif

#endif
