/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Scott Cohen, Alan Hindmarsh, Radu Serban, 
 *                  and Aaron Collier @ LLNL
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
 * Header file for the deprecated Scaled Preconditioned Iterative 
 * Linear Solver interface in KINSOL; these routines now just wrap 
 * the updated KINSOL generic linear solver interface in kinsol_ls.h.
 *-----------------------------------------------------------------*/

#ifndef _KINSPILS_H
#define _KINSPILS_H

#include <kinsol/kinsol_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  Function Types (typedefs for equivalent types in kinsol_ls.h)
  =================================================================*/

typedef KINLsPrecSetupFn KINSpilsPrecSetupFn;
typedef KINLsPrecSolveFn KINSpilsPrecSolveFn;
typedef KINLsJacTimesVecFn KINSpilsJacTimesVecFn;

/*=================================================================
  Exported Functions (wrappers for equivalent routines in kinsol_ls.h)
  =================================================================*/

int KINSpilsSetLinearSolver(void *kinmem, SUNLinearSolver LS)
{ return(KINSetLinearSolver(kinmem, LS, NULL)); }

int KINSpilsSetPreconditioner(void *kinmem, KINSpilsPrecSetupFn psetup,
                              KINSpilsPrecSolveFn psolve)
{ return(KINSetPreconditioner(kinmem, psetup, psolve)); }
  
int KINSpilsSetJacTimesVecFn(void *kinmem, KINSpilsJacTimesVecFn jtv)
{ return(KINSetJacTimesVecFn(kinmem, jtv)); }

int KINSpilsGetWorkSpace(void *kinmem, long int *lenrwLS, long int *leniwLS)
{ return(KINGetLinWorkSpace(kinmem, lenrwLS, leniwLS)); }
  
int KINSpilsGetNumPrecEvals(void *kinmem, long int *npevals)
{ return(KINGetNumPrecEvals(kinmem, npevals)); }
  
int KINSpilsGetNumPrecSolves(void *kinmem, long int *npsolves)
{ return(KINGetNumPrecSolves(kinmem, npsolves)); }
  
int KINSpilsGetNumLinIters(void *kinmem, long int *nliters)
{ return(KINGetNumLinIters(kinmem, nliters)); }
  
int KINSpilsGetNumConvFails(void *kinmem, long int *nlcfails)
{ return(KINGetNumLinConvFails(kinmem, nlcfails)); }
  
int KINSpilsGetNumJtimesEvals(void *kinmem, long int *njvevals)
{ return(KINGetNumJtimesEvals(kinmem, njvevals)); }
  
int KINSpilsGetNumFuncEvals(void *kinmem, long int *nfevals)
{ return(KINGetNumLinFuncEvals(kinmem, nfevals)); }
  
int KINSpilsGetLastFlag(void *kinmem, long int *flag)
{ return(KINGetLastLinFlag(kinmem, flag)); }
  
char *KINSpilsGetReturnFlagName(long int flag)
{ return(KINGetLinReturnFlagName(flag)); }


#ifdef __cplusplus
}
#endif

#endif
