/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Radu Serban @ LLNL
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
 * Header file for the deprecated direct linear solver interface in 
 * KINSOL; these routines now just wrap the updated KINSOL generic
 * linear solver interface in kinsol_ls.h.
 *-----------------------------------------------------------------*/

#ifndef _KINDLS_H
#define _KINDLS_H

#include <kinsol/kinsol_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  Function Types (typedefs for equivalent types in kinsol_ls.h)
  =================================================================*/

typedef KINLsJacFn KINDlsJacFn;
  
/*=================================================================
  Exported Functions (wrappers for equivalent routines in kinsol_ls.h)
  =================================================================*/
  
int KINDlsSetLinearSolver(void *kinmem, SUNLinearSolver LS, SUNMatrix A)
{ return(KINSetLinearSolver(kinmem, LS, A)); }

int KINDlsSetJacFn(void *kinmem, KINDlsJacFn jac)
{ return(KINSetJacFn(kinmem, jac)); }

int KINDlsGetWorkSpace(void *kinmem, long int *lenrw, long int *leniw)
{ return(KINGetLinWorkSpace(kinmem, lenrw, leniw)); }

int KINDlsGetNumJacEvals(void *kinmem, long int *njevals)
{ return(KINGetNumJacEvals(kinmem, njevals)); }

int KINDlsGetNumFuncEvals(void *kinmem, long int *nfevals)
{ return(KINGetNumLinFuncEvals(kinmem, nfevals)); }

int KINDlsGetLastFlag(void *kinmem, long int *flag)
{ return(KINGetLastLinFlag(kinmem, flag)); }

char *KINDlsGetReturnFlagName(long int flag)
{ return(KINGetLinReturnFlagName(flag)); }

#ifdef __cplusplus
}
#endif

#endif
