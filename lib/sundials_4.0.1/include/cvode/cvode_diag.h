/*
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the CVODE diagonal linear solver, CVDIAG.
 * -----------------------------------------------------------------
 */

#ifndef _CVDIAG_H
#define _CVDIAG_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : CVDiag
 * -----------------------------------------------------------------
 * A call to the CVDiag function links the main integrator with
 * the CVDIAG linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * The return value of CVDiag is one of:
 *    CVDIAG_SUCCESS   if successful
 *    CVDIAG_MEM_NULL  if the cvode memory was NULL
 *    CVDIAG_MEM_FAIL  if there was a memory allocation failure
 *    CVDIAG_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDiag(void *cvode_mem);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVDIAG linear solver
 * -----------------------------------------------------------------
 *
 * CVDiagGetWorkSpace returns the real and integer workspace used
 *                    by CVDIAG.
 * CVDiagGetNumRhsEvals returns the number of calls to the user
 *                      f routine due to finite difference Jacobian
 *                      evaluation.
 *                      Note: The number of diagonal approximate
 *                      Jacobians formed is equal to the number of
 *                      CVDiagSetup calls. This number is available
 *                      through CVodeGetNumLinSolvSetups.
 * CVDiagGetLastFlag returns the last error flag set by any of
 *                   the CVDIAG interface functions.
 *
 * The return value of CVDiagGet* is one of:
 *    CVDIAG_SUCCESS   if successful
 *    CVDIAG_MEM_NULL  if the cvode memory was NULL
 *    CVDIAG_LMEM_NULL if the cvdiag memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDiagGetWorkSpace(void *cvode_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int CVDiagGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS);
SUNDIALS_EXPORT int CVDiagGetLastFlag(void *cvode_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a CVDIAG return flag
 * -----------------------------------------------------------------
 */
  
SUNDIALS_EXPORT char *CVDiagGetReturnFlagName(long int flag);

/*
 * -----------------------------------------------------------------
 * CVDIAG return values
 * -----------------------------------------------------------------
 */

#define CVDIAG_SUCCESS          0
#define CVDIAG_MEM_NULL        -1
#define CVDIAG_LMEM_NULL       -2
#define CVDIAG_ILL_INPUT       -3
#define CVDIAG_MEM_FAIL        -4

/* Additional last_flag values */

#define CVDIAG_INV_FAIL        -5
#define CVDIAG_RHSFUNC_UNRECVR -6
#define CVDIAG_RHSFUNC_RECVR   -7

#ifdef __cplusplus
}
#endif

#endif
