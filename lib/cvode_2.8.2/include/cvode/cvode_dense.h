/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
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
 * This is the header file for the CVODE dense linear solver, CVDENSE.
 * -----------------------------------------------------------------
 */

#ifndef _CVDENSE_H
#define _CVDENSE_H

#include <cvode/cvode_direct.h>
#include <sundials/sundials_dense.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function: CVDense
 * -----------------------------------------------------------------
 * A call to the CVDense function links the main integrator with
 * the CVDENSE linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * N is the size of the ODE system.
 *
 * The return value of CVDense is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the cvode memory was NULL
 *    CVDLS_MEM_FAIL  if there was a memory allocation failure
 *    CVDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDense(void *cvode_mem, long int N);

#ifdef __cplusplus
}
#endif

#endif
