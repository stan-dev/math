/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
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
 * Header file for the CVODE dense linear solver CVLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _CVLAPACK_H
#define _CVLAPACK_H

#include <cvode/cvode_direct.h>
#include <sundials/sundials_lapack.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : CVLapackDense
 * -----------------------------------------------------------------
 * A call to the CVLapackDense function links the main integrator
 * with the CVLAPACK linear solver using dense Jacobians.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * N is the size of the ODE system.
 *
 * The return value of CVLapackDense is one of:
 *    CVLAPACK_SUCCESS   if successful
 *    CVLAPACK_MEM_NULL  if the CVODE memory was NULL
 *    CVLAPACK_MEM_FAIL  if there was a memory allocation failure
 *    CVLAPACK_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVLapackDense(void *cvode_mem, int N);

/*
 * -----------------------------------------------------------------
 * Function : CVLapackBand
 * -----------------------------------------------------------------
 * A call to the CVLapackBand function links the main integrator
 * with the CVLAPACK linear solver using banded Jacobians. 
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * N is the size of the ODE system.
 *
 * mupper is the upper bandwidth of the band Jacobian approximation.
 *
 * mlower is the lower bandwidth of the band Jacobian approximation.
 *
 * The return value of CVLapackBand is one of:
 *    CVLAPACK_SUCCESS   if successful
 *    CVLAPACK_MEM_NULL  if the CVODE memory was NULL
 *    CVLAPACK_MEM_FAIL  if there was a memory allocation failure
 *    CVLAPACK_ILL_INPUT if a required vector operation is missing or
 *                       if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVLapackBand(void *cvode_mem, int N, int mupper, int mlower);

#ifdef __cplusplus
}
#endif

#endif
