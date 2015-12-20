/*
 * -----------------------------------------------------------------
 * $Revision: 4357 $
 * $Date: 2015-02-09 13:22:31 -0800 (Mon, 09 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
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
 * This is the header file for the CVKLU linear solver module.
 * -----------------------------------------------------------------
 */

#ifndef _CVKLU_H
#define _CVKLU_H

#include "cvode/cvode_sparse.h"
#include "sundials/sundials_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : CVKLU
 * -----------------------------------------------------------------
 * A call to the CVKLU function links the main integrator      
 * with the CVKLU linear solver module.                        
 *                                                                
 * cv_mem is the pointer to integrator memory returned by        
 *     CVCreate.             
 *
 *                                                                
 * CVKLU returns:                                              
 *     CVSLU_SUCCESS   = 0  if successful                              
 *     CVSLU_LMEM_FAIL = -1 if there was a memory allocation failure   
 *     CVSLU_ILL_INPUT = -2 if NVECTOR found incompatible           
 *                                                                
 * NOTE: The KLU linear solver assumes a serial implementation  
 *       of the NVECTOR package. Therefore, CVKLU will first
 *       test for a compatible N_Vector internal representation
 *       by checking that the functions N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int CVKLU(void *cv_mem, int n, int nnz); 

/*
 * -----------------------------------------------------------------
 * CVKLUReInit
 * -----------------------------------------------------------------
 * This routine reinitializes memory and flags for a new factorization 
 * (symbolic and numeric) to be conducted at the next solver setup
 * call.  This routine is useful in the cases where the number of nonzeroes 
 * has changed or if the structure of the linear system has changed
 * which would require a new symbolic (and numeric factorization).
 *
 * The reinit_type argumenmt governs the level of reinitialization:
 *
 * reinit_type = 1: The Jacobian matrix will be destroyed and 
 *                  a new one will be allocated based on the nnz
 *                  value passed to this call. New symbolic and
 *                  numeric factorizations will be completed at the next
 *                  solver setup.
 *
 * reinit_type = 2: Only symbolic and numeric factorizations will be 
 *                  completed.  It is assumed that the Jacobian size
 *                  has not exceeded the size of nnz given in the prior
 *                  call to CVKLU.
 *
 * This routine assumes no other changes to solver use are necessary.
 *
 * The return value is CVSLS_SUCCESS = 0, CVSLS_MEM_NULL = -1,
 * CVSLS_LMEM_NULL = -2, CVSLS_ILL_INPUT = -3, or CVSLS_MEM_FAIL = -4.
 *
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int CVKLUReInit(void *cv_mem_v, int n, int nnz, 
				  int reinit_type);

/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * CVKLUSetOrdering sets the ordering used by KLU for reducing fill.
 * Options are: 0 for AMD, 1 for COLAMD, and 2 for the natural ordering.
 * The default used in CVODE is 1 for COLAMD.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int CVKLUSetOrdering(void *cv_mem, int ordering_choice); 


  
#ifdef __cplusplus
}
#endif

#endif
