/*
 * -----------------------------------------------------------------
 * $Rev $
 * $Date $
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
 * This is the header file for the CVSuperLUMT linear solver module.
 * -----------------------------------------------------------------
 */

#ifndef _CVSSUPERLUMT_H
#define _CVSSUPERLUMT_H

#include "cvodes/cvodes_sparse.h"
#include "sundials/sundials_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : CVSuperLUMT
 * -----------------------------------------------------------------
 * A call to the CVSuperLUMT function links the main integrator      
 * with the CVSuperLUMT linear solver module.                        
 *                                                                
 * cv_mem is the pointer to integrator memory returned by        
 *     CVCreate.             
 *
 *                                                                
 * CVSuperLUMT returns:                                              
 *     CVSLU_SUCCESS   = 0  if successful                              
 *     CVSLU_LMEM_FAIL = -1 if there was a memory allocation failure   
 *     CVSLU_ILL_INPUT = -2 if NVECTOR found incompatible           
 *                                                                
 * NOTE: The CVSuperLUMT linear solver assumes a serial implementation  
 *       of the NVECTOR package. Therefore, CVSuperLUMT will first
 *       test for a compatible N_Vector internal representation
 *       by checking that the functions N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int CVSuperLUMT(void *cv_mem, int num_threads,
				  int n, int nnz); 

/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * CVSuperLUMTSetOrdering sets the ordering used by CVSuperLUMT for 
 * reducing fill.
 * Options are: 
 * 0 for natural ordering
 * 1 for minimal degree ordering on A'*A
 * 2 for minimal degree ordering on A'+A
 * 3 for approximate minimal degree ordering for unsymmetric matrices
 * The default used in SUNDIALS is 3 for COLAMD.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int CVSuperLUMTSetOrdering(void *cv_mem, 
					     int ordering_choice); 


/*
 * -----------------------------------------------------------------
 * Function: CVSuperLUMTB
 * -----------------------------------------------------------------
 * CVSuperLUMTB links the main CVODES integrator with the CVSuperLUMT
 * linear solver for the backward integration.
 * The 'which' argument is the int returned by CVodeCreateB.
 * -----------------------------------------------------------------
 */
 
  SUNDIALS_EXPORT int CVSuperLUMTB(void *cvode_mem, int which, int num_threads, 
				   int nB, int nnzB);
 
 
/*
 * -----------------------------------------------------------------
 * Function: CVSuperLUMTSetOrderingB
 * -----------------------------------------------------------------
 * CVSuperLUMTSetOrderingB pulls off the memory block associated with the
 * which parameter and sets the ordering for the solver associated with that block.
 * The 'which' argument is the int returned by CVodeCreateB.
 * -----------------------------------------------------------------
 */
 
  SUNDIALS_EXPORT int CVSuperLUMTSetOrderingB(void *cvode_mem, int which, 
					      int ordering_choice);

  
#ifdef __cplusplus
}
#endif

#endif
