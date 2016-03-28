/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
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
 * This is the implementation file for an CVSLS linear solver.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "cvode_sparse_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * CVSlsSetSparseJacFn specifies the sparse Jacobian function.
 */
int CVSlsSetSparseJacFn(void *cvode_mem, CVSlsSparseJacFn jac)
{
  CVodeMem cv_mem;
  CVSlsMem cvsls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSLS", "CVSlsSetSparseJacFn", 
		    MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSLS_LMEM_NULL, "CVSLS", 
		    "CVSlsSetSparseJacFn", MSGSP_LMEM_NULL);
    return(CVSLS_LMEM_NULL);
  }
  cvsls_mem = (CVSlsMem) cv_mem->cv_lmem;

  cvsls_mem->s_jaceval = jac;

  return(CVSLS_SUCCESS);
}

/*
 * CVSlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int CVSlsGetNumJacEvals(void *cvode_mem, long int *njevals)
{
  CVodeMem cv_mem;
  CVSlsMem cvsls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSLS", "CVSlsGetNumJacEvals", MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSLS_LMEM_NULL, "CVSLS", 
		    "CVSlsGetNumJacEvals", MSGSP_LMEM_NULL);
    return(CVSLS_LMEM_NULL);
  }
  cvsls_mem = (CVSlsMem) cv_mem->cv_lmem;

  *njevals = cvsls_mem->s_nje;

  return(CVSLS_SUCCESS);
}

/*
 * CVSlsGetReturnFlagName returns the name associated with a CVSLS
 * return value.
 */
char *CVSlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CVSLS_SUCCESS:
    sprintf(name,"CVSLS_SUCCESS");
    break;   
  case CVSLS_MEM_NULL:
    sprintf(name,"CVSLS_MEM_NULL");
    break;
  case CVSLS_LMEM_NULL:
    sprintf(name,"CVSLS_LMEM_NULL");
    break;
  case CVSLS_ILL_INPUT:
    sprintf(name,"CVSLS_ILL_INPUT");
    break;
  case CVSLS_MEM_FAIL:
    sprintf(name,"CVSLS_MEM_FAIL");
    break;
  case CVSLS_JAC_NOSET:
    sprintf(name,"CVSLS_JAC_NOSET");
    break;
  case CVSLS_JACFUNC_UNRECVR:
    sprintf(name,"CVSLS_JACFUNC_UNRECVR");
    break;
  case CVSLS_JACFUNC_RECVR:
    sprintf(name,"CVSLS_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * CVSlsGetLastFlag returns the last flag set in a CVSLS function.
 */
int CVSlsGetLastFlag(void *cvode_mem, long int *flag)
{
  CVodeMem cv_mem;
  CVSlsMem cvsls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSLS", "CVSlsGetLastFlag", 
		    MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSLS_LMEM_NULL, "CVSLS", 
		    "CVSlsGetLastFlag", MSGSP_LMEM_NULL);
    return(CVSLS_LMEM_NULL);
  }
  cvsls_mem = (CVSlsMem) cv_mem->cv_lmem;

  *flag = cvsls_mem->s_last_flag;

  return(CVSLS_SUCCESS);
}

