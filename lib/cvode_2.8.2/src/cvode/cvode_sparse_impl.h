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
 * Implementation header file for the CVSLS linear solvers.
 * -----------------------------------------------------------------
 */

#ifndef _CVSPARSE_IMPL_H
#define _CVSPARSE_IMPL_H

#include "cvode/cvode_sparse.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 * C V S P A R S E    I N T E R N A L    C O N S T A N T S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * CVSLS solver constants
 * -----------------------------------------------------------------
 * CVS_MSBJ   maximum number of steps between Jacobian evaluations
 * CVS_DGMAX  maximum change in gamma between Jacobian evaluations
 * -----------------------------------------------------------------
 */
 
#define CVS_MSBJ  50
#define CVS_DGMAX RCONST(0.2)
 
/*
 * -----------------------------------------------------------------
 * Types : CVSlsMemRec, CVSlsMem                             
 * -----------------------------------------------------------------
 * CVSlsMem is pointer to a CVSlsMemRec structure.
 * -----------------------------------------------------------------
 */

typedef struct CVSlsMemRec {

  CVSlsSparseJacFn s_jaceval; /* user Jacobian evaluation routine 
				  to be called */
  void *s_jacdata;           /* J_data is passed to djac or bjac */

  long int s_nje;           /* nje = no. of calls to jac */

  long int s_last_flag;     /* last error return flag */

  int s_first_factorize;    /* flag telling whether the first 
			       factorization needs to happen */

  int s_nstlj;              /* time step of last Jacobian evaluation */

  SlsMat s_JacMat;          /* M = I - gamma * df/dy */

  SlsMat s_savedJ;          /* saved copy of Jacobian */

  void *s_solver_data;      /* structure for solver-specific data */
  

} *CVSlsMem;

/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */
  

/*
 * =================================================================
 * E R R O R   M E S S A G E S
 * =================================================================
 */

#define MSGSP_CVMEM_NULL "Integrator memory is NULL."
#define MSGSP_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGSP_MEM_FAIL "A memory request failed."
#define MSGSP_LMEM_NULL "Linear solver memory is NULL."
#define MSGSP_ILL_INPUT "Invalid input detected."
#define MSGSP_JAC_NOSET "Jacobian evaluation function has not been set."
#define MSGSP_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."
#define MSGSP_PACKAGE_FAIL "A call to an external package failed."

#ifdef __cplusplus
}
#endif

#endif
