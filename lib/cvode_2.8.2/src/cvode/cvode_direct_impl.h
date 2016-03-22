/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
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
 * Common implementation header file for the CVDLS linear solvers.
 * -----------------------------------------------------------------
 */

#ifndef _CVDLS_IMPL_H
#define _CVDLS_IMPL_H

#include "cvode/cvode_direct.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * CVDLS solver constants
 * -----------------------------------------------------------------
 * CVD_MSBJ   maximum number of steps between Jacobian evaluations
 * CVD_DGMAX  maximum change in gamma between Jacobian evaluations
 * -----------------------------------------------------------------
 */

#define CVD_MSBJ  50
#define CVD_DGMAX RCONST(0.2)

/*
 * -----------------------------------------------------------------
 * Types : CVDlsMemRec, CVDlsMem                             
 * -----------------------------------------------------------------
 * CVDlsMem is pointer to a CVDlsMemRec structure.
 * -----------------------------------------------------------------
 */

typedef struct CVDlsMemRec {

  int d_type;             /* SUNDIALS_DENSE or SUNDIALS_BAND              */

  long int d_n;           /* problem dimension                            */

  long int d_ml;          /* lower bandwidth of Jacobian                  */
  long int d_mu;          /* upper bandwidth of Jacobian                  */ 
  long int d_smu;         /* upper bandwith of M = MIN(N-1,d_mu+d_ml)     */

  booleantype d_jacDQ;    /* TRUE if using internal DQ Jacobian approx.   */
  CVDlsDenseJacFn d_djac; /* dense Jacobian routine to be called          */
  CVDlsBandJacFn d_bjac;  /* band Jacobian routine to be called           */
  void *d_J_data;         /* user data is passed to djac or bjac          */

  DlsMat d_M;             /* M = I - gamma * df/dy                        */
  DlsMat d_savedJ;        /* savedJ = old Jacobian                        */

  int *d_pivots;          /* pivots = int pivot array for PM = LU         */
  long int *d_lpivots;    /* lpivots = long int pivot array for PM = LU   */

  long int  d_nstlj;      /* nstlj = nst at last Jacobian eval.           */

  long int d_nje;         /* nje = no. of calls to jac                    */

  long int d_nfeDQ;       /* no. of calls to f due to DQ Jacobian approx. */

  long int d_last_flag;   /* last error return flag                       */
  
} *CVDlsMem;

/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

int cvDlsDenseDQJac(long int N, realtype t,
		    N_Vector y, N_Vector fy, 
		    DlsMat Jac, void *data,
		    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
int cvDlsBandDQJac(long int N, long int mupper, long int mlower,
		   realtype t, N_Vector y, N_Vector fy, 
		   DlsMat Jac, void *data,
		   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


/*
 * -----------------------------------------------------------------
 * Error Messages
 * -----------------------------------------------------------------
 */

#define MSGD_CVMEM_NULL "Integrator memory is NULL."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_MEM_FAIL "A memory request failed."
#define MSGD_LMEM_NULL "Linear solver memory is NULL."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
