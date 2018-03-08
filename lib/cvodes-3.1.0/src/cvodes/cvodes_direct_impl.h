/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 * -----------------------------------------------------------------
 * Implementation header file for the CVSDLS linear solver interface
 * -----------------------------------------------------------------
 */

#ifndef _CVSDLS_IMPL_H
#define _CVSDLS_IMPL_H

#include <cvodes/cvodes_direct.h>
#include "cvodes_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  CVSDLS Internal Constants
  =================================================================*/

/*-----------------------------------------------------------------
  CVDLS solver constants
  -----------------------------------------------------------------
  CVD_MSBJ   maximum number of steps between Jacobian evaluations
  CVD_DGMAX  maximum change in gamma between Jacobian evaluations
  -----------------------------------------------------------------*/

#define CVD_MSBJ  50
#define CVD_DGMAX RCONST(0.2)

  
/*=================================================================
  PART I:  Forward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  Types: CVDlsMemRec, CVDlsMem                             
  -----------------------------------------------------------------
  CVDlsMem is pointer to a CVDlsMemRec structure.
  -----------------------------------------------------------------*/

typedef struct CVDlsMemRec {

  booleantype jacDQ;    /* SUNTRUE if using internal DQ Jacobian approx. */
  CVDlsJacFn jac;       /* Jacobian routine to be called                 */
  void *J_data;         /* data pointer passed to jac                    */

  SUNMatrix A;          /* A = I - gamma * df/dy                         */
  SUNMatrix savedJ;     /* savedJ = old Jacobian                         */

  SUNLinearSolver LS;   /* generic direct linear solver object           */

  N_Vector x;           /* solution vector used by SUNLinearSolver       */
  
  long int nstlj;       /* nstlj = nst at last Jacobian eval.            */

  long int nje;         /* nje = no. of calls to jac                     */

  long int nfeDQ;       /* no. of calls to f due to DQ Jacobian approx.  */

  long int last_flag;   /* last error return flag                        */

} *CVDlsMem;

/*-----------------------------------------------------------------
  Prototypes of internal functions
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  PRIVATE FUNCTIONS 
  -----------------------------------------------------------------*/

/* difference-quotient Jacobian approximation routines */
int cvDlsDQJac(realtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
               void *data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
int cvDlsDenseDQJac(realtype t, N_Vector y, N_Vector fy, 
                    SUNMatrix Jac, CVodeMem cv_mem, N_Vector tmp1);
  
int cvDlsBandDQJac(realtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                   CVodeMem cv_mem, N_Vector tmp1, N_Vector tmp2);
  
/* generic linit/lsetup/lsolve/lfree interface routines for CVode to call */
int cvDlsInitialize(CVodeMem cv_mem);

int cvDlsSetup(CVodeMem cv_mem, int convfail, N_Vector y, 
               N_Vector fy, booleantype *jcurPtr, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3); 
  
int cvDlsSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
               N_Vector ycur, N_Vector fcur);

int cvDlsFree(CVodeMem cv_mem);

/* Auxilliary functions */
int cvDlsInitializeCounters(CVDlsMem cvdls_mem);


/*=================================================================
  PART II:  Backward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  Types : CVDlsMemRecB, CVDlsMemB       
  -----------------------------------------------------------------
  A CVDLS linear solver's specification function attaches such
  a structure to the lmemB filed of CVodeBMem
  -----------------------------------------------------------------*/

typedef struct CVDlsMemRecB {

  CVDlsJacFnB jacB;
  CVDlsJacFnBS jacBS;

} *CVDlsMemB;

/*-----------------------------------------------------------------
  Prototypes of internal functions
  -----------------------------------------------------------------*/

int cvDlsFreeB(CVodeBMem cvb_mem);

/*=================================================================
  Error Messages
  =================================================================*/

#define MSGD_CVMEM_NULL "Integrator memory is NULL."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_MEM_FAIL "A memory request failed."
#define MSGD_LMEM_NULL "Linear solver memory is NULL."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."
#define MSGD_MATCOPY_FAILED "The SUNMatCopy routine failed in an unrecoverable manner."
#define MSGD_MATZERO_FAILED "The SUNMatZero routine failed in an unrecoverable manner."
#define MSGD_MATSCALEADDI_FAILED "The SUNMatScaleAddI routine failed in an unrecoverable manner."

#define MSGD_NO_ADJ "Illegal attempt to call before calling CVodeAdjMalloc."
#define MSGD_BAD_WHICH "Illegal value for which."
#define MSGD_LMEMB_NULL "Linear solver memory is NULL for the backward integration."
#define MSGD_BAD_TINTERP "Bad t for interpolation."

#ifdef __cplusplus
}
#endif

#endif
