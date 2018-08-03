/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Radu Serban @ LLNL
 *-----------------------------------------------------------------
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
 *-----------------------------------------------------------------
 * Implementation header file for the IDADLS linear solver 
 * interface
 *-----------------------------------------------------------------*/

#ifndef _IDASDLS_IMPL_H
#define _IDASDLS_IMPL_H

#include <idas/idas_direct.h>
#include "idas_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  PART I:  Forward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  Types : IDADlsMemRec, IDADlsMem                             

  IDADlsMem is pointer to a IDADlsMemRec structure.
  -----------------------------------------------------------------*/
typedef struct IDADlsMemRec {

  booleantype jacDQ;    /* SUNTRUE if using internal DQ Jacobian approx.  */
  IDADlsJacFn jac;      /* dense Jacobian routine to be called            */
  void *J_data;         /* J_data is passed to jac                        */

  SUNLinearSolver LS;   /* generic direct linear solver object            */

  SUNMatrix J;          /* J = dF/dy + cj*dF/dy'                          */

  N_Vector x;           /* solution vector used by SUNLinearSolver        */
  
  long int nje;         /* nje = no. of calls to jac                      */

  long int nreDQ;       /* no. of calls to res due to DQ Jacobian approx. */

  long int last_flag;   /* last error return flag                         */
  
} *IDADlsMem;

/*---------------------------------------------------------------
  Prototypes of internal functions
  ---------------------------------------------------------------*/
  
/* difference-quotient Jacobian approximation routines */
int idaDlsDQJac(realtype tt, realtype c_j, N_Vector yy, 
                N_Vector yp, N_Vector rr, SUNMatrix Jac, 
                void *data, N_Vector tmp1, N_Vector tmp2, 
                N_Vector tmp3);
int idaDlsDenseDQJac(realtype tt, realtype c_j, N_Vector yy, 
                     N_Vector yp, N_Vector rr, SUNMatrix Jac,
                     IDAMem IDA_mem, N_Vector tmp1);
 int idaDlsBandDQJac(realtype tt, realtype c_j, N_Vector yy,
                     N_Vector yp, N_Vector rr, SUNMatrix Jac,
                    IDAMem IDA_mem, N_Vector tmp1,
                    N_Vector tmp2, N_Vector tmp3);

/* generic linit/lsetup/lsolve/lfree interface routines for IDA to call */
int idaDlsInitialize(IDAMem IDA_mem);

int idaDlsSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
                N_Vector resp, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3); 

int idaDlsSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                N_Vector ycur, N_Vector ypcur, N_Vector rescur);

int idaDlsFree(IDAMem IDA_mem);

/* Auxilliary functions */
int idaDlsInitializeCounters(IDADlsMem idadls_mem);


  
/*=================================================================
  PART II:  Backward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  Types : IDADlsMemRecB, IDADlsMemB       
  -----------------------------------------------------------------
  An IDADLS linear solver's specification function attaches such
  a structure to the lmemB filed of IDABMem
  -----------------------------------------------------------------*/
typedef struct IDADlsMemRecB {

  IDADlsJacFnB jacB;
  IDADlsJacFnBS jacBS;
  
} *IDADlsMemB;


/*-----------------------------------------------------------------
  Prototypes of internal functions
  -----------------------------------------------------------------*/

int idaDlsFreeB(IDABMem IDAB_mem);


/*=================================================================
  Error Messages
  =================================================================*/

#define MSGD_IDAMEM_NULL "Integrator memory is NULL."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_MEM_FAIL "A memory request failed."
#define MSGD_LMEM_NULL "Linear solver memory is NULL."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."
#define MSGD_MATZERO_FAILED "The SUNMatZero routine failed in an unrecoverable manner."

#define MSGD_CAMEM_NULL "idaadj_mem = NULL illegal."
#define MSGD_LMEMB_NULL "Linear solver memory is NULL for the backward integration."
#define MSGD_BAD_T "Bad t for interpolation."
#define MSGD_BAD_WHICH "Illegal value for which."
#define MSGD_NO_ADJ "Illegal attempt to call before calling IDAAdjInit."

#ifdef __cplusplus
}
#endif

#endif
