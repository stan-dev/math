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
 * Implementation header file for the scaled, preconditioned
 * linear solver interface.
 * -----------------------------------------------------------------
 */

#ifndef _CVSSPILS_IMPL_H
#define _CVSSPILS_IMPL_H

#include <cvodes/cvodes_spils.h>
#include "cvodes_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  PART I:  Forward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  Types : CVSpilsMemRec, CVSpilsMem
  -----------------------------------------------------------------
  The type CVSpilsMem is pointer to a CVSpilsMemRec.
  -----------------------------------------------------------------*/

typedef struct CVSpilsMemRec {

  realtype sqrtN;     /* sqrt(N)                                      */
  realtype eplifac;   /* eplifac = user specified or EPLIN_DEFAULT    */
  realtype deltar;    /* deltar = delt * tq4                          */
  realtype delta;     /* delta = deltar * sqrtN                       */

  booleantype jbad;   /* heuristic suggestion for pset/jtsetup        */
  long int nstlpre;   /* value of nst at the last pset call           */
  long int npe;       /* npe = total number of pset calls             */
  long int nli;       /* nli = total number of linear iterations      */
  long int nps;       /* nps = total number of psolve calls           */
  long int ncfl;      /* ncfl = total number of convergence failures  */
  long int njtsetup;  /* njtsetup = total number of calls to jtsetup  */
  long int njtimes;   /* njtimes = total number of calls to jtimes    */
  long int nfes;      /* nfeSG = total number of calls to f for     
                         difference quotient Jacobian-vector products */

  SUNLinearSolver LS; /* generic iterative linear solver object       */
  
  N_Vector ytemp;     /* temp vector passed to jtimes and psolve      */
  N_Vector x;         /* temp vector used by CVSpilsSolve             */
  N_Vector ycur;      /* CVODE current y vector in Newton Iteration   */
  N_Vector fcur;      /* fcur = f(tn, ycur)                           */

  /* Preconditioner computation
   * (a) user-provided:
   *     - P_data == user_data
   *     - pfree == NULL (the user dealocates memory for user_data)
   * (b) internal preconditioner module
   *     - P_data == cvode_mem
   *     - pfree == set by the prec. module and called in CVodeFree */
  CVSpilsPrecSetupFn pset;
  CVSpilsPrecSolveFn psolve;
  int (*pfree)(CVodeMem cv_mem);
  void *P_data;

  /* Jacobian times vector compuation
   * (a) jtimes function provided by the user:
   *     - j_data == user_data
   *     - jtimesDQ == SUNFALSE
   * (b) internal jtimes
   *     - j_data == cvode_mem
   *     - jtimesDQ == SUNTRUE  */
  booleantype jtimesDQ;
  CVSpilsJacTimesSetupFn jtsetup;
  CVSpilsJacTimesVecFn jtimes;
  void *j_data;

  long int last_flag;    /* last error flag returned by any function   */

} *CVSpilsMem;

/*-----------------------------------------------------------------
  Prototypes of internal functions
  -----------------------------------------------------------------*/

/* Interface routines called by system SUNLinearSolver */
int CVSpilsATimes(void *cv_mem, N_Vector v, N_Vector z);
int CVSpilsPSetup(void *cv_mem);
int CVSpilsPSolve(void *cv_mem, N_Vector r, N_Vector z,
                  realtype tol, int lr);

/* Difference quotient approximation for Jac times vector */
int CVSpilsDQJtimes(N_Vector v, N_Vector Jv, realtype t,
                    N_Vector y, N_Vector fy, void *data,
                    N_Vector work);

/* Generic linit/lsetup/lsolve/lfree interface routines for CVode to call */
int cvSpilsInitialize(CVodeMem cv_mem);
int cvSpilsSetup(CVodeMem cv_mem, int convfail, N_Vector y, 
                 N_Vector fy, booleantype *jcurPtr, 
                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3); 
int cvSpilsSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                 N_Vector ycur, N_Vector fcur);
int cvSpilsFree(CVodeMem cv_mem);

/* Auxilliary functions */
int cvSpilsInitializeCounters(CVSpilsMem cvspils_mem);


/*=================================================================
  PART II:  Backward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  Types : CVSpilsMemRecB, CVSpilsMemB       
  -----------------------------------------------------------------
  CVSpgmrB, CVSpbcgB, and CVSptfqmr attach such a structure to the 
  lmemB filed of CVodeBMem
  -----------------------------------------------------------------*/

typedef struct CVSpilsMemRecB {

  CVSpilsJacTimesSetupFnB jtsetupB;
  CVSpilsJacTimesSetupFnBS jtsetupBS;
  CVSpilsJacTimesVecFnB jtimesB;
  CVSpilsJacTimesVecFnBS jtimesBS;
  CVSpilsPrecSetupFnB psetB;
  CVSpilsPrecSetupFnBS psetBS;
  CVSpilsPrecSolveFnB psolveB;
  CVSpilsPrecSolveFnBS psolveBS;
  void *P_dataB;

} *CVSpilsMemB;


/*-----------------------------------------------------------------
  Prototypes of internal functions
  -----------------------------------------------------------------*/

int cvSpilsFreeB(CVodeBMem cvb_mem);

/*=================================================================
  Error Messages
  =================================================================*/

#define MSGS_CVMEM_NULL  "Integrator memory is NULL."
#define MSGS_MEM_FAIL    "A memory request failed."
#define MSGS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGS_BAD_LSTYPE  "Incompatible linear solver type."
#define MSGS_BAD_PRETYPE "Illegal value for pretype. Legal values are PREC_NONE, PREC_LEFT, PREC_RIGHT, and PREC_BOTH."
#define MSGS_PSOLVE_REQ  "pretype != PREC_NONE, but PSOLVE = NULL is illegal."
#define MSGS_LMEM_NULL   "Linear solver memory is NULL."
#define MSGS_BAD_GSTYPE  "Illegal value for gstype. Legal values are MODIFIED_GS and CLASSICAL_GS."
#define MSGS_BAD_EPLIN    "eplifac < 0 illegal."
  
#define MSGS_PSET_FAILED   "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTSETUP_FAILED "The Jacobian x vector setup routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."

#define MSGS_NO_ADJ      "Illegal attempt to call before calling CVodeAdjMalloc."
#define MSGS_BAD_WHICH   "Illegal value for which."
#define MSGS_LMEMB_NULL  "Linear solver memory is NULL for the backward integration."
#define MSGS_BAD_TINTERP "Bad t for interpolation."


#ifdef __cplusplus
}
#endif

#endif
