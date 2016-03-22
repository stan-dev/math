/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * Common implementation header file for the scaled, preconditioned
 * linear solver modules.
 * -----------------------------------------------------------------
 */

#ifndef _CVSPILS_IMPL_H
#define _CVSPILS_IMPL_H

#include <cvode/cvode_spils.h>
#include "cvode_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Types of iterative linear solvers */

#define SPILS_SPGMR   1
#define SPILS_SPBCG   2
#define SPILS_SPTFQMR 3

/*
 * -----------------------------------------------------------------
 * Types : CVSpilsMemRec, CVSpilsMem
 * -----------------------------------------------------------------
 * The type CVSpilsMem is pointer to a CVSpilsMemRec.
 * -----------------------------------------------------------------
 */

typedef struct CVSpilsMemRec {

  int s_type;           /* type of scaled preconditioned iterative LS   */

  int  s_pretype;       /* type of preconditioning                      */
  int  s_gstype;        /* type of Gram-Schmidt orthogonalization       */
  realtype s_sqrtN;     /* sqrt(N)                                      */
  realtype s_eplifac;   /* eplifac = user specified or EPLIN_DEFAULT    */
  realtype s_deltar;    /* deltar = delt * tq4                          */
  realtype s_delta;     /* delta = deltar * sqrtN                       */
  int  s_maxl;          /* maxl = maximum dimension of the Krylov space */

  long int s_nstlpre;   /* value of nst at the last pset call           */
  long int s_npe;       /* npe = total number of pset calls             */
  long int s_nli;       /* nli = total number of linear iterations      */
  long int s_nps;       /* nps = total number of psolve calls           */
  long int s_ncfl;      /* ncfl = total number of convergence failures  */
  long int s_njtimes;   /* njtimes = total number of calls to jtimes    */
  long int s_nfes;      /* nfeSG = total number of calls to f for     
                           difference quotient Jacobian-vector products */

  N_Vector s_ytemp;     /* temp vector passed to jtimes and psolve      */
  N_Vector s_x;         /* temp vector used by CVSpilsSolve             */
  N_Vector s_ycur;      /* CVODE current y vector in Newton Iteration   */
  N_Vector s_fcur;      /* fcur = f(tn, ycur)                           */

  void* s_spils_mem;    /* memory used by the generic solver            */

  /* Preconditioner computation
   * (a) user-provided:
   *     - P_data == user_data
   *     - pfree == NULL (the user dealocates memory for user_data)
   * (b) internal preconditioner module
   *     - P_data == cvode_mem
   *     - pfree == set by the prec. module and called in CVodeFree
   */
  CVSpilsPrecSetupFn s_pset;
  CVSpilsPrecSolveFn s_psolve;
  void (*s_pfree)(CVodeMem cv_mem);
  void *s_P_data;

  /* Jacobian times vector compuation
   * (a) jtimes function provided by the user:
   *     - j_data == user_data
   *     - jtimesDQ == FALSE
   * (b) internal jtimes
   *     - j_data == cvode_mem
   *     - jtimesDQ == TRUE
   */
  booleantype s_jtimesDQ;
  CVSpilsJacTimesVecFn s_jtimes;
  void *s_j_data;

  long int s_last_flag; /* last error flag returned by any function     */

} *CVSpilsMem;

/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

/* Atimes and PSolve routines called by generic solver */

int CVSpilsAtimes(void *cv_mem, N_Vector v, N_Vector z);

int CVSpilsPSolve(void *cv_mem, N_Vector r, N_Vector z, int lr);

/* Difference quotient approximation for Jac times vector */

int CVSpilsDQJtimes(N_Vector v, N_Vector Jv, realtype t,
                    N_Vector y, N_Vector fy, void *data,
                    N_Vector work);

/*
 * -----------------------------------------------------------------
 * Error Messages
 * -----------------------------------------------------------------
 */

#define MSGS_CVMEM_NULL  "Integrator memory is NULL."
#define MSGS_MEM_FAIL    "A memory request failed."
#define MSGS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGS_BAD_LSTYPE  "Incompatible linear solver type."
#define MSGS_BAD_PRETYPE "Illegal value for pretype. Legal values are PREC_NONE, PREC_LEFT, PREC_RIGHT, and PREC_BOTH."
#define MSGS_PSOLVE_REQ  "pretype != PREC_NONE, but PSOLVE = NULL is illegal."
#define MSGS_LMEM_NULL   "Linear solver memory is NULL."
#define MSGS_BAD_GSTYPE  "Illegal value for gstype. Legal values are MODIFIED_GS and CLASSICAL_GS."
#define MSGS_BAD_EPLIN   "eplifac < 0 illegal."

#define MSGS_PSET_FAILED "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."


#ifdef __cplusplus
}
#endif

#endif
