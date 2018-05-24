/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Alan C. Hindmarsh and Radu Serban @ LLNL
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
 * Implementation header file for the Scaled Preconditioned 
 * Iterative Linear Solver interface.
 *-----------------------------------------------------------------*/

#ifndef _IDASSPILS_IMPL_H
#define _IDASSPILS_IMPL_H

#include <idas/idas_spils.h>
#include "idas_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*-----------------------------------------------------------------
  Types : IDASpilsMemRec, IDASpilsMem                             
  -----------------------------------------------------------------*/

typedef struct IDASpilsMemRec {

  realtype sqrtN;    /* sqrt(N)                                      */
  realtype eplifac;  /* eplifac = linear convergence factor          */
  realtype dqincfac; /* dqincfac = optional increment factor in Jv   */
  realtype epslin;   /* SpgrmSolve tolerance parameter               */

  long int npe;      /* npe = total number of precond calls          */   
  long int nli;      /* nli = total number of linear iterations      */
  long int nps;      /* nps = total number of psolve calls           */
  long int ncfl;     /* ncfl = total number of convergence failures  */
  long int nres;     /* nres = total number of calls to res          */
  long int njtsetup; /* njtsetup = total number of calls to jtsetup  */
  long int njtimes;  /* njtimes = total number of calls to jtimes    */

  long int nst0;     /* nst0 = saved nst (for performance monitor)   */   
  long int nni0;     /* nni0 = saved nni (for performance monitor)   */   
  long int ncfn0;    /* ncfn0 = saved ncfn (for performance monitor) */   
  long int ncfl0;    /* ncfl0 = saved ncfl (for performance monitor) */   
  long int nwarn;    /* nwarn = no. of warnings (for perf. monitor)  */   

  N_Vector ytemp;    /* temp vector used by IDAAtimesDQ              */ 
  N_Vector yptemp;   /* temp vector used by IDAAtimesDQ              */ 
  N_Vector x;        /* temp vector used by the solve function       */
  N_Vector ycur;     /* current y vector in Newton iteration         */
  N_Vector ypcur;    /* current yp vector in Newton iteration        */
  N_Vector rcur;     /* rcur = F(tn, ycur, ypcur)                    */

  SUNLinearSolver LS; /* generic iterative linear solver object      */

  long int last_flag; /* last error return flag                      */

  /* Preconditioner computation
   * (a) user-provided:
   *     - pdata == user_data
   *     - pfree == NULL (the user dealocates memory)
   * (b) internal preconditioner module
   *     - pdata == ida_mem
   *     - pfree == set by the prec. module and called in IDASpilsFree */
  IDASpilsPrecSetupFn pset;
  IDASpilsPrecSolveFn psolve;
  int (*pfree)(IDAMem IDA_mem);
  void *pdata;
  
  /* Jacobian times vector compuation
   * (a) jtimes function provided by the user:
   *     - jdata == user_data
   *     - jtimesDQ == SUNFALSE
   * (b) internal jtimes
   *     - jdata == ida_mem
   *     - jtimesDQ == SUNTRUE */
  booleantype jtimesDQ;
  IDASpilsJacTimesSetupFn jtsetup;
  IDASpilsJacTimesVecFn jtimes;
  void *jdata;

} *IDASpilsMem;


/*-----------------------------------------------------------------
  Prototypes of internal functions
  -----------------------------------------------------------------*/

/* Interface routines called by system SUNLinearSolver */
int IDASpilsATimes(void *ida_mem, N_Vector v, N_Vector z);
int IDASpilsPSetup(void *ida_mem);
int IDASpilsPSolve(void *ida_mem, N_Vector r, N_Vector z,
                   realtype tol, int lr);

/* Difference quotient approximation for Jac times vector */
int IDASpilsDQJtimes(realtype tt, N_Vector yy, N_Vector yp,
                     N_Vector rr, N_Vector v, N_Vector Jv, 
                     realtype c_j, void *data, 
                     N_Vector work1, N_Vector work2);

/* Generic linit/lsetup/lsolve/lfree interface routines for IDA to call */
int idaSpilsInitialize(IDAMem IDA_mem);

int idaSpilsSetup(IDAMem IDA_mem, N_Vector y, N_Vector yp, N_Vector r, 
                  N_Vector vt1, N_Vector vt2, N_Vector vt3); 

int idaSpilsSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                  N_Vector ycur, N_Vector ypcur, N_Vector rescur);

int idaSpilsPerf(IDAMem IDA_mem, int perftask);

int idaSpilsFree(IDAMem IDA_mem);

/* Auxilliary functions */
int idaSpilsInitializeCounters(IDASpilsMem idaspils_mem);

  
/*---------------------------------------------------------------
  Error and Warning Messages
  ---------------------------------------------------------------*/
#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSGS_TIME "at t = %Lg, "
#define MSGS_FRMT "%Le."

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSGS_TIME "at t = %lg, "
#define MSGS_FRMT "%le."

#else

#define MSGS_TIME "at t = %g, "
#define MSGS_FRMT "%e."

#endif


/* Error Messages */
#define MSGS_IDAMEM_NULL   "Integrator memory is NULL."
#define MSGS_MEM_FAIL      "A memory request failed."
#define MSGS_BAD_NVECTOR   "A required vector operation is not implemented."
#define MSGS_BAD_LSTYPE    "Incompatible linear solver type."
#define MSGS_LMEM_NULL     "Linear solver memory is NULL."
#define MSGS_BAD_GSTYPE    "gstype has an illegal value."
#define MSGS_NEG_MAXRS     "maxrs < 0 illegal."
#define MSGS_NEG_EPLIFAC   "eplifac < 0.0 illegal."
#define MSGS_NEG_DQINCFAC  "dqincfac < 0.0 illegal."

#define MSGS_PSET_FAILED "The preconditioner setup routine failed in an unrecoverable manner."
#define MSGS_PSOLVE_FAILED "The preconditioner solve routine failed in an unrecoverable manner."
#define MSGS_JTSETUP_FAILED "The Jacobian x vector setup routine failed in an unrecoverable manner."
#define MSGS_JTIMES_FAILED "The Jacobian x vector routine failed in an unrecoverable manner."

/* Warning Messages */
#define MSGS_WARN  "Warning: " MSGS_TIME "poor iterative algorithm performance. "

#define MSGS_CFN_WARN  MSGS_WARN "Nonlinear convergence failure rate is " MSGS_FRMT
#define MSGS_CFL_WARN  MSGS_WARN "Linear convergence failure rate is " MSGS_FRMT

  
/*-----------------------------------------------------------------
  PART II - backward problems
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Types : IDASpilsMemRecB, IDASpilsMemB       

  IDASpilsSetLinearSolverB attaches such a structure to the lmemB 
  field of IDAadjMem
  -----------------------------------------------------------------*/

typedef struct IDASpilsMemRecB {

  IDASpilsJacTimesSetupFnB jtsetupB;
  IDASpilsJacTimesSetupFnBS jtsetupBS;
  IDASpilsJacTimesVecFnB jtimesB;
  IDASpilsJacTimesVecFnBS jtimesBS;
  IDASpilsPrecSetupFnB psetB;
  IDASpilsPrecSetupFnBS psetBS;
  IDASpilsPrecSolveFnB psolveB;
  IDASpilsPrecSolveFnBS psolveBS;
  void *P_dataB;

} *IDASpilsMemB;


/*-----------------------------------------------------------------
  Prototypes of internal functions
  -----------------------------------------------------------------*/

int idaSpilsFreeB(IDABMem IDAB_mem);

  
/*-----------------------------------------------------------------
  Error Messages 
  -----------------------------------------------------------------*/
#define MSGS_LMEMB_NULL "Linear solver memory is NULL for the backward integration."
#define MSGS_BAD_T      "Bad t for interpolation."
#define MSGS_BAD_WHICH  "Illegal value for which."
#define MSGS_NO_ADJ     "Illegal attempt to call before calling IDAAdjInit."

#define MSGS_LMEMB_NULL  "Linear solver memory is NULL for the backward integration."


#ifdef __cplusplus
}
#endif

#endif
