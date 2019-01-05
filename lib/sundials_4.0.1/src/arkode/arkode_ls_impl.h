/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and
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
 *---------------------------------------------------------------
 * Implementation header file for ARKode's linear solver interface.
 *--------------------------------------------------------------*/

#ifndef _ARKLS_IMPL_H
#define _ARKLS_IMPL_H

#include <arkode/arkode_ls.h>
#include "arkode_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
  ARKLS solver constants:

  ARKLS_MSBJ   default maximum number of steps between Jacobian /
               preconditioner evaluations

  ARKLS_EPLIN  default value for factor by which the tolerance
               on the nonlinear iteration is multiplied to get
               a tolerance on the linear iteration
  ---------------------------------------------------------------*/
#define ARKLS_MSBJ   50
#define ARKLS_EPLIN  RCONST(0.05)


/*---------------------------------------------------------------
  Types: ARKLsMemRec, ARKLsMem

  The type ARKLsMem is pointer to a ARKLsMemRec.
  ---------------------------------------------------------------*/
typedef struct ARKLsMemRec {

  /* Jacobian construction & storage */
  booleantype jacDQ;  /* SUNTRUE if using internal DQ Jacobian approx. */
  ARKLsJacFn jac;     /* Jacobian routine to be called                 */
  void *J_data;       /* user data is passed to jac                    */
  booleantype jbad;   /* heuristic suggestion for pset                 */

  /* Iterative solver tolerance */
  realtype sqrtN;     /* sqrt(N)                                       */
  realtype eplifac;   /* nonlinear -> linear tol scaling factor        */

  /* Linear solver, matrix and vector objects/pointers */
  SUNLinearSolver LS; /* generic linear solver object                  */
  SUNMatrix A;        /* A = M - gamma * df/dy                         */
  SUNMatrix savedJ;   /* savedJ = old Jacobian                         */
  N_Vector ytemp;     /* temp vector passed to jtimes and psolve       */
  N_Vector x;         /* solution vector used by SUNLinearSolver       */
  N_Vector ycur;      /* ptr to current y vector in ARKLs solve        */
  N_Vector fcur;      /* ptr to current fcur = fI(tcur, ycur)          */

  /* Statistics and associated parameters */
  long int msbj;      /* max num steps between jac/pset calls         */
  realtype tcur;      /* 'time' for current ARKLs solve               */
  long int nje;       /* no. of calls to jac                          */
  long int nfeDQ;     /* no. of calls to f due to DQ Jacobian or J*v
                         approximations                               */
  long int nstlj;     /* value of nst at the last jac/pset call       */
  long int npe;       /* npe = total number of pset calls             */
  long int nli;       /* nli = total number of linear iterations      */
  long int nps;       /* nps = total number of psolve calls           */
  long int ncfl;      /* ncfl = total number of convergence failures  */
  long int njtsetup;  /* njtsetup = total number of calls to jtsetup  */
  long int njtimes;   /* njtimes = total number of calls to jtimes    */

  /* Preconditioner computation
    (a) user-provided:
        - P_data == user_data
        - pfree == NULL (the user dealocates memory for user_data)
    (b) internal preconditioner module
        - P_data == arkode_mem
        - pfree == set by the prec. module and called in ARKodeFree  */
  ARKLsPrecSetupFn pset;
  ARKLsPrecSolveFn psolve;
  int (*pfree)(ARKodeMem ark_mem);
  void *P_data;

  /* Jacobian times vector computation
    (a) jtimes function provided by the user:
        - Jt_data == user_data
        - jtimesDQ == SUNFALSE
    (b) internal jtimes
        - Jt_data == arkode_mem
        - jtimesDQ == SUNTRUE   */
  booleantype jtimesDQ;
  ARKLsJacTimesSetupFn jtsetup;
  ARKLsJacTimesVecFn jtimes;
  void *Jt_data;

  long int last_flag; /* last error flag returned by any function */

} *ARKLsMem;


/*---------------------------------------------------------------
  Types: ARKLsMassMemRec, ARKLsMassMem

  The type ARKLsMassMem is pointer to a ARKLsMassMemRec.
  ---------------------------------------------------------------*/
typedef struct ARKLsMassMemRec {

  /* Mass matrix construction & storage */
  ARKLsMassFn mass;   /* user-provided mass matrix routine to call  */
  SUNMatrix M;        /* mass matrix structure                      */
  SUNMatrix M_lu;     /* mass matrix structure for LU decomposition */

  /* Iterative solver tolerance */
  realtype sqrtN;     /* sqrt(N)                                    */
  realtype eplifac;   /* nonlinear -> linear tol scaling factor     */

  /* Statistics and associated parameters */
  booleantype time_dependent;  /* flag whether M depends on t       */
  long int nmsetups;  /* total number of mass matrix-solver setups  */
  long int nmsolves;  /* total number of mass matrix-solver solves  */
  long int nmtsetup;  /* total number of calls to mtsetup           */
  long int nmtimes;   /* total number of calls to mtimes            */
  long int npe;       /* total number of pset calls                 */
  long int nli;       /* total number of linear iterations          */
  long int nps;       /* total number of psolve calls               */
  long int ncfl;      /* total number of convergence failures       */

  /* Linear solver, matrix and vector objects/pointers */
  SUNLinearSolver LS; /* generic linear solver object               */
  N_Vector x;         /* solution vector used by SUNLinearSolver    */
  N_Vector ycur;      /* ptr to ARKode current y vector             */

  /* Preconditioner computation
    (a) user-provided:
        - P_data == user_data
        - pfree == NULL (the user dealocates memory for user_data)
    (b) internal preconditioner module
        - P_data == arkode_mem
        - pfree == set by the prec. module and called in ARKodeFree */
  ARKLsMassPrecSetupFn pset;
  ARKLsMassPrecSolveFn psolve;
  int (*pfree)(ARKodeMem ark_mem);
  void *P_data;

  /* Mass matrix times vector setup and product routines, data */
  ARKLsMassTimesSetupFn mtsetup;
  ARKLsMassTimesVecFn mtimes;
  void *mt_data;

  long int last_flag; /* last error flag returned by any function  */

} *ARKLsMassMem;


/*---------------------------------------------------------------
  Prototypes of internal functions
  ---------------------------------------------------------------*/

/* Interface routines called by system SUNLinearSolver */
int arkLsATimes(void* arkode_mem, N_Vector v, N_Vector z);
int arkLsPSetup(void* arkode_mem);
int arkLsPSolve(void* arkode_mem, N_Vector r, N_Vector z,
                realtype tol, int lr);

/* Interface routines called by mass SUNLinearSolver */
int arkLsMTimes(void* arkode_mem, N_Vector v, N_Vector z);
int arkLsMPSetup(void* arkode_mem);
int arkLsMPSolve(void* arkode_mem, N_Vector r, N_Vector z,
                 realtype tol, int lr);

/* Difference quotient approximation for Jac times vector */
int arkLsDQJtimes(N_Vector v, N_Vector Jv, realtype t,
                  N_Vector y, N_Vector fy, void* data,
                  N_Vector work);

/* Difference-quotient Jacobian approximation routines */
int arkLsDQJac(realtype t, N_Vector y, N_Vector fy,
               SUNMatrix Jac, void* data, N_Vector tmp1,
               N_Vector tmp2, N_Vector tmp3);
int arkLsDenseDQJac(realtype t, N_Vector y, N_Vector fy,
                    SUNMatrix Jac, ARKodeMem ark_mem,
                    ARKLsMem arkls_mem, ARKRhsFn fi, N_Vector tmp1);
int arkLsBandDQJac(realtype t, N_Vector y, N_Vector fy,
                   SUNMatrix Jac, ARKodeMem ark_mem,
                   ARKLsMem arkls_mem, ARKRhsFn fi,
                   N_Vector tmp1, N_Vector tmp2);

/* Generic linit/lsetup/lsolve/lfree interface routines for ARKode to call */
int arkLsInitialize(void* arkode_mem);

int arkLsSetup(void* arkode_mem, int convfail, realtype tpred,
               N_Vector ypred, N_Vector fpred, booleantype* jcurPtr,
               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int arkLsSolve(void* arkode_mem, N_Vector b, realtype tcur,
               N_Vector ycur, N_Vector fcur, realtype eRnrm, int mnewt);

int arkLsFree(void* arkode_mem);

/* Generic minit/msetup/mmult/msolve/mfree routines for ARKode to call */
int arkLsMassInitialize(void* arkode_mem);

int arkLsMassSetup(void* arkode_mem, N_Vector vtemp1,
                   N_Vector vtemp2, N_Vector vtemp3);

int arkLsMassMult(void* arkode_mem, N_Vector v, N_Vector Mv);

int arkLsMassSolve(void* arkode_mem, N_Vector b, realtype nlscoef);

int arkLsMassFree(void* arkode_mem);

/* Auxilliary functions */
int arkLsInitializeCounters(ARKLsMem arkls_mem);

int arkLsInitializeMassCounters(ARKLsMassMem arkls_mem);

int arkLs_AccessLMem(void* arkode_mem, const char* fname,
                     ARKodeMem* ark_mem, ARKLsMem* arkls_mem);

int arkLs_AccessMassMem(void* arkode_mem, const char* fname,
                        ARKodeMem* ark_mem, ARKLsMassMem* arkls_mem);

/* Set/get routines called by time-stepper module */
int arkLSSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix A);

int arkLSSetMassLinearSolver(void* arkode_mem, SUNLinearSolver LS,
                             SUNMatrix M, booleantype time_dep);

int arkLSSetJacFn(void* arkode_mem, ARKLsJacFn jac);
int arkLSSetMassFn(void* arkode_mem, ARKLsMassFn mass);
int arkLSSetEpsLin(void* arkode_mem, realtype eplifac);
int arkLSSetMassEpsLin(void* arkode_mem, realtype eplifac);
int arkLSSetMaxStepsBetweenJac(void* arkode_mem, long int msbj);
int arkLSSetPreconditioner(void* arkode_mem, ARKLsPrecSetupFn psetup,
                           ARKLsPrecSolveFn psolve);
int arkLSSetMassPreconditioner(void* arkode_mem, ARKLsMassPrecSetupFn psetup,
                               ARKLsMassPrecSolveFn psolve);
int arkLSSetJacTimes(void* arkode_mem, ARKLsJacTimesSetupFn jtsetup,
                     ARKLsJacTimesVecFn jtimes);
int arkLSSetMassTimes(void* arkode_mem, ARKLsMassTimesSetupFn msetup,
                      ARKLsMassTimesVecFn mtimes, void* mtimes_data);

int arkLSGetWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS);
int arkLSGetNumJacEvals(void* arkode_mem, long int* njevals);
int arkLSGetNumPrecEvals(void* arkode_mem, long int* npevals);
int arkLSGetNumPrecSolves(void* arkode_mem, long int* npsolves);
int arkLSGetNumLinIters(void* arkode_mem, long int* nliters);
int arkLSGetNumConvFails(void* arkode_mem, long int* nlcfails);
int arkLSGetNumJTSetupEvals(void* arkode_mem, long int* njtsetups);
int arkLSGetNumJtimesEvals(void* arkode_mem, long int* njvevals);
int arkLSGetNumRhsEvals(void* arkode_mem, long int* nfevalsLS);
int arkLSGetLastFlag(void* arkode_mem, long int* flag);

int arkLSGetMassWorkSpace(void* arkode_mem, long int* lenrwMLS,
                          long int* leniwMLS);
int arkLSGetNumMassSetups(void* arkode_mem, long int* nmsetups);
int arkLSGetNumMassMult(void* arkode_mem, long int* nmvevals);
int arkLSGetNumMassSolves(void* arkode_mem, long int* nmsolves);
int arkLSGetNumMassPrecEvals(void* arkode_mem, long int* nmpevals);
int arkLSGetNumMassPrecSolves(void* arkode_mem, long int* nmpsolves);
int arkLSGetNumMassIters(void* arkode_mem, long int* nmiters);
int arkLSGetNumMassConvFails(void* arkode_mem, long int* nmcfails);
int arkLSGetNumMTSetups(void* arkode_mem, long int* nmtsetups);
int arkLSGetLastMassFlag(void* arkode_mem, long int* flag);

char* arkLSGetReturnFlagName(long int flag);

/*---------------------------------------------------------------
  Error Messages
  ---------------------------------------------------------------*/
#define MSG_LS_ARKMEM_NULL     "Integrator memory is NULL."
#define MSG_LS_MEM_FAIL        "A memory request failed."
#define MSG_LS_BAD_NVECTOR     "A required vector operation is not implemented."
#define MSG_LS_BAD_LSTYPE      "Incompatible linear solver type."
#define MSG_LS_LMEM_NULL       "Linear solver memory is NULL."
#define MSG_LS_MASSMEM_NULL    "Mass matrix solver memory is NULL."
#define MSG_LS_BAD_SIZES       "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."

#define MSG_LS_PSET_FAILED     "The preconditioner setup routine failed in an unrecoverable manner."
#define MSG_LS_PSOLVE_FAILED   "The preconditioner solve routine failed in an unrecoverable manner."
#define MSG_LS_JTSETUP_FAILED  "The Jacobian x vector setup routine failed in an unrecoverable manner."
#define MSG_LS_JTIMES_FAILED   "The Jacobian x vector routine failed in an unrecoverable manner."
#define MSG_LS_MTSETUP_FAILED  "The mass matrix x vector setup routine failed in an unrecoverable manner."
#define MSG_LS_MTIMES_FAILED   "The mass matrix x vector routine failed in an unrecoverable manner."

#define MSG_LS_JACFUNC_FAILED  "The Jacobian routine failed in an unrecoverable manner."
#define MSG_LS_MASSFUNC_FAILED "The mass matrix routine failed in an unrecoverable manner."
#define MSG_LS_SUNMAT_FAILED   "A SUNMatrix routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
