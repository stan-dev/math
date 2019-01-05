/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                Scott D. Cohen, Alan C. Hindmarsh, Radu Serban
 *                   and Dan Shumaker @ LLNL
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
 * Implementation header file for the main CVODE integrator.
 * -----------------------------------------------------------------
 */

#ifndef _CVODE_IMPL_H
#define _CVODE_IMPL_H

#include <stdarg.h>
#include <cvode/cvode.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 *   M A I N    I N T E G R A T O R    M E M O R Y    B L O C K
 * =================================================================
 */

/* Basic CVODE constants */

#define ADAMS_Q_MAX 12     /* max value of q for lmm == ADAMS     */
#define BDF_Q_MAX    5     /* max value of q for lmm == BDF       */
#define Q_MAX  ADAMS_Q_MAX /* max value of q for either lmm       */
#define L_MAX  (Q_MAX+1)   /* max value of L for either lmm       */
#define NUM_TESTS    5     /* number of error test quantities     */

#define HMIN_DEFAULT     RCONST(0.0)    /* hmin default value     */
#define HMAX_INV_DEFAULT RCONST(0.0)    /* hmax_inv default value */
#define MXHNIL_DEFAULT   10             /* mxhnil default value   */
#define MXSTEP_DEFAULT   500            /* mxstep default value   */

/* Return values for lower level routines used by CVode and functions
   provided to the nonlinear solver */

#define RHSFUNC_RECVR +9

/*
 * -----------------------------------------------------------------
 * Types : struct CVodeMemRec, CVodeMem
 * -----------------------------------------------------------------
 * The type CVodeMem is type pointer to struct CVodeMemRec.
 * This structure contains fields to keep track of problem state.
 * -----------------------------------------------------------------
 */

typedef struct CVodeMemRec {

  realtype cv_uround;    /* machine unit roundoff */

  /*--------------------------
    Problem Specification Data
    --------------------------*/

  CVRhsFn cv_f;              /* y' = f(t,y(t))                                */
  void *cv_user_data;        /* user pointer passed to f                      */
  int cv_lmm;                /* lmm = CV_ADAMS or CV_BDF                      */
  int cv_itol;               /* itol = CV_SS, CV_SV, CV_WF, CV_NN             */

  realtype cv_reltol;        /* relative tolerance                            */
  realtype cv_Sabstol;       /* scalar absolute tolerance                     */
  N_Vector cv_Vabstol;       /* vector absolute tolerance                     */
  booleantype cv_user_efun;  /* SUNTRUE if user sets efun                     */
  CVEwtFn cv_efun;           /* function to set ewt                           */
  void *cv_e_data;           /* user pointer passed to efun                   */

  booleantype cv_constraintsSet; /* constraints vector present:
                                    do constraints calc                       */

  /*-----------------------
    Nordsieck History Array
    -----------------------*/

  N_Vector cv_zn[L_MAX];  /* Nordsieck array, of size N x (q+1).
                             zn[j] is a vector of length N (j=0,...,q)
                             zn[j] = [1/factorial(j)] * h^j * (jth
                             derivative of the interpolating polynomial       */

  /*--------------------------
    other vectors of length N
    -------------------------*/

  N_Vector cv_ewt;     /* error weight vector                                 */
  N_Vector cv_y;       /* y is used as temporary storage by the solver
                          The memory is provided by the user to CVode
                          where the vector is named yout.                     */
  N_Vector cv_acor;    /* In the context of the solution of the nonlinear
                          equation, acor = y_n(m) - y_n(0). On return,
                          this vector is scaled to give the est. local err.   */
  N_Vector cv_tempv;   /* temporary storage vector                            */
  N_Vector cv_ftemp;   /* temporary storage vector                            */
  N_Vector cv_vtemp1;  /* temporary storage vector                            */
  N_Vector cv_vtemp2;  /* temporary storage vector                            */
  N_Vector cv_vtemp3;  /* temporary storage vector                            */

 N_Vector cv_mm;          /* mask vector in constraints tests                    */
 N_Vector cv_constraints; /* vector of inequality constraint options             */

  /*-----------------
    Tstop information
    -----------------*/

  booleantype cv_tstopset;
  realtype cv_tstop;

  /*---------
    Step Data
    ---------*/

  int cv_q;                    /* current order                               */
  int cv_qprime;               /* order to be used on the next step
                                  = q-1, q, or q+1                            */
  int cv_next_q;               /* order to be used on the next step           */
  int cv_qwait;                /* number of internal steps to wait before
                                  considering a change in q                   */
  int cv_L;                    /* L = q + 1                                   */

  realtype cv_hin;             /* initial step size                           */
  realtype cv_h;               /* current step size                           */
  realtype cv_hprime;          /* step size to be used on the next step       */
  realtype cv_next_h;          /* step size to be used on the next step       */
  realtype cv_eta;             /* eta = hprime / h                            */
  realtype cv_hscale;          /* value of h used in zn                       */
  realtype cv_tn;              /* current internal value of t                 */
  realtype cv_tretlast;        /* value of tret last returned by CVode        */

  realtype cv_tau[L_MAX+1];    /* array of previous q+1 successful step
                                  sizes indexed from 1 to q+1                 */
  realtype cv_tq[NUM_TESTS+1]; /* array of test quantities indexed from
                                  1 to NUM_TESTS(=5)                          */
  realtype cv_l[L_MAX];        /* coefficients of l(x) (degree q poly)        */

  realtype cv_rl1;              /* the scalar 1/l[1]                          */
  realtype cv_gamma;            /* gamma = h * rl1                            */
  realtype cv_gammap;           /* gamma at the last setup call               */
  realtype cv_gamrat;           /* gamma / gammap                             */

  realtype cv_crate;            /* estimated corrector convergence rate       */
  realtype cv_delp;             /* norm of previous nonlinear solver update   */
  realtype cv_acnrm;            /* | acor | wrms                              */
  realtype cv_nlscoef;          /* coeficient in nonlinear convergence test   */

  /*------
    Limits
    ------*/

  int cv_qmax;          /* q <= qmax                                          */
  long int cv_mxstep;   /* maximum number of internal steps for one user call */
  int cv_maxcor;        /* maximum number of corrector iterations for the
                           solution of the nonlinear equation                 */
  int cv_mxhnil;        /* maximum number of warning messages issued to the
                           user that t + h == t for the next internal step    */
  int cv_maxnef;        /* maximum number of error test failures              */
  int cv_maxncf;        /* maximum number of nonlinear convergence failures   */

  realtype cv_hmin;     /* |h| >= hmin                                        */
  realtype cv_hmax_inv; /* |h| <= 1/hmax_inv                                  */
  realtype cv_etamax;   /* eta <= etamax                                      */

  /*--------
    Counters
    --------*/

  long int cv_nst;         /* number of internal steps taken                  */
  long int cv_nfe;         /* number of f calls                               */
  long int cv_ncfn;        /* number of corrector convergence failures        */
  long int cv_netf;        /* number of error test failures                   */
  long int cv_nni;         /* number of Newton iterations performed           */
  long int cv_nsetups;     /* number of setup calls                           */
  int cv_nhnil;            /* number of messages issued to the user that
                              t + h == t for the next iternal step            */

  realtype cv_etaqm1;      /* ratio of new to old h for order q-1             */
  realtype cv_etaq;        /* ratio of new to old h for order q               */
  realtype cv_etaqp1;      /* ratio of new to old h for order q+1             */

  /*----------------------------
    Space requirements for CVODE
    ----------------------------*/

  sunindextype cv_lrw1;        /* no. of realtype words in 1 N_Vector             */
  sunindextype cv_liw1;        /* no. of integer words in 1 N_Vector              */
  long int cv_lrw;             /* no. of realtype words in CVODE work vectors     */
  long int cv_liw;             /* no. of integer words in CVODE work vectors      */

  /*---------------------
    Nonlinear Solver Data
    ---------------------*/

  SUNNonlinearSolver NLS;  /* Sundials generic nonlinear solver object */
  booleantype ownNLS;      /* flag indicating if CVODE created the nonlinear
                              solver object */
  int convfail;            /* flag to indicate when a Jacbian update may
                              be needed */

  /*------------------
    Linear Solver Data
    ------------------*/

  /* Linear Solver functions to be called */

  int (*cv_linit)(struct CVodeMemRec *cv_mem);

  int (*cv_lsetup)(struct CVodeMemRec *cv_mem, int convfail, N_Vector ypred,
		   N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
		   N_Vector vtemp2, N_Vector vtemp3);

  int (*cv_lsolve)(struct CVodeMemRec *cv_mem, N_Vector b, N_Vector weight,
		   N_Vector ycur, N_Vector fcur);

  int (*cv_lfree)(struct CVodeMemRec *cv_mem);

  /* Linear Solver specific memory */

  void *cv_lmem;

  /*------------
    Saved Values
    ------------*/

  int cv_qu;                   /* last successful q value used                */
  long int cv_nstlp;           /* step number of last setup call              */
  realtype cv_h0u;             /* actual initial stepsize                     */
  realtype cv_hu;              /* last successful h value used                */
  realtype cv_saved_tq5;       /* saved value of tq[5]                        */
  booleantype cv_jcur;         /* is Jacobian info. for lin. solver current?  */
  realtype cv_tolsf;           /* tolerance scale factor                      */
  int cv_qmax_alloc;           /* value of qmax used when allocating memory   */
  int cv_indx_acor;            /* index of the zn vector with saved acor      */

  booleantype cv_VabstolMallocDone;
  booleantype cv_MallocDone;  
  booleantype cv_constraintsMallocDone;

  /*-------------------------------------------
    Error handler function and error ouput file
    -------------------------------------------*/

  CVErrHandlerFn cv_ehfun;    /* error messages are handled by ehfun          */
  void *cv_eh_data;           /* data pointer passed to ehfun                 */
  FILE *cv_errfp;             /* CVODE error messages are sent to errfp       */

  /*-------------------------
    Stability Limit Detection
    -------------------------*/

  booleantype cv_sldeton;     /* is Stability Limit Detection on?             */
  realtype cv_ssdat[6][4];    /* scaled data array for STALD                  */
  int cv_nscon;               /* counter for STALD method                     */
  long int cv_nor;            /* counter for number of order reductions       */

  /*----------------
    Rootfinding Data
    ----------------*/

  CVRootFn cv_gfun;        /* function g for roots sought                     */
  int cv_nrtfn;            /* number of components of g                       */
  int *cv_iroots;          /* array for root information                      */
  int *cv_rootdir;         /* array specifying direction of zero-crossing     */
  realtype cv_tlo;         /* nearest endpoint of interval in root search     */
  realtype cv_thi;         /* farthest endpoint of interval in root search    */
  realtype cv_trout;       /* t value returned by rootfinding routine         */
  realtype *cv_glo;        /* saved array of g values at t = tlo              */
  realtype *cv_ghi;        /* saved array of g values at t = thi              */
  realtype *cv_grout;      /* array of g values at t = trout                  */
  realtype cv_toutc;       /* copy of tout (if NORMAL mode)                   */
  realtype cv_ttol;        /* tolerance on root location                      */
  int cv_taskc;            /* copy of parameter itask                         */
  int cv_irfnd;            /* flag showing whether last step had a root       */
  long int cv_nge;         /* counter for g evaluations                       */
  booleantype *cv_gactive; /* array with active/inactive event functions      */
  int cv_mxgnull;          /* number of warning messages about possible g==0  */

  /*-----------------------
    Fused Vector Operations
    -----------------------*/

  realtype cv_cvals[L_MAX]; /* array of scalars */
  N_Vector cv_Xvecs[L_MAX]; /* array of vectors */

} *CVodeMem;

/*
 * =================================================================
 *     I N T E R F A C E   T O    L I N E A R   S O L V E R S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Communication between CVODE and a CVODE Linear Solver
 * -----------------------------------------------------------------
 * convfail (input to cv_lsetup)
 *
 * CV_NO_FAILURES : Either this is the first cv_setup call for this
 *                  step, or the local error test failed on the
 *                  previous attempt at this step (but the Newton
 *                  iteration converged).
 *
 * CV_FAIL_BAD_J  : This value is passed to cv_lsetup if
 *
 *                  (a) The previous Newton corrector iteration
 *                      did not converge and the linear solver's
 *                      setup routine indicated that its Jacobian-
 *                      related data is not current
 *                                   or
 *                  (b) During the previous Newton corrector
 *                      iteration, the linear solver's solve routine
 *                      failed in a recoverable manner and the
 *                      linear solver's setup routine indicated that
 *                      its Jacobian-related data is not current.
 *
 * CV_FAIL_OTHER  : During the current internal step try, the
 *                  previous Newton iteration failed to converge
 *                  even though the linear solver was using current
 *                  Jacobian-related data.
 * -----------------------------------------------------------------
 */

/* Constants for convfail (input to cv_lsetup) */

#define CV_NO_FAILURES 0
#define CV_FAIL_BAD_J  1
#define CV_FAIL_OTHER  2

/*
 * -----------------------------------------------------------------
 * int (*cv_linit)(CVodeMem cv_mem);
 * -----------------------------------------------------------------
 * The purpose of cv_linit is to complete initializations for a
 * specific linear solver, such as counters and statistics.
 * An LInitFn should return 0 if it has successfully initialized the
 * CVODE linear solver and a negative value otherwise.
 * If an error does occur, an appropriate message should be sent to
 * the error handler function.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * int (*cv_lsetup)(CVodeMem cv_mem, int convfail, N_Vector ypred,
 *                 N_Vector fpred, booleantype *jcurPtr,
 *                 N_Vector vtemp1, N_Vector vtemp2,
 *                 N_Vector vtemp3);
 * -----------------------------------------------------------------
 * The job of cv_lsetup is to prepare the linear solver for
 * subsequent calls to cv_lsolve. It may recompute Jacobian-
 * related data is it deems necessary. Its parameters are as
 * follows:
 *
 * cv_mem - problem memory pointer of type CVodeMem. See the
 *          typedef earlier in this file.
 *
 * convfail - a flag to indicate any problem that occurred during
 *            the solution of the nonlinear equation on the
 *            current time step for which the linear solver is
 *            being used. This flag can be used to help decide
 *            whether the Jacobian data kept by a CVODE linear
 *            solver needs to be updated or not.
 *            Its possible values have been documented above.
 *
 * ypred - the predicted y vector for the current CVODE internal
 *         step.
 *
 * fpred - f(tn, ypred).
 *
 * jcurPtr - a pointer to a boolean to be filled in by cv_lsetup.
 *           The function should set *jcurPtr=SUNTRUE if its Jacobian
 *           data is current after the call and should set
 *           *jcurPtr=SUNFALSE if its Jacobian data is not current.
 *           Note: If cv_lsetup calls for re-evaluation of
 *           Jacobian data (based on convfail and CVODE state
 *           data), it should return *jcurPtr=SUNTRUE always;
 *           otherwise an infinite loop can result.
 *
 * vtemp1 - temporary N_Vector provided for use by cv_lsetup.
 *
 * vtemp3 - temporary N_Vector provided for use by cv_lsetup.
 *
 * vtemp3 - temporary N_Vector provided for use by cv_lsetup.
 *
 * The cv_lsetup routine should return 0 if successful, a positive
 * value for a recoverable error, and a negative value for an
 * unrecoverable error.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * int (*cv_lsolve)(CVodeMem cv_mem, N_Vector b, N_Vector weight,
 *                  N_Vector ycur, N_Vector fcur);
 * -----------------------------------------------------------------
 * cv_lsolve must solve the linear equation P x = b, where
 * P is some approximation to (I - gamma J), J = (df/dy)(tn,ycur)
 * and the RHS vector b is input. The N-vector ycur contains
 * the solver's current approximation to y(tn) and the vector
 * fcur contains the N_Vector f(tn,ycur). The solution is to be
 * returned in the vector b. cv_lsolve returns a positive value
 * for a recoverable error and a negative value for an
 * unrecoverable error. Success is indicated by a 0 return value.
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * int (*cv_lfree)(CVodeMem cv_mem);
 * -----------------------------------------------------------------
 * cv_lfree should free up any memory allocated by the linear
 * solver. This routine is called once a problem has been
 * completed and the linear solver is no longer needed.  It should
 * return 0 upon success, nonzero on failure.
 * -----------------------------------------------------------------
 */

/*
 * =================================================================
 *   C V O D E    I N T E R N A L   F U N C T I O N S
 * =================================================================
 */

/* Prototype of internal ewtSet function */

int cvEwtSet(N_Vector ycur, N_Vector weight, void *data);

/* High level error handler */

void cvProcessError(CVodeMem cv_mem,
		    int error_code, const char *module, const char *fname,
		    const char *msgfmt, ...);

/* Prototype of internal ErrHandler function */

void cvErrHandler(int error_code, const char *module, const char *function,
		  char *msg, void *data);

/* Nonlinear solver initializtion function */

int cvNlsInit(CVodeMem cv_mem);

/*
 * =================================================================
 *   C V O D E    E R R O R    M E S S A G E S
 * =================================================================
 */

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define MSG_TIME        "t = %Lg"
#define MSG_TIME_H      "t = %Lg and h = %Lg"
#define MSG_TIME_INT    "t = %Lg is not between tcur - hu = %Lg and tcur = %Lg."
#define MSG_TIME_TOUT   "tout = %Lg"
#define MSG_TIME_TSTOP  "tstop = %Lg"

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define MSG_TIME        "t = %lg"
#define MSG_TIME_H      "t = %lg and h = %lg"
#define MSG_TIME_INT    "t = %lg is not between tcur - hu = %lg and tcur = %lg."
#define MSG_TIME_TOUT   "tout = %lg"
#define MSG_TIME_TSTOP  "tstop = %lg"

#else

#define MSG_TIME        "t = %g"
#define MSG_TIME_H      "t = %g and h = %g"
#define MSG_TIME_INT    "t = %g is not between tcur - hu = %g and tcur = %g."
#define MSG_TIME_TOUT   "tout = %g"
#define MSG_TIME_TSTOP  "tstop = %g"

#endif

/* Initialization and I/O error messages */

#define MSGCV_NO_MEM "cvode_mem = NULL illegal."
#define MSGCV_CVMEM_FAIL "Allocation of cvode_mem failed."
#define MSGCV_MEM_FAIL "A memory request failed."
#define MSGCV_BAD_LMM  "Illegal value for lmm. The legal values are CV_ADAMS and CV_BDF."
#define MSGCV_NO_MALLOC "Attempt to call before CVodeInit."
#define MSGCV_NEG_MAXORD "maxord <= 0 illegal."
#define MSGCV_BAD_MAXORD  "Illegal attempt to increase maximum method order."
#define MSGCV_SET_SLDET  "Attempt to use stability limit detection with the CV_ADAMS method illegal."
#define MSGCV_NEG_HMIN "hmin < 0 illegal."
#define MSGCV_NEG_HMAX "hmax < 0 illegal."
#define MSGCV_BAD_HMIN_HMAX "Inconsistent step size limits: hmin > hmax."
#define MSGCV_BAD_RELTOL "reltol < 0 illegal."
#define MSGCV_BAD_ABSTOL "abstol has negative component(s) (illegal)."
#define MSGCV_NULL_ABSTOL "abstol = NULL illegal."
#define MSGCV_NULL_Y0 "y0 = NULL illegal."
#define MSGCV_Y0_FAIL_CONSTR "y0 fails to satisfy constraints."
#define MSGCV_NULL_F "f = NULL illegal."
#define MSGCV_NULL_G "g = NULL illegal."
#define MSGCV_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGCV_BAD_CONSTR "Illegal values in constraints vector."
#define MSGCV_BAD_K "Illegal value for k."
#define MSGCV_NULL_DKY "dky = NULL illegal."
#define MSGCV_BAD_T "Illegal value for t." MSG_TIME_INT
#define MSGCV_NO_ROOT "Rootfinding was not initialized."
#define MSGCV_NLS_INIT_FAIL "The nonlinear solver's init routine failed."

/* CVode Error Messages */

#define MSGCV_NO_TOLS "No integration tolerances have been specified."
#define MSGCV_LSOLVE_NULL "The linear solver's solve routine is NULL."
#define MSGCV_YOUT_NULL "yout = NULL illegal."
#define MSGCV_TRET_NULL "tret = NULL illegal."
#define MSGCV_BAD_EWT "Initial ewt has component(s) equal to zero (illegal)."
#define MSGCV_EWT_NOW_BAD "At " MSG_TIME ", a component of ewt has become <= 0."
#define MSGCV_BAD_ITASK "Illegal value for itask."
#define MSGCV_BAD_H0 "h0 and tout - t0 inconsistent."
#define MSGCV_BAD_TOUT "Trouble interpolating at " MSG_TIME_TOUT ". tout too far back in direction of integration"
#define MSGCV_EWT_FAIL "The user-provide EwtSet function failed."
#define MSGCV_EWT_NOW_FAIL "At " MSG_TIME ", the user-provide EwtSet function failed."
#define MSGCV_LINIT_FAIL "The linear solver's init routine failed."
#define MSGCV_HNIL_DONE "The above warning has been issued mxhnil times and will not be issued again for this problem."
#define MSGCV_TOO_CLOSE "tout too close to t0 to start integration."
#define MSGCV_MAX_STEPS "At " MSG_TIME ", mxstep steps taken before reaching tout."
#define MSGCV_TOO_MUCH_ACC "At " MSG_TIME ", too much accuracy requested."
#define MSGCV_HNIL "Internal " MSG_TIME_H " are such that t + h = t on the next step. The solver will continue anyway."
#define MSGCV_ERR_FAILS "At " MSG_TIME_H ", the error test failed repeatedly or with |h| = hmin."
#define MSGCV_CONV_FAILS "At " MSG_TIME_H ", the corrector convergence test failed repeatedly or with |h| = hmin."
#define MSGCV_SETUP_FAILED "At " MSG_TIME ", the setup routine failed in an unrecoverable manner."
#define MSGCV_SOLVE_FAILED "At " MSG_TIME ", the solve routine failed in an unrecoverable manner."
#define MSGCV_FAILED_CONSTR "At " MSG_TIME ", unable to satisfy inequality constraints."
#define MSGCV_RHSFUNC_FAILED "At " MSG_TIME ", the right-hand side routine failed in an unrecoverable manner."
#define MSGCV_RHSFUNC_UNREC "At " MSG_TIME ", the right-hand side failed in a recoverable manner, but no recovery is possible."
#define MSGCV_RHSFUNC_REPTD "At " MSG_TIME " repeated recoverable right-hand side function errors."
#define MSGCV_RHSFUNC_FIRST "The right-hand side routine failed at the first call."
#define MSGCV_RTFUNC_FAILED "At " MSG_TIME ", the rootfinding routine failed in an unrecoverable manner."
#define MSGCV_CLOSE_ROOTS "Root found at and very near " MSG_TIME "."
#define MSGCV_BAD_TSTOP "The value " MSG_TIME_TSTOP " is behind current " MSG_TIME " in the direction of integration."
#define MSGCV_INACTIVE_ROOTS "At the end of the first step, there are still some root functions identically 0. This warning will not be issued again."
#define MSGCV_NLS_SETUP_FAILED "At " MSG_TIME "the nonlinear solver setup failed unrecoverably."
#define MSGCV_NLS_INPUT_NULL "At " MSG_TIME "the nonlinear solver was passed a NULL input."

#ifdef __cplusplus
}
#endif

#endif
