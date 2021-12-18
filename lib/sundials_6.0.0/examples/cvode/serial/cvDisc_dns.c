/*
 * -----------------------------------------------------------------
 * Programmers: Radu Serban, Alan Hindmarsh, and Cody Balos @ LLNL
 * -----------------------------------------------------------------
 * Simple 1D example to illustrate integrating over discontinuities:
 *
 * A) Discontinuity in solution
 *       y' = -y   ; y(0) = 1    ; t = [0,1]
 *       y' = -y   ; y(1) = 1    ; t = [1,2]
 *
 * B) Discontinuity in RHS (y')
 *       y' = -y   ; y(0) = 1    ; t = [0,1]
 *       z' = -5*z ; z(1) = y(1) ; t = [1,2]
 *    This case is solved twice, first by explicitly treating the
 *    discontinuity point and secondly by letting the integrator
 *    deal with the discontinuity.
 * -----------------------------------------------------------------
 */

#include <stdio.h>

#include <cvode/cvode.h>               /* prototypes for CVODE functions and const */
#include <nvector/nvector_serial.h>    /* access to serial NVector                 */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix                */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver          */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype          */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* Problem Constants */
#define NEQ  1 /* number of equations */

#define RHS1 1
#define RHS2 2

/* User provided routine called by the solver to compute
 * the function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data);

/* Function for checking return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

int main()
{
  void *cvode_mem;
  SUNMatrix A;
  SUNLinearSolver LS;
  SUNContext sunctx;

  N_Vector y;
  int flag, retval;
  realtype reltol, abstol, t0, t1, t2, t;
  long int nst1, nst2, nst;

  reltol = RCONST(1.0e-3);
  abstol = RCONST(1.0e-4);

  t0 = RCONST(0.0);
  t1 = RCONST(1.0);
  t2 = RCONST(2.0);

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if(check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* Allocate the vector of initial conditions */
  y = N_VNew_Serial(NEQ, sunctx);

  /* Set initial condition */
  NV_Ith_S(y,0) = RCONST(1.0);

  /*
   * ------------------------------------------------------------
   *  Shared initialization and setup
   * ------------------------------------------------------------
   */

  /* Call CVodeCreate to create CVODE memory block and specify the
   * Backward Differentiaion Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Call CVodeInit to initialize integrator memory and specify the
   * user's right hand side function y'=f(t,y), the initial time T0
   * and the initial condiition vector y. */
  retval = CVodeInit(cvode_mem, f, t0, y);
  if (check_retval((void *)&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSStolerances to specify integration tolereances,
   * specifically the scalar relative and absolute tolerance. */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval((void *)&retval, "CVodeSStolerances", 1)) return(1);

  /* Provide RHS flag as user data which can be access in user provided routines */
  retval = CVodeSetUserData(cvode_mem, &flag);
  if (check_retval((void *)&retval, "CVodeSetUserData", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solver */
  A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense linear solver for use by CVode */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the linear solver and matrix to CVode by calling CVodeSetLinearSolver */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval((void *)&retval, "CVodeSetLinearSolver", 1)) return(1);

  /*
   * ---------------------------------------------------------------
   * Discontinuity in the solution
   *
   * 1) Integrate to the discontinuity
   * 2) Integrate from the discontinuity
   * ---------------------------------------------------------------
   */

  /* ---- Integrate to the discontinuity */

  printf("\nDiscontinuity in solution\n\n");

  /* set TSTOP (max time solution proceeds to) - this is not required */
  retval = CVodeSetStopTime(cvode_mem, t1);
  if (check_retval((void *)&retval, "CVodeSetStopTime", 1)) return(1);

  flag = RHS1; /* use -y for RHS */
  t = t0; /* set the integrator start time */

  printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));
  while (t<t1) {
    /* advance solver just one internal step */
    retval = CVode(cvode_mem, t1, y, &t, CV_ONE_STEP);
    if (check_retval((void *)&retval, "CVode", 1)) return(1);
    printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));
  }
  /* Get the number of steps the solver took to get to the discont. */
  retval = CVodeGetNumSteps(cvode_mem, &nst1);
  if (check_retval((void *)&retval, "CvodeGetNumSteps", 1)) return(1);

  /* ---- Integrate from the discontinuity */

  /* Include discontinuity */
  NV_Ith_S(y,0) = RCONST(1.0);

  /* Reinitialize the solver */
  retval = CVodeReInit(cvode_mem, t1, y);
  if (check_retval((void *)&retval, "CVodeReInit", 1)) return(1);

  /* set TSTOP (max time solution proceeds to) - this is not required */
  retval = CVodeSetStopTime(cvode_mem, t2);
  if (check_retval((void *)&retval, "CVodeSetStopTime", 1)) return(1);

  flag = RHS1; /* use -y for RHS */
  t = t1; /* set the integrator start time */

  printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));

  while (t<t2) {
    /* advance solver just one internal step */
    retval = CVode(cvode_mem, t2, y, &t, CV_ONE_STEP);
    if (check_retval((void *)&retval, "CVode", 1)) return(1);
    printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));
  }

  /* Get the number of steps the solver took after the discont. */
  retval = CVodeGetNumSteps(cvode_mem, &nst2);
  if (check_retval((void *)&retval, "CvodeGetNumSteps", 1)) return(1);

  /* Print statistics */
  nst = nst1 + nst2;
  printf("\nNumber of steps: %ld + %ld = %ld\n",nst1, nst2, nst);

  /*
   * ---------------------------------------------------------------
   * Discontinuity in RHS: Case 1 - explicit treatment
   * Note that it is not required to set TSTOP, but without it
   * we would have to find y(t1) to reinitialize the solver.
   * ---------------------------------------------------------------
   */

  printf("\nDiscontinuity in RHS: Case 1 - explicit treatment\n\n");

  /* Set initial condition */
  NV_Ith_S(y,0) = RCONST(1.0);

  /* Reinitialize the solver. CVodeReInit does not reallocate memory
   * so it can only be used when the new problem size is the same as
   * the problem size when CVodeCreate was called. */
  retval = CVodeReInit(cvode_mem, t0, y);
  if (check_retval((void *)&retval, "CVodeReInit", 1)) return(1);

  /* ---- Integrate to the discontinuity */

  /* Set TSTOP (max time solution proceeds to) to location of discont. */
  retval = CVodeSetStopTime(cvode_mem, t1);
  if (check_retval((void *)&retval, "CVodeSetStopTime", 1)) return(1);

  flag = RHS1; /* use -y for RHS */
  t = t0; /* set the integrator start time */

  printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));
  while (t<t1) {
    /* advance solver just one internal step */
    retval = CVode(cvode_mem, t1, y, &t, CV_ONE_STEP);
    if (check_retval((void *)&retval, "CVode", 1)) return(1);
    printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));
  }

  /* Get the number of steps the solver took to get to the discont. */
  retval = CVodeGetNumSteps(cvode_mem, &nst1);
  if (check_retval((void *)&retval, "CvodeGetNumSteps", 1)) return(1);

  /* If TSTOP was not set, we'd need to find y(t1): */
  /* CVodeGetDky(cvode_mem, t1, 0, y); */

  /* ---- Integrate from the discontinuity */

  /* Reinitialize solver */
  retval = CVodeReInit(cvode_mem, t1, y);

  /* set TSTOP (max time solution proceeds to) - this is not required */
  retval = CVodeSetStopTime(cvode_mem, t2);
  if (check_retval((void *)&retval, "CVodeSetStopTime", 1)) return(1);

  flag = RHS2; /* use -5y for RHS */
  t = t1; /* set the integrator start time */

  printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));

  while (t<t2) {
    /* advance solver just one internal step */
    retval = CVode(cvode_mem, t2, y, &t, CV_ONE_STEP);
    if (check_retval((void *)&retval, "CVode", 1)) return(1);
    printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));
  }

  /* Get the number of steps the solver took after the discont. */
  retval = CVodeGetNumSteps(cvode_mem, &nst2);
  if (check_retval((void *)&retval, "CvodeGetNumSteps", 1)) return(1);

  /* Print statistics */
  nst = nst1 + nst2;
  printf("\nNumber of steps: %ld + %ld = %ld\n",nst1, nst2, nst);


  /*
   * ---------------------------------------------------------------
   * Discontinuity in RHS: Case 2 - let CVODE deal with it
   * Note that here we MUST set TSTOP to ensure that the
   * change in the RHS happens at the appropriate time
   * ---------------------------------------------------------------
   */

  printf("\nDiscontinuity in RHS: Case 2 - let CVODE deal with it\n\n");

  /* Set initial condition */
  NV_Ith_S(y,0) = RCONST(1.0);

  /* Reinitialize the solver. CVodeReInit does not reallocate memory
   * so it can only be used when the new problem size is the same as
   * the problem size when CVodeCreate was called. */
  retval = CVodeReInit(cvode_mem, t0, y);
  if (check_retval((void *)&retval, "CVodeReInit", 1)) return(1);

  /* ---- Integrate to the discontinuity */

  /* Set TSTOP (max time solution proceeds to) to location of discont. */
  retval = CVodeSetStopTime(cvode_mem, t1);
  if (check_retval((void *)&retval, "CVodeSetStopTime", 1)) return(1);

  flag = RHS1; /* use -y for RHS */
  t = t0; /* set the integrator start time */

  printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));
  while (t<t1) {
    /* advance solver just one internal step */
    retval = CVode(cvode_mem, t1, y, &t, CV_ONE_STEP);
    if (check_retval((void *)&retval, "CVode", 1)) return(1);
    printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));
  }

  /* Get the number of steps the solver took to get to the discont. */
  retval = CVodeGetNumSteps(cvode_mem, &nst1);
  if (check_retval((void *)&retval, "CvodeGetNumSteps", 1)) return(1);

  /* ---- Integrate from the discontinuity */

  /* set TSTOP (max time solution proceeds to) - this is not required */
  retval = CVodeSetStopTime(cvode_mem, t2);
  if (check_retval((void *)&retval, "CVodeSetStopTime", 1)) return(1);

  flag = RHS2; /* use -5y for RHS */
  t = t1; /* set the integrator start time */

  printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));

  while (t<t2) {
    /* advance solver just one internal step */
    retval = CVode(cvode_mem, t2, y, &t, CV_ONE_STEP);
    if (check_retval((void *)&retval, "CVode", 1)) return(1);
    printf("%12.8"ESYM"  %12.8"ESYM"\n",t,NV_Ith_S(y,0));
  }

  /* Get the number of steps the solver took after the discont. */
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  if (check_retval((void *)&retval, "CvodeGetNumSteps", 1)) return(1);

  /* Print statistics */
  nst2 = nst - nst1;
  printf("\nNumber of steps: %ld + %ld = %ld\n",nst1, nst2, nst);

  /* Free memory */
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);

  return(0);
}

/*
 * RHS function
 * The form of the RHS function is controlled by the flag passed as f_data:
 *   flag = RHS1 -> y' = -y
 *   flag = RHS2 -> y' = -5*y
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  int *flag;

  flag = (int *) f_data;

  switch(*flag) {
  case RHS1:
    NV_Ith_S(ydot,0) = -NV_Ith_S(y,0);
    break;
  case RHS2:
    NV_Ith_S(ydot,0) = -5.0*NV_Ith_S(y,0);
    break;
  }

  return(0);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

