/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2020, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
 Routine to test the accuracy of the ARKodeInterp modules, using
 the problem
    [y'] = [ lambda*(y - sin(2*t)) + z - cos(3*t) + 2*cos(2*t) ]
    [z']   [     y - z - sin(2*t) + cos(3*t) - 3*sin(3*t)      ]
    [y(0); z(0)] = [0; 1]
 for various values of lambda<0.  This has analytical solution
    [y(t); z(t)] = [sin(2*t); cos(3*t)]
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#define NHVALS 9

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Private function to verify test results */
static int verify_results(realtype lambda, int rtype,
                          realtype *yrate, realtype *dyrate,
                          realtype *d2yrate, realtype *yrate2,
                          realtype *dyrate2, realtype *d2yrate2);

/* Main Program */
int main(int argc, char *argv[])
{
  /* initial declaration of all variables */
  realtype T0, Tf, lambda, t, t_test, hbase, rtol, atol;
  sunindextype NEQ;
  int flag, nttest, ideg, ih, itest, rtype;
  N_Vector y, ytest, dytest, d2ytest, d3ytest, d4ytest, d5ytest;
  N_Vector yerr, dyerr, d2yerr, d3yerr, d4yerr, d5yerr;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *arkode_mem;
  realtype hvals[NHVALS], yerrs[NHVALS], dyerrs[NHVALS], d2yerrs[NHVALS];
  realtype d3yerrs[NHVALS], d4yerrs[NHVALS], d5yerrs[NHVALS];
  realtype yrate[ARK_INTERP_MAX_DEGREE+1], dyrate[ARK_INTERP_MAX_DEGREE+1];
  realtype d2yrate[ARK_INTERP_MAX_DEGREE+1], d3yrate[ARK_INTERP_MAX_DEGREE+1];
  realtype d4yrate[ARK_INTERP_MAX_DEGREE+1], d5yrate[ARK_INTERP_MAX_DEGREE+1];
  realtype yferr[ARK_INTERP_MAX_DEGREE+1], dyferr[ARK_INTERP_MAX_DEGREE+1];
  realtype d2yferr[ARK_INTERP_MAX_DEGREE+1], d3yferr[ARK_INTERP_MAX_DEGREE+1];
  realtype d4yferr[ARK_INTERP_MAX_DEGREE+1], d5yferr[ARK_INTERP_MAX_DEGREE+1];
  realtype yrate2[ARK_INTERP_MAX_DEGREE+1], dyrate2[ARK_INTERP_MAX_DEGREE+1];
  realtype d2yrate2[ARK_INTERP_MAX_DEGREE+1], yferr2[ARK_INTERP_MAX_DEGREE+1];
  realtype dyferr2[ARK_INTERP_MAX_DEGREE+1], d2yferr2[ARK_INTERP_MAX_DEGREE+1];

  /* general problem parameters */
  T0 = RCONST(0.0);       /* initial time */
  Tf = RCONST(10.0);      /* final time */
  NEQ = 2;                /* number of dependent vars. */

  /* if an argument supplied, set lambda (otherwise use -100) */
  lambda = -RCONST(100.0);
  if (argc > 1)  lambda = strtod(argv[1], NULL);

  /* determine test configuration */
  if (sizeof(realtype) == 4) {
    rtype = 32;
  } else if (sizeof(realtype) == 8) {
    rtype = 64;
  } else {
    rtype = 128;
  }

  /* set IVP solver tolerances based on rtype */
  if (rtype == 32) {
    rtol = 1e-7;
    atol = 1e-7;
  } else if (rtype == 64) {
    rtol = 1e-14;
    atol = 1e-14;
  } else {
    rtol = 1e-19;
    atol = 1e-20;
  }

  /* initialize vector/solver objects to NULL */
  y = ytest = dytest = d2ytest = d3ytest = d4ytest = d5ytest = NULL;
  yerr = dyerr = d2yerr = d3yerr = d4yerr = d5yerr = NULL;
  A = NULL;
  LS = NULL;
  arkode_mem = NULL;

  /* Initial problem output */
  printf("\nARKode Interpolation module tester (rtype = %i):\n", rtype);

  /* Initialize vectors */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  ytest = N_VNew_Serial(NEQ);
  if (check_flag((void *)ytest, "N_VNew_Serial", 0)) return 1;
  dytest = N_VNew_Serial(NEQ);
  if (check_flag((void *)dytest, "N_VNew_Serial", 0)) return 1;
  d2ytest = N_VNew_Serial(NEQ);
  if (check_flag((void *)d2ytest, "N_VNew_Serial", 0)) return 1;
  d3ytest = N_VNew_Serial(NEQ);
  if (check_flag((void *)d3ytest, "N_VNew_Serial", 0)) return 1;
  d4ytest = N_VNew_Serial(NEQ);
  if (check_flag((void *)d4ytest, "N_VNew_Serial", 0)) return 1;
  d5ytest = N_VNew_Serial(NEQ);
  if (check_flag((void *)d5ytest, "N_VNew_Serial", 0)) return 1;
  yerr = N_VNew_Serial(NEQ);
  if (check_flag((void *)yerr, "N_VNew_Serial", 0)) return 1;
  dyerr = N_VNew_Serial(NEQ);
  if (check_flag((void *)dyerr, "N_VNew_Serial", 0)) return 1;
  d2yerr = N_VNew_Serial(NEQ);
  if (check_flag((void *)d2yerr, "N_VNew_Serial", 0)) return 1;
  d3yerr = N_VNew_Serial(NEQ);
  if (check_flag((void *)d3yerr, "N_VNew_Serial", 0)) return 1;
  d4yerr = N_VNew_Serial(NEQ);
  if (check_flag((void *)d4yerr, "N_VNew_Serial", 0)) return 1;
  d5yerr = N_VNew_Serial(NEQ);
  if (check_flag((void *)d5yerr, "N_VNew_Serial", 0)) return 1;

  /* initialize algebraic solver structures */
  A = SUNDenseMatrix(NEQ, NEQ);
  if (check_flag((void *) A, "SUNDenseMatrix", 0)) return 1;
  LS = SUNDenseLinearSolver(y, A);
  if (check_flag((void *) LS, "SUNDenseLinearSolver", 0)) return 1;

  /* test parameters */
  nttest = 500;
  hbase = RCONST(2.0);

  /*---- Part I: Hermite interpolation module ----*/
  printf("\nHermite interpolation module tests:\n");

  /* loop over dense output polynomial degrees */
  for (ideg=0; ideg<=ARK_INTERP_MAX_DEGREE; ideg++) {

    /* reset error/convergence arrays */
    for (ih=0; ih<NHVALS; ih++) {
      yerrs[ih]   = RCONST(0.0);
      dyerrs[ih]  = RCONST(0.0);
      d2yerrs[ih] = RCONST(0.0);
      d3yerrs[ih] = RCONST(0.0);
      d4yerrs[ih] = RCONST(0.0);
      d5yerrs[ih] = RCONST(0.0);
    }

    /* run tests for this dense output polynomial degree */
    for (ih=0; ih<NHVALS; ih++) {

      /* set test stepsize */
      hvals[ih] = RCONST(2.0)/pow(hbase,ih);

      /* Initialize y values */
      NV_Ith_S(y,0) = RCONST(0.0);
      NV_Ith_S(y,1) = RCONST(1.0);

      /* Create/initialize ARKStep module */
      arkode_mem = ARKStepCreate(NULL, f, T0, y);
      if (check_flag(arkode_mem, "ARKStepCreate", 0)) return 1;

      /* pass lambda to RHS routine */
      flag = ARKStepSetUserData(arkode_mem, &lambda);
      if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

      /* select Hermite interpolation module */
      flag = ARKStepSetInterpolantType(arkode_mem, ARK_INTERP_HERMITE);
      if (check_flag(&flag, "ARKStepSetInterpolantType", 1)) return 1;

      /* set dense output polynomial degree */
      flag = ARKStepSetInterpolantDegree(arkode_mem, ideg);
      if (check_flag(&flag, "ARKStepSetInterpolantDegree", 1)) return 1;

      /* set fixed time-stepping with desired time step size */
      flag = ARKStepSetFixedStep(arkode_mem, hvals[ih]);
      if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;

      /* set solver tolerances */
      flag = ARKStepSStolerances(arkode_mem, rtol, atol);
      if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

      /* indicate linearity of problem */
      flag = ARKStepSetLinear(arkode_mem, 0);
      if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;

      /* attach linear solver */
      flag = ARKStepSetLinearSolver(arkode_mem, LS, A);
      if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;

      /* increase maximum number of time steps */
      flag = ARKStepSetMaxNumSteps(arkode_mem, 100000);
      if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;

      /* set RK order to highest available value */
      flag = ARKStepSetOrder(arkode_mem, 5);
      if (check_flag(&flag, "ARKStepSetOrder", 1)) return 1;

      /* evolve to Tf to prepare interpolation structure */
      flag = ARKStepSetStopTime(arkode_mem, Tf);
      if (check_flag(&flag, "ARKStepSetStopTime", 1)) return 1;
      flag = ARKStepEvolve(arkode_mem, Tf, y, &t, ARK_NORMAL);
      if (check_flag(&flag, "ARKStepEvolve", 1)) return 1;

      /* loop over 100 evenly-spaced values within interior of preceding
         step to accumulate errors */
      for (itest=0; itest<nttest; itest++) {

        /* set test time */
        t_test = t - hvals[ih] + (itest+1)*hvals[ih]/(nttest + 2);

        /* call ARKStepGetDky to evaluate solution and derivatives at t_test */
        flag = ARKStepGetDky(arkode_mem, t_test, 0, ytest);
        if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;
        flag = ARKStepGetDky(arkode_mem, t_test, 1, dytest);
        if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;
        flag = ARKStepGetDky(arkode_mem, t_test, 2, d2ytest);
        if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;
        flag = ARKStepGetDky(arkode_mem, t_test, 3, d3ytest);
        if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;
        flag = ARKStepGetDky(arkode_mem, t_test, 4, d4ytest);
        if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;
        flag = ARKStepGetDky(arkode_mem, t_test, 5, d5ytest);
        if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;

        /* set error values */
        /*   y */
        NV_Ith_S(yerr,0) = SUNRabs(sin(RCONST(2.0)*t_test) - NV_Ith_S(ytest,0));
        NV_Ith_S(yerr,1) = SUNRabs(cos(RCONST(3.0)*t_test) - NV_Ith_S(ytest,1));
        /*   dy */
        NV_Ith_S(dyerr,0) = SUNRabs(RCONST(2.0)*cos(RCONST(2.0)*t_test) - NV_Ith_S(dytest,0));
        NV_Ith_S(dyerr,1) = SUNRabs(-RCONST(3.0)*sin(RCONST(3.0)*t_test) - NV_Ith_S(dytest,1));
        /*   d2y */
        NV_Ith_S(d2yerr,0) = SUNRabs(-RCONST(4.0)*sin(RCONST(2.0)*t_test) - NV_Ith_S(d2ytest,0));
        NV_Ith_S(d2yerr,1) = SUNRabs(-RCONST(9.0)*cos(RCONST(3.0)*t_test) - NV_Ith_S(d2ytest,1));
        /*   d3y */
        NV_Ith_S(d3yerr,0) = SUNRabs(-RCONST(8.0)*cos(RCONST(2.0)*t_test) - NV_Ith_S(d3ytest,0));
        NV_Ith_S(d3yerr,1) = SUNRabs(RCONST(27.0)*sin(RCONST(3.0)*t_test) - NV_Ith_S(d3ytest,1));
        /*   d4y */
        NV_Ith_S(d4yerr,0) = SUNRabs(RCONST(16.0)*sin(RCONST(2.0)*t_test) - NV_Ith_S(d4ytest,0));
        NV_Ith_S(d4yerr,1) = SUNRabs(RCONST(81.0)*cos(RCONST(3.0)*t_test) - NV_Ith_S(d4ytest,1));
        /*   d5y */
        NV_Ith_S(d5yerr,0) = SUNRabs(RCONST(32.0)*cos(RCONST(2.0)*t_test) - NV_Ith_S(d5ytest,0));
        NV_Ith_S(d5yerr,1) = SUNRabs(-RCONST(243.0)*sin(RCONST(3.0)*t_test) - NV_Ith_S(d5ytest,1));

        /* compute error norms (2-norm per test, max-norm over interval) */
        yerrs[ih] = SUNMAX(yerrs[ih], sqrt(N_VDotProd(yerr,yerr)));
        dyerrs[ih] = SUNMAX(dyerrs[ih], sqrt(N_VDotProd(dyerr,dyerr)));
        d2yerrs[ih] = SUNMAX(d2yerrs[ih], sqrt(N_VDotProd(d2yerr,d2yerr)));
        d3yerrs[ih] = SUNMAX(d3yerrs[ih], sqrt(N_VDotProd(d3yerr,d3yerr)));
        d4yerrs[ih] = SUNMAX(d4yerrs[ih], sqrt(N_VDotProd(d4yerr,d4yerr)));
        d5yerrs[ih] = SUNMAX(d5yerrs[ih], sqrt(N_VDotProd(d5yerr,d5yerr)));

      }  /* end itest loop */

      /* free ARKStep memory (to prepare for next call) */
      ARKStepFree(&arkode_mem);
      arkode_mem = NULL;

    }  /* end ih loop */

    printf("\nConvergence Test Results, dense output polynomial degree = %i:\n", ideg);
    printf("\n  Raw errors per h value:\n");
    printf("  ------------------------------------------------------------------------------\n");
    printf("    h         yerr       dyerr      d2yerr     d3yerr    d4yerr   d5yerr\n");
    printf("  ------------------------------------------------------------------------------\n");
    for (ih=0; ih<NHVALS; ih++)
      printf("    %.1e   %.2e   %.2e   %.2e   %.2e   %.2e   %.2e\n", (double) hvals[ih],
             (double) yerrs[ih], (double) dyerrs[ih], (double) d2yerrs[ih],
             (double) d3yerrs[ih], (double) d4yerrs[ih], (double) d5yerrs[ih]);
    printf("  ------------------------------------------------------------------------------\n");

    printf("  Estimated y convergence factors:\n  ");
    for (ih=1; ih<NHVALS; ih++)
      printf("  %.3f", log(yerrs[ih]/yerrs[ih-1])/log(hvals[ih]/hvals[ih-1]));
    yrate[ideg] = log(yerrs[NHVALS-1]/yerrs[0])/log(hvals[NHVALS-1]/hvals[0]);
    yferr[ideg] = yerrs[NHVALS-1];
    printf("  (%.3f, %.0e)\n\n", (double) yrate[ideg], (double) yferr[ideg]);

    printf("  Estimated dy convergence factors:\n  ");
    for (ih=1; ih<NHVALS; ih++)
      printf("  %.3f", log(dyerrs[ih]/dyerrs[ih-1])/log(hvals[ih]/hvals[ih-1]));
    dyrate[ideg] = log(dyerrs[NHVALS-1]/dyerrs[0])/log(hvals[NHVALS-1]/hvals[0]);
    dyferr[ideg] = dyerrs[NHVALS-1];
    printf("  (%.3f, %.0e)\n\n", (double) dyrate[ideg], (double) dyferr[ideg]);

    printf("  Estimated d2y convergence factors:\n  ");
    for (ih=1; ih<NHVALS; ih++)
      printf("  %.3f", log(d2yerrs[ih]/d2yerrs[ih-1])/log(hvals[ih]/hvals[ih-1]));
    d2yrate[ideg] = log(d2yerrs[NHVALS-1]/d2yerrs[0])/log(hvals[NHVALS-1]/hvals[0]);
    d2yferr[ideg] = d2yerrs[NHVALS-1];
    printf("  (%.3f, %.0e)\n\n", (double) d2yrate[ideg], (double) d2yferr[ideg]);

    printf("  Estimated d3y convergence factors:\n  ");
    for (ih=1; ih<NHVALS; ih++)
      printf("  %.3f", log(d3yerrs[ih]/d3yerrs[ih-1])/log(hvals[ih]/hvals[ih-1]));
    d3yrate[ideg] = log(d3yerrs[NHVALS-1]/d3yerrs[0])/log(hvals[NHVALS-1]/hvals[0]);
    d3yferr[ideg] = d3yerrs[NHVALS-1];
    printf("  (%.3f, %.0e)\n\n", (double) d3yrate[ideg], (double) d3yferr[ideg]);

    printf("  Estimated d4y convergence factors:\n  ");
    for (ih=1; ih<NHVALS; ih++)
      printf("  %.3f", log(d4yerrs[ih]/d4yerrs[ih-1])/log(hvals[ih]/hvals[ih-1]));
    d4yrate[ideg] = log(d4yerrs[NHVALS-1]/d4yerrs[0])/log(hvals[NHVALS-1]/hvals[0]);
    d4yferr[ideg] = d4yerrs[NHVALS-1];
    printf("  (%.3f, %.0e)\n\n", (double) d4yrate[ideg], (double) d4yferr[ideg]);

    printf("  Estimated d5y convergence factors:\n  ");
    for (ih=1; ih<NHVALS; ih++)
      printf("  %.3f", log(d5yerrs[ih]/d5yerrs[ih-1])/log(hvals[ih]/hvals[ih-1]));
    d5yrate[ideg] = log(d5yerrs[NHVALS-1]/d5yerrs[0])/log(hvals[NHVALS-1]/hvals[0]);
    d5yferr[ideg] = d5yerrs[NHVALS-1];
    printf("  (%.3f, %.0e)\n\n", (double) d5yrate[ideg], (double) d5yferr[ideg]);

  } /* end ideg loop */



  /*---- Part II: Lagrange interpolation module ----*/

  printf("\nLagrange interpolation module tests:\n");

  /* loop over dense output polynomial degrees */
  for (ideg=0; ideg<=ARK_INTERP_MAX_DEGREE; ideg++) {

    /* reset error/convergence arrays */
    for (ih=0; ih<NHVALS; ih++) {
      yerrs[ih] = RCONST(0.0);
      dyerrs[ih] = RCONST(0.0);
      d2yerrs[ih] = RCONST(0.0);
    }

    /* run tests for this dense output polynomial degree */
    for (ih=0; ih<NHVALS; ih++) {

      /* set test stepsize */
      hvals[ih] = RCONST(2.0)/pow(hbase,ih);

      /* Initialize y values */
      NV_Ith_S(y,0) = RCONST(0.0);
      NV_Ith_S(y,1) = RCONST(1.0);

      /* Create/initialize ARKStep module */
      arkode_mem = ARKStepCreate(NULL, f, T0, y);
      if (check_flag(arkode_mem, "ARKStepCreate", 0)) return 1;

      /* pass lambda to RHS routine */
      flag = ARKStepSetUserData(arkode_mem, &lambda);
      if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;

      /* select Lagrange interpolation module */
      flag = ARKStepSetInterpolantType(arkode_mem, ARK_INTERP_LAGRANGE);
      if (check_flag(&flag, "ARKStepSetInterpolantType", 1)) return 1;

      /* set dense output polynomial degree */
      flag = ARKStepSetInterpolantDegree(arkode_mem, ideg);
      if (check_flag(&flag, "ARKStepSetInterpolantDegree", 1)) return 1;

      /* set fixed time-stepping with desired time step size */
      flag = ARKStepSetFixedStep(arkode_mem, hvals[ih]);
      if (check_flag(&flag, "ARKStepSetFixedStep", 1)) return 1;

      /* set solver tolerances */
      flag = ARKStepSStolerances(arkode_mem, rtol, atol);
      if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

      /* indicate linearity of problem */
      flag = ARKStepSetLinear(arkode_mem, 0);
      if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;

      /* attach linear solver */
      flag = ARKStepSetLinearSolver(arkode_mem, LS, A);
      if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;

      /* increase maximum number of time steps */
      flag = ARKStepSetMaxNumSteps(arkode_mem, 100000);
      if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;

      /* set RK order to highest available value */
      flag = ARKStepSetOrder(arkode_mem, 5);
      if (check_flag(&flag, "ARKStepSetOrder", 1)) return 1;

      /* evolve to Tf to prepare interpolation structure */
      flag = ARKStepSetStopTime(arkode_mem, Tf);
      if (check_flag(&flag, "ARKStepSetStopTime", 1)) return 1;
      flag = ARKStepEvolve(arkode_mem, Tf, y, &t, ARK_NORMAL);
      if (check_flag(&flag, "ARKStepEvolve", 1)) return 1;

      /* loop over 100 evenly-spaced values within interior of this step to accumulate errors */
      for (itest=0; itest<nttest; itest++) {

        /* set test time */
        t_test = t - hvals[ih] + (itest+1)*hvals[ih]/(nttest + 2);

        /* call ARKStepGetDky to evaluate solution and derivatives at t_test */
        flag = ARKStepGetDky(arkode_mem, t_test, 0, ytest);
        if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;
        flag = ARKStepGetDky(arkode_mem, t_test, 1, dytest);
        if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;
        flag = ARKStepGetDky(arkode_mem, t_test, 2, d2ytest);
        if (check_flag(&flag, "ARKStepGetDky", 1)) return 1;

        /* set error values */
        /*   y */
        NV_Ith_S(yerr,0) = SUNRabs(sin(RCONST(2.0)*t_test) - NV_Ith_S(ytest,0));
        NV_Ith_S(yerr,1) = SUNRabs(cos(RCONST(3.0)*t_test) - NV_Ith_S(ytest,1));
        /*   dy */
        NV_Ith_S(dyerr,0) = SUNRabs(RCONST(2.0)*cos(RCONST(2.0)*t_test) - NV_Ith_S(dytest,0));
        NV_Ith_S(dyerr,1) = SUNRabs(-RCONST(3.0)*sin(RCONST(3.0)*t_test) - NV_Ith_S(dytest,1));
        /*   d2y */
        NV_Ith_S(d2yerr,0) = SUNRabs(-RCONST(4.0)*sin(RCONST(2.0)*t_test) - NV_Ith_S(d2ytest,0));
        NV_Ith_S(d2yerr,1) = SUNRabs(-RCONST(9.0)*cos(RCONST(3.0)*t_test) - NV_Ith_S(d2ytest,1));

        /* compute error norms (2-norm per test, max-norm over interval) */
        yerrs[ih] = SUNMAX(yerrs[ih], sqrt(N_VDotProd(yerr,yerr)));
        dyerrs[ih] = SUNMAX(dyerrs[ih], sqrt(N_VDotProd(dyerr,dyerr)));
        d2yerrs[ih] = SUNMAX(d2yerrs[ih], sqrt(N_VDotProd(d2yerr,d2yerr)));

      }  /* end itest loop */

      /* free ARKStep memory (to prepare for next call) */
      ARKStepFree(&arkode_mem);
      arkode_mem = NULL;

    }  /* end ih loop */

    printf("\nConvergence Test Results, dense output polynomial degree = %i:\n", ideg);
    printf("\n  Raw errors per h value:\n");
    printf("  ------------------------------------------\n");
    printf("    h         yerr       dyerr      d2yerr\n");
    printf("  ------------------------------------------\n");
    for (ih=0; ih<NHVALS; ih++)
      printf("    %.1e   %.2e   %.2e   %.2e\n", (double) hvals[ih],
             (double) yerrs[ih], (double) dyerrs[ih], (double) d2yerrs[ih]);
    printf("  ------------------------------------------\n");

    printf("  Estimated y convergence factors:\n  ");
    for (ih=1; ih<NHVALS; ih++)
      printf("  %.3f", log(yerrs[ih]/yerrs[ih-1])/log(hvals[ih]/hvals[ih-1]));
    yrate2[ideg] = log(yerrs[NHVALS-1]/yerrs[0])/log(hvals[NHVALS-1]/hvals[0]);
    yferr2[ideg] = yerrs[NHVALS-1];
    printf("  (%.3f, %.0e)\n\n", (double) yrate2[ideg], (double) yferr2[ideg]);

    printf("  Estimated dy convergence factors:\n  ");
    for (ih=1; ih<NHVALS; ih++)
      printf("  %.3f", log(dyerrs[ih]/dyerrs[ih-1])/log(hvals[ih]/hvals[ih-1]));
    dyrate2[ideg] = log(dyerrs[NHVALS-1]/dyerrs[0])/log(hvals[NHVALS-1]/hvals[0]);
    dyferr2[ideg] = dyerrs[NHVALS-1];
    printf("  (%.3f, %.0e)\n\n", (double) dyrate2[ideg], (double) dyferr2[ideg]);

    printf("  Estimated d2y convergence factors:\n  ");
    for (ih=1; ih<NHVALS; ih++)
      printf("  %.3f", log(d2yerrs[ih]/d2yerrs[ih-1])/log(hvals[ih]/hvals[ih-1]));
    d2yrate2[ideg] = log(d2yerrs[NHVALS-1]/d2yerrs[0])/log(hvals[NHVALS-1]/hvals[0]);
    d2yferr2[ideg] = d2yerrs[NHVALS-1];
    printf("  (%.3f, %.0e)\n\n", (double) d2yrate2[ideg], (double) d2yferr2[ideg]);

  } /* end ideg loop */


  /*---- Part III: output overall convergence rate tables ----*/

  printf("\n\nCumulative convergence rates (final errors):\n");
  printf("\n  Hermite Interpolation module:\n");
  printf("  ---------------------------------------------------------------------------------------------\n");
  printf("   degree |     y             dy           d2y           d3y         d4y         d5y\n");
  printf("  ---------------------------------------------------------------------------------------------\n");
  for (ideg=0; ideg<=ARK_INTERP_MAX_DEGREE; ideg++)
      printf("      %1d   |  %.1f (%.0e)   %.1f (%.0e)   %.1f (%.0e)   %.1f (%.0e)   %.1f (%.0e)   %.1f (%.0e)\n",
             ideg, (double) yrate[ideg], (double) yferr[ideg], (double) dyrate[ideg], (double) dyferr[ideg],
             (double) d2yrate[ideg], (double) d2yferr[ideg], (double) d3yrate[ideg], (double) d3yferr[ideg],
             (double) d4yrate[ideg], (double) d4yferr[ideg], (double) d5yrate[ideg], (double) d5yferr[ideg]);
  printf("  ---------------------------------------------------------------------------------------------\n");

  printf("\n  Lagrange Interpolation module:\n");
  printf("  ---------------------------------------------------\n");
  printf("   degree |     y             dy           d2y\n");
  printf("  ---------------------------------------------------\n");
  for (ideg=0; ideg<=ARK_INTERP_MAX_DEGREE; ideg++)
      printf("      %1d   |  %.1f (%.0e)   %.1f (%.0e)   %.1f (%.0e)\n", ideg,
             (double) yrate2[ideg], (double) yferr2[ideg], (double) dyrate2[ideg],
             (double) dyferr2[ideg], (double) d2yrate2[ideg], (double) d2yferr2[ideg]);
  printf("  ---------------------------------------------------\n");

  /* check validity of results */
  flag = verify_results(lambda, rtype, yrate, dyrate, d2yrate, yrate2, dyrate2, d2yrate2);

  /* Clean up and return 'flag' indicating success/failure of tests */
  N_VDestroy(y);
  N_VDestroy(ytest);
  N_VDestroy(dytest);
  N_VDestroy(d2ytest);
  N_VDestroy(d3ytest);
  N_VDestroy(d4ytest);
  N_VDestroy(d5ytest);
  N_VDestroy(yerr);
  N_VDestroy(dyerr);
  N_VDestroy(d2yerr);
  N_VDestroy(d3yerr);
  N_VDestroy(d4yerr);
  N_VDestroy(d5yerr);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  return(flag);
}


/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *lambda = (realtype *) user_data;
  NV_Ith_S(ydot,0) = (*lambda)*(NV_Ith_S(y,0) - sin(RCONST(2.0)*t))
    + NV_Ith_S(y,1) - cos(RCONST(3.0)*t) + RCONST(2.0)*cos(RCONST(2.0)*t);
  NV_Ith_S(ydot,1) = NV_Ith_S(y,0) - NV_Ith_S(y,1) - sin(RCONST(2.0)*t)
    + cos(RCONST(3.0)*t) - RCONST(3.0)*sin(RCONST(3.0)*t);
  return 0;  /* Return with success */
}



/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
   opt == 0 means SUNDIALS function allocates memory so check if
   returned NULL pointer
   opt == 1 means SUNDIALS function returns a flag so check if
   flag >= 0
   opt == 2 means function allocates memory so check if returned
   NULL pointer
*/
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


/* Verify convergence results...
   We only check convergence rates for y, y', and y'', for both the
   Hermite and Lagrange interpolation modules.  All of the benchmark
   results below were computed/entered manually, by running each
   configuration.  Any convergence rates corresponding to 'final'
   error values >1 are omitted from these checks, as non-convergence
   may occur differently on various architectures.

   We consider realtype precisions {32,64,128} and stiffness values
   lambda = {-1e2, -1e4, -1e6}.
*/
static int verify_results(realtype lambda, int rtype,
                          realtype *yrate, realtype *dyrate,
                          realtype *d2yrate, realtype *yrate2,
                          realtype *dyrate2, realtype *d2yrate2)
{
  int tests=0;
  int failure=0;

  /* single precision checks */
  if (rtype == 32) {
    if (SUNRabs(lambda) < 1.1e2) {        /* lambda = -1e2 */
      /* Hermite module convergence */
      failure += (yrate[0] < 0.8);  tests++;
      failure += (yrate[1] < 1.8);  tests++;
      failure += (yrate[2] < 2.2);  tests++;
      failure += (yrate[3] < 2.0);  tests++;
      failure += (yrate[4] < 2.7);  tests++;
      failure += (yrate[5] < 2.7);  tests++;
      failure += (dyrate[1] < 0.9);  tests++;
      failure += (dyrate[2] < 1.7);  tests++;
      failure += (dyrate[3] < 1.8);  tests++;
      failure += (dyrate[4] < 2.3);  tests++;
      failure += (dyrate[5] < 2.3);  tests++;
      failure += (d2yrate[2] < 0.7);  tests++;
      failure += (d2yrate[3] < 0.9);  tests++;
      failure += (d2yrate[4] < 1.3);  tests++;
      failure += (d2yrate[5] < 1.3);  tests++;
      /* Lagrange module convergence */
      failure += (yrate2[0] < 0.7);  tests++;
      failure += (yrate2[1] < 1.8);  tests++;
      failure += (yrate2[2] < 2.1);  tests++;
      failure += (yrate2[3] < 2.1);  tests++;
      failure += (yrate2[4] < 2.1);  tests++;
      failure += (yrate2[5] < 2.1);  tests++;
      failure += (dyrate2[1] < 0.9);  tests++;
      failure += (dyrate2[2] < 1.5);  tests++;
      failure += (dyrate2[3] < 1.9);  tests++;
      failure += (dyrate2[4] < 1.8);  tests++;
      failure += (dyrate2[5] < 1.8);  tests++;
      failure += (d2yrate2[2] < 0.6);  tests++;
      failure += (d2yrate2[3] < 1.2);  tests++;
      failure += (d2yrate2[4] < 0.9);  tests++;
      failure += (d2yrate2[5] < 0.9);  tests++;
    } else if (SUNRabs(lambda) < 1.1e4) { /* lambda = -1e4 */
      /* Hermite module convergence */
      failure += (yrate[0] < 0.8);  tests++;
      failure += (yrate[1] < 1.8);  tests++;
      failure += (yrate[2] < 2.2);  tests++;
      failure += (yrate[3] < 2.0);  tests++;
      failure += (yrate[4] < 3.2);  tests++;
      failure += (yrate[5] < 3.2);  tests++;
      failure += (dyrate[1] < 0.9);  tests++;
      failure += (dyrate[2] < 1.6);  tests++;
      failure += (dyrate[3] < 1.2);  tests++;
      failure += (dyrate[4] < 2.2);  tests++;
      failure += (dyrate[5] < 2.2);  tests++;
      failure += (d2yrate[2] < 0.7);  tests++;
      failure += (d2yrate[3] < 0.5);  tests++;
      /* Lagrange module convergence */
      failure += (yrate2[0] < 0.7);  tests++;
      failure += (yrate2[1] < 1.8);  tests++;
      failure += (yrate2[2] < 2.1);  tests++;
      failure += (yrate2[3] < 2.1);  tests++;
      failure += (yrate2[4] < 2.1);  tests++;
      failure += (yrate2[5] < 2.1);  tests++;
      failure += (dyrate2[1] < 0.9);  tests++;
      failure += (dyrate2[2] < 1.5);  tests++;
      failure += (dyrate2[3] < 1.6);  tests++;
      failure += (dyrate2[4] < 1.5);  tests++;
      failure += (dyrate2[5] < 1.5);  tests++;
      failure += (d2yrate2[2] < 0.6);  tests++;
      failure += (d2yrate2[3] < 0.8);  tests++;
      failure += (d2yrate2[4] < 0.7);  tests++;
      failure += (d2yrate2[5] < 0.7);  tests++;
    } else {                              /* lambda = -1e6 */
      /* Hermite module convergence */
      failure += (yrate[0] < 0.8);  tests++;
      failure += (yrate[1] < 1.7);  tests++;
      failure += (yrate[2] < 1.3);  tests++;
      failure += (yrate[3] < 1.1);  tests++;
      failure += (dyrate[1] < 0.9);  tests++;
      failure += (dyrate[2] < 0.4);  tests++;
      failure += (dyrate[3] < 0.1);  tests++;
      /* Lagrange module convergence */
      failure += (yrate2[0] < 0.7);  tests++;
      failure += (yrate2[1] < 1.7);  tests++;
      failure += (yrate2[2] < 1.7);  tests++;
      failure += (yrate2[3] < 1.7);  tests++;
      failure += (yrate2[4] < 1.6);  tests++;
      failure += (yrate2[5] < 1.6);  tests++;
      failure += (dyrate2[1] < 0.9);  tests++;
      failure += (dyrate2[2] < 1.1);  tests++;
      failure += (dyrate2[3] < 0.9);  tests++;
      failure += (dyrate2[4] < 0.8);  tests++;
      failure += (dyrate2[5] < 0.8);  tests++;
    }
  }

  /* double precision checks */
  if (rtype == 64) {
    if (SUNRabs(lambda) < 1.1e2) {        /* lambda = -1e2 */
      /* Hermite module convergence */
      failure += (yrate[0] < 0.8);  tests++;
      failure += (yrate[1] < 1.8);  tests++;
      failure += (yrate[2] < 2.6);  tests++;
      failure += (yrate[3] < 2.4);  tests++;
      failure += (yrate[4] < 3.1);  tests++;
      failure += (yrate[5] < 3.1);  tests++;
      failure += (dyrate[1] < 0.9);  tests++;
      failure += (dyrate[2] < 1.7);  tests++;
      failure += (dyrate[3] < 1.8);  tests++;
      failure += (dyrate[4] < 2.3);  tests++;
      failure += (dyrate[5] < 2.3);  tests++;
      failure += (d2yrate[2] < 0.7);  tests++;
      failure += (d2yrate[3] < 1.0);  tests++;
      failure += (d2yrate[4] < 1.3);  tests++;
      failure += (d2yrate[5] < 1.3);  tests++;
      /* Lagrange module convergence */
      failure += (yrate2[0] < 0.7);  tests++;
      failure += (yrate2[1] < 1.8);  tests++;
      failure += (yrate2[2] < 2.4);  tests++;
      failure += (yrate2[3] < 2.5);  tests++;
      failure += (yrate2[4] < 2.6);  tests++;
      failure += (yrate2[5] < 2.6);  tests++;
      failure += (dyrate2[1] < 0.9);  tests++;
      failure += (dyrate2[2] < 1.5);  tests++;
      failure += (dyrate2[3] < 2.4);  tests++;
      failure += (dyrate2[4] < 2.5);  tests++;
      failure += (dyrate2[5] < 2.5);  tests++;
      failure += (d2yrate2[2] < 0.6);  tests++;
      failure += (d2yrate2[3] < 1.6);  tests++;
      failure += (d2yrate2[4] < 2.0);  tests++;
      failure += (d2yrate2[5] < 2.0);  tests++;
    } else if (SUNRabs(lambda) < 1.1e4) { /* lambda = -1e4 */
      /* Hermite module convergence */
      failure += (yrate[0] < 0.8);  tests++;
      failure += (yrate[1] < 1.8);  tests++;
      failure += (yrate[2] < 2.7);  tests++;
      failure += (yrate[3] < 2.5);  tests++;
      failure += (yrate[4] < 3.5);  tests++;
      failure += (yrate[5] < 3.5);  tests++;
      failure += (dyrate[1] < 0.9);  tests++;
      failure += (dyrate[2] < 1.7);  tests++;
      failure += (dyrate[3] < 1.6);  tests++;
      failure += (dyrate[4] < 2.6);  tests++;
      failure += (dyrate[5] < 2.6);  tests++;
      failure += (d2yrate[2] < 0.7);  tests++;
      failure += (d2yrate[3] < 0.7);  tests++;
      /* Lagrange module convergence */
      failure += (yrate2[0] < 0.7);  tests++;
      failure += (yrate2[1] < 1.8);  tests++;
      failure += (yrate2[2] < 2.4);  tests++;
      failure += (yrate2[3] < 2.7);  tests++;
      failure += (yrate2[4] < 2.7);  tests++;
      failure += (yrate2[5] < 2.7);  tests++;
      failure += (dyrate2[1] < 0.9);  tests++;
      failure += (dyrate2[2] < 1.5);  tests++;
      failure += (dyrate2[3] < 2.5);  tests++;
      failure += (dyrate2[4] < 2.6);  tests++;
      failure += (dyrate2[5] < 2.6);  tests++;
      failure += (d2yrate2[2] < 0.6);  tests++;
      failure += (d2yrate2[3] < 1.6);  tests++;
      failure += (d2yrate2[4] < 2.0);  tests++;
      failure += (d2yrate2[5] < 2.0);  tests++;
    } else {                              /* lambda = -1e6 */
      /* Hermite module convergence */
      failure += (yrate[0] < 0.8);  tests++;
      failure += (yrate[1] < 1.8);  tests++;
      failure += (yrate[2] < 2.7);  tests++;
      failure += (yrate[3] < 2.4);  tests++;
      failure += (yrate[4] < 3.6);  tests++;
      failure += (yrate[5] < 3.6);  tests++;
      failure += (dyrate[1] < 0.9);  tests++;
      failure += (dyrate[2] < 1.6);  tests++;
      failure += (dyrate[3] < 1.4);  tests++;
      failure += (dyrate[4] < 2.6);  tests++;
      failure += (dyrate[5] < 2.6);  tests++;
      failure += (d2yrate[2] < 0.7);  tests++;
      failure += (d2yrate[3] < 0.6);  tests++;
      /* Lagrange module convergence */
      failure += (yrate2[0] < 0.7);  tests++;
      failure += (yrate2[1] < 1.8);  tests++;
      failure += (yrate2[2] < 2.4);  tests++;
      failure += (yrate2[3] < 2.6);  tests++;
      failure += (yrate2[4] < 2.7);  tests++;
      failure += (yrate2[5] < 2.7);  tests++;
      failure += (dyrate2[1] < 0.9);  tests++;
      failure += (dyrate2[2] < 1.5);  tests++;
      failure += (dyrate2[3] < 2.5);  tests++;
      failure += (dyrate2[4] < 2.6);  tests++;
      failure += (dyrate2[5] < 2.6);  tests++;
      failure += (d2yrate2[2] < 0.6);  tests++;
      failure += (d2yrate2[3] < 1.6);  tests++;
      failure += (d2yrate2[4] < 2.0);  tests++;
      failure += (d2yrate2[5] < 2.0);  tests++;
    }
  }

  /* extended precision checks */
  if (rtype == 128) {
    if (SUNRabs(lambda) < 1.1e2) {        /* lambda = -1e2 */
      /* Hermite module convergence */
      failure += (yrate[0] < 0.8);  tests++;
      failure += (yrate[1] < 1.8);  tests++;
      failure += (yrate[2] < 2.7);  tests++;
      failure += (yrate[3] < 2.7);  tests++;
      failure += (yrate[4] < 3.3);  tests++;
      failure += (yrate[5] < 3.3);  tests++;
      failure += (dyrate[1] < 0.9);  tests++;
      failure += (dyrate[2] < 1.7);  tests++;
      failure += (dyrate[3] < 2.1);  tests++;
      failure += (dyrate[4] < 2.6);  tests++;
      failure += (dyrate[5] < 2.6);  tests++;
      failure += (d2yrate[2] < 0.7);  tests++;
      failure += (d2yrate[3] < 1.3);  tests++;
      failure += (d2yrate[4] < 1.6);  tests++;
      failure += (d2yrate[5] < 1.6);  tests++;
      /* Lagrange module convergence */
      failure += (yrate2[0] < 0.7);  tests++;
      failure += (yrate2[1] < 1.8);  tests++;
      failure += (yrate2[2] < 2.5);  tests++;
      failure += (yrate2[3] < 2.8);  tests++;
      failure += (yrate2[4] < 2.8);  tests++;
      failure += (yrate2[5] < 2.8);  tests++;
      failure += (dyrate2[1] < 0.9);  tests++;
      failure += (dyrate2[2] < 1.5);  tests++;
      failure += (dyrate2[3] < 2.5);  tests++;
      failure += (dyrate2[4] < 2.8);  tests++;
      failure += (dyrate2[5] < 2.8);  tests++;
      failure += (d2yrate2[2] < 0.6);  tests++;
      failure += (d2yrate2[3] < 1.6);  tests++;
      failure += (d2yrate2[4] < 2.0);  tests++;
      failure += (d2yrate2[5] < 2.0);  tests++;
    } else if (SUNRabs(lambda) < 1.1e4) { /* lambda = -1e4 */
      /* Hermite module convergence */
      failure += (yrate[0] < 0.8);  tests++;
      failure += (yrate[1] < 1.8);  tests++;
      failure += (yrate[2] < 2.7);  tests++;
      failure += (yrate[3] < 2.5);  tests++;
      failure += (yrate[4] < 3.5);  tests++;
      failure += (yrate[5] < 3.5);  tests++;
      failure += (dyrate[1] < 0.9);  tests++;
      failure += (dyrate[2] < 1.7);  tests++;
      failure += (dyrate[3] < 1.6);  tests++;
      failure += (dyrate[4] < 2.6);  tests++;
      failure += (dyrate[5] < 2.6);  tests++;
      failure += (d2yrate[2] < 0.7);  tests++;
      failure += (d2yrate[3] < 0.7);  tests++;
      /* Lagrange module convergence */
      failure += (yrate2[0] < 0.7);  tests++;
      failure += (yrate2[1] < 1.8);  tests++;
      failure += (yrate2[2] < 2.4);  tests++;
      failure += (yrate2[3] < 2.7);  tests++;
      failure += (yrate2[4] < 2.7);  tests++;
      failure += (yrate2[5] < 2.7);  tests++;
      failure += (dyrate2[1] < 0.9);  tests++;
      failure += (dyrate2[2] < 1.5);  tests++;
      failure += (dyrate2[3] < 2.5);  tests++;
      failure += (dyrate2[4] < 2.6);  tests++;
      failure += (dyrate2[5] < 2.6);  tests++;
      failure += (d2yrate2[2] < 0.6);  tests++;
      failure += (d2yrate2[3] < 1.6);  tests++;
      failure += (d2yrate2[4] < 2.0);  tests++;
      failure += (d2yrate2[5] < 2.0);  tests++;
    } else {                              /* lambda = -1e6 */
      /* Hermite module convergence */
      failure += (yrate[0] < 0.8);  tests++;
      failure += (yrate[1] < 1.8);  tests++;
      failure += (yrate[2] < 2.5);  tests++;
      failure += (yrate[3] < 2.2);  tests++;
      failure += (yrate[4] < 3.4);  tests++;
      failure += (yrate[5] < 3.4);  tests++;
      failure += (dyrate[1] < 0.9);  tests++;
      failure += (dyrate[2] < 1.5);  tests++;
      failure += (dyrate[3] < 1.2);  tests++;
      failure += (dyrate[4] < 2.4);  tests++;
      failure += (dyrate[5] < 2.4);  tests++;
      failure += (d2yrate[2] < 0.6);  tests++;
      failure += (d2yrate[3] < 0.4);  tests++;
      /* Lagrange module convergence */
      failure += (yrate2[0] < 0.7);  tests++;
      failure += (yrate2[1] < 1.8);  tests++;
      failure += (yrate2[2] < 2.4);  tests++;
      failure += (yrate2[3] < 2.4);  tests++;
      failure += (yrate2[4] < 2.5);  tests++;
      failure += (yrate2[5] < 2.5);  tests++;
      failure += (dyrate2[1] < 0.9);  tests++;
      failure += (dyrate2[2] < 1.5);  tests++;
      failure += (dyrate2[3] < 2.3);  tests++;
      failure += (dyrate2[4] < 2.4);  tests++;
      failure += (dyrate2[5] < 2.4);  tests++;
      failure += (d2yrate2[2] < 0.6);  tests++;
      failure += (d2yrate2[3] < 1.6);  tests++;
      failure += (d2yrate2[4] < 2.0);  tests++;
      failure += (d2yrate2[5] < 2.0);  tests++;
    }
  }

  /* output result */
  printf("\nTotal failures: %i (out of %i)\n",failure,tests);
  return(0);
}


/*---- end of file ----*/
