/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem with analytical
 * solution,
 *     dy/dt = (t+1)*exp(-y)
 * for t in the interval [0.0, 10.0], with initial condition: y=0.
 * This has analytical solution
 *      y(t) = log(0.5*t^2 + t + 1)
 *
 * This program solves the problem with the ERK method.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 *-----------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode_erkstep.h>    /* prototypes for ERKStep fcts., consts */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>  /* def. of type 'realtype' */
#include <sundials/sundials_math.h>   /* def. of SUNRsqrt, etc. */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);     /* initial time */
  realtype Tf = RCONST(10.0);    /* final time */
  realtype dTout = RCONST(1.0);  /* time between outputs */
  sunindextype NEQ = 1;          /* number of dependent vars. */
  realtype reltol = 1.0e-6;      /* tolerances */
  realtype abstol = 1.0e-10;

  /* general problem variables */
  int flag;                      /* reusable error-checking flag */
  N_Vector y = NULL;             /* empty vector for storing solution */
  void *arkode_mem = NULL;       /* empty ARKode memory structure */
  FILE *UFID;
  realtype t, tout;
  long int nst, nst_a, nfe, netf;

  /* Initial problem output */
  printf("\nAnalytical ODE test problem:\n");
  printf("   reltol = %.1"ESYM"\n",  reltol);
  printf("   abstol = %.1"ESYM"\n\n",abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ);          /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  NV_Ith_S(y,0) = 0.0;             /* Specify initial condition */

  /* Call ERKStepCreate to initialize the ERK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y. */
  arkode_mem = ERKStepCreate(f, T0, y);
  if (check_flag((void *)arkode_mem, "ERKStepCreate", 0)) return 1;

  /* Specify tolerances */
  flag = ERKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ERKStepSStolerances", 1)) return 1;

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt","w");
  fprintf(UFID,"# t u\n");

  /* output initial condition to disk */
  fprintf(UFID," %.16"ESYM" %.16"ESYM"\n", T0, NV_Ith_S(y,0));

  /* Main time-stepping loop: calls ERKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  tout = T0+dTout;
  printf("        t           u\n");
  printf("   ---------------------\n");
  while (Tf - t > 1.0e-15) {

    flag = ERKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);       /* call integrator */
    if (check_flag(&flag, "ERKStepEvolve", 1)) break;
    printf("  %10.6"FSYM"  %10.6"FSYM"\n", t, NV_Ith_S(y,0));           /* access/print solution */
    fprintf(UFID," %.16"ESYM" %.16"ESYM"\n", t, NV_Ith_S(y,0));
    if (flag >= 0) {                                          /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                  /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }
  }
  printf("   ---------------------\n");
  fclose(UFID);

  /* Print some final statistics */
  flag = ERKStepGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ERKStepGetNumSteps", 1);
  flag = ERKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ERKStepGetNumStepAttempts", 1);
  flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
  check_flag(&flag, "ERKStepGetNumRhsEvals", 1);
  flag = ERKStepGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ERKStepGetNumErrTestFails", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals = %li\n", nfe);
  printf("   Total number of error test failures = %li\n\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy(y);               /* Free y vector */
  ERKStepFree(&arkode_mem);    /* Free integrator memory */
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  NV_Ith_S(ydot,0) = (t+1.0)*SUNRexp(-NV_Ith_S(y,0));
  return 0;
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


/*---- end of file ----*/
