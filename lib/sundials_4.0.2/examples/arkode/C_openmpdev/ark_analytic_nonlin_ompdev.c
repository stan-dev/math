/*-----------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ LLNL
 *---------------------------------------------------------------
 * This code is based on the serial code found in
 * ark_analytic_nonlin.c developed by Daniel R. Reynolds
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
#include <arkode/arkode_erkstep.h>     /* prototypes for ERKStep fcts., consts */
#include <nvector/nvector_openmpdev.h> /* OpenMPDEV N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>   /* def. of type 'realtype' */
#include <sundials/sundials_math.h>    /* def. of SUNRsqrt, etc. */

#ifdef _OPENMP
#include <omp.h>                       /* OpenMP functions */
#endif

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
  realtype *y_data = NULL;

  /* Initial problem output */
  printf("\nAnalytical ODE test problem:\n");
  printf("   reltol = %.1"ESYM"\n",  reltol);
  printf("   abstol = %.1"ESYM"\n\n",abstol);

  /* Initialize data structures */
  y = N_VNew_OpenMPDEV(NEQ);          /* Create OpenMPDEV vector for solution */
  if (check_flag((void *)y, "N_VNew_OpenMPDEV", 0)) return 1;
  y_data = N_VGetHostArrayPointer_OpenMPDEV(y);
  y_data[0] = 0.0;                                /* Specify initial condition */
  N_VCopyToDevice_OpenMPDEV(y);                   /* Copy to device */
  arkode_mem = ERKStepCreate(f, T0, y);     /* Create the solver memory */
  if (check_flag((void *)arkode_mem, "ERKStepCreate", 0)) return 1;

  /* Specify tolerances */
  flag = ERKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ERKStepSStolerances", 1)) return 1;

  /* Open output stream for results, output comment line */
  UFID = fopen("solution.txt","w");
  fprintf(UFID,"# t u\n");

  /* output initial condition to disk */
  N_VCopyFromDevice_OpenMPDEV(y);
  fprintf(UFID," %.16"ESYM" %.16"ESYM"\n", T0, y_data[0]);

  /* Main time-stepping loop: calls ERKStep to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  tout = T0+dTout;
  printf("        t           u\n");
  printf("   ---------------------\n");
  while (Tf - t > 1.0e-15) {

    flag = ERKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL); /* call integrator */
    if (check_flag(&flag, "ERKStep", 1)) break;
    N_VCopyFromDevice_OpenMPDEV(y);
    printf("  %10.6"FSYM"  %10.6"FSYM"\n", t, y_data[0]);      /* access/print solution */
    fprintf(UFID," %.16"ESYM" %.16"ESYM"\n", t, y_data[0]);
    if (flag >= 0) {                                           /* successful solve: update time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                   /* unsuccessful solve: break */
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
  ERKStepFree(&arkode_mem);     /* Free integrator memory */
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  int dev;
  realtype *y_data    = N_VGetDeviceArrayPointer_OpenMPDEV(y);
  realtype *ydot_data = N_VGetDeviceArrayPointer_OpenMPDEV(ydot);

  dev = omp_get_default_device();

#pragma omp target map(to:t) is_device_ptr(y_data, ydot_data) device(dev)
  {
    ydot_data[0] = (t+1.0)*SUNRexp(-1.0 * y_data[0]);
  }

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
