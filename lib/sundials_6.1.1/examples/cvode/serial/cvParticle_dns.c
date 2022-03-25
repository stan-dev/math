/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * Based on an example from Jean-Luc Fattebert @ ORNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example solves the equation for a particle moving conterclockwise with
 * velocity alpha on the unit circle in the xy-plane. The ODE system is given by
 *
 *   x' = -alpha * y
 *   y' =  alpha * x
 *
 * where x and y are subject to the constraint
 *
 *   x^2 + y^2 - 1 = 0
 *
 * with initial condition x = 1 and y = 0 at t = 0. The system has the analytic
 * solution
 *
 *  x(t) = cos(alpha * t)
 *  y(t) = sin(alpha * t)
 *
 * For a description of the command line options for this example run the
 * program with the --help flag.
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <cvode/cvode.h>               /* access to CVODE                 */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */

/* Precision specific formatting macros */
#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* Precision specific math function macros */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SIN(x)   (sin((x)))
#define COS(x)   (cos((x)))
#define SQRT(x)  (sqrt((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SIN(x)   (sinf((x)))
#define COS(x)   (cosf((x)))
#define SQRT(x)  (sqrtf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SIN(x)   (sinl((x)))
#define COS(x)   (cosl((x)))
#define SQRT(x)  (sqrtl((x)))
#endif

/* Problem Constants */
#define PI    RCONST(3.141592653589793238462643383279502884197169)
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

/* User-defined data structure */
typedef struct UserData_
{
  realtype alpha; /* particle velocity */

  int      orbits; /* number of orbits */
  realtype torbit; /* orbit time       */

  realtype rtol; /* integration tolerances */
  realtype atol;

  int proj;    /* enable/disable solution projection */
  int projerr; /* enable/disable error projection */

  int tstop; /* use tstop mode */
  int nout;  /* number of outputs per orbit */

} *UserData;

/* Functions provided to CVODE */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Proj(realtype t, N_Vector ycur, N_Vector corr, realtype epsProj,
                N_Vector err, void *user_data);

/* Utility functions */
static int InitUserData(int *argc, char ***argv, UserData udata);
static int PrintUserData(UserData udata);
static void InputHelp();
static int ComputeSolution(realtype t, N_Vector y, UserData udata);
static int ComputeError(realtype t, N_Vector y, N_Vector e, realtype *ec,
                        UserData udata);
static int WriteOutput(realtype t, N_Vector y, N_Vector e, realtype ec,
                       int screenfile, FILE *YFID, FILE *EFID);
static int PrintStats(void *cvode_mem);
static int check_retval(void *returnvalue, const char *funcname, int opt);


/* -----------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
  int      retval;          /* reusable return flag       */
  int      out      = 0;    /* output counter             */
  int      totalout = 0;    /* output counter             */
  realtype t        = ZERO; /* current integration time   */
  realtype dtout    = ZERO; /* output spacing             */
  realtype tout     = ZERO; /* next output time           */
  realtype ec       = ZERO; /* constraint error           */
  UserData udata    = NULL; /* user data structure        */

  SUNContext      sunctx     = NULL; /* SUNDIALS context     */
  void            *cvode_mem = NULL; /* CVODE memory         */
  N_Vector         y         = NULL; /* solution vector      */
  realtype        *ydata     = NULL; /* solution vector data */
  N_Vector         e         = NULL; /* error vector         */
  SUNMatrix        A         = NULL; /* Jacobian matrix      */
  SUNLinearSolver  LS        = NULL; /* linear solver        */

  FILE *YFID = NULL; /* solution output file */
  FILE *EFID = NULL; /* error output file    */

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if(check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* Allocate and initialize user data structure */
  udata = (UserData) malloc(sizeof *udata);
  if (check_retval((void *)udata, "malloc", 0)) return(1);

  retval = InitUserData(&argc, &argv, udata);
  if (check_retval(&retval, "InitUserData", 1)) return(1);

  /* Create serial vector to store the solution */
  y = N_VNew_Serial(2, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);

  /* Set initial contion */
  ydata    = N_VGetArrayPointer(y);
  ydata[0] = ONE;
  ydata[1] = ZERO;

  /* Create serial vector to store the solution error */
  e = N_VClone(y);
  if (check_retval((void *)y, "N_VClone", 0)) return(1);

  /* Set initial error */
  N_VConst(ZERO, e);

  /* Create CVODE memory */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Initialize CVODE */
  retval = CVodeInit(cvode_mem, f, t, y);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Attach user-defined data structure to CVODE */
  retval = CVodeSetUserData(cvode_mem, udata);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Set integration tolerances */
  retval = CVodeSStolerances(cvode_mem, udata->rtol, udata->atol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(2, 2, sunctx);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(y, A, sunctx);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver to CVODE */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* Set a user-supplied Jacobian function */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

  /* Set a user-supplied projection function */
  if (udata->proj)
  {
    retval = CVodeSetProjFn(cvode_mem, Proj);
    if(check_retval(&retval, "CVodeSetProjFn", 1)) return(1);

    retval = CVodeSetProjErrEst(cvode_mem, udata->projerr);
    if(check_retval(&retval, "CVodeSetProjErrEst", 1)) return(1);
  }

  /* Set max steps between outputs */
  retval = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return(1);

  /* Output problem setup */
  retval = PrintUserData(udata);
  if(check_retval(&retval, "PrintUserData", 1)) return(1);

  /* Output initial condition */
  printf("\n     t            x              y");
  printf("             err x          err y       err constr\n");
  WriteOutput(t, y, e, ec, 0, NULL, NULL);

  if (udata->nout > 0)
  {
    YFID = fopen("cvParticle_solution.txt","w");
    EFID = fopen("cvParticle_error.txt","w");
    WriteOutput(t, y, e, ec, 1, YFID, EFID);
  }

  /* Integrate in time and periodically output the solution and error */
  if (udata->nout > 0)
  {
    totalout = udata->orbits * udata->nout;
    dtout    = udata->torbit / udata->nout;
  }
  else
  {
    totalout = 1;
    dtout    = udata->torbit * udata->orbits;
  }
  tout = dtout;

  for (out = 0; out < totalout; out++)
  {
    /* Stop at output time (do not interpolate output) */
    if (udata->tstop || udata->nout == 0)
    {
      retval = CVodeSetStopTime(cvode_mem, tout);
      if (check_retval(&retval, "CVodeSetStopTime", 1)) return(1);
    }

    /* Advance in time */
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1)) break;

    /* Output solution and error */
    if (udata->nout > 0)
    {
      retval = ComputeError(t, y, e, &ec, udata);
      if (check_retval(&retval, "ComputeError", 1)) break;

      WriteOutput(t, y, e, ec, 1, YFID, EFID);
      if (check_retval(&retval, "WriteOutput", 1)) break;
    }

    /* Update output time */
    if (out < totalout - 1)
    {
      tout += dtout;
    }
    else
    {
      tout = udata->torbit * udata->orbits;
    }
  }

  /* Close output files */
  if (udata->nout > 0)
  {
    fclose(YFID);
    fclose(EFID);
  }

  /* Output final solution and error to screen */
  ComputeError(t, y, e, &ec, udata);
  if (check_retval(&retval, "ComputeError", 1)) return(1);

  WriteOutput(t, y, e, ec, 0, NULL, NULL);
  if (check_retval(&retval, "WriteOutput", 1)) return(1);

  /* Print some final statistics */
  PrintStats(cvode_mem);

  /* Free memory */
  N_VDestroy(y);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Functions provided to CVODE
 * ---------------------------------------------------------------------------*/


/* Compute the right-hand side function, y' = f(t,y) */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData  udata = (UserData) user_data;
  realtype *ydata = N_VGetArrayPointer(y);
  realtype *fdata = N_VGetArrayPointer(ydot);

  fdata[0] = -(udata->alpha) * ydata[1];
  fdata[1] =  (udata->alpha) * ydata[0];

  return(0);
}


/* Compute the Jacobian of the right-hand side function, J(t,y) = df/dy */
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData  udata = (UserData) user_data;
  realtype *Jdata = SUNDenseMatrix_Data(J);

  Jdata[0] =  ZERO;
  Jdata[1] = -(udata->alpha);
  Jdata[2] =  (udata->alpha);
  Jdata[3] =  ZERO;

  return(0);
}

/* Project the solution onto the constraint manifold */
static int Proj(realtype t, N_Vector ycur, N_Vector corr, realtype epsProj,
                N_Vector err, void *user_data)
{
  realtype *ydata = N_VGetArrayPointer(ycur);
  realtype *cdata = N_VGetArrayPointer(corr);
  realtype *edata = NULL;
  realtype  x = ydata[0];
  realtype  y = ydata[1];
  realtype  xp, yp, r;
  realtype  errxp, erryp;

  /* project onto the unit circle */
  r = SQRT(x * x + y * y);

  xp = x / r;
  yp = y / r;

  /* correction to the unprojected solution */
  cdata[0] = xp - x;
  cdata[1] = yp - y;

  /* project the error */
  if (err != NULL)
  {
    edata = N_VGetArrayPointer(err);

    errxp =  edata[0] * yp * yp - edata[1] * xp * yp;
    erryp = -edata[0] * xp * yp + edata[1] * xp * xp;

    edata[0] = errxp;
    edata[1] = erryp;
  }

  return(0);
}


/* -----------------------------------------------------------------------------
 * Private helper functions
 * ---------------------------------------------------------------------------*/

static int InitUserData(int *argc, char ***argv, UserData udata)
{
  int arg_idx = 1;

  /* set default values */
  udata->alpha = ONE;

  udata->orbits = 100;
  udata->torbit = (TWO * PI) / udata->alpha;

  udata->rtol = RCONST(1.0e-4);
  udata->atol = RCONST(1.0e-9);

  udata->proj    = 1;
  udata->projerr = 0;

  udata->tstop = 0;
  udata->nout  = 0;

  /* check for input args */
  while (arg_idx < (*argc))
  {
    if (strcmp((*argv)[arg_idx],"--alpha") == 0)
    {
      arg_idx++;
      udata->alpha = atof((*argv)[arg_idx++]);
      udata->torbit  = (TWO * PI) / udata->alpha;
    }
    else if (strcmp((*argv)[arg_idx],"--orbits") == 0)
    {
      arg_idx++;
      udata->orbits = atoi((*argv)[arg_idx++]);
    }
    else if (strcmp((*argv)[arg_idx],"--rtol") == 0)
    {
      arg_idx++;
      udata->rtol = atof((*argv)[arg_idx++]);
    }
    else if (strcmp((*argv)[arg_idx],"--atol") == 0)
    {
      arg_idx++;
      udata->atol = atof((*argv)[arg_idx++]);
    }
    else if (strcmp((*argv)[arg_idx],"--proj") == 0)
    {
      arg_idx++;
      udata->proj = atoi((*argv)[arg_idx++]);
    }
    else if (strcmp((*argv)[arg_idx],"--projerr") == 0)
    {
      arg_idx++;
      udata->projerr = atoi((*argv)[arg_idx++]);
    }
    else if (strcmp((*argv)[arg_idx],"--nout") == 0)
    {
      arg_idx++;
      udata->nout = atoi((*argv)[arg_idx++]);
    }
    else if (strcmp((*argv)[arg_idx],"--tstop") == 0)
    {
      arg_idx++;
      udata->tstop = 1;
    }
    else if (strcmp((*argv)[arg_idx],"--help") == 0 )
    {
      InputHelp();
      return(-1);
    }
    else
    {
      fprintf(stderr, "ERROR: Invalid input %s",(*argv)[arg_idx]);
      InputHelp();
      return(-1);
    }
  }

  /* If projection is disabled then disable error projection */
  if (!(udata->proj)) udata->projerr = 0;

  return(0);
}

static int PrintUserData(UserData udata)
{
  if (udata == NULL) return(-1);

  printf("\nParticle traveling on the unit circle example\n");
  printf("---------------------------------------------\n");
  printf("alpha      = %0.4" ESYM"\n", udata->alpha);
  printf("num orbits = %d\n", udata->orbits);
  printf("---------------------------------------------\n");
  printf("rtol       = %" GSYM"\n", udata->rtol);
  printf("atol       = %" GSYM"\n", udata->atol);
  printf("proj sol   = %d\n", udata->proj);
  printf("proj err   = %d\n", udata->projerr);
  printf("nout       = %d\n", udata->nout);
  printf("tstop      = %d\n", udata->tstop);
  printf("---------------------------------------------\n");

  return(0);
}


/* Print command line options */
static void InputHelp()
{
  printf("\nCommand line options:\n");
  printf("  --alpha <vel>      : particle velocity\n");
  printf("  --orbits <orbits>  : number of orbits to perform\n");
  printf("  --rtol <rtol>      : relative tolerance\n");
  printf("  --atol <atol>      : absoltue tolerance\n");
  printf("  --proj <1 or 0>    : enable (1) / disable (0) projection\n");
  printf("  --projerr <1 or 0> : enable (1) / disable (0) error projection\n");
  printf("  --nout <nout>      : outputs per period\n");
  printf("  --tstop            : stop at output time (do not interpolate)\n");
  return;
}


/* Compute the analytical solution */
static int ComputeSolution(realtype t, N_Vector y, UserData udata)
{
  realtype *ydata = N_VGetArrayPointer(y);

  ydata[0] = COS((udata->alpha) * t);
  ydata[1] = SIN((udata->alpha) * t);

  return(0);
}


/* Compute the error in the solution and constraint */
static int ComputeError(realtype t, N_Vector y, N_Vector e, realtype *ec,
                        UserData udata)
{
  realtype *ydata = N_VGetArrayPointer(y);
  int retval;

  /* solution error */
  retval = ComputeSolution(t, e, udata);
  if (check_retval(&retval, "ComputeSolution", 1)) return(1);
  N_VLinearSum(ONE, y, -ONE, e, e);

  /* constraint error */
  *ec = ydata[0] * ydata[0] + ydata[1] * ydata[1] - ONE;

  return(0);
}

/* Output the solution to the screen or disk */
static int WriteOutput(realtype t, N_Vector y, N_Vector e, realtype ec,
                       int screenfile, FILE* YFID, FILE* EFID)
{
  realtype *ydata = N_VGetArrayPointer(y);
  realtype *edata = N_VGetArrayPointer(e);

  if (screenfile == 0)
  {
    /* output solution and error to screen */
    printf("%0.4" ESYM" %14.6" ESYM" %14.6" ESYM" %14.6" ESYM" %14.6" ESYM" %14.6" ESYM"\n",
           t, ydata[0], ydata[1], edata[0], edata[1], ec);
  }
  else
  {
    /* check file pointers */
    if (YFID == NULL || EFID == NULL) return(1);

    /* output solution to disk */
    fprintf(YFID, "%24.16" ESYM" %24.16" ESYM" %24.16"ESYM"\n",
            t, ydata[0], ydata[1]);

    /* output error to disk */
    fprintf(EFID,
            "%24.16" ESYM" %24.16" ESYM" %24.16"ESYM" %24.16"ESYM"\n",
            t, edata[0], edata[1], ec);
  }

  return(0);
}


/* Print final statistics */
static int PrintStats(void *cvode_mem)
{
  int retval;
  long int nst, nfe, nsetups, nje, nni, ncfn, netf;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);

  printf("\nIntegration Statistics:\n");

  printf("Number of steps taken = %-6ld\n", nst);
  printf("Number of function evaluations = %-6ld\n", nfe);

  printf("Number of linear solver setups = %-6ld\n", nsetups);
  printf("Number of Jacobian evaluations = %-6ld\n", nje);

  printf("Number of nonlinear solver iterations = %-6ld\n", nni);
  printf("Number of convergence failures = %-6ld\n", ncfn);
  printf("Number of error test failures = %-6ld\n", netf);

  return(0);
}

/* Check function return value */
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Opt 0: Check if a NULL pointer was returned - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nERROR: %s() returned a NULL pointer\n\n",
            funcname);
    return(1);
  }
  /* Opt 1: Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int *) returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nERROR: %s() returned = %d\n\n",
              funcname, *retval);
      return(1);
    }
  }

  return(0);
}
