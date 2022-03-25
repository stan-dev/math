/* -----------------------------------------------------------------------------
 * Programmer(s): Radu Serban and David J. Gardner @ LLNL
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
 * This example solves a simple pendulum equation in Cartesian coordinates where
 * the pendulum bob has mass 1 and is suspended from the origin with a rod of
 * length 1. The governing equations are
 *
 * x'  = vx
 * y'  = vy
 * vx' = -x * T
 * vy' = -y * T - g
 *
 * with the constraints
 *
 * x^2 + y^2 - 1 = 0
 * x * vx + y * vy = 0
 *
 * where x and y are the pendulum bob position, vx and vy are the bob velocity
 * in the x and y directions respectively, T is the tension in the rod, and
 * g is acceleration due to gravity chosen such that the pendulum has period 2.
 * The initial condition at t = 0 is x = 1, y = 0, vx = 0, and vy = 0.
 *
 * A reference solution is computed using the pendulum equation in terms of the
 * angle between the x-axis and the pendulum rod i.e., theta in [0, -pi]. The
 * governing equations are
 *
 * theta'  = vtheta
 * vtheta' = -g * cos(theta)
 *
 * where theta is the angle from the x-axis, vtheta is the angular velocity, and
 * g the same acceleration due to gravity from above. The initial condition at
 * t = 0 is theta = 0 and vtheta = 0.
 *
 * The Cartesian formulation is run to a final time tf (default 30) with and
 * without projection for various integration tolerances. The error in the
 * position and velocity at tf compared to the reference solution, the error in
 * the position constraint equation, and various integrator statistics are
 * printed to the screen for each run.
 *
 * When projection is enabled a user-supplied function is used to project the
 * position, velocity, and error to the constraint manifold.
 *
 * Optional command line inputs may be used to change the final simulation time
 * (default 30), the initial tolerance (default 1e-5), the number of outputs
 * (default 1), or disable error projection. Use the option --help for a list
 * of the command line flags.
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <cvode/cvode.h>               /* access to CVODE                 */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNmatrix       */
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
#define ABS(x)   (fabs((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SIN(x)   (sinf((x)))
#define COS(x)   (cosf((x)))
#define SQRT(x)  (sqrtf((x)))
#define ABS(x)   (fabsf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SIN(x)   (sinl((x)))
#define COS(x)   (cosl((x)))
#define SQRT(x)  (sqrtl((x)))
#define ABS(x)   (fabsl((x)))
#endif

/* Problem Constants */
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define GRAV  RCONST(13.750371636040745654980191559621114395801712)

/* Functions provided to CVODE */
static int fref(realtype t, N_Vector yy, N_Vector fy, void *f_data);

static int f(realtype t, N_Vector yy, N_Vector fy, void *f_data);
static int proj(realtype t, N_Vector yy, N_Vector corr,
                realtype epsProj, N_Vector err, void *pdata);

/* Functions to integrate the Cartesian and reference solutions */
int GetSol(void *cvode_mem, N_Vector yy0, realtype rtol, realtype atol,
           realtype tf, int nout, booleantype proj, booleantype projerr,
           N_Vector yref);

int RefSol(realtype tf, N_Vector yref, int nout);

/* Utility functions */
static int ReadInputs(int *argc, char ***argv, realtype *rtol, realtype *atol,
                      realtype *tf, int *nout, booleantype *projerr);
static void InputHelp();
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* SUNDIALS context */
static SUNContext sunctx = NULL;

/* -----------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
  int         i;
  int         retval;                   /* reusable return flag    */
  int         nout    = 1;              /* number of outputs       */
  realtype    rtol    = RCONST(1.0e-5); /* base relative tolerance */
  realtype    atol    = RCONST(1.0e-5); /* base absolute tolerance */
  realtype    tf      = RCONST(30.0);   /* final integration time  */
  booleantype projerr = SUNTRUE;        /* enable error projection */

  void            *cvode_mem = NULL; /* CVODE memory              */
  N_Vector         yy0       = NULL; /* initial condition vector  */
  realtype        *yy0data   = NULL; /* vector data               */
  N_Vector         yref      = NULL; /* reference solution vector */
  SUNMatrix        A         = NULL; /* Jacobian matrix           */
  SUNLinearSolver  LS        = NULL; /* linear solver             */

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if(check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* Read command line inputs */
  retval = ReadInputs(&argc, &argv, &rtol, &atol, &tf, &nout, &projerr);
  if (check_retval(&retval, "ReadInputs", 1)) return(1);

  /* Compute reference solution */
  yref = N_VNew_Serial(4, sunctx);

  retval = RefSol(tf, yref, nout);
  if (check_retval(&retval, "RefSol", 1)) return(1);

  /* Create serial vector to store the initial condition */
  yy0 = N_VNew_Serial(4, sunctx);
  if (check_retval((void *)yy0, "N_VNew_Serial", 0)) return(1);

  /* Set the initial condition values */
  yy0data = N_VGetArrayPointer(yy0);

  yy0data[0] = ONE;  /* x  */
  yy0data[1] = ZERO; /* y  */
  yy0data[2] = ZERO; /* xd */
  yy0data[3] = ZERO; /* yd */

  /* Create CVODE memory */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Initialize CVODE */
  retval = CVodeInit(cvode_mem, f, ZERO, yy0);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(4, 4, sunctx);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy0, A, sunctx);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver to CVODE */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* Set a user-supplied projection function */
  retval = CVodeSetProjFn(cvode_mem, proj);
  if(check_retval(&retval, "CVodeSetProjFn", 1)) return(1);

  /* Set maximum number of steps between outputs */
  retval = CVodeSetMaxNumSteps(cvode_mem, 50000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return(1);

  /* Compute the solution with various tolerances */
  for (i = 0; i < 5; i++) {

    /* Output tolerance and output header for this run */
    printf("\n\nrtol = %8.2" ESYM", atol = %8.2" ESYM"\n", rtol, atol);
    printf("Project    x         y");
    printf("         x'        y'     |     g      |    ");
    printf("nst     rhs eval    setups (J eval)  |   cf   ef\n");

    /* Compute solution with projection */
    retval = GetSol(cvode_mem, yy0, rtol, atol, tf, nout, SUNTRUE, projerr, yref);
    if (check_retval(&retval, "GetSol", 1)) return(1);

    /* Compute solution without projection */
    retval = GetSol(cvode_mem, yy0, rtol, atol, tf, nout, SUNFALSE, SUNFALSE, yref);
    if (check_retval(&retval, "GetSol", 1)) return(1);

    /* Reduce tolerance for next run */
    rtol /= RCONST(10.0);
    atol /= RCONST(10.0);
  }

  /* Free memory */
  N_VDestroy_Serial(yref);
  N_VDestroy_Serial(yy0);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Functions to integrate the Cartesian and reference systems
 * ---------------------------------------------------------------------------*/


/* Compute the Cartesian system solution */
int GetSol(void *cvode_mem, N_Vector yy0, realtype rtol, realtype atol,
           realtype tf, int nout, booleantype proj, booleantype projerr,
           N_Vector yref)
{
  char      outname[100];  /* output file name */
  FILE     *FID    = NULL; /* output file      */
  N_Vector  yy     = NULL; /* solution vector  */
  realtype *yydata = NULL; /* vector data      */

  int      retval; /* reusable return flag */
  int      out;    /* output counter       */
  realtype dtout;  /* output frequency     */
  realtype tout;   /* output time          */
  realtype t;      /* return time          */
  realtype x, y;   /* position values      */
  realtype xd, yd; /* velocity values      */
  realtype g;      /* constraint value     */

  /* Integrator stats */
  long int nst, nfe, nsetups, nje, nfeLS, ncfn, netf;

  /* Enable or disable projection */
  if (proj)
  {
    printf("  YES   ");
    retval = CVodeSetProjFrequency(cvode_mem, 1);
    if(check_retval(&retval, "CVodeSetProjFrequency", 1)) return(1);

    /* Enable or disable error projection */
    retval = CVodeSetProjErrEst(cvode_mem, projerr);
    if(check_retval(&retval, "CVodeSetProjErrEst", 1)) return(1);
  }
  else
  {
    retval = CVodeSetProjFrequency(cvode_mem, 0);
    if(check_retval(&retval, "CVodeSetProjFrequency", 1)) return(1);
    printf("  NO    ");
  }

  /* Create vector to store the solution */
  yy = N_VNew_Serial(4, sunctx);

  /* Copy initial condition into solution vector */
  N_VScale(ONE, yy0, yy);

  /* Get pointer to vector data */
  yydata = N_VGetArrayPointer(yy);

  /* Reinitialize CVODE for this run */
  retval = CVodeReInit(cvode_mem, ZERO, yy0);
  if (check_retval(&retval, "CVodeReInit", 1))
  {
    N_VDestroy_Serial(yy);
    return(retval);
  }

  /* Set integration tolerances for this run */
  retval = CVodeSStolerances(cvode_mem, rtol, atol);
  if (check_retval(&retval, "CVodeSStolerances", 1))
  {
    N_VDestroy_Serial(yy);
    return(retval);
  }

  /* Open output file */
  if (proj)
  {
    sprintf(outname,
            "cvPendulum_dns_rtol_%03.2" ESYM"_atol_%03.2" ESYM"_proj.txt",
            rtol, atol);
  }
  else
  {
    sprintf(outname,
            "cvPendulum_dns_rtol_%03.2" ESYM"_atol_%03.2" ESYM".txt",
            rtol, atol);
  }
  FID = fopen(outname, "w");

  /* Output initial condition */
  fprintf(FID,
          "%24.16" ESYM" %24.16" ESYM" %24.16" ESYM" %24.16" ESYM" %24.16" ESYM"\n",
          ZERO, yydata[0], yydata[1], yydata[2], yydata[3]);

  /* Integrate to tf and peridoically output the solution */
  dtout = tf / nout;
  tout  = dtout;

  for (out = 0; out < nout; out++)
  {
    /* Set stop time (do not interpolate output) */
    retval = CVodeSetStopTime(cvode_mem, tout);
    if (check_retval(&retval, "CVodeSetStopTime", 1))
    {
      N_VDestroy_Serial(yy);
      fclose(FID);
      return(retval);
    }

    /* Integrate to tout */
    retval = CVode(cvode_mem, tout, yy, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1))
    {
      N_VDestroy_Serial(yy);
      fclose(FID);
      return(retval);
    }

    /* Write output */
    fprintf(FID,
            "%24.16" ESYM" %24.16" ESYM" %24.16" ESYM" %24.16" ESYM" %24.16" ESYM"\n",
            t, yydata[0], yydata[1], yydata[2], yydata[3]);

    /* Update output time */
    if (out < nout - 1)
    {
      tout += dtout;
    }
    else
    {
      tout = tf;
    }
  }

  /* Close output file */
  fclose(FID);

  /* Compute the constraint violation */
  x = yydata[0];
  y = yydata[1];
  g = ABS(x*x + y*y - ONE);

  /* Compute the absolute error compared to the reference solution */
  N_VLinearSum(ONE, yy, -ONE, yref, yy);
  N_VAbs(yy, yy);

  x  = yydata[0];
  y  = yydata[1];
  xd = yydata[2];
  yd = yydata[3];

  /* Output errors */
  printf("%8.2" ESYM"  %8.2" ESYM"  %8.2" ESYM"  %8.2" ESYM"  |  %8.2" ESYM"  |",
         x, y, xd, yd, g);

  /* Free solution vector */
  N_VDestroy_Serial(yy);

  /* Get integrator stats */
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  if (check_retval(&retval, "CVodeGetNumSteps", 1)) return(retval);

  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  if (check_retval(&retval, "CVodeGetNumFctEvals", 1)) return(retval);

  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  if (check_retval(&retval, "CVodeGetNumLinSolvSetups", 1)) return(retval);

  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  if (check_retval(&retval, "CVodeGetNumErrTestFails", 1)) return(retval);

  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  if (check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1)) return(retval);

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  if (check_retval(&retval, "CVodeGetNumJacEvals", 1)) return(retval);

  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  if (check_retval(&retval, "CVodeGetNumLinRhsEvals", 1)) return(retval);

  /* Output stats */
  printf(" %6ld   %6ld+%-4ld     %4ld (%3ld)     |  %3ld  %3ld\n",
         nst, nfe, nfeLS, nsetups, nje, ncfn, netf);

  return(0);
}


/* Compute the reference system solution */
int RefSol(realtype tf, N_Vector yref, int nout)
{
  FILE            *FID       = NULL; /* output file     */
  void            *cvode_mem = NULL; /* CVODE memory    */
  N_Vector         yy        = NULL; /* solution vector */
  realtype        *yydata    = NULL; /* vector data     */
  SUNMatrix        A         = NULL; /* Jacobian matrix */
  SUNLinearSolver  LS        = NULL; /* linear solver   */

  int      retval;                /* reusable return flag  */
  int      out;                   /* output counter        */
  realtype dtout;                 /* output frequency      */
  realtype tout;                  /* output time           */
  realtype t;                     /* return time           */
  realtype th, thd;               /* theta and theta dot   */
  realtype tol = RCONST(1.0e-14); /* integration tolerance */

  /* Create the solution vector */
  yy = N_VNew_Serial(2, sunctx);
  if (check_retval((void *)yy, "N_VNew_Serial", 0)) return(-1);

  /* Set the initial condition */
  yydata = N_VGetArrayPointer(yy);

  yydata[0] = ZERO; /* theta  */
  yydata[1] = ZERO; /* theta' */

  /* Create CVODE memory */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Initialize CVODE */
  retval = CVodeInit(cvode_mem, fref, ZERO, yy);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Set integration tolerances */
  retval = CVodeSStolerances(cvode_mem, tol, tol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(2, 2, sunctx);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A, sunctx);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver to CVODE */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* Set CVODE optional inputs */
  retval = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return(1);

  retval = CVodeSetStopTime(cvode_mem, tf);
  if (check_retval(&retval, "CVodeSetStopTime", 1)) return(1);

  /* Open output file */
  FID = fopen("cvPendulum_dns_ref.txt", "w");

  /* Output initial condition */
  th  = yydata[0];
  thd = yydata[1];
  fprintf(FID,
          "%24.16" ESYM" %24.16" ESYM" %24.16" ESYM" %24.16" ESYM" %24.16" ESYM"\n",
          ZERO, COS(th), SIN(th), -thd * SIN(th), thd * COS(th));

  /* Integrate to tf and periodically output the solution */
  dtout = tf / nout;
  tout  = dtout;

  for (out = 0; out < nout; out++)
  {
    /* Set stop time (do not interpolate output) */
    retval = CVodeSetStopTime(cvode_mem, tout);
    if (check_retval(&retval, "CVodeSetStopTime", 1))
    {
      N_VDestroy_Serial(yy);
      SUNMatDestroy(A);
      SUNLinSolFree(LS);
      CVodeFree(&cvode_mem);
      fclose(FID);
      return(retval);
    }

    /* Integrate to tout */
    retval = CVode(cvode_mem, tf, yy, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1))
    {
      N_VDestroy_Serial(yy);
      SUNMatDestroy(A);
      SUNLinSolFree(LS);
      CVodeFree(&cvode_mem);
      fclose(FID);
      return(retval);
    }

    /* Write output */
    th  = yydata[0];
    thd = yydata[1];
    fprintf(FID,
            "%24.16" ESYM" %24.16" ESYM" %24.16" ESYM" %24.16" ESYM" %24.16" ESYM"\n",
            t, COS(th), SIN(th), -thd * SIN(th), thd * COS(th));

    /* Update output time */
    if (out < nout - 1)
    {
      tout += dtout;
    }
    else
    {
      tout = tf;
    }
  }

  /* Close output file */
  fclose(FID);

  /* Get solution components */
  th  = yydata[0];
  thd = yydata[1];

  /* Convert to Cartesian reference solution */
  yydata = N_VGetArrayPointer(yref);

  yydata[0] = COS(th);
  yydata[1] = SIN(th);
  yydata[2] = -thd * SIN(th);
  yydata[3] =  thd * COS(th);

  /* Free memory */
  N_VDestroy_Serial(yy);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  CVodeFree(&cvode_mem);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Functions provided to CVODE
 * ---------------------------------------------------------------------------*/


/* ODE RHS function for the reference system */
static int fref(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  realtype *yydata = NULL; /* yy vector data */
  realtype *fydata = NULL; /* fy vector data */

  /* Get vector array pointers */
  yydata = N_VGetArrayPointer(yy);
  fydata = N_VGetArrayPointer(fy);

  fydata[0] = yydata[1];              /* theta'          */
  fydata[1] = -GRAV * COS(yydata[0]); /* -g * cos(theta) */
  return(0);
}


/* ODE RHS function for the Cartesian system */
static int f(realtype t, N_Vector yy, N_Vector fy, void *f_data)
{
  realtype *yydata = NULL; /* yy vector data */
  realtype *fydata = NULL; /* fy vector data */

  realtype x, y;   /* positions  */
  realtype xd, yd; /* velocities */
  realtype tmp;

  /* Get vector array pointers */
  yydata = N_VGetArrayPointer(yy);
  fydata = N_VGetArrayPointer(fy);

  /* Get vector components */
  x  = yydata[0];
  y  = yydata[1];
  xd = yydata[2];
  yd = yydata[3];

  /* Compute tension */
  tmp = xd * xd + yd * yd - GRAV * y;

  /* Compute RHS */
  fydata[0] = xd;
  fydata[1] = yd;
  fydata[2] = -x * tmp;
  fydata[3] = -y * tmp - GRAV;

  return(0);
}


/* Projection function */
static int proj(realtype t, N_Vector yy, N_Vector corr,
                realtype epsProj, N_Vector err, void *pdata)
{
  realtype *yydata = NULL; /* yy vector data   */
  realtype *cdata  = NULL; /* corr vector data */
  realtype *edata  = NULL; /* err vector data */

  realtype x, y, x_new, y_new;     /* positions  */
  realtype xd, yd, xd_new, yd_new; /* velocities */

  realtype e1, e2, e3, e4;
  realtype e1_new, e2_new, e3_new, e4_new;
  realtype R;

  /* Get vector array pointers */

  yydata = N_VGetArrayPointer(yy);
  cdata  = N_VGetArrayPointer(corr);

  /* Extract current solution */

  x  = yydata[0];
  y  = yydata[1];
  xd = yydata[2];
  yd = yydata[3];

  /* Project positions */

  R = SQRT(x * x + y * y);

  x_new = x / R;
  y_new = y / R;

  /* Project velocities
   *
   *        +-            -+  +-    -+
   *        |  y*y    -x*y |  |  xd  |
   *  P v = |              |  |      |
   *        | -x*y     x*x |  |  yd  |
   *        +-            -+  +-    -+
   */

  xd_new =   xd * y_new * y_new - yd * x_new * y_new;
  yd_new = - xd * x_new * y_new + yd * x_new * x_new;

  /* Return position and velocity corrections */

  cdata[0] = x_new  - x;
  cdata[1] = y_new  - y;
  cdata[2] = xd_new - xd;
  cdata[3] = yd_new - yd;

  /* Project error P * err */
  if (err != NULL)
  {
    edata = N_VGetArrayPointer(err);

    e1 = edata[0];
    e2 = edata[1];
    e3 = edata[2];
    e4 = edata[3];

    e1_new =  y_new * y_new * e1 - x_new * y_new * e2;
    e2_new = -x_new * y_new * e1 + x_new * x_new * e2;

    e3_new =  y_new * y_new * e3 - x_new * y_new * e4;
    e4_new = -x_new * y_new * e3 + x_new * x_new * e4;

    edata[0] = e1_new;
    edata[1] = e2_new;
    edata[2] = e3_new;
    edata[3] = e4_new;
  }

  return(0);
}


/* -----------------------------------------------------------------------------
 * Private helper functions
 * ---------------------------------------------------------------------------*/


/* Read command line unputs */
static int ReadInputs(int *argc, char ***argv, realtype *rtol, realtype *atol,
                      realtype *tf, int *nout, booleantype *projerr)
{
  int arg_idx = 1;

  /* check for input args */
  while (arg_idx < (*argc))
  {
    if (strcmp((*argv)[arg_idx],"--tol") == 0)
    {
      arg_idx++;
      *rtol = atof((*argv)[arg_idx++]);
      *atol = atof((*argv)[arg_idx++]);
    }
    else if (strcmp((*argv)[arg_idx],"--tf") == 0)
    {
      arg_idx++;
      *tf = atof((*argv)[arg_idx++]);
    }
    else if (strcmp((*argv)[arg_idx],"--nout") == 0)
    {
      arg_idx++;
      *nout = atoi((*argv)[arg_idx++]);
    }
    else if (strcmp((*argv)[arg_idx],"--noerrproj") == 0)
    {
      arg_idx++;
      *projerr = SUNFALSE;
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

  return(0);
}


/* Print command line options */
static void InputHelp()
{
  printf("\nCommand line options:\n");
  printf("  --tol <rtol> <atol> : relative and absolute tolerance\n");
  printf("  --tf <time>         : final simulation time\n");
  printf("  --nout <outputs>    : number of outputs\n");
  printf("  --noerrproj         : disable error projection\n");

  return;
}


/* Check function return value */
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Opt 0: Check if function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nERROR: %s() returned NULL pointer\n\n",
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
