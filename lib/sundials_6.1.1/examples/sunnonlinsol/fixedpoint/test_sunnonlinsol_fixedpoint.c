/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the testing routine to check the SUNNonlinearSolver fixed point
 * module. This test solves the nonlinear system
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * where the fixed point function is
 *
 * g1(x,y,z) = 1/3 cos((y-1)yz) + 1/6
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
 * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
 *
 * This system has the analytic solution x = 1/2, y = 1, z = -pi/6.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "nvector/nvector_serial.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

/* precision specific formatting macros */
#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* precision specific math function macros */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRsin(x) (sin((x)))
#define SUNRcos(x) (cos((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRsin(x) (sinf((x)))
#define SUNRcos(x) (cosf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRsin(x) (sinl((x)))
#define SUNRcos(x) (cosl((x)))
#endif

/* problem constants */
#define NEQ   3 /* number of equations */

#define ZERO         RCONST(0.0)             /* real 0.0  */
#define PTONE        RCONST(0.1)             /* real 0.1  */
#define HALF         RCONST(0.5)             /* real 0.5  */
#define PTNINE       RCONST(0.9)             /* real 0.9  */
#define ONE          RCONST(1.0)             /* real 1.0  */
#define ONEPTZEROSIX RCONST(1.06)            /* real 1.06 */
#define THREE        RCONST(3.0)             /* real 3.0  */
#define SIX          RCONST(6.0)             /* real 6.0  */
#define NINE         RCONST(9.0)             /* real 9.0  */
#define TEN          RCONST(10.0)            /* real 10.0 */
#define TWENTY       RCONST(20.0)            /* real 20.0 */
#define SIXTY        RCONST(60.0)            /* real 60.0 */
#define PI           RCONST(3.1415926535898) /* real pi   */

/* analytic solution */
#define XTRUE HALF
#define YTRUE ONE
#define ZTRUE -PI/SIX

/* Check the system solution */
static int check_ans(N_Vector ycur, realtype tol);

/* Check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

/* Nonlinear fixed point function */
static int FPFunction(N_Vector y, N_Vector f, void *mem);

/* Convergence test function */
static int ConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                    realtype tol, N_Vector ewt, void* mem);

/*
 * Proxy for integrator memory struct
 */

/* Integrator memory structure */
typedef struct IntegratorMemRec {
  N_Vector y0;
  N_Vector ycor;
  N_Vector ycur;
  N_Vector w;
} *IntegratorMem;

/* -----------------------------------------------------------------------------
 * Main testing routine
 * ---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  IntegratorMem      Imem    = NULL;
  int                retval  = 0;
  SUNNonlinearSolver NLS     = NULL;
  realtype           tol     = 100 * SUNRsqrt(UNIT_ROUNDOFF);
  int                mxiter  = 20;
  int                maa     = 0;           /* no acceleration */
  realtype           damping = RCONST(1.0); /* no damping      */
  long int           niters  = 0;
  realtype*          data    = NULL;
  SUNContext         sunctx     = NULL;

  /* Check if a acceleration/dampling values were provided */
  if (argc > 1) maa     = (long int) atoi(argv[1]);
  if (argc > 2) damping = (realtype) atof(argv[2]);

  /* Print problem description */
  printf("Solve the nonlinear system:\n");
  printf("    3x - cos((y-1)z) - 1/2 = 0\n");
  printf("    x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0\n");
  printf("    exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0\n");
  printf("Analytic solution:\n");
  printf("    x = %"GSYM"\n", XTRUE);
  printf("    y = %"GSYM"\n", YTRUE);
  printf("    z = %"GSYM"\n", ZTRUE);
  printf("Solution method: Anderson accelerated fixed point iteration.\n");
  printf("    tolerance = %"GSYM"\n", tol);
  printf("    max iters = %d\n", mxiter);
  printf("    accel vec = %d\n", maa);
  printf("    damping   = %"GSYM"\n", damping);

  /* create SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* create proxy for integrator memory */
  Imem = (IntegratorMem) malloc(sizeof(struct IntegratorMemRec));
  if (check_retval((void *)Imem, "Creating Integrator Memory", 0)) return(1);

  /* create vectors */
  Imem->y0 = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)Imem->y0, "N_VNew_Serial", 0)) return(1);

  Imem->ycor = N_VClone(Imem->y0);
  if (check_retval((void *)Imem->ycor, "N_VClone", 0)) return(1);

  Imem->ycur = N_VClone(Imem->y0);
  if (check_retval((void *)Imem->ycur, "N_VClone", 0)) return(1);

  Imem->w = N_VClone(Imem->y0);
  if (check_retval((void *)Imem->w, "N_VClone", 0)) return(1);

  /* set initial guess */
  data = N_VGetArrayPointer(Imem->y0);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return(1);

  data[0] =  PTONE;
  data[1] =  PTONE;
  data[2] = -PTONE;

  /* set inital correction */
  N_VConst(ZERO, Imem->ycor);

  /* set weights */
  N_VConst(ONE, Imem->w);

  /* create nonlinear solver */
  NLS = SUNNonlinSol_FixedPoint(Imem->y0, maa, sunctx);
  if (check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0)) return(1);

  /* set the nonlinear residual function */
  retval = SUNNonlinSolSetSysFn(NLS, FPFunction);
  if (check_retval(&retval, "SUNNonlinSolSetSysFn", 1)) return(1);

  /* set the convergence test function */
  retval = SUNNonlinSolSetConvTestFn(NLS, ConvTest, NULL);
  if (check_retval(&retval, "SUNNonlinSolSetConvTestFn", 1)) return(1);

  /* set the maximum number of nonlinear iterations */
  retval = SUNNonlinSolSetMaxIters(NLS, mxiter);
  if (check_retval(&retval, "SUNNonlinSolSetMaxIters", 1)) return(1);

  /* set the damping parameter */
  retval = SUNNonlinSolSetDamping_FixedPoint(NLS, damping);
  if (check_retval(&retval, "SUNNonlinSolSetDamping", 1)) return(1);

  /* solve the nonlinear system */
  retval = SUNNonlinSolSolve(NLS, Imem->y0, Imem->ycor, Imem->w, tol, SUNTRUE,
                             Imem);
  if (check_retval(&retval, "SUNNonlinSolSolve", 1)) return(1);

  /* update the initial guess with the final correction */
  N_VLinearSum(ONE, Imem->y0, ONE, Imem->ycor, Imem->ycur);

  /* get the number of linear iterations */
  retval = SUNNonlinSolGetNumIters(NLS, &niters);
  if (check_retval(&retval, "SUNNonlinSolGetNumIters", 1)) return(1);

  printf("Number of nonlinear iterations: %ld\n",niters);

  /* check solution */
  retval = check_ans(Imem->ycur, tol);

  /* Free vector, matrix, linear solver, and nonlinear solver */
  N_VDestroy(Imem->y0);
  N_VDestroy(Imem->ycor);
  N_VDestroy(Imem->ycur);
  N_VDestroy(Imem->w);
  SUNNonlinSolFree(NLS);
  free(Imem);
  SUNContext_Free(&sunctx);

  return(retval);
}

/* Proxy for integrator convergence test function */
int ConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del, realtype tol,
             N_Vector ewt, void* mem)
{
  realtype delnrm;

  /* compute the norm of the correction */
  delnrm = N_VMaxNorm(del);

  if (delnrm <= tol) return(SUN_NLS_SUCCESS);  /* success       */
  else               return(SUN_NLS_CONTINUE); /* not converged */
}


/* -----------------------------------------------------------------------------
 * Nonlinear system F(x,y,z):
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * Nonlinear fixed point function G(x,y,z):
 *
 * G1(x,y,z) = 1/3 cos((y-1)yz) + 1/6
 * G2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
 * G3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
 *
 * Corrector form g(x,y,z):
 *
 * g1(x,y,z) = 1/3 cos((y-1)yz) + 1/6 - x0
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9 - y0
 * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60 - z0
 *
 * ---------------------------------------------------------------------------*/
int FPFunction(N_Vector ycor, N_Vector gvec, void *mem)
{
  IntegratorMem Imem;
  realtype*     ydata = NULL;
  realtype*     gdata = NULL;
  realtype      x, y, z;

  if (mem == NULL) {
    printf("ERROR: Integrator memory is NULL");
    return(-1);
  }
  Imem = (IntegratorMem) mem;

  /* update state based on current correction */
  N_VLinearSum(ONE, Imem->y0, ONE, ycor, Imem->ycur);

  /* Get vector data arrays */
  ydata = N_VGetArrayPointer(Imem->ycur);
  if (check_retval((void*)ydata, "N_VGetArrayPointer", 0)) return(-1);

  gdata = N_VGetArrayPointer(gvec);
  if (check_retval((void*)gdata, "N_VGetArrayPointer", 0)) return(-1);

  /* get vector components */
  x = ydata[0];
  y = ydata[1];
  z = ydata[2];

  /* compute fixed point function */
  gdata[0] = (ONE/THREE) * SUNRcos((y-ONE)*z) + (ONE/SIX);
  gdata[1] = (ONE/NINE) * SUNRsqrt(x*x + SUNRsin(z) + ONEPTZEROSIX) + PTNINE;
  gdata[2] = -(ONE/TWENTY) * SUNRexp(-x*(y-ONE)) - (TEN * PI - THREE) / SIXTY;

  N_VLinearSum(ONE, gvec, -ONE, Imem->y0, gvec);

  return(0);
}

/* -----------------------------------------------------------------------------
 * Check the solution of the nonlinear system and return PASS or FAIL
 * ---------------------------------------------------------------------------*/
static int check_ans(N_Vector ycur, realtype tol)
{
  realtype* data = NULL;
  realtype  ex, ey, ez;

  /* Get vector data array */
  data = N_VGetArrayPointer(ycur);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return(1);

  /* print the solution */
  printf("Computed solution:\n");
  printf("    y1 = %"GSYM"\n", data[0]);
  printf("    y2 = %"GSYM"\n", data[1]);
  printf("    y3 = %"GSYM"\n", data[2]);

  /* solution error */
  ex = SUNRabs(data[0] - XTRUE);
  ey = SUNRabs(data[1] - YTRUE);
  ez = SUNRabs(data[2] - ZTRUE);

  /* print the solution error */
  printf("Solution error:\n");
  printf("    ex = %"GSYM"\n", ex);
  printf("    ey = %"GSYM"\n", ey);
  printf("    ez = %"GSYM"\n", ez);

  tol *= TEN;
  if (ex > tol || ey > tol || ez > tol) {
    printf("FAIL\n");
    return(1);
  }

  printf("PASS\n");
  return(0);
}

/* -----------------------------------------------------------------------------
 * Check function return value
 *   opt == 0 check if returned NULL pointer
 *   opt == 1 check if returned a non-zero value
 * ---------------------------------------------------------------------------*/
static int check_retval(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if the function returned a NULL pointer -- no memory allocated */
  if (opt == 0) {
    if (flagvalue == NULL) {
      fprintf(stderr, "\nERROR: %s() failed -- returned NULL\n\n", funcname);
      return(1);
    } else {
      return(0);
    }
  }

  /* Check if the function returned an non-zero value -- internal failure */
  if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag != 0) {
      fprintf(stderr, "\nERROR: %s() failed -- returned %d\n\n", funcname, *errflag);
      return(1);
    } else {
      return(0);
    }
  }

  /* if we make it here then opt was not 0 or 1 */
  fprintf(stderr, "\nERROR: check_retval failed -- Invalid opt value\n\n");
  return(1);
}
