/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This example solves the nonlinear system
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * using the accelerated fixed pointer solver in KINSOL. The nonlinear fixed
 * point function is
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

#include "kinsol/kinsol.h"           /* access to KINSOL func., consts. */
#include "nvector/nvector_serial.h"  /* access to serial N_Vector       */

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
#define ABS(x)  (fabs((x)))
#define SQRT(x) (sqrt((x)))
#define EXP(x)  (exp((x)))
#define SIN(x)  (sin((x)))
#define COS(x)  (cos((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define ABS(x)  (fabsf((x)))
#define SQRT(x) (sqrtf((x)))
#define EXP(x)  (expf((x)))
#define SIN(x)  (sinf((x)))
#define COS(x)  (cosf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define ABS(x)  (fabsl((x)))
#define SQRT(x) (sqrtl((x)))
#define EXP(x)  (expl((x)))
#define SIN(x)  (sinl((x)))
#define COS(x)  (cosl((x)))
#endif

/* problem constants */
#define NEQ 3 /* number of equations */

#define ZERO         RCONST(0.0)             /* real 0.0  */
#define PTONE        RCONST(0.1)             /* real 0.1  */
#define HALF         RCONST(0.5)             /* real 0.5  */
#define PTNINE       RCONST(0.9)             /* real 0.9  */
#define ONE          RCONST(1.0)             /* real 1.0  */
#define ONEPTZEROSIX RCONST(1.06)            /* real 1.06 */
#define ONEPTONE     RCONST(1.1)             /* real 1.1  */
#define THREE        RCONST(3.0)             /* real 3.0  */
#define FOUR         RCONST(4.0)             /* real 4.0  */
#define SIX          RCONST(6.0)             /* real 6.0  */
#define NINE         RCONST(9.0)             /* real 9.0  */
#define TEN          RCONST(10.0)            /* real 10.0 */
#define TWENTY       RCONST(20.0)            /* real 20.0 */
#define SIXTY        RCONST(60.0)            /* real 60.0 */
#define EIGHTYONE    RCONST(81.0)            /* real 81.0 */
#define PI           RCONST(3.1415926535898) /* real pi   */

/* analytic solution */
#define XTRUE HALF
#define YTRUE ONE
#define ZTRUE -PI/SIX

/* Nonlinear fixed point function */
static int FPFunction(N_Vector u, N_Vector f, void *user_data);

/* Check function return values */
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Check the system solution */
static int check_ans(N_Vector u, realtype tol);

/* -----------------------------------------------------------------------------
 * Main program
 * ---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int       retval  = 0;
  N_Vector  u       = NULL;
  N_Vector  scale   = NULL;
  realtype  tol     = 100 * SQRT(UNIT_ROUNDOFF);
  long int  mxiter  = 10;
  long int  maa     = 0;           /* no acceleration */
  realtype  damping = RCONST(1.0); /* no damping      */
  long int  nni, nfe;
  realtype* data;
  void*     kmem;

  /* Check if a acceleration/dampling values were provided */
  if (argc > 1) maa     = (long int) atoi(argv[1]);
  if (argc > 2) damping = (realtype) atof(argv[2]);

  /* -------------------------
   * Print problem description
   * ------------------------- */

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
  printf("    max iters = %ld\n", mxiter);
  printf("    accel vec = %ld\n", maa);
  printf("    damping   = %"GSYM"\n", damping);

  /* --------------------------------------
   * Create vectors for solution and scales
   * -------------------------------------- */

  u = N_VNew_Serial(NEQ);
  if (check_retval((void *)u, "N_VNew_Serial", 0)) return(1);

  scale = N_VClone(u);
  if (check_retval((void *)scale, "N_VClone", 0)) return(1);

  /* -----------------------------------------
   * Initialize and allocate memory for KINSOL
   * ----------------------------------------- */

  kmem = KINCreate();
  if (check_retval((void *)kmem, "KINCreate", 0)) return(1);

  /* Set number of prior residuals used in Anderson acceleration */
  retval = KINSetMAA(kmem, maa);

  retval = KINInit(kmem, FPFunction, u);
  if (check_retval(&retval, "KINInit", 1)) return(1);

  /* -------------------
   * Set optional inputs
   * ------------------- */

  /* Specify stopping tolerance based on residual */
  retval = KINSetFuncNormTol(kmem, tol);
  if (check_retval(&retval, "KINSetFuncNormTol", 1)) return(1);

  /* Set maximum number of iterations */
  retval = KINSetNumMaxIters(kmem, mxiter);
  if (check_retval(&retval, "KINSetNumMaxItersFuncNormTol", 1)) return(1);

  /* Set Anderson acceleration damping parameter */
  retval = KINSetDampingAA(kmem, damping);
  if (check_retval(&retval, "KINSetDampingAA", 1)) return(1);

  /* -------------
   * Initial guess
   * ------------- */

  /* Get vector data array */
  data = N_VGetArrayPointer(u);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return(1);

  data[0] =  PTONE;
  data[1] =  PTONE;
  data[2] = -PTONE;

  /* ----------------------------
   * Call KINSol to solve problem
   * ---------------------------- */

  /* No scaling used */
  N_VConst(ONE, scale);

  /* Call main solver */
  retval = KINSol(kmem,         /* KINSol memory block */
                  u,            /* initial guess on input; solution vector */
                  KIN_FP,       /* global strategy choice */
                  scale,        /* scaling vector, for the variable cc */
                  scale);       /* scaling vector for function values fval */
  if (check_retval(&retval, "KINSol", 1)) return(1);

  /* ------------------------------------
   * Get solver statistics
   * ------------------------------------ */

  /* get solver stats */
  retval = KINGetNumNonlinSolvIters(kmem, &nni);
  check_retval(&retval, "KINGetNumNonlinSolvIters", 1);

  retval = KINGetNumFuncEvals(kmem, &nfe);
  check_retval(&retval, "KINGetNumFuncEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("Number of nonlinear iterations: %6ld\n", nni);
  printf("Number of function evaluations: %6ld\n", nfe);

  /* ------------------------------------
   * Print solution and check error
   * ------------------------------------ */

  /* check solution */
  retval = check_ans(u, tol);

  /* -----------
   * Free memory
   * ----------- */

  N_VDestroy(u);
  N_VDestroy(scale);
  KINFree(&kmem);

  return(retval);
}

/* -----------------------------------------------------------------------------
 * Nonlinear system
 *
 * 3x - cos((y-1)z) - 1/2 = 0
 * x^2 - 81(y-0.9)^2 + sin(z) + 1.06 = 0
 * exp(-x(y-1)) + 20z + (10 pi - 3)/3 = 0
 *
 * Nonlinear fixed point function
 *
 * g1(x,y,z) = 1/3 cos((y-1)z) + 1/6
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) + 0.9
 * g3(x,y,z) = -1/20 exp(-x(y-1)) - (10 pi - 3) / 60
 *
 * ---------------------------------------------------------------------------*/
int FPFunction(N_Vector u, N_Vector g, void* user_data)
{
  realtype* udata = NULL;
  realtype* gdata = NULL;
  realtype  x, y, z;

  /* Get vector data arrays */
  udata = N_VGetArrayPointer(u);
  if (check_retval((void*)udata, "N_VGetArrayPointer", 0)) return(-1);

  gdata = N_VGetArrayPointer(g);
  if (check_retval((void*)gdata, "N_VGetArrayPointer", 0)) return(-1);

  x = udata[0];
  y = udata[1];
  z = udata[2];

  gdata[0] = (ONE/THREE) * COS((y-ONE)*z) + (ONE/SIX);
  gdata[1] = (ONE/NINE) * SQRT(x*x + SIN(z) + ONEPTZEROSIX) + PTNINE;
  gdata[2] = -(ONE/TWENTY) * EXP(-x*(y-ONE)) - (TEN * PI - THREE) / SIXTY;

  return(0);
}

/* -----------------------------------------------------------------------------
 * Check the solution of the nonlinear system and return PASS or FAIL
 * ---------------------------------------------------------------------------*/
static int check_ans(N_Vector u, realtype tol)
{
  realtype* data = NULL;
  realtype  ex, ey, ez;

  /* Get vector data array */
  data = N_VGetArrayPointer(u);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return(1);

  /* print the solution */
  printf("Computed solution:\n");
  printf("    x = %"GSYM"\n", data[0]);
  printf("    y = %"GSYM"\n", data[1]);
  printf("    z = %"GSYM"\n", data[2]);

  /* solution error */
  ex = ABS(data[0] - XTRUE);
  ey = ABS(data[1] - YTRUE);
  ez = ABS(data[2] - ZTRUE);

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
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if the function returned a NULL pointer -- no memory allocated */
  if (opt == 0) {
    if (returnvalue == NULL) {
      fprintf(stderr, "\nERROR: %s() failed -- returned NULL\n\n", funcname);
      return(1);
    } else {
      return(0);
    }
  }

  /* Check if the function returned an non-zero value -- internal failure */
  if (opt == 1) {
    errflag = (int *) returnvalue;
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
