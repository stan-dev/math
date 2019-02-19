/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the testing routine to check the SUNNonlinearSolver fixed point
 * module
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sundials/sundials_types.h"
#include "nvector/nvector_serial.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"

#define NEQ   3                /* number of equations        */
#define TOL   RCONST(1.0e-4)   /* nonlinear solver tolerance */
#define MAXIT 10               /* max nonlinear iterations   */

#define ZERO         RCONST(0.0)             /* real 0.0  */
#define PTONE        RCONST(0.1)             /* real 0.1  */
#define HALF         RCONST(0.5)             /* real 0.5  */
#define ONE          RCONST(1.0)             /* real 1.0  */
#define ONEPTZEROSIX RCONST(1.06)            /* real 1.06 */
#define THREE        RCONST(3.0)             /* real 3.0  */
#define FOUR         RCONST(4.0)             /* real 4.0  */
#define SIX          RCONST(6.0)             /* real 6.0  */
#define NINE         RCONST(9.0)             /* real 9.0  */
#define TEN          RCONST(10.0)            /* real 10.0 */
#define TWENTY       RCONST(20.0)            /* real 20.0 */
#define SIXTY        RCONST(60.0)            /* real 60.0 */
#define EIGHTYONE    RCONST(81.0)            /* real 81.0 */
#define PI           RCONST(3.1415926535898) /* real pi   */

/* approximate solution */
#define Y1 HALF
#define Y2 ZERO
#define Y3 -PI/SIX

/* Check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

/* Nonlinear fixed point function */
static int FPFunction(N_Vector y, N_Vector f, void *mem);

/* Convergence test function */
static int ConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                    realtype tol, N_Vector ewt, void* mem);

/* -----------------------------------------------------------------------------
 * Main testing routine
 * ---------------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int                retval = 0;
  N_Vector           x, y0, y, w;
  SUNNonlinearSolver NLS;
  long int           niters;

  /* create vectors */
  x  = N_VNew_Serial(NEQ);
  if (check_retval((void *)x, "N_VNew_Serial", 0)) return(1);
  
  y0 = N_VClone(x);
  if (check_retval((void *)y0, "N_VNew_Serial", 0)) return(1);

  y  = N_VClone(x);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);

  w  = N_VClone(x);
  if (check_retval((void *)w, "N_VNew_Serial", 0)) return(1);

  /* set initial guess */
  NV_Ith_S(y0,0) = PTONE;
  NV_Ith_S(y0,1) = PTONE;
  NV_Ith_S(y0,2) = -PTONE;

  /* set weights */
  NV_Ith_S(w,0) = ONE;
  NV_Ith_S(w,1) = ONE;
  NV_Ith_S(w,2) = ONE;

  /* create nonlinear solver */
  NLS = SUNNonlinSol_FixedPoint(y, 0);
  if (check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0)) return(1);

  /* set the nonlinear residual function */
  retval = SUNNonlinSolSetSysFn(NLS, FPFunction);
  if (check_retval(&retval, "SUNNonlinSolSetSysFn", 1)) return(1);

  /* set the convergence test function */
  retval = SUNNonlinSolSetConvTestFn(NLS, ConvTest);
  if (check_retval(&retval, "SUNNonlinSolSetConvTestFn", 1)) return(1);

  /* set the maximum number of nonlinear iterations */
  retval = SUNNonlinSolSetMaxIters(NLS, MAXIT);
  if (check_retval(&retval, "SUNNonlinSolSetMaxIters", 1)) return(1);

  /* solve the nonlinear system */
  retval = SUNNonlinSolSolve(NLS, y0, y, w, TOL, SUNTRUE, x);
  if (check_retval(&retval, "SUNNonlinSolSolve", 1)) return(1);

  /* print the solution */
  printf("Solution:\n");
  printf("y1 = %g\n",NV_Ith_S(y,0));
  printf("y2 = %g\n",NV_Ith_S(y,1));
  printf("y3 = %g\n",NV_Ith_S(y,2));

  /* print the solution error */
  printf("Solution Error:\n");
  printf("e1 = %g\n",NV_Ith_S(y,0) - Y1);
  printf("e2 = %g\n",NV_Ith_S(y,1) - Y2);
  printf("e3 = %g\n",NV_Ith_S(y,2) - Y3);

  /* get the number of linear iterations */
  retval = SUNNonlinSolGetNumIters(NLS, &niters);
  if (check_retval(&retval, "SUNNonlinSolGetNumIters", 1)) return(1);

  printf("Number of nonlinear iterations: %ld\n",niters);

  /* Free vector, matrix, linear solver, and nonlinear solver */
  N_VDestroy(x);
  N_VDestroy(y0);
  N_VDestroy(y);
  N_VDestroy(w);
  SUNNonlinSolFree(NLS);

  /* Print result */
  if (retval) {
    printf("FAIL\n");
  } else {
    printf("SUCCESS\n");
  }

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
 * Nonlinear system
 *
 * 3x - cos(yz) - 1/2 = 0
 * x^2 - 81(y+0.1)^2 + sin(z) + 1.06 = 0
 * exp(-xy) + 20z + (10 pi - 3)/3 = 0
 *
 * Nonlinear fixed point function
 *
 * g1(x,y,z) = 1/3 cos(yz) + 1/6
 * g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) - 0.1
 * g3(x,y,z) = -1/20 exp(-xy) - (10 pi - 3) / 60
 *
 * ---------------------------------------------------------------------------*/
int FPFunction(N_Vector y, N_Vector f, void *mem)
{
  realtype y1, y2, y3;

  y1 = NV_Ith_S(y,0);
  y2 = NV_Ith_S(y,1);
  y3 = NV_Ith_S(y,2);

  NV_Ith_S(f,0) = (ONE/THREE) * cos(y2*y3) + (ONE/SIX);
  NV_Ith_S(f,1) = (ONE/NINE) * sqrt(y1*y1 + sin(y3) + ONEPTZEROSIX) - PTONE;
  NV_Ith_S(f,2) = -(ONE/TWENTY) * exp(-y1*y2) - (TEN * PI - THREE) / SIXTY;

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
