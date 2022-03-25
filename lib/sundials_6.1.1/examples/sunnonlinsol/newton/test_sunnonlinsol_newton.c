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
 * This is the testing routine to check the SUNNonlinearSolver Newton module
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "sundials/sundials_types.h"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#define NEQ   3                /* number of equations        */
#define TOL   RCONST(1.0e-2)   /* nonlinear solver tolerance */
#define MAXIT 10               /* max nonlinear iterations   */

#define ZERO  RCONST(0.0)  /* real 0.0 */
#define HALF  RCONST(0.5)  /* real 0.5 */
#define ONE   RCONST(1.0)  /* real 1.0 */
#define TWO   RCONST(2.0)  /* real 2.0 */
#define THREE RCONST(3.0)  /* real 3.0 */
#define FOUR  RCONST(4.0)  /* real 4.0 */
#define SIX   RCONST(6.0)  /* real 6.0 */

/* approximate solution */
#define Y1 0.785196933062355226
#define Y2 0.496611392944656396
#define Y3 0.369922830745872357

/* Check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

/* Nonlinear residual function */
static int Res(N_Vector y, N_Vector f, void *mem);

/* Jacobian of the nonlinear residual */
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * Proxies for integrator memory struct and functions
 */

/* Integrator memory structure */
typedef struct IntegratorMemRec {
  N_Vector y0;
  N_Vector ycur;
  N_Vector ycor;
  N_Vector w;
  N_Vector x;
  SUNMatrix A;
  SUNLinearSolver LS;
} *IntegratorMem;

/* Linear solver setup interface function */
static int LSetup(booleantype jbad, booleantype* jcur, void* mem);

/* Linear solver solve interface function */
static int LSolve(N_Vector b, void* mem);

/* Convergence test function */
static int ConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                    realtype tol, N_Vector ewt, void* mem);

/* -----------------------------------------------------------------------------
 * Main testing routine
 * ---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

  IntegratorMem      Imem;       /* proxy for integrator memory */
  SUNNonlinearSolver NLS;        /* nonlinear solver object     */
  long int           niters;     /* number of nonlinear iters   */
  int                retval = 0; /* return value                */
  SUNContext         sunctx;

  /* create SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* create proxy for integrator memory */
  Imem = (IntegratorMem) malloc(sizeof(struct IntegratorMemRec));

  /* create vector */
  Imem->y0 = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)Imem->y0, "N_VNew_Serial", 0)) return(1);

  Imem->ycur = N_VClone(Imem->y0);
  if (check_retval((void *)Imem->ycur, "N_VClone", 0)) return(1);

  Imem->ycor = N_VClone(Imem->y0);
  if (check_retval((void *)Imem->ycor, "N_VClone", 0)) return(1);

  Imem->w = N_VClone(Imem->y0);
  if (check_retval((void *)Imem->w, "N_VClone", 0)) return(1);

  Imem->x = N_VClone(Imem->y0);
  if (check_retval((void *)Imem->x, "N_VClone", 0)) return(1);

  /* set initial guess for the state */
  NV_Ith_S(Imem->y0,0) = HALF;
  NV_Ith_S(Imem->y0,1) = HALF;
  NV_Ith_S(Imem->y0,2) = HALF;

  /* set initial guess for the correction */
  NV_Ith_S(Imem->ycor,0) = ZERO;
  NV_Ith_S(Imem->ycor,1) = ZERO;
  NV_Ith_S(Imem->ycor,2) = ZERO;

  /* set weights for norm */
  NV_Ith_S(Imem->w,0) = ONE;
  NV_Ith_S(Imem->w,1) = ONE;
  NV_Ith_S(Imem->w,2) = ONE;

  /* create dense matrix */
  Imem->A = SUNDenseMatrix(NEQ, NEQ, sunctx);
  if (check_retval((void *)Imem->A, "SUNDenseMatrix", 0)) return(1);

  /* create dense linear solver */
  Imem->LS = SUNLinSol_Dense(Imem->y0, Imem->A, sunctx);
  if (check_retval((void *)Imem->LS, "SUNLinSol_Dense", 0)) return(1);

  /* initialize the linear solver */
  retval = SUNLinSolInitialize(Imem->LS);
  if (check_retval(&retval, "SUNLinSolInitialize", 1)) return(1);

  /* create nonlinear solver */
  NLS = SUNNonlinSol_Newton(Imem->y0, sunctx);
  if (check_retval((void *)NLS, "SUNNonlinSol_Newton", 0)) return(1);

  /* set the nonlinear residual function */
  retval = SUNNonlinSolSetSysFn(NLS, Res);
  if (check_retval(&retval, "SUNNonlinSolSetSysFn", 1)) return(1);

  /* set the wrapper functions to linear solver setup and solve functions */
  retval = SUNNonlinSolSetLSetupFn(NLS, LSetup);
  if (check_retval(&retval, "SUNNonlinSolSetSetupFn", 1)) return(1);

  retval = SUNNonlinSolSetLSolveFn(NLS, LSolve);
  if (check_retval(&retval, "SUNNonlinSolSetSolveFn", 1)) return(1);

  retval = SUNNonlinSolSetConvTestFn(NLS, ConvTest, NULL);
  if (check_retval(&retval, "SUNNonlinSolSetConvTestFn", 1)) return(1);

  /* set the maximum number of nonlinear iterations */
  retval = SUNNonlinSolSetMaxIters(NLS, MAXIT);
  if (check_retval(&retval, "SUNNonlinSolSetMaxIters", 1)) return(1);

  /* solve the nonlinear system */
  retval = SUNNonlinSolSolve(NLS, Imem->y0, Imem->ycor, Imem->w, TOL, SUNTRUE,
                             Imem);
  if (check_retval(&retval, "SUNNonlinSolSolve", 1)) return(1);

  /* update the initial guess with the final correction */
  N_VLinearSum(ONE, Imem->y0, ONE, Imem->ycor, Imem->ycur);

  /* print the solution */
  printf("Solution:\n");
  printf("y1 = %"GSYM"\n",NV_Ith_S(Imem->ycur, 0));
  printf("y2 = %"GSYM"\n",NV_Ith_S(Imem->ycur, 1));
  printf("y3 = %"GSYM"\n",NV_Ith_S(Imem->ycur, 2));

  /* print the solution error */
  printf("Solution Error:\n");
  printf("e1 = %"GSYM"\n",NV_Ith_S(Imem->ycur, 0) - Y1);
  printf("e2 = %"GSYM"\n",NV_Ith_S(Imem->ycur, 1) - Y2);
  printf("e3 = %"GSYM"\n",NV_Ith_S(Imem->ycur, 2) - Y3);

  /* get the number of linear iterations */
  retval = SUNNonlinSolGetNumIters(NLS, &niters);
  if (check_retval(&retval, "SUNNonlinSolGetNumIters", 1)) return(1);

  printf("Number of nonlinear iterations: %ld\n",niters);

  /* Free vector, matrix, linear solver, and nonlinear solver */
  N_VDestroy(Imem->y0);
  N_VDestroy(Imem->ycur);
  N_VDestroy(Imem->ycor);
  N_VDestroy(Imem->w);
  N_VDestroy(Imem->x);
  SUNMatDestroy(Imem->A);
  SUNLinSolFree(Imem->LS);
  SUNNonlinSolFree(NLS);
  free(Imem);
  SUNContext_Free(&sunctx);

  /* Print result */
  if (retval) {
    printf("FAIL\n");
  } else {
    printf("SUCCESS\n");
  }

  return(retval);
}


/* Proxy for integrator lsetup function */
int LSetup(booleantype jbad, booleantype* jcur, void* mem)
{
  int retval;
  IntegratorMem Imem;

  if (mem == NULL) {
    printf("ERROR: Integrator memory is NULL");
    return(-1);
  }
  Imem = (IntegratorMem) mem;

  /* compute the Jacobian */
  retval = Jac(ZERO, Imem->ycur, NULL, Imem->A, NULL, NULL, NULL, NULL);
  if (retval != 0) return(retval);

  /* update Jacobian status */
  *jcur = SUNTRUE;

  /* setup the linear solver */
  retval = SUNLinSolSetup(Imem->LS, Imem->A);

  return(retval);
}


/* Proxy for integrator lsolve function */
int LSolve(N_Vector b, void* mem)
{
  int retval;
  IntegratorMem Imem;

  if (mem == NULL) {
    printf("ERROR: Integrator memory is NULL");
    return(-1);
  }
  Imem = (IntegratorMem) mem;

  retval = SUNLinSolSolve(Imem->LS, Imem->A, Imem->x, b, ZERO);
  N_VScale(ONE, Imem->x, b);

  return(retval);
}


/* Proxy for integrator convergence test function */
int ConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del, realtype tol,
             N_Vector ewt, void* mem)
{
  realtype delnrm;

  /* compute the norm of the correction */
  delnrm = N_VWrmsNorm(del, ewt);

  if (delnrm <= tol) return(SUN_NLS_SUCCESS);  /* success       */
  else               return(SUN_NLS_CONTINUE); /* not converged */
}


/* -----------------------------------------------------------------------------
 * Nonlinear residual function
 *
 * f1(x,y,z) = x^2 + y^2 + z^2 - 1 = 0
 * f2(x,y,z) = 2x^2 + y^2 - 4z     = 0
 * f3(x,y,z) = 3x^2 - 4y + z^2     = 0
 *
 * ---------------------------------------------------------------------------*/
int Res(N_Vector ycor, N_Vector f, void *mem)
{
  IntegratorMem Imem;
  realtype y1, y2, y3;

  if (mem == NULL) {
    printf("ERROR: Integrator memory is NULL");
    return(-1);
  }
  Imem = (IntegratorMem) mem;

  /* update state based on current correction */
  N_VLinearSum(ONE, Imem->y0, ONE, Imem->ycor, Imem->ycur);

  /* get vector components */
  y1 = NV_Ith_S(Imem->ycur,0);
  y2 = NV_Ith_S(Imem->ycur,1);
  y3 = NV_Ith_S(Imem->ycur,2);

  /* compute the residual function */
  NV_Ith_S(f,0) = y1*y1 + y2*y2 + y3*y3 - ONE;
  NV_Ith_S(f,1) = TWO * y1*y1 + y2*y2 - FOUR * y3;
  NV_Ith_S(f,2) = THREE * (y1*y1) - FOUR * y2 + y3*y3;

  /* return success */
  return(0);
}


/* -----------------------------------------------------------------------------
 * Jacobian of the nonlinear residual function
 *
 *            ( 2x  2y  2z )
 * J(x,y,z) = ( 4x  2y  -4 )
 *            ( 6x  -4  2z )
 *
 * ---------------------------------------------------------------------------*/
int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;

  y1 = NV_Ith_S(y,0);
  y2 = NV_Ith_S(y,1);
  y3 = NV_Ith_S(y,2);

  SM_ELEMENT_D(J,0,0) = TWO*y1;
  SM_ELEMENT_D(J,0,1) = TWO*y2;
  SM_ELEMENT_D(J,0,2) = TWO*y3;

  SM_ELEMENT_D(J,1,0) = FOUR*y1;
  SM_ELEMENT_D(J,1,1) = TWO*y2;
  SM_ELEMENT_D(J,1,2) = -FOUR;

  SM_ELEMENT_D(J,2,0) = SIX*y1;
  SM_ELEMENT_D(J,2,1) = -FOUR;
  SM_ELEMENT_D(J,2,2) = TWO*y3;

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
