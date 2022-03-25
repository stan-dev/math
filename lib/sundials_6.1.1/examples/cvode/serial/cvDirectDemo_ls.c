/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Demonstration program for CVODE - direct linear solvers.
 * Two separate problems are solved using both the CV_ADAMS and CV_BDF
 * linear multistep methods in combination with the
 * SUNNONLINSOL_FIXEDPOINT and SUNNONLINSOL_NEWTON nonlinear solver
 * modules:
 *
 * Problem 1: Van der Pol oscillator
 *   xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0.
 * This second-order ODE is converted to a first-order system by
 * defining y0 = x and y1 = xdot.
 * The NEWTON iteration cases use the following types of Jacobian
 * approximation: (1) dense, user-supplied, (2) dense, difference
 * quotient approximation, (3) diagonal approximation.
 *
 * Problem 2: ydot = A * y, where A is a banded lower triangular
 * matrix derived from 2-D advection PDE.
 * The NEWTON iteration cases use the following types of Jacobian
 * approximation: (1) band, user-supplied, (2) band, difference
 * quotient approximation, (3) diagonal approximation.
 *
 * For each problem, in the series of eight runs, CVodeInit is
 * called only once, for the first run, whereas CVodeReInit is
 * called for each of the remaining seven runs.
 *
 * Notes: This program demonstrates the usage of the sequential
 * macros NV_Ith_S, SM_ELEMENT_D, SM_COLUMN_B, and
 * SM_COLUMN_ELEMENT_B. The NV_Ith_S macro is used to reference the
 * components of an N_Vector. It works for any size N=NEQ, but
 * due to efficiency concerns it should only by used when the
 * problem size is small. The Problem 1 right hand side and
 * Jacobian functions f1 and Jac1 both use NV_Ith_S. The
 * N_VGetArrayPointer function gives the user access to the
 * memory used for the component storage of an N_Vector. In the
 * sequential case, the user may assume that this is one contiguous
 * array of reals. The N_VGetArrayPointer function
 * gives a more efficient means (than the NV_Ith_S macro) to
 * access the components of an N_Vector and should be used when the
 * problem size is large. The Problem 2 right hand side function f2
 * uses the N_VGetArrayPointer function. The SM_ELEMENT_D macro
 * used in Jac1 gives access to an element of a dense SUNMatrix. It
 * should be used only when the problem size is small (the
 * size of a Dense SUNMatrix is NEQ x NEQ) due to efficiency concerns. For
 * larger problem sizes, the macro SM_COLUMN_D can be used in order
 * to work directly with a column of a Dense SUNMatrix. The SM_COLUMN_B and
 * SM_COLUMN_ELEMENT_B allow efficient columnwise access to the elements
 * of a Banded SUNMatix. These macros are used in the Jac2 function.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvode/cvode.h>                          /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_serial.h>               /* access to serial N_Vector                    */
#include <sunmatrix/sunmatrix_dense.h>            /* access to dense SUNMatrix                    */
#include <sunlinsol/sunlinsol_dense.h>            /* access to dense SUNLinearSolver              */
#include <sunmatrix/sunmatrix_band.h>             /* access to band SUNMatrix                     */
#include <sunlinsol/sunlinsol_band.h>             /* access to band SUNLinearSolver               */
#include <cvode/cvode_diag.h>                     /* access to CVDIAG linear solver               */
#include "sunnonlinsol/sunnonlinsol_newton.h"     /* access to the newton SUNNonlinearSolver      */
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to the fixed point SUNNonlinearSolver */
#include <sundials/sundials_types.h>              /* definition of realtype                       */

/* helpful macros */

#ifndef SQR
#define SQR(A) ((A)*(A))
#endif

/* Shared Problem Constants */

#define ATOL RCONST(1.0e-6)
#define RTOL RCONST(0.0)

#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)
#define TWO    RCONST(2.0)
#define THIRTY RCONST(30.0)

/* Problem #1 Constants */

#define P1_NEQ        2
#define P1_ETA        RCONST(3.0)
#define P1_NOUT       4
#define P1_T0         RCONST(0.0)
#define P1_T1         RCONST(1.39283880203)
#define P1_DTOUT      RCONST(2.214773875)
#define P1_TOL_FACTOR RCONST(1.0e4)

/* Problem #2 Constants */

#define P2_MESHX      5
#define P2_MESHY      5
#define P2_NEQ        P2_MESHX*P2_MESHY
#define P2_ALPH1      RCONST(1.0)
#define P2_ALPH2      RCONST(1.0)
#define P2_NOUT       5
#define P2_ML         5
#define P2_MU         0
#define P2_T0         RCONST(0.0)
#define P2_T1         RCONST(0.01)
#define P2_TOUT_MULT  RCONST(10.0)
#define P2_TOL_FACTOR RCONST(1.0e3)

/* Linear Solver Options */

enum {FUNC, DENSE_USER, DENSE_DQ, DIAG, BAND_USER, BAND_DQ};

/* Private Helper Functions */

static int  Problem1(void);
static void PrintIntro1(void);
static void PrintHeader1(void);
static void PrintOutput1(realtype t, realtype y0, realtype y1, int qu, realtype hu);
static int  Problem2(void);
static void PrintIntro2(void);
static void PrintHeader2(void);
static void PrintOutput2(realtype t, realtype erm, int qu, realtype hu);
static realtype MaxError(N_Vector y, realtype t);
static int PrepareNextRun(SUNContext sunctx, void *cvode_mem, int lmm, int miter,
                          N_Vector y, SUNMatrix* A, sunindextype mu,
                          sunindextype ml, SUNLinearSolver* LS,
                          SUNNonlinearSolver* NLS);
static void PrintErrOutput(realtype tol_factor);
static void PrintFinalStats(void *cvode_mem, int miter, realtype ero);
static void PrintErrInfo(int nerr);

/* Functions Called by the Solver */

static int f1(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac1(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int f2(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac2(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Implementation */

int main(void)
{
  int nerr;

  nerr = Problem1();
  nerr += Problem2();
  PrintErrInfo(nerr);

  return(0);
}

static int Problem1(void)
{
  realtype reltol=RTOL, abstol=ATOL, t, tout, ero, er;
  int miter, retval, temp_retval, iout, nerr=0;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;
  void *cvode_mem;
  booleantype firstrun;
  int qu;
  realtype hu;
  SUNContext sunctx;

  y = NULL;
  A = NULL;
  LS = NULL;
  NLS = NULL;
  cvode_mem = NULL;

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

  y = N_VNew_Serial(P1_NEQ, sunctx);
  if(check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
  PrintIntro1();

  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  for (miter=FUNC; miter <= DIAG; miter++) {
    ero = ZERO;
    NV_Ith_S(y,0) = TWO;
    NV_Ith_S(y,1) = ZERO;

    firstrun = (miter==FUNC);
    if (firstrun) {

      /* initialize CVode */
      retval = CVodeInit(cvode_mem, f1, P1_T0, y);
      if(check_retval(&retval, "CVodeInit", 1)) return(1);

      /* set scalar tolerances */
      retval = CVodeSStolerances(cvode_mem, reltol, abstol);
      if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

    } else {

      /* reinitialize CVode */
      retval = CVodeReInit(cvode_mem, P1_T0, y);
      if(check_retval(&retval, "CVodeReInit", 1)) return(1);

    }

    retval = PrepareNextRun(sunctx, cvode_mem, CV_ADAMS, miter, y, &A, 0, 0, &LS, &NLS);
    if(check_retval(&retval, "PrepareNextRun", 1)) return(1);

    PrintHeader1();

    for(iout=1, tout=P1_T1; iout <= P1_NOUT; iout++, tout += P1_DTOUT) {
      retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      check_retval(&retval, "CVode", 1);
      temp_retval = CVodeGetLastOrder(cvode_mem, &qu);
      if(check_retval(&temp_retval, "CVodeGetLastOrder", 1)) ++nerr;
      temp_retval = CVodeGetLastStep(cvode_mem, &hu);
      if(check_retval(&temp_retval, "CVodeGetLastStep", 1)) ++nerr;
      PrintOutput1(t, NV_Ith_S(y,0), NV_Ith_S(y,1), qu, hu);
      if (retval != CV_SUCCESS) {
        nerr++;
        break;
      }
      if (iout%2 == 0) {
        er = fabs(NV_Ith_S(y,0)) / abstol;
        if (er > ero) ero = er;
        if (er > P1_TOL_FACTOR) {
          nerr++;
          PrintErrOutput(P1_TOL_FACTOR);
        }
      }
    }

    PrintFinalStats(cvode_mem, miter, ero);
  }

  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  NLS = NULL;
  LS = NULL;
  A = NULL;

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  for (miter=FUNC; miter <= DIAG; miter++) {
    ero = ZERO;
    NV_Ith_S(y,0) = TWO;
    NV_Ith_S(y,1) = ZERO;

    firstrun = (miter==FUNC);
    if (firstrun) {

      /* initialize CVode */
      retval = CVodeInit(cvode_mem, f1, P1_T0, y);
      if(check_retval(&retval, "CVodeInit", 1)) return(1);

      /* set scalar tolerances */
      retval = CVodeSStolerances(cvode_mem, reltol, abstol);
      if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

    } else {

      /* reinitialize CVode */
      retval = CVodeReInit(cvode_mem, P1_T0, y);
      if(check_retval(&retval, "CVodeReInit", 1)) return(1);

    }

    retval = PrepareNextRun(sunctx, cvode_mem, CV_BDF, miter, y, &A, 0, 0, &LS, &NLS);
    if(check_retval(&retval, "PrepareNextRun", 1)) return(1);

    PrintHeader1();

    for(iout=1, tout=P1_T1; iout <= P1_NOUT; iout++, tout += P1_DTOUT) {
      retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      check_retval(&retval, "CVode", 1);
      temp_retval = CVodeGetLastOrder(cvode_mem, &qu);
      if(check_retval(&temp_retval, "CVodeGetLastOrder", 1)) ++nerr;
      temp_retval = CVodeGetLastStep(cvode_mem, &hu);
      if(check_retval(&temp_retval, "CVodeGetLastStep", 1)) ++nerr;
      PrintOutput1(t, NV_Ith_S(y,0), NV_Ith_S(y,1), qu, hu);
      if (retval != CV_SUCCESS) {
        nerr++;
        break;
      }
      if (iout%2 == 0) {
        er = fabs(NV_Ith_S(y,0)) / abstol;
        if (er > ero) ero = er;
        if (er > P1_TOL_FACTOR) {
          nerr++;
          PrintErrOutput(P1_TOL_FACTOR);
        }
      }
    }

    PrintFinalStats(cvode_mem, miter, ero);
  }

  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  return(nerr);
}

static void PrintIntro1(void)
{
  printf("Demonstration program for CVODE package - direct linear solvers\n");
  printf("\n\n");
  printf("Problem 1: Van der Pol oscillator\n");
  printf(" xdotdot - 3*(1 - x^2)*xdot + x = 0, x(0) = 2, xdot(0) = 0\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" neq = %d,  reltol = %.2Lg,  abstol = %.2Lg",
         P1_NEQ, RTOL, ATOL);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf(" neq = %d,  reltol = %.2g,  abstol = %.2g",
         P1_NEQ, RTOL, ATOL);
#else
  printf(" neq = %d,  reltol = %.2g,  abstol = %.2g",
         P1_NEQ, RTOL, ATOL);
#endif
}

static void PrintHeader1(void)
{
  printf("\n     t           x              xdot         qu     hu \n");

  return;
}

static void PrintOutput1(realtype t, realtype y0, realtype y1, int qu, realtype hu)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%10.5Lf    %12.5Le   %12.5Le   %2d    %6.4Le\n", t, y0, y1, qu, hu);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%10.5f    %12.5e   %12.5e   %2d    %6.4e\n", t, y0, y1, qu, hu);
#else
  printf("%10.5f    %12.5e   %12.5e   %2d    %6.4e\n", t, y0, y1, qu, hu);
#endif

  return;
}

static int f1(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y0, y1;

  y0 = NV_Ith_S(y,0);
  y1 = NV_Ith_S(y,1);

  NV_Ith_S(ydot,0) = y1;
  NV_Ith_S(ydot,1) = (ONE - SQR(y0))* P1_ETA * y1 - y0;

  return(0);
}

static int Jac1(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y0, y1;

  y0 = NV_Ith_S(y,0);
  y1 = NV_Ith_S(y,1);

  SM_ELEMENT_D(J,0,1) = ONE;
  SM_ELEMENT_D(J,1,0) = -TWO * P1_ETA * y0 * y1 - ONE;
  SM_ELEMENT_D(J,1,1) = P1_ETA * (ONE - SQR(y0));

  return(0);
}

static int Problem2(void)
{
  realtype reltol=RTOL, abstol=ATOL, t, tout, er, erm, ero;
  int miter, retval, temp_retval, nerr=0;
  N_Vector y;
  SUNMatrix A;
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;
  void *cvode_mem;
  booleantype firstrun;
  int qu, iout;
  realtype hu;
  SUNContext sunctx;

  y = NULL;
  A = NULL;
  LS = NULL;
  NLS = NULL;
  cvode_mem = NULL;

  /* Create SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

  y = N_VNew_Serial(P2_NEQ, sunctx);
  if(check_retval((void *)y, "N_VNew_Serial", 0)) return(1);

  PrintIntro2();

  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  for (miter=FUNC; miter <= BAND_DQ; miter++) {
    if ((miter==DENSE_USER) || (miter==DENSE_DQ)) continue;
    ero = ZERO;
    N_VConst(ZERO, y);
    NV_Ith_S(y,0) = ONE;

    firstrun = (miter==FUNC);
    if (firstrun) {

      /* initialize CVode */
      retval = CVodeInit(cvode_mem, f2, P2_T0, y);
      if(check_retval(&retval, "CVodeInit", 1)) return(1);

      /* set scalar tolerances */
      retval = CVodeSStolerances(cvode_mem, reltol, abstol);
      if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

    } else {

      /* reinitialize CVode */
      retval = CVodeReInit(cvode_mem, P2_T0, y);
      if(check_retval(&retval, "CVodeReInit", 1)) return(1);

    }

    retval = PrepareNextRun(sunctx, cvode_mem, CV_ADAMS, miter, y, &A, P2_MU, P2_ML, &LS, &NLS);
    if(check_retval(&retval, "PrepareNextRun", 1)) return(1);

    PrintHeader2();

    for(iout=1, tout=P2_T1; iout <= P2_NOUT; iout++, tout*=P2_TOUT_MULT) {
      retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      check_retval(&retval, "CVode", 1);
      erm = MaxError(y, t);
      temp_retval = CVodeGetLastOrder(cvode_mem, &qu);
      if(check_retval(&temp_retval, "CVodeGetLastOrder", 1)) ++nerr;
      temp_retval = CVodeGetLastStep(cvode_mem, &hu);
      if(check_retval(&temp_retval, "CVodeGetLastStep", 1)) ++nerr;
      PrintOutput2(t, erm, qu, hu);
      if (retval != CV_SUCCESS) {
        nerr++;
        break;
      }
      er = erm / abstol;
        if (er > ero) ero = er;
        if (er > P2_TOL_FACTOR) {
          nerr++;
          PrintErrOutput(P2_TOL_FACTOR);
        }
    }

    PrintFinalStats(cvode_mem, miter, ero);
  }

  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  NLS = NULL;
  LS = NULL;
  A = NULL;

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  for (miter=FUNC; miter <= BAND_DQ; miter++) {
    if ((miter==DENSE_USER) || (miter==DENSE_DQ)) continue;
    ero = ZERO;
    N_VConst(ZERO, y);
    NV_Ith_S(y,0) = ONE;

    firstrun = (miter==FUNC);
    if (firstrun) {

      /* initialize CVode */
      retval = CVodeInit(cvode_mem, f2, P2_T0, y);
      if(check_retval(&retval, "CVodeInit", 1)) return(1);

      /* set scalar tolerances */
      retval = CVodeSStolerances(cvode_mem, reltol, abstol);
      if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

    } else {

      /* reinitialize CVode */
      retval = CVodeReInit(cvode_mem, P2_T0, y);
      if(check_retval(&retval, "CVodeReInit", 1)) return(1);

    }

    retval = PrepareNextRun(sunctx, cvode_mem, CV_BDF, miter, y, &A, P2_MU, P2_ML, &LS, &NLS);
    if(check_retval(&retval, "PrepareNextRun", 1)) return(1);

    PrintHeader2();

    for(iout=1, tout=P2_T1; iout <= P2_NOUT; iout++, tout*=P2_TOUT_MULT) {
      retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      check_retval(&retval, "CVode", 1);
      erm = MaxError(y, t);
      temp_retval = CVodeGetLastOrder(cvode_mem, &qu);
      if(check_retval(&temp_retval, "CVodeGetLastOrder", 1)) ++nerr;
      temp_retval = CVodeGetLastStep(cvode_mem, &hu);
      if(check_retval(&temp_retval, "CVodeGetLastStep", 1)) ++nerr;
      PrintOutput2(t, erm, qu, hu);
      if (retval != CV_SUCCESS) {
        nerr++;
        break;
      }
      er = erm / abstol;
        if (er > ero) ero = er;
        if (er > P2_TOL_FACTOR) {
          nerr++;
          PrintErrOutput(P2_TOL_FACTOR);
        }
    }

    PrintFinalStats(cvode_mem, miter, ero);
  }

  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  return(nerr);
}

static void PrintIntro2(void)
{
  printf("\n\n-------------------------------------------------------------");
  printf("\n-------------------------------------------------------------");
  printf("\n\nProblem 2: ydot = A * y, where A is a banded lower\n");
  printf("triangular matrix derived from 2-D advection PDE\n\n");
  printf(" neq = %d, ml = %d, mu = %d\n", P2_NEQ, P2_ML, P2_MU);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" itol = %s, reltol = %.2Lg, abstol = %.2Lg", "CV_SS", RTOL, ATOL);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf(" itol = %s, reltol = %.2g, abstol = %.2g", "CV_SS", RTOL, ATOL);
#else
  printf(" itol = %s, reltol = %.2g, abstol = %.2g", "CV_SS", RTOL, ATOL);
#endif
}

static void PrintHeader2(void)
{
  printf("\n      t        max.err      qu     hu \n");

  return;
}

static void PrintOutput2(realtype t, realtype erm, int qu, realtype hu)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%10.3Lf  %12.4Le   %2d   %12.4Le\n", t, erm, qu, hu);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%10.3f  %12.4e   %2d   %12.4e\n", t, erm, qu, hu);
#else
  printf("%10.3f  %12.4e   %2d   %12.4e\n", t, erm, qu, hu);
#endif

  return;
}

static int f2(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  sunindextype i, j, k;
  realtype d, *ydata, *dydata;

  ydata  = N_VGetArrayPointer(y);
  dydata = N_VGetArrayPointer(ydot);

  /*
     Excluding boundaries,

     ydot    = f    = -2 y    + alpha1 * y      + alpha2 * y
         i,j    i,j       i,j             i-1,j             i,j-1
  */

  for (j=0; j < P2_MESHY; j++) {
    for (i=0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      d = -TWO*ydata[k];
      if (i != 0) d += P2_ALPH1 * ydata[k-1];
      if (j != 0) d += P2_ALPH2 * ydata[k-P2_MESHX];
      dydata[k] = d;
    }
  }

  return(0);
}

static int Jac2(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int i, j, k;
  realtype *kthCol;

  /*
     The components of f(t,y) which depend on y    are
                                               i,j
     f    , f      , and f      :
      i,j    i+1,j        i,j+1

     f    = -2 y    + alpha1 * y      + alpha2 * y
      i,j       i,j             i-1,j             i,j-1

     f      = -2 y      + alpha1 * y    + alpha2 * y
      i+1,j       i+1,j             i,j             i+1,j-1

     f      = -2 y      + alpha1 * y        + alpha2 * y
      i,j+1       i,j+1             i-1,j+1             i,j
  */

  for (j=0; j < P2_MESHY; j++) {
    for (i=0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      kthCol = SM_COLUMN_B(J,k);
      SM_COLUMN_ELEMENT_B(kthCol,k,k) = -TWO;
      if (i != P2_MESHX-1) SM_COLUMN_ELEMENT_B(kthCol,k+1,k) = P2_ALPH1;
      if (j != P2_MESHY-1) SM_COLUMN_ELEMENT_B(kthCol,k+P2_MESHX,k) = P2_ALPH2;
    }
  }

  return(0);
}

static realtype MaxError(N_Vector y, realtype t)
{
  sunindextype i, j, k;
  realtype *ydata, er, ex=ZERO, yt, maxError=ZERO, ifact_inv, jfact_inv=ONE;

  if (t == ZERO) return(ZERO);

  ydata = N_VGetArrayPointer(y);
  if (t <= THIRTY) ex = exp(-TWO*t);

  for (j = 0; j < P2_MESHY; j++) {
    ifact_inv = ONE;
    for (i = 0; i < P2_MESHX; i++) {
      k = i + j * P2_MESHX;
      yt = pow(t, i+j) * ex * ifact_inv * jfact_inv;
      er = fabs(ydata[k] - yt);
      if (er > maxError) maxError = er;
      ifact_inv /= (i+1);
    }
    jfact_inv /= (j+1);
  }
  return(maxError);
}

static int PrepareNextRun(SUNContext sunctx, void *cvode_mem, int lmm, int miter,
                          N_Vector y, SUNMatrix* A, sunindextype mu,
                          sunindextype ml, SUNLinearSolver* LS,
                          SUNNonlinearSolver* NLS)
{
  int retval = CV_SUCCESS;

  if (*NLS)
    SUNNonlinSolFree(*NLS);
  if (*LS)
    SUNLinSolFree(*LS);
  if (*A)
    SUNMatDestroy(*A);

  printf("\n\n-------------------------------------------------------------");

  printf("\n\nLinear Multistep Method : ");
  if (lmm == CV_ADAMS) {
    printf("ADAMS\n");
  } else {
    printf("BDF\n");
  }

  printf("Iteration               : ");
  if (miter == FUNC) {
    printf("FIXEDPOINT\n");

    /* create fixed point nonlinear solver object */
    *NLS = SUNNonlinSol_FixedPoint(y, 0, sunctx);
    if(check_retval((void *)*NLS, "SUNNonlinSol_FixedPoint", 0)) return(1);

    /* attach nonlinear solver object to CVode */
    retval = CVodeSetNonlinearSolver(cvode_mem, *NLS);
    if(check_retval(&retval, "CVodeSetNonlinearSolver", 1)) return(1);

  } else {
    printf("NEWTON\n");

    /* create Newton nonlinear solver object */
    *NLS = SUNNonlinSol_Newton(y, sunctx);
    if(check_retval((void *)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* attach nonlinear solver object to CVode */
    retval = CVodeSetNonlinearSolver(cvode_mem, *NLS);
    if(check_retval(&retval, "CVodeSetNonlinearSolver", 1)) return(1);

    printf("Linear Solver           : ");

    switch(miter) {

    case DENSE_USER :
      printf("Dense, User-Supplied Jacobian\n");

      /* Create dense SUNMatrix for use in linear solves */
      *A = SUNDenseMatrix(P1_NEQ, P1_NEQ, sunctx);
      if(check_retval((void *)*A, "SUNDenseMatrix", 0)) return(1);

      /* Create dense SUNLinearSolver object for use by CVode */
      *LS = SUNLinSol_Dense(y, *A, sunctx);
      if(check_retval((void *)*LS, "SUNLinSol_Dense", 0)) return(1);

      /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
      retval = CVodeSetLinearSolver(cvode_mem, *LS, *A);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

      /* Set the user-supplied Jacobian routine Jac */
      retval = CVodeSetJacFn(cvode_mem, Jac1);
      if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
      break;

    case DENSE_DQ :
      printf("Dense, Difference Quotient Jacobian\n");

      /* Create dense SUNMatrix for use in linear solves */
      *A = SUNDenseMatrix(P1_NEQ, P1_NEQ, sunctx);
      if(check_retval((void *)*A, "SUNDenseMatrix", 0)) return(1);

      /* Create dense SUNLinearSolver object for use by CVode */
      *LS = SUNLinSol_Dense(y, *A, sunctx);
      if(check_retval((void *)*LS, "SUNLinSol_Dense", 0)) return(1);

      /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
      retval = CVodeSetLinearSolver(cvode_mem, *LS, *A);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

      /* Use a difference quotient Jacobian */
      retval = CVodeSetJacFn(cvode_mem, NULL);
      if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
      break;

    case DIAG :
      printf("Diagonal Jacobian\n");

      /* Call CVDiag to create/attach the CVODE-specific diagonal solver */
      retval = CVDiag(cvode_mem);
      if(check_retval(&retval, "CVDiag", 1)) return(1);
      break;

    case BAND_USER :
      printf("Band, User-Supplied Jacobian\n");

      /* Create band SUNMatrix for use in linear solves */
      *A = SUNBandMatrix(P2_NEQ, mu, ml, sunctx);
      if(check_retval((void *)*A, "SUNBandMatrix", 0)) return(1);

      /* Create banded SUNLinearSolver object for use by CVode */
      *LS = SUNLinSol_Band(y, *A, sunctx);
      if(check_retval((void *)*LS, "SUNLinSol_Band", 0)) return(1);

      /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
      retval = CVodeSetLinearSolver(cvode_mem, *LS, *A);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

      /* Set the user-supplied Jacobian routine Jac */
      retval = CVodeSetJacFn(cvode_mem, Jac2);
      if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
      break;

    case BAND_DQ  :
      printf("Band, Difference Quotient Jacobian\n");

      /* Create band SUNMatrix for use in linear solves */
      *A = SUNBandMatrix(P2_NEQ, mu, ml, sunctx);
      if(check_retval((void *)*A, "SUNBandMatrix", 0)) return(1);

      /* Create banded SUNLinearSolver object for use by CVode */
      *LS = SUNLinSol_Band(y, *A, sunctx);
      if(check_retval((void *)*LS, "SUNLinSol_Band", 0)) return(1);

      /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
      retval = CVodeSetLinearSolver(cvode_mem, *LS, *A);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

      /* Use a difference quotient Jacobian */
      retval = CVodeSetJacFn(cvode_mem, NULL);
      if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
      break;
    }
  }

  return(retval);
}

static void PrintErrOutput(realtype tol_factor)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("\n\n Error exceeds %Lg * tolerance \n\n", tol_factor);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#else
  printf("\n\n Error exceeds %g * tolerance \n\n", tol_factor);
#endif

  return;
}

static void PrintFinalStats(void *cvode_mem, int miter, realtype ero)
{
  long int lenrw, leniw, lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf, nje, nfeLS;
  int retval;

  retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_retval(&retval, "CVodeGetWorkSpace", 1);
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

  printf("\n Final statistics for this run:\n\n");
  printf(" CVode real workspace length              = %4ld \n",  lenrw);
  printf(" CVode integer workspace length           = %4ld \n",  leniw);
  printf(" Number of steps                          = %4ld \n",  nst);
  printf(" Number of f-s                            = %4ld \n",  nfe);
  printf(" Number of setups                         = %4ld \n",  nsetups);
  printf(" Number of nonlinear iterations           = %4ld \n",  nni);
  printf(" Number of nonlinear convergence failures = %4ld \n",  ncfn);
  printf(" Number of error test failures            = %4ld \n\n",netf);

  if (miter != FUNC) {
    if (miter != DIAG) {
      retval = CVodeGetNumJacEvals(cvode_mem, &nje);
      check_retval(&retval, "CVodeGetNumJacEvals", 1);
      retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
      check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);
      retval = CVodeGetLinWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
      check_retval(&retval, "CVodeGetLinWorkSpace", 1);
    } else {
      nje = nsetups;
      retval = CVDiagGetNumRhsEvals(cvode_mem, &nfeLS);
      check_retval(&retval, "CVDiagGetNumRhsEvals", 1);
      retval = CVDiagGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
      check_retval(&retval, "CVDiagGetWorkSpace", 1);
    }
    printf(" Linear solver real workspace length      = %4ld \n", lenrwLS);
    printf(" Linear solver integer workspace length   = %4ld \n", leniwLS);
    printf(" Number of Jacobian evaluations           = %4ld \n", nje);
    printf(" Number of f evals. in linear solver      = %4ld \n\n", nfeLS);
  }

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" Error overrun = %.3Lf \n", ero);
#else
  printf(" Error overrun = %.3f \n", ero);
#endif
}

static void PrintErrInfo(int nerr)
{
  printf("\n\n-------------------------------------------------------------");
  printf("\n-------------------------------------------------------------");
  printf("\n\n Number of errors encountered = %d \n", nerr);

  return;
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval < 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}
