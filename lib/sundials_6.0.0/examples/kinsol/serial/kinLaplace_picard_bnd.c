/* -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This example solves a 2D elliptic PDE
 *
 *    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u - 2.0
 *
 * subject to homogeneous Dirichlet boundary conditions.
 * The PDE is discretized on a uniform NX+2 by NY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving a system of size NEQ = NX*NY.
 * The nonlinear system is solved by KINSOL using the Picard
 * iteration and the SUNBAND linear solver.
 *
 * This file is strongly based on the kinLaplace_bnd.c file
 * developed by Radu Serban.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix        */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver  */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */
#include <sundials/sundials_math.h>    /* access to SUNRexp               */

/* Problem Constants */

#define NX   31             /* no. of points in x direction */
#define NY   31             /* no. of points in y direction */
#define NEQ  NX*NY          /* problem dimension */

#define SKIP 3              /* no. of points skipped for printing */

#define FTOL RCONST(1.e-12) /* function tolerance */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage.
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= NX, 1 <= j <= NY.
   The vdata array is obtained via the call vdata = N_VGetArrayPointer(v),
   where v is an N_Vector.
   The variables are ordered by the y index j, then by the x index i. */

#define IJth(vdata,i,j) (vdata[(j-1) + (i-1)*NY])

/* Private functions */

static int func(N_Vector u, N_Vector f, void *user_data);
static int jac(N_Vector u, N_Vector f, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2);
static void PrintOutput(N_Vector u);
static void PrintFinalStats(void *kmem);
static int check_retval(void *retvalvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main()
{
  SUNContext sunctx;
  realtype fnormtol, fnorm;
  N_Vector y, scale;
  int retval;
  void *kmem;
  SUNMatrix J;
  SUNLinearSolver LS;

  y = scale = NULL;
  kmem = NULL;
  J = NULL;
  LS = NULL;

  /* -------------------------
   * Print problem description
   * ------------------------- */

  printf("\n2D elliptic PDE on unit square\n");
  printf("   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0\n");
  printf(" + homogeneous Dirichlet boundary conditions\n\n");
  printf("Solution method: Anderson accelerated Picard iteration with band linear solver.\n");
  printf("Problem size: %2ld x %2ld = %4ld\n", (long int) NX, (long int) NY, (long int) NEQ);

  /* Create the SUNDIALS context that all SUNDIALS objects require */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* --------------------------------------
   * Create vectors for solution and scales
   * -------------------------------------- */

  y = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);

  scale = N_VNew_Serial(NEQ, sunctx);
  if (check_retval((void *)scale, "N_VNew_Serial", 0)) return(1);

  /* ----------------------------------------------------------------------------------
   * Initialize and allocate memory for KINSOL, set parametrs for Anderson acceleration
   * ---------------------------------------------------------------------------------- */

  kmem = KINCreate(sunctx);
  if (check_retval((void *)kmem, "KINCreate", 0)) return(1);

  /* y is used as a template */

  /* Use acceleration with up to 3 prior residuals */
  retval = KINSetMAA(kmem, 3);
  if (check_retval(&retval, "KINSetMAA", 1)) return(1);

  retval = KINInit(kmem, func, y);
  if (check_retval(&retval, "KINInit", 1)) return(1);

  /* -------------------
   * Set optional inputs
   * ------------------- */

  /* Specify stopping tolerance based on residual */

  fnormtol  = FTOL;
  retval = KINSetFuncNormTol(kmem, fnormtol);
  if (check_retval(&retval, "KINSetFuncNormTol", 1)) return(1);

  /* -------------------------
   * Create band SUNMatrix
   * ------------------------- */

  J = SUNBandMatrix(NEQ, NX, NX, sunctx);
  if(check_retval((void *)J, "SUNBandMatrix", 0)) return(1);

  /* ---------------------------
   * Create band SUNLinearSolver
   * --------------------------- */

  LS = SUNLinSol_Band(y, J, sunctx);
  if(check_retval((void *)LS, "SUNLinSol_Band", 0)) return(1);

  /* -------------------------
   * Attach band linear solver
   * ------------------------- */

  retval = KINSetLinearSolver(kmem, LS, J);
  if(check_retval(&retval, "KINSetLinearSolver", 1)) return(1);

  /* -------------------------
   * Set Jacobian function
   * ------------------------- */

  retval = KINSetJacFn(kmem, jac);
  if (check_retval(&retval, "KINSetJacFn", 1)) return(1);

  /* -------------
   * Initial guess
   * ------------- */

  N_VConst(ZERO, y);
  IJth(N_VGetArrayPointer(y), 2, 2) = ONE;

  /* ----------------------------
   * Call KINSol to solve problem
   * ---------------------------- */

  /* No scaling used */
  N_VConst(ONE,scale);

  /* Call main solver */
  retval = KINSol(kmem,           /* KINSol memory block */
                y,              /* initial guess on input; solution vector */
                KIN_PICARD,     /* global strategy choice */
                scale,          /* scaling vector, for the variable cc */
                scale);         /* scaling vector for function values fval */
  if (check_retval(&retval, "KINSol", 1)) return(1);


  /* ------------------------------------
   * Print solution and solver statistics
   * ------------------------------------ */

  /* Get scaled norm of the system function */

  retval = KINGetFuncNorm(kmem, &fnorm);
  if (check_retval(&retval, "KINGetfuncNorm", 1)) return(1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("\nComputed solution (||F|| = %Lg):\n\n",fnorm);
#else
  printf("\nComputed solution (||F|| = %g):\n\n",fnorm);
#endif
  PrintOutput(y);

  PrintFinalStats(kmem);

  /* -----------
   * Free memory
   * ----------- */

  N_VDestroy(y);
  N_VDestroy(scale);
  KINFree(&kmem);
  SUNLinSolFree(LS);
  SUNMatDestroy(J);
  SUNContext_Free(&sunctx);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * System function
 */

static int func(N_Vector u, N_Vector f, void *user_data)
{
  realtype dx, dy, hdiff, vdiff;
  realtype hdc, vdc;
  realtype uij, udn, uup, ult, urt;
  realtype *udata, *fdata;

  int i, j;

  dx = ONE/(NX+1);
  dy = ONE/(NY+1);
  hdc = ONE/(dx*dx);
  vdc = ONE/(dy*dy);

  udata = N_VGetArrayPointer(u);
  fdata = N_VGetArrayPointer(f);

  for (j=1; j <= NY; j++) {
    for (i=1; i <= NX; i++) {

      /* Extract u at x_i, y_j and four neighboring points */

      uij = IJth(udata, i, j);
      udn = (j == 1)  ? ZERO : IJth(udata, i, j-1);
      uup = (j == NY) ? ZERO : IJth(udata, i, j+1);
      ult = (i == 1)  ? ZERO : IJth(udata, i-1, j);
      urt = (i == NX) ? ZERO : IJth(udata, i+1, j);

      /* Evaluate diffusion components */

      hdiff = hdc*(ult - TWO*uij + urt);
      vdiff = vdc*(uup - TWO*uij + udn);

      /* Set residual at x_i, y_j */

      IJth(fdata, i, j) = hdiff + vdiff + uij - uij*uij*uij + 2.0;

    }
  }

  return(0);
}

/*
 * Jacobian function
 */

static int jac(N_Vector u, N_Vector f, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2)
{
  realtype dx, dy;
  realtype hdc, vdc;
  realtype *kthCol;

  int i, j, k;

  dx  = ONE/(NX+1);
  dy  = ONE/(NY+1);
  hdc = ONE/(dx*dx);
  vdc = ONE/(dy*dy);

  /*
     The components of f(t,u) which depend on u_{i,j} are
     f_{i,j}, f_{i-1,j}, f_{i+1,j}, f_{i,j+1}, and f_{i,j-1}.
     Thus, a column of the Jacobian will contain an entry from
     each of these equations exception the ones on the boundary.

     f_{i,j}   = hdc*(u_{i-1,j}  -2u_{i,j}  +u_{i+1,j})   + vdc*(u_{i,j-1}  -2u_{i,j}  +u_{i,j+1})
     f_{i-1,j} = hdc*(u_{i-2,j}  -2u_{i-1,j}+u_{i,j})     + vdc*(u_{i-1,j-1}-2u_{i-1,j}+u_{i-1,j+1})
     f_{i+1,j} = hdc*(u_{i,j}    -2u_{i+1,j}+u_{i+2,j})   + vdc*(u_{i+1,j-1}-2u_{i+1,j}+u_{i+1,j+1})
     f_{i,j-1} = hdc*(u_{i-1,j-1}-2u_{i,j-1}+u_{i+1,j-1}) + vdc*(u_{i,j-2}  -2u_{i,j-1}+u_{i,j})
     f_{i,j+1} = hdc*(u_{i-1,j+1}-2u_{i,j+1}+u_{i+1,j+1}) + vdc*(u_{i,j}    -2u_{i,j+1}+u_{i,j+2})
  */

  for (j=0; j <= NY-1; j++) {
    for (i=0; i <= NX-1; i++) {

      /* Evaluate diffusion coefficients */

      k = i + j*NX;
      kthCol = SUNBandMatrix_Column(J, k);
      SM_COLUMN_ELEMENT_B(kthCol,k,k) = -2.0*hdc - 2.0*vdc;
      if ( i != (NX-1) ) SM_COLUMN_ELEMENT_B(kthCol,k+1,k)  = hdc;
      if ( i != 0 )      SM_COLUMN_ELEMENT_B(kthCol,k-1,k)  = hdc;
      if ( j != (NY-1) ) SM_COLUMN_ELEMENT_B(kthCol,k+NX,k) = vdc;
      if ( j != 0 )      SM_COLUMN_ELEMENT_B(kthCol,k-NX,k) = vdc;
    }
  }

  return(0);
}

/*
 * Print solution at selected points
 */

static void PrintOutput(N_Vector u)
{
  int i, j;
  realtype dx, dy, x, y;
  realtype *udata;

  dx = ONE/(NX+1);
  dy = ONE/(NY+1);

  udata = N_VGetArrayPointer(u);

  printf("            ");
  for (i=1; i<=NX; i+= SKIP) {
    x = i*dx;
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("%-8.5Lf ", x);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("%-8.5f ", x);
#else
      printf("%-8.5f ", x);
#endif
  }
  printf("\n\n");

  for (j=1; j<=NY; j+= SKIP) {
    y = j*dy;
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("%-8.5Lf    ", y);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("%-8.5f    ", y);
#else
      printf("%-8.5f    ", y);
#endif
    for (i=1; i<=NX; i+= SKIP) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("%-8.5Lf ", IJth(udata,i,j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf("%-8.5f ", IJth(udata,i,j));
#else
      printf("%-8.5f ", IJth(udata,i,j));
#endif
    }
    printf("\n");
  }
}

/*
 * Print final statistics
 */

static void PrintFinalStats(void *kmem)
{
  long int nni, nfe, nje, nfeD;
  long int lenrwB, leniwB;
  int retval;

  /* Main solver statistics */

  retval = KINGetNumNonlinSolvIters(kmem, &nni);
  check_retval(&retval, "KINGetNumNonlinSolvIters", 1);
  retval = KINGetNumFuncEvals(kmem, &nfe);
  check_retval(&retval, "KINGetNumFuncEvals", 1);

  /* Band linear solver statistics */

  retval = KINGetNumJacEvals(kmem, &nje);
  check_retval(&retval, "KINGetNumJacEvals", 1);
  retval = KINGetNumLinFuncEvals(kmem, &nfeD);
  check_retval(&retval, "KINGetNumLinFuncEvals", 1);

  /* Band linear solver workspace size */

  retval = KINGetLinWorkSpace(kmem, &lenrwB, &leniwB);
  check_retval(&retval, "KINGetLinWorkSpace", 1);

  printf("\nFinal Statistics.. \n\n");
  printf("nni      = %6ld    nfe     = %6ld \n", nni, nfe);
  printf("nje      = %6ld    nfeB    = %6ld \n", nje, nfeD);
  printf("\n");
  printf("lenrwB   = %6ld    leniwB  = %6ld \n", lenrwB, leniwB);

}

/*
 * Check function return value...
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a retval so check if
 *             retval >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer
 */

static int check_retval(void *retvalvalue, const char *funcname, int opt)
{
  int *errretval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && retvalvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if retval < 0 */
  else if (opt == 1) {
    errretval = (int *) retvalvalue;
    if (*errretval < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *errretval);
      return(1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && retvalvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  return(0);
}
