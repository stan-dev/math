/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * Based on kinLaplace_picard_bnd.c by Carol S. Woodward @ LLNL
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
 * iteration and the SPGMR linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver */
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
static int jactimes(N_Vector v, N_Vector Jv, N_Vector u, booleantype *new_u,
                    void *user_data);
static void PrintOutput(N_Vector u);
static void PrintFinalStats(void *kmem);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main()
{
  realtype fnormtol, fnorm;
  N_Vector y, scale;
  int flag;
  void *kmem;
  SUNLinearSolver LS;
  FILE *infofp;

  y = scale = NULL;
  kmem = NULL;
  LS = NULL;

  /* -------------------------
   * Print problem description
   * ------------------------- */

  printf("\n2D elliptic PDE on unit square\n");
  printf("   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0\n");
  printf(" + homogeneous Dirichlet boundary conditions\n\n");
  printf("Solution method: Anderson accelerated Picard iteration with SPGMR linear solver.\n");
  printf("Problem size: %2ld x %2ld = %4ld\n", (long int) NX, (long int) NY, (long int) NEQ);

  /* --------------------------------------
   * Create vectors for solution and scales
   * -------------------------------------- */

  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  scale = N_VNew_Serial(NEQ);
  if (check_flag((void *)scale, "N_VNew_Serial", 0)) return(1);

  /* ----------------------------------------------------------------------------------
   * Initialize and allocate memory for KINSOL, set parametrs for Anderson acceleration
   * ---------------------------------------------------------------------------------- */

  kmem = KINCreate();
  if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

  /* y is used as a template */

  /* Use acceleration with up to 3 prior residuals */
  flag = KINSetMAA(kmem, 3);
  if (check_flag(&flag, "KINSetMAA", 1)) return(1);

  flag = KINInit(kmem, func, y);
  if (check_flag(&flag, "KINInit", 1)) return(1);

  /* -------------------
   * Set optional inputs
   * ------------------- */

  /* Specify stopping tolerance based on residual */

  fnormtol  = FTOL;
  flag = KINSetFuncNormTol(kmem, fnormtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);

  /* Set information file */

  infofp = fopen("KINSOL.log", "w");

  flag = KINSetInfoFile(kmem, infofp);
  if (check_flag(&flag, "KINSetInfoFile", 1)) return(1);

  flag = KINSetPrintLevel(kmem, 3);
  if (check_flag(&flag, "KINSetPrintLevel", 1)) return(1);

  /* ----------------------
   * Create SUNLinearSolver
   * ---------------------- */

  LS = SUNLinSol_SPGMR(y, PREC_NONE, 10);
  if(check_flag((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

  /* --------------------
   * Attach linear solver
   * -------------------- */

  flag = KINSetLinearSolver(kmem, LS, NULL);
  if(check_flag(&flag, "KINSetLinearSolver", 1)) return(1);

  /* ------------------------------------
   * Set Jacobian vector product function
   * ------------------------------------ */

  flag = KINSetJacTimesVecFn(kmem, jactimes);
  if (check_flag(&flag, "KINSetJacTimesVecFn", 1)) return(1);

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
  flag = KINSol(kmem,           /* KINSol memory block */
                y,              /* initial guess on input; solution vector */
                KIN_PICARD,     /* global strategy choice */
                scale,          /* scaling vector, for the variable cc */
                scale);         /* scaling vector for function values fval */
  if (check_flag(&flag, "KINSol", 1)) return(1);


  /* ------------------------------------
   * Print solution and solver statistics
   * ------------------------------------ */

  /* Get scaled norm of the system function */

  flag = KINGetFuncNorm(kmem, &fnorm);
  if (check_flag(&flag, "KINGetfuncNorm", 1)) return(1);

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

  fclose(infofp);
  N_VDestroy(y);
  N_VDestroy(scale);
  KINFree(&kmem);
  SUNLinSolFree(LS);

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
 * Jacobian vector product function
 */

static int jactimes(N_Vector v, N_Vector Jv, N_Vector u, booleantype *new_u,
                    void *user_data)
{
  realtype dx, dy, hdiff, vdiff;
  realtype hdc, vdc;
  realtype vij, vdn, vup, vlt, vrt;
  realtype *vdata, *Jvdata;

  int i, j;

  dx  = ONE/(NX+1);
  dy  = ONE/(NY+1);
  hdc = ONE/(dx*dx);
  vdc = ONE/(dy*dy);

  vdata  = N_VGetArrayPointer(v);
  Jvdata = N_VGetArrayPointer(Jv);

  for (j=1; j <= NY; j++) {
    for (i=1; i <= NX; i++) {

      /* Extract v at x_i, y_j and four neighboring points */

      vij = IJth(vdata, i, j);
      vdn = (j == 1)  ? ZERO : IJth(vdata, i, j-1);
      vup = (j == NY) ? ZERO : IJth(vdata, i, j+1);
      vlt = (i == 1)  ? ZERO : IJth(vdata, i-1, j);
      vrt = (i == NX) ? ZERO : IJth(vdata, i+1, j);

      /* Evaluate diffusion components */

      hdiff = hdc*(vlt - TWO*vij + vrt);
      vdiff = vdc*(vup - TWO*vij + vdn);

      /* Set Jv at x_i, y_j */

      IJth(Jvdata, i, j) = hdiff + vdiff;

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
  long int nni, nfe, nli, npe, nps, ncfl, nfeLS, njvevals;
  long int lenrw, leniw, lenrwLS, leniwLS;
  int flag;

  /* Main solver statistics */

  flag = KINGetNumNonlinSolvIters(kmem, &nni);
  check_flag(&flag, "KINGetNumNonlinSolvIters", 1);
  flag = KINGetNumFuncEvals(kmem, &nfe);
  check_flag(&flag, "KINGetNumFuncEvals", 1);

  /* Linear solver statistics */

  flag = KINGetNumLinIters(kmem, &nli);
  check_flag(&flag, "KINGetNumLinIters", 1);
  flag = KINGetNumLinFuncEvals(kmem, &nfeLS);
  check_flag(&flag, "KINGetNumLinFuncEvals", 1);
  flag = KINGetNumLinConvFails(kmem, &ncfl);
  check_flag(&flag, "KINGetNumLinConvFails", 1);
  flag = KINGetNumJtimesEvals(kmem, &njvevals);
  check_flag(&flag, "KINGetNumJtimesEvals", 1);
  flag = KINGetNumPrecEvals(kmem, &npe);
  check_flag(&flag, "KINGetNumPrecEvals", 1);
  flag = KINGetNumPrecSolves(kmem, &nps);
  check_flag(&flag, "KINGetNumPrecSolves", 1);

  /* Main solver workspace size */

  flag = KINGetWorkSpace(kmem, &lenrw, &leniw);
  check_flag(&flag, "KINGetWorkSpace", 1);

  /* Linear solver workspace size */

  flag = KINGetLinWorkSpace(kmem, &lenrwLS, &leniwLS);
  check_flag(&flag, "KINGetLinWorkSpace", 1);

  printf("\nFinal Statistics.. \n\n");
  printf("nni = %6ld  nli   = %6ld  ncfl = %6ld\n", nni, nli, ncfl);
  printf("nfe = %6ld  nfeLS = %6ld  njt  = %6ld\n", nfe, nfeLS, njvevals);
  printf("npe = %6ld  nps   = %6ld\n", npe, nps);
  printf("\n");
  printf("lenrw   = %6ld  leniw   = %6ld\n", lenrw, leniw);
  printf("lenrwLS = %6ld  leniwLS = %6ld\n", lenrwLS, leniwLS);

}

/*
 * Check function return value...
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  return(0);
}
