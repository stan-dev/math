/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, serial, banded.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the band solver and IDACalcIC.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform M x M
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = M^2. Here M = 10.
 *
 * The system is solved with IDA using the banded linear system
 * solver, half-bandwidths equal to M, and default
 * difference-quotient Jacobian. For purposes of illustration,
 * IDACalcIC is called to compute correct values at the boundary,
 * given incorrect values as input initial guesses. The constraints
 * u >= 0 are posed for all components. Output is taken at
 * t = 0, .01, .02, .04, ..., 10.24. (Output at t = 0 is for
 * IDACalcIC cost statistics only.)
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>                   /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver       */
#include <sundials/sundials_types.h>   /* definition of type realtype          */

/* Problem Constants */

#define NOUT  11
#define MGRID 10
#define NEQ   MGRID*MGRID
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define BVAL  RCONST(0.1)

/* Type: UserData */

typedef struct {
  sunindextype mm;
  realtype dx;
  realtype coeff;
} *UserData;

/* Prototypes of functions called by IDA */

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector resval, void *user_data);

/* Prototypes of private functions */

static void PrintHeader(realtype rtol, realtype atol);
static void PrintOutput(void *mem, realtype t, N_Vector u);
static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up,
                             N_Vector id, N_Vector res);

static int check_retval(void *returnvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(void)
{
  void *mem;
  UserData data;
  N_Vector uu, up, constraints, id, res;
  int retval, iout;
  long int netf, ncfn;
  sunindextype mu, ml;
  realtype rtol, atol, t0, t1, tout, tret;
  SUNMatrix A;
  SUNLinearSolver LS;

  mem = NULL;
  data = NULL;
  uu = up = constraints = id = res = NULL;
  A = NULL;
  LS = NULL;

  /* Create vectors uu, up, res, constraints, id. */
  uu = N_VNew_Serial(NEQ);
  if(check_retval((void *)uu, "N_VNew_Serial", 0)) return(1);
  up = N_VNew_Serial(NEQ);
  if(check_retval((void *)up, "N_VNew_Serial", 0)) return(1);
  res = N_VNew_Serial(NEQ);
  if(check_retval((void *)res, "N_VNew_Serial", 0)) return(1);
  constraints = N_VNew_Serial(NEQ);
  if(check_retval((void *)constraints, "N_VNew_Serial", 0)) return(1);
  id = N_VNew_Serial(NEQ);
  if(check_retval((void *)id, "N_VNew_Serial", 0)) return(1);

  /* Create and load problem data block. */
  data = (UserData) malloc(sizeof *data);
  if(check_retval((void *)data, "malloc", 2)) return(1);
  data->mm = MGRID;
  data->dx = ONE/(MGRID - ONE);
  data->coeff = ONE/( (data->dx) * (data->dx) );

  /* Initialize uu, up, id. */
  SetInitialProfile(data, uu, up, id, res);

  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst(ONE, constraints);

  /* Set remaining input parameters. */
  t0   = ZERO;
  t1   = RCONST(0.01);
  rtol = ZERO;
  atol = RCONST(1.0e-3);

  /* Call IDACreate and IDAMalloc to initialize solution */
  mem = IDACreate();
  if(check_retval((void *)mem, "IDACreate", 0)) return(1);

  retval = IDASetUserData(mem, data);
  if(check_retval(&retval, "IDASetUserData", 1)) return(1);

  /* Set which components are algebraic or differential */
  retval = IDASetId(mem, id);
  if(check_retval(&retval, "IDASetId", 1)) return(1);

  retval = IDASetConstraints(mem, constraints);
  if(check_retval(&retval, "IDASetConstraints", 1)) return(1);
  N_VDestroy(constraints);

  retval = IDAInit(mem, heatres, t0, uu, up);
  if(check_retval(&retval, "IDAInit", 1)) return(1);

  retval = IDASStolerances(mem, rtol, atol);
  if(check_retval(&retval, "IDASStolerances", 1)) return(1);

  /* Create banded SUNMatrix for use in linear solves */
  mu = MGRID; ml = MGRID;
  A = SUNBandMatrix(NEQ, mu, ml);
  if(check_retval((void *)A, "SUNBandMatrix", 0)) return(1);

  /* Create banded SUNLinearSolver object */
  LS = SUNLinSol_Band(uu, A);
  if(check_retval((void *)LS, "SUNLinSol_Band", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  /* Call IDACalcIC to correct the initial values. */

  retval = IDACalcIC(mem, IDA_YA_YDP_INIT, t1);
  if(check_retval(&retval, "IDACalcIC", 1)) return(1);

  /* Print output heading. */
  PrintHeader(rtol, atol);

  PrintOutput(mem, t0, uu);


  /* Loop over output times, call IDASolve, and print results. */

  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {

    retval = IDASolve(mem, tout, &tret, uu, up, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1)) return(1);

    PrintOutput(mem, tret, uu);

  }

  /* Print remaining counters and free memory. */
  retval = IDAGetNumErrTestFails(mem, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  printf("\n netf = %ld,   ncfn = %ld \n", netf, ncfn);

  IDAFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(uu);
  N_VDestroy(up);
  N_VDestroy(id);
  N_VDestroy(res);
  free(data);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * heatres: heat equation system residual function
 * This uses 5-point central differencing on the interior points, and
 * includes algebraic equations for the boundary values.
 * So for each interior point, the residual component has the form
 *    res_i = u'_i - (central difference)_i
 * while for each boundary point, it is res_i = u_i.
 */

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector resval,
            void *user_data)
{
  sunindextype mm, i, j, offset, loc;
  realtype *uv, *upv, *resv, coeff;
  UserData data;

  uv = N_VGetArrayPointer(uu); upv = N_VGetArrayPointer(up); resv = N_VGetArrayPointer(resval);

  data = (UserData)user_data;
  mm = data->mm;
  coeff = data->coeff;

  /* Initialize resval to uu, to take care of boundary equations. */
  N_VScale(ONE, uu, resval);

  /* Loop over interior points; set res = up - (central difference). */
  for (j = 1; j < mm-1; j++) {
    offset = mm*j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      resv[loc] = upv[loc] - coeff *
	  (uv[loc-1] + uv[loc+1] + uv[loc-mm] + uv[loc+mm] - RCONST(4.0)*uv[loc]);
    }
  }

  return(0);

}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * SetInitialProfile: routine to initialize u, up, and id vectors.
 */

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up,
                             N_Vector id, N_Vector res)
{
  realtype xfact, yfact, *udata, *updata, *iddata;
  sunindextype mm, mm1, i, j, offset, loc;

  mm = data->mm;
  mm1 = mm - 1;

  udata = N_VGetArrayPointer(uu);
  updata = N_VGetArrayPointer(up);
  iddata = N_VGetArrayPointer(id);

  /* Initialize id to 1's. */
  N_VConst(ONE, id);

  /* Initialize uu on all grid points. */
  for (j = 0; j < mm; j++) {
    yfact = data->dx * j;
    offset = mm*j;
    for (i = 0;i < mm; i++) {
      xfact = data->dx * i;
      loc = offset + i;
      udata[loc] = RCONST(16.0) * xfact * (ONE - xfact) * yfact * (ONE - yfact);
    }
  }

  /* Initialize up vector to 0. */
  N_VConst(ZERO, up);

  /* heatres sets res to negative of ODE RHS values at interior points. */
  heatres(ZERO, uu, up, res, data);

  /* Copy -res into up to get correct interior initial up values. */
  N_VScale(-ONE, res, up);

  /* Finally, set values of u, up, and id at boundary points. */
  for (j = 0; j < mm; j++) {
    offset = mm*j;
    for (i = 0;i < mm; i++) {
      loc = offset + i;
      if (j == 0 || j == mm1 || i == 0 || i == mm1 ) {
        udata[loc] = BVAL; updata[loc] = ZERO; iddata[loc] = ZERO; }
    }
  }

  return(0);

}

/*
 * Print first lines of output (problem description)
 */

static void PrintHeader(realtype rtol, realtype atol)
{
  printf("\nidaHeat2D_bnd: Heat equation, serial example problem for IDA\n");
  printf("          Discretized heat equation on 2D unit square.\n");
  printf("          Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("          Mesh dimensions: %d x %d", MGRID, MGRID);
  printf("        Total system size: %d\n\n", NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Constraints set to force all solution components >= 0. \n");
  printf("Linear solver: BAND, banded direct solver \n");
  printf("       difference quotient Jacobian, half-bandwidths = %d \n",MGRID);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("IDACalcIC called with input boundary values = %Lg \n",BVAL);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("IDACalcIC called with input boundary values = %g \n",BVAL);
#else
  printf("IDACalcIC called with input boundary values = %g \n",BVAL);
#endif
  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time       umax     k  nst  nni  nje   nre   nreLS    h      \n" );
  printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . \n");
}

/*
 * Print Output
 */

static void PrintOutput(void *mem, realtype t, N_Vector uu)
{
  int retval;
  realtype umax, hused;
  long int nst, nni, nje, nre, nreLS;
  int kused;

  umax = N_VMaxNorm(uu);

  retval = IDAGetLastOrder(mem, &kused);
  check_retval(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  check_retval(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_retval(&retval, "IDAGetNumResEvals", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_retval(&retval, "IDAGetLastStep", 1);
  retval = IDAGetNumJacEvals(mem, &nje);
  check_retval(&retval, "IDAGetNumJacEvals", 1);
  retval = IDAGetNumLinResEvals(mem, &nreLS);
  check_retval(&retval, "IDAGetNumLinResEvals", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2Le \n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused);
#else
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused);
#endif

}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if retval < 0 */
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1);
    }
  } else if (opt == 2 && returnvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  }

  return(0);
}
