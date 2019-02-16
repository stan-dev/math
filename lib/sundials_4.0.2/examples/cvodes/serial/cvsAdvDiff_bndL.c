/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
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
 * Example problem:
 *
 * The following is a simple example problem with a banded Jacobian,
 * with the program for its solution by CVODE.
 * The problem is the semi-discrete form of the advection-diffusion
 * equation in 2-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
 * interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
 * are posed, and the initial condition is
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX*MY.
 * This program solves the problem with the BDF method, Newton
 * iteration with the LAPACKBAND linear solver, and a user-supplied
 * Jacobian routine.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .1, .2, ..., 1.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvodes/cvodes.h>                   /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>          /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_band.h>        /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_lapackband.h>  /* access to band SUNLinearSolver       */
#include <sundials/sundials_types.h>         /* definition of type realtype          */
#include <sundials/sundials_math.h>          /* definition of ABS and EXP            */

/* Problem Constants */

#define XMAX  RCONST(2.0)    /* domain boundaries         */
#define YMAX  RCONST(1.0)
#define MX    10             /* mesh dimensions           */
#define MY    5
#define NEQ   MX*MY          /* number of equations       */
#define ATOL  RCONST(1.0e-5) /* scalar absolute tolerance */
#define T0    RCONST(0.0)    /* initial time              */
#define T1    RCONST(0.1)    /* first output time         */
#define DTOUT RCONST(0.1)    /* output time increment     */
#define NOUT  10             /* number of output times    */

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)
#define FIVE RCONST(5.0)

/* User-defined vector access macro IJth */

/* IJth is defined in order to isolate the translation from the
   mathematical 2-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. 
   IJth(vdata,i,j) references the element in the vdata array for
   u at mesh point (i,j), where 1 <= i <= MX, 1 <= j <= MY.
   The vdata array is obtained via the call vdata = N_VGetArrayPointer(v),
   where v is an N_Vector. 
   The variables are ordered by the y index j, then by the x index i. */

#define IJth(vdata,i,j) (vdata[(j-1) + (i-1)*MY])

/* Type : UserData (contains grid constants) */

typedef struct {
  realtype dx, dy, hdcoef, hacoef, vdcoef;
} *UserData;

/* Private Helper Functions */

static void SetIC(N_Vector u, UserData data);
static void PrintHeader(realtype reltol, realtype abstol, realtype umax);
static void PrintOutput(realtype t, realtype umax, long int nst);
static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);
static int Jac(realtype t, N_Vector u, N_Vector fu, 
               SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(void)
{
  realtype dx, dy, reltol, abstol, t, tout, umax;
  N_Vector u;
  UserData data;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, retval;
  long int nst;

  u = NULL;
  data = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Create a serial vector */

  u = N_VNew_Serial(NEQ);  /* Allocate u vector */
  if(check_retval((void*)u, "N_VNew_Serial", 0)) return(1);

  reltol = ZERO;  /* Set the tolerances */
  abstol = ATOL;

  data = (UserData) malloc(sizeof *data);  /* Allocate data memory */
  if(check_retval((void *)data, "malloc", 2)) return(1);
  dx = data->dx = XMAX/(MX+1);  /* Set grid coefficients in data */
  dy = data->dy = YMAX/(MY+1);
  data->hdcoef = ONE/(dx*dx);
  data->hacoef = HALF/(TWO*dx);
  data->vdcoef = ONE/(dy*dy);

  SetIC(u, data);  /* Initialize u vector */

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerance */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* Set the pointer to user-defined data */
  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Create banded SUNMatrix for use in linear solves -- since this will be factored, 
     set the storage bandwidth to be the sum of upper and lower bandwidths */
  A = SUNBandMatrix(NEQ, MY, MY);
  if(check_retval((void *)A, "SUNBandMatrix", 0)) return(1);

  /* Create SUNLinSol_LapackBand solver object for use by CVode */
  LS = SUNLinSol_LapackBand(u, A);
  if(check_retval((void *)LS, "SUNLinSol_LapackBand", 0)) return(1);
  
  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* Set the user-supplied Jacobian routine Jac */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

  /* In loop over output points: call CVode, print results, test for errors */

  umax = N_VMaxNorm(u);
  PrintHeader(reltol, abstol, umax);
  for(iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_retval(&retval, "CVode", 1)) break;
    umax = N_VMaxNorm(u);
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_retval(&retval, "CVodeGetNumSteps", 1);
    PrintOutput(t, umax, nst);
  }

  PrintFinalStats(cvode_mem);  /* Print some final statistics   */

  N_VDestroy(u);          /* Free the u vector */
  CVodeFree(&cvode_mem);  /* Free the integrator memory */
  free(data);             /* Free the user data */

  return(0);
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* f routine. Compute f(t,u). */

static int f(realtype t, N_Vector u,N_Vector udot, void *user_data)
{
  realtype uij, udn, uup, ult, urt, hordc, horac, verdc, hdiff, hadv, vdiff;
  realtype *udata, *dudata;
  int i, j;
  UserData data;

  udata  = N_VGetArrayPointer(u);
  dudata = N_VGetArrayPointer(udot);

  /* Extract needed constants from data */

  data = (UserData) user_data;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

  /* Loop over all grid points. */

  for (j=1; j <= MY; j++) {

    for (i=1; i <= MX; i++) {

      /* Extract u at x_i, y_j and four neighboring points */

      uij = IJth(udata, i, j);
      udn = (j == 1)  ? ZERO : IJth(udata, i, j-1);
      uup = (j == MY) ? ZERO : IJth(udata, i, j+1);
      ult = (i == 1)  ? ZERO : IJth(udata, i-1, j);
      urt = (i == MX) ? ZERO : IJth(udata, i+1, j);

      /* Set diffusion and advection terms and load into udot */

      hdiff = hordc*(ult - TWO*uij + urt);
      hadv = horac*(urt - ult);
      vdiff = verdc*(uup - TWO*uij + udn);
      IJth(dudata, i, j) = hdiff + hadv + vdiff;
    }
  }

  return(0);
}

/* Jacobian routine. Compute J(t,u). */

static int Jac(realtype t, N_Vector u, N_Vector fu, 
               SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int i, j, k;
  realtype *kthCol, hordc, horac, verdc;
  UserData data;
  
  /*
   * The components of f = udot that depend on u(i,j) are
   * f(i,j), f(i-1,j), f(i+1,j), f(i,j-1), f(i,j+1), with
   *   df(i,j)/du(i,j) = -2 (1/dx^2 + 1/dy^2)
   *   df(i-1,j)/du(i,j) = 1/dx^2 + .25/dx  (if i > 1)
   *   df(i+1,j)/du(i,j) = 1/dx^2 - .25/dx  (if i < MX)
   *   df(i,j-1)/du(i,j) = 1/dy^2           (if j > 1)
   *   df(i,j+1)/du(i,j) = 1/dy^2           (if j < MY)
   */

  data = (UserData) user_data;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

  /* set non-zero Jacobian entries */
  for (j=1; j <= MY; j++) {
    for (i=1; i <= MX; i++) {
      k = j-1 + (i-1)*MY;
      kthCol = SUNBandMatrix_Column(J,k);

      /* set the kth column of J */

      SM_COLUMN_ELEMENT_B(kthCol,k,k) = -TWO*(verdc+hordc);
      if (i != 1)  SM_COLUMN_ELEMENT_B(kthCol,k-MY,k) = hordc + horac;
      if (i != MX) SM_COLUMN_ELEMENT_B(kthCol,k+MY,k) = hordc - horac;
      if (j != 1)  SM_COLUMN_ELEMENT_B(kthCol,k-1,k)  = verdc;
      if (j != MY) SM_COLUMN_ELEMENT_B(kthCol,k+1,k)  = verdc;
    }
  }

  return(0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

/* Set initial conditions in u vector */

static void SetIC(N_Vector u, UserData data)
{
  int i, j;
  realtype x, y, dx, dy;
  realtype *udata;

  /* Extract needed constants from data */

  dx = data->dx;
  dy = data->dy;

  /* Set pointer to data array in vector u. */

  udata = N_VGetArrayPointer(u);

  /* Load initial profile into u vector */
  
  for (j=1; j <= MY; j++) {
    y = j*dy;
    for (i=1; i <= MX; i++) {
      x = i*dx;
      IJth(udata,i,j) = x*(XMAX - x)*y*(YMAX - y)*SUNRexp(FIVE*x*y);
    }
  }  
}

/* Print first lines of output (problem description) */

static void PrintHeader(realtype reltol, realtype abstol, realtype umax)
{
  printf("\n2-D Advection-Diffusion Equation\n");
  printf("Mesh dimensions = %d X %d\n", MX, MY);
  printf("Total system size = %d\n", NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters: reltol = %Lg   abstol = %Lg\n\n", reltol, abstol);
  printf("At t = %Lg      max.norm(u) =%14.6Le \n", T0, umax);
#else
  printf("Tolerance parameters: reltol = %g    abstol = %g\n\n", reltol, abstol);
  printf("At t = %g       max.norm(u) =%14.6e  \n", T0, umax);
#endif

  return;
}

/* Print current value */

static void PrintOutput(realtype t, realtype umax, long int nst)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %4.2Lf   max.norm(u) =%14.6Le   nst = %4ld\n", t, umax, nst);
#else
  printf("At t = %4.2f   max.norm(u) =%14.6e    nst = %4ld\n", t, umax, nst);
#endif

  return;
}

/* Get and print some final statistics */

static void PrintFinalStats(void *cvode_mem)
{
  int retval;
  long int nst, nfe, nsetups, netf, nni, ncfn, nje, nfeLS;

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
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %ld\n \n",
	 nni, ncfn, netf);

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
