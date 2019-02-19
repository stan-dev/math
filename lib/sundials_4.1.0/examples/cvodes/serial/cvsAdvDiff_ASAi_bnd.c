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
 * Adjoint sensitivity example problem:
 *
 * The following is a simple example problem with a banded Jacobian,
 * with the program for its solution by CVODES.
 * The problem is the semi-discrete form of the advection-diffusion
 * equation in 2-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
 * interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
 * are posed, and the initial condition is the following:
 *   u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX*MY.
 * This program solves the problem with the BDF method, Newton
 * iteration with the BAND linear solver, and a user-supplied
 * Jacobian routine.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .1, .2, ..., 1.
 * Run statistics (optional outputs) are printed at the end.
 *
 * Additionally, CVODES integrates backwards in time the
 * the semi-discrete form of the adjoint PDE:
 *   d(lambda)/dt = - d^2(lambda) / dx^2 + 0.5 d(lambda) / dx
 *                  - d^2(lambda) / dy^2 - 1.0
 * with homogeneous Dirichlet boundary conditions and final
 * conditions:
 *   lambda(x,y,t=t_final) = 0.0
 * whose solution at t = 0 represents the sensitivity of
 *   G = int_0^t_final int_x int _y u(t,x,y) dx dy dt
 * with respect to the initial conditions of the original problem.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Header files with a description of contents */

#include <cvodes/cvodes.h>             /* prototypes for CVODES fcts., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver       */
#include <sundials/sundials_types.h>   /* definition of type realtype          */
#include <sundials/sundials_math.h>    /* definition of SUNRabs and SUNRexp    */

/* Problem Constants */

#define XMAX  RCONST(2.0)   /* domain boundaries             */
#define YMAX  RCONST(1.0)
#define MX    40            /* mesh dimensions               */
#define MY    20
#define NEQ   MX*MY         /* number of equations           */
#define ATOL  RCONST(1.e-5)        
#define RTOLB RCONST(1.e-6)        
#define T0    RCONST(0.0)   /* initial time                  */
#define T1    RCONST(0.1)   /* first output time             */
#define DTOUT RCONST(0.1)   /* output time increment         */
#define NOUT  10            /* number of output times        */
#define TOUT  RCONST(1.0)   /* final time                    */
#define NSTEP 50            /* check point saved every NSTEP */

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

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

/* Type : UserData 
   contains grid constants */

typedef struct {
  realtype dx, dy, hdcoef, hacoef, vdcoef;
} *UserData;

/* Prototypes of user-supplied functions */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static int Jac(realtype t, N_Vector u, N_Vector fu, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3); 

static int fB(realtype tB, N_Vector u, N_Vector uB, N_Vector uBdot, void *user_dataB);

static int JacB(realtype tB, N_Vector u, N_Vector uB, N_Vector fuB, SUNMatrix JB,
                void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B); 

/* Prototypes of private functions */

static void SetIC(N_Vector u, UserData data);
static void PrintOutput(N_Vector uB, UserData data);
static int check_retval(void *returnvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;

  void *cvode_mem;
  SUNMatrix A, AB;
  SUNLinearSolver LS, LSB;

  realtype dx, dy, reltol, abstol, t;
  N_Vector u;

  int indexB;

  realtype reltolB, abstolB;
  N_Vector uB;
  
  int retval, ncheck;

  data = NULL;
  cvode_mem = NULL;
  u = uB = NULL;
  LS = LSB = NULL;
  A = AB = NULL;

  /* Allocate and initialize user data memory */

  data = (UserData) malloc(sizeof *data);
  if(check_retval((void *)data, "malloc", 2)) return(1);

  dx = data->dx = XMAX/(MX+1);
  dy = data->dy = YMAX/(MY+1);
  data->hdcoef = ONE/(dx*dx);
  data->hacoef = RCONST(1.5)/(TWO*dx);
  data->vdcoef = ONE/(dy*dy);

  /* Set the tolerances for the forward integration */
  reltol = ZERO;
  abstol = ATOL;

  /* Allocate u vector */
  u = N_VNew_Serial(NEQ);
  if(check_retval((void *)u, "N_VNew", 0)) return(1);

  /* Initialize u vector */
  SetIC(u, data);

  /* Create and allocate CVODES memory for forward run */

  printf("\nCreate and allocate CVODES memory for forward runs\n");

  cvode_mem = CVodeCreate(CV_BDF);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);

  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* Create banded SUNMatrix for the forward problem */
  A = SUNBandMatrix(NEQ, MY, MY);
  if(check_retval((void *)A, "SUNBandMatrix", 0)) return(1);

  /* Create banded SUNLinearSolver for the forward problem */
  LS = SUNLinSol_Band(u, A);
  if(check_retval((void *)LS, "SUNLinSol_Band", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* Set the user-supplied Jacobian routine for the forward problem */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

  /* Allocate global memory */

  printf("\nAllocate global memory\n");

  retval = CVodeAdjInit(cvode_mem, NSTEP, CV_HERMITE);
  if(check_retval(&retval, "CVodeAdjInit", 1)) return(1);

  /* Perform forward run */
  printf("\nForward integration\n");
  retval = CVodeF(cvode_mem, TOUT, u, &t, CV_NORMAL, &ncheck);
  if(check_retval(&retval, "CVodeF", 1)) return(1);

  printf("\nncheck = %d\n", ncheck);

  /* Set the tolerances for the backward integration */
  reltolB = RTOLB;
  abstolB = ATOL;

  /* Allocate uB */
  uB = N_VNew_Serial(NEQ);
  if(check_retval((void *)uB, "N_VNew", 0)) return(1);
  /* Initialize uB = 0 */
  N_VConst(ZERO, uB);

  /* Create and allocate CVODES memory for backward run */

  printf("\nCreate and allocate CVODES memory for backward run\n");

  retval = CVodeCreateB(cvode_mem, CV_BDF, &indexB);
  if(check_retval(&retval, "CVodeCreateB", 1)) return(1);

  retval = CVodeSetUserDataB(cvode_mem, indexB, data);
  if(check_retval(&retval, "CVodeSetUserDataB", 1)) return(1);

  retval = CVodeInitB(cvode_mem, indexB, fB, TOUT, uB);
  if(check_retval(&retval, "CVodeInitB", 1)) return(1);

  retval = CVodeSStolerancesB(cvode_mem, indexB, reltolB, abstolB);
  if(check_retval(&retval, "CVodeSStolerancesB", 1)) return(1);
 
  /* Create banded SUNMatrix for the backward problem */
  AB = SUNBandMatrix(NEQ, MY, MY);
  if(check_retval((void *)AB, "SUNBandMatrix", 0)) return(1);

  /* Create banded SUNLinearSolver for the backward problem */
  LSB = SUNLinSol_Band(uB, AB);
  if(check_retval((void *)LSB, "SUNLinSol_Band", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = CVodeSetLinearSolverB(cvode_mem, indexB, LSB, AB);
  if(check_retval(&retval, "CVodeSetLinearSolverB", 1)) return(1);

  /* Set the user-supplied Jacobian routine for the backward problem */
  retval = CVodeSetJacFnB(cvode_mem, indexB, JacB);
  if(check_retval(&retval, "CVodeSetJacFnB", 1)) return(1);

  /* Perform backward integration */
  printf("\nBackward integration\n");
  retval = CVodeB(cvode_mem, T0, CV_NORMAL);
  if(check_retval(&retval, "CVodeB", 1)) return(1);

  retval = CVodeGetB(cvode_mem, indexB, &t, uB);
  if(check_retval(&retval, "CVodeGetB", 1)) return(1);

  PrintOutput(uB, data);

  N_VDestroy(u);   /* Free the u vector                      */
  N_VDestroy(uB);  /* Free the uB vector                     */
  CVodeFree(&cvode_mem);  /* Free the CVODE problem memory          */
  SUNLinSolFree(LS);      /* Free the forward linear solver memory  */
  SUNMatDestroy(A);       /* Free the forward matrix memory         */
  SUNLinSolFree(LSB);     /* Free the backward linear solver memory */
  SUNMatDestroy(AB);      /* Free the backward matrix memory        */

  free(data);             /* Free the user data */

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. right-hand side of forward ODE.
 */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype uij, udn, uup, ult, urt, hordc, horac, verdc, hdiff, hadv, vdiff;
  realtype *udata, *dudata;
  int i, j;
  UserData data;

  udata = N_VGetArrayPointer(u);
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

/*
 * Jac function. Jacobian of forward ODE.
 */

static int Jac(realtype t, N_Vector u, N_Vector fu, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int i, j, k;
  realtype *kthCol, hordc, horac, verdc;
  UserData data;

  /*
    The components of f = udot that depend on u(i,j) are
    f(i,j), f(i-1,j), f(i+1,j), f(i,j-1), f(i,j+1), with
      df(i,j)/du(i,j) = -2 (1/dx^2 + 1/dy^2)
      df(i-1,j)/du(i,j) = 1/dx^2 + .25/dx  (if i > 1)
      df(i+1,j)/du(i,j) = 1/dx^2 - .25/dx  (if i < MX)
      df(i,j-1)/du(i,j) = 1/dy^2           (if j > 1)
      df(i,j+1)/du(i,j) = 1/dy^2           (if j < MY)
  */

  data = (UserData) user_data;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

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
 * fB function. Right-hand side of backward ODE.
 */

static int fB(realtype tB, N_Vector u, N_Vector uB, N_Vector uBdot, 
              void *user_dataB)
{
  UserData data;
  realtype *uBdata, *duBdata;
  realtype hordc, horac, verdc;
  realtype uBij, uBdn, uBup, uBlt, uBrt;
  realtype hdiffB, hadvB, vdiffB;
  int i, j;

  uBdata = N_VGetArrayPointer(uB);
  duBdata = N_VGetArrayPointer(uBdot);

  /* Extract needed constants from data */

  data = (UserData) user_dataB;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

  /* Loop over all grid points. */

  for (j=1; j <= MY; j++) {

    for (i=1; i <= MX; i++) {

      /* Extract u at x_i, y_j and four neighboring points */

      uBij = IJth(uBdata, i, j);
      uBdn = (j == 1)  ? ZERO : IJth(uBdata, i, j-1);
      uBup = (j == MY) ? ZERO : IJth(uBdata, i, j+1);
      uBlt = (i == 1)  ? ZERO : IJth(uBdata, i-1, j);
      uBrt = (i == MX) ? ZERO : IJth(uBdata, i+1, j);

      /* Set diffusion and advection terms and load into udot */

      hdiffB = hordc*(- uBlt + TWO*uBij - uBrt);
      hadvB  = horac*(uBrt - uBlt);
      vdiffB = verdc*(- uBup + TWO*uBij - uBdn);
      IJth(duBdata, i, j) = hdiffB + hadvB + vdiffB - ONE;
    }
  }

  return(0);
}

/*
 * JacB function. Jacobian of backward ODE
 */

static int JacB(realtype tB, N_Vector u, N_Vector uB, N_Vector fuB, SUNMatrix JB,
                void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  int i, j, k;
  realtype *kthCol, hordc, horac, verdc;
  UserData data;

  /* The Jacobian of the adjoint system is: JB = -J^T */

  data = (UserData) user_dataB;
  hordc = data->hdcoef;
  horac = data->hacoef;
  verdc = data->vdcoef;

  for (j=1; j <= MY; j++) {
    for (i=1; i <= MX; i++) {
      k = j-1 + (i-1)*MY;
      kthCol = SUNBandMatrix_Column(JB,k);

      /* set the kth column of J */

      SM_COLUMN_ELEMENT_B(kthCol,k,k) = TWO*(verdc+hordc);
      if (i != 1)  SM_COLUMN_ELEMENT_B(kthCol,k-MY,k) = - hordc + horac;
      if (i != MX) SM_COLUMN_ELEMENT_B(kthCol,k+MY,k) = - hordc - horac;
      if (j != 1)  SM_COLUMN_ELEMENT_B(kthCol,k-1,k)  = - verdc;
      if (j != MY) SM_COLUMN_ELEMENT_B(kthCol,k+1,k)  = - verdc;
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
 * Set initial conditions in u vector 
 */

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
      IJth(udata,i,j) = x*(XMAX - x)*y*(YMAX - y)*SUNRexp(RCONST(5.0)*x*y);
    }
  }  

}

/*
 * Print results after backward integration 
 */

static void PrintOutput(N_Vector uB, UserData data)
{
  realtype *uBdata, uBij, uBmax, x, y, dx, dy;
  int i, j;

  x = y = ZERO;

  dx = data->dx;
  dy = data->dy;

  uBdata = N_VGetArrayPointer(uB);

  uBmax = ZERO;
  for(j=1; j<= MY; j++) {
    for(i=1; i<=MX; i++) {
      uBij = IJth(uBdata, i, j);
      if (SUNRabs(uBij) > uBmax) {
        uBmax = uBij;
        x = i*dx;
        y = j*dy;
      }
    }
  }

  printf("\nMaximum sensitivity\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("  lambda max = %Le\n", uBmax);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("  lambda max = %e\n", uBmax);
#else
  printf("  lambda max = %e\n", uBmax);
#endif
  printf("at\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("  x = %Le\n  y = %Le\n", x, y);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("  x = %e\n  y = %e\n", x, y);
#else
  printf("  x = %e\n  y = %e\n", x, y);
#endif

}

/* 
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns an integer value so check if
 *             retval < 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
 */

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
