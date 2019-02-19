/* ------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ------------------------------------------------------------------
 * Based an example program by Rujeko Chinomona @ SMU.
 * ------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ------------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a simple 1D reaction-diffusion
 * equation,
 *
 *   y_t = k * y_xx + y^2 * (1-y)
 *
 * for t in [0, 3], x in [0, L] with boundary conditions,
 *
 *   y_x(0,t) = y_x(L,t) = 0
 *
 * and initial condition,
 *
 *   y(x,0) = (1 + exp(lambda*(x-1))^(-1),
 *
 * with parameter k = 1e-4/ep, lambda = 0.5*sqrt(2*ep*1e4),
 * ep = 1e-2, and L = 5.
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 * ----------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode_mristep.h>    /* prototypes for MRIStep fcts., consts */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h>  /* defs. of realtype, sunindextype, etc */
#include <sundials/sundials_math.h>   /* def. of SUNRsqrt, etc.               */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* user data structure */
typedef struct {
  sunindextype N;  /* number of intervals   */
  realtype dx;     /* mesh spacing          */
  realtype k;      /* diffusion coefficient */
  realtype lam;
} *UserData;

/* User-supplied Functions Called by the Solver */
static int fs(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int ff(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Private function to set initial condition */
static int SetInitialCondition(N_Vector y, UserData udata);

/* Private function to check function return values */
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Main Program */
int main() {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);     /* initial time */
  realtype Tf = RCONST(3.0);     /* final time */
  realtype dTout = RCONST(0.1);  /* time between outputs */
  int Nt = ceil(Tf/dTout);       /* number of output times */
  realtype hs = RCONST(0.001);   /* slow step size */
  realtype hf = RCONST(0.00002); /* fast step size */
  UserData udata = NULL;         /* user data */

  realtype *data;                /* array for solution output */
  realtype L = RCONST(5.0);      /* domain length */
  sunindextype N = 1001;         /* number of mesh points */
  realtype ep = RCONST(1e-2);
  sunindextype i;

  /* general problem variables */
  int retval;                    /* reusable error-checking flag */
  N_Vector y = NULL;             /* empty vector for storing solution */
  void *arkode_mem = NULL;       /* empty ARKode memory structure */
  FILE *FID, *UFID;
  realtype t, tout;
  int iout;
  long int nsts, nstf, nfs, nff;

  /*
   * Initialization
   */

  /* allocate and fill user data structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->N   = N;
  udata->dx  = L / (RCONST(1.0)*N - RCONST(1.0));
  udata->k   = RCONST(1e-4)/ep;
  udata->lam = RCONST(0.5)*SUNRsqrt(RCONST(2.0) * ep * RCONST(1e4));

  /* Initial problem output */
  printf("\n1D reaction-diffusion PDE test problem:\n");
  printf("  N = %li\n", (long int) udata->N);
  printf("  diffusion coefficient:  k = %"GSYM"\n", udata->k);

  /* Create and initialize serial vector for the solution */
  y = N_VNew_Serial(N);
  if (check_retval((void *) y, "N_VNew_Serial", 0)) return 1;

  retval = SetInitialCondition(y, udata);
  if (check_retval(&retval, "SetInitialCondition", 1)) return 1;

  /* Call MRIStepCreate to initialize the MRI timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = MRIStepCreate(fs, ff, T0, y);
  if (check_retval((void *) arkode_mem, "MRIStepCreate", 0)) return 1;

  /*
   * Set integrator options
   */

  /* Pass udata to user functions */
  retval = MRIStepSetUserData(arkode_mem, (void *) udata);
  if (check_retval(&retval, "MRIStepSetUserData", 1)) return 1;

  /* Specify slow and fast step sizes */
  retval = MRIStepSetFixedStep(arkode_mem, hs, hf);
  if (check_retval(&retval, "MRIStepSetFixedStep", 1)) return 1;

  /* Increase max num steps  */
  retval = MRIStepSetMaxNumSteps(arkode_mem, 10000);
  if (check_retval(&retval, "MRIStepSetMaxNumSteps", 1)) return 1;

  /*
   * Integrate ODE
   */

  /* output mesh to disk */
  FID=fopen("heat_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16"ESYM"\n", udata->dx*i);
  fclose(FID);

  /* Open output stream for results, access data array */
  UFID=fopen("heat1D.txt","w");
  data = N_VGetArrayPointer(y);

  /* output initial condition to disk */
  for (i=0; i<N; i++)  fprintf(UFID," %.16"ESYM"", data[i]);
  fprintf(UFID,"\n");

  /* Main time-stepping loop: calls MRIStepEvolve to perform the integration, then
     prints results. Stops when the final time has been reached */
  t = T0;
  dTout = (Tf-T0)/Nt;
  tout = T0+dTout;
  printf("        t      ||u||_rms\n");
  printf("   -------------------------\n");
  printf("  %10.6"FSYM"  %10.6"FSYM"\n", t, SUNRsqrt(N_VDotProd(y,y)/N));
  for (iout=0; iout<Nt; iout++) {

    /* call integrator */
    retval = MRIStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "MRIStepEvolve", 1)) break;

    /* print solution stats and output results to disk */
    printf("  %10.6"FSYM"  %10.6"FSYM"\n", t, SUNRsqrt(N_VDotProd(y,y)/N));
    for (i=0; i<N; i++)  fprintf(UFID," %.16"ESYM"", data[i]);
    fprintf(UFID,"\n");

    /* successful solve: update output time */
    tout += dTout;
    tout = (tout > Tf) ? Tf : tout;
  }
  printf("   -------------------------\n");
  fclose(UFID);

  /* Print some final statistics */
  retval = MRIStepGetNumSteps(arkode_mem, &nsts, &nstf);
  check_retval(&retval, "ARKStepGetNumSteps", 1);
  retval = MRIStepGetNumRhsEvals(arkode_mem, &nfs, &nff);
  check_retval(&retval, "ARKStepGetNumRhsEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Steps: nsts = %li, nstf = %li\n", nsts, nstf);
  printf("   Total RHS evals:  Fs = %li,  Ff = %li\n", nfs, nff);

  /* Clean up and return */
  N_VDestroy(y);            /* Free y vector */
  MRIStepFree(&arkode_mem); /* Free integrator memory */
  free(udata);              /* Free user data */

  return 0;
}

/* ------------------------------
 * Functions called by the solver
 * ------------------------------*/

/* ff routine to compute the fast portion of the ODE RHS. */
static int ff(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;    /* access problem data */
  sunindextype N = udata->N;                /* set variable shortcuts */
  realtype *Y=NULL, *Ydot=NULL;
  sunindextype i;

  /* access state array data */
  Y = N_VGetArrayPointer(y);
  if (check_retval((void *) Y, "N_VGetArrayPointer", 0)) return 1;

  /* access RHS array data */
  Ydot = N_VGetArrayPointer(ydot);
  if (check_retval((void *) Ydot, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over domain, computing reaction term */
  for (i = 0; i < N; i++)
    Ydot[i] = Y[i] * Y[i] * (RCONST(1.0) - Y[i]);

  /* Return with success */
  return 0;
}


/* fs routine to compute the slow portion of the ODE RHS. */
static int fs(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;    /* access problem data */
  sunindextype N = udata->N;                /* set variable shortcuts */
  realtype k  = udata->k;
  realtype dx = udata->dx;
  realtype *Y=NULL, *Ydot=NULL;
  realtype c1, c2;
  sunindextype i;

  /* access state array data */
  Y = N_VGetArrayPointer(y);
  if (check_retval((void *) Y, "N_VGetArrayPointer", 0)) return 1;

  /* access RHS array data */
  Ydot = N_VGetArrayPointer(ydot);
  if (check_retval((void *) Ydot, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over domain, computing diffusion term */
  c1 = k/dx/dx;
  c2 = RCONST(2.0)*k/dx/dx;

  /* left boundary condition */
  Ydot[0] = c2*(Y[1] - Y[0]);

  /* interior points */
  for (i=1; i<N-1; i++)
    Ydot[i] = c1*Y[i-1] - c2*Y[i] + c1*Y[i+1];

  /* right boundary condition */
  Ydot[N-1] = c2*(Y[N-2] - Y[N-1]);

  /* Return with success */
  return 0;
}

/* -----------------------------------------
 * Private function to set initial condition
 * -----------------------------------------*/

static int SetInitialCondition(N_Vector y, UserData user_data)
{
  UserData udata = (UserData) user_data;    /* access problem data */
  sunindextype N = udata->N;                /* set variable shortcuts */
  realtype lam = udata->lam;
  realtype dx = udata->dx;
  realtype *Y=NULL;
  sunindextype i;

  /* access state array data */
  Y = N_VGetArrayPointer(y);
  if (check_retval((void *) Y, "N_VGetArrayPointer", 0)) return -1;

  /* set initial condition */
  for (i = 0; i < N; i++)
    Y[i] = RCONST(1.0)/(1 + SUNRexp(lam*(i*dx-RCONST(1.0))));

  /* Return with success */
  return 0;
}

/* ------------------------------
 * Private helper functions
 * ------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval < 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


/*---- end of file ----*/
