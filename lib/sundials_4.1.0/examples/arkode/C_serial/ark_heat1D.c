/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a simple 1D heat equation,
 *    u_t = k*u_xx + f
 * for t in [0, 10], x in [0, 1], with initial conditions
 *    u(0,x) =  0
 * Dirichlet boundary conditions, i.e.
 *    u_t(t,0) = u_t(t,1) = 0,
 * and a point-source heating term,
 *    f = 1 for x=0.5.
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 *
 * This program solves the problem with either an ERK or DIRK
 * method.  For the DIRK method, we use a Newton iteration with
 * the SUNLinSol_PCG linear solver, and a user-supplied Jacobian-vector
 * product routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>    /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., macros */
#include <sunlinsol/sunlinsol_pcg.h>  /* access to PCG SUNLinearSolver        */
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
} *UserData;

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
               N_Vector fy, void *user_data, N_Vector tmp);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Main Program */
int main() {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);   /* initial time */
  realtype Tf = RCONST(1.0);   /* final time */
  int Nt = 10;                 /* total number of output times */
  realtype rtol = 1.e-6;       /* relative tolerance */
  realtype atol = 1.e-10;      /* absolute tolerance */
  UserData udata = NULL;
  realtype *data;
  sunindextype N = 201;        /* spatial mesh size */
  realtype k = 0.5;            /* heat conductivity */
  sunindextype i;

  /* general problem variables */
  int flag;                    /* reusable error-checking flag */
  N_Vector y = NULL;           /* empty vector for storing solution */
  SUNLinearSolver LS = NULL;   /* empty linear solver object */
  void *arkode_mem = NULL;     /* empty ARKode memory structure */
  FILE *FID, *UFID;
  realtype t, dTout, tout;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nlcf, nni, ncfn, netf;

  /* allocate and fill udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->N = N;
  udata->k = k;
  udata->dx = RCONST(1.0)/(1.0*N-1.0);     /* mesh spacing */

  /* Initial problem output */
  printf("\n1D Heat PDE test problem:\n");
  printf("  N = %li\n", (long int) udata->N);
  printf("  diffusion coefficient:  k = %"GSYM"\n", udata->k);

  /* Initialize data structures */
  y = N_VNew_Serial(N);            /* Create serial vector for solution */
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
  N_VConst(0.0, y);                /* Set initial conditions */

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y);
  if (check_flag((void *) arkode_mem, "ARKStepCreate", 0)) return 1;

  /* Set routines */
  flag = ARKStepSetUserData(arkode_mem, (void *) udata);   /* Pass udata to user functions */
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
  flag = ARKStepSetMaxNumSteps(arkode_mem, 10000);         /* Increase max num steps  */
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;
  flag = ARKStepSetPredictorMethod(arkode_mem, 1);         /* Specify maximum-order predictor */
  if (check_flag(&flag, "ARKStepSetPredictorMethod", 1)) return 1;
  flag = ARKStepSStolerances(arkode_mem, rtol, atol);      /* Specify tolerances */
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  /* Initialize PCG solver -- no preconditioning, with up to N iterations  */
  LS = SUNLinSol_PCG(y, 0, N);
  if (check_flag((void *)LS, "SUNLinSol_PCG", 0)) return 1;

  /* Linear solver interface -- set user-supplied J*v routine (no 'jtsetup' required) */
  flag = ARKStepSetLinearSolver(arkode_mem, LS, NULL);       /* Attach linear solver to ARKStep */
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;
  flag = ARKStepSetJacTimes(arkode_mem, NULL, Jac);     /* Set the Jacobian routine */
  if (check_flag(&flag, "ARKStepSetJacTimes", 1)) return 1;

  /* Specify linearly implicit RHS, with non-time-dependent Jacobian */
  flag = ARKStepSetLinear(arkode_mem, 0);
  if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;

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

  /* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  dTout = (Tf-T0)/Nt;
  tout = T0+dTout;
  printf("        t      ||u||_rms\n");
  printf("   -------------------------\n");
  printf("  %10.6"FSYM"  %10.6"FSYM"\n", t, SUNRsqrt(N_VDotProd(y,y)/N));
  for (iout=0; iout<Nt; iout++) {

    flag = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);         /* call integrator */
    if (check_flag(&flag, "ARKStepEvolve", 1)) break;
    printf("  %10.6"FSYM"  %10.6"FSYM"\n", t, SUNRsqrt(N_VDotProd(y,y)/N));   /* print solution stats */
    if (flag >= 0) {                                            /* successful solve: update output time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                                    /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }

    /* output results to disk */
    for (i=0; i<N; i++)  fprintf(UFID," %.16"ESYM"", data[i]);
    fprintf(UFID,"\n");
  }
  printf("   -------------------------\n");
  fclose(UFID);

  /* Print some final statistics */
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKStepGetNumSteps", 1);
  flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKStepGetNumStepAttempts", 1);
  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKStepGetNumRhsEvals", 1);
  flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1);
  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKStepGetNumErrTestFails", 1);
  flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1);
  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1);
  flag = ARKStepGetNumLinIters(arkode_mem, &nli);
  check_flag(&flag, "ARKStepGetNumLinIters", 1);
  flag = ARKStepGetNumJtimesEvals(arkode_mem, &nJv);
  check_flag(&flag, "ARKStepGetNumJtimesEvals", 1);
  flag = ARKStepGetNumLinConvFails(arkode_mem, &nlcf);
  check_flag(&flag, "ARKStepGetNumLinConvFails", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total linear iterations = %li\n", nli);
  printf("   Total number of Jacobian-vector products = %li\n", nJv);
  printf("   Total number of linear solver convergence failures = %li\n", nlcf);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy(y);               /* Free vectors */
  free(udata);                 /* Free user data */
  ARKStepFree(&arkode_mem);    /* Free integrator memory */
  SUNLinSolFree(LS);           /* Free linear solver */
  return 0;
}

/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;    /* access problem data */
  sunindextype N  = udata->N;                   /* set variable shortcuts */
  realtype k  = udata->k;
  realtype dx = udata->dx;
  realtype *Y=NULL, *Ydot=NULL;
  realtype c1, c2;
  sunindextype i, isource;

  Y = N_VGetArrayPointer(y);      /* access data arrays */
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return 1;
  Ydot = N_VGetArrayPointer(ydot);
  if (check_flag((void *) Ydot, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                      /* Initialize ydot to zero */

  /* iterate over domain, computing all equations */
  c1 = k/dx/dx;
  c2 = -RCONST(2.0)*k/dx/dx;
  isource = N/2;
  Ydot[0] = 0.0;                 /* left boundary condition */
  for (i=1; i<N-1; i++)
    Ydot[i] = c1*Y[i-1] + c2*Y[i] + c1*Y[i+1];
  Ydot[N-1] = 0.0;               /* right boundary condition */
  Ydot[isource] += 0.01/dx;      /* source term */

  return 0;                      /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
               N_Vector fy, void *user_data, N_Vector tmp)
{
  UserData udata = (UserData) user_data;     /* variable shortcuts */
  sunindextype N = udata->N;
  realtype k  = udata->k;
  realtype dx = udata->dx;
  realtype *V=NULL, *JV=NULL;
  realtype c1, c2;
  sunindextype i;

  V = N_VGetArrayPointer(v);       /* access data arrays */
  if (check_flag((void *) V, "N_VGetArrayPointer", 0)) return 1;
  JV = N_VGetArrayPointer(Jv);
  if (check_flag((void *) JV, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, Jv);                         /* initialize Jv product to zero */

  /* iterate over domain, computing all Jacobian-vector products */
  c1 = k/dx/dx;
  c2 = -RCONST(2.0)*k/dx/dx;
  JV[0] = 0.0;
  for (i=1; i<N-1; i++)
    JV[i] = c1*V[i-1] + c2*V[i] + c1*V[i+1];
  JV[N-1] = 0.0;

  return 0;                                  /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


/*---- end of file ----*/
