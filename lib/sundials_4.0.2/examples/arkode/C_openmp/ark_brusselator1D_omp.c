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
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is n PDE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    u_t = du*u_xx + a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx + w*u - v*u^2
 *    w_t = dw*w_xx + (b-w)/ep - w*u
 * for t in [0, 80], x in [0, 1], with initial conditions
 *    u(0,x) =  a  + 0.1*sin(pi*x)
 *    v(0,x) = b/a + 0.1*sin(pi*x)
 *    w(0,x) =  b  + 0.1*sin(pi*x),
 * and with stationary boundary conditions, i.e.
 *    u_t(t,0) = u_t(t,1) = 0,
 *    v_t(t,0) = v_t(t,1) = 0,
 *    w_t(t,0) = w_t(t,1) = 0.
 * Note: these can also be implemented as Dirichlet boundary
 * conditions with values identical to the initial conditions.
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 *
 * This program solves the problem with the DIRK method, using a
 * Newton iteration with the band linear solver, and a
 * user-supplied Jacobian routine.  This example uses the OpenMP
 * vector kernel, and employs OpenMP threading within the
 * right-hand side and Jacobian construction functions.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>    /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_openmp.h>   /* access to OpenMP N_Vector */
#include <sunmatrix/sunmatrix_band.h> /* access to band SUNMatrix */
#include <sunlinsol/sunlinsol_band.h> /* access to band SUNLinearSolver */
#include <sundials/sundials_types.h>  /* def. of type 'realtype' */
#include <sundials/sundials_math.h>   /* def. of SUNRsqrt, etc. */
#ifdef _OPENMP
#include <omp.h>                      /* OpenMP functions */
#endif

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* accessor macros between (x,v) location and 1D NVector array */
#define IDX(x,v) (3*(x)+v)

/* user data structure */
typedef struct {
  sunindextype N;  /* number of intervals      */
  int nthreads;    /* number of OpenMP threads */
  realtype dx;     /* mesh spacing             */
  realtype a;      /* constant forcing on u    */
  realtype b;      /* steady-state value of w  */
  realtype du;     /* diffusion coeff for u    */
  realtype dv;     /* diffusion coeff for v    */
  realtype dw;     /* diffusion coeff for w    */
  realtype ep;     /* stiffness parameter      */
} *UserData;


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy,
               SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private helper functions  */
static int LaplaceMatrix(realtype c, SUNMatrix Jac, UserData udata);
static int ReactionJac(realtype c, N_Vector y, SUNMatrix Jac, UserData udata);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Main Program */
int main(int argc, char *argv[])
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);    /* initial time */
  realtype Tf = RCONST(10.0);   /* final time */
  int Nt = 100;                 /* total number of output times */
  int Nvar = 3;                 /* number of solution fields */
  UserData udata = NULL;
  realtype *data;
  sunindextype N = 201;         /* spatial mesh size */
  realtype a = 0.6;             /* problem parameters */
  realtype b = 2.0;
  realtype du = 0.025;
  realtype dv = 0.025;
  realtype dw = 0.025;
  realtype ep = 1.0e-5;         /* stiffness parameter */
  realtype reltol = 1.0e-6;     /* tolerances */
  realtype abstol = 1.0e-10;
  sunindextype NEQ, i;

  /* general problem variables */
  int flag;                     /* reusable error-checking flag */
  N_Vector y = NULL;            /* empty vector for storing solution */
  N_Vector umask = NULL;        /* empty mask vectors for viewing solution components */
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  SUNMatrix A = NULL;           /* empty matrix for linear solver */
  SUNLinearSolver LS = NULL;    /* empty linear solver structure */
  void *arkode_mem = NULL;      /* empty ARKode memory structure */
  realtype pi, t, dTout, tout, u, v, w;
  FILE *FID, *UFID, *VFID, *WFID;
  int iout, num_threads;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  if (check_flag((void *) udata, "malloc", 2)) return 1;

  /* set the number of threads to use */
  num_threads = 1;                       /* default value */
#ifdef _OPENMP
  num_threads = omp_get_max_threads();   /* overwrite with OMP_NUM_THREADS environment variable */
#endif
  if (argc > 1)                          /* overwrite with command line value, if supplied */
    num_threads = strtol(argv[1], NULL, 0);

  /* store the inputs in the UserData structure */
  udata->N  = N;
  udata->a  = a;
  udata->b  = b;
  udata->du = du;
  udata->dv = dv;
  udata->dw = dw;
  udata->ep = ep;
  udata->nthreads = num_threads;

  /* set total allocated vector length */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D Brusselator PDE test problem:\n");
  printf("    N = %li,  NEQ = %li\n", (long int) udata->N, (long int) NEQ);
  printf("    num_threads = %i\n", num_threads);
  printf("    problem parameters:  a = %"GSYM",  b = %"GSYM",  ep = %"GSYM"\n",
         udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %"GSYM",  dv = %"GSYM",  dw = %"GSYM"\n",
         udata->du, udata->dv, udata->dw);
  printf("    reltol = %.1"ESYM",  abstol = %.1"ESYM"\n\n", reltol, abstol);

  /* Initialize vector data structures */
  y = N_VNew_OpenMP(NEQ, num_threads);      /* Create vector for solution */
  if (check_flag((void *)y, "N_VNew_OpenMP", 0)) return 1;
  udata->dx = RCONST(1.0)/(N-1);            /* set spatial mesh spacing */
  data = N_VGetArrayPointer(y);             /* Access data array for new NVector y */
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  umask = N_VNew_OpenMP(NEQ, num_threads);  /* Create vector masks */
  if (check_flag((void *)umask, "N_VNew_OpenMP", 0)) return 1;
  vmask = N_VNew_OpenMP(NEQ, num_threads);
  if (check_flag((void *)vmask, "N_VNew_OpenMP", 0)) return 1;
  wmask = N_VNew_OpenMP(NEQ, num_threads);
  if (check_flag((void *)wmask, "N_VNew_OpenMP", 0)) return 1;

  /* Set initial conditions into y */
  pi = RCONST(4.0)*atan(RCONST(1.0));
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*i*udata->dx);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*i*udata->dx);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*i*udata->dx);  /* w */
  }

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  data = N_VGetArrayPointer(umask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,0)] = RCONST(1.0);

  N_VConst(0.0, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,1)] = RCONST(1.0);

  N_VConst(0.0, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,2)] = RCONST(1.0);

  /* Initialize matrix and linear solver data structures */
  A = SUNBandMatrix(NEQ, 4, 4);
  if (check_flag((void *)A, "SUNBandMatrix", 0)) return 1;
  LS = SUNLinSol_Band(y, A);
  if (check_flag((void *)LS, "SUNLinSol_Band", 0)) return 1;

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y);
  if (check_flag((void *)arkode_mem, "ARKStepCreate", 0)) return 1;

  /* Set routines */
  flag = ARKStepSetUserData(arkode_mem, (void *) udata);     /* Pass udata to user functions */
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
  flag = ARKStepSStolerances(arkode_mem, reltol, abstol);    /* Specify tolerances */
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  /* Linear solver specification */
  flag = ARKStepSetLinearSolver(arkode_mem, LS, A);          /* Attach matrix and linear solver */
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;
  flag = ARKStepSetJacFn(arkode_mem, Jac);                   /* Set the Jacobian routine */
  if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;

  /* output spatial mesh to disk */
  FID=fopen("bruss_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16"ESYM"\n", udata->dx*i);
  fclose(FID);

  /* Open output stream for results, access data arrays */
  UFID=fopen("bruss_u.txt","w");
  VFID=fopen("bruss_v.txt","w");
  WFID=fopen("bruss_w.txt","w");

  /* output initial condition to disk */
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  fprintf(UFID," %.16"ESYM, data[IDX(i,0)]);
  for (i=0; i<N; i++)  fprintf(VFID," %.16"ESYM, data[IDX(i,1)]);
  for (i=0; i<N; i++)  fprintf(WFID," %.16"ESYM, data[IDX(i,2)]);
  fprintf(UFID,"\n");
  fprintf(VFID,"\n");
  fprintf(WFID,"\n");

  /* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  dTout = (Tf-T0)/Nt;
  tout = T0+dTout;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms\n");
  printf("   ----------------------------------------------\n");
  for (iout=0; iout<Nt; iout++) {

    flag = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);    /* call integrator */
    if (check_flag(&flag, "ARKStepEvolve", 1)) break;
    u = N_VWL2Norm(y,umask);                               /* access/print solution statistics */
    u = SUNRsqrt(u*u/N);
    v = N_VWL2Norm(y,vmask);
    v = SUNRsqrt(v*v/N);
    w = N_VWL2Norm(y,wmask);
    w = SUNRsqrt(w*w/N);
    printf("  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"\n", t, u, v, w);
    if (flag >= 0) {                                       /* successful solve: update output time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                               /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }

    /* output results to disk */
    for (i=0; i<N; i++)  fprintf(UFID," %.16"ESYM, data[IDX(i,0)]);
    for (i=0; i<N; i++)  fprintf(VFID," %.16"ESYM, data[IDX(i,1)]);
    for (i=0; i<N; i++)  fprintf(WFID," %.16"ESYM, data[IDX(i,2)]);
    fprintf(UFID,"\n");
    fprintf(VFID,"\n");
    fprintf(WFID,"\n");
  }
  printf("   ----------------------------------------------\n");
  fclose(UFID);
  fclose(VFID);
  fclose(WFID);

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
  flag = ARKStepGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKStepGetNumJacEvals", 1);
  flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);

  /* Clean up and return with successful completion */
  free(udata);                  /* Free user data */
  ARKStepFree(&arkode_mem);     /* Free integrator memory */
  SUNLinSolFree(LS);            /* Free linear solver */
  SUNMatDestroy(A);             /* Free matrix */
  N_VDestroy_OpenMP(y);         /* Free vectors */
  N_VDestroy_OpenMP(umask);
  N_VDestroy_OpenMP(vmask);
  N_VDestroy_OpenMP(wmask);
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;      /* access problem data */
  sunindextype N = udata->N;                  /* set variable shortcuts */
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  realtype dx = udata->dx;
  realtype *Ydata=NULL, *dYdata=NULL;
  realtype uconst, vconst, wconst, u, ul, ur, v, vl, vr, w, wl, wr;
  sunindextype i;

  /* clear out ydot (to be careful) */
  N_VConst(0.0, ydot);

  Ydata = N_VGetArrayPointer(y);     /* access data arrays */
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  dYdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *)dYdata, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                        /* initialize ydot to zero */

  /* iterate over domain, computing all equations */
  uconst = du/dx/dx;
  vconst = dv/dx/dx;
  wconst = dw/dx/dx;
#pragma omp parallel for default(shared) private(i,u,ul,ur,v,vl,vr,w,wl,wr) schedule(static) num_threads(udata->nthreads)
  for (i=1; i<N-1; i++) {

    /* set shortcuts */
    u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* u_t = du*u_xx + a - (w+1)*u + v*u^2 */
    dYdata[IDX(i,0)] = (ul - RCONST(2.0)*u + ur)*uconst + a - (w+RCONST(1.0))*u + v*u*u;

    /* v_t = dv*v_xx + w*u - v*u^2 */
    dYdata[IDX(i,1)] = (vl - RCONST(2.0)*v + vr)*vconst + w*u - v*u*u;

    /* w_t = dw*w_xx + (b-w)/ep - w*u */
    dYdata[IDX(i,2)] = (wl - RCONST(2.0)*w + wr)*wconst + (b-w)/ep - w*u;

  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  return 0;
}


/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(realtype t, N_Vector y, N_Vector fy,
               SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData udata = (UserData) user_data;  /* access problem data */
  SUNMatZero(J);                          /* Initialize Jacobian to zero */

  /* Fill in the Laplace matrix */
  if (LaplaceMatrix(RCONST(1.0), J, udata)) {
    printf("Jacobian calculation error in calling LaplaceMatrix!\n");
    return 1;
  }

  /* Add in the Jacobian of the reaction terms matrix */
  if (ReactionJac(RCONST(1.0), y, J, udata)) {
    printf("Jacobian calculation error in calling ReactionJac!\n");
    return 1;
  }

  return 0;
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Routine to compute the stiffness matrix from (L*y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int LaplaceMatrix(realtype c, SUNMatrix Jac, UserData udata)
{
  sunindextype N = udata->N;            /* set shortcuts */
  realtype dx = udata->dx;
  sunindextype i;
  realtype uconst = c*udata->du/dx/dx;
  realtype vconst = c*udata->dv/dx/dx;
  realtype wconst = c*udata->dw/dx/dx;

  /* iterate over intervals, filling in Jacobian entries */
#pragma omp parallel for default(shared) private(i) schedule(static) num_threads(udata->nthreads)
  for (i=1; i<N-1; i++) {

    /* Jacobian of (L*y) at this node */
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i-1,0)) += uconst;
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i-1,1)) += vconst;
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i-1,2)) += wconst;
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i,0)) -= RCONST(2.0)*uconst;
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i,1)) -= RCONST(2.0)*vconst;
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i,2)) -= RCONST(2.0)*wconst;
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i+1,0)) += uconst;
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i+1,1)) += vconst;
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i+1,2)) += wconst;
  }

  return 0;
}



/* Routine to compute the Jacobian matrix from R(y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int ReactionJac(realtype c, N_Vector y, SUNMatrix Jac, UserData udata)
{
  sunindextype N = udata->N;                   /* set shortcuts */
  realtype ep = udata->ep;
  sunindextype i;
  realtype u, v, w;
  realtype *Ydata = N_VGetArrayPointer(y);     /* access solution array */
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over nodes, filling in Jacobian entries */
#pragma omp parallel for default(shared) private(i,u,v,w) schedule(static) num_threads(udata->nthreads)
  for (i=1; i<N-1; i++) {

    /* set nodal value shortcuts (shifted index due to start at first interior node) */
    u = Ydata[IDX(i,0)];
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* all vars wrt u */
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i,0)) += c*(RCONST(2.0)*u*v-(w+RCONST(1.0)));
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i,0)) += c*(w - RCONST(2.0)*u*v);
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i,0)) += c*(-w);

    /* all vars wrt v */
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i,1)) += c*(u*u);
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i,1)) += c*(-u*u);

    /* all vars wrt w */
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i,2)) += c*(-u);
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i,2)) += c*(u);
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i,2)) += c*(-RCONST(1.0)/ep - u);

  }

  return 0;
}

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
