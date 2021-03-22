/*---------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ LLNL
 *---------------------------------------------------------------
 * Based on the serial code ark_heat1D_adapt.c developed
 * by Daniel R. Reynolds and parallelized with OpenMP 4.5
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
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
 * and a heating term of the form
 *    f = 2*exp(-200*(x-0.25)*(x-0.25))
 *        - exp(-400*(x-0.7)*(x-0.7))
 *        + exp(-500*(x-0.4)*(x-0.4))
 *        - 2*exp(-600*(x-0.55)*(x-0.55));
 *
 * The spatial derivatives are computed using a three-point
 * centered stencil (second order for a uniform mesh).  The data
 * is initially uniformly distributed over N points in the interval
 * [0, 1], but as the simulation proceeds the mesh is adapted.
 *
 * This program solves the problem with a DIRK method, solved with
 * a Newton iteration and SUNLinSol_PCG linear solver, with a
 * user-supplied Jacobian-vector product routine.
 *---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>     /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_openmpdev.h> /* OpenMPDEV N_Vector types, fcts., macros */
#include <sunlinsol/sunlinsol_pcg.h>   /* access to PCG SUNLinearSolver        */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype, etc */
#include <sundials/sundials_math.h>    /* def. of SUNRsqrt, etc.               */

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

/* constants */
#define ZERO        RCONST(0.0)
#define PT25        RCONST(0.25)
#define PT4         RCONST(0.4)
#define PT5         RCONST(0.5)
#define PT55        RCONST(0.55)
#define PT7         RCONST(0.7)
#define ONE         RCONST(1.0)
#define TWO         RCONST(2.0)
#define TWOHUNDRED  RCONST(200.0)
#define FOURHUNDRED RCONST(400.0)
#define FIVEHUNDRED RCONST(500.0)
#define SIXHUNDRED  RCONST(600.0)

/* user data structure */
typedef struct {
  sunindextype N;       /* current number of intervals */
  realtype *x_host;     /* current mesh on host */
  realtype *x_dev;      /* current mesh on device */
  realtype k;           /* diffusion coefficient */
  realtype refine_tol;  /* adaptivity tolerance */
} *UserData;

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
               N_Vector fy, void *user_data, N_Vector tmp);

/* Private function to check function return values */
realtype * adapt_mesh(N_Vector y, sunindextype *Nnew, UserData udata);
static int project(sunindextype Nold, realtype *xold, N_Vector yold,
                   sunindextype Nnew, realtype *xnew, N_Vector ynew);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/* Main Program */
int main() {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);         /* initial time */
  realtype Tf = RCONST(1.0);         /* final time */
  realtype rtol = RCONST(1.e-3);     /* relative tolerance */
  realtype atol = RCONST(1.e-10);    /* absolute tolerance */
  realtype hscale = RCONST(1.0);     /* time step change factor on resizes */
  UserData udata = NULL;
  realtype *data;
  sunindextype N = 21;               /* initial spatial mesh size */
  realtype refine = RCONST(3.0e-3);  /* adaptivity refinement tolerance */
  realtype k = RCONST(0.5);          /* heat conductivity */
  sunindextype i;
  long int nni, nni_tot=0, nli, nli_tot=0;
  int iout=0;

  /* general problem variables */
  int flag;                    /* reusable error-checking flag */
  N_Vector y  = NULL;          /* empty vector for storing solution */
  N_Vector y2 = NULL;          /* empty vector for storing solution */
  N_Vector yt = NULL;          /* empty vector for swapping */
  SUNLinearSolver LS = NULL;   /* empty linear solver object */
  void *arkode_mem = NULL;     /* empty ARKode memory structure */
  FILE *XFID, *UFID;
  realtype t, olddt, newdt;
  realtype *xnew_host = NULL;
  realtype *xnew_dev = NULL;
  sunindextype Nnew;
  int dev, host;

  /* get host and offloading device */
  dev  = omp_get_default_device();
  host = omp_get_initial_device();

  /* allocate and fill initial udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->N = N;
  udata->k = k;
  udata->refine_tol = refine;
  udata->x_host = malloc(N * sizeof(realtype));
  for (i=0; i<N; i++)  udata->x_host[i] = ONE*i/(N-1);
  udata->x_dev = omp_target_alloc(N * sizeof(realtype), dev);
  omp_target_memcpy(udata->x_dev, udata->x_host, N * sizeof(realtype), 0, 0, dev, host);

  /* Initial problem output */
  printf("\n1D adaptive Heat PDE test problem:\n");
  printf("  diffusion coefficient:  k = %"GSYM"\n", udata->k);
  printf("  initial N = %li\n", (long int) udata->N);

  /* Initialize data structures */
  y = N_VNew_OpenMPDEV(N);   /* Create initial OpenMPDEV vector for solution */
  if (check_flag((void *) y, "N_VNew_OpenMPDEV", 0)) return 1;
  N_VConst(ZERO, y);  /* Set initial conditions */

  /* output mesh to disk */
  XFID=fopen("heat_mesh.txt","w");

  /* output initial mesh to disk */
  for (i=0; i<udata->N; i++)  fprintf(XFID," %.16"ESYM, udata->x_host[i]);
  fprintf(XFID,"\n");

  /* Open output stream for results, access data array */
  UFID=fopen("heat1D.txt","w");

  /* output initial condition to disk */
  N_VCopyFromDevice_OpenMPDEV(y);
  data = N_VGetHostArrayPointer_OpenMPDEV(y);
  for (i=0; i<udata->N; i++)  fprintf(UFID," %.16"ESYM, data[i]);
  fprintf(UFID,"\n");

  /* Initialize the ARK timestepper */
  arkode_mem = ARKStepCreate(NULL, f, T0, y);
  if (check_flag((void *) arkode_mem, "ARKStepCreate", 0)) return 1;

  /* Set routines */
  flag = ARKStepSetUserData(arkode_mem, (void *) udata);   /* Pass udata to user functions */
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
  flag = ARKStepSetMaxNumSteps(arkode_mem, 10000);         /* Increase max num steps  */
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1)) return 1;
  flag = ARKStepSStolerances(arkode_mem, rtol, atol);      /* Specify tolerances */
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;
  flag = ARKStepSetAdaptivityMethod(arkode_mem, 2, 1, 0, NULL);  /* Set adaptivity method */
  if (check_flag(&flag, "ARKStepSetAdaptivityMethod", 1)) return 1;
  flag = ARKStepSetPredictorMethod(arkode_mem, 0);     /* Set predictor method */
  if (check_flag(&flag, "ARKStepSetPredictorMethod", 1)) return 1;

  /* Specify linearly implicit RHS, with time-dependent Jacobian */
  flag = ARKStepSetLinear(arkode_mem, 1);
  if (check_flag(&flag, "ARKStepSetLinear", 1)) return 1;

  /* Initialize PCG solver -- no preconditioning, with up to N iterations  */
  LS = SUNLinSol_PCG(y, 0, (int) N);
  if (check_flag((void *)LS, "SUNLinSol_PCG", 0)) return 1;

  /* Linear solver interface -- set user-supplied J*v routine (no 'jtsetup' required) */
  flag = ARKStepSetLinearSolver(arkode_mem, LS, NULL);        /* Attach linear solver to ARKStep */
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;
  flag = ARKStepSetJacTimes(arkode_mem, NULL, Jac);     /* Set the Jacobian routine */
  if (check_flag(&flag, "ARKStepSetJacTimes", 1)) return 1;

  /* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  olddt = ZERO;
  newdt = ZERO;
  printf("  iout          dt_old                 dt_new               ||u||_rms       N   NNI  NLI\n");
  printf(" ----------------------------------------------------------------------------------------\n");
  printf(" %4i  %19.15"ESYM"  %19.15"ESYM"  %19.15"ESYM"  %li   %2i  %3i\n",
         iout, olddt, newdt, SUNRsqrt(N_VDotProd(y,y)/udata->N),
         (long int) udata->N, 0, 0);
  while (t < Tf) {

    /* "set" routines */
    flag = ARKStepSetStopTime(arkode_mem, Tf);
    if (check_flag(&flag, "ARKStepSetStopTime", 1)) return 1;
    flag = ARKStepSetInitStep(arkode_mem, newdt);
    if (check_flag(&flag, "ARKStepSetInitStep", 1)) return 1;

    /* call integrator */
    flag = ARKStepEvolve(arkode_mem, Tf, y, &t, ARK_ONE_STEP);
    if (check_flag(&flag, "ARKStepEvolve", 1)) return 1;

    /* "get" routines */
    flag = ARKStepGetLastStep(arkode_mem, &olddt);
    if (check_flag(&flag, "ARKStepGetLastStep", 1)) return 1;
    flag = ARKStepGetCurrentStep(arkode_mem, &newdt);
    if (check_flag(&flag, "ARKStepGetCurrentStep", 1)) return 1;
    flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
    if (check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1)) return 1;
    flag = ARKStepGetNumLinIters(arkode_mem, &nli);
    if (check_flag(&flag, "ARKStepGetNumLinIters", 1)) return 1;

    /* print current solution stats */
    iout++;
    printf(" %4i  %19.15"ESYM"  %19.15"ESYM"  %19.15"ESYM"  %li   %2li  %3li\n",
           iout, olddt, newdt, SUNRsqrt(N_VDotProd(y,y)/udata->N),
           (long int) udata->N, nni, nli);
    nni_tot += nni;
    nli_tot += nli;

    /* output results and current mesh to disk */
    N_VCopyFromDevice_OpenMPDEV(y);
    data = N_VGetHostArrayPointer_OpenMPDEV(y);
    for (i=0; i<udata->N; i++)  fprintf(UFID," %.16"ESYM, data[i]);
    fprintf(UFID,"\n");
    for (i=0; i<udata->N; i++)  fprintf(XFID," %.16"ESYM, udata->x_host[i]);
    fprintf(XFID,"\n");

    /* adapt the spatial mesh */
    xnew_host = adapt_mesh(y, &Nnew, udata);
    if (check_flag(xnew_host, "ark_adapt", 0)) return 1;

    /* create N_Vector of new length */
    y2 = N_VNew_OpenMPDEV(Nnew);
    if (check_flag((void *) y2, "N_VNew_OpenMPDEV", 0)) return 1;

    /* copy new mesh from host array to device array */
    xnew_dev = omp_target_alloc(Nnew * sizeof(realtype), dev);
    omp_target_memcpy(xnew_dev, xnew_host, Nnew*sizeof(realtype), 0, 0, dev, host);

    /* project solution onto new mesh */
    flag = project(udata->N, udata->x_dev, y, Nnew, xnew_dev, y2);
    if (check_flag(&flag, "project", 1)) return 1;

    /* delete old vector, old mesh */
    N_VDestroy(y);
    free(udata->x_host);
    omp_target_free(udata->x_dev, dev);

    /* swap x and xnew so that new mesh is stored in udata structure */
    udata->x_host = xnew_host;
    xnew_host = NULL;
    udata->N = Nnew;   /* store size of new mesh */
    udata->x_dev = xnew_dev;
    xnew_dev = NULL;

    /* swap y and y2 so that y holds new solution */
    yt = y;
    y  = y2;
    y2 = yt;

    /* call ARKStepResize to notify integrator of change in mesh */
    flag = ARKStepResize(arkode_mem, y, hscale, t, NULL, NULL);
    if (check_flag(&flag, "ARKStepResize", 1)) return 1;

    /* destroy and re-allocate linear solver memory; reattach to ARKStep interface */
    SUNLinSolFree(LS);
    LS = SUNLinSol_PCG(y, 0, (int) N);
    if (check_flag((void *)LS, "SUNLinSol_PCG", 0)) return 1;
    flag = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
    if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;
    flag = ARKStepSetJacTimes(arkode_mem, NULL, Jac);
    if (check_flag(&flag, "ARKStepSetJacTimes", 1)) return 1;

  }
  printf(" ----------------------------------------------------------------------------------------\n");

  /* print some final statistics */
  printf(" Final solver statistics:\n");
  printf("   Total number of time steps = %i\n", iout);
  printf("   Total nonlinear iterations = %li\n", nni_tot);
  printf("   Total linear iterations    = %li\n\n", nli_tot);

  /* Clean up and return with successful completion */
  fclose(UFID);
  fclose(XFID);
  N_VDestroy(y);               /* Free vectors */
  free(udata->x_host);         /* Free user data */
  omp_target_free(udata->x_dev, dev);
  free(udata);
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
  sunindextype N = udata->N;                /* set variable shortcuts */
  realtype k  = udata->k;
  realtype *x = udata->x_dev;
  realtype *Y=NULL, *Ydot=NULL;
  realtype dxL, dxR;
  sunindextype i;
  int dev;

  dev = omp_get_default_device();

  /* access data arrays */
  Y = N_VGetDeviceArrayPointer_OpenMPDEV(y);
  if (check_flag((void *) Y, "N_VGetDeviceArrayPointer", 0)) return 1;
  Ydot = N_VGetDeviceArrayPointer_OpenMPDEV(ydot);
  if (check_flag((void *) Ydot, "N_VGetDeviceArrayPointer", 0)) return 1;

  /* Initialize ydot to zero - also handles boundary conditions */
  N_VConst(ZERO, ydot);

  /* iterate over domain interior, computing all equations */
#pragma omp target map(to:N) is_device_ptr(x, Ydot, Y) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i=1; i<N-1; i++) {        /* interior */
    dxL = x[i]-x[i-1];
    dxR = x[i+1]-x[i];
    Ydot[i] = Y[i-1]*k*TWO/(dxL*(dxL+dxR))
      - Y[i]*k*TWO/(dxL*dxR)
      + Y[i+1]*k*TWO/(dxR*(dxL+dxR))
      + TWO*SUNRexp(-TWOHUNDRED*(x[i]-PT25)*(x[i]-PT25)) /* source term */
      - SUNRexp(-FOURHUNDRED*(x[i]-PT7)*(x[i]-PT7))
      + SUNRexp(-FIVEHUNDRED*(x[i]-PT4)*(x[i]-PT4))
      - TWO*SUNRexp(-SIXHUNDRED*(x[i]-PT55)*(x[i]-PT55));
  }

  return 0;                      /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y,
               N_Vector fy, void *user_data, N_Vector tmp)
{
  UserData udata = (UserData) user_data;     /* variable shortcuts */
  sunindextype N = udata->N;
  realtype k  = udata->k;
  realtype *x = udata->x_dev;
  realtype *V=NULL, *JV=NULL;
  realtype dxL, dxR;
  sunindextype i;
  int dev;

  dev  = omp_get_default_device();

  /* access data arrays */
  V = N_VGetDeviceArrayPointer_OpenMPDEV(v);
  if (check_flag((void *) V, "N_VGetDeviceArrayPointer", 0)) return 1;
  JV = N_VGetDeviceArrayPointer_OpenMPDEV(Jv);
  if (check_flag((void *) JV, "N_VGetDeviceArrayPointer", 0)) return 1;

  /* initialize Jv product to zero - also handles boundary conditions */
  N_VConst(ZERO, Jv);

  /* iterate over domain, computing all Jacobian-vector products */
#pragma omp target map(to:N) is_device_ptr(x, JV, V) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i=1; i<N-1; i++) {
    dxL = x[i]-x[i-1];
    dxR = x[i+1]-x[i];
    JV[i] = V[i-1]*k*TWO/(dxL*(dxL+dxR))
          - V[i]*k*TWO/(dxL*dxR)
          + V[i+1]*k*TWO/(dxR*(dxL+dxR));
  }

  return 0;                                  /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Adapts the current mesh, using a simple adaptivity strategy of
   refining when an approximation of the scaled second-derivative is
   too large.  We only do this in one sweep, so no attempt is made to
   ensure the resulting mesh meets these same criteria after adaptivity:
      y [input] -- the current solution vector
      Nnew [output] -- the size of the new mesh
      udata [input] -- the current system information
   The return for this function is a pointer to the new mesh. */
realtype* adapt_mesh(N_Vector y, sunindextype *Nnew, UserData udata)
{
  sunindextype i, j;
  int *marks=NULL;
  realtype ydd, *xold=NULL, *Y=NULL, *xnew=NULL;
  sunindextype num_refine, N_new;

  /* Access current solution and mesh arrays */
  xold = udata->x_host;
  Y = N_VGetHostArrayPointer_OpenMPDEV(y);  /* assumes copy to host already done */
  if (check_flag((void *) Y, "N_VGetHostArrayPointer_OpenMPDEV", 0)) return NULL;

  /* create marking array */
  marks = calloc(udata->N-1, sizeof(int));

  /* perform marking:
      0 -> leave alone
      1 -> refine */
  for (i=1; i<udata->N-1; i++) {

    /* approximate scaled second-derivative */
    ydd = Y[i-1] - TWO*Y[i] + Y[i+1];

    /* check for refinement */
    if (fabs(ydd) > udata->refine_tol) {
      marks[i-1] = 1;
      marks[i] = 1;
    }

  }

  /* allocate new mesh */
  num_refine = 0;
  for (i=0; i<udata->N-1; i++)
    if (marks[i] == 1)   num_refine++;
  N_new = udata->N + num_refine;
  *Nnew = N_new;            /* Store new array length */
  xnew = malloc((N_new) * sizeof(realtype));


  /* fill new mesh */
  xnew[0] = xold[0];    /* store endpoints */
  xnew[N_new-1] = xold[udata->N-1];
  j=1;
  /* iterate over old intervals */
  for (i=0; i<udata->N-1; i++) {
    /* if mark is 0, reuse old interval */
    if (marks[i] == 0) {
      xnew[j++] = xold[i+1];
      continue;
    }

    /* if mark is 1, refine old interval */
    if (marks[i] == 1) {
      xnew[j++] = PT5*(xold[i]+xold[i+1]);
      xnew[j++] = xold[i+1];
      continue;
    }
  }

  /* verify that new mesh is legal */
  for (i=0; i<N_new-1; i++) {
    if (xnew[i+1] <= xnew[i]) {
      fprintf(stderr,"adapt_mesh error: illegal mesh created\n");
      free(xnew);
      return NULL;
    }
  }

  free(marks);              /* Delete marking array */
  return xnew;              /* Return with success */
}


/* Projects one vector onto another:
      Nold [input] -- the size of the old mesh
      xold [input] -- the old mesh
      yold [input] -- the vector defined over the old mesh
      Nnew [input] -- the size of the new mesh
      xnew [input] -- the new mesh
      ynew [output] -- the vector defined over the new mesh
                       (allocated prior to calling project) */
static int project(sunindextype Nold, realtype *xold, N_Vector yold,
                   sunindextype Nnew, realtype *xnew, N_Vector ynew)
{
  sunindextype iv, i, j;
  realtype *Yold=NULL, *Ynew=NULL;

  int dev = omp_get_default_device();

  /* Access data arrays */
  Yold = N_VGetDeviceArrayPointer_OpenMPDEV(yold);    /* access data arrays */
  if (check_flag((void *) Yold, "N_VGetDeviceArrayPointer_OpenMPDEV", 0)) return 1;
  Ynew = N_VGetDeviceArrayPointer_OpenMPDEV(ynew);
  if (check_flag((void *) Ynew, "N_VGetDeviceArrayPointer_OpenMPDEV", 0)) return 1;

  /* loop over new mesh, finding corresponding interval within old mesh,
     and perform piecewise linear interpolation from yold to ynew */
  iv=0;
#pragma omp target map(to:iv) is_device_ptr(Yold,Ynew,xnew,xold) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  {
    for (i=0; i<Nnew; i++) {

      /* find old interval, start with previous value since sorted */
      for (j=iv; j<Nold-1; j++) {
        if (xnew[i] >= xold[j] && xnew[i] <= xold[j+1]) {
          iv = j;
          break;
        }
        iv = Nold-1;     /* just in case it wasn't found above */
      }

      /* perform interpolation */
      Ynew[i] = Yold[iv]*(xnew[i]-xold[iv+1])/(xold[iv]-xold[iv+1])
              + Yold[iv+1]*(xnew[i]-xold[iv])/(xold[iv+1]-xold[iv]);
    }
  }

  return 0;            /* Return with success */
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
