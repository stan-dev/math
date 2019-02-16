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
 * The number of spatial points N, the parameters a, b, du, dv,
 * dw and ep, as well as the desired relative and absolute solver
 * tolerances, are provided in the input file
 * input_brusselator1D.txt.
 *
 * This program solves the problem with the DIRK method, using a
 * Newton iteration.  The inner linear systems are solved using
 * the SUNLinSol_KLU linear solver.
 *
 * 100 outputs are printed at equal intervals, and run statistics
 * are printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>       /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_serial.h>      /* serial N_Vector types, fcts., macros */
#include <sunmatrix/sunmatrix_sparse.h>  /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_klu.h>     /* access to KLU SUNLinearSolver        */
#include <sundials/sundials_types.h>     /* defs. of realtype, sunindextype, etc */
#include <sundials/sundials_math.h>      /* def. of SUNRsqrt, etc.               */

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

/* constants */
#define ONE (RCONST(1.0))
#define TWO (RCONST(2.0))

/* user data structure */
typedef struct {
  sunindextype N;  /* number of intervals     */
  realtype dx;     /* mesh spacing            */
  realtype a;      /* constant forcing on u   */
  realtype b;      /* steady-state value of w */
  realtype du;     /* diffusion coeff for u   */
  realtype dv;     /* diffusion coeff for v   */
  realtype dw;     /* diffusion coeff for w   */
  realtype ep;     /* stiffness parameter     */
  SUNMatrix R;     /* temporary storage       */
} *UserData;


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int LaplaceMatrix(SUNMatrix Jac, UserData udata);
static int ReactionJac(N_Vector y, SUNMatrix Jac, UserData udata);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const char *funcname, int opt);


/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);    /* initial time */
  realtype Tf = RCONST(10.0);   /* final time */
  int Nt = 10;                  /* total number of output times */
  int Nvar = 3;
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
  sunindextype NEQ, i, NNZ;

  /* general problem variables */
  int flag;                     /* reusable error-checking flag */
  N_Vector y = NULL;
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  SUNMatrix A = NULL;           /* empty matrix object for solver */
  SUNLinearSolver LS = NULL;    /* empty linear solver object */
  void *arkode_mem = NULL;
  realtype pi;
  FILE *FID, *UFID, *VFID, *WFID;
  realtype t = T0;
  realtype dTout = (Tf-T0)/Nt;
  realtype tout = T0+dTout;
  realtype u, v, w;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nni, ncfn, netf;

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  if (check_flag((void *) udata, "malloc", 2)) return 1;

  /* store the inputs in the UserData structure */
  udata->N  = N;
  udata->a  = a;
  udata->b  = b;
  udata->du = du;
  udata->dv = dv;
  udata->dw = dw;
  udata->ep = ep;
  udata->R  = NULL;

  /* set total allocated vector length */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D Brusselator PDE test problem (KLU solver):\n");
  printf("    N = %li,  NEQ = %li\n", (long int) udata->N, (long int) NEQ);
  printf("    problem parameters:  a = %"GSYM",  b = %"GSYM",  ep = %"GSYM"\n",
         udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %"GSYM",  dv = %"GSYM",  dw = %"GSYM"\n",
         udata->du, udata->dv, udata->dw);
  printf("    reltol = %.1"ESYM",  abstol = %.1"ESYM"\n\n", reltol, abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ);           /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  udata->dx = RCONST(1.0)/(N-1);    /* set spatial mesh spacing */
  data = N_VGetArrayPointer(y);     /* Access data array for new NVector y */
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  umask = N_VNew_Serial(NEQ);       /* Create serial vector masks */
  if (check_flag((void *)umask, "N_VNew_Serial", 0)) return 1;
  vmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)vmask, "N_VNew_Serial", 0)) return 1;
  wmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)wmask, "N_VNew_Serial", 0)) return 1;

  /* Set initial conditions into y */
  pi = RCONST(4.0)*atan(ONE);
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*i*udata->dx);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*i*udata->dx);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*i*udata->dx);  /* w */
  }

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  data = N_VGetArrayPointer(umask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,0)] = ONE;

  N_VConst(0.0, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,1)] = ONE;

  N_VConst(0.0, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,2)] = ONE;


  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y);
  if (check_flag((void *) arkode_mem, "ARKStepCreate", 0)) return 1;

  /* Set routines */
  flag = ARKStepSetUserData(arkode_mem, (void *) udata);     /* Pass udata to user functions */
  if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
  flag = ARKStepSStolerances(arkode_mem, reltol, abstol);    /* Specify tolerances */
  if (check_flag(&flag, "ARKStepSStolerances", 1)) return 1;

  /* Initialize sparse matrix data structure and KLU solver */
  NNZ = 5*NEQ;
  A = SUNSparseMatrix(NEQ, NEQ, NNZ, CSC_MAT);
  if (check_flag((void *)A, "SUNSparseMatrix", 0)) return 1;
  LS = SUNLinSol_KLU(y, A);
  if (check_flag((void *)LS, "SUNLinSol_KLU", 0)) return 1;

  /* Attach the matrix, linear solver, and Jacobian construction routine to ARKStep */
  flag = ARKStepSetLinearSolver(arkode_mem, LS, A);        /* Attach matrix and LS */
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return 1;
  flag = ARKStepSetJacFn(arkode_mem, Jac);                 /* Supply Jac routine */
  if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;

   /* output spatial mesh to disk */
  FID = fopen("bruss_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16"ESYM"\n", udata->dx*i);
  fclose(FID);

  /* Open output stream for results, access data arrays */
  UFID=fopen("bruss_u.txt","w");
  VFID=fopen("bruss_v.txt","w");
  WFID=fopen("bruss_w.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;

  /* output initial condition to disk */
  for (i=0; i<N; i++)  fprintf(UFID," %.16"ESYM"", data[IDX(i,0)]);
  for (i=0; i<N; i++)  fprintf(VFID," %.16"ESYM"", data[IDX(i,1)]);
  for (i=0; i<N; i++)  fprintf(WFID," %.16"ESYM"", data[IDX(i,2)]);
  fprintf(UFID,"\n");
  fprintf(VFID,"\n");
  fprintf(WFID,"\n");

  /* Main time-stepping loop: calls ARKStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t  = T0;
  dTout = Tf/Nt;
  tout = T0+dTout;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms\n");
  printf("   ----------------------------------------------\n");
  for (iout=0; iout<Nt; iout++) {

    flag = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);    /* call integrator */
    u = N_VWL2Norm(y,umask);
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
    for (i=0; i<N; i++)  fprintf(UFID," %.16"ESYM"", data[IDX(i,0)]);
    for (i=0; i<N; i++)  fprintf(VFID," %.16"ESYM"", data[IDX(i,1)]);
    for (i=0; i<N; i++)  fprintf(WFID," %.16"ESYM"", data[IDX(i,2)]);
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

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of nonlinear iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy(y);                /* Free vectors */
  N_VDestroy(umask);
  N_VDestroy(vmask);
  N_VDestroy(wmask);
  SUNMatDestroy(udata->R);      /* Free user data */
  free(udata);
  ARKStepFree(&arkode_mem);     /* Free integrator memory */
  SUNLinSolFree(LS);            /* Free linear solver */
  SUNMatDestroy(A);             /* Free A matrix */
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;      /* access problem data */
  sunindextype N = udata->N;                     /* set variable shortcuts */
  realtype a     = udata->a;
  realtype b     = udata->b;
  realtype ep    = udata->ep;
  realtype du    = udata->du;
  realtype dv    = udata->dv;
  realtype dw    = udata->dw;
  realtype dx    = udata->dx;
  realtype *Ydata=NULL, *dYdata=NULL;
  realtype uconst, vconst, wconst, u, ul, ur, v, vl, vr, w, wl, wr;
  sunindextype i;

  Ydata = N_VGetArrayPointer(y);     /* access data arrays */
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
  dYdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                        /* initialize ydot to zero */

  /* iterate over domain, computing all equations */
  uconst = du/dx/dx;
  vconst = dv/dx/dx;
  wconst = dw/dx/dx;
  for (i=1; i<N-1; i++) {

    /* set shortcuts */
    u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* u_t = du*u_xx + a - (w+1)*u + v*u^2 */
    dYdata[IDX(i,0)] = (ul - TWO*u + ur)*uconst + a - (w+ONE)*u + v*u*u;

    /* v_t = dv*v_xx + w*u - v*u^2 */
    dYdata[IDX(i,1)] = (vl - TWO*v + vr)*vconst + w*u - v*u*u;

    /* w_t = dw*w_xx + (b-w)/ep - w*u */
    dYdata[IDX(i,2)] = (wl - TWO*w + wr)*wconst + (b-w)/ep - w*u;

  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  return 0;
}


/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  /* problem data */
  UserData udata = (UserData) user_data;
  sunindextype N = udata->N;

  /* ensure that Jac is the correct size */
  if ( (SUNSparseMatrix_Rows(J) != N*3) ||
       (SUNSparseMatrix_Columns(J) != N*3) ) {
    printf("Jacobian calculation error: matrix is the wrong size!\n");
    return 1;
  }

  /* Fill in the Laplace matrix */
  if (LaplaceMatrix(J, udata)) {
    printf("Jacobian calculation error in calling LaplaceMatrix!\n");
    return 1;
  }

  /* Create matrix for Jacobian of the reaction terms if necessary */
  if (udata->R == NULL) {
    udata->R = SUNSparseMatrix(SUNSparseMatrix_Rows(J),
                               SUNSparseMatrix_Columns(J),
                               SUNSparseMatrix_NNZ(J), CSC_MAT);
    if (udata->R == NULL) {
      printf("Jacobian calculation error in allocating R matrix!\n");
      return 1;
    }
  }

  /* Compute the Jacobian of the reaction terms */
  if (ReactionJac(y, udata->R, udata)) {
    printf("Jacobian calculation error in calling ReactionJac!\n");
    return 1;
  }

  /* Add R to J */
  if (SUNMatScaleAdd(ONE, J, udata->R) != 0) {
    printf("Jacobian calculation error in adding sparse matrices!\n");
    return 1;
  }

  return 0;
}




/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Routine to compute the stiffness matrix from (L*y) */
static int LaplaceMatrix(SUNMatrix Lap, UserData udata)
{
  sunindextype N = udata->N;  /* set shortcuts */
  sunindextype i, nz=0;
  realtype uconst, uconst2, vconst, vconst2, wconst, wconst2;
  sunindextype *colptrs = SUNSparseMatrix_IndexPointers(Lap);
  sunindextype *rowvals = SUNSparseMatrix_IndexValues(Lap);
  realtype *data = SUNSparseMatrix_Data(Lap);

  /* clear out matrix */
  SUNMatZero(Lap);

  /* set first column to zero */
  colptrs[IDX(0,0)] = nz;
  colptrs[IDX(0,1)] = nz;
  colptrs[IDX(0,2)] = nz;

  /* iterate over nodes, filling in Laplacian entries depending on these */
  uconst  = (udata->du)/(udata->dx)/(udata->dx);
  uconst2 = -TWO*uconst;
  vconst  = (udata->dv)/(udata->dx)/(udata->dx);
  vconst2 = -TWO*vconst;
  wconst  = (udata->dw)/(udata->dx)/(udata->dx);
  wconst2 = -TWO*wconst;
  for (i=1; i<N-1; i++) {

    /* dependence on u at this node */
    colptrs[IDX(i,0)] = nz;
    if (i>1) {                /* node to left */
      data[nz] = uconst;
      rowvals[nz++] = IDX(i-1,0);
    }

    data[nz] = uconst2;  /* self */
    rowvals[nz++] = IDX(i,0);

    if (i<N-2) {              /* node to right */
      data[nz] = uconst;
      rowvals[nz++] = IDX(i+1,0);
    }

    /* dependence on v at this node */
    colptrs[IDX(i,1)] = nz;
    if (i>1) {                /* node to left */
      data[nz] = vconst;
      rowvals[nz++] = IDX(i-1,1);
    }

    data[nz] = vconst2;  /* self */
    rowvals[nz++] = IDX(i,1);

    if (i<N-2) {              /* node to right */
      data[nz] = vconst;
      rowvals[nz++] = IDX(i+1,1);
    }

    /* dependence on w at this node */
    colptrs[IDX(i,2)] = nz;
    if (i>1) {                /* node to left */
      data[nz] = wconst;
      rowvals[nz++] = IDX(i-1,2);
    }

    data[nz] = wconst2;  /* self */
    rowvals[nz++] = IDX(i,2);

    if (i<N-2) {              /* node to right */
      data[nz] = wconst;
      rowvals[nz++] = IDX(i+1,2);
    }

  }

  /* set last column to zero */
  colptrs[IDX(N-1,0)] = nz;
  colptrs[IDX(N-1,1)] = nz;
  colptrs[IDX(N-1,2)] = nz;

  /* end of data */
  colptrs[IDX(N-1,2)+1] = nz;

  return 0;
}



/* Routine to compute the Jacobian matrix from R(y) */
static int ReactionJac(N_Vector y, SUNMatrix Jac, UserData udata)
{
  sunindextype N = udata->N;                   /* set shortcuts */
  sunindextype i, nz=0;
  realtype u, v, w;
  realtype ep = udata->ep;
  sunindextype *colptrs = SUNSparseMatrix_IndexPointers(Jac);
  sunindextype *rowvals = SUNSparseMatrix_IndexValues(Jac);
  realtype *data = SUNSparseMatrix_Data(Jac);
  realtype *Ydata = N_VGetArrayPointer(y);     /* access solution array */
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;

  /* clear out matrix */
  SUNMatZero(Jac);

  /* set first matrix column to zero */
  colptrs[IDX(0,0)] = 0;
  colptrs[IDX(0,1)] = 0;
  colptrs[IDX(0,2)] = 0;

  /* iterate over interior nodes, filling in Jacobian entries */
  for (i=1; i<N-1; i++) {

    /* set nodal value shortcuts */
    u = Ydata[IDX(i,0)];
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* dependence on u at this node */
    colptrs[IDX(i,0)] = nz;

    rowvals[nz] = IDX(i,0);        /* fu wrt u */
    data[nz++] = TWO*u*v - w - ONE;

    rowvals[nz] = IDX(i,1);        /* fv wrt u */
    data[nz++] = w - TWO*u*v;

    rowvals[nz] = IDX(i,2);        /* fw wrt u */
    data[nz++] = -w;

    /* dependence on v at this node */
    colptrs[IDX(i,1)] = nz;

    rowvals[nz] = IDX(i,0);        /* fu wrt v */
    data[nz++] = u*u;

    rowvals[nz] = IDX(i,1);        /* fv wrt v */
    data[nz++] = -u*u;

    /* dependence on w at this node */
    colptrs[IDX(i,2)] = nz;

    rowvals[nz] = IDX(i,0);        /* fu wrt w */
    data[nz++] = -u;

    rowvals[nz] = IDX(i,1);        /* fv wrt w */
    data[nz++] = u;

    rowvals[nz] = IDX(i,2);        /* fw wrt w */
    data[nz++] = -ONE/ep - u;

  }

  /* set last matrix column to zero */
  colptrs[IDX(N-1,0)] = nz;
  colptrs[IDX(N-1,1)] = nz;
  colptrs[IDX(N-1,2)] = nz;

  /* end of data */
  colptrs[IDX(N-1,2)+1] = nz;

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
