/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * Based on example problem ark_brusselator1D_FEM by
 * Daniel R. Reynolds @ SMU.
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
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is a PDE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    u_t = du*u_xx + a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx + w*u - v*u^2
 *    w_t = dw*w_xx + (b-w)/ep - w*u
 * for t in [0, 80], x in [0, 1], with initial conditions
 *    u(0,x) =  a  + 0.1*sin(pi*x)
 *    v(0,x) = b/a + 0.1*sin(pi*x)
 *    w(0,x) =  b  + 0.1*sin(pi*x),
 * and with stationary boundary conditions, i.e.
 *    u_t(t,0) = u_t(t,1) = 0
 *    v_t(t,0) = v_t(t,1) = 0
 *    w_t(t,0) = w_t(t,1) = 0.
 *
 * Here, we use a piecewise linear Galerkin finite element
 * discretization in space, where all element-wise integrals are
 * computed using 3-node Gaussian quadrature (since we will have
 * quartic polynomials in the reaction terms for the u_t and v_t
 * equations, including the test function).  The time derivative
 * terms for this system will include a mass matrix, giving rise
 * to an ODE system of the form
 *      M y_t = L y + R(y),
 * where M is the block mass matrix for each component, L is
 * the block Laplace operator for each component, and R(y) is
 * a 3x3 block comprised of the nonlinear reaction terms for
 * each component.  Since it it highly inefficient to rewrite
 * this system as
 *      y_t = M^{-1}(L y + R(y)),
 * we solve this system using ARKStep, with a user-supplied mass
 * matrix.  We therefore provide functions to evaluate the ODE RHS
 *    f(t,y) = L y + R(y),
 * its Jacobian
 *    J(t,y) = L + dR/dy,
 * and the mass matrix, M.
 *
 * This program solves the problem with the DIRK method, using a
 * Newton iteration with the SuperLU_DIST SUNLinearSolver.
 *
 * 100 outputs are printed at equal time intervals, and run
 * statistics are printed at the end.
 *---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <arkode/arkode_arkstep.h>            /* prototypes for ARKStep fcts., consts     */
#include <nvector/nvector_serial.h>           /* serial N_Vector types, fcts., macros     */
#include <sunmatrix/sunmatrix_slunrloc.h>     /* access to SuperLU SLU_NR_loc SUNMatrix   */
#include <sunlinsol/sunlinsol_superludist.h>  /* access to SuperLU_DIST SUNLinearSolver   */
#include <sundials/sundials_types.h>          /* defs. of realtype, sunindextype, etc     */

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
/* [variables are grouped according to spatial location] */
#define IDX(x,v) (3*(x)+v)

/* constants */
#define ZERO (RCONST(0.0))
#define ONE  (RCONST(1.0))
#define TWO  (RCONST(2.0))
#define HALF (RCONST(0.5))

/* Gaussian quadrature nodes, weights and formula (3 node, 7th-order accurate) */
#define X1(xl,xr)   (HALF*(xl+xr) - HALF*(xr-xl)*RCONST(0.774596669241483377035853079956))
#define X2(xl,xr)   (HALF*(xl+xr))
#define X3(xl,xr)   (HALF*(xl+xr) + HALF*(xr-xl)*RCONST(0.774596669241483377035853079956))
#define W1          (RCONST(0.55555555555555555555555555555556))
#define W2          (RCONST(0.88888888888888888888888888888889))
#define W3          (RCONST(0.55555555555555555555555555555556))
#define Quad(f1,f2,f3,xl,xr) (HALF*(xr-xl)*(W1*f1 + W2*f2 + W3*f3))

/* evaluation macros for variables, basis functions and basis derivatives */
#define ChiL(xl,xr,x) ((xr-x)/(xr-xl))
#define ChiR(xl,xr,x) ((x-xl)/(xr-xl))
#define ChiL_x(xl,xr) (ONE/(xl-xr))
#define ChiR_x(xl,xr) (ONE/(xr-xl))
#define Eval(ul,ur,xl,xr,x) (ul*ChiL(xl,xr,x) + ur*ChiR(xl,xr,x))
#define Eval_x(ul,ur,xl,xr) (ul*ChiL_x(xl,xr) + ur*ChiR_x(xl,xr))


/* user data structure */
typedef struct {
  sunindextype N;   /* number of intervals     */
  realtype *x;      /* mesh node locations     */
  realtype a;       /* constant forcing on u   */
  realtype b;       /* steady-state value of w */
  realtype du;      /* diffusion coeff for u   */
  realtype dv;      /* diffusion coeff for v   */
  realtype dw;      /* diffusion coeff for w   */
  realtype ep;      /* stiffness parameter     */
  N_Vector tmp;     /* temporary vector        */
  SUNMatrix R;      /* temporary storage       */
} *UserData;


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f_diff(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f_rx(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int MassMatrix(realtype t, SUNMatrix M, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private helper functions  */
static int LaplaceMatrix(SUNMatrix Jac, UserData udata);
static int ReactionJac(N_Vector y, SUNMatrix Jac, UserData udata);

/* Private function to check function return values */
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Main Program */
int main(int argc, char *argv[]) {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);    /* initial time */
  realtype Tf = RCONST(10.0);   /* final time */
  int Nt = 100;                 /* total number of output times */
  int Nvar = 3;                 /* number of solution fields */
  UserData udata = NULL;
  realtype *data;
  realtype *Mdata, *Adata, *Rdata;           /* CSR matrix data for solver and mass matrix */
  sunindextype *Acolind, *Mcolind, *Rcolind; /* CSR matrix column indices */
  sunindextype *Arowptr, *Mrowptr, *Rrowptr; /* CSR matrix row pointers */
  sunindextype N = 201;         /* spatial mesh size */
  realtype a = RCONST(0.6);             /* problem parameters */
  realtype b = RCONST(2.0);
  realtype du = RCONST(0.025);
  realtype dv = RCONST(0.025);
  realtype dw = RCONST(0.025);
  realtype ep = RCONST(1.0e-5);         /* stiffness parameter */
  realtype reltol = RCONST(1.0e-6);     /* tolerances */
  realtype abstol = RCONST(1.0e-10);
  sunindextype i, NEQ, NNZ;
  int npes, my_pe;

  /* general problem variables */
  int retval;                 /* reusable error-checking retval */
  N_Vector y = NULL;
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  SUNMatrix A = NULL;         /* empty matrix object for solver */
  SUNMatrix M = NULL;         /* empty mass matrix object */
  SUNLinearSolver LS = NULL;  /* empty linear solver object */
  SUNLinearSolver MLS = NULL; /* empty mass matrix solver object */
  SuperMatrix Asuper;         /* empty SuperLU-DIST matrix object for solver */
  SuperMatrix Msuper;         /* empty SuperLU-DIST mass matrix object */
  SuperMatrix Rsuper;         /* empty SuperLU-DIST reaction matrix object */
  gridinfo_t grid;            /* SuperLU-DIST process grid */
  dLUstruct_t Alu, Mlu;        /* SuperLU-DIST dLUstruct_t */
  dScalePermstruct_t Ascaleperm, Mscaleperm; /* SuperLU-DIST dScalePermstruct_t */
  dSOLVEstruct_t Asolve, Msolve;             /* SuperLU-DIST dSOLVEstruct_t */
  SuperLUStat_t Astat, Mstat;               /* SuperLU-DIST SuperLUState_t */
  superlu_dist_options_t Aopts, Mopts;      /* SuperLU-DIST options struct */

  void *arkode_mem = NULL;
  FILE *FID, *UFID, *VFID, *WFID;
  realtype h, z, t, dTout, tout, u, v, w, pi;
  int iout;
  long int nst, nst_a, nfe, nfi, nsetups, nje, nni, ncfn;
  long int netf, nmset, nms, nMv;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_pe);

  /* This example only allows 1 MPI rank because we are demonstrating
   * SuperLU_DIST on one node. */
  if (npes != 1 && my_pe == 0) {
    printf("WARNING: This example only uses 1 MPI rank\n");
  }

  /* initialize SuperLU-DIST grid */
  superlu_gridinit(MPI_COMM_WORLD, 1, 1, &grid);
  if (grid.iam >= 1) {
    superlu_gridexit(&grid);
    MPI_Finalize();
    return(0);
  }

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->x = NULL;
  udata->tmp = NULL;
  udata->R = NULL;
  if (check_retval((void *)udata, "malloc", 2)) MPI_Abort(grid.comm, 1);

  /* store the inputs in the UserData structure */
  udata->N  = N;
  udata->a  = a;
  udata->b  = b;
  udata->du = du;
  udata->dv = dv;
  udata->dw = dw;
  udata->ep = ep;

  /* set total allocated vector length (N-1 intervals, Dirichlet end points) */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D FEM Brusselator PDE test problem:\n");
  printf("    N = %li,  NEQ = %li\n", (long int) udata->N, (long int) NEQ);
  printf("    problem parameters:  a = %" GSYM ",  b = %" GSYM ",  ep = %" GSYM "\n",
         udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %" GSYM ",  dv = %" GSYM ",  dw = %" GSYM "\n",
         udata->du, udata->dv, udata->dw);
  printf("    reltol = %.1" ESYM ",  abstol = %.1" ESYM "\n\n", reltol, abstol);

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ);           /* Create serial vector for solution */
  if (check_retval((void *)y, "N_VNew_Serial", 0)) MPI_Abort(grid.comm, 1);

  data = N_VGetArrayPointer(y);     /* Access data array for new NVector y */
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) MPI_Abort(grid.comm, 1);

  umask = N_VNew_Serial(NEQ);       /* Create serial vector masks */
  if (check_retval((void *)umask, "N_VNew_Serial", 0)) MPI_Abort(grid.comm, 1);

  vmask = N_VNew_Serial(NEQ);
  if (check_retval((void *)vmask, "N_VNew_Serial", 0)) MPI_Abort(grid.comm, 1);

  wmask = N_VNew_Serial(NEQ);
  if (check_retval((void *)wmask, "N_VNew_Serial", 0)) MPI_Abort(grid.comm, 1);

  udata->tmp = N_VNew_Serial(NEQ);  /* temporary N_Vector inside udata */
  if (check_retval((void *) udata->tmp, "N_VNew_Serial", 0)) MPI_Abort(grid.comm, 1);

  /* allocate and set up spatial mesh; this [arbitrarily] clusters
     more intervals near the end points of the interval */
  udata->x = (realtype *) malloc(N*sizeof(realtype));
  if (check_retval((void *)udata->x, "malloc", 2)) MPI_Abort(grid.comm, 1);
  h = RCONST(10.0)/(N-1);
  for (i=0; i<N; i++) {
    z = -RCONST(5.0) + h*i;
    udata->x[i] = HALF/atan(RCONST(5.0))*atan(z) + HALF;
  }

  /* Set initial conditions into y */
  pi = RCONST(4.0)*atan(RCONST(1.0));
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*udata->x[i]);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*udata->x[i]);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*udata->x[i]);  /* w */
  }

  /* Set mask array values for each solution component */
  N_VConst(ZERO, umask);
  data = N_VGetArrayPointer(umask);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) MPI_Abort(grid.comm, 1);
  for (i=0; i<N; i++)  data[IDX(i,0)] = ONE;

  N_VConst(ZERO, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) MPI_Abort(grid.comm, 1);
  for (i=0; i<N; i++)  data[IDX(i,1)] = ONE;

  N_VConst(ZERO, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) MPI_Abort(grid.comm, 1);
  for (i=0; i<N; i++)  data[IDX(i,2)] = ONE;


  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y), the inital time
     T0, and the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  arkode_mem = ARKStepCreate(NULL, f, T0, y);
  if (check_retval((void *)arkode_mem, "ARKStepCreate", 0)) MPI_Abort(grid.comm, 1);

  /* Set routines */

  /* Pass udata to user functions */
  retval = ARKStepSetUserData(arkode_mem, (void *) udata);
  if (check_retval(&retval, "ARKStepSetUserData", 1)) MPI_Abort(grid.comm, 1);

  /* Specify tolerances */
  retval = ARKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_retval(&retval, "ARKStepSStolerances", 1)) MPI_Abort(grid.comm, 1);

  /* Specify residual tolerance */
  retval = ARKStepResStolerance(arkode_mem, abstol);
  if (check_retval(&retval, "ARKStepResStolerance", 1)) MPI_Abort(grid.comm, 1);

  /* Initialize sparse matrix data structure and SuperLU_DIST solvers (system and mass) */
  NNZ = 15*NEQ;

  Adata = NULL;
  Adata = (realtype *) malloc(NNZ*sizeof(realtype));
  if (check_retval((void *)Adata, "malloc Adata", 2)) MPI_Abort(grid.comm, 1);

  Acolind = NULL;
  Acolind = (sunindextype *) malloc(NNZ*sizeof(sunindextype));
  if (check_retval((void *)Acolind, "malloc Acolind", 2)) MPI_Abort(grid.comm, 1);

  Arowptr = NULL;
  Arowptr = (sunindextype *) malloc((NEQ+1)*sizeof(sunindextype));
  if (check_retval((void *)Arowptr, "malloc Arowptr", 2)) MPI_Abort(grid.comm, 1);

  Mdata = NULL;
  Mdata = (realtype *) malloc(NNZ*sizeof(realtype));
  if (check_retval((void *)Mdata, "malloc Mdata", 2)) MPI_Abort(grid.comm, 1);

  Mcolind = NULL;
  Mcolind = (sunindextype *) malloc(NNZ*sizeof(sunindextype));
  if (check_retval((void *)Mcolind, "malloc Mcolind", 2)) MPI_Abort(grid.comm, 1);

  Mrowptr = NULL;
  Mrowptr = (sunindextype *) malloc((NEQ+1)*sizeof(sunindextype));
  if (check_retval((void *)Mrowptr, "malloc Mrowptr", 2)) MPI_Abort(grid.comm, 1);

  Rdata = NULL;
  Rdata = (realtype *) malloc(NNZ*sizeof(realtype));
  if (check_retval((void *)Rdata, "malloc Rdata", 2)) MPI_Abort(grid.comm, 1);

  Rcolind = NULL;
  Rcolind = (sunindextype *) malloc(NNZ*sizeof(sunindextype));
  if (check_retval((void *)Rcolind, "malloc Rcolind", 2)) MPI_Abort(grid.comm, 1);

  Rrowptr = NULL;
  Rrowptr = (sunindextype *) malloc((NEQ+1)*sizeof(sunindextype));
  if (check_retval((void *)Rrowptr, "malloc Rrowptr", 2)) MPI_Abort(grid.comm, 1);

  /* SuperLU_DIST structures */
  dCreate_CompRowLoc_Matrix_dist(&Asuper, NEQ, NEQ, NNZ, NEQ, 0,
                                 Adata, Acolind, Arowptr, SLU_NR_loc, SLU_D, SLU_GE);
  dScalePermstructInit(NEQ, NEQ, &Ascaleperm);
  dLUstructInit(NEQ, &Alu);
  PStatInit(&Astat);
  set_default_options_dist(&Aopts);
  Aopts.PrintStat = NO;

  dCreate_CompRowLoc_Matrix_dist(&Msuper, NEQ, NEQ, NNZ, NEQ, 0,
                                 Mdata, Mcolind, Mrowptr, SLU_NR_loc, SLU_D, SLU_GE);
  dScalePermstructInit(NEQ, NEQ, &Mscaleperm);
  dLUstructInit(NEQ, &Mlu);
  PStatInit(&Mstat);
  set_default_options_dist(&Mopts);
  Mopts.PrintStat = NO;

  dCreate_CompRowLoc_Matrix_dist(&Rsuper, NEQ, NEQ, NNZ, NEQ, 0,
                                 Rdata, Rcolind, Rrowptr, SLU_NR_loc, SLU_D, SLU_GE);

  /* SUNDIALS structures */
  A = SUNMatrix_SLUNRloc(&Asuper, &grid);
  if (check_retval((void *)A, "SUNMatrix_SLUNRloc", 0)) MPI_Abort(grid.comm, 1);

  LS = SUNLinSol_SuperLUDIST(y, A, &grid, &Alu, &Ascaleperm, &Asolve, &Astat, &Aopts);
  if (check_retval((void *)LS, "SUNLinSol_SuperLUDIST", 0)) MPI_Abort(grid.comm, 1);

  M = SUNMatrix_SLUNRloc(&Msuper, &grid);
  if (check_retval((void *)M, "SUNMatrix_SLUNRloc", 0)) MPI_Abort(grid.comm, 1);

  MLS = SUNLinSol_SuperLUDIST(y, M, &grid, &Mlu, &Mscaleperm, &Msolve, &Mstat, &Mopts);
  if (check_retval((void *)MLS, "SUNLinSol_SuperLUDIST", 0)) MPI_Abort(grid.comm, 1);

  udata->R = SUNMatrix_SLUNRloc(&Rsuper, &grid);
  if (check_retval((void *)udata->R, "SUNMatrix_SLUNRloc", 0)) MPI_Abort(grid.comm, 1);

  /* Attach the matrix, linear solver, and Jacobian construction routine to ARKStep */
  retval = ARKStepSetLinearSolver(arkode_mem, LS, A);
  if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) MPI_Abort(grid.comm, 1);

  /* Supply Jac routine */
  retval = ARKStepSetJacFn(arkode_mem, Jac);
  if (check_retval(&retval, "ARKStepSetJacFn", 1)) MPI_Abort(grid.comm, 1);

  /* Attach the mass matrix, linear solver and construction routines to ARKStep;
     notify ARKStep that the mass matrix is not time-dependent */
  retval = ARKStepSetMassLinearSolver(arkode_mem, MLS, M, SUNFALSE);
  if (check_retval(&retval, "ARKStepSetMassLinearSolver", 1)) MPI_Abort(grid.comm, 1);

  /* Supply M routine */
  retval = ARKStepSetMassFn(arkode_mem, MassMatrix);
  if (check_retval(&retval, "ARKStepSetMassFn", 1)) MPI_Abort(grid.comm, 1);

  /* output mesh to disk */
  FID=fopen("bruss_FEM_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16" ESYM "\n", udata->x[i]);
  fclose(FID);

  /* Open output stream for results, access data arrays */
  UFID = fopen("bruss_FEM_u.txt","w");
  VFID = fopen("bruss_FEM_v.txt","w");
  WFID = fopen("bruss_FEM_w.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) MPI_Abort(grid.comm, 1);

  /* output initial condition to disk */
  for (i=0; i<N; i++)  fprintf(UFID," %.16" ESYM, data[IDX(i,0)]);
  for (i=0; i<N; i++)  fprintf(VFID," %.16" ESYM, data[IDX(i,1)]);
  for (i=0; i<N; i++)  fprintf(WFID," %.16" ESYM, data[IDX(i,2)]);
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

    retval = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);  /* call integrator */
    if (check_retval(&retval, "ARKStepEvolve", 1)) break;
    u = N_VWL2Norm(y,umask);                               /* access/print solution statistics */
    u = sqrt(u*u/N);
    v = N_VWL2Norm(y,vmask);
    v = sqrt(v*v/N);
    w = N_VWL2Norm(y,wmask);
    w = sqrt(w*w/N);
    printf("  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "  %10.6" FSYM "\n", t, u, v, w);
    if (retval >= 0) {                                     /* successful solve: update output time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                               /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }

    /* output results to disk */
    for (i=0; i<N; i++)  fprintf(UFID," %.16" ESYM, data[IDX(i,0)]);
    for (i=0; i<N; i++)  fprintf(VFID," %.16" ESYM, data[IDX(i,1)]);
    for (i=0; i<N; i++)  fprintf(WFID," %.16" ESYM, data[IDX(i,2)]);
    fprintf(UFID,"\n");
    fprintf(VFID,"\n");
    fprintf(WFID,"\n");
  }
  printf("   ----------------------------------------------\n");
  fclose(UFID);
  fclose(VFID);
  fclose(WFID);

  /* Print some final statistics */
  retval = ARKStepGetNumSteps(arkode_mem, &nst);
  check_retval(&retval, "ARKStepGetNumSteps", 1);
  retval = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_retval(&retval, "ARKStepGetNumStepAttempts", 1);
  retval = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_retval(&retval, "ARKStepGetNumRhsEvals", 1);
  retval = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_retval(&retval, "ARKStepGetNumLinSolvSetups", 1);
  retval = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  check_retval(&retval, "ARKStepGetNumErrTestFails", 1);
  retval = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  check_retval(&retval, "ARKStepGetNumNonlinSolvIters", 1);
  retval = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_retval(&retval, "ARKStepGetNumNonlinSolvConvFails", 1);
  retval = ARKStepGetNumMassSetups(arkode_mem, &nmset);
  check_retval(&retval, "ARKStepGetNumMassSetups", 1);
  retval = ARKStepGetNumMassSolves(arkode_mem, &nms);
  check_retval(&retval, "ARKStepGetNumMassSolves", 1);
  retval = ARKStepGetNumMassMult(arkode_mem, &nMv);
  check_retval(&retval, "ARKStepGetNumMassMult", 1);
  retval = ARKStepGetNumJacEvals(arkode_mem, &nje);
  check_retval(&retval, "ARKStepGetNumJacEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total mass matrix setups = %li\n", nmset);
  printf("   Total mass matrix solves = %li\n", nms);
  printf("   Total mass times evals = %li\n", nMv);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy(y);                   /* Free vectors */
  N_VDestroy(umask);
  N_VDestroy(vmask);
  N_VDestroy(wmask);
  ARKStepFree(&arkode_mem);        /* Free integrator memory */
  SUNLinSolFree(LS);               /* Free linear solvers */
  SUNLinSolFree(MLS);
  SUNMatDestroy(A);                /* Free matrices */
  SUNMatDestroy(M);
  SUNMatDestroy(udata->R);         /* Free user data */
  N_VDestroy(udata->tmp);
  free(udata->x);
  free(udata);

  /* Free the SuperLU_DIST structures */
  PStatFree(&Astat); PStatFree(&Mstat);
  dScalePermstructFree(&Ascaleperm); dScalePermstructFree(&Mscaleperm);
  dLUstructFree(&Alu); dLUstructFree(&Mlu);
  Destroy_CompRowLoc_Matrix_dist(&Asuper);
  Destroy_CompRowLoc_Matrix_dist(&Msuper);
  Destroy_CompRowLoc_Matrix_dist(&Rsuper);
  superlu_gridexit(&grid);

  MPI_Finalize();

  return 0;
}


/*------------------------------
  Functions called by the solver
 *------------------------------*/


/* Routine to compute the ODE RHS function f(t,y), where system is of the form
        M y_t = f(t,y) := Ly + R(y)
   This routine only computes the f(t,y), leaving (M y_t) alone. */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* local data */
  int ier;

  /* clear out RHS (to be careful) */
  N_VConst(ZERO, ydot);

  /* add reaction terms to RHS */
  ier = f_rx(t, y, ydot, user_data);
  if (ier != 0)  return ier;

  /* add diffusion terms to RHS */
  ier = f_diff(t, y, ydot, user_data);
  if (ier != 0)  return ier;

  return 0;
}


/* Routine to compute the diffusion portion of the ODE RHS function f(t,y). */
static int f_diff(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  sunindextype N = udata->N;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;

  /* local variables */
  sunindextype i;
  realtype ul, ur, vl, vr, wl, wr;
  realtype xl, xr, f1;
  booleantype left, right;
  realtype *Ydata, *RHSdata;

  /* access data arrays */
  Ydata = N_VGetArrayPointer(y);
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return(1);
  RHSdata = N_VGetArrayPointer(ydot);
  if (check_retval((void *)RHSdata, "N_VGetArrayPointer", 0)) return(1);

  /* iterate over intervals, filling in residual function */
  for (i=0; i<N-1; i++) {

    /* set booleans to determine whether equations exist on the left/right */
    left  = (i==0)     ? SUNFALSE : SUNTRUE;
    right = (i==(N-2)) ? SUNFALSE : SUNTRUE;

    /* set nodal value shortcuts (interval index aligns with left node) */
    ul = Ydata[IDX(i,0)];
    vl = Ydata[IDX(i,1)];
    wl = Ydata[IDX(i,2)];
    ur = Ydata[IDX(i+1,0)];
    vr = Ydata[IDX(i+1,1)];
    wr = Ydata[IDX(i+1,2)];

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* evaluate L*y on this subinterval
       NOTE: all f values are the same since constant on interval */
    /*    left test function */
    if (left) {
      /*  u */
      f1 = -du * Eval_x(ul,ur,xl,xr) * ChiL_x(xl,xr);
      RHSdata[IDX(i,0)] += Quad(f1,f1,f1,xl,xr);

      /*  v */
      f1 = -dv * Eval_x(vl,vr,xl,xr) * ChiL_x(xl,xr);
      RHSdata[IDX(i,1)] += Quad(f1,f1,f1,xl,xr);

      /*  w */
      f1 = -dw * Eval_x(wl,wr,xl,xr) * ChiL_x(xl,xr);
      RHSdata[IDX(i,2)] += Quad(f1,f1,f1,xl,xr);
    }
    /*    right test function */
    if (right) {
      /*  u */
      f1 = -du * Eval_x(ul,ur,xl,xr) * ChiR_x(xl,xr);
      RHSdata[IDX(i+1,0)] += Quad(f1,f1,f1,xl,xr);

      /*  v */
      f1 = -dv * Eval_x(vl,vr,xl,xr) * ChiR_x(xl,xr);
      RHSdata[IDX(i+1,1)] += Quad(f1,f1,f1,xl,xr);

      /*  w */
      f1 = -dw * Eval_x(wl,wr,xl,xr) * ChiR_x(xl,xr);
      RHSdata[IDX(i+1,2)] += Quad(f1,f1,f1,xl,xr);
    }
  }

  return 0;
}



/* Routine to compute the reaction portion of the ODE RHS function f(t,y). */
static int f_rx(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  sunindextype N = udata->N;
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;

  /* local variables */
  sunindextype i;
  realtype ul, ur, vl, vr, wl, wr;
  realtype u, v, w, xl, xr, f1, f2, f3;
  booleantype left, right;
  realtype *Ydata, *RHSdata;

  /* access data arrays */
  Ydata = N_VGetArrayPointer(y);
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return(1);
  RHSdata = N_VGetArrayPointer(ydot);
  if (check_retval((void *)RHSdata, "N_VGetArrayPointer", 0)) return(1);

  /* iterate over intervals, filling in residual function */
  for (i=0; i<N-1; i++) {

    /* set booleans to determine whether equations exist on the left/right */
    left  = (i==0)     ? SUNFALSE : SUNTRUE;
    right = (i==(N-2)) ? SUNFALSE : SUNTRUE;

    /* set nodal value shortcuts (interval index aligns with left node) */
    ul = Ydata[IDX(i,0)];
    vl = Ydata[IDX(i,1)];
    wl = Ydata[IDX(i,2)];
    ur = Ydata[IDX(i+1,0)];
    vr = Ydata[IDX(i+1,1)];
    wr = Ydata[IDX(i+1,2)];

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* evaluate R(y) on this subinterval */
    /*    left test function */
    if (left) {
      /*  u */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (a - (w+ONE)*u + v*u*u) * ChiL(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (a - (w+ONE)*u + v*u*u) * ChiL(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (a - (w+ONE)*u + v*u*u) * ChiL(xl,xr,X3(xl,xr));
      RHSdata[IDX(i,0)] += Quad(f1,f2,f3,xl,xr);

      /*  v */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (w*u - v*u*u) * ChiL(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (w*u - v*u*u) * ChiL(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (w*u - v*u*u) * ChiL(xl,xr,X3(xl,xr));
      RHSdata[IDX(i,1)] += Quad(f1,f2,f3,xl,xr);

      /*  w */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = ((b-w)/ep - w*u) * ChiL(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = ((b-w)/ep - w*u) * ChiL(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = ((b-w)/ep - w*u) * ChiL(xl,xr,X3(xl,xr));
      RHSdata[IDX(i,2)] += Quad(f1,f2,f3,xl,xr);
    }
    /*    right test function */
    if (right) {
      /*  u */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (a - (w+ONE)*u + v*u*u) * ChiR(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (a - (w+ONE)*u + v*u*u) * ChiR(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (a - (w+ONE)*u + v*u*u) * ChiR(xl,xr,X3(xl,xr));
      RHSdata[IDX(i+1,0)] += Quad(f1,f2,f3,xl,xr);

      /*  v */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (w*u - v*u*u) * ChiR(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (w*u - v*u*u) * ChiR(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (w*u - v*u*u) * ChiR(xl,xr,X3(xl,xr));
      RHSdata[IDX(i+1,1)] += Quad(f1,f2,f3,xl,xr);

      /*  w */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = ((b-w)/ep - w*u) * ChiR(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = ((b-w)/ep - w*u) * ChiR(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = ((b-w)/ep - w*u) * ChiR(xl,xr,X3(xl,xr));
      RHSdata[IDX(i+1,2)] += Quad(f1,f2,f3,xl,xr);
    }
  }

  return 0;
}



/* Interface routine to compute the Jacobian of the full RHS function, f(y) */
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  /* temporary variables */
  int ier;
  UserData udata = (UserData) user_data;
  SuperMatrix *Jsuper  = SUNMatrix_SLUNRloc_SuperMatrix(J);
  sunindextype NEQ = 3*udata->N;

  /* ensure that Jac is the correct size */
  if (Jsuper->nrow != NEQ || Jsuper->ncol != NEQ) {
    printf("Jacobian calculation error: matrix is the wrong size!\n");
    return 1;
  }

  /* Fill in the Laplace matrix */
  ier = LaplaceMatrix(J, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling Laplace matrix = %i\n",ier);
    return 1;
  }

  if (udata->R == NULL) {
    printf("Jac: R matrix is NULL!\n");
    return 1;
  }

  /* Add in the Jacobian of the reaction terms matrix */
  ier = ReactionJac(y, udata->R, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling reaction Jacobian = %i\n",ier);
    return 1;
  }

  /* Add R to J */
  ier = SUNMatScaleAdd(ONE, J, udata->R);
  if (ier != 0) {
    printf("Jac: error in adding sparse matrices = %i!\n",ier);
    return 1;
  }

  return 0;
}



/* Routine to compute the mass matrix multiplying y_t. */
static int MassMatrix(realtype t, SUNMatrix M, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  /* user data structure */
  UserData udata = (UserData) user_data;

  /* set shortcuts */
  sunindextype N = udata->N, NEQ = 3*N;
  SuperMatrix *Msuper  = SUNMatrix_SLUNRloc_SuperMatrix(M);
  NRformat_loc *Mstore = (NRformat_loc *) Msuper->Store;
  sunindextype *rowptrs = Mstore->rowptr;
  sunindextype *colinds = Mstore->colind;
  realtype *Mdata = (realtype *) Mstore->nzval;
  realtype *Xdata = udata->x;

  /* local data */
  sunindextype i, nz=0;
  realtype xl, xc, xr, Ml, Mc, Mr, ChiL1, ChiL2, ChiL3, ChiR1, ChiR2, ChiR3;
  booleantype left, right, interior;

  /* check that vector/matrix dimensions match up */
  if ((Msuper->nrow != NEQ) || (Msuper->ncol != NEQ)) {
    printf("MassMatrix calculation error: matrix is wrong size!\n");
    return 1;
  }

  /* clear out mass matrix */
  SUNMatZero(M);

  /* initialize node value shortcuts */
  xl = ZERO;
  xc = ZERO;
  xr = ZERO;

  /* iterate through nodes, filling in matrix by rows */
  for (i=0; i<N; i++) {

    /* set booleans to determine whether intervals exist on the left/right */
    left  = (i==0)     ? SUNFALSE : SUNTRUE;
    right = (i==(N-1)) ? SUNFALSE : SUNTRUE;
    interior = (left && right) ? SUNTRUE : SUNFALSE;

    /* set nodal value shortcuts (interval index aligns with left node) */
    if (left) xl = Xdata[i-1];
    xc = Xdata[i];
    if (right) xr = Xdata[i+1];

    /* compute entries of all mass matrix rows at node ix */
    Ml = ZERO;
    Mc = ZERO;
    Mr = ZERO;

    /* first compute dependence on values to left and center */
    if (left) {
      ChiL1 = ChiL(xl, xc, X1(xl,xc));
      ChiL2 = ChiL(xl, xc, X2(xl,xc));
      ChiL3 = ChiL(xl, xc, X3(xl,xc));
      ChiR1 = ChiR(xl, xc, X1(xl,xc));
      ChiR2 = ChiR(xl, xc, X2(xl,xc));
      ChiR3 = ChiR(xl, xc, X3(xl,xc));
      Ml = Ml + Quad(ChiL1*ChiR1, ChiL2*ChiR2, ChiL3*ChiR3, xl, xc);
      Mc = Mc + Quad(ChiR1*ChiR1, ChiR2*ChiR2, ChiR3*ChiR3, xl, xc);
    }

    /* second compute dependence on values to center and right */
    if (right) {
      ChiL1 = ChiL(xc, xr, X1(xc,xr));
      ChiL2 = ChiL(xc, xr, X2(xc,xr));
      ChiL3 = ChiL(xc, xr, X3(xc,xr));
      ChiR1 = ChiR(xc, xr, X1(xc,xr));
      ChiR2 = ChiR(xc, xr, X2(xc,xr));
      ChiR3 = ChiR(xc, xr, X3(xc,xr));
      Mc = Mc + Quad(ChiL1*ChiL1, ChiL2*ChiL2, ChiL3*ChiL3, xc, xr);
      Mr = Mr + Quad(ChiL1*ChiR1, ChiL2*ChiR2, ChiL3*ChiR3, xc, xr);
    }

    /* insert mass matrix entires into CSR matrix structure */

    /* u row */
    rowptrs[IDX(i,0)] = nz;

    if (left) {
      Mdata[nz] = Ml;
      colinds[nz] = IDX(i-1,0);
      nz++;
      if (interior) {
        Mdata[nz]   = ZERO;
        Mdata[nz+1] = ZERO;
        colinds[nz]   = IDX(i-1,1);
        colinds[nz+1] = IDX(i-1,2);
        nz += 2;
      }
    }

    Mdata[nz] = Mc;
    colinds[nz] = IDX(i,0);
    nz++;
    if (interior) {
      Mdata[nz]   = ZERO;
      Mdata[nz+1] = ZERO;
      colinds[nz]   = IDX(i,1);
      colinds[nz+1] = IDX(i,2);
      nz += 2;
    }

    if (right) {
      Mdata[nz] = Mr;
      colinds[nz] = IDX(i+1,0);
      nz++;
      if (interior) {
        Mdata[nz]   = ZERO;
        Mdata[nz+1] = ZERO;
        colinds[nz]   = IDX(i+1,1);
        colinds[nz+1] = IDX(i+1,2);
        nz += 2;
      }
    }

    /* v row */
    rowptrs[IDX(i,1)] = nz;

    if (left) {
      if (interior) {
        Mdata[nz] = ZERO;
        colinds[nz] = IDX(i-1,0);
        nz++;
      }
      Mdata[nz] = Ml;
      colinds[nz] = IDX(i-1,1);
      nz++;
      if (interior) {
        Mdata[nz] = ZERO;
        colinds[nz+1] = IDX(i-1,2);
        nz++;
      }
    }

    if (interior) {
      Mdata[nz] = ZERO;
      colinds[nz] = IDX(i,0);
      nz++;
    }
    Mdata[nz] = Mc;
    colinds[nz] = IDX(i,1);
    nz++;
    if (interior) {
      Mdata[nz] = ZERO;
      colinds[nz] = IDX(i,2);
      nz++;
    }

    if (right) {
      if (interior) {
        Mdata[nz] = ZERO;
        colinds[nz] = IDX(i+1,0);
        nz++;
      }
      Mdata[nz] = Mr;
      colinds[nz] = IDX(i+1,1);
      nz++;
      if (interior) {
        Mdata[nz] = ZERO;
        colinds[nz] = IDX(i+1,2);
        nz++;
      }
    }

    /* w row */
    rowptrs[IDX(i,2)] = nz;

    if (left) {
      if (interior) {
        Mdata[nz]   = ZERO;
        Mdata[nz+1] = ZERO;
        colinds[nz]   = IDX(i-1,0);
        colinds[nz+1] = IDX(i-1,1);
        nz += 2;
      }
      Mdata[nz] = Ml;
      colinds[nz] = IDX(i-1,2);
      nz++;
    }

    if (interior) {
      Mdata[nz]   = ZERO;
      Mdata[nz+1] = ZERO;
      colinds[nz]   = IDX(i,0);
      colinds[nz+1] = IDX(i,1);
      nz += 2;
    }
    Mdata[nz] = Mc;
    colinds[nz] = IDX(i,2);
    nz++;

    if (right) {
      if (interior) {
        Mdata[nz]   = ZERO;
        Mdata[nz+1] = ZERO;
        colinds[nz]   = IDX(i+1,0);
        colinds[nz+1] = IDX(i+1,1);
        nz += 2;
      }
      Mdata[nz] = Mr;
      colinds[nz] = IDX(i+1,2);
      nz++;
    }

  }

  /* signal end of data */
  rowptrs[IDX(N-1,2)+1] = nz;

  return 0;
}


/*-------------------------------
 * Private helper functions
 *-------------------------------*/


/* Routine to compute the Laplace matrix */
static int LaplaceMatrix(SUNMatrix L, UserData udata)
{
  /* set shortcuts */
  SuperMatrix *Lsuper  = SUNMatrix_SLUNRloc_SuperMatrix(L);
  NRformat_loc *Lstore = (NRformat_loc *) Lsuper->Store;
  sunindextype N = udata->N;
  sunindextype *rowptrs = Lstore->rowptr;
  sunindextype *colinds = Lstore->colind;
  realtype *Ldata = (realtype *) Lstore->nzval;
  realtype *Xdata = udata->x;
  realtype du = udata->du, dv = udata->dv, dw = udata->dw;

  /* set local variables */
  sunindextype i, j, nz=0;
  realtype xl, xc, xr;
  realtype Lu[9], Lv[9], Lw[9];

  /* initialize all local variables to zero (to avoid uninitialized variable warnings) */
  xl = xc = xr = ZERO;

  /* clear out matrix */
  SUNMatZero(L);

  /* Dirichlet boundary at left */
  rowptrs[IDX(0,0)] = nz;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(0,0);
  nz++;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(1,0);
  nz++;

  rowptrs[IDX(0,1)] = nz;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(0,1);
  nz++;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(1,1);
  nz++;

  rowptrs[IDX(0,2)] = nz;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(0,2);
  nz++;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(1,2);
  nz++;

  /* iterate over columns, filling in Laplace matrix */
  for (i=1; i<(N-1); i++) {

    /* set nodal value shortcuts (interval index aligns with left node) */
    xl = Xdata[i-1];
    xc = Xdata[i];
    xr = Xdata[i+1];

    /* compute entries of all Jacobian rows at node i */
    for(j=0; j<9; j++) {
      Lu[j] = ZERO;
      Lv[j] = ZERO;
      Lw[j] = ZERO;
    }

    /* first compute dependence on values to left and center */

    /* compute diffusion Jacobian components */

    /* L_u = -du * u_x * ChiR_x */
    /*   dL_u/dul   */
    Lu[IDX(0,0)] = (-du) * Quad(ONE,ONE,ONE,xl,xc) * ChiL_x(xl,xc) * ChiR_x(xl,xc);
    /*   dL_u/duc   */
    Lu[IDX(0,1)] = (-du) * Quad(ONE,ONE,ONE,xl,xc) * ChiR_x(xl,xc) * ChiR_x(xl,xc);

    /* L_v = -dv * v_x * ChiR_x */
    /*   dL_v/dvl   */
    Lv[IDX(1,0)] = (-dv) * Quad(ONE,ONE,ONE,xl,xc) * ChiL_x(xl,xc) * ChiR_x(xl,xc);
    /*   dL_v/dvc   */
    Lv[IDX(1,1)] = (-dv) * Quad(ONE,ONE,ONE,xl,xc) * ChiR_x(xl,xc) * ChiR_x(xl,xc);

    /* L_w =  -dw * w_x * ChiR_x */
    /*   dL_w/dwl   */
    Lw[IDX(2,0)] = (-dw) * Quad(ONE,ONE,ONE,xl,xc) * ChiL_x(xl,xc) * ChiR_x(xl,xc);
    /*   dL_w/dwc   */
    Lw[IDX(2,1)] = (-dw) * Quad(ONE,ONE,ONE,xl,xc) * ChiR_x(xl,xc) * ChiR_x(xl,xc);


    /* second compute dependence on values to center and right */

    /* compute diffusion Jacobian components */

    /* L_u = -du * u_x * ChiL_x */
    /*    dL_u/duc    */
    Lu[IDX(0,1)] = Lu[IDX(0,1)] + (-du) * Quad(ONE,ONE,ONE,xc,xr) * ChiL_x(xc,xr) * ChiL_x(xc,xr);

    /*    dL_u/dur    */
    Lu[IDX(0,2)] = Lu[IDX(0,2)] + (-du) * Quad(ONE,ONE,ONE,xc,xr) * ChiL_x(xc,xr) * ChiR_x(xc,xr);

    /* L_v = -dv * v_x * ChiL_x */
    /*    dL_v/dvc    */
    Lv[IDX(1,1)] = Lv[IDX(1,1)] + (-dv) * Quad(ONE,ONE,ONE,xc,xr) * ChiL_x(xc,xr) * ChiL_x(xc,xr);

    /*    dL_v/dvr    */
    Lv[IDX(1,2)] = Lv[IDX(1,2)] + (-dv) * Quad(ONE,ONE,ONE,xc,xr) * ChiL_x(xc,xr) * ChiR_x(xc,xr);

    /* L_w =  -dw * w_x * ChiL_x */
    /*    dL_w/dwc    */
    Lw[IDX(2,1)] = Lw[IDX(2,1)] + (-dw) * Quad(ONE,ONE,ONE,xc,xr) * ChiL_x(xc,xr) * ChiL_x(xc,xr);

    /*    dL_w/dwr    */
    Lw[IDX(2,2)] = Lw[IDX(2,2)] + (-dw) * Quad(ONE,ONE,ONE,xc,xr) * ChiL_x(xc,xr) * ChiR_x(xc,xr);


    /* insert Jacobian entries into CSR matrix structure */

    /* Lu row */
    rowptrs[IDX(i,0)] = nz;

    Ldata[nz] = Lu[IDX(0,0)];
    Ldata[nz+1] = Lu[IDX(1,0)];
    Ldata[nz+2] = Lu[IDX(2,0)];
    colinds[nz] = IDX(i-1,0);
    colinds[nz+1] = IDX(i-1,1);
    colinds[nz+2] = IDX(i-1,2);
    nz += 3;

    Ldata[nz] = Lu[IDX(0,1)];
    Ldata[nz+1] = Lu[IDX(1,1)];
    Ldata[nz+2] = Lu[IDX(2,1)];
    colinds[nz] = IDX(i,0);
    colinds[nz+1] = IDX(i,1);
    colinds[nz+2] = IDX(i,2);
    nz += 3;

    Ldata[nz] = Lu[IDX(0,2)];
    Ldata[nz+1] = Lu[IDX(1,2)];
    Ldata[nz+2] = Lu[IDX(2,2)];
    colinds[nz] = IDX(i+1,0);
    colinds[nz+1] = IDX(i+1,1);
    colinds[nz+2] = IDX(i+1,2);
    nz += 3;

    /* Lv row */
    rowptrs[IDX(i,1)] = nz;

    Ldata[nz] = Lv[IDX(0,0)];
    Ldata[nz+1] = Lv[IDX(1,0)];
    Ldata[nz+2] = Lv[IDX(2,0)];
    colinds[nz] = IDX(i-1,0);
    colinds[nz+1] = IDX(i-1,1);
    colinds[nz+2] = IDX(i-1,2);
    nz += 3;

    Ldata[nz] = Lv[IDX(0,1)];
    Ldata[nz+1] = Lv[IDX(1,1)];
    Ldata[nz+2] = Lv[IDX(2,1)];
    colinds[nz] = IDX(i,0);
    colinds[nz+1] = IDX(i,1);
    colinds[nz+2] = IDX(i,2);
    nz += 3;

    Ldata[nz] = Lv[IDX(0,2)];
    Ldata[nz+1] = Lv[IDX(1,2)];
    Ldata[nz+2] = Lv[IDX(2,2)];
    colinds[nz] = IDX(i+1,0);
    colinds[nz+1] = IDX(i+1,1);
    colinds[nz+2] = IDX(i+1,2);
    nz += 3;

    /* Lw row */
    rowptrs[IDX(i,2)] = nz;

    Ldata[nz] = Lw[IDX(0,0)];
    Ldata[nz+1] = Lw[IDX(1,0)];
    Ldata[nz+2] = Lw[IDX(2,0)];
    colinds[nz] = IDX(i-1,0);
    colinds[nz+1] = IDX(i-1,1);
    colinds[nz+2] = IDX(i-1,2);
    nz += 3;

    Ldata[nz] = Lw[IDX(0,1)];
    Ldata[nz+1] = Lw[IDX(1,1)];
    Ldata[nz+2] = Lw[IDX(2,1)];
    colinds[nz] = IDX(i,0);
    colinds[nz+1] = IDX(i,1);
    colinds[nz+2] = IDX(i,2);
    nz += 3;

    Ldata[nz] = Lw[IDX(0,2)];
    Ldata[nz+1] = Lw[IDX(1,2)];
    Ldata[nz+2] = Lw[IDX(2,2)];
    colinds[nz] = IDX(i+1,0);
    colinds[nz+1] = IDX(i+1,1);
    colinds[nz+2] = IDX(i+1,2);
    nz += 3;

  }

  /* Dirichlet boundary at right */
  rowptrs[IDX(N-1,0)] = nz;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(N-2,0);
  nz++;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(N-1,0);
  nz++;

  rowptrs[IDX(N-1,1)] = nz;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(N-2,1);
  nz++;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(N-1,1);
  nz++;

  rowptrs[IDX(N-1,2)] = nz;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(N-2,2);
  nz++;

  Ldata[nz]   = ZERO;
  colinds[nz] = IDX(N-1,2);
  nz++;

  /* signal end of data */
  rowptrs[IDX(N-1,2)+1] = nz;

  return 0;
}



/* Routine to compute the Jacobian matrix from R(y) */
static int ReactionJac(N_Vector y, SUNMatrix Jac, UserData udata)
{
  /* set shortcuts */
  SuperMatrix *Jsuper  = SUNMatrix_SLUNRloc_SuperMatrix(Jac);
  NRformat_loc *Jstore = (NRformat_loc *) Jsuper->Store;
  sunindextype N = udata->N;
  sunindextype *rowptrs = Jstore->rowptr;
  sunindextype *colinds = Jstore->colind;
  realtype *Jdata = (realtype *) Jstore->nzval;
  realtype *Xdata = udata->x;

  /* set local variables */
  sunindextype i, j, nz=0;
  realtype ep = udata->ep;
  realtype ul, uc, ur, vl, vc, vr, wl, wc, wr, xl, xc, xr;
  realtype u1, u2, u3, v1, v2, v3, w1, w2, w3;
  realtype df1, df2, df3, dQdf1, dQdf2, dQdf3;
  realtype ChiL1, ChiL2, ChiL3, ChiR1, ChiR2, ChiR3;
  realtype Ju[9], Jv[9], Jw[9];

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_retval((void *) Ydata, "N_VGetArrayPointer", 0)) return(1);

  /* initialize all local variables to zero (to avoid uninitialized variable warnings) */
  ul = uc = ur = vl = vc = vr = wl = wc = wr = xl = xc = xr = 0.0;
  u1 = u2 = u3 = v1 = v2 = v3 = w1 = w2 = w3 = 0.0;
  df1 = df2 = df3 = dQdf1 = dQdf2 = dQdf3 = 0.0;
  ChiL1 = ChiL2 = ChiL3 = ChiR1 = ChiR2 = ChiR3 = 0.0;

  /* clear out matrix */
  SUNMatZero(Jac);

  /* Dirichlet boundary at left */
  rowptrs[IDX(0,0)] = nz;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(0,0);
  nz++;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(1,0);
  nz++;

  rowptrs[IDX(0,1)] = nz;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(0,1);
  nz++;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(1,1);
  nz++;

  rowptrs[IDX(0,2)] = nz;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(0,2);
  nz++;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(1,2);
  nz++;

  /* iterate over columns, filling in reaction Jacobian */
  for (i=1; i <= N-2; i++) {

    /* set nodal value shortcuts (interval index aligns with left node) */
    xl = Xdata[i-1];
    ul = Ydata[IDX(i-1,0)];
    vl = Ydata[IDX(i-1,1)];
    wl = Ydata[IDX(i-1,2)];
    xc = Xdata[i];
    uc = Ydata[IDX(i,0)];
    vc = Ydata[IDX(i,1)];
    wc = Ydata[IDX(i,2)];
    xr = Xdata[i+1];
    ur = Ydata[IDX(i+1,0)];
    vr = Ydata[IDX(i+1,1)];
    wr = Ydata[IDX(i+1,2)];

    /* compute entries of all Jacobian rows at node i */
    for(j=0; j<9; j++) {
      Ju[j] = ZERO;
      Jv[j] = ZERO;
      Jw[j] = ZERO;
    }

    /* first compute dependence on values to left and center */

    /* evaluate relevant variables in left subinterval */
    u1 = Eval(ul, uc, xl, xc, X1(xl,xc));
    v1 = Eval(vl, vc, xl, xc, X1(xl,xc));
    w1 = Eval(wl, wc, xl, xc, X1(xl,xc));
    u2 = Eval(ul, uc, xl, xc, X2(xl,xc));
    v2 = Eval(vl, vc, xl, xc, X2(xl,xc));
    w2 = Eval(wl, wc, xl, xc, X2(xl,xc));
    u3 = Eval(ul, uc, xl, xc, X3(xl,xc));
    v3 = Eval(vl, vc, xl, xc, X3(xl,xc));
    w3 = Eval(wl, wc, xl, xc, X3(xl,xc));

    dQdf1 = Quad(ONE, ZERO, ZERO, xl, xc);
    dQdf2 = Quad(ZERO, ONE, ZERO, xl, xc);
    dQdf3 = Quad(ZERO, ZERO, ONE, xl, xc);

    ChiL1 = ChiL(xl, xc, X1(xl,xc));
    ChiL2 = ChiL(xl, xc, X2(xl,xc));
    ChiL3 = ChiL(xl, xc, X3(xl,xc));
    ChiR1 = ChiR(xl, xc, X1(xl,xc));
    ChiR2 = ChiR(xl, xc, X2(xl,xc));
    ChiR3 = ChiR(xl, xc, X3(xl,xc));

    /* compute reaction Jacobian components */

    /* R_u = (a - (w+ONE)*u + v*u*u) */
    /*   dR_u/dul   */
    df1 = (-(w1+ONE) + TWO*v1*u1) * ChiL1 * ChiR1;
    df2 = (-(w2+ONE) + TWO*v2*u2) * ChiL2 * ChiR2;
    df3 = (-(w3+ONE) + TWO*v3*u3) * ChiL3 * ChiR3;
    Ju[IDX(0,0)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_u/duc   */
    df1 = (-(w1+ONE) + TWO*v1*u1) * ChiR1 * ChiR1;
    df2 = (-(w2+ONE) + TWO*v2*u2) * ChiR2 * ChiR2;
    df3 = (-(w3+ONE) + TWO*v3*u3) * ChiR3 * ChiR3;
    Ju[IDX(0,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_u/dvl   */
    df1 = (u1*u1) * ChiL1 * ChiR1;
    df2 = (u2*u2) * ChiL2 * ChiR2;
    df3 = (u3*u3) * ChiL3 * ChiR3;
    Ju[IDX(1,0)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_u/dvc   */
    df1 = (u1*u1) * ChiR1 * ChiR1;
    df2 = (u2*u2) * ChiR2 * ChiR2;
    df3 = (u3*u3) * ChiR3 * ChiR3;
    Ju[IDX(1,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_u/dwl   */
    df1 = (-u1) * ChiL1 * ChiR1;
    df2 = (-u2) * ChiL2 * ChiR2;
    df3 = (-u3) * ChiL3 * ChiR3;
    Ju[IDX(2,0)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_u/dwc   */
    df1 = (-u1) * ChiR1 * ChiR1;
    df2 = (-u2) * ChiR2 * ChiR2;
    df3 = (-u3) * ChiR3 * ChiR3;
    Ju[IDX(2,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;


    /* R_v = (w*u - v*u*u) */
    /*   dR_v/dul   */
    df1 = (w1 - TWO*v1*u1) * ChiL1 * ChiR1;
    df2 = (w2 - TWO*v2*u2) * ChiL2 * ChiR2;
    df3 = (w3 - TWO*v3*u3) * ChiL3 * ChiR3;
    Jv[IDX(0,0)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_v/duc   */
    df1 = (w1 - TWO*v1*u1) * ChiR1 * ChiR1;
    df2 = (w2 - TWO*v2*u2) * ChiR2 * ChiR2;
    df3 = (w3 - TWO*v3*u3) * ChiR3 * ChiR3;
    Jv[IDX(0,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_v/dvl   */
    df1 = (-u1*u1) * ChiL1 * ChiR1;
    df2 = (-u2*u2) * ChiL2 * ChiR2;
    df3 = (-u3*u3) * ChiL3 * ChiR3;
    Jv[IDX(1,0)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_v/dvc   */
    df1 = (-u1*u1) * ChiR1 * ChiR1;
    df2 = (-u2*u2) * ChiR2 * ChiR2;
    df3 = (-u3*u3) * ChiR3 * ChiR3;
    Jv[IDX(1,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_v/dwl   */
    df1 = (u1) * ChiL1 * ChiR1;
    df2 = (u2) * ChiL2 * ChiR2;
    df3 = (u3) * ChiL3 * ChiR3;
    Jv[IDX(2,0)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_v/dwc   */
    df1 = (u1) * ChiR1 * ChiR1;
    df2 = (u2) * ChiR2 * ChiR2;
    df3 = (u3) * ChiR3 * ChiR3;
    Jv[IDX(2,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;


    /* R_w = ((b-w)/ep - w*u) */
    /*   dR_w/dul   */
    df1 = (-w1) * ChiL1 * ChiR1;
    df2 = (-w2) * ChiL2 * ChiR2;
    df3 = (-w3) * ChiL3 * ChiR3;
    Jw[IDX(0,0)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_w/duc   */
    df1 = (-w1) * ChiR1 * ChiR1;
    df2 = (-w2) * ChiR2 * ChiR2;
    df3 = (-w3) * ChiR3 * ChiR3;
    Jw[IDX(0,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_w/dwl   */
    df1 = (-ONE/ep - u1) * ChiL1 * ChiR1;
    df2 = (-ONE/ep - u2) * ChiL2 * ChiR2;
    df3 = (-ONE/ep - u3) * ChiL3 * ChiR3;
    Jw[IDX(2,0)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*   dR_w/dwc   */
    df1 = (-ONE/ep - u1) * ChiR1 * ChiR1;
    df2 = (-ONE/ep - u2) * ChiR2 * ChiR2;
    df3 = (-ONE/ep - u3) * ChiR3 * ChiR3;
    Jw[IDX(2,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /* second compute dependence on values to center and right */

    /* evaluate relevant variables in right subinterval */
    u1 = Eval(uc, ur, xc, xr, X1(xc,xr));
    v1 = Eval(vc, vr, xc, xr, X1(xc,xr));
    w1 = Eval(wc, wr, xc, xr, X1(xc,xr));
    u2 = Eval(uc, ur, xc, xr, X2(xc,xr));
    v2 = Eval(vc, vr, xc, xr, X2(xc,xr));
    w2 = Eval(wc, wr, xc, xr, X2(xc,xr));
    u3 = Eval(uc, ur, xc, xr, X3(xc,xr));
    v3 = Eval(vc, vr, xc, xr, X3(xc,xr));
    w3 = Eval(wc, wr, xc, xr, X3(xc,xr));

    dQdf1 = Quad(ONE, ZERO, ZERO, xc, xr);
    dQdf2 = Quad(ZERO, ONE, ZERO, xc, xr);
    dQdf3 = Quad(ZERO, ZERO, ONE, xc, xr);

    ChiL1 = ChiL(xc, xr, X1(xc,xr));
    ChiL2 = ChiL(xc, xr, X2(xc,xr));
    ChiL3 = ChiL(xc, xr, X3(xc,xr));
    ChiR1 = ChiR(xc, xr, X1(xc,xr));
    ChiR2 = ChiR(xc, xr, X2(xc,xr));
    ChiR3 = ChiR(xc, xr, X3(xc,xr));

    /* compute reaction Jacobian components */

    /* R_u = (a - (w+ONE)*u + v*u*u) */
    /*    dR_u/duc    */
    df1 = (-(w1+ONE) + TWO*v1*u1) * ChiL1 * ChiL1;
    df2 = (-(w2+ONE) + TWO*v2*u2) * ChiL2 * ChiL2;
    df3 = (-(w3+ONE) + TWO*v3*u3) * ChiL3 * ChiL3;
    Ju[IDX(0,0)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_u/dur    */
    df1 = (-(w1+ONE) + TWO*v1*u1) * ChiL1 * ChiR1;
    df2 = (-(w2+ONE) + TWO*v2*u2) * ChiL2 * ChiR2;
    df3 = (-(w3+ONE) + TWO*v3*u3) * ChiL3 * ChiR3;
    Ju[IDX(0,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_u/dvc    */
    df1 = (u1*u1) * ChiL1 * ChiL1;
    df2 = (u2*u2) * ChiL2 * ChiL2;
    df3 = (u3*u3) * ChiL3 * ChiL3;
    Ju[IDX(1,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_u/dvr    */
    df1 = (u1*u1) * ChiL1 * ChiR1;
    df2 = (u2*u2) * ChiL2 * ChiR2;
    df3 = (u3*u3) * ChiL3 * ChiR3;
    Ju[IDX(1,2)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_u/dwc   */
    df1 = (-u1) * ChiL1 * ChiL1;
    df2 = (-u2) * ChiL2 * ChiL2;
    df3 = (-u3) * ChiL3 * ChiL3;
    Ju[IDX(2,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_u/dwr   */
    df1 = (-u1) * ChiL1 * ChiR1;
    df2 = (-u2) * ChiL2 * ChiR2;
    df3 = (-u3) * ChiL3 * ChiR3;
    Ju[IDX(2,2)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;


    /* R_v = (w*u - v*u*u) */
    /*    dR_v/duc     */
    df1 = (w1 - TWO*v1*u1) * ChiL1 * ChiL1;
    df2 = (w2 - TWO*v2*u2) * ChiL2 * ChiL2;
    df3 = (w3 - TWO*v3*u3) * ChiL3 * ChiL3;
    Jv[IDX(0,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_v/dur     */
    df1 = (w1 - TWO*v1*u1) * ChiL1 * ChiR1;
    df2 = (w2 - TWO*v2*u2) * ChiL2 * ChiR2;
    df3 = (w3 - TWO*v3*u3) * ChiL3 * ChiR3;
    Jv[IDX(0,2)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_v/dvc     */
    df1 = (-u1*u1) * ChiL1 * ChiL1;
    df2 = (-u2*u2) * ChiL2 * ChiL2;
    df3 = (-u3*u3) * ChiL3 * ChiL3;
    Jv[IDX(1,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_v/dvr    */
    df1 = (-u1*u1) * ChiL1 * ChiR1;
    df2 = (-u2*u2) * ChiL2 * ChiR2;
    df3 = (-u3*u3) * ChiL3 * ChiR3;
    Jv[IDX(1,2)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_v/dwc    */
    df1 = (u1) * ChiL1 * ChiL1;
    df2 = (u2) * ChiL2 * ChiL2;
    df3 = (u3) * ChiL3 * ChiL3;
    Jv[IDX(2,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_v/dwr    */
    df1 = (u1) * ChiL1 * ChiR1;
    df2 = (u2) * ChiL2 * ChiR2;
    df3 = (u3) * ChiL3 * ChiR3;
    Jv[IDX(2,2)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;


    /* R_w = ((b-w)/ep - w*u) */
    /*    dR_w/duc    */
    df1 = (-w1) * ChiL1 * ChiL1;
    df2 = (-w2) * ChiL2 * ChiL2;
    df3 = (-w3) * ChiL3 * ChiL3;
    Jw[IDX(0,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_w/dur    */
    df1 = (-w1) * ChiL1 * ChiR1;
    df2 = (-w2) * ChiL2 * ChiR2;
    df3 = (-w3) * ChiL3 * ChiR3;
    Jw[IDX(0,2)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_w/dwc    */
    df1 = (-ONE/ep - u1) * ChiL1 * ChiL1;
    df2 = (-ONE/ep - u2) * ChiL2 * ChiL2;
    df3 = (-ONE/ep - u3) * ChiL3 * ChiL3;
    Jw[IDX(2,1)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /*    dR_w/dwr    */
    df1 = (-ONE/ep - u1) * ChiL1 * ChiR1;
    df2 = (-ONE/ep - u2) * ChiL2 * ChiR2;
    df3 = (-ONE/ep - u3) * ChiL3 * ChiR3;
    Jw[IDX(2,2)] += dQdf1*df1 + dQdf2*df2 + dQdf3*df3;

    /* insert Jacobian entries into CSR matrix structure */

    /* Ju row */
    rowptrs[IDX(i,0)] = nz;

    Jdata[nz] = Ju[IDX(0,0)];
    Jdata[nz+1] = Ju[IDX(1,0)];
    Jdata[nz+2] = Ju[IDX(2,0)];
    colinds[nz] = IDX(i-1,0);
    colinds[nz+1] = IDX(i-1,1);
    colinds[nz+2] = IDX(i-1,2);
    nz += 3;

    Jdata[nz] = Ju[IDX(0,1)];
    Jdata[nz+1] = Ju[IDX(1,1)];
    Jdata[nz+2] = Ju[IDX(2,1)];
    colinds[nz] = IDX(i,0);
    colinds[nz+1] = IDX(i,1);
    colinds[nz+2] = IDX(i,2);
    nz += 3;

    Jdata[nz] = Ju[IDX(0,2)];
    Jdata[nz+1] = Ju[IDX(1,2)];
    Jdata[nz+2] = Ju[IDX(2,2)];
    colinds[nz] = IDX(i+1,0);
    colinds[nz+1] = IDX(i+1,1);
    colinds[nz+2] = IDX(i+1,2);
    nz += 3;

    /* Jv row */
    rowptrs[IDX(i,1)] = nz;

    Jdata[nz] = Jv[IDX(0,0)];
    Jdata[nz+1] = Jv[IDX(1,0)];
    Jdata[nz+2] = Jv[IDX(2,0)];
    colinds[nz] = IDX(i-1,0);
    colinds[nz+1] = IDX(i-1,1);
    colinds[nz+2] = IDX(i-1,2);
    nz += 3;

    Jdata[nz] = Jv[IDX(0,1)];
    Jdata[nz+1] = Jv[IDX(1,1)];
    Jdata[nz+2] = Jv[IDX(2,1)];
    colinds[nz] = IDX(i,0);
    colinds[nz+1] = IDX(i,1);
    colinds[nz+2] = IDX(i,2);
    nz += 3;

    Jdata[nz] = Jv[IDX(0,2)];
    Jdata[nz+1] = Jv[IDX(1,2)];
    Jdata[nz+2] = Jv[IDX(2,2)];
    colinds[nz] = IDX(i+1,0);
    colinds[nz+1] = IDX(i+1,1);
    colinds[nz+2] = IDX(i+1,2);
    nz += 3;

    /* Jw row */
    rowptrs[IDX(i,2)] = nz;

    Jdata[nz] = Jw[IDX(0,0)];
    Jdata[nz+1] = Jw[IDX(1,0)];
    Jdata[nz+2] = Jw[IDX(2,0)];
    colinds[nz] = IDX(i-1,0);
    colinds[nz+1] = IDX(i-1,1);
    colinds[nz+2] = IDX(i-1,2);
    nz += 3;

    Jdata[nz] = Jw[IDX(0,1)];
    Jdata[nz+1] = Jw[IDX(1,1)];
    Jdata[nz+2] = Jw[IDX(2,1)];
    colinds[nz] = IDX(i,0);
    colinds[nz+1] = IDX(i,1);
    colinds[nz+2] = IDX(i,2);
    nz += 3;

    Jdata[nz] = Jw[IDX(0,2)];
    Jdata[nz+1] = Jw[IDX(1,2)];
    Jdata[nz+2] = Jw[IDX(2,2)];
    colinds[nz] = IDX(i+1,0);
    colinds[nz+1] = IDX(i+1,1);
    colinds[nz+2] = IDX(i+1,2);
    nz += 3;

  }

  /* Dirichlet boundary at right */
  rowptrs[IDX(N-1,0)] = nz;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(N-2,0);
  nz++;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(N-1,0);
  nz++;

  rowptrs[IDX(N-1,1)] = nz;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(N-2,1);
  nz++;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(N-1,1);
  nz++;

  rowptrs[IDX(N-1,2)] = nz;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(N-2,2);
  nz++;

  Jdata[nz]   = ZERO;
  colinds[nz] = IDX(N-1,2);
  nz++;

  /* signal end of data */
  rowptrs[IDX(N-1,2)+1] = nz;

  return 0;
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a code so check if
             it is < 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) returnvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


/*---- end of file ----*/
