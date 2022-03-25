/* -----------------------------------------------------------------------------
 * Programmer(s): Rujeko Chinomona @SMU and @LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Example problem:
 *
 * The following test simulates a brusselator problem from chemical
 * kinetics.  This is n PDE system with 3 components, Y = [u,v,w],
 * satisfying the equations,
 *    u_t = du*u_xx - au*u_x +  a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx - av*v_x +  w*u - v*u^2
 *    w_t = dw*w_xx - aw*w_x + (b-w)/ep - w*u
 * for t in [0, 10], x in [0, 1], with initial conditions
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
 * We use parameters:
 * du = dv = dw = 0.01 (diffusion coefficients)
 * au = av = aw = -0.001 (advection coefficients - velocity)
 * a  = 0.6
 * b  = 2
 * ep = 0.01
 *
 * The spatial derivatives are computed using second-order
 * centered differences, with the data distributed over N points
 * on a uniform spatial grid.
 * Note: larger values of advection require advection schemes such as
 * upwinding not implemented here.
 *
 * This program solves the problem with multiple solvers listed below.
 * We select method to used based on solve_type input:
 * 0. MIS with third order dirk inner
 * 1. 5th order dirk method for reference solution
 * 2. MRI-GARK34a with erk inner
 * 3. MRI-GARK34a with dirk inner
 * 4. IMEX-MRI3b with erk inner
 * 5. IMEX-MRI3b with dirk inner
 * 6. IMEX-MRI4 with erk inner
 * 7. IMEX-MRI4 with dirk inner
 *
 *  We use Newton iteration with the SUNBAND linear solver and a user supplied
 * Jacobian routine for nonlinear solves.
 *
 * This program solves the problem with the MRI stepper. 10 outputs are printed
 * at equal intervals, and run statistics are printed at the end.
 * ---------------------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <arkode/arkode_mristep.h>    /* prototypes for MRIStep fcts., consts */
#include <arkode/arkode_arkstep.h>    /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_serial.h>   /* access to Serial N_Vector */
#include <sunmatrix/sunmatrix_band.h> /* access to band SUNMatrix */
#include <sunlinsol/sunlinsol_band.h> /* access to band SUNLinearSolver */
#include <sundials/sundials_types.h>  /* def. of type 'realtype' */
#include <sundials/sundials_math.h>   /* def. of SUNRsqrt, etc. */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

/* Define some constants */
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

/* accessor macros between (x,v) location and 1D NVector array */
#define IDX(x,v) (3*(x)+v)

/* user data structure */
typedef struct {
  sunindextype N;  /* number of intervals     */
  realtype dx;     /* mesh spacing            */
  realtype a;      /* constant forcing on u   */
  realtype b;      /* steady-state value of w */
  realtype pi;     /* value of pi             */
  realtype du;     /* diffusion coeff for u   */
  realtype dv;     /* diffusion coeff for v   */
  realtype dw;     /* diffusion coeff for w   */
  realtype au;     /* advection coeff for u   */
  realtype av;     /* advection coeff for v   */
  realtype aw;     /* advection coeff for w   */
  realtype ep;     /* stiffness parameter     */
} *UserData;

/* User-supplied Functions Called by the Solver */
static int ff(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fse(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fsi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fs(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jf(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jsi(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Js(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private helper functions  */
static int SetIC(N_Vector y, void *user_data);
static int AdvectionJac(realtype c, SUNMatrix Jac, UserData udata);
static int LaplaceMatrix(realtype c, SUNMatrix Jac, UserData udata);
static int ReactionJac(realtype c, N_Vector y, SUNMatrix Jac, UserData udata);

/* Private function to check function return values */
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Main Program */
int main(int argc, char *argv[])
{
  /* general problem parameters */
  realtype T0 = ZERO;    /* initial time                 */
  realtype Tf = RCONST(10.0);   /* final time                   */
  int Nt = 10;                  /* total number of output times */
  realtype dTout = (Tf-T0)/Nt;  /* time between outputs         */
  int Nvar = 3;                 /* number of solution fields    */
  sunindextype N = 101;         /* spatial mesh size            */
  realtype hs;                  /* slow step size       */
  realtype m = RCONST(10.0);    /* time-scale separation factor */
  int solve_type;               /* solver configuration */
  realtype dx = ONE/(N-1);      /* set spatial mesh spacing     */
  realtype a = 0.6;             /* problem parameters           */
  realtype b = 2.0;
  realtype pi = RCONST(4.0)*atan(ONE);
  realtype du = 0.01;
  realtype dv = 0.01;
  realtype dw = 0.01;
  realtype au = -0.001;
  realtype av = -0.001;
  realtype aw = -0.001;
  realtype ep = 1.0e-2;         /* stiffness parameter          */
  realtype reltol = 1.0e-12;    /* tolerances                   */
  realtype abstol = 1.0e-14;

  /* general problem variables */
  int retval;                               /* reusable return flag          */
  N_Vector y = NULL;                        /* empty solution vector         */
  void *arkode_mem = NULL;                  /* empty ARKode memory structure */
  void *inner_arkode_mem = NULL;            /* empty ARKode memory structure */
  MRIStepInnerStepper inner_stepper = NULL; /* inner stepper                 */
  ARKodeButcherTable B = NULL;              /* fast method Butcher table     */
  MRIStepCoupling C = NULL;                 /* slow coupling coefficients    */
  SUNMatrix Af = NULL;                      /* matrix for fast solver        */
  SUNLinearSolver LSf = NULL;               /* fast linear solver object     */
  SUNMatrix As = NULL;                      /* matrix for slow solver        */
  SUNLinearSolver LSs = NULL;               /* slow linear solver object     */
  booleantype implicit_slow;
  booleantype imex_slow = SUNFALSE;
  N_Vector umask = NULL;                    /* empty mask vectors            */
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  realtype t, tout;                         /* current/output time data      */
  realtype hf;                              /* fast time step                */
  realtype u, v, w;                         /* temp data values              */
  FILE *FID, *UFID, *VFID, *WFID;           /* output file pointers          */
  int iout;                                 /* output counter                */
  long int nsts, nstf;                      /* step stats                    */
  long int nfse, nfsi, nffe, nffi;          /* RHS stats                     */
  long int nnif, nncf, njef, nnis, nncs, njes;
  sunindextype NEQ;                         /* number of equations           */
  sunindextype i;                           /* counter                       */
  UserData udata = NULL;                    /* user data pointer             */
  realtype* data = NULL;                    /* array for vector data         */
  realtype  gamma, beta;
  char ofname[50]="";

  /* Create the SUNDIALS context object for this simulation. */
  SUNContext ctx = NULL;
  SUNContext_Create(NULL, &ctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return 1;

  /*
   * Initialization
   */

  /* Retrieve the command-line options: solve_type h  */
  if (argc < 2) {
    printf("ERROR: enter solve_type and hs \n");
    return(-1);
  }
  solve_type = (sunindextype) atol(argv[1]);
  hs = (realtype) atof(argv[2]);


  /* Check arguments for validity */
  /*   0 <= solve_type <= 7       */
  /*   h > 0                      */

  if ((solve_type < 0) || (solve_type > 7)) {
    printf("ERROR: solve_type be an integer in [0,7] \n");
    return(-1);
  }
  implicit_slow = SUNFALSE;
  if (solve_type > 1)
    implicit_slow = SUNTRUE;
  if (solve_type > 3)
    imex_slow = SUNTRUE;
  if (hs <= ZERO){
    printf("ERROR: hs must be in positive\n");
    return(-1);
  }
  hf = hs/m;
  NEQ = N*Nvar;


  /* Initial problem output */
  printf("\n1D Advection-Diffusion-Reaction (Brusselator) test problem:\n");
  printf("    time domain:  (%"GSYM",%"GSYM"]\n",T0,Tf);
  printf("    hs = %"GSYM"\n",hs);
  printf("    hf = %"GSYM"\n",hf);
  printf("    m  = %"GSYM"\n",m);
  printf("    N  = %li,  NEQ = %li\n", (long int) N, (long int) NEQ);
  printf("    dx = %"GSYM"\n",dx);
  printf("    problem parameters:  a = %"GSYM",  b = %"GSYM",  ep = %"GSYM"\n",
         a, b, ep);
  printf("    diffusion coefficients:  du = %"GSYM",  dv = %"GSYM",  dw = %"GSYM"\n",
         du, dv, dw);
  printf("    advection coefficients:  au = %"GSYM",  av = %"GSYM",  aw = %"GSYM"\n",
         au, av, aw);


  switch (solve_type) {
  case(0):
    /* reltol = SUNMAX(hs*hs*hs, 1e-10); */
    /* abstol = 1e-11; */
    printf("    solver: exp-3/dirk-3 (MIS / ESDIRK-3-3)\n\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n\n", reltol, abstol);
    break;
  case(1):
    reltol = SUNMAX(hs*hs*hs*hs*hs, 1e-14);
    abstol = 1e-14;
    printf("    solver: none/dirk-5 (no slow, default 5th order dirk fast)\n\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n\n", reltol, abstol);
    break;
  case(2):
    /* reltol = SUNMAX(hs*hs*hs, 1e-10); */
    /* abstol = 1e-11; */
    printf("    solver: dirk-3/exp-3 (MRI-GARK-ESDIRK34a / ERK-3-3) -- solve decoupled\n\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n\n", reltol, abstol);
    break;
  case(3):
    /* reltol = SUNMAX(hs*hs*hs, 1e-10); */
    /* abstol = 1e-11; */
    printf("    solver: dirk-3/dirk-3 (MRI-GARK-ESDIRK34a / ESDIRK-3-3) -- solve decoupled\n\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n\n", reltol, abstol);
    break;
  case(4):
    /* reltol = SUNMAX(hs*hs*hs, 1e-14); */
    /* abstol = 1e-14; */
    printf("    solver: ars343/exp-3 (IMEX-MRI3b / ERK-3-3) -- solve decoupled\n\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n\n", reltol, abstol);
    break;
  case(5):
    /* reltol = SUNMAX(hs*hs*hs, 1e-14); */
    /* abstol = 1e-14; */
    printf("    solver: ars343/dirk-3 (IMEX-MRI3b / ESDIRK-3-3) -- solve decoupled\n\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n\n", reltol, abstol);
    break;
  case(6):
    /* reltol = SUNMAX(hs*hs*hs*hs, 1e-14); */
    /* abstol = 1e-14; */
    printf("    solver: imexark4/exp-4 (IMEX-MRI4 / ERK-4-4) -- solve decoupled\n\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n\n", reltol, abstol);
    break;
  case(7):
    /* reltol = SUNMAX(hs*hs*hs*hs, 1e-14); */
    /* abstol = 1e-14; */
    printf("    solver: imexark4/dirk-4 (IMEX-MRI4 / CASH(5,3,4)-DIRK ) -- solve decoupled\n\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n\n", reltol, abstol);
    break;
  }

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  if (check_retval((void *) udata, "malloc", 2)) return 1;

  /* store the inputs in the UserData structure */
  udata->N  = N;
  udata->a  = a;
  udata->b  = b;
  udata->du = du;
  udata->dv = dv;
  udata->dw = dw;
  udata->au = au;
  udata->av = av;
  udata->aw = aw;
  udata->ep = ep;
  udata->pi = pi;
  udata->dx = dx;

  /* Create solution vector */
  y = N_VNew_Serial(NEQ, ctx);           /* Create vector for solution */
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return 1;

  /* Set initial condition */
  retval = SetIC(y, udata);
  if (check_retval(&retval, "SetIC", 1)) return 1;

  /* Create vector masks  */
  umask = N_VClone(y);
  if (check_retval((void *)umask, "N_VNew_Serial", 0)) return 1;

  vmask = N_VClone(y);
  if (check_retval((void *)vmask, "N_VNew_Serial", 0)) return 1;

  wmask = N_VClone(y);
  if (check_retval((void *)wmask, "N_VNew_Serial", 0)) return 1;

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  data = N_VGetArrayPointer(umask);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,0)] = RCONST(1.0);

  N_VConst(0.0, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,1)] = RCONST(1.0);

  N_VConst(0.0, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,2)] = RCONST(1.0);

  /*
   * Create the fast integrator and set options
   */

  /* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the inital time T0,
     and the initial dependent variable vector y. */
  switch (solve_type) {
  case(0):
  case(3):        /* esdirk-3-3 fast solver */
  case(5):
    inner_arkode_mem = ARKStepCreate(NULL, ff, T0, y, ctx);
    if (check_retval((void *) inner_arkode_mem, "ARKStepCreate", 0)) return 1;
    B = ARKodeButcherTable_Alloc(3, SUNFALSE);
    if (check_retval((void *)B, "ARKodeButcherTable_Alloc", 0)) return 1;
    beta  = SUNRsqrt(RCONST(3.0))/RCONST(6.0) + RCONST(0.5);
    gamma = (-ONE/RCONST(8.0))*(SUNRsqrt(RCONST(3.0))+ONE);
    B->A[1][0] = RCONST(4.0)*gamma+TWO*beta;
    B->A[1][1] = ONE-RCONST(4.0)*gamma-TWO*beta;
    B->A[2][0] = RCONST(0.5)-beta-gamma;
    B->A[2][1] = gamma;
    B->A[2][2] = beta;
    B->b[0] = ONE/RCONST(6.0);
    B->b[1] = ONE/RCONST(6.0);
    B->b[2] = TWO/RCONST(3.0);
    B->c[1] = ONE;
    B->c[2] = RCONST(0.5);
    B->q=3;
    retval = ARKStepSetTables(inner_arkode_mem, 3, 0, B, NULL);
    if (check_retval(&retval, "ARKStepSetTables", 1)) return 1;

    /* Initialize matrix and linear solver data structures */
    Af = SUNBandMatrix(NEQ, 4, 4, ctx);
    if (check_retval((void *)Af, "SUNBandMatrix", 0)) return 1;

    LSf = SUNLinSol_Band(y, Af, ctx);
    if (check_retval((void *)LSf, "SUNLinSol_Band", 0)) return 1;

    /* Specify fast tolerances */
    retval = ARKStepSStolerances(inner_arkode_mem, reltol, abstol);
    if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;

    /* Attach matrix and linear solver */
    retval = ARKStepSetLinearSolver(inner_arkode_mem, LSf, Af);
    if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return 1;

    /* Set max number of nonlinear iters */
    retval = ARKStepSetMaxNonlinIters(inner_arkode_mem, 10);
    if (check_retval(&retval, "ARKStepSetMaxNonlinIters", 1)) return 1;

    /* Set the Jacobian routine */
    retval = ARKStepSetJacFn(inner_arkode_mem, Jf);
    if (check_retval(&retval, "ARKStepSetJacFn", 1)) return 1;
    break;
  case(1):        /*dirk 5th order fast solver (full problem) */
    inner_arkode_mem = ARKStepCreate(NULL, f, T0, y, ctx);
    if (check_retval((void *) inner_arkode_mem, "ARKStepCreate", 0)) return 1;

    /* Set method order to use */
    retval = ARKStepSetOrder(inner_arkode_mem, 5);
    if (check_retval(&retval, "ARKStepSetOrder",1)) return 1;

    /* Initialize matrix and linear solver data structures */
    Af = SUNBandMatrix(NEQ, 4, 4, ctx);
    if (check_retval((void *)Af, "SUNBandMatrix", 0)) return 1;

    LSf = SUNLinSol_Band(y, Af, ctx);
    if (check_retval((void *)LSf, "SUNLinSol_Band", 0)) return 1;

    /* Specify fast tolerances */
    retval = ARKStepSStolerances(inner_arkode_mem, reltol, abstol);
    if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;

    /* Attach matrix and linear solver */
    retval = ARKStepSetLinearSolver(inner_arkode_mem, LSf, Af);
    if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return 1;

    /* Set the Jacobian routine */
    retval = ARKStepSetJacFn(inner_arkode_mem, Jac);
    if (check_retval(&retval, "ARKStepSetJacFn", 1)) return 1;
    break;
  case(2):           /* erk-3-3 fast solver */
  case(4):
    inner_arkode_mem = ARKStepCreate(ff, NULL, T0, y, ctx);
    if (check_retval((void *) inner_arkode_mem, "ARKStepCreate", 0)) return 1;
    B = ARKodeButcherTable_Alloc(3, SUNTRUE);
    if (check_retval((void *)B, "ARKodeButcherTable_Alloc", 0)) return 1;
    B->A[1][0] = RCONST(0.5);
    B->A[2][0] = -ONE;
    B->A[2][1] = TWO;
    B->b[0] = ONE/RCONST(6.0);
    B->b[1] = TWO/RCONST(3.0);
    B->b[2] = ONE/RCONST(6.0);
    B->d[1] = ONE;
    B->c[1] = RCONST(0.5);
    B->c[2] = ONE;
    B->q=3;
    B->p=2;
    retval = ARKStepSetTables(inner_arkode_mem, 3, 2, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) return 1;
    break;
  case(6):                 /* erk-4-4 fast solver */
    inner_arkode_mem = ARKStepCreate(ff, NULL, T0, y, ctx);
    if (check_retval((void *) inner_arkode_mem, "ARKStepCreate", 0)) return 1;
    B = ARKodeButcherTable_Alloc(4, SUNFALSE);
    if (check_retval((void *)B, "ARKodeButcherTable_Alloc", 0)) return 1;
    B->A[1][0] = RCONST(0.5);
    B->A[2][1] = RCONST(0.5);
    B->A[3][2] = ONE;
    B->b[0] = ONE/RCONST(6.0);
    B->b[1] = ONE/RCONST(3.0);
    B->b[2] = ONE/RCONST(3.0);
    B->b[3] = ONE/RCONST(6.0);
    B->c[1] = RCONST(0.5);
    B->c[2] = RCONST(0.5);
    B->c[3] = ONE;
    B->q=4;
    retval = ARKStepSetTables(inner_arkode_mem, 4, 0, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) return 1;
    break;
  case(7):                /* Cash(5,3,4)-SDIRK fast solver */
    inner_arkode_mem = ARKStepCreate(NULL, ff, T0, y, ctx);
    if (check_retval((void *) inner_arkode_mem, "ARKStepCreate", 0)) return 1;

    /* Set fast method */
    retval = ARKStepSetTableNum(inner_arkode_mem, ARKODE_CASH_5_3_4, -1);
    if (check_retval(&retval, "ARKStepSetTableNum", 1)) return 1;

    /* Initialize matrix and linear solver data structures */
    Af = SUNBandMatrix(NEQ, 4, 4, ctx);
    if (check_retval((void *)Af, "SUNBandMatrix", 0)) return 1;

    LSf = SUNLinSol_Band(y, Af, ctx);
    if (check_retval((void *)LSf, "SUNLinSol_Band", 0)) return 1;

    /* Specify fast tolerances */
    retval = ARKStepSStolerances(inner_arkode_mem, reltol, abstol);
    if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;

    /* Attach matrix and linear solver */
    retval = ARKStepSetLinearSolver(inner_arkode_mem, LSf, Af);
    if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return 1;

    /* Set max number of nonlinear iters */
    retval = ARKStepSetMaxNonlinIters(inner_arkode_mem, 10);
    if (check_retval(&retval, "ARKStepSetMaxNonlinIters", 1)) return 1;

    /* Set the Jacobian routine */
    retval = ARKStepSetJacFn(inner_arkode_mem, Jf);
    if (check_retval(&retval, "ARKStepSetJacFn", 1)) return 1;
    break;
 }

  /* Attach user data to fast integrator */
  retval = ARKStepSetUserData(inner_arkode_mem, (void *) udata);
  if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;

  /* Set the fast step size */
  retval = ARKStepSetFixedStep(inner_arkode_mem, hf);
  if (check_retval(&retval, "ARKStepSetFixedStep", 1)) return 1;

  /* Create inner stepper */
  retval = ARKStepCreateMRIStepInnerStepper(inner_arkode_mem,
                                            &inner_stepper);
  if (check_retval(&retval, "ARKStepCreateMRIStepInnerStepper", 1)) return 1;

  /*
   * Create the slow integrator and set options
   */

  /* Initialize the slow integrator. Specify the slow right-hand side
     function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the inital time
     T0, the initial dependent variable vector y, and the fast integrator. */
  switch (solve_type) {
  case(0):                    /* use MIS outer integrator default for MRIStep */
    arkode_mem = MRIStepCreate(fs, NULL, T0, y, inner_stepper, ctx);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;
    break;
  case(1):                    /* no slow dynamics (use ERK-2-2) */
    arkode_mem = MRIStepCreate(f0, NULL, T0, y, inner_stepper, ctx);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;
    B = ARKodeButcherTable_Alloc(2, SUNFALSE);
    if (check_retval((void *)B, "ARKodeButcherTable_Alloc", 0)) return 1;
    B->A[1][0] = TWO/RCONST(3.0);
    B->b[0] = RCONST(0.25);
    B->b[1] = RCONST(0.75);
    B->c[1] = TWO/RCONST(3.0);
    B->q=2;
    C = MRIStepCoupling_MIStoMRI(B, 2, 0);
    if (check_retval((void *)C, "MRIStepCoupling_MIStoMRI", 0)) return 1;
    retval = MRIStepSetCoupling(arkode_mem, C);
    if (check_retval(&retval, "MRIStepSetCoupling", 1)) return 1;
    break;
  case(2):
  case(3): /* MRI-GARK-ESDIRK34a, solve-decoupled slow solver */
    arkode_mem = MRIStepCreate(NULL, fs, T0, y, inner_stepper, ctx);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;

    C = MRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ESDIRK34a);
    if (check_retval((void *)C, "MRIStepCoupling_LoadTable", 0)) return 1;

    retval = MRIStepSetCoupling(arkode_mem, C);
    if (check_retval(&retval, "MRIStepSetCoupling", 1)) return 1;

    /* Initialize matrix and linear solver data structures */
    As = SUNBandMatrix(NEQ, 4, 4, ctx);
    if (check_retval((void *)As, "SUNBandMatrix", 0)) return 1;

    LSs = SUNLinSol_Band(y, As, ctx);
    if (check_retval((void *)LSs, "SUNLinSol_Band", 0)) return 1;

    /* Specify tolerances */
    retval = MRIStepSStolerances(arkode_mem, reltol, abstol);
    if (check_retval(&retval, "MRIStepSStolerances", 1)) return 1;

    /* Attach matrix and linear solver */
    retval = MRIStepSetLinearSolver(arkode_mem, LSs, As);
    if (check_retval(&retval, "MRIStepSetLinearSolver", 1)) return 1;

    /* Set the Jacobian routine */
    retval = MRIStepSetJacFn(arkode_mem, Js);
    if (check_retval(&retval, "MRIStepSetJacFn", 1)) return 1;
    break;
  case(4):
  case(5):        /* IMEX-MRI-GARK3b, solve-decoupled slow solver */
    arkode_mem = MRIStepCreate(fse, fsi, T0, y, inner_stepper, ctx);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;

    C = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK3b);
    if (check_retval((void *)C, "MRIStepCoupling_LoadTable", 0)) return 1;

    retval = MRIStepSetCoupling(arkode_mem, C);
    if (check_retval(&retval, "MRIStepSetCoupling", 1)) return 1;

    /* Initialize matrix and linear solver data structures */
    As = SUNBandMatrix(NEQ, 4, 4, ctx);
    if (check_retval((void *)As, "SUNBandMatrix", 0)) return 1;

    LSs = SUNLinSol_Band(y, As, ctx);
    if (check_retval((void *)LSs, "SUNLinSol_Band", 0)) return 1;

    /* Specify tolerances */
    retval = MRIStepSStolerances(arkode_mem, reltol, abstol);
    if (check_retval(&retval, "MRIStepSStolerances", 1)) return 1;

    /* Attach matrix and linear solver */
    retval = MRIStepSetLinearSolver(arkode_mem, LSs, As);
    if (check_retval(&retval, "MRIStepSetLinearSolver", 1)) return 1;

    /* Set the Jacobian routine */
    retval = MRIStepSetJacFn(arkode_mem, Jsi);
    if (check_retval(&retval, "MRIStepSetJacFn", 1)) return 1;
    break;
  case(6):
  case(7):     /* IMEX-MRI-GARK4, solve-decoupled slow solver */
    arkode_mem = MRIStepCreate(fse, fsi, T0, y, inner_stepper, ctx);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;

    C = MRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK4);
    if (check_retval((void *)C, "MRIStepCoupling_LoadTable", 0)) return 1;

    retval = MRIStepSetCoupling(arkode_mem, C);
    if (check_retval(&retval, "MRIStepSetCoupling", 1)) return 1;

    /* Initialize matrix and linear solver data structures */
    As = SUNBandMatrix(NEQ, 4, 4, ctx);
    if (check_retval((void *)As, "SUNBandMatrix", 0)) return 1;

    LSs = SUNLinSol_Band(y, As, ctx);
    if (check_retval((void *)LSs, "SUNLinSol_Band", 0)) return 1;

    /* Specify tolerances */
    retval = MRIStepSStolerances(arkode_mem, reltol, abstol);
    if (check_retval(&retval, "MRIStepSStolerances", 1)) return 1;

    /* Attach matrix and linear solver */
    retval = MRIStepSetLinearSolver(arkode_mem, LSs, As);
    if (check_retval(&retval, "MRIStepSetLinearSolver", 1)) return 1;

    /* Set the Jacobian routine */
    retval = MRIStepSetJacFn(arkode_mem, Jsi);
    if (check_retval(&retval, "MRIStepSetJacFn", 1)) return 1;
    break;
  }

  /* Pass udata to user functions */
  retval = MRIStepSetUserData(arkode_mem, (void *) udata);
  if (check_retval(&retval, "MRIStepSetUserData", 1)) return 1;

  /* Set the slow step size */
  retval = MRIStepSetFixedStep(arkode_mem, hs);
  if (check_retval(&retval, "MRIStepSetFixedStep", 1)) return 1;

  /* Set maximum number of steps taken by solver */
  retval = MRIStepSetMaxNumSteps(arkode_mem, 1000000);
  if (check_retval(&retval, "MRIStepSetMaxNumSteps", 1)) return 1;

  /*
   * Integrate ODE
   */

  /* output spatial mesh to disk */
  FID=fopen("bruss1D_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16"ESYM"\n", udata->dx*i);
  fclose(FID);

  /* Open output stream for results, access data arrays */
  strcpy(ofname, "bruss1D_");
  strcat(ofname, "u_");
  strcat(ofname, argv[1]);
  strcat(ofname, "_");
  strcat(ofname, argv[2]);
  strcat(ofname, ".txt");
  UFID=fopen(ofname,"w");

  strcpy(ofname, "bruss1D_");
  strcat(ofname, "v_");
  strcat(ofname, argv[1]);
  strcat(ofname, "_");
  strcat(ofname, argv[2]);
  strcat(ofname, ".txt");
  VFID=fopen(ofname,"w");

  strcpy(ofname, "bruss1D_");
  strcat(ofname, "w_");
  strcat(ofname, argv[1]);
  strcat(ofname, "_");
  strcat(ofname, argv[2]);
  strcat(ofname, ".txt");
  WFID=fopen(ofname,"w");

  /* output initial condition to disk */
  data = N_VGetArrayPointer(y);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return 1;

  for (i=0; i<N; i++)  fprintf(UFID," %.16"ESYM, data[IDX(i,0)]);
  for (i=0; i<N; i++)  fprintf(VFID," %.16"ESYM, data[IDX(i,1)]);
  for (i=0; i<N; i++)  fprintf(WFID," %.16"ESYM, data[IDX(i,2)]);
  fprintf(UFID,"\n");
  fprintf(VFID,"\n");
  fprintf(WFID,"\n");

  /* Main time-stepping loop: calls MRIStepEvolve to perform the integration, then
     prints results.  Stops when the final time has been reached */
  t = T0;
  tout = T0+dTout;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms\n");
  printf("   ----------------------------------------------\n");
  for (iout=0; iout<Nt; iout++) {

    /* call integrator */
    retval = MRIStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "MRIStepEvolve", 1)) break;

    /* access/print solution statistics */
    u = N_VWL2Norm(y,umask);
    u = SUNRsqrt(u*u/N);
    v = N_VWL2Norm(y,vmask);
    v = SUNRsqrt(v*v/N);
    w = N_VWL2Norm(y,wmask);
    w = SUNRsqrt(w*w/N);
    printf("  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"\n", t, u, v, w);

    /* successful solve: update output time */
    tout += dTout;
    tout = (tout > Tf) ? Tf : tout;

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

  /*
   * Finalize
   */

  /* Get some slow integrator statistics */
  retval = MRIStepGetNumSteps(arkode_mem, &nsts);
  check_retval(&retval, "MRIStepGetNumSteps", 1);
  retval = MRIStepGetNumRhsEvals(arkode_mem, &nfse, &nfsi);
  check_retval(&retval, "MRIStepGetNumRhsEvals", 1);

  /* Get some fast integrator statistics */
  retval = ARKStepGetNumSteps(inner_arkode_mem, &nstf);
  check_retval(&retval, "ARKStepGetNumSteps", 1);
  retval = ARKStepGetNumRhsEvals(inner_arkode_mem, &nffe, &nffi);
  check_retval(&retval, "ARKStepGetNumRhsEvals", 1);

  /* Print some final statistics */
  printf("\nFinal Solver Statistics:\n");
  printf("   Slow Steps: nsts = %li\n", nsts);
  printf("   Fast Steps: nstf = %li\n", nstf);
  if (imex_slow) {
    if ((solve_type==0) || (solve_type==1) || (solve_type==3) || (solve_type==5) || (solve_type==7)) {
      printf("   Total RHS evals:  Fse = %li, Fsi = %li,  Ff = %li\n", nfse, nfsi, nffi);
    } else {
      printf("   Total RHS evals:  Fse = %li, Fsi = %li,  Ff = %li\n", nfse, nfsi, nffe);
    }
  } else if (implicit_slow) {
    if ((solve_type==0) || (solve_type==1) || (solve_type==3) || (solve_type==5) || (solve_type==7)) {
      printf("   Total RHS evals:  Fs = %li,  Ff = %li\n", nfsi, nffi);
    } else {
      printf("   Total RHS evals:  Fs = %li,  Ff = %li\n", nfsi, nffe);
    }
  }
  else {
    if ((solve_type==0) || (solve_type==1) || (solve_type==3) || (solve_type==5) || (solve_type==7)) {
      printf("   Total RHS evals:  Fs = %li,  Ff = %li\n", nfse, nffi);
    } else {
      printf("   Total RHS evals:  Fs = %li,  Ff = %li\n", nfse, nffe);
    }
  }

  /* Get/print slow integrator decoupled implicit solver statistics */
  if (solve_type>1) {
    retval = MRIStepGetNonlinSolvStats(arkode_mem, &nnis, &nncs);
    check_retval(&retval, "MRIStepGetNonlinSolvStats", 1);
    retval = MRIStepGetNumJacEvals(arkode_mem, &njes);
    check_retval(&retval, "MRIStepGetNumJacEvals", 1);
    printf("   Slow Newton iters = %li\n", nnis);
    printf("   Slow Newton conv fails = %li\n", nncs);
    printf("   Slow Jacobian evals = %li\n", njes);
  }

  /* Get/print fast integrator implicit solver statistics */
  if ((solve_type==0) || (solve_type==1) || (solve_type==3) || (solve_type==5) || (solve_type==7)) {
    retval = ARKStepGetNonlinSolvStats(inner_arkode_mem, &nnif, &nncf);
    check_retval(&retval, "ARKStepGetNonlinSolvStats", 1);
    retval = ARKStepGetNumJacEvals(inner_arkode_mem, &njef);
    check_retval(&retval, "ARKStepGetNumJacEvals", 1);
    printf("   Fast Newton iters = %li\n", nnif);
    printf("   Fast Newton conv fails = %li\n", nncf);
    printf("   Fast Jacobian evals = %li\n", njef);
  }

  /* Clean up and return with successful completion */
  free(udata);                               /* Free user data             */
  ARKStepFree(&inner_arkode_mem);            /* Free integrator memory     */
  MRIStepInnerStepper_Free(&inner_stepper);  /* Free inner stepper         */
  MRIStepFree(&arkode_mem);                  /* Free integrator memory     */
  ARKodeButcherTable_Free(B);                /* Free Butcher table         */
  MRIStepCoupling_Free(C);                   /* Free coupling coefficients */
  SUNMatDestroy(Af);                         /* Free fast matrix           */
  SUNLinSolFree(LSf);                        /* Free fast linear solver    */
  SUNLinSolFree(LSs);                        /* Free slow linear solver    */
  SUNMatDestroy(As);                         /* Free slow matrix           */
  N_VDestroy(y);                             /* Free vectors               */
  N_VDestroy(umask);
  N_VDestroy(vmask);
  N_VDestroy(wmask);
  SUNContext_Free(&ctx);

  return 0;
}

/*------------------------------------
 * Functions called by the integrator
 *------------------------------------*/

/* ff routine to compute the fast portion of the ODE RHS. */
static int ff(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;      /* access problem data */
  sunindextype N = udata->N;                  /* set variable shortcuts */
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;
  realtype *Ydata=NULL, *dYdata=NULL;
  realtype u, v, w;
  sunindextype i;

  Ydata = N_VGetArrayPointer(y);     /* access data arrays */
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  dYdata = N_VGetArrayPointer(ydot);
  if (check_retval((void *)dYdata, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                        /* initialize ydot to zero */

  /* iterate over domain, computing all equations */
  for (i=1; i<N-1; i++) {
    /* set shortcuts */
    u = Ydata[IDX(i,0)];
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* Fill in ODE RHS for u */
    dYdata[IDX(i,0)] = a - (w+ONE)*u + v*u*u;

    /* Fill in ODE RHS for v */
    dYdata[IDX(i,1)] = w*u - v*u*u;

    /* Fill in ODE RHS for w */
    dYdata[IDX(i,2)] = (b-w)/ep - w*u;
  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  /* Return with success */
  return 0;
}

/* fse routine to compute the slow-explicit portion of the ODE RHS function. */
static int fse(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;      /* access problem data */
  sunindextype N = udata->N;                  /* set variable shortcuts */
  realtype au = udata->au;
  realtype av = udata->av;
  realtype aw = udata->aw;
  realtype dx = udata->dx;
  realtype *Ydata=NULL, *dYdata=NULL;
  realtype auconst, avconst, awconst, ul, ur, vl, vr, wl, wr;
  sunindextype i;

  Ydata = N_VGetArrayPointer(y);     /* access data arrays */
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  dYdata = N_VGetArrayPointer(ydot);
  if (check_retval((void *)dYdata, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                        /* initialize ydot to zero */

  /* iterate over domain, computing all equations */
  auconst = -au/RCONST(2.0)/dx;
  avconst = -av/RCONST(2.0)/dx;
  awconst = -aw/RCONST(2.0)/dx;
  for (i=1; i<N-1; i++) {
    /* set shortcuts */
    ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* Fill in ODE RHS for u */
    dYdata[IDX(i,0)] = (ur - ul)*auconst;

    /* Fill in ODE RHS for v */
    dYdata[IDX(i,1)] = (vr - vl)*avconst;

    /* Fill in ODE RHS for w */
    dYdata[IDX(i,2)] = (wr - wl)*awconst;
  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  /* Return with success */
  return 0;
}

/* fsi routine to compute the slow-implicit portion of the  ODE RHS. */
static int fsi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;      /* access problem data */
  sunindextype N = udata->N;                  /* set variable shortcuts */
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  realtype dx = udata->dx;
  realtype *Ydata=NULL, *dYdata=NULL;
  realtype duconst, dvconst, dwconst, u, ul, ur, v, vl, vr, w, wl, wr;
  sunindextype i;

  Ydata = N_VGetArrayPointer(y);     /* access data arrays */
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  dYdata = N_VGetArrayPointer(ydot);
  if (check_retval((void *)dYdata, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                        /* initialize ydot to zero */

  /* iterate over domain, computing all equations */
  duconst = du/dx/dx;
  dvconst = dv/dx/dx;
  dwconst = dw/dx/dx;
  for (i=1; i<N-1; i++) {
    /* set shortcuts */
    u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* Fill in ODE RHS for u */
    dYdata[IDX(i,0)] = (ul - RCONST(2.0)*u + ur)*duconst;

    /* Fill in ODE RHS for v */
    dYdata[IDX(i,1)] = (vl - RCONST(2.0)*v + vr)*dvconst;

    /* Fill in ODE RHS for w */
    dYdata[IDX(i,2)] = (wl - RCONST(2.0)*w + wr)*dwconst;
  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  /* Return with success */
  return 0;
}

/* fs routine to compute the slow portion of the ODE RHS. */
static int fs(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = (UserData) user_data;      /* access problem data */
  sunindextype N = udata->N;                  /* set variable shortcuts */
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  realtype au = udata->au;
  realtype av = udata->av;
  realtype aw = udata->aw;
  realtype dx = udata->dx;
  realtype *Ydata=NULL, *dYdata=NULL;
  realtype duconst, dvconst, dwconst, auconst, avconst, awconst, u, ul, ur, v, vl, vr, w, wl, wr;
  sunindextype i;

  Ydata = N_VGetArrayPointer(y);     /* access data arrays */
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  dYdata = N_VGetArrayPointer(ydot);
  if (check_retval((void *)dYdata, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                        /* initialize ydot to zero */

  /* iterate over domain, computing all equations */
  duconst = du/dx/dx;
  dvconst = dv/dx/dx;
  dwconst = dw/dx/dx;
  auconst = -au/TWO/dx;
  avconst = -av/TWO/dx;
  awconst = -aw/TWO/dx;
  for (i=1; i<N-1; i++) {
    /* set shortcuts */
    u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* Fill in ODE RHS for u */
    dYdata[IDX(i,0)] = (ul - TWO*u + ur)*duconst + (ur - ul)*auconst;

    /* Fill in ODE RHS for v */
    dYdata[IDX(i,1)] = (vl - TWO*v + vr)*dvconst + (vr - vl)*avconst;

    /* Fill in ODE RHS for w */
    dYdata[IDX(i,2)] = (wl - TWO*w + wr)*dwconst + (wr - wl)*awconst;
  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  /* Return with success */
  return 0;
}

/* f routine to compute the full ODE RHS function. */
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
  realtype au = udata->au;
  realtype av = udata->av;
  realtype aw = udata->aw;
  realtype dx = udata->dx;
  realtype *Ydata=NULL, *dYdata=NULL;
  realtype duconst, dvconst, dwconst, auconst, avconst, awconst, u, ul, ur, v, vl, vr, w, wl, wr;
  sunindextype i;

  Ydata = N_VGetArrayPointer(y);     /* access data arrays */
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  dYdata = N_VGetArrayPointer(ydot);
  if (check_retval((void *)dYdata, "N_VGetArrayPointer", 0)) return 1;
  N_VConst(0.0, ydot);                        /* initialize ydot to zero */

  /* iterate over domain, computing all equations */
  duconst = du/dx/dx;
  dvconst = dv/dx/dx;
  dwconst = dw/dx/dx;
  auconst = -au/RCONST(2.0)/dx;
  avconst = -av/RCONST(2.0)/dx;
  awconst = -aw/RCONST(2.0)/dx;
  for (i=1; i<N-1; i++) {
    /* set shortcuts */
    u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* Fill in ODE RHS for u */
    dYdata[IDX(i,0)] = (ul - RCONST(2.0)*u + ur)*duconst + (ur - ul)*auconst +  a - (w+RCONST(1.0))*u + v*u*u;

    /* Fill in ODE RHS for v */
    dYdata[IDX(i,1)] = (vl - RCONST(2.0)*v + vr)*dvconst + (vr - vl)*avconst +  w*u - v*u*u;

    /* Fill in ODE RHS for w */
    dYdata[IDX(i,2)] = (wl - RCONST(2.0)*w + wr)*dwconst + (wr - wl)*awconst + (b-w)/ep - w*u;
  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  /* Return with success */
  return 0;
}

/* Placeholder function of zeroes */
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(ZERO,ydot);
  return(0);
}

/* Jf routine to compute Jacobian of the fast portion of the ODE RHS */
static int Jf(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData udata = (UserData) user_data;     /* access problem data         */
  SUNMatZero(J);                             /* Initialize Jacobian to zero */

  /* Add in the Jacobian of the reaction terms matrix */
  ReactionJac(RCONST(1.0), y, J, udata);

  /* Return with success */
  return 0;
}

/* Jsi routine to compute the Jacobian of the slow-implicit portion of the ODE RHS. */
static int Jsi(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData udata = (UserData) user_data;     /* access problem data */
  SUNMatZero(J);                             /* Initialize Jacobian to zero */

  /* Fill in the Laplace matrix */
  LaplaceMatrix(RCONST(1.0), J, udata);

  /* Return with success */
  return 0;
}

/* Js routine to compute the Jacobian of the slow portion of ODE RHS. */
static int Js(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData udata = (UserData) user_data;     /* access problem data         */
  SUNMatZero(J);                             /* Initialize Jacobian to zero */

  /* Fill in the Laplace matrix */
  LaplaceMatrix(RCONST(1.0), J, udata);

  /* Add Jacobian of the advection terms  */
  AdvectionJac(RCONST(1.0), J, udata);

  /* Return with success */
  return 0;
}

/* Jac routine to compute the Jacobian of the full ODE RHS. */
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData udata = (UserData) user_data;     /* access problem data         */
  SUNMatZero(J);                             /* Initialize Jacobian to zero */

  /* Fill in the Laplace matrix */
  LaplaceMatrix(RCONST(1.0), J, udata);

  /* Add Jacobian of the advection terms  */
  AdvectionJac(RCONST(1.0), J, udata);

  /* Add in the Jacobian of the reaction terms matrix */
  ReactionJac(RCONST(1.0), y, J, udata);

  /* Return with success */
  return 0;
}




/*-------------------------------
 * Private helper functions
 *-------------------------------*/
/* Set the initial condition */
static int SetIC(N_Vector y, void *user_data)
{
  UserData     udata = (UserData) user_data;  /* access problem data    */
  sunindextype N     = udata->N;              /* set variable shortcuts */
  realtype     a     = udata->a;
  realtype     b     = udata->b;
  realtype     dx    = udata->dx;
  realtype     pi    = udata->pi;
  realtype*    data  = NULL;
  sunindextype i;

  /* Access data array from NVector y */
  data = N_VGetArrayPointer(y);

  /* Set initial conditions into y */
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*i*dx);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*i*dx);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*i*dx);  /* w */
  }

  /* Return  with success */
  return(0);
}

/* Routine to compute the Jacobian matrix from fse(t,y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int AdvectionJac(realtype c, SUNMatrix Jac, UserData udata)
{
  /* Set shortcuts */
  sunindextype N = udata->N;
  realtype dx    = udata->dx;
  realtype au = udata->au;
  realtype av = udata->av;
  realtype aw = udata->aw;
  sunindextype i;
  realtype auconst, avconst, awconst;
  auconst = -au/TWO/dx;
  avconst = -av/TWO/dx;
  awconst = -aw/TWO/dx;

  /* iterate over intervals, filling in Jacobian of (L*y) using SM_ELEMENT_B
     macro (see sunmatrix_band.h) */
  for (i=1; i<N-1; i++) {
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i-1,0)) += -c*auconst;
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i-1,1)) += -c*avconst;
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i-1,2)) += -c*awconst;
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i+1,0)) += c*auconst;
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i+1,1)) += c*avconst;
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i+1,2)) += c*awconst;
  }

  /* Return with success */
  return 0;
}

/* Routine to compute the stiffness matrix from (L*y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int LaplaceMatrix(realtype c, SUNMatrix Jac, UserData udata)
{
  sunindextype N = udata->N;           /* set shortcuts */
  realtype dx    = udata->dx;
  sunindextype i;

  /* iterate over intervals, filling in Jacobian of (L*y) using SM_ELEMENT_B
     macro (see sunmatrix_band.h) */
  for (i=1; i<N-1; i++) {
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i-1,0)) += c*udata->du/dx/dx;
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i-1,1)) += c*udata->dv/dx/dx;
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i-1,2)) += c*udata->dw/dx/dx;
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i,0)) += -c*RCONST(2.0)*udata->du/dx/dx;
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i,1)) += -c*RCONST(2.0)*udata->dv/dx/dx;
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i,2)) += -c*RCONST(2.0)*udata->dw/dx/dx;
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i+1,0)) += c*udata->du/dx/dx;
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i+1,1)) += c*udata->dv/dx/dx;
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i+1,2)) += c*udata->dw/dx/dx;
  }

  /* Return with success */
  return 0;
}

/* Routine to compute the Jacobian matrix from R(y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int ReactionJac(realtype c, N_Vector y, SUNMatrix Jac, UserData udata)
{
  sunindextype N = udata->N;                      /* set shortcuts */
  realtype ep = udata->ep;
  sunindextype i;
  realtype u, v, w;
  realtype *Ydata = N_VGetArrayPointer(y);     /* access solution array */
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over nodes, filling in Jacobian of reaction terms */
  for (i=1; i<N-1; i++) {

    u = Ydata[IDX(i,0)];                       /* set nodal value shortcuts */
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* all vars wrt u */
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i,0)) += c*(TWO*u*v-(w+ONE));
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i,0)) += c*(w - TWO*u*v);
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i,0)) += c*(-w);

    /* all vars wrt v */
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i,1)) += c*(u*u);
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i,1)) += c*(-u*u);

    /* all vars wrt w */
    SM_ELEMENT_B(Jac,IDX(i,0),IDX(i,2)) += c*(-u);
    SM_ELEMENT_B(Jac,IDX(i,1),IDX(i,2)) += c*(u);
    SM_ELEMENT_B(Jac,IDX(i,2),IDX(i,2)) += c*(-ONE/ep - u);

  }

  /* Return with success */
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
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *errvalue;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errvalue = (int *) returnvalue;
    if (*errvalue < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errvalue);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


/*---- end of file ----*/
