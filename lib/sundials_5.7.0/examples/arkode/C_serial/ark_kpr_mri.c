/* ----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * Multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:
 *
 *    [u]' = [ G  e ] [(-1+u^2-r)/(2u)] + [      r'(t)/(2u)        ]
 *    [v]    [ e -1 ] [(-2+v^2-s)/(2v)]   [ s'(t)/(2*sqrt(2+s(t))) ]
 *         = [ fs(t,u,v) ]
 *           [ ff(t,u,v) ]
 *
 * where r(t) = 0.5*cos(t),  s(t) = cos(w*t),  0 < t < 5.
 *
 * This problem has analytical solution given by
 *    u(t) = sqrt(1+r(t)),  v(t) = sqrt(2+s(t)).
 *
 * We use the parameters:
 *   e = 0.5 (fast/slow coupling strength) [default]
 *   G = -1e2 (stiffness at slow time scale) [default]
 *   w = 100  (time-scale separation factor) [default]
 *   hs = 0.01 (slow step size) [default]
 *
 * The stiffness of the slow time scale is essentially determined
 * by G, for |G| > 50 it is 'stiff' and ideally suited to a
 * multirate method that is implicit at the slow time scale.
 *
 * We select the MRI method to use based on an additional input,
 * solve_type; with options (slow type-order/fast type-order):
 * 0. exp-3/exp-3 (standard MIS) [default]
 * 1. none/exp-3 (no slow, explicit fast)
 * 2. none/dirk-3 (no slow, dirk fast)
 * 3. exp-3/none (explicit slow, no fast)
 * 4. dirk-2/none (dirk slow, no fast) -- solve-decoupled
 * 5. exp-4/exp-4 (MRI-GARK-ERK45a / ERK-4-4)
 * 6. exp-4/exp-3 (MRI-GARK-ERK45a / ERK-3-3)
 * 7. dirk-3/exp-3 (MRI-GARK-ESDIRK34a / ERK-3-3) -- solve decoupled
 *
 * The program should be run with arguments in the following order:
 *   $ a.out solve_type h G w e
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out solve_type h G w
 *   $ a.out solve_type h G
 *   $ a.out solve_type h
 *   $ a.out solve_type
 *   $ a.out
 * are acceptable.  We require:
 *   * 0 <= solve_type <= 7
 *   * 0 < h < 1/|G|
 *   * G < 0.0
 *   * w >= 1.0
 *
 * This program solves the problem with the MRI stepper. Outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 * ----------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode_mristep.h>      /* prototypes for MRIStep fcts., consts */
#include <arkode/arkode_arkstep.h>      /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_serial.h>     /* serial N_Vector type, fcts., macros  */
#include <sunmatrix/sunmatrix_dense.h>  /* dense matrix type, fcts., macros     */
#include <sunlinsol/sunlinsol_dense.h>  /* dense linear solver                  */
#include <sundials/sundials_math.h>     /* def. math fcns, 'realtype'           */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/* User-supplied functions called by the solver */
static int fs(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int ff(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Js(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jn(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static realtype r(realtype t, void *user_data);
static realtype s(realtype t, void *user_data);
static realtype rdot(realtype t, void *user_data);
static realtype sdot(realtype t, void *user_data);
static realtype utrue(realtype t, void *user_data);
static realtype vtrue(realtype t, void *user_data);
static int Ytrue(realtype t, N_Vector y, void *user_data);
static int check_retval(void *returnvalue, const char *funcname, int opt);


/* Main Program */
int main(int argc, char *argv[])
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);     /* initial time */
  realtype Tf = RCONST(5.0);     /* final time */
  realtype dTout = RCONST(0.1);  /* time between outputs */
  sunindextype NEQ = 2;          /* number of dependent vars. */
  int Nt = (int) ceil(Tf/dTout); /* number of output times */
  int solve_type = 0;            /* problem configuration type */
  realtype hs = RCONST(0.01);    /* slow step size */
  realtype e = RCONST(0.5);      /* fast/slow coupling strength */
  realtype G = RCONST(-100.0);   /* stiffness at slow time scale */
  realtype w = RCONST(100.0);    /* time-scale separation factor */
  realtype reltol = RCONST(0.01);
  realtype abstol = 1e-11;

  /* general problem variables */
  int retval;                    /* reusable error-checking flag */
  N_Vector y = NULL;             /* empty vector for the computed solution */
  void *arkode_mem = NULL;       /* empty ARKode memory structure */
  void *inner_arkode_mem = NULL; /* empty ARKode memory structure */
  ARKodeButcherTable B = NULL;   /* fast method Butcher table */
  MRIStepCoupling C = NULL;      /* slow coupling coefficients */
  SUNMatrix Af = NULL;           /* empty matrix for fast solver */
  SUNLinearSolver LSf = NULL;    /* empty fast linear solver object */
  SUNMatrix As = NULL;           /* empty matrix for slow solver */
  SUNLinearSolver LSs = NULL;    /* empty slow linear solver object */
  booleantype implicit_slow;
  FILE *UFID;
  realtype hf, gamma, beta, t, tout, rpar[3];
  realtype uerr, verr, uerrtot, verrtot, errtot;
  int iout;
  long int nsts, nstf, nfs, nff, nnif, nncf, njef, nnis, nncs, njes, tmp;

  /*
   * Initialization
   */

  /* Retrieve the command-line options: solve_type h G w e */
  if (argc > 1)  solve_type = (sunindextype) atol(argv[1]);
  if (argc > 2)  hs = (realtype) atof(argv[2]);
  if (argc > 3)  G = (realtype) atof(argv[3]);
  if (argc > 4)  w = (realtype) atof(argv[4]);
  if (argc > 5)  e = (realtype) atof(argv[5]);

  /* Check arguments for validity */
  /*   0 <= solve_type <= 7      */
  /*   G < 0.0                   */
  /*   h > 0                     */
  /*   h < 1/|G| (explicit slow) */
  /*   w >= 1.0                  */
  if ((solve_type < 0) || (solve_type > 7)) {
    printf("ERROR: solve_type be an integer in [0,7] \n");
    return(-1);
  }
  if (G >= ZERO) {
    printf("ERROR: G must be a negative real number\n");
    return(-1);
  }
  implicit_slow = SUNFALSE;
  if ((solve_type == 4) || (solve_type == 7))
    implicit_slow = SUNTRUE;
  if (hs <= ZERO){
    printf("ERROR: hs must be in positive\n");
    return(-1);
  }
  if ((hs > ONE/SUNRabs(G)) && (!implicit_slow)) {
    printf("ERROR: hs must be in (0, 1/|G|)\n");
    return(-1);
  }
  if (w < ONE) {
    printf("ERROR: w must be >= 1.0\n");
    return(-1);
  }
  rpar[0] = G;
  rpar[1] = w;
  rpar[2] = e;
  hf = hs/w;

  /* Initial problem output (and set implicit solver tolerances as needed) */
  printf("\nMultirate nonlinear Kvaerno-Prothero-Robinson test problem:\n");
  printf("    time domain:  (%"GSYM",%"GSYM"]\n",T0,Tf);
  printf("    hs = %"GSYM"\n",hs);
  printf("    hf = %"GSYM"\n",hf);
  printf("    G = %"GSYM"\n",G);
  printf("    w = %"GSYM"\n",w);
  printf("    e = %"GSYM"\n",e);
  switch (solve_type) {
  case(0):
    printf("    solver: exp-3/exp-3 (standard MIS)\n\n");
    break;
  case(1):
    printf("    solver: none/exp-3 (no slow, explicit fast)\n\n");
    break;
  case(2):
    reltol = SUNMAX(hs*hs*hs, 1e-10);
    abstol = 1e-11;
    printf("    solver: none/dirk-3 (no slow, dirk fast)\n\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n", reltol, abstol);
    break;
  case(3):
    printf("    solver: exp-3/none (explicit slow, no fast)\n");
    break;
  case(4):
    reltol = SUNMAX(hs*hs, 1e-10);
    abstol = 1e-11;
    printf("    solver: dirk-2/none (dirk slow, no fast)\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n", reltol, abstol);
    break;
  case(5):
    printf("    solver: exp-4/exp-4 (MRI-GARK-ERK45a / ERK-4-4)\n\n");
    break;
  case(6):
    printf("    solver: exp-4/exp-3 (MRI-GARK-ERK45a / ERK-3-3)\n\n");
    break;
  case(7):
    reltol = SUNMAX(hs*hs*hs, 1e-10);
    abstol = 1e-11;
    printf("    solver: dirk-3/exp-3 (MRI-GARK-ESDIRK34a / ERK-3-3) -- solve decoupled\n");
    printf("    reltol = %.2"ESYM",  abstol = %.2"ESYM"\n", reltol, abstol);
    break;
  }

  /* Create and initialize serial vector for the solution */
  y = N_VNew_Serial(NEQ);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return 1;
  retval = Ytrue(T0, y, rpar);
  if (check_retval(&retval, "Ytrue", 1)) return 1;

  /*
   * Create the fast integrator and set options
   */

  /* Initialize the fast integrator. Specify the fast right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, and the
     initial dependent variable vector y. */
  switch (solve_type) {
  case(0):
  case(6):
  case(7):  /* erk-3-3 fast solver */
    inner_arkode_mem = ARKStepCreate(ff, NULL, T0, y);
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
    B->q=2;
    retval = ARKStepSetTables(inner_arkode_mem, 3, 2, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) return 1;
    break;
  case(1):  /* erk-3-3 fast solver (full problem) */
    inner_arkode_mem = ARKStepCreate(fn, NULL, T0, y);
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
    B->q=2;
    retval = ARKStepSetTables(inner_arkode_mem, 3, 2, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) return 1;
    break;
  case(5):  /* erk-4-4 fast solver */
    inner_arkode_mem = ARKStepCreate(ff, NULL, T0, y);
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
  case(2):  /* esdirk-3-3 fast solver (full problem) */
    inner_arkode_mem = ARKStepCreate(NULL, fn, T0, y);
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
    Af = SUNDenseMatrix(NEQ, NEQ);
    if (check_retval((void *)Af, "SUNDenseMatrix", 0)) return 1;
    LSf = SUNLinSol_Dense(y, Af);
    if (check_retval((void *)LSf, "SUNLinSol_Dense", 0)) return 1;
    retval = ARKStepSetLinearSolver(inner_arkode_mem, LSf, Af);
    if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return 1;
    retval = ARKStepSetJacFn(inner_arkode_mem, Jn);
    if (check_retval(&retval, "ARKStepSetJacFn", 1)) return 1;
    retval = ARKStepSStolerances(inner_arkode_mem, reltol, abstol);
    if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
    break;
  case(3):  /* no fast dynamics ('evolve' explicitly w/ erk-3-3) */
  case(4):
    inner_arkode_mem = ARKStepCreate(f0, NULL, T0, y);
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
    B->q=2;
    retval = ARKStepSetTables(inner_arkode_mem, 3, 2, NULL, B);
    if (check_retval(&retval, "ARKStepSetTables", 1)) return 1;
  }

  /* Set the user data pointer */
  retval = ARKStepSetUserData(inner_arkode_mem, (void *) rpar);
  if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;

  /* Set the fast step size */
  retval = ARKStepSetFixedStep(inner_arkode_mem, hf);
  if (check_retval(&retval, "ARKStepSetFixedStep", 1)) return 1;

  /*
   * Create the slow integrator and set options
   */

  /* Initialize the slow integrator. Specify the slow right-hand side
     function in y'=fs(t,y)+ff(t,y), the inital time T0, the
     initial dependent variable vector y, and the fast integrator. */
  switch (solve_type) {
  case(0):  /* KW3 slow solver */
    arkode_mem = MRIStepCreate(fs, T0, y, MRISTEP_ARKSTEP, inner_arkode_mem);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;
    retval = MRIStepSetTableNum(arkode_mem, KNOTH_WOLKE_3_3);
    if (check_retval(&retval, "MRIStepSetTableNum", 1)) return 1;
    break;
  case(3):  /* KW3 slow solver (full problem) */
    arkode_mem = MRIStepCreate(fn, T0, y, MRISTEP_ARKSTEP, inner_arkode_mem);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;
    retval = MRIStepSetTableNum(arkode_mem, KNOTH_WOLKE_3_3);
    if (check_retval(&retval, "MRIStepSetTableNum", 1)) return 1;
    break;
  case(5):  /* MRI-GARK-ERK45a slow solver */
  case(6):
    arkode_mem = MRIStepCreate(fs, T0, y, MRISTEP_ARKSTEP, inner_arkode_mem);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;
    C = MRIStepCoupling_LoadTable(MRI_GARK_ERK45a);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) return 1;
    retval = MRIStepSetCoupling(arkode_mem, C);
    if (check_retval(&retval, "MRIStepSetCoupling", 1)) return 1;
    break;
  case(1):
  case(2):  /* no slow dynamics (use ERK-2-2) */
    arkode_mem = MRIStepCreate(f0, T0, y, MRISTEP_ARKSTEP, inner_arkode_mem);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;
    B = ARKodeButcherTable_Alloc(2, SUNFALSE);
    if (check_retval((void *)B, "ARKodeButcherTable_Alloc", 0)) return 1;
    B->A[1][0] = TWO/RCONST(3.0);
    B->b[0] = RCONST(0.25);
    B->b[1] = RCONST(0.75);
    B->c[1] = TWO/RCONST(3.0);
    B->q=2;
    retval = MRIStepSetTable(arkode_mem, 2, B);
    if (check_retval(&retval, "MRIStepSetTable", 1)) return 1;
    break;
  case(4):  /* dirk-2 (trapezoidal), solve-decoupled slow solver */
    arkode_mem = MRIStepCreate(fn, T0, y, MRISTEP_ARKSTEP, inner_arkode_mem);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;
    C = MRIStepCoupling_LoadTable(MRI_GARK_IRK21a);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) return 1;
    retval = MRIStepSetCoupling(arkode_mem, C);
    if (check_retval(&retval, "MRIStepSetCoupling", 1)) return 1;
    As = SUNDenseMatrix(NEQ, NEQ);
    if (check_retval((void *)As, "SUNDenseMatrix", 0)) return 1;
    LSs = SUNLinSol_Dense(y, As);
    if (check_retval((void *)LSs, "SUNLinSol_Dense", 0)) return 1;
    retval = MRIStepSetLinearSolver(arkode_mem, LSs, As);
    if (check_retval(&retval, "MRIStepSetLinearSolver", 1)) return 1;
    retval = MRIStepSetJacFn(arkode_mem, Jn);
    if (check_retval(&retval, "MRIStepSetJacFn", 1)) return 1;
    retval = MRIStepSStolerances(arkode_mem, reltol, abstol);
    if (check_retval(&retval, "MRIStepSStolerances", 1)) return 1;
    break;
  case(7):  /* MRI-GARK-ESDIRK34a, solve-decoupled slow solver */
    arkode_mem = MRIStepCreate(fs, T0, y, MRISTEP_ARKSTEP, inner_arkode_mem);
    if (check_retval((void *)arkode_mem, "MRIStepCreate", 0)) return 1;
    C = MRIStepCoupling_LoadTable(MRI_GARK_ESDIRK34a);
    if (check_retval((void*)C, "MRIStepCoupling_LoadTable", 1)) return 1;
    retval = MRIStepSetCoupling(arkode_mem, C);
    if (check_retval(&retval, "MRIStepSetCoupling", 1)) return 1;
    As = SUNDenseMatrix(NEQ, NEQ);
    if (check_retval((void *)As, "SUNDenseMatrix", 0)) return 1;
    LSs = SUNLinSol_Dense(y, As);
    if (check_retval((void *)LSs, "SUNLinSol_Dense", 0)) return 1;
    retval = MRIStepSetLinearSolver(arkode_mem, LSs, As);
    if (check_retval(&retval, "MRIStepSetLinearSolver", 1)) return 1;
    retval = MRIStepSetJacFn(arkode_mem, Js);
    if (check_retval(&retval, "MRIStepSetJacFn", 1)) return 1;
    retval = MRIStepSStolerances(arkode_mem, reltol, abstol);
    if (check_retval(&retval, "MRIStepSStolerances", 1)) return 1;
    break;
  }

  /* Set the user data pointer */
  retval = MRIStepSetUserData(arkode_mem, (void *) rpar);
  if (check_retval(&retval, "MRIStepSetUserData", 1)) return 1;

  /* Set the slow step size */
  retval = MRIStepSetFixedStep(arkode_mem, hs);
  if (check_retval(&retval, "MRIStepSetFixedStep", 1)) return 1;

  /*
   * Integrate ODE
   */

  /* Open output stream for results, output comment line */
  UFID = fopen("ark_kpr_mri_solution.txt","w");
  fprintf(UFID,"# t u v uerr verr\n");

  /* output initial condition to disk */
  fprintf(UFID," %.16"ESYM" %.16"ESYM" %.16"ESYM" %.16"ESYM" %.16"ESYM"\n",
          T0, NV_Ith_S(y,0), NV_Ith_S(y,1),
          SUNRabs(NV_Ith_S(y,0)-utrue(T0,rpar)),
          SUNRabs(NV_Ith_S(y,1)-vtrue(T0,rpar)));

  /* Main time-stepping loop: calls MRIStepEvolve to perform the
     integration, then prints results. Stops when the final time
     has been reached */
  t = T0;
  tout = T0+dTout;
  uerr = ZERO;
  verr = ZERO;
  uerrtot = ZERO;
  verrtot = ZERO;
  errtot  = ZERO;
  printf("        t           u           v       uerr      verr\n");
  printf("   ------------------------------------------------------\n");
  printf("  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"  %.2"ESYM"  %.2"ESYM"\n",
         t, NV_Ith_S(y,0), NV_Ith_S(y,1), uerr, verr);

  for (iout=0; iout<Nt; iout++) {

    /* call integrator */
    retval = MRIStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "MRIStepEvolve", 1)) break;

    /* access/print solution and error */
    uerr = SUNRabs(NV_Ith_S(y,0)-utrue(t,rpar));
    verr = SUNRabs(NV_Ith_S(y,1)-vtrue(t,rpar));
    printf("  %10.6"FSYM"  %10.6"FSYM"  %10.6"FSYM"  %.2"ESYM"  %.2"ESYM"\n",
           t, NV_Ith_S(y,0), NV_Ith_S(y,1), uerr, verr);
    fprintf(UFID," %.16"ESYM" %.16"ESYM" %.16"ESYM" %.16"ESYM" %.16"ESYM"\n",
            t, NV_Ith_S(y,0), NV_Ith_S(y,1), uerr, verr);
    uerrtot += uerr*uerr;
    verrtot += verr*verr;
    errtot += uerr*uerr + verr*verr;

    /* successful solve: update time */
    tout += dTout;
    tout = (tout > Tf) ? Tf : tout;
  }
  uerrtot = SUNRsqrt(uerrtot / Nt);
  verrtot = SUNRsqrt(verrtot / Nt);
  errtot = SUNRsqrt(errtot / Nt / 2);
  printf("   ------------------------------------------------------\n");
  fclose(UFID);

  /*
   * Finalize
   */

  /* Get some slow integrator statistics */
  retval = MRIStepGetNumSteps(arkode_mem, &nsts);
  check_retval(&retval, "MRIStepGetNumSteps", 1);
  retval = MRIStepGetNumRhsEvals(arkode_mem, &nfs);
  check_retval(&retval, "MRIStepGetNumRhsEvals", 1);

  /* Get some fast integrator statistics */
  retval = ARKStepGetNumSteps(inner_arkode_mem, &nstf);
  check_retval(&retval, "ARKStepGetNumSteps", 1);
  retval = ARKStepGetNumRhsEvals(inner_arkode_mem, &nff, &tmp);
  check_retval(&retval, "ARKStepGetNumRhsEvals", 1);

  /* Print some final statistics */
  printf("\nFinal Solver Statistics:\n");
  printf("   Steps: nsts = %li, nstf = %li\n", nsts, nstf);
  printf("   u error = %.3"ESYM", v error = %.3"ESYM", total error = %.3"ESYM"\n",
         uerrtot, verrtot, errtot);
  printf("   Total RHS evals:  Fs = %li,  Ff = %li\n", nfs, nff);

  /* Get/print slow integrator decoupled implicit solver statistics */
  if ((solve_type==4) || (solve_type==7)) {
    retval = MRIStepGetNonlinSolvStats(arkode_mem, &nnis, &nncs);
    check_retval(&retval, "MRIStepGetNonlinSolvStats", 1);
    retval = MRIStepGetNumJacEvals(arkode_mem, &njes);
    check_retval(&retval, "MRIStepGetNumJacEvals", 1);
    printf("   Slow Newton iters = %li\n", nnis);
    printf("   Slow Newton conv fails = %li\n", nncs);
    printf("   Slow Jacobian evals = %li\n", njes);
  }

  /* Get/print fast integrator implicit solver statistics */
  if (solve_type == 2) {
    retval = ARKStepGetNonlinSolvStats(inner_arkode_mem, &nnif, &nncf);
    check_retval(&retval, "ARKStepGetNonlinSolvStats", 1);
    retval = ARKStepGetNumJacEvals(inner_arkode_mem, &njef);
    check_retval(&retval, "ARKStepGetNumJacEvals", 1);
    printf("   Fast Newton iters = %li\n", nnif);
    printf("   Fast Newton conv fails = %li\n", nncf);
    printf("   Fast Jacobian evals = %li\n", njef);
  }

  /* Clean up and return */
  N_VDestroy(y);                  /* Free y vector */
  ARKodeButcherTable_Free(B);     /* Butcher table */
  MRIStepCoupling_Free(C);        /* free coupling coefficients */
  SUNMatDestroy(Af);              /* free fast matrix */
  SUNLinSolFree(LSf);             /* free fast linear solver */
  ARKStepFree(&inner_arkode_mem); /* Free fast integrator memory */
  MRIStepFree(&arkode_mem);       /* Free slow integrator memory */

  return 0;
}

/* ------------------------------
 * Functions called by the solver
 * ------------------------------*/

/* ff routine to compute the fast portion of the ODE RHS. */
static int ff(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rpar = (realtype *) user_data;
  const realtype e = rpar[2];
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype tmp1, tmp2;

  /* fill in the RHS function:
     [0  0]*[(-1+u^2-r(t))/(2*u)] + [         0          ]
     [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))] */
  tmp1 = (-ONE+u*u-r(t,rpar))/(TWO*u);
  tmp2 = (-TWO+v*v-s(t,rpar))/(TWO*v);
  NV_Ith_S(ydot,0) = ZERO;
  NV_Ith_S(ydot,1) = e*tmp1 - tmp2 + sdot(t,rpar)/(TWO*vtrue(t,rpar));

  /* Return with success */
  return 0;
}

/* fs routine to compute the slow portion of the ODE RHS. */
static int fs(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rpar = (realtype *) user_data;
  const realtype G = rpar[0];
  const realtype e = rpar[2];
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype tmp1, tmp2;

  /* fill in the RHS function:
     [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
     [0 0] [(-2+v^2-s(t))/(2*v)]    [      0      ] */
  tmp1 = (-ONE+u*u-r(t,rpar))/(TWO*u);
  tmp2 = (-TWO+v*v-s(t,rpar))/(TWO*v);
  NV_Ith_S(ydot,0) = G*tmp1 + e*tmp2 + rdot(t,rpar)/(TWO*u);
  NV_Ith_S(ydot,1) = ZERO;

  /* Return with success */
  return 0;
}

static int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rpar = (realtype *) user_data;
  const realtype G = rpar[0];
  const realtype e = rpar[2];
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype tmp1, tmp2;

  /* fill in the RHS function:
     [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
     [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))] */
  tmp1 = (-ONE+u*u-r(t,rpar))/(TWO*u);
  tmp2 = (-TWO+v*v-s(t,rpar))/(TWO*v);
  NV_Ith_S(ydot,0) = G*tmp1 + e*tmp2 + rdot(t,rpar)/(TWO*u);
  NV_Ith_S(ydot,1) = e*tmp1 - tmp2 + sdot(t,rpar)/(TWO*vtrue(t,rpar));

  /* Return with success */
  return 0;
}

static int f0(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(ZERO,ydot);
  return(0);
}

static int Js(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rpar = (realtype *) user_data;
  const realtype G = rpar[0];
  const realtype e = rpar[2];
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);

  /* fill in the Jacobian:
     [G/2 + (G*(1+r(t))+rdot(t))/(2*u^2)   e/2+e*(2+s(t))/(2*v^2)]
     [                 0                             0           ] */
  SM_ELEMENT_D(J,0,0) = G/TWO + (G*(ONE+r(t,rpar))+rdot(t,rpar))/(2*u*u);
  SM_ELEMENT_D(J,0,1) = e/TWO+e*(TWO+s(t,rpar))/(TWO*v*v);
  SM_ELEMENT_D(J,1,0) = ZERO;
  SM_ELEMENT_D(J,1,1) = ZERO;

  /* Return with success */
  return 0;
}

static int Jn(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rpar = (realtype *) user_data;
  const realtype G = rpar[0];
  const realtype e = rpar[2];
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);

  /* fill in the Jacobian:
     [G/2 + (G*(1+r(t))-rdot(t))/(2*u^2)     e/2+e*(2+s(t))/(2*v^2)]
     [e/2+e*(1+r(t))/(2*u^2)                -1/2 - (2+s(t))/(2*v^2)] */
  SM_ELEMENT_D(J,0,0) = G/TWO + (G*(ONE+r(t,rpar))-rdot(t,rpar))/(2*u*u);
  SM_ELEMENT_D(J,0,1) = e/TWO+e*(TWO+s(t,rpar))/(TWO*v*v);
  SM_ELEMENT_D(J,1,0) = e/TWO+e*(ONE+r(t,rpar))/(TWO*u*u);
  SM_ELEMENT_D(J,1,1) = -ONE/TWO - (TWO+s(t,rpar))/(TWO*v*v);

  /* Return with success */
  return 0;
}


/* ------------------------------
 * Private helper functions
 * ------------------------------*/

static realtype r(realtype t, void *user_data)
{
  return( RCONST(0.5)*cos(t) );
}
static realtype s(realtype t, void *user_data)
{
  realtype *rpar = (realtype *) user_data;
  return( cos(rpar[1]*t) );
}
static realtype rdot(realtype t, void *user_data)
{
  return( -RCONST(0.5)*sin(t) );
}
static realtype sdot(realtype t, void *user_data)
{
  realtype *rpar = (realtype *) user_data;
  return( -rpar[1]*sin(rpar[1]*t) );
}
static realtype utrue(realtype t, void *user_data)
{
  return( SUNRsqrt(ONE+r(t,user_data)) );
}
static realtype vtrue(realtype t, void *user_data)
{
  return( SUNRsqrt(TWO+s(t,user_data)) );
}
static int Ytrue(realtype t, N_Vector y, void *user_data)
{
  NV_Ith_S(y,0) = utrue(t,user_data);
  NV_Ith_S(y,1) = vtrue(t,user_data);
  return(0);
}


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
