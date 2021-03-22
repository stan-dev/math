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
 * Nonlinear Kvaerno-Prothero-Robinson ODE test problem with
 * possibly-time-dependent mass "matrix" and analytical solution,
 *
 *   M(t) * [u]' = M(t) * [ G  e ] [(u^2-r-1)/(2u)] + M(t) * [ r'(t)/(2u) ]
 *          [v]           [ e -1 ] [(v^2-s-2)/(2v)]          [ s'(t)/(2v) ]
 *
 * where r(t) = 0.5*cos(t),  s(t) = sin(t),  -3 < t < 7. This problem has
 * analytical solution given by
 *    u(t) = sqrt(1+r(t)),  v(t) = sqrt(2+s(t)).
 *
 * The mass matrix is either time-dependent,
 *
 *    M(t) = g*[cos(t) sin(t); -sin(t) cos(t)],
 *
 * or fixed throughout the simulation at M(pi/4), dependong on the
 * command-line argument TD (0=fixed, 1=time-dependent), with 1 the
 * default.
 *
 * We use the parameters: e = 0.5, G = -100 [default], g = 10 [default]
 *
 * The stiffness of the problem is essentially determined
 * by G, for |G| > 50 it is 'stiff' and ideally suited to an
 * implicit or IMEX method.
 *
 * The 'strength' of the mass matrix is proportional to the value of g.
 *
 * We use the ARKStep solver, with either an ARK (0), DIRK (1) or ERK (2)
 * method to solve this problem.  This defaults to using a 4th-order
 * ARK method.  For ARK and DIRK methods, we use either the Newton
 * nonlinear solver [default], or the accelerated fixed-point solver,
 * based on the N command-line option (0=Newton, 1=FP).
 *
 * By default, all runs use temporal adaptivity; however, if the 'ord'
 * command-line input is negative, we run with order |ord|, using a
 * variety of fixed step sizes, to estimate the attained order of
 * accuracy for the method.
 *
 * The program should be run with arguments in the following order:
 *   $ a.out RK ord N G TD g
 * Not all arguments are required, but these must be omitted from
 * end-to-beginning, i.e. any one of
 *   $ a.out RK ord N G TD
 *   $ a.out RK ord N G
 *   $ a.out RK ord N
 *   $ a.out RK ord
 *   $ a.out RK
 *   $ a.out
 * are acceptable.  We require:
 *   * 0 <= RK <= 2
 *   * G < 0.0
 *
 * When run with ord>=0 (temporal adaptivity enabled), outputs are
 * printed at equal intervals of 0.1 and run statistics are printed
 * at the end.
 *
 * When run with ord<0 (order of convergence tests), errors are
 * output for various fixed step sizes, and the overall estimated
 * convergence rate is printed at the end.
 * ----------------------------------------------------------------*/

// Header files
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <vector>
#include <arkode/arkode_arkstep.h>      // prototypes for ARKStep fcts., consts
#include <nvector/nvector_serial.h>     // serial N_Vector type, fcts., macros
#include <sunmatrix/sunmatrix_dense.h>  // dense matrix type, fcts., macros
#include <sunlinsol/sunlinsol_dense.h>  // dense linear solver
#include <sundials/sundials_math.h>     // def. math fcns, 'realtype'
#include <sunnonlinsol/sunnonlinsol_newton.h>      // Newton nonlinear solver
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>  // fixed-point nonlinear solver

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
#define PI4  RCONST(0.78539816339744830961566084581987572)

using namespace std;

// User data structure
struct UserData
{
  realtype G;
  realtype g;
  realtype e;
  booleantype M_timedep;
};

// User-supplied functions called by the solver
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Ji(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int Jn(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int MassMatrix(realtype t, SUNMatrix M, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private utility functions
static int adaptive_run(void *arkode_mem, N_Vector y, realtype T0, realtype Tf,
                        realtype dTout, realtype reltol, realtype abstol,
                        int rk_type, int nls_type, UserData &udata);
static int check_order(void *arkode_mem, N_Vector y, realtype T0,
                       realtype Tf, int order, int rk_type, UserData &udata);
static realtype r(realtype t);
static realtype s(realtype t);
static realtype rdot(realtype t);
static realtype sdot(realtype t);
static realtype utrue(realtype t);
static realtype vtrue(realtype t);
static int Ytrue(realtype t, N_Vector y);
static int check_retval(void *returnvalue, const char *funcname, int opt);

// Main Program
int main(int argc, char *argv[])
{
  // general problem parameters
  realtype T0 = RCONST(-3.0);       // initial time
  realtype Tf = RCONST(7.0);        // final time
  realtype dTout = RCONST(0.1);     // time between outputs
  sunindextype NEQ = 2;             // number of dependent vars.
  int rk_type = 0;                  // type of RK method [ARK=0, DIRK=1, ERK=2]
  int nls_type = 0;                 // type of nonlinear solver [Newton=0, FP=1]
  int order = 4;                    // order of accuracy for RK method
  booleantype adaptive = SUNTRUE;   // adaptive run vs convergence order
  realtype reltol = RCONST(1e-5);   // relative tolerance
  realtype abstol = RCONST(1e-11);  // absolute tolerance

  // general problem variables
  int retval;                    // reusable error-checking flag
  N_Vector y = NULL;             // empty vector for the computed solution
  void *arkode_mem = NULL;       // empty ARKStep memory structure
  SUNMatrix A = NULL;            // empty system matrix
  SUNMatrix M = NULL;            // empty mass matrix
  SUNLinearSolver LS = NULL;     // empty system linear solver object
  SUNLinearSolver MLS = NULL;    // empty mass linear solver object
  SUNNonlinearSolver NLS = NULL; // empty nonlinear solver object
  UserData udata;                // user-data structure
  udata.G = RCONST(-100.0);      // stiffness parameter
  udata.g = RCONST(10.0);        // mass matrix parameter
  udata.e = RCONST(0.5);         // coupling strength
  udata.M_timedep = SUNTRUE;     // time-dependence of mass matrix

  //
  // Initialization
  //

  // Retrieve the command-line options: RK ord N G TD g
  if (argc > 1)  rk_type = (int) atoi(argv[1]);
  if (argc > 2)  order = (int) atoi(argv[2]);
  if (argc > 3)  nls_type = (int) atoi(argv[3]);
  if (argc > 4)  udata.G = (realtype) atof(argv[4]);
  if (argc > 5)  udata.M_timedep = (int) atoi(argv[5]);
  if (argc > 6)  udata.g = (realtype) atof(argv[6]);

  // Check arguments for validity
  //   0 <= rk_type <= 2
  //   G < 0.0
  if ((rk_type < 0) || (rk_type > 2)) {
    cerr << "ERROR: RK type be an integer in [0,2] \n";
    return(-1);
  }
  if (udata.G >= ZERO) {
    cerr << "ERROR: G must be a negative real number\n";
    return(-1);
  }

  // Handle adaptive run vs order-of-convergence run
  if (order < 0) {
    adaptive = SUNFALSE;
    order = abs(order);
  }
  if (order == 0)  order = 4;

  // Initial problem output (and set implicit solver tolerances as needed)
  cout << "\nNonlinear Kvaerno-Prothero-Robinson test problem with mass matrix:\n";
  cout << "    time domain:  (" << T0 << "," << Tf << "]\n";
  if (udata.M_timedep) {
    cout << "    Time-dependent mass matrix\n";
  } else {
    cout << "    Fixed mass matrix\n";
  }
  if (rk_type == 0) {
    cout << "    ARK solver, order = " << order << endl;
  } else if (rk_type == 1) {
    cout << "    DIRK solver, order = " << order << endl;
  } else {
    cout << "    ERK solver, order = " << order << endl;
  }
  if (rk_type < 2) {
    if (nls_type == 0) {
      cout << "    Newton nonlinear solver\n";
    } else {
      cout << "    Fixed-point nonlinear solver\n";
    }
  }
  if (adaptive) {
    cout << "    Adaptive run with reltol = " << reltol << ", abstol = " << abstol << endl;
  } else {
    cout << "    Order-of-convergence run\n";
  }
  cout << "    G = " << udata.G << endl;
  cout << "    g = " << udata.g << endl;
  cout << "    e = " << udata.e << endl;

  //
  // Problem Setup
  //


  // Create and initialize serial vector for the solution
  y = N_VNew_Serial(NEQ);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return 1;
  retval = Ytrue(T0, y);
  if (check_retval(&retval, "Ytrue", 1)) return 1;

  // Initialize ARKStep. Specify the right-hand side function(s) for
  // M(t) * y' = fe(t,y) + fi(t,y), the inital time T0, and the
  // initial dependent variable vector y.
  if (rk_type == 0) {         // ARK method
    arkode_mem = ARKStepCreate(fe, fi, T0, y);
  } else if (rk_type == 1) {  // DIRK method
    arkode_mem = ARKStepCreate(NULL, fn, T0, y);
  } else {                    // ERK method
    arkode_mem = ARKStepCreate(fn, NULL, T0, y);
  }
  if (check_retval((void *) arkode_mem, "ARKStepCreate", 0)) return 1;

  // Initialize/attach nonlinear and linear solvers (if required)
  if (rk_type < 2) {

    if (nls_type == 0) {   // Newton

      NLS = SUNNonlinSol_Newton(y);
      if (check_retval((void *)NLS, "SUNNonlinSol_Newton", 0)) return 1;
      retval = ARKStepSetNonlinearSolver(arkode_mem, NLS);
      if (check_retval(&retval, "ARKStepSetNonlinearSolver", 1)) return(1);

      A = SUNDenseMatrix(NEQ, NEQ);
      if (check_retval((void *)A, "SUNDenseMatrix", 0)) return 1;
      LS = SUNLinSol_Dense(y, A);
      if (check_retval((void *)LS, "SUNLinSol_Dense", 0)) return 1;
      retval = ARKStepSetLinearSolver(arkode_mem, LS, A);
      if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return(1);
      if (rk_type == 0) {
        retval = ARKStepSetJacFn(arkode_mem, Ji);
      } else {
        retval = ARKStepSetJacFn(arkode_mem, Jn);
      }
      if (check_retval(&retval, "ARKStepSetJacFn", 1)) return 1;

    } else {               // Fixed-point

      NLS = SUNNonlinSol_FixedPoint(y,4);
      if (check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0)) return 1;
      retval = ARKStepSetNonlinearSolver(arkode_mem, NLS);
      if (check_retval(&retval, "ARKStepSetNonlinearSolver", 1)) return(1);

    }
  }

  // Set maximum stepsize for ERK run
  if (rk_type == 2) {
    retval = ARKStepSetMaxStep(arkode_mem, ONE/SUNRabs(udata.G));
    if (check_retval(&retval, "ARKStepSetMaxStep", 1)) return(1);
  }

  // Initialize/attach mass matrix solver
  M = SUNDenseMatrix(NEQ, NEQ);
  if (check_retval((void *)M, "SUNDenseMatrix", 0)) return 1;
  MLS = SUNLinSol_Dense(y, M);
  if (check_retval((void *)MLS, "SUNLinSol_Dense", 0)) return 1;
  retval = ARKStepSetMassLinearSolver(arkode_mem, MLS, M, SUNTRUE);
  if (check_retval(&retval, "ARKStepSetMassLinearSolver", 1)) return(1);
  retval = ARKStepSetMassFn(arkode_mem, MassMatrix);
  if (check_retval(&retval, "ARKStepSetMassFn", 1)) return(1);

  // Set desired solver order
  retval = ARKStepSetOrder(arkode_mem, order);
  if (check_retval(&retval, "ARKStepSetOrder", 1)) return 1;

  // Set the user data pointer
  retval = ARKStepSetUserData(arkode_mem, (void *) &udata);
  if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;

  // Integrate ODE, based on run type
  if (adaptive) {
    retval = adaptive_run(arkode_mem, y, T0, Tf, dTout,
                          reltol, abstol, rk_type, nls_type, udata);
    if (check_retval(&retval, "adaptive_run", 1)) return 1;
  } else {
    retval = check_order(arkode_mem, y, T0, Tf, order, rk_type, udata);
    if (check_retval(&retval, "check_order", 1)) return 1;
  }


  // Clean up and return
  ARKStepFree(&arkode_mem);       // Free integrator memory
  SUNLinSolFree(LS);              // free system linear solver
  SUNLinSolFree(MLS);             // free mass linear solver
  SUNNonlinSolFree(NLS);          // free nonlinear solver
  SUNMatDestroy(A);               // free system matrix
  SUNMatDestroy(M);               // free mass matrix
  N_VDestroy(y);                  // Free y vector
  return 0;
}

//------------------------------
// Functions called by the solver
//------------------------------

// fe routine to compute the explicit portion of the ODE RHS.
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData *udata = (UserData *) user_data;
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype gcos, gsin, tmp1, tmp2;

  // fill in the RHS function:
  //   g*[ cos(t) sin(t)] * [rdot(t)/(2u)]
  //     [-sin(t) cos(t)]   [sdot(t)/(2v)]
  gcos = (udata->M_timedep) ? udata->g*cos(t) : udata->g*cos(PI4);
  gsin = (udata->M_timedep) ? udata->g*sin(t) : udata->g*sin(PI4);
  tmp1 = rdot(t)/(TWO*u);
  tmp2 = sdot(t)/(TWO*v);
  NV_Ith_S(ydot,0) = gcos*tmp1 + gsin*tmp2;
  NV_Ith_S(ydot,1) = -gsin*tmp1 + gcos*tmp2;

  // Return with success
  return 0;
}

// fi routine to compute the implicit portion of the ODE RHS.
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData *udata = (UserData *) user_data;
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype gcos, gsin, tmp1, tmp2, tmp3, tmp4;

  // fill in the RHS function:
  //   g*[ cos(t) sin(t)]*[G  e]*[(-1+u^2-r(t))/(2*u)]
  //     [-sin(t) cos(t)] [e -1] [(-2+v^2-s(t))/(2*v)]
  gcos = (udata->M_timedep) ? udata->g*cos(t) : udata->g*cos(PI4);
  gsin = (udata->M_timedep) ? udata->g*sin(t) : udata->g*sin(PI4);
  tmp1 = (-ONE+u*u-r(t))/(TWO*u);
  tmp2 = (-TWO+v*v-s(t))/(TWO*v);
  tmp3 = udata->G*tmp1 + udata->e*tmp2;
  tmp4 = udata->e*tmp1 - tmp2;
  NV_Ith_S(ydot,0) = gcos*tmp3 + gsin*tmp4;
  NV_Ith_S(ydot,1) = -gsin*tmp3 + gcos*tmp4;

  // Return with success
  return 0;
}

static int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData *udata = (UserData *) user_data;
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype gcos, gsin, tmp1, tmp2, tmp3, tmp4;

  // fill in the RHS function:
  //   g*[ cos(t) sin(t)]*( [G  e]*[(-1+u^2-r(t))/(2*u)] + [rdot(t)/(2u)]
  //     [-sin(t) cos(t)] ( [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2v)]
  gcos = (udata->M_timedep) ? udata->g*cos(t) : udata->g*cos(PI4);
  gsin = (udata->M_timedep) ? udata->g*sin(t) : udata->g*sin(PI4);
  tmp1 = (-ONE+u*u-r(t))/(TWO*u);
  tmp2 = (-TWO+v*v-s(t))/(TWO*v);
  tmp3 = udata->G*tmp1 + udata->e*tmp2 + rdot(t)/(TWO*u);
  tmp4 = udata->e*tmp1 - tmp2 + sdot(t)/(TWO*v);
  NV_Ith_S(ydot,0) = gcos*tmp3 + gsin*tmp4;
  NV_Ith_S(ydot,1) = -gsin*tmp3 + gcos*tmp4;

  // Return with success
  return 0;
}

static int Ji(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData *udata = (UserData *) user_data;
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype gcos, gsin, t11, t12, t21, t22, m11, m12, m21, m22;

  // fill in the Jacobian:
  //   g*[ cos(t) sin(t)]*[G  e]*[1-(u^2-r(t)-1)/(2*u^2),  0]
  //     [-sin(t) cos(t)] [e -1] [0,  1-(v^2-s(t)-2)/(2*v^2)]
  gcos = (udata->M_timedep) ? udata->g*cos(t) : udata->g*cos(PI4);
  gsin = (udata->M_timedep) ? udata->g*sin(t) : udata->g*sin(PI4);
  t11 = ONE-(u*u-r(t)-ONE)/(TWO*u*u);
  t12 = ZERO;
  t22 = ONE-(v*v-s(t)-TWO)/(TWO*v*v);
  t21 = ZERO;
  m11 = udata->G*t11 + udata->e*t21;
  m12 = udata->G*t12 + udata->e*t22;
  m21 = udata->e*t11 - t21;
  m22 = udata->e*t12 - t22;
  SM_ELEMENT_D(J,0,0) = gcos*m11 + gsin*m21;
  SM_ELEMENT_D(J,0,1) = gcos*m12 + gsin*m22;
  SM_ELEMENT_D(J,1,0) = -gsin*m11 + gcos*m21;
  SM_ELEMENT_D(J,1,1) = -gsin*m12 + gcos*m22;

  // Return with success
  return 0;
}

static int Jn(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
              void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData *udata = (UserData *) user_data;
  const realtype u = NV_Ith_S(y,0);
  const realtype v = NV_Ith_S(y,1);
  realtype gcos, gsin, t11, t12, t21, t22, m11, m12, m21, m22;

  // fill in the Jacobian:
  //   g*[ cos(t) sin(t)]*( [G  e]*[1-(u^2-r(t)-1)/(2*u^2),  0] + [-r'(t)/(2*u^2),  0])
  //     [-sin(t) cos(t)] ( [e -1] [0,  1-(v^2-s(t)-2)/(2*v^2)]   [0,  -s'(t)/(2*v^2)])
  gcos = (udata->M_timedep) ? udata->g*cos(t) : udata->g*cos(PI4);
  gsin = (udata->M_timedep) ? udata->g*sin(t) : udata->g*sin(PI4);
  t11 = ONE-(u*u-r(t)-ONE)/(TWO*u*u);
  t12 = ZERO;
  t21 = ZERO;
  t22 = ONE-(v*v-s(t)-TWO)/(TWO*v*v);
  m11 = udata->G*t11 + udata->e*t21 - rdot(t)/(TWO*u*u);
  m12 = udata->G*t12 + udata->e*t22;
  m21 = udata->e*t11 - t21;
  m22 = udata->e*t12 - t22 - sdot(t)/(TWO*v*v);
  SM_ELEMENT_D(J,0,0) = gcos*m11 + gsin*m21;
  SM_ELEMENT_D(J,0,1) = gcos*m12 + gsin*m22;
  SM_ELEMENT_D(J,1,0) = -gsin*m11 + gcos*m21;
  SM_ELEMENT_D(J,1,1) = -gsin*m12 + gcos*m22;

  // Return with success
  return 0;
}

// Routine to compute the mass matrix multiplying y_t.
static int MassMatrix(realtype t, SUNMatrix M, void *user_data,
                      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData *udata = (UserData *) user_data;

  // fill in the mass matrix: g*[ cos(t) sin(t); -sin(t) cos(t)]
  SM_ELEMENT_D(M,0,0) = (udata->M_timedep) ? udata->g*cos(t) : udata->g*cos(PI4);
  SM_ELEMENT_D(M,0,1) = (udata->M_timedep) ? udata->g*sin(t) : udata->g*sin(PI4);
  SM_ELEMENT_D(M,1,0) = (udata->M_timedep) ? -udata->g*sin(t) : -udata->g*sin(PI4);
  SM_ELEMENT_D(M,1,1) = (udata->M_timedep) ? udata->g*cos(t) : udata->g*cos(PI4);

  return 0;
}




//------------------------------
// Private helper functions
//------------------------------

static int adaptive_run(void *arkode_mem, N_Vector y, realtype T0, realtype Tf,
                        realtype dTout, realtype reltol, realtype abstol,
                        int rk_type, int nls_type, UserData &udata)
{
  // reused variables
  int retval;

  // Set tolerances
  retval = ARKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
  retval = ARKStepResStolerance(arkode_mem, abstol);
  if (check_retval(&retval, "ARKStepResStolerance", 1)) return 1;

  // Open output stream for results, output comment line
  FILE *UFID = fopen("ark_kpr_Mt_solution.txt","w");
  fprintf(UFID, "# t u v uerr verr\n");

  // output initial condition to disk
  fprintf(UFID," %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM " %.16" ESYM"\n",
          T0, NV_Ith_S(y,0), NV_Ith_S(y,1),
          SUNRabs(NV_Ith_S(y,0)-utrue(T0)),
          SUNRabs(NV_Ith_S(y,1)-vtrue(T0)));

  // Main time-stepping loop: calls ARKStepEvolve to perform integration,
  // then prints results. Stops when the final time has been reached
  int Nt = (int) ceil((Tf-T0)/dTout);
  realtype t = T0;
  realtype tout = T0+dTout;
  realtype uerr = ZERO;
  realtype verr = ZERO;
  realtype uerrtot = ZERO;
  realtype verrtot = ZERO;
  realtype errtot  = ZERO;
  cout << "\n        t           u           v       uerr      verr\n";
  cout << "   ------------------------------------------------------\n";
  printf("  %10.6" FSYM"  %10.6" FSYM"  %10.6" FSYM"  %.2" ESYM"  %.2" ESYM"\n",
         t, NV_Ith_S(y,0), NV_Ith_S(y,1), uerr, verr);

  for (int iout=0; iout<Nt; iout++) {

    // call integrator
    retval = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "ARKStepEvolve", 1)) break;

    // access/print solution and error
    uerr = SUNRabs(NV_Ith_S(y,0)-utrue(t));
    verr = SUNRabs(NV_Ith_S(y,1)-vtrue(t));
    printf("  %10.6" FSYM"  %10.6" FSYM"  %10.6" FSYM"  %.2" ESYM"  %.2" ESYM"\n",
           t, NV_Ith_S(y,0), NV_Ith_S(y,1), uerr, verr);
    fprintf(UFID," %.16" ESYM" %.16" ESYM" %.16" ESYM" %.16" ESYM" %.16" ESYM"\n",
            t, NV_Ith_S(y,0), NV_Ith_S(y,1), uerr, verr);
    uerrtot += uerr*uerr;
    verrtot += verr*verr;
    errtot += uerr*uerr + verr*verr;

    // successful solve: update time
    tout += dTout;
    tout = (tout > Tf) ? Tf : tout;
  }
  uerrtot = SUNRsqrt(uerrtot / Nt);
  verrtot = SUNRsqrt(verrtot / Nt);
  errtot = SUNRsqrt(errtot / Nt / 2);
  cout << "   ------------------------------------------------------\n";
  fclose(UFID);

  // Get integrator statistics
  long int nst, nst_a, nfe, nfi, nni, nnc, nje, nsetups, netf, nmset, nms, nMv;
  retval = ARKStepGetNumSteps(arkode_mem, &nst);
  if (check_retval(&retval, "ARKStepGetNumSteps", 1)) return 1;
  retval = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  if (check_retval(&retval, "ARKStepGetNumStepAttempts", 1)) return 1;
  retval = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  if (check_retval(&retval, "ARKStepGetNumRhsEvals", 1)) return 1;
  retval = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  if (check_retval(&retval, "ARKStepGetNumErrTestFails", 1)) return 1;
  retval = ARKStepGetNumMassSetups(arkode_mem, &nmset);
  if (check_retval(&retval, "ARKStepGetNumMassSetups", 1)) return 1;
  retval = ARKStepGetNumMassSolves(arkode_mem, &nms);
  if (check_retval(&retval, "ARKStepGetNumMassSolves", 1)) return 1;
  retval = ARKStepGetNumMassMult(arkode_mem, &nMv);
  if (check_retval(&retval, "ARKStepGetNumMassMult", 1)) return 1;
  if (rk_type < 2) {
    retval = ARKStepGetNonlinSolvStats(arkode_mem, &nni, &nnc);
    if (check_retval(&retval, "ARKStepGetNonlinSolvStats", 1)) return 1;
    retval = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
    if (check_retval(&retval, "ARKStepGetNumLinSolvSetups", 1)) return 1;
    if (nls_type == 0) {
      retval = ARKStepGetNumJacEvals(arkode_mem, &nje);
      if (check_retval(&retval, "ARKStepGetNumJacEvals", 1)) return 1;
    } else {
      nje = 0;
    }
  }

  // Print some final statistics and return
  cout << "\nFinal Solver Statistics:\n";
  cout << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n";
  cout << "   Total number of error test failures = " << netf << endl;
  cout << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << endl;
  cout << "   Total mass matrix setups = " << nmset << endl;
  cout << "   Total mass matrix solves = " << nms << endl;
  cout << "   Total mass times evals = " << nMv << endl;
  if (rk_type < 2) {
    cout << "   Total number of Jacobian evaluations = " << nje << "\n";
    cout << "   Total linear solver setups = " << nsetups << endl;
    cout << "   Total number of Nonlinear iterations = " << nni << endl;
    cout << "   Total number of Nonlinear convergence failures = " << nnc << endl;
  }
  cout << "   Errors: u = " << uerrtot << ", v = " << verrtot
       << ", total = " << errtot << endl;
  return(0);
}

static int check_order(void *arkode_mem, N_Vector y, realtype T0,
                       realtype Tf, int order, int rk_type, UserData &udata)
{
  // local variables
  int retval;
  size_t Nout = 20;
  realtype reltol = RCONST(1.e-9);
  realtype abstol = RCONST(1.e-12);
  realtype a11, a12, a21, a22, b1, b2;
  a11 = a12 = a21 = a22 = b1 = b2 = ZERO;

  // Tighten implicit solver to accommodate fixed step sizes
  retval = ARKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;
  retval = ARKStepResStolerance(arkode_mem, abstol);
  if (check_retval(&retval, "ARKStepResStolerance", 1)) return(1);
  retval = ARKStepSetMaxNumSteps(arkode_mem, 1000000);
  if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return(1);
  if (rk_type < 2) {
    retval = ARKStepSetJacEvalFrequency(arkode_mem, 1);
    if (check_retval(&retval, "ARKStepSetJacEvalFrequency", 1)) return 1;
    retval = ARKStepSetLSetupFrequency(arkode_mem, 1);
    if (check_retval(&retval, "ARKStepSetLSetupFrequency", 1)) return 1;
    retval = ARKStepSetMaxNonlinIters(arkode_mem, 20);
    if (check_retval(&retval, "ARKStepSetMaxNonlinIters", 1)) return 1;
  }

  // Set array of fixed step sizes to use, storage for corresponding errors/orders
  realtype hmax = Tf - T0;
  if (rk_type == 2)  hmax = ONE/SUNRabs(udata.G);
  realtype Nmin = SUNMAX((realtype) Nout, (realtype) ceil((Tf-T0)/hmax));
  vector<realtype> hvals = {(Tf-T0)/Nmin, (Tf-T0)/2/Nmin, (Tf-T0)/4/Nmin, (Tf-T0)/8/Nmin,
                            (Tf-T0)/16/Nmin, (Tf-T0)/32/Nmin, (Tf-T0)/64/Nmin, (Tf-T0)/128/Nmin};
  vector<realtype> errs(hvals.size());
  vector<realtype> orders(hvals.size()-1);

  // Loop over fixed step sizes
  cout << "\n Fixed-step runs:\n";
  cout << " -----------------------------------------------------\n";
  for (size_t ih=0; ih<hvals.size(); ih++) {

    // Reset ARKStep for this run
    retval = Ytrue(T0, y);
    if (check_retval(&retval, "Ytrue", 1)) return 1;
    retval = ARKStepReset(arkode_mem, T0, y);
    if (check_retval(&retval, "ARKStepReset", 1)) return 1;
    retval = ARKStepSetFixedStep(arkode_mem, hvals[ih]);
    if (check_retval(&retval, "ARKStepSetFixedStep", 1)) return 1;

    // Main time-stepping loop: run for Nout periods, accumulating overall error
    realtype t = T0;
    realtype dTout = (Tf-T0)/((realtype) Nout);
    realtype tout = T0+dTout;
    realtype uerr = ZERO;
    realtype verr = ZERO;
    errs[ih] = ZERO;
    for (size_t iout=0; iout<Nout; iout++) {

      // call integrator and update output time
      retval = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
      if (check_retval(&retval, "ARKStepEvolve", 1)) break;
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;

      // accumulate solution error
      uerr = SUNRabs(NV_Ith_S(y,0)-utrue(t));
      verr = SUNRabs(NV_Ith_S(y,1)-vtrue(t));
      errs[ih] += uerr*uerr + verr*verr;
    }
    errs[ih] = SUNMIN(ONE, SUNRsqrt(errs[ih] / ((realtype) Nout * 2)));
    a11 += 1;
    a12 += log(hvals[ih]);
    a21 += log(hvals[ih]);
    a22 += (log(hvals[ih])*log(hvals[ih]));
    b1 += log(errs[ih]);
    b2 += (log(errs[ih])*log(hvals[ih]));
    if (ih > 0) {
      orders[ih-1] = log(errs[ih]/errs[ih-1])/log(hvals[ih]/hvals[ih-1]);
      printf("   h = %.3" ESYM ",  error = %.3" ESYM ",  order = %.2" FSYM "\n",
        hvals[ih], errs[ih], orders[ih-1]);
    } else {
      printf("   h = %.3" ESYM ",  error = %.3" ESYM "\n", hvals[ih], errs[ih]);
    }
  }
  cout << " -----------------------------------------------------\n";

  // Compute/print overall estimated convergence rate
  realtype ord_avg, ord_max, ord_est, det;
  ord_avg = ord_max = ord_est = ZERO;
  for (size_t ih=1; ih<hvals.size(); ih++) {
    ord_avg += orders[ih-1];
    ord_max = SUNMAX(ord_max,orders[ih-1]);
  }
  ord_avg = ord_avg/((realtype) hvals.size()-1);
  det = a11*a22-a12*a21;
  ord_est = (a11*b2 - a21*b1)/det;
  printf(" Order: max = %.2" FSYM ",  avg = %.2" FSYM ",  overall = %.2" FSYM,
         ord_max, ord_avg, ord_est);

  // clean up and return
  if (ord_max < (order-RCONST(0.5))) {
    cout << " [FAILURE]\n";
    return(1);
  } else {
    cout << " [SUCCESS]\n";
    return(0);
  }
}

static realtype r(realtype t)
{
  return( RCONST(0.5)*cos(t) );
}
static realtype s(realtype t)
{
  return( sin(t) );
}
static realtype rdot(realtype t)
{
  return( -RCONST(0.5)*sin(t) );
}
static realtype sdot(realtype t)
{
  return( cos(t) );
}
static realtype utrue(realtype t)
{
  return( SUNRsqrt(ONE+r(t)) );
}
static realtype vtrue(realtype t)
{
  return( SUNRsqrt(TWO+s(t)) );
}
static int Ytrue(realtype t, N_Vector y)
{
  NV_Ith_S(y,0) = utrue(t);
  NV_Ith_S(y,1) = vtrue(t);
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

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  // Check if retval < 0
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1; }}

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}


//---- end of file ----
