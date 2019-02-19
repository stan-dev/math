/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This example is based on idaHeat2D_kry
 *                   example by Allan Taylor, Alan Hindmarsh and
 *                   Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, serial, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver Spgmr.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y). The
 * PDE is treated with central differences on a uniform M x M grid.
 * The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = M^2. Here M = 10.
 *
 * The system is solved with IDA using the Krylov linear solver
 * SPGMR. The preconditioner uses the diagonal elements of the
 * Jacobian only. Routines for preconditioning, required by
 * SPGMR, are supplied here. The constraints u >= 0 are posed
 * for all components. Output is taken at t = 0, .01, .02, .04,
 * ..., 10.24. Two cases are run -- with the Gram-Schmidt type
 * being Modified in the first case, and Classical in the second.
 * The second run uses IDAReInit.
 * -----------------------------------------------------------------
 * When running this example for unit testing, set the environemnt
 * variables OMP_PROC_BIND=false and OMP_NUM_THREADS=2 to run
 * without thread binding and with two OpenMP threads.
 * -----------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

#include <ida/ida.h>                   /* prototypes for IDA methods           */
#include <nvector/trilinos/SundialsTpetraVectorInterface.hpp>
#include <nvector/nvector_trilinos.h>
#include <sunlinsol/sunlinsol_spgmr.h> /* access to spgmr SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* definition of type realtype          */

/* Problem Constants */

#define NOUT  11
#define MGRID 10
#define NEQ   MGRID*MGRID
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define FOUR  RCONST(4.0)

/* User data type */

struct UserData {
  sunindextype mm;  /* number of grid points in one dimension */
  sunindextype neq; /* number of equations */
  realtype dx;
  realtype coeff;
  N_Vector pp;  /* vector of prec. diag. elements */
};

typedef Sundials::TpetraVectorInterface::vector_type vector_type;
typedef vector_type::scalar_type scalar_type;
typedef vector_type::mag_type mag_type;
typedef vector_type::global_ordinal_type global_ordinal_type;
typedef vector_type::local_ordinal_type local_ordinal_type;
typedef vector_type::node_type::memory_space memory_space;
typedef vector_type::node_type::execution_space execution_space;
typedef vector_type::map_type map_type;


/* Prototypes for functions called by IDA */

int resHeat(realtype tres, N_Vector uu, N_Vector up,
            N_Vector resval, void *user_data);

int PsetupHeat(realtype tt,
               N_Vector uu, N_Vector up, N_Vector rr,
               realtype c_j, void *prec_data);

int PsolveHeat(realtype tt,
               N_Vector uu, N_Vector up, N_Vector rr,
               N_Vector rvec, N_Vector zvec,
               realtype c_j, realtype delta, void *prec_data);

/* Prototypes for private functions */

static int SetInitialProfile(UserData *data, N_Vector uu, N_Vector up,
                             N_Vector res);
static void PrintHeader(realtype rtol, realtype atol);
static void PrintOutput(void *mem, realtype t, N_Vector uu);
static int check_retval(void *flagvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  void *mem;
  N_Vector uu, up, constraints, res;
  int retval, iout;
  realtype rtol, atol, t0, t1, tout, tret;
  long int netf, ncfn, ncfl;
  SUNLinearSolver LS;

  mem = NULL;
  uu = up = constraints = res = NULL;
  LS = NULL;

  using Teuchos::RCP;
  using Teuchos::rcp;

  /* Assign parameters in the user data structure. */

  UserData *data = new UserData();
  data->pp = NULL;
  if(check_retval((void *)data, "malloc", 2)) return(1);

  data->mm  = MGRID;
  data->neq = data->mm * data->mm;
  data->dx = ONE/(data->mm-ONE);
  data->coeff = ONE/(data->dx * data->dx);

  /* Start an MPI session */
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  /* Create Tpetra communicator */
  auto comm = Tpetra::getDefaultComm();

  /* Make partitioning easy */
  const sunindextype local_length  = data->neq;
  const sunindextype global_length = local_length;

  /* Choose zero-based (C-style) indexing. */
  const sunindextype index_base = 0;

  /* Construct an MPI Map */
  Teuchos::RCP<const map_type> testMap =
    Teuchos::rcp(new map_type (global_length, index_base, comm,
                               Tpetra::GloballyDistributed));

  /* Construct a Tpetra vector and return refernce counting pointer to it. */
  Teuchos::RCP<vector_type> rcpuu =
    Teuchos::rcp(new vector_type(testMap));

  if (comm->getSize() != 1) {
    printf("Warning: This test case works only with one MPI rank!\n");
    return -1;
  }

  /* Allocate N-vectors and the user data structure objects. */

  uu = N_VMake_Trilinos(rcpuu);
  if(check_retval((void *)uu, "N_VMake_Trilinos", 0)) return(1);

  up = N_VClone(uu);
  if(check_retval((void *)up, "N_VClone", 0)) return(1);

  res = N_VClone(uu);
  if(check_retval((void *)res, "N_VClone", 0)) return(1);

  constraints = N_VClone(uu);
  if(check_retval((void *)constraints, "N_VClone", 0)) return(1);

  data->pp = N_VClone(uu);
  if(check_retval((void *)data->pp, "N_VClone", 0)) return(1);

  /* Initialize uu, up. */

  SetInitialProfile(data, uu, up, res);

  /* Set constraints to all 1's for nonnegative solution values. */

  N_VConst(ONE, constraints);

  /* Assign various parameters. */

  t0   = ZERO;
  t1   = RCONST(0.01);
  rtol = ZERO;
  atol = RCONST(1.0e-3);

  /* Call IDACreate and IDAMalloc to initialize solution */

  mem = IDACreate();
  if(check_retval((void *)mem, "IDACreate", 0)) return(1);

  retval = IDASetUserData(mem, data);
  if(check_retval(&retval, "IDASetUserData", 1)) return(1);

  retval = IDASetConstraints(mem, constraints);
  if(check_retval(&retval, "IDASetConstraints", 1)) return(1);
  N_VDestroy(constraints);

  retval = IDAInit(mem, resHeat, t0, uu, up);
  if(check_retval(&retval, "IDAInit", 1)) return(1);

  retval = IDASStolerances(mem, rtol, atol);
  if(check_retval(&retval, "IDASStolerances", 1)) return(1);

  /* Create the linear solver SUNSPGMR with left preconditioning
     and the default Krylov dimension */
  LS = SUNLinSol_SPGMR(uu, PREC_LEFT, 0);
  if(check_retval((void *)LS, "SUNSPGMR", 0)) return(1);

  /* IDA recommends allowing up to 5 restarts (default is 0) */
  retval = SUNLinSol_SPGMRSetMaxRestarts(LS, 5);
  if(check_retval(&retval, "SUNSPGMRSetMaxRestarts", 1)) return(1);

  /* Attach the linear sovler */
  retval = IDASetLinearSolver(mem, LS, NULL);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  /* Set the preconditioner solve and setup functions */
  retval = IDASetPreconditioner(mem, PsetupHeat, PsolveHeat);
  if(check_retval(&retval, "IDASetPreconditioner", 1)) return(1);

  /* Print output heading. */
  PrintHeader(rtol, atol);

  /*
   * -------------------------------------------------------------------------
   * CASE I
   * -------------------------------------------------------------------------
   */

  /* Print case number, output table heading, and initial line of table. */

  printf("\n\nCase 1: gsytpe = MODIFIED_GS\n");
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps\n" );
  printf("----------------------------------------------------------------------\n");

  /* Loop over output times, call IDASolve, and print results. */

  for (tout = t1,iout = 1; iout <= NOUT ; iout++, tout *= TWO) {
    retval = IDASolve(mem, tout, &tret, uu, up, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1)) return(1);
    PrintOutput(mem, tret, uu);
  }

  /* Print remaining counters. */

  retval = IDAGetNumErrTestFails(mem, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1);

  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1);

  retval = IDAGetNumLinConvFails(mem, &ncfl);
  check_retval(&retval, "IDAGetNumLinConvFails", 1);

  printf("\nError test failures            = %ld\n", netf);
  printf("Nonlinear convergence failures = %ld\n", ncfn);
  printf("Linear convergence failures    = %ld\n", ncfl);

  /*
   * -------------------------------------------------------------------------
   * CASE II
   * -------------------------------------------------------------------------
   */

  /* Re-initialize uu, up. */

  SetInitialProfile(data, uu, up, res);

  /* Re-initialize IDA and SPGMR */

  retval = IDAReInit(mem, t0, uu, up);
  if(check_retval(&retval, "IDAReInit", 1)) return(1);

  retval = SUNSPGMRSetGSType(LS, CLASSICAL_GS);
  if(check_retval(&retval, "SUNSPGMRSetGSType",1)) return(1);

  /* Print case number, output table heading, and initial line of table. */

  printf("\n\nCase 2: gstype = CLASSICAL_GS\n");
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps\n" );
  printf("----------------------------------------------------------------------\n");

  /* Loop over output times, call IDASolve, and print results. */

  for (tout = t1,iout = 1; iout <= NOUT ; iout++, tout *= TWO) {
    retval = IDASolve(mem, tout, &tret, uu, up, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1)) return(1);
    PrintOutput(mem, tret, uu);
  }

  /* Print remaining counters. */

  retval = IDAGetNumErrTestFails(mem, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1);

  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1);

  retval = IDAGetNumLinConvFails(mem, &ncfl);
  check_retval(&retval, "IDAGetNumLinConvFails", 1);

  printf("\nError test failures            = %ld\n", netf);
  printf("Nonlinear convergence failures = %ld\n", ncfn);
  printf("Linear convergence failures    = %ld\n", ncfl);

  /* Free Memory */

  IDAFree(&mem);
  SUNLinSolFree(LS);

  N_VDestroy(uu);
  N_VDestroy(up);
  N_VDestroy(res);

  N_VDestroy(data->pp);
  delete data;

  /* Free Tpetra vector */
  rcpuu = Teuchos::null;

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * resHeat: heat equation system residual function (user-supplied)
 * This uses 5-point central differencing on the interior points, and
 * includes algebraic equations for the boundary values.
 * So for each interior point, the residual component has the form
 *    res_i = u'_i - (central difference)_i
 * while for each boundary point, it is res_i = u_i.
 */
int resHeat(realtype tt,
            N_Vector uu, N_Vector up, N_Vector rr,
            void *user_data)
{
  UserData *data = reinterpret_cast<UserData*>(user_data);
  sunindextype mm;
  realtype coeff;

  Teuchos::RCP<vector_type> uutp = N_VGetVector_Trilinos(uu);
  Teuchos::RCP<vector_type> uptp = N_VGetVector_Trilinos(up);
  Teuchos::RCP<vector_type> rrtp = N_VGetVector_Trilinos(rr);

  const auto uu_2d = uutp->getLocalView<memory_space>();
  const auto uu_1d = Kokkos::subview (uu_2d, Kokkos::ALL(), 0);
  const auto up_2d = uptp->getLocalView<memory_space>();
  const auto up_1d = Kokkos::subview (up_2d, Kokkos::ALL(), 0);
  auto rr_2d = rrtp->getLocalView<memory_space>();
  auto rr_1d = Kokkos::subview (rr_2d, Kokkos::ALL(), 0);

  coeff = data->coeff;
  mm    = data->mm;

  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, mm*mm),
    KOKKOS_LAMBDA(const local_ordinal_type &loc)
    {
      sunindextype i = loc % mm;
      sunindextype j = loc / mm;
      if (j==0 || j==mm-1 || i==0 || i==mm-1) {
        /* Initialize rr to uu, to take care of boundary equations. */
        rr_1d(loc) = uu_1d(loc);
      } else {
        /* Loop over interior points; set res = up - (central difference). */
        realtype dif1 = uu_1d(loc-1)  + uu_1d(loc+1)  - TWO * uu_1d(loc);
        realtype dif2 = uu_1d(loc-mm) + uu_1d(loc+mm) - TWO * uu_1d(loc);
        rr_1d(loc) = up_1d(loc) - coeff * ( dif1 + dif2 );
      }
    }
  );

  return(0);
}


/*
 * PsetupHeat: setup for diagonal preconditioner for idaHeat2D_kry.
 *
 * The optional user-supplied functions PsetupHeat and
 * PsolveHeat together must define the left preconditoner
 * matrix P approximating the system Jacobian matrix
 *                   J = dF/du + cj*dF/du'
 * (where the DAE system is F(t,u,u') = 0), and solve the linear
 * systems P z = r.   This is done in this case by keeping only
 * the diagonal elements of the J matrix above, storing them as
 * inverses in a vector pp, when computed in PsetupHeat, for
 * subsequent use in PsolveHeat.
 *
 * In this instance, only cj and data (user data structure, with
 * pp etc.) are used from the PsetupdHeat argument list.
 */
int PsetupHeat(realtype tt,
               N_Vector uu, N_Vector up, N_Vector rr,
               realtype c_j, void *prec_data)
{
  UserData *data = reinterpret_cast<UserData*>(prec_data);
  sunindextype mm;

  mm = data->mm;
  realtype coeff = data->coeff;

  Teuchos::RCP<vector_type> pptp = N_VGetVector_Trilinos(data->pp);

  auto pp_2d = pptp->getLocalView<memory_space>();
  auto pp_1d = Kokkos::subview (pp_2d, Kokkos::ALL(), 0);


  Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mm*mm),
    KOKKOS_LAMBDA (const local_ordinal_type &loc)
    {
      sunindextype i = loc % mm;
      sunindextype j = loc / mm;
      if (j==0 || j==mm-1 || i==0 || i==mm-1) {
        /* Set ppv to one, to take care of boundary equations. */
        pp_1d(loc) = ONE;
      } else {
        /* Loop over interior points; ppv_i = 1/J_ii */
        pp_1d(loc) = ONE/(c_j + FOUR*coeff);
      }
    }
  );

  return(0);
}


/*
 * PsolveHeat: solve preconditioner linear system.
 * This routine multiplies the input vector rvec by the vector pp
 * containing the inverse diagonal Jacobian elements, returning the
 * result in zvec.
 */
int PsolveHeat(realtype tt,
               N_Vector uu, N_Vector up, N_Vector rr,
               N_Vector rvec, N_Vector zvec,
               realtype c_j, realtype delta, void *prec_data)
{
  UserData *data = reinterpret_cast<UserData*>(prec_data);
  N_VProd(data->pp, rvec, zvec);
  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * SetInitialProfile: routine to initialize u and up vectors.
 */
static int SetInitialProfile(UserData *data, N_Vector uu, N_Vector up,
                             N_Vector res)
{
  sunindextype mm, mm1, i, j;
  realtype xfact, yfact;

  mm = data->mm;

  Teuchos::RCP<vector_type> u = N_VGetVector_Trilinos(uu);
  Teuchos::RCP<vector_type> p = N_VGetVector_Trilinos(up);

  /* Sync the host with the device if needed */
  u->sync<Kokkos::HostSpace>();
  auto u_2d = u->getLocalView<Kokkos::HostSpace>();
  auto u_1d = Kokkos::subview(u_2d, Kokkos::ALL(), 0);

  u->modify<Kokkos::HostSpace>();

  /* Initialize uu on all grid points. */
  for (j = 0; j < mm; j++) {
    yfact = data->dx * j;
    for (i = 0; i < mm; i++) {
      xfact = data->dx * i;
      u_1d(mm*j + i) = RCONST(16.0) * xfact * (ONE - xfact) * yfact * (ONE - yfact);
    }
  }

  /* Sync the device with the host */
  u->sync<memory_space>();

  /* Initialize up vector to 0. */
  N_VConst(ZERO, up);

  /* resHeat sets res to negative of ODE RHS values at interior points. */
  resHeat(ZERO, uu, up, res, data);

  /* Copy -res into up to get correct interior initial up values. */
  N_VScale(-ONE, res, up);

  /* Set up at boundary points to zero. */
  auto p_2d = p->getLocalView<memory_space>();
  auto p_1d = Kokkos::subview(p_2d, Kokkos::ALL(), 0);
  mm1 = mm - 1;

  Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mm),
    KOKKOS_LAMBDA (const local_ordinal_type &loc)
    {
      sunindextype i = loc % mm;
      sunindextype j = loc / mm;
      if (j==0 || j==mm1 || i==0 || i==mm1) {
        p_1d(loc) = ZERO;
      }
    }
  );

  return(0);
}


/*
 * Print first lines of output (problem description)
 */
static void PrintHeader(realtype rtol, realtype atol)
{
  printf("\nidaHeat2D_kry_tpetra: Heat equation, serial example problem for IDA \n");
  printf("         Discretized heat equation on 2D unit square. \n");
  printf("         Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("         Mesh dimensions: %d x %d", MGRID, MGRID);
  printf("        Total system size: %d\n\n", NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Constraints set to force all solution components >= 0. \n");
  printf("Linear solver: SPGMR, preconditioner using diagonal elements. \n");
}


/*
 * PrintOutput: print max norm of solution and current solver statistics
 */
static void PrintOutput(void *mem, realtype t, N_Vector uu)
{
  realtype hused, umax;
  long int nst, nni, nje, nre, nreLS, nli, npe, nps;
  int kused, retval;

  umax = N_VMaxNorm(uu);

  retval = IDAGetLastOrder(mem, &kused);
  check_retval(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  check_retval(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_retval(&retval, "IDAGetNumResEvals", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_retval(&retval, "IDAGetLastStep", 1);
  retval = IDAGetNumJtimesEvals(mem, &nje);
  check_retval(&retval, "IDAGetNumJtimesEvals", 1);
  retval = IDAGetNumLinIters(mem, &nli);
  check_retval(&retval, "IDAGetNumLinIters", 1);
  retval = IDAGetNumLinResEvals(mem, &nreLS);
  check_retval(&retval, "IDAGetNumLinResEvals", 1);
  retval = IDAGetNumPrecEvals(mem, &npe);
  check_retval(&retval, "IDAGetPrecEvals", 1);
  retval = IDAGetNumPrecSolves(mem, &nps);
  check_retval(&retval, "IDAGetNumPrecSolves", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2Le  %3ld %3ld\n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#else
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#endif
}


/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */
static int check_retval(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  }

  return(0);
}
