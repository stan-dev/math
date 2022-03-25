/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * The following is a simple example problem based off of cvRoberts_klu.c. We
 * simulate a scenario where a set of independent ODEs are grouped together to
 * form a larger system. For simplicity, each set of ODEs is the same problem.
 * The problem is from chemical kinetics, and consists of the following three
 * rate equations:
 *
 *   dy1/dt = -.04 * y1 + 1.0e4 * y2 * y3
 *   dy2/dt = -dy1/dt - dy3/dt
 *   dy3/dt = 3.0e7 * y2 * y2
 *
 * on the interval from t = 0.0 to t = 4.0e10, with initial conditions:
 *
 *   y1 = 1.0, y2 = 0, and y3 = 0.
 *
 * The problem is stiff. By default this program solves the problem with the BDF
 * methods, a Newton iteration, user-supplied Jacobian routine, and since the
 * grouping of the independent systems results in a block diagonal linear
 * system, with the oneMKL SUNLinearSolver. Alternatively, the SPGMR linear
 * solver may be used with a Jacobi preconditioner. The problem uses a scalar
 * relative tolerance and a vector absolute tolerance. Output is printed in
 * decades from t = 0.1 to 1.0e6. Run statistics (optional outputs) are printed
 * at the end.
 *
 * The program takes three optional argument: the number of independent ODE
 * systems (default 100), a flag to use a direct (1, default) or iterative (0)
 * linear solver, and a flag to enable (1, default) or disable (0) solution
 * output:
 *
 *   ./cvRoberts_blockdiag_onemkl [number of groups] [solver type] [output]
 *
 * This problem is comparable to the cvRoberts_block_klu.c example.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <iostream>
#include <chrono>

#include <cvode/cvode.h>                      // access to CVODE fcts., consts.
#include <nvector/nvector_sycl.h>             // access the SYCL NVector
#include <sunmemory/sunmemory_sycl.h>         // access the SYCL Memory helper
#include <sunmatrix/sunmatrix_onemkldense.h>  // access the oneMKL SUNMatrix
#include <sunlinsol/sunlinsol_onemkldense.h>  // access the oneMKL SUNLinearSolver
#include <sunlinsol/sunlinsol_spgmr.h>        // access the GMRES SUNLinearSolver

using namespace std;

// Problem Constants

#define GROUPSIZE 3            // number of equations per group
#define Y1    RCONST(1.0)      // initial y components
#define Y2    RCONST(0.0)
#define Y3    RCONST(0.0)
#define RTOL  RCONST(1.0e-4)   // scalar relative tolerance
#define ATOL1 RCONST(1.0e-8)   // vector absolute tolerance components
#define ATOL2 RCONST(1.0e-14)
#define ATOL3 RCONST(1.0e-6)
#define T0    RCONST(0.0)      // initial time
#define T1    RCONST(0.1)      // first output time
#define TMULT RCONST(10.0)     // output time factor
#define NOUT  10               // number of output times

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)

// Functions Called by the Solver
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int PSolve(realtype t, N_Vector u, N_Vector f, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,
                  void *user_data);

// Private function to output results
static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3);

// Private function to print final statistics
static void PrintFinalStats(void *cvode_mem, bool direct);

// Private function to check function return values
static int check_retval(void *returnvalue, const char *funcname, int opt);

// User data structure
typedef struct
{
  sycl::queue* myQueue;
  int ngroups;
  int neq;
} UserData;

/* ---------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  // SUNDIALS simulation context
  sundials::Context sunctx;

  // return value flag
  int retval;

  // Parse command line arguments

  // Number of ODE groups
  int ngroups = 100;
  if (argc > 1) ngroups = atoi(argv[1]);

  // Use a direct or iterative linear sovler
  bool direct = true;
  if (argc > 2) direct = (atoi(argv[2])) ? true : false;

  // Write the solution to the screen
  bool output = true;
  if (argc > 3) output = (atoi(argv[3])) ? true : false;

  // Create an in-order GPU queue
  sycl::gpu_selector selector;
  sycl::queue myQueue(selector,
                      sycl::property_list{sycl::property::queue::in_order{}});

  sycl::device dev = myQueue.get_device();
  cout << "Running on "
       << (dev.get_info<sycl::info::device::name>())
       << endl;

  // Total number of equations
  sunindextype neq = ngroups * GROUPSIZE;

  // Set user data values
  UserData udata;
  udata.myQueue = &myQueue;
  udata.ngroups = ngroups;
  udata.neq     = neq;

  // Create the SYCL memory helper
  SUNMemoryHelper memhelper = SUNMemoryHelper_Sycl(sunctx);
  if (check_retval((void *)memhelper, "SUNMemoryHelper_Sycl", 0)) return 1;

  // Create SYCL vector for state and absolute tolerances
  N_Vector y = N_VNew_Sycl(neq, &myQueue, sunctx);
  if (check_retval((void *)y, "N_VNew", 0)) return 1;

  N_Vector abstol = N_VClone(y);
  if (check_retval((void *)abstol, "N_VClone", 0)) return 1;

  // Initialize y
  realtype* ydata = N_VGetArrayPointer(y);
  for (sunindextype groupj = 0; groupj < neq; groupj += GROUPSIZE)
  {
    ydata[groupj]     = Y1;
    ydata[groupj + 1] = Y2;
    ydata[groupj + 2] = Y3;
  }
  N_VCopyToDevice_Sycl(y);

  // Set the scalar relative tolerance
  realtype reltol = RTOL;

  // Set the vector absolute tolerance
  realtype* abstol_data = N_VGetArrayPointer(abstol);
  for (sunindextype groupj = 0; groupj < neq; groupj += GROUPSIZE)
  {
    abstol_data[groupj]     = ATOL1;
    abstol_data[groupj + 1] = ATOL2;
    abstol_data[groupj + 2] = ATOL3;
  }
  N_VCopyToDevice_Sycl(abstol);

  // Create CVODE with BDF methods
  void* cvode_mem = CVodeCreate(CV_BDF, NULL);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return 1;

  // Initialize CVODE
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return 1;

  // Call CVodeSetUserData to attach the user data structure
  retval = CVodeSetUserData(cvode_mem, &udata);
  if (check_retval(&retval, "CVodeSetUserData", 1)) return 1;

  // Set tolerances
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return 1;

  // Create and attach linear solver
  SUNMatrix       A  = NULL;
  SUNLinearSolver LS = NULL;

  if (direct)
  {
    // Create SUNMatrix for use in linear solves
    A = SUNMatrix_OneMklDenseBlock(ngroups, GROUPSIZE, GROUPSIZE,
                                   SUNMEMTYPE_DEVICE, memhelper, &myQueue,
                                   sunctx);
    if (check_retval((void *)A, "SUNMatrix_OneMklDenseBlock", 0)) return 1;

    // Create the SUNLinearSolver object for use by CVode
    LS = SUNLinSol_OneMklDense(y, A, sunctx);
    if (check_retval((void *)LS, "SUNLinSol_OneMklDense", 0)) return 1;

    // Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

    // Set the user-supplied Jacobian routine Jac
    retval = CVodeSetJacFn(cvode_mem, Jac);
    if (check_retval(&retval, "CVodeSetJacFn", 1)) return 1;
  }
  else
  {
    // Create SPGMR solver
    LS = SUNLinSol_SPGMR(y, SUN_PREC_RIGHT, 10, sunctx);
    if (check_retval(&retval, "SUNLinSol_SPGMR", 1)) return 1;

    // Call CVodeSetLinearSolver to attach the linear solver to CVode
    retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

    // Attach preconditioner
    retval = CVodeSetPreconditioner(cvode_mem, NULL, PSolve);
    if (check_retval(&retval, "CVodeSetPreconditioner", 1)) return 1;
  }

  // Loop over output times and print solution
  printf("\nGroup of independent 3-species kinetics problems\n");
  printf("  number of groups = %d\n", ngroups);
  if (direct)
    printf("  using direct linear solver\n");
  else
    printf("  using iterative linear solver\n");
  if (output)
    printf("  output enabled\n");
  else
    printf("  output disabled\n");

  int      iout = 0;
  realtype tout = T1;
  realtype t;

  // Start timer
  chrono::time_point<chrono::steady_clock> tstart = chrono::steady_clock::now();

  while(1)
  {
    // Evolve in time
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1)) break;

    if (output)
    {
      // Copy solution to host for output
      N_VCopyFromDevice_Sycl(y);
      for (sunindextype groupj = 0; groupj < ngroups; groupj += 10)
      {
        printf("group %ld: ", (long int) groupj);
        PrintOutput(t, ydata[GROUPSIZE * groupj],
                    ydata[1 + GROUPSIZE * groupj],
                    ydata[2 + GROUPSIZE * groupj]);
      }
    }

    // Update output counter and output time
    iout++;
    tout *= TMULT;

    // Stop after NOUT outputs
    if (iout == NOUT) break;
  }

  // Stop timer
  myQueue.wait();
  chrono::time_point<chrono::steady_clock> tstop = chrono::steady_clock::now();

  // Print some final statistics
  PrintFinalStats(cvode_mem, direct);

  // Print evoltuion time
  cout << "Evolution time: "
       << chrono::duration<double>(tstop - tstart).count() << endl;

  // Free objects and integrator
  N_VDestroy(y);
  N_VDestroy(abstol);
  SUNMatDestroy(A);
  SUNLinSolFree(LS);
  CVodeFree(&cvode_mem);
  SUNMemoryHelper_Destroy(memhelper);

  return 0;
}


/* ---------------------------------------------------------------------------
 * Functions called by the solver
 * ---------------------------------------------------------------------------*/


// Compute the right-hand side function, ydot = f(t, y).
static int f(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  UserData* udata    = (UserData*) user_data;
  realtype* ydata    = N_VGetDeviceArrayPointer(y);
  realtype* ydotdata = N_VGetDeviceArrayPointer(ydot);

  const size_t       ngroups = static_cast<size_t>(udata->ngroups);
  const sunindextype N       = GROUPSIZE;

  udata->myQueue->submit([&](sycl::handler& h)
  {
    h.parallel_for(sycl::range{ngroups}, [=](sycl::id<1> idx)
    {
      sunindextype groupj = idx[0];

      realtype y1 = ydata[N * groupj];
      realtype y2 = ydata[N * groupj + 1];
      realtype y3 = ydata[N * groupj + 2];

      realtype yd1 = RCONST(-0.04) * y1 + RCONST(1.0e4) * y2 * y3;
      realtype yd3 = RCONST(3.0e7) * y2 * y2;

      ydotdata[N * groupj]     = yd1;
      ydotdata[N * groupj + 1] = -yd1 - yd3;
      ydotdata[N * groupj + 2] = yd3;
    });
  });

  udata->myQueue->wait_and_throw();

  return 0;
}


// Compute the right-hand side Jacobian, J(t,y) = df/dy.
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData* udata = (UserData*) user_data;
  realtype* Jdata = SUNMatrix_OneMklDense_Data(J);
  realtype* ydata = N_VGetDeviceArrayPointer(y);

  const size_t       ngroups = static_cast<size_t>(udata->ngroups);
  const sunindextype N       = GROUPSIZE;
  const sunindextype NN      = GROUPSIZE * GROUPSIZE;

  udata->myQueue->submit([&](sycl::handler& h)
  {
    h.parallel_for(sycl::range{ngroups}, [=](sycl::id<1> idx)
    {
      sunindextype groupj = idx[0];

      // get y values
      realtype y2 = ydata[N * groupj + 1];
      realtype y3 = ydata[N * groupj + 2];

      // first col of block
      Jdata[NN * groupj]     = RCONST(-0.04);
      Jdata[NN * groupj + 1] = RCONST(0.04);
      Jdata[NN * groupj + 2] = ZERO;

      // second col of block
      Jdata[NN * groupj + 3] = RCONST(1.0e4) * y3;
      Jdata[NN * groupj + 4] = RCONST(-1.0e4) * y3 - RCONST(6.0e7) * y2;
      Jdata[NN * groupj + 5] = RCONST(6.0e7) * y2;

      // third col of block
      Jdata[NN * groupj + 6] = RCONST(1.0e4) * y2;
      Jdata[NN * groupj + 7] = RCONST(-1.0e4) * y2;
      Jdata[NN * groupj + 8] = ZERO;
    });
  });

  udata->myQueue->wait_and_throw();

  return 0;
}


static int PSolve(realtype t, N_Vector y, N_Vector f, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,
                  void *user_data)
{
  UserData* udata = (UserData*) user_data;
  realtype* ydata = N_VGetDeviceArrayPointer(y);
  realtype* rdata = N_VGetDeviceArrayPointer(r);
  realtype* zdata = N_VGetDeviceArrayPointer(z);

  const size_t       ngroups = static_cast<size_t>(udata->ngroups);
  const sunindextype N       = GROUPSIZE;

  udata->myQueue->submit([&](sycl::handler& h)
  {
    h.parallel_for(sycl::range{ngroups}, [=](sycl::id<1> idx)
    {
      sunindextype groupj = idx[0];
      sunindextype i0     = N * groupj;
      sunindextype i1     = N * groupj + 1;
      sunindextype i2     = N * groupj + 2;

      // Solve (I - gamma J) z = r
      //
      // [ 1 + a     -b      -c ] [ z0 ]   [ r0 ]
      // [  -a    1 + b + d   c ] [ z1 ] = [ r1 ]
      // [   0       -d       1 ] [ z2 ]   [ r2 ]

      // get y values
      realtype y2 = ydata[i1];
      realtype y3 = ydata[i2];

      // set matrix values
      realtype a = gamma * RCONST(0.04);
      realtype b = gamma * RCONST(1.0e4) * y3;
      realtype c = gamma * RCONST(1.0e4) * y2;
      realtype d = gamma * RCONST(6.0e7) * y2;

      // Initial Jacobi iteration with zero guess

      // z0 = r0 / (1 + a)
      zdata[i0] = rdata[i0] / (ONE + a);

      // z1 = r1 / (1 + b + d)
      zdata[i1] = rdata[i1] / (1 + b + d);

      // z2 = r2 + d
      zdata[i2] = rdata[i2];

      // Subsequent Jacobi iterations

      for (int i = 1; i < 10; ++i)
      {
        realtype z0 = zdata[i0];
        realtype z1 = zdata[i1];
        realtype z2 = zdata[i2];

        // z0 = (r0 + b * z1 + x * z2 / (1 + a)
        zdata[i0] = (rdata[i0] + b * z1 + c * z2) / (ONE + a);

        // z1 = r1 / (1 + b + d)
        zdata[i1] = (rdata[i1] + a * z0 - c * z2) / (1 + b + d);

        // z2 = r2 + d
        zdata[i2] = (rdata[i2] + d * z1);
      }
    });
  });

  return 0;
}

/* ---------------------------------------------------------------------------
 * Private helper functions
 * ---------------------------------------------------------------------------*/


// Output solution
static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}


// Get and print some final statistics
static void PrintFinalStats(void *cvode_mem, bool direct)
{
  int retval;

  printf("\nFinal Statistics:\n");

  // CVODE stats
  long int nst, nfe, netf;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);

  cout << "Time steps:       " << nst  << "\n";
  cout << "RHS evals:        " << nfe  << "\n";
  cout << "Error test fails: " << netf << "\n\n";

  // Nonlinear solver stats
  long int nni, ncfn;

  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  cout << "NLS iters: " << nni  << "\n";
  cout << "NLS fails: " << ncfn << "\n\n";

  // Linear solver stats
  if (direct)
  {
    long int nsetups, nje;
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
    retval = CVodeGetNumJacEvals(cvode_mem, &nje);
    check_retval(&retval, "CVodeGetNumJacEvals", 1);

    cout << "LS setups: " << nsetups << "\n";
    cout << "Jac evals: " << nje     << "\n\n";
  }
  else
  {
    long int nli, ncfl, nfeLS, nps;
    retval = CVodeGetNumLinIters(cvode_mem, &nli);
    check_retval(&retval, "CVodeGetNumLinIters", 1);
    retval = CVodeGetNumLinConvFails(cvode_mem, &ncfl);
    check_retval(&retval, "CVodeGetNumLinConvFails", 1);
    retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
    check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);
    retval = CVodeGetNumPrecSolves(cvode_mem, &nps);
    check_retval(&retval, "CVodeGetNumPrecSolves", 1);

    cout << "LS iters:     " << nli   << "\n";
    cout << "LS fails:     " << ncfl  << "\n";
    cout << "LS RHS evals: " << nfeLS << "\n";
    cout << "P solves:     " << nps   << "\n\n";
  }

  return;
}


/* Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer */
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  // Check if retval < 0
  else if (opt == 1)
  {
    retval = (int *) returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return 1;
    }
  }

  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1;
  }

  return 0;
}
