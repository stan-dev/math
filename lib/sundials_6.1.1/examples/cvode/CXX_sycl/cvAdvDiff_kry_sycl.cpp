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
 * Example problem:
 *
 * The following is a simple example problem with a banded Jacobian, with the
 * program for its solution by CVODE. The problem is the semi-discrete form of
 * the advection-diffusion equation in 2-D:
 *
 *   u_t = u_xx + u_yy + 0.5 u_x
 *
 * on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time interval 0 <= t <= 1.
 * Homogeneous Dirichlet boundary conditions are posed, and the initial
 * condition is
 *
 *   u(x,y,0) = x (2-x) y (1-y) exp(5xy).
 *
 * The PDE is discretized on a uniform MX+2 by MY+2 grid with central
 * differencing, and with boundary values eliminated, leaving an ODE system of
 * size NEQ = MX*MY.
 *
 * This program solves the problem with the BDF method, Newton iteration with
 * the GMRES linear solver, and a user-supplied Jacobian-vector product routine.
 * It uses scalar relative and absolute tolerances.
 *
 * Output is printed at t = .1, .2, ..., 1. Run statistics (optional outputs)
 * are printed at the end.
 * ---------------------------------------------------------------------------*/

#include <cstdlib>
#include <cmath>
#include <iostream>

#include <cvode/cvode.h>               // access CVODE fcts., consts.
#include <nvector/nvector_sycl.h>      // access the SYCL NVector
#include <sunlinsol/sunlinsol_spgmr.h> // access the SPGMR SUNLinearSolver

// Real Constants

#define ATOL  RCONST(1.0e-5) // scalar absolute tolerance
#define T0    RCONST(0.0)    // initial time
#define T1    RCONST(0.1)    // first output time
#define DTOUT RCONST(0.1)    // output time increment
#define NOUT  10             // number of output times

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)
#define FIVE RCONST(5.0)

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#if defined(SUNDIALS_INT64_T)
#define DSYM "ld"
#else
#define DSYM "d"
#endif

// User data sstructure contains model and discretization parameters
struct UserData
{
  sycl::queue* myQueue = NULL;              // SYCL queue
  sunindextype MX      = 10;                // interior nodes in the x-direction
  sunindextype MY      = 5;                 // interior nodes in the y-direction
  sunindextype NEQ     = MX * MY;           // number of equations
  realtype     xmax    = RCONST(2.0);       // x-domain boundary
  realtype     ymax    = RCONST(1.0);       // y-domain boundary
  realtype     dx      = xmax / (MX + 1);   // x-direction mesh spacing
  realtype     dy      = ymax / (MY + 1);   // y-directino mesh spacing
  realtype     hdcoef  = ONE / (dx * dx);   // x-diffusion
  realtype     vdcoef  = ONE / (dy * dy);   // y-diffusion
  realtype     hacoef  = HALF / (TWO * dx); // x-advection
};

// Functions Called by the Solver
static int f(realtype t, N_Vector u, N_Vector udot, void* user_data);
static int jtv(N_Vector v, N_Vector Jv, realtype t,
               N_Vector u, N_Vector fu,
               void* user_data, N_Vector tmp);

// Private Helper Functions
static void PrintHeader(realtype reltol, realtype abstol, realtype umax,
                        UserData* data);
static void PrintOutput(realtype t, realtype umax, long int nst);
static void PrintFinalStats(void* cvode_mem);

// Private function to check function return values
static int check_retval(void* returnvalue, const char *funcname, int opt);

/* ---------------------------------------------------------------------------
 * Main Program
 * ---------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
  // SUNDIALS simulation context
  sundials::Context sunctx;

  // return flag value
  int retval;

  // Create an in-order GPU queue
  sycl::gpu_selector selector;
  sycl::queue myQueue(selector,
                      sycl::property_list{sycl::property::queue::in_order{}});

  sycl::device dev = myQueue.get_device();
  std::cout << "Running on "
            << (dev.get_info<sycl::info::device::name>())
            << std::endl;

  // Create user data and set queue
  UserData data;
  data.myQueue = &myQueue;

  // Create a SYCL vector
  N_Vector u = N_VNew_Sycl(data.NEQ, &myQueue, sunctx);
  if (check_retval((void*)u, "N_VNew_Sycl", 0)) return 1;

  // Extract host pointer to solution vector data on the host
  realtype* udata = N_VGetArrayPointer(u);

  // Load initial profile into u vector
  for (sunindextype tid = 0; tid < data.NEQ; tid++)
  {
    sunindextype i = tid / data.MY; // x-node index
    sunindextype j = tid % data.MY; // y-node index

    realtype x = (i + 1) * data.dx; // x location
    realtype y = (j + 1) * data.dy; // y location

    udata[tid] =
      x * (data.xmax - x) * y * (data.ymax - y) * std::exp(FIVE * x * y);
  }

  // Copy data to device
  N_VCopyToDevice_Sycl(u);

  // Create CVODE and specify the Backward Differentiation Formula
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return 1;

  // Specify the right hand side function in f(t,u), initial condition (t0, u0)
  retval = CVodeInit(cvode_mem, f, T0, u);
  if (check_retval(&retval, "CVodeInit", 1)) return 1;

  // Specify the scalar relative tolerance and scalar absolute tolerance
  realtype reltol = ZERO;
  realtype abstol = ATOL;
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) return 1;

  // Set the pointer to user-defined data
  retval = CVodeSetUserData(cvode_mem, &data);
  if (check_retval(&retval, "CVodeSetUserData", 1)) return 1;

  // Create SPGMR solver without preconditioning and default Krylov dimension
  SUNLinearSolver LS = SUNLinSol_SPGMR(u, SUN_PREC_NONE, 0, sunctx);
  if (check_retval(&retval, "SUNLinSol_SPGMR", 1)) return 1;

  // Attach the linear sovler to CVODE
  retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

  // Set the Jacobian-times-vector function
  retval = CVodeSetJacTimes(cvode_mem, NULL, jtv);
  if (check_retval(&retval, "CVodeSetJacTimesVecFn", 1)) return 1;

  // In loop over output points: call CVODE, print results, test for errors
  realtype umax = N_VMaxNorm(u);
  PrintHeader(reltol, abstol, umax, &data);

  realtype tout = T1; // output time
  realtype t;         // CVODE return time
  long int nst;       // number of time steps

  for (int iout = 0; iout < NOUT; iout++)
  {
    // Advance in time
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1)) break;

    // Output status
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    if (check_retval(&retval, "CVodeGetNumSteps", 1)) break;

    umax = N_VMaxNorm(u);
    PrintOutput(t, umax, nst);

    // Update output time
    tout += DTOUT;
  }

  PrintFinalStats(cvode_mem);  // Print some final statistics

  N_VDestroy(u);          // Free the u vector
  CVodeFree(&cvode_mem);  // Free the integrator memory
  SUNLinSolFree(LS);      // Free linear solver memory

  return 0;
}


/* ---------------------------------------------------------------------------
 * Functions called by the solver
 * ---------------------------------------------------------------------------*/

// Compute the ODE right-hand side function f(t,u).
static int f(realtype t, N_Vector u, N_Vector udot, void* user_data)
{
  UserData* data = static_cast<UserData*>(user_data);

  // Extract needed constants from data
  const size_t   MX    = static_cast<size_t>(data->MX);
  const size_t   MY    = static_cast<size_t>(data->MY);
  const realtype hordc = data->hdcoef;
  const realtype horac = data->hacoef;
  const realtype verdc = data->vdcoef;

  // Extract pointers to vector data
  const realtype* udata  = N_VGetDeviceArrayPointer(u);
  realtype*       dudata = N_VGetDeviceArrayPointer(udot);

  data->myQueue->submit([&](sycl::handler& h)
  {
    h.parallel_for(sycl::range{MX, MY}, [=](sycl::id<2> idx)
    {
      sunindextype i   = idx[0];
      sunindextype j   = idx[1];
      sunindextype tid = i * MY + j;

      realtype uij = udata[tid];
      realtype udn = (j ==      0) ? ZERO : udata[tid - 1];
      realtype uup = (j == MY - 1) ? ZERO : udata[tid + 1];
      realtype ult = (i ==      0) ? ZERO : udata[tid - MY];
      realtype urt = (i == MX - 1) ? ZERO : udata[tid + MY];

      // Set diffusion and advection terms and load into udot
      realtype hdiff = hordc * (ult - TWO * uij + urt);
      realtype vdiff = verdc * (uup - TWO * uij + udn);
      realtype hadv  = horac * (urt - ult);

      dudata[tid] = hdiff + vdiff + hadv;
    });
  });

  return 0;
}


// Jacobian-times-vector routine.
static int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu,
               void* user_data, N_Vector tmp)
{
  UserData* data = static_cast<UserData*>(user_data);

  // Extract needed constants from data
  const size_t   MX    = static_cast<size_t>(data->MX);
  const size_t   MY    = static_cast<size_t>(data->MY);
  const realtype hordc = data->hdcoef;
  const realtype horac = data->hacoef;
  const realtype verdc = data->vdcoef;

  // Extract pointers to vector data
  const realtype *vdata  = N_VGetDeviceArrayPointer(v);
  realtype       *Jvdata = N_VGetDeviceArrayPointer(Jv);

  data->myQueue->submit([&](sycl::handler& h)
  {
    h.parallel_for(sycl::range{MX, MY}, [=](sycl::id<2> idx)
    {
      sunindextype i   = idx[0];
      sunindextype j   = idx[1];
      sunindextype tid = i * MY + j;

      // set the tid-th element of Jv
      Jvdata[tid] = -TWO * (verdc + hordc) * vdata[tid];

      if (i !=      0) Jvdata[tid] += (hordc - horac) * vdata[tid - MY];
      if (i != MX - 1) Jvdata[tid] += (hordc + horac) * vdata[tid + MY];
      if (j !=      0) Jvdata[tid] += verdc * vdata[tid - 1];
      if (j != MY - 1) Jvdata[tid] += verdc * vdata[tid + 1];
    });
  });

  return 0;
}

/* ---------------------------------------------------------------------------
 * Private helper functions
 * ---------------------------------------------------------------------------*/

// Print first lines of output (problem description)
static void PrintHeader(realtype reltol, realtype abstol, realtype umax,
                        UserData* data)
{
  std::cout << "\n2-D Advection-Diffusion Equation" << std::endl;
  std::cout << "Mesh dimensions = " << data->MX << " X " <<  data->MY << std::endl;
  std::cout << "Total system size = " << data->NEQ << std::endl;
  std::cout << "Tolerance parameters: reltol = " << reltol
            << "   abstol = " << abstol << std::endl << std::endl;
  std::cout << "At t = " << T0 << "      max.norm(u) = " << umax << std::endl;
  return;
}

// Print current value
static void PrintOutput(realtype t, realtype umax, long int nst)
{
  std::cout << "At t = " << t << "   max.norm(u) = "<< umax
            << "   nst = " << nst << std::endl;
  return;
}

// Get and print some final statistics
static void PrintFinalStats(void* cvode_mem)
{
  long lenrw, leniw ;
  long lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int retval;

  retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_retval(&retval, "CVodeGetWorkSpace", 1);
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  retval = CVodeGetLinWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
  check_retval(&retval, "CVodeGetLinWorkSpace", 1);
  retval = CVodeGetNumLinIters(cvode_mem, &nli);
  check_retval(&retval, "CVodeGetNumLinIters", 1);
  retval = CVodeGetNumPrecEvals(cvode_mem, &npe);
  check_retval(&retval, "CVodeGetNumPrecEvals", 1);
  retval = CVodeGetNumPrecSolves(cvode_mem, &nps);
  check_retval(&retval, "CVodeGetNumPrecSolves", 1);
  retval = CVodeGetNumLinConvFails(cvode_mem, &ncfl);
  check_retval(&retval, "CVodeGetNumLinConvFails", 1);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  std::cout << "\nFinal Statistics.. \n\n";
  std::cout << "lenrw   = " << lenrw   << "     leniw   = " << leniw   << "\n";
  std::cout << "lenrwLS = " << lenrwLS << "     leniwLS = " << leniwLS << "\n";
  std::cout << "nst     = " << nst     << "\n";
  std::cout << "nfe     = " << nfe     << "     nfeLS   = " << nfeLS << "\n";
  std::cout << "nni     = " << nni     << "     nli     = " << nli   << "\n";
  std::cout << "nsetups = " << nsetups << "     netf    = " << netf  << "\n";
  std::cout << "npe     = " << npe     << "     nps     = " << nps   << "\n";
  std::cout << "ncfn    = " << ncfn    << "     ncfl    = " << ncfl  << "\n\n";

  return;
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_retval(void* returnvalue, const char *funcname, int opt)
{
  int *retval;

  if (opt == 0 && returnvalue == NULL)
  {
    // Check if SUNDIALS function returned NULL pointer - no memory allocated
    std::cerr << "\nSUNDIALS_ERROR: " << funcname
              << " failed - returned NULL pointer\n\n";
    return 1;
  }
  else if (opt == 1)
  {
    // Check if retval < 0
    retval = static_cast<int*>(returnvalue);
    if (*retval < 0)
    {
      std::cerr << "\nSUNDIALS_ERROR: " << funcname
                << " failed with retval = " << *retval << "\n\n";
      return 1;
    }
  }
  else if (opt == 2 && returnvalue == NULL)
  {
    // Check if function returned NULL pointer - no memory allocated
    std::cerr << "\nMEMORY_ERROR: " << funcname
              << " failed - returned NULL pointer\n\n";
    return 1;
  }

  return 0;
}
