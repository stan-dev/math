/* -----------------------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ UIUC/LLNL
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
 * The following test implements a variation of the expectation-maximization
 * problem for mixture densities from [1] with performance and convergence
 * results presented in [2]. Here, we consider a mixture density composed of
 * three univariate normal densities with a mixture density given by
 *
 *   P(x) = \sum_{i=1}^{3} alpha_i * Z_i(x|mu_i, sigma_i)
 *
 * where
 *
 *   Z_i(x|mu_i, sigma_i) =
 *        1 / (sqrt(2*pi)*sigma_i) * e^{-(x-mu_i)^2 / (2*sigma_i^2)}
 *
 * Mixture proportions {alpha_i}_{i=1}^{3} are non-negative and sum to 1. The
 * mixture proportions and variances are assumed to be known and the means
 * {mu_i}_{i=1}^3 are estimated from a set of unlabeled samples {x_k}_{k=1}^N,
 * or samples of unknown origin. Determining the unknown means distribution
 * parameters is given by the following function for 1 <= i <= 3
 *
 *   G(mu_i) = mu_i =
 *      [ \sum_{k=1}^N x_k * (alpha_i * Z_i(x_k|mu_i, sigma_i)) / (P(x_k)) ] /
 *      [ \sum_{k=1}^N (alpha_i * Z_i(x_k|mu_i, sigma_i)) / (P(x_k)) ]
 *
 * with current mean estimations being applied alongside the known mixture
 * proportions and variances as the original test case,
 *
 *  alpha_1 = 0.3   alpha_2 = 0.3   alpha_3 = 0.4
 *
 * and
 *
 *  sigma_1 = sigma_2 = sigma_3 = 1.0.
 *
 * We generate 100,000 samples for the mean distribution set
 *
 *  mu_1 = 0   mu_2 = 0.5   mu_3 = 1.0
 *
 * corresponding to a poorly separated mixture and used the same AA parameter of
 * m=3 as in [1]. We estimate a single set of mean distribution parameters
 * redundantly for every entry in a global vector u.
 *
 * 1. H.F. Walker and P. Ni, "Anderson acceleration for fixed-point iterations",
 *    SIAM Journal on Numerical Analysis, 49 (2011), pp. 1715-1735.
 *
 * 2. S. Lockhart, D.J. Gardner, C.S. Woodward, S. Thomas and L.N. Olson,
 *    "Performance of Low Synchronization Orthogonliazation Methods in Anderson
 *    Accelerated Fixed Point Solvers." arXiv preprint arXiv:2110.09667 (2021).
 *
 * Several command line options are available to change the problem parameters
 * and KINSOL settings. Use the retval --help for more information.
 * ---------------------------------------------------------------------------*/

// Header file containing UserData and function declarations
#include "kin_em_mpicuda.hpp"

// -----------------------------------------------------------------------------
// Cuda Kernels
// -----------------------------------------------------------------------------

__global__
void PxKernel(realtype *mu, realtype *Px, realtype *x,
              realtype a1, realtype a2, realtype a3, realtype scale,
              sunindextype N)
{
  // Calculate all P(x_k) for each x value
  realtype val1, val2, val3;

  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < N) {
    val1 = x[tid] - mu[0];
    val2 = x[tid] - mu[1];
    val3 = x[tid] - mu[2];

    Px[tid] = a1 * scale * exp( -(val1 * val1)/TWO );
    Px[tid] += a2 * scale * exp( -(val2 * val2)/TWO );
    Px[tid] += a3 * scale * exp( -(val3 * val3)/TWO );
  }
}

__global__
void EMKernel(realtype *mu, realtype *mu_top, realtype *mu_bottom,
              realtype *x, realtype *Px,
              realtype a1, realtype a2, realtype a3, realtype scale,
              sunindextype N)
{
  realtype val1, val2, val3;
  realtype frac1, frac2, frac3;

  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < N) {
    val1 = x[tid] - mu[0];
    val2 = x[tid] - mu[1];
    val3 = x[tid] - mu[2];

    frac1 = a1 * scale * exp( -(val1 * val1)/TWO ) / Px[tid];
    frac2 = a2 * scale * exp( -(val2 * val2)/TWO ) / Px[tid];
    frac3 = a3 * scale * exp( -(val3 * val3)/TWO ) / Px[tid];

    atomicAdd(mu_top,     x[tid] * frac1);
    atomicAdd(mu_top + 1, x[tid] * frac2);
    atomicAdd(mu_top + 2, x[tid] * frac3);

    atomicAdd(mu_bottom,     frac1);
    atomicAdd(mu_bottom + 1, frac2);
    atomicAdd(mu_bottom + 2, frac3);
  }
}

__global__
void EMKernelFin(realtype *mu, realtype *mu_top, realtype *mu_bottom,
                 sunindextype localn)
{
  int tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < localn) {
    mu[3*tid]   = mu_top[0] / mu_bottom[0];
    mu[3*tid+1] = mu_top[1] / mu_bottom[1];
    mu[3*tid+2] = mu_top[2] / mu_bottom[2];
  }
}

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // Initialize MPI
  int retval = MPI_Init(&argc, &argv);
  if (check_retval(&retval, "MPI_Init", 1)) return 1;

  // Add scope so objects are destroyed before MPI_Finalize
  {
    // User data structure
    UserData *udata = NULL;

    // Timing variables
    double t1 = 0.0;
    double t2 = 0.0;

    // MPI communicator and process ID
    MPI_Comm comm_w = MPI_COMM_WORLD;
    int myid;

    retval = MPI_Comm_rank(comm_w, &myid);
    if (check_retval(&retval, "MPI_Comm_rank", 1)) return 1;

    // Set output process flag
    bool outproc = (myid == 0);

    // SUNDIALS context
    sundials::Context sunctx(&comm_w);

    // ------------------------------------------
    // Setup UserData and parallel decomposition
    // ------------------------------------------

    // Allocate and initialize user data structure with default values. The
    // defaults may be overwritten by command line inputs in ReadInputs below.
    udata = new UserData;
    retval = InitUserData(udata);
    if (check_retval(&retval, "InitUserData", 1)) return 1;

    // Parse command line inputs
    retval = ReadInputs(&argc, &argv, udata, outproc);
    if (retval != 0) return 1;

    // Output problem setup/options
    if (outproc)
    {
      retval = PrintUserData(udata);
      if (check_retval(&retval, "PrintUserData", 1)) return 1;
    }

    // --------------------------
    // Create MPI + Cuda vectors
    // --------------------------

    // Create vector for solution
    N_Vector ulocal = N_VNew_Cuda(3 * udata->nodes_loc, sunctx);
    if (check_retval((void *) ulocal, "N_VNew_Cuda", 0)) return 1;

    N_Vector u = N_VMake_MPIPlusX(udata->comm, ulocal, sunctx);
    if (check_retval((void *) u, "N_VMake_MPIPlusX", 0)) return 1;

    // Create vector for scaling initial value
    N_Vector scale = N_VClone(u);
    if (check_retval((void *) scale, "N_VClone", 0)) return 1;
    N_VConst(ONE, scale);

    // Set initial condition
    retval = SetStartGuess(u, udata);
    if (check_retval(&retval, "RandomVec", 1)) return 1;

    // Create vector true mu values
    udata->mu_true = N_VClone(u);
    if (check_retval((void *) (udata->mu_true), "N_VClone", 0)) return 1;

    // Create temporary vector for residual and error output
    udata->vtemp = N_VClone(u);
    if (check_retval((void *) (udata->vtemp), "N_VClone", 0)) return 1;

    // Create temporary vector for mu calculation
    udata->mu_bottom = N_VNew_Cuda(3, sunctx);
    if (check_retval((void *) (udata->mu_bottom), "N_VNewCuda", 0)) return 1;

    udata->mu_top = N_VNew_Cuda(3, sunctx);
    if (check_retval((void *) (udata->mu_top), "N_VNewCuda", 0)) return 1;

    // Create vector for samples
    udata->samples_local = N_VNew_Cuda(udata->num_samples, sunctx);
    if (check_retval((void *) udata->samples_local, "N_VNew_Cuda", 0)) return 1;

    // Clone samples for temporary vector
    udata->px = N_VClone(udata->samples_local);
    if (check_retval((void *) (udata->px), "N_VClone", 0)) return 1;

    // --------------
    // Setup Mus
    // --------------

    retval = SetMus(udata);
    if (check_retval(&retval, "SetMus", 1)) return 1;

    // --------------
    // Setup Samples
    // --------------

    retval = SetupSamples(udata);
    if (check_retval(&retval, "SetupSamples", 1)) return 1;

    // --------------
    // Setup KINSOL
    // --------------

    // Initialize KINSOL memory
    void* kin_mem = KINCreate(sunctx);
    if (check_retval((void *) kin_mem, "KINCreate", 0)) return 1;

    // Set number of prior residuals used in Anderson Accleration
    retval = KINSetMAA(kin_mem, udata->maa);
    if (check_retval(&retval, "KINSetMAA", 0)) return 1;

    // Set orthogonlization routine used in Anderson Accleration
    retval = KINSetOrthAA(kin_mem, udata->orthaa);
    if (check_retval(&retval, "KINSetOrthAA", 0)) return 1;

    // Set Fixed Point Function
    retval = KINInit(kin_mem, FPFunction, u);
    if (check_retval(&retval, "KINInit", 1)) return 1;

    // Specify tolerances
    retval = KINSetFuncNormTol(kin_mem, udata->rtol);
    if (check_retval(&retval, "KINSetFuncNormTol", 1)) return 1;

    // Set maximum number of iterations
    retval = KINSetNumMaxIters(kin_mem, udata->maxits);
    if (check_retval(&retval, "KINSetMaxNumIters", 1)) return 1;

    // Set Anderson Acceleration damping parameter
    retval = KINSetDampingAA(kin_mem, udata->damping);
    if (check_retval(&retval, "KINSetDampingAA", 1)) return 1;

    // Attach user data
    retval = KINSetUserData(kin_mem, (void *) udata);
    if (check_retval(&retval, "KINSetUserData", 1)) return 1;

    // Set debugging output file
    FILE* debugfp;
    if (udata->debug)
    {
      char fname[MXSTR];
      snprintf(fname, MXSTR, "kinsol_output_%06d.txt", myid);
      debugfp = fopen(fname,"w");

      retval = KINSetDebugFile(kin_mem, debugfp);
      if (check_retval(&retval, "KINSetDebugFile", 1)) return 1;
    }

    // ----------------------------
    // Call KINSol to solve problem
    // ----------------------------

    // No scaling used
    N_VConst(ONE, scale);

    if (udata->output > 1)
    {
      retval = OpenOutput(udata);
      if (check_retval(&retval, "OpenOutput", 1)) return 1;
    }

    // Start timer
    t1 = MPI_Wtime();

    // Call main solver
    retval = KINSol(kin_mem,        // KINSol memory block
                  u,              // inital guess on input; solution vector
                  KIN_FP,         // global strategy choice
                  scale,          // scaling vector, for the variable u
                  scale);         // scaling vector for function values fval
    if (check_retval(&retval, "KINSol", 1)) return(1);

    // Stop timer
    t2 = MPI_Wtime();

    // Update timer
    udata->totaltime = t2 - t1;

    // -----------------------
    // Get solver statistics
    // -----------------------

    if (udata->output > 0 && outproc)
    {
      cout << "Final statistics:" << endl;
      retval = OutputStats(kin_mem, udata);
      if (check_retval(&retval, "OutputStats", 1)) return 1;
    }
    if (udata->output > 1)
    {
      retval = CloseOutput(udata);
      if (check_retval(&retval, "CloseOutput", 1)) return 1;
    }

    // ------------------------------
    // Print timing
    // ------------------------------

    if (udata->timing)
    {
      retval = OutputTiming(udata);
      if (check_retval(&retval, "OutputTiming", 1)) return 1;
    }

    // ------------------------------
    // Free memory
    // ------------------------------

    if (udata->debug) fclose(debugfp);
    KINFree(&kin_mem);         // Free solver memory
    N_VDestroy(u);             // Free vectors
    N_VDestroy(scale);
    FreeUserData(udata);       // Free user data
    delete udata;
  }

  // Finalize MPI
  retval = MPI_Finalize();

  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the solver
// -----------------------------------------------------------------------------

// Fixed point function to compute G(u)
static int FPFunction(N_Vector u, N_Vector f, void *user_data)
{
  int retval;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Start timer
  double t1 = MPI_Wtime();

  // Call EM Algorithm
  retval = EM(u, f, user_data);
  if (check_retval(&retval, "EM", 1)) return -1;

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->fevaltime += t2 - t1;

  // Calculate and output residual and error history
  if (udata->output > 1)
  {
    retval = WriteOutput(u, f, udata);
    if (check_retval(&retval, "WriteOutput", 1)) return 1;
  }

  // Return success
  return 0;
}

// Setup mean distribution samples
static int SetupSamples(UserData *udata)
{
  sunindextype i, j, start, end;
  realtype mean, val;

  // Access problem data
  realtype *samples_local = N_VGetHostArrayPointer_Cuda(udata->samples_local);
  if (check_retval((void *) samples_local, "N_VGetHostArrayPointer_Cuda", 0)) return 1;

  realtype *mu_host = N_VGetHostArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(udata->mu_true));
  if (check_retval((void *) mu_host, "N_VGetHostArrayPointer_Cuda", 0)) return 1;

  realtype std_dev = ONE;

  for (i = 0; i < 3; i++) {
    // Set number of samples with this mean
    if (i == 0 || i == 1) {
      end = 3 * (udata->num_samples / 10);
      start = i * end;
      end += start;
    }
    else {
      end = 4 * (udata->num_samples / 10);
      start = 2 * (3 * (udata->num_samples / 10));
      end += start;
    }

    // Setup distribution parameters
    mean = mu_host[i];
    std::default_random_engine generator;
    std::normal_distribution<realtype> distribution(mean, std_dev);

    // Get samples
    for (j = start; j < end; j++) {
      val = distribution(generator);
      samples_local[j] = val;
    }
  }

  N_VCopyToDevice_Cuda(udata->samples_local);

  // Return success
  return 0;
}

// Fill the vector u with random data
static int SetMus(UserData *udata)
{
  sunindextype i;

  realtype *mu_host = N_VGetHostArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(udata->mu_true));
  if (check_retval((void *) mu_host, "N_VGetHostArrayPointer_Cuda", 0)) return 1;

  // Fill vectors with uniform random data in [-1,1]
  for (i = 0; i < udata->nodes_loc; i++)
  {
    mu_host[3*i]   = ZERO;
    mu_host[3*i+1] = HALF;
    mu_host[3*i+2] = ONE;
  }

  N_VCopyToDevice_Cuda(N_VGetLocalVector_MPIPlusX(udata->mu_true));

  // Return success
  return 0;
}

static int SetStartGuess(N_Vector u, UserData* udata)
{
  realtype *u_host = N_VGetHostArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(u));
  if (check_retval((void *) u_host, "N_VGetHostArrayPointer_Cuda", 0)) return 1;

  for (sunindextype i = 0; i < udata->nodes_loc; i++)
  {
    u_host[3 * i]     = RCONST(0.25);
    u_host[3 * i + 1] = RCONST(3.0);
    u_host[3 * i + 2] = RCONST(0.75);
  }

  N_VCopyToDevice_Cuda(N_VGetLocalVector_MPIPlusX(u));

  // Return success
  return 0;
}


static int EM(N_Vector u, N_Vector f, void *user_data)
{
  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Set grid and block sizes for kernel launch
  unsigned block = 256;
  unsigned grid1 = (udata->num_samples + block - 1) / block;
  unsigned grid2 = (udata->nodes_loc + block - 1) / block;

  // ---------
  // PX KERNEL
  // ---------

  // Scale value for functions
  realtype scale = ONE / sqrt(TWO * PI);

  // Get input device pointers
  realtype *u_dev  = N_VGetDeviceArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(u));
  realtype *x_dev  = N_VGetDeviceArrayPointer_Cuda(udata->samples_local);

  // Get output device pointer
  realtype *Px_dev = N_VGetDeviceArrayPointer_Cuda(udata->px);

  // Compute Px
  PxKernel<<<grid1, block>>>(u_dev, Px_dev, x_dev,
                             udata->alpha1, udata->alpha2, udata->alpha3, scale,
                             udata->num_samples);

  // ---------
  // EM KERNEL
  // ---------

  // Get output device pointers
  realtype *mu_bottom_dev = N_VGetDeviceArrayPointer_Cuda(udata->mu_bottom);
  realtype *mu_top_dev    = N_VGetDeviceArrayPointer_Cuda(udata->mu_top);

  // Initilaize output vectors to zero (for sum reduction)
  N_VConst(ZERO, udata->mu_bottom);
  N_VConst(ZERO, udata->mu_top);

  EMKernel<<<grid1, block>>>(u_dev, mu_top_dev, mu_bottom_dev, x_dev, Px_dev,
                             udata->alpha1, udata->alpha2, udata->alpha3, scale,
                             udata->num_samples);

  // ------------------
  // EM FINALIZE KERNEL
  // ------------------

  realtype *f_dev = N_VGetDeviceArrayPointer_Cuda(N_VGetLocalVector_MPIPlusX(f));

  EMKernelFin<<<grid2, block>>>(f_dev, mu_top_dev, mu_bottom_dev,
                                udata->nodes_loc);

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Initialize memory allocated within Userdata
static int InitUserData(UserData *udata)
{
  int retval;

  // Sigmas
  udata->sigma = ONE;

  // Alphas - mixture proportions
  udata->alpha1 = PTTHREE;
  udata->alpha2 = PTTHREE;
  udata->alpha3 = PTFOUR;

  // MPI variables
  udata->comm = MPI_COMM_WORLD;

  // Get the number of processes
  retval = MPI_Comm_size(udata->comm, &(udata->nprocs_w));
  if (retval != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << retval << endl;
    return -1;
  }

  // Get my rank
  retval = MPI_Comm_rank(udata->comm, &(udata->myid));
  if (retval != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << retval << endl;
    return -1;
  }

  // Local number of nodes
  udata->nodes_loc = 1;

  // Global total number of nodes
  udata->nodes = udata->nodes_loc * udata->nprocs_w;

  // Integrator settings
  udata->rtol        = RCONST(1.e-8);   // relative tolerance
  udata->maa         = 3;               // 3 vectors in Anderson Acceleration space
  udata->damping     = ONE;             // no damping for Anderson Acceleration
  udata->orthaa      = 0;               // use MGS for Anderson Acceleration
  udata->maxits      = 200;             // max number of fixed point iterations

  // Vectors
  udata->samples_local = NULL;
  udata->px            = NULL;
  udata->mu_bottom     = NULL;
  udata->mu_top        = NULL;
  udata->mu_true       = NULL;

  // Number samples
  udata->num_samples = 100000;

  // Output variables
  udata->output = 1;   // 0 = no output, 1 = stats output, 2 = output to disk
  udata->vtemp  = NULL;

  // Timing variables
  udata->timing       = false;
  udata->totaltime    = 0.0;
  udata->fevaltime    = 0.0;

  udata->debug = false;

  // Return success
  return 0;
}

// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{

  // Free samples vectors
  if (udata->samples_local)
  {
    N_VDestroy(udata->samples_local);
    udata->samples_local = NULL;
  }

  // Free temporary vectors
  if (udata->px)
  {
    N_VDestroy(udata->px);
    udata->px = NULL;
  }
  if (udata->mu_bottom)
  {
    N_VDestroy(udata->mu_bottom);
    udata->mu_bottom = NULL;
  }
  if (udata->mu_top)
  {
    N_VDestroy(udata->mu_top);
    udata->mu_top = NULL;
  }
  if (udata->mu_true)
  {
    N_VDestroy(udata->mu_true);
    udata->mu_true = NULL;
  }

  // Free error vector
  if (udata->vtemp)
  {
    N_VDestroy(udata->vtemp);
    udata->vtemp = NULL;
  }

  // Free MPI communicator
  udata->comm = MPI_COMM_NULL;

  // Return success
  return 0;
}

// Read command line inputs
static int ReadInputs(int *argc, char ***argv, UserData *udata, bool outproc)
{
  // Check for input args
  int arg_idx = 1;

  while (arg_idx < (*argc))
  {
    string arg = (*argv)[arg_idx++];

    // Mesh points
    if (arg == "--nodes_loc")
    {
      udata->nodes_loc = stoi((*argv)[arg_idx++]);
    }
    // Fixed Point settings
    else if (arg == "--rtol")
    {
      udata->rtol = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--maa")
    {
      udata->maa = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--damping")
    {
      udata->damping = stod((*argv)[arg_idx++]);
    }
    else if (arg == "--orthaa")
    {
      udata->orthaa = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--maxits")
    {
      udata->maxits = stoi((*argv)[arg_idx++]);
    }
    // Output settings
    else if (arg == "--output")
    {
      udata->output = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--timing")
    {
      udata->timing = true;
    }
    else if (arg == "--debug")
    {
      udata->debug = true;
    }
    // Help
    else if (arg == "--help")
    {
      if (outproc) InputHelp();
      return -1;
    }
    // Unknown input
    else
    {
      if (outproc)
      {
        cerr << "ERROR: Invalid input " << arg << endl;
        InputHelp();
      }
      return -1;
    }
  }

  // Recompute local number of nodes
  udata->nodes = udata->nodes_loc * udata->nprocs_w;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the solution error
static int SolutionError(N_Vector u_true, N_Vector u, N_Vector err,
                         UserData *udata)
{
  // Put true solution in error vector
  SetMus(udata);

  // Compute absolute error
  N_VLinearSum(ONE, u_true, -ONE, u, err);
  N_VAbs(err, err);

  return 0;
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --nodes                 : global number of values in vector" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --maa                   : size of Anderson Acceleration subspace" << endl;
  cout << "  --damping               : damping for Anderson Acceleration" << endl;
  cout << "  --orthaa                : orthogonalization routined used in Anderson Acceleration" << endl;
  cout << "  --maxits <iterations>   : max fixed point iterations" << endl;
  cout << "  --output                : output nonlinear solver statistics" << endl;
  cout << "  --timing                : print timing data" << endl;
  cout << "  --help                  : print this message and exit" << endl;
}

// Print user data
static int PrintUserData(UserData *udata)
{
  cout << endl;
  cout << "Expectation-Maximizaton Alg. for Mixture Densities Terms:" << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  nprocs             = " << udata->nprocs_w                << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  nodes              = " << udata->nodes                   << endl;
  cout << "  nodes_loc (proc 0) = " << udata->nodes_loc               << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  sigma              = {"
             << udata->sigma << ", " << udata->sigma << ", "
             << udata->sigma << "}" << endl;
  cout << "  alpha              = {"
             << udata->alpha1 << ", " << udata->alpha2 << ", "
             << udata->alpha3 << "}" << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  rtol               = " << udata->rtol                        << endl;
  cout << "  maa                = " << udata->maa                         << endl;
  cout << "  damping            = " << udata->damping                     << endl;
  cout << "  orthaa             = " << udata->orthaa                      << endl;
  cout << "  maxits             = " << udata->maxits                      << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << "  output             = " << udata->output                      << endl;
  cout << " ------------------------------------------------------- " << endl;
  cout << endl;

  return 0;
}

// Print nonlinear solver statistics
static int OutputStats(void *kinsol_mem, UserData* udata)
{
  int retval;

  // Get solver stats
  long int nfe, nni;
  retval = KINGetNumNonlinSolvIters(kinsol_mem, &nni);
  if (check_retval(&retval, "KINGetNumNonlinSolvIters", 1)) return(1);
  retval = KINGetNumFuncEvals(kinsol_mem, &nfe);
  if (check_retval(&retval, "KINGetNumFuncEvals", 1)) return(1);

  cout << setprecision(6);

  cout << "  Func evals       = " << nfe     << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << endl;

  return 0;
}

static int OutputTiming(UserData *udata)
{
  bool outproc = (udata->myid == 0);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(6);
  }

  double maxtime = 0.0;

  MPI_Reduce(&(udata->totaltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm);
  if (outproc)
  {
    cout << "  Total time                = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->fevaltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm);
  if (outproc)
  {
    cout << "  Function evaluation time  = " << maxtime << " sec" << endl;
  }

  return 0;
}

// Open residual and error output
static int OpenOutput(UserData *udata)
{
  bool outproc = (udata->myid == 0);

  if (outproc)
  {
    stringstream fname;

    // Open output stream for residual
    fname.str("");
    fname.clear();
    fname << "EM_res_m" << udata->maa << "_orth" << udata->orthaa
          << "_len" << udata->nodes_loc << ".txt";
    udata->rout.open(fname.str());

    udata->rout << scientific;
    udata->rout << setprecision(numeric_limits<realtype>::digits10);

    // Open output stream for error
    fname.str("");
    fname.clear();
    fname << "EM_err_m" << udata->maa << "_orth" << udata->orthaa
          << "_len" << udata->nodes_loc << ".txt";
    udata->eout.open(fname.str());

    udata->eout << scientific;
    udata->eout << setprecision(numeric_limits<realtype>::digits10);
  }

  return 0;
}

// Write residual and error out to file
static int WriteOutput(N_Vector u, N_Vector f, UserData *udata)
{
  int retval;
  bool outproc = (udata->myid == 0);

  // r = \|G(u) - u\|_inf
  N_VLinearSum(ONE, f, -ONE, u, udata->vtemp);
  realtype res = N_VMaxNorm(udata->vtemp);

  // e = \|u_exact - u\|_inf
  retval = SolutionError(udata->mu_true, u, udata->vtemp, udata);
  if (check_retval(&retval, "SolutionError", 1)) return 1;
  realtype err = N_VMaxNorm(udata->vtemp);

  if (outproc)
  {
    // Output residual
    udata->rout << res;
    udata->rout << endl;

    // Output error
    udata->eout << err;
    udata->eout << endl;
  }

  return 0;
}

// Close residual and error output files
static int CloseOutput(UserData *udata)
{
  bool outproc = (udata->myid == 0);

  if (outproc)
  {
    // Close output streams
    udata->rout.close();
    udata->eout.close();
  }

  return 0;
}

// Check function return value
static int check_retval(void *flagvalue, const string funcname, int opt)
{
  // Check if the function returned a NULL pointer
  if (opt == 0)
  {
    if (flagvalue == NULL)
    {
      cerr << endl << "ERROR: " << funcname << " returned NULL pointer" << endl
           << endl;
      return 1;
    }
  }
  // Check the function return value
  else if (opt == 1 || opt == 2)
  {
    int errflag = *((int *) flagvalue);
    if  ((opt == 1 && errflag < 0) || (opt == 2 && errflag != 0))
    {
      cerr << endl << "ERROR: " << funcname << " returned = "
           << errflag << endl << endl;
      return 1;
    }
  }
  else
  {
    cerr << endl << "ERROR: check_retval called with an invalid option value"
         << endl;
    return 1;
  }

  return 0;
}

//---- end of file ----
