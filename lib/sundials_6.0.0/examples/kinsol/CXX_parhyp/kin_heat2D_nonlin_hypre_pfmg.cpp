/* -----------------------------------------------------------------------------
 * Programmer(s): Shelby Lockhart @ LLNL
 *                David Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
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
 * The following test simulates a steady-state anisotropic 2D heat equation
 * with an additional nonlinear term,
 *
 *   b = kx u_xx + ky u_yy + c(u),
 *
 * for (x,y) in [0, 1]^2, with boundary conditions
 *
 *   u = 0,
 *
 * and the forcing term
 *
 *   b(x,y) = kx 2 pi^2 (cos^2(pi x) - sin^2(pi x)) sin^2(pi y)
 *            ky 2 pi^2 (cos^2(pi y) - sin^2(pi y)) sin^2(pi x)
 *            + c(u_exact).
 *
 * Under this setup, the problem has the analytical solution
 *
 *    u(x,y) = u_exact = sin^2(pi x) sin^2(pi y).
 *
 * The spatial derivatives are computed using second-order centered differences,
 * with the data distributed over nx * ny points on a uniform spatial grid. The
 * discretized nonlinear system is given by
 *
 *   b = A u + c(u)
 *
 * where A is the discretized Laplacian operator. Solving for u results in the
 * fixed point function
 *
 *   G(u) = A^{-1} (b - c(u)).
 *
 * This problem is solved via a fixed point iteration with Anderson
 * acceleration. The linear solve in each fixed-point iteration is performed the
 * PCG linear solver using hypre's PFMG preconditioner. Several command line
 * options are available to change the problem parameters and KINSOL settings.
 * Use the flag --help for more information.
 * ---------------------------------------------------------------------------*/

// Header file containing UserData struct, function declarations, and various
// nonlinear c(u) definitions
#include "kin_heat2D_nonlin_hypre_pfmg.hpp"

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

    // Setup parallel decomposition
    retval = SetupDecomp(comm_w, udata);
    if (check_retval(&retval, "SetupDecomp", 1)) return 1;

    // Output problem setup/options
    if (outproc)
    {
      retval = PrintUserData(udata);
      if (check_retval(&retval, "PrintUserData", 1)) return 1;
    }

    // ------------------------
    // Create parallel vectors
    // ------------------------

    // Create vector for solution
    N_Vector u = N_VNew_Parallel(udata->comm_c, udata->nodes_loc, udata->nodes,
                                 sunctx);
    if (check_retval((void *) u, "N_VNew_Parallel", 0)) return 1;

    // Create vector for scaling initial value
    N_Vector scale = N_VClone(u);
    if (check_retval((void *) scale, "N_VClone", 0)) return 1;
    N_VConst(ONE, scale);

    // Set initial condition
    N_VConst(ONE, u);

    // Create vector for error
    udata->e = N_VClone(u);
    if (check_retval((void *) (udata->e), "N_VClone", 0)) return 1;

    // Create vector for b
    udata->b = N_VClone(u);
    if (check_retval((void *) (udata->b), "N_VClone", 0)) return 1;

    // Create temp vector for FPFunction evaluation
    udata->vtemp = N_VClone(u);
    if (check_retval((void *) (udata->vtemp), "N_VClone", 0)) return 1;

    // ----------------------
    // Create rhs = b - c(u)
    // ----------------------
    retval = SetupRHS(udata);
    if (check_retval(&retval, "SetupRHS", 1)) return 1;

    // ---------------------
    // Create hypre objects
    // ---------------------

    retval = SetupHypre(udata);
    if (check_retval(&retval, "SetupHypre", 1)) return 1;

    // ---------------------
    // Create linear solver
    // ---------------------

    retval = SetupLS(u, udata, sunctx);
    if (check_retval(&retval, "SetupLS", 1)) return 1;

    // --------------
    // Setup KINSOL
    // --------------

    // Initialize KINSOL memory
    void* kin_mem = KINCreate(sunctx);
    if (check_retval((void *) kin_mem, "KINCreate", 0)) return 1;

    // Set number of prior residuals used in Anderson Acceleration
    retval = KINSetMAA(kin_mem, udata->maa);
    if (check_retval(&retval, "KINSetMAA", 1)) return 1;

    // Set orthogonalization routine used in Anderson Acceleration
    retval = KINSetOrthAA(kin_mem, udata->orthaa);
    if (check_retval(&retval, "KINSetOrthAA", 1)) return 1;

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
    retval = KINSol(kin_mem, // KINSol memory block
                    u,       // initial guess on input; solution vector
                    KIN_FP,  // global strategy choice
                    scale,   // scaling vector, for the variable u
                    scale);  // scaling vector, for the function values fval
    if (check_retval(&retval, "KINSol", 1)) return 1;

    // Stop timer
    t2 = MPI_Wtime();

    // Update timer
    udata->totaltime = t2 - t1;

    // ----------------------
    // Get solver statistics
    // ----------------------
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

      retval = WriteSolution(u, udata);
      if (check_retval(&retval, "WriteSolution", 1)) return 1;
    }

    // -------------------------
    // Calculate solution error
    // -------------------------

    // Output final error
    retval = SolutionError(u, udata->e, udata);
    if (check_retval(&retval, "SolutionError", 1)) return 1;

    realtype maxerr = N_VMaxNorm(udata->e);

    if (outproc)
    {
      cout << scientific;
      cout << setprecision(numeric_limits<realtype>::digits10);
      cout << "  Max error = " << maxerr << endl;
    }

    // ------------
    // Print timing
    // ------------
    if (udata->timing)
    {
      retval = OutputTiming(udata);
      if (check_retval(&retval, "OutputTiming", 1)) return 1;
    }

    // --------------------
    // Free memory
    // --------------------

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
// Setup the parallel decomposition
// -----------------------------------------------------------------------------

static int SetupDecomp(MPI_Comm comm_w, UserData *udata)
{
  int retval;

  // Check that this has not been called before
  if (udata->Erecv != NULL || udata->Wrecv != NULL ||
      udata->Srecv != NULL || udata->Nrecv != NULL)
  {
    cerr << "SetupDecomp error: parallel decomposition already set up" << endl;
    return -1;
  }

  // Get the number of processes
  retval = MPI_Comm_size(comm_w, &(udata->nprocs_w));
  if (retval != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << retval << endl;
    return -1;
  }

  // Check the processor grid
  if ((udata->npx * udata->npy) != udata->nprocs_w)
  {
    cerr << "Error: npx * npy != nproc" << endl;
    return -1;
  }

  // Set up 2D Cartesian communicator
  int dims[2];
  dims[0] = udata->npx;
  dims[1] = udata->npy;

  int periods[2];
  periods[0] = 0;
  periods[1] = 0;

  retval = MPI_Cart_create(comm_w, 2, dims, periods, 0, &(udata->comm_c));
  if (retval != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_create = " << retval << endl;
    return -1;
  }

  // Get my rank in the new Cartesian communicator
  retval = MPI_Comm_rank(udata->comm_c, &(udata->myid_c));
  if (retval != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << retval << endl;
    return -1;
  }

  // Get dimension of the Cartesian communicator and my coordinates
  int coords[2];
  retval = MPI_Cart_get(udata->comm_c, 2, dims, periods, coords);
  if (retval != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_get = " << retval << endl;
    return -1;
  }

  // Determine local extents in x-direction
  int idx         = coords[0];
  sunindextype qx = udata->nx / dims[0];
  sunindextype rx = udata->nx % dims[0];

  udata->is = qx * idx + (idx < rx ? idx : rx);
  udata->ie = udata->is + qx - 1 + (idx < rx ? 1 : 0);

  // Sanity check
  if (udata->ie > (udata->nx - 1))
  {
    cerr << "Error ie > nx - 1" << endl;
    return -1;
  }

  // Determine local extents in y-direction
  int idy         = coords[1];
  sunindextype qy = udata->ny / dims[1];
  sunindextype ry = udata->ny % dims[1];

  udata->js = qy * idy + (idy < ry ? idy : ry);
  udata->je = udata->js + qy - 1 + (idy < ry ? 1 : 0);

  // Sanity check
  if (udata->je > (udata->ny - 1))
  {
    cerr << "Error je > ny - 1" << endl;
    return -1;
  }

  // Number of local nodes
  udata->nx_loc = (udata->ie) - (udata->is) + 1;
  udata->ny_loc = (udata->je) - (udata->js) + 1;

  // Initialize global and local vector lengths
  udata->nodes     = udata->nx * udata->ny;
  udata->nodes_loc = udata->nx_loc * udata->ny_loc;

  // Determine if this proc has neighbors
  udata->HaveNbrW = (udata->is != 0);
  udata->HaveNbrE = (udata->ie != udata->nx-1);
  udata->HaveNbrS = (udata->js != 0);
  udata->HaveNbrN = (udata->je != udata->ny-1);

  // Allocate exchange buffers if necessary
  if (udata->HaveNbrW)
  {
    udata->Wrecv = new realtype[udata->ny_loc];
    udata->Wsend = new realtype[udata->ny_loc];
  }
  if (udata->HaveNbrE)
  {
    udata->Erecv = new realtype[udata->ny_loc];
    udata->Esend = new realtype[udata->ny_loc];
  }
  if (udata->HaveNbrS)
  {
    udata->Srecv = new realtype[udata->nx_loc];
    udata->Ssend = new realtype[udata->nx_loc];
  }
  if (udata->HaveNbrN)
  {
    udata->Nrecv = new realtype[udata->nx_loc];
    udata->Nsend = new realtype[udata->nx_loc];
  }

  // MPI neighborhood information
  int nbcoords[2];

  // West neighbor
  if (udata->HaveNbrW)
  {
    nbcoords[0] = coords[0]-1;
    nbcoords[1] = coords[1];
    retval = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipW));
    if (retval != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << retval << endl;
      return -1;
    }
  }

  // East neighbor
  if (udata->HaveNbrE)
  {
    nbcoords[0] = coords[0]+1;
    nbcoords[1] = coords[1];
    retval = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipE));
    if (retval != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << retval << endl;
      return -1;
    }
  }

  // South neighbor
  if (udata->HaveNbrS)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]-1;
    retval = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipS));
    if (retval != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << retval << endl;
      return -1;
    }
  }

  // North neighbor
  if (udata->HaveNbrN)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]+1;
    retval = MPI_Cart_rank(udata->comm_c, nbcoords, &(udata->ipN));
    if (retval != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << retval << endl;
      return -1;
    }
  }

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Functions called by the solver
// -----------------------------------------------------------------------------

// Fixed point function to compute G(u) =  A^{-1} (b - c(u))
static int FPFunction(N_Vector u, N_Vector f, void *user_data)
{
  int retval;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Start timer
  double t1 = MPI_Wtime();

  // ---------------------------
  // Calculate f = b - c(u)
  // ---------------------------

  // Calculate c(u)
  retval = c(u, udata->vtemp, user_data);
  if (check_retval(&retval, "c(u)", 1)) return 1;

  // f = b -  c(u)
  N_VLinearSum(ONE, udata->b, -ONE, udata->vtemp, f);

  // ---------------------------
  // Solve Au = b - c(u)
  // ---------------------------

  // Solve system Au = f, store solution in f
  retval = SUNLinSolSolve(udata->LS, NULL, f, f, udata->epslin);
  if (check_retval(&retval, "SUNLinSolSolve", 1)) return -1;

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

// Set nonlinear function c(u)
static int SetC(UserData *udata)
{
  if (0 < udata->c_int && udata->c_int < 18)
  {
    if (udata->c_int == 1) udata->c = c1;
    else if (udata->c_int == 2) udata->c = c2;
    else if (udata->c_int == 3) udata->c = c3;
    else if (udata->c_int == 4) udata->c = c4;
    else if (udata->c_int == 5) udata->c = c5;
    else if (udata->c_int == 6) udata->c = c6;
    else if (udata->c_int == 7) udata->c = c7;
    else if (udata->c_int == 8) udata->c = c8;
    else if (udata->c_int == 9) udata->c = c9;
    else if (udata->c_int == 10) udata->c = c10;
    else if (udata->c_int == 11) udata->c = c11;
    else if (udata->c_int == 12) udata->c = c12;
    else if (udata->c_int == 13) udata->c = c13;
    else if (udata->c_int == 14) udata->c = c14;
    else if (udata->c_int == 15) udata->c = c15;
    else if (udata->c_int == 16) udata->c = c16;
    else if (udata->c_int == 17) udata->c = c17;
  }
  else return 1;

  // Return success
  return 0;
}

// Create RHS = b - c(u)
static int SetupRHS(void *user_data)
{
  int          retval;
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // -------------------------------------------------------
  // Setup function c(u) for FPFunction based on user input
  // -------------------------------------------------------

  retval = SetC(udata);
  if (check_retval(&retval, "SetC", 1)) return 1;

  // ----------------------
  // Setup b for FPFunction
  // ----------------------

  // Access data array
  realtype *barray = N_VGetArrayPointer(udata->b);
  if (check_retval((void *) barray, "N_VGetArrayPointer", 0)) return 1;

  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, udata->b);

  // Iterate over subdomain and compute forcing term (b)
  realtype x, y;
  realtype sin_sqr_x, sin_sqr_y;
  realtype cos_sqr_x, cos_sqr_y;

  realtype bx = (udata->kx) * TWO * PI * PI;
  realtype by = (udata->ky) * TWO * PI * PI;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      x = (udata->is + i) * udata->dx;
      y = (udata->js + j) * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      cos_sqr_x = cos(PI * x) * cos(PI * x);
      cos_sqr_y = cos(PI * y) * cos(PI * y);

      barray[IDX(i,j,nx_loc)] =
        bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y
        + by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x;
    }
  }

  // Calculate c(u_exact) and add to forcing term (b)
  retval = Solution(udata->e, udata);
  if (check_retval(&retval, "rhs", 1)) return 1;

  retval = c(udata->e, udata->vtemp, user_data);
  if (check_retval(&retval, "c(u)", 1)) return 1;

  // b = Au_exact + c(u_xact)
  N_VLinearSum(ONE, udata->vtemp, ONE, udata->b, udata->b);

  // Return success
  return 0;
}

// c(u)
static int c(N_Vector u, N_Vector z, void *user_data)
{
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_retval((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  realtype *zarray = N_VGetArrayPointer(z);
  if (check_retval((void *) zarray, "N_VGetArrayPointer", 0)) return 1;

  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, z);

  // Iterate over subdomain and compute z = c(u)
  realtype u_val;

  for (j = jstart; j < jend; j++)
  {
    for (i = istart; i < iend; i++)
    {
      u_val = uarray[IDX(i,j,nx_loc)];

      zarray[IDX(i,j,nx_loc)] = udata->c(u_val);
    }
  }

  // Return success
  return 0;
}

// Create PCG Linear solver
static int SetupLS(N_Vector u, void *user_data, SUNContext sunctx)
{
  int retval;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  int prectype = PREC_RIGHT;

  // Create linear solver
  udata->LS = SUNLinSol_PCG(u, prectype, udata->liniters, sunctx);
  if (check_retval((void *) udata->LS, "SUNLinSol_PCG", 0)) return 1;

  // Set ATimes
  retval = SUNLinSolSetATimes(udata->LS, user_data, JTimes);
  if (check_retval(&retval, "SUNLinSolSetATimes", 1)) return 1;

  // Attach preconditioner
  retval = SUNLinSolSetPreconditioner(udata->LS, user_data, PSetup, PSolve);
  if (check_retval(&retval, "SUNLinSolSetPreconditioner", 1)) return 1;

  // Initialize solver
  retval = SUNLinSolInitialize(udata->LS);
  if (check_retval(&retval, "SUNLinSolInitialize", 1)) return 1;

  // Setup solver and preconditioner
  retval = SUNLinSolSetup(udata->LS, NULL);
  if (check_retval(&retval, "SUNLinSolSetup", 1)) return 1;

  // Return success
  return 0;
}

// Jacobian-vector product function
static int JTimes(void *user_data, N_Vector v, N_Vector Jv)
{
  int retval;

  // Start timer
  double t1 = MPI_Wtime();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Insert input N_Vector entries into HYPRE vector and assemble
  retval = HYPRE_StructVectorSetBoxValues(udata->vvec,
                                          udata->ilower, udata->iupper,
                                          N_VGetArrayPointer(v));
  if (retval != 0) return -1;

  retval = HYPRE_StructVectorAssemble(udata->vvec);
  if (retval != 0) return -1;

  // Initialize output HYPRE vector and assemble
  retval = HYPRE_StructVectorSetConstantValues(udata->Jvvec, ZERO);
  if (retval != 0) return -1;

  retval = HYPRE_StructVectorAssemble(udata->Jvvec);
  if (retval != 0) return -1;

  // Compute the matrix-vector product
  retval = HYPRE_StructMatrixMatvec(ONE,
                                    udata->Jmatrix,
                                    udata->vvec,
                                    ZERO,
                                    udata->Jvvec);
  if (retval != 0) return -1;

  // Extract matrix-vector product values
  retval = HYPRE_StructVectorGetBoxValues(udata->Jvvec,
                                          udata->ilower, udata->iupper,
                                          N_VGetArrayPointer(Jv));
  if (retval != 0) return -1;

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->jvtime += t2 - t1;

  // Return success
  return 0;
}

// Preconditioner setup routine
static int PSetup(void *user_data)
{
  int retval;

  // Start timer
  double t1 = MPI_Wtime();

  // Access problem data
  UserData *udata = (UserData *) user_data;

  // Assemble matrix
  retval = HYPRE_StructMatrixAssemble(udata->Jmatrix);
  if (retval != 0) return -1;

  // -----------
  // Setup PFMG
  // -----------

  // Set rhs/solution vectors as all zero for now
  retval = HYPRE_StructVectorSetConstantValues(udata->bvec, ZERO);
  if (retval != 0) return -1;

  retval = HYPRE_StructVectorAssemble(udata->bvec);
  if (retval != 0) return -1;

  retval = HYPRE_StructVectorSetConstantValues(udata->xvec, ZERO);
  if (retval != 0) return -1;

  retval = HYPRE_StructVectorAssemble(udata->xvec);
  if (retval != 0) return -1;

  // Free the existing preconditioner if necessary
  if (udata->precond) HYPRE_StructPFMGDestroy(udata->precond);

  // Create the new preconditioner
  retval = HYPRE_StructPFMGCreate(udata->comm_c, &(udata->precond));
  if (retval != 0) return -1;

  // Signal that the inital guess is zero
  retval = HYPRE_StructPFMGSetZeroGuess(udata->precond);
  if (retval != 0) return -1;

  // tol <= 0.0 means do the max number of iterations
  retval = HYPRE_StructPFMGSetTol(udata->precond, ZERO);
  if (retval != 0) return -1;

  // Use one v-cycle
  retval = HYPRE_StructPFMGSetMaxIter(udata->precond, 1);
  if (retval != 0) return -1;

  // Use non-Galerkin corase grid operator
  retval = HYPRE_StructPFMGSetRAPType(udata->precond, 1);
  if (retval != 0) return -1;

  // Set the relaxation type
  retval = HYPRE_StructPFMGSetRelaxType(udata->precond, udata->pfmg_relax);
  if (retval != 0) return -1;

  // Set the number of pre and post relaxation sweeps
  retval = HYPRE_StructPFMGSetNumPreRelax(udata->precond, udata->pfmg_nrelax);
  if (retval != 0) return -1;

  retval = HYPRE_StructPFMGSetNumPostRelax(udata->precond, udata->pfmg_nrelax);
  if (retval != 0) return -1;

  // Set up the solver
  retval = HYPRE_StructPFMGSetup(udata->precond, udata->Jmatrix,
                                 udata->bvec, udata->xvec);
  if (retval != 0) return -1;

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->psetuptime += t2 - t1;

  // Return success
  return 0;
}

// Preconditioner solve routine for Pz = r
static int PSolve(void *user_data, N_Vector r, N_Vector z,
                  realtype tol, int lr)
{
  int retval;

  // Start timer
  double t1 = MPI_Wtime();

  // Access user_data structure
  UserData *udata = (UserData *) user_data;

  // Insert rhs N_Vector entries into HYPRE vector b and assemble
  retval = HYPRE_StructVectorSetBoxValues(udata->bvec,
                                          udata->ilower, udata->iupper,
                                          N_VGetArrayPointer(r));
  if (retval != 0) return -1;

  retval = HYPRE_StructVectorAssemble(udata->bvec);
  if (retval != 0) return -1;

  // Set the initial guess into HYPRE vector x and assemble
  retval = HYPRE_StructVectorSetConstantValues(udata->xvec, ZERO);
  if (retval != 0) return -1;

  retval = HYPRE_StructVectorAssemble(udata->xvec);
  if (retval != 0) return -1;

  // Solve the linear system
  retval = HYPRE_StructPFMGSolve(udata->precond, udata->Jmatrix,
                                 udata->bvec, udata->xvec);

  // If a convergence error occured, clear the error and continue. For any
  // other error return with a recoverable error.
  if (retval == HYPRE_ERROR_CONV) HYPRE_ClearError(HYPRE_ERROR_CONV);
  else if (retval != 0) return 1;

  // Update precond statistics
  HYPRE_Int itmp;
  retval = HYPRE_StructPFMGGetNumIterations(udata->precond, &itmp);
  if (retval != 0) return -1;

  udata->pfmg_its += itmp;

  // Extract solution values
  retval = HYPRE_StructVectorGetBoxValues(udata->xvec,
                                          udata->ilower, udata->iupper,
                                          N_VGetArrayPointer(z));
  if (retval != 0) return -1;

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->psolvetime += t2 - t1;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Preconditioner helper functions
// -----------------------------------------------------------------------------

// Create hypre objects
static int SetupHypre(UserData *udata)
{
  int retval, result;

  // Check input
  if (udata == NULL) return -1;

  // Check if the grid or stencil have been created
  if ((udata->grid != NULL || udata->stencil != NULL))
  {
    cerr << "SetupHypre error: grid or stencil already exists" << endl;
    return -1;
  }

  // Check for valid 2D Cartesian MPI communicator
  retval = MPI_Topo_test(udata->comm_c, &result);
  if ((retval != MPI_SUCCESS) || (result != MPI_CART))
  {
    cerr << "SetupHypre error: communicator is not Cartesian" << endl;
    return -1;
  }

  retval = MPI_Cartdim_get(udata->comm_c, &result);
  if ((retval != MPI_SUCCESS) || (result != 2))
  {
    cerr << "SetupHypre error: communicator is not 2D" << endl;
    return -1;
  }

  // -----
  // Grid
  // -----

  // Create 2D grid object
  retval = HYPRE_StructGridCreate(udata->comm_c, 2, &(udata->grid));
  if (retval != 0) { FreeUserData(udata); return -1; }

  // Set grid extents (lower left and upper right corners)
  udata->ilower[0] = udata->is;
  udata->ilower[1] = udata->js;

  udata->iupper[0] = udata->ie;
  udata->iupper[1] = udata->je;

  retval = HYPRE_StructGridSetExtents(udata->grid, udata->ilower, udata->iupper);
  if (retval != 0) { FreeUserData(udata); return -1; }

  // Assemble the grid
  retval = HYPRE_StructGridAssemble(udata->grid);
  if (retval != 0) { FreeUserData(udata); return -1; }

  // --------
  // Stencil
  // --------

  // Create the 2D 5 point stencil object
  retval = HYPRE_StructStencilCreate(2, 5, &(udata->stencil));
  if (retval != 0) { FreeUserData(udata); return -1; }

  // Set the stencil entries (center, left, right, bottom, top)
  HYPRE_Int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

  for (int entry = 0; entry < 5; entry++)
  {
    retval = HYPRE_StructStencilSetElement(udata->stencil, entry, offsets[entry]);
    if (retval != 0) { FreeUserData(udata); return -1; }
  }

  // -----------
  // Work array
  // -----------

  udata->nwork = 5 * udata->nodes_loc;
  udata->work  = NULL;
  udata->work  = new HYPRE_Real[udata->nwork];
  if (udata->work == NULL) { FreeUserData(udata); return -1; }

  // ---------
  // x vector
  // ---------

  retval = HYPRE_StructVectorCreate(udata->comm_c, udata->grid, &(udata->xvec));
  if (retval != 0) { FreeUserData(udata); return -1; }

  retval = HYPRE_StructVectorInitialize(udata->xvec);
  if (retval != 0) { FreeUserData(udata); return -1; }

  // ---------
  // b vector
  // ---------

  retval = HYPRE_StructVectorCreate(udata->comm_c, udata->grid, &(udata->bvec));
  if (retval != 0) { FreeUserData(udata); return -1; }

  retval = HYPRE_StructVectorInitialize(udata->bvec);
  if (retval != 0) { FreeUserData(udata); return -1; }

  // ---------
  // v vector
  // ---------

  retval = HYPRE_StructVectorCreate(udata->comm_c, udata->grid, &(udata->vvec));
  if (retval != 0) { FreeUserData(udata); return -1; }

  retval = HYPRE_StructVectorInitialize(udata->vvec);
  if (retval != 0) { FreeUserData(udata); return -1; }

  // ----------
  // Jv vector
  // ----------

  retval = HYPRE_StructVectorCreate(udata->comm_c, udata->grid, &(udata->Jvvec));
  if (retval != 0) { FreeUserData(udata); return -1; }

  retval = HYPRE_StructVectorInitialize(udata->Jvvec);
  if (retval != 0) { FreeUserData(udata); return -1; }

  // ---------
  // J matrix
  // ---------

  retval = HYPRE_StructMatrixCreate(udata->comm_c, udata->grid, udata->stencil,
                                    &(udata->Jmatrix));
  if (retval != 0) { FreeUserData(udata); return -1; }

  retval = HYPRE_StructMatrixInitialize(udata->Jmatrix);
  if (retval != 0) { FreeUserData(udata); return -1; }

  // --------------------
  // PFMG preconditioner
  // --------------------

  // Note a new PFMG preconditioner must be created and attached each time the
  // linear system is updated. As such it is constructed in the preconditioner
  // setup function (if enabled).
  udata->precond = NULL;

  // --------------
  // Fill Jacobian
  // --------------

  retval = Jac(udata);
  if (retval != 0) { FreeUserData(udata); return -1; }

  retval = HYPRE_StructMatrixAssemble(udata->Jmatrix);
  if (retval != 0) { FreeUserData(udata); return -1; }

  return 0;
}

// Jac function to compute the ODE RHS function Jacobian, (df/dy)(t,y).
static int Jac(UserData *udata)
{
  // Shortcuts to hypre matrix and grid extents, work array, etc.
  HYPRE_StructMatrix Jmatrix = udata->Jmatrix;

  HYPRE_Int ilower[2];
  HYPRE_Int iupper[2];

  ilower[0] = udata->ilower[0];
  ilower[1] = udata->ilower[1];

  iupper[0] = udata->iupper[0];
  iupper[1] = udata->iupper[1];

  HYPRE_Int   nwork = udata->nwork;
  HYPRE_Real *work  = udata->work;

  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Matrix stencil: center, left, right, bottom, top
  HYPRE_Int entries[5] = {0, 1, 2, 3, 4};
  HYPRE_Int entry[1];

  // Grid extents for setting boundary entries
  HYPRE_Int bc_ilower[2];
  HYPRE_Int bc_iupper[2];

  // Loop counters
  HYPRE_Int idx, ix, iy;

  // hypre return retval
  int retval;

  // ----------
  // Compute J
  // ----------

  // Start timer
  double t1 = MPI_Wtime();

  // Only do work if the box is non-zero in size
  if ((ilower[0] <= iupper[0]) &&
      (ilower[1] <= iupper[1]))
  {
    // Jacobian values
    realtype cx = udata->kx / (udata->dx * udata->dx);
    realtype cy = udata->ky / (udata->dy * udata->dy);
    realtype cc = -TWO * (cx + cy);

    // --------------------------------
    // Set matrix values for all nodes
    // --------------------------------

    // Set the matrix interior entries (center, left, right, bottom, top)
    idx = 0;
    for (iy = 0; iy < ny_loc; iy++)
    {
      for (ix = 0; ix < nx_loc; ix++)
      {
        work[idx]     = cc;
        work[idx + 1] = cx;
        work[idx + 2] = cx;
        work[idx + 3] = cy;
        work[idx + 4] = cy;
        idx += 5;
      }
    }

    // Modify the matrix
    retval = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                            ilower, iupper,
                                            5, entries, work);
    if (retval != 0) return -1;

    // ----------------------------------------
    // Correct matrix values at boundary nodes
    // ----------------------------------------

    // Set the matrix boundary entries (center, left, right, bottom, top)
    if (ilower[1] == 0 ||
        iupper[1] == (udata->ny - 1) ||
        ilower[0] == 0 ||
        iupper[0] == (udata->nx - 1))
    {
      idx = 0;
      for (iy = 0; iy < ny_loc; iy++)
      {
        for (ix = 0; ix < nx_loc; ix++)
        {
          work[idx]     = ONE;
          work[idx + 1] = ZERO;
          work[idx + 2] = ZERO;
          work[idx + 3] = ZERO;
          work[idx + 4] = ZERO;
          idx += 5;
        }
      }
    }

    // Set cells on western boundary
    if (ilower[0] == 0)
    {
      // Grid cell on south-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];

      // Grid cell on north-west corner
      bc_iupper[0] = ilower[0];
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        retval = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                                bc_ilower, bc_iupper,
                                                5, entries, work);
        if (retval != 0) return -1;
      }
    }

    // Set cells on eastern boundary
    if (iupper[0] == (udata->nx - 1))
    {
      // Grid cell on south-east corner
      bc_ilower[0] = iupper[0];
      bc_ilower[1] = ilower[1];

      // Grid cell on north-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        retval = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                                bc_ilower, bc_iupper,
                                                5, entries, work);
        if (retval != 0) return -1;
      }
    }

    // Correct cells on southern boundary
    if (ilower[1] == 0)
    {
      // Grid cell on south-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = ilower[1];

      // Grid cell on south-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = ilower[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        retval = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                                bc_ilower, bc_iupper,
                                                5, entries, work);
        if (retval != 0) return -1;
      }
    }

    // Set cells on northern boundary
    if (iupper[1] == (udata->ny - 1))
    {
      // Grid cell on north-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = iupper[1];

      // Grid cell on north-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        retval = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                                bc_ilower, bc_iupper,
                                                5, entries, work);
        if (retval != 0) return -1;
      }
    }

    // -----------------------------------------------------------
    // Remove connections between the interior and boundary nodes
    // -----------------------------------------------------------

    // Zero out work array
    for (ix = 0; ix < nwork; ix++)
    {
      work[ix] = ZERO;
    }

    // Second column of nodes (depends on western boundary)
    if ((ilower[0] <= 1) && (iupper[0] >= 1))
    {
      // Remove western dependency
      entry[0] = 1;

      // Grid cell on south-west corner
      bc_ilower[0] = 1;
      bc_ilower[1] = ilower[1];

      // Grid cell on north-west corner
      bc_iupper[0] = 1;
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        retval = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                                bc_ilower, bc_iupper,
                                                1, entry, work);
        if (retval != 0) return -1;
      }
    }

    // Next to last column (depends on eastern boundary)
    if ((ilower[0] <= (udata->nx - 2)) &&
        (iupper[0] >= (udata->nx - 2)))
    {
      // Remove eastern dependency
      entry[0] = 2;

      // Grid cell on south-east corner
      bc_ilower[0] = udata->nx - 2;
      bc_ilower[1] = ilower[1];

      // Grid cell on north-east corner
      bc_iupper[0] = udata->nx - 2;
      bc_iupper[1] = iupper[1];

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        retval = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                                bc_ilower, bc_iupper,
                                                1, entry, work);
        if (retval != 0) return -1;
      }
    }

    // Second row of nodes (depends on southern boundary)
    if ((ilower[1] <= 1) && (iupper[1] >= 1))
    {
      // Remove southern dependency
      entry[0] = 3;

      // Grid cell on south-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = 1;

      // Grid cell on south-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = 1;

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        retval = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                                bc_ilower, bc_iupper,
                                                1, entry, work);
        if (retval != 0) return -1;
      }
    }

    // Next to last row of nodes (depends on northern boundary)
    if ((ilower[1] <= (udata->ny - 2)) &&
        (iupper[1] >= (udata->ny - 2)))
    {
      // Remove northern dependency
      entry[0] = 4;

      // Grid cell on north-west corner
      bc_ilower[0] = ilower[0];
      bc_ilower[1] = udata->ny - 2;

      // Grid cell on north-east corner
      bc_iupper[0] = iupper[0];
      bc_iupper[1] = udata->ny - 2;

      // Only do work if the box is non-zero in size
      if ((bc_ilower[0] <= bc_iupper[0]) && (bc_ilower[1] <= bc_iupper[1]))
      {
        // Modify the matrix
        retval = HYPRE_StructMatrixSetBoxValues(Jmatrix,
                                                bc_ilower, bc_iupper,
                                                1, entry, work);
        if (retval != 0) return -1;
      }
    }
  }

  // The matrix is assembled in hypre setup

  // Stop timer
  double t2 = MPI_Wtime();

  // Update timer
  udata->matfilltime += t2 - t1;

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// UserData and input functions
// -----------------------------------------------------------------------------

// Initialize memory allocated within Userdata
static int InitUserData(UserData *udata)
{
  // Diffusion coefficient
  udata->kx = ONE;
  udata->ky = ONE;

  // Upper bounds in x and y directions
  udata->xu = ONE;
  udata->yu = ONE;

  // Global number of nodes in the x and y directions
  udata->nx    = 64;
  udata->ny    = 64;
  udata->nodes = udata->nx * udata->ny;

  // Mesh spacing in the x and y directions
  udata->dx = udata->xu / (udata->nx - 1);
  udata->dy = udata->yu / (udata->ny - 1);

  // Locals number of nodes in the x and y directions (set in SetupDecomp)
  udata->nx_loc    = 0;
  udata->ny_loc    = 0;
  udata->nodes_loc = 0;

  // Global indices of this subdomain (set in SetupDecomp)
  udata->is = 0;
  udata->ie = 0;
  udata->js = 0;
  udata->je = 0;

  // MPI variables (set in SetupDecomp)
  udata->comm_c = MPI_COMM_NULL;

  udata->nprocs_w = 1;
  udata->npx      = 1;
  udata->npy      = 1;

  udata->myid_c = 0;

  // Flags denoting neighbors (set in SetupDecomp)
  udata->HaveNbrW = true;
  udata->HaveNbrE = true;
  udata->HaveNbrS = true;
  udata->HaveNbrN = true;

  // Exchange receive buffers (allocated in SetupDecomp)
  udata->Erecv = NULL;
  udata->Wrecv = NULL;
  udata->Nrecv = NULL;
  udata->Srecv = NULL;

  // Exchange send buffers (allocated in SetupDecomp)
  udata->Esend = NULL;
  udata->Wsend = NULL;
  udata->Nsend = NULL;
  udata->Ssend = NULL;

  // Neighbors IDs (set in SetupDecomp)
  udata->ipW = -1;
  udata->ipE = -1;
  udata->ipS = -1;
  udata->ipN = -1;

  // Fixed Point Solver settings
  udata->rtol        = RCONST(1.e-8);   // relative tolerance
  udata->maa         = 0;               // no Anderson Acceleration
  udata->damping     = ONE;             // no damping for Anderson Acceleration
  udata->orthaa      = 0;               // use MGS for Anderson Acceleration
  udata->maxits      = 100;             // max number of fixed point iterations


  // Linear solver and preconditioner options
  udata->lsinfo    = false;         // output residual history
  udata->liniters  = 20;            // max linear iterations
  udata->epslin    = RCONST(1.e-8); // use default (0.05)

  // Linear solver object
  udata->LS    = NULL;

  // c function
  udata->c     = NULL;
  udata->c_int = 1;

  // Vectors
  udata->b     = NULL;
  udata->vtemp = NULL;

  // hypre objects
  udata->grid    = NULL;
  udata->stencil = NULL;
  udata->Jmatrix = NULL;
  udata->bvec    = NULL;
  udata->xvec    = NULL;
  udata->vvec    = NULL;
  udata->Jvvec   = NULL;
  udata->precond = NULL;

  // hypre grid extents
  udata->ilower[0] = 0;
  udata->ilower[1] = 0;

  udata->iupper[0] = 0;
  udata->iupper[1] = 0;

  // hypre workspace
  udata->nwork = 0;
  udata->work  = NULL;

  // hypre counters
  udata->pfmg_its = 0;

  // hypre PFMG settings
  udata->pfmg_relax  = 2;
  udata->pfmg_nrelax = 2;

  // Output variables
  udata->output = 1;   // 0 = no output, 1 = stats output, 2 = output to disk
  udata->e      = NULL;

  // Timing variables
  udata->timing       = false;
  udata->totaltime    = 0.0;
  udata->fevaltime    = 0.0;
  udata->matfilltime  = 0.0;
  udata->jvtime       = 0.0;
  udata->psetuptime   = 0.0;
  udata->psolvetime   = 0.0;

  // Return success
  return 0;
}

// Free memory allocated within Userdata
static int FreeUserData(UserData *udata)
{
  // Free exchange buffers
  if (udata->Wrecv != NULL)  delete[] udata->Wrecv;
  if (udata->Wsend != NULL)  delete[] udata->Wsend;
  if (udata->Erecv != NULL)  delete[] udata->Erecv;
  if (udata->Esend != NULL)  delete[] udata->Esend;
  if (udata->Srecv != NULL)  delete[] udata->Srecv;
  if (udata->Ssend != NULL)  delete[] udata->Ssend;
  if (udata->Nrecv != NULL)  delete[] udata->Nrecv;
  if (udata->Nsend != NULL)  delete[] udata->Nsend;

  // Free Linear solver
  if (udata->LS != NULL) SUNLinSolFree_PCG(udata->LS);

  // Free b vector
  if (udata->b)
  {
    N_VDestroy(udata->b);
    udata->b = NULL;
  }

  // Free temporary vector
  if (udata->vtemp)
  {
    N_VDestroy(udata->vtemp);
    udata->vtemp = NULL;
  }

  // Free hypre preconditioner data
  if (udata->grid    != NULL) HYPRE_StructGridDestroy(udata->grid);
  if (udata->stencil != NULL) HYPRE_StructStencilDestroy(udata->stencil);
  if (udata->Jmatrix != NULL) HYPRE_StructMatrixDestroy(udata->Jmatrix);
  if (udata->bvec    != NULL) HYPRE_StructVectorDestroy(udata->bvec);
  if (udata->xvec    != NULL) HYPRE_StructVectorDestroy(udata->xvec);
  if (udata->vvec    != NULL) HYPRE_StructVectorDestroy(udata->vvec);
  if (udata->Jvvec   != NULL) HYPRE_StructVectorDestroy(udata->Jvvec);
  if (udata->precond != NULL) HYPRE_StructPFMGDestroy(udata->precond);
  if (udata->work    != NULL) delete[] udata->work;

  // Free MPI Cartesian communicator
  if (udata->comm_c != MPI_COMM_NULL)
    MPI_Comm_free(&(udata->comm_c));

  // Free error vector
  if (udata->e)
  {
    N_VDestroy(udata->e);
    udata->e = NULL;
  }

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
    if (arg == "--mesh")
    {
      udata->nx = stoi((*argv)[arg_idx++]);
      udata->ny = stoi((*argv)[arg_idx++]);
    }
    // MPI processes
    else if (arg == "--np")
    {
      udata->npx = stoi((*argv)[arg_idx++]);
      udata->npy = stoi((*argv)[arg_idx++]);
    }
    // Domain upper bounds
    else if (arg == "--domain")
    {
      udata->xu = stoi((*argv)[arg_idx++]);
      udata->yu = stoi((*argv)[arg_idx++]);
    }
    // Diffusion parameters
    else if (arg == "--k")
    {
      udata->kx = stod((*argv)[arg_idx++]);
      udata->ky = stod((*argv)[arg_idx++]);
    }
    // Solver settings
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
    // RHS settings
    else if (arg == "--c")
    {
      udata->c_int = stoi((*argv)[arg_idx++]);
    }
    // Linear solver settings
    else if (arg == "--lsinfo")
    {
      udata->lsinfo = true;
    }
    else if (arg == "--liniters")
    {
      udata->liniters = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--epslin")
    {
      udata->epslin = stod((*argv)[arg_idx++]);
    }
    // PFMG settings
    else if (arg == "--pfmg_relax")
    {
      udata->pfmg_relax = stoi((*argv)[arg_idx++]);
    }
    else if (arg == "--pfmg_nrelax")
    {
      udata->pfmg_nrelax = stoi((*argv)[arg_idx++]);
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

  // Recompute total number of nodes
  udata->nodes = udata->nx * udata->ny;

  // Recompute x and y mesh spacing
  udata->dx = (udata->xu) / (udata->nx - 1);
  udata->dy = (udata->yu) / (udata->ny - 1);

  // Return success
  return 0;
}

// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the exact solution
static int Solution(N_Vector u, UserData *udata)
{
  realtype x, y;
  realtype sin_sqr_x, sin_sqr_y;

  // Initialize u to zero (handles boundary conditions)
  N_VConst(ZERO, u);

  // Iterative over domain interior
  sunindextype istart = (udata->HaveNbrW) ? 0 : 1;
  sunindextype iend   = (udata->HaveNbrE) ? udata->nx_loc : udata->nx_loc - 1;

  sunindextype jstart = (udata->HaveNbrS) ? 0 : 1;
  sunindextype jend   = (udata->HaveNbrN) ? udata->ny_loc : udata->ny_loc - 1;

  realtype *uarray = N_VGetArrayPointer(u);
  if (check_retval((void *) uarray, "N_VGetArrayPointer", 0)) return 1;

  for (sunindextype j = jstart; j < jend; j++)
  {
    for (sunindextype i = istart; i < iend; i++)
    {
      x  = (udata->is + i) * udata->dx;
      y  = (udata->js + j) * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      uarray[IDX(i,j,udata->nx_loc)] = sin_sqr_x * sin_sqr_y;
    }
  }

  return 0;
}

// Compute the solution error
static int SolutionError(N_Vector u, N_Vector e, UserData *udata)
{
  // Compute true solution
  int retval = Solution(e, udata);
  if (retval != 0) return -1;

  // Compute absolute error
  N_VLinearSum(ONE, u, -ONE, e, e);
  N_VAbs(e, e);

  return 0;
}

// Print command line options
static void InputHelp()
{
  cout << endl;
  cout << "Command line options:" << endl;
  cout << "  --mesh <nx> <ny>        : mesh points in the x and y directions" << endl;
  cout << "  --np <npx> <npy>        : number of MPI processes in the x and y directions" << endl;
  cout << "  --domain <xu> <yu>      : domain upper bound in the x and y direction" << endl;
  cout << "  --k <kx> <ky>           : diffusion coefficients" << endl;
  cout << "  --rtol <rtol>           : relative tolerance" << endl;
  cout << "  --maa <maa>             : number of previous residuals for Anderson Acceleration" << endl;
  cout << "  --damping <damping>     : damping for Anderson Acceleration " << endl;
  cout << "  --orthaa <orthaa>       : orthogonalization routine used in Anderson Acceleration " << endl;
  cout << "  --c <c_int>             : nonlinear function parameter" << endl;
  cout << "  --lsinfo                : output residual history" << endl;
  cout << "  --liniters <iters>      : max number of iterations" << endl;
  cout << "  --epslin <factor>       : linear tolerance factor" << endl;
  cout << "  --pfmg_relax <types>    : relaxtion type in PFMG" << endl;
  cout << "  --pfmg_nrelax <iters>   : pre/post relaxtion sweeps in PFMG" << endl;
  cout << "  --output                : output nonlinear solver statistics" << endl;
  cout << "  --maxits <maxits>       : max fixed point iterations" << endl;
  cout << "  --timing                : print timing data" << endl;
  cout << "  --help                  : print this message and exit" << endl;
}

// Print user data
static int PrintUserData(UserData *udata)
{
  cout << endl;
  cout << "2D Stationary Heat PDE + Nonlinear term test problem:"  << endl;
  cout << " --------------------------------- "                    << endl;
  cout << "  nprocs         = " << udata->nprocs_w                 << endl;
  cout << "  npx            = " << udata->npx                      << endl;
  cout << "  npy            = " << udata->npy                      << endl;
  cout << " --------------------------------- "                    << endl;
  cout << "  kx             = " << udata->kx                       << endl;
  cout << "  ky             = " << udata->ky                       << endl;
  cout << "  xu             = " << udata->xu                       << endl;
  cout << "  yu             = " << udata->yu                       << endl;
  cout << "  nx             = " << udata->nx                       << endl;
  cout << "  ny             = " << udata->ny                       << endl;
  cout << "  nxl (proc 0)   = " << udata->nx_loc                   << endl;
  cout << "  nyl (proc 0)   = " << udata->ny_loc                   << endl;
  cout << "  dx             = " << udata->dx                       << endl;
  cout << "  dy             = " << udata->dy                       << endl;
  cout << " --------------------------------- "                    << endl;
  cout << "  rtol           = " << udata->rtol                     << endl;
  cout << "  maa            = " << udata->maa                      << endl;
  cout << "  damping        = " << udata->damping                  << endl;
  cout << "  orthaa         = " << udata->orthaa                   << endl;
  cout << "  maxits         = " << udata->maxits                   << endl;
  cout << " --------------------------------- "                    << endl;
  cout << "  c              = " << udata->c_int                    << endl;
  cout << " --------------------------------- "                    << endl;
  cout << "  linear solver  = PCG" << endl;
  cout << "  lin iters      = " << udata->liniters                 << endl;
  cout << "  eps lin        = " << udata->epslin                   << endl;
  cout << "  pfmg_relax     = " << udata->pfmg_relax               << endl;
  cout << "  pfmg_nrelax    = " << udata->pfmg_nrelax              << endl;
  cout << " --------------------------------- "                    << endl;
  cout << "  output         = " << udata->output                   << endl;
  cout << " --------------------------------- "                    << endl;
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
  if (check_retval(&retval, "KinGetNumNonLinSolvIters", 1)) return 1;
  retval = KINGetNumFuncEvals(kinsol_mem, &nfe);
  if (check_retval(&retval, "KinGetNumFuncEvals", 1)) return 1;

  cout << fixed;
  cout << setprecision(6);

  cout << "  Func evals       = " << nfe     << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << endl;

  return 0;
}

static int OutputTiming(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  if (outproc)
  {
    cout << scientific;
    cout << setprecision(6);
  }

  double maxtime = 0.0;

  MPI_Reduce(&(udata->totaltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Total time    = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->fevaltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  G(u) eval time = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->jvtime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  Jv time       = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->matfilltime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  MatFill time  = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->psetuptime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  PSetup time   = " << maxtime << " sec" << endl;
  }

  MPI_Reduce(&(udata->psolvetime), &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0,
             udata->comm_c);
  if (outproc)
  {
    cout << "  PSolve time   = " << maxtime << " sec" << endl;
    cout << endl;
  }

  return 0;
}

//Write solution to file to be plotted
static int WriteSolution(N_Vector u, UserData *udata)
{
  // Output problem information and open output streams
  // Each processor outputs subdomain information
  stringstream fname;
  fname << "heat2d_info." << setfill('0') << setw(5) << udata->myid_c
        << ".txt";

  ofstream dout;
  dout.open(fname.str());
  dout <<  "xu  " << udata->xu       << endl;
  dout <<  "yu  " << udata->yu       << endl;
  dout <<  "nx  " << udata->nx       << endl;
  dout <<  "ny  " << udata->ny       << endl;
  dout <<  "px  " << udata->npx      << endl;
  dout <<  "py  " << udata->npy      << endl;
  dout <<  "np  " << udata->nprocs_w << endl;
  dout <<  "is  " << udata->is       << endl;
  dout <<  "ie  " << udata->ie       << endl;
  dout <<  "js  " << udata->js       << endl;
  dout <<  "je  " << udata->je       << endl;
  dout <<  "nt  " << 1               << endl;
  dout.close();

  // Open output streams for solution
  fname.str("");
  fname.clear();
  fname << "heat2d_solution." << setfill('0') << setw(5) << udata->myid_c
        << ".txt";
  udata->uout.open(fname.str());

  udata->uout << scientific;
  udata->uout << setprecision(numeric_limits<realtype>::digits10);

  // Write solution and error to disk
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_retval((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  for (sunindextype i = 0; i < udata->nodes_loc; i++)
  {
    udata->uout << uarray[i] << " ";
  }
  udata->uout << endl;

  // Close output stream
  udata->uout.close();

  return 0;
}

static int OpenOutput(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  if (outproc)
  {
    stringstream fname;

    // Open output stream for residual
    fname.str("");
    fname.clear();
    fname << "heat2d_res_m" << udata->maa << "_orth" << udata->orthaa << ".txt";
    udata->rout.open(fname.str());

    udata->rout << scientific;
    udata->rout << setprecision(numeric_limits<realtype>::digits10);

    // Open output stream for error
    fname.str("");
    fname.clear();
    fname << "heat2d_err_m" << udata->maa << "_orth" << udata->orthaa << ".txt";
    udata->eout.open(fname.str());

    udata->eout << scientific;
    udata->eout << setprecision(numeric_limits<realtype>::digits10);
  }

  return 0;
}

static int WriteOutput(N_Vector u, N_Vector f, UserData *udata)
{
  int retval;
  bool outproc = (udata->myid_c == 0);

  // r = \|G(u) - u\|_2
  N_VLinearSum(ONE, f, -ONE, u, udata->e);
  realtype res = N_VDotProd(udata->e, udata->e);

  // e = \|u_exact - u\|_2
  retval = SolutionError(u, udata->e, udata);
  if (check_retval(&retval, "SolutionError", 1)) return 1;
  realtype err = N_VDotProd(udata->e, udata->e);

  if (outproc)
  {
    // Output residual
    udata->rout << sqrt(res);
    udata->rout << endl;

    // Output error
    udata->eout << sqrt(err);
    udata->eout << endl;
  }
  return 0;
}

static int CloseOutput(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

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
      cerr << endl << "ERROR: " << funcname << " returned "
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
