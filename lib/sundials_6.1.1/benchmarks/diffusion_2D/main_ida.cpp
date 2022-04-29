/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * IDA main for 2D diffusion benchmark problem
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"
#include "ida/ida.h"

struct UserOptions
{
  // Integrator settings
  realtype rtol        = RCONST(1.0e-5);   // relative tolerance
  realtype atol        = RCONST(1.0e-10);  // absolute tolerance
  int      maxsteps    = 0;                // max steps between outputs
  int      onestep     = 0;                // one step mode, number of steps

  // Linear solver and preconditioner settings
  bool     pcg      = true;   // use PCG (true) or GMRES (false)
  bool     prec     = true;   // preconditioner on/off
  bool     lsinfo   = false;  // output residual history
  int      liniters = 20;     // number of linear iterations
  realtype epslin   = ZERO;   // linear solver tolerance factor

  // Helper functions
  int parse_args(vector<string> &args, bool outproc);
  void help();
  void print();
};

// Print integration statistics and timings
static int OutputStats(void *ida_mem, UserData *udata);

// -----------------------------------------------------------------------------
// Main Program
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  // Reusable error-checking flag
  int flag;

  // Initialize MPI
  flag = MPI_Init(&argc, &argv);
  if (check_flag(&flag, "MPI_Init", 1)) return 1;

  // Create SUNDIALS context
  MPI_Comm    comm = MPI_COMM_WORLD;
  SUNContext  ctx  = NULL;
  SUNProfiler prof = NULL;

  flag = SUNContext_Create((void*) &comm, &ctx);
  if (check_flag(&flag, "SUNContextCreate", 1)) return 1;

  flag = SUNContext_GetProfiler(ctx, &prof);
  if (check_flag(&flag, "SUNContext_GetProfiler", 1)) return 1;

  // Add scope so objects are destroyed before MPI_Finalize
  {
    SUNDIALS_CXX_MARK_FUNCTION(prof);

    // MPI process ID
    int myid;
    flag = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (check_flag(&flag, "MPI_Comm_rank", 1)) return 1;

    bool outproc = (myid == 0);

    // --------------------------
    // Parse command line inputs
    // --------------------------

    UserData    udata(prof);
    UserOptions uopts;
    UserOutput  uout;

    vector<string> args(argv + 1, argv + argc);

    flag = udata.parse_args(args, outproc);
    if (check_flag(&flag, "UserData::parse_args", 1)) return 1;

    flag = uopts.parse_args(args, outproc);
    if (check_flag(&flag, "UserOptions::parse_args", 1)) return 1;

    flag = uout.parse_args(args, outproc);
    if (check_flag(&flag, "UserOutput::parse_args", 1)) return 1;

    // Check for unparsed inputs
    if (args.size() > 0)
    {
      if (find(args.begin(), args.end(), "--help") == args.end())
      {
        cerr << "ERROR: Unknown inputs: ";
        for (auto i = args.begin(); i != args.end(); ++i)
          cerr << *i << ' ';
        cerr << endl;
      }
      return 1;
    }

    // -----------------------------
    // Setup parallel decomposition
    // -----------------------------

    flag = udata.setup();
    if (check_flag(&flag, "UserData::setup", 1)) return 1;

    // Output problem setup/options
    FILE *diagfp = NULL;

    if (outproc)
    {
      udata.print();
      uopts.print();
      uout.print();

      // Open diagnostics output file
      if (uopts.lsinfo)
      {
        diagfp = fopen("diagnostics.txt", "w");
        if (check_flag((void *) diagfp, "fopen", 0)) return 1;
      }
    }

    // ---------------
    // Create vectors
    // ---------------

    // Create vector for solution
#if defined(USE_HIP)
    N_Vector u = N_VMake_MPIPlusX(udata.comm_c,
                                  N_VNew_Hip(udata.nodes_loc, ctx), ctx);
    if (check_flag((void *) u, "N_VMake_MPIPlusX", 0)) return 1;
#elif defined(USE_CUDA)
    N_Vector u = N_VMake_MPIPlusX(udata.comm_c,
                                  N_VNew_Cuda(udata.nodes_loc, ctx), ctx);
    if (check_flag((void *) u, "N_VMake_MPIPlusX", 0)) return 1;
#else
    N_Vector u = N_VNew_Parallel(udata.comm_c, udata.nodes_loc, udata.nodes,
                                 ctx);
    if (check_flag((void *) u, "N_VNew_Parallel", 0)) return 1;
#endif

    N_Vector up = N_VClone(u);
    if (check_flag((void *) up, "N_VClone", 0)) return 1;

    // Set initial condition
    flag = Solution(ZERO, u, &udata);
    if (check_flag(&flag, "Solution", 1)) return 1;

    flag = SolutionDerivative(ZERO, up, &udata);
    if (check_flag(&flag, "SolutionDerivative", 1)) return 1;

    // Create vector for error
    if (udata.forcing)
    {
      uout.error = N_VClone(u);
      if (check_flag((void *) (uout.error), "N_VClone", 0)) return 1;
    }

    // ---------------------
    // Create linear solver
    // ---------------------

    // Create linear solver
    SUNLinearSolver LS = NULL;

    int prectype = (uopts.prec) ? PREC_RIGHT : PREC_NONE;

    if (uopts.pcg)
    {
      LS = SUNLinSol_PCG(u, prectype, uopts.liniters, ctx);
      if (check_flag((void *) LS, "SUNLinSol_PCG", 0)) return 1;

      if (uopts.lsinfo && outproc)
      {
        flag = SUNLinSolSetPrintLevel_PCG(LS, 1);
        if (check_flag(&flag, "SUNLinSolSetPrintLevel_PCG", 1)) return(1);

        flag = SUNLinSolSetInfoFile_PCG(LS, diagfp);
        if (check_flag(&flag, "SUNLinSolSetInfoFile_PCG", 1)) return(1);
      }
    }
    else
    {
      LS = SUNLinSol_SPGMR(u, prectype, uopts.liniters, ctx);
      if (check_flag((void *) LS, "SUNLinSol_SPGMR", 0)) return 1;

      if (uopts.lsinfo && outproc)
      {
        flag = SUNLinSolSetPrintLevel_SPGMR(LS, 1);
        if (check_flag(&flag, "SUNLinSolSetPrintLevel_SPGMR", 1)) return(1);

        flag = SUNLinSolSetInfoFile_SPGMR(LS, diagfp);
        if (check_flag(&flag, "SUNLinSolSetInfoFile_SPGMR", 1)) return(1);
      }
    }

    // Allocate preconditioner workspace
    if (uopts.prec)
    {
      udata.diag = N_VClone(u);
      if (check_flag((void *) (udata.diag), "N_VClone", 0)) return 1;
    }

    // --------------
    // Setup IDA
    // --------------

    // Create integrator
    void *ida_mem = IDACreate(ctx);
    if (check_flag((void *) ida_mem, "IDACreate", 0)) return 1;

    // Initialize integrator
    flag = IDAInit(ida_mem, diffusion, ZERO, u, up);
    if (check_flag(&flag, "IDAInit", 1)) return 1;

    // Specify tolerances
    flag = IDASStolerances(ida_mem, uopts.rtol, uopts.atol);
    if (check_flag(&flag, "IDASStolerances", 1)) return 1;

    // Attach user data
    flag = IDASetUserData(ida_mem, (void *) &udata);
    if (check_flag(&flag, "IDASetUserData", 1)) return 1;

    // Attach linear solver
    flag = IDASetLinearSolver(ida_mem, LS, NULL);
    if (check_flag(&flag, "IDASetLinearSolver", 1)) return 1;

    if (uopts.prec)
    {
      // Attach preconditioner
      flag = IDASetPreconditioner(ida_mem, PSetup, PSolve);
      if (check_flag(&flag, "IDASetPreconditioner", 1)) return 1;
    }

    // Set linear solver tolerance factor
    flag = IDASetEpsLin(ida_mem, uopts.epslin);
    if (check_flag(&flag, "IDASetEpsLin", 1)) return 1;

    // Set max steps between outputs
    flag = IDASetMaxNumSteps(ida_mem, uopts.maxsteps);
    if (check_flag(&flag, "IDASetMaxNumSteps", 1)) return 1;

    // Set stopping time
    flag = IDASetStopTime(ida_mem, udata.tf);
    if (check_flag(&flag, "IDASetStopTime", 1)) return 1;

    // -----------------------
    // Loop over output times
    // -----------------------

    // Optionally run in one step mode for a fixed number of time steps (helpful
    // for debugging)
    int stepmode = IDA_NORMAL;

    if (uopts.onestep)
    {
      uout.nout = uopts.onestep;
      stepmode  = IDA_ONE_STEP;
    }

    realtype t     = ZERO;
    realtype dTout = udata.tf / uout.nout;
    realtype tout  = dTout;

    // Inital output
    flag = uout.open(&udata);
    if (check_flag(&flag, "UserOutput::open", 1)) return 1;

    flag = uout.write(t, u, &udata);
    if (check_flag(&flag, "UserOutput::write", 1)) return 1;

    for (int iout = 0; iout < uout.nout; iout++)
    {
      SUNDIALS_MARK_BEGIN(prof, "Evolve");

      // Evolve in time
      flag = IDASolve(ida_mem, tout, &t, u, up, stepmode);
      if (check_flag(&flag, "IDA", 1)) break;

      SUNDIALS_MARK_END(prof, "Evolve");

      // Output solution and error
      flag = uout.write(t, u, &udata);
      if (check_flag(&flag, "UserOutput::write", 1)) return 1;

      // Update output time
      tout += dTout;
      tout = (tout > udata.tf) ? udata.tf : tout;
    }

    // Close output
    flag = uout.close(&udata);
    if (check_flag(&flag, "UserOutput::close", 1)) return 1;

    // --------------
    // Final outputs
    // --------------

    // Print final integrator stats
    if (outproc)
    {
      cout << "Final integrator statistics:" << endl;
      flag = OutputStats(ida_mem, &udata);
      if (check_flag(&flag, "OutputStats", 1)) return 1;
    }

    // ---------
    // Clean up
    // ---------

    // Free MPI Cartesian communicator
    MPI_Comm_free(&(udata.comm_c));

    // Close diagnostics output file
    if (diagfp) fclose(diagfp);

    // Free integrator and linear solver
    IDAFree(&ida_mem);     // Free integrator memory
    SUNLinSolFree(LS);         // Free linear solver

    // Free vectors
#if defined(USE_HIP) || defined(USE_CUDA)
    N_VDestroy(N_VGetLocalVector_MPIPlusX(u));
#endif
    N_VDestroy(u);
  }
  // Close scope so objects are destroyed before MPI_Finalize

  // Stop timer and destroy context
    SUNContext_Free(&ctx);

  // Finalize MPI
  flag = MPI_Finalize();
  return 0;
}


// -----------------------------------------------------------------------------
// Output functions
// -----------------------------------------------------------------------------


// Print integrator statistics
static int OutputStats(void *ida_mem, UserData* udata)
{
  int flag;

  // Get integrator and solver stats
  long int nst, netf, nf, nni, ncfn, nli, nlcf, nsetups, nf_ls, nJv;
  flag = IDAGetNumSteps(ida_mem, &nst);
  if (check_flag(&flag, "IDAGetNumSteps", 1)) return -1;
  flag = IDAGetNumErrTestFails(ida_mem, &netf);
  if (check_flag(&flag, "IDAGetNumErrTestFails", 1)) return -1;
  flag = IDAGetNumResEvals(ida_mem, &nf);
  if (check_flag(&flag, "IDAGetNumResEvals", 1)) return -1;
  flag = IDAGetNumNonlinSolvIters(ida_mem, &nni);
  if (check_flag(&flag, "IDAGetNumNonlinSolvIters", 1)) return -1;
  flag = IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  if (check_flag(&flag, "IDAGetNumNonlinSolvConvFails", 1)) return -1;
  flag = IDAGetNumLinIters(ida_mem, &nli);
  if (check_flag(&flag, "IDAGetNumLinIters", 1)) return -1;
  flag = IDAGetNumLinConvFails(ida_mem, &nlcf);
  if (check_flag(&flag, "IDAGetNumLinConvFails", 1)) return -1;
  flag = IDAGetNumLinSolvSetups(ida_mem, &nsetups);
  if (check_flag(&flag, "IDAGetNumLinSolvSetups", 1)) return -1;
  flag = IDAGetNumLinResEvals(ida_mem, &nf_ls);
  if (check_flag(&flag, "IDAGetNumLinResEvals", 1)) return -1;
  flag = IDAGetNumJtimesEvals(ida_mem, &nJv);
  if (check_flag(&flag, "IDAGetNumJtimesEvals", 1)) return -1;

  cout << fixed;
  cout << setprecision(6);

  cout << "  Steps            = " << nst     << endl;
  cout << "  Error test fails = " << netf    << endl;
  cout << "  Res evals        = " << nf      << endl;
  cout << "  NLS iters        = " << nni     << endl;
  cout << "  NLS fails        = " << ncfn    << endl;
  cout << "  LS iters         = " << nli     << endl;
  cout << "  LS fails         = " << nlcf    << endl;
  cout << "  LS setups        = " << nsetups << endl;
  cout << "  LS Res evals     = " << nf_ls   << endl;
  cout << "  Jv products      = " << nJv     << endl;
  cout << endl;

  // Compute average nls iters per step attempt and ls iters per nls iter
  realtype avgnli = (realtype) nni / (realtype) nst;
  realtype avgli  = (realtype) nli / (realtype) nni;
  cout << "  Avg NLS iters per step    = " << avgnli << endl;
  cout << "  Avg LS iters per NLS iter = " << avgli  << endl;
  cout << endl;

  // Get preconditioner stats
  long int npe, nps;
  flag = IDAGetNumPrecEvals(ida_mem, &npe);
  if (check_flag(&flag, "IDAGetNumPrecEvals", 1)) return -1;
  flag = IDAGetNumPrecSolves(ida_mem, &nps);
  if (check_flag(&flag, "IDAGetNumPrecSolves", 1)) return -1;

  cout << "  Preconditioner setups = " << npe << endl;
  cout << "  Preconditioner solves = " << nps << endl;
  cout << endl;

  return 0;
}


// -----------------------------------------------------------------------------
// UserOptions Helper functions
// -----------------------------------------------------------------------------


int UserOptions::parse_args(vector<string> &args, bool outproc)
{
  vector<string>::iterator it;

  it = find(args.begin(), args.end(), "--help");
  if (it != args.end())
  {
    if (outproc) help();
    return 0;
  }

  it = find(args.begin(), args.end(), "--rtol");
  if (it != args.end())
  {
    rtol = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--atol");
  if (it != args.end())
  {
    atol = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--maxsteps");
  if (it != args.end())
  {
    maxsteps = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--onestep");
  if (it != args.end())
  {
    onestep = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--gmres");
  if (it != args.end())
  {
    pcg = false;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--lsinfo");
  if (it != args.end())
  {
    lsinfo = true;
    args.erase(it);
  }

  it = find(args.begin(), args.end(), "--liniters");
  if (it != args.end())
  {
    liniters = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--epslin");
  if (it != args.end())
  {
    epslin = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  return 0;
}


// Print command line options
void UserOptions::help()
{
  cout << endl;
  cout << "Integrator command line options:" << endl;
  cout << "  --rtol <rtol>       : relative tolerance" << endl;
  cout << "  --atol <atol>       : absoltue tolerance" << endl;
  cout << "  --gmres             : use GMRES linear solver" << endl;
  cout << "  --lsinfo            : output residual history" << endl;
  cout << "  --liniters <iters>  : max number of iterations" << endl;
  cout << "  --epslin <factor>   : linear tolerance factor" << endl;
  cout << "  --noprec            : disable preconditioner" << endl;
}


// Print user options
void UserOptions::print()
{
  cout << endl;
  cout << " Integrator options:" << endl;
  cout << " --------------------------------- " << endl;
  cout << " rtol        = " << rtol        << endl;
  cout << " atol        = " << atol        << endl;
  cout << " max steps   = " << maxsteps    << endl;
  cout << " --------------------------------- " << endl;

  cout << endl;
  cout << " Linear solver options:" << endl;
  cout << " --------------------------------- " << endl;
  cout << " PCG      = " << pcg      << endl;
  cout << " precond  = " << prec     << endl;
  cout << " LS info  = " << lsinfo   << endl;
  cout << " LS iters = " << liniters << endl;
  cout << " epslin   = " << epslin   << endl;
  cout << " --------------------------------- " << endl;
}

//---- end of file ----
