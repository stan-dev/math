/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
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
 * This demonstration problem simulates the advection and reaction of three
 * chemical species, u, v, and w, in a one dimensional domain. The reaction
 * mechanism is a variation of the Brusselator problem from chemical kinetics.
 * This is a PDE system with 3 components, Y = [u,v,w], satisfying the
 * equations,
 *
 *    u_t = -c * u_x + A - (w+1) * u + v * u^2
 *    v_t = -c * v_x + w * u - v * u^2
 *    w_t = -c * w_x + (B - w) / ep - w * u
 *
 * for t in [0, 10], x in [0, xmax] with periodic boundary conditions. The
 * initial condition is a Gaussian pertubation of the steady state
 * solution without advection
 *
 *    u(0,x) = k1 * A / k4 + p(x)
 *    v(0,x) = k2 * k4 * B / (k1 * k3 * A) + p(x)
 *    w(0,x) = 3.0 + p(x)
 *    p(x)   = alpha * e^( -(x - mu)^2 / (2*sigma^2) ).
 *
 * where alpha = 0.1, mu = xmax / 2.0, and sigma = xmax / 4.0.
 * The reaction rates are set so k_1 = k_2 = k_3 = k_4 = k, and k_5 = k_6 =
 * 1/5e-6. The spatial derivatives are discretized with first-order upwind
 * finite differences. An IMEX method is used to evolve the system in time with
 * the advection terms treated explicitly and the reaction terms implicitly. As
 * the reactions are purely local, the code uses a custom nonlinear solver to
 * exploit this locality so no parallel communication is needed in the implicit
 * solves. NOUT outputs are printed at equal intervals, and run statistics are
 * printed at the end.
 *
 * Command line options:
 *  --help           prints help message
 *  --printtime      print timing information
 *  --monitor        print solution information to screen (slower)
 *  --output-dir     the directory where all output files will be written
 *  --nout <int>     the number of output times
 *  --nx <int>       number of spatial mesh intervals
 *  --xmax <double>  maximum x value
 *  --explicit       use explicit method instead of IMEX
 *  --order <int>    method order
 *  --global-nls     use a global newton nonlinear solver instead of task-local
 *  --tf <double>    final time
 *  --A <double>     A parameter value
 *  --B <double>     B parameter value
 *  --k <double>     reaction rate
 *  --c <double>     advection speed
 *  --rtol <double>  relative tolerance
 *  --atol <double>  absolute tolerance
 * --------------------------------------------------------------------------*/

#include "ark_brusselator1D.h"

#ifdef USE_RAJA_NVEC
#define NVECTOR_ID_STRING "RAJA-unmanaged"
using EXEC_POLICY = RAJA::cuda_exec< 256, false >;
constexpr auto LocalNvector = N_VNew_Raja;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_Raja;

#elif USE_OMPDEV_NVEC
#define NVECTOR_ID_STRING "OpenMPDEV"
using EXEC_POLICY = RAJA::omp_target_parallel_for_exec< 256 >;
constexpr auto LocalNvector = N_VNew_OpenMPDEV;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_OpenMPDEV;

#elif USE_CUDAUVM_NVEC
#define NVECTOR_ID_STRING "CUDA-managed"
using EXEC_POLICY = RAJA::cuda_exec< 256, false >;
constexpr auto LocalNvector = N_VNewManaged_Cuda;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_Cuda;

#elif USE_CUDA_NVEC
#define NVECTOR_ID_STRING "CUDA-unmanaged"
using EXEC_POLICY = RAJA::cuda_exec< 256, false >;
constexpr auto LocalNvector = N_VNew_Cuda;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_Cuda;

#elif USE_HIP_NVEC
#define NVECTOR_ID_STRING "HIP-unmanaged"
using EXEC_POLICY = RAJA::hip_exec< 512, false >;
constexpr auto LocalNvector = N_VNew_Hip;
constexpr auto CopyVecFromDevice = N_VCopyFromDevice_Hip;

#else
#define NVECTOR_ID_STRING "Serial"
using EXEC_POLICY = RAJA::seq_exec;
constexpr auto LocalNvector = N_VNew_Serial;
#define CopyVecFromDevice(v)

#endif // USE_RAJA_NVEC

#ifdef USE_CUDA_OR_HIP
#define DEVICE_LAMBDA RAJA_DEVICE
#define GPU_SAFE_CALL(val) { gpuAssert((val), __FILE__, __LINE__, 1); }
#else
#define DEVICE_LAMBDA
#define GPU_SAFE_CALL
#endif


/* Main Program */
int main(int argc, char *argv[])
{
  /* Initialize MPI */
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  {
    /* general problem variables */
    N_Vector     y = NULL;      /* empty solution vector        */
    UserData     udata;         /* user data                    */
    UserOptions  uopt;          /* user options                 */
    int          retval;        /* reusable error-checking flag */
    int          i;             /* loop counter                 */
    FILE*        MFID;          /* mesh output file pointer     */
    char         fname[MXSTR];

    /* Start overall timing */
    udata.t_overall.start();

    /* Process input args and setup the problem */
    retval = SetupProblem(argc, argv, &udata, &uopt);
    if (check_retval(&retval, "SetupProblem", 1)) MPI_Abort(comm, 1);

    /* Create solution vector */
    y = N_VMake_MPIPlusX(udata.comm, LocalNvector(udata.NEQ));
    if (check_retval((void *) y, "N_VMake_MPIPlusX", 0)) MPI_Abort(comm, 1);

    /* Enabled fused vector ops */
    if (uopt.fused)
    {
      retval = EnableFusedVectorOps(y);
      if (check_retval(&retval, "EnableFusedVectorOps", 1)) MPI_Abort(comm, 1);
    }

    /* Set the initial condition */
    retval = SetIC(y, &udata);
    if (check_retval(&retval, "SetIC", 1)) MPI_Abort(comm, 1);

    /* Output spatial mesh to disk (add extra point for periodic BC) */
    if (udata.myid == 0 && uopt.nout > 0)
    {
      snprintf(fname, MXSTR, "%s/mesh.txt", uopt.outputdir);
      MFID = fopen(fname,"w");
      for (i = 0; i < udata.nx + 1; i++) fprintf(MFID,"  %.16e\n", udata.dx*i);
      fclose(MFID);
    }

    /* Integrate in time */
    if (uopt.expl) retval = EvolveProblemExplicit(y, &udata, &uopt);
    else           retval = EvolveProblemIMEX(y, &udata, &uopt);

    if (check_retval(&retval, "Evolve", 1)) MPI_Abort(comm, 1);

    /* End timing */
    udata.t_overall.stop();
    if (udata.myid == 0 && uopt.printtime)
    {
      printf("Setup wall clock time: %.4f seconds\n", udata.t_setup.total());
      printf("Communication wall clock time: %.4f seconds\n", udata.t_comm.total());
      printf("Advection wall clock time: %.4f seconds\n", udata.t_advec.total());
      printf("Reaction wall clock time: %.4f seconds\n", udata.t_react.total());
      if (uopt.global)
      {
        printf("PSolve wall clock time: %.4f seconds\n", udata.t_psolve.total());
      }
      else
      {
        printf("LSolve wall clock time: %.4f seconds\n", udata.t_lsolve.total());
      }
      printf("Total wall clock time: %.4f seconds\n", udata.t_overall.total());
    }

    /* Clean up */
    N_VDestroy(N_VGetLocalVector_MPIPlusX(y));
    N_VDestroy(y);
  }

  MPI_Finalize();
  return(0);
}


/* Setup ARKODE and evolve problem in time with IMEX method */
int EvolveProblemIMEX(N_Vector y, UserData* udata, UserOptions* uopt)
{
  void*              arkode_mem = NULL;  /* empty ARKODE memory structure    */
  SUNNonlinearSolver NLS = NULL;         /* empty nonlinear solver structure */
  SUNLinearSolver    LS  = NULL;         /* empty linear solver structure    */

  double   t, dtout, tout;    /* current/output time data     */
  int      retval;            /* reusable error-checking flag */
  int      iout;              /* output counter               */
  long int nst, nst_a, netf;  /* step stats                   */
  long int nfe, nfi;          /* RHS stats                    */
  long int nni, ncnf;         /* nonlinear solver stats       */
  long int nli, npsol;        /* linear solver stats          */
  FILE*    DFID = NULL;       /* diagnostics output file      */
  char     fname[MXSTR];

  /* Create the ARK timestepper module */
  arkode_mem = ARKStepCreate(Advection, Reaction, uopt->t0, y);
  if (check_retval((void*)arkode_mem, "ARKStepCreate", 0)) return 1;

  /* Select the method order */
  retval = ARKStepSetOrder(arkode_mem, uopt->order);
  if (check_retval(&retval, "ARKStepSetOrder", 1)) return 1;

  /* Attach user data */
  retval = ARKStepSetUserData(arkode_mem, (void*) udata);
  if (check_retval(&retval, "ARKStepSetUserData*", 1)) return 1;

  /* Specify tolerances */
  retval = ARKStepSStolerances(arkode_mem, uopt->rtol, uopt->atol);
  if (check_retval(&retval, "ARKStepSStolerances", 1)) return 1;

  /* Increase the max number of steps allowed between outputs */
  retval = ARKStepSetMaxNumSteps(arkode_mem, 100000);
  if (check_retval(&retval, "ARKStepSetMaxNumSteps", 1)) return 1;

  /* Open output file for integrator diagnostics */
  if (uopt->monitor)
  {
    sprintf(fname, "%s/diagnostics.%06d.txt", uopt->outputdir, udata->myid);
    DFID = fopen(fname, "w");

    retval = ARKStepSetDiagnostics(arkode_mem, DFID);
    if (check_retval(&retval, "ARKStepSetDiagnostics", 1)) return 1;
  }

  /* Create the (non)linear solver */
  if (uopt->global)
  {
    /* Create nonlinear solver */
    NLS = SUNNonlinSol_Newton(y);
    if (check_retval((void *)NLS, "SUNNonlinSol_Newton", 0)) return 1;

    /* Attach nonlinear solver */
    retval = ARKStepSetNonlinearSolver(arkode_mem, NLS);
    if (check_retval(&retval, "ARKStepSetNonlinearSolver", 1)) return 1;

    /* Create linear solver */
    LS = SUNLinSol_SPGMR(y, PREC_LEFT, 0);
    if (check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return 1;

    /* Attach linear solver */
    retval = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
    if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return 1;

    /* Attach preconditioner */
    retval = ARKStepSetPreconditioner(arkode_mem, NULL, PSolve);
    if (check_retval(&retval, "ARKStepSetPreconditioner", 1)) return 1;
  }
  else
  {
    /* The custom task-local nonlinear solver handles the linear solve
       as well, so we do not need a SUNLinearSolver. */
    NLS = TaskLocalNewton(y, DFID);
    if (check_retval((void *)NLS, "TaskLocalNewton", 0)) return 1;

    /* Attach nonlinear solver */
    retval = ARKStepSetNonlinearSolver(arkode_mem, NLS);
    if (check_retval(&retval, "ARKStepSetNonlinearSolver", 1)) return 1;
  }

  /* Output initial condition */
  if (udata->myid == 0 && uopt->monitor)
  {
    printf("\n          t         ||u||_rms   ||v||_rms   ||w||_rms\n");
    printf("   ----------------------------------------------------\n");
  }
  WriteOutput(uopt->t0, y, udata, uopt);

  /* Integrate to final time */
  t     = uopt->t0;
  dtout = (uopt->tf - uopt->t0);
  if (uopt->nout != 0)
    dtout /= uopt->nout;
  tout  = t + dtout;
  iout  = 0;

  do
  {
    /* Integrate to output time */
    retval = ARKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "ARKStepEvolve", 1)) break;

    /* Output state */
    WriteOutput(t, y, udata, uopt);

    /* Update output time */
    tout += dtout;
    tout = (tout > uopt->tf) ? uopt->tf : tout;

    iout++;
  } while (iout < uopt->nout);

  /* close output stream */
  if (uopt->monitor) fclose(DFID);

  /* Get final statistics */
  retval = ARKStepGetNumSteps(arkode_mem, &nst);
  check_retval(&retval, "ARKStepGetNumSteps", 1);
  retval = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_retval(&retval, "ARKStepGetNumStepAttempts", 1);
  retval = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_retval(&retval, "ARKStepGetNumRhsEvals", 1);
  retval = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  check_retval(&retval, "ARKStepGetNumErrTestFails", 1);
  retval = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  check_retval(&retval, "ARKStepGetNumNonlinSolvIters", 1);
  retval = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncnf);
  check_retval(&retval, "ARKStepGetNumNonlinSolvConvFails", 1);
  if (uopt->global)
  {
    retval = ARKStepGetNumLinIters(arkode_mem, &nli);
    check_retval(&retval, "ARKStepGetNumLinIters", 1);
    retval = ARKStepGetNumPrecSolves(arkode_mem, &npsol);
    check_retval(&retval, "ARKStepGetNumPrecSolves", 1);
  }

  /* Print final statistics */
  if (udata->myid == 0)
  {
    printf("\nFinal Solver Statistics (for processor 0):\n");
    printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
    printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe,
           nfi + udata->nnlfi);
    printf("   Total number of error test failures = %li\n", netf);
    printf("   Total number of nonlinear solver convergence failures = %li\n",
           ncnf);
    if (uopt->global)
    {
      printf("   Total number of nonlinear iterations = %li\n", nni);
      printf("   Total number of linear iterations = %li\n", nli);
      printf("   Total number of preconditioner solves = %li\n", npsol);
    }
  }

  /* Clean up */
  ARKStepFree(&arkode_mem);
  SUNNonlinSolFree(NLS);
  if (LS) SUNLinSolFree(LS);

  /* Return success */
  return(0);
}


/* Setup ARKODE and evolve problem in time explicitly */
int EvolveProblemExplicit(N_Vector y, UserData* udata, UserOptions* uopt)
{
  void*    arkode_mem = NULL; /* empty ARKODE memory structure */
  double   t, dtout, tout;    /* current/output time data      */
  int      retval;            /* reusable error-checking flag  */
  int      iout;              /* output counter                */
  long int nst, nst_a, netf;  /* step stats                    */
  long int nfe;               /* RHS stats                     */
  FILE*    DFID;              /* diagnostics output file       */
  char     fname[MXSTR];

  /* Create the ERK timestepper module */
  arkode_mem = ERKStepCreate(AdvectionReaction, uopt->t0, y);
  if (check_retval((void*)arkode_mem, "ERKStepCreate", 0)) return 1;

  /* Select the method order */
  retval = ERKStepSetOrder(arkode_mem, uopt->order);
  if (check_retval(&retval, "ERKStepSetOrder", 1)) return 1;

  /* Attach user data */
  retval = ERKStepSetUserData(arkode_mem, (void*) udata);
  if (check_retval(&retval, "ERKStepSetUserData", 1)) return 1;

  /* Specify tolerances */
  retval = ERKStepSStolerances(arkode_mem, uopt->rtol, uopt->atol);
  if (check_retval(&retval, "ERKStepSStolerances", 1)) return 1;

  /* Increase the max number of steps allowed between outputs */
  retval = ERKStepSetMaxNumSteps(arkode_mem, 1000000);
  if (check_retval(&retval, "ERKStepSetMaxNumSteps", 1)) return 1;

  /* Open output file for integrator diagnostics */
  if (uopt->monitor)
  {
    sprintf(fname, "%s/diagnostics.%06d.txt", uopt->outputdir, udata->myid);
    DFID = fopen(fname, "w");

    retval = ERKStepSetDiagnostics(arkode_mem, DFID);
    if (check_retval(&retval, "ERKStepSetDiagnostics", 1)) return 1;
  }

  /* Output initial condition */
  if (udata->myid == 0 && uopt->monitor)
  {
    printf("\n          t         ||u||_rms   ||v||_rms   ||w||_rms\n");
    printf("   ----------------------------------------------------\n");
  }
  WriteOutput(uopt->t0, y, udata, uopt);

  /* Integrate to final time */
  t     = uopt->t0;
  dtout = (uopt->tf - uopt->t0);
  if (uopt->nout != 0)
    dtout /= uopt->nout;
  tout  = t + dtout;
  iout  = 0;

  do
  {
    /* Integrate to output time */
    retval = ERKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_retval(&retval, "ERKStepEvolve", 1)) break;

    /* Output state */
    WriteOutput(t, y, udata, uopt);

    /* Update output time */
    tout += dtout;
    tout = (tout > uopt->tf) ? uopt->tf : tout;

    iout++;
  } while (iout < uopt->nout);

  /* close output stream */
  if (uopt->monitor) fclose(DFID);

  /* Get final statistics */
  retval = ERKStepGetNumSteps(arkode_mem, &nst);
  check_retval(&retval, "ERKStepGetNumSteps", 1);
  retval = ERKStepGetNumStepAttempts(arkode_mem, &nst_a);
  check_retval(&retval, "ERKStepGetNumStepAttempts", 1);
  retval = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
  check_retval(&retval, "ERKStepGetNumRhsEvals", 1);
  retval = ERKStepGetNumErrTestFails(arkode_mem, &netf);
  check_retval(&retval, "ERKStepGetNumErrTestFails", 1);

  /* Print final statistics */
  if (udata->myid == 0)
  {
    printf("\nFinal Solver Statistics (for processor 0):\n");
    printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
    printf("   Total RHS evals:  Fe = %li", nfe);
    printf("   Total number of error test failures = %li\n", netf);
  }

  /* Clean up */
  ERKStepFree(&arkode_mem);

  /* Return success */
  return(0);
}


/* Write time and solution to disk */
int WriteOutput(double t, N_Vector y, UserData* udata, UserOptions* uopt)
{
  int     nvar = udata->nvar;
  double  u, v, w, N;
  double* data = NULL;

  /* get left end point for output */
  ExchangeBCOnly(y, udata);

  /* get vector data array */
  data = N_VGetArrayPointer(y);
  if (check_retval((void *) data, "N_VGetArrayPointer", 0)) return 1;

  CopyVecFromDevice(N_VGetLocalVector_MPIPlusX(y));

  /* output current solution norm to screen */
  if (uopt->monitor)
  {
    N = (double) udata->nx;
    u = N_VWL2Norm(y,udata->umask);
    u = sqrt(u*u/N);
    v = N_VWL2Norm(y,udata->vmask);
    v = sqrt(v*v/N);
    w = N_VWL2Norm(y,udata->wmask);
    w = sqrt(w*w/N);
    if (udata->myid == 0)
      printf("     %10.6f   %10.6f   %10.6f   %10.6f\n", t, u, v, w);
  }

  if (uopt->nout > 0)
  {
      /* output the times to disk */
      if (udata->myid == 0 && udata->TFID)
        fprintf(udata->TFID," %.16e\n", t);

      /* output results to disk */
      for (int i = 0; i < udata->nxl; i++)
      {
        fprintf(udata->UFID," %.16e", data[IDX(nvar, i, 0)]);
        fprintf(udata->VFID," %.16e", data[IDX(nvar, i, 1)]);
        fprintf(udata->WFID," %.16e", data[IDX(nvar, i, 2)]);
      }

      /* we have one extra output because of the periodic BCs */
      if (udata->myid == (udata->nprocs - 1))
      {
        fprintf(udata->UFID," %.16e\n", udata->Erecv[IDX(nvar, 0, 0)]);
        fprintf(udata->VFID," %.16e\n", udata->Erecv[IDX(nvar, 0, 1)]);
        fprintf(udata->WFID," %.16e\n", udata->Erecv[IDX(nvar, 0, 2)]);
      }
      else
      {
        fprintf(udata->UFID,"\n");
        fprintf(udata->VFID,"\n");
        fprintf(udata->WFID,"\n");
      }
  }

  return(0);
}


/* Initial Condition Functions */
int SetIC(N_Vector y, UserData* udata)
{
  /* Variable shortcuts */
  const int    nvar = udata->nvar;
  const int    N    = udata->nxl;
  const double dx   = udata->dx;
  const double xmax = udata->xmax;
  const double A    = udata->A;
  const double B    = udata->B;
  const double k1   = udata->k1;
  const double k2   = udata->k2;
  const double k3   = udata->k3;
  const double k4   = udata->k4;
  const int    myid = udata->myid;

  /* Steady state solution */
  const double us = k1 * A / k4;
  const double vs = k2 * k4 * B / (k1 * k3 * A);
  const double ws = 3.0;

  /* Gaussian distribution defaults */
  const double mu    = xmax / 2.0;
  const double sigma = xmax / 4.0;
  const double alpha = 0.1;

  /* Gaussian perturbation of the steady state solution */
  double* data = GetVecData(y);
  RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, N),
    [=] DEVICE_LAMBDA (int i) {
    const double x = (myid * N + i) * dx;
    const double p = alpha * exp( -((x - mu) * (x - mu)) / (2.0 * sigma * sigma) );
    data[IDX(nvar, i, 0)] = us + p;
    data[IDX(nvar, i, 1)] = vs + p;
    data[IDX(nvar, i, 2)] = ws + p;
  });

  /* Return success */
  return(0);
}


/* --------------------------------------------------------------
 * Right Hand Side (RHS) Functions
 * --------------------------------------------------------------*/

/* Compute the advection term. */
int Advection(double t, N_Vector y, N_Vector ydot, void* user_data)
{
  /* access problem data */
  UserData* udata = (UserData*) user_data;

  /* time this function */
  udata->t_advec.start();

  /* set variable shortcuts */
  const int    nvar = udata->nvar;
  const int    N    = udata->nxl;
  const double dx   = udata->dx;
  const double c    = udata->c;
  double* Wrecv     = udata->Wrecv;
  double* Erecv     = udata->Erecv;

  /* local variables */
  double* Ydata  = NULL;
  double* dYdata = NULL;
  double  tmp;
  int     retval;

  /* set output to zero */
  N_VConst(0.0, ydot);

  /* access data arrays */
  Ydata = GetVecData(y);
  if (check_retval((void *)Ydata, "GetVecData", 0)) return(-1);

  dYdata = GetVecData(ydot);
  if (check_retval((void *)dYdata, "GetVecData", 0)) return(-1);

  /* begin exchanging boundary information */
  retval = ExchangeAllStart(y, udata);
  if (check_retval(&retval, "ExchangeAllStart", 1)) return(-1);

  /* finish exchanging boundary information */
  retval = ExchangeAllEnd(udata);
  if (check_retval(&retval, "ExchangeAllEnd", 1)) return(-1);

  /* iterate over domain interior, computing advection */
  tmp = -c / dx;

  if (c > 0.0)
  {
    /* right moving flow */
    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(1, N),
      [=] DEVICE_LAMBDA (int i) {
      dYdata[IDX(nvar, i, 0)] = tmp * (Ydata[IDX(nvar, i, 0)]
                                        - Ydata[IDX(nvar, i-1, 0)]);
      dYdata[IDX(nvar, i, 1)] = tmp * (Ydata[IDX(nvar, i, 1)]
                                        - Ydata[IDX(nvar, i-1, 1)]);
      dYdata[IDX(nvar, i, 2)] = tmp * (Ydata[IDX(nvar, i, 2)]
                                        - Ydata[IDX(nvar, i-1, 2)]);
    });
  }
  else if (c < 0.0)
  {
    /* left moving flow */
    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, N-1),
      [=] DEVICE_LAMBDA (int i) {
      dYdata[IDX(nvar, i, 0)] = tmp * (Ydata[IDX(nvar, i+1, 0)]
                                        - Ydata[IDX(nvar, i, 0)]);
      dYdata[IDX(nvar, i, 1)] = tmp * (Ydata[IDX(nvar, i+1, 1)]
                                        - Ydata[IDX(nvar, i, 1)]);
      dYdata[IDX(nvar, i, 2)] = tmp * (Ydata[IDX(nvar, i+1, 2)]
                                        - Ydata[IDX(nvar, i, 2)]);
    });
  }

  /* compute advection at local boundaries */
  if (c > 0.0)
  {
    /* right moving flow (left boundary) */
    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, nvar),
      [=] DEVICE_LAMBDA (int var) {
      dYdata[IDX(nvar, 0, var)] = tmp * (Ydata[IDX(nvar, 0, var)]
                                          - Wrecv[IDX(nvar, 0, var)]);
    });
  }
  else if (c < 0.0)
  {
    /* left moving flow (right boundary) */
    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, nvar),
      [=] DEVICE_LAMBDA (int var) {
      dYdata[IDX(nvar, N-1, var)] = tmp * (Erecv[IDX(nvar, 0, var)]
                                          - Ydata[IDX(nvar, N-1, var)]);
    });
  }

  /* stop timing */
  udata->t_advec.stop();

  /* return success */
  return(0);
}


/* Compute the reaction term. */
int Reaction(double t, N_Vector y, N_Vector ydot, void* user_data)
{
  /* access problem data */
  UserData* udata = (UserData*) user_data;

  /* time this function */
  udata->t_react.start();

  /* set variable shortcuts */
  const int    nvar = udata->nvar;
  const int    N    = udata->nxl;
  const double A    = udata->A;
  const double B    = udata->B;
  const double k1   = udata->k1;
  const double k2   = udata->k2;
  const double k3   = udata->k3;
  const double k4   = udata->k4;
  const double k5   = udata->k5;
  const double k6   = udata->k6;

  /* local variables */
  double* Ydata  = NULL;
  double* dYdata = NULL;

  /* access data arrays */
  Ydata = GetVecData(y);
  if (check_retval((void *)Ydata, "GetVecData", 0)) return(-1);

  dYdata = GetVecData(ydot);
  if (check_retval((void *)dYdata, "GetVecData", 0)) return(-1);

  /* iterate over domain, computing reactions */
  if (udata->uopt->expl)
  {
    /* when integrating explicitly, we add to ydot as we expect it
       to hold the advection term already */
    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, N),
      [=] DEVICE_LAMBDA (int i) {
      double u = Ydata[IDX(nvar, i ,0)];
      double v = Ydata[IDX(nvar, i, 1)];
      double w = Ydata[IDX(nvar, i, 2)];
      dYdata[IDX(nvar, i, 0)] += k1 * A - k2 * w * u + k3 * u * u * v - k4 * u;
      dYdata[IDX(nvar, i, 1)] += k2 * w * u - k3 * u * u * v;
      dYdata[IDX(nvar, i, 2)] += -k2 * w * u + k5 * B - k6 * w;
    });
  }
  else
  {
    /* set output to zero */
    N_VConst(0.0, ydot);

    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, N),
      [=] DEVICE_LAMBDA (int i) {
      double u = Ydata[IDX(nvar, i ,0)];
      double v = Ydata[IDX(nvar, i, 1)];
      double w = Ydata[IDX(nvar, i, 2)];
      dYdata[IDX(nvar, i, 0)] = k1 * A - k2 * w * u + k3 * u * u * v - k4 * u;
      dYdata[IDX(nvar, i, 1)] = k2 * w * u - k3 * u * u * v;
      dYdata[IDX(nvar, i, 2)] = -k2 * w * u + k5 * B - k6 * w;
    });
  }

  /* stop timing */
  udata->t_react.stop();

  /* return success */
  return(0);
}


/* Compute the RHS as Advection+Reaction. */
int AdvectionReaction(double t, N_Vector y, N_Vector ydot,
                      void *user_data)
{
  int retval;

  /* NOTE: The order in which Advection and Reaction are
           called is critical here. Advection must be
           computed first. */
  retval = Advection(t, y, ydot, user_data);
  if (check_retval((void *)&retval, "Advection", 1)) return(-1);

  retval = Reaction(t, y, ydot, user_data);
  if (check_retval((void *)&retval, "Reaction", 1)) return(-1);

  /* return success */
  return(0);
}


/* --------------------------------------------------------------
 * Linear solver functions
 * --------------------------------------------------------------*/

int SolveJacBlocks(N_Vector y, N_Vector x, N_Vector b,
                   double gamma, RAJA::RangeSegment blocks,
                   UserData* udata)
{
  /* shortcuts */
  int     nvar;
  double  k2, k3, k4, k6;
  double* Ydata;
  double* Bdata;
  double* Xdata;

  /* set shortcuts */
  nvar  = udata->nvar;
  k2    = udata->k2;
  k3    = udata->k3;
  k4    = udata->k4;
  k6    = udata->k6;
  Ydata = GetVecData(y);
  Bdata = GetVecData(b);
  Xdata = GetVecData(x);

  RAJA::forall< EXEC_POLICY >( blocks,
    [=] DEVICE_LAMBDA (int i) {

    /* and the corresponding vectors */
    double *b = &(Bdata[IDX(nvar, i, 0)]);
    double *x = &(Xdata[IDX(nvar, i, 0)]);

    /* shortcuts to u, v, w for the block */
    double u = Ydata[IDX(nvar, i, 0)];
    double v = Ydata[IDX(nvar, i, 1)];
    double w = Ydata[IDX(nvar, i, 2)];

    double A0, A1, A2, A3, A4, A5, A6, A7, A8;

    //
    // compute A = df/dz
    //

    /* 1st row: u, v, w */
    A0 = -k2 * w + 2.0 * k3 * u * v - k4;
    A1 =  k3 * u * u;
    A2 = -k2 * u;

    /* 2nd row: u, v, w */
    A3 =  k2 * w - 2.0 * k3 * u * v;
    A4 = -k3 * u * u;
    A5 =  k2 * u;

    /* 3rd row: u, v, w */
    A6 = -k2 * w;
    A7 =  0.0;
    A8 = -k2 * u - k6;

    //
    // compute A = I - gamma*J
    //

    A0 = 1. - (gamma * A0);
    A1 = -gamma * A1;
    A2 = -gamma * A2;
    A3 = -gamma * A3;
    A4 = 1. - (gamma * A4);
    A5 = -gamma * A5;
    A6 = -gamma * A6;
    A7 = -gamma * A7;
    A8 = 1. - (gamma * A8);

    //
    // compute x = A^{-1}b
    //

    double scratch_0 = A4*A8;
    double scratch_1 = A1*A5;
    double scratch_2 = A2*A7;
    double scratch_3 = A5*A7;
    double scratch_4 = A1*A8;
    double scratch_5 = A2*A4;
    double scratch_6 = 1.0/(A0*scratch_0 - A0*scratch_3 + A3*scratch_2 - A3*scratch_4 + A6*scratch_1 - A6*scratch_5);
    double scratch_7 = A2*A3;
    double scratch_8 = A6*b[0];
    double scratch_9 = A2*A6;
    double scratch_10 = A3*b[0];
    double scratch_11 = 1.0/A0;
    double scratch_12 = A1*scratch_11;
    double scratch_13 = (-A6*scratch_12 + A7)/(-A3*scratch_12 + A4);

    x[0] = scratch_6*(b[0]*scratch_0 - b[0]*scratch_3 + b[1]*scratch_2 - b[1]*scratch_4 + b[2]*scratch_1 - b[2]*scratch_5);
    x[1] = scratch_6*(-A0*A5*b[2] + A0*A8*b[1] + A5*scratch_8 - A8*scratch_10 - b[1]*scratch_9 + b[2]*scratch_7);
    x[2] = (-b[2] + scratch_11*scratch_8 + scratch_13*(b[1] - scratch_10*scratch_11))/(-A8 + scratch_11*scratch_9 + scratch_13*(A5 - scratch_11*scratch_7));
  });

  return(0);
}


/* --------------------------------------------------------------
 * (Non)linear system functions
 * --------------------------------------------------------------*/

int TaskLocalNlsResidual(N_Vector ycor, N_Vector F, void* arkode_mem)
{
  /* temporary variables */
  UserData* udata;
  int      retval;
  double   c[3];
  N_Vector X[3];

  /* nonlinear system data */
  N_Vector z, zpred, Fi, sdata;
  double   tcur, gamma;
  void     *user_data;

  retval = ARKStepGetNonlinearSystemData(arkode_mem, &tcur, &zpred, &z, &Fi,
                                         &gamma, &sdata, &user_data);
  if (check_retval((void *)&retval, "ARKStepGetNonlinearSystemData", 1))
    return(-1);

  udata = (UserData*) user_data;

  /* update 'z' value as stored predictor + current corrector */
  N_VLinearSum(1.0, N_VGetLocalVector_MPIPlusX(zpred),
               1.0, (ycor),
               N_VGetLocalVector_MPIPlusX(z));

  /* compute implicit RHS and save for later */
  retval = Reaction(tcur,
                    N_VGetLocalVector_MPIPlusX(z),
                    N_VGetLocalVector_MPIPlusX(Fi),
                    user_data);
  udata->nnlfi++; /* count calls to Fi as part of the nonlinear residual */
  if (retval < 0) return(-1);
  if (retval > 0) return(+1);

  /* update with y, sdata, and gamma * fy */
  X[0] = ycor;
  c[0] = 1.0;
  c[1] = -1.0;
  X[1] = N_VGetLocalVector_MPIPlusX(sdata);
  c[2] = -gamma;
  X[2] = N_VGetLocalVector_MPIPlusX(Fi);

  retval = N_VLinearCombination(3, c, X, F);
  if (retval != 0) return(-1);

  return(0);
}


int TaskLocalLSolve(N_Vector delta, void* arkode_mem)
{
  /* local variables */
  UserData* udata = NULL;
  int       retval;

  /* nonlinear system data */
  N_Vector z, zpred, Fi, sdata;
  double   tcur, gamma;
  void*    user_data = NULL;

  retval = ARKStepGetNonlinearSystemData(arkode_mem, &tcur, &zpred, &z, &Fi,
                                         &gamma, &sdata, &user_data);
  if (check_retval((void *)&retval, "ARKStepGetNonlinearSystemData", 1))
    return(-1);

  udata = (UserData*) user_data;

  /* time this function */
  udata->t_lsolve.start();

  /* set up I - gamma*J and solve */
  retval = SolveJacBlocks(z, delta, delta, gamma,
                          RAJA::RangeSegment(0, udata->nxl),
                          udata);

  /* stop timing */
  udata->t_lsolve.stop();

  return(retval);
}


SUNNonlinearSolver_Type TaskLocalNewton_GetType(SUNNonlinearSolver NLS)
{
  return SUNNONLINEARSOLVER_ROOTFIND;
}


int TaskLocalNewton_Initialize(SUNNonlinearSolver NLS)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  /* override default system and lsolve functions with local versions */
  SUNNonlinSolSetSysFn(LOCAL_NLS(NLS), TaskLocalNlsResidual);
  SUNNonlinSolSetLSolveFn(LOCAL_NLS(NLS), TaskLocalLSolve);

  return(SUNNonlinSolInitialize(LOCAL_NLS(NLS)));
}


int TaskLocalNewton_Solve(SUNNonlinearSolver NLS,
                          N_Vector y0, N_Vector ycor,
                          N_Vector w, double tol,
                          booleantype callLSetup, void* mem)
{
  /* local variables */
  MPI_Comm comm;
  int solve_status, recover, nonrecover;

  /* check that the inputs are non-null */
  if ((NLS  == NULL) ||
      (y0   == NULL) ||
      (ycor == NULL) ||
      (w    == NULL) ||
      (mem  == NULL))
    return SUN_NLS_MEM_NULL;

  /* shortcuts */
  comm = GET_NLS_CONTENT(NLS)->comm;

  /* each tasks solves the local nonlinear system */
  solve_status = SUNNonlinSolSolve(LOCAL_NLS(NLS),
                                   N_VGetLocalVector_MPIPlusX(y0),
                                   N_VGetLocalVector_MPIPlusX(ycor),
                                   N_VGetLocalVector_MPIPlusX(w),
                                   tol, callLSetup, mem);

  /* if any process had a nonrecoverable failure, return it */
  MPI_Allreduce(&solve_status, &nonrecover, 1, MPI_INT, MPI_MIN, comm);
  if (nonrecover < 0) return nonrecover;

  /* check if any process has a recoverable convergence failure */
  MPI_Allreduce(&solve_status, &recover, 1, MPI_INT, MPI_MAX, comm);
  if (recover == SUN_NLS_CONV_RECVR) GET_NLS_CONTENT(NLS)->ncnf++;

  /* return success (recover == 0) or a recoverable error code (recover > 0) */
  return recover;
}


int TaskLocalNewton_Free(SUNNonlinearSolver NLS)
{
  /* return if NLS is already free */
  if (NLS == NULL)
    return SUN_NLS_SUCCESS;

  /* free items from contents, then the generic structure */
  if (NLS->content)
  {
    SUNNonlinSolFree(LOCAL_NLS(NLS));
    free(NLS->content);
    NLS->content = NULL;
  }

  /* free the ops structure */
  if (NLS->ops)
  {
    free(NLS->ops);
    NLS->ops = NULL;
  }

  /* free the nonlinear solver */
  free(NLS);

  return SUN_NLS_SUCCESS;
}


int TaskLocalNewton_SetSysFn(SUNNonlinearSolver NLS,
                             SUNNonlinSolSysFn SysFn)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  return(SUNNonlinSolSetSysFn(LOCAL_NLS(NLS), SysFn));
}


int TaskLocalNewton_SetConvTestFn(SUNNonlinearSolver NLS,
                                  SUNNonlinSolConvTestFn CTestFn,
                                  void* ctest_data)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  return(SUNNonlinSolSetConvTestFn(LOCAL_NLS(NLS), CTestFn, ctest_data));
}


int TaskLocalNewton_GetNumConvFails(SUNNonlinearSolver NLS,
                                    long int *nconvfails)
{
  /* check that the nonlinear solver is non-null */
  if (NLS == NULL)
    return SUN_NLS_MEM_NULL;

  *nconvfails = GET_NLS_CONTENT(NLS)->ncnf;
  return(0);
}


SUNNonlinearSolver TaskLocalNewton(N_Vector y, FILE* DFID)
{
  SUNNonlinearSolver NLS;
  TaskLocalNewton_Content content;

  /* Check that the supplied N_Vector is non-NULL */
  if (y == NULL) return NULL;

  /* Check that the supplied N_Vector is an MPIPlusX */
  if (N_VGetVectorID(y) != SUNDIALS_NVEC_MPIPLUSX)
    return NULL;

  /* Create an empty nonlinear linear solver object */
  NLS = NULL;
  NLS = SUNNonlinSolNewEmpty();
  if (NLS == NULL) return NULL;

  /* Attach operations */
  NLS->ops->gettype         = TaskLocalNewton_GetType;
  NLS->ops->initialize      = TaskLocalNewton_Initialize;
  NLS->ops->solve           = TaskLocalNewton_Solve;
  NLS->ops->free            = TaskLocalNewton_Free;
  NLS->ops->setsysfn        = TaskLocalNewton_SetSysFn;
  NLS->ops->setctestfn      = TaskLocalNewton_SetConvTestFn;
  NLS->ops->getnumconvfails = TaskLocalNewton_GetNumConvFails;

  /* Create content */
  content = NULL;
  content = (TaskLocalNewton_Content) malloc(sizeof *content);
  if (content == NULL) { SUNNonlinSolFree(NLS); return NULL; }

  /* Initialize all components of content to 0/NULL */
  memset(content, 0, sizeof(*content));

  /* Attach content */
  NLS->content = content;

  /* Fill general content */
  content->comm = *((MPI_Comm*) N_VGetCommunicator(y));
  if (content->comm == NULL) { SUNNonlinSolFree(NLS); return NULL; }

  content->local_nls = SUNNonlinSol_Newton(N_VGetLocalVector_MPIPlusX(y));
  if (content->local_nls == NULL) { SUNNonlinSolFree(NLS); return NULL; }

  MPI_Comm_rank(content->comm, &content->myid);
  MPI_Comm_size(content->comm, &content->nprocs);

  content->ncnf = 0;

  /* Setup the local nonlinear solver monitoring */
  if (DFID != NULL)
  {
    SUNNonlinSolSetInfoFile_Newton(LOCAL_NLS(NLS), DFID);
    SUNNonlinSolSetPrintLevel_Newton(LOCAL_NLS(NLS), 1);
  }

  return NLS;
}


/* --------------------------------------------------------------
 * Preconditioner functions
 * --------------------------------------------------------------*/

/* Solves Pz = r */
int PSolve(double t, N_Vector y, N_Vector ydot, N_Vector r,
           N_Vector z, double gamma, double delta, int lr,
           void *user_data)
{
  /* local variables */
  UserData* udata = (UserData*) user_data;
  int       retval;

  /* time this function */
  udata->t_psolve.start();

  /* solve the task-local linear system Pz = r */
  retval = SolveJacBlocks(y, z, r, gamma,
                          RAJA::RangeSegment(0, udata->nxl),
                          udata);

  /* stop timing */
  udata->t_psolve.stop();

  return(retval);
}


/* --------------------------------------------------------------
 * Utility functions
 * --------------------------------------------------------------*/

/* Exchanges the periodic BCs only by sending the first
   mesh node to the last processor. */
int ExchangeBCOnly(N_Vector y, UserData* udata)
{
  int ierr;
  MPI_Status stat;
  MPI_Request reqR, reqS;

  /* shortcuts */
  int nvar  = udata->nvar;
  int myid  = udata->myid;
  int first = 0;
  int last  = udata->nprocs - 1;

  /* extract the data */
  double* Ydata = GetVecData(y);
  double* Wsend = udata->Wsend;

  /* open the East Irecv buffer */
  if (myid == last)
  {
    ierr = MPI_Irecv(udata->Erecv, nvar, MPI_DOUBLE, first,
                     MPI_ANY_TAG, udata->comm, &reqR);
  }

  /* send first mesh node to the last processor */
  if (myid == first)
  {
    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, nvar),
      [=] DEVICE_LAMBDA (int var) {
      Wsend[IDX(nvar, 0, var)] = Ydata[IDX(nvar, 0, var)];
    });
    ierr = MPI_Isend(udata->Wsend, nvar, MPI_DOUBLE,
                     last, 0, udata->comm, &reqS);
  }

  if (myid == last)
  {
    /* wait for exchange to finish */
    ierr = MPI_Wait(&reqR, &stat);
    if (ierr != MPI_SUCCESS)
    {
      fprintf(stderr, "\nERROR: error in MPI_Wait = %d\n", ierr);
      return -1;
    }
  }

  if (myid == first)
  {
    /* wait for exchange to finish */
    ierr = MPI_Wait(&reqS, &stat);
    if (ierr != MPI_SUCCESS)
    {
      fprintf(stderr, "\nERROR: error in MPI_Wait = %d\n", ierr);
      return -1;
    }
  }

  return(0);
}


/* Starts the exchange of the neighbor information */
int ExchangeAllStart(N_Vector y, UserData* udata)
{
  int retval;

  udata->t_comm.start();

  /* shortcuts */
  double c     = udata->c;
  int    N     = udata->nxl;
  int    nvar  = udata->nvar;
  int    myid  = udata->myid;
  int    first = 0;
  int    last  = udata->nprocs - 1;
  int    ipW   = (myid == first) ? last : udata->myid - 1; /* periodic BC */
  int    ipE   = (myid == last) ? first : udata->myid + 1; /* periodic BC */

  /* extract the data */
  double* Ydata = GetVecData(y);
  double* Esend = udata->Esend;
  double* Wsend = udata->Wsend;

  if (c > 0.0)
  {
    /* Right moving flow uses backward difference.
       Send from west to east (last processor sends to first) */
    retval = MPI_Irecv(udata->Wrecv, nvar, MPI_DOUBLE, ipW,
                       MPI_ANY_TAG, udata->comm, &udata->req[0]);
    if (retval != MPI_SUCCESS) MPI_Abort(udata->comm, 1);

    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, nvar),
      [=] DEVICE_LAMBDA (int var) {
      Esend[IDX(nvar, 0, var)] = Ydata[IDX(nvar, N - 1, var)];
    });

    retval = MPI_Isend(udata->Esend, nvar, MPI_DOUBLE, ipE,
                       0, udata->comm, &udata->req[1]);
    if (retval != MPI_SUCCESS) MPI_Abort(udata->comm, 1);
  }
  else if (c < 0.0)
  {
    /* Left moving flow uses forward difference.
       Send from east to west (first processor sends to last) */

    retval = MPI_Irecv(udata->Erecv, nvar, MPI_DOUBLE, ipE,
                       MPI_ANY_TAG, udata->comm, &udata->req[0]);
    if (retval != MPI_SUCCESS) MPI_Abort(udata->comm, 1);

    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, nvar),
      [=] DEVICE_LAMBDA (int var) {
      Wsend[IDX(nvar, 0, var)] = Ydata[IDX(nvar, 0, var)];
    });

    retval = MPI_Isend(udata->Wsend, nvar, MPI_DOUBLE, ipW,
                       0, udata->comm, &udata->req[1]);
    if (retval != MPI_SUCCESS) MPI_Abort(udata->comm, 1);
  }

  return(0);
}


/* Completes the exchange of the neighbor information */
int ExchangeAllEnd(UserData* udata)
{
  int ierr;
  MPI_Status stat[2];

  ierr = MPI_Waitall(2, udata->req, stat);
  if (ierr != MPI_SUCCESS)
  {
    fprintf(stderr, "\nERROR: error in MPI_Waitall = %d\n", ierr);
    return -1;
  }

  udata->t_comm.stop();

  return(0);
}


/* Get the vector data array pointer for the device
   if using the GPU, or host if not. */
double* GetVecData(N_Vector y)
{
#ifdef USE_GPU
  if (N_VGetVectorID(y) == SUNDIALS_NVEC_MPIPLUSX)
    return N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(y));
  else
    return N_VGetDeviceArrayPointer(y);
#else
  if (N_VGetVectorID(y) == SUNDIALS_NVEC_MPIPLUSX)
    return N_VGetArrayPointer(N_VGetLocalVector_MPIPlusX(y));
  else
    return N_VGetArrayPointer(y);
#endif
}


/* Turn on fused vector ops for y */
int EnableFusedVectorOps(N_Vector y)
{
  int retval = 0;

  retval = N_VEnableFusedOps_MPIPlusX(y, 1);
  if (check_retval(&retval, "N_VEnableFusedOps_MPIPlusX", 1)) return(-1);

#if defined(USE_CUDA_NVEC) || defined(USE_CUDAUVM_NVEC)
  retval = N_VEnableFusedOps_Cuda(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_Cuda", 1)) return(-1);
#elif defined(USE_HIP_NVEC)
  retval = N_VEnableFusedOps_Hip(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_Hip", 1)) return(-1);
#elif defined(USE_RAJA_NVEC)
  retval = N_VEnableFusedOps_Raja(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_Raja", 1)) return(-1);
#elif defined(USE_OMPDEV_NVEC)
  retval = N_VEnableFusedOps_OpenMPDEV(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_OpenMPDEV", 1)) return(-1);
#else
  retval = N_VEnableFusedOps_Serial(N_VGetLocalVector_MPIPlusX(y), 1);
  if (check_retval(&retval, "N_VEnableFusedOps_Serial", 1)) return(-1);
#endif

  return(0);
}


/* Parses the CLI arguments and sets up the problem */
int SetupProblem(int argc, char *argv[], UserData* udata, UserOptions* uopt)
{
  /* time this function */
  udata->t_setup.start();

  /* local variables */
  char fname[MXSTR]; /* output file name       */

  /* MPI variables */
  udata->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(udata->comm, &udata->myid);
  MPI_Comm_size(udata->comm, &udata->nprocs);

  /* default problem parameters */
  const int      nvar = 3;     /* number of solution fields               */
  const long int NX   = 100;   /* global spatial mesh size (NX intervals) */
  const double   xmax = 1.0;   /* maximum x value          */
  const double   A    = 1.0;   /* problem parameters                      */
  const double   B    = 3.5;
  const double   k    = 1.0;
  const double   c    = 0.01;

  /* set default user data values */
  udata->nvar = nvar;
  udata->nx   = NX;
  udata->nxl  = (int) NX / udata->nprocs;
  udata->NEQ  = (int) nvar * udata->nxl;
  udata->xmax = xmax;
  udata->dx   = xmax / (double) NX;
  udata->A    = A;
  udata->B    = B;
  udata->k1   = k;
  udata->k2   = k;
  udata->k3   = k;
  udata->k4   = k;
  udata->k5   = 1.0/5.0e-6;
  udata->k6   = 1.0/5.0e-6;
  udata->c    = c;
  udata->uopt = uopt;
  udata->TFID = NULL;
  udata->UFID = NULL;
  udata->VFID = NULL;
  udata->WFID = NULL;

  /* set default integrator options */
  uopt->order     = 3;            /* method order             */
  uopt->expl      = 0;            /* imex or explicit         */
  uopt->t0        = 0.0;          /* initial time             */
  uopt->tf        = 10.0;         /* final time               */
  uopt->rtol      = 1.0e-6;       /* relative tolerance       */
  uopt->atol      = 1.0e-9;       /* absolute tolerance       */
  uopt->global    = 0;            /* use global NLS           */
  uopt->fused     = 0;            /* use fused vector ops     */
  uopt->monitor   = 0;            /* print solution to screen */
  uopt->printtime = 0;            /* print timing             */
  uopt->nout      = 40;           /* number of output times   */
  uopt->outputdir = (char *) "."; /* output directory         */

  /* check for input args */
  if (argc > 1)
  {
    /* loop over input args and get value */
    for (int i = 1; i < argc; i++)
    {
      if (strcmp(argv[i],"--help") == 0)
      {
        InputError(argv[0]);
        return(-1);
      }
      else if (strcmp(argv[i],"--monitor") == 0)
      {
        uopt->monitor = 1;
      }
      else if (strcmp(argv[i],"--printtime") == 0)
      {
        uopt->printtime = 1;
      }
      else if (strcmp(argv[i],"--nout") == 0)
      {
        uopt->nout = atoi(argv[i+1]);
        i++;
      }
      else if (strcmp(argv[i],"--output-dir") == 0)
      {
        uopt->outputdir = argv[i+1];
        if (strlen(argv[i+1]) > MXSTR)
        {
          if (udata->myid == 0)
            fprintf(stderr, "ERROR: output directory string is too long\n");
          return(-1);
        }
        i++;
      }
      else if (strcmp(argv[i],"--nx") == 0)
      {
        udata->nx = atoi(argv[i+1]);
        i++;
      }
      else if (strcmp(argv[i],"--xmax") == 0)
      {
        udata->xmax = strtod(argv[i+1], NULL);
        i++;
      }
      else if (strcmp(argv[i],"--A") == 0)
      {
        udata->A = strtod(argv[i+1], NULL);
        i++;
      }
      else if (strcmp(argv[i],"--B") == 0)
      {
        udata->B = strtod(argv[i+1], NULL);
        i++;
      }
      else if (strcmp(argv[i],"--k") == 0)
      {
        udata->k1 = strtod(argv[i+1], NULL);
        udata->k2 = strtod(argv[i+1], NULL);
        udata->k3 = strtod(argv[i+1], NULL);
        udata->k4 = strtod(argv[i+1], NULL);
        i += 4;
      }
      else if (strcmp(argv[i],"--c") == 0)
      {
        udata->c = strtod(argv[i+1], NULL);
        i++;
      }
      else if (strcmp(argv[i],"--order") == 0)
      {
        uopt->order = atoi(argv[i+1]);
        i++;
      }
      else if (strcmp(argv[i],"--explicit") == 0)
      {
        uopt->expl = 1;
      }
      else if (strcmp(argv[i],"--global-nls") == 0)
      {
        uopt->global = 1;
      }
      else if (strcmp(argv[i],"--fused") == 0)
      {
        uopt->fused = 1;
      }
      else if (strcmp(argv[i],"--tf") == 0)
      {
        uopt->tf = strtod(argv[i+1], NULL);
        i++;
      }
      else if (strcmp(argv[i],"--rtol") == 0)
      {
        uopt->rtol = strtod(argv[i+1], NULL);
        i++;
      }
      else if (strcmp(argv[i],"--atol") == 0)
      {
        uopt->atol = strtod(argv[i+1], NULL);
        i++;
      }
      else
      {
        InputError(argv[0]);
        return(-1);
      }
    }
  }

  /* Setup the parallel decomposition */
  if (udata->nx % udata->nprocs)
  {
    fprintf(stderr, "\nERROR: The mesh size (nx = %li) must be divisible by the number of processors (%li)\n",
           (long int) udata->nx, (long int) udata->nprocs);
    return(-1);
  }
  udata->nxl = (int) udata->nx / udata->nprocs;
  udata->NEQ = udata->nvar * udata->nxl;
  udata->dx  = udata->xmax / (double) udata->nx;  /* nx is number of intervals */

  /* Create the MPI exchange buffers */
#if defined(USE_CUDA)
  GPU_SAFE_CALL( cudaHostAlloc(&(udata->Wsend), udata->nvar*sizeof(double), cudaHostAllocDefault) );
  GPU_SAFE_CALL( cudaHostAlloc(&(udata->Wrecv), udata->nvar*sizeof(double), cudaHostAllocDefault) );
  GPU_SAFE_CALL( cudaHostAlloc(&(udata->Esend), udata->nvar*sizeof(double), cudaHostAllocDefault) );
  GPU_SAFE_CALL( cudaHostAlloc(&(udata->Erecv), udata->nvar*sizeof(double), cudaHostAllocDefault) );
#elif defined(USE_HIP)
  GPU_SAFE_CALL( hipHostMalloc(&(udata->Wsend), udata->nvar*sizeof(double), hipHostMallocDefault) );
  GPU_SAFE_CALL( hipHostMalloc(&(udata->Wrecv), udata->nvar*sizeof(double), hipHostMallocDefault) );
  GPU_SAFE_CALL( hipHostMalloc(&(udata->Esend), udata->nvar*sizeof(double), hipHostMallocDefault) );
  GPU_SAFE_CALL( hipHostMalloc(&(udata->Erecv), udata->nvar*sizeof(double), hipHostMallocDefault) );
#elif defined(USE_OMPDEV_NVEC)
  int dev = omp_get_default_device();

  udata->Wsend = (double *) omp_target_alloc(udata->nvar*sizeof(double), dev);
  if (check_retval((void *)udata->Wsend, "omp_target_alloc", 0)) return 1;

  udata->Wrecv = (double *) omp_target_alloc(udata->nvar*sizeof(double), dev);
  if (check_retval((void *)udata->Wrecv, "omp_target_alloc", 0)) return 1;

  udata->Esend = (double *) omp_target_alloc(udata->nvar*sizeof(double), dev);
  if (check_retval((void *)udata->Esend, "omp_target_alloc", 0)) return 1;

  udata->Erecv = (double *) omp_target_alloc(udata->nvar*sizeof(double), dev);
  if (check_retval((void *)udata->Erecv, "omp_target_alloc", 0)) return 1;
#else
  udata->Wsend = (double *) calloc(udata->nvar, sizeof(double));
  if (check_retval((void *)udata->Wsend, "calloc", 0)) return 1;

  udata->Wrecv = (double *) calloc(udata->nvar, sizeof(double));
  if (check_retval((void *)udata->Wrecv, "calloc", 0)) return 1;

  udata->Esend = (double *) calloc(udata->nvar, sizeof(double));
  if (check_retval((void *)udata->Esend, "calloc", 0)) return 1;

  udata->Erecv = (double *) calloc(udata->nvar, sizeof(double));
  if (check_retval((void *)udata->Erecv, "calloc", 0)) return 1;
#endif

  /* Create the solution masks */
  udata->umask = N_VMake_MPIPlusX(udata->comm, LocalNvector(udata->NEQ));
  udata->vmask = N_VClone(udata->umask);
  udata->wmask = N_VClone(udata->umask);

  N_VConst(0.0, udata->umask);
  N_VConst(0.0, udata->vmask);
  N_VConst(0.0, udata->wmask);

  double* umask = GetVecData(udata->umask);
  double* vmask = GetVecData(udata->vmask);
  double* wmask = GetVecData(udata->wmask);

  RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, udata->nxl),
    [=] DEVICE_LAMBDA (int i) {
    umask[IDX(nvar, i, 0)] = 1.0;
    vmask[IDX(nvar, i, 1)] = 1.0;
    wmask[IDX(nvar, i, 2)] = 1.0;
  });

  /* Setup the linear system data structures */
  udata->nnlfi = 0;

  /* Open output files for results */
  if (uopt->nout > 0)
  {
      if (udata->myid == 0)
      {
        sprintf(fname, "%s/t.%06d.txt", uopt->outputdir, udata->myid);
        udata->TFID = fopen(fname, "w");
      }

      sprintf(fname, "%s/u.%06d.txt", uopt->outputdir, udata->myid);
      udata->UFID = fopen(fname, "w");

      sprintf(fname, "%s/v.%06d.txt", uopt->outputdir, udata->myid);
      udata->VFID = fopen(fname, "w");

      sprintf(fname, "%s/w.%06d.txt", uopt->outputdir, udata->myid);
      udata->WFID = fopen(fname, "w");
  }

  /* Print problem setup */
  if (udata->myid == 0)
  {
    printf("\n\t\t1D Advection-Reaction Test Problem\n\n");
    printf("Using the %s NVECTOR\n", NVECTOR_ID_STRING);
    printf("Number of Processors = %li\n", (long int) udata->nprocs);
    printf("Mesh Info:\n");
    printf("  NX = %li, NXL = %d, dx = %.6f, xmax = %.6f\n",
           udata->nx, udata->nxl, udata->dx, udata->xmax);
    printf("Problem Parameters:\n");
    printf("  A = %g\n", udata->A);
    printf("  B = %g\n", udata->B);
    printf("  k = %g\n", udata->k1);
    printf("  c = %g\n", udata->c);
    printf("Integrator Options:\n");
    printf("  order            = %d\n", uopt->order);
    printf("  method           = %s\n", uopt->expl ? "explicit" : "imex");
    printf("  fused vector ops = %d\n", uopt->fused);
    printf("  t0               = %g\n", uopt->t0);
    printf("  tf               = %g\n", uopt->tf);
    printf("  reltol           = %.1e\n", uopt->rtol);
    printf("  abstol           = %.1e\n", uopt->atol);
    printf("  nout             = %d\n", uopt->nout);
    if (uopt->expl)
      printf("  nonlinear solver = none\n");
    else
      printf("  nonlinear solver = %s\n", uopt->global ? "global" : "task local");
    printf("Output directory: %s\n", uopt->outputdir);
  }

  /* stop timing */
  udata->t_setup.stop();

  /* return success */
  return(0);
}


UserData::~UserData()
{
  /* free solution masks */
  N_VDestroy(N_VGetLocalVector_MPIPlusX(umask));
  N_VDestroy(umask);
  N_VDestroy(vmask);
  N_VDestroy(wmask);

  /* free exchange buffers */
#ifdef USE_CUDA_OR_HIP
  GPU_SAFE_CALL( GPU_PREFIX(FreeHost)(Wsend) );
  GPU_SAFE_CALL( GPU_PREFIX(FreeHost)(Wrecv) );
  GPU_SAFE_CALL( GPU_PREFIX(FreeHost)(Esend) );
  GPU_SAFE_CALL( GPU_PREFIX(FreeHost)(Erecv) );
#elif USE_OMPDEV_NVEC
  omp_target_free(Wsend);
  omp_target_free(Wrecv);
  omp_target_free(Esend);
  omp_target_free(Erecv);
#else
  delete Wsend;
  delete Wrecv;
  delete Esend;
  delete Erecv;
#endif

  /* close output streams */
  if (uopt->nout > 0)
  {
    fclose(UFID);
    fclose(VFID);
    fclose(WFID);
    if (myid == 0) fclose(TFID);
  }
}


void InputError(char *name)
{
  int myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0)
  {
    fprintf(stderr, "\nERROR: Invalid command line input\n");
    fprintf(stderr, "\nCommand line options for %s\n",name);
    fprintf(stderr, "  --help           prints this message\n");
    fprintf(stderr, "  --printtime      print timing information\n");
    fprintf(stderr, "  --monitor        print solution information to screen (slower)\n");
    fprintf(stderr, "  --output-dir     the directory where all output files will be written\n");
    fprintf(stderr, "  --nout <int>     number of output times\n");
    fprintf(stderr, "  --explicit       use an explicit method instead of IMEX\n");
    fprintf(stderr, "  --global-nls     use a global newton nonlinear solver instead of task-local (for IMEX only)\n");
    fprintf(stderr, "  --order <int>    the method order to use\n");
    fprintf(stderr, "  --nx <int>       number of mesh points\n");
    fprintf(stderr, "  --xmax <double>  maximum value of x (size of domain)\n");
    fprintf(stderr, "  --tf <double>    final time\n");
    fprintf(stderr, "  --A <double>     A parameter value\n");
    fprintf(stderr, "  --B <double>     B parameter value\n");
    fprintf(stderr, "  --k <double>     reaction rate\n");
    fprintf(stderr, "  --c <double>     advection speed\n");
    fprintf(stderr, "  --rtol <double>  relative tolerance\n");
    fprintf(stderr, "  --atol <double>  absolute tolerance\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
}


/* --------------------------------------------------------------
 * Function to check return values:
 *
 * opt == 0  means the function allocates memory and returns a
 *           pointer so check if a NULL pointer was returned
 * opt == 1  means the function returns an integer where a
 *           value < 0 indicates an error occured
 * --------------------------------------------------------------*/
int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *errvalue, myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (opt == 0 && returnvalue == NULL)
  {
    /* A NULL pointer was returned - no memory allocated */
    if (myid == 0)
      fprintf(stderr, "\nERROR: %s() failed - returned NULL pointer\n\n",
              funcname);
    return(1);
  }
  else if (opt == 1)
  {
    errvalue = (int *) returnvalue;

    /* A value < 0 was returned - function failed */
    if (*errvalue < 0)
    {
      if (myid == 0)
        fprintf(stderr, "\nERROR: %s() returned %d\n\n", funcname, *errvalue);
      return(1);
    }
  }

  /* return success */
  return(0);
}

/* function to sync the host and device */
void sync_device()
{
#ifdef USE_CUDA_OR_HIP
  GPU_PREFIX(DeviceSynchronize)();
#endif
}


#ifdef USE_CUDA_OR_HIP
void gpuAssert(GPU_PREFIX(Error_t) code, const char *file, int line, int abort)
{
   if (code != GPU_PREFIX(Success))
   {
      fprintf(stderr, "GPU ERROR: %s %s %d\n", GPU_PREFIX(GetErrorString)(code), file, line);
      if (abort) MPI_Abort(MPI_COMM_WORLD, -1);
   }
}
#endif

/*---- end of file ----*/
