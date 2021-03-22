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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "arkode/arkode_arkstep.h"             /* ARKStep                   */
#include "arkode/arkode_erkstep.h"             /* ERKStep                   */
#include "nvector/nvector_mpiplusx.h"          /* MPI+X N_Vector            */
#include "nvector/nvector_serial.h"            /* serial N_Vector           */
#include "sunmatrix/sunmatrix_dense.h"         /* dense SUNMatrix           */
#include "sunlinsol/sunlinsol_dense.h"         /* dense SUNLinearSolver     */
#include "sunlinsol/sunlinsol_spgmr.h"         /* GMRES SUNLinearSolver     */
#include "sunnonlinsol/sunnonlinsol_newton.h"  /* Newton SUNNonlinearSolver */


/* Maximum size of output directory string */
#define MXSTR 2048

/* Accessor macro:
   n = number of state variables
   i = mesh node index
   c = component */
#define IDX(n,i,c) ((n)*(i)+(c))


/*
 * User options structure
 */

typedef struct
{
  double t0;       /* initial time                 */
  double tf;       /* final time                   */
  double rtol;     /* relative tolerance           */
  double atol;     /* absolute tolerance           */
  int    order;    /* method order                 */
  int    explicit; /* imex method or explicit      */
  int    global;   /* use global nonlinear solve   */
  int    fused;    /* use fused vector ops         */
  int    nout;     /* number of outputs            */
  int    monitor;  /* print solution to screen     */
  int    printtime;/* print timing information     */
  FILE*  TFID;     /* time output file pointer     */
  FILE*  UFID;     /* solution output file pointer */
  FILE*  VFID;
  FILE*  WFID;
  char*  outputdir;
} *UserOptions;


/*
 * User data structure
 */

typedef struct
{
  /* MPI data */
  MPI_Comm    comm;
  int         myid;
  int         nprocs;
  MPI_Request reqS;
  MPI_Request reqR;
  double*     Wsend;
  double*     Esend;
  double*     Wrecv;
  double*     Erecv;

  /* data structures for per mesh node linear system */
  N_Vector        b_node;
  SUNMatrix       jac_node;
  SUNLinearSolver ls_node;
  long int        nnlfi;

  /* data structures for task-local preconditioner */
  SUNMatrix       pre;
  SUNLinearSolver prels;

  /* solution masks */
  N_Vector umask;
  N_Vector vmask;
  N_Vector wmask;

  /* problem paramaters */
  long long nvar; /* number of species            */
  long long nx;   /* number of intervals globally */
  long long nxl;  /* number of intervals locally  */
  long long NEQ;  /* number of equations locally  */
  double    dx;   /* mesh spacing                 */
  double    xmax; /* maximum x value              */
  double    A;    /* concentration of species A   */
  double    B;    /* w source rate                */
  double    k1;   /* reaction rates               */
  double    k2;
  double    k3;
  double    k4;
  double    k5;
  double    k6;
  double    c;    /* advection coefficient        */

  /* integrator options */
  UserOptions uopt;
} *UserData;


/*
 * Definitions for a custom task local SUNNonlinearSolver
 */

typedef struct
{
  int                myid;
  int                nprocs;
  long int           ncnf;
  MPI_Comm           comm;
  SUNNonlinearSolver local_nls;
} *TaskLocalNewton_Content;

/* Content accessor macors */
#define GET_NLS_CONTENT(NLS) ( (TaskLocalNewton_Content)(NLS->content) )
#define LOCAL_NLS(NLS)       ( GET_NLS_CONTENT(NLS)->local_nls )

/* SUNNonlinearSolver constructor */
SUNNonlinearSolver TaskLocalNewton(N_Vector y, FILE* DFID);


/*
 * RHS functions provided to the integrator
 */

static int Advection(double t, N_Vector y, N_Vector ydot, void *user_data);
static int Reaction(double t, N_Vector y, N_Vector ydot, void *user_data);
static int AdvectionReaction(double t, N_Vector y, N_Vector ydot,
                             void *user_data);

/*
 * Preconditioner functions (used only when using the global nonlinear solver)
 */

static int PSetup(double t, N_Vector y, N_Vector f, booleantype jok,
                  booleantype *jcurPtr, double gamma, void *user_data);
static int PSolve(double t, N_Vector y, N_Vector f, N_Vector r,
                  N_Vector z, double gamma, double delta, int lr,
                  void *user_data);


/*
 * Helper functions
 */

/* function that does ARKStep setup and evolves the solution */
static int EvolveProblemIMEX(N_Vector y, UserData udata, UserOptions uopt);

/* function that does ERKStep setup and evolves the solution */
static int EvolveProblemExplicit(N_Vector y, UserData udata, UserOptions uopt);

/* function to set initial condition */
static int SetIC(N_Vector y, UserData udata);

/* functions to exchange neighbor data */
static int ExchangeBCOnly(N_Vector y, UserData udata);
static int ExchangeAllStart(N_Vector y, UserData udata);
static int ExchangeAllEnd(UserData udata);

/* functions for processing command line args */
static int SetupProblem(int argc, char *argv[],
                        UserData udata, UserOptions uopt);
static void InputError(char *name);

/* function to write solution to disk */
static int WriteOutput(double t, N_Vector y, UserData udata, UserOptions uopt);

/* function to free the problem data */
static void FreeProblem(UserData udata, UserOptions uopt);

/* function to check sundials return values */
static int check_retval(void *returnvalue, const char *funcname, int opt);


/* Main Program */
int main(int argc, char *argv[])
{
  /* general problem variables */
  N_Vector     y = NULL;      /* empty solution vector        */
  UserData     udata;         /* user data                    */
  UserOptions  uopt;          /* user options                 */
  int          retval;        /* reusable error-checking flag */
  long long    i;             /* loop counter                 */
  FILE*        MFID;          /* mesh output file pointer     */
  char         fname[MXSTR];
  MPI_Comm     comm;
  double       starttime, endtime;

  /* Initialize MPI */
  comm = MPI_COMM_WORLD;
  retval = MPI_Init(&argc, &argv);
  if (check_retval(&retval, "MPI_Init", 1)) MPI_Abort(comm, 1);

  /* Start timing */
  starttime = MPI_Wtime();

  /* Allocate user data structure */
  udata = (UserData) malloc(sizeof(*udata));
  if (check_retval((void *) udata, "malloc", 0)) MPI_Abort(comm, 1);

  /* Allocate user options structure */
  uopt = (UserOptions) malloc(sizeof(*uopt));
  if (check_retval((void *) uopt, "malloc", 0)) MPI_Abort(comm, 1);

  /* Process input args and setup the problem */
  retval = SetupProblem(argc, argv, udata, uopt);
  if (check_retval(&retval, "SetupProblem", 1)) MPI_Abort(comm, 1);

  /* Create solution vector */
  y = N_VMake_MPIPlusX(udata->comm, N_VNew_Serial(udata->NEQ));
  if (check_retval((void *) y, "N_VMake_MPIPlusX", 0)) MPI_Abort(comm, 1);

  /* Enabled fused vector ops */
  if (uopt->fused)
  {
    retval = N_VEnableFusedOps_Serial(N_VGetLocalVector_MPIPlusX(y),
                                      uopt->fused);
    if (check_retval(&retval, "N_VEnableFusedOps_Serial", 1))
      MPI_Abort(comm, 1);

    retval = N_VEnableFusedOps_MPIManyVector(y, uopt->fused);
    if (check_retval(&retval, "N_VEnableFusedOps_MPIManyVector", 1))
      MPI_Abort(comm, 1);
  }

  /* Set the initial condition */
  retval = SetIC(y, udata);
  if (check_retval(&retval, "SetIC", 1)) MPI_Abort(comm, 1);

  /* Output spatial mesh to disk (add extra point for periodic BC) */
  if (udata->myid == 0 && uopt->nout > 0)
  {
    sprintf(fname, "%s/mesh.txt", uopt->outputdir);
    MFID = fopen(fname,"w");
    for (i = 0; i < udata->nx + 1; i++) fprintf(MFID,"  %.16e\n", udata->dx*i);
    fclose(MFID);
  }

  /* Integrate in time */
  if (uopt->explicit) retval = EvolveProblemExplicit(y, udata, uopt);
  else                retval = EvolveProblemIMEX(y, udata, uopt);

  if (check_retval(&retval, "Evolve", 1)) MPI_Abort(comm, 1);

  /* End timing */
  endtime = MPI_Wtime();
  if (udata->myid == 0 && uopt->printtime)
  {
    printf("\nTotal wall clock time: %.4f seconds\n", endtime-starttime);
  }

  /* Clean up */
  FreeProblem(udata, uopt);
  free(udata);
  free(uopt);
  N_VDestroy(N_VGetLocalVector_MPIPlusX(y));
  N_VDestroy(y);

  /* return success */
  MPI_Finalize();
  return 0;
}

/* Setup ARKODE and evolve problem in time with IMEX method*/
int EvolveProblemIMEX(N_Vector y, UserData udata, UserOptions uopt)
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
  long int nli, npre, npsol;  /* linear solver stats          */
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
  if (check_retval(&retval, "ARKStepSetUserData", 1)) return 1;

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
    LS  = SUNLinSol_SPGMR(y, PREC_LEFT, 0);
    if (check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return 1;

    /* Attach linear solver */
    retval = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
    if (check_retval(&retval, "ARKStepSetLinearSolver", 1)) return 1;

    /* Attach preconditioner */
    retval = ARKStepSetPreconditioner(arkode_mem, PSetup, PSolve);
    if (check_retval(&retval, "ARKStepSetPreconditioner", 1)) return 1;
  }
  else
  {
    /* The custom task-local nonlinear solver handles the linear solve
      as well, so we do not need a SUNLinearSolver */
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
    retval = ARKStepGetNumPrecEvals(arkode_mem, &npre);
    check_retval(&retval, "ARKStepGetNumPrecEvals", 1);
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
      printf("   Total number of preconditioner setups = %li\n", npre);
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
int EvolveProblemExplicit(N_Vector y, UserData udata, UserOptions uopt)
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
  return 0;
}


/* Write time and solution to disk */
int WriteOutput(double t, N_Vector y, UserData udata, UserOptions uopt)
{
  long long i;
  long long nvar = udata->nvar;
  double    u, v, w;
  double*   data = NULL;

  /* output current solution norm to screen */
  if(uopt->monitor)
  {
    u = N_VWL2Norm(y,udata->umask);
    u = sqrt(u*u/udata->nx);
    v = N_VWL2Norm(y,udata->vmask);
    v = sqrt(v*v/udata->nx);
    w = N_VWL2Norm(y,udata->wmask);
    w = sqrt(w*w/udata->nx);
    if (udata->myid == 0)
      printf("     %10.6f   %10.6f   %10.6f   %10.6f\n", t, u, v, w);
  }

  if (uopt->nout > 0)
  {
    /* get left end point for output */
    ExchangeBCOnly(y, udata);

    /* get vector data array */
    data = N_VGetArrayPointer(y);
    if (check_retval((void *) data, "N_VGetArrayPointer", 0)) return 1;

    /* output the times to disk */
    if (udata->myid == 0 && uopt->TFID)
      fprintf(uopt->TFID," %.16e\n", t);

    /* output results to disk */
    for (i = 0; i < udata->nxl; i++)
    {
      fprintf(uopt->UFID," %.16e", data[IDX(nvar, i, 0)]);
      fprintf(uopt->VFID," %.16e", data[IDX(nvar, i, 1)]);
      fprintf(uopt->WFID," %.16e", data[IDX(nvar, i, 2)]);
    }

    /* we have one extra output because of the periodic BCs */
    if (udata->myid == (udata->nprocs - 1))
    {
      fprintf(uopt->UFID," %.16e\n", udata->Erecv[IDX(nvar, 0, 0)]);
      fprintf(uopt->VFID," %.16e\n", udata->Erecv[IDX(nvar, 0, 1)]);
      fprintf(uopt->WFID," %.16e\n", udata->Erecv[IDX(nvar, 0, 2)]);
    }
    else
    {
      fprintf(uopt->UFID,"\n");
      fprintf(uopt->VFID,"\n");
      fprintf(uopt->WFID,"\n");
    }
  }

  return 0;
}


/* Initial Condition Functions */
int SetIC(N_Vector y, UserData udata)
{
  /* Variable shortcuts */
  long long nvar = udata->nvar;
  long long N    = udata->nxl;
  double    dx   = udata->dx;
  double    A    = udata->A;
  double    B    = udata->B;
  double    k1   = udata->k1;
  double    k2   = udata->k2;
  double    k3   = udata->k3;
  double    k4   = udata->k4;
  int       myid = udata->myid;

  /* Local variables */
  double*   data = NULL;
  double    x, us, vs, ws, p;
  long long i;

  /* Gaussian distribution defaults */
  double mu    = udata->xmax / 2.0;
  double sigma = udata->xmax / 4.0;
  double alpha = 0.1;

  /* Access data array from NVector y */
  data = N_VGetArrayPointer(y);

  /* Steady state solution */
  us = k1 * A / k4;
  vs = k2 * k4 * B / (k1 * k3 * A);
  ws = 3.0;

  /* Gaussian perturbation of the steady state solution */
  for (i = 0; i < N; i++)
  {
    x = (myid * N + i) * dx;
    p = alpha * exp( -((x - mu) * (x - mu)) / (2.0 * sigma * sigma) );
    data[IDX(nvar, i, 0)] = us + p;
    data[IDX(nvar, i, 1)] = vs + p;
    data[IDX(nvar, i, 2)] = ws + p;
  }

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
  UserData udata = (UserData) user_data;

  /* set variable shortcuts */
  long long nvar = udata->nvar;
  long long N    = udata->nxl;
  double    dx   = udata->dx;
  double    c    = udata->c;

  /* local variables */
  double*    Ydata  = NULL;
  double*    dYdata = NULL;
  double     tmp;
  long long  i, var;
  int        retval;

  /* set output to zero */
  N_VConst(0.0, ydot);

  /* access data arrays */
  Ydata = N_VGetArrayPointer(y);
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return(-1);

  dYdata = N_VGetArrayPointer(ydot);
  if (check_retval((void *)dYdata, "N_VGetArrayPointer", 0)) return(-1);

  /* begin exchanging boundary information */
  retval = ExchangeAllStart(y, udata);
  if (check_retval(&retval, "ExchangeAllStart", 1)) return(-1);

  /* iterate over domain interior, computing advection */
  tmp = -c / dx;

  if (c > 0.0)
  {
    /* right moving flow */
    for (i = 1; i < N; i++)
    {
      for (var = 0; var < nvar; var++)
      {
        dYdata[IDX(nvar, i, var)] = tmp * (Ydata[IDX(nvar, i, var)]
                                           - Ydata[IDX(nvar, i-1, var)]);
      }
    }
  }
  else if (c < 0.0)
  {
    /* left moving flow */
    for (i = 0; i < N - 1; i++)
    {
      for (var = 0; var < nvar; var++)
      {
        dYdata[IDX(nvar, i, var)] = tmp * (Ydata[IDX(nvar, i + 1, var)]
                                           - Ydata[IDX(nvar, i, var)]);
      }
    }
  }

  /* finish exchanging boundary information */
  retval = ExchangeAllEnd(udata);
  if (check_retval(&retval, "ExchangeAllEnd", 1)) return(-1);

  /* compute advection at local boundaries */
  if (c > 0.0)
  {
    /* right moving flow (left boundary) */
    for (var = 0; var < nvar; var++)
    {
      dYdata[IDX(nvar, 0, var)] = tmp * (Ydata[IDX(nvar, 0, var)]
                                         - udata->Wrecv[IDX(nvar, 0, var)]);
    }
  }
  else if (c < 0.0)
  {
    /* left moving flow (right boundary) */
    for (var = 0; var < nvar; var++)
    {
      dYdata[IDX(nvar, N - 1, var)] = tmp * (udata->Erecv[IDX(nvar, 0, var)]
                                             - Ydata[IDX(nvar, N - 1, var)]);
    }
  }

  /* return success */
  return(0);
}


/* Compute the reaction term. */
int Reaction(double t, N_Vector y, N_Vector ydot, void* user_data)
{
  /* access problem data */
  UserData udata = (UserData) user_data;

  /* set variable shortcuts */
  long long nvar = udata->nvar;
  long long N    = udata->nxl;
  double    A    = udata->A;
  double    B    = udata->B;
  double    k1   = udata->k1;
  double    k2   = udata->k2;
  double    k3   = udata->k3;
  double    k4   = udata->k4;
  double    k5   = udata->k5;
  double    k6   = udata->k6;

  /* local variables */
  double*   Ydata  = NULL;
  double*   dYdata = NULL;
  double    u, v, w;
  long long i;

  /* access data arrays */
  Ydata = N_VGetArrayPointer(y);
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return(-1);

  dYdata = N_VGetArrayPointer(ydot);
  if (check_retval((void *)dYdata, "N_VGetArrayPointer", 0)) return(-1);

  /* iterate over domain, computing reactions */
  if (udata->uopt->explicit)
  {
    /* when integrating explicitly, we add to ydot as we expect it
       to hold the advection term already */
    for (i = 0; i < N; i++)
    {
      u = Ydata[IDX(nvar, i ,0)];
      v = Ydata[IDX(nvar, i, 1)];
      w = Ydata[IDX(nvar, i, 2)];
      dYdata[IDX(nvar, i, 0)] += k1 * A - k2 * w * u + k3 * u * u * v - k4 * u;
      dYdata[IDX(nvar, i, 1)] += k2 * w * u - k3 * u * u * v;
      dYdata[IDX(nvar, i, 2)] += -k2 * w * u + k5 * B - k6 * w;
    }
  }
  else
  {
    /* set output to zero */
    N_VConst(0.0, ydot);

    for (i = 0; i < N; i++)
    {
      u = Ydata[IDX(nvar, i ,0)];
      v = Ydata[IDX(nvar, i, 1)];
      w = Ydata[IDX(nvar, i, 2)];
      dYdata[IDX(nvar, i, 0)] = k1 * A - k2 * w * u + k3 * u * u * v - k4 * u;
      dYdata[IDX(nvar, i, 1)] = k2 * w * u - k3 * u * u * v;
      dYdata[IDX(nvar, i, 2)] = -k2 * w * u + k5 * B - k6 * w;
    }
  }

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
 * (Non)linear system functions
 * --------------------------------------------------------------*/


int TaskLocalNlsResidual(N_Vector ycor, N_Vector F, void* arkode_mem)
{
  /* temporary variables */
  UserData udata;
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

  udata = (UserData) user_data;

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
  UserData  udata = NULL;
  double*   Zdata = NULL;
  double*   Bdata = NULL;
  double    u, v, w;
  long long i, j;
  int       retval;

  /* shortcuts */
  long long       nvar, N;
  double          k2, k3, k4, k6;
  N_Vector        b_node;
  SUNMatrix       Jac;
  SUNLinearSolver LS;

  /* nonlinear system data */
  N_Vector z, zpred, Fi, sdata;
  double   tcur, gamma;
  void*    user_data = NULL;

  retval = ARKStepGetNonlinearSystemData(arkode_mem, &tcur, &zpred, &z, &Fi,
                                         &gamma, &sdata, &user_data);
  if (check_retval((void *)&retval, "ARKStepGetNonlinearSystemData", 1))
    return(-1);

  udata = (UserData) user_data;

  /* set shortcuts */
  nvar   = udata->nvar;
  N      = udata->nxl;
  k2     = udata->k2;
  k3     = udata->k3;
  k4     = udata->k4;
  k6     = udata->k6;
  b_node = udata->b_node;
  Jac    = udata->jac_node;
  LS     = udata->ls_node;

  /* access solution array */
  Zdata = N_VGetArrayPointer(z);
  if (check_retval((void *)Zdata, "N_VGetArrayPointer", 0)) return(-1);

  Bdata = N_VGetArrayPointer(delta);
  if (check_retval((void *)Bdata, "N_VGetArrayPointer", 0)) return(-1);

  /* solve the linear system at each mesh node */
  for (i = 0; i < N; i++)
  {
    /* fill in Jacobian entries for this mesh node */

    /* set nodal value shortcuts */
    u = Zdata[IDX(nvar, i, 0)];
    v = Zdata[IDX(nvar, i, 1)];
    w = Zdata[IDX(nvar, i, 2)];

    /* all vars wrt u */
    SM_ELEMENT_D(Jac, 0, 0) = -k2 * w + 2.0 * k3 * u * v - k4;
    SM_ELEMENT_D(Jac, 1, 0) =  k2 * w - 2.0 * k3 * u * v;
    SM_ELEMENT_D(Jac, 2, 0) = -k2 * w;

    /* all vars wrt v */
    SM_ELEMENT_D(Jac, 0, 1) =  k3 * u * u;
    SM_ELEMENT_D(Jac, 1, 1) = -k3 * u * u;
    SM_ELEMENT_D(Jac, 2, 1) =  0.0;

    /* all vars wrt w */
    SM_ELEMENT_D(Jac, 0, 2) = -k2 * u;
    SM_ELEMENT_D(Jac, 1, 2) =  k2 * u;
    SM_ELEMENT_D(Jac, 2, 2) = -k2 * u - k6;

    /* I - gamma*J */
    SUNMatScaleAddI(-gamma, Jac);

    /* grab just the portion of the vector 'b' for this mesh node */
    for (j = 0; j < nvar; ++j)
      N_VGetArrayPointer(b_node)[j] = Bdata[IDX(nvar, i, j)];

    /* setup the linear system */
    retval = SUNLinSolSetup(LS, Jac);
    if (retval != 0) break;

    /* solve the linear system */
    retval = SUNLinSolSolve(LS, Jac, b_node, b_node, 0.0);
    if (retval != 0) break;

    /* set just the portion of the vector 'b' for this mesh node */
    for (j = 0; j < nvar; ++j)
      Bdata[IDX(nvar, i, j)] = N_VGetArrayPointer(b_node)[j];
  }

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

  return(GET_NLS_CONTENT(NLS)->ncnf);
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


/* Sets P = I - gamma * J */
int PSetup(double t, N_Vector y, N_Vector ydot, booleantype jok,
           booleantype *jcurPtr, double gamma, void *user_data)
{
  /* local variables */
  UserData  udata = (UserData) user_data;
  double    *Ydata;
  double    u, v, w;
  long long i, blocki;
  int       retval = 0;

  /* shortcuts */
  long long       nvar, N;
  double          k2, k3, k4, k6;
  SUNMatrix       P;
  SUNLinearSolver LS;

  /* set shortcuts */
  nvar = udata->nvar;
  N    = udata->nxl;
  k2   = udata->k2;
  k3   = udata->k3;
  k4   = udata->k4;
  k6   = udata->k6;
  P    = udata->pre;
  LS   = udata->prels;

  /* access solution data */
  Ydata = N_VGetArrayPointer(y);
  if (check_retval((void *)Ydata, "N_VGetArrayPointer", 0)) return(-1);

  if (!jok)
  {
    /* setup the block diagonal preconditioner matrix */
    for (i = 0; i < N; i++)
    {
      /* fill in Jacobian entries for this mesh node */
      blocki = nvar*i;

      /* set nodal value shortcuts */
      u = Ydata[IDX(nvar, i, 0)];
      v = Ydata[IDX(nvar, i, 1)];
      w = Ydata[IDX(nvar, i, 2)];

      /* all vars wrt u */
      SM_ELEMENT_D(P, blocki,   blocki) = -k2 * w + 2.0 * k3 * u * v - k4;
      SM_ELEMENT_D(P, blocki+1, blocki) =  k2 * w - 2.0 * k3 * u * v;
      SM_ELEMENT_D(P, blocki+2, blocki) = -k2 * w;

      /* all vars wrt v */
      SM_ELEMENT_D(P, blocki,   blocki+1) =  k3 * u * u;
      SM_ELEMENT_D(P, blocki+1, blocki+1) = -k3 * u * u;
      SM_ELEMENT_D(P, blocki+2, blocki+1) =  0.0;

      /* all vars wrt w */
      SM_ELEMENT_D(P, blocki,   blocki+2) = -k2 * u;
      SM_ELEMENT_D(P, blocki+1, blocki+2) =  k2 * u;
      SM_ELEMENT_D(P, blocki+2, blocki+2) = -k2 * u - k6;
    }
    SUNMatScaleAddI(-gamma, P);

    /* setup the linear system Pz = r */
    retval = SUNLinSolSetup(LS, P);
    if (retval) return(retval);

    /* indicate that J is now current */
    *jcurPtr = 1;
  }
  else
  {
    *jcurPtr = 0;
  }

  return(0);
}


/* Solves Pz = r */
int PSolve(double t, N_Vector y, N_Vector ydot, N_Vector r,
           N_Vector z, double gamma, double delta, int lr,
           void *user_data)
{
  /* local variables */
  UserData  udata = (UserData) user_data;
  int       retval;

  /* shortcuts */
  SUNMatrix       P;
  SUNLinearSolver LS;
  N_Vector        z_local, r_local;

  P       = udata->pre;
  LS      = udata->prels;
  z_local = N_VGetLocalVector_MPIPlusX(z);
  r_local = N_VGetLocalVector_MPIPlusX(r);

  /* solve the task-local linear system Pz = r */
  retval = SUNLinSolSolve(LS, P, z_local, r_local, 0.0);

  return(retval);
}


/* --------------------------------------------------------------
 * Utility functions
 * --------------------------------------------------------------*/

/* Exchanges the periodic BCs only by sending the first
   mesh node to the last processor. */
int ExchangeBCOnly(N_Vector y, UserData udata)
{
  int var, ierr;
  MPI_Status stat;
  MPI_Request reqS, reqR;

  /* shortcuts */
  int nvar  = udata->nvar;
  int myid  = udata->myid;
  int first = 0;
  int last  = udata->nprocs - 1;

  /* extract the data */
  double* Ydata = N_VGetArrayPointer(y);

  /* open the East Irecv buffer */
  if (myid == last)
  {
    ierr = MPI_Irecv(udata->Erecv, nvar, MPI_DOUBLE, first,
                     MPI_ANY_TAG, udata->comm, &reqR);
  }

  /* send first mesh node to the last processor */
  if (myid == first)
  {
    for (var = 0; var < nvar; var++)
      udata->Wsend[IDX(nvar, 0, var)] = Ydata[IDX(nvar, 0, var)];
    ierr = MPI_Isend(udata->Wsend, nvar, MPI_DOUBLE,
                     last, 0, udata->comm, &reqS);
  }

  /* wait for exchange to finish */
  if (myid == last)
  {
    ierr = MPI_Wait(&reqR, &stat);
    if (ierr != MPI_SUCCESS)
    {
      fprintf(stderr, "\nERROR: error in MPI_Wait = %d\n", ierr);
      return -1;
    }
  }

  if (myid == first)
  {
    ierr = MPI_Wait(&reqS, &stat);
    if (ierr != MPI_SUCCESS)
    {
      fprintf(stderr, "\nERROR: error in MPI_Wait = %d\n", ierr);
      return -1;
    }
  }

  return 0;
}


/* Starts the exchange of the neighbor information */
int ExchangeAllStart(N_Vector y, UserData udata)
{
  int var;
  int retval;

  /* shortcuts */
  double    c     = udata->c;
  long long N     = udata->nxl;
  int       nvar  = udata->nvar;
  int       myid  = udata->myid;
  int       first = 0;
  int       last  = udata->nprocs - 1;
  int       ipW   = (myid == first) ? last : udata->myid - 1; /* periodic BC */
  int       ipE   = (myid == last) ? first : udata->myid + 1; /* periodic BC */

  /* extract the data */
  double* Ydata = N_VGetArrayPointer(y);

  if (c > 0.0)
  {
    /* Right moving flow uses backward difference.
       Send from west to east (last processor sends to first) */

    retval = MPI_Irecv(udata->Wrecv, nvar, MPI_DOUBLE, ipW,
                       MPI_ANY_TAG, udata->comm, &udata->reqR);
    if (retval != MPI_SUCCESS) MPI_Abort(udata->comm, 1);

    for (var = 0; var < nvar; var++)
      udata->Esend[IDX(nvar, 0, var)] = Ydata[IDX(nvar, N - 1, var)];

    retval = MPI_Isend(udata->Esend, nvar, MPI_DOUBLE, ipE,
                       0, udata->comm, &udata->reqS);
    if (retval != MPI_SUCCESS) MPI_Abort(udata->comm, 1);
  }
  else if (c < 0.0)
  {
    /* Left moving flow uses forward difference.
       Send from east to west (first processor sends to last) */

    retval = MPI_Irecv(udata->Erecv, nvar, MPI_DOUBLE, ipE,
                       MPI_ANY_TAG, udata->comm, &udata->reqR);
    if (retval != MPI_SUCCESS) MPI_Abort(udata->comm, 1);

    for (var = 0; var < nvar; var++)
      udata->Wsend[IDX(nvar, 0, var)] = Ydata[IDX(nvar, 0, var)];

    retval = MPI_Isend(udata->Wsend, nvar, MPI_DOUBLE, ipW,
                       0, udata->comm, &udata->reqS);
    if (retval != MPI_SUCCESS) MPI_Abort(udata->comm, 1);
  }

  return 0;
}


/* Completes the exchange of the neighbor information */
int ExchangeAllEnd(UserData udata)
{
  int ierr;
  MPI_Status stat;

  /* wait for exchange to finish */
  if (udata->c < 0.0 || udata->c > 0.0)
  {
    ierr = MPI_Wait(&udata->reqR, &stat);
    if (ierr != MPI_SUCCESS)
    {
      fprintf(stderr, "\nERROR: error in MPI_Wait = %d\n", ierr);
      return -1;
    }

    ierr = MPI_Wait(&udata->reqS, &stat);
    if (ierr != MPI_SUCCESS)
    {
      fprintf(stderr, "\nERROR: error in MPI_Wait = %d\n", ierr);
      return -1;
    }
  }

  return 0;
}


int SetupProblem(int argc, char *argv[], UserData udata, UserOptions uopt)
{
  /* local variables */
  int      i, retval;
  double*  data = NULL;  /* data pointer           */
  N_Vector tmp  = NULL;  /* temporary local vector */
  char     fname[MXSTR]; /* output file name       */

  /* MPI variables */
  udata->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(udata->comm, &udata->myid);
  MPI_Comm_size(udata->comm, &udata->nprocs);

  /* default problem parameters */
  const long long nvar = 3;     /* number of solution fields               */
  const long long NX   = 100;   /* global spatial mesh size (NX intervals) */
  const double    xmax = 1.0;   /* maximum x value          */
  const double    A    = 1.0;   /* problem parameters                      */
  const double    B    = 3.5;
  const double    k    = 1.0;
  const double    c    = 0.01;

  /* set default user data values */
  udata->nvar  = nvar;
  udata->nx    = NX;
  udata->nxl   = NX / udata->nprocs;
  udata->NEQ   = nvar * udata->nxl;
  udata->xmax  = xmax;
  udata->dx    = xmax / NX;            /* NX is number of intervals */
  udata->A     = A;
  udata->B     = B;
  udata->k1    = k;
  udata->k2    = k;
  udata->k3    = k;
  udata->k4    = k;
  udata->k5    = 1.0/5.0e-6;
  udata->k6    = 1.0/5.0e-6;
  udata->c     = c;
  udata->Wsend = NULL;
  udata->Wrecv = NULL;
  udata->Esend = NULL;
  udata->Erecv = NULL;
  udata->uopt  = uopt;

  /* set default integrator options */
  uopt->order     = 3;       /* method order             */
  uopt->explicit  = 0;       /* imex or explicit         */
  uopt->t0        = 0.0;     /* initial time             */
  uopt->tf        = 10.0;    /* final time               */
  uopt->rtol      = 1.0e-6;  /* relative tolerance       */
  uopt->atol      = 1.0e-9;  /* absolute tolerance       */
  uopt->global    = 0;       /* use global NLS           */
  uopt->fused     = 0;       /* use fused vector ops     */
  uopt->monitor   = 0;       /* print solution to screen */
  uopt->printtime = 0;       /* print timing             */
  uopt->nout      = 40;      /* number of output times   */
  uopt->TFID      = NULL;
  uopt->UFID      = NULL;
  uopt->VFID      = NULL;
  uopt->WFID      = NULL;
  uopt->outputdir = ".";

  /* check for input args */
  if (argc > 1)
  {
    /* loop over input args and get value */
    for (i = 1; i < argc; i++)
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
        uopt->explicit = 1;
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
    fprintf(stderr,
            "\nERROR: The mesh size (nx = %li) must be divisible by the number of processors (%li)\n",
            (long int) udata->nx, (long int) udata->nprocs);
    return(-1);
  }
  udata->nxl = udata->nx / udata->nprocs;
  udata->NEQ = udata->nvar * udata->nxl;
  udata->dx  = udata->xmax / udata->nx;  /* nx is number of intervals */

  /* Create the MPI exchange buffers */
  udata->Wsend = (double *) calloc(udata->nvar, sizeof(double));
  if (check_retval((void *)udata->Wsend, "calloc", 0)) return 1;

  udata->Wrecv = (double *) calloc(udata->nvar, sizeof(double));
  if (check_retval((void *)udata->Wrecv, "calloc", 0)) return 1;

  udata->Esend = (double *) calloc(udata->nvar, sizeof(double));
  if (check_retval((void *)udata->Esend, "calloc", 0)) return 1;

  udata->Erecv = (double *) calloc(udata->nvar, sizeof(double));
  if (check_retval((void *)udata->Erecv, "calloc", 0)) return 1;

  /* Create the solution masks */
  udata->umask = N_VMake_MPIPlusX(udata->comm, N_VNew_Serial(udata->NEQ));
  if (uopt->fused)
  {
    retval = N_VEnableFusedOps_Serial(N_VGetLocalVector_MPIPlusX(udata->umask),
                                      uopt->fused);
    if (check_retval(&retval, "N_VEnableFusedOps_Serial", 1)) return 1;

    retval = N_VEnableFusedOps_MPIManyVector(udata->umask, uopt->fused);
    if (check_retval(&retval, "N_VEnableFusedOps_MPIManyVector", 1)) return 1;
  }

  N_VConst(0.0, udata->umask);
  data = N_VGetArrayPointer(udata->umask);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i = 0; i < udata->nxl; i++)  data[IDX(nvar, i, 0)] = 1.0;

  udata->vmask = N_VClone(udata->umask);
  N_VConst(0.0, udata->vmask);
  data = N_VGetArrayPointer(udata->vmask);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i = 0; i < udata->nxl; i++)  data[IDX(nvar, i ,1)] = 1.0;

  udata->wmask = N_VClone(udata->umask);
  N_VConst(0.0, udata->wmask);
  data = N_VGetArrayPointer(udata->wmask);
  if (check_retval((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i = 0; i < udata->nxl; i++)  data[IDX(nvar, i, 2)] = 1.0;

  /* Setup the linear system data structures */
  if (!uopt->explicit)
  {
    if (uopt->global)
    {
      /* Create MPI task-local data structures for preconditioning */
      udata->pre = SUNDenseMatrix(udata->NEQ, udata->NEQ);
      tmp = N_VNew_Serial(udata->NEQ);
      udata->prels = SUNLinSol_Dense(tmp, udata->pre);
      N_VDestroy(tmp);
    }
    else
    {
      /* Create MPI task-local data structures for mesh node solves */
      udata->b_node   = N_VNew_Serial(udata->nvar);
      udata->jac_node = SUNDenseMatrix(udata->nvar, udata->nvar);
      udata->ls_node  = SUNLinSol_Dense(udata->b_node, udata->jac_node);
    }
  }
  udata->nnlfi = 0;

  /* Open output files for results */
  if (uopt->nout > 0)
  {
      if (udata->myid == 0)
      {
        sprintf(fname, "%s/t.%06d.txt", uopt->outputdir, udata->myid);
        uopt->TFID = fopen(fname, "w");
      }

      sprintf(fname, "%s/u.%06d.txt", uopt->outputdir, udata->myid);
      uopt->UFID = fopen(fname, "w");

      sprintf(fname, "%s/v.%06d.txt", uopt->outputdir, udata->myid);
      uopt->VFID = fopen(fname, "w");

      sprintf(fname, "%s/w.%06d.txt", uopt->outputdir, udata->myid);
      uopt->WFID = fopen(fname, "w");
  }

  /* Print problem setup */
  if (udata->myid == 0)
  {
    printf("\n1D Advection-Reaction Test Problem\n\n");
    printf("Number of Processors = %li\n", (long int) udata->nprocs);
    printf("Mesh Info:\n");
    printf("  NX = %lli, NXL = %lli, dx = %.6f, xmax = %.6f\n",
           udata->nx, udata->nxl, udata->dx, udata->xmax);
    printf("Problem Parameters:\n");
    printf("  A = %g\n", udata->A);
    printf("  B = %g\n", udata->B);
    printf("  k = %g\n", udata->k1);
    printf("  c = %g\n", udata->c);
    printf("Integrator Options:\n");
    printf("  order            = %d\n", uopt->order);
    printf("  method           = %s\n", uopt->explicit ? "explicit" : "imex");
    printf("  fused vector ops = %d\n", uopt->fused);
    printf("  t0               = %g\n", uopt->t0);
    printf("  tf               = %g\n", uopt->tf);
    printf("  reltol           = %.1e\n", uopt->rtol);
    printf("  abstol           = %.1e\n", uopt->atol);
    printf("  nout             = %d\n", uopt->nout);
    if (uopt->explicit)
      printf("  nonlinear solver = none\n");
    else
      printf("  nonlinear solver = %s\n", uopt->global ? "global" : "task local");
    printf("Output directory: %s\n", uopt->outputdir);
  }


  /* return success */
  return(0);
}


void FreeProblem(UserData udata, UserOptions uopt)
{
  if (!uopt->explicit)
  {
    if (uopt->global)
    {
      /* free task-local preconditioner solve structures */
      SUNMatDestroy(udata->pre);
      SUNLinSolFree(udata->prels);
    }
    else
    {
      /* free task-local solve structures */
      N_VDestroy(udata->b_node);
      SUNMatDestroy(udata->jac_node);
      SUNLinSolFree(udata->ls_node);
    }
  }

  /* free solution masks */
  N_VDestroy(N_VGetLocalVector_MPIPlusX(udata->umask));
  N_VDestroy(udata->umask);
  N_VDestroy(udata->vmask);
  N_VDestroy(udata->wmask);

  /* free exchange buffers */
  free(udata->Wsend);
  free(udata->Wrecv);
  free(udata->Esend);
  free(udata->Erecv);

  /* close output streams */
  if (uopt->nout > 0)
  {
      fclose(uopt->UFID);
      fclose(uopt->VFID);
      fclose(uopt->WFID);
      if (udata->myid == 0) fclose(uopt->TFID);
  }
}


void InputError(char *name)
{
  int myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Barrier(MPI_COMM_WORLD);

  if (myid == 0)
  {
    fprintf(stderr, "\nERROR: Invalid command line input\n");
    fprintf(stderr, "\nCommand line options for %s\n",name);
    fprintf(stderr, "  --help           prints this message\n");
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

/*---- end of file ----*/
