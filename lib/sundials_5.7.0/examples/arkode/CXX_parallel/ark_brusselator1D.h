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
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <mpi.h>

#include "RAJA/RAJA.hpp"

#include "arkode/arkode_arkstep.h"                  /* ARKStep                   */
#include "arkode/arkode_erkstep.h"                  /* ERKStep                   */
#include "nvector/nvector_mpiplusx.h"               /* MPI+X N_Vector            */
#include "sunlinsol/sunlinsol_spgmr.h"              /* GMRES SUNLinearSolver     */
#include "sunnonlinsol/sunnonlinsol_newton.h"       /* Newton SUNNonlinearSolver */

#if defined(USE_RAJA_NVEC)
#include "nvector/nvector_raja.h"                   /* RAJA N_Vector             */
#elif defined(USE_OMPDEV_NVEC)
#include <omp.h>
#include "nvector/nvector_openmpdev.h"              /* OpenMPDEV N_Vector        */
#elif defined(USE_HIP_NVEC)
#include "nvector/nvector_hip.h"                    /* HIP N_Vector              */
#elif defined(USE_CUDA_NVEC) || defined(USE_CUDAUVM_NVEC)
#include "nvector/nvector_cuda.h"                   /* CUDA N_Vector             */
#else
#include "nvector/nvector_serial.h"                 /* serial N_Vector           */
#endif

#if defined(USE_CUDA_NVEC) || defined(USE_CUDAUVM_NVEC) || defined(USE_RAJA_NVEC)
#define USE_CUDA
#define USE_GPU
#elif defined(USE_OMPDEV_NVEC)
#define USE_GPU
#elif defined(USE_HIP_NVEC)
#define USE_HIP
#define USE_GPU
#endif

#if defined(USE_CUDA) || defined(USE_HIP)
#define USE_CUDA_OR_HIP
#endif

#if defined(USE_CUDA)
#define GPU_PREFIX(a) cuda##a
#elif defined(USE_HIP)
#define GPU_PREFIX(a) hip##a
#endif

/* Maximum size of output directory string */
#define MXSTR 2048

/* Accessor macro:
   n = number of state variables
   i = mesh node index
   c = component */
#define IDX(n,i,c) ((n)*(i)+(c))

/* function to sync the host and device */
static void sync_device();

/*
 * Simple timer class.
 */

class Timer
{
  public:
    Timer() : total_(0.0), start_(0.0), end_(0.0)  {}
    void start() { start_ = MPI_Wtime(); }
    void stop() {
      sync_device();
      end_ = MPI_Wtime();
      total_ += (end_-start_);
    }
    double total() const { return total_; }
  private:
    double total_;
    double start_;
    double end_;
};


/*
 * User options structure
 */

struct UserOptions
{
  double t0;       /* initial time                 */
  double tf;       /* final time                   */
  double rtol;     /* relative tolerance           */
  double atol;     /* absolute tolerance           */
  int    order;    /* method order                 */
  int    expl;     /* imex method or explicit      */
  int    global;   /* use global nonlinear solve   */
  int    fused;    /* use fused vector ops         */
  int    nout;     /* number of outputs            */
  int    monitor;  /* print solution to screen     */
  int    printtime;/* print timing information     */
  char*  outputdir;
};


/*
 * User data structure
 */

struct UserData
{
  /* MPI data */
  MPI_Comm    comm;
  int         myid;
  int         nprocs;
  MPI_Request req[2];
  double*     Wsend;
  double*     Esend;
  double*     Wrecv;
  double*     Erecv;

  /* file handles for output */
  FILE*  TFID;     /* time output file pointer     */
  FILE*  UFID;     /* solution output file pointer */
  FILE*  VFID;
  FILE*  WFID;

  /* solution masks */
  N_Vector umask;
  N_Vector vmask;
  N_Vector wmask;

  /* problem paramaters */
  int       nvar; /* number of species            */
  long int  nx;   /* number of intervals globally */
  int       nxl;  /* number of intervals locally  */
  int       NEQ;  /* number of equations locally  */
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

  /* count of implicit function evals by the task local nonlinear solver */
  long int nnlfi;

  /* integrator options */
  UserOptions* uopt;

  /* code timers */
  Timer t_overall;
  Timer t_setup;  /* initial setup */
  Timer t_comm;   /* communication in advection operator */
  Timer t_advec;  /* advection computation */
  Timer t_react;  /* reaction computation  */
  Timer t_lsolve; /* linear solve          */
  Timer t_psolve; /* preconditioner solve  */

  /* destructor frees the problem data */
  ~UserData();
};


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
 * Linear solver functions.
 */

int SolveJacBlocks(N_Vector y, N_Vector x, N_Vector b,
                   double gamma, RAJA::RangeSegment blocks,
                   UserData* udata);

/*
 * Preconditioner function (used only when using the global nonlinear solver)
 */

static int PSolve(double t, N_Vector y, N_Vector f, N_Vector r,
                  N_Vector z, double gamma, double delta, int lr,
                  void *user_data);


/*
 * Helper functions
 */

/* function that does ARKStep setup and evolves the solution */
static int EvolveProblemIMEX(N_Vector y, UserData* udata, UserOptions* uopt);

/* function that does ERKStep setup and evolves the solution */
static int EvolveProblemExplicit(N_Vector y, UserData* udata, UserOptions* uopt);

/* function to set initial condition */
static int SetIC(N_Vector y, UserData* udata);

/* functions to get the data pointers from the N_Vector */
static double* GetVecData(N_Vector y);

/* function to enable fused vector operations */
static int EnableFusedVectorOps(N_Vector y);

/* functions to exchange neighbor data */
static int ExchangeBCOnly(N_Vector y, UserData* udata);
static int ExchangeAllStart(N_Vector y, UserData* udata);
static int ExchangeAllEnd(UserData* udata);

/* functions for processing command line args */
static int SetupProblem(int argc, char *argv[],
                        UserData* udata, UserOptions* uopt);
static void InputError(char *name);

/* function to write solution to disk */
static int WriteOutput(double t, N_Vector y, UserData* udata, UserOptions* uopt);

/* function to check sundials return values */
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* function to check if GPU operation returned successfully */
#ifdef USE_CUDA_OR_HIP
static void gpuAssert(GPU_PREFIX(Error_t) code, const char *file, int line, int abort);
#endif

