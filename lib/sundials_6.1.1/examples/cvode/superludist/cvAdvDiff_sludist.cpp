/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * Based on cvAdvDiff_non_p.c by Scott D. Cohen, Alan C. Hindmarsh,
 * George Byrne, and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the program for
 * its solution by CVODE. The problem is the semi-discrete
 * form of the advection-diffusion equation in 1-D:
 *   du/dt = d^2 u / dx^2 + .5 du/dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is the following:
 *   u(x,t=0) = x(2-x)exp(2x) .
 * The PDE is discretized on a uniform grid of size MX+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX.
 * This program solves the problem with BDF and Newton Iteration.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .5, 1.0, ..., 5.
 * Run statistics (optional outputs) are printed at the end.
 *
 * This version uses MPI for user routines and the
 * SuperLU-DIST Linear Solver.
 *
 * Execute with Number of Processors = N,  with 1 <= N <= MX.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "superlu_ddefs.h"

#include <cvode/cvode.h>                          /* prototypes for CVODE fcts., consts.          */
#include <cvode/cvode_direct.h>                   /* CVODE direct linear solver interface         */
#include <nvector/nvector_parallel.h>             /* access to MPI-parallel N_Vector              */
#include <sunlinsol/sunlinsol_superludist.h>      /* access to the SuperLU-DIST SUNLinearSolver   */
#include <sunmatrix/sunmatrix_slunrloc.h>         /* access to the SuperLU SLU_NR_loc SUNMatrix   */
#include <sundials/sundials_types.h>              /* definition of type realtype                  */

#include <mpi.h> /* MPI constants and types */

/* Problem Constants */

#define XMAX  RCONST(2.0)    /* domain boundary           */
#define MX    10             /* mesh dimension            */
#define NEQ   MX             /* number of equations       */
#define ATOL  RCONST(1.0e-5) /* scalar absolute tolerance */
#define T0    ZERO           /* initial time              */
#define T1    RCONST(0.5)    /* first output time         */
#define DTOUT RCONST(0.5)    /* output time increment     */
#define NOUT  10             /* number of output times    */

#define ZERO  RCONST(0.0)
#define HALF  RCONST(0.5)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define FIVE  RCONST(5.0)

/* Type : UserData
   contains grid constants, parallel machine parameters, work array. */

typedef struct {
  realtype dx, hdcoef, hacoef;
  int npes, my_pe;
  MPI_Comm comm;
  realtype z[100];
} *UserData;

/* Private Helper Functions */

static void SetIC(N_Vector u, realtype dx, sunindextype my_length,
                  sunindextype my_base);

static void PrintIntro(int npes);

static void PrintData(realtype t, realtype umax, long int nst);

static void PrintFinalStats(void *cvode_mem);

/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);
static int Jac(realtype t, N_Vector u, N_Vector fu,
               SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  SUNContext sunctx;
  realtype dx, reltol, abstol, t, tout, umax;
  UserData data;
  void *cvode_mem;
  int iout, retval, my_pe, npes, nprow, npcol;
  sunindextype local_N, local_NNZ, nperpe, nrem, my_base;
  long int nst;

  gridinfo_t grid;
  dLUstruct_t LUstruct;
  dScalePermstruct_t scaleperm;
  dSOLVEstruct_t solve;
  SuperLUStat_t stat;
  superlu_dist_options_t options;
  SuperMatrix Asuper;
  sunindextype *rowptr, *colind;
  realtype *matdata;

  N_Vector u;
  SUNMatrix A;
  SUNLinearSolver LS;

  data      = NULL;
  cvode_mem = NULL;
  u         = NULL;
  A         = NULL;
  LS        = NULL;

  /* Get processor number, total number of pe's, and my_pe. */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_pe);

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if(check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* check for nprow and npcol arguments */
  if (argc < 2) {
    printf("ERROR: number of process rows and columns must be provided as arguments: ./cvAdvDiff <nprow> <npcol>\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* Initialize SuperLU-DIST process grid */
  nprow = atoi(argv[1]); npcol = atoi(argv[2]);
  superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
  /* Excess processes just exit */
  if (grid.iam >= nprow*npcol) {
    superlu_gridexit(&grid);
    MPI_Finalize();
    return(0);
  }
  npes = nprow*npcol;

  /* Set local vector length and matrix size. */
  nperpe    = NEQ/npes;
  nrem      = NEQ - npes*nperpe;
  local_N   = (my_pe < nrem) ? nperpe+1 : nperpe;
  my_base   = (my_pe < nrem) ? my_pe*local_N : my_pe*nperpe + nrem;

  data = (UserData) malloc(sizeof *data);  /* Allocate data memory */
  if(check_retval((void *)data, "malloc", 2, my_pe)) MPI_Abort(grid.comm, 1);

  data->comm = grid.comm;
  data->npes = npes;
  data->my_pe = my_pe;

  u = N_VNew_Parallel(grid.comm, local_N, NEQ, sunctx);  /* Allocate u vector */
  if(check_retval((void *)u, "N_VNew", 0, my_pe)) MPI_Abort(grid.comm, 1);

  reltol = ZERO;  /* Set the tolerances */
  abstol = ATOL;

  dx = data->dx = XMAX/((realtype)(MX+1));  /* Set grid coefficients in data */
  data->hdcoef = ONE/(dx*dx);
  data->hacoef = HALF/(TWO*dx);

  SetIC(u, dx, local_N, my_base);  /* Initialize u vector */

  if (!my_pe || (my_pe == (npes-1)))
    local_NNZ = 2 + 3*(local_N-1);
  else
    local_NNZ = 3*local_N;

  /* Create the SuperLU-DIST SuperMatrix which will be wrapped as A */
  matdata = (realtype *) malloc(local_NNZ*sizeof(realtype));
  colind  = (sunindextype *) calloc(local_NNZ, sizeof(sunindextype));
  rowptr  = (sunindextype *) calloc((local_N+1), sizeof(sunindextype));
  dCreate_CompRowLoc_Matrix_dist(&Asuper, NEQ, NEQ, local_NNZ, local_N, my_base, matdata, colind, rowptr,
                                 SLU_NR_loc, SLU_D, SLU_GE);

  /* Use the default SuperLU-DIST solver options */
  set_default_options_dist(&options);
  options.PrintStat = NO;

  /* Initialize SuperLU-DIST solver structures */
  dScalePermstructInit(NEQ, NEQ, &scaleperm);
  dLUstructInit(NEQ, &LUstruct);
  PStatInit(&stat);

  /* Call CVodeCreate to create the solver memory and specify the Adams-Moulton LMM */
  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0, my_pe)) MPI_Abort(grid.comm, 1);

  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1, my_pe)) MPI_Abort(grid.comm, 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1, my_pe)) MPI_Abort(grid.comm, 1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1, my_pe)) MPI_Abort(grid.comm, 1);

  /* create the SuperLU SLU_NR_loc SUNMatrix */
  A = SUNMatrix_SLUNRloc(&Asuper, &grid, sunctx);
  if (check_retval((void *)A, "SUNMatrix_SLUNRloc", 0, my_pe)) MPI_Abort(grid.comm, 1);

  /* create SuperLU-DIST linear solver object */
  LS = SUNLinSol_SuperLUDIST(u, A, &grid, &LUstruct, &scaleperm, &solve, &stat, &options, sunctx);
  if (check_retval((void *)LS, "SUNLinSol_SuperLUDIST", 0, my_pe)) MPI_Abort(grid.comm, 1);

  /* attach linear solver object to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1, my_pe)) MPI_Abort(grid.comm, 1);

  /* attach the Jacobian function */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if (check_retval(&retval, "CVodeSetJacFn", 1, my_pe)) MPI_Abort(grid.comm, 1);

  if (my_pe == 0) PrintIntro(npes);

  umax = N_VMaxNorm(u);

  if (my_pe == 0) {
    t = T0;
    PrintData(t, umax, 0);
  }

  /* In loop over output points, call CVode, print results, test for error */

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_retval(&retval, "CVode", 1, my_pe)) break;
    umax = N_VMaxNorm(u);
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_retval(&retval, "CVodeGetNumSteps", 1, my_pe);
    if (my_pe == 0) PrintData(t, umax, nst);
  }

  if (my_pe == 0)
    PrintFinalStats(cvode_mem);  /* Print some final statistics */

  /* Free SUNDIALS structures */
  N_VDestroy(u);                 /* Free the u vector */
  SUNLinSolFree(LS);             /* Free the linear solver */
  SUNMatDestroy(A);              /* Free the A matrix */
  CVodeFree(&cvode_mem);         /* Free the integrator memory */
  free(data);                    /* Free user data */
  SUNContext_Free(&sunctx);

  /* Free the SuperLU_DIST structures */
  PStatFree(&stat);
  dScalePermstructFree(&scaleperm);
  dLUstructFree(&LUstruct);
  Destroy_CompRowLoc_Matrix_dist(&Asuper);
  superlu_gridexit(&grid);

  MPI_Finalize();

  return(0);
}

/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */

static void SetIC(N_Vector u, realtype dx, sunindextype my_length,
                  sunindextype my_base)
{
  int i;
  sunindextype iglobal;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u. */
  udata = N_VGetArrayPointer(u);
  my_length = N_VGetLocalLength_Parallel(u);

  /* Load initial profile into u vector */
  for (i=1; i<=my_length; i++) {
    iglobal = my_base + i;
    x = iglobal*dx;
    udata[i-1] = x*(XMAX - x)*exp(TWO*x);
  }
}

/* Print problem introduction */

static void PrintIntro(int npes)
{
  printf("\n 1-D advection-diffusion equation, mesh size =%3d \n", MX);
  printf("\n Number of PEs = %3d \n\n", npes);

  return;
}

/* Print data */

static void PrintData(realtype t, realtype umax, long int nst)
{

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %4.2Lf  max.norm(u) =%14.6Le  nst =%4ld \n", t, umax, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %4.2f  max.norm(u) =%14.6e  nst =%4ld \n", t, umax, nst);
#else
  printf("At t = %4.2f  max.norm(u) =%14.6e  nst =%4ld \n", t, umax, nst);
#endif

  return;
}

/* Print some final statistics located in the iopt array */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nni, ncfn, netf;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1, 0);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1, 0);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1, 0);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1, 0);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1, 0);

  printf("\nFinal Statistics: \n\n");
  printf("nst = %-6ld  nfe  = %-6ld  ", nst, nfe);
  printf("nni = %-6ld  ncfn = %-6ld  netf = %ld\n \n", nni, ncfn, netf);
}

/***************** Function Called by the Solver ***********************/

/* f routine. Compute f(t,u). */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype ui, ult, urt, hordc, horac, hdiff, hadv;
  realtype *udata, *dudata, *z;
  int i;
  int npes, my_pe, my_length, my_pe_m1, my_pe_p1, last_pe;
  UserData data;
  MPI_Status status;
  MPI_Comm comm;

  udata = N_VGetArrayPointer(u);
  dudata = N_VGetArrayPointer(udot);

  /* Extract needed problem constants from data */
  data = (UserData) user_data;
  hordc = data->hdcoef;
  horac = data->hacoef;

  /* Extract parameters for parallel computation. */
  comm = data->comm;
  npes = data->npes;           /* Number of processes. */
  my_pe = data->my_pe;         /* Current process number. */
  my_length = N_VGetLocalLength_Parallel(u); /* Number of local elements of u. */
  z = data->z;

  /* Compute related parameters. */
  my_pe_m1 = my_pe - 1;
  my_pe_p1 = my_pe + 1;
  last_pe = npes - 1;

  /* Store local segment of u in the working array z. */
   for (i = 1; i <= my_length; i++)
     z[i] = udata[i - 1];

  /* Pass needed data to processes before and after current process. */
   if (my_pe != 0)
     MPI_Send(&z[1], 1, MPI_SUNREALTYPE, my_pe_m1, 0, comm);
   if (my_pe != last_pe)
     MPI_Send(&z[my_length], 1, MPI_SUNREALTYPE, my_pe_p1, 0, comm);

  /* Receive needed data from processes before and after current process. */
   if (my_pe != 0)
     MPI_Recv(&z[0], 1, MPI_SUNREALTYPE, my_pe_m1, 0, comm, &status);
   else z[0] = ZERO;
   if (my_pe != last_pe)
     MPI_Recv(&z[my_length+1], 1, MPI_SUNREALTYPE, my_pe_p1, 0, comm,
              &status);
   else z[my_length + 1] = ZERO;

  /* Loop over all grid points in current process. */
  for (i=1; i<=my_length; i++) {

    /* Extract u at x_i and two neighboring points */
    ui = z[i];
    ult = z[i-1];
    urt = z[i+1];

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - TWO*ui + urt);
    hadv = horac*(urt - ult);
    dudata[i-1] = hdiff + hadv;
  }

  return(0);
}

/* Jacobian routine. Compute J(t,u). */

static int Jac(realtype t, N_Vector u, N_Vector fu,
               SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunindextype i, j;
  realtype *nzval, hordc, horac;
  UserData data;
  SuperMatrix *Jsuper;
  NRformat_loc *Jstore;

  /*
   * The components of f = udot that depend on u(i) are
   * f(i), f(i-1), and f(i+1) with
   *   df(i)/du(i) = -2/dx^2
   *   df(i-1)/du(i) = 1/dx^2 + .25/dx  (if i > 0)
   *   df(i+1)/du(i) = 1/dx^2 - .25/dx  (if i < MX-1)
   */

  data = (UserData) user_data;
  hordc = data->hdcoef;
  horac = data->hacoef;

  Jsuper = SUNMatrix_SLUNRloc_SuperMatrix(J);
  Jstore = (NRformat_loc *) Jsuper->Store;
  nzval  = (realtype *) Jstore->nzval;
  sunindextype *colind = Jstore->colind;
  sunindextype *rowptr = Jstore->rowptr;

  /* set non-zero Jacobian entries */
  for (i=0; i < Jstore->m_loc; i++) {
    booleantype first_local_row;
    booleantype first_row, last_row;

    /* global row index */
    j = Jstore->fst_row + i;

    first_local_row = (i == 0);
    first_row       = (j == 0);
    last_row        = (j == (Jsuper->nrow-1));

    if (first_local_row) {
      rowptr[i] = 0;
      rowptr[i+1] = first_row ? 2 : 3;
    } else {
      rowptr[i+1] = (Jstore->fst_row == 0 || last_row) ? 3*i+2 :3*i+3;
    }

    /* local non-zero elements */
    if (first_row) {
      /* main diagonal */
      colind[rowptr[i]] = j;
      nzval[rowptr[i]]  = -TWO*hordc;

      colind[rowptr[i]+1] = j+1;
      nzval[rowptr[i]+1]  = hordc+horac;
    } else if (!first_row && !last_row) {
      colind[rowptr[i]] = j-1;
      nzval[rowptr[i]]  = hordc-horac;

      /* main diagonal */
      colind[rowptr[i]+1] = j;
      nzval[rowptr[i]+1]  = -TWO*hordc;

      colind[rowptr[i]+2] = j+1;
      nzval[rowptr[i]+2]  = hordc+horac;
    } else {
      colind[rowptr[i]] = j-1;
      nzval[rowptr[i]]  = hordc-horac;

      /* main diagonal */
      colind[rowptr[i]+1] = j;
      nzval[rowptr[i]+1]  = -TWO*hordc;
    }
  }

  return(0);
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed with retval = %d\n\n",
              id, funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1); }

  return(0);
}
