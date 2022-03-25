/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This examples is based on cvAdvDiff_non_p.c
 * by Cohen, Hindmarsh, Byrne, and Serban @ LLNL.
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
 * This program solves the problem with the option for nonstiff
 * systems: ADAMS method and the Anderson mixing method provided
 * by PETSc SNES using the SUNNonlinearSolver_PetscSNES module.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .5, 1.0, ..., 5.
 * Run statistics (optional outputs) are printed at the end.
 *
 * This version uses MPI for user routines.
 * Execute with Number of Processors = N,  with 1 <= N <= MX.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h> /* MPI constants and types */

#include <cvode/cvode.h>                          /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_petsc.h>                /* access to the PETSc N_Vector                 */
#include "sunnonlinsol/sunnonlinsol_petscsnes.h"  /* access to the fixed point SUNNonlinearSolver */
#include <sundials/sundials_types.h>              /* definition of type realtype                  */
#include <sundials/sundials_math.h>               /* definition of ABS and EXP                    */

/* Problem Constants */

#define ZERO  RCONST(0.0)

#define XMAX  RCONST(2.0)    /* domain boundary           */
#define MX    10             /* mesh dimension            */
#define NEQ   MX             /* number of equations       */
#define ATOL  RCONST(1.0e-5) /* scalar absolute tolerance */
#define T0    ZERO           /* initial time              */
#define T1    RCONST(0.5)    /* first output time         */
#define DTOUT RCONST(0.5)    /* output time increment     */
#define NOUT  10             /* number of output times    */

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

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  SUNContext sunctx;
  N_Vector u;
  SUNNonlinearSolver NLS;
  UserData data;
  void *cvode_mem;
  int iout, retval, my_pe, npes;
  long int nst;
  PetscErrorCode ptcerr;
  sunindextype local_N, nperpe, nrem, my_base;
  realtype dx, reltol, abstol, t, tout, umax;
  MPI_Comm comm;

  SNES snes;
  Vec uvec;

  u = NULL;
  data = NULL;
  cvode_mem = NULL;

  /* Get processor number, total number of pe's, and my_pe. */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  /* Initialize PETSc */
  retval = PetscInitialize(&argc, &argv, (char*)0, NULL);
  if(check_retval(&retval, "PetscInitialize", 1, my_pe)) MPI_Abort(comm, 1);

  /* Create SUNDIALS context */
  retval = SUNContext_Create(&comm, &sunctx);
  if(check_retval(&retval, "SUNContext_Create", 1, my_pe)) MPI_Abort(comm, 1);

  /* Set local vector length. */
  nperpe = NEQ/npes;
  nrem = NEQ - npes*nperpe;
  local_N = (my_pe < nrem) ? nperpe+1 : nperpe;
  my_base = (my_pe < nrem) ? my_pe*local_N : my_pe*nperpe + nrem;

  /* Create the vector */
  ptcerr = VecCreateMPI(comm, local_N, NEQ, &uvec); CHKERRQ(ptcerr);
  u = N_VMake_Petsc(uvec, sunctx);
  if(check_retval((void *)u, "N_VMake_Petsc", 0, my_pe)) MPI_Abort(comm, 1);

  /* Create the user data structure */
  data = (UserData) malloc(sizeof *data);  /* Allocate data memory */
  if(check_retval((void *)data, "malloc", 2, my_pe)) MPI_Abort(comm, 1);

  data->comm = comm;
  data->npes = npes;
  data->my_pe = my_pe;

  /* Set grid coefficients */
  dx = data->dx = XMAX/((realtype)(MX+1));
  data->hdcoef = RCONST(1.0)/(dx*dx);
  data->hacoef = RCONST(0.5)/(RCONST(2.0)*dx);

  /* Set the tolerances */
  reltol = ZERO;
  abstol = ATOL;

  /* Set the initial conditions */
  SetIC(u, dx, local_N, my_base);

  /* Call CVodeCreate to create the solver memory and specify the Adams-Moulton LMM */
  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0, my_pe)) MPI_Abort(comm, 1);

  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1, my_pe)) MPI_Abort(comm, 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1, my_pe)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_retval(&retval, "CVodeSStolerances", 1, my_pe)) return(1);

  /* Create SNES nonlinear solver object */
  ptcerr = SNESCreate(comm, &snes); CHKERRQ(ptcerr);
  NLS = SUNNonlinSol_PetscSNES(u, snes, sunctx); /* This will call SNESSetFunction appropriately */
  if(check_retval((void *)NLS, "SUNNonlinSol_PetscSNES", 0, my_pe)) return(1);

  /* Set SNES options */
  ptcerr = SNESSetType(snes, SNESANDERSON); CHKERRQ(ptcerr);

  /* Attach nonlinear solver object to CVode */
  retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if(check_retval(&retval, "CVodeSetNonlinearSolver", 1, my_pe)) return(1);

  /* Print the problem introduction header */
  if (my_pe == 0) PrintIntro(npes);

  /* Calculate and print initial maxnorm of u */
  umax = N_VMaxNorm(u);
  if (my_pe == 0) {
    t = T0;
    PrintData(t, umax, 0);
  }

  /* In loop over output points, call CVode, print results, test for error */
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_retval(&retval, "CVode", 1, my_pe)) break;
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    if(check_retval(&retval, "CVodeGetNumSteps", 1, my_pe)) break;
    umax = N_VMaxNorm(u);
    if (my_pe == 0) PrintData(t, umax, nst);
  }

  /* Print some final statistics */
  if (my_pe == 0) PrintFinalStats(cvode_mem);

  /* Free data structures */
  VecDestroy(&uvec);             /* Free the petsc u vector */
  SNESDestroy(&snes);            /* Free the petsc SNES context */
  N_VDestroy(u);                 /* Free the sundials u vector */
  CVodeFree(&cvode_mem);         /* Free the cvode integrator memory */
  SUNNonlinSolFree(NLS);         /* Free the sundials nonlinear solver */
  free(data);                    /* Free user data */
  SUNContext_Free(&sunctx);      /* Free context */

  MPI_Finalize();

  return(0);
}

/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */
static void SetIC(N_Vector u, realtype dx, sunindextype my_length,
                  sunindextype my_base)
{
  int i, iglobal;
  realtype x;
  Vec uvec;

  /* Set pointer to data array and get local length of u. */
  uvec = N_VGetVector_Petsc(u);

  /* Load initial profile into u vector */
  for (i=1; i<=my_length; i++) {
    iglobal = my_base + i; x = iglobal*dx;
    VecSetValue(uvec, iglobal-1, x*(XMAX - x)*SUNRexp(RCONST(2.0)*x), INSERT_VALUES);
  }

  VecAssemblyBegin(uvec);
  VecAssemblyEnd(uvec);
}

/* Print problem introduction */
static void PrintIntro(int npes)
{
  printf("\n 1-D advection-diffusion equation, mesh size =%3d \n", MX);
  printf("\n Number of PEs = %3d \n\n", npes);
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
  PetscErrorCode ptcerr;
  realtype ui, ult, urt, hordc, horac, hdiff, hadv;
  realtype *udata, *dudata, *z;
  sunindextype my_length;
  int i, npes, my_pe, my_pe_m1, my_pe_p1, last_pe;
  UserData data;
  MPI_Status status;
  MPI_Comm comm;

  ptcerr = VecGetArray(N_VGetVector_Petsc(u), &udata); CHKERRQ(ptcerr);
  ptcerr = VecGetArray(N_VGetVector_Petsc(udot), &dudata); CHKERRQ(ptcerr);

  /* Extract needed problem constants from data */
  data = (UserData) user_data;
  hordc = data->hdcoef;
  horac = data->hacoef;

  /* Extract parameters for parallel computation. */
  comm = data->comm;    /* MPI communicatopr */
  npes = data->npes;    /* Number of processes. */
  my_pe = data->my_pe;  /* Current process number. */
  z = data->z;          /* Work array */
  ptcerr = VecGetLocalSize(N_VGetVector_Petsc(u), &my_length); CHKERRQ(ptcerr); /* Number of local elements of u. */

  /* Compute related parameters. */
  my_pe_m1 = my_pe - 1;
  my_pe_p1 = my_pe + 1;
  last_pe = npes - 1;

  /* Store local segment of u in the working array z. */
   for (i = 1; i <= my_length; i++)
     z[i] = udata[i - 1];

  /* Pass needed data to processes before and after current process. */
   if (my_pe != 0)
     (void) MPI_Send(&z[1], 1, MPI_SUNREALTYPE, my_pe_m1, 0, comm);
   if (my_pe != last_pe)
     (void) MPI_Send(&z[my_length], 1, MPI_SUNREALTYPE, my_pe_p1, 0, comm);

  /* Receive needed data from processes before and after current process. */
   if (my_pe != 0)
     (void) MPI_Recv(&z[0], 1, MPI_SUNREALTYPE, my_pe_m1, 0, comm, &status);
   else z[0] = ZERO;
   if (my_pe != last_pe)
     (void) MPI_Recv(&z[my_length+1], 1, MPI_SUNREALTYPE, my_pe_p1, 0, comm,
                     &status);
   else z[my_length + 1] = ZERO;

  /* Loop over all grid points in current process. */
  for (i=1; i<=my_length; i++) {

    /* Extract u at x_i and two neighboring points */
    ui = z[i];
    ult = z[i-1];
    urt = z[i+1];

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - RCONST(2.0)*ui + urt);
    hadv = horac*(urt - ult);
    dudata[i-1] = hdiff + hadv;

  }

  ptcerr = VecRestoreArray(N_VGetVector_Petsc(u), &udata); CHKERRQ(ptcerr);
  ptcerr = VecRestoreArray(N_VGetVector_Petsc(udot), &dudata); CHKERRQ(ptcerr);

  return(0);
}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval < 0
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
