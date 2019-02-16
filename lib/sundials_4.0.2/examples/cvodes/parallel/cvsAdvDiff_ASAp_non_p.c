/*
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
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
 * its solution by CVODE. The problem is the semi-discrete form of
 * the advection-diffusion equation in 1-D:
 *   du/dt = p1 * d^2u / dx^2 + p2 * du / dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is:
 *   u(x,t=0) = x(2-x)exp(2x).
 * The nominal values of the two parameters are: p1=1.0, p2=0.5
 * The PDE is discretized on a uniform grid of size MX+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX.
 * This program solves the problem with the option for nonstiff
 * systems: ADAMS method and functional iteration.
 * It uses scalar relative and absolute tolerances.
 *
 * In addition to the solution, sensitivities with respect to p1
 * and p2 as well as with respect to initial conditions are
 * computed for the quantity:
 *    g(t, u, p) = int_x u(x,t) at t = 5
 * These sensitivities are obtained by solving the adjoint system:
 *    dv/dt = -p1 * d^2 v / dx^2 + p2 * dv / dx
 * with homogeneous Ditrichlet boundary conditions and the final
 * condition:
 *    v(x,t=5) = 1.0
 * Then, v(x, t=0) represents the sensitivity of g(5) with respect
 * to u(x, t=0) and the gradient of g(5) with respect to p1, p2 is
 *    (dg/dp)^T = [  int_t int_x (v * d^2u / dx^2) dx dt ]
 *                [  int_t int_x (v * du / dx) dx dt     ]
 *
 * This version uses MPI for user routines.
 * Execute with Number of Processors = N,  with 1 <= N <= MX.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvodes/cvodes.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to the fixed point SUNNonlinearSolver */

#include <mpi.h>

/* Problem Constants */

#define XMAX  RCONST(2.0)   /* domain boundary            */
#define MX    20            /* mesh dimension             */
#define NEQ   MX            /* number of equations        */
#define ATOL  RCONST(1.e-5) /* scalar absolute tolerance  */
#define T0    RCONST(0.0)   /* initial time               */
#define TOUT  RCONST(2.5)   /* output time increment      */

/* Adjoint Problem Constants */

#define NP    2            /* number of parameters       */
#define STEPS 200          /* steps between check points */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/* Type : UserData */

typedef struct {
  realtype p[2];            /* model parameters                         */
  realtype dx;              /* spatial discretization grid              */
  realtype hdcoef, hacoef;  /* diffusion and advection coefficients     */
  sunindextype local_N;
  int npes, my_pe;          /* total number of processes and current ID */
  sunindextype nperpe, nrem;
  MPI_Comm comm;            /* MPI communicator                         */
  realtype *z1, *z2;        /* work space                               */
} *UserData;

/* Prototypes of user-supplied funcitons */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);
static int fB(realtype t, N_Vector u, 
              N_Vector uB, N_Vector uBdot, void *user_dataB);

/* Prototypes of private functions */

static void SetIC(N_Vector u, realtype dx, sunindextype my_length, sunindextype my_base);
static void SetICback(N_Vector uB, sunindextype my_base);
static realtype Xintgr(realtype *z, sunindextype l, realtype dx);
static realtype Compute_g(N_Vector u, UserData data);
static void PrintOutput(realtype g_val, N_Vector uB, UserData data);
static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;

  void *cvode_mem;
  
  N_Vector u;
  realtype reltol, abstol;

  int indexB;
  N_Vector uB;

  realtype dx, t, g_val;
  int retval, my_pe, nprocs, npes, ncheck;
  sunindextype local_N=0, nperpe, nrem, my_base=-1;

  SUNNonlinearSolver NLS, NLSB;

  MPI_Comm comm;

  data = NULL;
  cvode_mem = NULL;
  u = uB = NULL;

  /*------------------------------------------------------
    Initialize MPI and get total number of pe's, and my_pe
    ------------------------------------------------------*/
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &my_pe);

  npes = nprocs - 1; /* pe's dedicated to PDE integration */

  if ( npes <= 0 ) {
    if (my_pe == npes)
      fprintf(stderr, "\nMPI_ERROR(%d): number of processes must be >= 2\n\n",
	      my_pe);
    MPI_Finalize();
    return(1);
  }

  /*-----------------------
    Set local vector length
    -----------------------*/
  nperpe = NEQ/npes;
  nrem = NEQ - npes*nperpe;
  if (my_pe < npes) {

    /* PDE vars. distributed to this proccess */
    local_N = (my_pe < nrem) ? nperpe+1 : nperpe;
    my_base = (my_pe < nrem) ? my_pe*local_N : my_pe*nperpe + nrem;

  } else {

    /* Make last process inactive for forward phase */
    local_N = 0;

  }

  /*-------------------------------------
    Allocate and load user data structure
    -------------------------------------*/
  data = (UserData) malloc(sizeof *data);
  if (check_retval((void *)data , "malloc", 2, my_pe)) MPI_Abort(comm, 1);
  data->p[0] = ONE;
  data->p[1] = RCONST(0.5);
  dx = data->dx = XMAX/((realtype)(MX+1));
  data->hdcoef = data->p[0]/(dx*dx);
  data->hacoef = data->p[1]/(TWO*dx);
  data->comm = comm;
  data->npes = npes;
  data->my_pe = my_pe;
  data->nperpe = nperpe;
  data->nrem = nrem;
  data->local_N = local_N;

  /*------------------------- 
    Forward integration phase
    -------------------------*/
  
  /* Set relative and absolute tolerances for forward phase */
  reltol = ZERO;
  abstol = ATOL;

  /* Allocate and initialize forward variables */
  u = N_VNew_Parallel(comm, local_N, NEQ);
  if (check_retval((void *)u, "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  SetIC(u, dx, local_N, my_base);

  /* Allocate CVODES memory for forward integration */
  cvode_mem = CVodeCreate(CV_ADAMS);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0, my_pe)) MPI_Abort(comm, 1);

  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval(&retval, "CVodeSetUserData", 1, my_pe)) MPI_Abort(comm, 1);

  retval = CVodeInit(cvode_mem, f, T0, u);
  if (check_retval(&retval, "CVodeInit", 1, my_pe)) MPI_Abort(comm, 1);

  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1, my_pe)) MPI_Abort(comm, 1);

  /* create fixed point nonlinear solver object */
  NLS = SUNNonlinSol_FixedPoint(u, 0);
  if(check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0, my_pe)) MPI_Abort(comm, 1);

  /* attach nonlinear solver object to CVode */
  retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if(check_retval(&retval, "CVodeSetNonlinearSolver", 1, my_pe)) MPI_Abort(comm, 1);

  /* Allocate combined forward/backward memory */
  retval = CVodeAdjInit(cvode_mem, STEPS, CV_HERMITE);
  if (check_retval(&retval, "CVadjInit", 1, my_pe)) MPI_Abort(comm, 1);

  /* Integrate to TOUT and collect check point information */
  retval = CVodeF(cvode_mem, TOUT, u, &t, CV_NORMAL, &ncheck);
  if (check_retval(&retval, "CVodeF", 1, my_pe)) MPI_Abort(comm, 1);

  /*---------------------------
    Compute and value of g(t_f)
    ---------------------------*/
  g_val = Compute_g(u, data);

  /*-------------------------- 
    Backward integration phase
    --------------------------*/

  if (my_pe == npes) {

    /* Activate last process for integration of the quadrature equations */
    local_N = NP;

  } else {

    /* Allocate work space */
    data->z1 = (realtype *)malloc(local_N*sizeof(realtype));
    if (check_retval((void *)data->z1, "malloc", 2, my_pe)) MPI_Abort(comm, 1);
    data->z2 = (realtype *)malloc(local_N*sizeof(realtype));
    if (check_retval((void *)data->z2, "malloc", 2, my_pe)) MPI_Abort(comm, 1);

  }

  /* Allocate and initialize backward variables */
  uB = N_VNew_Parallel(comm, local_N, NEQ+NP);
  if (check_retval((void *)uB, "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  SetICback(uB, my_base);

  /* Allocate CVODES memory for the backward integration */
  retval = CVodeCreateB(cvode_mem, CV_ADAMS, &indexB);
  if (check_retval(&retval, "CVodeCreateB", 1, my_pe)) MPI_Abort(comm, 1);

  retval = CVodeSetUserDataB(cvode_mem, indexB, data);
  if (check_retval(&retval, "CVodeSetUserDataB", 1, my_pe)) MPI_Abort(comm, 1);

  retval = CVodeInitB(cvode_mem, indexB, fB, TOUT, uB);
  if (check_retval(&retval, "CVodeInitB", 1, my_pe)) MPI_Abort(comm, 1);

  retval = CVodeSStolerancesB(cvode_mem, indexB, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerancesB", 1, my_pe)) MPI_Abort(comm, 1);

  /* create fixed point nonlinear solver object */
  NLSB = SUNNonlinSol_FixedPoint(uB, 0);
  if(check_retval((void *)NLSB, "SUNNonlinSol_FixedPoint", 0, my_pe)) MPI_Abort(comm, 1);

  /* attach nonlinear solver object to CVode */
  retval = CVodeSetNonlinearSolverB(cvode_mem, indexB, NLSB);
  if(check_retval(&retval, "CVodeSetNonlinearSolver", 1, my_pe)) MPI_Abort(comm, 1);

  /* Integrate to T0 */
  retval = CVodeB(cvode_mem, T0, CV_NORMAL);
  if (check_retval(&retval, "CVodeB", 1, my_pe)) MPI_Abort(comm, 1);

  retval = CVodeGetB(cvode_mem, indexB, &t, uB);
  if (check_retval(&retval, "CVodeGetB", 1, my_pe)) MPI_Abort(comm, 1);

  /* Print results (adjoint states and quadrature variables) */
  PrintOutput(g_val, uB, data);


  /* Free memory */
  N_VDestroy_Parallel(u);
  N_VDestroy_Parallel(uB);
  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);
  SUNNonlinSolFree(NLSB);

  if (my_pe != npes) {
    free(data->z1);
    free(data->z2);
  }
  free(data);

  MPI_Finalize();

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,u) for forward phase. 
 */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype uLeft, uRight, ui, ult, urt;
  realtype hordc, horac, hdiff, hadv;
  realtype *udata, *dudata;
  sunindextype i, my_length;
  int npes, my_pe, my_pe_m1, my_pe_p1, last_pe;
  UserData data;
  MPI_Status status;
  MPI_Comm comm;

  /* Extract MPI info. from data */
  data = (UserData) user_data;
  comm = data->comm;
  npes = data->npes;
  my_pe = data->my_pe;
  
  /* If this process is inactive, return now */
  if (my_pe == npes) return(0);

  /* Extract problem constants from data */
  hordc = data->hdcoef;
  horac = data->hacoef;

  /* Find related processes */
  my_pe_m1 = my_pe - 1;
  my_pe_p1 = my_pe + 1;
  last_pe = npes - 1;

  /* Obtain local arrays */
  udata = N_VGetArrayPointer_Parallel(u);
  dudata = N_VGetArrayPointer_Parallel(udot);
  my_length = N_VGetLocalLength_Parallel(u);

  /* Pass needed data to processes before and after current process. */
   if (my_pe != 0)
     MPI_Send(&udata[0], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm);
   if (my_pe != last_pe)
     MPI_Send(&udata[my_length-1], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm);   

  /* Receive needed data from processes before and after current process. */
   if (my_pe != 0)
     MPI_Recv(&uLeft, 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm, &status);
   else uLeft = ZERO;
   if (my_pe != last_pe)
     MPI_Recv(&uRight, 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm,
              &status);   
   else uRight = ZERO;

  /* Loop over all grid points in current process. */
  for (i=0; i<my_length; i++) {

    /* Extract u at x_i and two neighboring points */
    ui = udata[i];
    ult = (i==0) ? uLeft: udata[i-1];
    urt = (i==my_length-1) ? uRight : udata[i+1];

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - TWO*ui + urt);
    hadv = horac*(urt - ult);
    dudata[i] = hdiff + hadv;
  }

  return(0);
}

/*
 * fB routine. Compute right hand side of backward problem 
 */

static int fB(realtype t, N_Vector u, 
              N_Vector uB, N_Vector uBdot, void *user_dataB)
{
  realtype *uBdata, *duBdata, *udata;
  realtype uBLeft, uBRight, uBi, uBlt, uBrt;
  realtype uLeft, uRight, ui, ult, urt;
  realtype dx, hordc, horac, hdiff, hadv;
  realtype *z1, *z2, intgr1, intgr2;
  sunindextype i, my_length;
  int npes, my_pe, my_pe_m1, my_pe_p1, last_pe;
  UserData data;
  realtype data_in[2], data_out[2];
  MPI_Status status;
  MPI_Comm comm;

  /* Extract MPI info. from data */
  data = (UserData) user_dataB;
  comm = data->comm;
  npes = data->npes;
  my_pe = data->my_pe;

  if (my_pe == npes) { /* This process performs the quadratures */

    /* Obtain local arrays */
    duBdata = N_VGetArrayPointer_Parallel(uBdot);
    my_length = N_VGetLocalLength_Parallel(uB);

    /* Loop over all other processes and load right hand side of quadrature eqs. */
    duBdata[0] = ZERO;
    duBdata[1] = ZERO;
    for (i=0; i<npes; i++) {
      MPI_Recv(&intgr1, 1, PVEC_REAL_MPI_TYPE, i, 0, comm, &status); 
      duBdata[0] += intgr1;
      MPI_Recv(&intgr2, 1, PVEC_REAL_MPI_TYPE, i, 0, comm, &status); 
      duBdata[1] += intgr2;
    }

  } else { /* This process integrates part of the PDE */

    /* Extract problem constants and work arrays from data */
    dx    = data->dx;
    hordc = data->hdcoef;
    horac = data->hacoef;
    z1    = data->z1;
    z2    = data->z2;

    /* Obtain local arrays */
    uBdata = N_VGetArrayPointer_Parallel(uB);
    duBdata = N_VGetArrayPointer_Parallel(uBdot);
    udata = N_VGetArrayPointer_Parallel(u);
    my_length = N_VGetLocalLength_Parallel(uB);

    /* Compute related parameters. */
    my_pe_m1 = my_pe - 1;
    my_pe_p1 = my_pe + 1;
    last_pe  = npes - 1;

    /* Pass needed data to processes before and after current process. */
    if (my_pe != 0) {
      data_out[0] = udata[0];
      data_out[1] = uBdata[0];
    
      MPI_Send(data_out, 2, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm);
    }
    if (my_pe != last_pe) {
      data_out[0] = udata[my_length-1];
      data_out[1] = uBdata[my_length-1];

      MPI_Send(data_out, 2, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm);
    }
    
    /* Receive needed data from processes before and after current process. */
    if (my_pe != 0) {
      MPI_Recv(data_in, 2, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm, &status);
      
      uLeft = data_in[0];
      uBLeft = data_in[1];
    } else {
      uLeft = ZERO;
      uBLeft = ZERO;
    }
    if (my_pe != last_pe) {
      MPI_Recv(data_in, 2, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm, &status);

      uRight = data_in[0];
      uBRight = data_in[1];
    } else {
      uRight = ZERO;
      uBRight = ZERO;
    }

    /* Loop over all grid points in current process. */
    for (i=0; i<my_length; i++) {
      
      /* Extract uB at x_i and two neighboring points */
      uBi = uBdata[i];
      uBlt = (i==0) ? uBLeft: uBdata[i-1];
      uBrt = (i==my_length-1) ? uBRight : uBdata[i+1];
      
      /* Set diffusion and advection terms and load into udot */
      hdiff = hordc*(uBlt - TWO*uBi + uBrt);
      hadv = horac*(uBrt - uBlt);
      duBdata[i] = - hdiff + hadv;

      /* Extract u at x_i and two neighboring points */
      ui = udata[i];
      ult = (i==0) ? uLeft: udata[i-1];
      urt = (i==my_length-1) ? uRight : udata[i+1];

      /* Load integrands of the two space integrals */
      z1[i] = uBdata[i]*(ult - TWO*ui + urt)/(dx*dx);
      z2[i] = uBdata[i]*(urt - ult)/(TWO*dx);
    }

    /* Compute local integrals */
    intgr1 = Xintgr(z1, my_length, dx);
    intgr2 = Xintgr(z2, my_length, dx);

    /* Send local integrals to 'quadrature' process */
    MPI_Send(&intgr1, 1, PVEC_REAL_MPI_TYPE, npes, 0, comm);
    MPI_Send(&intgr2, 1, PVEC_REAL_MPI_TYPE, npes, 0, comm);

  }


  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/* 
 * Set initial conditions in u vector 
 */

static void SetIC(N_Vector u, realtype dx, sunindextype my_length, sunindextype my_base)
{
  int i;
  sunindextype iglobal;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u */
  udata = N_VGetArrayPointer_Parallel(u);
  my_length = N_VGetLocalLength_Parallel(u);

  /* Load initial profile into u vector */
  for (i=1; i<=my_length; i++) {
    iglobal = my_base + i;
    x = iglobal*dx;
    udata[i-1] = x*(XMAX - x)*SUNRexp(TWO*x);
  }  
}

/* 
 * Set final conditions in uB vector 
 */

static void SetICback(N_Vector uB, sunindextype my_base)
{
  int i;
  realtype *uBdata;
  sunindextype my_length;

  /* Set pointer to data array and get local length of uB */
  uBdata = N_VGetArrayPointer_Parallel(uB);
  my_length = N_VGetLocalLength_Parallel(uB);

  /* Set adjoint states to 1.0 and quadrature variables to 0.0 */
  if (my_base == -1) for (i=0; i<my_length; i++) uBdata[i] = ZERO;
  else               for (i=0; i<my_length; i++) uBdata[i] = ONE;
}

/*
 * Compute local value of the space integral int_x z(x) dx 
 */

static realtype Xintgr(realtype *z, sunindextype l, realtype dx)
{
  realtype my_intgr;
  sunindextype i;

  my_intgr = RCONST(0.5)*(z[0] + z[l-1]);
  for (i = 1; i < l-1; i++)
    my_intgr += z[i]; 
  my_intgr *= dx;

  return(my_intgr);
}

/*
 * Compute value of g(u) 
 */

static realtype Compute_g(N_Vector u, UserData data)
{
  realtype intgr, my_intgr, dx, *udata;
  sunindextype my_length;
  int npes, my_pe, i;
  MPI_Status status;
  MPI_Comm comm;

  /* Extract MPI info. from data */
  comm = data->comm;
  npes = data->npes;
  my_pe = data->my_pe;

  dx = data->dx;

  if (my_pe == npes) {  /* Loop over all other processes and sum */
    intgr = ZERO;
    for (i=0; i<npes; i++) {
      MPI_Recv(&my_intgr, 1, PVEC_REAL_MPI_TYPE, i, 0, comm, &status); 
      intgr += my_intgr;
    }
    return(intgr);
  } else {              /* Compute local portion of the integral */
    udata = N_VGetArrayPointer_Parallel(u);
    my_length = N_VGetLocalLength_Parallel(u);
    my_intgr = Xintgr(udata, my_length, dx);
    MPI_Send(&my_intgr, 1, PVEC_REAL_MPI_TYPE, npes, 0, comm);
    return(my_intgr);
  }
}

/* 
 * Print output after backward integration
 */

static void PrintOutput(realtype g_val, N_Vector uB, UserData data)
{
  MPI_Comm comm;
  MPI_Status status;
  int npes, my_pe;
  sunindextype i, Ni, indx, local_N, nperpe, nrem;
  realtype *uBdata;
  realtype *mu;

  comm = data->comm;
  npes = data->npes;
  my_pe = data->my_pe;
  local_N = data->local_N;
  nperpe = data->nperpe;
  nrem = data->nrem;

  uBdata = N_VGetArrayPointer_Parallel(uB);

  if (my_pe == npes) {

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("\ng(tf) = %8Le\n\n", g_val);
    printf("dgdp(tf)\n  [ 1]: %8Le\n  [ 2]: %8Le\n\n", -uBdata[0], -uBdata[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("\ng(tf) = %8e\n\n", g_val);
    printf("dgdp(tf)\n  [ 1]: %8e\n  [ 2]: %8e\n\n", -uBdata[0], -uBdata[1]);
#else
    printf("\ng(tf) = %8e\n\n", g_val);
    printf("dgdp(tf)\n  [ 1]: %8e\n  [ 2]: %8e\n\n", -uBdata[0], -uBdata[1]);
#endif

    mu = (realtype *)malloc(NEQ*sizeof(realtype));
    if (check_retval((void *)mu, "malloc", 2, my_pe)) MPI_Abort(comm, 1);

    indx = 0;
    for ( i = 0; i < npes; i++) {
      Ni = ( i < nrem ) ? nperpe+1 : nperpe;
      MPI_Recv(&mu[indx], Ni, PVEC_REAL_MPI_TYPE, i, 0, comm, &status);
      indx += Ni;
    }

    printf("mu(t0)\n");

#if defined(SUNDIALS_EXTENDED_PRECISION)
    for (i=0; i<NEQ; i++)
      printf("  [%2ld]: %8Le\n", (long int) i+1, mu[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    for (i=0; i<NEQ; i++)
      printf("  [%2ld]: %8e\n", (long int) i+1, mu[i]);
#else
    for (i=0; i<NEQ; i++)
      printf("  [%2ld]: %8e\n", (long int) i+1, mu[i]);
#endif

    free(mu);

  } else {

    MPI_Send(uBdata, local_N, PVEC_REAL_MPI_TYPE, npes, 0, comm);

  }

}

/* 
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns an integer value so check if
 *             retval < 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
 */

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
