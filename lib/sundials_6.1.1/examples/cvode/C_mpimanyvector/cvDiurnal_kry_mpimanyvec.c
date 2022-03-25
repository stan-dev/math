/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * Example problem (based on cvDiurnal_kry_p.c):
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(y) = Kv0*exp(y/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is posed on the square
 *   0 <= x <= 20,    30 <= y <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * mesh, with simple polynomial initial profiles.
 *
 * The problem is solved by CVODE on NPE processors, treated
 * as a rectangular process grid of size NPEX by NPEY, with
 * NPE = NPEX*NPEY. Each processor contains a subgrid of size MXSUB
 * by MYSUB of the (x,y) mesh.  Thus the actual mesh sizes are
 * MX = MXSUB*NPEX and MY = MYSUB*NPEY, and the ODE system size is
 * neq = 2*MX*MY.
 *
 * Each species is stored in its own parallel nvector, and are
 * combined together using the MPIManyVector module.
 *
 * The solution is done with the BDF/GMRES method (i.e. using the
 * SUNLinSol_SPGMR linear solver) and the block-diagonal part of the
 * Newton matrix as a left preconditioner. A copy of the
 * block-diagonal part of the Jacobian is saved and conditionally
 * reused within the preconditioner routine.
 *
 * Performance data and sampled solution values are printed at
 * selected output times, and all performance counters are printed
 * on completion.
 *
 * This version uses MPI for user routines.
 *
 * Execution: mpirun -np N cvDiurnal_kry_mpimanyvec
 * with N = NPEX*NPEY (see constants below).
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvode/cvode.h>                    /* prototypes for CVODE fcts., consts.      */
#include <nvector/nvector_mpimanyvector.h>  /* access to MPIManyVector (includes mpi.h) */
#include <nvector/nvector_parallel.h>       /* access to MPI-parallel N_Vector          */
#include <sunlinsol/sunlinsol_spgmr.h>      /* access to SPGMR SUNLinearSolver          */
#include <sundials/sundials_dense.h>        /* prototypes for small dense fcts.         */
#include <sundials/sundials_types.h>        /* definitions of realtype, booleantype     */

/* helpful macros */

#ifndef SQR
#define SQR(A) ((A)*(A))
#endif

/* Problem Constants */

#define NVARS        2                    /* number of species         */
#define KH           RCONST(4.0e-6)       /* horizontal diffusivity Kh */
#define VEL          RCONST(0.001)        /* advection velocity V      */
#define KV0          RCONST(1.0e-8)       /* coefficient in Kv(y)      */
#define Q1           RCONST(1.63e-16)     /* coefficients q1, q2, c3   */
#define Q2           RCONST(4.66e-16)
#define C3           RCONST(3.7e16)
#define A3           RCONST(22.62)     /* coefficient in expression for q3(t) */
#define A4           RCONST(7.601)     /* coefficient in expression for q4(t) */
#define C1_SCALE     RCONST(1.0e6)     /* coefficients in initial profiles    */
#define C2_SCALE     RCONST(1.0e12)

#define T0           RCONST(0.0)          /* initial time */
#define NOUT         12                   /* number of output times */
#define TWOHR        RCONST(7200.0)       /* number of seconds in two hours  */
#define HALFDAY      RCONST(4.32e4)       /* number of seconds in a half day */
#define PI       RCONST(3.1415926535898)  /* pi */

#define XMIN         RCONST(0.0)          /* grid boundaries in x  */
#define XMAX         RCONST(20.0)
#define YMIN         RCONST(30.0)         /* grid boundaries in y  */
#define YMAX         RCONST(50.0)

#define NPEX         2              /* no. PEs in x direction of PE array */
#define NPEY         2              /* no. PEs in y direction of PE array */
                                    /* Total no. PEs = NPEX*NPEY */
#define MXSUB        5              /* no. x points per subgrid */
#define MYSUB        5              /* no. y points per subgrid */

#define MX           (NPEX*MXSUB)   /* MX = number of x mesh points */
#define MY           (NPEY*MYSUB)   /* MY = number of y mesh points */
                                    /* Spatial mesh is MX by MY */
/* CVodeInit Constants */

#define RTOL    RCONST(1.0e-5)    /* scalar relative tolerance */
#define FLOOR   RCONST(100.0)     /* value of C1 or C2 at which tolerances */
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */


/* User-defined matrix accessor macro: IJth */

/* IJth is defined in order to write code which indexes into dense
   matrices with a (row,column) pair, where 1 <= row,column <= NVARS.

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NVARS. The small matrix routines in sundials_dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. */

#define IJth(a,i,j) (a[j-1][i-1])

/* Type : UserData
   contains problem constants, preconditioner blocks, pivot arrays,
   grid constants, and processor indices, as well as data needed
   for the preconditiner */

typedef struct {

  realtype q4, om, dx, dy, hdco, haco, vdco;
  realtype c1ext[(MXSUB+2)*(MYSUB+2)];
  realtype c2ext[(MXSUB+2)*(MYSUB+2)];
  realtype SendBufferE[NVARS*MYSUB];
  realtype SendBufferW[NVARS*MYSUB];
  realtype SendBufferN[NVARS*MXSUB];
  realtype SendBufferS[NVARS*MXSUB];
  realtype RecvBufferE[NVARS*MYSUB];
  realtype RecvBufferW[NVARS*MYSUB];
  realtype RecvBufferN[NVARS*MXSUB];
  realtype RecvBufferS[NVARS*MXSUB];
  int my_pe, isubx, isuby;
  int nvmxsub2;
  MPI_Comm comm;
  MPI_Request request[8];

  /* For preconditioner */
  realtype **P[MXSUB][MYSUB], **Jbd[MXSUB][MYSUB];
  sunindextype *pivot[MXSUB][MYSUB];

} *UserData;

/* Private Helper Functions */

static void InitUserData(int my_pe, MPI_Comm comm, UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector u, UserData data);
static void PrintOutput(void *cvode_mem, int my_pe, MPI_Comm comm,
                        N_Vector u, realtype t);
static void PrintFinalStats(void *cvode_mem);
static void BSend(realtype c1data[], realtype c2data[], UserData data);
static void BRecvPost(UserData data);
static void BRecvWait(UserData data);
static void PrepareExt(N_Vector u, UserData data);
static void fcalc(realtype t, N_Vector udot, UserData data);


/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static int Precond(realtype tn, N_Vector u, N_Vector fu,
                   booleantype jok, booleantype *jcurPtr,
                   realtype gamma, void *user_data);

static int PSolve(realtype tn, N_Vector u, N_Vector fu,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *user_data);


/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  SUNContext sunctx;
  realtype abstol, reltol, t, tout;
  N_Vector u, c[2];
  UserData data;
  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, retval, my_pe, npes;
  sunindextype neq, local_N;
  MPI_Comm comm;

  u = NULL;
  c[0] = c[1] = NULL;
  data = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Set per-species problem size, neq */
  neq = MX*MY;

  /* Get processor number and total number of pe's */
  retval = MPI_Init(&argc, &argv);
  if (retval != MPI_SUCCESS)  return(1);
  comm = MPI_COMM_WORLD;
  retval = MPI_Comm_size(comm, &npes);
  if (retval != MPI_SUCCESS)  MPI_Abort(comm,1);
  retval = MPI_Comm_rank(comm, &my_pe);
  if (retval != MPI_SUCCESS)  MPI_Abort(comm,1);

  if (npes != NPEX*NPEY) {
    if (my_pe == 0)
      fprintf(stderr, "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n\n",
	      npes,NPEX*NPEY);
    MPI_Finalize();
    MPI_Abort(comm, 1);
  }

  /* Allocate and load user data block; allocate preconditioner block */
  data = (UserData) malloc(sizeof *data);
  if (check_retval((void *)data, "malloc", 2, my_pe)) MPI_Abort(comm, 1);
  InitUserData(my_pe, comm, data);

  /* Set local length */
  local_N = MXSUB*MYSUB;

  /* Create the SUNDIALS context */
  retval = SUNContext_Create(&comm, &sunctx);
  if(check_retval(&retval, "SUNContext_Create", 1, my_pe)) MPI_Abort(comm, 1);

  /* Allocate c[0], c[1], u, and set initial values and tolerances */
  c[0] = N_VNew_Parallel(comm, local_N, neq, sunctx);
  if (check_retval((void *)c[0], "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  c[1] = N_VNew_Parallel(comm, local_N, neq, sunctx);
  if (check_retval((void *)c[1], "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  u = N_VNew_MPIManyVector(2, c, sunctx);
  if (check_retval((void *)u, "N_VNew_MPIManyVector", 0, my_pe)) MPI_Abort(comm, 1);
  SetInitialProfiles(u, data);
  abstol = ATOL; reltol = RTOL;

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0, my_pe)) MPI_Abort(comm, 1);

  /* Set the pointer to user-defined data */
  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval(&retval, "CVodeSetUserData", 1, my_pe)) MPI_Abort(comm, 1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1, my_pe)) MPI_Abort(comm, 1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1, my_pe)) MPI_Abort(comm, 1);

  /* Create SPGMR solver structure with left preconditioning
     and the default Krylov dimension maxl */
  LS = SUNLinSol_SPGMR(u, SUN_PREC_LEFT, 0, sunctx);
  if (check_retval((void *)LS, "SUNLinSol_SPGMR", 0, my_pe)) MPI_Abort(comm, 1);

  /* Attach SPGMR solver structure to CVode interface */
  retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1, my_pe)) MPI_Abort(comm, 1);

  /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
  retval = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
  if (check_retval(&retval, "CVodeSetPreconditioner", 1, my_pe)) MPI_Abort(comm, 1);

  if (my_pe == 0)
    printf("\n2-species diurnal advection-diffusion problem\n\n");

  PrintOutput(cvode_mem, my_pe, comm, u, T0);

  /* In loop over output points, call CVode, print results, test for error */
  for (iout=1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1, my_pe)) break;
    PrintOutput(cvode_mem, my_pe, comm, u, t);
  }

  /* Print final statistics */
  if (my_pe == 0) PrintFinalStats(cvode_mem);

  /* Free memory */
  N_VDestroy(u);
  N_VDestroy(c[0]);
  N_VDestroy(c[1]);
  FreeUserData(data);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNContext_Free(&sunctx);

  MPI_Finalize();

  return(0);
}


/*********************** Private Helper Functions ************************/


/* Load constants in data */

static void InitUserData(int my_pe, MPI_Comm comm, UserData data)
{
  int lx, ly;

  /* Set problem constants */
  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/((realtype)(MX-1));
  data->dy = (YMAX-YMIN)/((realtype)(MY-1));
  data->hdco = KH/SQR(data->dx);
  data->haco = VEL/(RCONST(2.0)*data->dx);
  data->vdco = (RCONST(1.0)/SQR(data->dy))*KV0;

  /* Set machine-related constants */
  data->comm = comm;
  data->my_pe = my_pe;

  /* isubx and isuby are the PE grid indices corresponding to my_pe */
  data->isuby = my_pe/NPEX;
  data->isubx = my_pe - (data->isuby)*NPEX;

  /* Set the sizes of a boundary x-line in u and uext */
  data->nvmxsub2 = MXSUB+2;

  /* Preconditioner-related fields */
  for (lx=0; lx<MXSUB; lx++) {
    for (ly=0; ly<MYSUB; ly++) {
      (data->P)[lx][ly] = SUNDlsMat_newDenseMat(NVARS, NVARS);
      (data->Jbd)[lx][ly] = SUNDlsMat_newDenseMat(NVARS, NVARS);
      (data->pivot)[lx][ly] = SUNDlsMat_newIndexArray(NVARS);
    }
  }
}

/* Free user data memory */

static void FreeUserData(UserData data)
{
  int lx, ly;

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      SUNDlsMat_destroyMat((data->P)[lx][ly]);
      SUNDlsMat_destroyMat((data->Jbd)[lx][ly]);
      SUNDlsMat_destroyArray((data->pivot)[lx][ly]);
    }
  }

  free(data);
}

/* Set initial conditions in u */

static void SetInitialProfiles(N_Vector u, UserData data)
{
  int lx, ly, jx, jy;
  sunindextype offset;
  realtype x, y, cx, cy, xmid, ymid;
  realtype *c1data, *c2data;

  /* Set pointer to data array in vector u */
  c1data = N_VGetSubvectorArrayPointer_MPIManyVector(u,0);
  c2data = N_VGetSubvectorArrayPointer_MPIManyVector(u,1);

  /* Load initial profiles of c1 and c2 into local u vector.
  Here lx and ly are local mesh point indices on the local subgrid,
  and jx and jy are the global mesh point indices. */
  offset = 0;
  xmid = RCONST(0.5)*(XMIN + XMAX);
  ymid = RCONST(0.5)*(YMIN + YMAX);
  for (ly = 0; ly < MYSUB; ly++) {
    jy = ly + (data->isuby)*MYSUB;
    y = YMIN + jy*(data->dy);
    cy = SQR(RCONST(0.1)*(y - ymid));
    cy = RCONST(1.0) - cy + RCONST(0.5)*SQR(cy);
    for (lx = 0; lx < MXSUB; lx++) {
      jx = lx + (data->isubx)*MXSUB;
      x = XMIN + jx*(data->dx);
      cx = SQR(RCONST(0.1)*(x - xmid));
      cx = RCONST(1.0) - cx + RCONST(0.5)*SQR(cx);
      c1data[offset] = C1_SCALE*cx*cy;
      c2data[offset] = C2_SCALE*cx*cy;
      offset++;
    }
  }
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(void *cvode_mem, int my_pe, MPI_Comm comm,
                        N_Vector u, realtype t)
{
  int qu, retval;
  realtype hu, *c1data, *c2data, tempu[2];
  int npelast;
  long int nst;
  MPI_Status status;

  npelast = NPEX*NPEY - 1;
  c1data = N_VGetSubvectorArrayPointer_MPIManyVector(u,0);
  c2data = N_VGetSubvectorArrayPointer_MPIManyVector(u,1);

  /* Send c1,c2 at top right mesh point to PE 0 */
  if (my_pe == npelast) {
    tempu[0] = c1data[MXSUB*MYSUB-1];
    tempu[1] = c2data[MXSUB*MYSUB-1];
    if (npelast != 0) {
      retval = MPI_Send(&(tempu[0]), 2, MPI_SUNREALTYPE, 0, 0, comm);
      if (check_retval(&retval, "MPI_Send", 3, my_pe))  MPI_Abort(comm,1);
    }
  }

  /* On PE 0, receive c1,c2 at top right, then print performance data
     and sampled solution values */
  if (my_pe == 0) {
    if (npelast != 0) {
      retval = MPI_Recv(&tempu[0], 2, MPI_SUNREALTYPE, npelast, 0, comm, &status);
      if (check_retval(&retval, "MPI_Recv", 3, my_pe))  MPI_Abort(comm,1);
    }
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_retval(&retval, "CVodeGetNumSteps", 1, my_pe);
    retval = CVodeGetLastOrder(cvode_mem, &qu);
    check_retval(&retval, "CVodeGetLastOrder", 1, my_pe);
    retval = CVodeGetLastStep(cvode_mem, &hu);
    check_retval(&retval, "CVodeGetLastStep", 1, my_pe);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("t = %.2Le   no. steps = %ld   order = %d   stepsize = %.2Le\n",
           t, nst, qu, hu);
    printf("At bottom left:  c1, c2 = %12.3Le %12.3Le \n", c1data[0], c2data[0]);
    printf("At top right:    c1, c2 = %12.3Le %12.3Le \n\n", tempu[0], tempu[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("t = %.2e   no. steps = %ld   order = %d   stepsize = %.2e\n",
           t, nst, qu, hu);
    printf("At bottom left:  c1, c2 = %12.3e %12.3e \n", c1data[0], c2data[0]);
    printf("At top right:    c1, c2 = %12.3e %12.3e \n\n", tempu[0], tempu[1]);
#else
    printf("t = %.2e   no. steps = %ld   order = %d   stepsize = %.2e\n",
           t, nst, qu, hu);
    printf("At bottom left:  c1, c2 = %12.3e %12.3e \n", c1data[0], c2data[0]);
    printf("At top right:    c1, c2 = %12.3e %12.3e \n\n", tempu[0], tempu[1]);
#endif
  }
}

/* Print final statistics contained in iopt */

static void PrintFinalStats(void *cvode_mem)
{
  long int lenrw, leniw;
  long int lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int retval;

  retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_retval(&retval, "CVodeGetWorkSpace", 1, 0);
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1, 0);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1, 0);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1, 0);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1, 0);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1, 0);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1, 0);

  retval = CVodeGetLinWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
  check_retval(&retval, "CVodeGetLinWorkSpace", 1, 0);
  retval = CVodeGetNumLinIters(cvode_mem, &nli);
  check_retval(&retval, "CVodeGetNumLinIters", 1, 0);
  retval = CVodeGetNumPrecEvals(cvode_mem, &npe);
  check_retval(&retval, "CVodeGetNumPrecEvals", 1, 0);
  retval = CVodeGetNumPrecSolves(cvode_mem, &nps);
  check_retval(&retval, "CVodeGetNumPrecSolves", 1, 0);
  retval = CVodeGetNumLinConvFails(cvode_mem, &ncfl);
  check_retval(&retval, "CVodeGetNumLinConvFails", 1, 0);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1, 0);

  printf("\nFinal Statistics: \n\n");
  printf("lenrw   = %5ld     leniw   = %5ld\n", lenrw, leniw);
  printf("lenrwls = %5ld     leniwls = %5ld\n", lenrwLS, leniwLS);
  printf("nst     = %5ld\n"                  , nst);
  printf("nfe     = %5ld     nfels   = %5ld\n"  , nfe, nfeLS);
  printf("nni     = %5ld     nli     = %5ld\n"  , nni, nli);
  printf("nsetups = %5ld     netf    = %5ld\n"  , nsetups, netf);
  printf("npe     = %5ld     nps     = %5ld\n"  , npe, nps);
  printf("ncfn    = %5ld     ncfl    = %5ld\n\n", ncfn, ncfl);
}

/* Routine to start receiving boundary data from neighboring PEs. */

static void BRecvPost(UserData data)
{
  int retval;

  /* If isuby > 0, receive data for bottom x-line */
  if (data->isuby > 0) {
    retval = MPI_Irecv(data->RecvBufferS, NVARS*MXSUB, MPI_SUNREALTYPE,
                       data->my_pe-NPEX, 0, data->comm, &(data->request[0]));
    if (check_retval(&retval, "MPI_Irecv", 3, data->my_pe))  MPI_Abort(data->comm,1);
  }

  /* If isuby < NPEY-1, receive data for top x-line */
  if (data->isuby < NPEY-1) {
    retval = MPI_Irecv(data->RecvBufferN, NVARS*MXSUB, MPI_SUNREALTYPE,
                       data->my_pe+NPEX, 1, data->comm, &(data->request[1]));
    if (check_retval(&retval, "MPI_Irecv", 3, data->my_pe))  MPI_Abort(data->comm,1);
  }

  /* If isubx > 0, receive data for left y-line */
  if (data->isubx > 0) {
    retval = MPI_Irecv(data->RecvBufferW, NVARS*MYSUB, MPI_SUNREALTYPE,
                       data->my_pe-1, 2, data->comm, &(data->request[2]));
    if (check_retval(&retval, "MPI_Irecv", 3, data->my_pe))  MPI_Abort(data->comm,1);
  }

  /* If isubx < NPEX-1, receive data for right y-line */
  if (data->isubx < NPEX-1) {
    retval = MPI_Irecv(data->RecvBufferE, NVARS*MYSUB, MPI_SUNREALTYPE,
                       data->my_pe+1, 3, data->comm, &(data->request[3]));
    if (check_retval(&retval, "MPI_Irecv", 3, data->my_pe))  MPI_Abort(data->comm,1);
  }
}

/* Routine to send boundary data to neighboring PEs */

static void BSend(realtype c1data[], realtype c2data[], UserData data)
{
  int ly, lx, retval;
  sunindextype offsetu;

  /* If isuby > 0, send data from bottom x-line of c1 and c2 */
  if (data->isuby > 0) {
    for (lx=0; lx<MXSUB; lx++)  data->SendBufferS[lx]       = c1data[lx];
    for (lx=0; lx<MXSUB; lx++)  data->SendBufferS[MXSUB+lx] = c2data[lx];

    retval = MPI_Isend(data->SendBufferS, NVARS*MXSUB, MPI_SUNREALTYPE,
                       data->my_pe-NPEX, 1, data->comm, &(data->request[4]));
    if (check_retval(&retval, "MPI_Send", 3, data->my_pe))  MPI_Abort(data->comm,1);
  }

  /* If isuby < NPEY-1, send data from top x-line of c1 and c2 */
  if (data->isuby < NPEY-1) {
    offsetu = (MYSUB-1)*MXSUB;
    for (lx=0; lx<MXSUB; lx++)  data->SendBufferN[lx]       = c1data[offsetu+lx];
    for (lx=0; lx<MXSUB; lx++)  data->SendBufferN[MXSUB+lx] = c2data[offsetu+lx];

    retval = MPI_Isend(data->SendBufferN, NVARS*MXSUB, MPI_SUNREALTYPE,
                       data->my_pe+NPEX, 0, data->comm, &(data->request[5]));
    if (check_retval(&retval, "MPI_Send", 3, data->my_pe))  MPI_Abort(data->comm,1);
  }

  /* If isubx > 0, send data from left y-line of c1 and c2 */
  if (data->isubx > 0) {
    for (ly=0; ly<MYSUB; ly++)  data->SendBufferW[ly]       = c1data[ly*MXSUB];
    for (ly=0; ly<MYSUB; ly++)  data->SendBufferW[MYSUB+ly] = c2data[ly*MXSUB];

    retval = MPI_Isend(data->SendBufferW, NVARS*MYSUB, MPI_SUNREALTYPE,
                       data->my_pe-1, 3, data->comm, &(data->request[6]));
    if (check_retval(&retval, "MPI_Send", 3, data->my_pe))  MPI_Abort(data->comm,1);
  }

  /* If isubx < NPEX-1, send data from right y-line of c1 and c2 */
  if (data->isubx < NPEX-1) {
    for (ly=0; ly<MYSUB; ly++)  data->SendBufferE[ly]       = c1data[(ly+1)*MXSUB-1];
    for (ly=0; ly<MYSUB; ly++)  data->SendBufferE[MYSUB+ly] = c2data[(ly+1)*MXSUB-1];

    retval = MPI_Isend(data->SendBufferE, NVARS*MYSUB, MPI_SUNREALTYPE,
                       data->my_pe+1, 2, data->comm, &(data->request[7]));
    if (check_retval(&retval, "MPI_Send", 3, data->my_pe))  MPI_Abort(data->comm,1);
  }
}

/* Routine to finish receiving boundary data from neighboring PEs. */

static void BRecvWait(UserData data)
{
  int lx, ly, retval;
  MPI_Status status;

  /* If isuby > 0, wait on communication for bottom x-line */
  if (data->isuby > 0) {
    retval = MPI_Wait(&(data->request[0]), &status);
    if (check_retval(&retval, "MPI_Wait", 3, data->my_pe))  MPI_Abort(data->comm,1);
    retval = MPI_Wait(&(data->request[4]), &status);
    if (check_retval(&retval, "MPI_Wait", 3, data->my_pe))  MPI_Abort(data->comm,1);

    /* Copy the receive buffer to c1ext and c2ext */
    for (lx=0; lx<MXSUB; lx++)  data->c1ext[1+lx] = data->RecvBufferS[lx];
    for (lx=0; lx<MXSUB; lx++)  data->c2ext[1+lx] = data->RecvBufferS[MXSUB+lx];
  }

  /* If isuby < NPEY-1, wait on communication for top x-line */
  if (data->isuby < NPEY-1) {
    retval = MPI_Wait(&(data->request[1]), &status);
    if (check_retval(&retval, "MPI_Wait", 3, data->my_pe))  MPI_Abort(data->comm,1);
    retval = MPI_Wait(&(data->request[5]), &status);
    if (check_retval(&retval, "MPI_Wait", 3, data->my_pe))  MPI_Abort(data->comm,1);

    /* Copy the receive buffer to c1ext and c2ext */
    for (lx=0; lx<MXSUB; lx++)  data->c1ext[(MYSUB+1)*(MXSUB+2)+1+lx] = data->RecvBufferN[lx];
    for (lx=0; lx<MXSUB; lx++)  data->c2ext[(MYSUB+1)*(MXSUB+2)+1+lx] = data->RecvBufferN[MXSUB+lx];
  }

  /* If isubx > 0, wait on communication for left y-line */
  if (data->isubx > 0) {
    retval = MPI_Wait(&(data->request[2]), &status);
    if (check_retval(&retval, "MPI_Wait", 3, data->my_pe))  MPI_Abort(data->comm,1);
    retval = MPI_Wait(&(data->request[6]), &status);
    if (check_retval(&retval, "MPI_Wait", 3, data->my_pe))  MPI_Abort(data->comm,1);

    /* Copy the receive buffer to c1ext and c2ext */
    for (ly=0; ly<MYSUB; ly++)  data->c1ext[(ly+1)*(MXSUB+2)] = data->RecvBufferW[ly];
    for (ly=0; ly<MYSUB; ly++)  data->c2ext[(ly+1)*(MXSUB+2)] = data->RecvBufferW[MYSUB+ly];
  }

  /* If isubx < NPEX-1, wait on communication for right y-line */
  if (data->isubx < NPEX-1) {
    retval = MPI_Wait(&(data->request[3]), &status);
    if (check_retval(&retval, "MPI_Wait", 3, data->my_pe))  MPI_Abort(data->comm,1);
    retval = MPI_Wait(&(data->request[7]), &status);
    if (check_retval(&retval, "MPI_Wait", 3, data->my_pe))  MPI_Abort(data->comm,1);

    /* Copy the receive buffer to c1ext and c2ext */
    for (ly=0; ly<MYSUB; ly++)  data->c1ext[(ly+2)*(MXSUB+2)-1] = data->RecvBufferE[ly];
    for (ly=0; ly<MYSUB; ly++)  data->c2ext[(ly+2)*(MXSUB+2)-1] = data->RecvBufferE[MYSUB+ly];
  }
}

/* PrepareExt routine.
   This routine performs all communication between processors of data needed to calculate f.
   It then copies the data from u to the extended work arrays c1ext and c2ext.
   Then to facilitate homogeneous Neumann boundary conditions, this copies data from the
   first interior mesh lines of c1,c2 to c1ext,c2ext for boundary PEs. */

static void PrepareExt(N_Vector u, UserData data)
{
  realtype *c1data, *c2data;
  int lx, ly;

  /* Access data arrays from u, and extended work arrays c1ext and c2ext */
  c1data = N_VGetSubvectorArrayPointer_MPIManyVector(u,0);
  c2data = N_VGetSubvectorArrayPointer_MPIManyVector(u,1);

  /* Start receiving boundary data from neighboring PEs */
  BRecvPost(data);

  /* Send data from boundary of local grid to neighboring PEs */
  BSend(c1data, c2data, data);

  /* Copy local segments of c1 and c2 vectors into the working extended arrays c1ext and c2ext */
  for (ly=0; ly<MYSUB; ly++)
    for (lx=0; lx<MXSUB; lx++)
      data->c1ext[(MXSUB+2)*(1+ly)+1+lx] = c1data[ly*MXSUB+lx];
  for (ly=0; ly<MYSUB; ly++)
    for (lx=0; lx<MXSUB; lx++)
      data->c2ext[(MXSUB+2)*(1+ly)+1+lx] = c2data[ly*MXSUB+lx];

  /* Copy data from the first interior mesh line of c1,c2 to c1ext,c2ext for boundary PEs */

  /* If isuby = 0, copy x-line 2 of c1,c2 to c1ext,c2ext */
  if (data->isuby == 0) {
    for (lx=0; lx<MXSUB; lx++)
      data->c1ext[1+lx] = c1data[MXSUB+lx];
    for (lx=0; lx<MXSUB; lx++)
      data->c2ext[1+lx] = c2data[MXSUB+lx];
  }

  /* If isuby = NPEY-1, copy x-line MYSUB-1 of c1,c2 to c1ext,c2ext */
  if (data->isuby == NPEY-1) {
    for (lx=0; lx<MXSUB; lx++)
      data->c1ext[(MYSUB+1)*(MXSUB+2)+1+lx] = c1data[(MYSUB-2)*MXSUB+lx];
    for (lx=0; lx<MXSUB; lx++)
      data->c2ext[(MYSUB+1)*(MXSUB+2)+1+lx] = c2data[(MYSUB-2)*MXSUB+lx];
  }

  /* If isubx = 0, copy y-line 2 of u to uext */
  if (data->isubx == 0) {
    for (ly=0; ly<MYSUB; ly++)
      data->c1ext[(ly+1)*(MXSUB+2)] = c1data[ly*MXSUB+1];
    for (ly=0; ly<MYSUB; ly++)
      data->c2ext[(ly+1)*(MXSUB+2)] = c2data[ly*MXSUB+1];
  }

  /* If isubx = NPEX-1, copy y-line MXSUB-1 of u to uext */
  if (data->isubx == NPEX-1) {
    for (ly=0; ly<MYSUB; ly++)
      data->c1ext[(ly+2)*(MXSUB+2)-1] = c1data[(ly+1)*MXSUB-2];
    for (ly=0; ly<MYSUB; ly++)
      data->c2ext[(ly+2)*(MXSUB+2)-1] = c2data[(ly+1)*MXSUB-2];
  }

  /* Finish receiving boundary data from neighboring PEs */
  BRecvWait(data);

}


static void fcalc(realtype t, N_Vector udot, UserData data)
{
  realtype *c1dot, *c2dot;
  realtype q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  realtype c1rt, c2rt, cydn, cyup, hord1, hord2, horad1, horad2;
  realtype qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup;
  realtype q4coef;
  int lx, ly, jy;
  sunindextype offset;

  /* Access output arrays */
  c1dot = N_VGetSubvectorArrayPointer_MPIManyVector(udot,0);
  c2dot = N_VGetSubvectorArrayPointer_MPIManyVector(udot,1);

  /* Set diurnal rate coefficients as functions of t, and save q4 in
  data block for use by preconditioner evaluation routine */
  s = sin((data->om)*t);
  if (s > RCONST(0.0)) {
    q3 = exp(-A3/s);
    q4coef = exp(-A4/s);
  } else {
    q3 = RCONST(0.0);
    q4coef = RCONST(0.0);
  }
  data->q4 = q4coef;

  /* Loop over all grid points in local subgrid */
  for (ly=0; ly<MYSUB; ly++) {

    /* get global index of y grid coordinate */
    jy = ly + (data->isuby)*MYSUB;

    /* Set vertical diffusion coefficients at jy +- 1/2 */
    ydn = YMIN + (jy - RCONST(0.5))*(data->dy);
    yup = ydn + data->dy;
    cydn = (data->vdco)*exp(RCONST(0.2)*ydn);
    cyup = (data->vdco)*exp(RCONST(0.2)*yup);

    /* Loop over x-direction */
    for (lx=0; lx<MXSUB; lx++) {

      /* Extract c1, c2 over 5-point stencil */
      offset = (lx+1) + (ly+1)*(MXSUB+2);
      c1 = data->c1ext[offset];
      c2 = data->c2ext[offset];
      c1dn = data->c1ext[offset-(MXSUB+2)];
      c2dn = data->c2ext[offset-(MXSUB+2)];
      c1up = data->c1ext[offset+(MXSUB+2)];
      c2up = data->c2ext[offset+(MXSUB+2)];
      c1lt = data->c1ext[offset-1];
      c2lt = data->c2ext[offset-1];
      c1rt = data->c1ext[offset+1];
      c2rt = data->c2ext[offset+1];

      /* Set kinetic rate terms */
      qq1 = Q1*c1*C3;
      qq2 = Q2*c1*c2;
      qq3 = q3*C3;
      qq4 = q4coef*c2;
      rkin1 = -qq1 - qq2 + RCONST(2.0)*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms */
      vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
      vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);

      /* Set horizontal diffusion and advection terms */
      hord1 = (data->hdco)*(c1rt - RCONST(2.0)*c1 + c1lt);
      hord2 = (data->hdco)*(c2rt - RCONST(2.0)*c2 + c2lt);
      horad1 = (data->haco)*(c1rt - c1lt);
      horad2 = (data->haco)*(c2rt - c2lt);

      /* Load all terms into dudata */
      c1dot[lx+ly*MXSUB] = vertd1 + hord1 + horad1 + rkin1;
      c2dot[lx+ly*MXSUB] = vertd2 + hord2 + horad2 + rkin2;
    }
  }
}


/***************** Functions Called by the Solver *************************/

/* f routine.  Evaluate f(t,y).  First call PrepareExt to do communication of
   subgrid boundary data and to copy data from u into uext.
   Then calculate f by a call to fcalc. */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  /* Call PrepareExt to set up work arrays for calculation of RHS
     (includes inter-processor communication) */
  PrepareExt(u, (UserData) user_data);

  /* Call fcalc to calculate all right-hand sides */
  fcalc(t, udot, (UserData) user_data);

  return(0);
}

/* Preconditioner setup routine. Generate and preprocess P. */
static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, void *user_data)
{
  realtype c1, c2, cydn, cyup, diag, ydn, yup;
  realtype **(*P)[MYSUB], **(*Jbd)[MYSUB];
  sunindextype retval;
  sunindextype *(*pivot)[MYSUB];
  sunindextype lx, ly, jy;
  realtype *c1data, *c2data, **a, **j;
  UserData data;

  /* Make local copies of pointers in user_data, pointer to u's data,
     and PE index pair */
  data = (UserData) user_data;
  P = data->P;
  Jbd = data->Jbd;
  pivot = data->pivot;
  c1data = N_VGetSubvectorArrayPointer_MPIManyVector(u,0);
  c2data = N_VGetSubvectorArrayPointer_MPIManyVector(u,1);

  /* jok = SUNTRUE: Copy Jbd to P, and update jcurPtr */
  if (jok) {
    for (ly=0; ly<MYSUB; ly++)
      for (lx=0; lx<MXSUB; lx++)
        SUNDlsMat_denseCopy(Jbd[lx][ly], P[lx][ly], NVARS, NVARS);
    *jcurPtr = SUNFALSE;
  }

  /* jok = SUNFALSE: Generate Jbd from scratch and copy to P, and update jcurPtr */
  else {

    /* Compute 2x2 diagonal Jacobian blocks (using q4 values
       computed on the last f call).  Load into P. */
    for (ly=0; ly<MYSUB; ly++) {
      jy = ly + (data->isuby)*MYSUB;
      ydn = YMIN + (jy - RCONST(0.5))*(data->dy);
      yup = ydn + (data->dy);
      cydn = (data->vdco)*exp(RCONST(0.2)*ydn);
      cyup = (data->vdco)*exp(RCONST(0.2)*yup);
      diag = -(cydn + cyup + RCONST(2.0)*(data->hdco));
      for (lx = 0; lx < MXSUB; lx++) {
        c1 = c1data[lx+ly*MXSUB];
        c2 = c2data[lx+ly*MXSUB];
        j = Jbd[lx][ly];
        a = P[lx][ly];
        IJth(j,1,1) = (-Q1*C3 - Q2*c2) + diag;
        IJth(j,1,2) = -Q2*c1 + data->q4;
        IJth(j,2,1) = Q1*C3 - Q2*c2;
        IJth(j,2,2) = (-Q2*c1 - data->q4) + diag;
        SUNDlsMat_denseCopy(j, a, NVARS, NVARS);
      }
    }

    *jcurPtr = SUNTRUE;

  }

  /* Scale by -gamma */
  for (lx=0; lx<MXSUB; lx++)
    for (ly=0; ly<MYSUB; ly++)
      SUNDlsMat_denseScale(-gamma, P[lx][ly], NVARS, NVARS);

  /* Add identity matrix and do LU decompositions on blocks in place */
  for (lx=0; lx<MXSUB; lx++) {
    for (ly=0; ly<MYSUB; ly++) {
      SUNDlsMat_denseAddIdentity(P[lx][ly], NVARS);
      retval = SUNDlsMat_denseGETRF(P[lx][ly], NVARS, NVARS, pivot[lx][ly]);
      if (retval != 0) return(1);
    }
  }

  return(0);
}

/* Preconditioner solve routine */
static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  realtype **(*P)[MYSUB];
  sunindextype *(*pivot)[MYSUB];
  int lx, ly;
  realtype *z1data, *z2data, v[2];
  UserData data;

  /* Extract the P and pivot arrays from user_data */
  data = (UserData) user_data;
  P = data->P;
  pivot = data->pivot;

  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z.
     First copy vector r to z. */
  N_VScale(RCONST(1.0), r, z);
  z1data = N_VGetSubvectorArrayPointer_MPIManyVector(z,0);
  z2data = N_VGetSubvectorArrayPointer_MPIManyVector(z,1);
  for (lx=0; lx<MXSUB; lx++) {
    for (ly=0; ly<MYSUB; ly++) {
      v[0] = z1data[lx + ly*MXSUB];
      v[1] = z2data[lx + ly*MXSUB];
      SUNDlsMat_denseGETRS(P[lx][ly], NVARS, pivot[lx][ly], v);
    }
  }

  return(0);
}


/*********************** Private Helper Function ************************/

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval < 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer
     opt == 3 means MPI function, so check for MPI_SUCCESS */

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

  /* Check if retval != MPI_SUCCESS */
  else if (opt == 3) {
    retval = (int *) returnvalue;
    if (*retval != MPI_SUCCESS) {
    fprintf(stderr, "\nMPI_ERROR(%d): %s() failed with retval = %d\n\n",
	    id, funcname, *retval);
    return(1); }}

  return(0);
}
