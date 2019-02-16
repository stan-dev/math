/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds, Jean Sexton @ SMU
 *                Slaven Peles @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Example problem:
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
 * The problem is solved by ARKODE on NPE processors, treated
 * as a rectangular process grid of size NPEX by NPEY, with
 * NPE = NPEX*NPEY. Each processor contains a subgrid of size MXSUB
 * by MYSUB of the (x,y) mesh.  Thus the actual mesh sizes are
 * MX = MXSUB*NPEX and MY = MYSUB*NPEY, and the ODE system size is
 * neq = 2*MX*MY.
 *
 * The solution is done with the DIRK/GMRES method (i.e. using the
 * SUNLinSol_SPGMR linear solver) and the block-diagonal part of the
 * Newton matrix as a left preconditioner. A copy of the
 * block-diagonal part of the Jacobian is saved and conditionally
 * reused within the preconditioner routine.
 *
 * Performance data and sampled solution values are printed at
 * selected output times, and all performance counters are printed
 * on completion.
 *
 * This example uses Hypre vector and MPI parallelization. User is
 * expected to be familiar with the Hypre library.
 *
 * Execution: mpiexec -n N ark_diurnal_kry_ph  with N = NPEX*NPEY
 * (see constants below).
 *---------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode_arkstep.h>     /* prototypes for ARKStep fcts., consts */
#include <nvector/nvector_parhyp.h>    /* declaration of N_Vector  */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver  */
#include <sundials/sundials_dense.h>   /* prototypes for small dense fcts. */
#include <sundials/sundials_types.h>   /* definitions of realtype, booleantype */
#include <sundials/sundials_math.h>    /* definition of macros SUNSQR and EXP */
#include <mpi.h>                       /* MPI constants and types */

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>

/* Problem Constants */
#define NVARS        2                    /* number of species         */
#define KH           RCONST(4.0e-6)       /* horizontal diffusivity Kh */
#define VEL          RCONST(0.001)        /* advection velocity V      */
#define KV0          RCONST(1.0e-8)       /* coefficient in Kv(y)      */
#define Q1           RCONST(1.63e-16)     /* coefficients q1, q2, c3   */
#define Q2           RCONST(4.66e-16)
#define C3           RCONST(3.7e16)
#define A3           RCONST(22.62)        /* coefficient in expression for q3(t) */
#define A4           RCONST(7.601)        /* coefficient in expression for q4(t) */
#define C1_SCALE     RCONST(1.0e6)        /* coefficients in initial profiles    */
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

#define NPEX         2                    /* no. PEs in x direction of PE array */
#define NPEY         2                    /* no. PEs in y direction of PE array */
                                          /* Total no. PEs = NPEX*NPEY */
#define MXSUB        5                    /* no. x points per subgrid */
#define MYSUB        5                    /* no. y points per subgrid */

#define MX           (NPEX*MXSUB)         /* MX = number of x mesh points */
#define MY           (NPEY*MYSUB)         /* MY = number of y mesh points */
                                          /* Spatial mesh is MX by MY */

/* initialization constants */
#define RTOL    RCONST(1.0e-5)            /* scalar relative tolerance */
#define FLOOR   RCONST(100.0)             /* value of C1 or C2 at which tolerances */
                                          /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)              /* scalar absolute tolerance */


/* User-defined matrix accessor macro: IJth

   IJth is defined in order to write code which indexes into dense
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
  realtype uext[NVARS*(MXSUB+2)*(MYSUB+2)];
  int my_pe, isubx, isuby;
  int nvmxsub, nvmxsub2;
  MPI_Comm comm;

  /* For preconditioner */
  realtype **P[MXSUB][MYSUB], **Jbd[MXSUB][MYSUB];
  sunindextype *pivot[MXSUB][MYSUB];

} *UserData;

/* Private Helper Functions */
static void InitUserData(int my_pe, MPI_Comm comm, UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(HYPRE_IJVector Uij, UserData data,
                               sunindextype local_length,
                               sunindextype my_base);
static void PrintOutput(void *arkode_mem, int my_pe, MPI_Comm comm,
                        N_Vector u, realtype t);
static void PrintFinalStats(void *arkode_mem);
static void BSend(MPI_Comm comm,
                  int my_pe, int isubx, int isuby,
                  sunindextype dsizex, sunindextype dsizey,
                  realtype udata[]);
static void BRecvPost(MPI_Comm comm, MPI_Request request[],
                      int my_pe, int isubx, int isuby,
                      sunindextype dsizex, sunindextype dsizey,
                      realtype uext[], realtype buffer[]);
static void BRecvWait(MPI_Request request[],
                      int isubx, int isuby,
                      sunindextype dsizex, realtype uext[],
                      realtype buffer[]);
static void ucomm(realtype t, N_Vector u, UserData data);
static void fcalc(realtype t, realtype udata[], realtype dudata[],
                  UserData data);


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
static int check_flag(void *flagvalue, const char *funcname, int opt, int id);


/***************************** Main Program ******************************/
int main(int argc, char *argv[])
{
  realtype abstol, reltol, t, tout;
  N_Vector u;
  UserData data;
  SUNLinearSolver LS;
  void *arkode_mem;
  int iout, flag, my_pe, npes;
  MPI_Comm comm;

  HYPRE_Int local_N;
  HYPRE_ParVector Upar; /* Declare HYPRE parallel vector */
  HYPRE_IJVector  Uij;  /* Declare "IJ" interface to HYPRE vector */

  u = NULL;
  data = NULL;
  LS = NULL;
  arkode_mem = NULL;

  /* Get processor number and total number of pe's */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  if (npes != NPEX*NPEY) {
    if (my_pe == 0)
      fprintf(stderr, "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n\n",
              npes,NPEX*NPEY);
    MPI_Finalize();
    return(1);
  }

  /* Set local length */
  local_N = NVARS*MXSUB*MYSUB;

  /* Allocate hypre vector */
  HYPRE_IJVectorCreate(comm, my_pe*local_N, (my_pe + 1)*local_N - 1, &Uij);
  HYPRE_IJVectorSetObjectType(Uij, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(Uij);

  /* Allocate and load user data block; allocate preconditioner block */
  data = (UserData) malloc(sizeof *data);
  if (check_flag((void *)data, "malloc", 2, my_pe)) MPI_Abort(comm, 1);
  InitUserData(my_pe, comm, data);

  /* Set initial values and allocate u */
  SetInitialProfiles(Uij, data, local_N, my_pe*local_N);
  HYPRE_IJVectorAssemble(Uij);
  HYPRE_IJVectorGetObject(Uij, (void**) &Upar);

  u = N_VMake_ParHyp(Upar);  /* Create wrapper u around hypre vector */
  if (check_flag((void *)u, "N_VNew", 0, my_pe)) MPI_Abort(comm, 1);

  /* Set tolerances */
  abstol = ATOL; reltol = RTOL;

  /* Create SPGMR solver structure -- use left preconditioning
     and the default Krylov dimension maxl */
  LS = SUNLinSol_SPGMR(u, PREC_LEFT, 0);
  if (check_flag((void *)LS, "SUNLinSol_SPGMR", 0, my_pe)) MPI_Abort(comm, 1);

  /* Call ARKStepCreate to initialize the integrator memory and specify the
     user's right hand side function in u'=fi(t,u) [here fe is NULL],
     the inital time T0, and the initial dependent variable vector u. */
  arkode_mem = ARKStepCreate(NULL, f, T0, u);
  if (check_flag((void *)arkode_mem, "ARKStepCreate", 0, my_pe)) MPI_Abort(comm, 1);

  /* Set the pointer to user-defined data */
  flag = ARKStepSetUserData(arkode_mem, data);
  if (check_flag(&flag, "ARKStepSetUserData", 1, my_pe)) MPI_Abort(comm, 1);

  /* Call ARKStepSetMaxNumSteps to increase default */
  flag = ARKStepSetMaxNumSteps(arkode_mem, 10000);
  if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1, my_pe)) return(1);

  /* Call ARKStepSStolerances to specify the scalar relative tolerance
     and scalar absolute tolerances */
  flag = ARKStepSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKStepSStolerances", 1, my_pe)) return(1);

  /* Attach SPGMR solver structure to ARKStep interface */
  flag = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
  if (check_flag(&flag, "ARKStepSetLinearSolver", 1, my_pe)) MPI_Abort(comm, 1);

  /* Set preconditioner setup and solve routines Precond and PSolve,
     and the pointer to the user-defined block data */
  flag = ARKStepSetPreconditioner(arkode_mem, Precond, PSolve);
  if (check_flag(&flag, "ARKStepSetPreconditioner", 1, my_pe)) MPI_Abort(comm, 1);

  /* Print heading */
  if (my_pe == 0)
    printf("\n2-species diurnal advection-diffusion problem\n\n");

  /* In loop over output points, call ARKStepEvolve, print results, test for error */
  for (iout=1, tout=TWOHR; iout<=NOUT; iout++, tout+=TWOHR) {
    flag = ARKStepEvolve(arkode_mem, tout, u, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKStepEvolve", 1, my_pe)) break;
    PrintOutput(arkode_mem, my_pe, comm, u, t);
  }

  /* Print final statistics */
  if (my_pe == 0) PrintFinalStats(arkode_mem);

  /* Free memory */
  N_VDestroy(u);              /* Free hypre vector wrapper */
  HYPRE_IJVectorDestroy(Uij); /* Free the underlying hypre vector */
  FreeUserData(data);
  ARKStepFree(&arkode_mem);
  SUNLinSolFree(LS);
  MPI_Finalize();
  return(0);
}

/*********************** Private Helper Functions ************************/

/* Load constants in data */
static void InitUserData(int my_pe, MPI_Comm comm, UserData data)
{
  int isubx, isuby;
  int lx, ly;

  /* Set problem constants */
  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/((realtype)(MX-1));
  data->dy = (YMAX-YMIN)/((realtype)(MY-1));
  data->hdco = KH/SUNSQR(data->dx);
  data->haco = VEL/(RCONST(2.0)*data->dx);
  data->vdco = (RCONST(1.0)/SUNSQR(data->dy))*KV0;

  /* Set machine-related constants */
  data->comm = comm;
  data->my_pe = my_pe;

  /* isubx and isuby are the PE grid indices corresponding to my_pe */
  isuby = my_pe/NPEX;
  isubx = my_pe - isuby*NPEX;
  data->isubx = isubx;
  data->isuby = isuby;

  /* Set the sizes of a boundary x-line in u and uext */
  data->nvmxsub = NVARS*MXSUB;
  data->nvmxsub2 = NVARS*(MXSUB+2);

  /* Preconditioner-related fields */
  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      (data->P)[lx][ly] = newDenseMat(NVARS, NVARS);
      (data->Jbd)[lx][ly] = newDenseMat(NVARS, NVARS);
      (data->pivot)[lx][ly] = newIndexArray(NVARS);
    }
  }
}

/* Free user data memory */
static void FreeUserData(UserData data)
{
  int lx, ly;
  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      destroyMat((data->P)[lx][ly]);
      destroyMat((data->Jbd)[lx][ly]);
      destroyArray((data->pivot)[lx][ly]);
    }
  }
  free(data);
}

/* Set initial conditions in u */
static void SetInitialProfiles(HYPRE_IJVector Uij, UserData data,
                               sunindextype local_length,
                               sunindextype my_base)
{
  int isubx, isuby, lx, ly, jx, jy;
  sunindextype offset;
  realtype dx, dy, x, y, cx, cy, xmid, ymid;
  realtype *udata;
  HYPRE_Int *iglobal;

  /* Set pointer to data array in vector u */
  udata   = (realtype*) malloc(local_length*sizeof(realtype));
  iglobal = (HYPRE_Int*) malloc(local_length*sizeof(HYPRE_Int));


  /* Get mesh spacings, and subgrid indices for this PE */
  dx = data->dx;         dy = data->dy;
  isubx = data->isubx;   isuby = data->isuby;

  /* Load initial profiles of c1 and c2 into local u vector.
  Here lx and ly are local mesh point indices on the local subgrid,
  and jx and jy are the global mesh point indices. */
  offset = 0;
  xmid = RCONST(0.5)*(XMIN + XMAX);
  ymid = RCONST(0.5)*(YMIN + YMAX);
  for (ly = 0; ly < MYSUB; ly++) {
    jy = ly + isuby*MYSUB;
    y = YMIN + jy*dy;
    cy = SUNSQR(RCONST(0.1)*(y - ymid));
    cy = RCONST(1.0) - cy + RCONST(0.5)*SUNSQR(cy);
    for (lx = 0; lx < MXSUB; lx++) {
      jx = lx + isubx*MXSUB;
      x = XMIN + jx*dx;
      cx = SUNSQR(RCONST(0.1)*(x - xmid));
      cx = RCONST(1.0) - cx + RCONST(0.5)*SUNSQR(cx);
      iglobal[offset] = my_base + offset;
      udata[offset++] = C1_SCALE*cx*cy;
      iglobal[offset] = my_base + offset;
      udata[offset++] = C2_SCALE*cx*cy;
    }
  }
  HYPRE_IJVectorSetValues(Uij, local_length, iglobal, udata);
  free(iglobal);
  free(udata);
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */
static void PrintOutput(void *arkode_mem, int my_pe, MPI_Comm comm,
                        N_Vector u, realtype t)
{
  int flag;
  realtype hu, *udata, tempu[2];
  int npelast;
  sunindextype i0, i1;
  long int nst;
  MPI_Status status;
  HYPRE_ParVector uhyp;

  npelast = NPEX*NPEY - 1;

  uhyp  = N_VGetVector_ParHyp(u);
  udata = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));

  /* Send c1,c2 at top right mesh point to PE 0 */
  if (my_pe == npelast) {
    i0 = NVARS*MXSUB*MYSUB - 2;
    i1 = i0 + 1;
    if (npelast != 0)
      MPI_Send(&udata[i0], 2, PVEC_REAL_MPI_TYPE, 0, 0, comm);
    else {
      tempu[0] = udata[i0];
      tempu[1] = udata[i1];
    }
  }

  /* On PE 0, receive c1,c2 at top right, then print performance data
     and sampled solution values */
  if (my_pe == 0) {
    if (npelast != 0)
      MPI_Recv(&tempu[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);
    flag = ARKStepGetNumSteps(arkode_mem, &nst);
    check_flag(&flag, "ARKStepGetNumSteps", 1, my_pe);
    flag = ARKStepGetLastStep(arkode_mem, &hu);
    check_flag(&flag, "ARKStepGetLastStep", 1, my_pe);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("t = %.2Le   no. steps = %ld   stepsize = %.2Le\n",
           t, nst, hu);
    printf("At bottom left:  c1, c2 = %12.3Le %12.3Le \n", udata[0], udata[1]);
    printf("At top right:    c1, c2 = %12.3Le %12.3Le \n\n", tempu[0], tempu[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("t = %.2e   no. steps = %ld   stepsize = %.2e\n",
           t, nst, hu);
    printf("At bottom left:  c1, c2 = %12.3e %12.3e \n", udata[0], udata[1]);
    printf("At top right:    c1, c2 = %12.3e %12.3e \n\n", tempu[0], tempu[1]);
#else
    printf("t = %.2e   no. steps = %ld   stepsize = %.2e\n",
           t, nst, hu);
    printf("At bottom left:  c1, c2 = %12.3e %12.3e \n", udata[0], udata[1]);
    printf("At top right:    c1, c2 = %12.3e %12.3e \n\n", tempu[0], tempu[1]);
#endif
  }
}

/* Print final statistics contained in iopt */
static void PrintFinalStats(void *arkode_mem)
{
  long int lenrw, leniw;
  long int lenrwLS, leniwLS;
  long int nst, nfe, nfi, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int flag;

  flag = ARKStepGetWorkSpace(arkode_mem, &lenrw, &leniw);
  check_flag(&flag, "ARKStepGetWorkSpace", 1, 0);
  flag = ARKStepGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKStepGetNumSteps", 1, 0);
  flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKStepGetNumRhsEvals", 1, 0);
  flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1, 0);
  flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKStepGetNumErrTestFails", 1, 0);
  flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1, 0);
  flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1, 0);

  flag = ARKStepGetLinWorkSpace(arkode_mem, &lenrwLS, &leniwLS);
  check_flag(&flag, "ARKStepGetLinWorkSpace", 1, 0);
  flag = ARKStepGetNumLinIters(arkode_mem, &nli);
  check_flag(&flag, "ARKStepGetNumLinIters", 1, 0);
  flag = ARKStepGetNumPrecEvals(arkode_mem, &npe);
  check_flag(&flag, "ARKStepGetNumPrecEvals", 1, 0);
  flag = ARKStepGetNumPrecSolves(arkode_mem, &nps);
  check_flag(&flag, "ARKStepGetNumPrecSolves", 1, 0);
  flag = ARKStepGetNumLinConvFails(arkode_mem, &ncfl);
  check_flag(&flag, "ARKStepGetNumLinConvFails", 1, 0);
  flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1, 0);

  printf("\nFinal Statistics: \n\n");
  printf("lenrw   = %5ld     leniw   = %5ld\n", lenrw, leniw);
  printf("lenrwls = %5ld     leniwls = %5ld\n", lenrwLS, leniwLS);
  printf("nst     = %5ld     nfe     = %5ld\n", nst, nfe);
  printf("nfi     = %5ld     nfels   = %5ld\n", nfi, nfeLS);
  printf("nni     = %5ld     nli     = %5ld\n", nni, nli);
  printf("nsetups = %5ld     netf    = %5ld\n", nsetups, netf);
  printf("npe     = %5ld     nps     = %5ld\n", npe, nps);
  printf("ncfn    = %5ld     ncfl    = %5ld\n\n", ncfn, ncfl);
}

/* Routine to send boundary data to neighboring PEs */
static void BSend(MPI_Comm comm, int my_pe, int isubx,
                  int isuby, sunindextype dsizex,
                  sunindextype dsizey, realtype udata[])
{
  int i, ly;
  sunindextype offsetu, offsetbuf;
  realtype bufleft[NVARS*MYSUB], bufright[NVARS*MYSUB];

  /* If isuby > 0, send data from bottom x-line of u */
  if (isuby != 0)
    MPI_Send(&udata[0], dsizex, PVEC_REAL_MPI_TYPE, my_pe-NPEX, 0, comm);

  /* If isuby < NPEY-1, send data from top x-line of u */
  if (isuby != NPEY-1) {
    offsetu = (MYSUB-1)*dsizex;
    MPI_Send(&udata[offsetu], dsizex, PVEC_REAL_MPI_TYPE, my_pe+NPEX, 0, comm);
  }

  /* If isubx > 0, send data from left y-line of u (via bufleft) */
  if (isubx != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetu = ly*dsizex;
      for (i = 0; i < NVARS; i++)
        bufleft[offsetbuf+i] = udata[offsetu+i];
    }
    MPI_Send(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe-1, 0, comm);
  }

  /* If isubx < NPEX-1, send data from right y-line of u (via bufright) */
  if (isubx != NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetu = offsetbuf*MXSUB + (MXSUB-1)*NVARS;
      for (i = 0; i < NVARS; i++)
        bufright[offsetbuf+i] = udata[offsetu+i];
    }
    MPI_Send(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe+1, 0, comm);
  }
}

/* Routine to start receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */

static void BRecvPost(MPI_Comm comm, MPI_Request request[],
                      int my_pe, int isubx, int isuby,
                      sunindextype dsizex, sunindextype dsizey,
                      realtype uext[], realtype buffer[])
{
  sunindextype offsetue;
  /* Have bufleft and bufright use the same buffer */
  realtype *bufleft = buffer, *bufright = buffer+NVARS*MYSUB;

  /* If isuby > 0, receive data for bottom x-line of uext */
  if (isuby != 0)
    MPI_Irecv(&uext[NVARS], dsizex, PVEC_REAL_MPI_TYPE,
              my_pe-NPEX, 0, comm, &request[0]);

  /* If isuby < NPEY-1, receive data for top x-line of uext */
  if (isuby != NPEY-1) {
    offsetue = NVARS*(1 + (MYSUB+1)*(MXSUB+2));
    MPI_Irecv(&uext[offsetue], dsizex, PVEC_REAL_MPI_TYPE,
                                         my_pe+NPEX, 0, comm, &request[1]);
  }

  /* If isubx > 0, receive data for left y-line of uext (via bufleft) */
  if (isubx != 0) {
    MPI_Irecv(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE,
                                         my_pe-1, 0, comm, &request[2]);
  }

  /* If isubx < NPEX-1, receive data for right y-line of uext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Irecv(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE,
                                         my_pe+1, 0, comm, &request[3]);
  }
}

/* Routine to finish receiving boundary data from neighboring PEs.
   Notes:
   1) buffer should be able to hold 2*NVARS*MYSUB realtype entries, should be
   passed to both the BRecvPost and BRecvWait functions, and should not
   be manipulated between the two calls.
   2) request should have 4 entries, and should be passed in both calls also. */

static void BRecvWait(MPI_Request request[],
                      int isubx, int isuby,
                      sunindextype dsizex, realtype uext[],
                      realtype buffer[])
{
  int i, ly;
  sunindextype dsizex2, offsetue, offsetbuf;
  realtype *bufleft = buffer, *bufright = buffer+NVARS*MYSUB;
  MPI_Status status;

  dsizex2 = dsizex + 2*NVARS;

  /* If isuby > 0, receive data for bottom x-line of uext */
  if (isuby != 0)
    MPI_Wait(&request[0],&status);

  /* If isuby < NPEY-1, receive data for top x-line of uext */
  if (isuby != NPEY-1)
    MPI_Wait(&request[1],&status);

  /* If isubx > 0, receive data for left y-line of uext (via bufleft) */
  if (isubx != 0) {
    MPI_Wait(&request[2],&status);

    /* Copy the buffer to uext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetue = (ly+1)*dsizex2;
      for (i = 0; i < NVARS; i++)
        uext[offsetue+i] = bufleft[offsetbuf+i];
    }
  }

  /* If isubx < NPEX-1, receive data for right y-line of uext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Wait(&request[3],&status);

    /* Copy the buffer to uext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NVARS;
      offsetue = (ly+2)*dsizex2 - NVARS;
      for (i = 0; i < NVARS; i++)
        uext[offsetue+i] = bufright[offsetbuf+i];
    }
  }
}

/* ucomm routine.  This routine performs all communication
   between processors of data needed to calculate f. */

static void ucomm(realtype t, N_Vector u, UserData data)
{

  realtype *udata, *uext, buffer[2*NVARS*MYSUB];
  MPI_Comm comm;
  int my_pe, isubx, isuby;
  sunindextype nvmxsub, nvmysub;
  MPI_Request request[4];
  HYPRE_ParVector uhyp;

  uhyp  = N_VGetVector_ParHyp(u);
  udata = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));

  /* Get comm, my_pe, subgrid indices, data sizes, extended array uext */
  comm = data->comm;  my_pe = data->my_pe;
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub;
  nvmysub = NVARS*MYSUB;
  uext = data->uext;

  /* Start receiving boundary data from neighboring PEs */
  BRecvPost(comm, request, my_pe, isubx, isuby, nvmxsub, nvmysub, uext, buffer);

  /* Send data from boundary of local grid to neighboring PEs */
  BSend(comm, my_pe, isubx, isuby, nvmxsub, nvmysub, udata);

  /* Finish receiving boundary data from neighboring PEs */
  BRecvWait(request, isubx, isuby, nvmxsub, uext, buffer);
}


/* fcalc routine. Compute f(t,y).  This routine assumes that communication
   between processors of data needed to calculate f has already been done,
   and this data is in the work array uext. */

static void fcalc(realtype t, realtype udata[],
                  realtype dudata[], UserData data)
{
  realtype *uext;
  realtype q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  realtype c1rt, c2rt, cydn, cyup, hord1, hord2, horad1, horad2;
  realtype qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup;
  realtype q4coef, dely, verdco, hordco, horaco;
  int i, lx, ly, jy;
  int isubx, isuby;
  sunindextype nvmxsub, nvmxsub2, offsetu, offsetue;

  /* Get subgrid indices, data sizes, extended work array uext */
  isubx = data->isubx;   isuby = data->isuby;
  nvmxsub = data->nvmxsub; nvmxsub2 = data->nvmxsub2;
  uext = data->uext;

  /* Copy local segment of u vector into the working extended array uext */
  offsetu = 0;
  offsetue = nvmxsub2 + NVARS;
  for (ly = 0; ly < MYSUB; ly++) {
    for (i = 0; i < nvmxsub; i++) uext[offsetue+i] = udata[offsetu+i];
    offsetu = offsetu + nvmxsub;
    offsetue = offsetue + nvmxsub2;
  }

  /* To facilitate homogeneous Neumann boundary conditions, when this is
  a boundary PE, copy data from the first interior mesh line of u to uext */

  /* If isuby = 0, copy x-line 2 of u to uext */
  if (isuby == 0) {
    for (i = 0; i < nvmxsub; i++) uext[NVARS+i] = udata[nvmxsub+i];
  }

  /* If isuby = NPEY-1, copy x-line MYSUB-1 of u to uext */
  if (isuby == NPEY-1) {
    offsetu = (MYSUB-2)*nvmxsub;
    offsetue = (MYSUB+1)*nvmxsub2 + NVARS;
    for (i = 0; i < nvmxsub; i++) uext[offsetue+i] = udata[offsetu+i];
  }

  /* If isubx = 0, copy y-line 2 of u to uext */
  if (isubx == 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = ly*nvmxsub + NVARS;
      offsetue = (ly+1)*nvmxsub2;
      for (i = 0; i < NVARS; i++) uext[offsetue+i] = udata[offsetu+i];
    }
  }

  /* If isubx = NPEX-1, copy y-line MXSUB-1 of u to uext */
  if (isubx == NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetu = (ly+1)*nvmxsub - 2*NVARS;
      offsetue = (ly+2)*nvmxsub2 - NVARS;
      for (i = 0; i < NVARS; i++) uext[offsetue+i] = udata[offsetu+i];
    }
  }

  /* Make local copies of problem variables, for efficiency */
  dely = data->dy;
  verdco = data->vdco;
  hordco  = data->hdco;
  horaco  = data->haco;

  /* Set diurnal rate coefficients as functions of t, and save q4 in
  data block for use by preconditioner evaluation routine */
  s = sin((data->om)*t);
  if (s > RCONST(0.0)) {
    q3 = SUNRexp(-A3/s);
    q4coef = SUNRexp(-A4/s);
  } else {
    q3 = RCONST(0.0);
    q4coef = RCONST(0.0);
  }
  data->q4 = q4coef;

  /* Loop over all grid points in local subgrid */
  for (ly = 0; ly < MYSUB; ly++) {

    jy = ly + isuby*MYSUB;

    /* Set vertical diffusion coefficients at jy +- 1/2 */
    ydn = YMIN + (jy - RCONST(0.5))*dely;
    yup = ydn + dely;
    cydn = verdco*SUNRexp(RCONST(0.2)*ydn);
    cyup = verdco*SUNRexp(RCONST(0.2)*yup);
    for (lx = 0; lx < MXSUB; lx++) {

      /* Extract c1 and c2, and set kinetic rate terms */
      offsetue = (lx+1)*NVARS + (ly+1)*nvmxsub2;
      c1 = uext[offsetue];
      c2 = uext[offsetue+1];
      qq1 = Q1*c1*C3;
      qq2 = Q2*c1*c2;
      qq3 = q3*C3;
      qq4 = q4coef*c2;
      rkin1 = -qq1 - qq2 + RCONST(2.0)*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms */
      c1dn = uext[offsetue-nvmxsub2];
      c2dn = uext[offsetue-nvmxsub2+1];
      c1up = uext[offsetue+nvmxsub2];
      c2up = uext[offsetue+nvmxsub2+1];
      vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
      vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);

      /* Set horizontal diffusion and advection terms */
      c1lt = uext[offsetue-2];
      c2lt = uext[offsetue-1];
      c1rt = uext[offsetue+2];
      c2rt = uext[offsetue+3];
      hord1 = hordco*(c1rt - RCONST(2.0)*c1 + c1lt);
      hord2 = hordco*(c2rt - RCONST(2.0)*c2 + c2lt);
      horad1 = horaco*(c1rt - c1lt);
      horad2 = horaco*(c2rt - c2lt);

      /* Load all terms into dudata */
      offsetu = lx*NVARS + ly*nvmxsub;
      dudata[offsetu]   = vertd1 + hord1 + horad1 + rkin1;
      dudata[offsetu+1] = vertd2 + hord2 + horad2 + rkin2;
    }
  }
}


/***************** Functions Called by the Solver *************************/

/* f routine.  Evaluate f(t,y).  First call ucomm to do communication of
   subgrid boundary data into uext.  Then calculate f by a call to fcalc. */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype *udata, *udotdata;
  UserData data;
  HYPRE_ParVector uhyp;
  HYPRE_ParVector udothyp;

  /* Extract hypre vectors */
  uhyp  = N_VGetVector_ParHyp(u);
  udothyp  = N_VGetVector_ParHyp(udot);

  /* Access hypre vectors local data */
  udata = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));
  udotdata = hypre_VectorData(hypre_ParVectorLocalVector(udothyp));

  data = (UserData) user_data;

  /* Call ucomm to do inter-processor communication */
  ucomm(t, u, data);

  /* Call fcalc to calculate all right-hand sides */
  fcalc(t, udata, udotdata, data);

  return(0);
}


/* Preconditioner setup routine. Generate and preprocess P. */

static int Precond(realtype tn, N_Vector u, N_Vector fu,
                   booleantype jok, booleantype *jcurPtr,
                   realtype gamma, void *user_data)
{
  realtype c1, c2, cydn, cyup, diag, ydn, yup, q4coef, dely, verdco, hordco;
  realtype **(*P)[MYSUB], **(*Jbd)[MYSUB];
  int nvmxsub, ier, offset;
  sunindextype *(*pivot)[MYSUB];
  int lx, ly, jy, isuby;
  realtype *udata, **a, **j;
  HYPRE_ParVector uhyp;
  UserData data;

  /* Make local copies of pointers in user_data, pointer to u's data,
     and PE index pair */
  data = (UserData) user_data;
  P = data->P;
  Jbd = data->Jbd;
  pivot = data->pivot;
  isuby = data->isuby;
  nvmxsub = data->nvmxsub;

  uhyp  = N_VGetVector_ParHyp(u);
  udata = hypre_VectorData(hypre_ParVectorLocalVector(uhyp));

  if (jok) {
  /* jok = SUNTRUE: Copy Jbd to P */
    for (ly = 0; ly < MYSUB; ly++)
      for (lx = 0; lx < MXSUB; lx++)
        denseCopy(Jbd[lx][ly], P[lx][ly], NVARS, NVARS);

  *jcurPtr = SUNFALSE;

  }

  else {

    /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */

    /* Make local copies of problem variables, for efficiency */
    q4coef = data->q4;
    dely = data->dy;
    verdco = data->vdco;
    hordco  = data->hdco;

    /* Compute 2x2 diagonal Jacobian blocks (using q4 values
     computed on the last f call).  Load into P. */
    for (ly = 0; ly < MYSUB; ly++) {
      jy = ly + isuby*MYSUB;
      ydn = YMIN + (jy - RCONST(0.5))*dely;
      yup = ydn + dely;
      cydn = verdco*SUNRexp(RCONST(0.2)*ydn);
      cyup = verdco*SUNRexp(RCONST(0.2)*yup);
      diag = -(cydn + cyup + RCONST(2.0)*hordco);
      for (lx = 0; lx < MXSUB; lx++) {
        offset = lx*NVARS + ly*nvmxsub;
        c1 = udata[offset];
        c2 = udata[offset+1];
        j = Jbd[lx][ly];
        a = P[lx][ly];
        IJth(j,1,1) = (-Q1*C3 - Q2*c2) + diag;
        IJth(j,1,2) = -Q2*c1 + q4coef;
        IJth(j,2,1) = Q1*C3 - Q2*c2;
        IJth(j,2,2) = (-Q2*c1 - q4coef) + diag;
        denseCopy(j, a, NVARS, NVARS);
      }
    }

    *jcurPtr = SUNTRUE;

  }

  /* Scale by -gamma */
  for (ly = 0; ly < MYSUB; ly++)
    for (lx = 0; lx < MXSUB; lx++)
      denseScale(-gamma, P[lx][ly], NVARS, NVARS);

  /* Add identity matrix and do LU decompositions on blocks in place */
  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      denseAddIdentity(P[lx][ly], NVARS);
      ier = denseGETRF(P[lx][ly], NVARS, NVARS, pivot[lx][ly]);
      if (ier != 0) return(1);
    }
  }

  return(0);
}


/* Preconditioner solve routine */
static int PSolve(realtype tn, N_Vector u, N_Vector fu,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *user_data)
{
  realtype **(*P)[MYSUB];
  int nvmxsub;
  sunindextype *(*pivot)[MYSUB];
  int lx, ly;
  realtype *zdata, *v;
  HYPRE_ParVector zhyp;
  UserData data;

  /* Extract the P and pivot arrays from user_data */
  data = (UserData) user_data;
  P = data->P;
  pivot = data->pivot;

  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z.
     First copy vector r to z. */
  N_VScale(RCONST(1.0), r, z);
  nvmxsub = data->nvmxsub;
  zhyp  = N_VGetVector_ParHyp(z); /* extract hypre vector */
  zdata = hypre_VectorData(hypre_ParVectorLocalVector(zhyp));

  for (lx = 0; lx < MXSUB; lx++) {
    for (ly = 0; ly < MYSUB; ly++) {
      v = &(zdata[lx*NVARS + ly*nvmxsub]);
      denseGETRS(P[lx][ly], NVARS, pivot[lx][ly], v);
    }
  }

  return(0);
}


/*********************** Private Helper Function ************************/

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag(void *flagvalue, const char *funcname, int opt, int id)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1);
  } else if (opt == 1) { /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed with flag = %d\n\n",
              id, funcname, *errflag);
      return(1);
    }
  } else if (opt == 2 && flagvalue == NULL) { /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1);
  }

  return(0);
}
