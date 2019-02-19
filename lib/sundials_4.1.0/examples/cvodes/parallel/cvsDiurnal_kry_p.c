/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds
 *                S. D. Cohen, A. C. Hindmarsh, M. R. Wittman, and
 *                Radu Serban  @ LLNL
 * -----------------------------------------------------------------
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
 * The problem is solved by CVODES on NPE processors, treated
 * as a rectangular process grid of size NPEX by NPEY, with
 * NPE = NPEX*NPEY. Each processor contains a subgrid of size MXSUB
 * by MYSUB of the (x,y) mesh.  Thus the actual mesh sizes are
 * MX = MXSUB*NPEX and MY = MYSUB*NPEY, and the ODE system size is
 * neq = 2*MX*MY.
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
 * Execution: mpirun -np N cvsDiurnal_kry_p   with N = NPEX*NPEY
 * (see constants below).
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvodes/cvodes.h>              /* prototypes for CVODE fcts. */
#include <sunlinsol/sunlinsol_spgmr.h>  /* prototypes and constants for SUNLinSol_SPGMR solver */
#include <nvector/nvector_parallel.h>   /* definition N_Vector  */
#include <sundials/sundials_dense.h>    /* prototypes for small dense matrix fcts. */
#include <sundials/sundials_types.h>    /* definitions of realtype, booleantype */
#include <sundials/sundials_math.h>     /* definition of macros SUNSQR and EXP */

#include <mpi.h>                       /* MPI constants and types */

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
  realtype uext[NVARS*(MXSUB+2)*(MYSUB+2)];
  int my_pe, isubx, isuby;
  sunindextype nvmxsub, nvmxsub2;
  MPI_Comm comm;

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
                  N_Vector r, N_Vector z, realtype gamma, 
                  realtype delta, int lr, void *user_data);


/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);


/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  realtype abstol, reltol, t, tout;
  N_Vector u;
  UserData data;
  SUNLinearSolver LS;
  void *cvode_mem;
  int iout, retval, my_pe, npes;
  sunindextype neq, local_N;
  MPI_Comm comm;

  u = NULL;
  data = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Set problem size neq */
  neq = NVARS*MX*MY;

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

  /* Allocate and load user data block; allocate preconditioner block */
  data = (UserData) malloc(sizeof *data);
  if (check_retval((void *)data, "malloc", 2, my_pe)) MPI_Abort(comm, 1);
  InitUserData(my_pe, comm, data);

  /* Allocate u, and set initial values and tolerances */ 
  u = N_VNew_Parallel(comm, local_N, neq);
  if (check_retval((void *)u, "N_VNew", 0, my_pe)) MPI_Abort(comm, 1);
  SetInitialProfiles(u, data);
  abstol = ATOL; reltol = RTOL;

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
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

  /* Create SPGMR solver structure -- use left preconditioning 
     and the default Krylov dimension maxl */
  LS = SUNLinSol_SPGMR(u, PREC_LEFT, 0);
  if (check_retval((void *)LS, "SUNLinSol_SPGMR", 0, my_pe)) MPI_Abort(comm, 1);
  
  /* Attach the linear solver */
  retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1, my_pe)) MPI_Abort(comm, 1);

  /* Set preconditioner setup and solve routines Precond and PSolve, 
     and the pointer to the user-defined block data */
  retval = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
  if (check_retval(&retval, "CVodeSetPreconditioner", 1, my_pe)) MPI_Abort(comm, 1);

  if (my_pe == 0)
    printf("\n2-species diurnal advection-diffusion problem\n\n");

  /* In loop over output points, call CVode, print results, test for error */
  for (iout=1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1, my_pe)) break;
    PrintOutput(cvode_mem, my_pe, comm, u, t);
  }

  /* Print final statistics */  
  if (my_pe == 0) PrintFinalStats(cvode_mem);

  /* Free memory */
  N_VDestroy_Parallel(u);
  FreeUserData(data);
  CVodeFree(&cvode_mem);
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

static void SetInitialProfiles(N_Vector u, UserData data)
{
  int isubx, isuby, lx, ly, jx, jy;
  sunindextype offset;
  realtype dx, dy, x, y, cx, cy, xmid, ymid;
  realtype *udata;

  /* Set pointer to data array in vector u */
  udata = N_VGetArrayPointer_Parallel(u);

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
      udata[offset  ] = C1_SCALE*cx*cy; 
      udata[offset+1] = C2_SCALE*cx*cy;
      offset = offset + 2;
    }
  }
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(void *cvode_mem, int my_pe, MPI_Comm comm,
                        N_Vector u, realtype t)
{
  int qu, retval;
  realtype hu, *udata, tempu[2];
  int npelast;
  long int  nst;
  sunindextype i0, i1;
  MPI_Status status;

  npelast = NPEX*NPEY - 1;
  udata = N_VGetArrayPointer_Parallel(u);

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
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_retval(&retval, "CVodeGetNumSteps", 1, my_pe);
    retval = CVodeGetLastOrder(cvode_mem, &qu);
    check_retval(&retval, "CVodeGetLastOrder", 1, my_pe);
    retval = CVodeGetLastStep(cvode_mem, &hu);
    check_retval(&retval, "CVodeGetLastStep", 1, my_pe);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("t = %.2Le   no. steps = %ld   order = %d   stepsize = %.2Le\n",
           t, nst, qu, hu);
    printf("At bottom left:  c1, c2 = %12.3Le %12.3Le \n", udata[0], udata[1]);
    printf("At top right:    c1, c2 = %12.3Le %12.3Le \n\n", tempu[0], tempu[1]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("t = %.2e   no. steps = %ld   order = %d   stepsize = %.2e\n",
           t, nst, qu, hu);
    printf("At bottom left:  c1, c2 = %12.3e %12.3e \n", udata[0], udata[1]);
    printf("At top right:    c1, c2 = %12.3e %12.3e \n\n", tempu[0], tempu[1]);
#else
    printf("t = %.2e   no. steps = %ld   order = %d   stepsize = %.2e\n",
           t, nst, qu, hu);
    printf("At bottom left:  c1, c2 = %12.3e %12.3e \n", udata[0], udata[1]);
    printf("At top right:    c1, c2 = %12.3e %12.3e \n\n", tempu[0], tempu[1]);
#endif
  }
}

/* Print final statistics contained in iopt */

static void PrintFinalStats(void *cvode_mem)
{
  long int lenrw, leniw ;
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
  printf("lenrw   = %5ld     leniw   = %5ld\n"  , lenrw, leniw);
  printf("lenrwls = %5ld     leniwls = %5ld\n"  , lenrwLS, leniwLS);
  printf("nst     = %5ld\n"                     , nst);
  printf("nfe     = %5ld     nfels   = %5ld\n"  , nfe, nfeLS);
  printf("nni     = %5ld     nli     = %5ld\n"  , nni, nli);
  printf("nsetups = %5ld     netf    = %5ld\n"  , nsetups, netf);
  printf("npe     = %5ld     nps     = %5ld\n"  , npe, nps);
  printf("ncfn    = %5ld     ncfl    = %5ld\n\n", ncfn, ncfl); 
}
 
/* Routine to send boundary data to neighboring PEs */

static void BSend(MPI_Comm comm, 
                  int my_pe, int isubx, int isuby, 
                  sunindextype dsizex, sunindextype dsizey,
                  realtype udata[])
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

  udata = N_VGetArrayPointer_Parallel(u);

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
  realtype *udata, *dudata;
  UserData data;

  udata = N_VGetArrayPointer_Parallel(u);
  dudata = N_VGetArrayPointer_Parallel(udot);
  data = (UserData) user_data;

  /* Call ucomm to do inter-processor communication */
  ucomm(t, u, data);

  /* Call fcalc to calculate all right-hand sides */
  fcalc(t, udata, dudata, data);

  return(0);
}

/* Preconditioner setup routine. Generate and preprocess P. */
static int Precond(realtype tn, N_Vector u, N_Vector fu,
                   booleantype jok, booleantype *jcurPtr, 
                   realtype gamma, void *user_data)
{
  realtype c1, c2, cydn, cyup, diag, ydn, yup, q4coef, dely, verdco, hordco;
  realtype **(*P)[MYSUB], **(*Jbd)[MYSUB];
  sunindextype nvmxsub, offset;
  sunindextype *(*pivot)[MYSUB], retval;
  int lx, ly, jy, isuby;
  realtype *udata, **a, **j;
  UserData data;

  /* Make local copies of pointers in user_data, pointer to u's data,
     and PE index pair */
  data = (UserData) user_data;
  P = data->P;
  Jbd = data->Jbd;
  pivot = data->pivot;
  udata = N_VGetArrayPointer_Parallel(u);
  isuby = data->isuby;
  nvmxsub = data->nvmxsub;

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
      retval = denseGETRF(P[lx][ly], NVARS, NVARS, pivot[lx][ly]);
      if (retval != 0) return(1);
    }
  }

  return(0);
}

/* Preconditioner solve routine */
static int PSolve(realtype tn, N_Vector u, N_Vector fu, 
                  N_Vector r, N_Vector z, realtype gamma, 
                  realtype delta, int lr, void *user_data)
{
  realtype **(*P)[MYSUB];
  sunindextype nvmxsub;
  sunindextype *(*pivot)[MYSUB];
  int lx, ly;
  realtype *zdata, *v;
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
  zdata = N_VGetArrayPointer_Parallel(z);

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
