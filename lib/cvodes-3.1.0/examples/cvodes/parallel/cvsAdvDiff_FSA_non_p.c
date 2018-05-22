/*
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, George D. Byrne,
 *                and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the program for
 * its solution by CVODES. The problem is the semi-discrete form of
 * the advection-diffusion equation in 1-D:
 *   du/dt = q1 * d^2 u / dx^2 + q2 * du/dx
 * on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
 * Homogeneous Dirichlet boundary conditions are posed, and the
 * initial condition is:
 *   u(x,y,t=0) = x(2-x)exp(2x).
 * The PDE is discretized on a uniform grid of size MX+2 with
 * central differencing, and with boundary values eliminated,
 * leaving an ODE system of size NEQ = MX.
 * This program solves the problem with the option for nonstiff
 * systems: ADAMS method and functional iteration.
 * It uses scalar relative and absolute tolerances.
 * Output is printed at t = .5, 1.0, ..., 5.
 * Run statistics (optional outputs) are printed at the end.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters q1 and q2.
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on FULL or PARTIAL,
 * respectively).
 *
 * Execution:
 *
 * Note: This version uses MPI for user routines, and the CVODES
 *       solver. In what follows, N is the number of processors,
 *       N = NPEX*NPEY (see constants below) and it is assumed that
 *       the MPI script mpirun is used to run a parallel
 *       application.
 * If no sensitivities are desired:
 *    % mpirun -np N cvsAdvDiff_FSA_non_p -nosensi
 * If sensitivities are to be computed:
 *    % mpirun -np N cvsAdvDiff_FSA_non_p -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <cvodes/cvodes.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <mpi.h>

/* Problem Constants */
#define XMAX  RCONST(2.0)   /* domain boundary           */
#define MX    10            /* mesh dimension            */
#define NEQ   MX            /* number of equations       */
#define ATOL  RCONST(1.e-5) /* scalar absolute tolerance */
#define T0    RCONST(0.0)   /* initial time              */
#define T1    RCONST(0.5)   /* first output time         */
#define DTOUT RCONST(0.5)   /* output time increment     */
#define NOUT  10            /* number of output times    */

#define NP    2
#define NS    2

#define ZERO  RCONST(0.0)

/* Type : UserData 
   contains problem parameters, grid constants, work array. */

typedef struct {
  realtype *p;
  realtype dx;
  int npes, my_pe;
  MPI_Comm comm;
  realtype z[100];
} *UserData;


/* Prototypes of user-supplied functins */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

/* Prototypes of private functions */

static void ProcessArgs(int argc, char *argv[], int my_pe,
                        booleantype *sensi, int *sensi_meth, booleantype *err_con);
static void WrongArgs(int my_pe, char *name);
static void SetIC(N_Vector u, realtype dx, sunindextype my_length, sunindextype my_base);
static void PrintOutput(void *cvode_mem, int my_pe, realtype t, N_Vector u);
static void PrintOutputS(int my_pe, N_Vector *uS);
static void PrintFinalStats(void *cvode_mem, booleantype sensi); 
static int check_flag(void *flagvalue, const char *funcname, int opt, int id);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  realtype dx, reltol, abstol, t, tout;
  N_Vector u;
  UserData data;
  void *cvode_mem;
  int iout, flag, my_pe, npes;
  sunindextype local_N, nperpe, nrem, my_base;

  realtype *pbar;
  int is, *plist;
  N_Vector *uS;
  booleantype sensi, err_con;
  int sensi_meth;

  MPI_Comm comm;

  u = NULL;
  data = NULL;
  cvode_mem = NULL;
  pbar = NULL;
  plist = NULL;
  uS = NULL;

  /* Get processor number, total number of pe's, and my_pe. */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  /* Process arguments */
  ProcessArgs(argc, argv, my_pe, &sensi, &sensi_meth, &err_con);

  /* Set local vector length. */
  nperpe = NEQ/npes;
  nrem = NEQ - npes*nperpe;
  local_N = (my_pe < nrem) ? nperpe+1 : nperpe;
  my_base = (my_pe < nrem) ? my_pe*local_N : my_pe*nperpe + nrem;

  /* USER DATA STRUCTURE */
  data = (UserData) malloc(sizeof *data); /* Allocate data memory */
  data->p = NULL;
  if(check_flag((void *)data, "malloc", 2, my_pe)) MPI_Abort(comm, 1);
  data->comm = comm;
  data->npes = npes;
  data->my_pe = my_pe;
  data->p = (realtype *) malloc(NP * sizeof(realtype));
  if(check_flag((void *)data->p, "malloc", 2, my_pe)) MPI_Abort(comm, 1);
  dx = data->dx = XMAX/((realtype)(MX+1));
  data->p[0] = RCONST(1.0);
  data->p[1] = RCONST(0.5);

  /* INITIAL STATES */
  u = N_VNew_Parallel(comm, local_N, NEQ);    /* Allocate u vector */
  if(check_flag((void *)u, "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  SetIC(u, dx, local_N, my_base);    /* Initialize u vector */

  /* TOLERANCES */
  reltol = ZERO;                /* Set the tolerances */
  abstol = ATOL;

  /* CVODE_CREATE & CVODE_MALLOC */
  cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  if(check_flag((void *)cvode_mem, "CVodeCreate", 0, my_pe)) MPI_Abort(comm, 1);

  flag = CVodeSetUserData(cvode_mem, data);
  if(check_flag(&flag, "CVodeSetUserData", 1, my_pe)) MPI_Abort(comm, 1);

  flag = CVodeInit(cvode_mem, f, T0, u);
  if(check_flag(&flag, "CVodeInit", 1, my_pe)) MPI_Abort(comm, 1);
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_flag(&flag, "CVodeSStolerances", 1, my_pe)) MPI_Abort(comm, 1);

 
  if (my_pe == 0) {
    printf("\n1-D advection-diffusion equation, mesh size =%3d \n", MX);
    printf("\nNumber of PEs = %3d \n",npes);
  }

  if(sensi) {

    plist = (int *) malloc(NS * sizeof(int));
    if(check_flag((void *)plist, "malloc", 2, my_pe)) MPI_Abort(comm, 1);
    for(is=0; is<NS; is++)
      plist[is] = is; /* sensitivity w.r.t. i-th parameter */

    pbar  = (realtype *) malloc(NS * sizeof(realtype));
    if(check_flag((void *)pbar, "malloc", 2, my_pe)) MPI_Abort(comm, 1);
    for(is=0; is<NS; is++) pbar[is] = data->p[plist[is]];

    uS = N_VCloneVectorArray_Parallel(NS, u);
    if(check_flag((void *)uS, "N_VCloneVectorArray_Parallel", 0, my_pe)) 
      MPI_Abort(comm, 1);
    for(is=0;is<NS;is++)
      N_VConst(ZERO,uS[is]);

    flag = CVodeSensInit1(cvode_mem, NS, sensi_meth, NULL, uS);
    if(check_flag(&flag, "CVodeSensInit1", 1, my_pe)) MPI_Abort(comm, 1);

    flag = CVodeSensEEtolerances(cvode_mem);
    if(check_flag(&flag, "CVodeSensEEtolerances", 1, my_pe)) MPI_Abort(comm, 1);

    flag = CVodeSetSensErrCon(cvode_mem, err_con);
    if(check_flag(&flag, "CVodeSetSensErrCon", 1, my_pe)) MPI_Abort(comm, 1);

    flag = CVodeSetSensDQMethod(cvode_mem, CV_CENTERED, ZERO);
    if(check_flag(&flag, "CVodeSetSensDQMethod", 1, my_pe)) MPI_Abort(comm, 1);

    flag = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);
    if(check_flag(&flag, "CVodeSetSensParams", 1, my_pe)) MPI_Abort(comm, 1);

    if(my_pe == 0) {
      printf("Sensitivity: YES ");
      if(sensi_meth == CV_SIMULTANEOUS)   
        printf("( SIMULTANEOUS +");
      else 
        if(sensi_meth == CV_STAGGERED) printf("( STAGGERED +");
        else                           printf("( STAGGERED1 +");   
      if(err_con) printf(" FULL ERROR CONTROL )");
      else        printf(" PARTIAL ERROR CONTROL )");
    }

  } else {

    if(my_pe == 0) printf("Sensitivity: NO ");

  }

  /* In loop over output points, call CVode, print results, test for error */

  if(my_pe == 0) {
    printf("\n\n");
    printf("============================================================\n");
    printf("     T     Q       H      NST                    Max norm   \n");
    printf("============================================================\n");
  }

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {

    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_flag(&flag, "CVode", 1, my_pe)) break;
    PrintOutput(cvode_mem, my_pe, t, u);
    if (sensi) {
      flag = CVodeGetSens(cvode_mem, &t, uS);
      if(check_flag(&flag, "CVodeGetSens", 1, my_pe)) break;
      PrintOutputS(my_pe, uS);
    }
    if (my_pe == 0)
      printf("------------------------------------------------------------\n");

  }

  /* Print final statistics */
  if (my_pe == 0) 
    PrintFinalStats(cvode_mem, sensi);

  /* Free memory */
  N_VDestroy(u);                   /* Free the u vector              */
  if (sensi) 
    N_VDestroyVectorArray(uS, NS); /* Free the uS vectors            */
  free(data->p);                   /* Free the p vector              */
  free(data);                      /* Free block of UserData         */
  CVodeFree(&cvode_mem);           /* Free the CVODES problem memory */
  free(pbar);
  if(sensi) free(plist);

  MPI_Finalize();

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,u). 
 */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype ui, ult, urt, hordc, horac, hdiff, hadv;
  realtype *udata, *dudata, *z;
  realtype dx;
  int i;
  int npes, my_pe, my_length, my_pe_m1, my_pe_p1, last_pe, my_last;
  UserData data;
  MPI_Status status;
  MPI_Comm comm;

  udata = N_VGetArrayPointer_Parallel(u);
  dudata = N_VGetArrayPointer_Parallel(udot);

  /* Extract needed problem constants from data */
  data  = (UserData) user_data;
  dx    = data->dx; 
  hordc = data->p[0]/(dx*dx);
  horac = data->p[1]/(RCONST(2.0)*dx);

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
  my_last = my_length - 1;

  /* Store local segment of u in the working array z. */
   for (i = 1; i <= my_length; i++)
     z[i] = udata[i - 1];

  /* Pass needed data to processes before and after current process. */
   if (my_pe != 0)
     MPI_Send(&z[1], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm);
   if (my_pe != last_pe)
     MPI_Send(&z[my_length], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm);   

  /* Receive needed data from processes before and after current process. */
   if (my_pe != 0)
     MPI_Recv(&z[0], 1, PVEC_REAL_MPI_TYPE, my_pe_m1, 0, comm, &status);
   else z[0] = ZERO;
   if (my_pe != last_pe)
     MPI_Recv(&z[my_length+1], 1, PVEC_REAL_MPI_TYPE, my_pe_p1, 0, comm,
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

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/* 
 * Process and verify arguments to cvsfwdnonx_p.
 */

static void ProcessArgs(int argc, char *argv[], int my_pe,
                        booleantype *sensi, int *sensi_meth, booleantype *err_con)
{
  *sensi = SUNFALSE;
  *sensi_meth = -1;
  *err_con = SUNFALSE;

  if (argc < 2) WrongArgs(my_pe, argv[0]);

  if (strcmp(argv[1],"-nosensi") == 0)
    *sensi = SUNFALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    *sensi = SUNTRUE;
  else
    WrongArgs(my_pe, argv[0]);
  
  if (*sensi) {

    if (argc != 4)
      WrongArgs(my_pe, argv[0]);

    if (strcmp(argv[2],"sim") == 0)
      *sensi_meth = CV_SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      *sensi_meth = CV_STAGGERED;
    else if (strcmp(argv[2],"stg1") == 0)
      *sensi_meth = CV_STAGGERED1;
    else 
      WrongArgs(my_pe, argv[0]);

    if (strcmp(argv[3],"t") == 0)
      *err_con = SUNTRUE;
    else if (strcmp(argv[3],"f") == 0)
      *err_con = SUNFALSE;
    else
      WrongArgs(my_pe, argv[0]);
  }

}

static void WrongArgs(int my_pe, char *name)
{
  if (my_pe == 0) {
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = t or f\n");
  }  
  MPI_Finalize();
  exit(0);
}

/*
 * Set initial conditions in u vector 
 */

static void SetIC(N_Vector u, realtype dx, sunindextype my_length, 
                  sunindextype my_base)
{
  int i;
  sunindextype iglobal;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u. */
  udata = N_VGetArrayPointer_Parallel(u);
  my_length = N_VGetLocalLength_Parallel(u);

  /* Load initial profile into u vector */
  for (i=1; i<=my_length; i++) {
    iglobal = my_base + i;
    x = iglobal*dx;
    udata[i-1] = x*(XMAX - x)*SUNRexp(2.0*x);
  }  
}

/*
 * Print current t, step count, order, stepsize, and max norm of solution  
 */

static void PrintOutput(void *cvode_mem, int my_pe, realtype t, N_Vector u)
{
  long int nst;
  int qu, flag;
  realtype hu, umax;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1, my_pe);
  flag = CVodeGetLastOrder(cvode_mem, &qu);
  check_flag(&flag, "CVodeGetLastOrder", 1, my_pe);
  flag = CVodeGetLastStep(cvode_mem, &hu);
  check_flag(&flag, "CVodeGetLastStep", 1, my_pe);

  umax = N_VMaxNorm(u);

  if (my_pe == 0) {

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%8.3Le %2d  %8.3Le %5ld\n", t,qu,hu,nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%8.3e %2d  %8.3e %5ld\n", t,qu,hu,nst);
#else
    printf("%8.3e %2d  %8.3e %5ld\n", t,qu,hu,nst);
#endif

    printf("                                Solution       ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%12.4Le \n", umax);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%12.4e \n", umax);
#else
    printf("%12.4e \n", umax);
#endif

  }  

}

/*
 * Print max norm of sensitivities 
 */

static void PrintOutputS(int my_pe, N_Vector *uS)
{
  realtype smax;

  smax = N_VMaxNorm(uS[0]);
  if (my_pe == 0) {
    printf("                                Sensitivity 1  ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%12.4Le \n", smax);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%12.4e \n", smax);
#else
    printf("%12.4e \n", smax);
#endif
  }

  smax = N_VMaxNorm(uS[1]);
  if (my_pe == 0) {
    printf("                                Sensitivity 2  ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%12.4Le \n", smax);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%12.4e \n", smax);
#else
    printf("%12.4e \n", smax);
#endif
  }

}

/*
 * Print some final statistics located in the iopt array 
 */

static void PrintFinalStats(void *cvode_mem, booleantype sensi) 
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1, 0);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1, 0);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1, 0);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1, 0);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1, 0);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1, 0);

  if (sensi) {
    flag = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    check_flag(&flag, "CVodeGetSensNumRhsEvals", 1, 0);
    flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    check_flag(&flag, "CVodeGetNumRhsEvalsSens", 1, 0);
    flag = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    check_flag(&flag, "CVodeGetSensNumLinSolvSetups", 1, 0);
    flag = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
    check_flag(&flag, "CVodeGetSensNumErrTestFails", 1, 0);
    flag = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
    check_flag(&flag, "CVodeGetSensNumNonlinSolvIters", 1, 0);
    flag = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
    check_flag(&flag, "CVodeGetSensNumNonlinSolvConvFails", 1, 0);
  }

  printf("\nFinal Statistics\n\n");
  printf("nst     = %5ld\n\n", nst);
  printf("nfe     = %5ld\n",   nfe);
  printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  if(sensi) {
    printf("\n");
    printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }

}

/* 
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt, int id)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n",
	    id, funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed with flag = %d\n\n",
	      id, funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n",
	    id, funcname);
    return(1); }

  return(0);
}
