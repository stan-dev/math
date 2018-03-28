/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, George D. Byrne,
 *              and Radu Serban @ LLNL
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
 * If no sensitivities are desired:
 *    % cvsAdvDiff_FSA_non -nosensi
 * If sensitivities are to be computed:
 *    % cvsAdvDiff_FSA_non -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

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
} *UserData;

/* Functions Called by the CVODES Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

/* Private Helper Functions */

static void ProcessArgs(int argc, char *argv[],
                        booleantype *sensi, int *sensi_meth,
			booleantype *err_con);
static void WrongArgs(char *name);
static void SetIC(N_Vector u, realtype dx);
static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
static void PrintOutputS(N_Vector *uS);
static void PrintFinalStats(void *cvode_mem, booleantype sensi);

static int check_flag(void *flagvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  void *cvode_mem;
  UserData data;
  realtype dx, reltol, abstol, t, tout;
  N_Vector u;
  int iout, flag;

  realtype *pbar;
  int is, *plist;
  N_Vector *uS;
  booleantype sensi, err_con;
  int sensi_meth;

  cvode_mem = NULL;
  data = NULL;
  u = NULL;
  pbar = NULL;
  plist = NULL;
  uS = NULL;

  /* Process arguments */
  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);

  /* Set user data */
  data = (UserData) malloc(sizeof *data); /* Allocate data memory */
  if(check_flag((void *)data, "malloc", 2)) return(1);
  data->p = (realtype *) malloc(NP * sizeof(realtype));
  dx = data->dx = XMAX/((realtype)(MX+1));
  data->p[0] = RCONST(1.0);
  data->p[1] = RCONST(0.5);

  /* Allocate and set initial states */
  u = N_VNew_Serial(NEQ);
  if(check_flag((void *)u, "N_VNew_Serial", 0)) return(1);
  SetIC(u, dx);

  /* Set integration tolerances */
  reltol = ZERO;
  abstol = ATOL;

  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  flag = CVodeSetUserData(cvode_mem, data);
  if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Allocate CVODES memory */
  flag = CVodeInit(cvode_mem, f, T0, u);
  if(check_flag(&flag, "CVodeInit", 1)) return(1);

  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  printf("\n1-D advection-diffusion equation, mesh size =%3d\n", MX);

  /* Sensitivity-related settings */
  if(sensi) {

    plist = (int *) malloc(NS * sizeof(int));
    if(check_flag((void *)plist, "malloc", 2)) return(1);
    for(is=0; is<NS; is++) plist[is] = is;

    pbar  = (realtype *) malloc(NS * sizeof(realtype));
    if(check_flag((void *)pbar, "malloc", 2)) return(1);
    for(is=0; is<NS; is++) pbar[is] = data->p[plist[is]];

    uS = N_VCloneVectorArray(NS, u);
    if(check_flag((void *)uS, "N_VCloneVectorArray", 0)) return(1);
    for(is=0;is<NS;is++)
      N_VConst(ZERO, uS[is]);

    flag = CVodeSensInit1(cvode_mem, NS, sensi_meth, NULL, uS);
    if(check_flag(&flag, "CVodeSensInit1", 1)) return(1);

    flag = CVodeSensEEtolerances(cvode_mem);
    if(check_flag(&flag, "CVodeSensEEtolerances", 1)) return(1);

    flag = CVodeSetSensErrCon(cvode_mem, err_con);
    if(check_flag(&flag, "CVodeSetSensErrCon", 1)) return(1);

    flag = CVodeSetSensDQMethod(cvode_mem, CV_CENTERED, ZERO);
    if(check_flag(&flag, "CVodeSetSensDQMethod", 1)) return(1);

    flag = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);
    if(check_flag(&flag, "CVodeSetSensParams", 1)) return(1);

    printf("Sensitivity: YES ");
    if(sensi_meth == CV_SIMULTANEOUS)   
      printf("( SIMULTANEOUS +");
    else 
      if(sensi_meth == CV_STAGGERED) printf("( STAGGERED +");
      else                           printf("( STAGGERED1 +");   
    if(err_con) printf(" FULL ERROR CONTROL )");
    else        printf(" PARTIAL ERROR CONTROL )");

  } else {

    printf("Sensitivity: NO ");

  }

  /* In loop over output points, call CVode, print results, test for error */

  printf("\n\n");
  printf("============================================================\n");
  printf("     T     Q       H      NST                    Max norm   \n");
  printf("============================================================\n");

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_flag(&flag, "CVode", 1)) break;
    PrintOutput(cvode_mem, t, u);
    if (sensi) {
      flag = CVodeGetSens(cvode_mem, &t, uS);
      if(check_flag(&flag, "CVodeGetSens", 1)) break;
      PrintOutputS(uS);
    } 
    printf("------------------------------------------------------------\n");
  }

  /* Print final statistics */
  PrintFinalStats(cvode_mem, sensi);

  /* Free memory */
  N_VDestroy(u);
  if (sensi) {
    N_VDestroyVectorArray(uS, NS);
    free(plist);
    free(pbar);
  }
  free(data->p);
  free(data);
  CVodeFree(&cvode_mem);

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
  realtype dx;
  realtype *udata, *dudata;
  int i;
  UserData data;

  udata = N_VGetArrayPointer(u);
  dudata = N_VGetArrayPointer(udot);

  /* Extract needed problem constants from data */
  data = (UserData) user_data;
  dx    = data->dx;
  hordc = data->p[0]/(dx*dx);
  horac = data->p[1]/(RCONST(2.0)*dx);

  /* Loop over all grid points. */
  for (i=0; i<NEQ; i++) {

    /* Extract u at x_i and two neighboring points */
    ui = udata[i];
    if(i!=0) 
      ult = udata[i-1];
    else
      ult = ZERO;
    if(i!=NEQ-1)
      urt = udata[i+1];
    else
      urt = ZERO;

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - RCONST(2.0)*ui + urt);
    hadv = horac*(urt - ult);
    dudata[i] = hdiff + hadv;
  }

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/* 
 * Process and verify arguments.
 */

static void ProcessArgs(int argc, char *argv[], 
                        booleantype *sensi, int *sensi_meth, booleantype *err_con)
{
  *sensi = SUNFALSE;
  *sensi_meth = -1;
  *err_con = SUNFALSE;

  if (argc < 2) WrongArgs(argv[0]);

  if (strcmp(argv[1],"-nosensi") == 0)
    *sensi = SUNFALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    *sensi = SUNTRUE;
  else
    WrongArgs(argv[0]);
  
  if (*sensi) {

    if (argc != 4)
      WrongArgs(argv[0]);

    if (strcmp(argv[2],"sim") == 0)
      *sensi_meth = CV_SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      *sensi_meth = CV_STAGGERED;
    else if (strcmp(argv[2],"stg1") == 0)
      *sensi_meth = CV_STAGGERED1;
    else 
      WrongArgs(argv[0]);

    if (strcmp(argv[3],"t") == 0)
      *err_con = SUNTRUE;
    else if (strcmp(argv[3],"f") == 0)
      *err_con = SUNFALSE;
    else
      WrongArgs(argv[0]);
  }

}

static void WrongArgs(char *name)
{
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = t or f\n");
    
    exit(0);
}

/* 
 * Set initial conditions in u vector.
 */

static void SetIC(N_Vector u, realtype dx)
{
  int i;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u. */
  udata = N_VGetArrayPointer(u);

  /* Load initial profile into u vector */
  for (i=0; i<NEQ; i++) {
    x = (i+1)*dx;
    udata[i] = x*(XMAX - x)*SUNRexp(RCONST(2.0)*x);
  }  
}

/*
 * Print current t, step count, order, stepsize, and max norm of solution  
 */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector u)
{
  long int nst;
  int qu, flag;
  realtype hu;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetLastOrder(cvode_mem, &qu);
  check_flag(&flag, "CVodeGetLastOrder", 1);
  flag = CVodeGetLastStep(cvode_mem, &hu);
  check_flag(&flag, "CVodeGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.3Le %2d  %8.3Le %5ld\n", t, qu, hu ,nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu ,nst);
#else
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu ,nst);
#endif

  printf("                                Solution       ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le \n", N_VMaxNorm(u));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e \n", N_VMaxNorm(u));
#else
  printf("%12.4e \n", N_VMaxNorm(u));
#endif
}

/*
 * Print max norm of sensitivities 
 */

static void PrintOutputS(N_Vector *uS)
{
  printf("                                Sensitivity 1  ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le \n", N_VMaxNorm(uS[0]));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e \n", N_VMaxNorm(uS[0]));
#else
  printf("%12.4e \n", N_VMaxNorm(uS[0]));
#endif

  printf("                                Sensitivity 2  ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le \n", N_VMaxNorm(uS[1]));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e \n", N_VMaxNorm(uS[1]));
#else
  printf("%12.4e \n", N_VMaxNorm(uS[1]));
#endif
}


/*
 * Print some final statistics located in the CVODES memory
 */

static void PrintFinalStats(void *cvode_mem, booleantype sensi)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  if (sensi) {
    flag = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    check_flag(&flag, "CVodeGetSensNumRhsEvals", 1);
    flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    check_flag(&flag, "CVodeGetNumRhsEvalsSens", 1);
    flag = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    check_flag(&flag, "CVodeGetSensNumLinSolvSetups", 1);
    flag = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
    check_flag(&flag, "CVodeGetSensNumErrTestFails", 1);
    flag = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
    check_flag(&flag, "CVodeGetSensNumNonlinSolvIters", 1);
    flag = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
    check_flag(&flag, "CVodeGetSensNumNonlinSolvConvFails", 1);
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

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1); }

  return(0);
}
