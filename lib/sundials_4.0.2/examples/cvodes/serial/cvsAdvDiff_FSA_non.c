/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, George D. Byrne,
 *              and Radu Serban @ LLNL
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
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to the fixed point SUNNonlinearSolver */

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
static void PrintFinalStats(void *cvode_mem, booleantype sensi,
                            booleantype err_con, int sensi_meth);

static int check_retval(void *returnvalue, const char *funcname, int opt);

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
  int iout, retval;

  realtype *pbar;
  int is, *plist;
  N_Vector *uS;
  booleantype sensi, err_con;
  int sensi_meth;

  SUNNonlinearSolver NLS, NLSsens;

  cvode_mem = NULL;
  data = NULL;
  u = NULL;
  pbar = NULL;
  plist = NULL;
  uS = NULL;
  NLS = NULL;
  NLSsens = NULL;

  /* Process arguments */
  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);

  /* Set user data */
  data = (UserData) malloc(sizeof *data); /* Allocate data memory */
  if(check_retval((void *)data, "malloc", 2)) return(1);
  data->p = (realtype *) malloc(NP * sizeof(realtype));
  dx = data->dx = XMAX/((realtype)(MX+1));
  data->p[0] = RCONST(1.0);
  data->p[1] = RCONST(0.5);

  /* Allocate and set initial states */
  u = N_VNew_Serial(NEQ);
  if(check_retval((void *)u, "N_VNew_Serial", 0)) return(1);
  SetIC(u, dx);

  /* Set integration tolerances */
  reltol = ZERO;
  abstol = ATOL;

  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_ADAMS);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Allocate CVODES memory */
  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);

  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* create fixed point nonlinear solver object */
  NLS = SUNNonlinSol_FixedPoint(u, 0);
  if(check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0)) return(1);

  /* attach nonlinear solver object to CVode */
  retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if(check_retval(&retval, "CVodeSetNonlinearSolver", 1)) return(1);

  printf("\n1-D advection-diffusion equation, mesh size =%3d\n", MX);

  /* Sensitivity-related settings */
  if(sensi) {

    plist = (int *) malloc(NS * sizeof(int));
    if(check_retval((void *)plist, "malloc", 2)) return(1);
    for(is=0; is<NS; is++) plist[is] = is;

    pbar  = (realtype *) malloc(NS * sizeof(realtype));
    if(check_retval((void *)pbar, "malloc", 2)) return(1);
    for(is=0; is<NS; is++) pbar[is] = data->p[plist[is]];

    uS = N_VCloneVectorArray(NS, u);
    if(check_retval((void *)uS, "N_VCloneVectorArray", 0)) return(1);
    for(is=0;is<NS;is++)
      N_VConst(ZERO, uS[is]);

    retval = CVodeSensInit1(cvode_mem, NS, sensi_meth, NULL, uS);
    if(check_retval(&retval, "CVodeSensInit1", 1)) return(1);

    retval = CVodeSensEEtolerances(cvode_mem);
    if(check_retval(&retval, "CVodeSensEEtolerances", 1)) return(1);

    retval = CVodeSetSensErrCon(cvode_mem, err_con);
    if(check_retval(&retval, "CVodeSetSensErrCon", 1)) return(1);

    retval = CVodeSetSensDQMethod(cvode_mem, CV_CENTERED, ZERO);
    if(check_retval(&retval, "CVodeSetSensDQMethod", 1)) return(1);

    retval = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);
    if(check_retval(&retval, "CVodeSetSensParams", 1)) return(1);

    /* create sensitivity fixed point nonlinear solver object */
    if (sensi_meth == CV_SIMULTANEOUS)
      NLSsens = SUNNonlinSol_FixedPointSens(NS+1, u, 0);
    else if(sensi_meth == CV_STAGGERED)
      NLSsens = SUNNonlinSol_FixedPointSens(NS, u, 0);
    else
      NLSsens = SUNNonlinSol_FixedPoint(u, 0);
    if(check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0)) return(1);

    /* attach nonlinear solver object to CVode */
    if (sensi_meth == CV_SIMULTANEOUS)
      retval = CVodeSetNonlinearSolverSensSim(cvode_mem, NLSsens);
    else if(sensi_meth == CV_STAGGERED)
      retval = CVodeSetNonlinearSolverSensStg(cvode_mem, NLSsens);
    else
      retval = CVodeSetNonlinearSolverSensStg1(cvode_mem, NLSsens);
    if(check_retval(&retval, "CVodeSetNonlinearSolver", 1)) return(1);

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
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_retval(&retval, "CVode", 1)) break;
    PrintOutput(cvode_mem, t, u);
    if (sensi) {
      retval = CVodeGetSens(cvode_mem, &t, uS);
      if(check_retval(&retval, "CVodeGetSens", 1)) break;
      PrintOutputS(uS);
    } 
    printf("------------------------------------------------------------\n");
  }

  /* Print final statistics */
  PrintFinalStats(cvode_mem, sensi, err_con, sensi_meth);

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
  SUNNonlinSolFree(NLS);
  if (sensi) SUNNonlinSolFree(NLSsens);

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
  int qu, retval;
  realtype hu;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetLastOrder(cvode_mem, &qu);
  check_retval(&retval, "CVodeGetLastOrder", 1);
  retval = CVodeGetLastStep(cvode_mem, &hu);
  check_retval(&retval, "CVodeGetLastStep", 1);

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

static void PrintFinalStats(void *cvode_mem, booleantype sensi,
                            booleantype err_con, int sensi_meth)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  if (sensi) {
    retval = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    check_retval(&retval, "CVodeGetSensNumRhsEvals", 1);
    retval = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    check_retval(&retval, "CVodeGetNumRhsEvalsSens", 1);
    retval = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    check_retval(&retval, "CVodeGetSensNumLinSolvSetups", 1);
    if (err_con) {
      retval = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
      check_retval(&retval, "CVodeGetSensNumErrTestFails", 1);
    } else {
      netfS = 0;
    }
    if ((sensi_meth == CV_STAGGERED) || (sensi_meth == CV_STAGGERED1)) {
      retval = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
      check_retval(&retval, "CVodeGetSensNumNonlinSolvIters", 1);
      retval = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
      check_retval(&retval, "CVodeGetSensNumNonlinSolvConvFails", 1);
    } else {
      nniS = 0;
      ncfnS = 0;
    }
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
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n", 
              funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1); }

  return(0);
}
