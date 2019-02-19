/*
 * -----------------------------------------------------------------
 * Programmers: Radu Serban and Alan Hindmarsh, and Cody Balos @ LLNL
 * -----------------------------------------------------------------
 * Modification of the cvsRoberts_FSA_dns to illustrate switching
 * on and off sensitivity computations.
 *
 * Example problem (from cvsRoberts_FSA_dns):
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES for Forward Sensitivity
 * Analysis. The problem is from chemical kinetics, and consists
 * of the following three rate equations:
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions y1 = 1.0, y2 = y3 = 0. The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the DENSE linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 * Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 *------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cvodes/cvodes.h>             /* prototypes for CVODE functions and const */
#include <nvector/nvector_serial.h>    /* access to serial NVector                 */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix                */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver          */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype          */

/* Problem Constants */
#define MXSTEPS 2000            /* max number of steps */
#define NEQ     3               /* number of equations */
#define T0      RCONST(0.0)     /* initial time        */
#define T1      RCONST(4.0e10)  /* first output time   */

#define ZERO    RCONST(0.0)

/* Type : UserData */
typedef struct {
  booleantype sensi;     /* turn on (T) or off (F) sensitivity analysis    */
  booleantype errconS;   /* full (T) or partial error control (F)          */
  booleantype fsDQ;      /* user provided r.h.s sensitivity analysis (T/F) */
  int meth;              /* sensitivity method                             */
  realtype p[3];         /* sensitivity variables                          */
} *UserData;

/* User provided routine called by the solver to compute
 * the function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *udata);

/* User provided routine called by the solver to
 * approximate the Jacobian J(t,y).  */
static int Jac(realtype t, N_Vector y, N_Vector fy,
               SUNMatrix J, void *udata,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* User provided routine called by the solver to compute
 * r.h.s. sensititivy. */
static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot,
              int iS, N_Vector yS, N_Vector ySdot,
              void *udata, N_Vector tmp1, N_Vector tmp2);

/* Prototypes of private functions */
static int runCVode(void *cvode_mem, N_Vector y, N_Vector *yS, UserData data);
static void PrintHeader(UserData data);
static int PrintFinalStats(void *cvode_mem, UserData data);
static int check_retval(void *returnvalue, const char *funcname, int opt);


/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;
  void *cvode_mem;

  realtype reltol;
  N_Vector y0, y, abstol;

  int Ns;
  realtype *pbar;
  int is, *plist, retval;
  N_Vector *yS0, *yS;

  SUNMatrix A;
  SUNLinearSolver LS;

  /* Allocate user data structure */
  data = (UserData) malloc(sizeof *data);

  /* Initialize sensitivity variables (reaction rates for this problem) */
  data->p[0] = RCONST(0.04);
  data->p[1] = RCONST(1.0e4);
  data->p[2] = RCONST(3.0e7);

  /* Allocate vectors */
  y0 = N_VNew_Serial(NEQ);      /* initial conditions  */
  y = N_VNew_Serial(NEQ);       /* solution vector     */
  abstol = N_VNew_Serial(NEQ);  /* absolute tolerances */

  /* Set initial conditions */
  NV_Ith_S(y0,0) = RCONST(1.0);
  NV_Ith_S(y0,1) = RCONST(0.0);
  NV_Ith_S(y0,2) = RCONST(0.0);

  /* Set integration tolerances */
  reltol = RCONST(1e-6);
  NV_Ith_S(abstol,0) = RCONST(1e-8);
  NV_Ith_S(abstol,1) = RCONST(1e-14);
  NV_Ith_S(abstol,2) = RCONST(1e-6);

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function y'=f(t,y), the initial time T0, and
   * the intiial dependenet variable vector y0. */
  retval = CVodeInit(cvode_mem, f, T0, y0);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSVTolerances to specify the scalar relative tolerance
   * and vector absolute tolereance */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

  /* Call CVodeSetUserData so the sensitivity params can be accessed
   * from user provided routines. */
  retval = CVodeSetUserData(cvode_mem, data);
  if (check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeSetMaxNumSteps to set the maximum number of steps the
   * solver will take in an attempt to reach the next output time. */
  retval = CVodeSetMaxNumSteps(cvode_mem, MXSTEPS);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solvers */
  A = SUNDenseMatrix(NEQ, NEQ);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(y, A);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* Specifiy the Jacobian approximation routine to be used */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if (check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

  /* Sensitivity-related settings */
  data->sensi   = SUNTRUE;          /* sensitivity ON                */
  data->meth    = CV_SIMULTANEOUS;  /* simultaneous corrector method */
  data->errconS = SUNTRUE;          /* full error control            */
  data->fsDQ    = SUNFALSE;         /* user-provided sensitvity RHS  */

  Ns = 3;

  pbar = (realtype *) malloc(Ns * sizeof(realtype));
  pbar[0] = data->p[0];
  pbar[1] = data->p[1];
  pbar[2] = data->p[2];

  plist = (int *) malloc(Ns * sizeof(int));
  for (is=0; is<Ns; is++) plist[is] = is;

  yS0 = N_VCloneVectorArray_Serial(Ns, y);
  for (is=0;is<Ns;is++) N_VConst(ZERO, yS0[is]);

  yS = N_VCloneVectorArray_Serial(Ns, y);

  retval = CVodeSensInit1(cvode_mem, Ns, data->meth, fS, yS0);
  if (check_retval(&retval, "CVodeSensInit1", 1)) return(1);

  retval = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);
  if (check_retval(&retval, "CVodeSetSensParams", 1)) return(1);

  /*
    Sensitivities are enabled
    Set full error control
    Set user-provided sensitivity RHS
    Run CVODES
  */

  retval = CVodeSensEEtolerances(cvode_mem);
  if (check_retval(&retval, "CVodeSensEEtolerances", 1)) return(1);

  retval = CVodeSetSensErrCon(cvode_mem, data->errconS);
  if (check_retval(&retval, "CVodeSetSensErrCon", 1)) return(1);

  retval = runCVode(cvode_mem, y, yS, data);
  if (check_retval(&retval, "runCVode", 1)) return(1);

  /*
    Change parameters
    Toggle sensitivities OFF
    Reinitialize and run CVODES
  */

  data->p[0] = RCONST(0.05);
  data->p[1] = RCONST(2.0e4);
  data->p[2] = RCONST(2.9e7);

  data->sensi = SUNFALSE;

  retval = CVodeReInit(cvode_mem, T0, y0);
  if (check_retval(&retval, "CVodeReInit", 1)) return(1);

  retval = CVodeSensToggleOff(cvode_mem);
  if (check_retval(&retval, "CVodeSensToggleOff", 1)) return(1);

  retval = runCVode(cvode_mem, y, yS, data);
  if (check_retval(&retval, "runCVode", 1)) return(1);

  /*
    Change parameters
    Switch to internal DQ sensitivity RHS function
    Toggle sensitivities ON (reinitialize sensitivities)
    Reinitialize and run CVODES
  */

  data->p[0] = RCONST(0.06);
  data->p[1] = RCONST(3.0e4);
  data->p[2] = RCONST(2.8e7);

  data->sensi = SUNTRUE;
  data->fsDQ  = SUNTRUE;

  retval = CVodeReInit(cvode_mem, T0, y0);
  if (check_retval(&retval, "CVodeReInit", 1)) return(1);

  CVodeSensFree(cvode_mem);
  retval = CVodeSensInit1(cvode_mem, Ns, data->meth, NULL, yS0);
  if (check_retval(&retval, "CVodeSensInit1", 1)) return(1);

  retval = runCVode(cvode_mem, y, yS, data);
  if (check_retval(&retval, "runCVode", 1)) return(1);

  /*
    Switch to partial error control
    Switch back to user-provided sensitivity RHS
    Toggle sensitivities ON (reinitialize sensitivities)
    Change method to staggered
    Reinitialize and run CVODES
  */

  data->sensi   = SUNTRUE;
  data->errconS = SUNFALSE;
  data->fsDQ    = SUNFALSE;
  data->meth    = CV_STAGGERED;

  retval = CVodeReInit(cvode_mem, T0, y0);
  if (check_retval(&retval, "CVodeReInit", 1)) return(1);

  retval = CVodeSetSensErrCon(cvode_mem, data->errconS);
  if (check_retval(&retval, "CVodeSetSensErrCon", 1)) return(1);

  CVodeSensFree(cvode_mem);
  retval = CVodeSensInit1(cvode_mem, Ns, data->meth, fS, yS0);
  if (check_retval(&retval, "CVodeSensInit1", 1)) return(1);

  retval = runCVode(cvode_mem, y, yS, data);
  if (check_retval(&retval, "runCVode", 1)) return(1);

  /*
    Free sensitivity-related memory
    (CVodeSensToggle is not needed, as CVodeSensFree toggles sensitivities OFF)
    Reinitialize and run CVODES
  */

  data->sensi = SUNFALSE;

  CVodeSensFree(cvode_mem);

  retval = CVodeReInit(cvode_mem, T0, y0);
  if (check_retval(&retval, "CVodeReInit", 1)) return(1);

  retval = runCVode(cvode_mem, y, yS, data);
  if (check_retval(&retval, "runCVode", 1)) return(1);

  /* Free memory */

  N_VDestroy(y0);                 /* Free y0 vector         */
  N_VDestroy(y);                  /* Free y vector          */
  N_VDestroy(abstol);             /* Free abstol vector     */
  N_VDestroyVectorArray(yS0, Ns); /* Free yS0 vector        */
  N_VDestroyVectorArray(yS, Ns);  /* Free yS vector         */
  free(plist);                    /* Free plist             */
  free(pbar);                     /* Free pbar              */
  free(data);                     /* Free user data         */
  CVodeFree(&cvode_mem);          /* Free integrator memory */
  SUNLinSolFree(LS);              /* Free solver memory     */
  SUNMatDestroy(A);               /* Free the matrix memory */

  return(0);

}

/*
 * Runs integrator and prints final statistics when complete.
 */

static int runCVode(void *cvode_mem, N_Vector y, N_Vector *yS, UserData data)
{
  realtype t;
  int retval;

  /* Print header for current run */
  PrintHeader(data);

  /* Call CVode in CV_NORMAL mode */
  retval = CVode(cvode_mem, T1, y, &t, CV_NORMAL);
  if (retval != 0) return(retval);
  
  /* Print final statistics */
  retval = PrintFinalStats(cvode_mem, data);
  printf("\n");

  return(retval);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY THE SOLVER
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y).
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *udata)
{
  realtype y1, y2, y3, yd1, yd3;
  UserData data;
  realtype p1, p2, p3;

  y1 = NV_Ith_S(y,0); y2 = NV_Ith_S(y,1); y3 = NV_Ith_S(y,2);
  data = (UserData) udata;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  yd1 = NV_Ith_S(ydot,0) = -p1*y1 + p2*y2*y3;
  yd3 = NV_Ith_S(ydot,2) = p3*y2*y2;
        NV_Ith_S(ydot,1) = -yd1 - yd3;

  return(0);
}


/*
 * Jacobian routine. Compute J(t,y).
 */

static int Jac(realtype t, N_Vector y, N_Vector fy,
               SUNMatrix J, void *udata,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y2, y3;
  UserData data;
  realtype p1, p2, p3;

  y2 = NV_Ith_S(y,1); y3 = NV_Ith_S(y,2);
  data = (UserData) udata;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  SM_ELEMENT_D(J,0,0) = -p1;  SM_ELEMENT_D(J,0,1) = p2*y3;          SM_ELEMENT_D(J,0,2) = p2*y2;
  SM_ELEMENT_D(J,1,0) =  p1;  SM_ELEMENT_D(J,1,1) = -p2*y3-2*p3*y2; SM_ELEMENT_D(J,1,2) = -p2*y2;
                      SM_ELEMENT_D(J,2,1) = 2*p3*y2;

  return(0);
}

/*
 * fS routine. Compute sensitivity r.h.s.
 */

static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot,
              int iS, N_Vector yS, N_Vector ySdot,
              void *udata, N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;
  realtype s1, s2, s3;
  realtype sd1, sd2, sd3;

  data = (UserData) udata;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  y1 = NV_Ith_S(y,0);  y2 = NV_Ith_S(y,1);  y3 = NV_Ith_S(y,2);
  s1 = NV_Ith_S(yS,0); s2 = NV_Ith_S(yS,1); s3 = NV_Ith_S(yS,2);

  sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
  sd3 = 2*p3*y2*s2;
  sd2 = -sd1-sd3;

  switch (iS) {
  case 0:
    sd1 += -y1;
    sd2 +=  y1;
    break;
  case 1:
    sd1 +=  y2*y3;
    sd2 += -y2*y3;
    break;
  case 2:
    sd2 += -y2*y2;
    sd3 +=  y2*y2;
    break;
  }

  NV_Ith_S(ySdot,0) = sd1;
  NV_Ith_S(ySdot,1) = sd2;
  NV_Ith_S(ySdot,2) = sd3;

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

static void PrintHeader(UserData data)
{
  /* Print sensitivity control retvals */
  printf("Sensitivity: ");
  if (data->sensi) {
    printf("YES (");
    switch (data->meth) {
    case CV_SIMULTANEOUS:
      printf("SIMULTANEOUS + ");
      break;
    case CV_STAGGERED:
      printf("STAGGERED + ");
      break;
    case CV_STAGGERED1:
      printf("STAGGERED-1 + ");
      break;
    }
    if (data->errconS) printf("FULL ERROR CONTROL + ");
    else         printf("PARTIAL ERROR CONTROL + ");
    if (data->fsDQ)    printf("DQ sensitivity RHS)\n");
    else         printf("user-provided sensitivity RHS)\n");
  } else {
    printf("NO\n");
  }

  /* Print current problem parameters */
  printf("Parameters: [%8.4e  %8.4e  %8.4e]\n",data->p[0], data->p[1], data->p[2]);
}

/*
 * Print some final statistics from the CVODES memory.
 */

static int PrintFinalStats(void *cvode_mem, UserData data)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int njeD, nfeD;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

  if (data->sensi) {
    retval = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    retval = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    retval = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    if (data->errconS)
      retval = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
    else
      netfS = 0;
    if (data->meth == CV_STAGGERED) {
      retval = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
      retval = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
    } else {
      nniS = 0;
      ncfnS = 0;
    }
  }

  retval = CVodeGetNumJacEvals(cvode_mem, &njeD);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeD);

  printf("Run statistics:\n");

  printf("   nst     = %5ld\n", nst);
  printf("   nfe     = %5ld\n",   nfe);
  printf("   netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("   nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  printf("   njeD    = %5ld    nfeD     = %5ld\n", njeD, nfeD);

  if(data->sensi) {
    printf("   -----------------------------------\n");  /* simultaneous corrector method */
    printf("   nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("   netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("   nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }

  return(retval);
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
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

