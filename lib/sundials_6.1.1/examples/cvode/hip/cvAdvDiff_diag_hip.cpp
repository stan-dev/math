/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This program solves the problem with the ADAMS integration method,
 * and with Newton iteration using diagonal approximate Jacobians.
 * It can use scalar (default) relative and absolute tolerances or a
 * vector of absolute tolerances (controlled by a runtime argument).
 * The constraint u_i >= 0 is posed for all components.
 * Output is printed at t = .5, 1.0, ..., 5.
 * Run statistics (optional outputs) are printed at the end.
 *
 * ./cvAdvDiff_diag_hip [0 (scalar atol) | 1 (vector atol)]
 *                       [0 (unfused) | 1 (fused)]
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <hip/hip_runtime.h>

#include <cvode/cvode.h>              /* prototypes for CVODE fcts., consts.  */
#include <cvode/cvode_diag.h>         /* prototypes for CVODE diagonal solver */
#include <nvector/nvector_hip.h>      /* access to hip N_Vector               */
#include <sundials/sundials_types.h>  /* definition of type realtype          */

/* Problem Constants */

#define ZERO  RCONST(0.0)

#define XMAX  RCONST(2.0)    /* domain boundary           */
#define MX    10             /* mesh dimension            */
#define NEQ   MX             /* number of equations       */
#define ATOL  RCONST(1e-10)  /* scalar absolute tolerance */
#define T0    ZERO           /* initial time              */
#define T1    RCONST(0.5)    /* first output time         */
#define DTOUT RCONST(0.5)    /* output time increment     */
#define NOUT  10             /* number of output times    */

/* Type : UserData
   contains mesh spacing and problem parameters. */

typedef struct {
  realtype dx;
  realtype hdcoef;
  realtype hacoef;
} *UserData;

/* Private Helper Functions */

static void SetIC(N_Vector u, realtype dx);

static void PrintIntro(int toltype, int usefused);

static void PrintData(realtype t, realtype umax, long int nst);

static void PrintFinalStats(void *cvode_mem);

/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

/***************************** Main Program ******************************/

int main(int argc, char *argv[])
{
  sundials::Context sunctx;
  realtype dx, reltol, abstol, t, tout, umax;
  N_Vector u;
  UserData data;
  void *cvode_mem;
  int iout, retval, toltype, usefused;
  long int nst;

  u = NULL;
  data = NULL;
  cvode_mem = NULL;
  toltype = 0;
  usefused = 0;

  if (argc >= 2) {
    /* use vector or scalar atol? */
    toltype = atoi(argv[1]);
    /* use fused operations? */
    if (argc == 3)
      usefused = atoi(argv[2]);
  }

  data = (UserData) malloc(sizeof *data);  /* Allocate data memory */
  if(check_retval((void *)data, "malloc", 2)) return 1;

  u = N_VNew_Hip(NEQ, sunctx);  /* Allocate u vector */
  if(check_retval((void *)u, "N_VNew", 0)) return 1;

  reltol = ZERO;  /* Set the tolerances */
  abstol = ATOL;

  dx = data->dx = XMAX/((realtype)(MX+1));  /* Set grid coefficients in data */
  data->hdcoef = RCONST(1.0)/(dx*dx);
  data->hacoef = RCONST(0.5)/(RCONST(2.0)*dx);

  SetIC(u, dx);  /* Initialize u vector */

  /* Call CVodeCreate to create the solver memory and specify the
   * Adams-Moulton LMM */
  cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return 1;

  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return 1;

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */

  if (toltype == 0) {
    retval = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_retval(&retval, "CVodeSStolerances", 1)) return(1);
  } else {
    N_Vector vabstol = N_VClone_Hip(u);
    if (check_retval(&vabstol, "N_VClone_Hip", 0)) return(1);
    N_VConst(abstol, vabstol);
    retval = CVodeSVtolerances(cvode_mem, reltol, vabstol);
    if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);
    N_VDestroy(vabstol);
  }

  /* Call CVDiag to create and attach CVODE-specific diagonal linear solver */
  retval = CVDiag(cvode_mem);
  if(check_retval(&retval, "CVDiag", 1)) return(1);

  /* Tell CVode to use fused kernels if they are available. */
  retval = CVodeSetUseIntegratorFusedKernels(cvode_mem, usefused);
  check_retval(&retval, "CVodeSetUseIntegratorFusedKernels", 1);

  PrintIntro(toltype, usefused);

  umax = N_VMaxNorm(u);

  t = T0;
  PrintData(t, umax, 0);

  /* In loop over output points, call CVode, print results, test for error */

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
    retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if(check_retval(&retval, "CVode", 1)) break;
    umax = N_VMaxNorm(u);
    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_retval(&retval, "CVodeGetNumSteps", 1);
    PrintData(t, umax, nst);
  }

  PrintFinalStats(cvode_mem);  /* Print some final statistics */

  N_VDestroy(u);                 /* Free the u vector */
  CVodeFree(&cvode_mem);         /* Free the integrator memory */
  free(data);                    /* Free user data */

  return(0);
}

/************************ Private Helper Functions ***********************/

/* Set initial conditions in u vector */

static void SetIC(N_Vector u, realtype dx)
{
  int i;
  sunindextype N;
  realtype x;
  realtype *udata;

  /* Set pointer to data array and get local length of u. */
  udata = N_VGetHostArrayPointer_Hip(u);
  N = N_VGetLength(u);

  /* Load initial profile into u vector */
  for (i=1; i<=N; i++) {
    x = i*dx;
    udata[i-1] = x*(XMAX - x)*exp(RCONST(2.0)*x);
  }
  N_VCopyToDevice_Hip(u);
}

/* Print problem introduction */

static void PrintIntro(int toltype, int usefused)
{
  printf("\n 1-D advection-diffusion equation, mesh size =%3d \n", MX);
  printf("\n Diagonal linear solver CVDiag \n");
  if (usefused)
    printf(" Using fused CVODE kernels \n");
  if (toltype == 0)
    printf(" Using scalar ATOL\n");
  else
    printf(" Using vector ATOL\n");
  printf("\n");

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
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  printf("\nFinal Statistics: \n\n");
  printf("nst = %-6ld  nfe  = %-6ld  ", nst, nfe);
  printf("nni = %-6ld  ncfn = %-6ld  netf = %ld\n \n", nni, ncfn, netf);
}

 /***************** Function Called by the Solver ***********************/

 /* f routine. Compute f(t,u). */

__global__
static void f_kernel(sunindextype N,
                     realtype hordc, realtype horac,
                     const realtype* u, realtype* udot)
{
  sunindextype i = blockDim.x*blockIdx.x + threadIdx.x;
  realtype ui, ult, urt, hdiff, hadv;

  if (i < N) {
    /* Extract u at x_i and two neighboring points */
    ui = u[i];
    ult = (i == 0) ? ZERO : u[i-1];
    urt = (i == N-1) ? ZERO : u[i+1];

    /* Set diffusion and advection terms and load into udot */
    hdiff = hordc*(ult - RCONST(2.0)*ui + urt);
    hadv = horac*(urt - ult);
    udot[i] = hdiff + hadv;
  }
}

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype hordc, horac;
  realtype *udata, *dudata;
  sunindextype N;
  size_t grid, block;
  UserData data;
  hipError_t cuerr;

  udata = N_VGetDeviceArrayPointer_Hip(u);
  dudata = N_VGetDeviceArrayPointer_Hip(udot);

  /* Extract needed problem constants from data */
  data = (UserData) user_data;
  hordc = data->hdcoef;
  horac = data->hacoef;

  /* Extract parameters for parallel computation. */
  N = N_VGetLength(u); /* Number of elements of u. */

  block = 64;
  grid  = (block + N - 1)/block;
  f_kernel<<<grid, block>>>(N, hordc, horac, udata, dudata);

  hipDeviceSynchronize();
  cuerr = hipGetLastError();
  if (cuerr != hipSuccess) {
    fprintf(stderr, "ERROR in f: f_kernel --> %s\n", hipGetErrorString(cuerr));
    return(-1);
  }

  return(0);
}

 /* Check function return value...
      opt == 0 means SUNDIALS function allocates memory so check if
               returned NULL pointer
      opt == 1 means SUNDIALS function returns an integer value so check if
               retval < 0
      opt == 2 means function allocates memory so check if returned
               NULL pointer */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n", funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1); }

  return(0);
}
