/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
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
 * This example loops through the available iterative linear solvers:
 * SPGMR, SPBCG and SPTFQMR.
 *
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
 * 10 x 10 mesh, with simple polynomial initial profiles.
 * The problem is solved with CVODE, with the BDF/GMRES,
 * BDF/Bi-CGStab, and BDF/TFQMR methods (i.e. using the SUNLinSol_SPGMR,
 * SUNLinSol_SPBCGS and SUNLinSol_SPTFQMR linear solvers) and the
 * block-diagonal part of the Newton matrix as a left preconditioner.
 * A copy of the block-diagonal part of the Jacobian is saved and
 * conditionally reused within the Precond routine.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvode/cvode.h>                 /* main integrator header file       */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver   */
#include <sunlinsol/sunlinsol_spbcgs.h>  /* access to SPBCGS SUNLinearSolver  */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to SPTFQMR SUNLinearSolver */
#include <nvector/nvector_serial.h>      /* serial N_Vector types, fct. and macros */
#include <sundials/sundials_dense.h>     /* use generic DENSE solver in preconditioning */
#include <sundials/sundials_types.h>     /* definition of realtype */
#include <sundials/sundials_math.h>      /* contains the macros ABS, SUNSQR, and EXP */

/* Problem Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

#define NUM_SPECIES  2                 /* number of species         */
#define KH           RCONST(4.0e-6)    /* horizontal diffusivity Kh */
#define VEL          RCONST(0.001)     /* advection velocity V      */
#define KV0          RCONST(1.0e-8)    /* coefficient in Kv(y)      */
#define Q1           RCONST(1.63e-16)  /* coefficients q1, q2, c3   */ 
#define Q2           RCONST(4.66e-16)
#define C3           RCONST(3.7e16)
#define A3           RCONST(22.62)     /* coefficient in expression for q3(t) */
#define A4           RCONST(7.601)     /* coefficient in expression for q4(t) */
#define C1_SCALE     RCONST(1.0e6)     /* coefficients in initial profiles    */
#define C2_SCALE     RCONST(1.0e12)

#define T0           ZERO                 /* initial time */
#define NOUT         12                   /* number of output times */
#define TWOHR        RCONST(7200.0)       /* number of seconds in two hours  */
#define HALFDAY      RCONST(4.32e4)       /* number of seconds in a half day */
#define PI       RCONST(3.1415926535898)  /* pi */ 

#define XMIN         ZERO                 /* grid boundaries in x  */
#define XMAX         RCONST(20.0)           
#define YMIN         RCONST(30.0)         /* grid boundaries in y  */
#define YMAX         RCONST(50.0)
#define XMID         RCONST(10.0)         /* grid midpoints in x,y */          
#define YMID         RCONST(40.0)

#define MX           10             /* MX = number of x mesh points */
#define MY           10             /* MY = number of y mesh points */
#define NSMX         20             /* NSMX = NUM_SPECIES*MX */
#define MM           (MX*MY)        /* MM = MX*MY */

/* CVodeInit Constants */

#define RTOL    RCONST(1.0e-5)    /* scalar relative tolerance */
#define FLOOR   RCONST(100.0)     /* value of C1 or C2 at which tolerances */
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */
#define NEQ     (NUM_SPECIES*MM)  /* NEQ = number of equations */

/* Linear Solver Loop Constants */

#define USE_SPGMR   0
#define USE_SPBCG   1
#define USE_SPTFQMR 2

/* User-defined vector and matrix accessor macros: IJKth, IJth */

/* IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.   
   
   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MY-1. The vdata array is obtained via
   the call vdata = N_VGetArrayPointer(v), where v is an N_Vector. 
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.

   IJth(a,i,j) references the (i,j)th entry of the matrix realtype **a,
   where 1 <= i,j <= NUM_SPECIES. The small matrix routines in 
   sundials_dense.h work with matrices stored by column in a 2-dimensional 
   array. In C, arrays are indexed starting at 0, not 1. */

#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMX])
#define IJth(a,i,j)        (a[j-1][i-1])

/* Type : UserData 
   contains preconditioner blocks, pivot arrays, and problem constants */

typedef struct {
  realtype **P[MX][MY], **Jbd[MX][MY];
  sunindextype *pivot[MX][MY];
  realtype q4, om, dx, dy, hdco, haco, vdco;
} *UserData;

/* Private Helper Functions */

static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector u, realtype dx, realtype dy);
static void PrintOutput(void *cvode_mem, N_Vector u, realtype t);
static void PrintFinalStats(void *cvode_mem, int linsolver);
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Functions Called by the Solver */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data);

static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, void *user_data);

static int PSolve(realtype tn, N_Vector u, N_Vector fu,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *user_data);


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(void)
{
  realtype abstol, reltol, t, tout;
  N_Vector u;
  UserData data;
  SUNLinearSolver LS;
  void *cvode_mem;
  int linsolver, iout, retval;

  u = NULL;
  data = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Allocate memory, and set problem data, initial values, tolerances */ 
  u = N_VNew_Serial(NEQ);
  if(check_retval((void *)u, "N_VNew_Serial", 0)) return(1);
  data = AllocUserData();
  if(check_retval((void *)data, "AllocUserData", 2)) return(1);
  InitUserData(data);
  SetInitialProfiles(u, data->dx, data->dy);
  abstol=ATOL; 
  reltol=RTOL;

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Set the pointer to user-defined data */
  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in u'=f(t,u), the inital time T0, and
   * the initial dependent variable vector u. */
  retval = CVodeInit(cvode_mem, f, T0, u);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* START: Loop through SPGMR, SPBCG and SPTFQMR linear solver modules */
  for (linsolver = 0; linsolver < 3; ++linsolver) {

    if (linsolver != 0) {

      /* Re-initialize user data */
      InitUserData(data);
      SetInitialProfiles(u, data->dx, data->dy);

    /* Re-initialize CVode for the solution of the same problem, but
       using a different linear solver module */
      retval = CVodeReInit(cvode_mem, T0, u);
      if (check_retval(&retval, "CVodeReInit", 1)) return(1);

    }

    /* Free previous linear solver and attach a new linear solver module */
    SUNLinSolFree(LS);

    switch(linsolver) {

    /* (a) SPGMR */
    case(USE_SPGMR):

      /* Print header */
      printf(" -------");
      printf(" \n| SPGMR |\n");
      printf(" -------\n");

      /* Call SUNLinSol_SPGMR to specify the linear solver SPGMR with
         left preconditioning and the default maximum Krylov dimension */
      LS = SUNLinSol_SPGMR(u, PREC_LEFT, 0);
      if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

      retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

      break;

    /* (b) SPBCG */
    case(USE_SPBCG):

      /* Print header */
      printf(" -------");
      printf(" \n| SPBCGS |\n");
      printf(" -------\n");

      /* Call SUNLinSol_SPBCGS to specify the linear solver SPBCGS with
         left preconditioning and the default maximum Krylov dimension */
      LS = SUNLinSol_SPBCGS(u, PREC_LEFT, 0);
      if(check_retval((void *)LS, "SUNLinSol_SPBCGS", 0)) return(1);

      retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

      break;

    /* (c) SPTFQMR */
    case(USE_SPTFQMR):

      /* Print header */
      printf(" ---------");
      printf(" \n| SPTFQMR |\n");
      printf(" ---------\n");

      /* Call SUNLinSol_SPTFQMR to specify the linear solver SPTFQMR with
         left preconditioning and the default maximum Krylov dimension */
      LS = SUNLinSol_SPTFQMR(u, PREC_LEFT, 0);
      if(check_retval((void *)LS, "SUNLinSol_SPTFQMR", 0)) return(1);

      retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
      if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

      break;

    }


    /* Set preconditioner setup and solve routines Precond and PSolve,
       and the pointer to the user-defined block data */
    retval = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
    if(check_retval(&retval, "CVodeSetPreconditioner", 1)) return(1);

    /* In loop over output points, call CVode, print results, test for error */
    printf(" \n2-species diurnal advection-diffusion problem\n\n");
    for (iout=1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
      retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      PrintOutput(cvode_mem, u, t);
      if(check_retval(&retval, "CVode", 1)) break;
    }

    PrintFinalStats(cvode_mem, linsolver);

  }  /* END: Loop through SPGMR, SPBCG and SPTFQMR linear solver modules */

  /* Free memory */
  N_VDestroy(u);
  FreeUserData(data);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);

  return(0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

/* Allocate memory for data structure of type UserData */

static UserData AllocUserData(void)
{
  int jx, jy;
  UserData data;

  data = (UserData) malloc(sizeof *data);

  for (jx=0; jx < MX; jx++) {
    for (jy=0; jy < MY; jy++) {
      (data->P)[jx][jy] = newDenseMat(NUM_SPECIES, NUM_SPECIES);
      (data->Jbd)[jx][jy] = newDenseMat(NUM_SPECIES, NUM_SPECIES);
      (data->pivot)[jx][jy] = newIndexArray(NUM_SPECIES);
    }
  }

  return(data);
}

/* Load problem constants in data */

static void InitUserData(UserData data)
{
  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/(MX-1);
  data->dy = (YMAX-YMIN)/(MY-1);
  data->hdco = KH/SUNSQR(data->dx);
  data->haco = VEL/(TWO*data->dx);
  data->vdco = (ONE/SUNSQR(data->dy))*KV0;
}

/* Free data memory */

static void FreeUserData(UserData data)
{
  int jx, jy;

  for (jx=0; jx < MX; jx++) {
    for (jy=0; jy < MY; jy++) {
      destroyMat((data->P)[jx][jy]);
      destroyMat((data->Jbd)[jx][jy]);
      destroyArray((data->pivot)[jx][jy]);
    }
  }

  free(data);
}

/* Set initial conditions in u */

static void SetInitialProfiles(N_Vector u, realtype dx, realtype dy)
{
  int jx, jy;
  realtype x, y, cx, cy;
  realtype *udata;

  /* Set pointer to data array in vector u. */

  udata = N_VGetArrayPointer(u);

  /* Load initial profiles of c1 and c2 into u vector */

  for (jy=0; jy < MY; jy++) {
    y = YMIN + jy*dy;
    cy = SUNSQR(RCONST(0.1)*(y - YMID));
    cy = ONE - cy + RCONST(0.5)*SUNSQR(cy);
    for (jx=0; jx < MX; jx++) {
      x = XMIN + jx*dx;
      cx = SUNSQR(RCONST(0.1)*(x - XMID));
      cx = ONE - cx + RCONST(0.5)*SUNSQR(cx);
      IJKth(udata,1,jx,jy) = C1_SCALE*cx*cy; 
      IJKth(udata,2,jx,jy) = C2_SCALE*cx*cy;
    }
  }
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(void *cvode_mem, N_Vector u, realtype t)
{
  long int nst;
  int qu, retval;
  realtype hu, *udata;
  int mxh = MX/2 - 1, myh = MY/2 - 1, mx1 = MX - 1, my1 = MY - 1;

  udata = N_VGetArrayPointer(u);

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetLastOrder(cvode_mem, &qu);
  check_retval(&retval, "CVodeGetLastOrder", 1);
  retval = CVodeGetLastStep(cvode_mem, &hu);
  check_retval(&retval, "CVodeGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("t = %.2Le   no. steps = %ld   order = %d   stepsize = %.2Le\n",
         t, nst, qu, hu);
  printf("c1 (bot.left/middle/top rt.) = %12.3Le  %12.3Le  %12.3Le\n",
         IJKth(udata,1,0,0), IJKth(udata,1,mxh,myh), IJKth(udata,1,mx1,my1));
  printf("c2 (bot.left/middle/top rt.) = %12.3Le  %12.3Le  %12.3Le\n\n",
         IJKth(udata,2,0,0), IJKth(udata,2,mxh,myh), IJKth(udata,2,mx1,my1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("t = %.2e   no. steps = %ld   order = %d   stepsize = %.2e\n",
         t, nst, qu, hu);
  printf("c1 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n",
         IJKth(udata,1,0,0), IJKth(udata,1,mxh,myh), IJKth(udata,1,mx1,my1));
  printf("c2 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n\n",
         IJKth(udata,2,0,0), IJKth(udata,2,mxh,myh), IJKth(udata,2,mx1,my1));
#else
  printf("t = %.2e   no. steps = %ld   order = %d   stepsize = %.2e\n",
         t, nst, qu, hu);
  printf("c1 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n",
         IJKth(udata,1,0,0), IJKth(udata,1,mxh,myh), IJKth(udata,1,mx1,my1));
  printf("c2 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n\n",
         IJKth(udata,2,0,0), IJKth(udata,2,mxh,myh), IJKth(udata,2,mx1,my1));
#endif
}

/* Get and print final statistics */

static void PrintFinalStats(void *cvode_mem, int linsolver)
{
  long int lenrw, leniw ;
  long int lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int retval;

  retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  check_retval(&retval, "CVodeGetWorkSpace", 1);
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

  retval = CVodeGetLinWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
  check_retval(&retval, "CVodeGetLinWorkSpace", 1);
  retval = CVodeGetNumLinIters(cvode_mem, &nli);
  check_retval(&retval, "CVodeGetNumLinIters", 1);
  retval = CVodeGetNumPrecEvals(cvode_mem, &npe);
  check_retval(&retval, "CVodeGetNumPrecEvals", 1);
  retval = CVodeGetNumPrecSolves(cvode_mem, &nps);
  check_retval(&retval, "CVodeGetNumPrecSolves", 1);
  retval = CVodeGetNumLinConvFails(cvode_mem, &ncfl);
  check_retval(&retval, "CVodeGetNumLinConvFails", 1);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  printf("\nFinal Statistics.. \n\n");
  printf("lenrw   = %5ld     leniw   = %5ld\n"  , lenrw, leniw);
  printf("lenrwLS = %5ld     leniwLS = %5ld\n"  , lenrwLS, leniwLS);
  printf("nst     = %5ld\n"                     , nst);
  printf("nfe     = %5ld     nfeLS   = %5ld\n"  , nfe, nfeLS);
  printf("nni     = %5ld     nli     = %5ld\n"  , nni, nli);
  printf("nsetups = %5ld     netf    = %5ld\n"  , nsetups, netf);
  printf("npe     = %5ld     nps     = %5ld\n"  , npe, nps);
  printf("ncfn    = %5ld     ncfl    = %5ld\n\n", ncfn, ncfl);

  if (linsolver < 2)
    printf("======================================================================\n\n");
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

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* f routine. Compute RHS function f(t,u). */

static int f(realtype t, N_Vector u, N_Vector udot, void *user_data)
{
  realtype q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  realtype c1rt, c2rt, cydn, cyup, hord1, hord2, horad1, horad2;
  realtype qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup;
  realtype q4coef, dely, verdco, hordco, horaco;
  realtype *udata, *dudata;
  int jx, jy, idn, iup, ileft, iright;
  UserData data;

  data   = (UserData) user_data;
  udata  = N_VGetArrayPointer(u);
  dudata = N_VGetArrayPointer(udot);

  /* Set diurnal rate coefficients. */

  s = sin(data->om*t);
  if (s > ZERO) {
    q3 = SUNRexp(-A3/s);
    data->q4 = SUNRexp(-A4/s);
  } else {
      q3 = ZERO;
      data->q4 = ZERO;
  }

  /* Make local copies of problem variables, for efficiency. */

  q4coef = data->q4;
  dely = data->dy;
  verdco = data->vdco;
  hordco  = data->hdco;
  horaco  = data->haco;

  /* Loop over all grid points. */

  for (jy=0; jy < MY; jy++) {

    /* Set vertical diffusion coefficients at jy +- 1/2 */

    ydn = YMIN + (jy - RCONST(0.5))*dely;
    yup = ydn + dely;
    cydn = verdco*SUNRexp(RCONST(0.2)*ydn);
    cyup = verdco*SUNRexp(RCONST(0.2)*yup);
    idn = (jy == 0) ? 1 : -1;
    iup = (jy == MY-1) ? -1 : 1;
    for (jx=0; jx < MX; jx++) {

      /* Extract c1 and c2, and set kinetic rate terms. */

      c1 = IJKth(udata,1,jx,jy); 
      c2 = IJKth(udata,2,jx,jy);
      qq1 = Q1*c1*C3;
      qq2 = Q2*c1*c2;
      qq3 = q3*C3;
      qq4 = q4coef*c2;
      rkin1 = -qq1 - qq2 + TWO*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms. */

      c1dn = IJKth(udata,1,jx,jy+idn);
      c2dn = IJKth(udata,2,jx,jy+idn);
      c1up = IJKth(udata,1,jx,jy+iup);
      c2up = IJKth(udata,2,jx,jy+iup);
      vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
      vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);

      /* Set horizontal diffusion and advection terms. */

      ileft = (jx == 0) ? 1 : -1;
      iright =(jx == MX-1) ? -1 : 1;
      c1lt = IJKth(udata,1,jx+ileft,jy); 
      c2lt = IJKth(udata,2,jx+ileft,jy);
      c1rt = IJKth(udata,1,jx+iright,jy);
      c2rt = IJKth(udata,2,jx+iright,jy);
      hord1 = hordco*(c1rt - TWO*c1 + c1lt);
      hord2 = hordco*(c2rt - TWO*c2 + c2lt);
      horad1 = horaco*(c1rt - c1lt);
      horad2 = horaco*(c2rt - c2lt);

      /* Load all terms into udot. */

      IJKth(dudata, 1, jx, jy) = vertd1 + hord1 + horad1 + rkin1; 
      IJKth(dudata, 2, jx, jy) = vertd2 + hord2 + horad2 + rkin2;
    }
  }

  return(0);
}

/* Preconditioner setup routine. Generate and preprocess P. */

static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, void *user_data)
{
  realtype c1, c2, cydn, cyup, diag, ydn, yup, q4coef, dely, verdco, hordco;
  realtype **(*P)[MY], **(*Jbd)[MY];
  sunindextype *(*pivot)[MY], retval;
  int jx, jy;
  realtype *udata, **a, **j;
  UserData data;
  
  /* Make local copies of pointers in user_data, and of pointer to u's data */
  
  data = (UserData) user_data;
  P = data->P;
  Jbd = data->Jbd;
  pivot = data->pivot;
  udata = N_VGetArrayPointer(u);
  
  if (jok) {
    
    /* jok = SUNTRUE: Copy Jbd to P */
    
    for (jy=0; jy < MY; jy++)
      for (jx=0; jx < MX; jx++)
        denseCopy(Jbd[jx][jy], P[jx][jy], NUM_SPECIES, NUM_SPECIES);
    
    *jcurPtr = SUNFALSE;
    
  }
  
  else {
    /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */
    
    /* Make local copies of problem variables, for efficiency. */
    
    q4coef = data->q4;
    dely = data->dy;
    verdco = data->vdco;
    hordco  = data->hdco;
    
    /* Compute 2x2 diagonal Jacobian blocks (using q4 values 
       computed on the last f call).  Load into P. */
    
    for (jy=0; jy < MY; jy++) {
      ydn = YMIN + (jy - RCONST(0.5))*dely;
      yup = ydn + dely;
      cydn = verdco*SUNRexp(RCONST(0.2)*ydn);
      cyup = verdco*SUNRexp(RCONST(0.2)*yup);
      diag = -(cydn + cyup + TWO*hordco);
      for (jx=0; jx < MX; jx++) {
        c1 = IJKth(udata,1,jx,jy);
        c2 = IJKth(udata,2,jx,jy);
        j = Jbd[jx][jy];
        a = P[jx][jy];
        IJth(j,1,1) = (-Q1*C3 - Q2*c2) + diag;
        IJth(j,1,2) = -Q2*c1 + q4coef;
        IJth(j,2,1) = Q1*C3 - Q2*c2;
        IJth(j,2,2) = (-Q2*c1 - q4coef) + diag;
        denseCopy(j, a, NUM_SPECIES, NUM_SPECIES);
      }
    }
    
    *jcurPtr = SUNTRUE;
    
  }
  
  /* Scale by -gamma */
  
  for (jy=0; jy < MY; jy++)
    for (jx=0; jx < MX; jx++)
      denseScale(-gamma, P[jx][jy], NUM_SPECIES, NUM_SPECIES);
  
  /* Add identity matrix and do LU decompositions on blocks in place. */
  
  for (jx=0; jx < MX; jx++) {
    for (jy=0; jy < MY; jy++) {
      denseAddIdentity(P[jx][jy], NUM_SPECIES);
      retval =denseGETRF(P[jx][jy], NUM_SPECIES, NUM_SPECIES, pivot[jx][jy]);
      if (retval != 0) return(1);
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
  realtype **(*P)[MY];
  sunindextype *(*pivot)[MY];
  int jx, jy;
  realtype *zdata, *v;
  UserData data;

  /* Extract the P and pivot arrays from user_data. */

  data = (UserData) user_data;
  P = data->P;
  pivot = data->pivot;
  zdata = N_VGetArrayPointer(z);
  
  N_VScale(ONE, r, z);
  
  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  
  for (jx=0; jx < MX; jx++) {
    for (jy=0; jy < MY; jy++) {
      v = &(IJKth(zdata, 1, jx, jy));
      denseGETRS(P[jx][jy], NUM_SPECIES, pivot[jx][jy], v);
    }
  }

  return(0);
}
