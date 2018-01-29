/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 *
 * This example loops through the available iterative linear solvers:
 * SPGMR, SPBCG and SPTFQMR.
 *
 * Example problem for IDA: 2D heat equation, serial, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version loops through the Krylov solvers Spgmr, Spbcg
 * and Sptfqmr.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y). The
 * PDE is treated with central differences on a uniform M x M grid.
 * The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = M^2. Here M = 10.
 *
 * The system is solved with IDA using the following Krylov
 * linear solvers: SPGMR, SPBCG and SPTFQMR. The
 * preconditioner uses the diagonal elements of the Jacobian only.
 * Routines for preconditioning, required by SP*, are supplied
 * here. The constraints u >= 0 are posed for all components. Output
 * is taken at t = 0, .01, .02, .04,..., 10.24.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <idas/idas.h>                   /* main integrator header file       */
#include <idas/idas_spils.h>             /* access to IDASpils interface      */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver   */
#include <sunlinsol/sunlinsol_spbcgs.h>  /* access to SPBCGS SUNLinearSolver  */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to SPTFQMR SUNLinearSolver */
#include <nvector/nvector_serial.h>      /* serial N_Vector types, fct. and macros */
#include <sundials/sundials_types.h>     /* definition of realtype */

/* Problem Constants */

#define NOUT  11
#define MGRID 10
#define NEQ   MGRID*MGRID
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define FOUR  RCONST(4.0)

/* Linear Solver Loop Constants */

#define USE_SPGMR   0
#define USE_SPBCG   1
#define USE_SPTFQMR 2

/* User data type */

typedef struct {  
  sunindextype mm;  /* number of grid points */
  realtype dx;
  realtype coeff;
  N_Vector pp;  /* vector of prec. diag. elements */
} *UserData;

/* Prototypes for functions called by IDA */

int resHeat(realtype tres, N_Vector uu, N_Vector up,
            N_Vector resval, void *user_data);

int PsetupHeat(realtype tt, 
               N_Vector uu, N_Vector up, N_Vector rr, 
               realtype c_j, void *user_data);

int PsolveHeat(realtype tt, 
               N_Vector uu, N_Vector up, N_Vector rr, 
               N_Vector rvec, N_Vector zvec, 
               realtype c_j, realtype delta, void *user_data);

/* Prototypes for private functions */

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector res);
static void PrintHeader(realtype rtol, realtype atol, int linsolver);
static void PrintOutput(void *mem, realtype t, N_Vector uu, int linsolver);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(void)
{
  void *mem;
  UserData data;
  N_Vector uu, up, constraints, res;
  int ier, iout, linsolver;
  realtype rtol, atol, t0, t1, tout, tret;
  long int netf, ncfn, ncfl;
  SUNLinearSolver LS;

  mem = NULL;
  data = NULL;
  uu = up = constraints = res = NULL;
  LS = NULL;

  /* Allocate N-vectors and the user data structure. */

  uu = N_VNew_Serial(NEQ);
  if(check_flag((void *)uu, "N_VNew_Serial", 0)) return(1);

  up = N_VNew_Serial(NEQ);
  if(check_flag((void *)up, "N_VNew_Serial", 0)) return(1);

  res = N_VNew_Serial(NEQ);
  if(check_flag((void *)res, "N_VNew_Serial", 0)) return(1);

  constraints = N_VNew_Serial(NEQ);
  if(check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);

  data = (UserData) malloc(sizeof *data);
  data->pp = NULL;
  if(check_flag((void *)data, "malloc", 2)) return(1);

  /* Assign parameters in the user data structure. */

  data->mm  = MGRID;
  data->dx = ONE/(MGRID-ONE);
  data->coeff = ONE/(data->dx * data->dx);
  data->pp = N_VNew_Serial(NEQ);
  if(check_flag((void *)data->pp, "N_VNew_Serial", 0)) return(1);

  /* Initialize uu, up. */

  SetInitialProfile(data, uu, up, res);

  /* Set constraints to all 1's for nonnegative solution values. */

  N_VConst(ONE, constraints);

  /* Assign various parameters. */

  t0   = ZERO;
  t1   = RCONST(0.01);
  rtol = ZERO;
  atol = RCONST(1.0e-3); 

  /* Call IDACreate and IDAMalloc to initialize solution */

  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);

  ier = IDASetUserData(mem, data);
  if(check_flag(&ier, "IDASetUserData", 1)) return(1);

  ier = IDASetConstraints(mem, constraints);
  if(check_flag(&ier, "IDASetConstraints", 1)) return(1);
  N_VDestroy(constraints);

  ier = IDAInit(mem, resHeat, t0, uu, up);
  if(check_flag(&ier, "IDAInit", 1)) return(1);

  ier = IDASStolerances(mem, rtol, atol);
  if(check_flag(&ier, "IDASStolerances", 1)) return(1);

  /* START: Loop through SPGMR, SPBCG and SPTFQMR linear solver modules */
  for (linsolver = 0; linsolver < 3; ++linsolver) {

    if (linsolver != 0) {

      /* Re-initialize uu, up. */
      SetInitialProfile(data, uu, up, res);

      /* Re-initialize IDA */
      ier = IDAReInit(mem, t0, uu, up);
      if (check_flag(&ier, "IDAReInit", 1)) return(1);

    }

    /* Attach a linear solver module */
    switch(linsolver) {

    /* (a) SPGMR */
    case(USE_SPGMR):

      /* Print header */
      printf(" -------");
      printf(" \n| SPGMR |\n");
      printf(" -------\n");

      /* Call SUNSPGMR to specify the linear solver SPGMR with
         left preconditioning and the default maximum Krylov dimension */
      LS = SUNSPGMR(uu, PREC_LEFT, 0);
      if(check_flag((void *)LS, "SUNSPGMR", 0)) return(1);

      /* Attach the linear solver */
      ier = IDASpilsSetLinearSolver(mem, LS);
      if(check_flag(&ier, "IDASpilsSetLinearSolver", 1)) return 1;

      break;

    /* (b) SPBCG */
    case(USE_SPBCG):

      /* Print header */
      printf(" -------");
      printf(" \n| SPBCGS |\n");
      printf(" -------\n");

      /* Call SUNSPBCGS to specify the linear solver SPBCGS with
         left preconditioning and the default maximum Krylov dimension */
      LS = SUNSPBCGS(uu, PREC_LEFT, 0);
      if(check_flag((void *)LS, "SUNSPBCGS", 0)) return(1);

      /* Attach the linear solver */
      ier = IDASpilsSetLinearSolver(mem, LS);
      if(check_flag(&ier, "IDASpilsSetLinearSolver", 1)) return 1;

      break;

    /* (c) SPTFQMR */
    case(USE_SPTFQMR):

      /* Print header */
      printf(" ---------");
      printf(" \n| SPTFQMR |\n");
      printf(" ---------\n");

      /* Call SUNSPTFQMR to specify the linear solver SPTFQMR with
         left preconditioning and the default maximum Krylov dimension */
      LS = SUNSPTFQMR(uu, PREC_LEFT, 0);
      if(check_flag((void *)LS, "SUNSPTFQMR", 0)) return(1);

      /* Attach the linear solver */
      ier = IDASpilsSetLinearSolver(mem, LS);
      if(check_flag(&ier, "IDASpilsSetLinearSolver", 1)) return 1;

      break;

    }

    /* Specify preconditioner */
    ier = IDASpilsSetPreconditioner(mem, PsetupHeat, PsolveHeat);
    if(check_flag(&ier, "IDASpilsSetPreconditioner", 1)) return(1);

    /* Print output heading. */
    PrintHeader(rtol, atol, linsolver);

    /* Print output table heading, and initial line of table. */

    printf("\n   Output Summary (umax = max-norm of solution) \n\n");
    printf("  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps\n" );
    printf("----------------------------------------------------------------------\n");

    /* Loop over output times, call IDASolve, and print results. */

    for (tout = t1,iout = 1; iout <= NOUT ; iout++, tout *= TWO) {
      ier = IDASolve(mem, tout, &tret, uu, up, IDA_NORMAL);
      if(check_flag(&ier, "IDASolve", 1)) return(1);
      PrintOutput(mem, tret, uu, linsolver);
    }

    /* Print remaining counters. */
    ier = IDAGetNumErrTestFails(mem, &netf);
    check_flag(&ier, "IDAGetNumErrTestFails", 1);

    ier = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
    check_flag(&ier, "IDAGetNumNonlinSolvConvFails", 1);

    ier = IDASpilsGetNumConvFails(mem, &ncfl);
    check_flag(&ier, "IDASpilsGetNumConvFails", 1);

    printf("\nError test failures            = %ld\n", netf);
    printf("Nonlinear convergence failures = %ld\n", ncfn);
    printf("Linear convergence failures    = %ld\n", ncfl);

    if (linsolver < 2)
      printf("\n======================================================================\n\n");

  } /* END: Loop through SPGMR, SPBCG and SPTFQMR linear solver modules */

  /* Free Memory */

  IDAFree(&mem);
  SUNLinSolFree(LS);

  N_VDestroy(uu);
  N_VDestroy(up);
  N_VDestroy(res);

  N_VDestroy(data->pp);
  free(data);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * resHeat: heat equation system residual function (user-supplied)      
 * This uses 5-point central differencing on the interior points, and   
 * includes algebraic equations for the boundary values.                
 * So for each interior point, the residual component has the form      
 *    res_i = u'_i - (central difference)_i                             
 * while for each boundary point, it is res_i = u_i.                     
 */

int resHeat(realtype tt, 
            N_Vector uu, N_Vector up, N_Vector rr, 
            void *user_data)
{
  sunindextype i, j, offset, loc, mm;
  realtype *uu_data, *up_data, *rr_data, coeff, dif1, dif2;
  UserData data;
  
  uu_data = N_VGetArrayPointer(uu); 
  up_data = N_VGetArrayPointer(up); 
  rr_data = N_VGetArrayPointer(rr);

  data = (UserData) user_data;
  
  coeff = data->coeff;
  mm    = data->mm;
  
  /* Initialize rr to uu, to take care of boundary equations. */
  N_VScale(ONE, uu, rr);
  
  /* Loop over interior points; set res = up - (central difference). */
  for (j = 1; j < MGRID-1; j++) {
    offset = mm*j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      dif1 = uu_data[loc-1]  + uu_data[loc+1]  - TWO * uu_data[loc];
      dif2 = uu_data[loc-mm] + uu_data[loc+mm] - TWO * uu_data[loc];
      rr_data[loc]= up_data[loc] - coeff * ( dif1 + dif2 );
    }
  }

  return(0);
}

/*
 * PsetupHeat: setup for diagonal preconditioner.   
 *                                                                 
 * The optional user-supplied functions PsetupHeat and          
 * PsolveHeat together must define the left preconditoner        
 * matrix P approximating the system Jacobian matrix               
 *                   J = dF/du + cj*dF/du'                         
 * (where the DAE system is F(t,u,u') = 0), and solve the linear   
 * systems P z = r.   This is done in this case by keeping only    
 * the diagonal elements of the J matrix above, storing them as    
 * inverses in a vector pp, when computed in PsetupHeat, for    
 * subsequent use in PsolveHeat.                                 
 *                                                                 
 * In this instance, only cj and data (user data structure, with    
 * pp etc.) are used from the PsetupdHeat argument list.         
 */
  
int PsetupHeat(realtype tt, 
               N_Vector uu, N_Vector up, N_Vector rr, 
               realtype c_j, void *user_data)
{
  
  sunindextype i, j, offset, loc, mm;
  realtype *ppv, pelinv;
  UserData data;
  
  data = (UserData) user_data;
  ppv = N_VGetArrayPointer(data->pp);
  mm = data->mm;

  /* Initialize the entire vector to 1., then set the interior points to the
     correct value for preconditioning. */
  N_VConst(ONE,data->pp);
  
  /* Compute the inverse of the preconditioner diagonal elements. */
  pelinv = ONE/(c_j + FOUR*data->coeff); 
  
  for (j = 1; j < mm-1; j++) {
    offset = mm * j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      ppv[loc] = pelinv;
    }
  }
  
  return(0);  
}

/*
 * PsolveHeat: solve preconditioner linear system.              
 * This routine multiplies the input vector rvec by the vector pp 
 * containing the inverse diagonal Jacobian elements (previously  
 * computed in PrecondHeateq), returning the result in zvec.      
 */

int PsolveHeat(realtype tt, 
               N_Vector uu, N_Vector up, N_Vector rr, 
               N_Vector rvec, N_Vector zvec, 
               realtype c_j, realtype delta, void *user_data)
{
  UserData data;
  data = (UserData) user_data;
  N_VProd(data->pp, rvec, zvec);
  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * SetInitialProfile: routine to initialize u and up vectors.
 */

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up, 
                             N_Vector res)
{
  sunindextype mm, mm1, i, j, offset, loc;
  realtype xfact, yfact, *udata, *updata;

  mm = data->mm;

  udata = N_VGetArrayPointer(uu);
  updata = N_VGetArrayPointer(up);

  /* Initialize uu on all grid points. */ 
  mm1 = mm - 1;
  for (j = 0; j < mm; j++) {
    yfact = data->dx * j;
    offset = mm*j;
    for (i = 0;i < mm; i++) {
      xfact = data->dx * i;
      loc = offset + i;
      udata[loc] = RCONST(16.0) * xfact * (ONE - xfact) * yfact * (ONE - yfact);
    }
  }
  
  /* Initialize up vector to 0. */
  N_VConst(ZERO, up);

  /* resHeat sets res to negative of ODE RHS values at interior points. */
  resHeat(ZERO, uu, up, res, data);

  /* Copy -res into up to get correct interior initial up values. */
  N_VScale(-ONE, res, up);

  /* Set up at boundary points to zero. */
  for (j = 0; j < mm; j++) {
    offset = mm*j;
    for (i = 0; i < mm; i++) {
      loc = offset + i;
      if (j == 0 || j == mm1 || i == 0 || i == mm1 ) updata[loc] = ZERO;
    }
  }
  
  return(0);
  }

/*
 * Print first lines of output (problem description)
 */

static void PrintHeader(realtype rtol, realtype atol, int linsolver)
{
  printf("\nidasKrylovDemo_ls: Heat equation, serial example problem for IDA\n");
  printf("                   Discretized heat equation on 2D unit square.\n");
  printf("                   Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("                   Mesh dimensions: %d x %d", MGRID, MGRID);
  printf("       Total system size: %d\n\n", NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION) 
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Constraints set to force all solution components >= 0. \n");

  switch(linsolver) {

  case(USE_SPGMR):
    printf("Linear solver: SPGMR, preconditioner using diagonal elements. \n");
    break;

  case(USE_SPBCG):
    printf("Linear solver: SPBCG, preconditioner using diagonal elements. \n");
    break;

  case(USE_SPTFQMR):
    printf("Linear solver: SPTFQMR, preconditioner using diagonal elements. \n");
    break;
  }
}

/*
 * PrintOutput: print max norm of solution and current solver statistics
 */

static void PrintOutput(void *mem, realtype t, N_Vector uu, int linsolver)
{
  realtype hused, umax;
  long int nst, nni, nje, nre, nreLS, nli, npe, nps;
  int kused, ier;
  
  umax = N_VMaxNorm(uu);

  ier = IDAGetLastOrder(mem, &kused);
  check_flag(&ier, "IDAGetLastOrder", 1);
  ier = IDAGetNumSteps(mem, &nst);
  check_flag(&ier, "IDAGetNumSteps", 1);
  ier = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&ier, "IDAGetNumNonlinSolvIters", 1);
  ier = IDAGetNumResEvals(mem, &nre);
  check_flag(&ier, "IDAGetNumResEvals", 1);
  ier = IDAGetLastStep(mem, &hused);
  check_flag(&ier, "IDAGetLastStep", 1);

  ier = IDASpilsGetNumJtimesEvals(mem, &nje);
  check_flag(&ier, "IDASpilsGetNumJtimesEvals", 1);
  ier = IDASpilsGetNumLinIters(mem, &nli);
  check_flag(&ier, "IDASpilsGetNumLinIters", 1);
  ier = IDASpilsGetNumResEvals(mem, &nreLS);
  check_flag(&ier, "IDASpilsGetNumResEvals", 1);
  ier = IDASpilsGetNumPrecEvals(mem, &npe);
  check_flag(&ier, "IDASpilsGetPrecEvals", 1);
  ier = IDASpilsGetNumPrecSolves(mem, &nps);
  check_flag(&ier, "IDASpilsGetNumPrecSolves", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2Le  %3ld %3ld\n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#elif defined(SUNDIALS_DOUBLE_PRECISION) 
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#else
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
         t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#endif
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
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  }

  return(0);
}
