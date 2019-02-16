/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Cosmin Petra @ LLNL
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
 * Adjoint sensitivity example problem
 *
 * This IVP is a stiff system of 6 non-linear DAEs of index 1. The 
 * problem originates from Akzo Nobel Central research in Arnhern, 
 * The Netherlands, and describes a chemical process in which 2 
 * species are mixed, while carbon dioxide is continuously added.
 * See http://pitagora.dm.uniba.it/~testset/report/chemakzo.pdf  
 * -----------------------------------------------------------------*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include <idas/idas.h>                 /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* defs. of SUNRabs, SUNRexp, etc.      */

/* Accessor macros */
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component */

/* Problem Constants */
#define NEQ 6
#define T0  RCONST(0.0)
#define T1  RCONST(1e-8)  /* first time for output */

#define TF  RCONST(180.0) /* Final time. */
#define NF  25            /* Total number of outputs. */ 

#define RTOL  RCONST(1.0e-08)
#define ATOL  RCONST(1.0e-10)
#define RTOLQ RCONST(1.0e-10)
#define ATOLQ RCONST(1.0e-12)

#define ZERO  RCONST(0.0)
#define HALF  RCONST(0.5)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)


typedef struct {
  realtype k1, k2, k3, k4;
  realtype K, klA, Ks, pCO2, H;
} *UserData;

static int res(realtype t, N_Vector yy, N_Vector yd, N_Vector resval, void *userdata);

static int rhsQ(realtype t, N_Vector yy, N_Vector yp, 
              N_Vector qdot, void *user_data);

static void PrintHeader(realtype rtol, realtype avtol, N_Vector y);
static void PrintOutput(void *mem, realtype t, N_Vector y);
static int PrintFinalStats(void *mem);
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Main program */
int main()
{
  UserData data;
  void *mem;
  N_Vector yy, yp, rr, q;
  int retval;
  realtype time, tout, incr;
  int nout;
  SUNMatrix A;
  SUNLinearSolver LS;

  /* Consistent IC for  y, y'. */
  const realtype y01 = RCONST(0.444);
  const realtype y02 = RCONST(0.00123);
  const realtype y03 = RCONST(0.0);
  const realtype y04 = RCONST(0.007);
  const realtype y05 = RCONST(0.0);

  mem = NULL;
  yy = yp = NULL;
  A = NULL;
  LS = NULL;

  /* Allocate user data. */
  data = (UserData) malloc(sizeof(*data));

  /* Fill user's data with the appropriate values for coefficients. */
  data->k1 = RCONST(18.7);
  data->k2 = RCONST(0.58);
  data->k3 = RCONST(0.09);
  data->k4 = RCONST(0.42);
  data->K = RCONST(34.4);
  data->klA = RCONST(3.3);
  data->Ks = RCONST(115.83);
  data->pCO2 = RCONST(0.9);
  data->H = RCONST(737.0);

  /* Allocate N-vectors. */
  yy = N_VNew_Serial(NEQ);
  if (check_retval((void *)yy, "N_VNew_Serial", 0)) return(1);
  yp = N_VNew_Serial(NEQ);
  if (check_retval((void *)yp, "N_VNew_Serial", 0)) return(1);

  /* Set IC */
  Ith(yy,1) = y01;
  Ith(yy,2) = y02;
  Ith(yy,3) = y03;
  Ith(yy,4) = y04;
  Ith(yy,5) = y05;
  Ith(yy,6) = data->Ks * y01 * y04;

  /* Get y' = - res(t0, y, 0) */
  N_VConst(ZERO, yp);

  rr = N_VNew_Serial(NEQ);
  res(T0, yy, yp, rr, data);
  N_VScale(-ONE, rr, yp);
  N_VDestroy(rr);
  
 /* Create and initialize q0 for quadratures. */
  q = N_VNew_Serial(1);
  if (check_retval((void *)q, "N_VNew_Serial", 0)) return(1);
  Ith(q,1) = ZERO;

  /* Call IDACreate and IDAInit to initialize IDA memory */
  mem = IDACreate();
  if(check_retval((void *)mem, "IDACreate", 0)) return(1);

  retval = IDAInit(mem, res, T0, yy, yp);
  if(check_retval(&retval, "IDAInit", 1)) return(1);


  /* Set tolerances. */
  retval = IDASStolerances(mem, RTOL, ATOL);
  if(check_retval(&retval, "IDASStolerances", 1)) return(1);

  /* Attach user data. */
  retval = IDASetUserData(mem, data);
  if(check_retval(&retval, "IDASetUserData", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  /* Initialize QUADRATURE(S). */
  retval = IDAQuadInit(mem, rhsQ, q);
  if (check_retval(&retval, "IDAQuadInit", 1)) return(1);

  /* Set tolerances and error control for quadratures. */
  retval = IDAQuadSStolerances(mem, RTOLQ, ATOLQ);
  if (check_retval(&retval, "IDAQuadSStolerances", 1)) return(1);

  retval = IDASetQuadErrCon(mem, SUNTRUE);
  if (check_retval(&retval, "IDASetQuadErrCon", 1)) return(1);

  PrintHeader(RTOL, ATOL, yy);
  /* Print initial states */
  PrintOutput(mem,0.0,yy);

  tout = T1; nout = 0;
  incr = SUNRpowerR(TF/T1,ONE/NF);
 
  /* FORWARD run. */
  while (1) {

    retval = IDASolve(mem, tout, &time, yy, yp, IDA_NORMAL);
    if (check_retval(&retval, "IDASolve", 1)) return(1);

    PrintOutput(mem, time, yy);

    nout++;
    tout *= incr;

    if (nout>NF) break;
  }

  retval = IDAGetQuad(mem, &time, q);
  if (check_retval(&retval, "IDAGetQuad", 1)) return(1);

  printf("\n--------------------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("G:          %24.16Lf \n",Ith(q,1));
#else
  printf("G:          %24.16f \n",Ith(q,1));
#endif  
  printf("--------------------------------------------------------\n\n");

  retval = PrintFinalStats(mem);
  if (check_retval(&retval, "PrintFinalStats", 1)) return(1);
  
  IDAFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(yy);
  N_VDestroy(yp);
  N_VDestroy(q);
  free(data);

  return(0);
}


static int res(realtype t, N_Vector yy, N_Vector yd, N_Vector resval, void *userdata)
{
  UserData data;
  realtype k1, k2, k3, k4;
  realtype K, klA, Ks, pCO2, H;

  realtype y1, y2, y3, y4, y5, y6;
  realtype yd1, yd2, yd3, yd4, yd5;

  realtype r1, r2, r3, r4, r5, Fin;

  data = (UserData) userdata;
  k1 = data->k1;
  k2 = data->k2;
  k3 = data->k3;
  k4 = data->k4;
  K = data->K;
  klA = data->klA;
  Ks = data->Ks;
  pCO2 = data->pCO2;
  H = data->H;

  y1 = Ith(yy,1);
  y2 = Ith(yy,2);
  y3 = Ith(yy,3);
  y4 = Ith(yy,4);
  y5 = Ith(yy,5);
  y6 = Ith(yy,6);

  yd1 = Ith(yd,1);
  yd2 = Ith(yd,2);
  yd3 = Ith(yd,3);
  yd4 = Ith(yd,4);
  yd5 = Ith(yd,5);

  r1 = k1 * SUNRpowerI(y1,4) * SUNRsqrt(y2);
  r2 = k2 * y3 * y4;
  r3 = k2/K * y1 * y5;
  r4 = k3 * y1 * y4 * y4;
  r5 = k4 * y6 * y6 * SUNRsqrt(y2);
  Fin = klA * ( pCO2/H - y2 );

  Ith(resval,1) = yd1 + TWO*r1 - r2 + r3 + r4;
  Ith(resval,2) = yd2 + HALF*r1 + r4 + HALF*r5 - Fin;
  Ith(resval,3) = yd3 - r1 + r2 - r3;
  Ith(resval,4) = yd4 + r2 - r3 + TWO*r4;
  Ith(resval,5) = yd5 - r2 + r3 - r5;
  Ith(resval,6) = Ks*y1*y4 - y6;

  return(0);
}

/* 
 * rhsQ routine. Computes quadrature(t,y). 
 */
 
static int rhsQ(realtype t, N_Vector yy, N_Vector yp, N_Vector qdot, void *user_data)
{
  Ith(qdot,1) = Ith(yy,1);  

  return(0);
}

static void PrintHeader(realtype rtol, realtype avtol, N_Vector y)
{
  printf("\nidasAkzoNob_dns: Akzo Nobel chemical kinetics DAE serial example problem for IDAS\n");
  printf("Linear solver: DENSE, Jacobian is computed by IDAS.\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n",
         rtol, avtol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n",
         rtol, avtol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n",
         rtol, avtol);
#endif
  printf("---------------------------------------------------------------------------------\n");
  printf("   t        y1        y2       y3       y4       y5");
  printf("      y6    | nst  k      h\n");
  printf("---------------------------------------------------------------------------------\n");
}


static void PrintOutput(void *mem, realtype t, N_Vector y)
{
  realtype *yval;
  int retval, kused;
  long int nst;
  realtype hused;

  yval  = N_VGetArrayPointer(y);

  retval = IDAGetLastOrder(mem, &kused);
  check_retval(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_retval(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.2Le %8.2Le %8.2Le %8.2Le %8.2Le %8.2Le %8.2Le | %3ld  %1d %8.2Le\n", 
         t, yval[0], yval[1], yval[2], yval[3], yval[4], yval[5], nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e | %3ld  %1d %8.2e\n", 
         t, yval[0], yval[1], yval[2], yval[3], yval[4], yval[5], nst, kused, hused);
#else
  printf("%8.2e %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e | %3ld  %1d %8.2e\n", 
         t, yval[0], yval[1], yval[2], yval[3], yval[4], yval[5], nst, kused, hused);
#endif
}


static int PrintFinalStats(void *mem)
{
  int retval;
  long int nst, nni, nje, nre, nreLS, netf, ncfn;

  retval = IDAGetNumSteps(mem, &nst);
  retval = IDAGetNumResEvals(mem, &nre);
  retval = IDAGetNumJacEvals(mem, &nje);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  retval = IDAGetNumErrTestFails(mem, &netf);
  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  retval = IDAGetNumLinResEvals(mem, &nreLS);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);

  return(retval);
}


/* 
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns an integer value so check if
 *             retval < 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
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
