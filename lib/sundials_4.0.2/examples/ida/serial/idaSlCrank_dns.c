/* -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
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
 * Simulation of a slider-crank mechanism modelled with 3 generalized
 * coordinates: crank angle, connecting bar angle, and slider location.
 * The mechanism moves under the action of a constant horizontal force
 * applied to the connecting rod and a spring-damper connecting the crank
 * and connecting rod.
 *
 * The equations of motion are formulated as a system of stabilized
 * index-2 DAEs (Gear-Gupta-Leimkuhler formulation).
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>                   /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* defs. of SUNRabs, SUNRexp, etc.      */

/* Problem Constants */

#define NEQ   10

#define TEND  RCONST(10.0)
#define NOUT  41


#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)
#define FOUR RCONST(4.0)

typedef struct {
  realtype a;
  realtype J1, J2, m2;
  realtype k, c, l0;
  realtype F;
} *UserData;

int ressc(realtype tres, N_Vector yy, N_Vector yp,
           N_Vector resval, void *user_data);

void setIC(N_Vector yy, N_Vector yp, UserData data);
void force(N_Vector yy, realtype *Q, UserData data);

/* Prototypes of private functions */
static void PrintHeader(realtype rtol, realtype atol, N_Vector y);
static int PrintOutput(void *mem, realtype t, N_Vector y);
static int PrintFinalStats(void *mem);
static int check_retval(void *returnvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

int main(void)
{
  UserData data;

  void *mem;
  N_Vector yy, yp, id;
  realtype rtol, atol;
  realtype t0, tf, tout, dt, tret;
  int retval, iout;
  SUNMatrix A;
  SUNLinearSolver LS;

  A = NULL;
  LS = NULL;

  /* User data */

  data = (UserData) malloc(sizeof *data);

  data->a = 0.5;   /* half-length of crank */
  data->J1 = 1.0;  /* crank moment of inertia */
  data->m2 = 1.0;  /* mass of connecting rod */
  data->J2 = 2.0;  /* moment of inertia of connecting rod */
  data->k = 1.0;   /* spring constant */
  data->c = 1.0;   /* damper constant */
  data->l0 = 1.0;  /* spring free length */
  data->F = 1.0;   /* external constant force */

  /* Create N_Vectors */
  yy = N_VNew_Serial(NEQ);
  yp = N_VNew_Serial(NEQ);
  id = N_VNew_Serial(NEQ);

  /* Consistent IC */
  setIC(yy, yp, data);

  /* ID array */
  N_VConst(ONE, id);
  NV_Ith_S(id,6) = ZERO;
  NV_Ith_S(id,7) = ZERO;
  NV_Ith_S(id,8) = ZERO;
  NV_Ith_S(id,9) = ZERO;

  /* Tolerances */
  rtol = RCONST(1.0e-6);
  atol = RCONST(1.0e-6);

  /* Integration limits */
  t0 = ZERO;
  tf = TEND;
  dt = (tf-t0)/(NOUT-1);

  /* IDA initialization */
  mem = IDACreate();
  retval = IDAInit(mem, ressc, t0, yy, yp);
  retval = IDASStolerances(mem, rtol, atol);
  retval = IDASetUserData(mem, data);
  retval = IDASetId(mem, id);
  retval = IDASetSuppressAlg(mem, SUNTRUE);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  PrintHeader(rtol, atol, yy);

  /* In loop, call IDASolve, print results, and test for error. */

  retval = PrintOutput(mem,t0,yy);

  tout = dt;
  for (iout=1; iout<NOUT; iout++) {
    tout = iout*dt;
    retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
    if (retval < 0) break;

    retval = PrintOutput(mem,tret,yy);

  }

  retval = PrintFinalStats(mem);

  /* Free memory */

  free(data);
  IDAFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(yy);
  N_VDestroy(yp);
  N_VDestroy(id);

  return(0);

}

void setIC(N_Vector yy, N_Vector yp, UserData data)
{
  realtype pi;
  realtype a, J1, m2, J2;
  realtype q, p, x;
  realtype Q[3];

  N_VConst(ZERO, yy);
  N_VConst(ZERO, yp);

  pi = FOUR*atan(ONE);

  a = data->a;
  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  q = pi/TWO;
  p = asin(-a);
  x = cos(p);

  NV_Ith_S(yy,0) = q;
  NV_Ith_S(yy,1) = x;
  NV_Ith_S(yy,2) = p;

  force(yy, Q, data);

  NV_Ith_S(yp,3) = Q[0]/J1;
  NV_Ith_S(yp,4) = Q[1]/m2;
  NV_Ith_S(yp,5) = Q[2]/J2;

}

void force(N_Vector yy, realtype *Q, UserData data)
{
  realtype a, k, c, l0, F;
  realtype q, x, p;
  realtype qd, xd, pd;
  realtype s1, c1, s2, c2, s21, c21;
  realtype l2, l, ld;
  realtype f, fl;

  a = data->a;
  k = data->k;
  c = data->c;
  l0 = data->l0;
  F = data->F;

  q = NV_Ith_S(yy,0);
  x = NV_Ith_S(yy,1);
  p = NV_Ith_S(yy,2);

  qd = NV_Ith_S(yy,3);
  xd = NV_Ith_S(yy,4);
  pd = NV_Ith_S(yy,5);

  s1 = sin(q);
  c1 = cos(q);
  s2 = sin(p);
  c2 = cos(p);
  s21 = s2*c1 - c2*s1;
  c21 = c2*c1 + s2*s1;

  l2 = x*x - x*(c2+a*c1) + (ONE + a*a)/FOUR + a*c21/TWO;
  l = SUNRsqrt(l2);
  ld = TWO*x*xd - xd*(c2+a*c1) + x*(s2*pd+a*s1*qd) - a*s21*(pd-qd)/TWO;
  ld /= TWO*l;

  f = k*(l-l0) + c*ld;
  fl = f/l;

  Q[0] = - fl * a * (s21/TWO + x*s1) / TWO;
  Q[1] = fl * (c2/TWO - x + a*c1/TWO) + F;
  Q[2] = - fl * (x*s2 - a*s21/TWO) / TWO - F*s2;

}

int ressc(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
  UserData data;
  realtype Q[3];
  realtype a, J1, m2, J2;
  realtype *yval, *ypval, *rval;
  realtype q, x, p;
  realtype qd, xd, pd;
  realtype lam1, lam2, mu1, mu2;
  realtype s1, c1, s2, c2;

  data = (UserData) user_data;

  a  = data->a;
  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  yval = N_VGetArrayPointer(yy);
  ypval = N_VGetArrayPointer(yp);
  rval = N_VGetArrayPointer(rr);

  q = yval[0];
  x = yval[1];
  p = yval[2];

  qd = yval[3];
  xd = yval[4];
  pd = yval[5];

  lam1 = yval[6];
  lam2 = yval[7];

  mu1 = yval[8];
  mu2 = yval[9];

  s1 = sin(q);
  c1 = cos(q);
  s2 = sin(p);
  c2 = cos(p);

  force(yy, Q, data);

  rval[0] = ypval[0] - qd + a*s1*mu1 - a*c1*mu2;
  rval[1] = ypval[1] - xd + mu1;
  rval[2] = ypval[2] - pd + s2*mu1 - c2*mu2;

  rval[3] = J1*ypval[3] - Q[0] + a*s1*lam1 - a*c1*lam2;
  rval[4] = m2*ypval[4] - Q[1] + lam1;
  rval[5] = J2*ypval[5] - Q[2] + s2*lam1 - c2*lam2;

  rval[6] = x - c2 - a*c1;
  rval[7] = -s2 - a*s1;

  rval[8] = a*s1*qd + xd + s2*pd;
  rval[9] = -a*c1*qd - c2*pd;

  return(0);
}

static void PrintHeader(realtype rtol, realtype atol, N_Vector y)
{
  printf("\nidaSlCrank_dns: Slider-Crank DAE serial example problem for IDA\n");
  printf("Linear solver: DENSE, Jacobian is computed by IDA.\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n",
         rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n",
         rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n",
         rtol, atol);
#endif
  printf("-----------------------------------------------------------------------\n");
  printf("  t            y1          y2           y3");
  printf("      | nst  k      h\n");
  printf("-----------------------------------------------------------------------\n");
}

static int PrintOutput(void *mem, realtype t, N_Vector y)
{
  realtype *yval;
  int retval, kused;
  long int nst;
  realtype hused;

  yval  = N_VGetArrayPointer(y);

  retval = IDAGetLastOrder(mem, &kused);
  retval = IDAGetNumSteps(mem, &nst);
  retval = IDAGetLastStep(mem, &hused);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%10.4Le %12.4Le %12.4Le %12.4Le %3ld  %1d %12.4Le\n",
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#else
  printf("%10.4e %12.4e %12.4e %12.4e %3ld  %1d %12.4e\n",
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#endif

  return(retval);
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
    return(1);
  } else if (opt == 1) {
    /* Check if retval < 0 */
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return(1);
    }
  } else if (opt == 2 && returnvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  }

  return(0);
}
