/* -----------------------------------------------------------------
 * Programmer: Radu Serban and Cosmin Petra @ LLNL
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
 *
 * IDAS also computes the average kinetic energy as the quadrature:
 *   G = int_t0^tend g(t,y,p) dt, 
 * where
 *   g(t,y,p) = 0.5*J1*v1^2 + 0.5*J2*v3^2 + 0.5*m2*v2^2
 *              
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <idas/idas.h>                 /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* defs. of SUNRabs, SUNRexp, etc.      */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i= 1..NEQ */

/* Problem Constants */

#define NEQ   10

#define TBEGIN  RCONST(0.0)
#define TEND    RCONST(10.0)

#define NOUT  25

#define RTOLF   RCONST(1.0e-06)
#define ATOLF   RCONST(1.0e-07)

#define RTOLQ   RCONST(1.0e-06)
#define ATOLQ   RCONST(1.0e-08)

#define ZERO     RCONST(0.00)
#define QUARTER  RCONST(0.25)
#define HALF     RCONST(0.50)
#define ONE      RCONST(1.00)
#define TWO      RCONST(2.00)
#define FOUR     RCONST(4.00)

typedef struct {
  realtype a;
  realtype J1, J2, m1, m2;
  realtype l0;
  realtype params[2];
  realtype F;
} *UserData;

static int ressc(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data);
static int rhsQ(realtype t, N_Vector yy, N_Vector yp, N_Vector qdot, void *user_data);

static void setIC(N_Vector yy, N_Vector yp, UserData data);
static void force(N_Vector yy, realtype *Q, UserData data);

/* Prototypes of private functions */
static void PrintHeader(realtype rtol, realtype avtol, N_Vector y);
static void PrintOutput(void *mem, realtype t, N_Vector y);
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
  N_Vector yy, yp, id, q;
  realtype tret, tout;
  int retval;
  SUNMatrix A;
  SUNLinearSolver LS;

  A = NULL;
  LS = NULL;

  id = N_VNew_Serial(NEQ);
  yy = N_VNew_Serial(NEQ);
  yp = N_VNew_Serial(NEQ);
  q = N_VNew_Serial(1);

  data = (UserData) malloc(sizeof *data);

  data->a = 0.5;   /* half-length of crank */
  data->J1 = 1.0;  /* crank moment of inertia */
  data->m2 = 1.0;  /* mass of connecting rod */
  data->m1 = 1.0;
  data->J2 = 2.0;  /* moment of inertia of connecting rod */
  data->params[0] = 1.0;   /* spring constant */
  data->params[1] = 1.0;   /* damper constant */
  data->l0 = 1.0;  /* spring free length */
  data->F = 1.0;   /* external constant force */

  N_VConst(ONE, id);
  NV_Ith_S(id, 9) = ZERO;
  NV_Ith_S(id, 8) = ZERO;
  NV_Ith_S(id, 7) = ZERO;
  NV_Ith_S(id, 6) = ZERO;
  
  /* Consistent IC*/
  setIC(yy, yp, data);

  /* IDAS initialization */
  mem = IDACreate();
  retval = IDAInit(mem, ressc, TBEGIN, yy, yp);
  retval = IDASStolerances(mem, RTOLF, ATOLF);
  retval = IDASetUserData(mem, data);
  retval = IDASetId(mem, id);
  retval = IDASetSuppressAlg(mem, SUNTRUE);
  retval = IDASetMaxNumSteps(mem, 20000);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  N_VConst(ZERO, q);
  retval = IDAQuadInit(mem, rhsQ, q);
  retval = IDAQuadSStolerances(mem, RTOLQ, ATOLQ);
  retval = IDASetQuadErrCon(mem, SUNTRUE);

  PrintHeader(RTOLF, ATOLF, yy);

  /* Print initial states */
  PrintOutput(mem,0.0,yy);

  /* Perform forward run */
  tout = TEND/NOUT;

  while (1) {

    retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
    if (check_retval(&retval, "IDASolve", 1)) return(1);

    PrintOutput(mem,tret,yy);

    tout += TEND/NOUT;
    
    if (tret > TEND) break;
  }

  retval = PrintFinalStats(mem);
  if (check_retval(&retval, "PrintFinalStats", 1)) return(1);

  IDAGetQuad(mem, &tret, q);
  printf("--------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("  G = %24.16Lf\n", Ith(q,1));
#else
  printf("  G = %24.16f\n", Ith(q,1));
#endif  
  printf("--------------------------------------------\n\n");
  
  IDAFree(&mem);

  /* Free memory */

  free(data);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(id);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(q);

  return(0);  
}

static void setIC(N_Vector yy, N_Vector yp, UserData data)
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

static void force(N_Vector yy, realtype *Q, UserData data)
{
  realtype a, k, c, l0, F;
  realtype q, x, p;
  realtype qd, xd, pd;  
  realtype s1, c1, s2, c2, s21, c21;
  realtype l2, l, ld;
  realtype f, fl;

  a = data->a;
  k = data->params[0];
  c = data->params[1];
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

static int ressc(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
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

static int rhsQ(realtype t, N_Vector yy, N_Vector yp, N_Vector qdot, void *user_data)
{
  realtype v1, v2, v3;
  realtype J1, m2, J2;
  UserData data;
  
  data = (UserData) user_data;
  J1 = data->J1;
  m2 = data->m2;
  J2 = data->J2;

  v1 = Ith(yy,4); 
  v2 = Ith(yy,5); 
  v3 = Ith(yy,6);

  Ith(qdot,1) = HALF*(J1*v1*v1 + m2*v2*v2 + J2*v3*v3);

  return(0);
}

static void PrintHeader(realtype rtol, realtype avtol, N_Vector y)
{
  printf("\nidasSlCrank_dns: Slider-Crank DAE serial example problem for IDAS\n");
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
  printf("-----------------------------------------------------------------------\n");
  printf("  t         y1          y2           y3");
  printf("      | nst  k      h\n");
  printf("-----------------------------------------------------------------------\n");
}

static void PrintOutput(void *mem, realtype t, N_Vector y)
{
  realtype *yval;
  int retval, kused;
  long int nst;
  realtype hused;

  yval  = N_VGetArrayPointer_Serial(y);

  retval = IDAGetLastOrder(mem, &kused);
  check_retval(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_retval(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%5.2Lf %12.4Le %12.4Le %12.4Le | %3ld  %1d %12.4Le\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%5.2f %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#else
  printf("%5.2f %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
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
