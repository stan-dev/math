/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Adjoint sensitivity example problem.
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES. The problem is from chemical
 * kinetics, and consists of the following three rate equations.
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The reaction rates are:
 * p1=0.04, p2=1e4, and p3=3e7. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVODE dense linear solver, and a user-supplied
 * Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 * Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * 
 * Optionally, CVODES can compute sensitivities with respect to
 * the problem parameters p1, p2, and p3 of the following quantity:
 *   G = int_t0^t1 g(t,p,y) dt
 * where
 *   g(t,p,y) = y3
 *        
 * The gradient dG/dp is obtained as:
 *   dG/dp = int_t0^t1 (g_p - lambda^T f_p ) dt - lambda^T(t0)*y0_p
 *         = - xi^T(t0) - lambda^T(t0)*y0_p
 * where lambda and xi are solutions of:
 *   d(lambda)/dt = - (f_y)^T * lambda - (g_y)^T
 *   lambda(t1) = 0
 * and
 *   d(xi)/dt = - (f_p)^T * lambda + (g_p)^T
 *   xi(t1) = 0
 * 
 * During the backward integration, CVODES also evaluates G as
 *   G = - phi(t0)
 * where
 *   d(phi)/dt = g(t,y,p)
 *   phi(t1) = 0
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvodes/cvodes_direct.h>      /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* defs. of SUNRabs, SUNRexp, etc.      */

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component, i=1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* (i,j)-th matrix el., i,j=1..NEQ */

/* Problem Constants */

#define NEQ      3             /* number of equations                  */

#define RTOL     RCONST(1e-6)  /* scalar relative tolerance            */

#define ATOL1    RCONST(1e-8)  /* vector absolute tolerance components */
#define ATOL2    RCONST(1e-14)
#define ATOL3    RCONST(1e-6)

#define ATOLl    RCONST(1e-8)  /* absolute tolerance for adjoint vars. */
#define ATOLq    RCONST(1e-6)  /* absolute tolerance for quadratures   */

#define T0       RCONST(0.0)   /* initial time                         */
#define TOUT     RCONST(4e7)   /* final time                           */

#define TB1      RCONST(4e7)   /* starting point for adjoint problem   */
#define TB2      RCONST(50.0)  /* starting point for adjoint problem   */
#define TBout1   RCONST(40.0)  /* intermediate t for adjoint problem   */

#define STEPS    150           /* number of steps between check points */

#define NP       3             /* number of problem parameters         */

#define ZERO     RCONST(0.0)


/* Type : UserData */

typedef struct {
  realtype p[3];
} *UserData;

/* Prototypes of user-supplied functions */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data);
static int ewt(N_Vector y, N_Vector w, void *user_data);

static int fB(realtype t, N_Vector y, 
              N_Vector yB, N_Vector yBdot, void *user_dataB);
static int JacB(realtype t, N_Vector y, N_Vector yB, N_Vector fyB, SUNMatrix JB,
                void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
static int fQB(realtype t, N_Vector y, N_Vector yB, 
               N_Vector qBdot, void *user_dataB);


/* Prototypes of private functions */

static void PrintHead(realtype tB0);
static void PrintOutput(realtype tfinal, N_Vector y, N_Vector yB, N_Vector qB);
static void PrintOutput1(realtype time, realtype t, N_Vector y, N_Vector yB);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;

  SUNMatrix A, AB;
  SUNLinearSolver LS, LSB;
  void *cvode_mem;

  realtype reltolQ, abstolQ;
  N_Vector y, q;

  int steps;

  int indexB;

  realtype reltolB, abstolB, abstolQB;
  N_Vector yB, qB;

  realtype time;
  int flag, ncheck;

  long int nst, nstB;

  CVadjCheckPointRec *ckpnt;

  data = NULL;
  A = AB = NULL;
  LS = LSB = NULL;
  cvode_mem = NULL;
  ckpnt = NULL;
  y = yB = qB = NULL;

  /* Print problem description */
  printf("\nAdjoint Sensitivity Example for Chemical Kinetics\n");
  printf("-------------------------------------------------\n\n");
  printf("ODE: dy1/dt = -p1*y1 + p2*y2*y3\n");
  printf("     dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2\n");
  printf("     dy3/dt =  p3*(y2)^2\n\n");
  printf("Find dG/dp for\n");
  printf("     G = int_t0^tB0 g(t,p,y) dt\n");
  printf("     g(t,p,y) = y3\n\n\n");

  /* User data structure */
  data = (UserData) malloc(sizeof *data);
  if (check_flag((void *)data, "malloc", 2)) return(1);
  data->p[0] = RCONST(0.04);
  data->p[1] = RCONST(1.0e4);
  data->p[2] = RCONST(3.0e7);

  /* Initialize y */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  Ith(y,1) = RCONST(1.0);
  Ith(y,2) = ZERO;
  Ith(y,3) = ZERO;

  /* Initialize q */
  q = N_VNew_Serial(1);
  if (check_flag((void *)q, "N_VNew_Serial", 0)) return(1);
  Ith(q,1) = ZERO;

  /* Set the scalar realtive and absolute tolerances reltolQ and abstolQ */
  reltolQ = RTOL;
  abstolQ = ATOLq;

  /* Create and allocate CVODES memory for forward run */
  printf("Create and allocate CVODES memory for forward runs\n");

  /* Call CVodeCreate to create the solver memory and specify the 
     Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the initial time T0, and
     the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Call CVodeWFtolerances to specify a user-supplied function ewt that sets
     the multiplicative error weights w_i for use in the weighted RMS norm */
  flag = CVodeWFtolerances(cvode_mem, ewt);
  if (check_flag(&flag, "CVodeWFtolerances", 1)) return(1);

  /* Attach user data */
  flag = CVodeSetUserData(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ);
  if (check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNDenseLinearSolver(y, A);
  if (check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

  /* Attach the matrix and linear solver */
  flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
  if (check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

  /* Set the user-supplied Jacobian routine Jac */
  flag = CVDlsSetJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetJacFn", 1)) return(1);

  /* Call CVodeQuadInit to allocate initernal memory and initialize
     quadrature integration*/
  flag = CVodeQuadInit(cvode_mem, fQ, q);
  if (check_flag(&flag, "CVodeQuadInit", 1)) return(1);

  /* Call CVodeSetQuadErrCon to specify whether or not the quadrature variables
     are to be used in the step size control mechanism within CVODES. Call
     CVodeQuadSStolerances or CVodeQuadSVtolerances to specify the integration
     tolerances for the quadrature variables. */
  flag = CVodeSetQuadErrCon(cvode_mem, SUNTRUE);
  if (check_flag(&flag, "CVodeSetQuadErrCon", 1)) return(1);

  /* Call CVodeQuadSStolerances to specify scalar relative and absolute
     tolerances. */
  flag = CVodeQuadSStolerances(cvode_mem, reltolQ, abstolQ);
  if (check_flag(&flag, "CVodeQuadSStolerances", 1)) return(1);

  /* Allocate global memory */

  /* Call CVodeAdjInit to update CVODES memory block by allocting the internal 
     memory needed for backward integration.*/
  steps = STEPS; /* no. of integration steps between two consecutive ckeckpoints*/
  flag = CVodeAdjInit(cvode_mem, steps, CV_HERMITE);
  /*
  flag = CVodeAdjInit(cvode_mem, steps, CV_POLYNOMIAL);
  */
  if (check_flag(&flag, "CVodeAdjInit", 1)) return(1);

  /* Perform forward run */
  printf("Forward integration ... ");

  /* Call CVodeF to integrate the forward problem over an interval in time and
     saves checkpointing data */
  flag = CVodeF(cvode_mem, TOUT, y, &time, CV_NORMAL, &ncheck);
  if (check_flag(&flag, "CVodeF", 1)) return(1);
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  if (check_flag(&flag, "CVodeGetNumSteps", 1)) return(1);

  printf("done ( nst = %ld )\n",nst);
  printf("\nncheck = %d\n\n", ncheck);

  flag = CVodeGetQuad(cvode_mem, &time, q);
  if (check_flag(&flag, "CVodeGetQuad", 1)) return(1);

  printf("--------------------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("G:          %12.4Le \n",Ith(q,1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("G:          %12.4e \n",Ith(q,1));
#else
  printf("G:          %12.4e \n",Ith(q,1));
#endif
  printf("--------------------------------------------------------\n\n");

  /* Test check point linked list 
     (uncomment next block to print check point information) */
  
  /*
  {
    int i;
    
    printf("\nList of Check Points (ncheck = %d)\n\n", ncheck);
    ckpnt = (CVadjCheckPointRec *) malloc ( (ncheck+1)*sizeof(CVadjCheckPointRec));
    CVodeGetAdjCheckPointsInfo(cvode_mem, ckpnt);
    for (i=0;i<=ncheck;i++) {
      printf("Address:       %p\n",ckpnt[i].my_addr);
      printf("Next:          %p\n",ckpnt[i].next_addr);
      printf("Time interval: %le  %le\n",ckpnt[i].t0, ckpnt[i].t1);
      printf("Step number:   %ld\n",ckpnt[i].nstep);
      printf("Order:         %d\n",ckpnt[i].order);
      printf("Step size:     %le\n",ckpnt[i].step);
      printf("\n");
    }
    
  }
  */
  
  /* Initialize yB */
  yB = N_VNew_Serial(NEQ);
  if (check_flag((void *)yB, "N_VNew_Serial", 0)) return(1);
  Ith(yB,1) = ZERO;
  Ith(yB,2) = ZERO;
  Ith(yB,3) = ZERO;

  /* Initialize qB */
  qB = N_VNew_Serial(NP);
  if (check_flag((void *)qB, "N_VNew", 0)) return(1);
  Ith(qB,1) = ZERO;
  Ith(qB,2) = ZERO;
  Ith(qB,3) = ZERO;

  /* Set the scalar relative tolerance reltolB */
  reltolB = RTOL;               

  /* Set the scalar absolute tolerance abstolB */
  abstolB = ATOLl;

  /* Set the scalar absolute tolerance abstolQB */
  abstolQB = ATOLq;

  /* Create and allocate CVODES memory for backward run */
  printf("Create and allocate CVODES memory for backward run\n");

  /* Call CVodeCreateB to specify the solution method for the backward 
     problem. */
  flag = CVodeCreateB(cvode_mem, CV_BDF, CV_NEWTON, &indexB);
  if (check_flag(&flag, "CVodeCreateB", 1)) return(1);

  /* Call CVodeInitB to allocate internal memory and initialize the 
     backward problem. */
  flag = CVodeInitB(cvode_mem, indexB, fB, TB1, yB);
  if (check_flag(&flag, "CVodeInitB", 1)) return(1);

  /* Set the scalar relative and absolute tolerances. */
  flag = CVodeSStolerancesB(cvode_mem, indexB, reltolB, abstolB);
  if (check_flag(&flag, "CVodeSStolerancesB", 1)) return(1);

  /* Attach the user data for backward problem. */
  flag = CVodeSetUserDataB(cvode_mem, indexB, data);
  if (check_flag(&flag, "CVodeSetUserDataB", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solves */
  AB = SUNDenseMatrix(NEQ, NEQ);
  if (check_flag((void *)AB, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LSB = SUNDenseLinearSolver(yB, AB);
  if (check_flag((void *)LSB, "SUNDenseLinearSolver", 0)) return(1);

  /* Attach the matrix and linear solver */
  flag = CVDlsSetLinearSolverB(cvode_mem, indexB, LSB, AB);
  if (check_flag(&flag, "CVDlsSetLinearSolverB", 1)) return(1);

  /* Set the user-supplied Jacobian routine JacB */
  flag = CVDlsSetJacFnB(cvode_mem, indexB, JacB);
  if (check_flag(&flag, "CVDlsSetJacFnB", 1)) return(1);

  /* Call CVodeQuadInitB to allocate internal memory and initialize backward
     quadrature integration. */
  flag = CVodeQuadInitB(cvode_mem, indexB, fQB, qB);
  if (check_flag(&flag, "CVodeQuadInitB", 1)) return(1);

  /* Call CVodeSetQuadErrCon to specify whether or not the quadrature variables
     are to be used in the step size control mechanism within CVODES. Call
     CVodeQuadSStolerances or CVodeQuadSVtolerances to specify the integration
     tolerances for the quadrature variables. */
  flag = CVodeSetQuadErrConB(cvode_mem, indexB, SUNTRUE);
  if (check_flag(&flag, "CVodeSetQuadErrConB", 1)) return(1);

  /* Call CVodeQuadSStolerancesB to specify the scalar relative and absolute tolerances
     for the backward problem. */
  flag = CVodeQuadSStolerancesB(cvode_mem, indexB, reltolB, abstolQB);
  if (check_flag(&flag, "CVodeQuadSStolerancesB", 1)) return(1);

  /* Backward Integration */

  PrintHead(TB1);

  /* First get results at t = TBout1 */

  /* Call CVodeB to integrate the backward ODE problem. */
  flag = CVodeB(cvode_mem, TBout1, CV_NORMAL);
  if (check_flag(&flag, "CVodeB", 1)) return(1);

  /* Call CVodeGetB to get yB of the backward ODE problem. */
  flag = CVodeGetB(cvode_mem, indexB, &time, yB);
  if (check_flag(&flag, "CVodeGetB", 1)) return(1);

  /* Call CVodeGetAdjY to get the interpolated value of the forward solution
     y during a backward integration. */
  flag = CVodeGetAdjY(cvode_mem, TBout1, y);
  if (check_flag(&flag, "CVodeGetAdjY", 1)) return(1);

  PrintOutput1(time, TBout1, y, yB);

  /* Then at t = T0 */

  flag = CVodeB(cvode_mem, T0, CV_NORMAL);
  if (check_flag(&flag, "CVodeB", 1)) return(1);
  CVodeGetNumSteps(CVodeGetAdjCVodeBmem(cvode_mem, indexB), &nstB);
  printf("Done ( nst = %ld )\n", nstB);

  flag = CVodeGetB(cvode_mem, indexB, &time, yB);
  if (check_flag(&flag, "CVodeGetB", 1)) return(1);

  /* Call CVodeGetQuadB to get the quadrature solution vector after a 
     successful return from CVodeB. */
  flag = CVodeGetQuadB(cvode_mem, indexB, &time, qB);
  if (check_flag(&flag, "CVodeGetQuadB", 1)) return(1);

  flag = CVodeGetAdjY(cvode_mem, T0, y);
  if (check_flag(&flag, "CVodeGetAdjY", 1)) return(1);

  PrintOutput(time, y, yB, qB);

  /* Reinitialize backward phase (new tB0) */

  Ith(yB,1) = ZERO;
  Ith(yB,2) = ZERO;
  Ith(yB,3) = ZERO;

  Ith(qB,1) = ZERO;
  Ith(qB,2) = ZERO;
  Ith(qB,3) = ZERO;

  printf("Re-initialize CVODES memory for backward run\n");

  flag = CVodeReInitB(cvode_mem, indexB, TB2, yB);
  if (check_flag(&flag, "CVodeReInitB", 1)) return(1);

  flag = CVodeQuadReInitB(cvode_mem, indexB, qB); 
  if (check_flag(&flag, "CVodeQuadReInitB", 1)) return(1);

  PrintHead(TB2);

  /* First get results at t = TBout1 */

  flag = CVodeB(cvode_mem, TBout1, CV_NORMAL);
  if (check_flag(&flag, "CVodeB", 1)) return(1);

  flag = CVodeGetB(cvode_mem, indexB, &time, yB);
  if (check_flag(&flag, "CVodeGetB", 1)) return(1);

  flag = CVodeGetAdjY(cvode_mem, TBout1, y);
  if (check_flag(&flag, "CVodeGetAdjY", 1)) return(1);

  PrintOutput1(time, TBout1, y, yB);

  /* Then at t = T0 */

  flag = CVodeB(cvode_mem, T0, CV_NORMAL);
  if (check_flag(&flag, "CVodeB", 1)) return(1);
  CVodeGetNumSteps(CVodeGetAdjCVodeBmem(cvode_mem, indexB), &nstB);
  printf("Done ( nst = %ld )\n", nstB);

  flag = CVodeGetB(cvode_mem, indexB, &time, yB);
  if (check_flag(&flag, "CVodeGetB", 1)) return(1);

  flag = CVodeGetQuadB(cvode_mem, indexB, &time, qB);
  if (check_flag(&flag, "CVodeGetQuadB", 1)) return(1);

  flag = CVodeGetAdjY(cvode_mem, T0, y);
  if (check_flag(&flag, "CVodeGetAdjY", 1)) return(1);

  PrintOutput(time, y, yB, qB);

  /* Free memory */
  printf("Free memory\n\n");

  CVodeFree(&cvode_mem);
  N_VDestroy(y); 
  N_VDestroy(q);
  N_VDestroy(yB);
  N_VDestroy(qB);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNLinSolFree(LSB);
  SUNMatDestroy(AB);

  if (ckpnt != NULL) free(ckpnt);
  free(data);

  return(0);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y). 
*/

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, y3, yd1, yd3;
  UserData data;
  realtype p1, p2, p3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) user_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  yd1 = Ith(ydot,1) = -p1*y1 + p2*y2*y3;
  yd3 = Ith(ydot,3) = p3*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;

  return(0);
}

/* 
 * Jacobian routine. Compute J(t,y). 
*/

static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;
  UserData data;
  realtype p1, p2, p3;
 
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) user_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
 
  IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
  IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
  IJth(J,3,1) = ZERO; IJth(J,3,2) = 2*p3*y2;        IJth(J,3,3) = ZERO;

  return(0);
}

/* 
 * fQ routine. Compute fQ(t,y). 
*/

static int fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data)
{
  Ith(qdot,1) = Ith(y,3);  

  return(0);
}
 
/*
 * EwtSet function. Computes the error weights at the current solution.
 */

static int ewt(N_Vector y, N_Vector w, void *user_data)
{
  int i;
  realtype yy, ww, rtol, atol[3];

  rtol    = RTOL;
  atol[0] = ATOL1;
  atol[1] = ATOL2;
  atol[2] = ATOL3;

  for (i=1; i<=3; i++) {
    yy = Ith(y,i);
    ww = rtol * SUNRabs(yy) + atol[i-1];
    if (ww <= 0.0) return (-1);
    Ith(w,i) = 1.0/ww;
  }

  return(0);
}

/* 
 * fB routine. Compute fB(t,y,yB). 
*/

static int fB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_dataB)
{
  UserData data;
  realtype y1, y2, y3;
  realtype p1, p2, p3;
  realtype l1, l2, l3;
  realtype l21, l32, y23;
  
  data = (UserData) user_dataB;

  /* The p vector */
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  /* The y vector */
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  
  /* The lambda vector */
  l1 = Ith(yB,1); l2 = Ith(yB,2); l3 = Ith(yB,3);

  /* Temporary variables */
  l21 = l2-l1;
  l32 = l3-l2;
  y23 = y2*y3;

  /* Load yBdot */
  Ith(yBdot,1) = - p1*l21;
  Ith(yBdot,2) = p2*y3*l21 - RCONST(2.0)*p3*y2*l32;
  Ith(yBdot,3) = p2*y2*l21 - RCONST(1.0);

  return(0);
}

/* 
 * JacB routine. Compute JB(t,y,yB). 
 */

static int JacB(realtype t, N_Vector y, N_Vector yB, N_Vector fyB, SUNMatrix JB,
                void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  UserData data;
  realtype y1, y2, y3;
  realtype p1, p2, p3;
  
  data = (UserData) user_dataB;

  /* The p vector */
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  /* The y vector */
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  /* Load JB */
  IJth(JB,1,1) = p1;     IJth(JB,1,2) = -p1;             IJth(JB,1,3) = ZERO;
  IJth(JB,2,1) = -p2*y3; IJth(JB,2,2) = p2*y3+2.0*p3*y2; IJth(JB,2,3) = RCONST(-2.0)*p3*y2;
  IJth(JB,3,1) = -p2*y2; IJth(JB,3,2) = p2*y2;           IJth(JB,3,3) = ZERO;

  return(0);
}

/*
 * fQB routine. Compute integrand for quadratures 
*/

static int fQB(realtype t, N_Vector y, N_Vector yB, 
               N_Vector qBdot, void *user_dataB)
{
  realtype y1, y2, y3;
  realtype l1, l2, l3;
  realtype l21, l32, y23;

  /* The y vector */
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  
  /* The lambda vector */
  l1 = Ith(yB,1); l2 = Ith(yB,2); l3 = Ith(yB,3);

  /* Temporary variables */
  l21 = l2-l1;
  l32 = l3-l2;
  y23 = y2*y3;

  Ith(qBdot,1) = y1*l21;
  Ith(qBdot,2) = - y23*l21;
  Ith(qBdot,3) = y2*y2*l32;

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Print heading for backward integration
 */

static void PrintHead(realtype tB0)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Backward integration from tB0 = %12.4Le\n\n",tB0);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Backward integration from tB0 = %12.4e\n\n",tB0);
#else
  printf("Backward integration from tB0 = %12.4e\n\n",tB0);
#endif
}

/*
 * Print intermediate results during backward integration
 */

static void PrintOutput1(realtype time, realtype t, N_Vector y, N_Vector yB)
{
  printf("--------------------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("returned t: %12.4Le\n",time);
  printf("tout:       %12.4Le\n",t);
  printf("lambda(t):  %12.4Le %12.4Le %12.4Le\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
  printf("y(t):       %12.4Le %12.4Le %12.4Le\n", 
         Ith(y,1), Ith(y,2), Ith(y,3));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("returned t: %12.4e\n",time);
  printf("tout:       %12.4e\n",t);
  printf("lambda(t):  %12.4e %12.4e %12.4e\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
  printf("y(t):       %12.4e %12.4e %12.4e\n", 
         Ith(y,1), Ith(y,2), Ith(y,3));
#else
  printf("returned t: %12.4e\n",time);
  printf("tout:       %12.4e\n",t);
  printf("lambda(t):  %12.4e %12.4e %12.4e\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
  printf("y(t)      : %12.4e %12.4e %12.4e\n", 
         Ith(y,1), Ith(y,2), Ith(y,3));
#endif
  printf("--------------------------------------------------------\n\n");
}

/*
 * Print final results of backward integration
 */

static void PrintOutput(realtype tfinal, N_Vector y, N_Vector yB, N_Vector qB)
{
  printf("--------------------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("returned t: %12.4Le\n",tfinal);
  printf("lambda(t0): %12.4Le %12.4Le %12.4Le\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
  printf("y(t0):      %12.4Le %12.4Le %12.4Le\n", 
         Ith(y,1), Ith(y,2), Ith(y,3));
  printf("dG/dp:      %12.4Le %12.4Le %12.4Le\n", 
         -Ith(qB,1), -Ith(qB,2), -Ith(qB,3));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("returned t: %12.4e\n",tfinal);
  printf("lambda(t0): %12.4e %12.4e %12.4e\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
  printf("y(t0):      %12.4e %12.4e %12.4e\n", 
         Ith(y,1), Ith(y,2), Ith(y,3));
  printf("dG/dp:      %12.4e %12.4e %12.4e\n", 
         -Ith(qB,1), -Ith(qB,2), -Ith(qB,3));
#else
  printf("returned t: %12.4e\n",tfinal);
  printf("lambda(t0): %12.4e %12.4e %12.4e\n", 
         Ith(yB,1), Ith(yB,2), Ith(yB,3));
  printf("y(t0)     : %12.4e %12.4e %12.4e\n", 
         Ith(y,1), Ith(y,2), Ith(y,3));
  printf("dG/dp:      %12.4e %12.4e %12.4e\n", 
         -Ith(qB,1), -Ith(qB,2), -Ith(qB,3));
#endif
  printf("--------------------------------------------------------\n\n");
}

/* 
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
