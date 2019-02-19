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
 *
 * Hessian using adjoint sensitivity example problem. 
 * 
 * This simple example problem for IDAS, due to Robertson, 
 * is from chemical kinetics, and consists of the following three 
 * equations:
 *
 *   [ y1' + p1 * y1 - p2 * y2 * y3              = 0 
 *   [ y2' - p1 * y1 + p2 * y2 * y3 + p3 * y2^2  = 0 
 *   [ y1 + y2 + y3 -1                               = 0 
 * 
 *        [1]        [-p1]
 *   y(0)=[0]  y'(0)=[ p1]   p1 = 0.04   p2 = 1e4   p3 = 1e07   
 *        [0]        [ 0 ]
 *
 *       80
 *      / 
 *  G = | 0.5 * (y1^2 + y2^2 + y3^2) dt
 *      /
 *      0
 * Compute the gradient (using FSA and ASA) and Hessian (FSA over ASA)
 * of G with respect to parameters p1 and p2.
 *
 * Reference: D.B. Ozyurt and P.I. Barton, SISC 26(5) 1725-1743, 2005.
 *
 * Error handling was suppressed for code readibility reasons.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <idas/idas.h>                 /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* defs. of SUNRabs, SUNRexp, etc.      */

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i= 1..NEQ */

/* Problem Constants */
#define NEQ      3             /* number of equations                  */
#define NP       2             /* number of sensitivities              */

#define T0       RCONST(0.0)   /* Initial time. */
#define TF       RCONST(80.0)  /* Final time. */

/* Tolerances */
#define RTOL     RCONST(1e-08) /* scalar relative tolerance            */
#define ATOL     RCONST(1e-10) /* vector absolute tolerance components */
#define RTOLA    RCONST(1e-08) /* for adjoint integration              */
#define ATOLA    RCONST(1e-08) /* for adjoint integration              */

/* Parameters */
#define P1 RCONST(0.04)
#define P2 RCONST(1.0e4)
#define P3 RCONST(3.0e7)

/* Predefined consts */
#define HALF RCONST(0.5)
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/* User defined struct */
typedef struct { 
  realtype p[3];
} *UserData;

/* residual for forward problem */
static int res(realtype t, N_Vector yy, N_Vector yp, 
               N_Vector resval, void *user_data);

static int resS(int Ns, realtype t, 
                N_Vector yy, N_Vector yp, N_Vector resval,
                N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int rhsQ(realtype t, N_Vector yy, N_Vector yp, N_Vector qdot, void *user_data);

static int rhsQS(int Ns, realtype t, N_Vector yy, N_Vector yp, 
                 N_Vector *yyS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS,
                 void *user_data,  N_Vector yytmp, N_Vector yptmp, N_Vector tmpQS);

static int resBS1(realtype tt, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS,
                  N_Vector yyB, N_Vector ypB, N_Vector resvalBQ, void *user_dataB);


static int rhsQBS1(realtype tt, N_Vector yy, N_Vector yp,
                   N_Vector *yyS, N_Vector *ypS, N_Vector yyB, N_Vector ypB,
                   N_Vector rhsBQS, void *user_dataB);

static int resBS2(realtype tt, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS,
                  N_Vector yyB, N_Vector ypB, N_Vector resvalBQ, void *user_dataB);


static int rhsQBS2(realtype tt, N_Vector yy, N_Vector yp,
                   N_Vector *yyS, N_Vector *ypS, N_Vector yyB, N_Vector ypB,
                   N_Vector rhsBQS, void *user_dataB);

static int check_retval(void *returnvalue, const char *funcname, int opt);


int main(int argc, char *argv[])
{
  N_Vector yy, yp, q, *yyS, *ypS, *qS; 
  N_Vector yyB1, ypB1, qB1, yyB2, ypB2, qB2;
  void *ida_mem;
  UserData data;
  realtype time, ti, tf;
  int retval, nckp, indexB1, indexB2;
  realtype G, Gm, Gp, dp1, dp2, grdG_fwd[2], grdG_bck[2], grdG_cntr[2], H11, H22;
  realtype rtolFD, atolFD;
  SUNMatrix A, AB1, AB2;
  SUNLinearSolver LS, LSB1, LSB2;

  /* Print problem description */
  printf("\nAdjoint Sensitivity Example for Chemical Kinetics\n");
  printf("---------------------------------------------------------\n");
  printf("DAE: dy1/dt + p1*y1 - p2*y2*y3 = 0\n");
  printf("     dy2/dt - p1*y1 + p2*y2*y3 + p3*(y2)^2 = 0\n");
  printf("               y1  +  y2  +  y3 = 0\n\n");
  printf("Find dG/dp and d^2G/dp^2, where p=[p1,p2] for\n");
  printf("     G = int_t0^tB0 g(t,p,y) dt\n");
  printf("     g(t,p,y) = y3\n\n\n");

  /* Alocate and initialize user data. */
  data = (UserData) malloc(sizeof(*data));
  data->p[0] = P1; data->p[1] = P2; data->p[2] = P3;

  /* Consistent IC */
  yy = N_VNew_Serial(NEQ);
  yp = N_VNew_Serial(NEQ);
  Ith(yy,1) = ONE; Ith(yy,2) = ZERO; Ith(yy,3) = ZERO;
  Ith(yp,1) = -P1; Ith(yp,2) = P1; Ith(yp,3) = 0;

  q = N_VNew_Serial(1);
  N_VConst(ZERO, q);

  yyS = N_VCloneVectorArray(NP, yy);
  ypS = N_VCloneVectorArray(NP, yp);
  N_VConst(ZERO, yyS[0]); N_VConst(ZERO, yyS[1]);
  N_VConst(ZERO, ypS[0]); N_VConst(ZERO, ypS[1]);

  qS = N_VCloneVectorArray(NP, q);
  N_VConst(ZERO, qS[0]);

  ida_mem = IDACreate();

  ti = T0;
  retval = IDAInit(ida_mem, res, ti, yy, yp);

  /* Forward problem's setup. */
  retval = IDASStolerances(ida_mem, RTOL, ATOL);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(ida_mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  retval = IDASetUserData(ida_mem, data);
  retval = IDASetMaxNumSteps(ida_mem, 1500);

  /* Quadrature's setup. */
  retval = IDAQuadInit(ida_mem, rhsQ, q);
  retval = IDAQuadSStolerances(ida_mem, RTOL, ATOL);
  retval = IDASetQuadErrCon(ida_mem, SUNTRUE);

  /* Sensitivity's setup. */
  retval = IDASensInit(ida_mem, NP, IDA_SIMULTANEOUS, resS, yyS, ypS);
  retval = IDASensEEtolerances(ida_mem);
  retval = IDASetSensErrCon(ida_mem, SUNTRUE);

  /* Setup of quadrature's sensitivities */
  retval = IDAQuadSensInit(ida_mem, rhsQS, qS);
  retval = IDAQuadSensEEtolerances(ida_mem);
  retval = IDASetQuadSensErrCon(ida_mem, SUNTRUE); 
  
  /* Initialize ASA. */
  retval = IDAAdjInit(ida_mem, 100, IDA_HERMITE);

  printf("---------------------------------------------------------\n");
  printf("Forward integration\n");
  printf("---------------------------------------------------------\n\n");

  tf = TF;
  retval = IDASolveF(ida_mem, tf, &time, yy, yp, IDA_NORMAL, &nckp);

  IDAGetQuad(ida_mem, &time, q);
  G = Ith(q,1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("     G:    %12.4Le\n", G);
#else
  printf("     G:    %12.4e\n", G);
#endif  

  /* Sensitivities are needed for IC of backward problems. */
  IDAGetSensDky(ida_mem, tf, 0, yyS);
  IDAGetSensDky(ida_mem, tf, 1, ypS);

  IDAGetQuadSens(ida_mem, &time, qS);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le\n", Ith(qS[0],1), Ith(qS[1],1));
#else
  printf("   dG/dp:  %12.4e %12.4e\n", Ith(qS[0],1), Ith(qS[1],1));
#endif  
  printf("\n");
  /******************************
  * BACKWARD PROBLEM #1
  *******************************/

  /* Consistent IC. */
  yyB1 = N_VNew_Serial(2*NEQ);
  ypB1 = N_VNew_Serial(2*NEQ);

  N_VConst(ZERO, yyB1);
  Ith(yyB1,3) = Ith(yy,3);
  Ith(yyB1,6) = Ith(yyS[0], 3);

  N_VConst(ZERO, ypB1);
  Ith(ypB1,1) = Ith(yy,3)-Ith(yy,1);
  Ith(ypB1,2) = Ith(yy,3)-Ith(yy,2);
  Ith(ypB1,4) = Ith(yyS[0],3) - Ith(yyS[0],1);
  Ith(ypB1,5) = Ith(yyS[0],3) - Ith(yyS[0],2);
  
  qB1 = N_VNew_Serial(2*NP);
  N_VConst(ZERO, qB1);

  retval = IDACreateB(ida_mem, &indexB1);
  retval = IDAInitBS(ida_mem, indexB1, resBS1, tf, yyB1, ypB1);
  retval = IDASStolerancesB(ida_mem, indexB1, RTOLA, ATOLA);   
  retval = IDASetUserDataB(ida_mem, indexB1, data);
  retval = IDASetMaxNumStepsB(ida_mem, indexB1, 5000);

  /* Create dense SUNMatrix for use in linear solves */
  AB1 = SUNDenseMatrix(2*NEQ, 2*NEQ);
  if(check_retval((void *)AB1, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LSB1 = SUNLinSol_Dense(yyB1, AB1);
  if(check_retval((void *)LSB1, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolverB(ida_mem, indexB1, LSB1, AB1);
  if(check_retval(&retval, "IDASetLinearSolverB", 1)) return(1);

  retval = IDAQuadInitBS(ida_mem, indexB1, rhsQBS1, qB1);

  /******************************
  * BACKWARD PROBLEM #2  
  *******************************/

  /* Consistent IC. */
  yyB2 = N_VNew_Serial(2*NEQ);
  ypB2 = N_VNew_Serial(2*NEQ);

  N_VConst(ZERO, yyB2);
  Ith(yyB2,3) = Ith(yy,3);
  Ith(yyB2,6) = Ith(yyS[1],3);

  N_VConst(ZERO, ypB2);
  Ith(ypB2,1) = Ith(yy,3)-Ith(yy,1);
  Ith(ypB2,2) = Ith(yy,3)-Ith(yy,2);
  Ith(ypB2,4) = Ith(yyS[1],3) - Ith(yyS[1],1);
  Ith(ypB2,5) = Ith(yyS[1],3) - Ith(yyS[1],2);
  
  qB2 = N_VNew_Serial(2*NP);
  N_VConst(ZERO, qB2);

  retval = IDACreateB(ida_mem, &indexB2);
  retval = IDAInitBS(ida_mem, indexB2, resBS2, tf, yyB2, ypB2);
  retval = IDASStolerancesB(ida_mem, indexB2, RTOLA, ATOLA);   
  retval = IDASetUserDataB(ida_mem, indexB2, data);
  retval = IDASetMaxNumStepsB(ida_mem, indexB2, 2500);

  /* Create dense SUNMatrix for use in linear solves */
  AB2 = SUNDenseMatrix(2*NEQ, 2*NEQ);
  if(check_retval((void *)AB2, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LSB2 = SUNLinSol_Dense(yyB2, AB2);
  if(check_retval((void *)LSB2, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolverB(ida_mem, indexB2, LSB2, AB2);
  if(check_retval(&retval, "IDASetLinearSolverB", 1)) return(1);

  retval = IDAQuadInitBS(ida_mem, indexB2, rhsQBS2, qB2);

  /* Integrate backward problems. */
  printf("---------------------------------------------------------\n");
  printf("Backward integration \n");
  printf("---------------------------------------------------------\n\n"); 

  retval = IDASolveB(ida_mem, ti, IDA_NORMAL);

  retval = IDAGetB(ida_mem, indexB1, &time, yyB1, ypB1);
  /* 
     retval = IDAGetNumSteps(IDAGetAdjIDABmem(ida_mem, indexB1), &nst);
     printf("at time=%g \tpb 1 Num steps:%d\n", time, nst); 
     retval = IDAGetNumSteps(IDAGetAdjIDABmem(ida_mem, indexB2), &nst);
     printf("at time=%g \tpb 2 Num steps:%d\n\n", time, nst); 
  */

  retval = IDAGetQuadB(ida_mem, indexB1, &time, qB1);
  retval = IDAGetQuadB(ida_mem, indexB2, &time, qB2);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le   (from backward pb. 1)\n", Ith(qB1,1), Ith(qB1,2));
  printf("   dG/dp:  %12.4Le %12.4Le   (from backward pb. 2)\n", Ith(qB2,1), Ith(qB2,2));
#else
  printf("   dG/dp:  %12.4e %12.4e   (from backward pb. 1)\n", Ith(qB1,1), Ith(qB1,2));
  printf("   dG/dp:  %12.4e %12.4e   (from backward pb. 2)\n", Ith(qB2,1), Ith(qB2,2));
#endif  

  printf("\n");
  printf("   H = d2G/dp2:\n");
  printf("        (1)            (2)\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("  %12.4Le  %12.4Le\n", Ith(qB1,3), Ith(qB2,3));
  printf("  %12.4Le  %12.4Le\n", Ith(qB1,4), Ith(qB2,4));
#else
  printf("  %12.4e  %12.4e\n", Ith(qB1,3), Ith(qB2,3));
  printf("  %12.4e  %12.4e\n", Ith(qB1,4), Ith(qB2,4));
#endif  

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNLinSolFree(LSB1);
  SUNMatDestroy(AB1);
  SUNLinSolFree(LSB2);
  SUNMatDestroy(AB2);


  /*********************************
  * Use Finite Differences to verify
  **********************************/

  /* Perturbations are of different magnitudes as p1 and p2 are. */
  dp1 = RCONST(1.0e-3);
  dp2 = RCONST(2.5e+2);

  printf("\n");
  printf("---------------------------------------------------------\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Finite Differences ( dp1=%6.1Le and dp2 = %6.1Le )\n", dp1, dp2);
#else
  printf("Finite Differences ( dp1=%6.1e and dp2 = %6.1e )\n", dp1, dp2);
#endif  
  printf("---------------------------------------------------------\n\n");

  ida_mem = IDACreate();

  /********************
  * Forward FD for p1
  ********************/
  data->p[0] += dp1;

  Ith(yy,1) = ONE; Ith(yy,2) = ZERO; Ith(yy,3) = ZERO;
  Ith(yp,1) = -data->p[0]; Ith(yp,2) = -Ith(yp,1); Ith(yp,3) = 0;
  N_VConst(ZERO, q);
  ti = T0;
  tf = TF;

  retval = IDAInit(ida_mem, res, ti, yy, yp);

  rtolFD = RCONST(1.0e-12);
  atolFD = RCONST(1.0e-14);

  retval = IDASStolerances(ida_mem, rtolFD, atolFD);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ);
  if(check_retval((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object */
  LS = SUNLinSol_Dense(yy, A);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(ida_mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  retval = IDASetUserData(ida_mem, data);
  retval = IDASetMaxNumSteps(ida_mem, 10000);

  retval = IDAQuadInit(ida_mem, rhsQ, q);
  retval = IDAQuadSStolerances(ida_mem, rtolFD, atolFD);
  retval = IDASetQuadErrCon(ida_mem, SUNTRUE);

  retval = IDASolve(ida_mem, tf, &time, yy, yp, IDA_NORMAL);
  retval = IDAGetQuad(ida_mem, &time, q);
  Gp = Ith(q,1);

  /********************
  * Backward FD for p1
  ********************/
  data->p[0] -= 2*dp1;

  Ith(yy,1) = ONE; Ith(yy,2) = ZERO; Ith(yy,3) = ZERO;
  Ith(yp,1) = -data->p[0]; Ith(yp,2) = -Ith(yp,1); Ith(yp,3) = 0;
  N_VConst(ZERO, q);
  
  retval = IDAReInit(ida_mem, ti, yy, yp);
  retval = IDAQuadReInit(ida_mem, q);

  retval = IDASolve(ida_mem, tf, &time, yy, yp, IDA_NORMAL);
  retval = IDAGetQuad(ida_mem, &time, q);
  Gm = Ith(q,1);

  /* Compute FD for p1. */
  grdG_fwd[0] = (Gp-G)/dp1;
  grdG_bck[0] = (G-Gm)/dp1;
  grdG_cntr[0] = (Gp-Gm)/(2.0*dp1);
  H11 = (Gp - 2.0*G + Gm) / (dp1*dp1);

  /********************
  * Forward FD for p2
  ********************/
  /*restore p1*/
  data->p[0] += dp1; 
  data->p[1] += dp2;

  Ith(yy,1) = ONE; Ith(yy,2) = ZERO; Ith(yy,3) = ZERO;
  Ith(yp,1) = -data->p[0]; Ith(yp,2) = -Ith(yp,1); Ith(yp,3) = 0;
  N_VConst(ZERO, q);

  retval = IDAReInit(ida_mem, ti, yy, yp);
  retval = IDAQuadReInit(ida_mem, q);

  retval = IDASolve(ida_mem, tf, &time, yy, yp, IDA_NORMAL);
  retval = IDAGetQuad(ida_mem, &time, q);
  Gp = Ith(q,1);

  /********************
  * Backward FD for p2
  ********************/
  data->p[1] -= 2*dp2;

  Ith(yy,1) = ONE; Ith(yy,2) = ZERO; Ith(yy,3) = ZERO;
  Ith(yp,1) = -data->p[0]; Ith(yp,2) = -Ith(yp,1); Ith(yp,3) = 0;
  N_VConst(ZERO, q);
  
  retval = IDAReInit(ida_mem, ti, yy, yp);
  retval = IDAQuadReInit(ida_mem, q);

  retval = IDASolve(ida_mem, tf, &time, yy, yp, IDA_NORMAL);
  retval = IDAGetQuad(ida_mem, &time, q);
  Gm = Ith(q,1);

  /* Compute FD for p2. */
  grdG_fwd[1] = (Gp-G)/dp2;
  grdG_bck[1] = (G-Gm)/dp2;
  grdG_cntr[1] = (Gp-Gm)/(2.0*dp2);
  H22 = (Gp - 2.0*G + Gm) / (dp2*dp2);


  printf("\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le  %12.4Le   (fwd FD)\n",  grdG_fwd[0],  grdG_fwd[1]);
  printf("           %12.4Le  %12.4Le   (bck FD)\n",  grdG_bck[0],  grdG_bck[1]);
  printf("           %12.4Le  %12.4Le   (cntr FD)\n", grdG_cntr[0], grdG_cntr[1]);
  printf("\n");
  printf("  H(1,1):  %12.4Le\n", H11);
  printf("  H(2,2):  %12.4Le\n", H22);
#else
  printf("   dG/dp:  %12.4e  %12.4e   (fwd FD)\n",  grdG_fwd[0],  grdG_fwd[1]);
  printf("           %12.4e  %12.4e   (bck FD)\n",  grdG_bck[0],  grdG_bck[1]);
  printf("           %12.4e  %12.4e   (cntr FD)\n", grdG_cntr[0], grdG_cntr[1]);
  printf("\n");
  printf("  H(1,1):  %12.4e\n", H11);
  printf("  H(2,2):  %12.4e\n", H22);
#endif  

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);

  N_VDestroy(yyB1);
  N_VDestroy(ypB1);
  N_VDestroy(qB1);

  N_VDestroy(yyB2);
  N_VDestroy(ypB2);
  N_VDestroy(qB2);

  N_VDestroy(yy);
  N_VDestroy(yp);
  N_VDestroy(q);
  N_VDestroyVectorArray(yyS, NP);
  N_VDestroyVectorArray(ypS, NP);
  N_VDestroyVectorArray(qS, NP);

  free(data);
  return 0;
}



static int res(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
  realtype y1, y2, y3, yp1, yp2, *rval;
  UserData data;
  realtype p1, p2, p3;

  y1  = Ith(yy,1); y2  = Ith(yy,2); y3  = Ith(yy,3); 
  yp1 = Ith(yp,1); yp2 = Ith(yp,2);
  rval = N_VGetArrayPointer(rr);

  data = (UserData) user_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  rval[0] = p1*y1-p2*y2*y3;
  rval[1] = -rval[0] + p3*y2*y2 + yp2;
  rval[0]+= yp1;
  rval[2] = y1+y2+y3-1;

  return(0);
}

static int resS(int Ns, realtype t, 
                N_Vector yy, N_Vector yp, N_Vector resval,
                N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                void *user_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;
  realtype s1, s2, s3;
  realtype sd1, sd2;
  realtype rs1, rs2, rs3;
  int is;

  data = (UserData) user_data;
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  y1 = Ith(yy,1);
  y2 = Ith(yy,2);
  y3 = Ith(yy,3);

  for (is=0; is<NP; is++) {

    s1 = Ith(yyS[is],1);
    s2 = Ith(yyS[is],2);
    s3 = Ith(yyS[is],3);

    sd1 = Ith(ypS[is],1);
    sd2 = Ith(ypS[is],2);

    rs1 = sd1 + p1*s1 - p2*y3*s2 - p2*y2*s3;
    rs2 = sd2 - p1*s1 + p2*y3*s2 + p2*y2*s3 + TWO*p3*y2*s2;
    rs3 = s1 + s2 + s3;

    switch (is) {
    case 0:
      rs1 += y1;
      rs2 -= y1;
      break;
    case 1:
      rs1 -= y2*y3;
      rs2 += y2*y3;
      break;
    }
  
    Ith(resvalS[is],1) = rs1;
    Ith(resvalS[is],2) = rs2;
    Ith(resvalS[is],3) = rs3;

  }

  return(0);
}

static int rhsQ(realtype t, N_Vector yy, N_Vector yp, N_Vector qdot, void *user_data)
{

  realtype y1, y2, y3;

  y1 = Ith(yy,1); y2 = Ith(yy,2); y3 = Ith(yy,3); 
  Ith(qdot,1) = HALF*(y1*y1+y2*y2+y3*y3);

  return(0);
}

static int rhsQS(int Ns, realtype t,
                 N_Vector yy, N_Vector yp, 
                 N_Vector *yyS, N_Vector *ypS, 
                 N_Vector rrQ, N_Vector *rhsvalQS,
                 void *user_data,
                 N_Vector yytmp, N_Vector yptmp, N_Vector tmpQS)
{

  realtype y1, y2, y3;
  realtype s1, s2, s3;

  y1 = Ith(yy,1); 
  y2 = Ith(yy,2); 
  y3 = Ith(yy,3);

  /* 1st sensitivity RHS */
  s1 = Ith(yyS[0],1);
  s2 = Ith(yyS[0],2);
  s3 = Ith(yyS[0],3);
  Ith(rhsvalQS[0],1) = y1*s1 + y2*s2 + y3*s3;

  /* 2nd sensitivity RHS */
  s1 = Ith(yyS[1],1);
  s2 = Ith(yyS[1],2);
  s3 = Ith(yyS[1],3);
  Ith(rhsvalQS[1],1) = y1*s1 + y2*s2 + y3*s3;

  return(0);
}

/* Residuals for adjoint model. */
static int resBS1(realtype tt, 
                  N_Vector yy, N_Vector yp, 
                  N_Vector *yyS, N_Vector *ypS,
                  N_Vector yyB, N_Vector ypB,
                  N_Vector rrBS, void *user_dataB)

{
  UserData data;
  realtype y1, y2, y3;
  realtype p1, p2, p3;
  realtype l1, l2, l3, m1, m2, m3;
  realtype lp1, lp2, mp1, mp2;
  realtype s1, s2, s3;
  realtype l21;
  
  data = (UserData) user_dataB;

  /* The parameters. */
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  /* The y vector. */
  y1 = Ith(yy,1); y2 = Ith(yy,2); y3 = Ith(yy,3);

  /* The lambda vector. */
  l1 = Ith(yyB,1); l2 = Ith(yyB,2); l3 = Ith(yyB,3);
  /* The mu vector. */
  m1 = Ith(yyB,4); m2 = Ith(yyB,5); m3 = Ith(yyB,6);

  /* The lambda dot vector. */
  lp1 = Ith(ypB,1); lp2 = Ith(ypB,2);
  /* The mu dot vector. */
  mp1 = Ith(ypB,4); mp2 = Ith(ypB,5);

  /* The sensitivity with respect to p1 */
  s1 = Ith(yyS[0],1); s2 = Ith(yyS[0],2); s3 = Ith(yyS[0],3);

  /* Temporary variables */
  l21 = l2-l1;

  Ith(rrBS,1) = lp1 + p1*l21 - l3 + y1;
  Ith(rrBS,2) = lp2 - p2*y3*l21 - TWO*p3*y2*l2 - l3 + y2;
  Ith(rrBS,3) = -p2*y2*l21 - l3 + y3;

  Ith(rrBS,4) = mp1 + p1*(-m1+m2) - m3 + l21 + s1;
  Ith(rrBS,5) = mp2 + p2*y3*m1 - (p2*y3+TWO*p3*y2)*m2 - m3 + p2*s3*l1 - (TWO*p3*s2+p2*s3)*l2 + s2;
  Ith(rrBS,6) = p2*y2*(m1-m2) - m3 - p2*s2*l21 + s3;

  return(0);
}

static int rhsQBS1(realtype tt, 
                 N_Vector yy, N_Vector yp,
                 N_Vector *yyS, N_Vector *ypS,
                 N_Vector yyB, N_Vector ypB,
                 N_Vector rhsBQS, void *user_dataB)
{
  realtype y1, y2, y3;
  realtype l1, l2, m1, m2;
  realtype s1, s2, s3;
  realtype l21;

  /* The y vector */
  y1 = Ith(yy,1); y2 = Ith(yy,2); y3 = Ith(yy,3);
  
  /* The lambda vector. */
  l1 = Ith(yyB,1); l2 = Ith(yyB,2);

  /* The mu vector. */
  m1 = Ith(yyB,4); m2 = Ith(yyB,5);

  /* The sensitivity with respect to p1 */
  s1 = Ith(yyS[0],1); s2 = Ith(yyS[0],2); s3 = Ith(yyS[0],3);
  
  /* Temporary variables */
  l21 = l2-l1;

  Ith(rhsBQS,1) = -y1*l21;
  Ith(rhsBQS,2) = y2*y3*l21;

  Ith(rhsBQS,3) = y1*(m1-m2) - s1*l21; 
  Ith(rhsBQS,4) = y2*y3*(m2-m1) + (y3*s2+y2*s3)*l21;

  return(0);
}

static int resBS2(realtype tt, 
                  N_Vector yy, N_Vector yp, 
                  N_Vector *yyS, N_Vector *ypS,
                  N_Vector yyB, N_Vector ypB,
                  N_Vector rrBS, void *user_dataB)

{
  UserData data;
  realtype y1, y2, y3;
  realtype p1, p2, p3;
  realtype l1, l2, l3, m1, m2, m3;
  realtype lp1, lp2, mp1, mp2;
  realtype s1, s2, s3;
  realtype l21;
  
  data = (UserData) user_dataB;

  /* The parameters. */
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  /* The y vector. */
  y1 = Ith(yy,1); y2 = Ith(yy,2); y3 = Ith(yy,3);

  /* The lambda vector. */
  l1 = Ith(yyB,1); l2 = Ith(yyB,2); l3 = Ith(yyB,3);
  /* The mu vector. */
  m1 = Ith(yyB,4); m2 = Ith(yyB,5); m3 = Ith(yyB,6);

  /* The lambda dot vector. */
  lp1 = Ith(ypB,1); lp2 = Ith(ypB,2);

  /* The mu dot vector. */
  mp1 = Ith(ypB,4); mp2 = Ith(ypB,5);

  /* The sensitivity with respect to p2 */
  s1 = Ith(yyS[1],1); s2 = Ith(yyS[1],2); s3 = Ith(yyS[1],3);

  /* Temporary variables */
  l21 = l2-l1;

  Ith(rrBS,1) = lp1 + p1*l21 - l3 + y1;
  Ith(rrBS,2) = lp2 - p2*y3*l21 - TWO*p3*y2*l2 - l3 + y2;
  Ith(rrBS,3) = -p2*y2*l21 - l3 + y3;

  Ith(rrBS,4) = mp1 + p1*(-m1+m2) - m3 + s1;
  Ith(rrBS,5) = mp2 + p2*y3*m1 - (p2*y3+TWO*p3*y2)*m2 - m3 + (y3+p2*s3)*l1 - (y3+TWO*p3*s2+p2*s3)*l2 + s2;
  Ith(rrBS,6) = p2*y2*(m1-m2) - m3 - (y2+p2*s2)*l21 + s3;

  return(0);
}

static int rhsQBS2(realtype tt, 
                 N_Vector yy, N_Vector yp,
                 N_Vector *yyS, N_Vector *ypS,
                 N_Vector yyB, N_Vector ypB,
                 N_Vector rhsBQS, void *user_dataB)
{
  realtype y1, y2, y3;
  realtype l1, l2, m1, m2;
  realtype s1, s2, s3;
  realtype l21;

  /* The y vector */
  y1 = Ith(yy,1); y2 = Ith(yy,2); y3 = Ith(yy,3);
  
  /* The lambda vector. */
  l1 = Ith(yyB,1); l2 = Ith(yyB,2);

  /* The mu vector. */
  m1 = Ith(yyB,4); m2 = Ith(yyB,5);

  /* The sensitivity with respect to p2 */
  s1 = Ith(yyS[1],1); s2 = Ith(yyS[1],2); s3 = Ith(yyS[1],3);
  
  /* Temporary variables */
  l21 = l2-l1;

  Ith(rhsBQS,1) = -y1*l21;
  Ith(rhsBQS,2) =  y2*y3*l21;

  Ith(rhsBQS,3) = y1*(m1-m2) - s1*l21; 
  Ith(rhsBQS,4) = y2*y3*(m2-m1) + (y3*s2+y2*s3)*l21;

  return(0);
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

