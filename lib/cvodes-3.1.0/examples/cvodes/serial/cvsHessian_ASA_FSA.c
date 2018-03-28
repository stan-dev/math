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
 *
 * Hessian through adjoint sensitivity example problem.
 *
 *        [ - p1 * y1^2 - y3 ]           [ 1 ]
 *   y' = [    - y2          ]    y(0) = [ 1 ]
 *        [ -p2^2 * y2 * y3  ]           [ 1 ]
 *
 *   p1 = 1.0
 *   p2 = 2.0
 *
 *           2
 *          /
 *   G(p) = |  0.5 * ( y1^2 + y2^2 + y3^2 ) dt
 *          /
 *          0
 *
 * Compute the gradient (ASA) and Hessian (FSA over ASA) of G(p).
 *
 * See D.B. Ozyurt and P.I. Barton, SISC 26(5) 1725-1743, 2005.
 *
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes.h>              /* prototypes for CVODES fcts., consts. */
#include <nvector/nvector_serial.h>     /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>  /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_dense.h>  /* access to band SUNLinearSolver       */
#include <cvodes/cvodes_direct.h>       /* access to CVDls interface            */
#include <sundials/sundials_math.h>     /* definition of SUNRabs, SUNRexp, etc. */

#define Ith(v,i)    NV_Ith_S(v,i-1)

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

typedef struct {
  realtype p1, p2;
} *UserData;

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data);
static int fS(int Ns, realtype t,
              N_Vector y, N_Vector ydot,
              N_Vector *yS, N_Vector *ySdot,
              void *user_data,
              N_Vector tmp1, N_Vector tmp2);
static int fQS(int Ns, realtype t,
               N_Vector y, N_Vector *yS, 
               N_Vector yQdot, N_Vector *yQSdot,
               void *user_data,
               N_Vector tmp, N_Vector tmpQ);

static int fB1(realtype t, N_Vector y, N_Vector *yS, 
               N_Vector yB, N_Vector yBdot, void *user_dataB);
static int fQB1(realtype t, N_Vector y, N_Vector *yS, 
                N_Vector yB, N_Vector qBdot, void *user_dataB);


static int fB2(realtype t, N_Vector y, N_Vector *yS, 
               N_Vector yB, N_Vector yBdot, void *user_dataB);
static int fQB2(realtype t, N_Vector y, N_Vector *yS,
                N_Vector yB, N_Vector qBdot, void *user_dataB);

void PrintFwdStats(void *cvode_mem);
void PrintBckStats(void *cvode_mem, int idx);

/* Private function to check function return values */

static int check_flag(void *flagvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  UserData data;

  SUNMatrix A, AB1, AB2;
  SUNLinearSolver LS, LSB1, LSB2;
  void *cvode_mem;

  sunindextype Neq, Np2;
  int Np;

  realtype t0, tf;

  realtype reltol;
  realtype abstol, abstolQ, abstolB, abstolQB;

  N_Vector y, yQ;
  N_Vector *yS, *yQS;
  N_Vector yB1, yB2;
  N_Vector yQB1, yQB2;

  int steps, ncheck;
  int indexB1, indexB2;

  int flag;
  realtype time;

  realtype dp;
  realtype G, Gp, Gm;
  realtype grdG_fwd[2], grdG_bck[2], grdG_cntr[2];
  realtype H11, H22;

  data = NULL;
  y = yQ = NULL;
  yB1 = yB2 = NULL;
  yQB1 = yQB2 = NULL;
  A = AB1 = AB2 = NULL;
  LS = LSB1 = LSB2 = NULL;
  cvode_mem = NULL;

  /* User data structure */

  data = (UserData) malloc(sizeof *data);
  data->p1 = RCONST(1.0);
  data->p2 = RCONST(2.0);

  /* Problem size, integration interval, and tolerances */

  Neq = 3;
  Np  = 2;
  Np2 = 2*Np;

  t0 = 0.0;
  tf = 2.0;

  reltol = 1.0e-8;

  abstol = 1.0e-8;
  abstolQ = 1.0e-8;

  abstolB = 1.0e-8;
  abstolQB = 1.0e-8;

  /* Initializations for forward problem */

  y = N_VNew_Serial(Neq);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  N_VConst(ONE, y);

  yQ = N_VNew_Serial(1);
  if (check_flag((void *)yQ, "N_VNew_Serial", 0)) return(1);
  N_VConst(ZERO, yQ);

  yS = N_VCloneVectorArray(Np, y);
  if (check_flag((void *)yS, "N_VCloneVectorArray", 0)) return(1);
  N_VConst(ZERO, yS[0]);
  N_VConst(ZERO, yS[1]);

  yQS = N_VCloneVectorArray(Np, yQ);
  if (check_flag((void *)yQS, "N_VCloneVectorArray", 0)) return(1);
  N_VConst(ZERO, yQS[0]);
  N_VConst(ZERO, yQS[1]);

  /* Create and initialize forward problem */

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if(check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  flag = CVodeInit(cvode_mem, f, t0, y);
  if(check_flag(&flag, "CVodeInit", 1)) return(1);

  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  flag = CVodeSetUserData(cvode_mem, data);
  if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Create a dense SUNMatrix */
  A = SUNDenseMatrix(Neq, Neq);
  if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create banded SUNLinearSolver for the forward problem */
  LS = SUNDenseLinearSolver(y, A);
  if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

  /* Attach the matrix and linear solver */
  flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
  if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

  flag = CVodeQuadInit(cvode_mem, fQ, yQ);
  if(check_flag(&flag, "CVodeQuadInit", 1)) return(1);

  flag = CVodeQuadSStolerances(cvode_mem, reltol, abstolQ);
  if(check_flag(&flag, "CVodeQuadSStolerances", 1)) return(1);

  flag = CVodeSetQuadErrCon(cvode_mem, SUNTRUE);
  if(check_flag(&flag, "CVodeSetQuadErrCon", 1)) return(1);

  flag = CVodeSensInit(cvode_mem, Np, CV_SIMULTANEOUS, fS, yS);
  if(check_flag(&flag, "CVodeSensInit", 1)) return(1);

  flag = CVodeSensEEtolerances(cvode_mem);
  if(check_flag(&flag, "CVodeSensEEtolerances", 1)) return(1);

  flag = CVodeSetSensErrCon(cvode_mem, SUNTRUE);
  if(check_flag(&flag, "CVodeSetSensErrCon", 1)) return(1);

  flag = CVodeQuadSensInit(cvode_mem, fQS, yQS);
  if(check_flag(&flag, "CVodeQuadSensInit", 1)) return(1);

  flag = CVodeQuadSensEEtolerances(cvode_mem);
  if(check_flag(&flag, "CVodeQuadSensEEtolerances", 1)) return(1);

  flag = CVodeSetQuadSensErrCon(cvode_mem, SUNTRUE);
  if(check_flag(&flag, "CVodeSetQuadSensErrCon", 1)) return(1);

  /* Initialize ASA */

  steps = 100;
  flag = CVodeAdjInit(cvode_mem, steps, CV_POLYNOMIAL);
  if(check_flag(&flag, "CVodeAdjInit", 1)) return(1);

  /* Forward integration */

  printf("-------------------\n");
  printf("Forward integration\n");
  printf("-------------------\n\n");

  flag = CVodeF(cvode_mem, tf, y, &time, CV_NORMAL, &ncheck);
  if(check_flag(&flag, "CVodeF", 1)) return(1);

  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  if(check_flag(&flag, "CVodeGetQuad", 1)) return(1);

  G = Ith(yQ,1);

  flag = CVodeGetSens(cvode_mem, &time, yS);
  if(check_flag(&flag, "CVodeGetSens", 1)) return(1);

  flag = CVodeGetQuadSens(cvode_mem, &time, yQS);
  if(check_flag(&flag, "CVodeGetQuadSens", 1)) return(1);

  printf("ncheck = %d\n", ncheck);
  printf("\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("     y:    %12.4Le %12.4Le %12.4Le", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:    %12.4Le\n", Ith(yQ,1));
  printf("\n");
  printf("     yS1:  %12.4Le %12.4Le %12.4Le\n", Ith(yS[0],1), Ith(yS[0],2), Ith(yS[0],3));
  printf("     yS2:  %12.4Le %12.4Le %12.4Le\n", Ith(yS[1],1), Ith(yS[1],2), Ith(yS[1],3));
  printf("\n");
  printf("   dG/dp:  %12.4Le %12.4Le\n", Ith(yQS[0],1), Ith(yQS[1],1));
#else
  printf("     y:    %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:    %12.4e\n", Ith(yQ,1));
  printf("\n");
  printf("     yS1:  %12.4e %12.4e %12.4e\n", Ith(yS[0],1), Ith(yS[0],2), Ith(yS[0],3));
  printf("     yS2:  %12.4e %12.4e %12.4e\n", Ith(yS[1],1), Ith(yS[1],2), Ith(yS[1],3));
  printf("\n");
  printf("   dG/dp:  %12.4e %12.4e\n", Ith(yQS[0],1), Ith(yQS[1],1));
#endif
  printf("\n");

  printf("Final Statistics for forward pb.\n");
  printf("--------------------------------\n");
  PrintFwdStats(cvode_mem);


  /* Initializations for backward problems */

  yB1 = N_VNew_Serial(2*Neq);
  if (check_flag((void *)yB1, "N_VNew_Serial", 0)) return(1);
  N_VConst(ZERO, yB1);

  yQB1 = N_VNew_Serial(Np2);
  if (check_flag((void *)yQB1, "N_VNew_Serial", 0)) return(1);
  N_VConst(ZERO, yQB1);

  yB2 = N_VNew_Serial(2*Neq);
  if (check_flag((void *)yB2, "N_VNew_Serial", 0)) return(1);
  N_VConst(ZERO, yB2);

  yQB2 = N_VNew_Serial(Np2);
  if (check_flag((void *)yQB2, "N_VNew_Serial", 0)) return(1);
  N_VConst(ZERO, yQB2);

  /* Create and initialize backward problems (one for each column of the Hessian) */

  /* -------------------------
     First backward problem
     -------------------------*/

  flag = CVodeCreateB(cvode_mem, CV_BDF, CV_NEWTON, &indexB1);
  if(check_flag(&flag, "CVodeCreateB", 1)) return(1);

  flag = CVodeInitBS(cvode_mem, indexB1, fB1, tf, yB1);
  if(check_flag(&flag, "CVodeInitBS", 1)) return(1);

  flag = CVodeSStolerancesB(cvode_mem, indexB1, reltol, abstolB);
  if(check_flag(&flag, "CVodeSStolerancesB", 1)) return(1);

  flag = CVodeSetUserDataB(cvode_mem, indexB1, data);
  if(check_flag(&flag, "CVodeSetUserDataB", 1)) return(1);

  flag = CVodeQuadInitBS(cvode_mem, indexB1, fQB1, yQB1);
  if(check_flag(&flag, "CVodeQuadInitBS", 1)) return(1);

  flag = CVodeQuadSStolerancesB(cvode_mem, indexB1, reltol, abstolQB);
  if(check_flag(&flag, "CVodeQuadSStolerancesB", 1)) return(1);

  flag = CVodeSetQuadErrConB(cvode_mem, indexB1, SUNTRUE);
  if(check_flag(&flag, "CVodeSetQuadErrConB", 1)) return(1);

  /* Create a dense SUNMatrix */
  AB1 = SUNDenseMatrix(2*Neq, 2*Neq);
  if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create banded SUNLinearSolver for the forward problem */
  LSB1 = SUNDenseLinearSolver(yB1, AB1);
  if(check_flag((void *)LSB1, "SUNDenseLinearSolver", 0)) return(1);

  /* Attach the matrix and linear solver */
  flag = CVDlsSetLinearSolverB(cvode_mem, indexB1, LSB1, AB1);
  if(check_flag(&flag, "CVDlsSetLinearSolverB", 1)) return(1);

  /* -------------------------
     Second backward problem
     -------------------------*/

  flag = CVodeCreateB(cvode_mem, CV_BDF, CV_NEWTON, &indexB2);
  if(check_flag(&flag, "CVodeCreateB", 1)) return(1);

  flag = CVodeInitBS(cvode_mem, indexB2, fB2, tf, yB2);
  if(check_flag(&flag, "CVodeInitBS", 1)) return(1);

  flag = CVodeSStolerancesB(cvode_mem, indexB2, reltol, abstolB);
  if(check_flag(&flag, "CVodeSStolerancesB", 1)) return(1);

  flag = CVodeSetUserDataB(cvode_mem, indexB2, data);
  if(check_flag(&flag, "CVodeSetUserDataB", 1)) return(1);

  flag = CVodeQuadInitBS(cvode_mem, indexB2, fQB2, yQB2);
  if(check_flag(&flag, "CVodeQuadInitBS", 1)) return(1);

  flag = CVodeQuadSStolerancesB(cvode_mem, indexB2, reltol, abstolQB);
  if(check_flag(&flag, "CVodeQuadSStolerancesB", 1)) return(1);

  flag = CVodeSetQuadErrConB(cvode_mem, indexB2, SUNTRUE);
  if(check_flag(&flag, "CVodeSetQuadErrConB", 1)) return(1);

  /* Create a dense SUNMatrix */
  AB2 = SUNDenseMatrix(2*Neq, 2*Neq);
  if(check_flag((void *)AB2, "SUNDenseMatrix", 0)) return(1);

  /* Create banded SUNLinearSolver for the forward problem */
  LSB2 = SUNDenseLinearSolver(yB2, AB2);
  if(check_flag((void *)LSB2, "SUNDenseLinearSolver", 0)) return(1);

  /* Attach the matrix and linear solver */
  flag = CVDlsSetLinearSolverB(cvode_mem, indexB2, LSB2, AB2);
  if(check_flag(&flag, "CVDlsSetLinearSolverB", 1)) return(1);

  /* Backward integration */

  printf("---------------------------------------------\n");
  printf("Backward integration ... (2 adjoint problems)\n");
  printf("---------------------------------------------\n\n");

  flag = CVodeB(cvode_mem, t0, CV_NORMAL);
  if(check_flag(&flag, "CVodeB", 1)) return(1);

  flag = CVodeGetB(cvode_mem, indexB1, &time, yB1);
  if(check_flag(&flag, "CVodeGetB", 1)) return(1);

  flag = CVodeGetQuadB(cvode_mem, indexB1, &time, yQB1);
  if(check_flag(&flag, "CVodeGetQuadB", 1)) return(1);

  flag = CVodeGetB(cvode_mem, indexB2, &time, yB2);
  if(check_flag(&flag, "CVodeGetB", 1)) return(1);

  flag = CVodeGetQuadB(cvode_mem, indexB2, &time, yQB2);
  if(check_flag(&flag, "CVodeGetQuadB", 1)) return(1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le   (from backward pb. 1)\n", -Ith(yQB1,1), -Ith(yQB1,2));
  printf("           %12.4Le %12.4Le   (from backward pb. 2)\n", -Ith(yQB2,1), -Ith(yQB2,2));
#else
  printf("   dG/dp:  %12.4e %12.4e   (from backward pb. 1)\n", -Ith(yQB1,1), -Ith(yQB1,2));
  printf("           %12.4e %12.4e   (from backward pb. 2)\n", -Ith(yQB2,1), -Ith(yQB2,2));
#endif  
  printf("\n");
  printf("   H = d2G/dp2:\n");
  printf("        (1)            (2)\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("  %12.4Le   %12.4Le\n", -Ith(yQB1,3) , -Ith(yQB2,3));
  printf("  %12.4Le   %12.4Le\n", -Ith(yQB1,4) , -Ith(yQB2,4));
#else
  printf("  %12.4e   %12.4e\n", -Ith(yQB1,3) , -Ith(yQB2,3));
  printf("  %12.4e   %12.4e\n", -Ith(yQB1,4) , -Ith(yQB2,4));
#endif
  printf("\n");

  printf("Final Statistics for backward pb. 1\n");
  printf("-----------------------------------\n");
  PrintBckStats(cvode_mem, indexB1);

  printf("Final Statistics for backward pb. 2\n");
  printf("-----------------------------------\n");
  PrintBckStats(cvode_mem, indexB2);

  /* Free memory */

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNLinSolFree(LSB1);
  SUNMatDestroy(AB1);
  SUNLinSolFree(LSB2);
  SUNMatDestroy(AB2);

  /* Finite difference tests */

  dp = RCONST(1.0e-2);

  printf("-----------------------\n");
  printf("Finite Difference tests\n");
  printf("-----------------------\n\n");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("del_p = %Lg\n\n",dp);
#else
  printf("del_p = %g\n\n",dp);
#endif

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);

  flag = CVodeInit(cvode_mem, f, t0, y);
  if(check_flag(&flag, "CVodeInit", 1)) return(1);

  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  flag = CVodeSetUserData(cvode_mem, data);
  if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Create a dense SUNMatrix */
  A = SUNDenseMatrix(Neq, Neq);
  if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create banded SUNLinearSolver for the forward problem */
  LS = SUNDenseLinearSolver(y, A);
  if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

  /* Attach the matrix and linear solver */
  flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
  if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

  flag = CVodeQuadInit(cvode_mem, fQ, yQ);
  if(check_flag(&flag, "CVodeQuadInit", 1)) return(1);

  flag = CVodeQuadSStolerances(cvode_mem, reltol, abstolQ);
  if(check_flag(&flag, "CVodeQuadSStolerances", 1)) return(1);

  flag = CVodeSetQuadErrCon(cvode_mem, SUNTRUE);
  if(check_flag(&flag, "CVodeSetQuadErrCon", 1)) return(1);

  data->p1 += dp;

  flag = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  if(check_flag(&flag, "CVode", 1)) return(1);

  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  if(check_flag(&flag, "CVodeGetQuad", 1)) return(1);

  Gp = Ith(yQ,1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("p1+  y:   %12.4Le %12.4Le %12.4Le", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4Le\n",Ith(yQ,1));
#else
  printf("p1+  y:   %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4e\n",Ith(yQ,1));
#endif  
  data->p1 -= 2.0*dp;

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);

  CVodeReInit(cvode_mem, t0, y);
  CVodeQuadReInit(cvode_mem, yQ);

  flag = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  if(check_flag(&flag, "CVode", 1)) return(1);

  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  if(check_flag(&flag, "CVodeGetQuad", 1)) return(1);

  Gm = Ith(yQ,1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("p1-  y:   %12.4Le %12.4Le %12.4Le", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4Le\n",Ith(yQ,1));
#else
  printf("p1-  y:   %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4e\n",Ith(yQ,1));
#endif
  data->p1 += dp;

  grdG_fwd[0] = (Gp-G)/dp;
  grdG_bck[0] = (G-Gm)/dp;
  grdG_cntr[0] = (Gp-Gm)/(2.0*dp);
  H11 = (Gp - 2.0*G + Gm) / (dp*dp);

  data->p2 += dp;

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);

  CVodeReInit(cvode_mem, t0, y);
  CVodeQuadReInit(cvode_mem, yQ);

  flag = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  if(check_flag(&flag, "CVode", 1)) return(1);

  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  if(check_flag(&flag, "CVodeGetQuad", 1)) return(1);

  Gp = Ith(yQ,1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("p2+  y:   %12.4Le %12.4Le %12.4Le", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4Le\n",Ith(yQ,1));
#else
  printf("p2+  y:   %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4e\n",Ith(yQ,1));
#endif
  data->p2 -= 2.0*dp;

  N_VConst(ONE, y);
  N_VConst(ZERO, yQ);

  CVodeReInit(cvode_mem, t0, y);
  CVodeQuadReInit(cvode_mem, yQ);

  flag = CVode(cvode_mem, tf, y, &time, CV_NORMAL);
  if(check_flag(&flag, "CVode", 1)) return(1);

  flag = CVodeGetQuad(cvode_mem, &time, yQ);
  if(check_flag(&flag, "CVodeGetQuad", 1)) return(1);

  Gm = Ith(yQ,1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("p2-  y:   %12.4Le %12.4Le %12.4Le", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4Le\n",Ith(yQ,1));
#else
  printf("p2-  y:   %12.4e %12.4e %12.4e", Ith(y,1), Ith(y,2), Ith(y,3));
  printf("     G:   %12.4e\n",Ith(yQ,1));
#endif  
  data->p2 += dp;

  grdG_fwd[1] = (Gp-G)/dp;
  grdG_bck[1] = (G-Gm)/dp;
  grdG_cntr[1] = (Gp-Gm)/(2.0*dp);
  H22 = (Gp - 2.0*G + Gm) / (dp*dp);

  printf("\n");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("   dG/dp:  %12.4Le %12.4Le   (fwd FD)\n", grdG_fwd[0], grdG_fwd[1]);
  printf("           %12.4Le %12.4Le   (bck FD)\n", grdG_bck[0], grdG_bck[1]);
  printf("           %12.4Le %12.4Le   (cntr FD)\n", grdG_cntr[0], grdG_cntr[1]);
  printf("\n");
  printf("  H(1,1):  %12.4Le\n", H11);
  printf("  H(2,2):  %12.4Le\n", H22);
#else
  printf("   dG/dp:  %12.4e %12.4e   (fwd FD)\n", grdG_fwd[0], grdG_fwd[1]);
  printf("           %12.4e %12.4e   (bck FD)\n", grdG_bck[0], grdG_bck[1]);
  printf("           %12.4e %12.4e   (cntr FD)\n", grdG_cntr[0], grdG_cntr[1]);
  printf("\n");
  printf("  H(1,1):  %12.4e\n", H11);
  printf("  H(2,2):  %12.4e\n", H22);
#endif  

  /* Free memory */

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);

  N_VDestroy(y);
  N_VDestroy(yQ);

  N_VDestroyVectorArray(yS, Np);
  N_VDestroyVectorArray(yQS, Np);

  N_VDestroy(yB1);
  N_VDestroy(yQB1);
  N_VDestroy(yB2);
  N_VDestroy(yQB2);

  free(data);

  return(0);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */


static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, y3;
  UserData data;
  realtype p1, p2;

  data = (UserData) user_data;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  Ith(ydot,1) = -p1*y1*y1 - y3;
  Ith(ydot,2) = -y2;
  Ith(ydot,3) = -p2*p2*y2*y3;

  return(0);
}

static int fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data)
{
  realtype y1, y2, y3;

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  Ith(qdot,1) = 0.5 * ( y1*y1 + y2*y2 + y3*y3 );

  return(0);
}

static int fS(int Ns, realtype t,
              N_Vector y, N_Vector ydot,
              N_Vector *yS, N_Vector *ySdot,
              void *user_data,
              N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  realtype y1, y2, y3;
  realtype s1, s2, s3;
  realtype fys1, fys2, fys3;
  realtype p1, p2;

  data = (UserData) user_data;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  /* 1st sensitivity RHS */

  s1 = Ith(yS[0],1);
  s2 = Ith(yS[0],2);
  s3 = Ith(yS[0],3);

  fys1 = - 2.0*p1*y1 * s1 - s3;
  fys2 = - s2;
  fys3 = - p2*p2*y3 * s2 - p2*p2*y2 * s3;

  Ith(ySdot[0],1) = fys1 - y1*y1;
  Ith(ySdot[0],2) = fys2;
  Ith(ySdot[0],3) = fys3;

  /* 2nd sensitivity RHS */

  s1 = Ith(yS[1],1);
  s2 = Ith(yS[1],2);
  s3 = Ith(yS[1],3);

  fys1 = - 2.0*p1*y1 * s1 - s3;
  fys2 = - s2;
  fys3 = - p2*p2*y3 * s2 - p2*p2*y2 * s3;

  Ith(ySdot[1],1) = fys1;
  Ith(ySdot[1],2) = fys2;
  Ith(ySdot[1],3) = fys3 - 2.0*p2*y2*y3;

  return(0);
}

static int fQS(int Ns, realtype t,
               N_Vector y, N_Vector *yS, 
               N_Vector yQdot, N_Vector *yQSdot,
               void *user_data,
               N_Vector tmp, N_Vector tmpQ)
{
  realtype y1, y2, y3;
  realtype s1, s2, s3;

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);


  /* 1st sensitivity RHS */

  s1 = Ith(yS[0],1);
  s2 = Ith(yS[0],2);
  s3 = Ith(yS[0],3);

  Ith(yQSdot[0],1) = y1*s1 + y2*s2 + y3*s3;


  /* 1st sensitivity RHS */

  s1 = Ith(yS[1],1);
  s2 = Ith(yS[1],2);
  s3 = Ith(yS[1],3);

  Ith(yQSdot[1],1) = y1*s1 + y2*s2 + y3*s3;

  return(0);
}

static int fB1(realtype t, N_Vector y, N_Vector *yS, 
               N_Vector yB, N_Vector yBdot, void *user_dataB)
{
  UserData data;
  realtype p1, p2;
  realtype y1, y2, y3;  /* solution */
  realtype s1, s2, s3;  /* sensitivity 1 */
  realtype l1, l2, l3;  /* lambda */
  realtype m1, m2, m3;  /* mu */
  
  data = (UserData) user_dataB;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);
  
  s1 = Ith(yS[0],1); 
  s2 = Ith(yS[0],2); 
  s3 = Ith(yS[0],3);

  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 
  l3 = Ith(yB,3);

  m1 = Ith(yB,4); 
  m2 = Ith(yB,5); 
  m3 = Ith(yB,6);

  
  Ith(yBdot,1) = 2.0*p1*y1 * l1     - y1;
  Ith(yBdot,2) = l2 + p2*p2*y3 * l3 - y2;
  Ith(yBdot,3) = l1 + p2*p2*y2 * l3 - y3;

  Ith(yBdot,4) = 2.0*p1*y1 * m1     + l1 * 2.0*(y1 + p1*s1) - s1;
  Ith(yBdot,5) = m2 + p2*p2*y3 * m3 + l3 * p2*p2*s3         - s2;
  Ith(yBdot,6) = m1 + p2*p2*y2 * m3 + l3 * p2*p2*s2         - s3;

  return(0);
}

static int fQB1(realtype t, N_Vector y, N_Vector *yS,
                N_Vector yB, N_Vector qBdot, void *user_dataB)
{
  UserData data;
  realtype p1, p2;
  realtype y1, y2, y3;  /* solution */
  realtype s1, s2, s3;  /* sensitivity 1 */
  realtype l1, l2, l3;  /* lambda */
  realtype m1, m2, m3;  /* mu */
  
  data = (UserData) user_dataB;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);
  
  s1 = Ith(yS[0],1); 
  s2 = Ith(yS[0],2); 
  s3 = Ith(yS[0],3);
  
  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 
  l3 = Ith(yB,3);

  m1 = Ith(yB,4); 
  m2 = Ith(yB,5); 
  m3 = Ith(yB,6);

  Ith(qBdot,1) = -y1*y1 * l1;
  Ith(qBdot,2) = -2.0*p2*y2*y3 * l3;

  Ith(qBdot,3) = -y1*y1 * m1        - l1 * 2.0*y1*s1;
  Ith(qBdot,4) = -2.0*p2*y2*y3 * m3 - l3 * 2.0*(p2*y3*s2 + p2*y2*s3);

  return(0);
}




static int fB2(realtype t, N_Vector y, N_Vector *yS, 
               N_Vector yB, N_Vector yBdot, void *user_dataB)
{
  UserData data;
  realtype p1, p2;
  realtype y1, y2, y3;  /* solution */
  realtype s1, s2, s3;  /* sensitivity 2 */
  realtype l1, l2, l3;  /* lambda */
  realtype m1, m2, m3;  /* mu */

  data = (UserData) user_dataB;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);
  
  s1 = Ith(yS[1],1); 
  s2 = Ith(yS[1],2); 
  s3 = Ith(yS[1],3);
  
  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 
  l3 = Ith(yB,3);

  m1 = Ith(yB,4); 
  m2 = Ith(yB,5); 
  m3 = Ith(yB,6);

  Ith(yBdot,1) = 2.0*p1*y1 * l1     - y1;
  Ith(yBdot,2) = l2 + p2*p2*y3 * l3 - y2;
  Ith(yBdot,3) = l1 + p2*p2*y2 * l3 - y3;

  Ith(yBdot,4) = 2.0*p1*y1 * m1     + l1 * 2.0*p1*s1              - s1;
  Ith(yBdot,5) = m2 + p2*p2*y3 * m3 + l3 * (2.0*p2*y3 + p2*p2*s3) - s2;
  Ith(yBdot,6) = m1 + p2*p2*y2 * m3 + l3 * (2.0*p2*y2 + p2*p2*s2) - s3;


  return(0);
}


static int fQB2(realtype t, N_Vector y, N_Vector *yS,
                N_Vector yB, N_Vector qBdot, void *user_dataB)
{
  UserData data;
  realtype p1, p2;
  realtype y1, y2, y3;  /* solution */
  realtype s1, s2, s3;  /* sensitivity 2 */
  realtype l1, l2, l3;  /* lambda */
  realtype m1, m2, m3;  /* mu */
  
  data = (UserData) user_dataB;
  p1 = data->p1; 
  p2 = data->p2; 

  y1 = Ith(y,1); 
  y2 = Ith(y,2); 
  y3 = Ith(y,3);

  s1 = Ith(yS[1],1); 
  s2 = Ith(yS[1],2); 
  s3 = Ith(yS[1],3);

  l1 = Ith(yB,1); 
  l2 = Ith(yB,2); 
  l3 = Ith(yB,3);

  m1 = Ith(yB,4); 
  m2 = Ith(yB,5); 
  m3 = Ith(yB,6);

  Ith(qBdot,1) = -y1*y1 * l1;
  Ith(qBdot,2) = -2.0*p2*y2*y3 * l3;

  Ith(qBdot,3) = -y1*y1 * m1        - l1 * 2.0*y1*s1;
  Ith(qBdot,4) = -2.0*p2*y2*y3 * m3 - l3 * 2.0*(p2*y3*s2 + p2*y2*s3 + y2*y3);

  return(0);
}


/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

void PrintFwdStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nfQe, netfQ;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int nfQSe, netfQS;

  int qlast, qcur;
  realtype h0u, hlast, hcur, tcur;

  int flag;


  flag = CVodeGetIntegratorStats(cvode_mem, &nst, &nfe, &nsetups, &netf, 
                                 &qlast, &qcur,
                                 &h0u, &hlast, &hcur,
                                 &tcur);

  flag = CVodeGetNonlinSolvStats(cvode_mem, &nni, &ncfn);

  flag = CVodeGetQuadStats(cvode_mem, &nfQe, &netfQ);

  flag = CVodeGetSensStats(cvode_mem, &nfSe, &nfeS, &netfS, &nsetupsS);

  flag = CVodeGetSensNonlinSolvStats(cvode_mem, &nniS, &ncfnS);

  flag = CVodeGetQuadSensStats(cvode_mem, &nfQSe, &netfQS);


  printf(" Number steps: %5ld\n\n", nst);
  printf(" Function evaluations:\n");
  printf("  f:        %5ld\n  fQ:       %5ld\n  fS:       %5ld\n  fQS:      %5ld\n",
         nfe, nfQe, nfSe, nfQSe);
  printf(" Error test failures:\n");
  printf("  netf:     %5ld\n  netfQ:    %5ld\n  netfS:    %5ld\n  netfQS:   %5ld\n",
         netf, netfQ, netfS, netfQS);
  printf(" Linear solver setups:\n");
  printf("  nsetups:  %5ld\n  nsetupsS: %5ld\n", nsetups, nsetupsS);
  printf(" Nonlinear iterations:\n");
  printf("  nni:      %5ld\n  nniS:     %5ld\n", nni, nniS);
  printf(" Convergence failures:\n");
  printf("  ncfn:     %5ld\n  ncfnS:    %5ld\n", ncfn, ncfnS);

  printf("\n");

}


void PrintBckStats(void *cvode_mem, int idx)
{
  void *cvode_mem_bck;

  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nfQe, netfQ;

  int qlast, qcur;
  realtype h0u, hlast, hcur, tcur;

  int flag;

  cvode_mem_bck = CVodeGetAdjCVodeBmem(cvode_mem, idx);

  flag = CVodeGetIntegratorStats(cvode_mem_bck, &nst, &nfe, &nsetups, &netf, 
                                 &qlast, &qcur,
                                 &h0u, &hlast, &hcur,
                                 &tcur);

  flag = CVodeGetNonlinSolvStats(cvode_mem_bck, &nni, &ncfn);

  flag = CVodeGetQuadStats(cvode_mem_bck, &nfQe, &netfQ);

  printf(" Number steps: %5ld\n\n", nst);
  printf(" Function evaluations:\n");
  printf("  f:        %5ld\n  fQ:       %5ld\n", nfe, nfQe);
  printf(" Error test failures:\n");
  printf("  netf:     %5ld\n  netfQ:    %5ld\n", netf, netfQ);
  printf(" Linear solver setups:\n");
  printf("  nsetups:  %5ld\n", nsetups);
  printf(" Nonlinear iterations:\n");
  printf("  nni:      %5ld\n", nni);
  printf(" Convergence failures:\n");
  printf("  ncfn:     %5ld\n", ncfn);

  printf("\n");


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

