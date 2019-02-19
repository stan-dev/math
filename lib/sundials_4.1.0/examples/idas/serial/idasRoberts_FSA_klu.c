/* -----------------------------------------------------------------
 * Programmer(s): Ting Yan @ SMU
 *      Based on idasRoberts_FSA_dns.c and modified to use KLU
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
 * Example problem:
 *
 * This simple example problem for IDA, due to Robertson, 
 * is from chemical kinetics, and consists of the following three 
 * equations:
 *
 *      dy1/dt = -p1*y1 + p2*y2*y3
 *      dy2/dt = p1*y1 - p2*y2*y3 - p3*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7
 *
 * Optionally, IDAS can compute sensitivities with respect to the
 * problem parameters p1, p2, and p3.
 * The sensitivity right hand side is given analytically through the
 * user routine fS (of type SensRhs1Fn).
 * Any of two sensitivity methods (SIMULTANEOUS and STAGGERED can be
 * used and sensitivities may be included in the error test or not 
 *(error control set on SUNTRUE or SUNFALSE, respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % idasRoberts_FSA_klu -nosensi
 * If sensitivities are to be computed:
 *    % idasRoberts_FSA_klu -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Header files with a description of contents used */

#include <idas/idas.h>                     /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>        /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_sparse.h>    /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_klu.h>       /* access to KLU linear solver          */
#include <sundials/sundials_types.h>       /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>        /* defs. of SUNRabs, SUNRexp, etc.      */

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i=1..NEQ */

/* Problem Constants */

#define NEQ   3             /* number of equations  */
#define T0    RCONST(0.0)   /* initial time */
#define T1    RCONST(0.4)   /* first output time */
#define TMULT RCONST(10.0)  /* output time factor */
#define NOUT  12            /* number of output times */

#define NP    3             /* number of problem parameters */
#define NS    3             /* number of sensitivities computed */

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)

/* Type : UserData */

typedef struct {
  realtype p[3];           /* problem parameters */
  realtype coef;
} *UserData;

/* Prototypes of functions by IDAS */

static int res(realtype t, N_Vector y, N_Vector yp, N_Vector resval, void *user_data);

static int Jac(realtype t, realtype cj, 
               N_Vector yy, N_Vector yp, N_Vector resvec, 
               SUNMatrix JJ, void *user_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int resS(int Ns, realtype t, 
                N_Vector y, N_Vector yp, N_Vector resval,
                N_Vector *yyS, N_Vector *ypS, N_Vector *resvalS,
                void *user_data, 
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int rhsQ(realtype tres, N_Vector yy, N_Vector yp,
                 N_Vector rrQ, void *user_data);

/* Prototypes of private functions */

static void ProcessArgs(int argc, char *argv[],
                        booleantype *sensi, int *sensi_meth, 
                        booleantype *err_con);
static void WrongArgs(char *name);

static void PrintIC(N_Vector y, N_Vector yp);
static void PrintSensIC(N_Vector y, N_Vector yp, N_Vector* yS, N_Vector* ypS); 

static void PrintOutput(void *ida_mem, realtype t, N_Vector u);
static void PrintSensOutput(N_Vector *uS);

static void PrintFinalStats(void *ida_mem, booleantype sensi);

static int check_retval(void *returnvalue, char *funcname, int opt);
/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  void *ida_mem;
  SUNMatrix A;
  SUNLinearSolver LS;
  UserData data;
  realtype reltol, t, tout;
  N_Vector y, yp, abstol, id;
  int iout, retval, nnz;

  realtype pbar[NS];
  int is; 
  N_Vector *yS, *ypS;
  booleantype sensi, err_con;
  int sensi_meth;

  N_Vector yQ, *yQS;

  ida_mem = NULL;
  data    = NULL;
  y       = NULL;
  yS      = NULL;
  ypS     = NULL;
  A       = NULL;
  LS      = NULL;

  /* Process arguments */
  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);

  /* User data structure */
  data = (UserData) malloc(sizeof *data);
  if (check_retval((void *)data, "malloc", 2)) return(1);
  data->p[0] = RCONST(0.040);
  data->p[1] = RCONST(1.0e4);
  data->p[2] = RCONST(3.0e7);
  data->coef = 0.5;

  /* Initial conditions */
  y = N_VNew_Serial(NEQ);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) return(1);

  Ith(y,1) = ONE;
  Ith(y,2) = ZERO;
  Ith(y,3) = ZERO;

  yp = N_VNew_Serial(NEQ);
  if(check_retval((void *)yp, "N_VNew_Serial", 0)) return(1);

  /* These initial conditions are NOT consistent. See IDACalcIC below. */
  Ith(yp,1) = RCONST(0.1);
  Ith(yp,2) = ZERO;
  Ith(yp,3) = ZERO;  

  /* Create IDAS object */
  ida_mem = IDACreate();
  if (check_retval((void *)ida_mem, "IDACreate", 0)) return(1);

  /* Allocate space for IDAS */
  retval = IDAInit(ida_mem, res, T0, y, yp);
  if (check_retval(&retval, "IDAInit", 1)) return(1);

  /* Specify scalar relative tol. and vector absolute tol. */
  reltol = RCONST(1.0e-6);
  abstol = N_VNew_Serial(NEQ);
  Ith(abstol,1) = RCONST(1.0e-8);
  Ith(abstol,2) = RCONST(1.0e-14);
  Ith(abstol,3) = RCONST(1.0e-6);
  retval = IDASVtolerances(ida_mem, reltol, abstol);
  if (check_retval(&retval, "IDASVtolerances", 1)) return(1);

  /* Set ID vector */
  id = N_VNew_Serial(NEQ);
  Ith(id,1) = 1.0;
  Ith(id,2) = 1.0;
  Ith(id,3) = 0.0;
  retval = IDASetId(ida_mem, id);
  if (check_retval(&retval, "IDASetId", 1)) return(1);

  /* Attach user data */
  retval = IDASetUserData(ida_mem, data);
  if (check_retval(&retval, "IDASetUserData", 1)) return(1);

  /* Create sparse SUNMatrix for use in linear solves */
  nnz = NEQ * NEQ;
  A = SUNSparseMatrix(NEQ, NEQ, nnz, CSC_MAT);
  if(check_retval((void *)A, "SUNSparseMatrix", 0)) return(1);

  /* Create KLU SUNLinearSolver object (one thread) */
  LS = SUNLinSol_KLU(y, A);
  if(check_retval((void *)LS, "SUNLinSol_KLU", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(ida_mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  /* Set the user-supplied Jacobian routine */
  retval = IDASetJacFn(ida_mem, Jac);
  if(check_retval(&retval, "IDASetJacFn", 1)) return(1);

  printf("\n3-species chemical kinetics problem\n");

  /* Sensitivity-related settings */
  if (sensi) {

    pbar[0] = data->p[0];
    pbar[1] = data->p[1];
    pbar[2] = data->p[2];

    yS = N_VCloneVectorArray(NS, y);
    if (check_retval((void *)yS, "N_VCloneVectorArray", 0)) return(1);
    for (is=0;is<NS;is++) N_VConst(ZERO, yS[is]);
    
    ypS = N_VCloneVectorArray(NS, y);
    if (check_retval((void *)ypS, "N_VCloneVectorArray", 0)) return(1);
    for (is=0;is<NS;is++) N_VConst(ZERO, ypS[is]);

    /* 
    * Only non-zero sensitivity I.C. are ypS[0]: 
    * - Ith(ypS[0],1) = -ONE;
    * - Ith(ypS[0],2) =  ONE;
    *
    * They are not set. IDACalcIC also computes consistent IC for sensitivities.
    */
    
    retval = IDASensInit(ida_mem, NS, sensi_meth, resS, yS, ypS);
    if(check_retval(&retval, "IDASensInit", 1)) return(1);

    retval = IDASensEEtolerances(ida_mem);
    if(check_retval(&retval, "IDASensEEtolerances", 1)) return(1);

    retval = IDASetSensErrCon(ida_mem, err_con);
    if (check_retval(&retval, "IDASetSensErrCon", 1)) return(1);

    retval = IDASetSensParams(ida_mem, data->p, pbar, NULL);
    if (check_retval(&retval, "IDASetSensParams", 1)) return(1);

    printf("Sensitivity: YES ");
    if(sensi_meth == IDA_SIMULTANEOUS)   
      printf("( SIMULTANEOUS +");
    else 
      printf("( STAGGERED +");   
    if(err_con) printf(" FULL ERROR CONTROL )");
    else        printf(" PARTIAL ERROR CONTROL )");

  } else {

    printf("Sensitivity: NO ");

  }

  /*----------------------------------------------------------
   *               Q U A D R A T U R E S
   * ---------------------------------------------------------*/
  yQ = N_VNew_Serial(2);

  Ith(yQ,1) = 0;
  Ith(yQ,2) = 0;

  IDAQuadInit(ida_mem, rhsQ, yQ);

  yQS = N_VCloneVectorArray(NS, yQ);
  for (is=0;is<NS;is++) N_VConst(ZERO, yQS[is]);  

  IDAQuadSensInit(ida_mem, NULL, yQS);

  /* Call IDACalcIC to compute consistent initial conditions. If sensitivity is
     enabled, this function also try to find consistent IC for the sensitivities. */

  retval = IDACalcIC(ida_mem, IDA_YA_YDP_INIT, T1);;
  if (check_retval(&retval, "IDACalcIC", 1)) return(1);

  retval = IDAGetConsistentIC(ida_mem, y, yp);
  if (check_retval(&retval, "IDAGetConsistentIC", 1)) return(1);

  PrintIC(y, yp);

  if(sensi) {
      IDAGetSensConsistentIC(ida_mem, yS, ypS);
      PrintSensIC(y, yp, yS, ypS); 
    }
      
  /* In loop over output points, call IDA, print results, test for error */
  
  printf("\n\n");
  printf("===========================================");
  printf("============================\n");
  printf("     T     Q       H      NST           y1");
  printf("           y2           y3    \n");
  printf("===========================================");
  printf("============================\n");

  for (iout=1, tout=T1; iout <= NOUT; iout++, tout *= TMULT) {
    
    retval = IDASolve(ida_mem, tout, &t, y, yp, IDA_NORMAL);
    if (check_retval(&retval, "IDASolve", 1)) break;

    PrintOutput(ida_mem, t, y);
    
    if (sensi) {
      retval = IDAGetSens(ida_mem, &t, yS);
      if (check_retval(&retval, "IDAGetSens", 1)) break;
      PrintSensOutput(yS);
    } 
    printf("-----------------------------------------");
    printf("------------------------------\n");

  }

  printf("\nQuadrature:\n");
  IDAGetQuad(ida_mem, &t, yQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("G:      %10.4Le\n", Ith(yQ,1));
#else
  printf("G:      %10.4e\n", Ith(yQ,1));
#endif  

  if(sensi) {
    IDAGetQuadSens(ida_mem, &t, yQS);
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("\nSensitivities at t=%Lg:\n",t);
    printf("dG/dp1: %11.4Le\n", Ith(yQS[0], 1));
    printf("dG/dp1: %11.4Le\n", Ith(yQS[1], 1));
    printf("dG/dp1: %11.4Le\n", Ith(yQS[2], 1));
#else
    printf("\nSensitivities at t=%g:\n",t);
    printf("dG/dp1: %11.4e\n", Ith(yQS[0], 1));
    printf("dG/dp1: %11.4e\n", Ith(yQS[1], 1));
    printf("dG/dp1: %11.4e\n", Ith(yQS[2], 1));
#endif    
  }

  /* Print final statistics */
  PrintFinalStats(ida_mem, sensi);

  /* Free memory */
  N_VDestroy(y);
  if (sensi) {
    N_VDestroyVectorArray(yS, NS);
  }
  free(data);
  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(yQ);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDAS
 *--------------------------------------------------------------------
 */

/*
 * Residual routine. Compute F(t,y,y',p). 
 */
static int res(realtype t, N_Vector yy, N_Vector yp, N_Vector resval, void *user_data)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;
  realtype yp1, yp2;

  data = (UserData) user_data;
  p1 = data->p[0];
  p2 = data->p[1];
  p3 = data->p[2];

  y1 = Ith(yy,1);
  y2 = Ith(yy,2);
  y3 = Ith(yy,3);

  yp1 = Ith(yp,1);
  yp2 = Ith(yp,2);

  Ith(resval,1) = yp1 + p1*y1 - p2*y2*y3;
  Ith(resval,2) = yp2 - p1*y1 + p2*y2*y3 + p3*y2*y2;
  Ith(resval,3) = y1 + y2 + y3 - ONE;

  return(0);
}


/* 
 * Jacobian routine. Compute J(t,y). 
*/

static int Jac(realtype t, realtype cj,
               N_Vector yy, N_Vector yp, N_Vector resvec, 
               SUNMatrix JJ, void *user_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *yval;
  sunindextype *colptrs = SUNSparseMatrix_IndexPointers(JJ);
  sunindextype *rowvals = SUNSparseMatrix_IndexValues(JJ);
  realtype *data = SUNSparseMatrix_Data(JJ);

  UserData userdata;
  realtype p1, p2, p3;
 
  yval = N_VGetArrayPointer(yy);

  userdata = (UserData) user_data;
  p1 = userdata->p[0]; p2 = userdata->p[1]; p3 = userdata->p[2];

  SUNMatZero(JJ);

  colptrs[0] = 0;
  colptrs[1] = 3;
  colptrs[2] = 6;
  colptrs[3] = 9;

  data[0] = p1+cj;
  rowvals[0] = 0;
  data[1] = -p1;
  rowvals[1] = 1;
  data[2] = ONE;
  rowvals[2] = 2;

  data[3] = -p2*yval[2];
  rowvals[3] = 0;
  data[4] = p2*yval[2]+2*p3*yval[1]+cj;
  rowvals[4] = 1;
  data[5] = ONE;
  rowvals[5] = 2;

  data[6] = -p2*yval[1];
  rowvals[6] = 0;
  data[7] = p2*yval[1];
  rowvals[7] = 1;
  data[8] = ONE;
  rowvals[8] = 2;

  return(0);
}
 
/* 
 * resS routine. Compute sensitivity r.h.s. 
 */

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

  for (is=0; is<NS; is++) {

    s1 = Ith(yyS[is],1);
    s2 = Ith(yyS[is],2);
    s3 = Ith(yyS[is],3);

    sd1 = Ith(ypS[is],1);
    sd2 = Ith(ypS[is],2);

    rs1 = sd1 + p1*s1 - p2*y3*s2 - p2*y2*s3;
    rs2 = sd2 - p1*s1 + p2*y3*s2 + p2*y2*s3 + 2*p3*y2*s2;
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
    case 2:
      rs2 += y2*y2;
      break;
    }
  
    Ith(resvalS[is],1) = rs1;
    Ith(resvalS[is],2) = rs2;
    Ith(resvalS[is],3) = rs3;

  }

  return(0);
}

static int rhsQ(realtype t, N_Vector y, N_Vector yp, 
              N_Vector ypQ, void* user_data)
{
  UserData data;

  data = (UserData) user_data;

  Ith(ypQ,1) = Ith(y,3);
  
  Ith(ypQ,2) = data->coef*( Ith(y,1)*Ith(y,1)+ 
                            Ith(y,2)*Ith(y,2)+ 
                            Ith(y,3)*Ith(y,3) );
                    
  return(0);
}


/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Process and verify arguments to idasfwddenx.
 */

static void ProcessArgs(int argc, char *argv[], 
                        booleantype *sensi, int *sensi_meth, booleantype *err_con)
{
  *sensi = SUNFALSE;
  *sensi_meth = -1;
  *err_con = SUNFALSE;

  if (argc < 2) WrongArgs(argv[0]);

  if (strcmp(argv[1],"-nosensi") == 0)
    *sensi = SUNFALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    *sensi = SUNTRUE;
  else
    WrongArgs(argv[0]);
  
  if (*sensi) {

    if (argc != 4)
      WrongArgs(argv[0]);

    if (strcmp(argv[2],"sim") == 0)
      *sensi_meth = IDA_SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      *sensi_meth = IDA_STAGGERED;
    else 
      WrongArgs(argv[0]);

    if (strcmp(argv[3],"t") == 0)
      *err_con = SUNTRUE;
    else if (strcmp(argv[3],"f") == 0)
      *err_con = SUNFALSE;
    else
      WrongArgs(argv[0]);
  }

}

static void WrongArgs(char *name)
{
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim or stg\n");
    printf("         err_con    = t or f\n");
    
    exit(0);
}


static void PrintIC(N_Vector y, N_Vector yp)
{
  realtype* data;

  data = N_VGetArrayPointer(y);
  printf("\n\nConsistent IC:\n");
  printf("\ty = ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", data[0], data[1], data[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", data[0], data[1], data[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", data[0], data[1], data[2]);
#endif

  data = N_VGetArrayPointer(yp);
  printf("\typ= ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", data[0], data[1], data[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", data[0], data[1], data[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", data[0], data[1], data[2]);
#endif

}
static void PrintSensIC(N_Vector y, N_Vector yp, N_Vector* yS, N_Vector* ypS)
{
  realtype *sdata;

  sdata = N_VGetArrayPointer(yS[0]);
  printf("                  Sensitivity 1  ");

  printf("\n\ts1 = ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
  sdata = N_VGetArrayPointer(ypS[0]);
  printf("\ts1'= ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif


  printf("                  Sensitivity 2  ");
  sdata = N_VGetArrayPointer(yS[1]);
  printf("\n\ts2 = ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
  sdata = N_VGetArrayPointer(ypS[1]);
  printf("\ts2'= ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif


  printf("                  Sensitivity 3  ");
  sdata = N_VGetArrayPointer(yS[2]);
  printf("\n\ts3 = ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
  sdata = N_VGetArrayPointer(ypS[2]);
  printf("\ts3'= ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif


}

/*
 * Print current t, step count, order, stepsize, and solution.
 */

static void PrintOutput(void *ida_mem, realtype t, N_Vector u)
{
  long int nst;
  int qu, retval;
  realtype hu, *udata;
  
  udata = N_VGetArrayPointer(u);

  retval = IDAGetNumSteps(ida_mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastOrder(ida_mem, &qu);
  check_retval(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetLastStep(ida_mem, &hu);
  check_retval(&retval, "IDAGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.3Le %2d  %8.3Le %5ld\n", t, qu, hu, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#else
  printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#endif

  printf("                  Solution       ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", udata[0], udata[1], udata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#endif

}

/* 
 * Print sensitivities.
*/

static void PrintSensOutput(N_Vector *uS)
{
  realtype *sdata;

  sdata = N_VGetArrayPointer(uS[0]);
  printf("                  Sensitivity 1  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
  
  sdata = N_VGetArrayPointer(uS[1]);
  printf("                  Sensitivity 2  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif

  sdata = N_VGetArrayPointer(uS[2]);
  printf("                  Sensitivity 3  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
  printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
}

/* 
 * Print some final statistics from the IDAS memory.
 */

static void PrintFinalStats(void *ida_mem, booleantype sensi)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  int retval;

  retval = IDAGetNumSteps(ida_mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(ida_mem, &nfe);
  check_retval(&retval, "IDAGetNumRhsEvals", 1);
  retval = IDAGetNumLinSolvSetups(ida_mem, &nsetups);
  check_retval(&retval, "IDAGetNumLinSolvSetups", 1);
  retval = IDAGetNumErrTestFails(ida_mem, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvIters(ida_mem, &nni);
  check_retval(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1);

  if (sensi) {
    retval = IDAGetSensNumResEvals(ida_mem, &nfSe);
    check_retval(&retval, "IDAGetSensNumRhsEvals", 1);
    retval = IDAGetNumResEvalsSens(ida_mem, &nfeS);
    check_retval(&retval, "IDAGetNumResEvalsSens", 1);
    retval = IDAGetSensNumLinSolvSetups(ida_mem, &nsetupsS);
    check_retval(&retval, "IDAGetSensNumLinSolvSetups", 1);
    retval = IDAGetSensNumErrTestFails(ida_mem, &netfS);
    check_retval(&retval, "IDAGetSensNumErrTestFails", 1);
    retval = IDAGetSensNumNonlinSolvIters(ida_mem, &nniS);
    check_retval(&retval, "IDAGetSensNumNonlinSolvIters", 1);
    retval = IDAGetSensNumNonlinSolvConvFails(ida_mem, &ncfnS);
    check_retval(&retval, "IDAGetSensNumNonlinSolvConvFails", 1);
  }

  printf("\nFinal Statistics\n\n");
  printf("nst     = %5ld\n\n", nst);
  printf("nfe     = %5ld\n",   nfe);
  printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  if(sensi) {
    printf("\n");
    printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }

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

static int check_retval(void *returnvalue, char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
