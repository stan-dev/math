/* -----------------------------------------------------------------
 * Programmer(s): Chris Nguyen @ LLNL
 *                based on idaHeat2D_bnd.c and idaRoberts_klu.c
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
 * Example problem for IDA: 2D heat equation, serial, sparse.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the KLU solver and IDACalcIC.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform MGRID x MGRID
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = MGRID^2. Here MGRID = 10.
 *
 * The system is solved with IDA using the sparse linear system
 * solver and a user supplied Jacobian.
 * For purposes of illustration,
 * IDACalcIC is called to compute correct values at the boundary,
 * given incorrect values as input initial guesses. The constraints
 * u >= 0 are posed for all components. Output is taken at
 * t = 0, .01, .02, .04, ..., 10.24. (Output at t = 0 is for
 * IDACalcIC cost statistics only.)
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>                       /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>        /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_sparse.h>    /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_klu.h>       /* access to KLU linear solver          */
#include <sundials/sundials_types.h>       /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>        /* defs. of SUNRabs, SUNRexp, etc.      */

/* Problem Constants */

#define NOUT  11
#define MGRID 10
#define NEQ   MGRID*MGRID
#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)
#define BVAL  RCONST(0.0)
#define TOTAL 4*MGRID+8*(MGRID-2)+(MGRID-4)*(MGRID+4*(MGRID-2)) /* total num of nonzero elements */

/* Type: UserData */

typedef struct {
  sunindextype mm;
  realtype dx;
  realtype coeff;
} *UserData;

/* Prototypes of functions called by IDA */

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector resval, void *user_data);

int jacHeat(realtype tt,  realtype cj,
	    N_Vector yy, N_Vector yp, N_Vector resvec,
	    SUNMatrix JJ, void *user_data,
	    N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Exact same setup as jacHeat. Function needed for special case MGRID=3  */
int jacHeat3(realtype tt,  realtype cj,
             N_Vector yy, N_Vector yp, N_Vector resvec,
             SUNMatrix JJ, void *user_data,
             N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Prototypes of private functions */

static void PrintHeader(realtype rtol, realtype atol);
static void PrintOutput(void *mem, realtype t, N_Vector u);
static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up,
                             N_Vector id, N_Vector res);

static int check_retval(void *returnvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(void)
{
  void *mem;
  UserData data;
  N_Vector uu, up, constraints, id, res;
  int retval, iout;
  long int netf, ncfn;
  realtype rtol, atol, t0, t1, tout, tret;
  SUNMatrix A;
  SUNLinearSolver LS;
  sunindextype nnz;

  mem = NULL;
  data = NULL;
  uu = up = constraints = id = res = NULL;
  A = NULL;
  LS = NULL;

  /* Create vectors uu, up, res, constraints, id. */
  uu = N_VNew_Serial(NEQ);
  if(check_retval((void *)uu, "N_VNew_Serial", 0)) return(1);
  up = N_VNew_Serial(NEQ);
  if(check_retval((void *)up, "N_VNew_Serial", 0)) return(1);
  res = N_VNew_Serial(NEQ);
  if(check_retval((void *)res, "N_VNew_Serial", 0)) return(1);
  constraints = N_VNew_Serial(NEQ);
  if(check_retval((void *)constraints, "N_VNew_Serial", 0)) return(1);
  id = N_VNew_Serial(NEQ);
  if(check_retval((void *)id, "N_VNew_Serial", 0)) return(1);

  /* Create and load problem data block. */
  data = (UserData) malloc(sizeof *data);
  if(check_retval((void *)data, "malloc", 2)) return(1);
  data->mm = MGRID;
  data->dx = ONE/(MGRID - ONE);
  data->coeff = ONE/( (data->dx) * (data->dx) );

  /* Initialize uu, up, id. */
  SetInitialProfile(data, uu, up, id, res);

  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst(ONE, constraints);

  /* Set remaining input parameters. */
  t0   = ZERO;
  t1   = RCONST(0.01);
  rtol = ZERO;
  atol = RCONST(1.0e-8);

  /* Call IDACreate and IDAMalloc to initialize solution */
  mem = IDACreate();
  if(check_retval((void *)mem, "IDACreate", 0)) return(1);

  retval = IDASetUserData(mem, data);
  if(check_retval(&retval, "IDASetUserData", 1)) return(1);

  /* Set which components are algebraic or differential */
  retval = IDASetId(mem, id);
  if(check_retval(&retval, "IDASetId", 1)) return(1);

  retval = IDASetConstraints(mem, constraints);
  if(check_retval(&retval, "IDASetConstraints", 1)) return(1);
  N_VDestroy(constraints);

  retval = IDAInit(mem, heatres, t0, uu, up);
  if(check_retval(&retval, "IDAInit", 1)) return(1);

  retval = IDASStolerances(mem, rtol, atol);
  if(check_retval(&retval, "IDASStolerances", 1)) return(1);

  /* Create sparse SUNMatrix for use in linear solves */
  nnz = NEQ*NEQ;
  A = SUNSparseMatrix(NEQ, NEQ, nnz, CSC_MAT);
  if(check_retval((void*)A, "SUNSparseMtarix", 0)) return(1);

  /* Create KLU SUNLinearSolver object */
  LS = SUNLinSol_KLU(uu, A);
  if(check_retval((void *)LS, "SUNLinSol_KLU", 0)) return(1);

  /* Attach the matrix and linear solver */
  retval = IDASetLinearSolver(mem, LS, A);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  /* Set the user-supplied Jacobian routine */
  if(MGRID >= 4){
    retval = IDASetJacFn(mem, jacHeat);
  } else if(MGRID == 3) {
    retval = IDASetJacFn(mem, jacHeat3);
  } else {
    /* MGRID<=2 is pure boundary points, nothing to solve */
    printf("MGRID size is too small to run.\n");
    return(1);
  }
  if(check_retval(&retval, "IDASetJacFn", 1)) return(1);

  /* Call IDACalcIC to correct the initial values. */

  retval = IDACalcIC(mem, IDA_YA_YDP_INIT, t1);
  if(check_retval(&retval, "IDACalcIC", 1)) return(1);

  /* Print output heading. */
  PrintHeader(rtol, atol);

  PrintOutput(mem, t0, uu);


  /* Loop over output times, call IDASolve, and print results. */

  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {

    retval = IDASolve(mem, tout, &tret, uu, up, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1)) return(1);

    PrintOutput(mem, tret, uu);

  }

  /* Print remaining counters and free memory. */
  retval = IDAGetNumErrTestFails(mem, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  printf("\n netf = %ld,   ncfn = %ld \n", netf, ncfn);

  IDAFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(uu);
  N_VDestroy(up);
  N_VDestroy(id);
  N_VDestroy(res);
  free(data);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * heatres: heat equation system residual function
 * This uses 5-point central differencing on the interior points, and
 * includes algebraic equations for the boundary values.
 * So for each interior point, the residual component has the form
 *    res_i = u'_i - (central difference)_i
 * while for each boundary point, it is res_i = u_i.
 */

int heatres(realtype tres, N_Vector uu, N_Vector up, N_Vector resval,
            void *user_data)
{
  sunindextype mm, i, j, offset, loc;
  realtype *uv, *upv, *resv, coeff;
  UserData data;

  uv = N_VGetArrayPointer(uu); upv = N_VGetArrayPointer(up); resv = N_VGetArrayPointer(resval);

  data = (UserData)user_data;
  mm = data->mm;
  coeff = data->coeff;

  /* Initialize resval to uu, to take care of boundary equations. */
  N_VScale(ZERO, uu, resval);

  /* Loop over interior points; set res = up - (central difference). */
  for (j = 1; j < mm-1; j++) {
    offset = mm*j;
    for (i = 1; i < mm-1; i++) {
      loc = offset + i;
      resv[loc] = upv[loc] - coeff *
	  (uv[loc-1] + uv[loc+1] + uv[loc-mm] + uv[loc+mm] - RCONST(4.0)*uv[loc]);
    }
  }

  return(0);

}

/* Jacobian matrix setup for MGRID=3  */
int jacHeat3(realtype tt,  realtype cj,
             N_Vector yy, N_Vector yp, N_Vector resvec,
             SUNMatrix JJ, void *user_data,
             N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype dx =  ONE/(MGRID - ONE);
  realtype beta = RCONST(4.0)/(dx*dx) + cj;

  sunindextype *colptrs = SUNSparseMatrix_IndexPointers(JJ);
  sunindextype *rowvals = SUNSparseMatrix_IndexValues(JJ);
  realtype *data = SUNSparseMatrix_Data(JJ);

  SUNMatZero(JJ);

  /*
   * set up number of elements in each column
   */
  colptrs[0]  = 0;
  colptrs[1]  = 1;
  colptrs[2]  = 3;
  colptrs[3]  = 4;
  colptrs[4]  = 6;
  colptrs[5]  = 7;
  colptrs[6]  = 9;
  colptrs[7]  = 10;
  colptrs[8]  = 12;
  colptrs[9]  = 13;

  /*
   * set up data and row values stored
   */

  data[0] = ONE;
  rowvals[0] = 0;
  data[1] = ONE;
  rowvals[1] = 1;
  data[2] = -ONE/(dx*dx);
  rowvals[2] = 4;
  data[3] = ONE;
  rowvals[3] = 2;
  data[4] = ONE;
  rowvals[4] = 3;
  data[5] = -ONE/(dx*dx);
  rowvals[5] = 4;
  data[6] = beta;
  rowvals[6] = 4;
  data[7] = -ONE/(dx*dx);
  rowvals[7] = 4;
  data[8] = ONE;
  rowvals[8] = 5;
  data[9] = ONE;
  rowvals[9] = 6;
  data[10] = -ONE/(dx*dx);
  rowvals[10] = 4;
  data[11] = ONE;
  rowvals[11] = 7;
  data[12] = ONE;
  rowvals[12] = 8;

  return(0);
}

/* Jacobian matrix setup for MGRID>=4  */
int jacHeat(realtype tt,  realtype cj,
            N_Vector yy, N_Vector yp, N_Vector resvec,
            SUNMatrix JJ, void *user_data,
            N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype dx =  ONE/(MGRID - ONE);
  realtype beta = RCONST(4.0)/(dx*dx) + cj;
  int i,j, repeat=0;

  sunindextype *colptrs = SUNSparseMatrix_IndexPointers(JJ);
  sunindextype *rowvals = SUNSparseMatrix_IndexValues(JJ);
  realtype *data = SUNSparseMatrix_Data(JJ);

  SUNMatZero(JJ);

  /*
   *-----------------------------------------------
   * set up number of elements in each column
   *-----------------------------------------------
   */

  /**** first column block ****/
  colptrs[0] = 0;
  colptrs[1] = 1;
  /* count by twos in the middle  */
  for(i=2;i<MGRID;i++) colptrs[i] = (colptrs[i-1])+2;
  colptrs[MGRID] = 2*MGRID-2;

  /**** second column block ****/
  colptrs[MGRID+1] = 2*MGRID;
  colptrs[MGRID+2] = 2*MGRID+3;
  /* count by fours in the middle */
  for(i=0;i<MGRID-4;i++) colptrs[MGRID+3+i] = (colptrs[MGRID+3+i-1])+4;
  colptrs[2*MGRID-1] = 2*MGRID+4*(MGRID-2)-2;
  colptrs[2*MGRID] = 2*MGRID+4*(MGRID-2);

  /**** repeated (MGRID-4 times) middle column blocks ****/
  for(i=0;i<MGRID-4;i++){
    colptrs[2*MGRID+1+repeat]   = (colptrs[2*MGRID+1+repeat-1])+2;
    colptrs[2*MGRID+1+repeat+1] = (colptrs[2*MGRID+1+repeat])+4;

    /* count by fives in the middle */
    for(j=0;j<MGRID-4;j++) colptrs[2*MGRID+1+repeat+2+j] =
			  (colptrs[2*MGRID+1+repeat+1+j])+5;

    colptrs[2*MGRID+1+repeat+(MGRID-4)+2] = (colptrs[2*MGRID+1+repeat+(MGRID-4)+1])+4;
    colptrs[2*MGRID+1+repeat+(MGRID-4)+3] = (colptrs[2*MGRID+1+repeat+(MGRID-4)+2])+2;

    repeat+=MGRID; /* shift that accounts for accumulated number of columns */
  }

  /**** last-1 column block ****/
  colptrs[MGRID*MGRID-2*MGRID+1] = TOTAL-2*MGRID-4*(MGRID-2)+2;
  colptrs[MGRID*MGRID-2*MGRID+2] = TOTAL-2*MGRID-4*(MGRID-2)+5;
  /* count by fours in the middle */
  for(i=0;i<MGRID-4;i++) colptrs[MGRID*MGRID-2*MGRID+3+i] =
			(colptrs[MGRID*MGRID-2*MGRID+3+i-1])+4;
  colptrs[MGRID*MGRID-MGRID-1] = TOTAL-2*MGRID;
  colptrs[MGRID*MGRID-MGRID]   = TOTAL-2*MGRID+2;

  /**** last column block ****/
  colptrs[MGRID*MGRID-MGRID+1] = TOTAL-MGRID-(MGRID-2)+1;
  /* count by twos in the middle */
  for(i=0;i<MGRID-2;i++) colptrs[MGRID*MGRID-MGRID+2+i] =
			(colptrs[MGRID*MGRID-MGRID+2+i-1])+2;
  colptrs[MGRID*MGRID-1] = TOTAL-1;
  colptrs[MGRID*MGRID]   = TOTAL;


  /*
   *-----------------------------------------------
   * set up data stored
   *-----------------------------------------------
   */

  /**** first column block ****/
  data[0] = ONE;
  /* alternating pattern in data, separate loop for each pattern  */
  for(i=1;i<MGRID+(MGRID-2)  ;i+=2) data[i] = ONE;
  for(i=2;i<MGRID+(MGRID-2)-1;i+=2) data[i] = -ONE/(dx*dx);

  /**** second column block ****/
  data[MGRID+MGRID-2] = ONE;
  data[MGRID+MGRID-1] = -ONE/(dx*dx);
  data[MGRID+MGRID]   = beta;
  data[MGRID+MGRID+1] = -ONE/(dx*dx);
  data[MGRID+MGRID+2] = -ONE/(dx*dx);
  /* middle data elements */
  for(i=0;i<(MGRID-4);i++) data[MGRID+MGRID+3+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) data[MGRID+MGRID+4+4*i] = beta;
  for(i=0;i<(MGRID-4);i++) data[MGRID+MGRID+5+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) data[MGRID+MGRID+6+4*i] = -ONE/(dx*dx);
  data[2*MGRID+4*(MGRID-2)-5] = -ONE/(dx*dx);
  data[2*MGRID+4*(MGRID-2)-4] = beta;
  data[2*MGRID+4*(MGRID-2)-3] = -ONE/(dx*dx);
  data[2*MGRID+4*(MGRID-2)-2] = -ONE/(dx*dx);
  data[2*MGRID+4*(MGRID-2)-1] = ONE;

  /**** repeated (MGRID-4 times) middle column blocks ****/
  repeat=0;
  for(i=0;i<MGRID-4;i++){
    data[2*MGRID+4*(MGRID-2)+repeat]   = ONE;
    data[2*MGRID+4*(MGRID-2)+repeat+1] = -ONE/(dx*dx);

    data[2*MGRID+4*(MGRID-2)+repeat+2] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+3] = beta;
    data[2*MGRID+4*(MGRID-2)+repeat+4] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+5] = -ONE/(dx*dx);

    /* 5 in 5*j chosen since there are 5 elements in each column */
    /* this column loops MGRID-4 times within the outer loop */
    for(j=0;j<MGRID-4;j++){
      data[2*MGRID+4*(MGRID-2)+repeat+6+5*j]  = -ONE/(dx*dx);
      data[2*MGRID+4*(MGRID-2)+repeat+7+5*j]  = -ONE/(dx*dx);
      data[2*MGRID+4*(MGRID-2)+repeat+8+5*j]  = beta;
      data[2*MGRID+4*(MGRID-2)+repeat+9+5*j]  = -ONE/(dx*dx);
      data[2*MGRID+4*(MGRID-2)+repeat+10+5*j] = -ONE/(dx*dx);
    }

    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+6] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+7] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+8] = beta;
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+9] = -ONE/(dx*dx);

    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+10] = -ONE/(dx*dx);
    data[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+11] = ONE;

    repeat+=MGRID+4*(MGRID-2); /* shift that accounts for accumulated columns and elements */
  }

  /**** last-1 column block ****/
  data[TOTAL-6*(MGRID-2)-4] = ONE;
  data[TOTAL-6*(MGRID-2)-3] = -ONE/(dx*dx);
  data[TOTAL-6*(MGRID-2)-2] = -ONE/(dx*dx);
  data[TOTAL-6*(MGRID-2)-1] = beta;
  data[TOTAL-6*(MGRID-2)  ] = -ONE/(dx*dx);
  /* middle data elements */
  for(i=0;i<(MGRID-4);i++) data[TOTAL-6*(MGRID-2)+1+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) data[TOTAL-6*(MGRID-2)+2+4*i] = -ONE/(dx*dx);
  for(i=0;i<(MGRID-4);i++) data[TOTAL-6*(MGRID-2)+3+4*i] = beta;
  for(i=0;i<(MGRID-4);i++) data[TOTAL-6*(MGRID-2)+4+4*i] = -ONE/(dx*dx);
  data[TOTAL-2*(MGRID-2)-7] = -ONE/(dx*dx);
  data[TOTAL-2*(MGRID-2)-6] = -ONE/(dx*dx);
  data[TOTAL-2*(MGRID-2)-5] = beta;
  data[TOTAL-2*(MGRID-2)-4] = -ONE/(dx*dx);
  data[TOTAL-2*(MGRID-2)-3] = ONE;

  /**** last column block ****/
  data[TOTAL-2*(MGRID-2)-2] = ONE;
  /* alternating pattern in data, separate loop for each pattern  */
  for(i=TOTAL-2*(MGRID-2)-1;i<TOTAL-2;i+=2) data[i] = -ONE/(dx*dx);
  for(i=TOTAL-2*(MGRID-2)  ;i<TOTAL-1;i+=2) data[i] = ONE;
  data[TOTAL-1] = ONE;

  /*
   *-----------------------------------------------
   * row values
   *-----------------------------------------------
   */

  /**** first block ****/
  rowvals[0] = 0;
  /* alternating pattern in data, separate loop for each pattern */
  for(i=1;i<MGRID+(MGRID-2)  ;i+=2) rowvals[i] = (i+1)/2;
  for(i=2;i<MGRID+(MGRID-2)-1;i+=2) rowvals[i] = i/2+MGRID; /* i+1 unnecessary here */

  /**** second column block ****/
  rowvals[MGRID+MGRID-2] = MGRID;
  rowvals[MGRID+MGRID-1] = MGRID+1;
  rowvals[MGRID+MGRID]   = MGRID+1;
  rowvals[MGRID+MGRID+1] = MGRID+2;
  rowvals[MGRID+MGRID+2] = 2*MGRID+1;
  /* middle row values */
  for(i=0;i<(MGRID-4);i++) rowvals[MGRID+MGRID+3+4*i] = MGRID+1+i;
  for(i=0;i<(MGRID-4);i++) rowvals[MGRID+MGRID+4+4*i] = MGRID+2+i;
  for(i=0;i<(MGRID-4);i++) rowvals[MGRID+MGRID+5+4*i] = MGRID+3+i;
  for(i=0;i<(MGRID-4);i++) rowvals[MGRID+MGRID+6+4*i] = 2*MGRID+2+i;
  rowvals[2*MGRID+4*(MGRID-2)-5] = MGRID+(MGRID-2)-1;
  rowvals[2*MGRID+4*(MGRID-2)-4] = MGRID+(MGRID-2); /* starting from here, add two diag patterns */
  rowvals[2*MGRID+4*(MGRID-2)-3] = 2*MGRID+(MGRID-2);
  rowvals[2*MGRID+4*(MGRID-2)-2] = MGRID+(MGRID-2);
  rowvals[2*MGRID+4*(MGRID-2)-1] = MGRID+(MGRID-2)+1;

  /**** repeated (MGRID-4 times) middle column blocks ****/
  repeat=0;
  for(i=0;i<MGRID-4;i++){
    rowvals[2*MGRID+4*(MGRID-2)+repeat]   = MGRID+(MGRID-2)+2+MGRID*i;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+1] = MGRID+(MGRID-2)+2+MGRID*i+1;

    rowvals[2*MGRID+4*(MGRID-2)+repeat+2] = MGRID+(MGRID-2)+2+MGRID*i+1-MGRID;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+3] = MGRID+(MGRID-2)+2+MGRID*i+1;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+4] = MGRID+(MGRID-2)+2+MGRID*i+2; /* *this */
    rowvals[2*MGRID+4*(MGRID-2)+repeat+5] = MGRID+(MGRID-2)+2+MGRID*i+1+MGRID;

    /* 5 in 5*j chosen since there are 5 elements in each column */
    /* column repeats MGRID-4 times within the outer loop */
    for(j=0;j<MGRID-4;j++){
      rowvals[2*MGRID+4*(MGRID-2)+repeat+6+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+1-MGRID+1+j;
      rowvals[2*MGRID+4*(MGRID-2)+repeat+7+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+1+j;
      rowvals[2*MGRID+4*(MGRID-2)+repeat+8+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+2+j;
      rowvals[2*MGRID+4*(MGRID-2)+repeat+9+5*j]  = MGRID+(MGRID-2)+2+MGRID*i+2+1+j;
      rowvals[2*MGRID+4*(MGRID-2)+repeat+10+5*j] = MGRID+(MGRID-2)+2+MGRID*i+1+MGRID+1+j;
    }

    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+6] = MGRID+(MGRID-2)+2+MGRID*i-2;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+7] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID-1;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+8] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID; /* *this+MGRID */
    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+9] = MGRID+(MGRID-2)+2+MGRID*i-2+2*MGRID;

    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+10] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID;
    rowvals[2*MGRID+4*(MGRID-2)+repeat+(MGRID-4)*5+11] = MGRID+(MGRID-2)+2+MGRID*i-2+MGRID+1;

    repeat+=MGRID+4*(MGRID-2); /* shift that accounts for accumulated columns and elements */
  }


  /**** last-1 column block ****/
  rowvals[TOTAL-6*(MGRID-2)-4] = MGRID*MGRID-1-2*(MGRID-1)-1;
  rowvals[TOTAL-6*(MGRID-2)-3] = MGRID*MGRID-1-2*(MGRID-1); /* starting with this as base */
  rowvals[TOTAL-6*(MGRID-2)-2] = MGRID*MGRID-1-2*(MGRID-1)-MGRID;
  rowvals[TOTAL-6*(MGRID-2)-1] = MGRID*MGRID-1-2*(MGRID-1);
  rowvals[TOTAL-6*(MGRID-2)  ] = MGRID*MGRID-1-2*(MGRID-1)+1;
  /* middle row values */
  for(i=0;i<(MGRID-4);i++) rowvals[TOTAL-6*(MGRID-2)+1+4*i] = MGRID*MGRID-1-2*(MGRID-1)-MGRID+1+i;
  for(i=0;i<(MGRID-4);i++) rowvals[TOTAL-6*(MGRID-2)+2+4*i] = MGRID*MGRID-1-2*(MGRID-1)+i;
  for(i=0;i<(MGRID-4);i++) rowvals[TOTAL-6*(MGRID-2)+3+4*i] = MGRID*MGRID-1-2*(MGRID-1)+1+i;/*copied above*/
  for(i=0;i<(MGRID-4);i++) rowvals[TOTAL-6*(MGRID-2)+4+4*i] = MGRID*MGRID-1-2*(MGRID-1)+2+i;
  rowvals[TOTAL-2*(MGRID-2)-7] = MGRID*MGRID-2*MGRID-2;
  rowvals[TOTAL-2*(MGRID-2)-6] = MGRID*MGRID-MGRID-3;
  rowvals[TOTAL-2*(MGRID-2)-5] = MGRID*MGRID-MGRID-2;
  rowvals[TOTAL-2*(MGRID-2)-4] = MGRID*MGRID-MGRID-2;
  rowvals[TOTAL-2*(MGRID-2)-3] = MGRID*MGRID-MGRID-1;

  /* last column block */
  rowvals[TOTAL-2*(MGRID-2)-2] = MGRID*MGRID-MGRID;
  /* alternating pattern in data, separate loop for each pattern  */
  for(i=0;i<(MGRID-2);i++) rowvals[TOTAL-2*(MGRID-2)-1+2*i] = MGRID*MGRID-2*MGRID+1+i;
  for(i=0;i<(MGRID-2);i++) rowvals[TOTAL-2*(MGRID-2)  +2*i] = MGRID*MGRID-MGRID+1+i;
  rowvals[TOTAL-1] = MGRID*MGRID-1;

  return(0);
}


/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * SetInitialProfile: routine to initialize u, up, and id vectors.
 */

static int SetInitialProfile(UserData data, N_Vector uu, N_Vector up,
                             N_Vector id, N_Vector res)
{
  realtype xfact, yfact, *udata, *updata, *iddata;
  sunindextype mm, mm1, i, j, offset, loc;

  mm = data->mm;
  mm1 = mm - 1;

  udata = N_VGetArrayPointer(uu);
  updata = N_VGetArrayPointer(up);
  iddata = N_VGetArrayPointer(id);

  /* Initialize id to 1's. */
  N_VConst(ONE, id);

  /* Initialize uu on all grid points. */
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

  /* heatres sets res to negative of ODE RHS values at interior points. */
  heatres(ZERO, uu, up, res, data);

  /* Copy -res into up to get correct interior initial up values. */
  N_VScale(-ONE, res, up);

  /* Finally, set values of u, up, and id at boundary points. */
  for (j = 0; j < mm; j++) {
    offset = mm*j;
    for (i = 0;i < mm; i++) {
      loc = offset + i;
      if (j == 0 || j == mm1 || i == 0 || i == mm1 ) {
        udata[loc] = BVAL; updata[loc] = ZERO; iddata[loc] = ZERO; }
    }
  }

  return(0);

}

/*
 * Print first lines of output (problem description)
 */

static void PrintHeader(realtype rtol, realtype atol)
{
  printf("\nidaHeat2D_klu: Heat equation, serial example problem for IDA\n");
  printf("          Discretized heat equation on 2D unit square.\n");
  printf("          Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("          Mesh dimensions: %d x %d", MGRID, MGRID);
  printf("        Total system size: %d\n\n", NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Constraints set to force all solution components >= 0. \n");
  printf("Linear solver: KLU, sparse direct solver \n");
  printf("       difference quotient Jacobian\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("IDACalcIC called with input boundary values = %Lg \n",BVAL);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("IDACalcIC called with input boundary values = %g \n",BVAL);
#else
  printf("IDACalcIC called with input boundary values = %g \n",BVAL);
#endif
  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time       umax     k  nst  nni  nje   nre     h       \n" );
  printf(" .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . \n");
}

/*
 * Print Output
 */

static void PrintOutput(void *mem, realtype t, N_Vector uu)
{
  int retval;
  realtype umax, hused;
  long int nst, nni, nje, nre;
  int kused;

  umax = N_VMaxNorm(uu);

  retval = IDAGetLastOrder(mem, &kused);
  check_retval(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  check_retval(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_retval(&retval, "IDAGetNumResEvals", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_retval(&retval, "IDAGetLastStep", 1);
  retval = IDAGetNumJacEvals(mem, &nje);
  check_retval(&retval, "IDAGetNumJacEvals", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %4ld  %9.2Le \n",
         t, umax, kused, nst, nni, nje, nre, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, nre, hused);
#else
  printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %9.2e \n",
         t, umax, kused, nst, nni, nje, nre, hused);
#endif

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
