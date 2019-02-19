/* --------------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * --------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * --------------------------------------------------------------------
 * Demonstration program for CVODES - Krylov linear solver.
 * ODE system from ns-species interaction PDE in 2 dimensions.
 * 
 * This program solves a stiff ODE system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion on
 * the unit square in two dimensions. The dependent variable vector is:
 *
 *        1   2        ns
 *  c = (c , c , ..., c  )
 *
 * and the PDEs are as follows:
 *
 *    i               i      i
 *  dc /dt  =  d(i)*(c    + c   )  +  f (x,y,c)  (i=1,...,ns)
 *                    xx     yy        i   
 *
 * where
 *
 *                 i          ns         j
 *  f (x,y,c)  =  c *(b(i) + sum a(i,j)*c )
 *   i                       j=1                                         
 *                                                                       
 * The number of species is ns = 2*np, with the first np being prey
 * and the last np being predators. The coefficients a(i,j), b(i),
 * d(i) are:
 *
 *  a(i,i) = -a  (all i)
 *  a(i,j) = -g  (i <= np, j > np)
 *  a(i,j) =  e  (i > np, j <= np)
 *  b(i) =  b*(1 + alpha*x*y)  (i <= np)
 *  b(i) = -b*(1 + alpha*x*y)  (i > np)
 *  d(i) = Dprey  (i <= np)
 *  d(i) = Dpred  (i > np)
 *
 * The spatial domain is the unit square. The final time is 10.
 * The boundary conditions are: normal derivative = 0.
 * A polynomial in x and y is used to set the initial conditions.
 *
 * The PDEs are discretized by central differencing on an MX by MY mesh.
 *
 * The resulting ODE system is stiff.
 *
 * The ODE system is solved using Newton iteration and the SUNLinSol_SPGMR
 * linear solver (scaled preconditioned GMRES).
 *
 * The preconditioner matrix used is the product of two matrices:
 * (1) A matrix, only defined implicitly, based on a fixed number
 * of Gauss-Seidel iterations using the diffusion terms only.
 * (2) A block-diagonal matrix based on the partial derivatives
 * of the interaction terms f only, using block-grouping (computing
 * only a subset of the ns by ns blocks).
 *
 * Four different runs are made for this problem.
 * The product preconditoner is applied on the left and on the
 * right. In each case, both the modified and classical Gram-Schmidt
 * options are tested.
 * In the series of runs, CVodeInit, SUNLinSol_SPGMR, and 
 * CVDlsSetLinearSolver are called only for the first run, whereas 
 * CVodeReInit, SUNLinSol_SPGMRSetPrecType, and SUNSLinSol_PGMRSetGSType
 * are called for each of the remaining three runs.
 *
 * A problem description, performance statistics at selected output
 * times, and final statistics are written to standard output.
 * On the first run, solution values are also printed at output
 * times. Error and warning messages are written to standard error,
 * but there should be no such messages.
 *
 * Note: This program requires the dense linear solver functions
 * newDenseMat, newIndexArray, denseAddIdentity, denseGETRF, denseGETRS, 
 * destroyMat and destroyArray.
 *
 * Note: This program assumes the sequential implementation for the
 * type N_Vector and uses the N_VGetArrayPointer function to gain
 * access to the contiguous array of components of an N_Vector.
 * --------------------------------------------------------------------
 * Reference: Peter N. Brown and Alan C. Hindmarsh, Reduced Storage
 * Matrix Methods in Stiff ODE Systems, J. Appl. Math. & Comp., 31
 * (1989), pp. 40-91.  Also available as Lawrence Livermore National
 * Laboratory Report UCRL-95088, Rev. 1, June 1987.
 * --------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvodes/cvodes.h>              /* main integrator header file                 */
#include <sunlinsol/sunlinsol_spgmr.h>  /* access to SPGMR SUNLinearSolver             */
#include <nvector/nvector_serial.h>     /* serial N_Vector types, fct. and macros      */
#include <sundials/sundials_dense.h>    /* use generic DENSE solver in preconditioning */
#include <sundials/sundials_types.h>    /* definition of realtype                      */
#include <sundials/sundials_math.h>     /* contains the macros ABS and SUNSQR          */

/* Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* Problem Specification Constants */

#define AA    ONE               /* AA = a */
#define EE    RCONST(1.0e4)     /* EE = e */
#define GG    RCONST(0.5e-6)    /* GG = g */
#define BB    ONE               /* BB = b */
#define DPREY ONE
#define DPRED RCONST(0.5)
#define ALPH  ONE
#define NP    3
#define NS    (2*NP)

/* Method Constants */

#define MX    6
#define MY    6
#define MXNS  (MX*NS)
#define AX    ONE
#define AY    ONE
#define DX    (AX/(realtype)(MX-1))
#define DY    (AY/(realtype)(MY-1))
#define MP    NS
#define MQ    (MX*MY)
#define MXMP  (MX*MP)
#define NGX   2
#define NGY   2
#define NGRP  (NGX*NGY)
#define ITMAX 5

/* CVodeInit Constants */

#define NEQ  (NS*MX*MY)
#define T0   ZERO
#define RTOL RCONST(1.0e-5)
#define ATOL RCONST(1.0e-5)

/* Spgmr/CVLS Constants */

#define MAXL 0     /* => use default = MIN(NEQ, 5)            */
#define DELT ZERO  /* => use default = 0.05                   */

/* Output Constants */

#define T1        RCONST(1.0e-8)
#define TOUT_MULT RCONST(10.0)
#define DTOUT     ONE
#define NOUT      18

/* Note: The value for species i at mesh point (j,k) is stored in */
/* component number (i-1) + j*NS + k*NS*MX of an N_Vector,        */
/* where 1 <= i <= NS, 0 <= j < MX, 0 <= k < MY.                  */

/* Structure for user data */

typedef struct {
  realtype **P[NGRP];
  sunindextype *pivot[NGRP];
  int ns, mxns;
  int mp, mq, mx, my, ngrp, ngx, ngy, mxmp;
  int jgx[NGX+1], jgy[NGY+1], jigx[MX], jigy[MY];
  int jxr[NGX], jyr[NGY];
  realtype acoef[NS][NS], bcoef[NS], diff[NS];
  realtype cox[NS], coy[NS], dx, dy, srur;
  realtype fsave[NEQ];
  N_Vector tmp;
  N_Vector rewt;
  void *cvode_mem;
} *WebData;

/* Private Helper Functions */

static WebData AllocUserData(void);
static void InitUserData(WebData wdata);
static void SetGroups(int m, int ng, int jg[], int jig[], int jr[]);
static void CInit(N_Vector c, WebData wdata);
static void PrintIntro(void);
static void PrintHeader(int jpre, int gstype);
static void PrintAllSpecies(N_Vector c, int ns, int mxns, realtype t);
static void PrintOutput(void *cvode_mem, realtype t);
static void PrintFinalStats(void *cvode_mem);
static void FreeUserData(WebData wdata);
static void WebRates(realtype x, realtype y, realtype t, realtype c[],
		     realtype rate[], WebData wdata);
static void fblock (realtype t, realtype cdata[], int jx, int jy,
		    realtype cdotdata[], WebData wdata);
static void GSIter(realtype gamma, N_Vector z, N_Vector x, WebData wdata);

/* Small Vector Kernels */

static void v_inc_by_prod(realtype u[], realtype v[], realtype w[], int n);
static void v_sum_prods(realtype u[], realtype p[], realtype q[], realtype v[],
                        realtype w[], int n);
static void v_prod(realtype u[], realtype v[], realtype w[], int n);
static void v_zero(realtype u[], int n);

/* Functions Called By The Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Precond(realtype tn, N_Vector c, N_Vector fc, booleantype jok, 
                   booleantype *jcurPtr, realtype gamma, void *user_data);

static int PSolve(realtype tn, N_Vector c, N_Vector fc, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Implementation */

int main()
{
  realtype abstol=ATOL, reltol=RTOL, t, tout;
  N_Vector c;
  WebData wdata;
  SUNLinearSolver LS;
  void *cvode_mem;
  booleantype firstrun;
  int jpre, gstype, retval;
  int ns, mxns, iout;

  c = NULL;
  wdata = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Initializations */
  c = N_VNew_Serial(NEQ);
  if(check_retval((void *)c, "N_VNew_Serial", 0)) return(1);
  wdata = AllocUserData();
  if(check_retval((void *)wdata, "AllocUserData", 2)) return(1);
  InitUserData(wdata);
  ns = wdata->ns;
  mxns = wdata->mxns;

  /* Print problem description */
  PrintIntro();

  /* Loop over jpre and gstype (four cases) */
  for (jpre = PREC_LEFT; jpre <= PREC_RIGHT; jpre++) {
    for (gstype = MODIFIED_GS; gstype <= CLASSICAL_GS; gstype++) {

      /* Initialize c and print heading */
      CInit(c, wdata);
      PrintHeader(jpre, gstype);

      /* Call CVodeInit or CVodeReInit, then SUNLinSol_SPGMR to set up problem */

      firstrun = (jpre == PREC_LEFT) && (gstype == MODIFIED_GS);
      if (firstrun) {
        cvode_mem = CVodeCreate(CV_BDF);
        if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

        wdata->cvode_mem = cvode_mem;

        retval = CVodeSetUserData(cvode_mem, wdata);
        if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

        retval = CVodeInit(cvode_mem, f, T0, c);
        if(check_retval(&retval, "CVodeInit", 1)) return(1);

        retval = CVodeSStolerances(cvode_mem, reltol, abstol);
        if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

        LS = SUNLinSol_SPGMR(c, jpre, MAXL);
        if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

        retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
        if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

        retval = SUNLinSol_SPGMRSetGSType(LS, gstype);
        if(check_retval(&retval, "SUNLinSol_SPGMRSetGSType", 1)) return(1);

        retval = CVodeSetEpsLin(cvode_mem, DELT);
        if(check_retval(&retval, "CVodeSetEpsLin", 1)) return(1);

        retval = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
        if(check_retval(&retval, "CVodeSetPreconditioner", 1)) return(1);

      } else {

        retval = CVodeReInit(cvode_mem, T0, c);
        if(check_retval(&retval, "CVodeReInit", 1)) return(1);

        retval = SUNLinSol_SPGMRSetPrecType(LS, jpre);
        if(check_retval(&retval, "SUNLinSol_SPGMRSetPrecType", 1)) return(1);
        retval = SUNLinSol_SPGMRSetGSType(LS, gstype);
        if(check_retval(&retval, "SUNLinSol_SPGMRSetGSType", 1)) return(1);

      }
      
      /* Print initial values */
      if (firstrun) PrintAllSpecies(c, ns, mxns, T0);
      
      /* Loop over output points, call CVode, print sample solution values. */
      tout = T1;
      for (iout = 1; iout <= NOUT; iout++) {
        retval = CVode(cvode_mem, tout, c, &t, CV_NORMAL);
        PrintOutput(cvode_mem, t);
        if (firstrun && (iout % 3 == 0)) PrintAllSpecies(c, ns, mxns, t);
        if(check_retval(&retval, "CVode", 1)) break;
        if (tout > RCONST(0.9)) tout += DTOUT; else tout *= TOUT_MULT; 
      }
      
      /* Print final statistics, and loop for next case */
      PrintFinalStats(cvode_mem);
      
    }
  }

  /* Free all memory */
  CVodeFree(&cvode_mem);
  N_VDestroy(c);
  SUNLinSolFree(LS);
  FreeUserData(wdata);

  return(0);
}

static WebData AllocUserData(void)
{
  int i, ngrp = NGRP;
  sunindextype ns = NS;
  WebData wdata;
  
  wdata = (WebData) malloc(sizeof *wdata);
  for(i=0; i < ngrp; i++) {
    (wdata->P)[i] = newDenseMat(ns, ns);
    (wdata->pivot)[i] = newIndexArray(ns);
  }
  wdata->rewt = N_VNew_Serial(NEQ);
  wdata->tmp = N_VNew_Serial(NEQ);
  return(wdata);
}

static void InitUserData(WebData wdata)
{
  int i, j, ns;
  realtype *bcoef, *diff, *cox, *coy, dx, dy;
  realtype (*acoef)[NS];
  
  acoef = wdata->acoef;
  bcoef = wdata->bcoef;
  diff = wdata->diff;
  cox = wdata->cox;
  coy = wdata->coy;
  ns = wdata->ns = NS;
  
  for (j = 0; j < NS; j++) { for (i = 0; i < NS; i++) acoef[i][j] = 0.; }
  for (j = 0; j < NP; j++) {
    for (i = 0; i < NP; i++) {
      acoef[NP+i][j] = EE;
      acoef[i][NP+j] = -GG;
    }
    acoef[j][j] = -AA;
    acoef[NP+j][NP+j] = -AA;
    bcoef[j] = BB;
    bcoef[NP+j] = -BB;
    diff[j] = DPREY;
    diff[NP+j] = DPRED;
  }

  /* Set remaining problem parameters */

  wdata->mxns = MXNS;
  dx = wdata->dx = DX;
  dy = wdata->dy = DY;
  for (i = 0; i < ns; i++) {
    cox[i] = diff[i]/SUNSQR(dx);
    coy[i] = diff[i]/SUNSQR(dy);
  }

  /* Set remaining method parameters */

  wdata->mp = MP;
  wdata->mq = MQ;
  wdata->mx = MX;
  wdata->my = MY;
  wdata->srur = SUNRsqrt(UNIT_ROUNDOFF);
  wdata->mxmp = MXMP;
  wdata->ngrp = NGRP;
  wdata->ngx = NGX;
  wdata->ngy = NGY;
  SetGroups(MX, NGX, wdata->jgx, wdata->jigx, wdata->jxr);
  SetGroups(MY, NGY, wdata->jgy, wdata->jigy, wdata->jyr);
}

/*
 This routine sets arrays jg, jig, and jr describing
 a uniform partition of (0,1,2,...,m-1) into ng groups.
 The arrays set are:
   jg    = length ng+1 array of group boundaries.
           Group ig has indices j = jg[ig],...,jg[ig+1]-1.
   jig   = length m array of group indices vs node index.
           Node index j is in group jig[j].
   jr    = length ng array of indices representing the groups.
           The index for group ig is j = jr[ig].
*/
static void SetGroups(int m, int ng, int jg[], int jig[], int jr[])
{
  int ig, j, len1, mper, ngm1;

  mper = m/ng; /* does integer division */
  for (ig=0; ig < ng; ig++) jg[ig] = ig*mper;
  jg[ng] = m;
  
  ngm1 = ng - 1;
  len1 = ngm1*mper;
  for (j = 0; j < len1; j++) jig[j] = j/mper;
  for (j = len1; j < m; j++) jig[j] = ngm1;

  for (ig = 0; ig < ngm1; ig++) jr[ig] = ((2*ig+1)*mper-1)/2;
  jr[ngm1] = (ngm1*mper+m-1)/2;
}

/* This routine computes and loads the vector of initial values. */
static void CInit(N_Vector c, WebData wdata)
{
  int jx, jy, ns, mxns, ioff, iyoff, i, ici;
  realtype argx, argy, x, y, dx, dy, x_factor, y_factor, *cdata;
  
  cdata = N_VGetArrayPointer(c);
  ns = wdata->ns;
  mxns = wdata->mxns;
  dx = wdata->dx;
  dy = wdata->dy;

  x_factor = RCONST(4.0)/SUNSQR(AX);
  y_factor = RCONST(4.0)/SUNSQR(AY);
  for (jy = 0; jy < MY; jy++) {
    y = jy*dy;
    argy = SUNSQR(y_factor*y*(AY-y));
    iyoff = mxns*jy;
    for (jx = 0; jx < MX; jx++) {
      x = jx*dx;
      argx = SUNSQR(x_factor*x*(AX-x));
      ioff = iyoff + ns*jx;
      for (i = 1; i <= ns; i++) {
        ici = ioff + i-1;
        cdata[ici] = RCONST(10.0) + i*argx*argy;
      }
    }
  }
}

static void PrintIntro(void)
{
  printf("\n\nDemonstration program for CVODES - SPGMR linear solver\n\n");
  printf("Food web problem with ns species, ns = %d\n", NS);
  printf("Predator-prey interaction and diffusion on a 2-D square\n\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Matrix parameters: a = %.2Lg   e = %.2Lg   g = %.2Lg\n",
         AA, EE, GG);
  printf("b parameter = %.2Lg\n", BB);
  printf("Diffusion coefficients: Dprey = %.2Lg   Dpred = %.2Lg\n",
         DPREY, DPRED);
  printf("Rate parameter alpha = %.2Lg\n\n", ALPH);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Matrix parameters: a = %.2g   e = %.2g   g = %.2g\n",
         AA, EE, GG);
  printf("b parameter = %.2g\n", BB);
  printf("Diffusion coefficients: Dprey = %.2g   Dpred = %.2g\n",
         DPREY, DPRED);
  printf("Rate parameter alpha = %.2g\n\n", ALPH);
#else
  printf("Matrix parameters: a = %.2g   e = %.2g   g = %.2g\n",
         AA, EE, GG);
  printf("b parameter = %.2g\n", BB);
  printf("Diffusion coefficients: Dprey = %.2g   Dpred = %.2g\n",
         DPREY, DPRED);
  printf("Rate parameter alpha = %.2g\n\n", ALPH);
#endif
  printf("Mesh dimensions (mx,my) are %d, %d.  ", MX, MY);
  printf("Total system size is neq = %d \n\n", NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerances: reltol = %.2Lg, abstol = %.2Lg \n\n",
         RTOL, ATOL);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerances: reltol = %.2g, abstol = %.2g \n\n",
         RTOL, ATOL);
#else
  printf("Tolerances: reltol = %.2g, abstol = %.2g \n\n",
         RTOL, ATOL);
#endif
  printf("Preconditioning uses a product of:\n");
  printf("  (1) Gauss-Seidel iterations with ");
  printf("itmax = %d iterations, and\n", ITMAX);
  printf("  (2) interaction-only block-diagonal matrix ");
  printf("with block-grouping\n");
  printf("  Number of diagonal block groups = ngrp = %d", NGRP);
  printf("  (ngx by ngy, ngx = %d, ngy = %d)\n", NGX, NGY);
  printf("\n\n--------------------------------------------------------------");
  printf("--------------\n");
}

static void PrintHeader(int jpre, int gstype)
{
  if(jpre == PREC_LEFT)
    printf("\n\nPreconditioner type is           jpre = %s\n", "PREC_LEFT");
  else
    printf("\n\nPreconditioner type is           jpre = %s\n", "PREC_RIGHT");

  if(gstype == MODIFIED_GS)
    printf("\nGram-Schmidt method type is    gstype = %s\n\n\n", "MODIFIED_GS");
  else
    printf("\nGram-Schmidt method type is    gstype = %s\n\n\n", "CLASSICAL_GS");
}

static void PrintAllSpecies(N_Vector c, int ns, int mxns, realtype t)
{
  int i, jx ,jy;
  realtype *cdata;
  
  cdata = N_VGetArrayPointer(c);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("c values at t = %Lg:\n\n", t);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("c values at t = %g:\n\n", t);
#else
  printf("c values at t = %g:\n\n", t);
#endif
  for (i=1; i <= ns; i++) {
    printf("Species %d\n", i);
    for (jy=MY-1; jy >= 0; jy--) {
      for (jx=0; jx < MX; jx++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%-10.6Lg", cdata[(i-1) + jx*ns + jy*mxns]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%-10.6g", cdata[(i-1) + jx*ns + jy*mxns]);
#else
        printf("%-10.6g", cdata[(i-1) + jx*ns + jy*mxns]);
#endif
      }
      printf("\n");
    }
    printf("\n");
  }
}

static void PrintOutput(void *cvode_mem, realtype t)
{
  long int nst, nfe, nni;
  int qu, retval;
  realtype hu;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetLastOrder(cvode_mem, &qu);
  check_retval(&retval, "CVodeGetLastOrder", 1);
  retval = CVodeGetLastStep(cvode_mem, &hu);
  check_retval(&retval, "CVodeGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("t = %10.2Le  nst = %ld  nfe = %ld  nni = %ld", t, nst, nfe, nni);
  printf("  qu = %d  hu = %11.2Le\n\n", qu, hu);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("t = %10.2e  nst = %ld  nfe = %ld  nni = %ld", t, nst, nfe, nni);
  printf("  qu = %d  hu = %11.2e\n\n", qu, hu);
#else
  printf("t = %10.2e  nst = %ld  nfe = %ld  nni = %ld", t, nst, nfe, nni);
  printf("  qu = %d  hu = %11.2e\n\n", qu, hu);
#endif
}

static void PrintFinalStats(void *cvode_mem)
{
  long int lenrw, leniw ;
  long int lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int retval;
  realtype avdim;
  
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

  printf("\n\n Final statistics for this run:\n\n");
  printf(" CVode real workspace length           = %4ld \n", lenrw);
  printf(" CVode integer workspace length        = %4ld \n", leniw);
  printf(" CVLS real workspace length            = %4ld \n", lenrwLS);
  printf(" CVLS integer workspace length         = %4ld \n", leniwLS);
  printf(" Number of steps                       = %4ld \n", nst);
  printf(" Number of f-s                         = %4ld \n", nfe);
  printf(" Number of f-s (SPGMR)                 = %4ld \n", nfeLS);
  printf(" Number of f-s (TOTAL)                 = %4ld \n", nfe + nfeLS);
  printf(" Number of setups                      = %4ld \n", nsetups);
  printf(" Number of nonlinear iterations        = %4ld \n", nni);
  printf(" Number of linear iterations           = %4ld \n", nli);
  printf(" Number of preconditioner evaluations  = %4ld \n", npe);
  printf(" Number of preconditioner solves       = %4ld \n", nps);
  printf(" Number of error test failures         = %4ld \n", netf);
  printf(" Number of nonlinear conv. failures    = %4ld \n", ncfn);
  printf(" Number of linear convergence failures = %4ld \n", ncfl);
  avdim = (nni > 0) ? ((realtype)nli)/((realtype)nni) : ZERO;
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" Average Krylov subspace dimension     = %.3Lf \n", avdim);
#else
  printf(" Average Krylov subspace dimension     = %.3f \n", avdim);
#endif
  printf("\n\n--------------------------------------------------------------");
  printf("--------------\n");
  printf(    "--------------------------------------------------------------");
  printf("--------------\n");
}

static void FreeUserData(WebData wdata)
{
  int i, ngrp;

  ngrp = wdata->ngrp;
  for(i=0; i < ngrp; i++) {
    destroyMat((wdata->P)[i]);
    destroyArray((wdata->pivot)[i]);
  }
  N_VDestroy(wdata->rewt);
  N_VDestroy(wdata->tmp);
  free(wdata);
}

/*
 This routine computes the right-hand side of the ODE system and
 returns it in cdot. The interaction rates are computed by calls to WebRates,
 and these are saved in fsave for use in preconditioning.
*/
static int f(realtype t, N_Vector c, N_Vector cdot,void *user_data)
{
  int i, ic, ici, idxl, idxu, jx, ns, mxns, iyoff, jy, idyu, idyl;
  realtype dcxli, dcxui, dcyli, dcyui, x, y, *cox, *coy, *fsave, dx, dy;
  realtype *cdata, *cdotdata;
  WebData wdata;
  
  wdata = (WebData) user_data;
  cdata = N_VGetArrayPointer(c);
  cdotdata = N_VGetArrayPointer(cdot);
  
  mxns = wdata->mxns;
  ns = wdata->ns;
  fsave = wdata->fsave;
  cox = wdata->cox;
  coy = wdata->coy;
  mxns = wdata->mxns;
  dx = wdata->dx;
  dy = wdata->dy;
  
  for (jy = 0; jy < MY; jy++) {
    y = jy*dy;
    iyoff = mxns*jy;
    idyu = (jy == MY-1) ? -mxns : mxns;
    idyl = (jy == 0) ? -mxns : mxns;
    for (jx = 0; jx < MX; jx++) {
      x = jx*dx;
      ic = iyoff + ns*jx;
      /* Get interaction rates at one point (x,y). */
      WebRates(x, y, t, cdata+ic, fsave+ic, wdata);
      idxu = (jx == MX-1) ? -ns : ns;
      idxl = (jx == 0) ? -ns : ns;
      for (i = 1; i <= ns; i++) {
        ici = ic + i-1;
        /* Do differencing in y. */
        dcyli = cdata[ici] - cdata[ici-idyl];
        dcyui = cdata[ici+idyu] - cdata[ici];
        /* Do differencing in x. */
        dcxli = cdata[ici] - cdata[ici-idxl];
        dcxui = cdata[ici+idxu] - cdata[ici];
        /* Collect terms and load cdot elements. */
        cdotdata[ici] = coy[i-1]*(dcyui - dcyli) + cox[i-1]*(dcxui - dcxli) +
          fsave[ici];
      }
    }
  }

  return(0);
}

/*
  This routine computes the interaction rates for the species
  c_1, ... ,c_ns (stored in c[0],...,c[ns-1]), at one spatial point 
  and at time t.
*/
static void WebRates(realtype x, realtype y, realtype t, realtype c[],
                     realtype rate[], WebData wdata)
{
  int i, j, ns;
  realtype fac, *bcoef;
  realtype (*acoef)[NS];
  
  ns = wdata->ns;
  acoef = wdata->acoef;
  bcoef = wdata->bcoef;
  
  for (i = 0; i < ns; i++)
    rate[i] = ZERO;
  
  for (j = 0; j < ns; j++) 
    for (i = 0; i < ns; i++) 
      rate[i] += c[j] * acoef[i][j];
  
  fac = ONE + ALPH*x*y;
  for (i = 0; i < ns; i++) 
    rate[i] = c[i]*(bcoef[i]*fac + rate[i]);
}

/*
 This routine generates the block-diagonal part of the Jacobian
 corresponding to the interaction rates, multiplies by -gamma, adds
 the identity matrix, and calls denseGETRF to do the LU decomposition of
 each diagonal block. The computation of the diagonal blocks uses
 the preset block and grouping information. One block per group is
 computed. The Jacobian elements are generated by difference
 quotients using calls to the routine fblock.

 This routine can be regarded as a prototype for the general case
 of a block-diagonal preconditioner. The blocks are of size mp, and
 there are ngrp=ngx*ngy blocks computed in the block-grouping scheme.
*/ 
static int Precond(realtype t, N_Vector c, N_Vector fc, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, void *user_data)
{
  realtype ***P;
  sunindextype **pivot;
  int i, if0, if00, ig, igx, igy, j, jj, jx, jy;
  int *jxr, *jyr, ngrp, ngx, ngy, mxmp, retval;
  sunindextype mp, denseretval;
  realtype uround, fac, r, r0, save, srur;
  realtype *f1, *fsave, *cdata, *rewtdata;
  WebData wdata;
  void *cvode_mem;
  N_Vector rewt;
  
  wdata = (WebData) user_data;
  cvode_mem = wdata->cvode_mem;
  cdata = N_VGetArrayPointer(c);
  rewt = wdata->rewt;
  retval = CVodeGetErrWeights(cvode_mem, rewt);
  if(check_retval(&retval, "CVodeGetErrWeights", 1)) return(1);
  rewtdata = N_VGetArrayPointer(rewt);

  uround = UNIT_ROUNDOFF;

  P = wdata->P;
  pivot = wdata->pivot;
  jxr = wdata->jxr;
  jyr = wdata->jyr;
  mp = wdata->mp;
  srur = wdata->srur;
  ngrp = wdata->ngrp;
  ngx = wdata->ngx;
  ngy = wdata->ngy;
  mxmp = wdata->mxmp;
  fsave = wdata->fsave;
  
  /* Make mp calls to fblock to approximate each diagonal block of Jacobian.
     Here, fsave contains the base value of the rate vector and 
     r0 is a minimum increment factor for the difference quotient. */
  
  f1 = N_VGetArrayPointer(wdata->tmp);
  
  fac = N_VWrmsNorm (fc, rewt);
  r0 = RCONST(1000.0)*SUNRabs(gamma)*uround*NEQ*fac;
  if (r0 == ZERO) r0 = ONE;
  
  for (igy = 0; igy < ngy; igy++) {
    jy = jyr[igy];
    if00 = jy*mxmp;
    for (igx = 0; igx < ngx; igx++) { 
      jx = jxr[igx];
      if0 = if00 + jx*mp;
      ig = igx + igy*ngx; 
      /* Generate ig-th diagonal block */
      for (j = 0; j < mp; j++) {
        /* Generate the jth column as a difference quotient */
        jj = if0 + j; 
        save = cdata[jj];
        r = SUNMAX(srur*SUNRabs(save),r0/rewtdata[jj]);
        cdata[jj] += r;
        fac = -gamma/r;
        fblock (t, cdata, jx, jy, f1, wdata);
        for (i = 0; i < mp; i++) {
          P[ig][j][i] = (f1[i] - fsave[if0+i])*fac;
        }
        cdata[jj] = save;
      }
    }
  }
  
  /* Add identity matrix and do LU decompositions on blocks. */
  
  for (ig = 0; ig < ngrp; ig++) {
    denseAddIdentity(P[ig], mp);
    denseretval = denseGETRF(P[ig], mp, mp, pivot[ig]);
    if (denseretval != 0) return(1);
  }
  
  *jcurPtr = SUNTRUE;
  return(0);
}

/*
  This routine computes one block of the interaction terms of the
  system, namely block (jx,jy), for use in preconditioning.
  Here jx and jy count from 0.
*/
static void fblock(realtype t, realtype cdata[], int jx, int jy,
                   realtype cdotdata[], WebData wdata)
{
  int iblok, ic;
  realtype x, y;
  
  iblok = jx + jy*(wdata->mx);
  y = jy*(wdata->dy);
  x = jx*(wdata->dx);
  ic = (wdata->ns)*(iblok);
  WebRates(x, y, t, cdata+ic, cdotdata, wdata);
}

/*
  This routine applies two inverse preconditioner matrices
  to the vector r, using the interaction-only block-diagonal Jacobian
  with block-grouping, denoted Jr, and Gauss-Seidel applied to the
  diffusion contribution to the Jacobian, denoted Jd.
  It first calls GSIter for a Gauss-Seidel approximation to
  ((I - gamma*Jd)-inverse)*r, and stores the result in z.
  Then it computes ((I - gamma*Jr)-inverse)*z, using LU factors of the
  blocks in P, and pivot information in pivot, and returns the result in z.
*/
static int PSolve(realtype tn, N_Vector c, N_Vector fc, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  realtype   ***P;
  sunindextype **pivot;
  int jx, jy, igx, igy, iv, ig, *jigx, *jigy, mx, my, ngx;
  sunindextype mp;
  WebData wdata;
  
  wdata = (WebData) user_data;
  
  N_VScale(ONE, r, z);
  
  /* call GSIter for Gauss-Seidel iterations */
  
  GSIter(gamma, z, wdata->tmp, wdata);
  
  /* Do backsolves for inverse of block-diagonal preconditioner factor */
  
  P = wdata->P;
  pivot = wdata->pivot;
  mx = wdata->mx;
  my = wdata->my;
  ngx = wdata->ngx;
  mp = wdata->mp;
  jigx = wdata->jigx;
  jigy = wdata->jigy;
  
  iv = 0;
  for (jy = 0; jy < my; jy++) {
    igy = jigy[jy];
    for (jx = 0; jx < mx; jx++) {
      igx = jigx[jx];
      ig = igx + igy*ngx;
      denseGETRS(P[ig], mp, pivot[ig], &(N_VGetArrayPointer(z)[iv]));
      iv += mp;
    }
  }
  
  return(0);
}

/*
  This routine performs ITMAX=5 Gauss-Seidel iterations to compute an
  approximation to (P-inverse)*z, where P = I - gamma*Jd, and
  Jd represents the diffusion contributions to the Jacobian.
  The answer is stored in z on return, and x is a temporary vector.
  The dimensions below assume a global constant NS >= ns.
  Some inner loops of length ns are implemented with the small
  vector kernels v_sum_prods, v_prod, v_inc_by_prod.
*/
static void GSIter(realtype gamma, N_Vector z, N_Vector x, WebData wdata)
{
  int jx, jy, mx, my, x_loc, y_loc;
  int ns, mxns, i, iyoff, ic, iter;
  realtype beta[NS], beta2[NS], cof1[NS], gam[NS], gam2[NS];
  realtype temp, *cox, *coy, *xd, *zd;
  
  xd = N_VGetArrayPointer(x);
  zd = N_VGetArrayPointer(z);
  ns = wdata->ns;
  mx = wdata->mx;
  my = wdata->my;
  mxns = wdata->mxns;
  cox = wdata->cox;
  coy = wdata->coy;
  
  /* Write matrix as P = D - L - U.
     Load local arrays beta, beta2, gam, gam2, and cof1. */
  
  for (i = 0; i < ns; i++) {
    temp = ONE/(ONE + RCONST(2.0)*gamma*(cox[i] + coy[i]));
    beta[i] = gamma*cox[i]*temp;
    beta2[i] = RCONST(2.0)*beta[i];
    gam[i] = gamma*coy[i]*temp;
    gam2[i] = RCONST(2.0)*gam[i];
    cof1[i] = temp;
  }
  
  /* Begin iteration loop.
     Load vector x with (D-inverse)*z for first iteration. */
  
  for (jy = 0; jy < my; jy++) {
    iyoff = mxns*jy;
    for (jx = 0; jx < mx; jx++) {
      ic = iyoff + ns*jx;
      v_prod(xd+ic, cof1, zd+ic, ns); /* x[ic+i] = cof1[i]z[ic+i] */
    }
  }
  N_VConst(ZERO, z);
  
  /* Looping point for iterations. */
  
  for (iter=1; iter <= ITMAX; iter++) {
    
    /* Calculate (D-inverse)*U*x if not the first iteration. */
    
    if (iter > 1) {
      for (jy=0; jy < my; jy++) {
        iyoff = mxns*jy;
        for (jx=0; jx < mx; jx++) { /* order of loops matters */
          ic = iyoff + ns*jx;
          x_loc = (jx == 0) ? 0 : ((jx == mx-1) ? 2 : 1);
          y_loc = (jy == 0) ? 0 : ((jy == my-1) ? 2 : 1);
          switch (3*y_loc+x_loc) {
          case 0 : 
            /* jx == 0, jy == 0 */
            /* x[ic+i] = beta2[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] */
            v_sum_prods(xd+ic, beta2, xd+ic+ns, gam2, xd+ic+mxns, ns);
            break;
          case 1 : 
            /* 1 <= jx <= mx-2, jy == 0 */
            /* x[ic+i] = beta[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] */
            v_sum_prods(xd+ic, beta, xd+ic+ns, gam2, xd+ic+mxns, ns);
            break;
          case 2 : 
            /* jx == mx-1, jy == 0 */
            /* x[ic+i] = gam2[i]x[ic+mxns+i] */
            v_prod(xd+ic, gam2, xd+ic+mxns, ns);
            break;
          case 3 : 
            /* jx == 0, 1 <= jy <= my-2 */
            /* x[ic+i] = beta2[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] */
            v_sum_prods(xd+ic, beta2, xd+ic+ns, gam, xd+ic+mxns, ns);
            break;
          case 4 : 
            /* 1 <= jx <= mx-2, 1 <= jy <= my-2 */
            /* x[ic+i] = beta[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] */
            v_sum_prods(xd+ic, beta, xd+ic+ns, gam, xd+ic+mxns, ns);
            break;
          case 5 : 
            /* jx == mx-1, 1 <= jy <= my-2 */
            /* x[ic+i] = gam[i]x[ic+mxns+i] */
            v_prod(xd+ic, gam, xd+ic+mxns, ns);
            break;
          case 6 : 
            /* jx == 0, jy == my-1 */
            /* x[ic+i] = beta2[i]x[ic+ns+i] */
            v_prod(xd+ic, beta2, xd+ic+ns, ns);
            break;
          case 7 : 
            /* 1 <= jx <= mx-2, jy == my-1 */
            /* x[ic+i] = beta[i]x[ic+ns+i] */
            v_prod(xd+ic, beta, xd+ic+ns, ns);
            break;
          case 8 : 
            /* jx == mx-1, jy == my-1 */
            /* x[ic+i] = 0.0 */
            v_zero(xd+ic, ns);
            break;
          }
        }
      }
    }  /* end if (iter > 1) */
    
    /* Overwrite x with [(I - (D-inverse)*L)-inverse]*x. */
    
    for (jy=0; jy < my; jy++) {
      iyoff = mxns*jy;
      for (jx=0; jx < mx; jx++) { /* order of loops matters */
        ic = iyoff + ns*jx;
        x_loc = (jx == 0) ? 0 : ((jx == mx-1) ? 2 : 1);
        y_loc = (jy == 0) ? 0 : ((jy == my-1) ? 2 : 1);
        switch (3*y_loc+x_loc) {
        case 0 : 
          /* jx == 0, jy == 0 */
          break;
        case 1 : 
          /* 1 <= jx <= mx-2, jy == 0 */
          /* x[ic+i] += beta[i]x[ic-ns+i] */
          v_inc_by_prod(xd+ic, beta, xd+ic-ns, ns);
          break;
        case 2 : 
          /* jx == mx-1, jy == 0 */
          /* x[ic+i] += beta2[i]x[ic-ns+i] */
          v_inc_by_prod(xd+ic, beta2, xd+ic-ns, ns);
          break;
        case 3 : 
          /* jx == 0, 1 <= jy <= my-2 */
          /* x[ic+i] += gam[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, gam, xd+ic-mxns, ns);
          break;
        case 4 : 
          /* 1 <= jx <= mx-2, 1 <= jy <= my-2 */
          /* x[ic+i] += beta[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, beta, xd+ic-ns, ns);
          v_inc_by_prod(xd+ic, gam, xd+ic-mxns, ns);
          break;
        case 5 : 
          /* jx == mx-1, 1 <= jy <= my-2 */
          /* x[ic+i] += beta2[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, beta2, xd+ic-ns, ns);
          v_inc_by_prod(xd+ic, gam, xd+ic-mxns, ns);
          break;
        case 6 : 
          /* jx == 0, jy == my-1 */
          /* x[ic+i] += gam2[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, gam2, xd+ic-mxns, ns);
          break;
        case 7 : 
          /* 1 <= jx <= mx-2, jy == my-1 */
          /* x[ic+i] += beta[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, beta, xd+ic-ns, ns);
          v_inc_by_prod(xd+ic, gam2, xd+ic-mxns, ns);
          break;
        case 8 : 
          /* jx == mx-1, jy == my-1 */
          /* x[ic+i] += beta2[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, beta2, xd+ic-ns, ns);
          v_inc_by_prod(xd+ic, gam2, xd+ic-mxns, ns);
          break;
        }
      }
    }
    
    /* Add increment x to z : z <- z+x */
    
    N_VLinearSum(ONE, z, ONE, x, z);
    
  }
}

static void v_inc_by_prod(realtype u[], realtype v[], realtype w[], int n)
{
  int i;  
  for (i=0; i < n; i++) u[i] += v[i]*w[i];
}

static void v_sum_prods(realtype u[], realtype p[], realtype q[],
                        realtype v[], realtype w[], int n)
{
  int i;  
  for (i=0; i < n; i++) u[i] = p[i]*q[i] + v[i]*w[i];
}

static void v_prod(realtype u[], realtype v[], realtype w[], int n)
{ 
  int i;
  for (i=0; i < n; i++) u[i] = v[i]*w[i];
}

static void v_zero(realtype u[], int n)
{
  int i;  
  for (i=0; i < n; i++) u[i] = ZERO;
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
