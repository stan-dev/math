/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
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
 * This program solves a stiff ODE system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion on
 * the unit square in two dimensions. The dependent variable vector
 * is the following:
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
 * The PDEs are discretized by central differencing on an MX by
 * MY mesh. The resulting ODE system is stiff.
 *
 * The ODE system is solved by CVODES using Newton iteration and
 * the SUNLinSol_SPGMR linear solver (scaled preconditioned GMRES).
 *
 * The preconditioner matrix used is the product of two matrices:
 * (1) A matrix, only defined implicitly, based on a fixed number
 * of Gauss-Seidel iterations using the diffusion terms only.
 * (2) A block-diagonal matrix based on the partial derivatives
 * of the interaction terms f only, using block-grouping (computing
 * only a subset of the ns by ns blocks).
 *
 * Additionally, CVODES can integrate backwards in time the
 * the semi-discrete form of the adjoint PDE:
 *   d(lambda)/dt = - D^T ( lambda_xx + lambda_yy )
 *                  - F_c^T lambda
 * with homogeneous Neumann boundary conditions and final conditions
 *   lambda(x,y,t=t_final) = - g_c^T(t_final)
 * whose solution at t = 0 represents the sensitivity of
 *   int_x int _y g(t_final,c) dx dy dt
 * with respect to the initial conditions of the original problem.
 *
 * In this example,
 *   g(t,c) = c(ISPEC), with ISPEC defined below.
 * -----------------------------------------------------------------
 * Reference:  Peter N. Brown and Alan C. Hindmarsh, Reduced Storage
 * Matrix Methods in Stiff ODE Systems, J. Appl. Math. & Comp., 31
 * (1989), pp. 40-91.  Also available as Lawrence Livermore National
 * Laboratory Report UCRL-95088, Rev. 1, June 1987.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cvodes/cvodes.h>             /* main integrator header file          */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver      */
#include <sundials/sundials_dense.h>   /* use generic dense solver in precond. */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* contains the macros ABS, SUNSQR, EXP */

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

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

#define MX    20
#define MY    20
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

#define NEQ   (NS*MX*MY)
#define T0    ZERO
#define RTOL  RCONST(1.0e-5)
#define ATOL  RCONST(1.0e-5)

/* Output Constants */

#define TOUT RCONST(10.0)

/* Note: The value for species i at mesh point (j,k) is stored in */
/* component number (i-1) + j*NS + k*NS*MX of an N_Vector,        */
/* where 1 <= i <= NS, 0 <= j < MX, 0 <= k < MY.                  */

/* Structure for user data */

typedef struct {
  realtype **P[NGRP];
  sunindextype *pivot[NGRP];
  int ns, mxns, mp, mq, mx, my, ngrp, ngx, ngy, mxmp;
  int jgx[NGX+1], jgy[NGY+1], jigx[MX], jigy[MY];
  int jxr[NGX], jyr[NGY];
  realtype acoef[NS][NS], bcoef[NS], diff[NS];
  realtype cox[NS], coy[NS], dx, dy, srur;
  realtype fsave[NEQ];
  realtype fBsave[NEQ];
  N_Vector rewt;
  N_Vector vtemp;
  void *cvode_mem;
  int indexB;
} *WebData;

/* Adjoint calculation constants */
/* g = int_x int_y c(ISPEC) dy dx at t = Tfinal */
#define NSTEPS 80  /* check points every NSTEPS steps */
#define ISPEC  6   /* species # in objective */

/* Prototypes for user-supplied functions */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Precond(realtype t, N_Vector c, N_Vector fc,
                   booleantype jok, booleantype *jcurPtr, 
                   realtype gamma, void *user_data);

static int PSolve(realtype t, N_Vector c, N_Vector fc,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *user_data);

static int fB(realtype t, N_Vector c, N_Vector cB, 
               N_Vector cBdot, void *user_data);

static int PrecondB(realtype t, N_Vector c, 
                    N_Vector cB, N_Vector fcB, booleantype jok, 
                    booleantype *jcurPtr, realtype gamma,
                    void *user_data);

static int PSolveB(realtype t, N_Vector c, 
                   N_Vector cB, N_Vector fcB, 
                   N_Vector r, N_Vector z,
                   realtype gamma, realtype delta, 
                   int lr, void *user_data);

/* Prototypes for private functions */

static WebData AllocUserData(void);
static void InitUserData(WebData wdata);
static void SetGroups(int m, int ng, int jg[], int jig[], int jr[]);
static void CInit(N_Vector c, WebData wdata);
static void CbInit(N_Vector c, int is, WebData wdata);
static void PrintOutput(N_Vector c, int ns, int mxns, WebData wdata);
static void FreeUserData(WebData wdata);
static void WebRates(realtype x, realtype y, realtype t, realtype c[], realtype rate[],
                     WebData wdata);
static void WebRatesB(realtype x, realtype y, realtype t, realtype c[], realtype cB[], 
                      realtype rate[], realtype rateB[], WebData wdata);
static void fblock (realtype t, realtype cdata[], int jx, int jy, realtype cdotdata[],
                    WebData wdata);
static void GSIter(realtype gamma, N_Vector z, N_Vector x, WebData wdata);
static realtype doubleIntgr(N_Vector c, int i, WebData wdata);
static int check_retval(void *returnvalue, const char *funcname, int opt);

/* Small Vector Kernels */

static void v_inc_by_prod(realtype u[], realtype v[], realtype w[], int n);
static void v_sum_prods(realtype u[], realtype p[], realtype q[], realtype v[], 
                        realtype w[], int n);
static void v_prod(realtype u[], realtype v[], realtype w[], int n);
static void v_zero(realtype u[], int n);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  realtype abstol=ATOL, reltol=RTOL, t;
  N_Vector c;
  WebData wdata;
  void *cvode_mem;
  SUNLinearSolver LS, LSB;

  int retval, ncheck;

  int indexB;
  
  realtype reltolB=RTOL, abstolB=ATOL;
  N_Vector cB;

  c = cB = NULL;
  wdata = NULL;
  cvode_mem = NULL;
  LS = LSB = NULL;

  /* Allocate and initialize user data */

  wdata = AllocUserData();
  if(check_retval((void *)wdata, "AllocUserData", 2)) return(1);
  InitUserData(wdata);

  /* Set-up forward problem */

  /* Initializations */
  c = N_VNew_Serial(NEQ);
  if(check_retval((void *)c, "N_VNew_Serial", 0)) return(1);
  CInit(c, wdata);

  /* Call CVodeCreate/CVodeInit for forward run */
  printf("\nCreate and allocate CVODES memory for forward run\n");
  cvode_mem = CVodeCreate(CV_BDF);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  wdata->cvode_mem = cvode_mem; /* Used in Precond */
  retval = CVodeSetUserData(cvode_mem, wdata);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);
  retval = CVodeInit(cvode_mem, f, T0, c);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* Create SUNLinSol_SPGMR linear solver for forward run */
  LS = SUNLinSol_SPGMR(c, PREC_LEFT, 0);
  if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

  /* Attach the linear sovler */
  retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

  /* Set the preconditioner solve and setup functions */
  retval = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
  if(check_retval(&retval, "CVodeSetPreconditioner", 1)) return(1);

  /* Set-up adjoint calculations */

  printf("\nAllocate global memory\n");
  retval = CVodeAdjInit(cvode_mem, NSTEPS, CV_HERMITE);
  if(check_retval(&retval, "CVodeAdjInit", 1)) return(1);

  /* Perform forward run */

  printf("\nForward integration\n");
  retval = CVodeF(cvode_mem, TOUT, c, &t, CV_NORMAL, &ncheck);
  if(check_retval(&retval, "CVodeF", 1)) return(1);

  printf("\nncheck = %d\n", ncheck);


#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("\n   g = int_x int_y c%d(Tfinal,x,y) dx dy = %Lf \n\n", 
         ISPEC, doubleIntgr(c,ISPEC,wdata));
#else
  printf("\n   g = int_x int_y c%d(Tfinal,x,y) dx dy = %f \n\n", 
         ISPEC, doubleIntgr(c,ISPEC,wdata));
#endif

  /* Set-up backward problem */

  /* Allocate cB */
  cB = N_VNew_Serial(NEQ);
  if(check_retval((void *)cB, "N_VNew_Serial", 0)) return(1);
  /* Initialize cB = 0 */
  CbInit(cB, ISPEC, wdata);

  /* Create and allocate CVODES memory for backward run */
  printf("\nCreate and allocate CVODES memory for backward run\n");
  retval = CVodeCreateB(cvode_mem, CV_BDF, &indexB);
  if(check_retval(&retval, "CVodeCreateB", 1)) return(1);
  retval = CVodeSetUserDataB(cvode_mem, indexB, wdata);
  if(check_retval(&retval, "CVodeSetUserDataB", 1)) return(1);
  retval = CVodeInitB(cvode_mem, indexB, fB, TOUT, cB);
  if(check_retval(&retval, "CVodeInitB", 1)) return(1);
  retval = CVodeSStolerancesB(cvode_mem, indexB, reltolB, abstolB);
  if(check_retval(&retval, "CVodeSStolerancesB", 1)) return(1);

  wdata->indexB = indexB;

  /* Create SUNLinSol_SPGMR linear solver for backward run */
  LSB = SUNLinSol_SPGMR(cB, PREC_LEFT, 0);
  if(check_retval((void *)LSB, "SUNLinSol_SPGMR", 0)) return(1);

  /* Attach the linear sovler */
  retval = CVodeSetLinearSolverB(cvode_mem, indexB, LSB, NULL);
  if (check_retval(&retval, "CVodeSetLinearSolverB", 1)) return 1;

  /* Set the preconditioner solve and setup functions */
  retval = CVodeSetPreconditionerB(cvode_mem, indexB, PrecondB, PSolveB);
  if(check_retval(&retval, "CVodeSetPreconditionerB", 1)) return(1);

  /* Perform backward integration */

  printf("\nBackward integration\n");
  retval = CVodeB(cvode_mem, T0, CV_NORMAL);
  if(check_retval(&retval, "CVodeB", 1)) return(1);

  retval = CVodeGetB(cvode_mem, indexB, &t, cB);
  if(check_retval(&retval, "CVodeGetB", 1)) return(1);

  PrintOutput(cB, NS, MXNS, wdata);

  /* Free all memory */
  CVodeFree(&cvode_mem);

  N_VDestroy(c);
  N_VDestroy(cB);
  SUNLinSolFree(LS);
  SUNLinSolFree(LSB);

  FreeUserData(wdata);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * This routine computes the right-hand side of the ODE system and
 * returns it in cdot. The interaction rates are computed by calls to WebRates,
 * and these are saved in fsave for use in preconditioning.
 */

static int f(realtype t, N_Vector c, N_Vector cdot, void *user_data)
{
  int i, ic, ici, idxl, idxu, idyl, idyu, iyoff, jx, jy, ns, mxns;
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
        cdotdata[ici] = coy[i-1]*(dcyui - dcyli) + 
                        cox[i-1]*(dcxui - dcxli) +
                        fsave[ici];
      }
    }
  }

  return(0);
}

/*
 * This routine generates the block-diagonal part of the Jacobian
 * corresponding to the interaction rates, multiplies by -gamma, adds
 * the identity matrix, and calls denseGETRF to do the LU decomposition of
 * each diagonal block. The computation of the diagonal blocks uses
 * the preset block and grouping information. One block per group is
 * computed. The Jacobian elements are generated by difference
 * quotients using calls to the routine fblock.
 *
 * This routine can be regarded as a prototype for the general case
 * of a block-diagonal preconditioner. The blocks are of size mp, and
 * there are ngrp=ngx*ngy blocks computed in the block-grouping scheme.
 */
 
static int Precond(realtype t, N_Vector c, N_Vector fc,
                   booleantype jok, booleantype *jcurPtr, 
                   realtype gamma, void *user_data)
{
  realtype ***P;
  sunindextype **pivot;
  int i, if0, if00, ig, igx, igy, j, jj, jx, jy;
  int *jxr, *jyr, ngrp, ngx, ngy, mxmp, retval;
  sunindextype mp, denseretval;
  realtype uround, fac, r, r0, save, srur;
  realtype *f1, *fsave, *cdata, *rewtdata;
  WebData wdata;
  N_Vector rewt;

  wdata = (WebData) user_data;
  rewt = wdata->rewt;
  retval = CVodeGetErrWeights(wdata->cvode_mem, rewt);
  if(check_retval(&retval, "CVodeGetErrWeights", 1)) return(1);

  cdata = N_VGetArrayPointer(c);
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

  f1 = N_VGetArrayPointer(wdata->vtemp);

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
 * This routine applies two inverse preconditioner matrices
 * to the vector r, using the interaction-only block-diagonal Jacobian
 * with block-grouping, denoted Jr, and Gauss-Seidel applied to the
 * diffusion contribution to the Jacobian, denoted Jd.
 * It first calls GSIter for a Gauss-Seidel approximation to
 * ((I - gamma*Jd)-inverse)*r, and stores the result in z.
 * Then it computes ((I - gamma*Jr)-inverse)*z, using LU factors of the
 * blocks in P, and pivot information in pivot, and returns the result in z.
 */

static int PSolve(realtype t, N_Vector c, N_Vector fc,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *user_data)
{
  realtype   ***P;
  sunindextype **pivot;
  int jx, jy, igx, igy, iv, ig, *jigx, *jigy, mx, my, ngx, mp;
  WebData wdata;

  wdata = (WebData) user_data;

  N_VScale(ONE, r, z);

  /* call GSIter for Gauss-Seidel iterations */

  GSIter(gamma, z, wdata->vtemp, wdata);

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
 * This routine computes the right-hand side of the adjoint ODE system and
 * returns it in cBdot. The interaction rates are computed by calls to WebRates,
 * and these are saved in fsave for use in preconditioning. The adjoint 
 * interaction rates are computed by calls to WebRatesB.
 */

static int fB(realtype t, N_Vector c, N_Vector cB, 
              N_Vector cBdot, void *user_data)
{
  int i, ic, ici, idxl, idxu, idyl, idyu, iyoff, jx, jy, ns, mxns;
  realtype dcxli, dcxui, dcyli, dcyui, x, y, *cox, *coy, *fsave, *fBsave, dx, dy;
  realtype *cdata, *cBdata, *cBdotdata;
  WebData wdata;

  wdata = (WebData) user_data;
  cdata = N_VGetArrayPointer(c);
  cBdata = N_VGetArrayPointer(cB);
  cBdotdata = N_VGetArrayPointer(cBdot);

  mxns = wdata->mxns;
  ns = wdata->ns;
  fsave = wdata->fsave;
  fBsave = wdata->fBsave;
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
      WebRatesB(x, y, t, cdata+ic, cBdata+ic, fsave+ic, fBsave+ic, wdata);
      idxu = (jx == MX-1) ? -ns : ns;
      idxl = (jx == 0) ? -ns : ns;
      for (i = 1; i <= ns; i++) {
        ici = ic + i-1;
        /* Do differencing in y. */
        dcyli = cBdata[ici] - cBdata[ici-idyl];
        dcyui = cBdata[ici+idyu] - cBdata[ici];
        /* Do differencing in x. */
        dcxli = cBdata[ici] - cBdata[ici-idxl];
        dcxui = cBdata[ici+idxu] - cBdata[ici];
        /* Collect terms and load cdot elements. */
        cBdotdata[ici] = - coy[i-1]*(dcyui - dcyli) 
                         - cox[i-1]*(dcxui - dcxli)
	                 - fBsave[ici];
      }
    }
  }

  return(0);
}

/*
 * Preconditioner setup function for the backward problem
 */

static int PrecondB(realtype t, N_Vector c, 
                    N_Vector cB, N_Vector fcB, booleantype jok, 
                    booleantype *jcurPtr, realtype gamma,
                    void *user_data)
{
  realtype ***P;
  sunindextype **pivot;
  sunindextype denseretval;
  int i, if0, if00, ig, igx, igy, j, jj, jx, jy;
  int *jxr, *jyr, mp, ngrp, ngx, ngy, mxmp, retval;
  realtype uround, fac, r, r0, save, srur;
  realtype *f1, *fsave, *cdata, *rewtdata;
  void *cvode_mem;
  WebData wdata;
  N_Vector rewt;

  wdata = (WebData) user_data;
  cvode_mem = CVodeGetAdjCVodeBmem(wdata->cvode_mem, wdata->indexB);
  if(check_retval((void *)cvode_mem, "CVadjGetCVodeBmem", 0)) return(1);
  rewt = wdata->rewt;
  retval = CVodeGetErrWeights(cvode_mem, rewt);
  if(check_retval(&retval, "CVodeGetErrWeights", 1)) return(1);

  cdata = N_VGetArrayPointer(c);
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

  f1 = N_VGetArrayPointer(wdata->vtemp);

  fac = N_VWrmsNorm (fcB, rewt);
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
        fac = gamma/r;
        fblock (t, cdata, jx, jy, f1, wdata);
        for (i = 0; i < mp; i++) {
          P[ig][i][j] = (f1[i] - fsave[if0+i])*fac;
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
 * Preconditioner solve function for the backward problem
 */

static int PSolveB(realtype t, N_Vector c, 
                   N_Vector cB, N_Vector fcB, 
                   N_Vector r, N_Vector z,
                   realtype gamma, realtype delta, 
                   int lr, void *user_data)
{
  realtype ***P;
  sunindextype **pivot;
  int jx, jy, igx, igy, iv, ig, *jigx, *jigy, mx, my, ngx;
  sunindextype mp;
  WebData wdata;

  wdata = (WebData) user_data;

  N_VScale(ONE, r, z);

  /* call GSIter for Gauss-Seidel iterations (same routine but with gamma=-gamma) */

  GSIter(-gamma, z, wdata->vtemp, wdata);

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
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Allocate space for user data structure
 */

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
  wdata->rewt  = N_VNew_Serial(NEQ);
  wdata->vtemp = N_VNew_Serial(NEQ);

  return(wdata);
}

/*
 * Initialize user data structure
 */

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

  for (j = 0; j < NS; j++) { for (i = 0; i < NS; i++) acoef[i][j] = ZERO; }
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
 * This routine sets arrays jg, jig, and jr describing
 * a uniform partition of (0,1,2,...,m-1) into ng groups.
 * The arrays set are:
 *   jg    = length ng+1 array of group boundaries.
 *           Group ig has indices j = jg[ig],...,jg[ig+1]-1.
 *   jig   = length m array of group indices vs node index.
 *           Node index j is in group jig[j].
 *   jr    = length ng array of indices representing the groups.
 *           The index for group ig is j = jr[ig].
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

/*
 * This routine computes and loads the vector of initial values. 
 */

static void CInit(N_Vector c, WebData wdata)
{
  int i, ici, ioff, iyoff, jx, jy, ns, mxns;
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

/*
 * This function computes and loads the final values for the adjoint variables
 */

static void CbInit(N_Vector c, int is, WebData wdata)
{
  int i, ici, ioff, iyoff, jx, jy, ns, mxns;
  realtype *cdata;

  realtype gu[NS];

  cdata = N_VGetArrayPointer(c);
  ns = wdata->ns;
  mxns = wdata->mxns;

  for ( i = 1; i <= ns; i++ ) gu[i-1] = ZERO; 
  gu[ISPEC-1] = ONE;

  for (jy = 0; jy < MY; jy++) {
    iyoff = mxns*jy;
    for (jx = 0; jx < MX; jx++) {
      ioff = iyoff + ns*jx;
      for (i = 1; i <= ns; i++) {
        ici = ioff + i-1;
        cdata[ici] = gu[i-1];
      }
    }
  }
}

/*
 * This routine computes the interaction rates for the species
 * c_1, ... ,c_ns (stored in c[0],...,c[ns-1]), at one spatial point 
 * and at time t.
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
 * This routine computes the interaction rates for the backward problem
 */

static void WebRatesB(realtype x, realtype y, realtype t, realtype c[], realtype cB[], 
                      realtype rate[], realtype rateB[], WebData wdata)
{
  int i, j, ns;
  realtype fac, *bcoef;
  realtype (*acoef)[NS];

  ns = wdata->ns;
  acoef = wdata->acoef;
  bcoef = wdata->bcoef;

  fac = ONE + ALPH*x*y;

  for (i = 0; i < ns; i++)
    rate[i] = bcoef[i]*fac;
  
  for (j = 0; j < ns; j++) 
    for (i = 0; i < ns; i++)
      rate[i] += acoef[i][j]*c[j];

  for (i = 0; i < ns; i++) {
    rateB[i] = cB[i]*rate[i];
    rate[i] = c[i]*rate[i];
  }

  for (j = 0; j < ns; j++) 
    for (i = 0; i < ns; i++)
      rateB[i] += acoef[j][i]*c[j]*cB[j];

}

/*
 * This routine computes one block of the interaction terms of the
 * system, namely block (jx,jy), for use in preconditioning.
 * Here jx and jy count from 0.
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
 * This routine performs ITMAX=5 Gauss-Seidel iterations to compute an
 * approximation to (P-inverse)*z, where P = I - gamma*Jd, and
 * Jd represents the diffusion contributions to the Jacobian.
 * The answer is stored in z on return, and x is a temporary vector.
 * The dimensions below assume a global constant NS >= ns.
 * Some inner loops of length ns are implemented with the small
 * vector kernels v_sum_prods, v_prod, v_inc_by_prod.
 */

static void GSIter(realtype gamma, N_Vector z, N_Vector x, WebData wdata)
{
  int i, ic, iter, iyoff, jx, jy, ns, mxns, mx, my, x_loc, y_loc;
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
    temp = ONE/(ONE + TWO*gamma*(cox[i] + coy[i]));
    beta[i] = gamma*cox[i]*temp;
    beta2[i] = TWO*beta[i];
    gam[i] = gamma*coy[i]*temp;
    gam2[i] = TWO*gam[i];
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
          case 0 : /* jx == 0, jy == 0 */
            /* x[ic+i] = beta2[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] */
            v_sum_prods(xd+ic, beta2, xd+ic+ns, gam2, xd+ic+mxns, ns);
            break;
          case 1 : /* 1 <= jx <= mx-2, jy == 0 */
            /* x[ic+i] = beta[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] */
            v_sum_prods(xd+ic, beta, xd+ic+ns, gam2, xd+ic+mxns, ns);
            break;
          case 2 : /* jx == mx-1, jy == 0 */
            /* x[ic+i] = gam2[i]x[ic+mxns+i] */
            v_prod(xd+ic, gam2, xd+ic+mxns, ns);
            break;
          case 3 : /* jx == 0, 1 <= jy <= my-2 */
            /* x[ic+i] = beta2[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] */
            v_sum_prods(xd+ic, beta2, xd+ic+ns, gam, xd+ic+mxns, ns);
            break;
          case 4 : /* 1 <= jx <= mx-2, 1 <= jy <= my-2 */
            /* x[ic+i] = beta[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] */
            v_sum_prods(xd+ic, beta, xd+ic+ns, gam, xd+ic+mxns, ns);
            break;
          case 5 : /* jx == mx-1, 1 <= jy <= my-2 */
            /* x[ic+i] = gam[i]x[ic+mxns+i] */
            v_prod(xd+ic, gam, xd+ic+mxns, ns);
            break;
          case 6 : /* jx == 0, jy == my-1 */
            /* x[ic+i] = beta2[i]x[ic+ns+i] */
            v_prod(xd+ic, beta2, xd+ic+ns, ns);
            break;
          case 7 : /* 1 <= jx <= mx-2, jy == my-1 */
            /* x[ic+i] = beta[i]x[ic+ns+i] */
            v_prod(xd+ic, beta, xd+ic+ns, ns);
            break;
          case 8 : /* jx == mx-1, jy == my-1 */
            /* x[ic+i] = ZERO */
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
        case 0 : /* jx == 0, jy == 0 */
          break;
        case 1 : /* 1 <= jx <= mx-2, jy == 0 */
          /* x[ic+i] += beta[i]x[ic-ns+i] */
          v_inc_by_prod(xd+ic, beta, xd+ic-ns, ns);
          break;
        case 2 : /* jx == mx-1, jy == 0 */
          /* x[ic+i] += beta2[i]x[ic-ns+i] */
          v_inc_by_prod(xd+ic, beta2, xd+ic-ns, ns);
          break;
        case 3 : /* jx == 0, 1 <= jy <= my-2 */
          /* x[ic+i] += gam[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, gam, xd+ic-mxns, ns);
          break;
        case 4 : /* 1 <= jx <= mx-2, 1 <= jy <= my-2 */
          /* x[ic+i] += beta[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, beta, xd+ic-ns, ns);
          v_inc_by_prod(xd+ic, gam, xd+ic-mxns, ns);
          break;
        case 5 : /* jx == mx-1, 1 <= jy <= my-2 */
          /* x[ic+i] += beta2[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, beta2, xd+ic-ns, ns);
          v_inc_by_prod(xd+ic, gam, xd+ic-mxns, ns);
          break;
        case 6 : /* jx == 0, jy == my-1 */
          /* x[ic+i] += gam2[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, gam2, xd+ic-mxns, ns);
          break;
        case 7 : /* 1 <= jx <= mx-2, jy == my-1 */
          /* x[ic+i] += beta[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] */
          v_inc_by_prod(xd+ic, beta, xd+ic-ns, ns);
          v_inc_by_prod(xd+ic, gam2, xd+ic-mxns, ns);
          break;
        case 8 : /* jx == mx-1, jy == my-1 */
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

/*
 * Print maximum sensitivity of G for each species
 */

static void PrintOutput(N_Vector cB, int ns, int mxns, WebData wdata)
{
  int i, jx, jy;
  realtype *cdata, cij, cmax, x, y;

  x = y = ZERO;

  cdata = N_VGetArrayPointer(cB);

  for (i=1; i <= ns; i++) {

    cmax = ZERO;

    for (jy=MY-1; jy >= 0; jy--) {
      for (jx=0; jx < MX; jx++) {
        cij = cdata[(i-1) + jx*ns + jy*mxns];
        if (SUNRabs(cij) > cmax) {
          cmax = cij;
          x = jx * wdata->dx;
          y = jy * wdata->dy;
        }
      }
    }

    printf("\nMaximum sensitivity with respect to I.C. of species %d\n", i);
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("  mu max = %Le\n",cmax);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("  mu max = %e\n",cmax);
#else
    printf("  mu max = %e\n",cmax);
#endif
    printf("at\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("  x = %Le\n  y = %Le\n", x, y);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("  x = %e\n  y = %e\n", x, y);
#else
    printf("  x = %e\n  y = %e\n", x, y);
#endif

  }

}

/*
 * Compute double space integral
 */

static realtype doubleIntgr(N_Vector c, int i, WebData wdata)
{
  realtype *cdata;
  int ns, mx, my, mxns;
  realtype dx, dy;
  realtype intgr_xy, intgr_x;
  int jx, jy;

  cdata = N_VGetArrayPointer(c);

  ns   = wdata->ns;
  mx   = wdata->mx;
  my   = wdata->my;
  mxns = wdata->mxns;
  dx   = wdata->dx;
  dy   = wdata->dy;

  jy = 0;
  intgr_x = cdata[(i-1)+jy*mxns];
  for (jx = 1; jx < mx-1; jx++) {
    intgr_x += TWO*cdata[(i-1) + jx*ns + jy*mxns]; 
  }
  intgr_x += cdata[(i-1)+(mx-1)*ns+jy*mxns];
  intgr_x *= RCONST(0.5)*dx;
  
  intgr_xy = intgr_x;
  
  for (jy = 1; jy < my-1; jy++) {
    
    intgr_x = cdata[(i-1)+jy*mxns];
    for (jx = 1; jx < mx-1; jx++) {
      intgr_x += TWO*cdata[(i-1) + jx*ns + jy*mxns]; 
    }
    intgr_x += cdata[(i-1)+(mx-1)*ns+jy*mxns];
    intgr_x *= RCONST(0.5)*dx;
    
    intgr_xy += TWO*intgr_x;

  }
  
  jy = my-1;
  intgr_x = cdata[(i-1)+jy*mxns];
  for (jx = 1; jx < mx-1; jx++) {
    intgr_x += TWO*cdata[(i-1) + jx*ns + jy*mxns]; 
  }
  intgr_x += cdata[(i-1)+(mx-1)*ns+jy*mxns];
  intgr_x *= RCONST(0.5)*dx;
  
  intgr_xy += intgr_x;
  
  intgr_xy *= RCONST(0.5)*dy;

  return(intgr_xy);
}

/*
 * Free space allocated for the user data structure
 */

static void FreeUserData(WebData wdata)
{
  int i, ngrp;

  ngrp = wdata->ngrp;
  for(i=0; i < ngrp; i++) {
    destroyMat((wdata->P)[i]);
    destroyArray((wdata->pivot)[i]);
  }
  N_VDestroy(wdata->rewt);
  N_VDestroy(wdata->vtemp);
  free(wdata);
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
