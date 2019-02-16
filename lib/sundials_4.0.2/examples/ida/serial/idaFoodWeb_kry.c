/* -----------------------------------------------------------------
 * Programmer(s): Ting Yan @ SMU
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
 * Example program for IDA: Food web problem, OpenMP, GMRES, 
 * user-supplied preconditioner
 *
 * This example program uses the SPGMR as the linear
 * solver, and IDACalcIC for initial condition calculation.
 *
 * The mathematical problem solved in this example is a DAE system
 * that arises from a system of partial differential equations after
 * spatial discretization. The PDE system is a food web population
 * model, with predator-prey interaction and diffusion on the unit
 * square in two dimensions. The dependent variable vector is:
 *
 *         1   2         ns
 *   c = (c , c ,  ..., c  ) , ns = 2 * np
 *
 * and the PDE's are as follows:
 *
 *     i             i      i
 *   dc /dt = d(i)*(c    + c  )  +  R (x,y,c)   (i = 1,...,np)
 *                   xx     yy       i
 *
 *              i      i
 *   0 = d(i)*(c    + c  )  +  R (x,y,c)   (i = np+1,...,ns)
 *              xx     yy       i
 *
 *   where the reaction terms R are:
 *
 *                   i             ns         j
 *   R  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
 *    i                           j=1
 *
 * The number of species is ns = 2 * np, with the first np being
 * prey and the last np being predators. The coefficients a(i,j),
 * b(i), d(i) are:
 *
 *  a(i,i) = -AA   (all i)
 *  a(i,j) = -GG   (i <= np , j >  np)
 *  a(i,j) =  EE   (i >  np, j <= np)
 *  all other a(i,j) = 0
 *  b(i) = BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y)) (i <= np)
 *  b(i) =-BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y)) (i  > np)
 *  d(i) = DPREY   (i <= np)
 *  d(i) = DPRED   (i > np)
 *
 * The various scalar parameters required are set using '#define'
 * statements or directly in routine InitUserData. In this program,
 * np = 1, ns = 2. The boundary conditions are homogeneous Neumann:
 * normal derivative = 0.
 *
 * A polynomial in x and y is used to set the initial values of the
 * first np variables (the prey variables) at each x,y location,
 * while initial values for the remaining (predator) variables are
 * set to a flat value, which is corrected by IDACalcIC.
 *
 * The PDEs are discretized by central differencing on a MX by MY
 * mesh.
 *
 * The DAE system is solved by IDA using the SPGMR linear solver.
 * Output is printed at t = 0, .001, .01, .1, .4, .7, 1.
 * -----------------------------------------------------------------
 * References:
 * [1] Peter N. Brown and Alan C. Hindmarsh,
 *     Reduced Storage Matrix Methods in Stiff ODE systems, Journal
 *     of Applied Mathematics and Computation, Vol. 31 (May 1989),
 *     pp. 40-91.
 *
 * [2] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
 *     Using Krylov Methods in the Solution of Large-Scale
 *     Differential-Algebraic Systems, SIAM J. Sci. Comput., 15
 *     (1994), pp. 1467-1488.
 *
 * [3] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
 *     Consistent Initial Condition Calculation for Differential-
 *     Algebraic Systems, SIAM J. Sci. Comput., 19 (1998),
 *     pp. 1495-1512.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>                   /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to spgmr SUNLinearSolver      */
#include <sundials/sundials_dense.h>   /* use generic dense solver in precond. */
#include <sundials/sundials_types.h>   /* definition of type realtype          */
#include <sundials/sundials_math.h>    /* macros SUNRabs, SUNRsqrt, etc.       */

/* Problem Constants. */

#define NPREY       1              /* No. of prey (= no. of predators). */
#define NUM_SPECIES 2*NPREY

#define PI          RCONST(3.1415926535898)
#define FOURPI      (RCONST(4.0)*PI)

#define MX          20             /* MX = number of x mesh points      */
#define MY          20             /* MY = number of y mesh points      */
#define NSMX        (NUM_SPECIES * MX)
#define NEQ         (NUM_SPECIES*MX*MY)
#define AA          RCONST(1.0)    /* Coefficient in above eqns. for a  */
#define EE          RCONST(10000.) /* Coefficient in above eqns. for a  */
#define GG          RCONST(0.5e-6) /* Coefficient in above eqns. for a  */
#define BB          RCONST(1.0)    /* Coefficient in above eqns. for b  */
#define DPREY       RCONST(1.0)    /* Coefficient in above eqns. for d  */
#define DPRED       RCONST(0.05)   /* Coefficient in above eqns. for d  */
#define ALPHA       RCONST(50.)    /* Coefficient alpha in above eqns.  */
#define BETA        RCONST(1000.)  /* Coefficient beta in above eqns.   */
#define AX          RCONST(1.0)    /* Total range of x variable         */
#define AY          RCONST(1.0)    /* Total range of y variable         */
#define RTOL        RCONST(1.e-5)  /* Relative tolerance                */
#define ATOL        RCONST(1.e-5)  /* Absolute tolerance                */
#define NOUT        6              /* Number of output times            */
#define TMULT       RCONST(10.0)   /* Multiplier for tout values        */
#define TADD        RCONST(0.3)    /* Increment for tout values         */
#define ZERO        RCONST(0.)
#define ONE         RCONST(1.0)

/*
 * User-defined vector and accessor macro: IJ_Vptr.
 * IJ_Vptr is defined in order to express the underlying 3-D structure of
 * the dependent variable vector from its underlying 1-D storage (an N_Vector).
 * IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to
 * species index is = 0, x-index ix = i, and y-index jy = j.
 */

#define IJ_Vptr(vv,i,j) (&NV_Ith_S(vv, (i)*NUM_SPECIES + (j)*NSMX))

/* Type: UserData.  Contains problem constants, etc. */

typedef struct {
  sunindextype Neq, ns, np, mx, my;
  realtype dx, dy, **acoef;
  realtype cox[NUM_SPECIES], coy[NUM_SPECIES], bcoef[NUM_SPECIES];
  realtype **PP[MX][MY];
  sunindextype *pivot[MX][MY];
  N_Vector rates;
  N_Vector ewt;
  void *ida_mem;
} *UserData;

/* Prototypes for functions called by the IDA Solver. */

static int resweb(realtype time, N_Vector cc, N_Vector cp, N_Vector resval,
                  void *user_data);

static int Precond(realtype tt,
		   N_Vector cc, N_Vector cp, N_Vector rr,
		   realtype cj, void *user_data);

static int PSolve(realtype tt,
		  N_Vector cc, N_Vector cp, N_Vector rr,
		  N_Vector rvec, N_Vector zvec,
		  realtype cj, realtype delta, void *user_data);

/* Prototypes for private Helper Functions. */

static void InitUserData(UserData webdata);
static void SetInitialProfiles(N_Vector cc, N_Vector cp, N_Vector id,
                               UserData webdata);
static void PrintHeader(sunindextype maxl, realtype rtol, realtype atol);
static void PrintOutput(void *mem, N_Vector c, realtype t);
static void PrintFinalStats(void *mem);
static void Fweb(realtype tcalc, N_Vector cc, N_Vector crate, UserData webdata);
static void WebRates(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy,
                     UserData webdata);
static realtype dotprod(sunindextype size, realtype *x1, realtype *x2);
static int check_retval(void *returnvalue, char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main()
{
  void *mem;
  UserData webdata;
  N_Vector cc, cp, id;
  int iout, jx, jy, retval;
  sunindextype maxl;
  realtype rtol, atol, t0, tout, tret;
  SUNLinearSolver LS;

  mem = NULL;
  webdata = NULL;
  cc = cp = id = NULL;
  LS = NULL;

  /* Allocate and initialize user data block webdata. */

  webdata = (UserData) malloc(sizeof *webdata);
  webdata->rates = N_VNew_Serial(NEQ);
  webdata->acoef = newDenseMat(NUM_SPECIES, NUM_SPECIES);
  webdata->ewt = N_VNew_Serial(NEQ);
  for (jx = 0; jx < MX; jx++) {
    for (jy = 0; jy < MY; jy++) {
      (webdata->pivot)[jx][jy] = newIndexArray(NUM_SPECIES);
      (webdata->PP)[jx][jy] = newDenseMat(NUM_SPECIES, NUM_SPECIES);
    }
  }

  InitUserData(webdata);

  /* Allocate N-vectors and initialize cc, cp, and id. */

  cc  = N_VNew_Serial(NEQ);
  if(check_retval((void *)cc, "N_VNew_Serial", 0)) return(1);

  cp  = N_VNew_Serial(NEQ);
  if(check_retval((void *)cp, "N_VNew_Serial", 0)) return(1);

  id  = N_VNew_Serial(NEQ);
  if(check_retval((void *)id, "N_VNew_Serial", 0)) return(1);

  SetInitialProfiles(cc, cp, id, webdata);

  /* Set remaining inputs to IDAMalloc. */

  t0 = ZERO;
  rtol = RTOL;
  atol = ATOL;

  /* Call IDACreate and IDAMalloc to initialize IDA. */

  mem = IDACreate();
  if(check_retval((void *)mem, "IDACreate", 0)) return(1);

  retval = IDASetUserData(mem, webdata);
  if(check_retval(&retval, "IDASetUserData", 1)) return(1);

  retval = IDASetId(mem, id);
  if(check_retval(&retval, "IDASetId", 1)) return(1);

  retval = IDAInit(mem, resweb, t0, cc, cp);
  if(check_retval(&retval, "IDAInit", 1)) return(1);

  retval = IDASStolerances(mem, rtol, atol);
  if(check_retval(&retval, "IDASStolerances", 1)) return(1);

  webdata->ida_mem = mem;

  /* Create the linear solver SUNLinSol_SPGMR with left preconditioning
     and maximum Krylov dimension maxl */
  maxl = 16;
  LS = SUNLinSol_SPGMR(cc, PREC_LEFT, maxl);
  if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

  /* IDA recommends allowing up to 5 restarts (default is 0) */
  retval = SUNLinSol_SPGMRSetMaxRestarts(LS, 5);
  if(check_retval(&retval, "SUNLinSol_SPGMRSetMaxRestarts", 1)) return(1);

  /* Attach the linear sovler */
  retval = IDASetLinearSolver(mem, LS, NULL);
  if(check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

  /* Set the preconditioner solve and setup functions */
  retval = IDASetPreconditioner(mem, Precond, PSolve);
  if(check_retval(&retval, "IDASetPreconditioner", 1)) return(1);

  /* Call IDACalcIC (with default options) to correct the initial values. */

  tout = RCONST(0.001);
  retval = IDACalcIC(mem, IDA_YA_YDP_INIT, tout);
  if(check_retval(&retval, "IDACalcIC", 1)) return(1);

  /* Print heading, basic parameters, and initial values. */

  PrintHeader(maxl, rtol, atol);
  PrintOutput(mem, cc, ZERO);

  /* Loop over iout, call IDASolve (normal mode), print selected output. */

  for (iout = 1; iout <= NOUT; iout++) {

    retval = IDASolve(mem, tout, &tret, cc, cp, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1)) return(retval);

    PrintOutput(mem, cc, tret);

    if (iout < 3) tout *= TMULT; else tout += TADD;

  }

  /* Print final statistics and free memory. */

  PrintFinalStats(mem);

  /* Free memory */

  IDAFree(&mem);
  SUNLinSolFree(LS);

  N_VDestroy(cc);
  N_VDestroy(cp);
  N_VDestroy(id);


  destroyMat(webdata->acoef);
  N_VDestroy(webdata->rates);
  N_VDestroy(webdata->ewt);
  for (jx = 0; jx < MX; jx++) {
    for (jy = 0; jy < MY; jy ++) {
      destroyArray((webdata->pivot)[jx][jy]);
      destroyMat((webdata->PP)[jx][jy]);
    }
  }
  free(webdata);

  return(0);
}

/* Define lines for readability in later routines */

#define acoef  (webdata->acoef)
#define bcoef  (webdata->bcoef)
#define cox    (webdata->cox)
#define coy    (webdata->coy)

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * resweb: System residual function for predator-prey system.
 * This routine calls Fweb to get all the right-hand sides of the
 * equations, then loads the residual vector accordingly,
 * using cp in the case of prey species.
 */

static int resweb(realtype tt, N_Vector cc, N_Vector cp,
                  N_Vector res,  void *user_data)
{
  sunindextype jx, jy, is, yloc, loc, np;
  realtype *resv, *cpv;
  UserData webdata;

  webdata = (UserData)user_data;

  cpv = N_VGetArrayPointer(cp);
  resv = N_VGetArrayPointer(res);;
  np = webdata->np;

  /* Call Fweb to set res to vector of right-hand sides. */
  Fweb(tt, cc, res, webdata);

  /* Loop over all grid points, setting residual values appropriately
     for differential or algebraic components.                        */

  for (jy = 0; jy < MY; jy++) {
    yloc = NSMX * jy;
    for (jx = 0; jx < MX; jx++) {
      loc = yloc + NUM_SPECIES * jx;
      for (is = 0; is < NUM_SPECIES; is++) {
        if (is < np)
          resv[loc+is] = cpv[loc+is] - resv[loc+is];
        else
          resv[loc+is] = -resv[loc+is];
      }
    }
  }

  return(0);

}


static int Precond(realtype tt,
		   N_Vector cc, N_Vector cp, N_Vector rr,
		   realtype cj, void *user_data)
{
  int retval;
  realtype uround, xx, yy, del_x, del_y;
  realtype **Pxy, *ratesxy, *Pxycol, *cxy, *cpxy, *ewtxy, cctmp;
  realtype inc, fac, sqru, perturb_rates[NUM_SPECIES];
  int is, js, jx, jy, ret;
  void *mem;
  N_Vector ewt;
  realtype hh;
  UserData webdata;

  webdata = (UserData) user_data;
  del_x = webdata->dx;
  del_y = webdata->dy;

  uround = UNIT_ROUNDOFF;
  sqru = SUNRsqrt(uround);

  mem = webdata->ida_mem;
  ewt = webdata->ewt;
  retval = IDAGetErrWeights(mem, ewt);
  if(check_retval(&retval, "IDAGetErrWeights", 1)) return(1);
  retval = IDAGetCurrentStep(mem, &hh);
  if(check_retval(&retval, "IDAGetCurrentStep", 1)) return(1);

  for (jy = 0; jy < MY; jy++) {
    yy = jy * del_y;

    for (jx = 0; jx < MX; jx++) {
      xx = jx * del_x;
      Pxy = (webdata->PP)[jx][jy];
      cxy = IJ_Vptr(cc, jx, jy);
      cpxy = IJ_Vptr(cp, jx, jy);
      ewtxy = IJ_Vptr(ewt, jx, jy);
      ratesxy = IJ_Vptr((webdata->rates), jx, jy);

      for (js = 0; js < NUM_SPECIES; js++) {
	inc = sqru*(SUNMAX(SUNRabs(cxy[js]), SUNMAX(hh*SUNRabs(cpxy[js]), ONE/ewtxy[js])));
	cctmp = cxy[js];
	cxy[js] += inc;
	fac = -ONE/inc;

	WebRates(xx, yy, cxy, perturb_rates, webdata);

	Pxycol = Pxy[js];

	for (is = 0; is < NUM_SPECIES; is++)
	  Pxycol[is] = (perturb_rates[is] - ratesxy[is])*fac;

	if (js < 1) Pxycol[js] += cj;

	cxy[js] = cctmp;
      }

      ret = denseGETRF(Pxy, NUM_SPECIES, NUM_SPECIES, (webdata->pivot)[jx][jy]);

      if (ret != 0) return(1);
    }
  }

  return(0);

}


static int PSolve(realtype tt,
		  N_Vector cc, N_Vector cp, N_Vector rr,
		  N_Vector rvec, N_Vector zvec,
		  realtype cj, realtype dalta,
		  void *user_data)
{
  realtype **Pxy, *zxy;
  sunindextype *pivot;
  int jx, jy;
  UserData webdata;

  webdata = (UserData) user_data;

  N_VScale(ONE, rvec, zvec);

  for (jx = 0; jx < MX; jx++) {
    for (jy = 0; jy <MY; jy++) {

      zxy = IJ_Vptr(zvec, jx, jy);
      Pxy = (webdata->PP)[jx][jy];
      pivot = (webdata->pivot)[jx][jy];
      denseGETRS(Pxy, NUM_SPECIES, pivot, zxy);
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
 * InitUserData: Load problem constants in webdata (of type UserData).
 */

static void InitUserData(UserData webdata)
{
  int i, j, np;
  realtype *a1,*a2, *a3, *a4, dx2, dy2;

  webdata->mx = MX;
  webdata->my = MY;
  webdata->ns = NUM_SPECIES;
  webdata->np = NPREY;
  webdata->dx = AX/(MX-1);
  webdata->dy = AY/(MY-1);
  webdata->Neq= NEQ;

  /* Set up the coefficients a and b, and others found in the equations. */
  np = webdata->np;
  dx2 = (webdata->dx)*(webdata->dx); dy2 = (webdata->dy)*(webdata->dy);

  for (i = 0; i < np; i++) {
    a1 = &(acoef[i][np]);
    a2 = &(acoef[i+np][0]);
    a3 = &(acoef[i][0]);
    a4 = &(acoef[i+np][np]);
    /*  Fill in the portion of acoef in the four quadrants, row by row. */
    for (j = 0; j < np; j++) {
      *a1++ =  -GG;
      *a2++ =   EE;
      *a3++ = ZERO;
      *a4++ = ZERO;
    }

    /* Reset the diagonal elements of acoef to -AA. */
    acoef[i][i] = -AA; acoef[i+np][i+np] = -AA;

    /* Set coefficients for b and diffusion terms. */
    bcoef[i] = BB; bcoef[i+np] = -BB;
    cox[i] = DPREY/dx2; cox[i+np] = DPRED/dx2;
    coy[i] = DPREY/dy2; coy[i+np] = DPRED/dy2;
  }

}

/*
 * SetInitialProfiles: Set initial conditions in cc, cp, and id.
 * A polynomial profile is used for the prey cc values, and a constant
 * (1.0e5) is loaded as the initial guess for the predator cc values.
 * The id values are set to 1 for the prey and 0 for the predators.
 * The prey cp values are set according to the given system, and
 * the predator cp values are set to zero.
 */

static void SetInitialProfiles(N_Vector cc, N_Vector cp, N_Vector id,
                               UserData webdata)
{
  sunindextype loc, yloc, is, jx, jy, np;
  realtype xx, yy, xyfactor;
  realtype *ccv, *cpv, *idv;

  ccv = N_VGetArrayPointer(cc);
  cpv = N_VGetArrayPointer(cp) ;
  idv = N_VGetArrayPointer(id);
  np = webdata->np;

  /* Loop over grid, load cc values and id values. */
  for (jy = 0; jy < MY; jy++) {
    yy = jy * webdata->dy;
    yloc = NSMX * jy;
    for (jx = 0; jx < MX; jx++) {
      xx = jx * webdata->dx;
      xyfactor = RCONST(16.0)*xx*(ONE-xx)*yy*(ONE-yy);
      xyfactor *= xyfactor;
      loc = yloc + NUM_SPECIES*jx;

      for (is = 0; is < NUM_SPECIES; is++) {
        if (is < np) {
          ccv[loc+is] = RCONST(10.0) + (realtype)(is+1) * xyfactor;
          idv[loc+is] = ONE;
        }
        else {
	  ccv[loc+is] = RCONST(1.0e5);
          idv[loc+is] = ZERO;
        }
      }
    }
  }

  /* Set c' for the prey by calling the function Fweb. */
  Fweb(ZERO, cc, cp, webdata);

  /* Set c' for predators to 0. */
  for (jy = 0; jy < MY; jy++) {
    yloc = NSMX * jy;
    for (jx = 0; jx < MX; jx++) {
      loc = yloc + NUM_SPECIES * jx;
      for (is = np; is < NUM_SPECIES; is++) {
        cpv[loc+is] = ZERO;
      }
    }
  }
}

/*
 * Print first lines of output (problem description)
 */

static void PrintHeader(sunindextype maxl, realtype rtol, realtype atol)
{
  printf("\nidaFoodWeb_kry: Predator-prey DAE serial example problem for IDA \n\n");
  printf("Number of species ns: %d", NUM_SPECIES);
  printf("     Mesh dimensions: %d x %d", MX, MY);
  printf("     System size: %d\n", NEQ);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Linear solver: SPGMR,  SPGMR parameters maxl = %ld\n",(long int) maxl);
  printf("CalcIC called to correct initial predator concentrations.\n\n");
  printf("-----------------------------------------------------------\n");
  printf("  t        bottom-left  top-right");
  printf("    | nst  k      h\n");
  printf("-----------------------------------------------------------\n\n");

}

/*
 * PrintOutput: Print output values at output time t = tt.
 * Selected run statistics are printed.  Then values of the concentrations
 * are printed for the bottom left and top right grid points only.
 */

static void PrintOutput(void *mem, N_Vector c, realtype t)
{
  int i, kused, retval;
  long int nst;
  realtype *c_bl, *c_tr, hused;

  retval = IDAGetLastOrder(mem, &kused);
  check_retval(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_retval(&retval, "IDAGetLastStep", 1);

  c_bl = IJ_Vptr(c,0,0);
  c_tr = IJ_Vptr(c,MX-1,MY-1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.2Le %12.4Le %12.4Le   | %3ld  %1d %12.4Le\n",
         t, c_bl[0], c_tr[0], nst, kused, hused);
  for (i=1;i<NUM_SPECIES;i++)
    printf("         %12.4Le %12.4Le   |\n",c_bl[i],c_tr[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%8.2e %12.4e %12.4e   | %3ld  %1d %12.4e\n",
         t, c_bl[0], c_tr[0], nst, kused, hused);
  for (i=1;i<NUM_SPECIES;i++)
    printf("         %12.4e %12.4e   |\n",c_bl[i],c_tr[i]);
#else
  printf("%8.2e %12.4e %12.4e   | %3ld  %1d %12.4e\n",
         t, c_bl[0], c_tr[0], nst, kused, hused);
  for (i=1;i<NUM_SPECIES;i++)
    printf("         %12.4e %12.4e   |\n",c_bl[i],c_tr[i]);
#endif

  printf("\n");
}

/*
 * PrintFinalStats: Print final run data contained in iopt.
 */

static void PrintFinalStats(void *mem)
{
  long int nst, nre, sli, netf, nps, npevals, nrevalsLS;
  int retval;

  retval = IDAGetNumSteps(mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumLinIters(mem, &sli);
  check_retval(&retval, "IDAGetNumLinIters", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_retval(&retval, "IDAGetNumResEvals", 1);
  retval = IDAGetNumErrTestFails(mem, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumPrecSolves(mem, &nps);
  check_retval(&retval, "IDAGetNumPrecSolves", 1);
  retval = IDAGetNumPrecEvals(mem, &npevals);
  check_retval(&retval, "IDAGetNumPrecEvals", 1);
  retval = IDAGetNumLinResEvals(mem, &nrevalsLS);
  check_retval(&retval, "IDAGetNumLinResEvals", 1);

  printf("-----------------------------------------------------------\n");
  printf("Final run statistics: \n\n");
  printf("Number of steps                       = %ld\n", nst);
  printf("Number of residual evaluations        = %ld\n", nre);
  printf("Number of Preconditioner evaluations  = %ld\n", npevals);
  printf("Number of linear iterations           = %ld\n", sli);
  printf("Number of error test failures         = %ld\n", netf);
  printf("Number of precond solve fun called    = %ld\n", nps);

}

/*
 * Fweb: Rate function for the food-web problem.
 * This routine computes the right-hand sides of the system equations,
 * consisting of the diffusion term and interaction term.
 * The interaction term is computed by the function WebRates.
 */

static void Fweb(realtype tcalc, N_Vector cc, N_Vector crate,
                 UserData webdata)
{
  sunindextype jx, jy, is, idyu, idyl, idxu, idxl;
  realtype xx, yy, *cxy, *ratesxy, *cratexy, dcyli, dcyui, dcxli, dcxui;

  /* Loop over grid points, evaluate interaction vector (length ns),
     form diffusion difference terms, and load crate.                    */

  for (jy = 0; jy < MY; jy++) {
    yy = (webdata->dy) * jy ;
    idyu = (jy!=MY-1) ? NSMX : -NSMX;
    idyl = (jy!= 0  ) ? NSMX : -NSMX;

    for (jx = 0; jx < MX; jx++) {
      xx = (webdata->dx) * jx;
      idxu = (jx!= MX-1) ?  NUM_SPECIES : -NUM_SPECIES;
      idxl = (jx!=  0  ) ?  NUM_SPECIES : -NUM_SPECIES;
      cxy = IJ_Vptr(cc,jx,jy);
      ratesxy = IJ_Vptr(webdata->rates,jx,jy);
      cratexy = IJ_Vptr(crate,jx,jy);

      /* Get interaction vector at this grid point. */
      WebRates(xx, yy, cxy, ratesxy, webdata);

      /* Loop over species, do differencing, load crate segment. */
      for (is = 0; is < NUM_SPECIES; is++) {

        /* Differencing in y. */
        dcyli = *(cxy+is) - *(cxy - idyl + is) ;
        dcyui = *(cxy + idyu + is) - *(cxy+is);

        /* Differencing in x. */
        dcxli = *(cxy+is) - *(cxy - idxl + is);
        dcxui = *(cxy + idxu +is) - *(cxy+is);

        /* Compute the crate values at (xx,yy). */
        cratexy[is] = coy[is] * (dcyui - dcyli) +
          cox[is] * (dcxui - dcxli) + ratesxy[is];

      } /* End is loop */
    } /* End of jx loop */
  } /* End of jy loop */

}

/*
 * WebRates: Evaluate reaction rates at a given spatial point.
 * At a given (x,y), evaluate the array of ns reaction terms R.
 */

static void WebRates(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy,
                     UserData webdata)
{
  int is;
  realtype fac;

  for (is = 0; is < NUM_SPECIES; is++)
    ratesxy[is] = dotprod(NUM_SPECIES, cxy, acoef[is]);

  fac = ONE + ALPHA*xx*yy + BETA*sin(FOURPI*xx)*sin(FOURPI*yy);

  for (is = 0; is < NUM_SPECIES; is++)
    ratesxy[is] = cxy[is]*( bcoef[is]*fac + ratesxy[is] );

}

/*
 * dotprod: dot product routine for realtype arrays, for use by WebRates.
 */

static realtype dotprod(sunindextype size, realtype *x1, realtype *x2)
{
  sunindextype i;
  realtype *xx1, *xx2, temp = ZERO;

  xx1 = x1; xx2 = x2;
  for (i = 0; i < size; i++) temp += (*xx1++) * (*xx2++);
  return(temp);

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

static int check_retval(void *returnvalue, char *funcname, int opt)
{
  int *retval;

  if (opt == 0 && returnvalue == NULL) {
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
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
