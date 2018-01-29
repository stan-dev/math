/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds and Ting Yan @ SMU
 * -----------------------------------------------------------------
 * Example program for IDAS: Food web problem, OpenMP, GMRES, 
 * user-supplied preconditioner
 *
 * This example program uses SUNSPGMR as the linear 
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
 * The DAE system is solved by IDAS using the SUNSPGMR linear solver.
 * Output is printed at t = 0, .001, .01, .1, .4, .7, 1.
 *
 * Optionally, we can set the number of threads from environment 
 * variable or command line. To check the current value for number
 * of threads from environment:
 *      % echo $OMP_NUM_THREADS
 *
 * Execution:
 *
 * To use the default value for the number of threads from 
 * the OMP_NUM_THREADS environment value:
 *      % ./idasFoodWeb_kry_omp 
 * To specify the number of threads at the command line, use
 *      % ./idasFoodWeb_kry_omp num_threads
 * where num_threads is the desired number of threads. 
 *
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
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <idas/idas.h>
#include <idas/idas_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <nvector/nvector_openmp.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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

#define IJ_Vptr(vv,i,j) (&NV_Ith_OMP(vv, (i)*NUM_SPECIES + (j)*NSMX))

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
  int nthreads;
} *UserData;

/* Prototypes for functions called by the IDA Solver. */

static int resweb(realtype time, N_Vector cc, N_Vector cp, N_Vector resval, 
                  void *user_data);

static int Precond(realtype tt, N_Vector cc, N_Vector cp,
                   N_Vector rr, realtype cj, void *user_data);

static int PSolve(realtype tt, N_Vector cc, N_Vector cp,
                  N_Vector rr, N_Vector rvec, N_Vector zvec,
		  realtype cj, realtype delta, void *user_data);

/* Prototypes for private Helper Functions. */

static void InitUserData(UserData webdata);
static void SetInitialProfiles(N_Vector cc, N_Vector cp, N_Vector id,
                               UserData webdata);
static void PrintHeader(sunindextype maxl, realtype rtol, realtype atol);
static void PrintOutput(void *ida_mem, N_Vector c, realtype t);
static void PrintFinalStats(void *ida_mem);
static void Fweb(realtype tcalc, N_Vector cc, N_Vector crate, UserData webdata);
static void WebRates(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy, 
                     UserData webdata);
static realtype dotprod(sunindextype size, realtype *x1, realtype *x2);
static int check_flag(void *flagvalue, char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{ 
  void *ida_mem;
  SUNLinearSolver LS;
  UserData webdata;
  N_Vector cc, cp, id;
  int iout, jx, jy, flag;
  sunindextype maxl;
  realtype rtol, atol, t0, tout, tret;
  int num_threads;

  ida_mem = NULL;
  LS = NULL;
  webdata = NULL;
  cc = cp = id = NULL;

  /* Set the number of threads to use */
  num_threads = 1;   /* default value */
#ifdef _OPENMP
  num_threads = omp_get_max_threads();    /* overwrite with OMP_NUM_THREADS */
#endif
  if (argc > 1)
    num_threads = strtol(argv[1], NULL, 0);

  /* Allocate and initialize user data block webdata. */

  webdata = (UserData) malloc(sizeof *webdata);
  webdata->rates = N_VNew_OpenMP(NEQ, num_threads);
  webdata->acoef = newDenseMat(NUM_SPECIES, NUM_SPECIES);
  webdata->ewt = N_VNew_OpenMP(NEQ, num_threads);
  for (jx = 0; jx < MX; jx++) {
    for (jy = 0; jy < MY; jy++) {
      (webdata->pivot)[jx][jy] = newIndexArray(NUM_SPECIES);
      (webdata->PP)[jx][jy] = newDenseMat(NUM_SPECIES, NUM_SPECIES);
    }
  }
  webdata->nthreads = num_threads;

  InitUserData(webdata);

  /* Allocate N-vectors and initialize cc, cp, and id. */

  cc  = N_VNew_OpenMP(NEQ, num_threads);
  if(check_flag((void *)cc, "N_VNew_OpenMP", 0)) return(1);

  cp  = N_VNew_OpenMP(NEQ, num_threads);
  if(check_flag((void *)cp, "N_VNew_OpenMP", 0)) return(1);

  id  = N_VNew_OpenMP(NEQ, num_threads);
  if(check_flag((void *)id, "N_VNew_OpenMP", 0)) return(1);
  
  SetInitialProfiles(cc, cp, id, webdata);
  
  /* Set remaining inputs to IDAMalloc. */
  
  t0 = ZERO;
  rtol = RTOL; 
  atol = ATOL;

  /* Call IDACreate and IDAMalloc to initialize IDA. */
  
  ida_mem = IDACreate();
  if(check_flag((void *)ida_mem, "IDACreate", 0)) return(1);

  flag = IDASetUserData(ida_mem, webdata);
  if(check_flag(&flag, "IDASetUserData", 1)) return(1);

  flag = IDASetId(ida_mem, id);
  if(check_flag(&flag, "IDASetId", 1)) return(1);

  flag = IDAInit(ida_mem, resweb, t0, cc, cp);
  if(check_flag(&flag, "IDAInit", 1)) return(1);

  flag = IDASStolerances(ida_mem, rtol, atol);
  if(check_flag(&flag, "IDASStolerances", 1)) return(1);

  webdata->ida_mem = ida_mem;

  /* Create SUNSPGMR linear solver, attach to IDA, and set 
     preconditioning routines. */

  maxl = 16;                               /* max dimension of the Krylov subspace */
  LS = SUNSPGMR(cc, PREC_LEFT, maxl);      /* IDA only allows left preconditioning */
  if(check_flag((void *)LS, "SUNSPGMR", 0)) return(1);

  flag = IDASpilsSetLinearSolver(ida_mem, LS);
  if(check_flag(&flag, "IDASpilsSetLinearSolver", 1)) return(1);

  flag = IDASpilsSetPreconditioner(ida_mem, Precond, PSolve);
  if(check_flag(&flag, "IDASpilsSetPreconditioner", 1)) return(1);

  /* Call IDACalcIC (with default options) to correct the initial values. */

  tout = RCONST(0.001);
  flag = IDACalcIC(ida_mem, IDA_YA_YDP_INIT, tout);
  if(check_flag(&flag, "IDACalcIC", 1)) return(1);
  
  /* Print heading, basic parameters, and initial values. */

  PrintHeader(maxl, rtol, atol);
  PrintOutput(ida_mem, cc, ZERO);
  
  /* Loop over iout, call IDASolve (normal mode), print selected output. */
  
  for (iout = 1; iout <= NOUT; iout++) {
    
    flag = IDASolve(ida_mem, tout, &tret, cc, cp, IDA_NORMAL);
    if(check_flag(&flag, "IDASolve", 1)) return(flag);
    
    PrintOutput(ida_mem, cc, tret);
    
    if (iout < 3) tout *= TMULT; else tout += TADD;
    
  }
  
  /* Print final statistics and free memory. */  
  
  PrintFinalStats(ida_mem);
  printf("num_threads = %i\n\n", num_threads);

  /* Free memory */

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);

  N_VDestroy_OpenMP(cc);
  N_VDestroy_OpenMP(cp);
  N_VDestroy_OpenMP(id);


  destroyMat(webdata->acoef);
  N_VDestroy_OpenMP(webdata->rates);
  N_VDestroy_OpenMP(webdata->ewt);
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
  
  cpv = NV_DATA_OMP(cp);
  resv = NV_DATA_OMP(res);
  np = webdata->np;
  
  /* Call Fweb to set res to vector of right-hand sides. */
  Fweb(tt, cc, res, webdata);
  
  /* Loop over all grid points, setting residual values appropriately
     for differential or algebraic components.                        */
#pragma omp parallel for default(shared) private(jy, jx, is, yloc, loc) schedule(static) num_threads(webdata->nthreads)  
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


static int Precond(realtype tt, N_Vector cc, N_Vector cp,
                   N_Vector rr, realtype cj, void *user_data)
{
  int flag;
  realtype uround, xx, yy, del_x, del_y;
  realtype **Pxy, *ratesxy, *Pxycol, *cxy, *cpxy, *ewtxy, cctmp;
  realtype inc, fac, sqru, perturb_rates[NUM_SPECIES];
  int is, js, jx, jy, ret;
  void *ida_mem;
  N_Vector ewt;
  realtype hh;
  UserData webdata;

  webdata = (UserData) user_data;
  del_x = webdata->dx;
  del_y = webdata->dy;

  uround = UNIT_ROUNDOFF;
  sqru = SUNRsqrt(uround);

  ida_mem = webdata->ida_mem;
  ewt = webdata->ewt;
  flag = IDAGetErrWeights(ida_mem, ewt);
  if(check_flag(&flag, "IDAGetErrWeights", 1)) return(1);
  flag = IDAGetCurrentStep(ida_mem, &hh);
  if(check_flag(&flag, "IDAGetCurrentStep", 1)) return(1);

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


static int PSolve(realtype tt, N_Vector cc, N_Vector cp,
                  N_Vector rr, N_Vector rvec, N_Vector zvec,
		  realtype cj, realtype dalta, void *user_data) 
{
  realtype **Pxy, *zxy;
  sunindextype *pivot;
  int jx, jy;
  UserData webdata;
  
  webdata = (UserData) user_data;

  N_VScale(ONE, rvec, zvec);

#pragma omp parallel for collapse(2) default(shared) private(jx, jy, zxy, Pxy, pivot) schedule(static) num_threads(webdata->nthreads)
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
  
  ccv = NV_DATA_OMP(cc);
  cpv = NV_DATA_OMP(cp);
  idv = NV_DATA_OMP(id);
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
  printf("\nidasFoodWeb_kry_omp: Predator-prey DAE OpenMP example problem using Krylov solver for IDAS \n\n");
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
  printf("Linear solver: SUNSPGMR, maxl = %ld\n",(long int) maxl);
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

static void PrintOutput(void *ida_mem, N_Vector c, realtype t)
{
  int i, kused, flag;
  long int nst;
  realtype *c_bl, *c_tr, hused;

  flag = IDAGetLastOrder(ida_mem, &kused);
  check_flag(&flag, "IDAGetLastOrder", 1);
  flag = IDAGetNumSteps(ida_mem, &nst);
  check_flag(&flag, "IDAGetNumSteps", 1);
  flag = IDAGetLastStep(ida_mem, &hused);
  check_flag(&flag, "IDAGetLastStep", 1);
  
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

static void PrintFinalStats(void *ida_mem)
{ 
  long int nst, nre, sli, netf, nps, npevals, nrevalsLS;
  int flag;

  flag = IDAGetNumSteps(ida_mem, &nst);
  check_flag(&flag, "IDAGetNumSteps", 1);
  flag = IDASpilsGetNumLinIters(ida_mem, &sli);
  check_flag(&flag, "IDAGetNumNonlinSolvIters", 1);
  flag = IDAGetNumResEvals(ida_mem, &nre);
  check_flag(&flag, "IDAGetNumResEvals", 1);
  flag = IDAGetNumErrTestFails(ida_mem, &netf);
  check_flag(&flag, "IDAGetNumErrTestFails", 1);
  flag = IDASpilsGetNumPrecSolves(ida_mem, &nps);
  check_flag(&flag, "IDAGetNumNonlinSolvConvFails", 1);
  flag = IDASpilsGetNumPrecEvals(ida_mem, &npevals);
  check_flag(&flag, "IDADlsGetNumJacEvals", 1);
  flag = IDASpilsGetNumResEvals(ida_mem, &nrevalsLS);
  check_flag(&flag, "IDADlsGetNumResEvals", 1);

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
#pragma omp parallel for default(shared) private(is, dcyli, dcyui, dcxli, dcxui) schedule(static) num_threads(webdata->nthreads)
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
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  if (opt == 0 && flagvalue == NULL) {
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1); 
  }

  return(0);
}
