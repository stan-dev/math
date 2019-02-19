/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example (serial):
 *
 * This example solves a nonlinear system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion
 * on the unit square in two dimensions. The dependent variable
 * vector is the following:
 *
 *       1   2         ns
 * c = (c , c ,  ..., c  )     (denoted by the variable cc)
 *
 * and the PDE's are as follows:
 *
 *                    i       i
 *         0 = d(i)*(c     + c    )  +  f  (x,y,c)   (i=1,...,ns)
 *                    xx      yy         i
 *
 *   where
 *
 *                   i             ns         j
 *   f  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
 *    i                           j=1
 *
 * The number of species is ns = 2 * np, with the first np being
 * prey and the last np being predators. The number np is both the
 * number of prey and predator species. The coefficients a(i,j),
 * b(i), d(i) are:
 *
 *   a(i,i) = -AA   (all i)
 *   a(i,j) = -GG   (i <= np , j >  np)
 *   a(i,j) =  EE   (i >  np,  j <= np)
 *   b(i) = BB * (1 + alpha * x * y)   (i <= np)
 *   b(i) =-BB * (1 + alpha * x * y)   (i >  np)
 *   d(i) = DPREY   (i <= np)
 *   d(i) = DPRED   ( i > np)
 *
 * The various scalar parameters are set using define's or in
 * routine InitUserData.
 *
 * The boundary conditions are: normal derivative = 0, and the
 * initial guess is constant in x and y, but the final solution
 * is not.
 *
 * The PDEs are discretized by central differencing on an MX by
 * MY mesh.
 *
 * The nonlinear system is solved by KINSOL using the method
 * specified in the local variable globalstrat.
 *
 * The preconditioner matrix is a block-diagonal matrix based on
 * the partial derivatives of the interaction terms f only.
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 * References:
 *
 * 1. Peter N. Brown and Youcef Saad,
 *    Hybrid Krylov Methods for Nonlinear Systems of Equations
 *    LLNL report UCRL-97645, November 1987.
 *
 * 2. Peter N. Brown and Alan C. Hindmarsh,
 *    Reduced Storage Matrix Methods in Stiff ODE systems,
 *    Lawrence Livermore National Laboratory Report  UCRL-95088,
 *    Rev. 1, June 1987, and  Journal of Applied Mathematics and
 *    Computation, Vol. 31 (May 1989), pp. 40-91. (Presents a
 *    description of the time-dependent version of this test
 *    problem.)
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts.      */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver      */
#include <sundials/sundials_dense.h>   /* use generic dense solver in precond. */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* access to SUNMAX, SUNRabs, SUNRsqrt  */

/* Problem Constants */

#define NUM_SPECIES     6  /* must equal 2*(number of prey or predators)
                              number of prey = number of predators       */ 

#define PI       RCONST(3.1415926535898)   /* pi */ 

#define MX          8              /* MX = number of x mesh points */
#define MY          8              /* MY = number of y mesh points */
#define NSMX        (NUM_SPECIES * MX)
#define NEQ         (NSMX * MY)    /* number of equations in the system */
#define AA          RCONST(1.0)    /* value of coefficient AA in above eqns */
#define EE          RCONST(10000.) /* value of coefficient EE in above eqns */
#define GG          RCONST(0.5e-6) /* value of coefficient GG in above eqns */
#define BB          RCONST(1.0)    /* value of coefficient BB in above eqns */
#define DPREY       RCONST(1.0)    /* value of coefficient dprey above */
#define DPRED       RCONST(0.5)    /* value of coefficient dpred above */
#define ALPHA       RCONST(1.0)    /* value of coefficient alpha above */
#define AX          RCONST(1.0)    /* total range of x variable */
#define AY          RCONST(1.0)    /* total range of y variable */
#define FTOL        RCONST(1.e-7)  /* ftol tolerance */
#define STOL        RCONST(1.e-13) /* stol tolerance */
#define THOUSAND    RCONST(1000.0) /* one thousand */
#define ZERO        RCONST(0.0)    /* 0. */
#define ONE         RCONST(1.0)    /* 1. */
#define TWO         RCONST(2.0)    /* 2. */
#define PREYIN      RCONST(1.0)    /* initial guess for prey concentrations. */
#define PREDIN      RCONST(30000.0)/* initial guess for predator concs.      */

/* User-defined vector access macro: IJ_Vptr */

/* IJ_Vptr is defined in order to translate from the underlying 3D structure
   of the dependent variable vector to the 1D storage scheme for an N-vector.
   IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to 
   indices is = 0, jx = i, jy = j.    */

#define IJ_Vptr(vv,i,j)   (&NV_Ith_S(vv, i*NUM_SPECIES + j*NSMX))

/* Type : UserData 
   contains preconditioner blocks, pivot arrays, and problem constants */

typedef struct {
  realtype **P[MX][MY];
  sunindextype *pivot[MX][MY];
  realtype **acoef, *bcoef;
  N_Vector rates;
  realtype *cox, *coy;
  realtype ax, ay, dx, dy;
  realtype uround, sqruround;
  sunindextype mx, my, ns, np;
} *UserData;

/* Functions Called by the KINSOL Solver */

static int func(N_Vector cc, N_Vector fval, void *user_data);

static int PrecSetupBD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       void *user_data);

static int PrecSolveBD(N_Vector cc, N_Vector cscale, 
                       N_Vector fval, N_Vector fscale, 
                       N_Vector vv, void *user_data);

/* Private Helper Functions */

static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector cc, N_Vector sc);
static void PrintHeader(int globalstrategy, int maxl, int maxlrst, 
                        realtype fnormtol, realtype scsteptol);
static void PrintOutput(N_Vector cc);
static void PrintFinalStats(void *kmem);
static void WebRate(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy, 
                    void *user_data);
static realtype DotProd(sunindextype size, realtype *x1, realtype *x2);
static int check_flag(void *flagvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(void)
{
  int globalstrategy;
  realtype fnormtol, scsteptol;
  N_Vector cc, sc, constraints;
  UserData data;
  int flag, maxl, maxlrst;
  void *kmem;
  SUNLinearSolver LS;

  cc = sc = constraints = NULL;
  kmem = NULL;
  LS = NULL;
  data = NULL;

  /* Allocate memory, and set problem data, initial values, tolerances */ 
  globalstrategy = KIN_NONE;

  data = AllocUserData();
  if (check_flag((void *)data, "AllocUserData", 2)) return(1);
  InitUserData(data);

  /* Create serial vectors of length NEQ */
  cc = N_VNew_Serial(NEQ);
  if (check_flag((void *)cc, "N_VNew_Serial", 0)) return(1);
  sc = N_VNew_Serial(NEQ);
  if (check_flag((void *)sc, "N_VNew_Serial", 0)) return(1);
  data->rates = N_VNew_Serial(NEQ);
  if (check_flag((void *)data->rates, "N_VNew_Serial", 0)) return(1);

  constraints = N_VNew_Serial(NEQ);
  if (check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);
  N_VConst(TWO, constraints);

  SetInitialProfiles(cc, sc);

  fnormtol=FTOL; scsteptol=STOL;

  /* Call KINCreate/KINInit to initialize KINSOL.
     A pointer to KINSOL problem memory is returned and stored in kmem. */
  kmem = KINCreate();
  if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

  /* Vector cc passed as template vector. */
  flag = KINInit(kmem, func, cc);
  if (check_flag(&flag, "KINInit", 1)) return(1);

  flag = KINSetUserData(kmem, data);
  if (check_flag(&flag, "KINSetUserData", 1)) return(1);
  flag = KINSetConstraints(kmem, constraints);
  if (check_flag(&flag, "KINSetConstraints", 1)) return(1);
  flag = KINSetFuncNormTol(kmem, fnormtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);
  flag = KINSetScaledStepTol(kmem, scsteptol);
  if (check_flag(&flag, "KINSetScaledStepTol", 1)) return(1);

  /* We no longer need the constraints vector since KINSetConstraints
     creates a private copy for KINSOL to use. */
  N_VDestroy_Serial(constraints);


  /* Create SUNLinSol_SPGMR object with right preconditioning and the 
     maximum Krylov dimension maxl */
  maxl = 15; 
  LS = SUNLinSol_SPGMR(cc, PREC_RIGHT, maxl);
  if(check_flag((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

  /* Attach the linear solver to KINSOL */
  flag = KINSetLinearSolver(kmem, LS, NULL);
  if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

  /* Set the maximum number of restarts */
  maxlrst = 2;
  flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);
  if (check_flag(&flag, "SUNLinSol_SPGMRSetMaxRestarts", 1)) return(1);

  /* Specify the preconditioner setup and solve routines */
  flag = KINSetPreconditioner(kmem, PrecSetupBD, PrecSolveBD);
  if (check_flag(&flag, "KINSetPreconditioner", 1)) return(1);

  /* Print out the problem size, solution parameters, initial guess. */
  PrintHeader(globalstrategy, maxl, maxlrst, fnormtol, scsteptol);

  /* Call KINSol and print output concentration profile */
  flag = KINSol(kmem,           /* KINSol memory block */
                cc,             /* initial guess on input; solution vector */
                globalstrategy, /* global strategy choice */
                sc,             /* scaling vector for the variable cc */
                sc);            /* scaling vector for function values fval */
  if (check_flag(&flag, "KINSol", 1)) return(1);

  printf("\n\nComputed equilibrium species concentrations:\n");
  PrintOutput(cc);

  /* Print final statistics and free memory */  
  PrintFinalStats(kmem);

  N_VDestroy_Serial(cc);
  N_VDestroy_Serial(sc);
  KINFree(&kmem);
  SUNLinSolFree(LS);
  FreeUserData(data);

  return(0);
}

/* Readability definitions used in other routines below */

#define acoef  (data->acoef)
#define bcoef  (data->bcoef)
#define cox    (data->cox)
#define coy    (data->coy)

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

/* 
 * System function for predator-prey system 
 */

static int func(N_Vector cc, N_Vector fval, void *user_data)
{
  realtype xx, yy, delx, dely, *cxy, *rxy, *fxy, dcyli, dcyui, dcxli, dcxri;
  sunindextype jx, jy, is, idyu, idyl, idxr, idxl;
  UserData data;
  
  data = (UserData)user_data;
  delx = data->dx;
  dely = data->dy;
  
  /* Loop over all mesh points, evaluating rate array at each point*/
  for (jy = 0; jy < MY; jy++) {
    
    yy = dely*jy;

    /* Set lower/upper index shifts, special at boundaries. */
    idyl = (jy != 0   ) ? NSMX : -NSMX;
    idyu = (jy != MY-1) ? NSMX : -NSMX;
    
    for (jx = 0; jx < MX; jx++) {

      xx = delx*jx;

      /* Set left/right index shifts, special at boundaries. */
      idxl = (jx !=  0  ) ?  NUM_SPECIES : -NUM_SPECIES;
      idxr = (jx != MX-1) ?  NUM_SPECIES : -NUM_SPECIES;

      cxy = IJ_Vptr(cc,jx,jy);
      rxy = IJ_Vptr(data->rates,jx,jy);
      fxy = IJ_Vptr(fval,jx,jy);

      /* Get species interaction rate array at (xx,yy) */
      WebRate(xx, yy, cxy, rxy, user_data);
      
      for(is = 0; is < NUM_SPECIES; is++) {
        
        /* Differencing in x direction */
        dcyli = *(cxy+is) - *(cxy - idyl + is) ;
        dcyui = *(cxy + idyu + is) - *(cxy+is);
        
        /* Differencing in y direction */
        dcxli = *(cxy+is) - *(cxy - idxl + is);
        dcxri = *(cxy + idxr +is) - *(cxy+is);
        
        /* Compute the total rate value at (xx,yy) */
        fxy[is] = (coy)[is] * (dcyui - dcyli) +
          (cox)[is] * (dcxri - dcxli) + rxy[is];
        
      } /* end of is loop */
      
    } /* end of jx loop */
    
  } /* end of jy loop */

  return(0);
}

/*
 * Preconditioner setup routine. Generate and preprocess P. 
 */

static int PrecSetupBD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       void *user_data)
{
  realtype r, r0, uround, sqruround, xx, yy, delx, dely, csave, fac;
  realtype *cxy, *scxy, **Pxy, *ratesxy, *Pxycol, perturb_rates[NUM_SPECIES];
  sunindextype i, j, jx, jy, ret;
  UserData data;
  
  data = (UserData) user_data;
  delx = data->dx;
  dely = data->dy;
  
  uround = data->uround;
  sqruround = data->sqruround;
  fac = N_VWL2Norm(fval, fscale);
  r0 = THOUSAND * uround * fac * NEQ;
  if(r0 == ZERO) r0 = ONE;
  
  /* Loop over spatial points; get size NUM_SPECIES Jacobian block at each */
  for (jy = 0; jy < MY; jy++) {
    yy = jy*dely;
    
    for (jx = 0; jx < MX; jx++) {
      xx = jx*delx;
      Pxy = (data->P)[jx][jy];
      cxy = IJ_Vptr(cc,jx,jy);
      scxy= IJ_Vptr(cscale,jx,jy);
      ratesxy = IJ_Vptr((data->rates),jx,jy);

      /* Compute difference quotients of interaction rate fn. */
      for (j = 0; j < NUM_SPECIES; j++) {
        
        csave = cxy[j];  /* Save the j,jx,jy element of cc */
        r = SUNMAX(sqruround*SUNRabs(csave), r0/scxy[j]);
        cxy[j] += r; /* Perturb the j,jx,jy element of cc */
        fac = ONE/r;
        
        WebRate(xx, yy, cxy, perturb_rates, data);
        
        /* Restore j,jx,jy element of cc */
        cxy[j] = csave;
        
        /* Load the j-th column of difference quotients */
        Pxycol = Pxy[j];
        for (i = 0; i < NUM_SPECIES; i++)
          Pxycol[i] = (perturb_rates[i] - ratesxy[i]) * fac;
        
      } /* end of j loop */
      
      /* Do LU decomposition of size NUM_SPECIES preconditioner block */
      ret = denseGETRF(Pxy, NUM_SPECIES, NUM_SPECIES, (data->pivot)[jx][jy]);
      if (ret != 0) return(1);
      
    } /* end of jx loop */
    
  } /* end of jy loop */
  
  return(0);
}

/*
 * Preconditioner solve routine 
 */

static int PrecSolveBD(N_Vector cc, N_Vector cscale, 
                       N_Vector fval, N_Vector fscale, 
                       N_Vector vv, void *user_data)
{
  realtype **Pxy, *vxy;
  sunindextype *piv, jx, jy;
  UserData data;
  
  data = (UserData)user_data;
  
  for (jx=0; jx<MX; jx++) {
    
    for (jy=0; jy<MY; jy++) {
      
      /* For each (jx,jy), solve a linear system of size NUM_SPECIES.
         vxy is the address of the corresponding portion of the vector vv;
         Pxy is the address of the corresponding block of the matrix P;
         piv is the address of the corresponding block of the array pivot. */
      vxy = IJ_Vptr(vv,jx,jy);
      Pxy = (data->P)[jx][jy];
      piv = (data->pivot)[jx][jy];
      denseGETRS(Pxy, NUM_SPECIES, piv, vxy);
      
    } /* end of jy loop */
    
  } /* end of jx loop */
  
  return(0);
}

/*
 * Interaction rate function routine 
 */

static void WebRate(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy, 
                    void *user_data)
{
  sunindextype i;
  realtype fac;
  UserData data;
  
  data = (UserData)user_data;
  
  for (i = 0; i<NUM_SPECIES; i++)
    ratesxy[i] = DotProd(NUM_SPECIES, cxy, acoef[i]);
  
  fac = ONE + ALPHA * xx * yy;
  
  for (i = 0; i < NUM_SPECIES; i++)
    ratesxy[i] = cxy[i] * ( bcoef[i] * fac + ratesxy[i] );
}

/*
 * Dot product routine for realtype arrays 
 */

static realtype DotProd(sunindextype size, realtype *x1, realtype *x2)
{
  sunindextype i;
  realtype *xx1, *xx2, temp = ZERO;
  
  xx1 = x1; xx2 = x2;
  for (i = 0; i < size; i++) temp += (*xx1++) * (*xx2++);

  return(temp);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Allocate memory for data structure of type UserData 
 */

static UserData AllocUserData(void)
{
  int jx, jy;
  UserData data;

  data = (UserData) malloc(sizeof *data);
  
  for (jx=0; jx < MX; jx++) {
    for (jy=0; jy < MY; jy++) {
      (data->P)[jx][jy] = newDenseMat(NUM_SPECIES, NUM_SPECIES);
      (data->pivot)[jx][jy] = newIndexArray(NUM_SPECIES);
    }
  }

  acoef = newDenseMat(NUM_SPECIES, NUM_SPECIES);
  bcoef = (realtype *)malloc(NUM_SPECIES * sizeof(realtype));
  cox   = (realtype *)malloc(NUM_SPECIES * sizeof(realtype));
  coy   = (realtype *)malloc(NUM_SPECIES * sizeof(realtype));
  
  return(data);
}

/*
 * Load problem constants in data 
 */

static void InitUserData(UserData data)
{
  sunindextype i, j, np;
  realtype *a1,*a2, *a3, *a4, dx2, dy2;

  data->mx = MX;
  data->my = MY;
  data->ns = NUM_SPECIES;
  data->np = NUM_SPECIES/2;
  data->ax = AX;
  data->ay = AY;
  data->dx = (data->ax)/(MX-1);
  data->dy = (data->ay)/(MY-1);
  data->uround = UNIT_ROUNDOFF;
  data->sqruround = SUNRsqrt(data->uround);

  /* Set up the coefficients a and b plus others found in the equations */
  np = data->np;

  dx2=(data->dx)*(data->dx); dy2=(data->dy)*(data->dy);

  for (i = 0; i < np; i++) {
    a1= &(acoef[i][np]);
    a2= &(acoef[i+np][0]);
    a3= &(acoef[i][0]);
    a4= &(acoef[i+np][np]);

    /*  Fill in the portion of acoef in the four quadrants, row by row */
    for (j = 0; j < np; j++) {
      *a1++ =  -GG;
      *a2++ =   EE;
      *a3++ = ZERO;
      *a4++ = ZERO;
    }

    /* and then change the diagonal elements of acoef to -AA */
    acoef[i][i]=-AA;
    acoef[i+np][i+np] = -AA;

    bcoef[i] = BB;
    bcoef[i+np] = -BB;

    cox[i]=DPREY/dx2;
    cox[i+np]=DPRED/dx2;

    coy[i]=DPREY/dy2;
    coy[i+np]=DPRED/dy2;
  }
}

/*
 * Free data memory 
 */

static void FreeUserData(UserData data)
{
  int jx, jy;

  for (jx=0; jx < MX; jx++) {
    for (jy=0; jy < MY; jy++) {
      destroyMat((data->P)[jx][jy]);
      destroyArray((data->pivot)[jx][jy]);
    }
  }

  destroyMat(acoef);
  free(bcoef);
  free(cox);
  free(coy);
  N_VDestroy_Serial(data->rates);
  free(data);
}

/*
 * Set initial conditions in cc 
 */

static void SetInitialProfiles(N_Vector cc, N_Vector sc)
{
  int i, jx, jy;
  realtype *cloc, *sloc;
  realtype  ctemp[NUM_SPECIES], stemp[NUM_SPECIES];
  
  /* Initialize arrays ctemp and stemp used in the loading process */
  for (i = 0; i < NUM_SPECIES/2; i++) {
    ctemp[i] = PREYIN;
    stemp[i] = ONE;
  }
  for (i = NUM_SPECIES/2; i < NUM_SPECIES; i++) {
    ctemp[i] = PREDIN;
    stemp[i] = RCONST(0.00001);
  }

  /* Load initial profiles into cc and sc vector from ctemp and stemp. */
  for (jy = 0; jy < MY; jy++) {
    for (jx = 0; jx < MX; jx++) {
      cloc = IJ_Vptr(cc,jx,jy);
      sloc = IJ_Vptr(sc,jx,jy);
      for (i = 0; i < NUM_SPECIES; i++) {
        cloc[i] = ctemp[i];
        sloc[i] = stemp[i];
      }
    }
  }
}

/* 
 * Print first lines of output (problem description)
 */

static void PrintHeader(int globalstrategy, int maxl, int maxlrst, 
                        realtype fnormtol, realtype scsteptol)
{
  printf("\nPredator-prey test problem --  KINSol (serial version)\n\n");
  printf("Mesh dimensions = %d X %d\n", MX, MY);
  printf("Number of species = %d\n", NUM_SPECIES);
  printf("Total system size = %d\n\n", NEQ);
  printf("Flag globalstrategy = %d (0 = None, 1 = Linesearch)\n",
         globalstrategy);
  printf("Linear solver is SPGMR with maxl = %d, maxlrst = %d\n",
         maxl, maxlrst);
  printf("Preconditioning uses interaction-only block-diagonal matrix\n");
  printf("Positivity constraints imposed on all components \n");
#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf("Tolerance parameters:  fnormtol = %Lg   scsteptol = %Lg\n",
         fnormtol, scsteptol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  fnormtol = %g   scsteptol = %g\n",
         fnormtol, scsteptol);
#else
  printf("Tolerance parameters:  fnormtol = %g   scsteptol = %g\n",
         fnormtol, scsteptol);
#endif

  printf("\nInitial profile of concentration\n");
#if defined(SUNDIALS_EXTENDED_PRECISION) 
  printf("At all mesh points:  %Lg %Lg %Lg   %Lg %Lg %Lg\n", 
         PREYIN, PREYIN, PREYIN,
         PREDIN, PREDIN, PREDIN);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At all mesh points:  %g %g %g   %g %g %g\n", 
         PREYIN, PREYIN, PREYIN,
         PREDIN, PREDIN, PREDIN);
#else
  printf("At all mesh points:  %g %g %g   %g %g %g\n", 
         PREYIN, PREYIN, PREYIN,
         PREDIN, PREDIN, PREDIN);
#endif
}

/*
 * Print sampled values of current cc 
 */

static void PrintOutput(N_Vector cc)
{
  int is, jx, jy;
  realtype *ct;
  
  jy = 0; jx = 0;
  ct = IJ_Vptr(cc,jx,jy);
  printf("\nAt bottom left:");

  /* Print out lines with up to 6 values per line */
  for (is = 0; is < NUM_SPECIES; is++){
    if ((is%6)*6 == is) printf("\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %Lg",ct[is]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %g",ct[is]);
#else
    printf(" %g",ct[is]);
#endif
  }
  
  jy = MY-1; jx = MX-1;
  ct = IJ_Vptr(cc,jx,jy);
  printf("\n\nAt top right:");

  /* Print out lines with up to 6 values per line */
  for (is = 0; is < NUM_SPECIES; is++) {
    if ((is%6)*6 == is) printf("\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %Lg",ct[is]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %g",ct[is]);
#else
    printf(" %g",ct[is]);
#endif
  }
  printf("\n\n");
}

/*
 * Print final statistics contained in iopt 
 */

static void PrintFinalStats(void *kmem)
{
  long int nni, nfe, nli, npe, nps, ncfl, nfeSG;
  int flag;

  flag = KINGetNumNonlinSolvIters(kmem, &nni);
  check_flag(&flag, "KINGetNumNonlinSolvIters", 1);
  flag = KINGetNumFuncEvals(kmem, &nfe);
  check_flag(&flag, "KINGetNumFuncEvals", 1);
  flag = KINGetNumLinIters(kmem, &nli);
  check_flag(&flag, "KINGetNumLinIters", 1);
  flag = KINGetNumPrecEvals(kmem, &npe);
  check_flag(&flag, "KINGetNumPrecEvals", 1);
  flag = KINGetNumPrecSolves(kmem, &nps);
  check_flag(&flag, "KINGetNumPrecSolves", 1);
  flag = KINGetNumLinConvFails(kmem, &ncfl);
  check_flag(&flag, "KINGetNumLinConvFails", 1);
  flag = KINGetNumLinFuncEvals(kmem, &nfeSG);
  check_flag(&flag, "KINGetNumLinFuncEvals", 1);

  printf("Final Statistics.. \n");
  printf("nni    = %5ld    nli   = %5ld\n", nni, nli);
  printf("nfe    = %5ld    nfeSG = %5ld\n", nfe, nfeSG);
  printf("nps    = %5ld    npe   = %5ld     ncfl  = %5ld\n", nps, npe, ncfl);

}

/*
 * Check function return value...
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
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  return(0);
}
