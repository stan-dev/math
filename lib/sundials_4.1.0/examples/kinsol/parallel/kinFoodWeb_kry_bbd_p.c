/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
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
 * Example problem for KINSOL (parallel machine case) using the BBD
 * preconditioner.
 *
 * This example solves a nonlinear system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion on
 * the unit square in two dimensions. The dependent variable vector
 * is the following:
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
 * The preconditioner matrix is a band-block-diagonal matrix
 * using the KINBBDPRE module. The half-bandwidths are as follows:
 *
 *   Difference quotient half-bandwidths mldq = mudq = 2*ns - 1
 *   Retained banded blocks have half-bandwidths mlkeep = mukeep = ns.
 *
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
 * --------------------------------------------------------------------------
 *  Run command line: mpirun -np N -machinefile machines kinFoodWeb_kry_bbd_p
 *  where N = NPEX * NPEY is the number of processors.
 * --------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts.      */
#include <nvector/nvector_parallel.h>  /* access to MPI parallel N_Vector      */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver      */
#include <kinsol/kinsol_bbdpre.h>      /* access to BBD preconditioner         */
#include <sundials/sundials_dense.h>   /* use generic dense solver in precond. */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* access to SUNMAX, SUNRabs, SUNRsqrt  */

#include <mpi.h>

/* Problem Constants */

#define NUM_SPECIES     6  /* must equal 2*(number of prey or predators)
                              number of prey = number of predators       */ 

#define PI       RCONST(3.1415926535898)   /* pi */ 

#define NPEX        2            /* number of processors in the x-direction  */
#define NPEY        2            /* number of processors in the y-direction  */
#define MXSUB       10           /* number of x mesh points per subgrid      */
#define MYSUB       10           /* number of y mesh points per subgrid      */
#define MX          (NPEX*MXSUB) /* number of mesh points in the x-direction */
#define MY          (NPEY*MYSUB) /* number of mesh points in the y-direction */
#define NSMXSUB     (NUM_SPECIES * MXSUB)
#define NSMXSUB2    (NUM_SPECIES * (MXSUB+2))
#define NEQ         (NUM_SPECIES*MX*MY)  /* number of equations in the system */
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
#define PREYIN      RCONST(1.0)    /* initial guess for prey concentrations. */
#define PREDIN      RCONST(30000.0)/* initial guess for predator concs.      */

/* User-defined vector access macro: IJ_Vptr */

/* IJ_Vptr is defined in order to translate from the underlying 3D structure
   of the dependent variable vector to the 1D storage scheme for an N-vector.
   IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to 
   indices is = 0, jx = i, jy = j.    */

#define IJ_Vptr(vv,i,j)   (&NV_Ith_P(vv, i*NUM_SPECIES + j*NSMXSUB))

/* Type : UserData 
   contains problem constants and extended array */

typedef struct {
  realtype **acoef, *bcoef;
  N_Vector rates;
  realtype *cox, *coy;
  realtype ax, ay, dx, dy;
  sunindextype Nlocal;
  int mx, my, ns, np;
  realtype cext[NUM_SPECIES * (MXSUB+2)*(MYSUB+2)];
  int my_pe, isubx, isuby, nsmxsub, nsmxsub2;
  MPI_Comm comm;
} *UserData;

/* Functions called by the KINSOL Solver */

static int func(N_Vector cc, N_Vector fval, void *user_data);

static int ccomm(sunindextype Nlocal, N_Vector cc, void *data);

static int func_local(sunindextype Nlocal, N_Vector cc, N_Vector fval, void *user_data);

/* Private Helper Functions */

static UserData AllocUserData(void);
static void InitUserData(int my_pe, sunindextype Nlocal, MPI_Comm comm, UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector cc, N_Vector sc);
static void PrintHeader(int globalstrategy, int maxl, int maxlrst,
                        sunindextype mudq, sunindextype mldq,
                        sunindextype mukeep, sunindextype mlkeep,
                        realtype fnormtol, realtype scsteptol);
static void PrintOutput(int my_pe, MPI_Comm comm, N_Vector cc);
static void PrintFinalStats(void *kmem);
static void WebRate(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy,
                    void *user_data);
static realtype DotProd(int size, realtype *x1, realtype *x2);
static void BSend(MPI_Comm comm, int my_pe, int isubx,
                  int isuby, int dsizex, int dsizey,
                  realtype *cdata);
static void BRecvPost(MPI_Comm comm, MPI_Request request[], int my_pe,
                      int isubx, int isuby,
                      int dsizex, int dsizey,
                      realtype *cext, realtype *buffer);
static void BRecvWait(MPI_Request request[], int isubx,
                      int isuby, int dsizex, realtype *cext,
                      realtype *buffer);
static int check_flag(void *flagvalue, const char *funcname, int opt, int id);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  int globalstrategy;
  sunindextype Nlocal;
  realtype fnormtol, scsteptol, dq_rel_uu;
  N_Vector cc, sc, constraints;
  UserData data;
  int flag, maxl, maxlrst;
  int my_pe, npes, npelast = NPEX*NPEY-1;
  void *kmem;
  SUNLinearSolver LS;
  sunindextype mudq, mldq, mukeep, mlkeep;
  MPI_Comm comm;

  cc = sc = constraints = NULL;
  kmem = NULL;
  LS = NULL;
  data = NULL;

  /* Get processor number and total number of pe's */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &my_pe);

  if (npes != NPEX*NPEY) {
    if (my_pe == 0)
      fprintf(stderr, "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n",
	          npes,NPEX*NPEY);
    MPI_Finalize();
    return(1);
  }

  /* Allocate memory, and set problem data, initial values, tolerances */ 

  /* Set local vector length */
  Nlocal = NUM_SPECIES*MXSUB*MYSUB;

  /* Allocate and initialize user data block */
  data = AllocUserData();
  if (check_flag((void *)data, "AllocUserData", 2, my_pe)) MPI_Abort(comm, 1);
  InitUserData(my_pe, Nlocal, comm, data);

  /* Set global strategy flag */
  globalstrategy = KIN_NONE;

  /* Allocate and initialize vectors */
  cc = N_VNew_Parallel(comm, Nlocal, NEQ);
  if (check_flag((void *)cc, "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  sc = N_VNew_Parallel(comm, Nlocal, NEQ);
  if (check_flag((void *)sc, "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  data->rates = N_VNew_Parallel(comm, Nlocal, NEQ);
  if (check_flag((void *)data->rates, "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  constraints = N_VNew_Parallel(comm, Nlocal, NEQ);
  if (check_flag((void *)constraints, "N_VNew_Parallel", 0, my_pe)) MPI_Abort(comm, 1);
  N_VConst(ZERO, constraints);

  SetInitialProfiles(cc, sc);

  fnormtol = FTOL; scsteptol = STOL;

  /* Call KINCreate/KINInit to initialize KINSOL:
     nvSpec  points to machine environment data
     A pointer to KINSOL problem memory is returned and stored in kmem. */
  kmem = KINCreate();
  if (check_flag((void *)kmem, "KINCreate", 0, my_pe)) MPI_Abort(comm, 1);

  /* Vector cc passed as template vector. */
  flag = KINInit(kmem, func, cc);
  if (check_flag(&flag, "KINInit", 1, my_pe)) MPI_Abort(comm, 1);
  flag = KINSetUserData(kmem, data);
  if (check_flag(&flag, "KINSetUserData", 1, my_pe)) MPI_Abort(comm, 1);
  flag = KINSetConstraints(kmem, constraints);
  if (check_flag(&flag, "KINSetConstraints", 1, my_pe)) MPI_Abort(comm, 1);
  flag = KINSetFuncNormTol(kmem, fnormtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1, my_pe)) MPI_Abort(comm, 1);
  flag = KINSetScaledStepTol(kmem, scsteptol);
  if (check_flag(&flag, "KINSetScaledStepTol", 1, my_pe)) MPI_Abort(comm, 1);

  /* We no longer need the constraints vector since KINSetConstraints
     creates a private copy for KINSOL to use. */
  N_VDestroy_Parallel(constraints);

  /* Create SUNLinSol_SPGMR object with right preconditioning and the 
     maximum Krylov dimension maxl */
  maxl = 20; 
  LS = SUNLinSol_SPGMR(cc, PREC_RIGHT, maxl);
  if(check_flag((void *)LS, "SUNLinSol_SPGMR", 0, my_pe)) MPI_Abort(comm, 1);

  /* Attach the linear solver to KINSOL */
  flag = KINSetLinearSolver(kmem, LS, NULL);
  if (check_flag(&flag, "KINSetLinearSolver", 1, my_pe)) MPI_Abort(comm, 1);

  /* Set the maximum number of restarts */
  maxlrst = 2;
  flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);
  if (check_flag(&flag, "SUNLinSol_SPGMRSetMaxRestarts", 1, my_pe)) MPI_Abort(comm, 1);
  
  /* Call KINBBDPrecInit to initialize and allocate memory for the
     band-block-diagonal preconditioner, and specify the local and
     communication functions func_local and gcomm=NULL (all communication
     needed for the func_local is already done in func). */
  dq_rel_uu = ZERO;
  mudq = mldq = 2*NUM_SPECIES - 1;
  mukeep = mlkeep = NUM_SPECIES;

  /* Initialize BBD preconditioner */
  flag = KINBBDPrecInit(kmem, Nlocal, mudq, mldq, mukeep, mlkeep,
                        dq_rel_uu, func_local, NULL);
  if (check_flag(&flag, "KINBBDPrecInit", 1, my_pe)) MPI_Abort(comm, 1);

  /* Print out the problem size, solution parameters, initial guess. */
  if (my_pe == 0)
    PrintHeader(globalstrategy, maxl, maxlrst, mudq, mldq, mukeep,
		mlkeep, fnormtol, scsteptol);

  /* Call KINSol and print output concentration profile */
  flag = KINSol(kmem,           /* KINSol memory block */
                cc,             /* initial guess on input; solution vector */
                globalstrategy, /* global strategy choice */
                sc,             /* scaling vector for the variable cc */
                sc);            /* scaling vector for function values fval */
  if (check_flag(&flag, "KINSol", 1, my_pe)) MPI_Abort(comm, 1);

  if (my_pe == 0) printf("\n\nComputed equilibrium species concentrations:\n");
  if (my_pe == 0 || my_pe == npelast) PrintOutput(my_pe, comm, cc);

  /* Print final statistics and free memory */
  if (my_pe == 0) PrintFinalStats(kmem);

  N_VDestroy_Parallel(cc);
  N_VDestroy_Parallel(sc);
  KINFree(&kmem);
  SUNLinSolFree(LS);
  FreeUserData(data);

  MPI_Finalize();

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
 * ccomm routine.  This routine performs all communication 
 * between processors of data needed to calculate f. 
 */

static int ccomm(sunindextype Nlocal, N_Vector cc, void *userdata)
{

  realtype *cdata, *cext, buffer[2*NUM_SPECIES*MYSUB];
  UserData data;
  MPI_Comm comm;
  int my_pe, isubx, isuby, nsmxsub, nsmysub;
  MPI_Request request[4];

  /* Get comm, my_pe, subgrid indices, data sizes, extended array cext */
  data = (UserData) userdata;
  comm = data->comm;  my_pe = data->my_pe;
  isubx = data->isubx;   isuby = data->isuby;
  nsmxsub = data->nsmxsub;
  nsmysub = NUM_SPECIES*MYSUB;
  cext = data->cext;

  cdata = N_VGetArrayPointer_Parallel(cc);

  /* Start receiving boundary data from neighboring PEs */
  BRecvPost(comm, request, my_pe, isubx, isuby, nsmxsub, nsmysub, cext, buffer);

  /* Send data from boundary of local grid to neighboring PEs */
  BSend(comm, my_pe, isubx, isuby, nsmxsub, nsmysub, cdata);

  /* Finish receiving boundary data from neighboring PEs */
  BRecvWait(request, isubx, isuby, nsmxsub, cext, buffer);

  return(0);
}

/*
 * System function for predator-prey system - calculation part 
 */

static int func_local(sunindextype Nlocal, N_Vector cc, N_Vector fval, void *user_data)
{
  realtype xx, yy, *cxy, *rxy, *fxy, dcydi, dcyui, dcxli, dcxri;
  realtype *cext, dely, delx, *cdata;
  int i, jx, jy, is, ly;
  int isubx, isuby, nsmxsub, nsmxsub2;
  int shifty, offsetc, offsetce, offsetcl, offsetcr, offsetcd, offsetcu;
  UserData data;
  
  data = (UserData)user_data;
  cdata = N_VGetArrayPointer_Parallel(cc);
  
  /* Get subgrid indices, data sizes, extended work array cext */
  isubx = data->isubx;   isuby = data->isuby;
  nsmxsub = data->nsmxsub; nsmxsub2 = data->nsmxsub2;
  cext = data->cext;
  
  /* Copy local segment of cc vector into the working extended array cext */
  offsetc = 0;
  offsetce = nsmxsub2 + NUM_SPECIES;
  for (ly = 0; ly < MYSUB; ly++) {
    for (i = 0; i < nsmxsub; i++) cext[offsetce+i] = cdata[offsetc+i];
    offsetc = offsetc + nsmxsub;
    offsetce = offsetce + nsmxsub2;
  }
  
  /* To facilitate homogeneous Neumann boundary conditions, when this is a
     boundary PE, copy data from the first interior mesh line of cc to cext */
  
  /* If isuby = 0, copy x-line 2 of cc to cext */
  if (isuby == 0) {
    for (i = 0; i < nsmxsub; i++) cext[NUM_SPECIES+i] = cdata[nsmxsub+i];
  }
  
  /* If isuby = NPEY-1, copy x-line MYSUB-1 of cc to cext */
  if (isuby == NPEY-1) {
    offsetc = (MYSUB-2)*nsmxsub;
    offsetce = (MYSUB+1)*nsmxsub2 + NUM_SPECIES;
    for (i = 0; i < nsmxsub; i++) cext[offsetce+i] = cdata[offsetc+i];
  }
  
  /* If isubx = 0, copy y-line 2 of cc to cext */
  if (isubx == 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetc = ly*nsmxsub + NUM_SPECIES;
      offsetce = (ly+1)*nsmxsub2;
      for (i = 0; i < NUM_SPECIES; i++) cext[offsetce+i] = cdata[offsetc+i];
    }
  }
  
  /* If isubx = NPEX-1, copy y-line MXSUB-1 of cc to cext */
  if (isubx == NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetc = (ly+1)*nsmxsub - 2*NUM_SPECIES;
      offsetce = (ly+2)*nsmxsub2 - NUM_SPECIES;
      for (i = 0; i < NUM_SPECIES; i++) cext[offsetce+i] = cdata[offsetc+i];
    }
  }

  /* Loop over all mesh points, evaluating rate arra at each point */
  delx = data->dx;
  dely = data->dy;
  shifty = (MXSUB+2)*NUM_SPECIES;

  for (jy = 0; jy < MYSUB; jy++) {

    yy = dely*(jy + isuby * MYSUB);

    for (jx = 0; jx < MXSUB; jx++) {

      xx = delx * (jx + isubx * MXSUB);
      cxy = IJ_Vptr(cc,jx,jy);
      rxy = IJ_Vptr(data->rates,jx,jy);
      fxy = IJ_Vptr(fval,jx,jy);
      
      WebRate(xx, yy, cxy, rxy, user_data);

      offsetc = (jx+1)*NUM_SPECIES + (jy+1)*NSMXSUB2;
      offsetcd = offsetc - shifty;
      offsetcu = offsetc + shifty;
      offsetcl = offsetc - NUM_SPECIES;
      offsetcr = offsetc + NUM_SPECIES;
      
      for (is = 0; is < NUM_SPECIES; is++) {
        
        /* differencing in x */
        dcydi = cext[offsetc+is]  - cext[offsetcd+is];
        dcyui = cext[offsetcu+is] - cext[offsetc+is];
        
        /* differencing in y */
        dcxli = cext[offsetc+is]  - cext[offsetcl+is];
        dcxri = cext[offsetcr+is] - cext[offsetc+is];
        
        /* compute the value at xx , yy */
        fxy[is] = (coy)[is] * (dcyui - dcydi) +
          (cox)[is] * (dcxri - dcxli) + rxy[is];
        
      } /* end of is loop */
      
    } /* end of jx loop */
    
  } /* end of jy loop */

  return(0);
}

/*
 * System function routine.  Evaluate f(cc).  First call ccomm to do 
 * communication of subgrid boundary data into cext.  Then calculate f
 * by a call to func_local. 
 */

static int func(N_Vector cc, N_Vector fval, void *user_data)
{
  UserData data;

  data = (UserData) user_data;

  /* Call ccomm to do inter-processor communicaiton */
  ccomm(data->Nlocal, cc, data);

  /* Call func_local to calculate all right-hand sides */
  func_local(data->Nlocal, cc, fval, data);

  return(0);
}

/*
 * Interaction rate function routine 
 */

static void WebRate(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy, 
                    void *user_data)
{
  int i;
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

static realtype DotProd(int size, realtype *x1, realtype *x2)
{
  int i;
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
  UserData data;
  
  data = (UserData) malloc(sizeof *data);
  
  acoef = newDenseMat(NUM_SPECIES, NUM_SPECIES);
  bcoef = (realtype *)malloc(NUM_SPECIES * sizeof(realtype));
  cox   = (realtype *)malloc(NUM_SPECIES * sizeof(realtype));
  coy   = (realtype *)malloc(NUM_SPECIES * sizeof(realtype));
  
  return(data);
}

/* 
 * Load problem constants in data 
 */

static void InitUserData(int my_pe, sunindextype Nlocal, MPI_Comm comm, UserData data)
{
  int i, j, np;
  realtype *a1,*a2, *a3, *a4, dx2, dy2;

  data->mx = MX;
  data->my = MY;
  data->ns = NUM_SPECIES;
  data->np = NUM_SPECIES/2;
  data->ax = AX;
  data->ay = AY;
  data->dx = (data->ax)/(MX-1);
  data->dy = (data->ay)/(MY-1);
  data->my_pe = my_pe;
  data->Nlocal = Nlocal;
  data->comm = comm;
  data->isuby = my_pe/NPEX;
  data->isubx = my_pe - data->isuby*NPEX;
  data->nsmxsub = NUM_SPECIES * MXSUB;
  data->nsmxsub2 = NUM_SPECIES * (MXSUB+2);
  
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

  destroyMat(acoef);
  free(bcoef);
  free(cox); free(coy);
  N_VDestroy_Parallel(data->rates);

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
  for (jy = 0; jy < MYSUB; jy++) {
    for (jx=0; jx < MXSUB; jx++) {
      cloc = IJ_Vptr(cc,jx,jy);
      sloc = IJ_Vptr(sc,jx,jy);
      for (i = 0; i < NUM_SPECIES; i++){
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
                        sunindextype mudq, sunindextype mldq,
			sunindextype mukeep, sunindextype mlkeep,
                        realtype fnormtol, realtype scsteptol)
{
    printf("\nPredator-prey test problem--  KINSol (parallel-BBD version)\n\n");
    printf("Mesh dimensions = %d X %d\n", MX, MY);
    printf("Number of species = %d\n", NUM_SPECIES);
    printf("Total system size = %d\n\n", NEQ);
    printf("Subgrid dimensions = %d X %d\n", MXSUB, MYSUB);
    printf("Processor array is %d X %d\n\n", NPEX, NPEY);
    printf("Flag globalstrategy = %d (0 = None, 1 = Linesearch)\n",
           globalstrategy);
    printf("Linear solver is SPGMR with maxl = %d, maxlrst = %d\n",
           maxl, maxlrst);
    printf("Preconditioning uses band-block-diagonal matrix from KINBBDPRE\n");
    printf("  Difference quotient half-bandwidths: mudq = %ld, mldq = %ld\n",
	   (long int) mudq, (long int) mldq);
    printf("  Retained band block half-bandwidths: mukeep = %ld, mlkeep = %ld\n",
	   (long int) mukeep, (long int) mlkeep);
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
           PREYIN,PREYIN,PREYIN, PREDIN,PREDIN,PREDIN);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("At all mesh points:  %g %g %g   %g %g %g\n",
           PREYIN,PREYIN,PREYIN, PREDIN,PREDIN,PREDIN);
#else
    printf("At all mesh points:  %g %g %g   %g %g %g\n", PREYIN,PREYIN,PREYIN,
           PREDIN,PREDIN,PREDIN);
#endif
}

/*
 * Print sample of current cc values 
 */

static void PrintOutput(int my_pe, MPI_Comm comm, N_Vector cc)
{
  int is, i0, npelast;
  realtype  *ct, tempc[NUM_SPECIES];
  MPI_Status status;
  
  npelast = NPEX*NPEY - 1;
  
  ct = N_VGetArrayPointer_Parallel(cc);
  
  /* Send the cc values (for all species) at the top right mesh point to PE 0 */
  if (my_pe == npelast) {
    i0 = NUM_SPECIES*(MXSUB*MYSUB-1);
    if (npelast!=0)
      MPI_Send(&ct[i0],NUM_SPECIES,PVEC_REAL_MPI_TYPE,0,0,comm);
    else  /* single processor case */
      for (is = 0; is < NUM_SPECIES; is++) tempc[is]=ct[i0+is];   
  }
  
  /* On PE 0, receive the cc values at top right, then print performance data 
     and sampled solution values */
  if (my_pe == 0) {
    
    if (npelast != 0)
      MPI_Recv(&tempc[0],NUM_SPECIES,PVEC_REAL_MPI_TYPE,npelast,0,comm,&status);
    
    printf("\nAt bottom left:");
    for (is = 0; is < NUM_SPECIES; is++){
      if ((is%6)*6== is) printf("\n");
#if defined(SUNDIALS_EXTENDED_PRECISION) 
      printf(" %Lg",ct[is]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf(" %g",ct[is]);
#else
      printf(" %g",ct[is]);
#endif
    }
    
    printf("\n\nAt top right:");
    for (is = 0; is < NUM_SPECIES; is++) {
      if ((is%6)*6 == is) printf("\n");
#if defined(SUNDIALS_EXTENDED_PRECISION) 
      printf(" %Lg",tempc[is]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      printf(" %g",tempc[is]);
#else
      printf(" %g",tempc[is]);
#endif
    }
    printf("\n\n");
  } 
}

/*
 * Print final statistics contained in iopt 
 */

static void PrintFinalStats(void *kmem)
{
  long int nni, nfe, nli, npe, nps, ncfl, nfeSG;
  int flag;
  
  flag = KINGetNumNonlinSolvIters(kmem, &nni);
  check_flag(&flag, "KINGetNumNonlinSolvIters", 1, 0);
  flag = KINGetNumFuncEvals(kmem, &nfe);
  check_flag(&flag, "KINGetNumFuncEvals", 1, 0);
  flag = KINGetNumLinIters(kmem, &nli);
  check_flag(&flag, "KINGetNumLinIters", 1, 0);
  flag = KINGetNumPrecEvals(kmem, &npe);
  check_flag(&flag, "KINGetNumPrecEvals", 1, 0);
  flag = KINGetNumPrecSolves(kmem, &nps);
  check_flag(&flag, "KINGetNumPrecSolves", 1, 0);
  flag = KINGetNumLinConvFails(kmem, &ncfl);
  check_flag(&flag, "KINGetNumLinConvFails", 1, 0);
  flag = KINGetNumLinFuncEvals(kmem, &nfeSG);
  check_flag(&flag, "KINGetNumLinFuncEvals", 1, 0);

  printf("Final Statistics.. \n");
  printf("nni    = %5ld    nli   = %5ld\n", nni, nli);
  printf("nfe    = %5ld    nfeSG = %5ld\n", nfe, nfeSG);
  printf("nps    = %5ld    npe   = %5ld     ncfl  = %5ld\n", nps, npe, ncfl);
  
}

/*
 * Routine to send boundary data to neighboring PEs 
 */

static void BSend(MPI_Comm comm, int my_pe, 
                  int isubx, int isuby, 
                  int dsizex, int dsizey, realtype *cdata)
{
  int i, ly;
  int offsetc, offsetbuf;
  realtype bufleft[NUM_SPECIES*MYSUB], bufright[NUM_SPECIES*MYSUB];

  /* If isuby > 0, send data from bottom x-line of u */
  if (isuby != 0)
    MPI_Send(&cdata[0], dsizex, PVEC_REAL_MPI_TYPE, my_pe-NPEX, 0, comm);

  /* If isuby < NPEY-1, send data from top x-line of u */
  if (isuby != NPEY-1) {
    offsetc = (MYSUB-1)*dsizex;
    MPI_Send(&cdata[offsetc], dsizex, PVEC_REAL_MPI_TYPE, my_pe+NPEX, 0, comm);
  }

  /* If isubx > 0, send data from left y-line of u (via bufleft) */
  if (isubx != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetc = ly*dsizex;
      for (i = 0; i < NUM_SPECIES; i++)
        bufleft[offsetbuf+i] = cdata[offsetc+i];
    }
    MPI_Send(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe-1, 0, comm);   
  }

  /* If isubx < NPEX-1, send data from right y-line of u (via bufright) */
  if (isubx != NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetc = offsetbuf*MXSUB + (MXSUB-1)*NUM_SPECIES;
      for (i = 0; i < NUM_SPECIES; i++)
        bufright[offsetbuf+i] = cdata[offsetc+i];
    }
    MPI_Send(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe+1, 0, comm);   
  }
}

/*
 * Routine to start receiving boundary data from neighboring PEs.
 *  Notes:
 *  1) buffer should be able to hold 2*NUM_SPECIES*MYSUB realtype entries,
 *     should be passed to both the BRecvPost and BRecvWait functions, and
 *     should not be manipulated between the two calls.
 *  2) request should have 4 entries, and should be passed in both calls also. 
 */

static void BRecvPost(MPI_Comm comm, MPI_Request request[], int my_pe,
                      int isubx, int isuby,
                      int dsizex, int dsizey,
                      realtype *cext, realtype *buffer)
{
  int offsetce;

  /* Have bufleft and bufright use the same buffer */
  realtype *bufleft = buffer, *bufright = buffer+NUM_SPECIES*MYSUB;
  
  /* If isuby > 0, receive data for bottom x-line of cext */
  if (isuby != 0)
    MPI_Irecv(&cext[NUM_SPECIES], dsizex, PVEC_REAL_MPI_TYPE,
              my_pe-NPEX, 0, comm, &request[0]);
  
  /* If isuby < NPEY-1, receive data for top x-line of cext */
  if (isuby != NPEY-1) {
    offsetce = NUM_SPECIES*(1 + (MYSUB+1)*(MXSUB+2));
    MPI_Irecv(&cext[offsetce], dsizex, PVEC_REAL_MPI_TYPE,
              my_pe+NPEX, 0, comm, &request[1]);
  }
  
  /* If isubx > 0, receive data for left y-line of cext (via bufleft) */
  if (isubx != 0) {
    MPI_Irecv(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE,
              my_pe-1, 0, comm, &request[2]);
  }
  
  /* If isubx < NPEX-1, receive data for right y-line of cext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Irecv(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE,
              my_pe+1, 0, comm, &request[3]);
  } 
}

/*
 * Routine to finish receiving boundary data from neighboring PEs.
 *  Notes:
 *  1) buffer should be able to hold 2*NUM_SPECIES*MYSUB realtype entries,
 *     should be passed to both the BRecvPost and BRecvWait functions, and
 *     should not be manipulated between the two calls.
 *  2) request should have 4 entries, and should be passed in both calls also. 
 */

static void BRecvWait(MPI_Request request[], int isubx,
                      int isuby, int dsizex, realtype *cext,
                      realtype *buffer)
{
  int i, ly;
  int dsizex2, offsetce, offsetbuf;
  realtype *bufleft = buffer, *bufright = buffer+NUM_SPECIES*MYSUB;
  MPI_Status status;
  
  dsizex2 = dsizex + 2*NUM_SPECIES;
  
  /* If isuby > 0, receive data for bottom x-line of cext */
  if (isuby != 0)
    MPI_Wait(&request[0],&status);
  
  /* If isuby < NPEY-1, receive data for top x-line of cext */
  if (isuby != NPEY-1)
    MPI_Wait(&request[1],&status);

  /* If isubx > 0, receive data for left y-line of cext (via bufleft) */
  if (isubx != 0) {
    MPI_Wait(&request[2],&status);
    
    /* Copy the buffer to cext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetce = (ly+1)*dsizex2;
      for (i = 0; i < NUM_SPECIES; i++)
        cext[offsetce+i] = bufleft[offsetbuf+i];
    }
  }
  
  /* If isubx < NPEX-1, receive data for right y-line of cext (via bufright) */
  if (isubx != NPEX-1) {
    MPI_Wait(&request[3],&status);
    
    /* Copy the buffer to cext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetce = (ly+2)*dsizex2 - NUM_SPECIES;
      for (i = 0; i < NUM_SPECIES; i++)
        cext[offsetce+i] = bufright[offsetbuf+i];
    }
  }  
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

static int check_flag(void *flagvalue, const char *funcname, int opt, int id)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n",
	    id, funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR(%d): %s() failed with flag = %d\n\n",
	      id, funcname, *errflag);
      return(1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n",
	    id, funcname);
    return(1);
  }

  return(0);
}
