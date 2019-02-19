/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *         Allan Taylor, Alan Hindmarsh and Radu Serban @ LLNL
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
 * Example program for IDA: Food web, parallel, GMRES, IDABBD
 * preconditioner.
 *
 * This example program for IDA uses SUNLinSol_SPGMR as the linear solver.
 * It is written for a parallel computer system and uses the
 * IDABBDPRE band-block-diagonal preconditioner module for the
 * SUNLinSol_SPGMR package. 
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
 *   0 = d(i)*(c    + c  )  +  R  (x,y,c)   (i = np+1,...,ns)
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
 *   a(i,i) = -AA  (all i)
 *   a(i,j) = -GG  (i <= np , j >  np)
 *   a(i,j) =  EE  (i >  np,  j <= np)
 *   all other a(i,j) = 0
 *   b(i) = BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y))  (i <= np)
 *   b(i) =-BB*(1+ alpha * x*y + beta*sin(4 pi x)*sin(4 pi y))  (i  > np)
 *   d(i) = DPREY  (i <= np)
 *   d(i) = DPRED  (i > np)
 *
 * Note: The above equations are written in 1-based indices,
 * whereas the code has 0-based indices, being written in C.
 *
 * The various scalar parameters required are set using '#define'
 * statements or directly in routine InitUserData. In this program,
 * np = 1, ns = 2. The boundary conditions are homogeneous Neumann:
 * normal derivative  =  0.
 *
 * A polynomial in x and y is used to set the initial values of the
 * first np variables (the prey variables) at each x,y location,
 * while initial values for the remaining (predator) variables are
 * set to a flat value, which is corrected by IDACalcIC.
 *
 * The PDEs are discretized by central differencing on a MX by MY
 * mesh, and so the system size Neq is the product
 * MX * MY * NUM_SPECIES. The system is actually implemented on
 * submeshes, processor by processor, with an MXSUB by MYSUB mesh
 * on each of NPEX * NPEY processors.
 *
 * The DAE system is solved by IDA using the SUNLinSol_SPGMR linear solver,
 * in conjunction with the preconditioner module IDABBDPRE. The
 * preconditioner uses a 5-diagonal band-block-diagonal
 * approximation (half-bandwidths = 2). Output is printed at
 * t = 0, .001, .01, .1, .4, .7, 1.
 * -----------------------------------------------------------------
 * References:
 * [1] Peter N. Brown and Alan C. Hindmarsh,
 *     Reduced Storage Matrix Methods in Stiff ODE systems,
 *     Journal of Applied Mathematics and Computation, Vol. 31
 *     (May 1989), pp. 40-91.
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

#include <ida/ida.h>
#include <ida/ida_bbdpre.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_math.h>

#include <mpi.h>

/* Problem Constants */

#define NPREY       1        /* Number of prey (= number of predators). */
#define NUM_SPECIES 2*NPREY

#define PI          RCONST(3.1415926535898) /* pi */ 
#define FOURPI      (RCONST(4.0)*PI)        /* 4 pi */

#define MXSUB       10    /* Number of x mesh points per processor subgrid */
#define MYSUB       10    /* Number of y mesh points per processor subgrid */
#define NPEX        2     /* Number of subgrids in the x direction */
#define NPEY        2     /* Number of subgrids in the y direction */
#define MX          (MXSUB*NPEX)      /* MX = number of x mesh points */
#define MY          (MYSUB*NPEY)      /* MY = number of y mesh points */
#define NSMXSUB     (NUM_SPECIES * MXSUB)
#define NEQ         (NUM_SPECIES*MX*MY) /* Number of equations in system */
#define AA          RCONST(1.0)    /* Coefficient in above eqns. for a */
#define EE          RCONST(10000.) /* Coefficient in above eqns. for a */
#define GG          RCONST(0.5e-6) /* Coefficient in above eqns. for a */
#define BB          RCONST(1.0)    /* Coefficient in above eqns. for b */
#define DPREY       RCONST(1.0)    /* Coefficient in above eqns. for d */
#define DPRED       RCONST(0.05)   /* Coefficient in above eqns. for d */
#define ALPHA       RCONST(50.)    /* Coefficient alpha in above eqns. */
#define BETA        RCONST(1000.)  /* Coefficient beta in above eqns. */
#define AX          RCONST(1.0)    /* Total range of x variable */
#define AY          RCONST(1.0)    /* Total range of y variable */
#define RTOL        RCONST(1.e-5)  /*  rtol tolerance */
#define ATOL        RCONST(1.e-5)  /*  atol tolerance */
#define ZERO        RCONST(0.)     /* 0. */
#define ONE         RCONST(1.0)    /* 1. */
#define NOUT        6  
#define TMULT       RCONST(10.0)   /* Multiplier for tout values */
#define TADD        RCONST(0.3)    /* Increment for tout values */

/* User-defined vector accessor macro IJ_Vptr. */

/*
 * IJ_Vptr is defined in order to express the underlying 3-d structure of the 
 * dependent variable vector from its underlying 1-d storage (an N_Vector).
 * IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to 
 * species index is = 0, x-index ix = i, and y-index jy = j.                
 */

#define IJ_Vptr(vv,i,j) (&NV_Ith_P(vv, (i)*NUM_SPECIES + (j)*NSMXSUB ))

/* Type: UserData.  Contains problem constants, preconditioner data, etc. */

typedef struct {
  int ns, np, thispe, npes, ixsub, jysub, npex, npey;
  int mxsub, mysub, nsmxsub, nsmxsub2;
  realtype dx, dy, **acoef;
  realtype cox[NUM_SPECIES], coy[NUM_SPECIES], bcoef[NUM_SPECIES],
    rhs[NUM_SPECIES], cext[(MXSUB+2)*(MYSUB+2)*NUM_SPECIES];
  MPI_Comm comm;
  N_Vector rates;
  sunindextype n_local;
} *UserData;

/* Prototypes for functions called by the IDA Solver. */

static int resweb(realtype tt, N_Vector cc, N_Vector cp,
                  N_Vector rr, void *user_data);

static int reslocal(sunindextype Nlocal, realtype tt, 
                    N_Vector cc, N_Vector cp, N_Vector res, 
                    void *user_data);

static int rescomm(sunindextype Nlocal, realtype tt,
                   N_Vector cc, N_Vector cp, 
                   void *user_data);

/* Prototypes for supporting functions */

static void BSend(MPI_Comm comm, int thispe, int ixsub, int jysub,
                  int dsizex, int dsizey, realtype carray[]);

static void BRecvPost(MPI_Comm comm, MPI_Request request[], int thispe,
                      int ixsub, int jysub,
                      int dsizex, int dsizey,
                      realtype cext[], realtype buffer[]);

static void BRecvWait(MPI_Request request[], int ixsub, int jysub,
                      int dsizex, realtype cext[], realtype buffer[]);

static void WebRates(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy, 
                     UserData webdata);

static realtype dotprod(int size, realtype *x1, realtype *x2);

/* Prototypes for private functions */

static void InitUserData(UserData webdata, int thispe, int npes, 
                         MPI_Comm comm);

static void SetInitialProfiles(N_Vector cc, N_Vector cp, N_Vector id,
                               N_Vector scrtch, UserData webdata);

static void PrintHeader(sunindextype SystemSize, int maxl, 
                        sunindextype mudq, sunindextype mldq, 
                        sunindextype mukeep, sunindextype mlkeep,
                        realtype rtol, realtype atol);

static void PrintOutput(void *ida_mem, N_Vector cc, realtype time,
                        UserData webdata, MPI_Comm comm);

static void PrintFinalStats(void *ida_mem);

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  MPI_Comm comm;
  void *ida_mem;
  SUNLinearSolver LS;
  UserData webdata;
  sunindextype SystemSize, local_N, mudq, mldq, mukeep, mlkeep;
  realtype rtol, atol, t0, tout, tret;
  N_Vector cc, cp, res, id;
  int thispe, npes, maxl, iout, retval;

  cc = cp = res = id = NULL;
  webdata = NULL;
  LS = NULL;
  ida_mem = NULL;

  /* Set communicator, and get processor number and total number of PE's. */

  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &thispe);
  MPI_Comm_size(comm, &npes);

  if (npes != NPEX*NPEY) {
    if (thispe == 0)
      fprintf(stderr, 
              "\nMPI_ERROR(0): npes = %d not equal to NPEX*NPEY = %d\n", 
              npes, NPEX*NPEY);
    MPI_Finalize();
    return(1); 
  }
  
  /* Set local length (local_N) and global length (SystemSize). */

  local_N = MXSUB*MYSUB*NUM_SPECIES;
  SystemSize = NEQ;

  /* Set up user data block webdata. */

  webdata = (UserData) malloc(sizeof *webdata);
  webdata->rates = N_VNew_Parallel(comm, local_N, SystemSize);
  webdata->acoef = newDenseMat(NUM_SPECIES, NUM_SPECIES);

  InitUserData(webdata, thispe, npes, comm);
  
  /* Create needed vectors, and load initial values.
     The vector res is used temporarily only.        */
  
  cc  = N_VNew_Parallel(comm, local_N, SystemSize);
  if(check_retval((void *)cc, "N_VNew_Parallel", 0, thispe)) MPI_Abort(comm, 1);

  cp  = N_VNew_Parallel(comm, local_N, SystemSize);
  if(check_retval((void *)cp, "N_VNew_Parallel", 0, thispe)) MPI_Abort(comm, 1);

  res = N_VNew_Parallel(comm, local_N, SystemSize);
  if(check_retval((void *)res, "N_VNew_Parallel", 0, thispe)) MPI_Abort(comm, 1);

  id  = N_VNew_Parallel(comm, local_N, SystemSize);
  if(check_retval((void *)id, "N_VNew_Parallel", 0, thispe)) MPI_Abort(comm, 1);
  
  SetInitialProfiles(cc, cp, id, res, webdata);
  
  N_VDestroy_Parallel(res);
  
  /* Set remaining inputs to IDAMalloc. */
  
  t0 = ZERO;
  rtol = RTOL; 
  atol = ATOL;
  
  /* Call IDACreate and IDAMalloc to initialize solution */

  ida_mem = IDACreate();
  if(check_retval((void *)ida_mem, "IDACreate", 0, thispe)) MPI_Abort(comm, 1);

  retval = IDASetUserData(ida_mem, webdata);
  if(check_retval(&retval, "IDASetUserData", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetId(ida_mem, id);
  if(check_retval(&retval, "IDASetId", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDAInit(ida_mem, resweb, t0, cc, cp);
  if(check_retval(&retval, "IDAInit", 1, thispe)) MPI_Abort(comm, 1);
  
  retval = IDASStolerances(ida_mem, rtol, atol);
  if(check_retval(&retval, "IDASStolerances", 1, thispe)) MPI_Abort(comm, 1);

  /* Call SUNLinSol_SPGMR and IDASetLinearSolver to specify the linear solver */

  maxl = 16;
  LS = SUNLinSol_SPGMR(cc, PREC_LEFT, maxl);
  if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0, thispe)) MPI_Abort(comm, 1);
  retval = IDASetLinearSolver(ida_mem, LS, NULL);
  if(check_retval(&retval, "IDASetLinearSolver", 1, thispe)) MPI_Abort(comm, 1);

  /* Call IDABBDPrecInit to initialize the band-block-diagonal preconditioner.
     The half-bandwidths for the difference quotient evaluation are exact
     for the system Jacobian, but only a 5-diagonal band matrix is retained. */
  
  mudq = mldq = NSMXSUB;
  mukeep = mlkeep = 2;
  retval = IDABBDPrecInit(ida_mem, local_N, mudq, mldq, mukeep, mlkeep, 
                          ZERO, reslocal, NULL);
  if(check_retval(&retval, "IDABBDPrecInit", 1, thispe)) MPI_Abort(comm, 1);
  
  /* Call IDACalcIC (with default options) to correct the initial values. */
  
  tout = RCONST(0.001);
  retval = IDACalcIC(ida_mem, IDA_YA_YDP_INIT, tout);
  if(check_retval(&retval, "IDACalcIC", 1, thispe)) MPI_Abort(comm, 1);
  
  /* On PE 0, print heading, basic parameters, initial values. */
 
  if (thispe == 0) PrintHeader(SystemSize, maxl, 
                               mudq, mldq, mukeep, mlkeep,
                               rtol, atol);
  PrintOutput(ida_mem, cc, t0, webdata, comm);

  /* Call IDA in tout loop, normal mode, and print selected output. */
  
  for (iout = 1; iout <= NOUT; iout++) {
    
    retval = IDASolve(ida_mem, tout, &tret, cc, cp, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1, thispe)) MPI_Abort(comm, 1);
    
    PrintOutput(ida_mem, cc, tret, webdata, comm);
    
    if (iout < 3) tout *= TMULT; 
    else          tout += TADD;

  }
  
  /* On PE 0, print final set of statistics. */
  
  if (thispe == 0)  PrintFinalStats(ida_mem);

  /* Free memory. */

  N_VDestroy_Parallel(cc);
  N_VDestroy_Parallel(cp);
  N_VDestroy_Parallel(id);

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);

  destroyMat(webdata->acoef);
  N_VDestroy_Parallel(webdata->rates);
  free(webdata);

  MPI_Finalize();

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

static void InitUserData(UserData webdata, int thispe, int npes, 
                         MPI_Comm comm)
{
  int i, j, np;
  realtype *a1,*a2, *a3, *a4, dx2, dy2, **acoef, *bcoef, *cox, *coy;

  webdata->jysub = thispe / NPEX;
  webdata->ixsub = thispe - (webdata->jysub)*NPEX;
  webdata->mxsub = MXSUB;
  webdata->mysub = MYSUB;
  webdata->npex = NPEX;
  webdata->npey = NPEY;
  webdata->ns = NUM_SPECIES;
  webdata->np = NPREY;
  webdata->dx = AX/(MX-1);
  webdata->dy = AY/(MY-1);
  webdata->thispe = thispe;
  webdata->npes   = npes;
  webdata->nsmxsub = MXSUB * NUM_SPECIES;
  webdata->nsmxsub2 = (MXSUB+2)*NUM_SPECIES;
  webdata->comm = comm;
  webdata->n_local = MXSUB*MYSUB*NUM_SPECIES;

  /* Set up the coefficients a and b plus others found in the equations. */

  np = webdata->np;
  dx2 = (webdata->dx)*(webdata->dx); 
  dy2 = (webdata->dy)*(webdata->dy);

  acoef = webdata->acoef;
  bcoef = webdata->bcoef;
  cox = webdata->cox;
  coy = webdata->coy;

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
                               N_Vector res, UserData webdata)
{
  int ixsub, jysub, mxsub, mysub, np, ix, jy, is;
  realtype *cxy, *idxy, *cpxy, dx, dy, xx, yy, xyfactor;

  ixsub = webdata->ixsub;
  jysub = webdata->jysub;
  mxsub = webdata->mxsub;
  mysub = webdata->mxsub;
  dx = webdata->dx;
  dy = webdata->dy;
  np = webdata->np;

  /* Loop over grid, load cc values and id values. */
  for (jy = 0; jy < mysub; jy++) {
    yy = (jy + jysub*mysub) * dy;
    for (ix = 0; ix < mxsub; ix++) {
      xx = (ix + ixsub*mxsub) * dx;
      xyfactor = 16.*xx*(1. - xx)*yy*(1. - yy);
      xyfactor *= xyfactor;

      cxy = IJ_Vptr(cc,ix,jy); 
      idxy = IJ_Vptr(id,ix,jy); 
      for (is = 0; is < NUM_SPECIES; is++) {
	if (is < np)
           {cxy[is] = RCONST(10.0)+(realtype)(is+1)*xyfactor; idxy[is] = ONE;}
        else { cxy[is] = 1.0e5; idxy[is] = ZERO; }
      }
    }
  }

  /* Set c' for the prey by calling the residual function with cp = 0. */
  
  N_VConst(ZERO, cp);
  resweb(ZERO, cc, cp, res, webdata);
  N_VScale(-ONE, res, cp);
  
  /* Set c' for predators to 0. */
  
  for (jy = 0; jy < mysub; jy++) {
    for (ix = 0; ix < mxsub; ix++) {
      cpxy = IJ_Vptr(cp,ix,jy); 
      for (is = np; is < NUM_SPECIES; is++) cpxy[is] = ZERO;
    }
  }
}

/*
 * Print first lines of output (problem description)
 * and table headerr
 */

static void PrintHeader(sunindextype SystemSize, int maxl, 
                        sunindextype mudq, sunindextype mldq, 
                        sunindextype mukeep, sunindextype mlkeep,
                        realtype rtol, realtype atol)
{
  printf("\nidaFoodWeb_kry_bbd_p: Predator-prey DAE parallel example problem for IDA \n\n");
  printf("Number of species ns: %d", NUM_SPECIES);
  printf("     Mesh dimensions: %d x %d", MX, MY);
  printf("     Total system size: %ld\n",(long int) SystemSize);
  printf("Subgrid dimensions: %d x %d", MXSUB, MYSUB);
  printf("     Processor array: %d x %d\n", NPEX, NPEY);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Linear solver: SUNLinSol_SPGMR     Max. Krylov dimension maxl: %d\n", maxl);
  printf("Preconditioner: band-block-diagonal (IDABBDPRE), with parameters\n");
  printf("     mudq = %ld,  mldq = %ld,  mukeep = %ld,  mlkeep = %ld\n",
         (long int) mudq, (long int) mldq, (long int) mukeep, (long int) mlkeep);
  printf("CalcIC called to correct initial predator concentrations \n\n");
  printf("-----------------------------------------------------------\n");
  printf("  t        bottom-left  top-right");
  printf("    | nst  k      h\n");
  printf("-----------------------------------------------------------\n\n");
}


/*
 * PrintOutput: Print output values at output time t = tt.
 * Selected run statistics are printed.  Then values of c1 and c2
 * are printed for the bottom left and top right grid points only.
 */

static void PrintOutput(void *ida_mem, N_Vector cc, realtype tt,
                        UserData webdata, MPI_Comm comm)
{
  MPI_Status status;
  realtype *cdata, clast[2], hused;
  long int nst;
  int i, kused, retval, thispe, npelast, ilast;;

  thispe = webdata->thispe; 
  npelast = webdata->npes - 1;
  cdata = N_VGetArrayPointer_Parallel(cc);
  
  /* Send conc. at top right mesh point from PE npes-1 to PE 0. */
  if (thispe == npelast) {
    ilast = NUM_SPECIES*MXSUB*MYSUB - 2;
    if (npelast != 0)
      MPI_Send(&cdata[ilast], 2, PVEC_REAL_MPI_TYPE, 0, 0, comm);
    else { clast[0] = cdata[ilast]; clast[1] = cdata[ilast+1]; }
  }
  
  /* On PE 0, receive conc. at top right from PE npes - 1.
     Then print performance data and sampled solution values. */
  
  if (thispe == 0) {
    
    if (npelast != 0)
      MPI_Recv(&clast[0], 2, PVEC_REAL_MPI_TYPE, npelast, 0, comm, &status);
    
    retval = IDAGetLastOrder(ida_mem, &kused);
    check_retval(&retval, "IDAGetLastOrder", 1, thispe);
    retval = IDAGetNumSteps(ida_mem, &nst);
    check_retval(&retval, "IDAGetNumSteps", 1, thispe);
    retval = IDAGetLastStep(ida_mem, &hused);
    check_retval(&retval, "IDAGetLastStep", 1, thispe);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%8.2Le %12.4Le %12.4Le   | %3ld  %1d %12.4Le\n", 
         tt, cdata[0], clast[0], nst, kused, hused);
    for (i=1;i<NUM_SPECIES;i++)
      printf("         %12.4Le %12.4Le   |\n",cdata[i],clast[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%8.2e %12.4e %12.4e   | %3ld  %1d %12.4e\n", 
         tt, cdata[0], clast[0], nst, kused, hused);
    for (i=1;i<NUM_SPECIES;i++)
      printf("         %12.4e %12.4e   |\n",cdata[i],clast[i]);
#else
    printf("%8.2e %12.4e %12.4e   | %3ld  %1d %12.4e\n", 
         tt, cdata[0], clast[0], nst, kused, hused);
    for (i=1;i<NUM_SPECIES;i++)
      printf("         %12.4e %12.4e   |\n",cdata[i],clast[i]);
#endif
    printf("\n");

  }

}

/*
 * PrintFinalStats: Print final run data contained in iopt.              
 */

static void PrintFinalStats(void *ida_mem)
{
  long int nst, nre, nreLS, netf, ncfn, nni, ncfl, nli, npe, nps, nge;
  int retval;

  retval = IDAGetNumSteps(ida_mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1, 0);
  retval = IDAGetNumResEvals(ida_mem, &nre);
  check_retval(&retval, "IDAGetNumResEvals", 1, 0);
  retval = IDAGetNumErrTestFails(ida_mem, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1, 0);
  retval = IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1, 0);
  retval = IDAGetNumNonlinSolvIters(ida_mem, &nni);
  check_retval(&retval, "IDAGetNumNonlinSolvIters", 1, 0);

  retval = IDAGetNumLinConvFails(ida_mem, &ncfl);
  check_retval(&retval, "IDAGetNumLinConvFails", 1, 0);
  retval = IDAGetNumLinIters(ida_mem, &nli);
  check_retval(&retval, "IDAGetNumLinIters", 1, 0);
  retval = IDAGetNumPrecEvals(ida_mem, &npe);
  check_retval(&retval, "IDAGetNumPrecEvals", 1, 0);
  retval = IDAGetNumPrecSolves(ida_mem, &nps);
  check_retval(&retval, "IDAGetNumPrecSolves", 1, 0);
  retval = IDAGetNumLinResEvals(ida_mem, &nreLS);
  check_retval(&retval, "IDAGetNumLinResEvals", 1, 0);

  retval = IDABBDPrecGetNumGfnEvals(ida_mem, &nge);
  check_retval(&retval, "IDABBDPrecGetNumGfnEvals", 1, 0);

  printf("-----------------------------------------------------------\n");
  printf("\nFinal statistics: \n\n");

  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n\n", ncfn);

  printf("Number of linear iterations        = %ld\n", nli);
  printf("Number of linear conv. failures    = %ld\n\n", ncfl);

  printf("Number of preconditioner setups    = %ld\n", npe);
  printf("Number of preconditioner solves    = %ld\n", nps);
  printf("Number of local residual evals.    = %ld\n", nge);

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

static int check_retval(void *returnvalue, const char *funcname, int opt, int id)
{
  int *retval;

  if (opt == 0 && returnvalue == NULL) {
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n", 
            id, funcname);
    return(1); 
  } else if (opt == 1) {
    /* Check if retval < 0 */
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR(%d): %s() failed with retval = %d\n\n", 
              id, funcname, *retval);
      return(1); 
    }
  } else if (opt == 2 && returnvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n", 
            id, funcname);
    return(1); 
  }

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA & SUPPORTING FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * resweb: System residual function for predator-prey system.
 * To compute the residual function F, this routine calls:
 * rescomm, for needed communication, and then
 * reslocal, for computation of the residuals on this processor.      
 */

static int resweb(realtype tt, 
                  N_Vector cc, N_Vector cp, N_Vector rr, 
                  void *user_data)
{
  int retval;
  UserData webdata;
  sunindextype Nlocal;
  
  webdata = (UserData) user_data;
  
  Nlocal = webdata->n_local;

  /* Call rescomm to do inter-processor communication. */
  retval = rescomm(Nlocal, tt, cc, cp, user_data);
  
  /* Call reslocal to calculate the local portion of residual vector. */
  retval = reslocal(Nlocal, tt, cc, cp, rr, user_data);
  
  return(retval);
}

/*
 * rescomm: Communication routine in support of resweb.
 * This routine performs all inter-processor communication of components
 * of the cc vector needed to calculate F, namely the components at all
 * interior subgrid boundaries (ghost cell data).  It loads this data
 * into a work array cext (the local portion of c, extended).
 * The message-passing uses blocking sends, non-blocking receives,
 * and receive-waiting, in routines BRecvPost, BSend, BRecvWait.         
 */

static int rescomm(sunindextype Nlocal, realtype tt, 
                   N_Vector cc, N_Vector cp,
                   void *user_data)
{

  UserData webdata;
  realtype *cdata, *cext, buffer[2*NUM_SPECIES*MYSUB];
  int thispe, ixsub, jysub, nsmxsub, nsmysub;
  MPI_Comm comm;
  MPI_Request request[4];
  
  webdata = (UserData) user_data;
  cdata = N_VGetArrayPointer_Parallel(cc);
  
  /* Get comm, thispe, subgrid indices, data sizes, extended array cext. */
  
  comm = webdata->comm;     
  thispe = webdata->thispe;

  ixsub = webdata->ixsub;   
  jysub = webdata->jysub;
  cext = webdata->cext;
  nsmxsub = webdata->nsmxsub; 
  nsmysub = (webdata->ns)*(webdata->mysub);
  
  /* Start receiving boundary data from neighboring PEs. */

  BRecvPost(comm, request, thispe, ixsub, jysub, nsmxsub, nsmysub, 
            cext, buffer);
  
  /* Send data from boundary of local grid to neighboring PEs. */
  
  BSend(comm, thispe, ixsub, jysub, nsmxsub, nsmysub, cdata);
  
  /* Finish receiving boundary data from neighboring PEs. */
  
  BRecvWait(request, ixsub, jysub, nsmxsub, cext, buffer);
  
  return(0);
}

/*
 * BRecvPost: Start receiving boundary data from neighboring PEs.
 * (1) buffer should be able to hold 2*NUM_SPECIES*MYSUB realtype entries,
 *     should be passed to both the BRecvPost and BRecvWait functions, and
 *     should not be manipulated between the two calls.
 * (2) request should have 4 entries, and is also passed in both calls. 
 */

static void BRecvPost(MPI_Comm comm, MPI_Request request[], int my_pe,
                      int ixsub, int jysub,
                      int dsizex, int dsizey,
                      realtype cext[], realtype buffer[])
{
  int offsetce;
  /* Have bufleft and bufright use the same buffer. */
  realtype *bufleft = buffer, *bufright = buffer+NUM_SPECIES*MYSUB;

  /* If jysub > 0, receive data for bottom x-line of cext. */
  if (jysub != 0)
    MPI_Irecv(&cext[NUM_SPECIES], dsizex, PVEC_REAL_MPI_TYPE,
              my_pe-NPEX, 0, comm, &request[0]);
  
  /* If jysub < NPEY-1, receive data for top x-line of cext. */
  if (jysub != NPEY-1) {
    offsetce = NUM_SPECIES*(1 + (MYSUB+1)*(MXSUB+2));
    MPI_Irecv(&cext[offsetce], dsizex, PVEC_REAL_MPI_TYPE,
              my_pe+NPEX, 0, comm, &request[1]);
  }
  
  /* If ixsub > 0, receive data for left y-line of cext (via bufleft). */
  if (ixsub != 0) {
    MPI_Irecv(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE,
              my_pe-1, 0, comm, &request[2]);
  }
  
  /* If ixsub < NPEX-1, receive data for right y-line of cext (via bufright). */
  if (ixsub != NPEX-1) {
    MPI_Irecv(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE,
              my_pe+1, 0, comm, &request[3]);
  }
  
}

/*
 * BRecvWait: Finish receiving boundary data from neighboring PEs.
 * (1) buffer should be able to hold 2*NUM_SPECIES*MYSUB realtype entries,
 *     should be passed to both the BRecvPost and BRecvWait functions, and
 *     should not be manipulated between the two calls.
 * (2) request should have 4 entries, and is also passed in both calls.  
 */

static void BRecvWait(MPI_Request request[], int ixsub, int jysub,
                      int dsizex, realtype cext[], realtype buffer[])
{
  int i;
  int ly, dsizex2, offsetce, offsetbuf;
  realtype *bufleft = buffer, *bufright = buffer+NUM_SPECIES*MYSUB;
  MPI_Status status;
  
  dsizex2 = dsizex + 2*NUM_SPECIES;
  
  /* If jysub > 0, receive data for bottom x-line of cext. */
  if (jysub != 0)
    MPI_Wait(&request[0],&status);
  
  /* If jysub < NPEY-1, receive data for top x-line of cext. */
  if (jysub != NPEY-1)
    MPI_Wait(&request[1],&status);

  /* If ixsub > 0, receive data for left y-line of cext (via bufleft). */
  if (ixsub != 0) {
    MPI_Wait(&request[2],&status);

    /* Copy the buffer to cext */
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetce = (ly+1)*dsizex2;
      for (i = 0; i < NUM_SPECIES; i++)
        cext[offsetce+i] = bufleft[offsetbuf+i];
    }
  }
  
  /* If ixsub < NPEX-1, receive data for right y-line of cext (via bufright). */
  if (ixsub != NPEX-1) {
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
 * BSend: Send boundary data to neighboring PEs.
 * This routine sends components of cc from internal subgrid boundaries
 * to the appropriate neighbor PEs.                                      
 */
 
static void BSend(MPI_Comm comm, int my_pe, int ixsub, int jysub,
                  int dsizex, int dsizey, realtype cdata[])
{
  int i;
  int ly, offsetc, offsetbuf;
  realtype bufleft[NUM_SPECIES*MYSUB], bufright[NUM_SPECIES*MYSUB];

  /* If jysub > 0, send data from bottom x-line of cc. */

  if (jysub != 0)
    MPI_Send(&cdata[0], dsizex, PVEC_REAL_MPI_TYPE, my_pe-NPEX, 0, comm);

  /* If jysub < NPEY-1, send data from top x-line of cc. */

  if (jysub != NPEY-1) {
    offsetc = (MYSUB-1)*dsizex;
    MPI_Send(&cdata[offsetc], dsizex, PVEC_REAL_MPI_TYPE, my_pe+NPEX, 0, comm);
  }

  /* If ixsub > 0, send data from left y-line of cc (via bufleft). */

  if (ixsub != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetc = ly*dsizex;
      for (i = 0; i < NUM_SPECIES; i++)
        bufleft[offsetbuf+i] = cdata[offsetc+i];
    }
    MPI_Send(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe-1, 0, comm);   
  }
  
  /* If ixsub < NPEX-1, send data from right y-line of cc (via bufright). */
  
  if (ixsub != NPEX-1) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetc = offsetbuf*MXSUB + (MXSUB-1)*NUM_SPECIES;
      for (i = 0; i < NUM_SPECIES; i++)
        bufright[offsetbuf+i] = cdata[offsetc+i];
    }
    MPI_Send(&bufright[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe+1, 0, comm);   
  }
}
 
/* Define lines are for ease of readability in the following functions. */

#define mxsub      (webdata->mxsub)
#define mysub      (webdata->mysub)
#define npex       (webdata->npex)
#define npey       (webdata->npey)
#define ixsub      (webdata->ixsub)
#define jysub      (webdata->jysub)
#define nsmxsub    (webdata->nsmxsub)
#define nsmxsub2   (webdata->nsmxsub2)
#define np         (webdata->np)
#define dx         (webdata->dx)
#define dy         (webdata->dy)
#define cox        (webdata->cox)
#define coy        (webdata->coy)
#define rhs        (webdata->rhs)
#define cext       (webdata->cext)
#define rates      (webdata->rates)
#define ns         (webdata->ns)
#define acoef      (webdata->acoef)
#define bcoef      (webdata->bcoef)

/*
 * reslocal: Compute res = F(t,cc,cp).
 * This routine assumes that all inter-processor communication of data
 * needed to calculate F has already been done.  Components at interior
 * subgrid boundaries are assumed to be in the work array cext.
 * The local portion of the cc vector is first copied into cext.
 * The exterior Neumann boundary conditions are explicitly handled here
 * by copying data from the first interior mesh line to the ghost cell
 * locations in cext.  Then the reaction and diffusion terms are
 * evaluated in terms of the cext array, and the residuals are formed.
 * The reaction terms are saved separately in the vector webdata->rates
 * for use by the preconditioner setup routine.                          
 */

static int reslocal(sunindextype Nlocal, realtype tt, 
                    N_Vector cc, N_Vector cp, N_Vector rr,
                    void *user_data)
{
  realtype *cdata, *ratesxy, *cpxy, *resxy,
    xx, yy, dcyli, dcyui, dcxli, dcxui;
  int ix, jy, is, i, locc, ylocce, locce;
  UserData webdata;
  
  webdata = (UserData) user_data;
  
  /* Get data pointers, subgrid data, array sizes, work array cext. */
  
  cdata = N_VGetArrayPointer_Parallel(cc);
  
  /* Copy local segment of cc vector into the working extended array cext. */
  
  locc = 0;
  locce = nsmxsub2 + NUM_SPECIES;
  for (jy = 0; jy < mysub; jy++) {
    for (i = 0; i < nsmxsub; i++) cext[locce+i] = cdata[locc+i];
    locc = locc + nsmxsub;
    locce = locce + nsmxsub2;
  }

  /* To facilitate homogeneous Neumann boundary conditions, when this is
     a boundary PE, copy data from the first interior mesh line of cc to cext. */
  
  /* If jysub = 0, copy x-line 2 of cc to cext. */
  if (jysub == 0)
    { for (i = 0; i < nsmxsub; i++) cext[NUM_SPECIES+i] = cdata[nsmxsub+i]; }
  
  /* If jysub = npey-1, copy x-line mysub-1 of cc to cext. */
  if (jysub == npey-1) {
    locc = (mysub-2)*nsmxsub;
    locce = (mysub+1)*nsmxsub2 + NUM_SPECIES;
    for (i = 0; i < nsmxsub; i++) cext[locce+i] = cdata[locc+i];
  }
  
  /* If ixsub = 0, copy y-line 2 of cc to cext. */
  if (ixsub == 0) {
    for (jy = 0; jy < mysub; jy++) {
      locc = jy*nsmxsub + NUM_SPECIES;
      locce = (jy+1)*nsmxsub2;
      for (i = 0; i < NUM_SPECIES; i++) cext[locce+i] = cdata[locc+i];
    }
  }
  
  /* If ixsub = npex-1, copy y-line mxsub-1 of cc to cext. */
  if (ixsub == npex-1) {
    for (jy = 0; jy < mysub; jy++) {
      locc  = (jy+1)*nsmxsub - 2*NUM_SPECIES;
      locce = (jy+2)*nsmxsub2 - NUM_SPECIES;
      for (i = 0; i < NUM_SPECIES; i++) cext[locce+i] = cdata[locc+i];
    }
  }

  /* Loop over all grid points, setting local array rates to right-hand sides.
     Then set rr values appropriately for prey/predator components of F. */

  for (jy = 0; jy < mysub; jy++) {
    ylocce = (jy+1)*nsmxsub2;
    yy     = (jy+jysub*mysub)*dy;

    for (ix = 0; ix < mxsub; ix++) {
      locce = ylocce + (ix+1)*NUM_SPECIES;
      xx = (ix + ixsub*mxsub)*dx;

      ratesxy = IJ_Vptr(rates,ix,jy);
      WebRates(xx, yy, &(cext[locce]), ratesxy, webdata);

      resxy = IJ_Vptr(rr,ix,jy); 
      cpxy = IJ_Vptr(cp,ix,jy); 
      
      for (is = 0; is < NUM_SPECIES; is++) {
        dcyli = cext[locce+is]          - cext[locce+is-nsmxsub2];
        dcyui = cext[locce+is+nsmxsub2] - cext[locce+is];

        dcxli = cext[locce+is]             - cext[locce+is-NUM_SPECIES];
        dcxui = cext[locce+is+NUM_SPECIES] - cext[locce+is];

        rhs[is] = cox[is]*(dcxui-dcxli) + coy[is]*(dcyui-dcyli) + ratesxy[is];

        if (is < np) resxy[is] = cpxy[is] - rhs[is];
        else         resxy[is] =          - rhs[is];

      }
    }
  }
  
  return(0);
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

static realtype dotprod(int size, realtype *x1, realtype *x2)
{
  int i;
  realtype *xx1, *xx2, temp = ZERO;
  
  xx1 = x1; 
  xx2 = x2;
  for (i = 0; i < size; i++) 
    temp += (*xx1++) * (*xx2++);

  return(temp);
}

