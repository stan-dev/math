/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Cosmin Petra and Radu Serban @ LLNL
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
 * Example program for IDA: Brusselator, parallel, GMRES, IDABBD
 * preconditioner, FSA.
 *
 * This example program for IDAS uses SUNLinSol_SPGMR as the linear solver.
 * It is written for a parallel computer system and uses the
 * IDABBDPRE band-block-diagonal preconditioner module for the
 * IDALS interface.
 *
 * The mathematical problem solved in this example is a DAE system
 * that arises from a system of partial differential equations after
 * spatial discretization.
 * 
 * The PDE system is a two-species time-dependent PDE known as
 * Brusselator PDE and models a chemically reacting system.
 *
 *                  
 *  du/dt = eps1(u  + u  ) + u^2 v -(B+1)u + A  
 *                xx   yy
 *                                              domain [0,L]X[0,L]
 *  dv/dt = eps2(v  + v  ) - u^2 v + Bu
 *                xx   yy
 *
 *  B.C. Neumann
 *  I.C  u(x,y,t0) = u0(x,y) =  1  - 0.5*cos(pi*y/L) 
 *       v(x,y,t0) = v0(x,y) = 3.5 - 2.5*cos(pi*x/L) 
 *
 * The PDEs are discretized by central differencing on a MX by MY
 * mesh, and so the system size Neq is the product MX*MY*NUM_SPECIES. 
 * The system is actually implemented on submeshes, processor by 
 * processor, with an MXSUB by MYSUB mesh on each of NPEX * NPEY 
 * processors.
 *
 * The average of the solution u at final time is also computed.
 *            / /
 *        g = | | u(x,y,tf) dx dy
 *            / /
 * Also the sensitivities of g with respect to parameters eps1 and 
 * eps2 are computed.
 *                  / /
 *       dg/d eps = | | u  (x,y,tf)  dx dy  
 *                  / /  eps
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <idas/idas.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <idas/idas_bbdpre.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <mpi.h>

/* Problem Constants */
#define NUM_SPECIES 2 
#define ctL         RCONST(1.0)    /* Domain =[0,L]^2 */
#define ctA         RCONST(1.0)
#define ctB         RCONST(3.4)
#define ctEps       RCONST(2.0e-3)

#define NS          2

#define PI          RCONST(3.1415926535898) /* pi */ 

#define MXSUB       41    /* Number of x mesh points per processor subgrid */
#define MYSUB       41    /* Number of y mesh points per processor subgrid */
#define NPEX        2     /* Number of subgrids in the x direction */
#define NPEY        2     /* Number of subgrids in the y direction */
#define MX          (MXSUB*NPEX)      /* MX = number of x mesh points */
#define MY          (MYSUB*NPEY)      /* MY = number of y mesh points */
#define NSMXSUB     (NUM_SPECIES * MXSUB)
#define NEQ         (NUM_SPECIES*MX*MY) /* Number of equations in system */


#define RTOL        RCONST(1.e-5)  /*  rtol tolerance */
#define ATOL        RCONST(1.e-5)  /*  atol tolerance */
#define NOUT        6  
#define TMULT       RCONST(10.0)   /* Multiplier for tout values */
#define TADD        RCONST(0.3)    /* Increment for tout values */

#define ZERO        RCONST(0.0)
#define HALF        RCONST(0.5)
#define ONE         RCONST(1.0)


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
  int ns, thispe, npes, ixsub, jysub, npex, npey;
  int mxsub, mysub, nsmxsub, nsmxsub2;
  realtype A, B, L, eps[NUM_SPECIES];
  realtype dx, dy;
  realtype cox[NUM_SPECIES], coy[NUM_SPECIES];
  realtype gridext[(MXSUB+2)*(MYSUB+2)*NUM_SPECIES];
  realtype rhs[NUM_SPECIES];
  MPI_Comm comm;
  realtype rates[2];
  sunindextype n_local;
} *UserData;

/* Prototypes for functions called by the IDA Solver. */
static int res(realtype tt, N_Vector uv, N_Vector uvp,
               N_Vector rr, void *user_data);

static int reslocal(sunindextype Nlocal, realtype tt, 
                    N_Vector uv, N_Vector uvp, N_Vector res, 
                    void *user_data);

static int rescomm(sunindextype Nlocal, realtype tt,
                   N_Vector uv, N_Vector uvp, void *user_data);

/* Integrate over spatial domain. */
static int integr(MPI_Comm comm, N_Vector uv, void *user_data, realtype *intval);

/* Prototypes for supporting functions */
static void BSend(MPI_Comm comm, int thispe, int ixsub, int jysub,
                  int dsizex, int dsizey, realtype carray[]);

static void BRecvPost(MPI_Comm comm, MPI_Request request[], int thispe,
                      int ixsub, int jysub,
                      int dsizex, int dsizey,
                      realtype cext[], realtype buffer[]);

static void BRecvWait(MPI_Request request[], int ixsub, int jysub,
                      int dsizex, realtype cext[], realtype buffer[]);

static void ReactRates(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy, 
                     UserData data);

/* Prototypes for private functions */
static void InitUserData(UserData data, int thispe, int npes, MPI_Comm comm);

static void SetInitialProfiles(N_Vector uv, N_Vector uvp, N_Vector id,
                               N_Vector resid, UserData data);

static void PrintHeader(sunindextype SystemSize, int maxl, 
                        sunindextype mudq, sunindextype mldq,
                        sunindextype mukeep, sunindextype mlkeep,
                        realtype rtol, realtype atol);

static void PrintOutput(void *mem, N_Vector uv, realtype time,
                        UserData data, MPI_Comm comm);

static void PrintSol(void* mem, N_Vector uv, N_Vector uvp, UserData data,
                     MPI_Comm comm);

static void PrintFinalStats(void *mem);

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
  UserData data;
  sunindextype SystemSize, local_N, mudq, mldq, mukeep, mlkeep;
  realtype rtol, atol, t0, tout, tret;
  N_Vector uv, uvp, resid, id, *uvS, *uvpS;
  int thispe, npes, maxl, iout, retval;
  realtype pbar[NS];
  int is;
  realtype intval;

  uv = uvp = resid = id = NULL;
  uvS = uvpS = NULL;
  data = NULL;
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

  /* Set up user data block data. */
  data = (UserData) malloc(sizeof *data);

  InitUserData(data, thispe, npes, comm);
  
  /* Create needed vectors, and load initial values.
     The vector resid is used temporarily only.        */
  
  uv  = N_VNew_Parallel(comm, local_N, SystemSize);
  if(check_retval((void *)uv, "N_VNew_Parallel", 0, thispe)) MPI_Abort(comm, 1);

  uvp  = N_VNew_Parallel(comm, local_N, SystemSize);
  if(check_retval((void *)uvp, "N_VNew_Parallel", 0, thispe)) MPI_Abort(comm, 1);

  resid = N_VNew_Parallel(comm, local_N, SystemSize);
  if(check_retval((void *)resid, "N_VNew_Parallel", 0, thispe)) MPI_Abort(comm, 1);

  id  = N_VNew_Parallel(comm, local_N, SystemSize);
  if(check_retval((void *)id, "N_VNew_Parallel", 0, thispe)) MPI_Abort(comm, 1);

  uvS = N_VCloneVectorArray_Parallel(NS, uv);
  if (check_retval((void *)uvS, "N_VCloneVectorArray_Parallel", 0, thispe)) 
    MPI_Abort(comm, 1);
  for (is=0;is<NS;is++) N_VConst(ZERO, uvS[is]);
    
  uvpS = N_VCloneVectorArray_Parallel(NS, uv);
  if (check_retval((void *)uvpS, "N_VCloneVectorArray_Parallel", 0, thispe))  
    MPI_Abort(comm, 1);
  for (is=0;is<NS;is++) N_VConst(ZERO, uvpS[is]);

  SetInitialProfiles(uv, uvp, id, resid, data);

  /* Set remaining inputs to IDAS. */
  t0 = ZERO;
  rtol = RTOL; 
  atol = ATOL;
  
  /* Call IDACreate and IDAInit to initialize solution */
  ida_mem = IDACreate();
  if(check_retval((void *)ida_mem, "IDACreate", 0, thispe)) MPI_Abort(comm, 1);

  retval = IDASetUserData(ida_mem, data);
  if(check_retval(&retval, "IDASetUserData", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetId(ida_mem, id);
  if(check_retval(&retval, "IDASetId", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDAInit(ida_mem, res, t0, uv, uvp);
  if(check_retval(&retval, "IDAInit", 1, thispe)) MPI_Abort(comm, 1);
  
  retval = IDASStolerances(ida_mem, rtol, atol);
  if(check_retval(&retval, "IDASStolerances", 1, thispe)) MPI_Abort(comm, 1);


  /* Enable forward sensitivity analysis. */
  retval = IDASensInit(ida_mem, NS, IDA_SIMULTANEOUS, NULL, uvS, uvpS);
  if(check_retval(&retval, "IDASensInit", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASensEEtolerances(ida_mem);
  if(check_retval(&retval, "IDASensEEtolerances", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetSensErrCon(ida_mem, SUNTRUE);
  if (check_retval(&retval, "IDASetSensErrCon", 1, thispe)) MPI_Abort(comm, 1);

  pbar[0] = data->eps[0];
  pbar[1] = data->eps[1];
  retval = IDASetSensParams(ida_mem, data->eps, pbar, NULL);
  if (check_retval(&retval, "IDASetSensParams", 1, thispe)) MPI_Abort(comm, 1);

  /* Call SUNLinSol_SPGMR and IDASetLinearSolver to specify the IDAS linear solver */
  maxl = 16;                                      /* max dimension of the Krylov subspace */
  LS = SUNLinSol_SPGMR(uv, PREC_LEFT, maxl);      /* IDA only allows left preconditioning */
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
  PrintOutput(ida_mem, uv, t0, data, comm);


  /* Call IDAS in tout loop, normal mode, and print selected output. */
  for (iout = 1; iout <= NOUT; iout++) {
    
    retval = IDASolve(ida_mem, tout, &tret, uv, uvp, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1, thispe)) MPI_Abort(comm, 1);

    PrintOutput(ida_mem, uv, tret, data, comm);

    if (iout < 3) tout *= TMULT;
    else          tout += TADD;

  }
  /* Print each PE's portion of the solution in a separate file. */
  /* PrintSol(ida_mem, uv, uvp, data, comm); */


  /* On PE 0, print final set of statistics. */  
  if (thispe == 0) {
    PrintFinalStats(ida_mem);
  }

  /* calculate integral of u over domain. */
  integr(comm, uv, data, &intval);
  if (thispe == 0) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("\n\nThe average of u on the domain:\ng = %Lg\n", intval);
#else
    printf("\n\nThe average of u on the domain:\ng = %g\n", intval);
#endif    
  }

  /* integrate the sensitivities of u over domain. */
  IDAGetSens(ida_mem, &tret, uvS);
  if (thispe == 0)
    printf("\nSensitivities of g:\n");

  for (is=0; is<NS; is++) {
    integr(comm, uvS[is], data, &intval);
    if (thispe == 0) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      printf("w.r.t. eps%d = %14.10Lf\n", is, intval);
#else
      printf("w.r.t. eps%d = %14.10f\n", is, intval);
#endif      
    }
  }

  /* Free memory. */
  N_VDestroy_Parallel(uv);
  N_VDestroy_Parallel(uvp);
  N_VDestroy_Parallel(id);
  N_VDestroy_Parallel(resid);
  N_VDestroyVectorArray_Parallel(uvS, NS);
  N_VDestroyVectorArray_Parallel(uvpS, NS);
  IDAFree(&ida_mem);
  SUNLinSolFree(LS);

  free(data);

  MPI_Finalize();

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * InitUserData: Load problem constants in data (of type UserData).   
 */

static void InitUserData(UserData data, int thispe, int npes, 
                         MPI_Comm comm)
{
  data->jysub = thispe / NPEX;
  data->ixsub = thispe - (data->jysub)*NPEX;
  data->mxsub = MXSUB;
  data->mysub = MYSUB;
  data->npex = NPEX;
  data->npey = NPEY;
  data->ns = NUM_SPECIES;
  data->dx = ctL/(MX-1);
  data->dy = ctL/(MY-1);
  data->thispe = thispe;
  data->npes   = npes;
  data->nsmxsub = MXSUB * NUM_SPECIES;
  data->nsmxsub2 = (MXSUB+2)*NUM_SPECIES;
  data->comm = comm;
  data->n_local = MXSUB*MYSUB*NUM_SPECIES;

  data->A = ctA;
  data->B = ctB;
  data->L = ctL;
  data->eps[0] = data->eps[1] = ctEps;
}

/*
 * SetInitialProfiles: Set initial conditions in uv, uvp, and id.
 */

static void SetInitialProfiles(N_Vector uv, N_Vector uvp, N_Vector id,
                               N_Vector resid, UserData data)
{
  int ixsub, jysub, mxsub, mysub, ix, jy;
  realtype *idxy, dx, dy, x, y, *uvxy, *uvxy1, L, npex, npey;

  ixsub = data->ixsub;
  jysub = data->jysub;
  mxsub = data->mxsub;
  mysub = data->mysub;
  npex = data->npex;
  npey = data->npey;
  dx = data->dx;
  dy = data->dy;
  L = data->L;

  /* Loop over grid, load uv values and id values. */
  for (jy = 0; jy < mysub; jy++) {
    y = (jy + jysub*mysub) * dy;
    for (ix = 0; ix < mxsub; ix++) {

      x = (ix + ixsub*mxsub) * dx;
      uvxy = IJ_Vptr(uv,ix,jy); 

      uvxy[0] = RCONST(1.0) - HALF*cos(PI*y/L); 
      uvxy[1] = RCONST(3.5) - RCONST(2.5)*cos(PI*x/L);
    }
  }
 
  N_VConst(ONE, id);

  if (jysub == 0) {
    for (ix=0; ix<mxsub; ix++) {
      idxy = IJ_Vptr(id,ix,0);
      idxy[0] = idxy[1] = ZERO;

      uvxy  = IJ_Vptr(uv,ix,0);
      uvxy1 = IJ_Vptr(uv,ix,1);
      uvxy[0] = uvxy1[0];
      uvxy[1] = uvxy1[1];
    }
  }

  if (ixsub == npex-1) {
    for (jy = 0; jy < mysub; jy++) {
      idxy = IJ_Vptr(id,mxsub-1,jy);
      idxy[0] = idxy[1] = ZERO;

      uvxy  = IJ_Vptr(uv,mxsub-1,jy);
      uvxy1 = IJ_Vptr(uv,mxsub-2,jy);
      uvxy[0] = uvxy1[0];
      uvxy[1] = uvxy1[1];

    }
  }


  if (ixsub == 0) {
    for (jy = 0; jy < mysub; jy++) {
      idxy = IJ_Vptr(id,0,jy);
      idxy[0] = idxy[1] = ZERO;

      uvxy  = IJ_Vptr(uv,0,jy);
      uvxy1 = IJ_Vptr(uv,1,jy);
      uvxy[0] = uvxy1[0];
      uvxy[1] = uvxy1[1];

    }    
  }

  if (jysub == npey-1) {
    for (ix=0; ix<mxsub; ix++) {
      idxy = IJ_Vptr(id,ix,jysub);
      idxy[0] = idxy[1] = ZERO;

      uvxy  = IJ_Vptr(uv,ix,mysub-1);
      uvxy1 = IJ_Vptr(uv,ix,mysub-2);
      uvxy[0] = uvxy1[0];
      uvxy[1] = uvxy1[1];

    }
  }

  /* Derivative found by calling the residual function with uvp = 0. */
  N_VConst(ZERO, uvp);
  res(ZERO, uv, uvp, resid, data);
  N_VScale(-ONE, resid, uvp);


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
  printf("\n Brusselator PDE -  DAE parallel example problem for IDA \n\n");
  printf("Number of species ns: %d", NUM_SPECIES);
  printf("     Mesh dimensions: %d x %d\n", MX, MY);
  printf("Total system size: %ld\n",(long int) SystemSize);
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
  printf("CalcIC called to correct initial concentrations \n\n");
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

static void PrintOutput(void *ida_mem, N_Vector uv, realtype tt,
                        UserData data, MPI_Comm comm)
{
  MPI_Status status;
  realtype *cdata, clast[2], hused;
  long int nst;
  int i, kused, retval, thispe, npelast, ilast;;

  thispe = data->thispe; 
  npelast = data->npes - 1;
  cdata = N_VGetArrayPointer_Parallel(uv);
  
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
 * PrintSol the PE's portion of the solution to a file.
 */
static void PrintSol(void* ida_mem, N_Vector uv, N_Vector uvp, 
                     UserData data, MPI_Comm comm)
{
  FILE* fout;
  realtype *uvxy;
  int ix, jy, mxsub, mysub, thispe;
  char szFilename[128];

  thispe = data->thispe;
  sprintf(szFilename, "ysol%d.txt", thispe);

  fout = fopen(szFilename, "w+");
  if (fout==NULL) {
    printf("PE[% 2d] is unable to write solution to disk!\n", thispe);
    return;
  }

  mxsub = data->mxsub;
  mysub = data->mysub;

  for (jy=0; jy<mysub; jy++) {
    for (ix=0; ix<mxsub; ix++) {
    
      uvxy  = IJ_Vptr(uv, ix, jy);
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(fout, "%Lg\n%Lg\n", uvxy[0], uvxy[1]);
#else
      fprintf(fout, "%g\n%g\n", uvxy[0], uvxy[1]);
#endif      
    }
  }    
  fclose(fout);
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
 * res: System residual function 
 *
 * To compute the residual function F, this routine calls:
 * rescomm, for needed communication, and then
 * reslocal, for computation of the residuals on this processor.      
 */

static int res(realtype tt, N_Vector uv, N_Vector uvp,
               N_Vector rr, void *user_data)
{
  int retval;
  UserData data;
  sunindextype Nlocal;
  
  data = (UserData) user_data;
  
  Nlocal = data->n_local;

  /* Call rescomm to do inter-processor communication. */
  retval = rescomm(Nlocal, tt, uv, uvp, user_data);
  
  /* Call reslocal to calculate the local portion of residual vector. */
  retval = reslocal(Nlocal, tt, uv, uvp, rr, user_data);
  
  return(retval);
}

/*
 * rescomm: Communication routine in support of resweb.
 * This routine performs all inter-processor communication of components
 * of the uv vector needed to calculate F, namely the components at all
 * interior subgrid boundaries (ghost cell data).  It loads this data
 * into a work array cext (the local portion of c, extended).
 * The message-passing uses blocking sends, non-blocking receives,
 * and receive-waiting, in routines BRecvPost, BSend, BRecvWait.         
 */
static int rescomm(sunindextype Nlocal, realtype tt, 
                   N_Vector uv, N_Vector uvp, void *user_data)
{

  UserData data;
  realtype *cdata, *gridext, buffer[2*NUM_SPECIES*MYSUB];
  int thispe, ixsub, jysub, nsmxsub, nsmysub;
  MPI_Comm comm;
  MPI_Request request[4];
  
  data = (UserData) user_data;
  cdata = N_VGetArrayPointer_Parallel(uv);

  /* Get comm, thispe, subgrid indices, data sizes, extended array cext. */  
  comm = data->comm;     
  thispe = data->thispe;

  ixsub = data->ixsub;   
  jysub = data->jysub;
  gridext = data->gridext;
  nsmxsub = data->nsmxsub; 
  nsmysub = (data->ns)*(data->mysub);
  
  /* Start receiving boundary data from neighboring PEs. */
  BRecvPost(comm, request, thispe, ixsub, jysub, nsmxsub, nsmysub, 
            gridext, buffer);
  
  /* Send data from boundary of local grid to neighboring PEs. */
  BSend(comm, thispe, ixsub, jysub, nsmxsub, nsmysub, cdata);
  
  /* Finish receiving boundary data from neighboring PEs. */
  BRecvWait(request, ixsub, jysub, nsmxsub, gridext, buffer);
  
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
                      int ixsub, int jysub, int dsizex, int dsizey,
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
 * This routine sends components of uv from internal subgrid boundaries
 * to the appropriate neighbor PEs.                                      
 */
 
static void BSend(MPI_Comm comm, int my_pe, int ixsub, int jysub,
                  int dsizex, int dsizey, realtype cdata[])
{
  int i;
  int ly, offsetc, offsetbuf;
  realtype bufleft[NUM_SPECIES*MYSUB], bufright[NUM_SPECIES*MYSUB];

  /* If jysub > 0, send data from bottom x-line of uv. */

  if (jysub != 0)
    MPI_Send(&cdata[0], dsizex, PVEC_REAL_MPI_TYPE, my_pe-NPEX, 0, comm);

  /* If jysub < NPEY-1, send data from top x-line of uv. */

  if (jysub != NPEY-1) {
    offsetc = (MYSUB-1)*dsizex;
    MPI_Send(&cdata[offsetc], dsizex, PVEC_REAL_MPI_TYPE, my_pe+NPEX, 0, comm);
  }

  /* If ixsub > 0, send data from left y-line of uv (via bufleft). */

  if (ixsub != 0) {
    for (ly = 0; ly < MYSUB; ly++) {
      offsetbuf = ly*NUM_SPECIES;
      offsetc = ly*dsizex;
      for (i = 0; i < NUM_SPECIES; i++)
        bufleft[offsetbuf+i] = cdata[offsetc+i];
    }
    MPI_Send(&bufleft[0], dsizey, PVEC_REAL_MPI_TYPE, my_pe-1, 0, comm);   
  }
  
  /* If ixsub < NPEX-1, send data from right y-line of uv (via bufright). */
  
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

#define mxsub      (data->mxsub)
#define mysub      (data->mysub)
#define npex       (data->npex)
#define npey       (data->npey)
#define ixsub      (data->ixsub)
#define jysub      (data->jysub)
#define nsmxsub    (data->nsmxsub)
#define nsmxsub2   (data->nsmxsub2)
#define dx         (data->dx)
#define dy         (data->dy)
#define cox        (data->cox)
#define coy        (data->coy)
#define gridext    (data->gridext)
#define eps        (data->eps)
#define ns         (data->ns)

/*
 * reslocal: Compute res = F(t,uv,uvp).
 * This routine assumes that all inter-processor communication of data
 * needed to calculate F has already been done.  Components at interior
 * subgrid boundaries are assumed to be in the work array cext.
 * The local portion of the uv vector is first copied into cext.
 * The exterior Neumann boundary conditions are explicitly handled here
 * by copying data from the first interior mesh line to the ghost cell
 * locations in cext.  Then the reaction and diffusion terms are
 * evaluated in terms of the cext array, and the residuals are formed.
 * The reaction terms are saved separately in the vector data->rates
 * for use by the preconditioner setup routine.                          
 */

static int reslocal(sunindextype Nlocal, realtype tt, N_Vector uv, 
                    N_Vector uvp, N_Vector rr, void *user_data)
{
  realtype *uvdata, *uvpxy, *resxy, xx, yy, dcyli, dcyui, dcxli, dcxui, dx2, dy2;
  realtype ixend, ixstart, jystart, jyend;
  int ix, jy, is, i, locc, ylocce, locce;
  realtype rates[2];
  UserData data;
  
  data = (UserData) user_data;

  /* Get data pointers, subgrid data, array sizes, work array cext. */
  uvdata = N_VGetArrayPointer_Parallel(uv);

  dx2 = dx * dx;
  dy2 = dy * dy;

  /* Copy local segment of uv vector into the working extended array gridext. */
  locc = 0;
  locce = nsmxsub2 + NUM_SPECIES;
  for (jy = 0; jy < mysub; jy++) {
    for (i = 0; i < nsmxsub; i++) gridext[locce+i] = uvdata[locc+i];
    locc = locc + nsmxsub;
    locce = locce + nsmxsub2;
  }

  /* To facilitate homogeneous Neumann boundary conditions, when this is
     a boundary PE, copy data from the first interior mesh line of uv to gridext. */
  
  /* If jysub = 0, copy x-line 2 of uv to gridext. */
  if (jysub == 0)
    { for (i = 0; i < nsmxsub; i++) gridext[NUM_SPECIES+i] = uvdata[nsmxsub+i]; }
  

  /* If jysub = npey-1, copy x-line mysub-1 of uv to gridext. */
  if (jysub == npey-1) {
    locc = (mysub-2)*nsmxsub;
    locce = (mysub+1)*nsmxsub2 + NUM_SPECIES;
    for (i = 0; i < nsmxsub; i++) gridext[locce+i] = uvdata[locc+i];
  }
  

  /* If ixsub = 0, copy y-line 2 of uv to gridext. */
  if (ixsub == 0) {
    for (jy = 0; jy < mysub; jy++) {
      locc = jy*nsmxsub + NUM_SPECIES;
      locce = (jy+1)*nsmxsub2;
      for (i = 0; i < NUM_SPECIES; i++) gridext[locce+i] = uvdata[locc+i];
    }
  }
  

  /* If ixsub = npex-1, copy y-line mxsub-1 of uv to gridext. */
  if (ixsub == npex-1) {
    for (jy = 0; jy < mysub; jy++) {
      locc  = (jy+1)*nsmxsub - 2*NUM_SPECIES;
      locce = (jy+2)*nsmxsub2 - NUM_SPECIES;
      for (i = 0; i < NUM_SPECIES; i++) gridext[locce+i] = uvdata[locc+i];
    }
  }
  
  /* Loop over all grid points, setting local array rates to right-hand sides.
     Then set rr values appropriately (ODE in the interior and DAE on the boundary)*/
  ixend = ixstart = jystart = jyend = 0;  

  if (jysub==0)      jystart = 1;
  if (jysub==npey-1) jyend   = 1;
  if (ixsub==0)      ixstart = 1;
  if (ixsub==npex-1) ixend   = 1;

  for (jy = jystart; jy < mysub-jyend; jy++) { 
    ylocce = (jy+1)*nsmxsub2;
    yy     = (jy+jysub*mysub)*dy;

    for (ix = ixstart; ix < mxsub-ixend; ix++) {
      locce = ylocce + (ix+1)*NUM_SPECIES;
      xx = (ix + ixsub*mxsub)*dx;

      ReactRates(xx, yy, &(gridext[locce]), rates, data);      
      
      resxy = IJ_Vptr(rr,ix,jy); 
      uvpxy = IJ_Vptr(uvp,ix,jy); 
      
      for (is = 0; is < NUM_SPECIES; is++) {
        dcyli = gridext[locce+is]          - gridext[locce+is-nsmxsub2];
        dcyui = gridext[locce+is+nsmxsub2] - gridext[locce+is];
        
        dcxli = gridext[locce+is]             - gridext[locce+is-NUM_SPECIES];
        dcxui = gridext[locce+is+NUM_SPECIES] - gridext[locce+is];
        
        resxy[is] = uvpxy[is]-eps[is]*((dcxui-dcxli)/dx2+(dcyui-dcyli)/dy2)-rates[is];
      }
    }
  }

  /* Algebraic equation correspoding to boundary mesh point. */
  if (jysub==0) {
    for (ix=0; ix<mxsub; ix++) {
      
      locce = nsmxsub2 + NUM_SPECIES * (ix+1);
      resxy = IJ_Vptr(rr,ix,0); 

      for (is=0; is<NUM_SPECIES; is++)
        resxy[is] = gridext[locce+is+nsmxsub2] - gridext[locce+is];
    }
  }

  if (ixsub==npex-1) {
    for(jy=0; jy<mysub; jy++) {
      locce = (jy+1)*nsmxsub2 + nsmxsub2-NUM_SPECIES;
      resxy = IJ_Vptr(rr,mxsub-1,jy); 

      for (is=0; is<NUM_SPECIES; is++)
        resxy[is] = gridext[locce+is-NUM_SPECIES] - gridext[locce+is];
    }
  }
  
  if (ixsub==0) {
    for (jy=0; jy<mysub; jy++) {
      locce = (jy+1)*nsmxsub2 + NUM_SPECIES;
      resxy = IJ_Vptr(rr,0,jy); 

      for (is=0; is<NUM_SPECIES; is++)
        resxy[is] = gridext[locce+is-NUM_SPECIES] - gridext[locce+is];
    }    
  }

  if (jysub==npey-1) {
    for(ix=0; ix<mxsub; ix++) {
      locce = nsmxsub2*mysub + (ix+1)*NUM_SPECIES;
      resxy = IJ_Vptr(rr,ix, mysub-1);

      for (is=0; is<NUM_SPECIES; is++)
        resxy[is] = gridext[locce+is-nsmxsub2] - gridext[locce+is];
    }
  }
  return(0);
}

/* 
 * ReactRates: Evaluate reaction rates at a given spatial point.   
 * At a given (x,y), evaluate the array of ns reaction terms R.          
 */

static void ReactRates(realtype xx, realtype yy, realtype *uvval, 
                       realtype *rates, UserData data)
{
  realtype A, B;

  A = data->A; B = data->B;

  rates[0] = uvval[0]*uvval[0]*uvval[1];
  rates[1] = - rates[0];

  rates[0] += A-(B+1)*uvval[0];
  rates[1] += B*uvval[0];
}

/* Integrate over the spatial domain. Each process computes the integral on its
   grid. Then processes call MPI_REDUCE to compute sum of the local values. */
static int integr(MPI_Comm comm, N_Vector uv, void *user_data, realtype *intval)
{
  int ix, jy;
  int retval;
  realtype *uvdata;
  UserData data;
  
  realtype buf[2];
  
  data = (UserData) user_data;

  /* compute the integral on the (local) grid */
  uvdata = N_VGetArrayPointer_Parallel(uv);

  *intval = 0;

  for (jy=1; jy<mysub; jy++)
    for (ix=1; ix<mxsub; ix++) 
       /* consider only u */
      *intval += uvdata[ix*NUM_SPECIES + jy*NSMXSUB];
  *intval *= (dx*dy);

  buf[0] = *intval;

  /* Sum local values and get the result on all processors. */
  retval = MPI_Allreduce(buf, buf+1, 1, MPI_DOUBLE, MPI_SUM, comm);

  *intval = buf[1];
  return(retval);
}
