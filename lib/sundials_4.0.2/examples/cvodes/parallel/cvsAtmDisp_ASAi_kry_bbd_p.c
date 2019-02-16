/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Lukas Jager and Radu Serban @ LLNL
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
 * Parallel Krylov adjoint sensitivity example problem.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_bbdpre.h>
#include <sunlinsol/sunlinsol_spgmr.h> 
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <mpi.h>

/*
 *------------------------------------------------------------------
 * Constants
 *------------------------------------------------------------------
 */

#ifdef USE3D
#define DIM 3
#else
#define DIM 2
#endif

/* Domain definition */

#define XMIN RCONST(0.0)
#define XMAX RCONST(20.0)
#define MX   80    /* no. of divisions in x dir. */
#define NPX  2     /* no. of procs. in x dir.    */

#define YMIN RCONST(0.0)
#define YMAX RCONST(20.0)
#define MY   80    /* no. of divisions in y dir. */
#define NPY  4     /* no. of procs. in y dir.    */

#ifdef USE3D
#define ZMIN RCONST(0.0)
#define ZMAX RCONST(20.0)
#define MZ   40    /* no. of divisions in z dir. */
#define NPZ  2     /* no. of procs. in z dir.    */
#endif

/* Parameters for source Gaussians */

#define G1_AMPL   RCONST(1.0)
#define G1_SIGMA  RCONST(1.7) 
#define G1_X      RCONST(4.0)
#define G1_Y      RCONST(8.0)
#ifdef USE3D
#define G1_Z      RCONST(8.0)
#endif

#define G2_AMPL   RCONST(0.8)
#define G2_SIGMA  RCONST(3.0)
#define G2_X      RCONST(16.0)
#define G2_Y      RCONST(12.0)
#ifdef USE3D
#define G2_Z      RCONST(12.0)
#endif

#define G_MIN     RCONST(1.0e-5)

/* Diffusion coeff., max. velocity, domain width in y dir. */

#define DIFF_COEF RCONST(1.0)
#define V_MAX     RCONST(1.0)
#define L         ((YMAX-YMIN)/RCONST(2.0))
#define V_COEFF   V_MAX/L/L

/* Initial and final times */

#define ti    RCONST(0.0)
#define tf    RCONST(10.0)

/* Integration tolerances */

#define RTOL    RCONST(1.0e-8) /* states */
#define ATOL    RCONST(1.0e-6)

#define RTOL_Q  RCONST(1.0e-8) /* forward quadrature */
#define ATOL_Q  RCONST(1.0e-6)

#define RTOL_B  RCONST(1.0e-8) /* adjoint variables */
#define ATOL_B  RCONST(1.0e-6)

#define RTOL_QB RCONST(1.0e-8) /* backward quadratures */
#define ATOL_QB RCONST(1.0e-6)

/* Steps between check points */

#define STEPS 200

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/*
 *------------------------------------------------------------------
 * Macros
 *------------------------------------------------------------------
 */

#define FOR_DIM for(dim=0; dim<DIM; dim++)

/* IJth:     (i[0],i[1],i[2])-th vector component                       */
/* IJth_ext: (i[0],i[1],i[2])-th vector component in the extended array */

#ifdef USE3D
#define IJth(y,i)     ( y[(i[0])+(l_m[0]*((i[1])+(i[2])*l_m[1]))] )
#define IJth_ext(y,i) ( y[(i[0]+1)+((l_m[0]+2)*((i[1]+1)+(i[2]+1)*(l_m[1]+2)))] )
#else
#define IJth(y,i)     (y[i[0]+(i[1])*l_m[0]])
#define IJth_ext(y,i) (y[ (i[0]+1) + (i[1]+1) * (l_m[0]+2)])
#endif

/*
 *------------------------------------------------------------------
 * Type definition: ProblemData 
 *------------------------------------------------------------------
 */

typedef struct {
  /* Domain */
  realtype xmin[DIM];  /* "left" boundaries */  
  realtype xmax[DIM];  /* "right" boundaries */
  int m[DIM];          /* number of grid points */
  realtype dx[DIM];    /* grid spacing */
  realtype dOmega;     /* differential volume */

  /* Parallel stuff */
  MPI_Comm comm;       /* MPI communicator */
  int myId;            /* process id */ 
  int npes;            /* total number of processes */
  int num_procs[DIM];  /* number of processes in each direction */
  int nbr_left[DIM];   /* MPI ID of "left" neighbor */
  int nbr_right[DIM];  /* MPI ID of "right" neighbor */
  int m_start[DIM];    /* "left" index in the global domain */
  int l_m[DIM];        /* number of local grid points */ 
  realtype *y_ext;     /* extended data array */
  realtype *buf_send;  /* Send buffer */
  realtype *buf_recv;  /* Receive buffer */
  int buf_size;        /* Buffer size */

  /* Source */
  N_Vector p;          /* Source parameters */ 

} *ProblemData;

/*
 *------------------------------------------------------------------
 * Interface functions to CVODES
 *------------------------------------------------------------------
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f_local(sunindextype Nlocal, realtype t, N_Vector y, 
                   N_Vector ydot, void *user_data);

static int fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data);


static int fB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, 
              void *user_dataB);
static int fB_local(sunindextype NlocalB, realtype t, 
                    N_Vector y, N_Vector yB, N_Vector yBdot, 
                    void *user_dataB);

static int fQB(realtype t, N_Vector y, N_Vector yB, 
               N_Vector qBdot, void *user_dataB);

/*
 *------------------------------------------------------------------
 * Private functions
 *------------------------------------------------------------------
 */

static void SetData(ProblemData d, MPI_Comm comm, int npes, int myId,
                    sunindextype *neq, sunindextype *l_neq);
static void SetSource(ProblemData d);
static void f_comm(sunindextype Nlocal, realtype t, N_Vector y, void *user_data);
static void Load_yext(realtype *src, ProblemData d);
static void PrintHeader();
static int PrintFinalStats(void *cvode_mem);
static void OutputGradient(int myId, N_Vector qB, ProblemData d);
static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

/*
 *------------------------------------------------------------------
 * Main program
 *------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  ProblemData d;

  MPI_Comm comm;
  int npes, npes_needed;
  int myId;
 
  sunindextype neq, l_neq;

  void *cvode_mem;
  SUNLinearSolver LS;
  N_Vector y, q;
  realtype abstol, reltol, abstolQ, reltolQ;
  sunindextype mudq, mldq, mukeep, mlkeep;

  int indexB;
  N_Vector yB, qB;
  realtype abstolB, reltolB, abstolQB, reltolQB;
  sunindextype mudqB, mldqB, mukeepB, mlkeepB;

  realtype tret, *qdata, G;

  int ncheckpnt, retval;

  booleantype output;

  /* Initialize MPI and set Ids */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myId);

  /* Check number of processes */
  npes_needed = NPX * NPY;
#ifdef USE3D
  npes_needed *= NPZ;
#endif
  MPI_Comm_size(comm, &npes);
  if (npes_needed != npes) {
    if (myId == 0)
      fprintf(stderr,"I need %d processes but I only got %d\n",
              npes_needed, npes);
    MPI_Abort(comm, EXIT_FAILURE);
  }

  /* Test if matlab output is requested */
  if (argc > 1) output = SUNTRUE;
  else          output = SUNFALSE;

  /* Allocate and set problem data structure */
  d = (ProblemData) malloc(sizeof *d);
  SetData(d, comm, npes, myId, &neq, &l_neq);
  
  if (myId == 0) PrintHeader();

  /*-------------------------- 
    Forward integration phase
    --------------------------*/

  /* Allocate space for y and set it with the I.C. */
  y = N_VNew_Parallel(comm, l_neq, neq);
  if(check_retval(y, "N_VNew_Parallel", 0, myId)) MPI_Abort(comm, 1);
  N_VConst(ZERO, y);
  
  /* Allocate and initialize qB (local contribution to cost) */
  q = N_VNew_Parallel(comm, 1, npes); 
  if(check_retval(q, "N_VNew_Parallel", 0, myId)) MPI_Abort(comm, 1);
  N_VConst(ZERO, q);

  /* Create CVODES object, attach user data, and allocate space */
  cvode_mem = CVodeCreate(CV_BDF);
  if(check_retval(cvode_mem, "CVodeCreate", 0, myId)) MPI_Abort(comm, 1);

  retval = CVodeSetUserData(cvode_mem, d);
  if(check_retval(&retval, "CVodeSetUserData", 1, myId)) MPI_Abort(comm, 1);
  
  retval = CVodeInit(cvode_mem, f, ti, y);
  if(check_retval(&retval, "CVodeInit", 1, myId)) MPI_Abort(comm, 1);
  
  abstol = ATOL;  
  reltol = RTOL;   
  
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_retval(&retval, "CVodeSStolerances", 1, myId)) MPI_Abort(comm, 1);

  /* create and attach linear solver */
  LS = SUNLinSol_SPGMR(y, PREC_LEFT, 0);
  if(check_retval(LS, "SUNLinSol_SPGMR", 0, myId)) MPI_Abort(comm, 1);

  retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1, myId)) MPI_Abort(comm, 1);
  
  /* Attach preconditioner and linear solver modules */
  mudq = mldq = d->l_m[0]+1;
  mukeep = mlkeep = 2;  
  retval = CVBBDPrecInit(cvode_mem, l_neq, mudq, mldq, 
                       mukeep, mlkeep, ZERO,
                       f_local, NULL);
  if(check_retval(&retval, "CVBBDPrecInit", 1, myId)) MPI_Abort(comm, 1);
  
  /* Initialize quadrature calculations */
  abstolQ = ATOL_Q;
  reltolQ = RTOL_Q;
  
  retval = CVodeQuadInit(cvode_mem, fQ, q);
  if(check_retval(&retval, "CVodeQuadInit", 1, myId)) MPI_Abort(comm, 1);
  
  retval = CVodeQuadSStolerances(cvode_mem, reltolQ, abstolQ);
  if(check_retval(&retval, "CVodeQuadSStolerances", 1, myId)) MPI_Abort(comm, 1);
  
  retval = CVodeSetQuadErrCon(cvode_mem, SUNTRUE);
  if(check_retval(&retval, "CVodesSetQuadErrCon", 1, myId)) MPI_Abort(comm, 1);

  /* Allocate space for the adjoint calculation */
  retval = CVodeAdjInit(cvode_mem, STEPS, CV_HERMITE);
  if(check_retval(&retval, "CVodeAdjInit", 1, myId)) MPI_Abort(comm, 1);

  /* Integrate forward in time while storing check points */
  if (myId == 0) printf("Begin forward integration... ");
  retval = CVodeF(cvode_mem, tf, y, &tret, CV_NORMAL, &ncheckpnt);
  if(check_retval(&retval, "CVodeF", 1, myId)) MPI_Abort(comm, 1);
  if (myId == 0) printf("done. ");

   /* Extract quadratures */
  retval = CVodeGetQuad(cvode_mem, &tret, q);
  if(check_retval(&retval, "CVodeGetQuad", 1, myId)) MPI_Abort(comm, 1);
  
  qdata = N_VGetArrayPointer_Parallel(q);
  MPI_Allreduce(&qdata[0], &G, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  if (myId == 0) printf("  G = %Le\n",G);
#else
  if (myId == 0) printf("  G = %e\n",G);
#endif

  /* Print statistics for forward run */
  if (myId == 0) {
    retval = PrintFinalStats(cvode_mem);
    if(check_retval(&retval, "PrintFinalStats", 1, myId)) MPI_Abort(comm, 1);
  } 

  /*-------------------------- 
    Backward integration phase
    --------------------------*/
 
  /* Allocate and initialize yB */
  yB = N_VNew_Parallel(comm, l_neq, neq); 
  if(check_retval(yB, "N_VNew_Parallel", 0, myId)) MPI_Abort(comm, 1);
  N_VConst(ZERO, yB);

  /* Allocate and initialize qB (gradient) */
  qB = N_VNew_Parallel(comm, l_neq, neq); 
  if(check_retval(qB, "N_VNew_Parallel", 0, myId)) MPI_Abort(comm, 1);
  N_VConst(ZERO, qB);

  /* Create and allocate backward CVODE memory */
  retval = CVodeCreateB(cvode_mem, CV_BDF, &indexB);
  if(check_retval(&retval, "CVodeCreateB", 1, myId)) MPI_Abort(comm, 1);
  
  retval = CVodeSetUserDataB(cvode_mem, indexB, d);
  if(check_retval(&retval, "CVodeSetUserDataB", 1, myId)) MPI_Abort(comm, 1);
  
  retval = CVodeInitB(cvode_mem, indexB, fB, tf, yB);
  if(check_retval(&retval, "CVodeInitB", 1, myId)) MPI_Abort(comm, 1);
  
  abstolB = ATOL_B;  
  reltolB = RTOL_B;

  retval = CVodeSStolerancesB(cvode_mem, indexB, reltolB, abstolB);
  if(check_retval(&retval, "CVodeSStolerancesB", 1, myId)) MPI_Abort(comm, 1);

  /* Attach preconditioner and linear solver modules */
  retval = CVodeSetLinearSolverB(cvode_mem, indexB, LS, NULL);
  if(check_retval(&retval, "CVodeSetLinearSolverB", 1, myId)) MPI_Abort(comm, 1);
  
  mudqB = mldqB = d->l_m[0]+1;
  mukeepB = mlkeepB = 2;  
  
  retval = CVBBDPrecInitB(cvode_mem, indexB, l_neq, mudqB, mldqB, 
                        mukeepB, mlkeepB, ZERO, fB_local, NULL);
  if(check_retval(&retval, "CVBBDPrecInitB", 1, myId)) MPI_Abort(comm, 1);

  /* Initialize quadrature calculations */
  abstolQB = ATOL_QB;
  reltolQB = RTOL_QB;
  
  retval = CVodeQuadInitB(cvode_mem, indexB, fQB, qB);
  if(check_retval(&retval, "CVodeQuadInitB", 1, myId)) MPI_Abort(comm, 1);
  
  retval = CVodeQuadSStolerancesB(cvode_mem, indexB, reltolQB, abstolQB);
  if(check_retval(&retval, "CVodeQuadSStolerancesB", 1, myId)) MPI_Abort(comm, 1);
  
  retval = CVodeSetQuadErrConB(cvode_mem, indexB, SUNTRUE);
  if(check_retval(&retval, "CVodeSetQuadErrConB", 1, myId)) MPI_Abort(comm, 1);

  /* Integrate backwards */
  if (myId == 0) printf("Begin backward integration... ");
  retval = CVodeB(cvode_mem, ti, CV_NORMAL);
  if(check_retval(&retval, "CVodeB", 1, myId)) MPI_Abort(comm, 1);
  if (myId == 0) printf("done.\n");
  
  /* Extract solution */
  retval = CVodeGetB(cvode_mem, indexB, &tret, yB);
  if(check_retval(&retval, "CVodeGetB", 1, myId)) MPI_Abort(comm, 1);

  /* Extract quadratures */
  retval = CVodeGetQuadB(cvode_mem, indexB, &tret, qB);
  if(check_retval(&retval, "CVodeGetQuadB", 1, myId)) MPI_Abort(comm, 1);

  /* Print statistics for backward run */
  if (myId == 0) {
    retval = PrintFinalStats(CVodeGetAdjCVodeBmem(cvode_mem, indexB));
    if(check_retval(&retval, "PrintFinalStats", 1, myId)) MPI_Abort(comm, 1);
  }

  /* Process 0 collects the gradient components and prints them */
  if (output) {
    OutputGradient(myId, qB, d);
    if (myId == 0) printf("Wrote matlab file 'grad.m'.\n");
  }

  /* Free memory */

  N_VDestroy_Parallel(y);
  N_VDestroy_Parallel(q);
  N_VDestroy_Parallel(qB);
  N_VDestroy_Parallel(yB);

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);

  MPI_Finalize();

  return(0);
}

/*
 *------------------------------------------------------------------
 * SetData:
 * Allocate space for the ProblemData structure.
 * Set fields in the ProblemData structure.
 * Return local and global problem dimensions.
 *
 * SetSource:
 * Instantiates the source parameters for a combination of two
 * Gaussian sources.
 *------------------------------------------------------------------
 */

static void SetData(ProblemData d, MPI_Comm comm, int npes, int myId,
                    sunindextype *neq, sunindextype *l_neq)
{
  int n[DIM], nd[DIM];
  int dim, size;

  /* Set MPI communicator, id, and total number of processes */

  d->comm = comm;
  d->myId = myId;
  d->npes = npes;

  /* Set domain boundaries */

  d->xmin[0] = XMIN;
  d->xmax[0] = XMAX;
  d->m[0]    = MX;

  d->xmin[1] = YMIN;
  d->xmax[1] = YMAX;
  d->m[1]    = MY;

#ifdef USE3D
  d->xmin[2] = ZMIN;
  d->xmax[2] = ZMAX;
  d->m[2]    = MZ;
#endif

  /* Calculate grid spacing and differential volume */

  d->dOmega = ONE;
  FOR_DIM {
    d->dx[dim] = (d->xmax[dim] - d->xmin[dim]) / d->m[dim];
    d->m[dim] +=1;
    d->dOmega *= d->dx[dim];
  }

  /* Set partitioning */

  d->num_procs[0] = NPX;
  n[0] = NPX; 
  nd[0] = d->m[0] / NPX;

  d->num_procs[1] = NPY;
  n[1] = NPY; 
  nd[1] = d->m[1] / NPY;

#ifdef USE3D
  d->num_procs[2] = NPZ;
  n[2] = NPZ; 
  nd[2] = d->m[2] / NPZ;
#endif
  
  /* Compute the neighbors */

  d->nbr_left[0]  = (myId%n[0]) == 0                ? myId : myId-1;
  d->nbr_right[0] = (myId%n[0]) == n[0]-1           ? myId : myId+1;

  d->nbr_left[1]  = (myId/n[0])%n[1] == 0           ? myId : myId-n[0];
  d->nbr_right[1] = (myId/n[0])%n[1] == n[1]-1      ? myId : myId+n[0];

#ifdef USE3D
  d->nbr_left[2]  = (myId/n[0]/n[1])%n[2] == 0      ? myId : myId-n[0]*n[1];
  d->nbr_right[2] = (myId/n[0]/n[1])%n[2] == n[2]-1 ? myId : myId+n[0]*n[1];
#endif
 
  /* Compute the local subdomains 
     m_start: left border in global index space 
     l_m:     length of the subdomain */
  
  d->m_start[0] = (myId%n[0])*nd[0];
  d->l_m[0]     = d->nbr_right[0] == myId ? d->m[0] - d->m_start[0] : nd[0];

  d->m_start[1] = ((myId/n[0])%n[1])*nd[1];
  d->l_m[1]     = d->nbr_right[1] == myId ? d->m[1] - d->m_start[1] : nd[1];

#ifdef USE3D
  d->m_start[2] = (myId/n[0]/n[1])*nd[2];
  d->l_m[2]     = d->nbr_right[2] == myId ? d->m[2] - d->m_start[2] : nd[2];
#endif

  /* Allocate memory for the y_ext array 
     (local solution + data from neighbors) */

  size = 1;
  FOR_DIM size *= d->l_m[dim]+2;
  d->y_ext = (realtype *) malloc( size*sizeof(realtype));

  /* Initialize Buffer field.
     Size of buffer is checked when needed */

  d->buf_send = NULL;
  d->buf_recv = NULL;
  d->buf_size = 0;   

  /* Allocate space for the source parameters */

  *neq = 1; *l_neq = 1;
  FOR_DIM {*neq *= d->m[dim];  *l_neq *= d->l_m[dim];}
  d->p = N_VNew_Parallel(comm, *l_neq, *neq);

  /* Initialize the parameters for a source with Gaussian profile */

  SetSource(d);

}

static void SetSource(ProblemData d)
{
  int *l_m, *m_start;
  realtype *xmin, *dx;
  realtype x[DIM], g, *pdata;
  int i[DIM];

  l_m  = d->l_m;
  m_start = d->m_start;
  xmin = d->xmin;
  dx = d->dx;


  pdata = N_VGetArrayPointer_Parallel(d->p);

  for(i[0]=0; i[0]<l_m[0]; i[0]++) {
    x[0] = xmin[0] + (m_start[0]+i[0]) * dx[0];
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {
      x[1] = xmin[1] + (m_start[1]+i[1]) * dx[1];
#ifdef USE3D
      for(i[2]=0; i[2]<l_m[2]; i[2]++) {
        x[2] = xmin[2] + (m_start[2]+i[2]) * dx[2];
        
        g = G1_AMPL 
          * SUNRexp( -SUNSQR(G1_X-x[0])/SUNSQR(G1_SIGMA) )
          * SUNRexp( -SUNSQR(G1_Y-x[1])/SUNSQR(G1_SIGMA) )
          * SUNRexp( -SUNSQR(G1_Z-x[2])/SUNSQR(G1_SIGMA) );
        
        g += G2_AMPL 
          * SUNRexp( -SUNSQR(G2_X-x[0])/SUNSQR(G2_SIGMA) )
          * SUNRexp( -SUNSQR(G2_Y-x[1])/SUNSQR(G2_SIGMA) )
          * SUNRexp( -SUNSQR(G2_Z-x[2])/SUNSQR(G2_SIGMA) );
        
        if( g < G_MIN ) g = ZERO;

        IJth(pdata, i) = g;
      }
#else
      g = G1_AMPL 
        * SUNRexp( -SUNSQR(G1_X-x[0])/SUNSQR(G1_SIGMA) )
        * SUNRexp( -SUNSQR(G1_Y-x[1])/SUNSQR(G1_SIGMA) );

      g += G2_AMPL 
        * SUNRexp( -SUNSQR(G2_X-x[0])/SUNSQR(G2_SIGMA) )
        * SUNRexp( -SUNSQR(G2_Y-x[1])/SUNSQR(G2_SIGMA) );
      
      if( g < G_MIN ) g = ZERO;

      IJth(pdata, i) = g;
#endif 
    }
  }
}

/*
 *------------------------------------------------------------------
 * f_comm: 
 * Function for inter-process communication
 * Used both for the forward and backward phase.
 *------------------------------------------------------------------
 */

static void f_comm(sunindextype N_local, realtype t, N_Vector y, void *user_data)
{
  int id, n[DIM], proc_cond[DIM], nbr[DIM][2];
  ProblemData d;
  realtype *yextdata, *ydata;
  int l_m[DIM], dim;
  int c, i[DIM], l[DIM-1];
  realtype *buf_send, *buf_recv;
  MPI_Status stat;
  MPI_Comm comm;
  int dir, size = 1, small = INT_MAX;

  d  = (ProblemData) user_data;
  comm = d->comm;
  id = d->myId;
  
  /* extract data from domain*/
  FOR_DIM {
    n[dim] = d->num_procs[dim];
    l_m[dim] = d->l_m[dim];
  }
  yextdata = d->y_ext;
  ydata    = N_VGetArrayPointer_Parallel(y);
  
  /* Calculate required buffer size */
  FOR_DIM {
    size *= l_m[dim];
    if( l_m[dim] < small) small = l_m[dim];
  }
  size /= small;
  
  /* Adjust buffer size if necessary */
  if( d->buf_size < size ) {
    d->buf_send = (realtype*) realloc( d->buf_send, size * sizeof(realtype));
    d->buf_recv = (realtype*) realloc( d->buf_recv, size * sizeof(realtype));
    d->buf_size = size;
  }

  buf_send = d->buf_send;
  buf_recv = d->buf_recv;
  
  /* Compute the communication pattern; who sends first? */
  /* if proc_cond==1 , process sends first in this dimension */
  proc_cond[0] = (id%n[0])%2;
  proc_cond[1] = ((id/n[0])%n[1])%2;
#ifdef USE3D
  proc_cond[2] = (id/n[0]/n[1])%2;
#endif

  /* Compute the actual communication pattern */
  /* nbr[dim][0] is first proc to communicate with in dimension dim */
  /* nbr[dim][1] the second one */
  FOR_DIM {
    nbr[dim][proc_cond[dim]]  = d->nbr_left[dim];
    nbr[dim][!proc_cond[dim]] = d->nbr_right[dim];
  }
  
  /* Communication: loop over dimension and direction (left/right) */
  FOR_DIM {

    for (dir=0; dir<=1; dir++) {

      /* If subdomain at boundary, no communication in this direction */

      if (id != nbr[dim][dir]) {
        c=0;
        /* Compute the index of the boundary (right or left) */
        i[dim] = (dir ^ proc_cond[dim]) ? (l_m[dim]-1) : 0;
        /* Loop over all other dimensions and copy data into buf_send */
        l[0]=(dim+1)%DIM;
#ifdef USE3D
        l[1]=(dim+2)%DIM;
        for(i[l[1]]=0; i[l[1]]<l_m[l[1]]; i[l[1]]++) 
#endif
          for(i[l[0]]=0; i[l[0]]<l_m[l[0]]; i[l[0]]++) 
            buf_send[c++] = IJth(ydata, i);
	  
        if ( proc_cond[dim] ) {
          /* Send buf_send and receive into buf_recv */
          MPI_Send(buf_send, c, PVEC_REAL_MPI_TYPE, nbr[dim][dir], 0, comm);
          MPI_Recv(buf_recv, c, PVEC_REAL_MPI_TYPE, nbr[dim][dir], 0, comm, &stat);
        } else {
          /* Receive into buf_recv and send buf_send*/
          MPI_Recv(buf_recv, c, PVEC_REAL_MPI_TYPE, nbr[dim][dir], 0, comm, &stat);
          MPI_Send(buf_send, c, PVEC_REAL_MPI_TYPE, nbr[dim][dir], 0, comm);
        }

        c=0;

        /* Compute the index of the boundary (right or left) in yextdata */
        i[dim] = (dir ^ proc_cond[dim]) ? l_m[dim] : -1;

        /* Loop over all other dimensions and copy data into yextdata */
#ifdef USE3D
        for(i[l[1]]=0; i[l[1]]<l_m[l[1]]; i[l[1]]++)
#endif
          for(i[l[0]]=0; i[l[0]]<l_m[l[0]]; i[l[0]]++)
            IJth_ext(yextdata, i) = buf_recv[c++];
      }
    } /* end loop over direction */
  } /* end loop over dimension */ 
}

/*
 *------------------------------------------------------------------
 * f and f_local:
 * Forward phase ODE right-hand side
 *------------------------------------------------------------------
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  ProblemData d;
  sunindextype l_neq=1;
  int dim;

  d = (ProblemData) user_data;
  FOR_DIM l_neq *= d->l_m[dim];
  
  /* Do all inter-processor communication */
  f_comm(l_neq, t, y, user_data);

  /* Compute right-hand side locally */
  f_local(l_neq, t, y, ydot, user_data);

  return(0);
}

static int f_local(sunindextype Nlocal, realtype t, N_Vector y, 
                   N_Vector ydot, void *user_data)
{
  realtype *Ydata, *dydata, *pdata;
  realtype dx[DIM], c, v[DIM], cl[DIM], cr[DIM];
  realtype adv[DIM], diff[DIM];
  realtype xmin[DIM], x[DIM], x1;
  int i[DIM], l_m[DIM], m_start[DIM], nbr_left[DIM], nbr_right[DIM], id;
  ProblemData d;
  int dim;

  d = (ProblemData) user_data;

  /* Extract stuff from data structure */
  id = d->myId;
  FOR_DIM {
    xmin[dim]      = d->xmin[dim];
    l_m[dim]       = d->l_m[dim];
    m_start[dim]   = d->m_start[dim];
    dx[dim]        = d->dx[dim];
    nbr_left[dim]  = d->nbr_left[dim];
    nbr_right[dim] = d->nbr_right[dim];
  } 

  /* Get pointers to vector data */
  dydata = N_VGetArrayPointer_Parallel(ydot);
  pdata  = N_VGetArrayPointer_Parallel(d->p);

  /* Copy local segment of y to y_ext */
  Load_yext(N_VGetArrayPointer_Parallel(y), d);
  Ydata = d->y_ext;

  /* Velocity components in x1 and x2 directions (Poiseuille profile) */
  v[1] = ZERO;
#ifdef USE3D
  v[2] = ZERO;
#endif

  /* Local domain is [xmin+(m_start+1)*dx, xmin+(m_start+1+l_m-1)*dx] */
#ifdef USE3D
  for(i[2]=0; i[2]<l_m[2]; i[2]++) {

    x[2] = xmin[2] + (m_start[2]+i[2])*dx[2];
#endif    
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {

      x[1] = xmin[1] + (m_start[1]+i[1])*dx[1];

      /* Velocity component in x0 direction (Poiseuille profile) */
      x1 = x[1] - xmin[1] - L;
      v[0] = V_COEFF * (L + x1) * (L - x1);

      for(i[0]=0; i[0]<l_m[0]; i[0]++) {

        x[0] = xmin[0] + (m_start[0]+i[0])*dx[0];

        c  = IJth_ext(Ydata, i);	       

        /* Source term*/
        IJth(dydata, i) = IJth(pdata, i);

        FOR_DIM {
          i[dim]+=1;
          cr[dim] = IJth_ext(Ydata, i);
          i[dim]-=2;
          cl[dim] = IJth_ext(Ydata, i);
          i[dim]+=1;

          /* Boundary conditions for the state variables */
          if( i[dim]==l_m[dim]-1 && nbr_right[dim]==id)
            cr[dim] = cl[dim];
          else if( i[dim]==0 && nbr_left[dim]==id )
            cl[dim] = cr[dim];

          adv[dim]  = v[dim] * (cr[dim]-cl[dim]) / (TWO*dx[dim]);
          diff[dim] = DIFF_COEF * (cr[dim]-TWO*c+cl[dim]) / SUNSQR(dx[dim]);

          IJth(dydata, i) += (diff[dim] - adv[dim]);
        } 
      }
    }
#ifdef USE3D
  }
#endif

  return(0);
}

/*
 *------------------------------------------------------------------
 * fQ:
 * Right-hand side of quadrature equations on forward integration.
 * The only quadrature on this phase computes the local contribution
 * to the function G.
 *------------------------------------------------------------------
 */

static int fQ(realtype t, N_Vector y, N_Vector qdot, void *user_data)
{
  ProblemData d;
  realtype *dqdata;

  d = (ProblemData) user_data;

  dqdata = N_VGetArrayPointer_Parallel(qdot);

  dqdata[0] = N_VDotProd_Parallel(y,y);
  dqdata[0] *= RCONST(0.5) * (d->dOmega);

  return(0);
}

/*
 *------------------------------------------------------------------
 * fB and fB_local:
 * Backward phase ODE right-hand side (the discretized adjoint PDE)
 *------------------------------------------------------------------
 */

static int fB(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, 
              void *user_dataB)
{
  ProblemData d;
  sunindextype l_neq=1;
  int dim;

  d = (ProblemData) user_dataB;
  FOR_DIM l_neq *= d->l_m[dim];
  
  /* Do all inter-processor communication */
  f_comm(l_neq, t, yB, user_dataB);

  /* Compute right-hand side locally */
  fB_local(l_neq, t, y, yB, yBdot, user_dataB);

  return(0);
}

static int fB_local(sunindextype NlocalB, realtype t, 
                    N_Vector y, N_Vector yB, N_Vector dyB, 
                    void *user_dataB)
{
  realtype *YBdata, *dyBdata, *ydata;
  realtype dx[DIM], c, v[DIM], cl[DIM], cr[DIM];
  realtype adv[DIM], diff[DIM];
  realtype xmin[DIM], x[DIM], x1;
  int i[DIM], l_m[DIM], m_start[DIM], nbr_left[DIM], nbr_right[DIM], id;
  ProblemData d;
  int dim;
  
  d = (ProblemData) user_dataB;

  /* Extract stuff from data structure */
  id = d->myId;
  FOR_DIM {
    xmin[dim]      = d->xmin[dim];
    l_m[dim]       = d->l_m[dim];
    m_start[dim]   = d->m_start[dim];
    dx[dim]        = d->dx[dim];
    nbr_left[dim]  = d->nbr_left[dim];
    nbr_right[dim] = d->nbr_right[dim];
  }
 
  dyBdata = N_VGetArrayPointer_Parallel(dyB);
  ydata   = N_VGetArrayPointer_Parallel(y);

  /* Copy local segment of yB to y_ext */
  Load_yext(N_VGetArrayPointer_Parallel(yB), d);
  YBdata = d->y_ext;

  /* Velocity components in x1 and x2 directions (Poiseuille profile) */
  v[1] = ZERO;
#ifdef USE3D
  v[2] = ZERO;
#endif
 
  /* local domain is [xmin+(m_start)*dx, xmin+(m_start+l_m-1)*dx] */
#ifdef USE3D
  for(i[2]=0; i[2]<l_m[2]; i[2]++) {

    x[2] = xmin[2] + (m_start[2]+i[2])*dx[2];
#endif
    
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {
      
      x[1] = xmin[1] + (m_start[1]+i[1])*dx[1];
	  
      /* Velocity component in x0 direction (Poiseuille profile) */
      x1 = x[1] - xmin[1] - L;
      v[0] = V_COEFF * (L + x1) * (L - x1);

      for(i[0]=0; i[0]<l_m[0]; i[0]++) {

        x[0] = xmin[0] + (m_start[0]+i[0])*dx[0];
        
        c  = IJth_ext(YBdata, i);	       
        
        /* Source term for adjoint PDE */
        IJth(dyBdata, i) = -IJth(ydata, i);
        
        FOR_DIM {
          
          i[dim]+=1;
          cr[dim] = IJth_ext(YBdata, i);
          i[dim]-=2;
          cl[dim] = IJth_ext(YBdata, i);
          i[dim]+=1;

          /* Boundary conditions for the adjoint variables */
          if( i[dim]==l_m[dim]-1 && nbr_right[dim]==id)
	    cr[dim] = cl[dim]-(TWO*dx[dim]*v[dim]/DIFF_COEF)*c;
          else if( i[dim]==0 && nbr_left[dim]==id )
	      cl[dim] = cr[dim]+(TWO*dx[dim]*v[dim]/DIFF_COEF)*c;
		  
          adv[dim]  = v[dim] * (cr[dim]-cl[dim]) / (TWO*dx[dim]);
          diff[dim] = DIFF_COEF * (cr[dim]-TWO*c+cl[dim]) / SUNSQR(dx[dim]);
          
          IJth(dyBdata, i) -= (diff[dim] + adv[dim]);
        } 
      }
    }
#ifdef USE3D
  }
#endif

  return(0);
}

/*
 *------------------------------------------------------------------
 * fQB:
 * Right-hand side of quadrature equations on backward integration
 * The i-th component of the gradient is nothing but int_t yB_i dt
 *------------------------------------------------------------------
 */

static int fQB(realtype t, N_Vector y, N_Vector yB, N_Vector qBdot, 
               void *user_dataB)
{
  ProblemData d;

  d = (ProblemData) user_dataB;

  N_VScale_Parallel(-(d->dOmega), yB, qBdot);

  return(0);
}

/*
 *------------------------------------------------------------------
 * Load_yext: 
 * copies data from src (y or yB) into y_ext, which already contains
 * data from neighboring processes.
 *------------------------------------------------------------------
 */

static void Load_yext(realtype *src, ProblemData d)
{
  int i[DIM], l_m[DIM], dim;

  FOR_DIM l_m[dim] = d->l_m[dim];
     
  /* copy local segment */
#ifdef USE3D
  for  (i[2]=0; i[2]<l_m[2]; i[2]++)
#endif
    for(i[1]=0; i[1]<l_m[1]; i[1]++)
      for(i[0]=0; i[0]<l_m[0]; i[0]++)
	IJth_ext(d->y_ext, i) = IJth(src, i);
}

/*
 *------------------------------------------------------------------
 * PrintHeader:
 * Print first lins of output (problem description)
 *------------------------------------------------------------------
 */

static void PrintHeader()
{
    printf("\nParallel Krylov adjoint sensitivity analysis example\n");
    printf("%1dD Advection diffusion PDE with homogeneous Neumann B.C.\n",DIM);
    printf("Computes gradient of G = int_t_Omega ( c_i^2 ) dt dOmega\n");
    printf("with respect to the source values at each grid point.\n\n");

    printf("Domain:\n");

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("   %Lf < x < %Lf   mx = %d  npe_x = %d \n",XMIN,XMAX,MX,NPX);
    printf("   %Lf < y < %Lf   my = %d  npe_y = %d \n",YMIN,YMAX,MY,NPY);
#else
    printf("   %f < x < %f   mx = %d  npe_x = %d \n",XMIN,XMAX,MX,NPX);
    printf("   %f < y < %f   my = %d  npe_y = %d \n",YMIN,YMAX,MY,NPY);
#endif

#ifdef USE3D
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("   %Lf < z < %Lf   mz = %d  npe_z = %d \n",ZMIN,ZMAX,MZ,NPZ);
#else
    printf("   %f < z < %f   mz = %d  npe_z = %d \n",ZMIN,ZMAX,MZ,NPZ);
#endif
#endif

    printf("\n");
  }

/*
 *------------------------------------------------------------------
 * PrintFinalStats:
 * Print final statistics contained in cvode_mem
 *------------------------------------------------------------------
 */

static int PrintFinalStats(void *cvode_mem)
{
  long int lenrw, leniw;
  long int lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int retval;

  retval = CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
  if(check_retval(&retval, "CVodeGetWorkSpace", 1, 0)) return(-1);
  
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  if(check_retval(&retval, "CVodeGetNumSteps", 1, 0)) return(-1);
  
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  if(check_retval(&retval, "CVodeGetNumRhsEvals", 1, 0)) return(-1);
  
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  if(check_retval(&retval, "CVodeGetNumLinSolvSetups", 1, 0)) return(-1);
  
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  if(check_retval(&retval, "CVodeGetNumErrTestFails", 1, 0)) return(-1);
  
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  if(check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1, 0)) return(-1);
  
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  if(check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1, 0)) return(-1);

  retval = CVodeGetLinWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
  if(check_retval(&retval, "CVodeGetLinWorkSpace", 1, 0)) return(-1);
  
  retval = CVodeGetNumLinIters(cvode_mem, &nli);
  if(check_retval(&retval, "CVodeGetNumLinIters", 1, 0)) return(-1);
  
  retval = CVodeGetNumPrecEvals(cvode_mem, &npe);
  if(check_retval(&retval, "CVodeGetNumPrecEvals", 1, 0)) return(-1);
  
  retval = CVodeGetNumPrecSolves(cvode_mem, &nps);
  if(check_retval(&retval, "CVodeGetNumPrecSolves", 1, 0)) return(-1);
  
  retval = CVodeGetNumLinConvFails(cvode_mem, &ncfl);
  if(check_retval(&retval, "CVodeGetNumLinConvFails", 1, 0)) return(-1);
  
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  if(check_retval(&retval, "CVodeGetNumLinRhsEvals", 1, 0)) return(-1);

  printf("\nFinal Statistics.. \n\n");
  printf("lenrw   = %6ld     leniw = %6ld\n", lenrw, leniw);
  printf("llrw    = %6ld     lliw  = %6ld\n", lenrwLS, leniwLS);
  printf("nst     = %6ld\n"                  , nst);
  printf("nfe     = %6ld     nfel  = %6ld\n"  , nfe, nfeLS);
  printf("nni     = %6ld     nli   = %6ld\n"  , nni, nli);
  printf("nsetups = %6ld     netf  = %6ld\n"  , nsetups, netf);
  printf("npe     = %6ld     nps   = %6ld\n"  , npe, nps);
  printf("ncfn    = %6ld     ncfl  = %6ld\n\n", ncfn, ncfl); 

  return(0);
}

/*
 *------------------------------------------------------------------
 * OutputGradient:
 * Generate matlab m files for visualization
 * One file gradXXXX.m from each process + a driver grad.m
 *------------------------------------------------------------------
 */

static void OutputGradient(int myId, N_Vector qB, ProblemData d)
{
  FILE *fid;
  char filename[20];
  int *l_m, *m_start, i[DIM],ip;
  realtype *xmin, *dx;
  realtype x[DIM], *pdata, p, *qBdata, g;

  sprintf(filename,"grad%03d.m",myId);
  fid = fopen(filename,"w");

  l_m  = d->l_m;
  m_start = d->m_start;
  xmin = d->xmin;
  dx = d->dx;

  qBdata = N_VGetArrayPointer_Parallel(qB);
  pdata  = N_VGetArrayPointer_Parallel(d->p);

  /* Write matlab files with solutions from each process */

  /*   Allocate Matlab storage for data */

  fprintf(fid,"x%d = zeros(%d,1); \n",  myId, l_m[0]);
  fprintf(fid,"y%d = zeros(%d,1); \n",  myId, l_m[1]);
#ifdef USE3D
  fprintf(fid,"z%d = zeros(%d,1); \n",  myId, l_m[2]);
  fprintf(fid,"p%d = zeros(%d,%d,%d); \n", myId, l_m[1], l_m[0], l_m[2]);
  fprintf(fid,"g%d = zeros(%d,%d,%d); \n", myId, l_m[1], l_m[0], l_m[2]);
#else
  fprintf(fid,"p%d = zeros(%d,%d); \n", myId, l_m[1], l_m[0]);
  fprintf(fid,"g%d = zeros(%d,%d); \n", myId, l_m[1], l_m[0]);
#endif

  /*   Write mesh information */

  for(i[0]=0; i[0]<l_m[0]; i[0]++) {
    x[0] = xmin[0] + (m_start[0]+i[0]) * dx[0];
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(fid,"x%d(%d,1) = %Le; \n",  myId, i[0]+1, x[0]);
#else
    fprintf(fid,"x%d(%d,1) = %e; \n",  myId, i[0]+1, x[0]);
#endif
  }

  for(i[1]=0; i[1]<l_m[1]; i[1]++) {
    x[1] = xmin[1] + (m_start[1]+i[1]) * dx[1];
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(fid,"y%d(%d,1) = %Le; \n",  myId, i[1]+1, x[1]);
#else
    fprintf(fid,"y%d(%d,1) = %e; \n",  myId, i[1]+1, x[1]);
#endif
  }

#ifdef USE3D
  for(i[2]=0; i[2]<l_m[2]; i[2]++) {
    x[2] = xmin[2] + (m_start[2]+i[2]) * dx[2];
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(fid,"z%d(%d,1) = %Le; \n",  myId, i[2]+1, x[2]);
#else
    fprintf(fid,"z%d(%d,1) = %e; \n",  myId, i[2]+1, x[2]);
#endif
  }
#endif

  /*   Write solution data */

  for(i[0]=0; i[0]<l_m[0]; i[0]++) {
    x[0] = xmin[0] + (m_start[0]+i[0]) * dx[0];
    for(i[1]=0; i[1]<l_m[1]; i[1]++) {
      x[1] = xmin[1] + (m_start[1]+i[1]) * dx[1];
#ifdef USE3D
      for(i[2]=0; i[2]<l_m[2]; i[2]++) {
        x[2] = xmin[2] + (m_start[2]+i[2]) * dx[2];
        g = IJth(qBdata, i);
        p = IJth(pdata, i);
#if defined(SUNDIALS_EXTENDED_PRECISION)
        fprintf(fid,"p%d(%d,%d,%d) = %Le; \n", myId, i[1]+1, i[0]+1, i[2]+1, p);
        fprintf(fid,"g%d(%d,%d,%d) = %Le; \n", myId, i[1]+1, i[0]+1, i[2]+1, g);
#else
        fprintf(fid,"p%d(%d,%d,%d) = %e; \n", myId, i[1]+1, i[0]+1, i[2]+1, p);
        fprintf(fid,"g%d(%d,%d,%d) = %e; \n", myId, i[1]+1, i[0]+1, i[2]+1, g);
#endif
      }
#else
      g = IJth(qBdata, i);
      p = IJth(pdata, i);
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(fid,"p%d(%d,%d) = %Le; \n", myId, i[1]+1, i[0]+1, p);
      fprintf(fid,"g%d(%d,%d) = %Le; \n", myId, i[1]+1, i[0]+1, g);
#else
      fprintf(fid,"p%d(%d,%d) = %e; \n", myId, i[1]+1, i[0]+1, p);
      fprintf(fid,"g%d(%d,%d) = %e; \n", myId, i[1]+1, i[0]+1, g);
#endif
#endif 
    }
  }
  fclose(fid);


  /* Write matlab driver */

  if (myId == 0) {

    fid = fopen("grad.m","w");

#ifdef USE3D
    fprintf(fid,"clear;\n");
    fprintf(fid,"figure(1);\nhold on\n");
    fprintf(fid,"figure(2);\nhold on\n");
    fprintf(fid,"trans = 0.7;\n");
    fprintf(fid,"ecol  = 'none';\n");
    fprintf(fid,"glev1 = 0.4;\n");
    fprintf(fid,"glev2 = 0.25;\n");
    fprintf(fid,"gcol1 = 'blue';\n");
    fprintf(fid,"gcol2 = 'green';\n");
    fprintf(fid,"gtrans = 0.5;\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(fid,"xp=[%Lf %Lf];\n",G1_X,G2_X);
    fprintf(fid,"yp=[%Lf];\n",G2_Y);
    fprintf(fid,"zp=[%Lf];\n",G1_Z);
#else
    fprintf(fid,"xp=[%f %f];\n",G1_X,G2_X);
    fprintf(fid,"yp=[%f];\n",G2_Y);
    fprintf(fid,"zp=[];\n");
#endif
    for (ip=0; ip<d->npes; ip++) {
      fprintf(fid,"\ngrad%03d;\n",ip);
      fprintf(fid,"figure(1)\n");
      fprintf(fid,"[X,Y,Z]=meshgrid(x%d,y%d,z%d);\n",ip,ip,ip);
      fprintf(fid,"s%d=slice(X,Y,Z,p%d,xp,yp,zp);\n",ip,ip);
      fprintf(fid,"for i = 1:length(s%d)\n",ip);
      fprintf(fid,"  set(s%d(i),'FaceAlpha',trans);\n",ip);
      fprintf(fid,"  set(s%d(i),'EdgeColor',ecol);\n",ip);
      fprintf(fid,"end\n");

      fprintf(fid,"\nfigure(2)\n");
      fprintf(fid,"p=patch(isosurface(X,Y,Z,g%d,glev1));\n",ip);
      fprintf(fid,"p.FaceColor = gcol1;\n");
      fprintf(fid,"p.EdgeColor = ecol;\n");
      fprintf(fid,"p=patch(isosurface(X,Y,Z,g%d,glev2));\n",ip);
      fprintf(fid,"p.FaceColor = gcol2;\n");
      fprintf(fid,"p.EdgeColor = ecol;\n");
      fprintf(fid,"p.FaceAlpha = gtrans;\n");
      fprintf(fid,"clear x%d y%d z%d p%d g%d;\n",ip,ip,ip,ip,ip);
    }
    
    fprintf(fid,"\nfigure(1)\n");
    fprintf(fid,"view(3)\n");
    fprintf(fid,"shading interp\naxis equal\n");
    fprintf(fid,"hold off\n");
    fprintf(fid,"xlabel('x')\n");
    fprintf(fid,"ylabel('y')\n");
    fprintf(fid,"zlabel('z')\n");
    fprintf(fid,"axis([%f, %f, %f, %f, %f, %f])\n",XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX);
    fprintf(fid,"print('cvsadjkryx_p3Dcf','-depsc')\n");
    fprintf(fid,"savefig('cvsadjkryx_p3Dcf.fig')\n");

    fprintf(fid,"\nfigure(2)\n");
    fprintf(fid,"view(3)\n");
    fprintf(fid,"axis equal\n");
    fprintf(fid,"hold off\n");
    fprintf(fid,"xlabel('x')\n");
    fprintf(fid,"ylabel('y')\n");
    fprintf(fid,"zlabel('z')\n");
    fprintf(fid,"axis([%f, %f, %f, %f, %f, %f])\n",XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX);
    fprintf(fid,"camlight\n");
    fprintf(fid,"lighting gouraud\n");
    fprintf(fid,"print('cvsadjkryx_p3Dgrad','-depsc')\n");
    fprintf(fid,"savefig('cvsadjkryx_p3Dgrad.fig')\n");
#else
    fprintf(fid,"clear;\n");
    fprintf(fid,"figure('units','normalized','position',[.1 .1 .6 .4])\n");
    fprintf(fid,"trans = 0.7;\n");
    fprintf(fid,"ecol  = 'none';\n");

    for (ip=0; ip<d->npes; ip++) {

      fprintf(fid,"\ngrad%03d;\n",ip);

      fprintf(fid,"\nax(1) = subplot(1,2,1);\n");
      fprintf(fid,"s = surf(x%d,y%d,g%d);\n",ip,ip,ip);
      fprintf(fid,"set(s, 'FaceAlpha', trans);\n");
      fprintf(fid,"set(s, 'EdgeColor', ecol);\n");
      fprintf(fid,"hold on\n");
      fprintf(fid,"axis tight\n");
      fprintf(fid,"box on\n");
      fprintf(fid,"colorbar('Position', [0.5 0.1 0.025 0.8])\n");
      
      fprintf(fid,"\nax(2) = subplot(1,2,2);\n");
      fprintf(fid,"s = surf(x%d,y%d,p%d);\n",ip,ip,ip);
      fprintf(fid,"set(s, 'CData', g%d);\n",ip);
      fprintf(fid,"set(s, 'FaceAlpha', trans);\n");
      fprintf(fid,"set(s, 'EdgeColor', ecol);\n");
      fprintf(fid,"hold on\n");
      fprintf(fid,"axis tight\n");
      fprintf(fid,"box on\n");

      fprintf(fid,"clear x%d y%d p%d g%d;\n",ip,ip,ip,ip);
    }

    fprintf(fid,"\nax(1) = subplot(1,2,1);\n");
    fprintf(fid,"pos = get(ax(1), 'Position');\n");
    fprintf(fid,"set(ax(1), 'Position', [pos(1)-0.02 pos(2) pos(3) pos(4)]);\n");
    fprintf(fid,"xlabel('x'), ylabel('y')\n");
    fprintf(fid,"hold off\n");

    fprintf(fid,"\nax(2) = subplot(1,2,2);\n");
    fprintf(fid,"pos = get(ax(2), 'Position');\n");
    fprintf(fid,"set(ax(2), 'Position', [pos(1)+0.02 pos(2) pos(3) pos(4)]);\n");
    fprintf(fid,"xlabel('x'), ylabel('y')\n");
    fprintf(fid,"hold off\n");

    fprintf(fid,"\nfig = gcf;\n");
    fprintf(fid,"fig.PaperPositionMode = 'auto';\n");
    fprintf(fid,"print('cvsadjkryx_p2D','-depsc','-r0')\n");
    fprintf(fid,"savefig('cvsadjkryx_p2D.fig')\n");
#endif
    fclose(fid);
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

static int check_retval(void *returnvalue, const char *funcname, int opt, int id)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n",
	    id, funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR(%d): %s() failed with retval = %d\n\n",
	      id, funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n",
	    id, funcname);
    return(1); }

  return(0);
}
