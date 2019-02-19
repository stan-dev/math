/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen and Alan C. Hindmarsh and
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
 * Example problem:
 *
 * An ODE system is generated from the following 2-species diurnal
 * kinetics advection-diffusion PDE system in 2 space dimensions:
 *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)
 *                 + Ri(c1,c2,t)      for i = 1,2,   where
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
 *   Kv(z) = Kv0*exp(z/5) ,
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
 * vary diurnally. The problem is posed on the square
 *   0 <= x <= 20,    30 <= z <= 50   (all in km),
 * with homogeneous Neumann boundary conditions, and for time t in
 *   0 <= t <= 86400 sec (1 day).
 * The PDE system is treated by central differences on a uniform
 * 10 x 10 mesh, with simple polynomial initial profiles.
 * The problem is solved with CVODES, with the BDF/GMRES method
 * (i.e. using the SUNLinSol_SPGMR linear solver) and the block-diagonal
 * part of the Newton matrix as a left preconditioner. A copy of
 * the block-diagonal part of the Jacobian is saved and
 * conditionally reused within the Precond routine.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters q1 and q2.
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on FULL or PARTIAL,
 * respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % cvsDiurnal_FSA_kry -nosensi
 * If sensitivities are to be computed:
 *    % cvsDiurnal_FSA_kry -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cvodes/cvodes.h>             /* main CVODES header file              */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver      */
#include <sundials/sundials_dense.h>   /* use generic dense solver in precond. */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* contains the macros ABS, SUNSQR, EXP */

/* Problem Constants */

#define NUM_SPECIES  2                /* number of species */
#define C1_SCALE     RCONST(1.0e6)    /* coefficients in initial profiles */
#define C2_SCALE     RCONST(1.0e12)

#define T0           RCONST(0.0)      /* initial time */
#define NOUT         12               /* number of output times */
#define TWOHR        RCONST(7200.0)   /* number of seconds in two hours  */
#define HALFDAY      RCONST(4.32e4)   /* number of seconds in a half day */
#define PI           RCONST(3.1415926535898)   /* pi */ 

#define XMIN         RCONST(0.0)      /* grid boundaries in x  */
#define XMAX         RCONST(20.0)           
#define ZMIN         RCONST(30.0)     /* grid boundaries in z  */
#define ZMAX         RCONST(50.0)
#define XMID         RCONST(10.0)     /* grid midpoints in x,z */          
#define ZMID         RCONST(40.0)

#define MX           15               /* MX = number of x mesh points */
#define MZ           15               /* MZ = number of z mesh points */
#define NSMX         NUM_SPECIES*MX   /* NSMX = NUM_SPECIES*MX */
#define MM           (MX*MZ)          /* MM = MX*MZ */

/* CVodeInit Constants */
#define RTOL         RCONST(1.0e-5)   /* scalar relative tolerance */
#define FLOOR        RCONST(100.0)    /* value of C1 or C2 at which tolerances */
                                      /* change from relative to absolute      */
#define ATOL         (RTOL*FLOOR)     /* scalar absolute tolerance */
#define NEQ          (NUM_SPECIES*MM) /* NEQ = number of equations */

/* Sensitivity Constants */
#define NP           8
#define NS           2

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* User-defined vector and matrix accessor macros: IJKth, IJth */

/* IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into small dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.   
   
   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MZ-1. The vdata array is obtained via
   the call vdata = N_VGetArrayPointer(v), where v is an N_Vector. 
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.

   IJth(a,i,j) references the (i,j)th entry of the small matrix realtype **a,
   where 1 <= i,j <= NUM_SPECIES. The small matrix routines in dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. */

#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMX])
#define IJth(a,i,j)        (a[j-1][i-1])

/* Type : UserData 
   contains preconditioner blocks, pivot arrays, 
   problem parameters, and problem constants     */

typedef struct {
  realtype *p;
  realtype **P[MX][MZ], **Jbd[MX][MZ];
  sunindextype *pivot[MX][MZ];
  realtype q4, om, dx, dz, hdco, haco, vdco;
} *UserData;


/* Prototypes of user-supplied functions */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Precond(realtype tn, N_Vector y, N_Vector fy, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, void *user_data);

static int PSolve(realtype tn, N_Vector y, N_Vector fy,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *user_data);

/* Prototypes of private functions */

static void ProcessArgs(int argc, char *argv[],
                        booleantype *sensi, int *sensi_meth, booleantype *err_con);
static void WrongArgs(char *name);
static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data);
static void SetInitialProfiles(N_Vector y, realtype dx, realtype dz);
static void PrintOutput(void *cvode_mem, realtype t, N_Vector y);
static void PrintOutputS(N_Vector *uS);
static void PrintFinalStats(void *cvode_mem, booleantype sensi,
                            booleantype err_con, int sensi_meth);
static int check_retval(void *returnvalue, const char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  void *cvode_mem;
  SUNLinearSolver LS;
  UserData data;
  realtype abstol, reltol, t, tout;
  N_Vector y;
  int iout, retval;

  realtype *pbar;
  int is, *plist;
  N_Vector *uS;
  booleantype sensi, err_con;
  int sensi_meth;

  pbar = NULL;
  plist = NULL;
  uS = NULL;
  y = NULL;
  data = NULL;
  cvode_mem = NULL;
  LS = NULL;

  /* Process arguments */
  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);

  /* Problem parameters */
  data = AllocUserData();
  if(check_retval((void *)data, "AllocUserData", 2)) return(1);
  InitUserData(data);

  /* Initial states */
  y = N_VNew_Serial(NEQ);
  if(check_retval((void *)y, "N_VNew_Serial", 0)) return(1);
  SetInitialProfiles(y, data->dx, data->dz);
  
  /* Tolerances */
  abstol=ATOL; 
  reltol=RTOL;

  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_BDF);
  if(check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  retval = CVodeSetUserData(cvode_mem, data);
  if(check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  retval = CVodeSetMaxNumSteps(cvode_mem, 2000);
  if(check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return(1);

  /* Allocate CVODES memory */
  retval = CVodeInit(cvode_mem, f, T0, y);
  if(check_retval(&retval, "CVodeInit", 1)) return(1);

  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  if(check_retval(&retval, "CVodeSStolerances", 1)) return(1);

  /* Create the SUNLinSol_SPGMR linear solver with left
     preconditioning and the default Krylov dimension */
  LS = SUNLinSol_SPGMR(y, PREC_LEFT, 0);
  if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

  /* Attach the linear sovler */
  retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
  if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

  /* Set the preconditioner solve and setup functions */
  retval = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
  if(check_retval(&retval, "CVodeSetPreconditioner", 1)) return(1);

  printf("\n2-species diurnal advection-diffusion problem\n");

  /* Forward sensitivity analysis */
  if(sensi) {

    plist = (int *) malloc(NS * sizeof(int));
    if(check_retval((void *)plist, "malloc", 2)) return(1);
    for(is=0; is<NS; is++) plist[is] = is;

    pbar = (realtype *) malloc(NS * sizeof(realtype));
    if(check_retval((void *)pbar, "malloc", 2)) return(1);
    for(is=0; is<NS; is++) pbar[is] = data->p[plist[is]];

    uS = N_VCloneVectorArray(NS, y);
    if(check_retval((void *)uS, "N_VCloneVectorArray", 0)) return(1);
    for(is=0;is<NS;is++)
      N_VConst(ZERO,uS[is]);

    retval = CVodeSensInit1(cvode_mem, NS, sensi_meth, NULL, uS);
    if(check_retval(&retval, "CVodeSensInit", 1)) return(1);

    retval = CVodeSensEEtolerances(cvode_mem);
    if(check_retval(&retval, "CVodeSensEEtolerances", 1)) return(1);

    retval = CVodeSetSensErrCon(cvode_mem, err_con);
    if(check_retval(&retval, "CVodeSetSensErrCon", 1)) return(1);

    retval = CVodeSetSensDQMethod(cvode_mem, CV_CENTERED, ZERO);
    if(check_retval(&retval, "CVodeSetSensDQMethod", 1)) return(1);

    retval = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);
    if(check_retval(&retval, "CVodeSetSensParams", 1)) return(1);

    printf("Sensitivity: YES ");
    if(sensi_meth == CV_SIMULTANEOUS)   
      printf("( SIMULTANEOUS +");
    else 
      if(sensi_meth == CV_STAGGERED) printf("( STAGGERED +");
      else                           printf("( STAGGERED1 +");   
    if(err_con) printf(" FULL ERROR CONTROL )");
    else        printf(" PARTIAL ERROR CONTROL )");
    
  } else {

    printf("Sensitivity: NO ");

  }

  /* In loop over output points, call CVode, print results, test for error */

  printf("\n\n");
  printf("========================================================================\n");
  printf("     T     Q       H      NST                    Bottom left  Top right \n");
  printf("========================================================================\n");

  for (iout=1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if(check_retval(&retval, "CVode", 1)) break;
    PrintOutput(cvode_mem, t, y);
    if (sensi) {
      retval = CVodeGetSens(cvode_mem, &t, uS);
      if(check_retval(&retval, "CVodeGetSens", 1)) break;
      PrintOutputS(uS);
    }
    
    printf("------------------------------------------------------------------------\n");

  }

  /* Print final statistics */
  PrintFinalStats(cvode_mem, sensi, err_con, sensi_meth);

  /* Free memory */
  N_VDestroy(y);
  if (sensi) {
    N_VDestroyVectorArray(uS, NS);
    free(pbar);
    free(plist);
  }
  FreeUserData(data);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y). 
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  realtype c1rt, c2rt, czdn, czup, hord1, hord2, horad1, horad2;
  realtype qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, zdn, zup;
  realtype q4coef, delz, verdco, hordco, horaco;
  realtype *ydata, *dydata;
  int jx, jz, idn, iup, ileft, iright;
  UserData data;
  realtype Q1, Q2, C3, A3, A4;

  data = (UserData) user_data;
  ydata = N_VGetArrayPointer(y);
  dydata = N_VGetArrayPointer(ydot);

  /* Load problem coefficients and parameters */

  Q1 = data->p[0];
  Q2 = data->p[1];
  C3 = data->p[2];
  A3 = data->p[3];
  A4 = data->p[4];

  /* Set diurnal rate coefficients. */

  s = sin(data->om*t);
  if (s > ZERO) {
    q3 = SUNRexp(-A3/s);
    data->q4 = SUNRexp(-A4/s);
  } else {
    q3 = ZERO;
    data->q4 = ZERO;
  }

  /* Make local copies of problem variables, for efficiency. */

  q4coef = data->q4;
  delz = data->dz;
  verdco = data->vdco;
  hordco  = data->hdco;
  horaco  = data->haco;

  /* Loop over all grid points. */

  for (jz=0; jz < MZ; jz++) {

    /* Set vertical diffusion coefficients at jz +- 1/2 */

    zdn = ZMIN + (jz - RCONST(0.5))*delz;
    zup = zdn + delz;
    czdn = verdco*SUNRexp(RCONST(0.2)*zdn);
    czup = verdco*SUNRexp(RCONST(0.2)*zup);
    idn = (jz == 0) ? 1 : -1;
    iup = (jz == MZ-1) ? -1 : 1;
    for (jx=0; jx < MX; jx++) {

      /* Extract c1 and c2, and set kinetic rate terms. */

      c1 = IJKth(ydata,1,jx,jz); 
      c2 = IJKth(ydata,2,jx,jz);
      qq1 = Q1*c1*C3;
      qq2 = Q2*c1*c2;
      qq3 = q3*C3;
      qq4 = q4coef*c2;
      rkin1 = -qq1 - qq2 + RCONST(2.0)*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms. */

      c1dn = IJKth(ydata,1,jx,jz+idn);
      c2dn = IJKth(ydata,2,jx,jz+idn);
      c1up = IJKth(ydata,1,jx,jz+iup);
      c2up = IJKth(ydata,2,jx,jz+iup);
      vertd1 = czup*(c1up - c1) - czdn*(c1 - c1dn);
      vertd2 = czup*(c2up - c2) - czdn*(c2 - c2dn);

      /* Set horizontal diffusion and advection terms. */

      ileft = (jx == 0) ? 1 : -1;
      iright =(jx == MX-1) ? -1 : 1;
      c1lt = IJKth(ydata,1,jx+ileft,jz); 
      c2lt = IJKth(ydata,2,jx+ileft,jz);
      c1rt = IJKth(ydata,1,jx+iright,jz);
      c2rt = IJKth(ydata,2,jx+iright,jz);
      hord1 = hordco*(c1rt - RCONST(2.0)*c1 + c1lt);
      hord2 = hordco*(c2rt - RCONST(2.0)*c2 + c2lt);
      horad1 = horaco*(c1rt - c1lt);
      horad2 = horaco*(c2rt - c2lt);

      /* Load all terms into ydot. */

      IJKth(dydata, 1, jx, jz) = vertd1 + hord1 + horad1 + rkin1; 
      IJKth(dydata, 2, jx, jz) = vertd2 + hord2 + horad2 + rkin2;
    }
  }

  return(0);
}

/*
 * Preconditioner setup routine. Generate and preprocess P. 
 */

static int Precond(realtype tn, N_Vector y, N_Vector fy, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, void *user_data)
{
  realtype c1, c2, czdn, czup, diag, zdn, zup, q4coef, delz, verdco, hordco;
  realtype **(*P)[MZ], **(*Jbd)[MZ];
  sunindextype *(*pivot)[MZ];
  int retval, jx, jz;
  realtype *ydata, **a, **j;
  UserData data;
  realtype Q1, Q2, C3;

  /* Make local copies of pointers in user_data, and of pointer to y's data */
  data = (UserData) user_data;
  P = data->P;
  Jbd = data->Jbd;
  pivot = data->pivot;
  ydata = N_VGetArrayPointer(y);

  /* Load problem coefficients and parameters */
  Q1 = data->p[0];
  Q2 = data->p[1];
  C3 = data->p[2];

  if (jok) {

  /* jok = SUNTRUE: Copy Jbd to P */

    for (jz=0; jz < MZ; jz++)
      for (jx=0; jx < MX; jx++)
        denseCopy(Jbd[jx][jz], P[jx][jz], NUM_SPECIES, NUM_SPECIES);

  *jcurPtr = SUNFALSE;

  }

  else {
  /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */

  /* Make local copies of problem variables, for efficiency. */

  q4coef = data->q4;
  delz = data->dz;
  verdco = data->vdco;
  hordco  = data->hdco;

  /* Compute 2x2 diagonal Jacobian blocks (using q4 values 
     computed on the last f call).  Load into P. */

    for (jz=0; jz < MZ; jz++) {
      zdn = ZMIN + (jz - RCONST(0.5))*delz;
      zup = zdn + delz;
      czdn = verdco*SUNRexp(RCONST(0.2)*zdn);
      czup = verdco*SUNRexp(RCONST(0.2)*zup);
      diag = -(czdn + czup + RCONST(2.0)*hordco);
      for (jx=0; jx < MX; jx++) {
        c1 = IJKth(ydata,1,jx,jz);
        c2 = IJKth(ydata,2,jx,jz);
        j = Jbd[jx][jz];
        a = P[jx][jz];
        IJth(j,1,1) = (-Q1*C3 - Q2*c2) + diag;
        IJth(j,1,2) = -Q2*c1 + q4coef;
        IJth(j,2,1) = Q1*C3 - Q2*c2;
        IJth(j,2,2) = (-Q2*c1 - q4coef) + diag;
        denseCopy(j, a, NUM_SPECIES, NUM_SPECIES);
      }
    }

  *jcurPtr = SUNTRUE;

  }

  /* Scale by -gamma */

    for (jz=0; jz < MZ; jz++)
      for (jx=0; jx < MX; jx++)
        denseScale(-gamma, P[jx][jz], NUM_SPECIES, NUM_SPECIES);

  /* Add identity matrix and do LU decompositions on blocks in place. */

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      denseAddIdentity(P[jx][jz], NUM_SPECIES);
      retval = denseGETRF(P[jx][jz], NUM_SPECIES, NUM_SPECIES, pivot[jx][jz]);
      if (retval != 0) return(1);
    }
  }

  return(0);
}

/*
 * Preconditioner solve routine 
 */

static int PSolve(realtype tn, N_Vector y, N_Vector fy,
                  N_Vector r, N_Vector z,
                  realtype gamma, realtype delta,
                  int lr, void *user_data)
{
  realtype **(*P)[MZ];
  sunindextype *(*pivot)[MZ];
  int jx, jz;
  realtype *zdata, *v;
  UserData data;

  /* Extract the P and pivot arrays from user_data. */

  data = (UserData) user_data;
  P = data->P;
  pivot = data->pivot;
  zdata = N_VGetArrayPointer(z);

  N_VScale(ONE, r, z);

  /* Solve the block-diagonal system Px = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      v = &(IJKth(zdata, 1, jx, jz));
      denseGETRS(P[jx][jz], NUM_SPECIES, pivot[jx][jz], v);
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
 * Process and verify arguments to cvsfwdkryx.
 */

static void ProcessArgs(int argc, char *argv[], 
                        booleantype *sensi, int *sensi_meth, booleantype *err_con)
{
  *sensi = SUNFALSE;
  *sensi_meth = -1;
  *err_con = SUNFALSE;

  if (argc < 2) WrongArgs(argv[0]);

  if (strcmp(argv[1],"-nosensi") == 0)
    *sensi = SUNFALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    *sensi = SUNTRUE;
  else
    WrongArgs(argv[0]);
  
  if (*sensi) {

    if (argc != 4)
      WrongArgs(argv[0]);

    if (strcmp(argv[2],"sim") == 0)
      *sensi_meth = CV_SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      *sensi_meth = CV_STAGGERED;
    else if (strcmp(argv[2],"stg1") == 0)
      *sensi_meth = CV_STAGGERED1;
    else 
      WrongArgs(argv[0]);

    if (strcmp(argv[3],"t") == 0)
      *err_con = SUNTRUE;
    else if (strcmp(argv[3],"f") == 0)
      *err_con = SUNFALSE;
    else
      WrongArgs(argv[0]);
  }

}

static void WrongArgs(char *name)
{
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = t or f\n");
    
    exit(0);
}

/*
 * Allocate memory for data structure of type UserData 
 */

static UserData AllocUserData(void)
{
  int jx, jz;
  UserData data;

  data = (UserData) malloc(sizeof *data);

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      (data->P)[jx][jz] = newDenseMat(NUM_SPECIES, NUM_SPECIES);
      (data->Jbd)[jx][jz] = newDenseMat(NUM_SPECIES, NUM_SPECIES);
      (data->pivot)[jx][jz] = newIndexArray(NUM_SPECIES);
    }
  }

  data->p = (realtype *) malloc(NP*sizeof(realtype));

  return(data);
}

/*
 * Load problem constants in data 
 */

static void InitUserData(UserData data)
{
  realtype Q1, Q2, C3, A3, A4, KH, VEL, KV0;

  /* Set problem parameters */
  Q1 = RCONST(1.63e-16); /* Q1  coefficients q1, q2, c3             */
  Q2 = RCONST(4.66e-16); /* Q2                                      */
  C3 = RCONST(3.7e16);   /* C3                                      */
  A3 = RCONST(22.62);    /* A3  coefficient in expression for q3(t) */
  A4 = RCONST(7.601);    /* A4  coefficient in expression for q4(t) */
  KH = RCONST(4.0e-6);   /* KH  horizontal diffusivity Kh           */ 
  VEL = RCONST(0.001);   /* VEL advection velocity V                */
  KV0 = RCONST(1.0e-8);  /* KV0 coefficient in Kv(z)                */  

  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/(MX-1);
  data->dz = (ZMAX-ZMIN)/(MZ-1);
  data->hdco = KH/SUNSQR(data->dx);
  data->haco = VEL/(RCONST(2.0)*data->dx);
  data->vdco = (ONE/SUNSQR(data->dz))*KV0;

  data->p[0] = Q1;
  data->p[1] = Q2;
  data->p[2] = C3;
  data->p[3] = A3;
  data->p[4] = A4;
  data->p[5] = KH;
  data->p[6] = VEL;
  data->p[7] = KV0;
}

/*
 * Free user data memory 
 */

static void FreeUserData(UserData data)
{
  int jx, jz;

  for (jx=0; jx < MX; jx++) {
    for (jz=0; jz < MZ; jz++) {
      destroyMat((data->P)[jx][jz]);
      destroyMat((data->Jbd)[jx][jz]);
      destroyArray((data->pivot)[jx][jz]);
    }
  }

  free(data->p);

  free(data);
}

/*
 * Set initial conditions in y 
 */

static void SetInitialProfiles(N_Vector y, realtype dx, realtype dz)
{
  int jx, jz;
  realtype x, z, cx, cz;
  realtype *ydata;

  /* Set pointer to data array in vector y. */

  ydata = N_VGetArrayPointer(y);

  /* Load initial profiles of c1 and c2 into y vector */

  for (jz=0; jz < MZ; jz++) {
    z = ZMIN + jz*dz;
    cz = SUNSQR(RCONST(0.1)*(z - ZMID));
    cz = ONE - cz + RCONST(0.5)*SUNSQR(cz);
    for (jx=0; jx < MX; jx++) {
      x = XMIN + jx*dx;
      cx = SUNSQR(RCONST(0.1)*(x - XMID));
      cx = ONE - cx + RCONST(0.5)*SUNSQR(cx);
      IJKth(ydata,1,jx,jz) = C1_SCALE*cx*cz; 
      IJKth(ydata,2,jx,jz) = C2_SCALE*cx*cz;
    }
  }
}

/*
 * Print current t, step count, order, stepsize, and sampled c1,c2 values 
 */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector y)
{  
  long int nst;
  int qu, retval;
  realtype hu;
  realtype *ydata;

  ydata = N_VGetArrayPointer(y);

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetLastOrder(cvode_mem, &qu);
  check_retval(&retval, "CVodeGetLastOrder", 1);
  retval = CVodeGetLastStep(cvode_mem, &hu);
  check_retval(&retval, "CVodeGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%8.3Le %2d  %8.3Le %5ld\n", t,qu,hu,nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%8.3e %2d  %8.3e %5ld\n", t,qu,hu,nst);
#else
  printf("%8.3e %2d  %8.3e %5ld\n", t,qu,hu,nst);
#endif

  printf("                                Solution       ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le \n", IJKth(ydata,1,0,0), IJKth(ydata,1,MX-1,MZ-1)); 
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e \n", IJKth(ydata,1,0,0), IJKth(ydata,1,MX-1,MZ-1)); 
#else
  printf("%12.4e %12.4e \n", IJKth(ydata,1,0,0), IJKth(ydata,1,MX-1,MZ-1)); 
#endif
  printf("                                               ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le \n", IJKth(ydata,2,0,0), IJKth(ydata,2,MX-1,MZ-1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e \n", IJKth(ydata,2,0,0), IJKth(ydata,2,MX-1,MZ-1));
#else
  printf("%12.4e %12.4e \n", IJKth(ydata,2,0,0), IJKth(ydata,2,MX-1,MZ-1));
#endif
}

/*
 * Print sampled sensitivities 
 */

static void PrintOutputS(N_Vector *uS)
{
  realtype *sdata;

  sdata = N_VGetArrayPointer(uS[0]);

  printf("                                ----------------------------------------\n"); 
  printf("                                Sensitivity 1  ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
#else
  printf("%12.4e %12.4e \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
#endif
  printf("                                               ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));
#else
  printf("%12.4e %12.4e \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));
#endif

  sdata = N_VGetArrayPointer(uS[1]);

  printf("                                ----------------------------------------\n"); 
  printf("                                Sensitivity 2  ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
#else
  printf("%12.4e %12.4e \n", IJKth(sdata,1,0,0), IJKth(sdata,1,MX-1,MZ-1)); 
#endif
  printf("                                               ");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%12.4Le %12.4Le \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%12.4e %12.4e \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));
#else
  printf("%12.4e %12.4e \n", IJKth(sdata,2,0,0), IJKth(sdata,2,MX-1,MZ-1));
#endif
}

/*
 * Print final statistics contained in iopt 
 */

static void PrintFinalStats(void *cvode_mem, booleantype sensi,
                            booleantype err_con, int sensi_meth)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int nli, ncfl, npe, nps;
  int retval;

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

  if (sensi) {
    retval = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    check_retval(&retval, "CVodeGetSensNumRhsEvals", 1);
    retval = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    check_retval(&retval, "CVodeGetNumRhsEvalsSens", 1);
    retval = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    check_retval(&retval, "CVodeGetSensNumLinSolvSetups", 1);
    if (err_con) {
      retval = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
      check_retval(&retval, "CVodeGetSensNumErrTestFails", 1);
    } else {
      netfS = 0;
    }
    if ((sensi_meth == CV_STAGGERED) || (sensi_meth == CV_STAGGERED1)) {
      retval = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
      check_retval(&retval, "CVodeGetSensNumNonlinSolvIters", 1);
      retval = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
      check_retval(&retval, "CVodeGetSensNumNonlinSolvConvFails", 1);
    } else {
      nniS = 0;
      ncfnS = 0;
    }
  }

  retval = CVodeGetNumLinIters(cvode_mem, &nli);
  check_retval(&retval, "CVodeGetNumLinIters", 1);
  retval = CVodeGetNumLinConvFails(cvode_mem, &ncfl);
  check_retval(&retval, "CVodeGetNumLinConvFails", 1);
  retval = CVodeGetNumPrecEvals(cvode_mem, &npe);
  check_retval(&retval, "CVodeGetNumPrecEvals", 1);
  retval = CVodeGetNumPrecSolves(cvode_mem, &nps);
  check_retval(&retval, "CVodeGetNumPrecSolves", 1);

  printf("\nFinal Statistics\n\n");
  printf("nst     = %5ld\n\n", nst);
  printf("nfe     = %5ld\n",   nfe);
  printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  if(sensi) {
    printf("\n");
    printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }

  printf("\n");
  printf("nli     = %5ld    ncfl     = %5ld\n", nli, ncfl);
  printf("npe     = %5ld    nps      = %5ld\n", npe, nps);

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
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n", funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1); }

  return(0);
}
