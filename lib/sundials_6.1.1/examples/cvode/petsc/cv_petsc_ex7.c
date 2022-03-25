/*-----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *-----------------------------------------------------------------
 * Acknowledgement: This example is based on the PETSc TS ex7.c
 *-----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------*/

static char help[] = "CVODE example based on PETSc TS ex7.c: Nonlinear, time-dependent PDE in 2d.\n";

#include <mpi.h>

/*
   Include "petscdmda.h" so that we can use distributed arrays (DMDAs).
*/
#include <petscdm.h>
#include <petscdmda.h>

/*
   Include "cvode.h" for access to the CVODE BDF integrator. Include
   "sunnonlinsol_petscsnes.h" for access to the SUNNonlinearSolver
   wrapper for PETSc SNES.
*/
#include <cvode/cvode.h>
#include <nvector/nvector_petsc.h>
#include <sunnonlinsol/sunnonlinsol_petscsnes.h>

/*
   User-defined routines in PETSc TS format
*/
extern PetscErrorCode FormFunction(DM,PetscReal,Vec,Vec,void*);
extern PetscErrorCode FormInitialSolution(DM,Vec);
extern PetscErrorCode MySNESMonitor(SNES,PetscInt,PetscReal,PetscViewerAndFormat*);

/*
   User-defined routines in CVODE format
*/

/* f - computes f(t,x); this interfaces FormFunction to the CVODE expected format */
extern int f(PetscReal t, N_Vector x, N_Vector xdot, void *ptr);
extern PetscErrorCode MyCVodeMonitor(long int,PetscReal,Vec,void*);

/* private helper function for checking return value from SUNDIALS calls */
static int check_retval(void *value, const char *funcname, int opt);

int main(int argc,char **argv)
{
  MPI_Comm           comm = PETSC_COMM_WORLD;

  /* SUNDIALS data structures */
  SUNContext         sunctx;
  void*              cvode_mem;        /* integrator memory */
  N_Vector           nvecx;
  SUNNonlinearSolver NLS;
  long int           nsteps = 0;

  /* PETSc data structures */
  SNES                 snes;
  Vec                  x,r;            /* solution, residual vectors */
  Mat                  Jmf;
  PetscErrorCode       ierr;
  DM                   da;
  PetscViewerAndFormat *vf;
  PetscReal            T0, t, tf;

  /* set start and stop time */
  T0 = 0.;
  t  = 0.;
  tf = 0.0005;

  printf("%s\n",help);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
  ierr = SUNContext_Create(&comm, &sunctx);if (ierr) return ierr;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,8,8,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(da);CHKERRQ(ierr);
  ierr = DMSetUp(da);CHKERRQ(ierr);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA; then duplicate for remaining
     vectors that are the same types
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create N_Vector wrapper of petsc vector
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  nvecx = N_VMake_Petsc(x, sunctx);
  if (check_retval((void *)nvecx, "N_VMake_Petsc", 0)) return 1;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create CVODE integrator
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return 1;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver and set its options.

     Set Jacobian matrix data structure and default Jacobian evaluation
     routine. User can override with:
     -snes_mf : matrix-free Newton-Krylov method with no preconditioning
                (unless user explicitly sets preconditioner)
     -snes_mf_operator : form preconditioning matrix as set by the user,
                         but use matrix-free approx for Jacobian-vector
                         products within Newton-Krylov method

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* create SUNNonlinearSolver object which interfaces to SNES */
  NLS = SUNNonlinSol_PetscSNES(nvecx,snes,sunctx); /* This will call SNESSetFunction appropriately */
  if (check_retval((void *)NLS,"SUNNonlinSol_PetscSNES",0)) return 1;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);CHKERRQ(ierr);
  ierr = SNESMonitorSet(snes,(PetscErrorCode (*)(SNES,PetscInt,PetscReal,void*))MySNESMonitor,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);CHKERRQ(ierr);

  /* use matrix free */
  ierr = MatCreateSNESMF(snes,&Jmf);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,Jmf,Jmf,MatMFFDComputeJacobian,0);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions and integrator options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = FormInitialSolution(da,x);CHKERRQ(ierr);
  ierr = CVodeInit(cvode_mem,f,T0,nvecx);
  if (check_retval(&ierr,"CVodeInit",1)) return 1;

  /* provide the DM context as user data so we can access it in the RHS */
  ierr = CVodeSetUserData(cvode_mem,(void *)da);
  if (check_retval(&ierr,"CVodeSetUserData",1)) return 1;

  /* use the PETSc TS default tolerances */
  ierr = CVodeSStolerances(cvode_mem,1e-4,1e-4);
  if (check_retval(&ierr,"CVodeSStolerances",1)) return 1;

  /* set the max order to 1 for Backward Euler */
  ierr = CVodeSetMaxOrd(cvode_mem,1);
  if (check_retval(&ierr,"CVodeSetMaxOrd",1)) return 1;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set the nonlinear solver
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = CVodeSetNonlinearSolver(cvode_mem,NLS);
  if (check_retval(&ierr,"CVodeSetNonlinearSolver",1)) return 1;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  while (t<tf) {
     /* CV_ONE_STEP mode causes CVODE to return after every time step.
        We use it here to demonstrate how to print monitoring information
        at every time step. */
     MyCVodeMonitor(nsteps,t,x,NULL);
     ierr = CVode(cvode_mem,tf,nvecx,&t,CV_ONE_STEP);
     if (check_retval(&ierr,"CVode",1)) break;
     ierr = CVodeGetNumSteps(cvode_mem, &nsteps);
     if (check_retval(&ierr,"CVodeGetNumSteps",1)) break;
  }
  MyCVodeMonitor(nsteps,t,x,NULL);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  CVodeFree(&cvode_mem);
  N_VDestroy(nvecx);
  ierr = SUNNonlinSolFree(NLS);
  if (check_retval(&ierr,"SUNonlinSolFree",0)) return 1;

  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = MatDestroy(&Jmf);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  ierr = SUNContext_Free(&sunctx);
  ierr = PetscFinalize();
  return ierr;
}

/* ------------------------------------------------------------------- */

/* f - computes f(t,x); this is in the CVODE expected format */
int f(PetscReal t, N_Vector x, N_Vector xdot, void *ptr)
{
   PetscErrorCode ierr;
   ierr = FormFunction((DM)ptr,t,N_VGetVector_Petsc(x),N_VGetVector_Petsc(xdot),NULL);
   return ierr;
}

/*
   FormFunction - Evaluates nonlinear function, F(x).

   Input Parameters:
.  DM - the DM context
.  X - input vector
.  ptr - optional user-defined context, as set by SNESSetFunction()

   Output Parameter:
.  F - function vector
 */
PetscErrorCode FormFunction(DM da,PetscReal ftime,Vec X,Vec F,void *ptr)
{
  PetscErrorCode ierr;
  PetscInt       i,j,Mx,My,xs,ys,xm,ym;
  PetscReal      two = 2.0,hx,hy,sx,sy;
  PetscScalar    u,uxx,uyy,**x,**f;
  Vec            localX;

  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);

  hx = 1.0/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = 1.0/(PetscReal)(My-1); sy = 1.0/(hy*hy);

  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  /*
     Get pointers to vector data
  */
  ierr = DMDAVecGetArrayRead(da,localX,&x);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,F,&f);CHKERRQ(ierr);

  /*
     Get local grid boundaries
  */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  /*
     Compute function over the locally owned part of the grid
  */
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        f[j][i] = x[j][i];
        continue;
      }
      u   = x[j][i];
      uxx = (two*u - x[j][i-1] - x[j][i+1])*sx;
      uyy = (two*u - x[j-1][i] - x[j+1][i])*sy;
      /*      f[j][i] = -(uxx + uyy); */
      f[j][i] = -u*(uxx + uyy) - (4.0 - 1.0)*((x[j][i+1] - x[j][i-1])*(x[j][i+1] - x[j][i-1])*.25*sx +
                                              (x[j+1][i] - x[j-1][i])*(x[j+1][i] - x[j-1][i])*.25*sy);
    }
  }

  /*
     Restore vectors
  */
  ierr = DMDAVecRestoreArrayRead(da,localX,&x);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da,F,&f);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localX);CHKERRQ(ierr);
  ierr = PetscLogFlops(11.0*ym*xm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
PetscErrorCode FormInitialSolution(DM da,Vec U)
{
  PetscErrorCode ierr;
  PetscInt       i,j,xs,ys,xm,ym,Mx,My;
  PetscScalar    **u;
  PetscReal      hx,hy,x,y,r;

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);

  hx = 1.0/(PetscReal)(Mx-1);
  hy = 1.0/(PetscReal)(My-1);

  /*
     Get pointers to vector data
  */
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);

  /*
     Get local grid boundaries
  */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);

  /*
     Compute function over the locally owned part of the grid
  */
  for (j=ys; j<ys+ym; j++) {
    y = j*hy;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx;
      r = PetscSqrtReal((x-.5)*(x-.5) + (y-.5)*(y-.5));
      if (r < .125) u[j][i] = PetscExpReal(-30.0*r*r*r);
      else          u[j][i] = 0.0;
    }
  }

  /*
     Restore vectors
  */
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MyCVodeMonitor(long int step,PetscReal ptime,Vec v,void *ctx)
{
  PetscErrorCode ierr;
  PetscReal      norm;

  PetscFunctionBeginUser;
  ierr = VecNorm(v,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"timestep %D time %g norm %g\n",step,(double)ptime,(double)norm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
   MySNESMonitor - illustrate how to set user-defined monitoring routine for SNES.
   Input Parameters:
     snes - the SNES context
     its - iteration number
     fnorm - 2-norm function value (may be estimated)
     ctx - optional user-defined context for private data for the
         monitor routine, as set by SNESMonitorSet()
 */
PetscErrorCode MySNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,PetscViewerAndFormat *vf)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = SNESMonitorDefaultShort(snes,its,fnorm,vf);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void *value, const char *funcname, int opt)
{
  int *errretval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && value == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if retval < 0 */
  else if (opt == 1) {
    errretval = (int *) value;
    if (*errretval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *errretval);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && value == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}
