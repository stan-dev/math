/*-----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *-----------------------------------------------------------------
 * Acknowledgement: This example is based on the PETSc TS ex25.c.
 *-----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *----------------------------------------------------------------
   u_t - alpha u_xx = A + u^2 v - (B+1) u
   v_t - alpha v_xx = B u - u^2 v
   0 < x < 1;
   A = 1, B = 3, alpha = 1/50

   Initial conditions:
   u(x,0) = 1 + sin(2 pi x)
   v(x,0) = 3

   Boundary conditions:
   u(0,t) = u(1,t) = 1
   v(0,t) = v(1,t) = 3
 -----------------------------------------------------------------*/

static const char help[] = "ARKode example based on PETSc TS ex25.c.\nTime-dependent Brusselator reaction-diffusion PDE in 1d. Demonstrates IMEX methods.\n";

#include <petscdm.h>
#include <petscdmda.h>

#include <arkode/arkode_arkstep.h>
#include <nvector/nvector_petsc.h>
#include <sunnonlinsol/sunnonlinsol_petscsnes.h>

typedef struct {
  PetscScalar u,v;
} Field;

typedef struct _User *User;
struct _User {
  PetscReal A,B;                /* Reaction coefficients */
  PetscReal alpha;              /* Diffusion coefficient */
  PetscReal uleft,uright;       /* Dirichlet boundary conditions */
  PetscReal vleft,vright;       /* Dirichlet boundary conditions */
  DM da;                        /* PETSc DM for the problem */
  void* arkode_mem;             /* ARKode memory structure */
};

static PetscErrorCode FormRHSFunction(PetscReal,Vec,Vec,void*);
static PetscErrorCode FormIFunction(PetscReal,Vec,Vec,Vec,void*);
static PetscErrorCode FormIJacobian(SNES,Vec,Mat,Mat,void*);
static PetscErrorCode FormInitialSolution(Vec,void*);

/* User-supplied Functions called by ARKStep */
static int f_I(PetscReal,N_Vector,N_Vector,void*);
static int f_E(PetscReal,N_Vector,N_Vector,void*);

/* Private function to check function return values */
static int check_retval(void *retvalvalue, const char *funcname, int opt);

int main(int argc, char **argv)
{
  long int           steps=0;

  /* SUNDIALS data structures */
  void              *arkode_mem; /* integrator memory */
  N_Vector           nvecx;      /* SUNDIALS N_Vector wrapper of X */
  SUNNonlinearSolver NLS;        /* SUNDIALS nonlinear solver */

  /* PETSc data structures */
  SNES              snes;       /* nonlinear solver */
  Vec               X;          /* solution, residual vectors */
  Field             *xdata;     /* underlying x data */
  Mat               J;          /* Jacobian matrix */
  PetscInt          mx;
  PetscErrorCode    ierr;
  PetscReal         T0,t,tstop,dtout,ftime,hx,dt;
  PetscReal         rtol,atol;
  struct _User      user;       /* user-defined work context */


  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /* Solution start and end time */
  T0    = 0.0;
  ftime = 10.0;
  tstop = T0;
  t     = T0;
  dtout = 1.0;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,11,2,2,NULL,&user.da);CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.da);CHKERRQ(ierr);
  ierr = DMSetUp(user.da);CHKERRQ(ierr);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA;
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = DMCreateGlobalVector(user.da,&X);CHKERRQ(ierr);
  nvecx = N_VMake_Petsc(X);
  if (check_retval((void *)nvecx,"N_VMake_Petsc",0)) return 1;

  /* Initialize user application context */
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Advection-reaction options","");
  {
    rtol        = 1e-4;
    atol        = 1e-4;
    user.A      = 1;
    user.B      = 3;
    user.alpha  = 0.02;
    user.uleft  = 1;
    user.uright = 1;
    user.vleft  = 3;
    user.vright = 3;
    ierr        = PetscOptionsReal("-rtol","Relative tolerance","",rtol,&rtol,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsReal("-atol","Absolute tolerance","",atol,&atol,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsReal("-A","Reaction rate","",user.A,&user.A,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsReal("-B","Reaction rate","",user.B,&user.B,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsReal("-alpha","Diffusion coefficient","",user.alpha,&user.alpha,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsReal("-uleft","Dirichlet boundary condition","",user.uleft,&user.uleft,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsReal("-uright","Dirichlet boundary condition","",user.uright,&user.uright,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsReal("-vleft","Dirichlet boundary condition","",user.vleft,&user.vleft,NULL);CHKERRQ(ierr);
    ierr        = PetscOptionsReal("-vright","Dirichlet boundary condition","",user.vright,&user.vright,NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = FormInitialSolution(X,&user);CHKERRQ(ierr);
  ierr = VecGetSize(X,&mx);CHKERRQ(ierr);

  hx = 1.0/(PetscReal)(mx/2-1);
  dt = 0.4 * PetscSqr(hx) / user.alpha; /* Diffusive stability limit */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create ARKStep time stepper
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* Call ARKStepCreate to initialize the ARK timestepper module and
     specify the right-hand side function in y'=f(t,y),the inital time
     T0,and the initial dependent variable vector y. */
  arkode_mem = ARKStepCreate(f_E,f_I,T0,nvecx);
  if (check_retval((void *)arkode_mem,"ARKStepCreate",0)) return 1;

  /* Store the arkode mem in the user data so we can access it in the Jacobian routine */
  user.arkode_mem = arkode_mem;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create the nonlinear solver
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* Create SUNNonlinearSolver object which interfaces to SNES */
  NLS = SUNNonlinSol_PetscSNES(nvecx,snes); /* this will call SNESSetFunction appropriately */
  if (check_retval((void *)NLS,"SUNNonlinSol_PetscSNES",0)) return 1;

  /* Set the Jacobian routine */
  ierr = DMSetMatType(user.da,MATAIJ);CHKERRQ(ierr);
  ierr = DMCreateMatrix(user.da,&J);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,FormIJacobian,&user);CHKERRQ(ierr);

  /* Allow SNES/KSP/PC runtime options */
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set ARKStep options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = ARKStepSStolerances(arkode_mem,rtol,atol);
  if (check_retval(&ierr,"ARKStepSStolerances",1)) return 1;

  ierr = ARKStepSetOrder(arkode_mem,3);
  if (check_retval(&ierr,"ARKStepSetOrder",1)) return 1;

  ierr = ARKStepSetUserData(arkode_mem,(void *) &user);
  if (check_retval(&ierr,"ARKStepSetUserData",1)) return 1;

  ierr = ARKStepSetNonlinearSolver(arkode_mem,NLS);
  if (check_retval(&ierr,"ARKStepSetNonlinearSolver",1)) return 1;

  ierr = ARKStepSetAdaptivityMethod(arkode_mem,2,1,0,NULL);
  if (check_retval(&ierr,"ARKStepSetAdaptivity",1)) return 1;

  ierr = ARKStepSetInitStep(arkode_mem,dt);
  if (check_retval(&ierr,"ARKStepSetInitStep",1)) return 1;


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Perform the integration
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  printf("%s\n",help);

  /* Extract underlying data of x for printing */
  ierr = DMDAVecGetArrayRead(user.da,X,&xdata);CHKERRQ(ierr);

  /* Print out the solution every dt */
  while(tstop<=ftime) {
    printf("%ld TS dt %.6e time %.6f\n",steps,dt,t);

    /* Advance time */
    tstop+=dtout;
    ierr = ARKStepSetStopTime(arkode_mem,tstop);
    if (check_retval(&ierr,"ARKStepSetStopTime",1)) return 1;

    /* Evolve solution in time */
    ierr = ARKStepEvolve(arkode_mem,ftime,nvecx,&t,ARK_NORMAL);
    if (check_retval(&ierr,"ARKStepEvolve",1)) return 1;

    /* Get statistics */
    ierr = ARKStepGetCurrentStep(arkode_mem,&dt);
    if (check_retval(&ierr,"ARKStepGetCurrntStep",1)) return 1;
    ierr = ARKStepGetNumSteps(arkode_mem,&steps);
    if (check_retval(&ierr,"ARKStepGetNumSteps",1)) return 1;
  }

  printf("Converged at time %g after %ld steps.\n",t,steps);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* Free SUNDIALS data structures */
  N_VDestroy(nvecx);           /* Free x nvector         */
  SUNNonlinSolFree(NLS);       /* Free nonlinear solver  */
  ARKStepFree(&arkode_mem);    /* Free integrator memory */

  /* Free petsc data structures */
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  User provided functions in ARKStep format
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Implicit component of RHS */
static int f_I(PetscReal t,N_Vector x,N_Vector xdot,void *ptr)
{
  PetscErrorCode ierr;
  ierr = FormIFunction(t,N_VGetVector_Petsc(x),NULL,N_VGetVector_Petsc(xdot),ptr);
  return ierr;
}

/* Explicit component of RHS */
static int f_E(PetscReal t,N_Vector x,N_Vector xdot,void *ptr)
{
  PetscErrorCode ierr;
  ierr = FormRHSFunction(t,N_VGetVector_Petsc(x),N_VGetVector_Petsc(xdot),ptr);
  return ierr;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  User provided functions in Petsc TS format (minus the TS context)
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

static PetscErrorCode FormIFunction(PetscReal t,Vec X,Vec Xdot,Vec F,void *ptr)
{
  User           user = (User)ptr;
  DMDALocalInfo  info;
  PetscInt       i;
  Field          *x,*f;
  PetscReal      hx;
  Vec            Xloc;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da,&info);CHKERRQ(ierr);
  hx   = 1.0/(PetscReal)(info.mx-1);

  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(), DMGlobalToLocalEnd().
     By placing code between these two statements,computations can be
     done while messages are in transition.
  */
  ierr = DMGetLocalVector(user->da,&Xloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(user->da,X,INSERT_VALUES,Xloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(user->da,X,INSERT_VALUES,Xloc);CHKERRQ(ierr);

  /* Get pointers to vector data */
  ierr = DMDAVecGetArrayRead(user->da,Xloc,&x);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da,F,&f);CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */
  for (i=info.xs; i<info.xs+info.xm; i++) {
    if (i == 0) {
      f[i].u = 0.;
      f[i].v = 0.;
    } else if (i == info.mx-1) {
      f[i].u = 0.;
      f[i].v = 0.;
    } else {
      f[i].u = user->alpha * (x[i-1].u - 2.*x[i].u + x[i+1].u) / (hx*hx);
      f[i].v = user->alpha * (x[i-1].v - 2.*x[i].v + x[i+1].v) / (hx*hx);
    }
  }

  /* Restore vectors */
  ierr = DMDAVecRestoreArrayRead(user->da,Xloc,&x);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da,F,&f);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(user->da,&Xloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode FormRHSFunction(PetscReal t,Vec X,Vec F,void *ptr)
{
  User           user = (User)ptr;
  DMDALocalInfo  info;
  PetscInt       i;
  Field          *x,*f;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da,&info);CHKERRQ(ierr);

  /* Get pointers to vector data */
  ierr = DMDAVecGetArrayRead(user->da,X,&x);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da,F,&f);CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */
  for (i=info.xs; i<info.xs+info.xm; i++) {
    PetscScalar u = x[i].u,v = x[i].v;
    if (i == 0) {
      f[i].u = 0.;
      f[i].v = 0.;
    } else if (i == info.mx-1) {
      f[i].u = 0.;
      f[i].v = 0.;
    } else {
      f[i].u = user->A + u*u*v - (user->B+1)*u;
      f[i].v = user->B*u - u*u*v;
    }
  }

  /* Restore vectors */
  ierr = DMDAVecRestoreArrayRead(user->da,X,&x);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da,F,&f);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*
  IJacobian - Compute J = I - gamma * df_i/dx
*/
PetscErrorCode FormIJacobian(SNES snes,Vec X,Mat J,Mat Jpre,void *ptr)
{
  User           user = (User)ptr;
  PetscErrorCode ierr;
  DMDALocalInfo  info;
  PetscInt       i;
  PetscReal      hx,gamma,a;
  Field          *x;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da,&info);CHKERRQ(ierr);
  hx   = 1.0/(PetscReal)(info.mx-1);

  /* Get current gamma value from ARKode */
  ierr = ARKStepGetCurrentGamma(user->arkode_mem,&gamma);

  /* Get pointers to vector data */
  ierr = DMDAVecGetArrayRead(user->da,X,&x);CHKERRQ(ierr);

  /* Set shortcut value */
  a = user->alpha/hx/hx * gamma;

  /* Compute function over the locally owned part of the grid */
  for (i=info.xs; i<info.xs+info.xm; i++) {
    if (i == 0 || i == info.mx-1) {
      const PetscInt    row        = i,col = i;
      const PetscScalar vals[2][2] = {{1.,0},{0,1.}};
      ierr = MatSetValuesBlocked(Jpre,1,&row,1,&col,&vals[0][0],INSERT_VALUES);CHKERRQ(ierr);
    } else {
      const PetscInt    row           = i,col[] = {i-1,i,i+1};
      const PetscScalar dxxL          = -a,dxx0 = 2.*a,dxxR = -a;
      const PetscScalar vals[2][3][2] = {{{dxxL,0},{1.+dxx0,0},{dxxR,0}},
                                         {{0,dxxL},{0,1.+dxx0},{0,dxxR}}};
      ierr = MatSetValuesBlocked(Jpre,1,&row,3,col,&vals[0][0][0],INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  /* Restore vectors */
  ierr = DMDAVecRestoreArrayRead(user->da,X,&x);CHKERRQ(ierr);

  /* Assemble matrix */
  ierr = MatAssemblyBegin(Jpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Jpre,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (J != Jpre) {
    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  /* Add I */
  ierr = MatShift(Jpre,1.0);

  PetscFunctionReturn(0);
}

PetscErrorCode FormInitialSolution(Vec X,void *ctx)
{
  User           user = (User)ctx;
  PetscInt       i;
  DMDALocalInfo  info;
  Field          *x;
  PetscReal      hx;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMDAGetLocalInfo(user->da,&info);CHKERRQ(ierr);
  hx   = 1.0/(PetscReal)(info.mx-1);

  /* Get pointers to vector data */
  ierr = DMDAVecGetArray(user->da,X,&x);CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */
  for (i=info.xs; i<info.xs+info.xm; i++) {
    PetscReal xi = i*hx;
    x[i].u = user->uleft*(1.-xi) + user->uright*xi + PetscSinReal(2.*PETSC_PI*xi);
    x[i].v = user->vleft*(1.-xi) + user->vright*xi;
  }
  ierr = DMDAVecRestoreArray(user->da,X,&x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a retval so check if
             retval >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer
*/
static int check_retval(void *value,const char *funcname,int opt)
{
  int *errretval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && value == NULL) {
    fprintf(stderr,"\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  /* Check if retval < 0 */
  else if (opt == 1) {
    errretval = (int *) value;
    if (*errretval < 0) {
      fprintf(stderr,"\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname,*errretval);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && value == NULL) {
    fprintf(stderr,"\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return 1; }

  return 0;
}
