/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * Based on PETSc TS example 15 and idaHeat2D_petsc_spgmr.c by
 * Slaven Peles @ LLNL and Daniel Reynolds @ SMU.
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, PETSc vector, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver SUNLinSol_SPGMR.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform MX x MY
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = MX * MY. Here MX = MY = 10.
 *
 * The system is actually implemented on submeshes, processor by
 * processor, with an MXSUB by MYSUB mesh on each of NPEX * NPEY
 * processors.
 *
 * The system is solved with IDA using the SUNNonlinearSolver
 * interface to PETSc SNES. By default it uses the SNES "newtonls"
 * solver with  the KSP GMRES linear solver. The '-pre' runtime
 * option can be used to toggle on/off user-supplied routines for
 * preconditioning. The preconditioner uses the diagonal elements
 * of the Jacobian only. The '-jac' runtime option can be used to
 * toggle on/off a user-supplied Jacobian routine.
 *
 * Local error testing on the boundary values is suppressed.
 * Output is taken at t = 0, .01, .02, .04, ..., 10.24.
 *
 * The example also uses the PETSc data management functions to generate
 * the mesh and handle MPI communication.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ida/ida.h>
#include <nvector/nvector_petsc.h>
#include <sunnonlinsol/sunnonlinsol_petscsnes.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <mpi.h>
#include <petscdm.h>
#include <petscdmda.h>

#if defined(SUNDIALS_INT64_T)
#define DSYM "ld"
#else
#define DSYM "d"
#endif

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define NOUT         11             /* Number of output times */

#define NPEX         2              /* No. PEs in x direction of PE array */
#define NPEY         2              /* No. PEs in y direction of PE array */
                                    /* Total no. PEs = NPEX*NPEY */
#define MXSUB        5              /* No. x points per subgrid */
#define MYSUB        5              /* No. y points per subgrid */

#define MX           (NPEX*MXSUB)   /* MX = number of x mesh points */
#define MY           (NPEY*MYSUB)   /* MY = number of y mesh points */
                                    /* Spatial mesh is MX by MY */

typedef struct {
  N_Vector pp;      /* vector of diagonal preconditioner elements */
  DM       da;      /* PETSc data management object */
  void*    ida_mem; /* IDA integrator memory */
} *UserData;

/* User-supplied residual function */

int resHeat(realtype tt, N_Vector uu, N_Vector up, N_Vector rr, void *user_data);
int jacHeat(SNES snes, Vec x, Mat Jpre, Mat J, void *user_data);

/* User-supplied preconditioner routines */

int PsetupHeat(PC pc);

int PsolveHeat(PC pc, Vec x, Vec y);

/* Private function to check function return values */

static int SetInitialProfile(N_Vector uu, N_Vector up, N_Vector id,
                             N_Vector res, void *user_data);

static void PrintHeader(sunindextype Neq, realtype rtol, realtype atol);

static void PrintOutput(int id, void *ida_mem, realtype t, N_Vector uu, SNES snes);

static void PrintFinalStats(void *ida_mem, SNES snes);

static int check_retval(void *returnvalue, const char *funcname, int opt, int id);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  MPI_Comm comm;
  int iout, thispe, retval, npes;
  realtype rtol, atol, t0, t1, tout, tret;
  sunindextype Neq;
  PetscBool pre, jac;
  /* declare SUNDIALS data structures */
  void *ida_mem;
  SUNNonlinearSolver NLS;
  N_Vector uu, up, constraints, id, res;
  UserData data;
  /* declare PETSc data structures */
  PetscErrorCode ierr;
  SNES snes;
  KSP ksp;
  PC pc;
  Mat J;
  Vec uvec;

  /* Get processor number and total number of pe's. */
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &thispe);

  if (npes != NPEX*NPEY) {
    if (thispe == 0)
      fprintf(stderr,
              "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n",
              npes,NPEX*NPEY);
    MPI_Finalize();
    return(1);
  }

  /* Set global vector length Neq. */
  Neq = MX * MY;

  /* Initialize PETSc */
  ierr = PetscInitialize(&argc, &argv, (char*)0, NULL);

  /* Initialize user application context */
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Heat2D options", "");
  {
    jac = PETSC_FALSE;
    pre = PETSC_FALSE;
    ierr = PetscOptionsBool("-jac", "Utilize user-supplied Jacobian", "", jac, &jac, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-pre", "Utilize user-supplied preconditioner", "", pre, &pre, NULL); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  /* Allocate and initialize the data structure and N-vectors. */
  data = (UserData) malloc(sizeof *data);
  if(check_retval((void *)data, "malloc", 2, thispe))
    MPI_Abort(comm, 1);
  data->pp = NULL;
  data->da = NULL;

  /* Create the object that will manage the communication of 2D data */
  ierr = DMDACreate2d(comm,
                      DM_BOUNDARY_NONE,  /* NONE, PERIODIC, GHOSTED */
                      DM_BOUNDARY_NONE,
                      DMDA_STENCIL_STAR, /* STAR, BOX */
                      MX,
                      MY,
                      NPEX,              /* Set numbers or use PETSC_DECIDE */
                      NPEY,
                      1,                 /* degrees of freedom per node */
                      1,                 /* stencil width */
                      NULL,              /* number of nodes per cell in x */
                      NULL,              /* number of nodes per cell in y */
                      &(data->da));
  CHKERRQ(ierr);
  DMSetFromOptions(data->da);
  DMSetUp(data->da);

  /* Create PETSc vector */
  ierr = DMCreateGlobalVector(data->da, &uvec);
  CHKERRQ(ierr);

  /* Make N_Vector wrapper for uvec */
  uu = N_VMake_Petsc(uvec);
  if(check_retval((void *)uu, "N_VNew_Petsc", 0, thispe))
    MPI_Abort(comm, 1);

  up = N_VClone(uu);
  if(check_retval((void *)up, "N_VClone", 0, thispe))
    MPI_Abort(comm, 1);

  res = N_VClone(uu);
  if(check_retval((void *)res, "N_VClone", 0, thispe))
    MPI_Abort(comm, 1);

  constraints = N_VClone(uu);
  if(check_retval((void *)constraints, "N_VClone", 0, thispe))
    MPI_Abort(comm, 1);

  id = N_VClone(uu);
  if(check_retval((void *)id, "N_VClone", 0, thispe))
    MPI_Abort(comm, 1);

  /* An N-vector to hold preconditioner. */
  data->pp = N_VClone(uu);
  if(check_retval((void *)data->pp, "N_VClone", 0, thispe))
    MPI_Abort(comm, 1);

  /* Initialize the uu, up, id, and res profiles. */
  SetInitialProfile(uu, up, id, res, data);

  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst(ONE, constraints);
  t0 = ZERO; t1 = RCONST(0.01);

  /* Scalar relative and absolute tolerance. */
  rtol = ZERO;
  atol = RCONST(1.0e-4);

  /*
   * Call IDACreate and IDAInit, set integration tolerances, then set optional inputs
   */

  ida_mem = IDACreate();
  if(check_retval((void *)ida_mem, "IDACreate", 0, thispe)) MPI_Abort(comm, 1);
  data->ida_mem = ida_mem;

  retval = IDAInit(ida_mem, resHeat, t0, uu, up);
  if(check_retval(&retval, "IDAInit", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASStolerances(ida_mem, rtol, atol);
  if(check_retval(&retval, "IDASStolerances", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetUserData(ida_mem, data);
  if(check_retval(&retval, "IDASetUserData", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetSuppressAlg(ida_mem, SUNTRUE);
  if(check_retval(&retval, "IDASetSuppressAlg", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetId(ida_mem, id);
  if(check_retval(&retval, "IDASetId", 1, thispe)) MPI_Abort(comm, 1);

  retval = IDASetConstraints(ida_mem, constraints);
  if(check_retval(&retval, "IDASetConstraints", 1, thispe)) MPI_Abort(comm, 1);
  N_VDestroy_Petsc(constraints);

  /*
   * Create SNES context, then wrap the SNES context in a SUNNonlinsol_PetscSNES
   * object and set options. Finally, tell IDA to use it.
   */

  ierr = SNESCreate(comm, &snes); CHKERRQ(ierr);

  /* Wrap the SNES context in a SUNNonlinsol_PetscSNES object */
  NLS = SUNNonlinSol_PetscSNES(uu, snes);  /* This will call SNESSetFunction appropriately */
  if(check_retval((void*)NLS, "SUNNonlinsol_PetscSNES", 0, thispe)) MPI_Abort(comm, 1);

  ierr = SNESSetDM(snes, data->da); CHKERRQ(ierr);

  if (jac) {
    ierr = DMSetMatType(data->da, MATAIJ);
    ierr = DMCreateMatrix(data->da, &J); CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, J, J, jacHeat, data); CHKERRQ(ierr);
  } else {
    /* Tell SNES to use the matrix-free finite difference scheme for both
      the preconditioner matrix and the Jacobian. */
    ierr = MatCreateSNESMF(snes, &J); CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, J, J, MatMFFDComputeJacobian, NULL); CHKERRQ(ierr);
  }

  /* Setup the preconditioner if it is turned on */
  if (pre) {
    SNESGetKSP(snes, &ksp);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCSHELL);
    PCShellSetSetUp(pc, PsetupHeat);
    PCShellSetApply(pc, PsolveHeat);
    PCShellSetContext(pc, (void *) data);
  }

  /* Set SNES/KSP/PC runtime options, these will override the defaults set */
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  /* Tell ida to use the SUNNonlinsol_PetscSNES nonlinear solver */
  retval = IDASetNonlinearSolver(ida_mem, NLS);
  if(check_retval(&retval, "IDASetNonlinearSolver", 0, thispe)) MPI_Abort(comm, 1);

  /*
   * Solve the problem, printing output at the desired points in time
   */

  /* Print output heading (on processor 0 only) and intial solution  */
  if (thispe == 0) PrintHeader(Neq, rtol, atol);
  PrintOutput(thispe, ida_mem, t0, uu, snes);

  /* Loop over tout, call IDASolve, print output. */
  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {
    retval = IDASolve(ida_mem, tout, &tret, uu, up, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1, thispe)) MPI_Abort(comm, 1);
    PrintOutput(thispe, ida_mem, tret, uu, snes);
  }

  /* Print remaining counters. */
  if (thispe == 0) PrintFinalStats(ida_mem, snes);

  /*
   * Free memory
   */

  IDAFree(&ida_mem);

  N_VDestroy_Petsc(id);
  N_VDestroy_Petsc(res);
  N_VDestroy_Petsc(up);
  N_VDestroy_Petsc(uu);

  N_VDestroy_Petsc(data->pp);
  DMDestroy(&data->da);
  free(data);

  VecDestroy(&uvec);
  MatDestroy(&J);
  SNESDestroy(&snes);

  /*
   * Finalize and exit
   */

  ierr = PetscFinalize();
  MPI_Finalize();

  return(0);

}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * resHeat: heat equation system residual function
 * This uses 5-point central differencing on the interior points, and
 * includes algebraic equations for the boundary values.
 * So for each interior point, the residual component has the form
 *    res_i = u'_i - (central difference)_i
 * while for each boundary point, it is res_i = u_i.
 *
 */

int resHeat(realtype tt, N_Vector uu, N_Vector up, N_Vector rr,
            void *user_data)
{
  PetscErrorCode ierr;
  UserData       data = (UserData) user_data;
  DM             da   = (DM) data->da;
  PetscInt       i, j, Mx, My, xs, ys, xm, ym;
  PetscReal      hx, hy, sx, sy;
  PetscScalar    u, uxx, uyy, **uarray, **f, **udot;
  Vec localU;
  Vec U    = N_VGetVector_Petsc(uu);
  Vec Udot = N_VGetVector_Petsc(up);
  Vec F    = N_VGetVector_Petsc(rr);

  PetscFunctionBeginUser;
  ierr = DMGetLocalVector(da, &localU);
  CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,
                     PETSC_IGNORE,
                     &Mx,          /* Get global grid x size */
                     &My,          /* Get global grid y size */
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE);

  hx = 1.0/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = 1.0/(PetscReal)(My-1); sy = 1.0/(hy*hy);

  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  ierr = DMGlobalToLocalBegin(da, U, INSERT_VALUES, localU);
  CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da, U, INSERT_VALUES, localU);
  CHKERRQ(ierr);

  /* Get pointers to vector data */
  ierr = DMDAVecGetArrayRead(da, localU, &uarray);
  CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, F, &f);
  CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, Udot, &udot);
  CHKERRQ(ierr);

  /* Get local grid boundaries */
  ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      /* Boundary conditions */
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        f[j][i] = uarray[j][i]; /* F = U */
      } else { /* Interior */
        u = uarray[j][i];
        /* 5-point stencil */
        uxx = (-2.0*u + uarray[j][i-1] + uarray[j][i+1]);
        uyy = (-2.0*u + uarray[j-1][i] + uarray[j+1][i]);
        f[j][i] = udot[j][i] - (uxx*sx + uyy*sy);
      }
    }
  }

  /* Restore vectors */
  ierr = DMDAVecRestoreArrayRead(da, localU, &uarray);
  CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, F, &f);
  CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, Udot, &udot);
  CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da, &localU);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*
 * jacHeat: Heat equation system Jacobian matrix.
 *
 * The optional user-supplied functions jacHeat provides Jacobian
 * matrix
 *        J = dF/du + cj*dF/du'
 * where the DAE system is F(t,u,u') = 0.
 *
 */
int jacHeat(SNES snes, Vec x, Mat Jpre, Mat J, void *user_data)
{
  PetscErrorCode ierr;
  PetscInt       i, j, Mx , My, xs, ys, xm, ym, nc;
  UserData       data = (UserData) user_data;
  DM             da   = (DM) data->da;
  MatStencil     col[5], row;
  PetscScalar    vals[5], hx, hy, sx, sy, cj;

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,
                     PETSC_IGNORE,
                     &Mx,
                     &My,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE);
  ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL); CHKERRQ(ierr);

  hx = 1.0/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = 1.0/(PetscReal)(My-1); sy = 1.0/(hy*hy);

  ierr = IDAGetCurrentCj(data->ida_mem, &cj);

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      nc    = 0;
      row.j = j; row.i = i;
      if ((i == 0 || i == Mx-1 || j == 0 || j == My-1)) { /* Dirichlet BC */
        col[nc].j = j; col[nc].i = i; vals[nc++] = 1.0;
      } else {   /* Interior */
        col[nc].j = j-1; col[nc].i = i;   vals[nc++] = -sy;
        col[nc].j = j;   col[nc].i = i-1; vals[nc++] = -sx;
        col[nc].j = j;   col[nc].i = i;   vals[nc++] = 2.0*(sx + sy) + cj;
        col[nc].j = j;   col[nc].i = i+1; vals[nc++] = -sx;
        col[nc].j = j+1; col[nc].i = i;   vals[nc++] = -sy;
      }
      ierr = MatSetValuesStencil(Jpre, 1, &row, nc, col, vals, INSERT_VALUES); CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  if (J != Jpre) {
    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*
 * PsetupHeat: setup for diagonal preconditioner for heatsk.
 *
 * The optional user-supplied functions PsetupHeat and
 * PsolveHeat together must define the left preconditoner
 * matrix P approximating the system Jacobian matrix
 *                   J = dF/du + cj*dF/du'
 * (where the DAE system is F(t,u,u') = 0), and solve the linear
 * systems P z = r.   This is done in this case by keeping only
 * the diagonal elements of the J matrix above, storing them as
 * inverses in a vector pp, when computed in PsetupHeat, for
 * subsequent use in PsolveHeat.
 *
 * In this instance, only cj and data (user data structure, with
 * pp etc.) are used from the PsetupHeat argument list.
 *
 */

int PsetupHeat(PC pc)
{
  PetscErrorCode ierr;
  PetscInt    i, j, Mx, My, xs, ys, xm, ym;
  PetscReal   hx,hy,sx,sy,c_j;
  PetscScalar pelinv;
  PetscScalar **ppv;
  UserData data;

  /* Get the user data out of the PC object */
  PCShellGetContext(pc, (void **) &data);

  DM da = (DM) data->da;
  Vec ppvec = N_VGetVector_Petsc((data->pp));
  IDAGetCurrentCj(data->ida_mem, &c_j);

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,
                     PETSC_IGNORE,
                     &Mx,          /* Get global grid x size */
                     &My,          /* Get global grid y size */
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE);
  ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
  CHKERRQ(ierr);

  hx = 1.0/(PetscReal)(Mx-1); sx = 1.0/(hx*hx);
  hy = 1.0/(PetscReal)(My-1); sy = 1.0/(hy*hy);

  pelinv = ONE/(2.0*(sx + sy) + c_j);

  ierr = DMDAVecGetArray(da, ppvec, &ppv);
  CHKERRQ(ierr);

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        ppv[j][i] = ONE;    /* Boundary */
      } else {
        ppv[j][i] = pelinv; /* Interior */
      }
    }
  }

  ierr = DMDAVecRestoreArray(da, ppvec, &ppv);
  CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/*
 * PsolveHeat: solve preconditioner linear system.
 * This routine multiplies the input vector rvec by the vector pp
 * containing the inverse diagonal Jacobian elements (previously
 * computed in PsetupHeat), returning the result in zvec.
 */

int PsolveHeat(PC pc, Vec r, Vec z)
{
  UserData data;

  PetscFunctionBeginUser;

  /* Get the user data out of the PC object */
  PCShellGetContext(pc, (void **) &data);
  /* Multiply r by pp and store in z */
  VecPointwiseMult(z, N_VGetVector_Petsc(data->pp), r);

  PetscFunctionReturn(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */


/*
 * SetInitialProfile sets the initial values for the problem.
 */

static int SetInitialProfile(N_Vector uu, N_Vector up, N_Vector id,
                             N_Vector res, void *user_data)
{
  UserData       data = (UserData) user_data;
  DM             da   = data->da;
  PetscErrorCode ierr;
  PetscInt       i, j, xs, ys, xm, ym, Mx, My;
  PetscScalar    **u;
  PetscScalar    **iddat;
  PetscReal      hx, hy, x, y;
  Vec U     = N_VGetVector_Petsc(uu);
  Vec idvec = N_VGetVector_Petsc(id);

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(da,
                     PETSC_IGNORE,
                     &Mx,
                     &My,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE,
                     PETSC_IGNORE);

  hx = 1.0/(PetscReal)(Mx-1);
  hy = 1.0/(PetscReal)(My-1);

  /* Get pointers to vector data */
  ierr = DMDAVecGetArray(da, U, &u);
  CHKERRQ(ierr);

  /* Get pointers to differentiable variable IDs */
  ierr = DMDAVecGetArray(da, idvec, &iddat);
  CHKERRQ(ierr);

  /* Get local grid boundaries */
  ierr = DMDAGetCorners(da,&xs,&ys,NULL,&xm,&ym,NULL);
  CHKERRQ(ierr);

  /* Compute function over the locally owned part of the grid */
  for (j=ys; j<ys+ym; j++) {
    y = j*hy;
    for (i=xs; i<xs+xm; i++) {
      x = i*hx;
      u[j][i] = 16.0 * x*(1.0 - x) * y*(1.0 - y);
      if (i == 0 || j == 0 || i == Mx-1 || j == My-1) {
        iddat[j][i] = 0.0; /* algebraic variables on the boundary */
      } else {
        iddat[j][i] = 1.0; /* differential variables in the interior */
      }
    }
  }

  /* Restore vectors */
  ierr = DMDAVecRestoreArray(da, U, &u);
  CHKERRQ(ierr);

   /* Restore vectors */
  ierr = DMDAVecRestoreArray(da, idvec, &iddat);
  CHKERRQ(ierr);

  /* Initialize up. */
  N_VConst(ZERO, up);    /* Initially set up = 0. */

  /* resHeat sets res to negative of ODE RHS values at interior points. */
  resHeat(ZERO, uu, up, res, data);

  /* Copy -res into up to get correct initial up values. */
  N_VScale(-ONE, res, up);

  PetscFunctionReturn(0);
}


/*
 * Print first lines of output and table heading
 */

static void PrintHeader(sunindextype Neq, realtype rtol, realtype atol)
{
  printf("\nidaHeat2D_kry_petscsnes: Heat equation, parallel example problem for IDA\n");
  printf("                         Discretized heat equation on 2D unit square.\n");
  printf("                         Zero boundary conditions, polynomial initial conditions.\n");
  printf("                         Mesh dimensions: %d x %d", MX, MY);
  printf("\tTotal system size: %ld\n", (long int) Neq);
  printf("                         Subgrid dimensions: %d x %d", MXSUB, MYSUB);
  printf("\tProcessor array: %d x %d\n", NPEX, NPEY);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("                         Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("                         Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("                         Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("\nConstraints set to force all solution components >= 0.\n");
  printf("SUPPRESSALG = SUNTRUE to suppress local error testing on all boundary components.\n");
  printf("Linear solver: PETSc KSP GMRES.\n");
  printf("Preconditioner: user-supplied, diagonal elements only.\n");
  printf("This example uses the PETSc SNES nonlinear solver interface.\n");
  printf("Use the SNES command line options to override defaults.\n");

  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time     umax       k  nst  nni  nli   nre     h      \n");
  printf("--------------------------------------------------------\n");
}

/*
 * PrintOutput: print max norm of solution and current solver statistics
 */

static void PrintOutput(int id, void *ida_mem, realtype t, N_Vector uu, SNES snes)
{
  realtype hused, umax;
  long int nst, nni, njve, nre, nreLS, npe, nps;
  int kused;
  PetscInt nli;
  int retval;
  KSP ksp;

  nst = nni = njve = nre = nreLS = npe = nps = -99;
  nli = -99;
  umax = N_VMaxNorm(uu);
  SNESGetKSP(snes, &ksp);

  if (id == 0) {

    retval = IDAGetLastOrder(ida_mem, &kused);
    check_retval(&retval, "IDAGetLastOrder", 1, id);
    retval = IDAGetNumSteps(ida_mem, &nst);
    check_retval(&retval, "IDAGetNumSteps", 1, id);
    retval = IDAGetLastStep(ida_mem, &hused);
    check_retval(&retval, "IDAGetLastStep", 1, id);
    retval = IDAGetNumResEvals(ida_mem, &nre);
    check_retval(&retval, "IDAGetNumResEvals", 1, id);
    retval = IDAGetNumNonlinSolvIters(ida_mem, &nni);
    check_retval(&retval, "IDAGetNumNonlinSolvIters", 1, id);
    retval = KSPGetTotalIterations(ksp, &nli); CHKERRV(retval);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3"DSYM"  %4ld  %9.2Le \n",
           t, umax, kused, nst, nni, nli, nre, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3"DSYM"  %4ld   %9.2e  \n",
           t, umax, kused, nst, nni, nli, nre, hused);
#else
    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3"DSYM"  %4ld   %9.2e  \n",
           t, umax, kused, nst, nni, nli, nre, hused);
#endif

  }
}

/*
 * Print some final integrator statistics
 */

static void PrintFinalStats(void *ida_mem, SNES snes)
{
  long int netf, ncfn;
  PetscInt ncfl;

  netf = ncfn = ncfl = -99;

  IDAGetNumErrTestFails(ida_mem, &netf);
  IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  SNESGetLinearSolveFailures(snes, &ncfl);

  printf("\nError test failures            = %ld\n", netf);
  printf("Nonlinear convergence failures = %ld\n", ncfn);
  printf("Linear convergence failures    = %"DSYM"\n", ncfl);
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
