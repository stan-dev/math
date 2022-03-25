/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the testing routine to check the SUNNonlinearSolver PetscSNES module
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "sundials/sundials_types.h"
#include "nvector/nvector_petsc.h"
#include "sunnonlinsol/sunnonlinsol_petscsnes.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif

#define NEQ   3                /* number of equations        */
#define TOL   RCONST(1.0e-2)   /* nonlinear solver tolerance */
#define MAXIT 10               /* max nonlinear iterations   */

#define ZERO  RCONST(0.0)  /* real 0.0 */
#define HALF  RCONST(0.5)  /* real 0.5 */
#define ONE   RCONST(1.0)  /* real 1.0 */
#define TWO   RCONST(2.0)  /* real 2.0 */
#define THREE RCONST(3.0)  /* real 3.0 */
#define FOUR  RCONST(4.0)  /* real 4.0 */
#define SIX   RCONST(6.0)  /* real 6.0 */

/* approximate solution */
#define Y1 0.785196933062355226
#define Y2 0.496611392944656396
#define Y3 0.369922830745872357

/* Check function return values */
static int check_retval(void *flagvalue, const char *funcname, int opt);

/* Nonlinear residual function */
static int Res(N_Vector y, N_Vector f, void *mem);

/* Jacobian of the nonlinear residual */
int Jac(SNES snes, Vec y, Mat J, Mat Jpre, void *ctx);

/* -----------------------------------------------------------------------------
 * Main testing routine
 * ---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int                retval = 0;
  N_Vector           y, y0, w;
  SUNNonlinearSolver NLS;
  long int           niters;
  SUNContext         sunctx;

  SNES snes;
  Vec X, Y0, Y, W;
  Mat J;

  retval = PetscInitializeNoArguments();CHKERRQ(retval);

  /* create SUNDIALS context */
  retval = SUNContext_Create(NULL, &sunctx);
  if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

  /* create vector */
  VecCreate(PETSC_COMM_WORLD, &X);
  VecSetSizes(X, PETSC_DECIDE, NEQ);
  VecSetFromOptions(X);
  VecDuplicate(X, &Y0);
  VecDuplicate(X, &Y);
  VecDuplicate(X, &W);

  /* create nvector wrappers */
  y0 = N_VMake_Petsc(Y0, sunctx);
  y  = N_VMake_Petsc(Y, sunctx);
  w  = N_VMake_Petsc(W, sunctx);

  /* create Jacobian matrix */
  MatCreate(PETSC_COMM_WORLD, &J);
  MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, NEQ, NEQ);
  MatSetFromOptions(J);
  MatSetUp(J);

  /* set initial guess */
  VecSet(Y0, ZERO);
  VecSet(Y, HALF);
  VecSet(W, ONE);

  /* create SNES context */
  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESSetJacobian(snes, J, J, Jac, NULL);
  /* set the maximum number of nonlinear iterations */
  SNESSetTolerances(snes,
                    TOL,
                    PETSC_DEFAULT,
                    PETSC_DEFAULT,
                    MAXIT,
                    PETSC_DEFAULT);
  SNESSetFromOptions(snes);

  /* create nonlinear solver */
  NLS = SUNNonlinSol_PetscSNES(y, snes, sunctx);
  if (check_retval((void *)NLS, "SUNNonlinSol_PetscSNES", 0)) return(1);

  /* set the nonlinear residual function */
  retval = SUNNonlinSolSetSysFn(NLS, Res);
  if (check_retval(&retval, "SUNNonlinSolSetSysFn", 1)) return(1);

  /* solve the nonlinear system */
  retval = SUNNonlinSolSolve(NLS, y0, y, w, TOL, SUNFALSE, NULL);
  if (check_retval(&retval, "SUNNonlinSolSolve", 1)) return(1);

  /* get the solution */
  realtype yvals[3];
  sunindextype indc[3] = {0, 1, 2};
  VecGetValues(Y, 3, indc, yvals);

  /* print the solution */
  printf("Solution:\n");
  printf("y1 = %"GSYM"\n", yvals[0]);
  printf("y2 = %"GSYM"\n", yvals[1]);
  printf("y3 = %"GSYM"\n", yvals[2]);

  /* print the solution error */
  printf("Solution Error:\n");
  printf("e1 = %"GSYM"\n", yvals[0] - Y1);
  printf("e2 = %"GSYM"\n", yvals[1] - Y2);
  printf("e3 = %"GSYM"\n", yvals[2] - Y3);

  /* get the number of linear iterations */
  retval = SUNNonlinSolGetNumIters(NLS, &niters);
  if (check_retval(&retval, "SUNNonlinSolGetNumIters", 1)) return(1);

  printf("Number of nonlinear iterations: %ld\n",niters);

  /* Free vector, matrix, and nonlinear solver */
  VecDestroy(&X);
  VecDestroy(&Y0);
  VecDestroy(&Y);
  VecDestroy(&W);
  MatDestroy(&J);
  SNESDestroy(&snes);
  N_VDestroy(y);
  N_VDestroy(y0);
  N_VDestroy(w);
  SUNNonlinSolFree(NLS);
  SUNContext_Free(&sunctx);

  /* Print result */
  if (retval) {
    printf("FAIL\n");
  } else {
    printf("SUCCESS\n");
  }

  return(retval);
}


/* -----------------------------------------------------------------------------
 * Nonlinear residual function
 *
 * f1(x,y,z) = x^2 + y^2 + z^2 - 1 = 0
 * f2(x,y,z) = 2x^2 + y^2 - 4z     = 0
 * f3(x,y,z) = 3x^2 - 4y + z^2     = 0
 *
 * ---------------------------------------------------------------------------*/
int Res(N_Vector y, N_Vector f, void *mem)
{
  Vec yvec, fvec;
  realtype vals[3];
  realtype y1, y2, y3;

  yvec = N_VGetVector_Petsc(y);
  fvec = N_VGetVector_Petsc(f);

  /* set vector indices */
  sunindextype indc[3] = {0, 1, 2};

  /* get y vector values */
  VecGetValues(yvec, 3, indc, vals);
  y1 = vals[0]; y2 = vals[1]; y3 = vals[2];

  /* set f vector values */
  vals[0] = y1*y1 + y2*y2 + y3*y3 - ONE;
  vals[1] = TWO * y1*y1 + y2*y2 - FOUR * y3;
  vals[2] = THREE * (y1*y1) - FOUR * y2 + y3*y3;
  VecSetValues(fvec, 3, indc, vals, INSERT_VALUES);

  /* assemble the f vector */
  VecAssemblyBegin(fvec);
  VecAssemblyEnd(fvec);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Jacobian of the nonlinear residual function
 *
 *            ( 2x  2y  2z )
 * J(x,y,z) = ( 4x  2y  -4 )
 *            ( 6x  -4  2z )
 *
 * ---------------------------------------------------------------------------*/
int Jac(SNES snes, Vec y, Mat J, Mat Jpre, void *ctx)
{
  realtype y1, y2, y3;
  realtype yvals[3];

  /* set vector indices */
  sunindextype indc[3] = { 0, 1, 2 };

  /* get y vector values */
  VecGetValues(y, 3, indc, yvals);
  y1 = yvals[0]; y2 = yvals[1]; y3 = yvals[2];

  /* set the Jacobian values */
  realtype jvals[3][3] =  { { TWO*y1,  TWO*y2,  TWO*y3 },
                            { FOUR*y1, TWO*y2, -FOUR   },
                            { SIX*y1, -FOUR,    TWO*y3 } };
  MatSetValues(J, 3, indc, 3, indc, &jvals[0][0], INSERT_VALUES);

  /* assemble the matrix */
  if (J != Jpre) {
    MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY);
  }
  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Check function return value
 *   opt == 0 check if returned NULL pointer
 *   opt == 1 check if returned a non-zero value
 * ---------------------------------------------------------------------------*/
static int check_retval(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if the function returned a NULL pointer -- no memory allocated */
  if (opt == 0) {
    if (flagvalue == NULL) {
      fprintf(stderr, "\nERROR: %s() failed -- returned NULL\n\n", funcname);
      return(1);
    } else {
      return(0);
    }
  }

  /* Check if the function returned an non-zero value -- internal failure */
  if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag != 0) {
      fprintf(stderr, "\nERROR: %s() failed -- returned %d\n\n", funcname, *errflag);
      return(1);
    } else {
      return(0);
    }
  }

  /* if we make it here then opt was not 0 or 1 */
  fprintf(stderr, "\nERROR: check_retval failed -- Invalid opt value\n\n");
  return(1);
}
