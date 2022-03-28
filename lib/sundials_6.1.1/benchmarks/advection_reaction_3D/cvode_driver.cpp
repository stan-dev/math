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
 * ---------------------------------------------------------------------------*/

#include "cvode/cvode.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"
#include "advection_reaction_3D.hpp"
#include "rhs3D.hpp"


/* Setup CVODE and evolve problem in time with BDF method */
int EvolveProblemBDF(N_Vector y, UserData* udata, UserOptions* uopt)
{
  void*              cvode_mem = NULL;   /* empty CVODE memory structure    */
  SUNNonlinearSolver NLS = NULL;         /* empty nonlinear solver structure */
  SUNLinearSolver    LS  = NULL;         /* empty linear solver structure    */

  realtype t, dtout, tout;    /* current/output time data     */
  int      retval;            /* reusable error-checking flag */
  int      iout;              /* output counter               */
  long int nst, netf;         /* step stats                   */
  long int nfi;               /* RHS stats                    */
  long int nni, ncnf;         /* nonlinear solver stats       */
  long int nli, npsol;        /* linear solver stats          */

  /* Additively split methods should not add the advection and reaction terms */
  udata->add_reactions = true;

  /* Create CVode */
  cvode_mem = CVodeCreate(CV_BDF, udata->ctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0, udata->myid)) return 1;

  /* Initialize CVode */
  retval = CVodeInit(cvode_mem, AdvectionReaction, uopt->t0, y);
  if (check_retval((void*)cvode_mem, "CVodeInit", 0, udata->myid)) return 1;

  /* Attach user data */
  retval = CVodeSetUserData(cvode_mem, (void*) udata);
  if (check_retval(&retval, "CVodeSetUserData*", 1, udata->myid)) return 1;

  /* Specify tolerances */
  retval = CVodeSStolerances(cvode_mem, uopt->rtol, uopt->atol);
  if (check_retval(&retval, "CVodeSStolerances", 1, udata->myid)) return 1;

  /* Increase the max number of steps allowed between outputs */
  retval = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1, udata->myid)) return 1;

  /* Create the (non)linear solver */
  if (uopt->nls == "newton")
  {
    /* Create nonlinear solver */
    NLS = SUNNonlinSol_Newton(y, udata->ctx);
    if (check_retval((void *)NLS, "SUNNonlinSol_Newton", 0, udata->myid)) return 1;

    /* Attach nonlinear solver */
    retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
    if (check_retval(&retval, "CVodeSetNonlinearSolver", 1, udata->myid)) return 1;

    /* Create linear solver */
    LS = uopt->precond ? SUNLinSol_SPGMR(y, PREC_LEFT, 0, udata->ctx) : SUNLinSol_SPGMR(y, PREC_NONE, 0, udata->ctx);
    if (check_retval((void *)LS, "SUNLinSol_SPGMR", 0, udata->myid)) return 1;

    /* Attach linear solver */
    retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1, udata->myid)) return 1;

    /* Attach preconditioner */
    retval = CVodeSetPreconditioner(cvode_mem, NULL, PSolve);
    if (check_retval(&retval, "CVodeSetPreconditioner", 1, udata->myid)) return 1;
  }
  else if (uopt->nls == "fixedpoint")
  {
    /* Create nonlinear solver */
    NLS = SUNNonlinSol_FixedPoint(y, uopt->fpaccel, udata->ctx);
    if (check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0, udata->myid)) return 1;

    /* Attach nonlinear solver */
    retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
    if (check_retval(&retval, "CVodeSetNonlinearSolver", 1, udata->myid)) return 1;
  }
  else
  {
    fprintf(stderr, "\nERROR: CV-BDF method is not compatible with the nls option provided\n");
    return 1;
  }

  /* Output initial condition */
  if (uopt->nout > 0)
  {
    if (udata->myid == 0)
    {
      printf("\n          t         ||u||_rms   ||v||_rms   ||w||_rms\n");
      printf("   ----------------------------------------------------\n");
    }
    WriteOutput(uopt->t0, y, udata, uopt);
  }

  /* Integrate to final time */
  t     = uopt->t0;
  dtout = (uopt->tf - uopt->t0);
  if (uopt->nout != 0)
    dtout /= uopt->nout;
  tout  = t + dtout;
  iout  = 0;

  do
  {
    /* Integrate to output time */
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1, udata->myid)) break;

    /* Output state */
    if (uopt->nout > 0) WriteOutput(t, y, udata, uopt);

    /* Update output time */
    tout += dtout;
    tout = (tout > uopt->tf) ? uopt->tf : tout;

    iout++;
  } while (iout < uopt->nout);

  /* Get final statistics */
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1, udata->myid);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfi);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1, udata->myid);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1, udata->myid);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1, udata->myid);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncnf);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1, udata->myid);
  if (uopt->nls == "newton")
  {
    retval = CVodeGetNumLinIters(cvode_mem, &nli);
    check_retval(&retval, "CVodeGetNumLinIters", 1, udata->myid);
    retval = CVodeGetNumPrecSolves(cvode_mem, &npsol);
    check_retval(&retval, "CVodeGetNumPrecSolves", 1, udata->myid);
  }

  /* Print final statistics */
  if (udata->myid == 0)
  {
    printf("\nFinal Solver Statistics (for processor 0):\n");
    printf("   Internal solver steps = %li\n", nst);
    printf("   Total RHS evals: %li\n", nfi + udata->nnlfi);
    printf("   Total number of error test failures = %li\n", netf);
    printf("   Total number of nonlinear solver convergence failures = %li\n",
           ncnf);
    printf("   Total number of nonlinear iterations = %li\n", nni);
    if (uopt->nls == "newton")
    {
      printf("   Total number of linear iterations = %li\n", nli);
      printf("   Total number of preconditioner solves = %li\n", npsol);
    }
  }

  /* Clean up */
  CVodeFree(&cvode_mem);
  if (NLS) SUNNonlinSolFree(NLS);
  if (LS)  SUNLinSolFree(LS);

  /* Return success */
  return(0);
}


/* Setup CVODE and evolve problem in time with Adams method */
int EvolveProblemAdams(N_Vector y, UserData* udata, UserOptions* uopt)
{
  void*              cvode_mem = NULL;   /* empty CVODE memory structure    */
  SUNNonlinearSolver NLS = NULL;         /* empty nonlinear solver structure */

  realtype t, dtout, tout;    /* current/output time data     */
  int      retval;            /* reusable error-checking flag */
  int      iout;              /* output counter               */
  long int nst, netf;         /* step stats                   */
  long int nfi;               /* RHS stats                    */
  long int nni, ncnf;         /* nonlinear solver stats       */

  /* Additively split methods should not add the advection and reaction terms */
  udata->add_reactions = true;

  /* Create CVode */
  cvode_mem = CVodeCreate(CV_ADAMS, udata->ctx);
  if (check_retval((void*)cvode_mem, "CVodeCreate", 0, udata->myid)) return 1;

  /* Initialize CVode */
  retval = CVodeInit(cvode_mem, AdvectionReaction, uopt->t0, y);
  if (check_retval((void*)cvode_mem, "CVodeInit", 0, udata->myid)) return 1;

  /* Attach user data */
  retval = CVodeSetUserData(cvode_mem, (void*) udata);
  if (check_retval(&retval, "CVodeSetUserData*", 1, udata->myid)) return 1;

  /* Specify tolerances */
  retval = CVodeSStolerances(cvode_mem, uopt->rtol, uopt->atol);
  if (check_retval(&retval, "CVodeSStolerances", 1, udata->myid)) return 1;

  /* Increase the max number of steps allowed between outputs */
  retval = CVodeSetMaxNumSteps(cvode_mem, 100000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1, udata->myid)) return 1;

  /* Create nonlinear solver */
  NLS = SUNNonlinSol_FixedPoint(y, uopt->fpaccel, udata->ctx);
  if (check_retval((void *)NLS, "SUNNonlinSol_FixedPoint", 0, udata->myid)) return 1;

  /* Attach nonlinear solver */
  retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
  if (check_retval(&retval, "CVodeSetNonlinearSolver", 1, udata->myid)) return 1;

  /* Output initial condition */
  if (uopt->nout > 0)
  {
    if (udata->myid == 0)
    {
      printf("\n          t         ||u||_rms   ||v||_rms   ||w||_rms\n");
      printf("   ----------------------------------------------------\n");
    }
    WriteOutput(uopt->t0, y, udata, uopt);
  }

  /* Integrate to final time */
  t     = uopt->t0;
  dtout = (uopt->tf - uopt->t0);
  if (uopt->nout != 0)
    dtout /= uopt->nout;
  tout  = t + dtout;
  iout  = 0;

  do
  {
    /* Integrate to output time */
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_retval(&retval, "CVode", 1, udata->myid)) break;

    /* Output state */
    if (uopt->nout > 0) WriteOutput(t, y, udata, uopt);

    /* Update output time */
    tout += dtout;
    tout = (tout > uopt->tf) ? uopt->tf : tout;

    iout++;
  } while (iout < uopt->nout);

  /* Get final statistics */
  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1, udata->myid);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfi);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1, udata->myid);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1, udata->myid);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1, udata->myid);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncnf);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1, udata->myid);

  /* Print final statistics */
  if (udata->myid == 0)
  {
    printf("\nFinal Solver Statistics (for processor 0):\n");
    printf("   Internal solver steps = %li\n", nst);
    printf("   Total RHS evals: %li\n", nfi + udata->nnlfi);
    printf("   Total number of error test failures = %li\n", netf);
    printf("   Total number of nonlinear solver convergence failures = %li\n",
           ncnf);
  }

  /* Clean up */
  CVodeFree(&cvode_mem);
  SUNNonlinSolFree(NLS);

  /* Return success */
  return(0);
}
