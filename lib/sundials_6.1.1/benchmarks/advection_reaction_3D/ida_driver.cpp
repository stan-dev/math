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

#include "ida/ida.h"
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"
#include "advection_reaction_3D.hpp"
#include "rhs3D.hpp"


/* Initial condition function */
int SetICDot(N_Vector y, N_Vector yp, UserData* udata)
{
  int retval;

  retval = AdvectionReaction(0, y, yp, (void*)udata);
  if (check_retval(&retval, "AdvectionReaction", 1, udata->myid)) return 1;

  /* Return success */
  return(0);
}


/* Setup IDA and evolve problem in time with BDF method */
int EvolveDAEProblem(N_Vector y, UserData* udata, UserOptions* uopt)
{
  void*              ida_mem = NULL;  /* empty IDA memory structure       */
  SUNNonlinearSolver NLS = NULL;      /* empty nonlinear solver structure */
  SUNLinearSolver    LS  = NULL;      /* empty linear solver structure    */
  N_Vector           yp  = NULL;      /* empty vector structure           */

  realtype t, dtout, tout;    /* current/output time data     */
  int      retval;            /* reusable error-checking flag */
  int      iout;              /* output counter               */
  long int nst, netf;         /* step stats                   */
  long int nfi;               /* RHS stats                    */
  long int nni, ncnf;         /* nonlinear solver stats       */
  long int nli, npsol;        /* linear solver stats          */

  /* Additively split methods should not add the advection and reaction terms */
  udata->add_reactions = true;

  /* Create ydot' vector */
  yp = N_VClone(y);
  if (check_retval((void*)yp, "N_VClone", 0, udata->myid)) return 1;

  /* Create IDA */
  ida_mem = IDACreate(udata->ctx);
  if (check_retval((void*)ida_mem, "IDACreate", 0, udata->myid)) return 1;

  /* Initialize IDA */
  retval = IDAInit(ida_mem, AdvectionReactionResidual, uopt->t0, y, yp);
  if (check_retval(&retval, "IDAInit", 1, udata->myid)) return 1;

  /* Attach user data */
  retval = IDASetUserData(ida_mem, (void*) udata);
  if (check_retval(&retval, "IDASetUserData*", 1, udata->myid)) return 1;

  /* Specify tolerances */
  retval = IDASStolerances(ida_mem, uopt->rtol, uopt->atol);
  if (check_retval(&retval, "IDASStolerances", 1, udata->myid)) return 1;

  /* Increase the max number of steps allowed between outputs */
  retval = IDASetMaxNumSteps(ida_mem, 100000);
  if (check_retval(&retval, "IDASetMaxNumSteps", 1, udata->myid)) return 1;

  /* Increase the max number of ETF allowed between outputs */
  retval = IDASetMaxErrTestFails(ida_mem, 25);
  if (check_retval(&retval, "IDASetMaxErrTestFails", 1, udata->myid)) return 1;

  /* Create the (non)linear solver */
  if (uopt->nls == "newton")
  {
    /* Create nonlinear solver */
    NLS = SUNNonlinSol_Newton(y, udata->ctx);
    if (check_retval((void *)NLS, "SUNNonlinSol_Newton", 0, udata->myid)) return 1;

    /* Attach nonlinear solver */
    retval = IDASetNonlinearSolver(ida_mem, NLS);
    if (check_retval(&retval, "IDASetNonlinearSolver", 1, udata->myid)) return 1;

    /* Create linear solver */
    LS = uopt->precond ? SUNLinSol_SPGMR(y, PREC_LEFT, 0, udata->ctx) : SUNLinSol_SPGMR(y, PREC_NONE, 0, udata->ctx);
    if (check_retval((void *)LS, "SUNLinSol_SPGMR", 0, udata->myid)) return 1;

    /* Attach linear solver */
    retval = IDASetLinearSolver(ida_mem, LS, NULL);
    if (check_retval(&retval, "IDASetLinearSolver", 1, udata->myid)) return 1;

    // /* Attach preconditioner */
    retval = IDASetPreconditioner(ida_mem, NULL, PSolveRes);
    if (check_retval(&retval, "IDASetPreconditioner", 1, udata->myid)) return 1;
  }
  else
  {
    fprintf(stderr, "\nERROR: IDA method is not compatible with the nls option provided\n");
    return 1;
  }

  /* Set ydot' initial condition */
  retval = SetICDot(y, yp, udata);
  if (check_retval(&retval, "SetICDot", 1, udata->myid)) return 1;

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
    retval = IDASolve(ida_mem, tout, &t, y, yp, IDA_NORMAL);
    if (check_retval(&retval, "IDA", 1, udata->myid)) break;

    /* Output state */
    if(uopt->nout > 0) WriteOutput(t, y, udata, uopt);

    /* Update output time */
    tout += dtout;
    tout = (tout > uopt->tf) ? uopt->tf : tout;

    iout++;
  } while (iout < uopt->nout);

  /* Get final statistics */
  retval = IDAGetNumSteps(ida_mem, &nst);
  check_retval(&retval, "IDAGetNumSteps", 1, udata->myid);
  retval = IDAGetNumResEvals(ida_mem, &nfi);
  check_retval(&retval, "IDAGetNumResEvals", 1, udata->myid);
  retval = IDAGetNumErrTestFails(ida_mem, &netf);
  check_retval(&retval, "IDAGetNumErrTestFails", 1, udata->myid);
  retval = IDAGetNumNonlinSolvIters(ida_mem, &nni);
  check_retval(&retval, "IDAGetNumNonlinSolvIters", 1, udata->myid);
  retval = IDAGetNumNonlinSolvConvFails(ida_mem, &ncnf);
  check_retval(&retval, "IDAGetNumNonlinSolvConvFails", 1, udata->myid);
  if (uopt->nls == "newton")
  {
    retval = IDAGetNumLinIters(ida_mem, &nli);
    check_retval(&retval, "IDAGetNumLinIters", 1, udata->myid);
    retval = IDAGetNumPrecSolves(ida_mem, &npsol);
    check_retval(&retval, "IDAGetNumPrecSolves", 1, udata->myid);
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
  IDAFree(&ida_mem);
  if (yp) N_VDestroy(yp);
  if (NLS) SUNNonlinSolFree(NLS);
  if (LS)  SUNLinSolFree(LS);

  /* Return success */
  return(0);
}
