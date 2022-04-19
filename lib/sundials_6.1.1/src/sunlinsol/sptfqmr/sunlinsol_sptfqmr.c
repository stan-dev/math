/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on sundials_sptfqmr.c code, written by Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the SPTFQMR implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_sptfqmr.h>
#include <sundials/sundials_math.h>

#include "sundials_debug.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * SPTFQMR solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define SPTFQMR_CONTENT(S)  ( (SUNLinearSolverContent_SPTFQMR)(S->content) )
#define LASTFLAG(S)         ( SPTFQMR_CONTENT(S)->last_flag )

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SPTFQMR linear solver
 */

SUNLinearSolver SUNLinSol_SPTFQMR(N_Vector y, int pretype, int maxl, SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_SPTFQMR content;

  /* check for legal pretype and maxl values; if illegal use defaults */
  if ((pretype != SUN_PREC_NONE)  && (pretype != SUN_PREC_LEFT) &&
      (pretype != SUN_PREC_RIGHT) && (pretype != SUN_PREC_BOTH))
    pretype = SUN_PREC_NONE;
  if (maxl <= 0)
    maxl = SUNSPTFQMR_MAXL_DEFAULT;

  /* check that the supplied N_Vector supports all requisite operations */
  if ( (y->ops->nvclone == NULL) || (y->ops->nvdestroy == NULL) ||
       (y->ops->nvlinearsum == NULL) || (y->ops->nvconst == NULL) ||
       (y->ops->nvprod == NULL) || (y->ops->nvdiv == NULL) ||
       (y->ops->nvscale == NULL) || (y->ops->nvdotprod == NULL) )
    return(NULL);

  /* Create linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype           = SUNLinSolGetType_SPTFQMR;
  S->ops->getid             = SUNLinSolGetID_SPTFQMR;
  S->ops->setatimes         = SUNLinSolSetATimes_SPTFQMR;
  S->ops->setpreconditioner = SUNLinSolSetPreconditioner_SPTFQMR;
  S->ops->setscalingvectors = SUNLinSolSetScalingVectors_SPTFQMR;
  S->ops->setzeroguess      = SUNLinSolSetZeroGuess_SPTFQMR;
  S->ops->initialize        = SUNLinSolInitialize_SPTFQMR;
  S->ops->setup             = SUNLinSolSetup_SPTFQMR;
  S->ops->solve             = SUNLinSolSolve_SPTFQMR;
  S->ops->numiters          = SUNLinSolNumIters_SPTFQMR;
  S->ops->resnorm           = SUNLinSolResNorm_SPTFQMR;
  S->ops->resid             = SUNLinSolResid_SPTFQMR;
  S->ops->lastflag          = SUNLinSolLastFlag_SPTFQMR;
  S->ops->space             = SUNLinSolSpace_SPTFQMR;
  S->ops->free              = SUNLinSolFree_SPTFQMR;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SPTFQMR) malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag   = 0;
  content->maxl        = maxl;
  content->pretype     = pretype;
  content->zeroguess   = SUNFALSE;
  content->numiters    = 0;
  content->resnorm     = ZERO;
  content->r_star      = NULL;
  content->q           = NULL;
  content->d           = NULL;
  content->v           = NULL;
  content->p           = NULL;
  content->r           = NULL;
  content->u           = NULL;
  content->vtemp1      = NULL;
  content->vtemp2      = NULL;
  content->vtemp3      = NULL;
  content->s1          = NULL;
  content->s2          = NULL;
  content->ATimes      = NULL;
  content->ATData      = NULL;
  content->Psetup      = NULL;
  content->Psolve      = NULL;
  content->PData       = NULL;
  content->print_level = 0;
  content->info_file   = stdout;

  /* Allocate content */
  content->r_star = N_VClone(y);
  if (content->r_star == NULL) { SUNLinSolFree(S); return(NULL); }

  content->q = N_VClone(y);
  if (content->q == NULL) { SUNLinSolFree(S); return(NULL); }

  content->d = N_VClone(y);
  if (content->d == NULL) { SUNLinSolFree(S); return(NULL); }

  content->v = N_VClone(y);
  if (content->v == NULL) { SUNLinSolFree(S); return(NULL); }

  content->p = N_VClone(y);
  if (content->p == NULL) { SUNLinSolFree(S); return(NULL); }

  content->r = N_VCloneVectorArray(2, y);
  if (content->r == NULL) { SUNLinSolFree(S); return(NULL); }

  content->u = N_VClone(y);
  if (content->u == NULL) { SUNLinSolFree(S); return(NULL); }

  content->vtemp1 = N_VClone(y);
  if (content->vtemp1 == NULL) { SUNLinSolFree(S); return(NULL); }

  content->vtemp2 = N_VClone(y);
  if (content->vtemp2 == NULL) { SUNLinSolFree(S); return(NULL); }

  content->vtemp3 = N_VClone(y);
  if (content->vtemp3 == NULL) { SUNLinSolFree(S); return(NULL); }

  return(S);
}


/* ----------------------------------------------------------------------------
 * Function to set the type of preconditioning for SPTFQMR to use
 */

int SUNLinSol_SPTFQMRSetPrecType(SUNLinearSolver S, int pretype)
{
  /* Check for legal pretype */
  if ((pretype != SUN_PREC_NONE)  && (pretype != SUN_PREC_LEFT) &&
      (pretype != SUN_PREC_RIGHT) && (pretype != SUN_PREC_BOTH)) {
    return(SUNLS_ILL_INPUT);
  }

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SUNLS_MEM_NULL);

  /* Set pretype */
  SPTFQMR_CONTENT(S)->pretype = pretype;
  return(SUNLS_SUCCESS);
}


/* ----------------------------------------------------------------------------
 * Function to set the maximum number of iterations for SPTFQMR to use
 */

int SUNLinSol_SPTFQMRSetMaxl(SUNLinearSolver S, int maxl)
{
  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SUNLS_MEM_NULL);

  /* Check for legal pretype */
  if (maxl <= 0)
    maxl = SUNSPTFQMR_MAXL_DEFAULT;

  /* Set pretype */
  SPTFQMR_CONTENT(S)->maxl = maxl;
  return(SUNLS_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SPTFQMR(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_ITERATIVE);
}


SUNLinearSolver_ID SUNLinSolGetID_SPTFQMR(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_SPTFQMR);
}


int SUNLinSolInitialize_SPTFQMR(SUNLinearSolver S)
{
  SUNLinearSolverContent_SPTFQMR content;

  /* set shortcut to SPTFQMR memory structure */
  if (S == NULL) return(SUNLS_MEM_NULL);
  content = SPTFQMR_CONTENT(S);

  /* ensure valid options */
  if (content->maxl <= 0)
    content->maxl = SUNSPTFQMR_MAXL_DEFAULT;

  if (content->ATimes == NULL) {
    LASTFLAG(S) = SUNLS_ATIMES_NULL;
    return(LASTFLAG(S));
  }

  if ( (content->pretype != SUN_PREC_LEFT) &&
       (content->pretype != SUN_PREC_RIGHT) &&
       (content->pretype != SUN_PREC_BOTH) )
    content->pretype = SUN_PREC_NONE;

  if ((content->pretype != SUN_PREC_NONE) && (content->Psolve == NULL)) {
    LASTFLAG(S) = SUNLS_PSOLVE_NULL;
    return(LASTFLAG(S));
  }

  /* no additional memory to allocate */

  /* return with success */
  content->last_flag = SUNLS_SUCCESS;
  return(SUNLS_SUCCESS);
}


int SUNLinSolSetATimes_SPTFQMR(SUNLinearSolver S, void* ATData,
                               SUNATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  if (S == NULL) return(SUNLS_MEM_NULL);
  SPTFQMR_CONTENT(S)->ATimes = ATimes;
  SPTFQMR_CONTENT(S)->ATData = ATData;
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSetPreconditioner_SPTFQMR(SUNLinearSolver S, void* PData,
                                       SUNPSetupFn Psetup, SUNPSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  if (S == NULL) return(SUNLS_MEM_NULL);
  SPTFQMR_CONTENT(S)->Psetup = Psetup;
  SPTFQMR_CONTENT(S)->Psolve = Psolve;
  SPTFQMR_CONTENT(S)->PData = PData;
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSetScalingVectors_SPTFQMR(SUNLinearSolver S,
                                       N_Vector s1,
                                       N_Vector s2)
{
  /* set N_Vector pointers to integrator-supplied scaling vectors,
     and return with success */
  if (S == NULL) return(SUNLS_MEM_NULL);
  SPTFQMR_CONTENT(S)->s1 = s1;
  SPTFQMR_CONTENT(S)->s2 = s2;
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSetZeroGuess_SPTFQMR(SUNLinearSolver S, booleantype onoff)
{
  /* set flag indicating a zero initial guess */
  if (S == NULL) return(SUNLS_MEM_NULL);
  SPTFQMR_CONTENT(S)->zeroguess = onoff;
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSetup_SPTFQMR(SUNLinearSolver S, SUNMatrix A)
{
  int ier;
  SUNPSetupFn Psetup;
  void* PData;

  /* Set shortcuts to SPTFQMR memory structures */
  if (S == NULL) return(SUNLS_MEM_NULL);
  Psetup = SPTFQMR_CONTENT(S)->Psetup;
  PData = SPTFQMR_CONTENT(S)->PData;

  /* no solver-specific setup is required, but if user-supplied
     Psetup routine exists, call that here */
  if (Psetup != NULL) {
    ier = Psetup(PData);
    if (ier != 0) {
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_PSET_FAIL_UNREC : SUNLS_PSET_FAIL_REC;
      return(LASTFLAG(S));
    }
  }

  /* return with success */
  return(SUNLS_SUCCESS);
}


int SUNLinSolSolve_SPTFQMR(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                           N_Vector b, realtype delta)
{
  /* local data and shortcut variables */
  realtype alpha, tau, eta, beta, c, sigma, v_bar, omega;
  realtype rho[2];
  realtype r_init_norm, r_curr_norm;
  realtype temp_val;
  booleantype preOnLeft, preOnRight, scale_x, scale_b, converged, b_ok;
  booleantype *zeroguess;
  int n, m, ier, l_max;
  void *A_data, *P_data;
  SUNATimesFn atimes;
  SUNPSolveFn psolve;
  realtype *res_norm;
  int *nli;
  N_Vector sx, sb, r_star, q, d, v, p, *r, u, vtemp1, vtemp2, vtemp3;

  /* local variables for fused vector operations */
  realtype cv[3];
  N_Vector Xv[3];

  /* Make local shorcuts to solver variables. */
  if (S == NULL) return(SUNLS_MEM_NULL);
  l_max        = SPTFQMR_CONTENT(S)->maxl;
  r_star       = SPTFQMR_CONTENT(S)->r_star;
  q            = SPTFQMR_CONTENT(S)->q;
  d            = SPTFQMR_CONTENT(S)->d;
  v            = SPTFQMR_CONTENT(S)->v;
  p            = SPTFQMR_CONTENT(S)->p;
  r            = SPTFQMR_CONTENT(S)->r;
  u            = SPTFQMR_CONTENT(S)->u;
  vtemp1       = SPTFQMR_CONTENT(S)->vtemp1;
  vtemp2       = SPTFQMR_CONTENT(S)->vtemp2;
  vtemp3       = SPTFQMR_CONTENT(S)->vtemp3;
  sb           = SPTFQMR_CONTENT(S)->s1;
  sx           = SPTFQMR_CONTENT(S)->s2;
  A_data       = SPTFQMR_CONTENT(S)->ATData;
  P_data       = SPTFQMR_CONTENT(S)->PData;
  atimes       = SPTFQMR_CONTENT(S)->ATimes;
  psolve       = SPTFQMR_CONTENT(S)->Psolve;
  zeroguess    = &(SPTFQMR_CONTENT(S)->zeroguess);
  nli          = &(SPTFQMR_CONTENT(S)->numiters);
  res_norm     = &(SPTFQMR_CONTENT(S)->resnorm);

  /* Initialize counters and convergence flag */
  temp_val = r_curr_norm = -ONE;
  *nli = 0;
  converged = SUNFALSE;
  b_ok = SUNFALSE;

  /* set booleantype flags for internal solver options */
  preOnLeft  = ( (SPTFQMR_CONTENT(S)->pretype == SUN_PREC_LEFT) ||
                 (SPTFQMR_CONTENT(S)->pretype == SUN_PREC_BOTH) );
  preOnRight = ( (SPTFQMR_CONTENT(S)->pretype == SUN_PREC_RIGHT) ||
                 (SPTFQMR_CONTENT(S)->pretype == SUN_PREC_BOTH) );
  scale_x = (sx != NULL);
  scale_b = (sb != NULL);

  /* Check for unsupported use case */
  if (preOnRight && !(*zeroguess)) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return(SUNLS_ILL_INPUT);
  }

#ifdef SUNDIALS_BUILD_WITH_MONITORING
  if (SPTFQMR_CONTENT(S)->print_level && SPTFQMR_CONTENT(S)->info_file)
    STAN_SUNDIALS_FPRINTF(SPTFQMR_CONTENT(S)->info_file, "SUNLINSOL_SPTFQMR:\n");
#endif

  /* Check if Atimes function has been set */
  if (atimes == NULL) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_ATIMES_NULL;
    return(LASTFLAG(S));
  }

  /* If preconditioning, check if psolve has been set */
  if ((preOnLeft || preOnRight) && psolve == NULL) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_PSOLVE_NULL;
    return(LASTFLAG(S));
  }

  /* Set r_star to initial (unscaled) residual r_star = r_0 = b - A*x_0 */
  /* NOTE: if x == 0 then just set residual to b and continue */
  if (*zeroguess) {
    N_VScale(ONE, b, r_star);
  } else {
    ier = atimes(A_data, x, r_star);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }
    N_VLinearSum(ONE, b, -ONE, r_star, r_star);
  }

  /* Apply left preconditioner and b-scaling to r_star (or really just r_0) */
  if (preOnLeft) {
    ier = psolve(P_data, r_star, vtemp1, delta, SUN_PREC_LEFT);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
      return(LASTFLAG(S));
    }
  }
  else N_VScale(ONE, r_star, vtemp1);
  if (scale_b) N_VProd(sb, vtemp1, r_star);
  else N_VScale(ONE, vtemp1, r_star);

  /* Initialize rho[0] */
  /* NOTE: initialized here to reduce number of computations - avoid need
           to compute r_star^T*r_star twice, and avoid needlessly squaring
           values */
  rho[0] = N_VDotProd(r_star, r_star);

  /* Compute norm of initial residual (r_0) to see if we really need
     to do anything */
  *res_norm = r_init_norm = SUNRsqrt(rho[0]);

#ifdef SUNDIALS_BUILD_WITH_MONITORING
  /* print initial residual */
  if (SPTFQMR_CONTENT(S)->print_level && SPTFQMR_CONTENT(S)->info_file)
  {
    STAN_SUNDIALS_FPRINTF(SPTFQMR_CONTENT(S)->info_file,
            SUNLS_MSG_RESIDUAL,
            (long int) 0, *res_norm);
  }
#endif

  if (r_init_norm <= delta) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_SUCCESS;
    return(LASTFLAG(S));
  }

  /* Set v = A*r_0 (preconditioned and scaled) */
  if (scale_x) N_VDiv(r_star, sx, vtemp1);
  else N_VScale(ONE, r_star, vtemp1);
  if (preOnRight) {
    N_VScale(ONE, vtemp1, v);
    ier = psolve(P_data, v, vtemp1, delta, SUN_PREC_RIGHT);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
      return(LASTFLAG(S));
    }
  }
  ier = atimes(A_data, vtemp1, v);
  if (ier != 0) {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = (ier < 0) ?
      SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
    return(LASTFLAG(S));
  }
  if (preOnLeft) {
    ier = psolve(P_data, v, vtemp1, delta, SUN_PREC_LEFT);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
      return(LASTFLAG(S));
    }
  }
  else N_VScale(ONE, v, vtemp1);
  if (scale_b) N_VProd(sb, vtemp1, v);
  else N_VScale(ONE, vtemp1, v);

  /* Initialize remaining variables */
  N_VScale(ONE, r_star, r[0]);
  N_VScale(ONE, r_star, u);
  N_VScale(ONE, r_star, p);
  N_VConst(ZERO, d);

  /* Set x = sx x if non-zero guess */
  if (scale_x && !(*zeroguess)) N_VProd(sx, x, x);

  tau = r_init_norm;
  v_bar = eta = ZERO;

  /* START outer loop */
  for (n = 0; n < l_max; ++n) {

    /* Increment linear iteration counter */
    (*nli)++;

    /* sigma = r_star^T*v */
    sigma = N_VDotProd(r_star, v);

    /* alpha = rho[0]/sigma */
    alpha = rho[0]/sigma;

    /* q = u-alpha*v */
    N_VLinearSum(ONE, u, -alpha, v, q);

    /* r[1] = r[0]-alpha*A*(u+q) */
    N_VLinearSum(ONE, u, ONE, q, r[1]);
    if (scale_x) N_VDiv(r[1], sx, r[1]);
    if (preOnRight) {
      N_VScale(ONE, r[1], vtemp1);
      ier = psolve(P_data, vtemp1, r[1], delta, SUN_PREC_RIGHT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
        return(LASTFLAG(S));
      }
    }
    ier = atimes(A_data, r[1], vtemp1);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }
    if (preOnLeft) {
      ier = psolve(P_data, vtemp1, r[1], delta, SUN_PREC_LEFT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
        return(LASTFLAG(S));
      }
    }
    else N_VScale(ONE, vtemp1, r[1]);
    if (scale_b) N_VProd(sb, r[1], vtemp1);
    else N_VScale(ONE, r[1], vtemp1);
    N_VLinearSum(ONE, r[0], -alpha, vtemp1, r[1]);

    /* START inner loop */
    for (m = 0; m < 2; ++m) {

      /* d = [*]+(v_bar^2*eta/alpha)*d */
      /* NOTES:
       *   (1) [*] = u if m == 0, and q if m == 1
       *   (2) using temp_val reduces the number of required computations
       *       if the inner loop is executed twice
       */
      if (m == 0) {
        temp_val = SUNRsqrt(N_VDotProd(r[1], r[1]));
        omega = SUNRsqrt(SUNRsqrt(N_VDotProd(r[0], r[0]))*temp_val);
        N_VLinearSum(ONE, u, SUNSQR(v_bar)*eta/alpha, d, d);
      }
      else {
        omega = temp_val;
        N_VLinearSum(ONE, q, SUNSQR(v_bar)*eta/alpha, d, d);
      }

      /* v_bar = omega/tau */
      v_bar = omega/tau;

      /* c = (1+v_bar^2)^(-1/2) */
      c = ONE / SUNRsqrt(ONE+SUNSQR(v_bar));

      /* tau = tau*v_bar*c */
      tau = tau*v_bar*c;

      /* eta = c^2*alpha */
      eta = SUNSQR(c)*alpha;

      /* x = x+eta*d */
      if (n == 0 && m == 0 && *zeroguess)
        N_VScale(eta, d, x);
      else
        N_VLinearSum(ONE, x, eta, d, x);

      /* Check for convergence... */
      /* NOTE: just use approximation to norm of residual, if possible */
      *res_norm = r_curr_norm = tau*SUNRsqrt(m+1);

#ifdef SUNDIALS_BUILD_WITH_MONITORING
      /* print current iteration number and the residual */
      if (SPTFQMR_CONTENT(S)->print_level && SPTFQMR_CONTENT(S)->info_file)
      {
        STAN_SUNDIALS_FPRINTF(SPTFQMR_CONTENT(S)->info_file,
                SUNLS_MSG_RESIDUAL,
                (long int) *nli, *res_norm);
      }
#endif

      /* Exit inner loop if iteration has converged based upon approximation
         to norm of current residual */
      if (r_curr_norm <= delta) {
        converged = SUNTRUE;
        break;
      }

      /* Decide if actual norm of residual vector should be computed */
      /* NOTES:
       *   (1) if r_curr_norm > delta, then check if actual residual norm
       *       is OK (recall we first compute an approximation)
       *   (2) if r_curr_norm >= r_init_norm and m == 1 and n == l_max, then
       *       compute actual residual norm to see if the iteration can be
       *       saved
       *   (3) the scaled and preconditioned right-hand side of the given
       *       linear system (denoted by b) is only computed once, and the
       *       result is stored in vtemp3 so it can be reused - reduces the
       *       number of psovles if using left preconditioning
       */
      if ((r_curr_norm > delta) ||
          (r_curr_norm >= r_init_norm && m == 1 && n == l_max)) {

        /* Compute norm of residual ||b-A*x||_2 (preconditioned and scaled) */
        if (scale_x) N_VDiv(x, sx, vtemp1);
        else N_VScale(ONE, x, vtemp1);
        if (preOnRight) {
          ier = psolve(P_data, vtemp1, vtemp2, delta, SUN_PREC_RIGHT);
          if (ier != 0) {
            *zeroguess  = SUNFALSE;
            LASTFLAG(S) = (ier < 0) ?
              SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_UNREC;
            return(LASTFLAG(S));
          }
          N_VScale(ONE, vtemp2, vtemp1);
        }
        ier = atimes(A_data, vtemp1, vtemp2);
        if (ier != 0) {
          *zeroguess  = SUNFALSE;
          LASTFLAG(S) = (ier < 0) ?
            SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
          return(LASTFLAG(S));
        }
        if (preOnLeft) {
          ier = psolve(P_data, vtemp2, vtemp1, delta, SUN_PREC_LEFT);
          if (ier != 0) {
            *zeroguess  = SUNFALSE;
            LASTFLAG(S) = (ier < 0) ?
              SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
            return(LASTFLAG(S));
          }
        }
        else N_VScale(ONE, vtemp2, vtemp1);
        if (scale_b) N_VProd(sb, vtemp1, vtemp2);
        else N_VScale(ONE, vtemp1, vtemp2);
        /* Only precondition and scale b once (result saved for reuse) */
        if (!b_ok) {
          b_ok = SUNTRUE;
          if (preOnLeft) {
            ier = psolve(P_data, b, vtemp3, delta, SUN_PREC_LEFT);
            if (ier != 0) {
              *zeroguess  = SUNFALSE;
              LASTFLAG(S) = (ier < 0) ?
                SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
              return(LASTFLAG(S));
            }
          }
          else N_VScale(ONE, b, vtemp3);
          if (scale_b) N_VProd(sb, vtemp3, vtemp3);
        }
        N_VLinearSum(ONE, vtemp3, -ONE, vtemp2, vtemp1);
        *res_norm = r_curr_norm = SUNRsqrt(N_VDotProd(vtemp1, vtemp1));

        /* Exit inner loop if inequality condition is satisfied
           (meaning exit if we have converged) */
        if (r_curr_norm <= delta) {
          converged = SUNTRUE;
          break;
        }

      }

    }  /* END inner loop */

    /* If converged, then exit outer loop as well */
    if (converged == SUNTRUE) break;

    /* rho[1] = r_star^T*r_[1] */
    rho[1] = N_VDotProd(r_star, r[1]);

    /* beta = rho[1]/rho[0] */
    beta = rho[1]/rho[0];

    /* u = r[1]+beta*q */
    N_VLinearSum(ONE, r[1], beta, q, u);

    /* p = u+beta*(q+beta*p) = beta*beta*p + beta*q + u */
    cv[0] = SUNSQR(beta);
    Xv[0] = p;

    cv[1] = beta;
    Xv[1] = q;

    cv[2] = ONE;
    Xv[2] = u;

    ier = N_VLinearCombination(3, cv, Xv, p);
    if (ier != SUNLS_SUCCESS) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = SUNLS_VECTOROP_ERR;
      return(SUNLS_VECTOROP_ERR);
    }

    /* v = A*p */
    if (scale_x) N_VDiv(p, sx, vtemp1);
    else N_VScale(ONE, p, vtemp1);
    if (preOnRight) {
      N_VScale(ONE, vtemp1, v);
      ier = psolve(P_data, v, vtemp1, delta, SUN_PREC_RIGHT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
        return(LASTFLAG(S));
      }
    }
    ier = atimes(A_data, vtemp1, v);
    if (ier != 0) {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (ier < 0) ?
        SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
      return(LASTFLAG(S));
    }
    if (preOnLeft) {
      ier = psolve(P_data, v, vtemp1, delta, SUN_PREC_LEFT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_REC;
        return(LASTFLAG(S));
      }
    }
    else N_VScale(ONE, v, vtemp1);
    if (scale_b) N_VProd(sb, vtemp1, v);
    else N_VScale(ONE, vtemp1, v);

    /* Shift variable values */
    /* NOTE: reduces storage requirements */
    N_VScale(ONE, r[1], r[0]);
    rho[0] = rho[1];

  }  /* END outer loop */

  /* Determine return value */
  /* If iteration converged or residual was reduced, then return current iterate (x) */
  if ((converged == SUNTRUE) || (r_curr_norm < r_init_norm)) {
    if (scale_x) N_VDiv(x, sx, x);
    if (preOnRight) {
      ier = psolve(P_data, x, vtemp1, delta, SUN_PREC_RIGHT);
      if (ier != 0) {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (ier < 0) ?
          SUNLS_PSOLVE_FAIL_UNREC : SUNLS_PSOLVE_FAIL_UNREC;
        return(LASTFLAG(S));
      }
      N_VScale(ONE, vtemp1, x);
    }

    *zeroguess = SUNFALSE;
    if (converged == SUNTRUE)
      LASTFLAG(S) = SUNLS_SUCCESS;
    else
      LASTFLAG(S) = SUNLS_RES_REDUCED;
    return(LASTFLAG(S));
  }
  /* Otherwise, return error code */
  else {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_CONV_FAIL;
    return(LASTFLAG(S));
  }
}


int SUNLinSolNumIters_SPTFQMR(SUNLinearSolver S)
{
  /* return the stored 'numiters' value */
  if (S == NULL) return(-1);
  return (SPTFQMR_CONTENT(S)->numiters);
}


realtype SUNLinSolResNorm_SPTFQMR(SUNLinearSolver S)
{
  /* return the stored 'resnorm' value */
  if (S == NULL) return(-ONE);
  return (SPTFQMR_CONTENT(S)->resnorm);
}


N_Vector SUNLinSolResid_SPTFQMR(SUNLinearSolver S)
{
  /* return the stored 'vtemp1' vector */
  return (SPTFQMR_CONTENT(S)->vtemp1);
}


sunindextype SUNLinSolLastFlag_SPTFQMR(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  if (S == NULL) return(-1);
  return (LASTFLAG(S));
}


int SUNLinSolSpace_SPTFQMR(SUNLinearSolver S,
                           long int *lenrwLS,
                           long int *leniwLS)
{
  sunindextype liw1, lrw1;
  if (SPTFQMR_CONTENT(S)->vtemp1->ops->nvspace)
    N_VSpace(SPTFQMR_CONTENT(S)->vtemp1, &lrw1, &liw1);
  else
    lrw1 = liw1 = 0;
  *lenrwLS = lrw1*11;
  *leniwLS = liw1*11;
  return(SUNLS_SUCCESS);
}

int SUNLinSolFree_SPTFQMR(SUNLinearSolver S)
{
  if (S == NULL) return(SUNLS_SUCCESS);

  if (S->content) {
    /* delete items from within the content structure */
    if (SPTFQMR_CONTENT(S)->r_star) {
      N_VDestroy(SPTFQMR_CONTENT(S)->r_star);
      SPTFQMR_CONTENT(S)->r_star = NULL;
    }
    if (SPTFQMR_CONTENT(S)->q) {
      N_VDestroy(SPTFQMR_CONTENT(S)->q);
      SPTFQMR_CONTENT(S)->q = NULL;
    }
    if (SPTFQMR_CONTENT(S)->d) {
      N_VDestroy(SPTFQMR_CONTENT(S)->d);
      SPTFQMR_CONTENT(S)->d = NULL;
    }
    if (SPTFQMR_CONTENT(S)->v) {
      N_VDestroy(SPTFQMR_CONTENT(S)->v);
      SPTFQMR_CONTENT(S)->v = NULL;
    }
    if (SPTFQMR_CONTENT(S)->p) {
      N_VDestroy(SPTFQMR_CONTENT(S)->p);
      SPTFQMR_CONTENT(S)->p = NULL;
    }
    if (SPTFQMR_CONTENT(S)->r) {
      N_VDestroyVectorArray(SPTFQMR_CONTENT(S)->r, 2);
      SPTFQMR_CONTENT(S)->r = NULL;
    }
    if (SPTFQMR_CONTENT(S)->u) {
      N_VDestroy(SPTFQMR_CONTENT(S)->u);
      SPTFQMR_CONTENT(S)->u = NULL;
    }
    if (SPTFQMR_CONTENT(S)->vtemp1) {
      N_VDestroy(SPTFQMR_CONTENT(S)->vtemp1);
      SPTFQMR_CONTENT(S)->vtemp1 = NULL;
    }
    if (SPTFQMR_CONTENT(S)->vtemp2) {
      N_VDestroy(SPTFQMR_CONTENT(S)->vtemp2);
      SPTFQMR_CONTENT(S)->vtemp2 = NULL;
    }
    if (SPTFQMR_CONTENT(S)->vtemp3) {
      N_VDestroy(SPTFQMR_CONTENT(S)->vtemp3);
      SPTFQMR_CONTENT(S)->vtemp3 = NULL;
    }
    free(S->content); S->content = NULL;
  }
  if (S->ops) { free(S->ops); S->ops = NULL; }
  free(S); S = NULL;
  return(SUNLS_SUCCESS);
}


int SUNLinSolSetInfoFile_SPTFQMR(SUNLinearSolver S,
                                 FILE* info_file)
{
#ifdef SUNDIALS_BUILD_WITH_MONITORING
  /* check that the linear solver is non-null */
  if (S == NULL)
    return(SUNLS_MEM_NULL);

  SPTFQMR_CONTENT(S)->info_file = info_file;

  return(SUNLS_SUCCESS);
#else
  SUNDIALS_DEBUG_PRINT("ERROR in SUNLinSolSetInfoFile_SPTFQMR: SUNDIALS was not built with monitoring\n");
  return(SUNLS_ILL_INPUT);
#endif
}


int SUNLinSolSetPrintLevel_SPTFQMR(SUNLinearSolver S,
                                   int print_level)
{
#ifdef SUNDIALS_BUILD_WITH_MONITORING
  /* check that the linear solver is non-null */
  if (S == NULL)
    return(SUNLS_MEM_NULL);

  /* check for valid print level */
  if (print_level < 0 || print_level > 1)
    return(SUNLS_ILL_INPUT);

  SPTFQMR_CONTENT(S)->print_level = print_level;

  return(SUNLS_SUCCESS);
#else
  SUNDIALS_DEBUG_PRINT("ERROR in SUNLinSolSetPrintLevel_SPTFQMR: SUNDIALS was not built with monitoring\n");
  return(SUNLS_ILL_INPUT);
#endif
}
