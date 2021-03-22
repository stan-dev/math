/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the interface between MRIStep and the
 * SUNNonlinearSolver object
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

/* constants */
#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)


/*===============================================================
  Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepSetNonlinearSolver:

  This routine attaches a SUNNonlinearSolver object to the MRIStep
  module.
  ---------------------------------------------------------------*/
int MRIStepSetNonlinearSolver(void *arkode_mem, SUNNonlinearSolver NLS)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepSetNonlinearSolver",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return immediately if NLS input is NULL */
  if (NLS == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepSetNonlinearSolver",
                    "The NLS input must be non-NULL");
    return(ARK_ILL_INPUT);
  }

  /* check for required nonlinear solver functions */
  if ( (NLS->ops->gettype    == NULL) ||
       (NLS->ops->solve      == NULL) ||
       (NLS->ops->setsysfn   == NULL) ) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "MRIStepSetNonlinearSolver",
                    "NLS does not support required operations");
    return(ARK_ILL_INPUT);
  }

  /* free any existing nonlinear solver */
  if ((step_mem->NLS != NULL) && (step_mem->ownNLS))
    retval = SUNNonlinSolFree(step_mem->NLS);

  /* set SUNNonlinearSolver pointer */
  step_mem->NLS = NLS;
  step_mem->ownNLS = SUNFALSE;

  /* set the nonlinear residual/fixed-point function, based on solver type */
  if (SUNNonlinSolGetType(NLS) == SUNNONLINEARSOLVER_ROOTFIND) {
    retval = SUNNonlinSolSetSysFn(step_mem->NLS, mriStep_NlsResidual);
  } else if (SUNNonlinSolGetType(NLS) ==  SUNNONLINEARSOLVER_FIXEDPOINT) {
    retval = SUNNonlinSolSetSysFn(step_mem->NLS, mriStep_NlsFPFunction);
  } else {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepSetNonlinearSolver",
                    "Invalid nonlinear solver type");
    return(ARK_ILL_INPUT);
  }
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepSetNonlinearSolver",
                    "Setting nonlinear system function failed");
    return(ARK_ILL_INPUT);
  }

  /* set convergence test function */
  retval = SUNNonlinSolSetConvTestFn(step_mem->NLS, mriStep_NlsConvTest,
                                     arkode_mem);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepSetNonlinearSolver",
                    "Setting convergence test function failed");
    return(ARK_ILL_INPUT);
  }

  /* set default nonlinear iterations */
  retval = SUNNonlinSolSetMaxIters(step_mem->NLS, step_mem->maxcor);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "MRIStepSetNonlinearSolver",
                    "Setting maximum number of nonlinear iterations failed");
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  MRIStepGetNonlinearSystemData:

  This routine provides access to the relevant data needed to
  compute the nonlinear system function.
  ---------------------------------------------------------------*/
int MRIStepGetNonlinearSystemData(void *arkode_mem, realtype *tcur,
                                  N_Vector *zpred, N_Vector *z,
                                  N_Vector *F, realtype *gamma,
                                  N_Vector *sdata, void **user_data)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "MRIStepGetNonlinearSystemData",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  *tcur      = ark_mem->tcur;
  *zpred     = step_mem->zpred;
  *z         = ark_mem->ycur;
  *F         = step_mem->F[step_mem->istage];
  *gamma     = step_mem->gamma;
  *sdata     = step_mem->sdata;
  *user_data = ark_mem->user_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Utility routines called by MRIStep
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  mriStep_NlsInit:

  This routine attaches the linear solver 'setup' and 'solve'
  routines to the nonlinear solver object, and then initializes
  the nonlinear solver object itself.  This should only be
  called at the start of a simulation, after a re-init, or after
  a re-size.
  ---------------------------------------------------------------*/
int mriStep_NlsInit(ARKodeMem ark_mem)
{
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "mriStep_NlsInit", MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* reset counters */
  step_mem->nls_iters = 0;

  /* set the linear solver setup wrapper function */
  if (step_mem->lsetup)
    retval = SUNNonlinSolSetLSetupFn(step_mem->NLS, mriStep_NlsLSetup);
  else
    retval = SUNNonlinSolSetLSetupFn(step_mem->NLS, NULL);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_NlsInit",
                    "Setting the linear solver setup function failed");
    return(ARK_NLS_INIT_FAIL);
  }

  /* set the linear solver solve wrapper function */
  if (step_mem->lsolve)
    retval = SUNNonlinSolSetLSolveFn(step_mem->NLS, mriStep_NlsLSolve);
  else
    retval = SUNNonlinSolSetLSolveFn(step_mem->NLS, NULL);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_NlsInit",
                    "Setting linear solver solve function failed");
    return(ARK_NLS_INIT_FAIL);
  }

  /* initialize nonlinear solver */
  retval = SUNNonlinSolInitialize(step_mem->NLS);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::MRIStep",
                    "mriStep_NlsInit", MSG_NLS_INIT_FAIL);
    return(ARK_NLS_INIT_FAIL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_Nls

  This routine attempts to solve the nonlinear system associated
  with a single solve-decoupled implicit stage. It calls the
  supplied SUNNonlinearSolver object to perform the solve.

  Upon entry, the predicted solution is held in step_mem->zpred,
  which is never changed throughout this routine.  If an initial
  attempt at solving the nonlinear system fails (e.g. due to a
  stale Jacobian), this allows for new attempts at the solution.

  Upon a successful solve, the solution is held in ark_mem->ycur.
  ---------------------------------------------------------------*/
int mriStep_Nls(ARKodeMem ark_mem, int nflag)
{
  ARKodeMRIStepMem step_mem;
  booleantype callLSetup;
  long int nls_iters_inc;
  int retval;

  /* access ARKodeMRIStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::MRIStep",
                    "mriStep_Nls", MSG_MRISTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeMRIStepMem) ark_mem->step_mem;

  /* If a linear solver 'setup' is supplied, set various flags for
     determining whether it should be called */
  if (step_mem->lsetup) {

    /* Set interface 'convfail' flag for use inside lsetup */
    if (step_mem->linear) {
      step_mem->convfail = (nflag == FIRST_CALL) ? ARK_NO_FAILURES : ARK_FAIL_OTHER;
    } else {
      step_mem->convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
        ARK_NO_FAILURES : ARK_FAIL_OTHER;
    }

    /* Decide whether to recommend call to lsetup within nonlinear solver */
    callLSetup = (ark_mem->firststage) || (step_mem->msbp < 0) ||
      (SUNRabs(step_mem->gamrat-ONE) > step_mem->dgmax);
    if (step_mem->linear) {   /* linearly-implicit problem */
      callLSetup = callLSetup || (step_mem->linear_timedep);
    } else {                  /* nonlinearly-implicit problem */
      callLSetup = callLSetup ||
        (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
        (ark_mem->nst >= step_mem->nstlp + abs(step_mem->msbp));
    }
  } else {
    step_mem->crate = ONE;
    callLSetup = SUNFALSE;
  }

  /* set a zero guess for correction */
  N_VConst(ZERO, step_mem->zcor);

  /* Reset the stored residual norm (for iterative linear solvers) */
  step_mem->eRNrm = RCONST(0.1) * step_mem->nlscoef;

  /* solve the nonlinear system for the actual correction */
  retval = SUNNonlinSolSolve(step_mem->NLS, step_mem->zpred, step_mem->zcor,
                             ark_mem->ewt, step_mem->nlscoef, callLSetup, ark_mem);

  /* apply the correction to construct ycur */
  N_VLinearSum(ONE, step_mem->zcor, ONE, step_mem->zpred, ark_mem->ycur);

  /* increment counter */
  nls_iters_inc = 0;
  (void) SUNNonlinSolGetNumIters(step_mem->NLS, &(nls_iters_inc));
  step_mem->nls_iters += nls_iters_inc;

  /* on successful solve, reset the jcur flag */
  if (retval == ARK_SUCCESS)  step_mem->jcur = SUNFALSE;

  /* if convergence failure, return ARKode::CONV_FAIL */
  if (retval == SUN_NLS_CONV_RECVR) return(CONV_FAIL);

  return(retval);
}


/*---------------------------------------------------------------
  Interface routines supplied to the SUNNonlinearSolver module
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  mriStep_NlsLSetup:

  This routine wraps the ARKode linear solver interface 'setup'
  routine for use by the nonlinear solver object.
  ---------------------------------------------------------------*/
int mriStep_NlsLSetup(booleantype jbad, booleantype* jcur, void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_NlsLSetup",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* update convfail based on jbad flag */
  if (jbad)  step_mem->convfail = ARK_FAIL_BAD_J;

  /* Use ARKode's tempv1, tempv2 and tempv3 as
     temporary vectors for the linear solver setup routine */
  step_mem->nsetups++;
  retval = step_mem->lsetup(ark_mem, step_mem->convfail, ark_mem->tcur,
                            ark_mem->ycur, step_mem->F[step_mem->istage],
                            &(step_mem->jcur), ark_mem->tempv1,
                            ark_mem->tempv2, ark_mem->tempv3);

  /* update Jacobian status */
  *jcur = step_mem->jcur;

  /* update flags and 'gamma' values for last lsetup call */
  ark_mem->firststage = SUNFALSE;
  step_mem->gamrat = step_mem->crate = ONE;
  step_mem->gammap = step_mem->gamma;
  step_mem->nstlp  = ark_mem->nst;

  if (retval < 0) return(ARK_LSETUP_FAIL);
  if (retval > 0) return(CONV_FAIL);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_NlsLSolve:

  This routine wraps the ARKode linear solver interface 'solve'
  routine for use by the nonlinear solver object.
  ---------------------------------------------------------------*/
int mriStep_NlsLSolve(N_Vector b, void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval, nonlin_iter;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_NlsLSolve",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* retrieve nonlinear solver iteration from module */
  retval = SUNNonlinSolGetCurIter(step_mem->NLS, &nonlin_iter);
  if (retval != SUN_NLS_SUCCESS)
    return(ARK_NLS_OP_ERR);

  /* call linear solver interface, and handle return value */
  retval = step_mem->lsolve(ark_mem, b, ark_mem->tcur,
                            ark_mem->ycur, step_mem->F[step_mem->istage],
                            step_mem->eRNrm, nonlin_iter);

  if (retval < 0) return(ARK_LSOLVE_FAIL);
  if (retval > 0) return(CONV_FAIL);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_NlsResidual:

  This routine evaluates the nonlinear residual for this
  solve-decoupled implicit MRI stage.  It assumes that any data
  from previous time steps/stages is contained in step_mem, and
  merely combines this old data with the current implicit ODE
  RHS vector to compute the nonlinear residual r.

  At the ith stage, we compute the residual vector:
     r = zc - gamma*F(z) - sdata
  where the current stage solution is z = zp + zc,
     gamma = h*A(i,i),
     zc is stored in the input, zcor, and
     sdata is the old solution/stage data stored in step_mem->sdata.
  Hence we really just compute:
     z = zp + zc (stored in ark_mem->ycur)
     F(z) (stored step_mem->F[step_mem->istage])
     r = zc - gamma*F(z) - step_mem->sdata
  ---------------------------------------------------------------*/
int mriStep_NlsResidual(N_Vector zcor, N_Vector r, void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;
  realtype c[3];
  N_Vector X[3];

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_NlsResidual",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* update 'ycur' value as stored predictor + current corrector */
  N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, ark_mem->ycur);

  /* compute slow RHS and save for later */
  retval = step_mem->fs(ark_mem->tcur, ark_mem->ycur,
                        step_mem->F[step_mem->istage],
                        ark_mem->user_data);
  step_mem->nfs++;
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  /* compute residual: zcor - gamma*F - sdata */
  c[0] = ONE;
  X[0] = zcor;
  c[1] = -ONE;
  X[1] = step_mem->sdata;
  c[2] = -step_mem->gamma;
  X[2] = step_mem->F[step_mem->istage];
  retval = N_VLinearCombination(3, c, X, r);
  if (retval != 0)  return(ARK_VECTOROP_ERR);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_NlsFPFunction:

  This routine evaluates the fixed point iteration function for
  this solve-decoupled implicit MRI stage.  It assumes that any
  data from previous time steps/stages is contained in step_mem,
  and merely combines this old data with the current guess and
  current slow RHS vector to compute the iteration function g.

  At the ith stage, the new stage solution z=(zc+zp) should solve:
     zc = g(zc) := gamma*F(z) + sdata
  where
     gamma = h*A(i,i),
     zp is the predicted stage solution,
     zc is stored in the input, zcor, and
     sdata is the old solution/stage data stored in step_mem->sdata.
  So we really just compute:
     z = zp + zc (stored in ark_mem->ycur)
     F(z) (store in step_mem->F[step_mem->istage])
     g = gamma*F(z) + step_mem->sdata
  ---------------------------------------------------------------*/
int mriStep_NlsFPFunction(N_Vector zcor, N_Vector g, void* arkode_mem)
{
  /* temporary variables */
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  int retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_NlsFPFunction",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* update 'ycur' value as stored predictor + current corrector */
  N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, ark_mem->ycur);

  /* compute slow RHS and save for later */
  retval = step_mem->fs(ark_mem->tcur, ark_mem->ycur,
                        step_mem->F[step_mem->istage],
                        ark_mem->user_data);
  step_mem->nfs++;
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  /* combine parts:  g = gamma*F(z) + sdata */
  N_VLinearSum(step_mem->gamma, step_mem->F[step_mem->istage],
               ONE, step_mem->sdata, g);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  mriStep_NlsConvTest:

  This routine provides the nonlinear solver convergence test for
  this solve-decoupled implicit MRI stage.  We have two modes.

  Standard:
      delnorm = ||del||_WRMS
      if (m==0) crate = 1
      if (m>0)  crate = max(crdown*crate, delnorm/delp)
      dcon = min(crate, ONE) * del / nlscoef
      if (dcon<=1)  return convergence
      if ((m >= 2) && (del > rdiv*delp))  return divergence

  Linearly-implicit mode:
      if the user specifies that the problem is linearly
      implicit, then we just declare 'success' no matter what
      is provided.
  ---------------------------------------------------------------*/
int mriStep_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                        realtype tol, N_Vector ewt, void* arkode_mem)
{
  /* temporary variables */
  ARKodeMem ark_mem;
  ARKodeMRIStepMem step_mem;
  realtype delnrm, dcon;
  int m, retval;

  /* access ARKodeMRIStepMem structure */
  retval = mriStep_AccessStepMem(arkode_mem, "mriStep_NlsConvTest",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* if the problem is linearly implicit, just return success */
  if (step_mem->linear)
    return(SUN_NLS_SUCCESS);

  /* compute the norm of the correction */
  delnrm = N_VWrmsNorm(del, ewt);

  /* get the current nonlinear solver iteration count */
  retval = SUNNonlinSolGetCurIter(NLS, &m);
  if (retval != ARK_SUCCESS)  return(ARK_MEM_NULL);

  /* update the stored estimate of the convergence rate (assumes linear convergence) */
  if (m > 0)
    step_mem->crate = SUNMAX(step_mem->crdown*step_mem->crate, delnrm/step_mem->delp);

  /* compute our scaled error norm for testing convergence */
  dcon = SUNMIN(step_mem->crate, ONE) * delnrm / tol;

  /* check for convergence; if so return with success */
  if (dcon <= ONE)  return(SUN_NLS_SUCCESS);

  /* check for divergence */
  if ((m >= 1) && (delnrm > step_mem->rdiv*step_mem->delp))
    return(SUN_NLS_CONV_RECVR);

  /* save norm of correction for next iteration */
  step_mem->delp = delnrm;

  /* return with flag that there is more work to do */
  return(SUN_NLS_CONTINUE);
}


/*===============================================================
  EOF
  ===============================================================*/
