/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * This is the interface between ARKStep and the
 * SUNNonlinearSolver object
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_arkstep_impl.h"
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
  ARKStepSetNonlinearSolver:

  This routine attaches a SUNNonlinearSolver object to the ARKStep
  module.
  ---------------------------------------------------------------*/
int ARKStepSetNonlinearSolver(void *arkode_mem, SUNNonlinearSolver NLS)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "ARKStepSetNonlinearSolver",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return immediately if NLS input is NULL */
  if (NLS == NULL) {
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetNonlinearSolver",
                    "The NLS input must be non-NULL");
    return(ARK_ILL_INPUT);
  }

  /* check for required nonlinear solver functions */
  if ( (NLS->ops->gettype    == NULL) ||
       (NLS->ops->initialize == NULL) ||
       (NLS->ops->solve      == NULL) ||
       (NLS->ops->free       == NULL) ||
       (NLS->ops->setsysfn   == NULL) ) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "ARKStepSetNonlinearSolver",
                    "NLS does not support required operations");
    return(ARK_ILL_INPUT);
  }

  /* free any existing nonlinear solver */
  if (step_mem->NLS)
    retval = SUNNonlinSolFree(step_mem->NLS);

  /* set SUNNonlinearSolver pointer */
  step_mem->NLS = NLS;
  step_mem->ownNLS = SUNFALSE;

  /* set the nonlinear residual/fixed-point function, based on solver type */
  if (SUNNonlinSolGetType(NLS) == SUNNONLINEARSOLVER_ROOTFIND) {
    retval = SUNNonlinSolSetSysFn(step_mem->NLS, arkStep_NlsResidual);
  } else if (SUNNonlinSolGetType(NLS) ==  SUNNONLINEARSOLVER_FIXEDPOINT) {
    retval = SUNNonlinSolSetSysFn(step_mem->NLS, arkStep_NlsFPFunction);
  } else {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetNonlinearSolver",
                    "Invalid nonlinear solver type");
    return(ARK_ILL_INPUT);
  }
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetNonlinearSolver",
                    "Setting nonlinear system function failed");
    return(ARK_ILL_INPUT);
  }

  /* set convergence test function */
  retval = SUNNonlinSolSetConvTestFn(step_mem->NLS, arkStep_NlsConvTest);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetNonlinearSolver",
                    "Setting convergence test function failed");
    return(ARK_ILL_INPUT);
  }

  /* set default nonlinear iterations */
  retval = SUNNonlinSolSetMaxIters(step_mem->NLS, step_mem->maxcor);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "ARKStepSetNonlinearSolver",
                    "Setting maximum number of nonlinear iterations failed");
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  Utility routines called by ARKStep 
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  arkStep_NlsInit:

  This routine attaches the linear solver 'setup' and 'solve'
  routines to the nonlinear solver object, and then initializes
  the nonlinear solver object itself.  This should only be
  called at the start of a simulation, after a re-init, or after 
  a re-size.
  ---------------------------------------------------------------*/
int arkStep_NlsInit(ARKodeMem ark_mem)
{
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_NlsInit", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

  /* set the linear solver setup wrapper function */
  if (step_mem->lsetup)
    retval = SUNNonlinSolSetLSetupFn(step_mem->NLS, arkStep_NlsLSetup);
  else
    retval = SUNNonlinSolSetLSetupFn(step_mem->NLS, NULL);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_NlsInit",
                    "Setting the linear solver setup function failed");
    return(ARK_NLS_INIT_FAIL);
  }

  /* set the linear solver solve wrapper function */
  if (step_mem->lsolve)
    retval = SUNNonlinSolSetLSolveFn(step_mem->NLS, arkStep_NlsLSolve);
  else
    retval = SUNNonlinSolSetLSolveFn(step_mem->NLS, NULL);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_NlsInit",
                    "Setting linear solver solve function failed");
    return(ARK_NLS_INIT_FAIL);
  }

  /* initialize nonlinear solver */
  retval = SUNNonlinSolInitialize(step_mem->NLS);
  if (retval != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode::ARKStep",
                    "arkStep_NlsInit", MSG_NLS_INIT_FAIL);
    return(ARK_NLS_INIT_FAIL);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_Nls

  This routine attempts to solve the nonlinear system associated
  with a single implicit step of the linear multistep method.
  It calls the supplied SUNNonlinearSolver object to perform the
  solve.

  Upon entry, the predicted solution is held in step_mem->zpred;
  this array is never changed throughout this routine.  If an
  initial attempt at solving the nonlinear system fails (e.g. due
  to a stale Jacobian), this allows for new attempts at the
  solution.

  Upon a successful solve, the solution is held in ark_mem->ycur.
  ---------------------------------------------------------------*/
int arkStep_Nls(ARKodeMem ark_mem, int nflag)
{
  ARKodeARKStepMem step_mem;
  booleantype callLSetup;
  N_Vector zcor0;
  int retval;

  /* access ARKodeARKStepMem structure */
  if (ark_mem->step_mem==NULL) {
    arkProcessError(ark_mem, ARK_MEM_NULL, "ARKode::ARKStep",
                    "arkStep_Nls", MSG_ARKSTEP_NO_MEM);
    return(ARK_MEM_NULL);
  }
  step_mem = (ARKodeARKStepMem) ark_mem->step_mem;

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
  
  /* call nonlinear solver based on method type:
     FP methods solve for the updated solution directly, but
     Newton uses predictor-corrector form */
  if (SUNNonlinSolGetType(step_mem->NLS) == SUNNONLINEARSOLVER_FIXEDPOINT) {
    
    /* solve the nonlinear system, place solution directly in ycur */
    retval = SUNNonlinSolSolve(step_mem->NLS, step_mem->zpred, ark_mem->ycur, ark_mem->ewt,
                               step_mem->nlscoef, callLSetup, ark_mem);
    
  } else {
    
    /* set a zero guess for correction */
    zcor0 = ark_mem->tempv4;
    N_VConst(ZERO, zcor0);
    
    /* Reset the stored residual norm (for iterative linear solvers) */
    step_mem->eRNrm = RCONST(0.1) * step_mem->nlscoef;
    
    /* solve the nonlinear system for the actual correction */
    retval = SUNNonlinSolSolve(step_mem->NLS, zcor0, step_mem->zcor, ark_mem->ewt,
                               step_mem->nlscoef, callLSetup, ark_mem);
    
    /* apply the correction to construct ycur */
    N_VLinearSum(ONE, step_mem->zcor, ONE, step_mem->zpred, ark_mem->ycur);
    
  }
  
  /* on successful solve, reset the jcur flag */
  if (retval == ARK_SUCCESS)  step_mem->jcur = SUNFALSE;
  
  return(retval);
}


/*---------------------------------------------------------------
  Interface routines supplied to SUNNonlinearSolver module
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  arkStep_NlsLSetup:

  This routine wraps the ARKode linear solver interface 'setup'
  routine for use by the nonlinear solver object.
  ---------------------------------------------------------------*/
int arkStep_NlsLSetup(N_Vector zcor, N_Vector res, booleantype jbad,
                      booleantype* jcur, void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_NlsLSetup",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* update convfail based on jbad flag */
  if (jbad)
    step_mem->convfail = ARK_FAIL_BAD_J;

  /* Use ARKode's tempv1, tempv2 and tempv3 as
     temporary vectors for the linear solver setup routine */
  step_mem->nsetups++;
  retval = step_mem->lsetup(ark_mem, step_mem->convfail, ark_mem->tcur,
                            ark_mem->ycur, step_mem->Fi[step_mem->istage],
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
  arkStep_NlsLSolve:

  This routine wraps the ARKode linear solver interface 'solve'
  routine for use by the nonlinear solver object.
  ---------------------------------------------------------------*/
int arkStep_NlsLSolve(N_Vector zcor, N_Vector b, void* arkode_mem)
{
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval, nonlin_iter;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_NlsLSolve",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* retrieve nonlinear solver iteration from module */
  retval = SUNNonlinSolGetCurIter(step_mem->NLS, &nonlin_iter);
  if (retval != SUN_NLS_SUCCESS)
    return(ARK_NLS_OP_ERR);

  /* call linear solver interface, and handle return value */
  retval = step_mem->lsolve(ark_mem, b, ark_mem->tcur,
                            ark_mem->ycur, step_mem->Fi[step_mem->istage],
                            step_mem->eRNrm, nonlin_iter);

  if (retval < 0) return(ARK_LSOLVE_FAIL);
  if (retval > 0) return(retval);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_NlsResidual:

  This routine evaluates the nonlinear residual for the additive 
  Runge-Kutta method.  It assumes that any data from previous 
  time steps/stages is contained in step_mem, and merely combines 
  this old data with the current implicit ODE RHS vector to 
  compute the nonlinear residual r.

  At the ith stage, we compute the residual vector:
     r = M*z - M*yn - h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                    - h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
     r = M*zp + M*zc - M*yn - h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
                            - h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
     r = (M*zc - gamma*Fi(z)) - (M*yn - M*zp + data)
  where the current stage solution z = zp + zc, and where
     zc is stored in the input, zcor
     (M*yn-M*zp+data) is stored in step_mem->sdata,
  so we really just compute:
     z = zp + zc (stored in ark_mem->ycur)
     Fi(z) (stored step_mem->Fi[step_mem->istage])
     r = M*zc - gamma*Fi(z) - step_mem->sdata
  ---------------------------------------------------------------*/
int arkStep_NlsResidual(N_Vector zcor, N_Vector r, void* arkode_mem)
{
  /* temporary variables */
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  int retval;
  realtype c[3];
  N_Vector X[3];

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_NlsResidual",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* update 'ycur' value as stored predictor + current corrector */
  N_VLinearSum(ONE, step_mem->zpred, ONE, zcor, ark_mem->ycur);

  /* compute implicit RHS and save for later */
  retval = step_mem->fi(ark_mem->tcur, ark_mem->ycur,
                        step_mem->Fi[step_mem->istage],
                        ark_mem->user_data);
  step_mem->nfi++;
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  /* put M*zcor in r */
  if (step_mem->mass_mem != NULL) {
    retval = step_mem->mmult((void *) ark_mem, zcor, r);
    if (retval != ARK_SUCCESS)  return (ARK_MASSMULT_FAIL);
    X[0] = r;
  } else {
    X[0] = zcor;
  }

  /* update with My, sdata and gamma*fy */
  c[0] = ONE;
  c[1] = -ONE;
  X[1] = step_mem->sdata;
  c[2] = -step_mem->gamma;
  X[2] = step_mem->Fi[step_mem->istage];
  retval = N_VLinearCombination(3, c, X, r);
  if (retval != 0)  return(ARK_VECTOROP_ERR);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_NlsFPFunction:

  This routine evaluates the fixed point iteration function for
  the additive Runge-Kutta method.  It assumes that any data from
  previous time steps/stages is contained in step_mem, and
  merely combines this old data with the current guess and
  implicit ODE RHS vector to compute the iteration function g.

  At the ith stage, the new stage solution z should solve:
     z = yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j)
            + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
  <=>
     z = yn + gamma*Fi(z) + h*sum_{j=0}^{i-1} ( Ae(i,j)*Fe(j)
                                              + Ai(i,j)*Fi(j) )
  <=>
     z = yn + gamma*Fi(z) + data
  Our fixed-point problem is z=g(z), so the FP function is just:
     g(z) = yn + gamma*Fi(z) + data
  <=>
     g(z) = zp - zp + gamma*Fi(z) + yn + data
  <=>
     g(z) = zp + gamma*Fi(z) + (yn - zp + data)
  where the current nonlinear guess is z = zp + zc, and where
     z is stored in the input, z,
     zp is stored in step_mem->zpred,
     (yn-zp+data) is stored in step_mem->sdata,
  so we really just compute:
     Fi(z) (store in step_mem->Fi[step_mem->istage])
     g = zp + gamma*Fi(z) + step_mem->sdata
  ---------------------------------------------------------------*/
int arkStep_NlsFPFunction(N_Vector z, N_Vector g, void* arkode_mem)
{
  /* temporary variables */
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  realtype c[3];
  N_Vector X[3];
  int retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_NlsFPFunction",
                                 &ark_mem, &step_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* compute implicit RHS and save for later */
  retval = step_mem->fi(ark_mem->tcur, z,
                        step_mem->Fi[step_mem->istage],
                        ark_mem->user_data);
  step_mem->nfi++;
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  /* combine parts:  g = zpred + sdata + gamma*Fi(z) */
  c[0] = ONE;
  X[0] = step_mem->zpred;
  c[1] = ONE;
  X[1] = step_mem->sdata;
  c[2] = step_mem->gamma;
  X[2] = step_mem->Fi[step_mem->istage];
  retval = N_VLinearCombination(3, c, X, g);
  if (retval != 0)  return(ARK_VECTOROP_ERR);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkStep_NlsConvTest:

  This routine provides the nonlinear solver convergence test for
  the additive Runge-Kutta method.  We have two modes.

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
int arkStep_NlsConvTest(SUNNonlinearSolver NLS, N_Vector y, N_Vector del,
                        realtype tol, N_Vector ewt, void* arkode_mem)
{
  /* temporary variables */
  ARKodeMem ark_mem;
  ARKodeARKStepMem step_mem;
  realtype delnrm, dcon;
  int m, retval;

  /* access ARKodeARKStepMem structure */
  retval = arkStep_AccessStepMem(arkode_mem, "arkStep_NlsConvTest",
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
