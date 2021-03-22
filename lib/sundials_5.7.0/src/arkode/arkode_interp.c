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
 * This is the implementation file for ARKode's temporal
 * interpolation utility.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "arkode_impl.h"
#include "arkode_interp_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif


/*---------------------------------------------------------------
  Section I: generic ARKInterp functions provided by all
  interpolation modules
  ---------------------------------------------------------------*/

int arkInterpResize(void* arkode_mem, ARKInterp interp,
                    ARKVecResizeFn resize, void *resize_data,
                    sunindextype lrw_diff, sunindextype liw_diff,
                    N_Vector tmpl)
{
  if (interp == NULL)  return(ARK_SUCCESS);
  return((int) interp->ops->resize(arkode_mem, interp,
                                   resize, resize_data,
                                   lrw_diff, liw_diff, tmpl));
}

void arkInterpFree(void* arkode_mem, ARKInterp interp)
{
  if (interp == NULL)  return;
  interp->ops->free(arkode_mem, interp);
  return;
}

void arkInterpPrintMem(ARKInterp interp, FILE *outfile)
{
  if (interp == NULL)  return;
  interp->ops->print(interp, outfile);
  return;
}

int arkInterpSetDegree(void* arkode_mem, ARKInterp interp,
                       int degree)
{
  if (interp == NULL)  return(ARK_SUCCESS);
  return((int) interp->ops->setdegree(arkode_mem, interp, degree));
}

int arkInterpInit(void* arkode_mem, ARKInterp interp,
                  realtype tnew)
{
  if (interp == NULL)  return(ARK_SUCCESS);
  return((int) interp->ops->init(arkode_mem, interp, tnew));
}

int arkInterpUpdate(void* arkode_mem, ARKInterp interp, realtype tnew)
{
  if (interp == NULL)  return(ARK_SUCCESS);
  return((int) interp->ops->update(arkode_mem, interp, tnew));
}

int arkInterpEvaluate(void* arkode_mem, ARKInterp interp,
                      realtype tau, int d, int order, N_Vector yout)
{
  if (interp == NULL)  return(ARK_SUCCESS);
  return((int) interp->ops->evaluate(arkode_mem, interp,
                                     tau, d, order, yout));
}



/*---------------------------------------------------------------
  Section II: Hermite interpolation module implementation
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  arkInterpCreate_Hermite:

  This routine creates an ARKInterp structure, through
  cloning an input template N_Vector.  This returns a non-NULL
  structure if no errors occurred, or a NULL value otherwise.
  ---------------------------------------------------------------*/
ARKInterp arkInterpCreate_Hermite(void* arkode_mem, int degree)
{
  ARKInterp interp;
  ARKInterpContent_Hermite content;
  ARKInterpOps ops;
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for valid degree */
  if (degree < 0 || degree > ARK_INTERP_MAX_DEGREE) return(NULL);

  /* allocate overall structure */
  interp = NULL;
  interp = (ARKInterp) malloc(sizeof *interp);
  if (interp == NULL)  return(NULL);

  /* allocate ops structure and set entries */
  ops = NULL;
  ops = (ARKInterpOps) malloc(sizeof *ops);
  if (ops == NULL) { free(interp); return(NULL); }
  ops->resize    = arkInterpResize_Hermite;
  ops->free      = arkInterpFree_Hermite;
  ops->print     = arkInterpPrintMem_Hermite;
  ops->setdegree = arkInterpSetDegree_Hermite;
  ops->init      = arkInterpInit_Hermite;
  ops->update    = arkInterpUpdate_Hermite;
  ops->evaluate  = arkInterpEvaluate_Hermite;

  /* create content, and initialize everything to zero/NULL */
  content = NULL;
  content = (ARKInterpContent_Hermite) malloc(sizeof *content);
  if (content == NULL) { free(ops); free(interp); return(NULL); }
  memset(content, 0, sizeof(struct _ARKInterpContent_Hermite));

  /* attach ops and content structures to overall structure */
  interp->ops = ops;
  interp->content = content;

  /* fill content */

  /* initialize local N_Vectors to NULL */
  content->fold = NULL;
  content->yold = NULL;
  content->fa   = NULL;
  content->fb   = NULL;

  /* set maximum interpolant degree */
  content->degree = SUNMIN(ARK_INTERP_MAX_DEGREE, degree);

  /* set ynew and fnew pointers to ark_mem->yn and ark_mem->fn, respectively */
  content->ynew = ark_mem->yn;
  content->fnew = ark_mem->fn;

  /* update workspace sizes */
  ark_mem->lrw += 2;
  ark_mem->liw += 5;

  /* initialize time values */
  content->told = ark_mem->tcur;
  content->tnew = ark_mem->tcur;
  content->h    = RCONST(0.0);

  return(interp);
}


/*---------------------------------------------------------------
  arkInterpResize_Hermite:

  This routine resizes the internal vectors.
  ---------------------------------------------------------------*/
int arkInterpResize_Hermite(void* arkode_mem, ARKInterp interp,
                            ARKVecResizeFn resize, void *resize_data,
                            sunindextype lrw_diff, sunindextype liw_diff,
                            N_Vector y0)
{
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* resize vectors */
  if (interp == NULL)  return(ARK_SUCCESS);

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, y0, &HINT_FOLD(interp)))
    return(ARK_MEM_FAIL);

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, y0, &HINT_YOLD(interp)))
    return(ARK_MEM_FAIL);

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, y0, &HINT_FA(interp)))
    return(ARK_MEM_FAIL);

  if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                    liw_diff, y0, &HINT_FB(interp)))
    return(ARK_MEM_FAIL);

  /* update ynew and fnew pointers */
  HINT_YNEW(interp) = ark_mem->yn;
  HINT_FNEW(interp) = ark_mem->fn;

  /* reinitialize time values */
  HINT_TOLD(interp) = ark_mem->tcur;
  HINT_TNEW(interp) = ark_mem->tcur;
  HINT_H(interp)    = RCONST(0.0);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkInterpFree_Hermite:

  This routine frees the Hermite ARKInterp structure.
  ---------------------------------------------------------------*/
void arkInterpFree_Hermite(void* arkode_mem, ARKInterp interp)
{
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return;
  ark_mem = (ARKodeMem) arkode_mem;

  /* if interpolation structure is NULL, just return */
  if (interp == NULL) return;

  /* free content */
  if (interp->content != NULL) {
    if (HINT_FOLD(interp) != NULL) {
      arkFreeVec(ark_mem, &(HINT_FOLD(interp)));
      HINT_FOLD(interp) = NULL;
    }
    if (HINT_YOLD(interp) != NULL) {
      arkFreeVec(ark_mem, &(HINT_YOLD(interp)));
      HINT_YOLD(interp) = NULL;
    }
    if (HINT_FA(interp) != NULL) {
      arkFreeVec(ark_mem, &(HINT_FA(interp)));
      HINT_FA(interp) = NULL;
    }
    if (HINT_FB(interp) != NULL) {
      arkFreeVec(ark_mem, &(HINT_FB(interp)));
      HINT_FB(interp) = NULL;
    }

    /* update work space sizes */
    ark_mem->lrw -= 2;
    ark_mem->liw -= 5;

    free(interp->content);
    interp->content = NULL;
  }

  /* free ops and interpolation structures */
  if (interp->ops) {
    free(interp->ops);
    interp->ops = NULL;
  }
  free(interp);
  interp = NULL;

  return;
}


/*---------------------------------------------------------------
  arkInterpPrintMem_Hermite

  This routine outputs the Hermite temporal interpolation memory
  structure to a specified file pointer.
  ---------------------------------------------------------------*/
void arkInterpPrintMem_Hermite(ARKInterp interp, FILE *outfile)
{
  if (interp != NULL) {
    fprintf(outfile, "arkode_interp (Hermite): degree = %d\n", HINT_DEGREE(interp));
    fprintf(outfile, "arkode_interp (Hermite): told = %"RSYM"\n", HINT_TOLD(interp));
    fprintf(outfile, "arkode_interp (Hermite): tnew = %"RSYM"\n", HINT_TNEW(interp));
    fprintf(outfile, "arkode_interp (Hermite): h = %"RSYM"\n", HINT_H(interp));
#ifdef SUNDIALS_DEBUG_PRINTVEC
    fprintf(outfile, "arkode_interp (Hermite): fold:\n");
    N_VPrintFile(HINT_FOLD(interp), outfile);
    fprintf(outfile, "arkode_interp (Hermite): fnew:\n");
    N_VPrintFile(HINT_FNEW(interp), outfile);
    fprintf(outfile, "arkode_interp (Hermite): yold:\n");
    N_VPrintFile(HINT_YOLD(interp), outfile);
    fprintf(outfile, "arkode_interp (Hermite): ynew:\n");
    N_VPrintFile(HINT_YNEW(interp), outfile);
    fprintf(outfile, "arkode_interp (Hermite): fa:\n");
    N_VPrintFile(HINT_FA(interp), outfile);
    fprintf(outfile, "arkode_interp (Hermite): fb:\n");
    N_VPrintFile(HINT_FB(interp), outfile);
#endif
  }
}


/*---------------------------------------------------------------
  arkInterpSetDegree_Hermite

  This routine sets a supplied interpolation degree.  If the
  argument is positive, then we require that
      0 <= degree <= ARK_INTERP_MAX_DEGREE
  and use this value as the user-specified (or default) degree.

  If the argument is negative, then we assume that this has been
  called by a time-step module to limit the interpolant degree
  based on the temporal method order. In this case we set the
  Hermite polynomial degree to be the minimum of (-degree),
  ARK_INTERP_MAX_DEGREE, and the previously-set value [i.e., in
  case the user has already specified use of a lower-degree
  polynomial].

  Return values:
    ARK_MEM_NULL -- if either arkode_mem or interp are NULL
    ARK_ILL_INPUT -- if the input is outside of allowable bounds
    ARK_INTERP_FAIL -- if the interpolation module has already
       been initialized,
    ARK_SUCCESS -- successful completion.
  ---------------------------------------------------------------*/
int arkInterpSetDegree_Hermite(void* arkode_mem, ARKInterp interp,
                               int degree)
{
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* if this degree is already stored, just return */
  if (abs(degree) == HINT_DEGREE(interp))  return(ARK_SUCCESS);

  /* on positive degree, check for allowable value and overwrite stored degree */
  if (degree >= 0) {
    if (degree > ARK_INTERP_MAX_DEGREE) {
      arkProcessError(ark_mem, ARK_INTERP_FAIL, "ARKode",
                      "arkInterpSetDegree_Hermite",
                      "Illegal degree specified.");
      return(ARK_ILL_INPUT);
    }

    HINT_DEGREE(interp) = degree;
    return(ARK_SUCCESS);
  }

  /* on negative degree, check for allowable value and update stored degree */
  degree = -degree;
  if (degree > ARK_INTERP_MAX_DEGREE) degree = ARK_INTERP_MAX_DEGREE;
  HINT_DEGREE(interp) = SUNMIN(HINT_DEGREE(interp), degree);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkInterpInit_Hermite

  This routine performs the following steps:
  1. Sets tnew and told to the input time
  2. Allocates any missing/needed N_Vector storage (for reinit)
  3. Copies ark_mem->yn into yold
  4. Calls the full RHS routine to fill fnew
  5. Copies fnew into fold
  ---------------------------------------------------------------*/
int arkInterpInit_Hermite(void* arkode_mem, ARKInterp interp,
                          realtype tnew)
{
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* initialize time values */
  HINT_TOLD(interp) = tnew;
  HINT_TNEW(interp) = tnew;
  HINT_H(interp)    = RCONST(0.0);

  /* allocate vectors based on interpolant degree */
  if (HINT_FOLD(interp) == NULL)
    if (!arkAllocVec(ark_mem, ark_mem->yn, &(HINT_FOLD(interp)))) {
      arkInterpFree(ark_mem, interp); return(ARK_MEM_FAIL);
    }
  if (HINT_YOLD(interp) == NULL)
    if (!arkAllocVec(ark_mem, ark_mem->yn, &(HINT_YOLD(interp)))) {
      arkInterpFree(ark_mem, interp); return(ARK_MEM_FAIL);
    }
  if ((HINT_DEGREE(interp) > 3) && (HINT_FA(interp) == NULL)) {
    if (!arkAllocVec(ark_mem, ark_mem->yn, &(HINT_FA(interp)))) {
      arkInterpFree(ark_mem, interp); return(ARK_MEM_FAIL);
    }
  }
  if ((HINT_DEGREE(interp) > 4) && (HINT_FB(interp) == NULL)) {
    if (!arkAllocVec(ark_mem, ark_mem->yn, &(HINT_FB(interp)))) {
      arkInterpFree(ark_mem, interp); return(ARK_MEM_FAIL);
    }
  }

  /* copy current solution into yold */
  N_VScale(ONE, ark_mem->yn, HINT_YOLD(interp));

  /* copy fnew into fold */
  N_VScale(ONE, HINT_FNEW(interp), HINT_FOLD(interp));

  /* signal that fullrhs is required after each step */
  ark_mem->call_fullrhs = SUNTRUE;

  /* return with success */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkInterpUpdate_Hermite

  This routine copies ynew into yold, and fnew into fold, so that
  yold and fold contain the previous values.
  ---------------------------------------------------------------*/
int arkInterpUpdate_Hermite(void* arkode_mem, ARKInterp interp, realtype tnew)
{
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* copy ynew and fnew into yold and fold, respectively */
  N_VScale(ONE, HINT_YNEW(interp), HINT_YOLD(interp));
  N_VScale(ONE, HINT_FNEW(interp), HINT_FOLD(interp));

  /* update time values */
  HINT_TOLD(interp) = HINT_TNEW(interp);
  HINT_TNEW(interp) = tnew;
  HINT_H(interp)    = ark_mem->h;

  /* return with success */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkInterpEvaluate_Hermite

  This routine evaluates a temporal interpolation/extrapolation
  based on the data in the interpolation structure:
     yold = y(told)
     ynew = y(tnew)
     fold = f(told, yold)
     fnew = f(told, ynew)
  This typically consists of using a cubic Hermite interpolating
  formula with this data.  If greater polynomial degree than 3 is
  requested, then we can bootstrap up to a 5th-order interpolant.
  For lower order interpolants than cubic, we use:
     {yold,ynew,fnew} for quadratic
     {yold,ynew} for linear
     {0.5*(yold+ynew)} for constant.

  Derivatives have lower accuracy than the interpolant
  itself, losing one order per derivative.  We will provide
  derivatives up to d = min(5,q).

  The input 'tau' specifies the time at which to evaluate the Hermite
  polynomial.  The formula for tau is defined using the
  most-recently-completed solution interval [told,tnew], and is
  given by:
               t = tnew + tau*(tnew-told),
  where h = tnew-told, i.e. values -1<tau<0 provide interpolation,
  other values result in extrapolation.
  ---------------------------------------------------------------*/
int arkInterpEvaluate_Hermite(void* arkode_mem, ARKInterp interp,
                              realtype tau, int d, int order, N_Vector yout)
{
  /* local variables */
  int q, retval;
  realtype tval, a0, a1, tau2, tau3, tau4, tau5;
  realtype h, h2, h3, h4, h5;
  realtype a[6];
  N_Vector X[6];
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* set constants */
  tau2 = tau*tau;
  tau3 = tau*tau2;
  tau4 = tau*tau3;
  tau5 = tau*tau4;

  h  = HINT_H(interp);
  h2 = h*h;
  h3 = h*h2;
  h4 = h*h3;
  h5 = h*h4;

  /* determine polynomial order q */
  q = SUNMAX(order, 0);                 /* respect lower bound  */
  q = SUNMIN(q, HINT_DEGREE(interp));   /* respect max possible */

  /* error on illegal d */
  if (d < 0) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInterpEvaluate_Hermite",
                    "Requested illegal derivative.");
    return (ARK_ILL_INPUT);
  }

  /* if d is too high, just return zeros */
  if (d > q) {
    N_VConst(ZERO, yout);
    return(ARK_SUCCESS);
  }

  /* build polynomial based on order */
  switch (q) {

  case(0):    /* constant interpolant, yout = 0.5*(yn+yp) */
    N_VLinearSum(HALF, HINT_YOLD(interp), HALF, HINT_YNEW(interp), yout);
    break;

  case(1):    /* linear interpolant */
    if (d == 0) {
      a0 = -tau;
      a1 = ONE+tau;
    } else {  /* d=1 */
      a0 = -ONE/h;
      a1 =  ONE/h;
    }
    N_VLinearSum(a0, HINT_YOLD(interp), a1, HINT_YNEW(interp), yout);
    break;

  case(2):    /* quadratic interpolant */
    if (d == 0) {
      a[0] = tau2;
      a[1] = ONE - tau2;
      a[2] = h*(tau2 + tau);
    } else if (d == 1) {
      a[0] = TWO*tau/h;
      a[1] = -TWO*tau/h;
      a[2] = (ONE + TWO*tau);
    } else {  /* d == 2 */
      a[0] = TWO/h/h;
      a[1] = -TWO/h/h;
      a[2] = TWO/h;
    }
    X[0] = HINT_YOLD(interp);
    X[1] = HINT_YNEW(interp);
    X[2] = HINT_FNEW(interp);
    retval = N_VLinearCombination(3, a, X, yout);
    if (retval != 0)  return(ARK_VECTOROP_ERR);
    break;

  case(3):    /* cubic interpolant */
    if (d == 0) {
      a[0] = THREE*tau2 + TWO*tau3;
      a[1] = ONE - THREE*tau2 - TWO*tau3;
      a[2] = h*(tau2 + tau3);
      a[3] = h*(tau + TWO*tau2 + tau3);
    } else if (d == 1) {
      a[0] = SIX*(tau + tau2)/h;
      a[1] = -SIX*(tau + tau2)/h;
      a[2] = TWO*tau + THREE*tau2;
      a[3] = ONE + FOUR*tau + THREE*tau2;
    } else if (d == 2) {
      a[0] = SIX*(ONE + TWO*tau)/h2;
      a[1] = -SIX*(ONE + TWO*tau)/h2;
      a[2] = (TWO + SIX*tau)/h;
      a[3] = (FOUR + SIX*tau)/h;
    } else {  /* d == 3 */
      a[0] = TWELVE/h3;
      a[1] = -TWELVE/h3;
      a[2] = SIX/h2;
      a[3] = SIX/h2;
    }
    X[0] = HINT_YOLD(interp);
    X[1] = HINT_YNEW(interp);
    X[2] = HINT_FOLD(interp);
    X[3] = HINT_FNEW(interp);
    retval = N_VLinearCombination(4, a, X, yout);
    if (retval != 0) return(ARK_VECTOROP_ERR);
   break;

  case(4):    /* quartic interpolant */

    /* first, evaluate cubic interpolant at tau=-1/3 */
    tval = -ONE/THREE;
    retval = arkInterpEvaluate(arkode_mem, interp, tval, 0, 3, yout);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* second, evaluate RHS at tau=-1/3, storing the result in fa */
    tval = HINT_TNEW(interp) - h/THREE;
    retval = ark_mem->step_fullrhs(ark_mem, tval, yout, HINT_FA(interp), 2);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* evaluate desired function */
    if (d == 0) {
      a[0] = -SIX*tau2 - RCONST(16.0)*tau3 - RCONST(9.0)*tau4;
      a[1] = ONE + SIX*tau2 + RCONST(16.0)*tau3 + RCONST(9.0)*tau4;
      a[2] = h*FOURTH*(-FIVE*tau2 - RCONST(14.0)*tau3 - RCONST(9.0)*tau4);
      a[3] = h*(tau + TWO*tau2 + tau3);
      a[4] = h*RCONST(27.0)*FOURTH*(-tau4 - TWO*tau3 - tau2);
    } else if (d == 1) {
      a[0] = (-TWELVE*tau - RCONST(48.0)*tau2 - RCONST(36.0)*tau3)/h;
      a[1] = (TWELVE*tau + RCONST(48.0)*tau2 + RCONST(36.0)*tau3)/h;
      a[2] = HALF*(-FIVE*tau - RCONST(21.0)*tau2 - RCONST(18.0)*tau3);
      a[3] = (ONE + FOUR*tau + THREE*tau2);
      a[4] = -RCONST(27.0)*HALF*(TWO*tau3 + THREE*tau2 + tau);
    } else if (d == 2) {
      a[0] = (-TWELVE - RCONST(96.0)*tau - RCONST(108.0)*tau2)/h2;
      a[1] = (TWELVE + RCONST(96.0)*tau + RCONST(108.0)*tau2)/h2;
      a[2] = (-FIVE*HALF - RCONST(21.0)*tau - RCONST(27.0)*tau2)/h;
      a[3] = (FOUR + SIX*tau)/h;
      a[4] = (-RCONST(27.0)*HALF - RCONST(81.0)*tau - RCONST(81.0)*tau2)/h;
    } else if (d == 3) {
      a[0] = (-RCONST(96.0) - RCONST(216.0)*tau)/h3;
      a[1] = (RCONST(96.0) + RCONST(216.0)*tau)/h3;
      a[2] = (-RCONST(21.0) - RCONST(54.0)*tau)/h2;
      a[3] = SIX/h2;
      a[4] = (-RCONST(81.0) - RCONST(162.0)*tau)/h2;
    } else {  /* d == 4 */
      a[0] = -RCONST(216.0)/h4;
      a[1] = RCONST(216.0)/h4;
      a[2] = -RCONST(54.0)/h3;
      a[3] = ZERO;
      a[4] = -RCONST(162.0)/h3;
    }
    X[0] = HINT_YOLD(interp);
    X[1] = HINT_YNEW(interp);
    X[2] = HINT_FOLD(interp);
    X[3] = HINT_FNEW(interp);
    X[4] = HINT_FA(interp);
    retval = N_VLinearCombination(5, a, X, yout);
    if (retval != 0) return(ARK_VECTOROP_ERR);
    break;

  case(5):    /* quintic interpolant */

    /* first, evaluate quartic interpolant at tau=-1/3 */
    tval = -ONE/THREE;
    retval = arkInterpEvaluate(arkode_mem, interp, tval, 0, 4, yout);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* second, evaluate RHS at tau=-1/3, storing the result in fa */
    tval = HINT_TNEW(interp) - h/THREE;
    retval = ark_mem->step_fullrhs(ark_mem, tval, yout, HINT_FA(interp), 2);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* third, evaluate quartic interpolant at tau=-2/3 */
    tval = -TWO/THREE;
    retval = arkInterpEvaluate(arkode_mem, interp, tval, 0, 4, yout);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* fourth, evaluate RHS at tau=-2/3, storing the result in fb */
    tval = HINT_TNEW(interp) - h*TWO/THREE;
    retval = ark_mem->step_fullrhs(ark_mem, tval, yout, HINT_FB(interp), 2);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* evaluate desired function */
    if (d == 0) {
      a[0] = RCONST(54.0)*tau5 + RCONST(135.0)*tau4 + RCONST(110.0)*tau3 + RCONST(30.0)*tau2;
      a[1] = ONE - a[0];
      a[2] = h/FOUR*(RCONST(27.0)*tau5 + RCONST(63.0)*tau4 + RCONST(49.0)*tau3 + RCONST(13.0)*tau2);
      a[3] = h/FOUR*(RCONST(27.0)*tau5 + RCONST(72.0)*tau4 + RCONST(67.0)*tau3 + RCONST(26.0)*tau2 + FOUR*tau);
      a[4] = h/FOUR*(RCONST(81.0)*tau5 + RCONST(189.0)*tau4 + RCONST(135.0)*tau3 + RCONST(27.0)*tau2);
      a[5] = h/FOUR*(RCONST(81.0)*tau5 + RCONST(216.0)*tau4 + RCONST(189.0)*tau3 + RCONST(54.0)*tau2);
    } else if (d == 1) {
      a[0] = (RCONST(270.0)*tau4 + RCONST(540.0)*tau3 + RCONST(330.0)*tau2 + RCONST(60.0)*tau)/h;
      a[1] = -a[0];
      a[2] = (RCONST(135.0)*tau4 + RCONST(252.0)*tau3 + RCONST(147.0)*tau2 + RCONST(26.0)*tau)/FOUR;
      a[3] = (RCONST(135.0)*tau4 + RCONST(288.0)*tau3 + RCONST(201.0)*tau2 + RCONST(52.0)*tau + FOUR)/FOUR;
      a[4] = (RCONST(405.0)*tau4 + RCONST(4.0)*189*tau3 + RCONST(405.0)*tau2 + RCONST(54.0)*tau)/FOUR;
      a[5] = (RCONST(405.0)*tau4 + RCONST(864.0)*tau3 + RCONST(567.0)*tau2 + RCONST(108.0)*tau)/FOUR;
    } else if (d == 2) {
      a[0] = (RCONST(1080.0)*tau3 + RCONST(1620.0)*tau2 + RCONST(660.0)*tau + RCONST(60.0))/h2;
      a[1] = -a[0];
      a[2] = (RCONST(270.0)*tau3 + RCONST(378.0)*tau2 + RCONST(147.0)*tau + RCONST(13.0))/(TWO*h);
      a[3] = (RCONST(270.0)*tau3 + RCONST(432.0)*tau2 + RCONST(201.0)*tau + RCONST(26.0))/(TWO*h);
      a[4] = (RCONST(810.0)*tau3 + RCONST(1134.0)*tau2 + RCONST(405.0)*tau + RCONST(27.0))/(TWO*h);
      a[5] = (RCONST(810.0)*tau3 + RCONST(1296.0)*tau2 + RCONST(567.0)*tau + RCONST(54.0))/(TWO*h);
    } else if (d == 3) {
      a[0] = (RCONST(3240.0)*tau2 + RCONST(3240.0)*tau + RCONST(660.0))/h3;
      a[1] = -a[0];
      a[2] = (RCONST(810.0)*tau2 + RCONST(756.0)*tau + RCONST(147.0))/(TWO*h2);
      a[3] = (RCONST(810.0)*tau2 + RCONST(864.0)*tau + RCONST(201.0))/(TWO*h2);
      a[4] = (RCONST(2430.0)*tau2 + RCONST(2268.0)*tau + RCONST(405.0))/(TWO*h2);
      a[5] = (RCONST(2430.0)*tau2 + RCONST(2592.0)*tau + RCONST(567.0))/(TWO*h2);
    } else if (d == 4) {
      a[0] = (RCONST(6480.0)*tau + RCONST(3240.0))/h4;
      a[1] = -a[0];
      a[2] = (RCONST(810.0)*tau + RCONST(378.0))/h3;
      a[3] = (RCONST(810.0)*tau + RCONST(432.0))/h3;
      a[4] = (RCONST(2430.0)*tau + RCONST(1134.0))/h3;
      a[5] = (RCONST(2430.0)*tau + RCONST(1296.0))/h3;
    } else {  /* d == 5 */
      a[0] = RCONST(6480.0)/h5;
      a[1] = -a[0];
      a[2] = RCONST(810.0)/h4;
      a[3] = a[2];
      a[4] = RCONST(2430.0)/h4;
      a[5] = a[4];
    }
    X[0] = HINT_YOLD(interp);
    X[1] = HINT_YNEW(interp);
    X[2] = HINT_FOLD(interp);
    X[3] = HINT_FNEW(interp);
    X[4] = HINT_FA(interp);
    X[5] = HINT_FB(interp);
    retval = N_VLinearCombination(6, a, X, yout);
    if (retval != 0) return(ARK_VECTOROP_ERR);
    break;

  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInterpEvaluate_Hermite",
                    "Illegal polynomial order");
    return (ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}




/*---------------------------------------------------------------
  Section III: Lagrange interpolation module implementation
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  arkInterpCreate_Lagrange:

  This routine creates an ARKInterp structure, through
  cloning an input template N_Vector.  This returns a non-NULL
  structure if no errors occurred, or a NULL value otherwise.
  ---------------------------------------------------------------*/
ARKInterp arkInterpCreate_Lagrange(void* arkode_mem, int degree)
{
  ARKInterp interp;
  ARKInterpContent_Lagrange content;
  ARKInterpOps ops;
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for valid degree */
  if (degree < 0 || degree > ARK_INTERP_MAX_DEGREE) return(NULL);

  /* allocate overall structure */
  interp = NULL;
  interp = (ARKInterp) malloc(sizeof *interp);
  if (interp == NULL)  return(NULL);

  /* allocate ops structure and set entries */
  ops = NULL;
  ops = (ARKInterpOps) malloc(sizeof *ops);
  if (ops == NULL) { free(interp); return(NULL); }
  ops->resize    = arkInterpResize_Lagrange;
  ops->free      = arkInterpFree_Lagrange;
  ops->print     = arkInterpPrintMem_Lagrange;
  ops->setdegree = arkInterpSetDegree_Lagrange;
  ops->init      = arkInterpInit_Lagrange;
  ops->update    = arkInterpUpdate_Lagrange;
  ops->evaluate  = arkInterpEvaluate_Lagrange;

  /* create content, and initialize everything to zero/NULL */
  content = NULL;
  content = (ARKInterpContent_Lagrange) malloc(sizeof *content);
  if (content == NULL) { free(ops); free(interp); return(NULL); }
  memset(content, 0, sizeof(struct _ARKInterpContent_Lagrange));

  /* attach ops and content structures to overall structure */
  interp->ops = ops;
  interp->content = content;

  /* fill content */

  /* maximum/current history length */
  content->nmax      = SUNMIN(degree+1, ARK_INTERP_MAX_DEGREE+1);  /* respect maximum possible */
  content->nmaxalloc = 0;
  content->nhist     = 0;

  /* initialize time/solution history arrays to NULL */
  content->thist = NULL;
  content->yhist = NULL;

  /* initial t roundoff value */
  content->tround = FUZZ_FACTOR*ark_mem->uround;

  /* update workspace sizes */
  ark_mem->lrw += content->nmax + 1;
  ark_mem->liw += content->nmax + 2;

  return(interp);
}


/*---------------------------------------------------------------
  arkInterpResize_Lagrange:

  This routine resizes the internal vectors.
  ---------------------------------------------------------------*/
int arkInterpResize_Lagrange(void* arkode_mem, ARKInterp I,
                             ARKVecResizeFn resize, void *resize_data,
                             sunindextype lrw_diff, sunindextype liw_diff,
                             N_Vector y0)
{
  int i;
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* resize vectors */
  if (I == NULL)  return(ARK_SUCCESS);
  if (LINT_YHIST(I) != NULL) {
    for (i=0; i<LINT_NMAXALLOC(I); i++) {
      if (!arkResizeVec(ark_mem, resize, resize_data, lrw_diff,
                        liw_diff, y0, &(LINT_YJ(I,i))))
        return(ARK_MEM_FAIL);
    }
  }

  /* reset active history length */
  LINT_NHIST(I) = 0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkInterpFree_Lagrange:

  This routine frees the Lagrange ARKInterp structure.
  ---------------------------------------------------------------*/
void arkInterpFree_Lagrange(void* arkode_mem, ARKInterp I)
{
  int i;
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return;
  ark_mem = (ARKodeMem) arkode_mem;

  /* if interpolation structure is NULL, just return */
  if (I == NULL) return;

  /* free content */
  if (I->content != NULL) {
    if (LINT_YHIST(I) != NULL) {
      for (i=0; i<LINT_NMAX(I); i++) {
        if (LINT_YJ(I,i) != NULL) {
          arkFreeVec(ark_mem, &(LINT_YJ(I,i)));
          LINT_YJ(I,i) = NULL;
        }
      }
      free(LINT_YHIST(I));
      LINT_YHIST(I) = NULL;
    }
    if (LINT_THIST(I) != NULL) {
      free(LINT_THIST(I));
      LINT_THIST(I) = NULL;
    }

    /* update work space sizes */
    ark_mem->lrw -= (LINT_NMAX(I) + 1);
    ark_mem->liw -= (LINT_NMAX(I) + 2);

    free(I->content);
    I->content = NULL;
  }

  /* free ops and interpolation structures */
  if (I->ops) {
    free(I->ops);
    I->ops = NULL;
  }
  free(I);
  I = NULL;

  return;
}


/*---------------------------------------------------------------
  arkInterpPrintMem_Lagrange

  This routine outputs the Lagrange temporal interpolation memory
  structure to a specified file pointer.
  ---------------------------------------------------------------*/
void arkInterpPrintMem_Lagrange(ARKInterp I, FILE *outfile)
{
  int i;
  if (I != NULL) {
    fprintf(outfile, "arkode_interp (Lagrange): nmax = %i\n", LINT_NMAX(I));
    fprintf(outfile, "arkode_interp (Lagrange): nhist = %i\n", LINT_NHIST(I));
    if (LINT_THIST(I) != NULL) {
      fprintf(outfile, "arkode_interp (Lagrange): thist =");
      for (i=0; i<LINT_NMAX(I); i++)
        fprintf(outfile, "  %"RSYM,LINT_TJ(I,i));
      fprintf(outfile, "\n");
    }
    if (LINT_YHIST(I) != NULL) {
      fprintf(outfile, "arkode_interp (Lagrange): yhist ptrs =");
      for (i=0; i<LINT_NMAX(I); i++)
        fprintf(outfile, "  %p",(void*) LINT_YJ(I,i));
      fprintf(outfile, "\n");
    }
#ifdef SUNDIALS_DEBUG_PRINTVEC
    if (LINT_YHIST(I) != NULL) {
      for (i=0; i<LINT_NMAX(I); i++) {
        fprintf(outfile, "arkode_interp (Lagrange): yhist[%i]:\n",i);
        N_VPrintFile(LINT_YJ(I,i), outfile);
      }
    }
#endif
  }
}


/*---------------------------------------------------------------
  arkInterpSetDegree_Lagrange

  This routine sets a supplied interpolation degree.  If the
  argument is positive, then we require that
      0 <= degree <= ARK_INTERP_MAX_DEGREE
  and use this value as the user-specified (or default) degree.

  If the argument is negative, then we assume that this has been
  called by a time-step module to limit the interpolant degree
  based on the temporal method order. In this case we set the
  Lagrange polynomial degree to be the minimum of (-degree),
  ARK_INTERP_MAX_DEGREE, and the previously-set value [i.e., in
  case the user has already specified use of a lower-degree
  polynomial].

  Return values:
    ARK_MEM_NULL -- if either arkode_mem or interp are NULL
    ARK_ILL_INPUT -- if the input is outside of allowable bounds
    ARK_INTERP_FAIL -- if the interpolation module has already
       been initialized,
    ARK_SUCCESS -- successful completion.
  ---------------------------------------------------------------*/
int arkInterpSetDegree_Lagrange(void* arkode_mem, ARKInterp I,
                                int degree)
{
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* if this degree is already stored, just return */
  if (abs(degree)+1 == LINT_NMAX(I))  return(ARK_SUCCESS);

  /* on positive degree, check for allowable value and overwrite stored degree */
  if (degree >= 0) {
    if (degree > ARK_INTERP_MAX_DEGREE) {
      arkProcessError(ark_mem, ARK_INTERP_FAIL, "ARKode",
                      "arkInterpSetDegree_Lagrange",
                      "Illegal degree specified.");
      return(ARK_ILL_INPUT);
    }

    LINT_NMAX(I) = degree+1;
    return(ARK_SUCCESS);
  }

  /* on negative degree, check for allowable value and update stored degree */
  degree = -degree;
  if (degree > ARK_INTERP_MAX_DEGREE) degree = ARK_INTERP_MAX_DEGREE;
  LINT_NMAX(I) = SUNMIN(LINT_NMAX(I), degree+1);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkInterpInit_Lagrange

  This routine performs the following steps:
  1. allocates any missing/needed (t,y) history arrays
  2. zeros out stored (t,y) history
  3. copies current (t,y) from main ARKode memory into history
  4. updates the 'active' history counter to 1
  ---------------------------------------------------------------*/
int arkInterpInit_Lagrange(void* arkode_mem, ARKInterp I,
                           realtype tnew)
{
  int i;
  ARKodeMem ark_mem;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* check if storage has increased since the last init */
  if (LINT_NMAX(I) > LINT_NMAXALLOC(I)) {
    if (LINT_THIST(I) != NULL) {
      free(LINT_THIST(I));
      LINT_THIST(I) = NULL;
    }
    if (LINT_YHIST(I) != NULL) {
      for (i=0; i<LINT_NMAX(I); i++) {
        if (LINT_YJ(I,i) != NULL) {
          arkFreeVec(ark_mem, &(LINT_YJ(I,i)));
          LINT_YJ(I,i) = NULL;
        }
      }
      free(LINT_YHIST(I));
      LINT_YHIST(I) = NULL;
    }
  }

  /* allocate storage for time and solution histories */
  if (LINT_THIST(I) == NULL) {
    LINT_THIST(I) = (realtype*) malloc(LINT_NMAX(I) * sizeof(realtype));
    if (LINT_THIST(I) == NULL) {
      arkInterpFree(ark_mem, I); return(ARK_MEM_FAIL);
    }
  }

  /* solution history allocation */
  if (LINT_YHIST(I) == NULL) {
    LINT_YHIST(I) = (N_Vector*) malloc(LINT_NMAX(I) * sizeof(N_Vector));
    if (LINT_YHIST(I) == NULL) {
      arkInterpFree(ark_mem, I); return(ARK_MEM_FAIL);
    }
    for (i=0; i<LINT_NMAX(I); i++) {
      LINT_YJ(I,i) = NULL;
      if (!arkAllocVec(ark_mem, ark_mem->yn, &(LINT_YJ(I,i)))) {
        arkInterpFree(ark_mem, I); return(ARK_MEM_FAIL);
      }
    }
  }

  /* update allocated size if necesary */
  if (LINT_NMAX(I) > LINT_NMAXALLOC(I))
    LINT_NMAXALLOC(I) = LINT_NMAX(I);

  /* zero out history (to be safe) */
  for (i=0; i<LINT_NMAXALLOC(I); i++)
    LINT_TJ(I,i) = RCONST(0.0);
  if (N_VConstVectorArray(LINT_NMAXALLOC(I), RCONST(0.0), LINT_YHIST(I)))
    return(ARK_VECTOROP_ERR);

  /* set current time and state as first entries of (t,y) history, update counter */
  LINT_TJ(I,0) = tnew;
  N_VScale(ONE, ark_mem->yn, LINT_YJ(I,0));
  LINT_NHIST(I) = 1;

  /* return with success */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkInterpUpdate_Lagrange

  If the current time is 'different enough' from the stored
  values, then this routine performs the following steps:
  1. shifts the t-history array values, and prepends the current
     time
  2. shifts the y-history pointers, and copies the current state
     into the first history vector
  Otherwise it just returns with success.
  ---------------------------------------------------------------*/
int arkInterpUpdate_Lagrange(void* arkode_mem, ARKInterp I, realtype tnew)
{
  int i;
  ARKodeMem ark_mem;
  realtype tdiff;
  N_Vector ytmp;
  int nhist, nmax;
  realtype *thist;
  N_Vector *yhist;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* set readability shortcuts */
  nhist = LINT_NHIST(I);
  nmax = LINT_NMAX(I);
  thist = LINT_THIST(I);
  yhist = LINT_YHIST(I);

  /* update t roundoff value */
  LINT_TROUND(I) = FUZZ_FACTOR*ark_mem->uround *
    (SUNRabs(ark_mem->tcur) + SUNRabs(ark_mem->h));

  /* determine if tnew differs sufficiently from stored values */
  tdiff = SUNRabs(tnew - thist[0]);
  for (i=1; i<nhist; i++)
    tdiff = SUNMIN(tdiff, SUNRabs(tnew - thist[i]));
  if (tdiff <= LINT_TROUND(I))  return(ARK_SUCCESS);

  /* shift (t,y) history arrays by one */
  ytmp = yhist[nmax-1];
  for (i=nmax-1; i>0; i--) {
    thist[i] = thist[i-1];
    yhist[i] = yhist[i-1];
  }
  yhist[0] = ytmp;

  /* copy tnew and ycur into first entry of history arrays */
  thist[0] = tnew;
  N_VScale(ONE, ark_mem->ycur, yhist[0]);

  /* update 'nhist' (first few steps) */
  LINT_NHIST(I) = nhist = SUNMIN(nhist+1, nmax);

  /* return with success */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkInterpEvaluate_Lagrange

  This routine evaluates a temporal interpolation/extrapolation
  based on the stored solution data in the interpolation structure.

  Derivatives have lower accuracy than the interpolant
  itself, losing one order per derivative.  This module can provide
  up to 3rd derivatives.

  The input 'tau' specifies the time at which to evaluate the
  Lagrange polynomial.  The formula for tau is defined using the
  most-recently-completed solution interval [t1,t0], and is
  given by:
               t = t0 + tau*(t0-t1),
  here t0 and t1 are the 2 most-recent entries in the 'thist'
  array within the interpolation structure.  Thus values
            -(nhist-1) <= tau < = 0
  provide interpolation, others result in extrapolation (assuming
  fixed step sizes, otherwise the stated lower bound is only
  approximate).
  ---------------------------------------------------------------*/
int arkInterpEvaluate_Lagrange(void* arkode_mem, ARKInterp I,
                               realtype tau, int deriv, int degree,
                               N_Vector yout)
{
  /* local variables */
  int q, retval, i, j;
  realtype tval;
  realtype a[6];
  N_Vector X[6];
  ARKodeMem ark_mem;
  int nhist;
  realtype *thist;
  N_Vector *yhist;

  /* access ARKodeMem structure */
  if (arkode_mem == NULL)  return(ARK_MEM_NULL);
  ark_mem = (ARKodeMem) arkode_mem;

  /* set readability shortcuts */
  nhist = LINT_NHIST(I);
  thist = LINT_THIST(I);
  yhist = LINT_YHIST(I);

  /* determine polynomial degree q */
  q = SUNMAX(degree, 0);    /* respect lower bound */
  q = SUNMIN(q, nhist-1);   /* respect max possible */

  /* error on illegal deriv */
  if ((deriv < 0) || (deriv > 3)) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode",
                    "arkInterpEvaluate_Lagrange",
                    "Requested illegal derivative.");
    return (ARK_ILL_INPUT);
  }

  /* if deriv is too high, just return zeros */
  if (deriv > q) {
    N_VConst(ZERO, yout);
    return(ARK_SUCCESS);
  }

  /* if constant interpolant is requested, just return ynew */
  if (q == 0) {
    N_VScale(ONE, yhist[0], yout);
    return(ARK_SUCCESS);
  }

  /* convert from tau back to t (both tnew and told are valid since q>0 => NHIST>1) */
  tval = thist[0] + tau*(thist[0]-thist[1]);

  /* linear interpolant */
  if (q == 1) {
    if (deriv == 0) {
      a[0] = LBasis(I,0,tval);
      a[1] = LBasis(I,1,tval);
    } else {  /* deriv == 1 */
      a[0] = LBasisD(I,0,tval);
      a[1] = LBasisD(I,1,tval);
    }
    N_VLinearSum(a[0], yhist[0], a[1], yhist[1], yout);
    return(ARK_SUCCESS);
  }

  /* higher-degree interpolant */
  /*    initialize arguments for N_VLinearCombination */
  for (i=0; i<q+1; i++) {
    a[i] = ZERO;
    X[i] = yhist[i];
  }

  /*    construct linear combination coefficients based on derivative requested */
  switch (deriv) {

  case(0):    /* p(t) */
    for (j=0; j<q+1; j++)  a[j] = LBasis(I,j,tval);
    break;

  case(1):    /* p'(t) */
    for (j=0; j<q+1; j++)  a[j] = LBasisD(I,j,tval);
    break;

  case(2):    /* p''(t) */
    for (j=0; j<q+1; j++)  a[j] = LBasisD2(I,j,tval);
    break;

  case(3):    /* p'''(t) */
    for (j=0; j<q+1; j++)  a[j] = LBasisD3(I,j,tval);
    break;
  }

  /*    call N_VLinearCombination to evaluate the result, and return */
  retval = N_VLinearCombination(q+1, a, X, yout);
  if (retval != 0)  return(ARK_VECTOROP_ERR);

  return(ARK_SUCCESS);
}


/* Lagrange utility routines (basis functions and their derivatives) */
realtype LBasis(ARKInterp I, int j, realtype t)
{
  int k;
  realtype p = ONE;
  for (k=0; k<LINT_NHIST(I); k++) {
    if (k == j) continue;
    p *= (t-LINT_TJ(I,k))/(LINT_TJ(I,j)-LINT_TJ(I,k));
  }
  return(p);
}


realtype LBasisD(ARKInterp I, int j, realtype t)
{
  int i, k;
  realtype p, q;
  p = ZERO;
  for (i=0; i<LINT_NHIST(I); i++) {
    if (i == j) continue;
    q = ONE;
    for (k=0; k<LINT_NHIST(I); k++) {
      if (k == j) continue;
      if (k == i) continue;
      q *= (t-LINT_TJ(I,k))/(LINT_TJ(I,j)-LINT_TJ(I,k));
    }
    p += q/(LINT_TJ(I,j)-LINT_TJ(I,i));
  }

  return(p);
}


realtype LBasisD2(ARKInterp I, int j, realtype t)
{
  int i, k, l;
  realtype p, q, r;
  p = ZERO;
  for (l=0; l<LINT_NHIST(I); l++) {
    if (l == j) continue;
    q = ZERO;
    for (i=0; i<LINT_NHIST(I); i++) {
      if (i == j) continue;
      if (i == l) continue;
      r = ONE;
      for (k=0; k<LINT_NHIST(I); k++) {
        if (k == j) continue;
        if (k == i) continue;
        if (k == l) continue;
        r *= (t-LINT_TJ(I,k))/(LINT_TJ(I,j)-LINT_TJ(I,k));
      }
      q += r/(LINT_TJ(I,j)-LINT_TJ(I,i));
    }
    p += q/(LINT_TJ(I,j)-LINT_TJ(I,l));
  }

  return(p);
}


realtype LBasisD3(ARKInterp I, int j, realtype t)
{
  int i, k, l, m;
  realtype p, q, r, s;
  p = ZERO;
  for (m=0; m<LINT_NHIST(I); m++) {
    if (m == j) continue;
    q = ZERO;
    for (l=0; l<LINT_NHIST(I); l++) {
      if (l == j) continue;
      if (l == m) continue;
      r = ZERO;
      for (i=0; i<LINT_NHIST(I); i++) {
        if (i == j) continue;
        if (i == m) continue;
        if (i == l) continue;
        s = ONE;
        for (k=0; k<LINT_NHIST(I); k++) {
          if (k == j) continue;
          if (k == m) continue;
          if (k == l) continue;
          if (k == i) continue;
          s *= (t-LINT_TJ(I,k))/(LINT_TJ(I,j)-LINT_TJ(I,k));
        }
        r += s/(LINT_TJ(I,j)-LINT_TJ(I,i));
      }
      q += r/(LINT_TJ(I,j)-LINT_TJ(I,l));
    }
    p += q/(LINT_TJ(I,j)-LINT_TJ(I,m));
  }

  return(p);
}


/*===============================================================
  EOF
  ===============================================================*/
