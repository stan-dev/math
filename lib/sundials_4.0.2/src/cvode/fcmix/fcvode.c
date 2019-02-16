/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *     Alan C. Hindmarsh, Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the CVODE package.  See fcvode.h for usage.
 * NOTE: some routines are necessarily stored elsewhere to avoid
 * linking problems.  Therefore, see the othe C files in this folder
 * for all the options available.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fcvode.h"                    /* actual function names, prototypes, global vars.*/ 
#include "cvode_impl.h"                /* definition of CVodeMem type                    */
#include <sundials/sundials_matrix.h>
#include <cvode/cvode_ls.h>
#include <cvode/cvode_diag.h>


/***************************************************************************/

/* Definitions for global variables shared amongst various routines */

void *CV_cvodemem;
long int *CV_iout;
realtype *CV_rout;
int CV_nrtfn;
int CV_ls;

/***************************************************************************/

/* private constant(s) */
#define ZERO RCONST(0.0)

/***************************************************************************/

/* Prototypes of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_FUN(realtype*,     /* T    */
                      realtype*,     /* Y    */
                      realtype*,     /* YDOT */
                      long int*,     /* IPAR */
                      realtype*,     /* RPAR */
                      int*);         /* IER  */
#ifdef __cplusplus
}
#endif

/**************************************************************************/

void FCV_MALLOC(realtype *t0, realtype *y0,
                int *meth, int *iatol,
                realtype *rtol, realtype *atol,
                long int *iout, realtype *rout,
                long int *ipar, realtype *rpar,
                int *ier)
{
  int lmm;
  N_Vector Vatol;
  FCVUserData CV_userdata;

  *ier = 0;

  /* Check for required vector operations */
  if(F2C_CVODE_vec->ops->nvgetarraypointer == NULL ||
     F2C_CVODE_vec->ops->nvsetarraypointer == NULL) {
    *ier = -1;
    fprintf(stderr, "A required vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize all pointers to NULL */
  CV_cvodemem = NULL;
  Vatol = NULL;
  FCVNullNonlinSol();

  /* initialize global constants to disable each option */
  CV_nrtfn = 0;
  CV_ls = -1;
  
  /* Create CVODE object */
  lmm = (*meth == 1) ? CV_ADAMS : CV_BDF;
  CV_cvodemem = CVodeCreate(lmm);
  if (CV_cvodemem == NULL) {
    *ier = -1;
    return;
  }
  
  /* Set and attach user data */
  CV_userdata = NULL;
  CV_userdata = (FCVUserData) malloc(sizeof *CV_userdata);
  if (CV_userdata == NULL) {
    *ier = -1;
    return;
  }
  CV_userdata->rpar = rpar;
  CV_userdata->ipar = ipar;

  *ier = CVodeSetUserData(CV_cvodemem, CV_userdata);
  if(*ier != CV_SUCCESS) {
    free(CV_userdata); CV_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Set data in F2C_CVODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_CVODE_vec);

  /* Call CVodeInit */
  *ier = CVodeInit(CV_cvodemem, FCVf, *t0, F2C_CVODE_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_CVODE_vec);

  /* On failure, exit */
  if(*ier != CV_SUCCESS) {
    free(CV_userdata); CV_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Set tolerances */
  switch (*iatol) {
  case 1:
    *ier = CVodeSStolerances(CV_cvodemem, *rtol, *atol); 
    break;
  case 2:
    Vatol = NULL;
    Vatol = N_VCloneEmpty(F2C_CVODE_vec);
    if (Vatol == NULL) {
      free(CV_userdata); CV_userdata = NULL;
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    *ier = CVodeSVtolerances(CV_cvodemem, *rtol, Vatol);
    N_VDestroy(Vatol);
    break;
  }

  /* On failure, exit */
  if(*ier != CV_SUCCESS) {
    free(CV_userdata); CV_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Grab optional output arrays and store them in global variables */
  CV_iout = iout;
  CV_rout = rout;

  /* Store the unit roundoff in rout for user access */
  CV_rout[5] = UNIT_ROUNDOFF;

  return;
}

/***************************************************************************/

void FCV_REINIT(realtype *t0, realtype *y0, 
                int *iatol, realtype *rtol, realtype *atol, 
                int *ier)
{
  N_Vector Vatol;

  *ier = 0;

  /* Initialize all pointers to NULL */
  Vatol = NULL;

  /* Set data in F2C_CVODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_CVODE_vec);

  /* Call CVReInit */
  *ier = CVodeReInit(CV_cvodemem, *t0, F2C_CVODE_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_CVODE_vec);

  /* On failure, exit */
  if (*ier != CV_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Set tolerances */
  switch (*iatol) {
  case 1:
    *ier = CVodeSStolerances(CV_cvodemem, *rtol, *atol); 
    break;
  case 2:
    Vatol = NULL;
    Vatol = N_VCloneEmpty(F2C_CVODE_vec);
    if (Vatol == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    *ier = CVodeSVtolerances(CV_cvodemem, *rtol, Vatol);
    N_VDestroy(Vatol);
    break;
  }

  /* On failure, exit */
  if (*ier != CV_SUCCESS) {
    *ier = -1;
    return;
  }

  return;
}

/***************************************************************************/

void FCV_SETIIN(char key_name[], long int *ival, int *ier)
{
  if (!strncmp(key_name,"MAX_ORD",7))
    *ier = CVodeSetMaxOrd(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NSTEPS",10))
    *ier = CVodeSetMaxNumSteps(CV_cvodemem, (long int) *ival);
  else if (!strncmp(key_name,"MAX_ERRFAIL",11))
    *ier = CVodeSetMaxErrTestFails(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NITERS",10))
    *ier = CVodeSetMaxNonlinIters(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"MAX_CONVFAIL",12))
    *ier = CVodeSetMaxConvFails(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"HNIL_WARNS",10))
    *ier = CVodeSetMaxHnilWarns(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"STAB_LIM",8))
    *ier = CVodeSetStabLimDet(CV_cvodemem, (booleantype) *ival);
  else {
    *ier = -99;
    fprintf(stderr, "FCVSETIIN: Unrecognized key.\n\n");
  }

}

/***************************************************************************/

void FCV_SETRIN(char key_name[], realtype *rval, int *ier)
{
  if (!strncmp(key_name,"INIT_STEP",9))
    *ier = CVodeSetInitStep(CV_cvodemem, *rval);
  else if (!strncmp(key_name,"MAX_STEP",8))
    *ier = CVodeSetMaxStep(CV_cvodemem, *rval);
  else if (!strncmp(key_name,"MIN_STEP",8))
    *ier = CVodeSetMinStep(CV_cvodemem, *rval);
  else if (!strncmp(key_name,"STOP_TIME",9))
    *ier = CVodeSetStopTime(CV_cvodemem, *rval);
  else if (!strncmp(key_name,"NLCONV_COEF",11))
    *ier = CVodeSetNonlinConvCoef(CV_cvodemem, *rval);
  else {
    *ier = -99;
    fprintf(stderr, "FCVSETRIN: Unrecognized key.\n\n");
  }

}

/***************************************************************************/

void FCV_SETVIN(char key_name[], realtype *vval, int *ier)
{
  N_Vector Vec;

  *ier = 0;

  if (!strncmp(key_name,"CONSTR_VEC",10)) {
    Vec = NULL;
    Vec = N_VCloneEmpty(F2C_CVODE_vec);
    if (Vec == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(vval, Vec);
    CVodeSetConstraints(CV_cvodemem, Vec);
    N_VDestroy(Vec);
  }
  else {
    *ier = -99;
    fprintf(stderr, "FCVSETVIN: Unrecognized key. \n\n");
  }

}

/***************************************************************************/

void FCV_LSINIT(int *ier) {
  if ( (CV_cvodemem == NULL) || (F2C_CVODE_linsol == NULL) ) {
    *ier = -1;
    return;
  }
  *ier = CVodeSetLinearSolver(CV_cvodemem, F2C_CVODE_linsol, 
                              F2C_CVODE_matrix);
  CV_ls = CV_LS_STD;
  return;
}

/***************************************************************************/

/* ---DEPRECATED--- */
void FCV_DLSINIT(int *ier)
{ FCV_LSINIT(ier); }

/***************************************************************************/

/* ---DEPRECATED--- */
void FCV_SPILSINIT(int *ier)
{ FCV_LSINIT(ier); }

/*=============================================================*/

/* ---DEPRECATED--- */
void FCV_SPILSSETEPSLIN(realtype *eplifac, int *ier)
{ FCV_LSSETEPSLIN(eplifac, ier); }

void FCV_LSSETEPSLIN(realtype *eplifac, int *ier)
{ *ier = CVodeSetEpsLin(CV_cvodemem, *eplifac); }

/***************************************************************************/

void FCV_DIAG(int *ier)
{
  if (CV_cvodemem == NULL) {
    *ier = -1;
    return;
  }
  *ier = CVDiag(CV_cvodemem);
  CV_ls = CV_LS_DIAG;
  return;
}

/***************************************************************************/

void FCV_CVODE(realtype *tout, realtype *t, realtype *y, int *itask, int *ier)
{
  /* 
     tout          is the t value where output is desired
     F2C_CVODE_vec is the N_Vector containing the solution on return
     t             is the returned independent variable value
     itask         is the task indicator (1 = CV_NORMAL, 2 = CV_ONE_STEP, 
                                          3 = CV_NORMAL_TSTOP, 4 = CV_ONE_STEP_TSTOP) 
  */

  int qu, qcur;

  N_VSetArrayPointer(y, F2C_CVODE_vec);

  *ier = CVode(CV_cvodemem, *tout, F2C_CVODE_vec, t, *itask);

  N_VSetArrayPointer(NULL, F2C_CVODE_vec);

  /* Load optional outputs in iout & rout */
  CVodeGetWorkSpace(CV_cvodemem,
                    &CV_iout[0],                          /* LENRW   */
                    &CV_iout[1]);                         /* LENIW   */
  CVodeGetIntegratorStats(CV_cvodemem, 
                          &CV_iout[2],                    /* NST     */
                          &CV_iout[3],                    /* NFE     */ 
                          &CV_iout[7],                    /* NSETUPS */ 
                          &CV_iout[4],                    /* NETF    */ 
                          &qu,                            /* QU      */
                          &qcur,                          /* QCUR    */
                          &CV_rout[0],                    /* H0U     */
                          &CV_rout[1],                    /* HU      */ 
                          &CV_rout[2],                    /* HCUR    */ 
                          &CV_rout[3]);                   /* TCUR    */ 
  CV_iout[8] = (long int) qu;
  CV_iout[9] = (long int) qcur;
  CVodeGetTolScaleFactor(CV_cvodemem, 
                         &CV_rout[4]);                    /* TOLSFAC */
  CVodeGetNonlinSolvStats(CV_cvodemem,
                          &CV_iout[6],                    /* NNI     */
                          &CV_iout[5]);                   /* NCFN    */
  CVodeGetNumStabLimOrderReds(CV_cvodemem, &CV_iout[10]); /* NOR     */
  
  /* Root finding is on */
  if (CV_nrtfn != 0)
    CVodeGetNumGEvals(CV_cvodemem, &CV_iout[11]);         /* NGE     */
  
  switch(CV_ls) {
  case CV_LS_STD:
    CVodeGetLinWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);   /* LENRWLS,LENIWLS */
    CVodeGetLastLinFlag(CV_cvodemem, &CV_iout[14]);                  /* LSTF */
    CVodeGetNumLinRhsEvals(CV_cvodemem, &CV_iout[15]);               /* NFELS */
    CVodeGetNumJacEvals(CV_cvodemem, &CV_iout[16]);                  /* NJE */
    CVodeGetNumJTSetupEvals(CV_cvodemem, &CV_iout[17]);              /* NJTS */
    CVodeGetNumJtimesEvals(CV_cvodemem, &CV_iout[18]);               /* NJTV */
    CVodeGetNumPrecEvals(CV_cvodemem, &CV_iout[19]);                 /* NPE */
    CVodeGetNumPrecSolves(CV_cvodemem, &CV_iout[20]);                /* NPS */
    CVodeGetNumLinIters(CV_cvodemem, &CV_iout[21]);                  /* NLI */
    CVodeGetNumLinConvFails(CV_cvodemem, &CV_iout[22]);              /* NCFL */
    break;
  case CV_LS_DIAG:
    CVDiagGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);  /* LENRWLS,LENIWLS */
    CVDiagGetLastFlag(CV_cvodemem, &CV_iout[14]);                 /* LSTF */
    CVDiagGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);              /* NFELS */
  }
}

/***************************************************************************/

void FCV_DKY (realtype *t, int *k, realtype *dky, int *ier)
{
  /* 
     t             is the t value where output is desired
     k             is the derivative order
     F2C_CVODE_vec is the N_Vector containing the solution derivative on return 
  */

  realtype *f2c_data = N_VGetArrayPointer(F2C_CVODE_vec);
  N_VSetArrayPointer(dky, F2C_CVODE_vec);

  *ier = 0;
  *ier = CVodeGetDky(CV_cvodemem, *t, *k, F2C_CVODE_vec);

  N_VSetArrayPointer(f2c_data, F2C_CVODE_vec);

}

/*************************************************/

void FCV_GETERRWEIGHTS(realtype *eweight, int *ier)
{
  /* Attach user data to vector */
  realtype *f2c_data = N_VGetArrayPointer(F2C_CVODE_vec);
  N_VSetArrayPointer(eweight, F2C_CVODE_vec);

  *ier = 0;
  *ier = CVodeGetErrWeights(CV_cvodemem, F2C_CVODE_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(f2c_data, F2C_CVODE_vec);

  return;
}

/*************************************************/

void FCV_GETESTLOCALERR(realtype *ele, int *ier)
{
  /* Attach user data to vector */
  realtype *f2c_data = N_VGetArrayPointer(F2C_CVODE_vec);
  N_VSetArrayPointer(ele, F2C_CVODE_vec);

  *ier = 0;
  *ier = CVodeGetEstLocalErrors(CV_cvodemem, F2C_CVODE_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(f2c_data, F2C_CVODE_vec);

  return;
}

/***************************************************************************/

void FCV_FREE ()
{
  CVodeMem cv_mem;

  cv_mem = (CVodeMem) CV_cvodemem;

  if (cv_mem->cv_lfree)
    cv_mem->cv_lfree(cv_mem);
  cv_mem->cv_lmem = NULL;
  
  free(cv_mem->cv_user_data); cv_mem->cv_user_data = NULL;

  CVodeFree(&CV_cvodemem);

  N_VSetArrayPointer(NULL, F2C_CVODE_vec);
  N_VDestroy(F2C_CVODE_vec);
  if (F2C_CVODE_matrix)
    SUNMatDestroy(F2C_CVODE_matrix);
  if (F2C_CVODE_linsol)
    SUNLinSolFree(F2C_CVODE_linsol);
  /* already freed by CVodeFree */
  if (F2C_CVODE_nonlinsol)
    F2C_CVODE_nonlinsol = NULL;
  return;
}

/***************************************************************************/

/* 
 * C function CVf to interface between CVODE and a Fortran subroutine FCVFUN.
 * Addresses of t, y, and ydot are passed to CVFUN, using the
 * routine N_VGetArrayPointer from the NVECTOR module.
 * Auxiliary data is assumed to be communicated by Common. 
 */

int FCVf(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  int ier;
  realtype *ydata, *dydata;
  FCVUserData CV_userdata;

  ydata  = N_VGetArrayPointer(y);
  dydata = N_VGetArrayPointer(ydot);

  CV_userdata = (FCVUserData) user_data;

  FCV_FUN(&t, ydata, dydata, CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}

/* Fortran interface to C routine CVodeSetNonlinearSolver; see 
   fcvode.h for further details */
void FCV_NLSINIT(int *ier) {
  if ( (CV_cvodemem == NULL) || (F2C_CVODE_nonlinsol == NULL) ) {
    *ier = -1;
    return;
  }
  if (((CVodeMem) CV_cvodemem)->cv_lsolve == NULL) {
    FCVNullMatrix();
    FCVNullLinsol();
  }

  *ier = CVodeSetNonlinearSolver(CV_cvodemem, F2C_CVODE_nonlinsol);
  return;
}
