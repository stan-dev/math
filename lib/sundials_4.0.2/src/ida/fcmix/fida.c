/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Aaron Collier and Radu Serban @ LLNL
 *-----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *-----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the IDA package. See fida.h for usage.
 * NOTE: Some routines are necessarily stored elsewhere to avoid
 * linking problems.
 *-----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fida.h"            /* function names, prototypes, global variables   */
#include "ida_impl.h"        /* definition of IDAMem type                      */
#include <ida/ida_ls.h>      /* prototypes for IDALS interface routines        */

/*************************************************/

/* Definitions for global variables shared amongst various routines */

N_Vector F2C_IDA_ypvec, F2C_IDA_ewtvec;

void *IDA_idamem;
long int *IDA_iout;
realtype *IDA_rout;
int IDA_nrtfn;

/*************************************************/

/* private constant(s) */
#define ZERO RCONST(0.0)

/*************************************************/

/* Prototype of user-supplied Fortran routine (IDAResFn) */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_RESFUN(realtype*,    /* T    */
                          realtype*,    /* Y    */
                          realtype*,    /* YP   */
                          realtype*,    /* R    */
                          long int*,    /* IPAR */
                          realtype*,    /* RPAR */
                          int*);        /* IER  */

#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_MALLOC(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 long int *iout, realtype *rout, 
                 long int *ipar, realtype *rpar, int *ier)
{
  N_Vector Vatol;
  FIDAUserData IDA_userdata;

  *ier = 0;

  /* Check for required vector operations */
  if ((F2C_IDA_vec->ops->nvgetarraypointer == NULL) ||
      (F2C_IDA_vec->ops->nvsetarraypointer == NULL)) {
    *ier = -1;
    fprintf(stderr, "A required vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize all pointers to NULL */
  IDA_idamem = NULL;
  Vatol = NULL;
  F2C_IDA_ypvec = F2C_IDA_ewtvec = NULL;
  FIDANullNonlinSol();

  /* Create IDA object */
  IDA_idamem = IDACreate();
  if (IDA_idamem == NULL) {
    *ier = -1;
    return;
  }

  /* Set and attach user data */
  IDA_userdata = NULL;
  IDA_userdata = (FIDAUserData) malloc(sizeof *IDA_userdata);
  if (IDA_userdata == NULL) {
    *ier = -1;
    return;
  }
  IDA_userdata->rpar = rpar;
  IDA_userdata->ipar = ipar;

  *ier = IDASetUserData(IDA_idamem, IDA_userdata);
  if(*ier != IDA_SUCCESS) {
    free(IDA_userdata); IDA_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Attach user's yy0 to F2C_IDA_vec */
  N_VSetArrayPointer(yy0, F2C_IDA_vec);

  /* Create F2C_IDA_ypvec and attach user's yp0 to it */
  F2C_IDA_ypvec = NULL;
  F2C_IDA_ypvec = N_VCloneEmpty(F2C_IDA_vec);
  if (F2C_IDA_ypvec == NULL) {
    free(IDA_userdata); IDA_userdata = NULL;
    *ier = -1;
  }
  N_VSetArrayPointer(yp0, F2C_IDA_ypvec);

  /* Call IDAInit */
  *ier = IDAInit(IDA_idamem, FIDAresfn, *t0, F2C_IDA_vec, F2C_IDA_ypvec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  /* On failure, clean-up and exit */
  if (*ier != IDA_SUCCESS) {
    N_VDestroy(F2C_IDA_ypvec);
    free(IDA_userdata); IDA_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Set tolerances */
  switch (*iatol) {
  case 1:
    *ier = IDASStolerances(IDA_idamem, *rtol, *atol);
    break;
  case 2:
    Vatol = NULL;
    Vatol= N_VCloneEmpty(F2C_IDA_vec);
    if (Vatol == NULL) {
      free(IDA_userdata); IDA_userdata = NULL;
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    *ier = IDASVtolerances(IDA_idamem, *rtol, Vatol);
    N_VDestroy(Vatol);
    break;
  }

  /* On failure, clean-up and exit */
  if (*ier != IDA_SUCCESS) {
    free(IDA_userdata); IDA_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Grab optional output arrays and store them in global variables */
  IDA_iout = iout;
  IDA_rout = rout;

  /* Store the unit roundoff in rout for user access */
  IDA_rout[5] = UNIT_ROUNDOFF;

  /* Set F2C_IDA_ewtvec on NULL */
  F2C_IDA_ewtvec = NULL;

  return;
}

/*************************************************/

void FIDA_REINIT(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 int *ier)
{
  N_Vector Vatol;

  *ier = 0;

  /* Initialize all pointers to NULL */
  Vatol = NULL;

  /* Attach user's yy0 to F2C_IDA_vec */
  N_VSetArrayPointer(yy0, F2C_IDA_vec);

  /* Attach user's yp0 to F2C_IDA_ypvec */
  N_VSetArrayPointer(yp0, F2C_IDA_ypvec);

  /* Call IDAReInit */
  *ier = IDAReInit(IDA_idamem, *t0, F2C_IDA_vec, F2C_IDA_ypvec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  /* On failure, exit */
  if (*ier != IDA_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Set tolerances */
  switch (*iatol) {
  case 1:
    *ier = IDASStolerances(IDA_idamem, *rtol, *atol);
    break;
  case 2:
    Vatol = NULL;
    Vatol= N_VCloneEmpty(F2C_IDA_vec);
    if (Vatol == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    *ier = IDASVtolerances(IDA_idamem, *rtol, Vatol);
    N_VDestroy(Vatol);
    break;
  }

  /* On failure, exit */
  if (*ier != IDA_SUCCESS) {
    *ier = -1;
    return;
  }

  return;
}

/*************************************************/

void FIDA_SETIIN(char key_name[], long int *ival, int *ier)
{
  if (!strncmp(key_name,"MAX_ORD",7))
    *ier = IDASetMaxOrd(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NSTEPS",10))
    *ier = IDASetMaxNumSteps(IDA_idamem, (long int) *ival);
  else if (!strncmp(key_name,"MAX_ERRFAIL",11))
    *ier = IDASetMaxErrTestFails(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NITERS",10))
    *ier = IDASetMaxNonlinIters(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_CONVFAIL",12))
    *ier = IDASetMaxConvFails(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"SUPPRESS_ALG",12))
    *ier = IDASetSuppressAlg(IDA_idamem, (booleantype) *ival);
  else if (!strncmp(key_name,"MAX_NSTEPS_IC",13))
    *ier = IDASetMaxNumStepsIC(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NITERS_IC",13)) 
    *ier = IDASetMaxNumItersIC(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NJE_IC",10))
    *ier = IDASetMaxNumJacsIC(IDA_idamem, (int) *ival);
  else if (!strncmp(key_name,"LS_OFF_IC",9))
    *ier = IDASetLineSearchOffIC(IDA_idamem, (booleantype) *ival);
  else {
    *ier = -99;
    fprintf(stderr, "FIDASETIIN: Unrecognized key.\n\n");
  }

}

/***************************************************************************/

void FIDA_SETRIN(char key_name[], realtype *rval, int *ier)
{

  if (!strncmp(key_name,"INIT_STEP",9))
    *ier = IDASetInitStep(IDA_idamem, *rval);
  else if (!strncmp(key_name,"MAX_STEP",8))
    *ier = IDASetMaxStep(IDA_idamem, *rval);
  else if (!strncmp(key_name,"STOP_TIME",9))
    *ier = IDASetStopTime(IDA_idamem, *rval);
  else if (!strncmp(key_name,"NLCONV_COEF_IC",14))
    *ier = IDASetNonlinConvCoefIC(IDA_idamem, *rval);
  else if (!strncmp(key_name,"NLCONV_COEF",11))
    *ier = IDASetNonlinConvCoef(IDA_idamem, *rval);
  else if (!strncmp(key_name,"STEP_TOL_IC",11))
    *ier = IDASetStepToleranceIC(IDA_idamem, *rval);
  else {
    *ier = -99;
    fprintf(stderr, "FIDASETRIN: Unrecognized key.\n\n");
  }

}

/*************************************************/

void FIDA_SETVIN(char key_name[], realtype *vval, int *ier)
{
  N_Vector Vec;

  *ier = 0;

  if (!strncmp(key_name,"ID_VEC",6)) {
    Vec = NULL;
    Vec = N_VCloneEmpty(F2C_IDA_vec);
    if (Vec == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(vval, Vec);
    IDASetId(IDA_idamem, Vec);
    N_VDestroy(Vec);
  } else if (!strncmp(key_name,"CONSTR_VEC",10)) {
    Vec = NULL;
    Vec = N_VCloneEmpty(F2C_IDA_vec);
    if (Vec == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(vval, Vec);
    IDASetConstraints(IDA_idamem, Vec);
    N_VDestroy(Vec);
  } else {
    *ier = -99;
    fprintf(stderr, "FIDASETVIN: Unrecognized key.\n\n");
  }

}

/*************************************************/

void FIDA_TOLREINIT(int *iatol, realtype *rtol, realtype *atol, int *ier)
{
  N_Vector Vatol=NULL;

  *ier = 0;

  if (*iatol == 1) {
    *ier = IDASStolerances(IDA_idamem, *rtol, *atol);
  } else {
    Vatol = NULL;
    Vatol = N_VCloneEmpty(F2C_IDA_vec);
    if (Vatol == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    *ier = IDASVtolerances(IDA_idamem, *rtol, Vatol);
    N_VDestroy(Vatol);
  }

  return;
}

/*************************************************/

void FIDA_CALCIC(int *icopt, realtype *tout1, int *ier)
{
  *ier = 0;
  *ier = IDACalcIC(IDA_idamem, *icopt, *tout1);
  return;
}

/*************************************************/

/* Fortran interface to C routine IDASetLinearSolver; see 
   fida.h for further details */
void FIDA_LSINIT(int *ier) {
  if ( (IDA_idamem == NULL) || (F2C_IDA_linsol == NULL) ) {
    *ier = -1;
    return;
  }
  *ier = IDASetLinearSolver(IDA_idamem, F2C_IDA_linsol, 
                            F2C_IDA_matrix);
  return;
}


/*************************************************/

/*** DEPRECATED ***/
void FIDA_DLSINIT(int *ier)
{ FIDA_LSINIT(ier); }

/*************************************************/

/*** DEPRECATED ***/
void FIDA_SPILSINIT(int *ier) {
  FIDANullMatrix();
  FIDA_LSINIT(ier);
}

/*************************************************/

/* Fortran interfaces to C "set" routines for the IDALS solver; 
   see fida.h for further details */
void FIDA_LSSETEPSLIN(realtype *eplifac, int *ier) {
  *ier = IDASetEpsLin(IDA_idamem, *eplifac);
  return;
}

void FIDA_LSSETINCREMENTFACTOR(realtype *dqincfac, int *ier) {
  *ier = IDASetIncrementFactor(IDA_idamem, *dqincfac);
  return;
}

/*** DEPRECATED ***/
void FIDA_SPILSSETEPSLIN(realtype *eplifac, int *ier)
{ FIDA_LSSETEPSLIN(eplifac, ier); }

/*** DEPRECATED ***/
void FIDA_SPILSSETINCREMENTFACTOR(realtype *dqincfac, int *ier) 
{ FIDA_LSSETINCREMENTFACTOR(dqincfac, ier); }


/*************************************************/

void FIDA_SOLVE(realtype *tout, realtype *tret, realtype *yret,
		realtype *ypret, int *itask, int *ier)
{
  int klast, kcur;

  *ier = 0;

  /* Attach user data to vectors */
  N_VSetArrayPointer(yret, F2C_IDA_vec);
  N_VSetArrayPointer(ypret, F2C_IDA_ypvec);

  *ier = IDASolve(IDA_idamem, *tout, tret, F2C_IDA_vec, F2C_IDA_ypvec, *itask);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);

  /* Set optional outputs */

  IDAGetWorkSpace(IDA_idamem,
                  &IDA_iout[0],                 /* LENRW */
                  &IDA_iout[1]);                /* LENIW */

  IDAGetIntegratorStats(IDA_idamem,
                        &IDA_iout[2],           /* NST */
                        &IDA_iout[3],           /* NRE */
                        &IDA_iout[7],           /* NSETUPS */
                        &IDA_iout[4],           /* NETF */
                        &klast,                 /* KLAST */
                        &kcur,                  /* KCUR */
                        &IDA_rout[0],           /* HINUSED */
                        &IDA_rout[1],           /* HLAST */
                        &IDA_rout[2],           /* HCUR */
                        &IDA_rout[3]);          /* TCUR */
  IDA_iout[8] = (long int) klast;
  IDA_iout[9] = (long int) kcur;
  IDAGetNonlinSolvStats(IDA_idamem,
                        &IDA_iout[6],           /* NNI */
                        &IDA_iout[5]);          /* NCFN */
  IDAGetNumBacktrackOps(IDA_idamem, 
                        &IDA_iout[10]);         /* NBCKTRK */
  IDAGetTolScaleFactor(IDA_idamem,    
                       &IDA_rout[4]);           /* TOLSFAC */
  
  /* Root finding is on */
  if (IDA_nrtfn != 0)
    IDAGetNumGEvals(IDA_idamem, &IDA_iout[11]); /* NGE */

  /* Linear solver optional outputs */
  IDAGetLinWorkSpace(IDA_idamem, &IDA_iout[12], &IDA_iout[13]);   /* LENRWLS, LENIWLS */
  IDAGetLastLinFlag(IDA_idamem, &IDA_iout[14]);                   /* LSTF */
  IDAGetNumLinResEvals(IDA_idamem, &IDA_iout[15]);                /* NRE */
  IDAGetNumJacEvals(IDA_idamem, &IDA_iout[16]);                   /* NJE */
  IDAGetNumJTSetupEvals(IDA_idamem, &IDA_iout[17]);               /* NJTS */
  IDAGetNumJtimesEvals(IDA_idamem, &IDA_iout[18]);                /* NJT */
  IDAGetNumPrecEvals(IDA_idamem, &IDA_iout[19]);                  /* NPE */
  IDAGetNumPrecSolves(IDA_idamem, &IDA_iout[20]);                 /* NPS */
  IDAGetNumLinIters(IDA_idamem, &IDA_iout[21]);                   /* NLI */
  IDAGetNumLinConvFails(IDA_idamem, &IDA_iout[22]);               /* NCFL */

  return;
}

/*************************************************/

void FIDA_GETDKY(realtype *t, int *k, realtype *dky, int *ier)
{
  /* Store existing F2C_IDA_vec data pointer */
  realtype *f2c_data = N_VGetArrayPointer(F2C_IDA_vec);

  /* Attach user data to vectors */
  N_VSetArrayPointer(dky, F2C_IDA_vec);

  *ier = 0;
  *ier = IDAGetDky(IDA_idamem, *t, *k, F2C_IDA_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(f2c_data, F2C_IDA_vec);

  return;
}

/*************************************************/

void FIDA_GETERRWEIGHTS(realtype *eweight, int *ier)
{
  /* Store existing F2C_IDA_vec data pointer */
  realtype *f2c_data = N_VGetArrayPointer(F2C_IDA_vec);

  /* Attach user data to vector */
  N_VSetArrayPointer(eweight, F2C_IDA_vec);

  *ier = 0;
  *ier = IDAGetErrWeights(IDA_idamem, F2C_IDA_vec);

  /* Reset data pointer */
  N_VSetArrayPointer(f2c_data, F2C_IDA_vec);

  return;
}

/*************************************************/

void FIDA_GETESTLOCALERR(realtype *ele, int *ier)
{
  /* Store existing F2C_IDA_vec data pointer */
  realtype *f2c_data = N_VGetArrayPointer(F2C_IDA_vec);

  /* Attach user data to vector */
  N_VSetArrayPointer(ele, F2C_IDA_vec);

  *ier = 0;
  *ier = IDAGetEstLocalErrors(IDA_idamem, F2C_IDA_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(f2c_data, F2C_IDA_vec);

  return;
}

/*************************************************/

void FIDA_FREE(void)
{
  IDAMem ida_mem;

  ida_mem = (IDAMem) IDA_idamem;

  /* free IDALS interface */
  if (ida_mem->ida_lfree)
    ida_mem->ida_lfree(ida_mem);
  ida_mem->ida_lmem = NULL;

  /* free user_data structure */
  if (ida_mem->ida_user_data)
    free(ida_mem->ida_user_data);
  ida_mem->ida_user_data = NULL;

  /* free main integrator memory structure */
  IDAFree(&IDA_idamem);

  /* free interface vectors / matrices / linear solvers / nonlinear solver */
  N_VSetArrayPointer(NULL, F2C_IDA_vec);
  N_VDestroy(F2C_IDA_vec);
  N_VSetArrayPointer(NULL, F2C_IDA_ypvec);
  N_VDestroy(F2C_IDA_ypvec);
  if (F2C_IDA_ewtvec != NULL) {
    N_VSetArrayPointer(NULL, F2C_IDA_ewtvec);
    N_VDestroy(F2C_IDA_ewtvec);
  }
  if (F2C_IDA_matrix)
    SUNMatDestroy(F2C_IDA_matrix);
  if (F2C_IDA_linsol)
    SUNLinSolFree(F2C_IDA_linsol);
  /* already freed by IDAFree */
  if (F2C_IDA_nonlinsol)
    F2C_IDA_nonlinsol = NULL;
  return;
}

/*************************************************/

int FIDAresfn(realtype t, N_Vector yy, N_Vector yp,
	      N_Vector rr, void *user_data)
{
  int ier;
  realtype *yy_data, *yp_data, *rr_data;
  FIDAUserData IDA_userdata;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_RESFUN(&t, yy_data, yp_data, rr_data, 
              IDA_userdata->ipar, IDA_userdata->rpar, &ier);

  return(ier);
}

/*************************************************/

/* Fortran interface to C routine IDASetNonlinearSolver; see 
   fida.h for further details */
void FIDA_NLSINIT(int *ier) {
  if ( (IDA_idamem == NULL) || (F2C_IDA_nonlinsol == NULL) ) {
    *ier = -1;
    return;
  }
  *ier = IDASetNonlinearSolver(IDA_idamem, F2C_IDA_nonlinsol);
  return;
}
