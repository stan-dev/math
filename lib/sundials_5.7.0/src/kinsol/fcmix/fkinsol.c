/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 *                David J. Gardner @ LLNL
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
 * This is the implementation file for the Fortran interface to
 * the KINSOL package. See fkinsol.h for usage.
 *
 * Note: Some routines are necessarily stored elsewhere to avoid
 * linking problems. See also, therefore, fkinpreco.c, fkinjtimes.c,
 * and fkinbbd.c.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fkinsol.h"     /* prototypes of interfaces and global vars. */
#include "kinsol_impl.h" /* definition of KINMem type                 */

#include <kinsol/kinsol_ls.h> /* KINLS interface routine prototypes   */

/*------------------------------------------------------------------
  definitions of global variables shared amongst various routines
  ------------------------------------------------------------------*/

void *KIN_kinmem;
long int *KIN_iout;
realtype *KIN_rout;

/*------------------------------------------------------------------
  private constants
  ------------------------------------------------------------------*/

#define ZERO RCONST(0.0)

/*------------------------------------------------------------------
  prototype of user-supplied fortran routine
  ------------------------------------------------------------------*/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_FUN(realtype*, realtype*, int*);

#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------
  Function : FKIN_CREATE
  ------------------------------------------------------------------*/

void FKIN_CREATE(int *ier)
{
  
  *ier = 0;
  /* check for required vector operations */
  if ((F2C_KINSOL_vec->ops->nvgetarraypointer == NULL) ||
      (F2C_KINSOL_vec->ops->nvsetarraypointer == NULL)) {
    *ier = -1;
    STAN_SUNDIALS_FPRINTF(stderr, "FKINCREATE: A required vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize pointers to NULL */
  KIN_kinmem = NULL;

  /* Create KINSOL object */
  KIN_kinmem = KINCreate();
  if (KIN_kinmem == NULL) {
    *ier = -1;
    return;
  }
}

/*------------------------------------------------------------------
  Function : FKIN_INIT
  ------------------------------------------------------------------*/

void FKIN_INIT(long int *iout, realtype *rout, int *ier)
{
  
  /* Call KINInit */
  *ier = 0;
  *ier = KINInit(KIN_kinmem, FKINfunc, F2C_KINSOL_vec);

  /* On failure, exit */
  if (*ier != KIN_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Grab optional output arrays and store them in global variables */
  KIN_iout = iout;
  KIN_rout = rout;

  return;
}

/*------------------------------------------------------------------
  Function : FKIN_MALLOC
  ------------------------------------------------------------------*/

void FKIN_MALLOC(long int *iout, realtype *rout, int *ier)
{
  
  /* check for required vector operations */
  if ((F2C_KINSOL_vec->ops->nvgetarraypointer == NULL) ||
      (F2C_KINSOL_vec->ops->nvsetarraypointer == NULL)) {
    *ier = -1;
    STAN_SUNDIALS_FPRINTF(stderr, "A required vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize pointers to NULL */
  KIN_kinmem = NULL;

  /* Create KINSOL object */
  KIN_kinmem = KINCreate();
  if (KIN_kinmem == NULL) {
    *ier = -1;
    return;
  }

  /* Call KINInit */
  *ier = 0;
  *ier = KINInit(KIN_kinmem, FKINfunc, F2C_KINSOL_vec);

  /* On failure, exit */
  if (*ier != KIN_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Grab optional output arrays and store them in global variables */
  KIN_iout = iout;
  KIN_rout = rout;

  return;
}

/*------------------------------------------------------------------
  Function : FKIN_SETIIN
  ------------------------------------------------------------------*/

void FKIN_SETIIN(char key_name[], long int *ival, int *ier)
{
  if (!strncmp(key_name,"PRNT_LEVEL",10))
    *ier = KINSetPrintLevel(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NITERS",10))
    *ier = KINSetNumMaxIters(KIN_kinmem, (long int) *ival);
  else if (!strncmp(key_name,"ETA_FORM",8))
    *ier = KINSetEtaForm(KIN_kinmem, (int) *ival);
  else if (!strncmp(key_name,"MAA",3))
    *ier = KINSetMAA(KIN_kinmem, (long int) *ival);
  else if (!strncmp(key_name,"MAX_SETUPS",10))
    *ier = KINSetMaxSetupCalls(KIN_kinmem, (long int) *ival);
  else if (!strncmp(key_name,"MAX_SP_SETUPS",13))
    *ier = KINSetMaxSubSetupCalls(KIN_kinmem, (long int) *ival);
  else if (!strncmp(key_name,"NO_INIT_SETUP",13))
    *ier = KINSetNoInitSetup(KIN_kinmem, (booleantype) *ival);
  else if (!strncmp(key_name,"NO_MIN_EPS",10))
    *ier = KINSetNoMinEps(KIN_kinmem, (booleantype) *ival);
  else if (!strncmp(key_name,"NO_RES_MON",10))
    *ier = KINSetNoResMon(KIN_kinmem, (booleantype) *ival);
  else {
    *ier = -99;
    STAN_SUNDIALS_FPRINTF(stderr, "FKINSETIIN: Unrecognized key.\n\n");
  }

}

/*------------------------------------------------------------------
  Function : FKIN_SETRIN
  ------------------------------------------------------------------*/

void FKIN_SETRIN(char key_name[], realtype *rval, int *ier)
{

  if (!strncmp(key_name,"FNORM_TOL",9))
    *ier = KINSetFuncNormTol(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"SSTEP_TOL",9))
    *ier = KINSetScaledStepTol(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"MAX_STEP",8))
    *ier = KINSetMaxNewtonStep(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"RERR_FUNC",9))
    *ier = KINSetRelErrFunc(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"ETA_CONST",9))
    *ier = KINSetEtaConstValue(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"ETA_PARAMS",10))
    *ier = KINSetEtaParams(KIN_kinmem, rval[0], rval[1]);
  else if (!strncmp(key_name,"RMON_CONST",10))
    *ier = KINSetResMonConstValue(KIN_kinmem, *rval);
  else if (!strncmp(key_name,"RMON_PARAMS",11))
    *ier = KINSetResMonParams(KIN_kinmem, rval[0], rval[1]);
  else {
    *ier = -99;
    STAN_SUNDIALS_FPRINTF(stderr, "FKINSETRIN: Unrecognized key.\n\n");
  }

}

/*------------------------------------------------------------------
  Function : FKIN_SETVIN
  ------------------------------------------------------------------*/

void FKIN_SETVIN(char key_name[], realtype *vval, int *ier)
{
  N_Vector Vec;

  if (!strncmp(key_name,"CONSTR_VEC",10)) {
    Vec = NULL;
    Vec = N_VCloneEmpty(F2C_KINSOL_vec);
    if (Vec == NULL) {
      *ier = -1;
      return;
    }
    *ier = 0;
    N_VSetArrayPointer(vval, Vec);
    KINSetConstraints(KIN_kinmem, Vec);
    N_VDestroy(Vec);
  } else {
    *ier = -99;
    STAN_SUNDIALS_FPRINTF(stderr, "FKINSETVIN: Unrecognized key.\n\n");
  }

}

/*------------------------------------------------------------------
  Function : FKIN_LSINIT
  ------------------------------------------------------------------*/

/* Fortran interface to C routine KINSetLinearSolver */
void FKIN_LSINIT(int *ier) {
  if ( (KIN_kinmem == NULL) || (F2C_KINSOL_linsol == NULL) ) {
    *ier = -1;
    return;
  }
  *ier = KINSetLinearSolver(KIN_kinmem, F2C_KINSOL_linsol,
                            F2C_KINSOL_matrix);
  return;
}

/*------------------------------------------------------------------
  Function : FKIN_DLSINIT -- DEPRECATED
  ------------------------------------------------------------------*/

void FKIN_DLSINIT(int *ier)
{ FKIN_LSINIT(ier); }

/*------------------------------------------------------------------
  Function : FKIN_SPILSINIT -- DEPRECATED
  ------------------------------------------------------------------*/

void FKIN_SPILSINIT(int *ier)
{ FKIN_LSINIT(ier); }

/*------------------------------------------------------------------
  Function : FKIN_SOL
  ------------------------------------------------------------------*/

void FKIN_SOL(realtype *uu, int *globalstrategy, 
              realtype *uscale , realtype *fscale, int *ier)

{
  N_Vector uuvec, uscalevec, fscalevec;

  *ier = 0;
  uuvec = uscalevec = fscalevec = NULL;

  uuvec = F2C_KINSOL_vec;
  N_VSetArrayPointer(uu, uuvec);

  uscalevec = NULL;
  uscalevec = N_VCloneEmpty(F2C_KINSOL_vec);
  if (uscalevec == NULL) {
    *ier = -4;  /* KIN_MEM_FAIL */
    return;
  }
  N_VSetArrayPointer(uscale, uscalevec);

  fscalevec = NULL;
  fscalevec = N_VCloneEmpty(F2C_KINSOL_vec);
  if (fscalevec == NULL) {
    N_VDestroy(uscalevec);
    *ier = -4;  /* KIN_MEM_FAIL */
    return;
  }
  N_VSetArrayPointer(fscale, fscalevec);

  /* If using the fixed-point solver, initialize F2C_KINSOL_linsol 
     and F2C_KINSOL_matrix to NULL */
  if (*globalstrategy == KIN_FP) {
    FKINNullMatrix();
    FKINNullLinsol();
  }
  
  /* Call main solver function */
  *ier = KINSol(KIN_kinmem, uuvec, *globalstrategy, uscalevec, fscalevec);

  N_VSetArrayPointer(NULL, uuvec);

  N_VSetArrayPointer(NULL, uscalevec);
  N_VDestroy(uscalevec);

  N_VSetArrayPointer(NULL, fscalevec);
  N_VDestroy(fscalevec);

  /* load optional outputs into iout[] and rout[] */
  KINGetWorkSpace(KIN_kinmem, &KIN_iout[0], &KIN_iout[1]);   /* LENRW & LENIW */
  KINGetNumNonlinSolvIters(KIN_kinmem, &KIN_iout[2]);        /* NNI */
  KINGetNumFuncEvals(KIN_kinmem, &KIN_iout[3]);              /* NFE */
  KINGetNumBetaCondFails(KIN_kinmem, &KIN_iout[4]);          /* NBCF */
  KINGetNumBacktrackOps(KIN_kinmem, &KIN_iout[5]);           /* NBCKTRK */

  KINGetFuncNorm(KIN_kinmem, &KIN_rout[0]);                  /* FNORM */
  KINGetStepLength(KIN_kinmem, &KIN_rout[1]);                /* SSTEP */

  KINGetLinWorkSpace(KIN_kinmem, &KIN_iout[6], &KIN_iout[7]); /* LRW & LIW */
  KINGetLastLinFlag(KIN_kinmem, &KIN_iout[8]);                /* LSTF */
  KINGetNumLinFuncEvals(KIN_kinmem, &KIN_iout[9]);            /* NFE */
  KINGetNumJacEvals(KIN_kinmem, &KIN_iout[10]);               /* NJE */    
  KINGetNumJtimesEvals(KIN_kinmem, &KIN_iout[11]);            /* NJT */
  KINGetNumPrecEvals(KIN_kinmem, &KIN_iout[12]);              /* NPE */
  KINGetNumPrecSolves(KIN_kinmem, &KIN_iout[13]);             /* NPS */
  KINGetNumLinIters(KIN_kinmem, &KIN_iout[14]);               /* NLI */
  KINGetNumLinConvFails(KIN_kinmem, &KIN_iout[15]);           /* NCFL */

  return;
}

/*------------------------------------------------------------------
  Function : FKIN_FREE
  ------------------------------------------------------------------*/

void FKIN_FREE(void)
{
  KINMem kin_mem;

  kin_mem = (KINMem) KIN_kinmem;

  /* free LS interface */
  if (kin_mem->kin_lfree)
    kin_mem->kin_lfree(kin_mem);
  kin_mem->kin_lmem = NULL;

  /* free user_data structure */
  if (kin_mem->kin_user_data)
    free(kin_mem->kin_user_data);
  kin_mem->kin_user_data = NULL;

  /* free main solver memory structure */
  KINFree(&KIN_kinmem);

  /* free interface vectors / matrices / linear solvers */
  N_VSetArrayPointer(NULL, F2C_KINSOL_vec);
  N_VDestroy(F2C_KINSOL_vec);
  if (F2C_KINSOL_matrix)
    SUNMatDestroy(F2C_KINSOL_matrix);
  if (F2C_KINSOL_linsol)
    SUNLinSolFree(F2C_KINSOL_linsol);

  return;
}


/*------------------------------------------------------------------
  Function : FKINfunc
  ------------------------------------------------------------------
  The C function FKINfunc acts as an interface between KINSOL and
  the Fortran user-supplied subroutine FKFUN. Addresses of the
  data uu and fdata are passed to FKFUN, using the routine
  N_VGetArrayPointer from the NVECTOR module. The data in the
  returned N_Vector fval is set using N_VSetArrayPointer. Auxiliary
  data is assumed to be communicated by 'Common'.
  ------------------------------------------------------------------*/

int FKINfunc(N_Vector uu, N_Vector fval, void *user_data)
{
  realtype *udata, *fdata;
  int ier;

  udata = N_VGetArrayPointer(uu);
  fdata = N_VGetArrayPointer(fval);

  FK_FUN(udata, fdata, &ier);

  return(ier);
}
