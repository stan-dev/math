/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the ARKODE package.  See farkode.h for usage.
 * NOTE: some routines are necessarily stored elsewhere to avoid
 * linking problems.  Therefore, see also the other C files in
 * this folder for all of the available options.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <sundials/sundials_matrix.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_arkstep.h>

/*=============================================================*/

/* Constants and default values (in case of illegal inputs) */
#define  ABSTOL  RCONST(1.0e-9)
#define  RELTOL  RCONST(1.0e-4)
#define  ZERO    RCONST(0.0)

/*=============================================================*/

/* Definitions for global variables shared between Fortran/C
   interface routines */
void     *ARK_arkodemem;
long int *ARK_iout;
realtype *ARK_rout;
int       ARK_nrtfn;
int       ARK_ls;
int       ARK_mass_ls;

/*=============================================================*/

/* Prototypes of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_IMP_FUN(realtype *T, realtype *Y, realtype *YDOT,
                           long int *IPAR, realtype *RPAR, int *IER);
  extern void FARK_EXP_FUN(realtype *T, realtype *Y, realtype *YDOT,
                           long int *IPAR, realtype *RPAR, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface routine to initialize ARKStep memory
   structure; functions as an all-in-one interface to the C
   routines ARKStepCreate, ARKStepSetUserData, and
   ARKStepSStolerances (or ARKStepSVtolerances); see farkode.h
   for further details */
void FARK_MALLOC(realtype *t0, realtype *y0, int *imex,
                 int *iatol, realtype *rtol, realtype *atol,
                 long int *iout, realtype *rout,
                 long int *ipar, realtype *rpar, int *ier) {

  N_Vector Vatol;
  FARKUserData ARK_userdata;
  realtype reltol, abstol;

  *ier = 0;

  /* Check for required vector operations */
  if(F2C_ARKODE_vec->ops->nvgetarraypointer == NULL) {
    *ier = -1;
    fprintf(stderr, "Error: getarraypointer vector operation is not implemented.\n\n");
    return;
  }
  if(F2C_ARKODE_vec->ops->nvsetarraypointer == NULL) {
    *ier = -1;
    fprintf(stderr, "Error: setarraypointer vector operation is not implemented.\n\n");
    return;
  }
  if(F2C_ARKODE_vec->ops->nvcloneempty == NULL) {
    *ier = -1;
    fprintf(stderr, "Error: cloneempty vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize all pointers to NULL */
  ARK_arkodemem = NULL;
  Vatol = NULL;

  /* initialize global constants to disable each option */
  ARK_nrtfn = 0;
  ARK_ls = SUNFALSE;
  ARK_mass_ls = SUNFALSE;

  /* Set data in F2C_ARKODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_ARKODE_vec);

  /* Call ARKStepCreate based on imex argument */
  switch (*imex) {
  case 0:  /* purely implicit */
    ARK_arkodemem = ARKStepCreate(NULL, FARKfi, *t0, F2C_ARKODE_vec);
    break;
  case 1:  /* purely explicit */
    ARK_arkodemem = ARKStepCreate(FARKfe, NULL, *t0, F2C_ARKODE_vec);
    FARKNullMatrix();
    FARKNullLinsol();
    FARKNullNonlinsol();
    break;
  case 2:  /* imex */
    ARK_arkodemem = ARKStepCreate(FARKfe, FARKfi, *t0, F2C_ARKODE_vec);
    break;
  }
  if (ARK_arkodemem == NULL) {
    *ier = -1;
    return;
  }

  /* Set and attach user data */
  ARK_userdata = NULL;
  ARK_userdata = (FARKUserData) malloc(sizeof *ARK_userdata);
  if (ARK_userdata == NULL) {
    *ier = -1;
    return;
  }
  ARK_userdata->rpar = rpar;
  ARK_userdata->ipar = ipar;
  *ier = ARKStepSetUserData(ARK_arkodemem, ARK_userdata);
  if(*ier != ARK_SUCCESS) {
    free(ARK_userdata); ARK_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);

  /* Set tolerances -- if <= 0, keep as defaults */
  reltol = RELTOL;
  abstol = ABSTOL;
  if (*rtol > ZERO)  reltol = *rtol;
  switch (*iatol) {
  case 1:
    if (*atol > ZERO)  abstol = *atol;
    *ier = ARKStepSStolerances(ARK_arkodemem, reltol, abstol);
    break;
  case 2:
    Vatol = N_VCloneEmpty(F2C_ARKODE_vec);
    if (Vatol == NULL) {
      free(ARK_userdata);
      ARK_userdata = NULL;
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    if (N_VMin(Vatol) <= ZERO)  N_VConst(abstol, Vatol);
    *ier = ARKStepSVtolerances(ARK_arkodemem, reltol, Vatol);
    N_VDestroy(Vatol);
    break;
  }

  /* On failure, exit */
  if(*ier != ARK_SUCCESS) {
    free(ARK_userdata);
    ARK_userdata = NULL;
    *ier = -1;
    return;
  }

  /* store pointers to optional output arrays in global vars */
  ARK_iout = iout;
  ARK_rout = rout;

  /* Store the unit roundoff in rout for user access */
  ARK_rout[5] = UNIT_ROUNDOFF;

  return;
}

/*=============================================================*/

/* Fortran interface routine to re-initialize ARKStep memory
   structure; functions as an all-in-one interface to the C
   routines ARKStepReInit and ARKStepSStolerances (or
   ARKStepSVtolerances); see farkode.h for further details */
void FARK_REINIT(realtype *t0, realtype *y0, int *imex, int *iatol,
                 realtype *rtol, realtype *atol, int *ier) {

  N_Vector Vatol;
  realtype reltol, abstol;
  *ier = 0;

  /* Initialize all pointers to NULL */
  Vatol = NULL;

  /* Set data in F2C_ARKODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_ARKODE_vec);

  /* Call ARKStepReInit based on imex argument */
  switch (*imex) {
  case 0:  /* purely implicit */
    *ier = ARKStepReInit(ARK_arkodemem, NULL, FARKfi,
                         *t0, F2C_ARKODE_vec);
    break;
  case 1:  /* purely explicit */
    *ier = ARKStepReInit(ARK_arkodemem, FARKfe, NULL,
                         *t0, F2C_ARKODE_vec);
    break;
  case 2:  /* imex */
    *ier = ARKStepReInit(ARK_arkodemem, FARKfe, FARKfi,
                         *t0, F2C_ARKODE_vec);
    break;
  }

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);

  /* On failure, exit */
  if (*ier != ARK_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Set tolerances */
  reltol = RELTOL;
  abstol = ABSTOL;
  if (*rtol > ZERO)  reltol = *rtol;
  switch (*iatol) {
  case 1:
    if (*atol > ZERO)  abstol = *atol;
    *ier = ARKStepSStolerances(ARK_arkodemem, reltol, abstol);
    break;
  case 2:
    Vatol = N_VCloneEmpty(F2C_ARKODE_vec);
    if (Vatol == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    if (N_VMin(Vatol) <= ZERO)  N_VConst(abstol, Vatol);
    *ier = ARKStepSVtolerances(ARK_arkodemem, reltol, Vatol);
    N_VDestroy(Vatol);
    break;
  }

  /* On failure, exit */
  if (*ier != ARK_SUCCESS) {
    *ier = -1;
    return;
  }

  return;
}

/*=============================================================*/

/* Fortran interface routine to re-initialize ARKStep memory
   structure for a problem with a new size but similar time
   scale; functions as an all-in-one interface to the C
   routines ARKStepResize (and potentially ARKStepSVtolerances);
   see farkode.h for further details */
void FARK_RESIZE(realtype *t0, realtype *y0, realtype *hscale,
                 int *itol, realtype *rtol, realtype *atol, int *ier) {

  *ier = 0;

  /* Set data in F2C_ARKODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_ARKODE_vec);

  /* Call ARKStepResize (currently does not allow Fortran
     user-supplied vector resize function) */
  *ier = ARKStepResize(ARK_arkodemem, F2C_ARKODE_vec, *hscale,
                       *t0, NULL, NULL);

  /* Reset data pointer */
  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);

  /* On failure, exit */
  if (*ier != ARK_SUCCESS) {
    *ier = -1;
    return;
  }

  /* Set tolerances, based on itol argument */
  if (*itol) {
    N_Vector Vatol = NULL;
    Vatol = N_VCloneEmpty(F2C_ARKODE_vec);
    if (Vatol == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    *ier = ARKStepSVtolerances(ARK_arkodemem, *rtol, Vatol);
    N_VDestroy(Vatol);
  }

  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetDefaults; see
   farkode.h for further details */
void FARK_SETDEFAULTS(int *ier) {
  *ier += ARKStepSetDefaults(ARK_arkodemem);
  return;
}

/*=============================================================*/

/* Fortran interface to C "set" routines having integer
   arguments; see farkode.h for further details */
void FARK_SETIIN(char key_name[], long int *ival, int *ier) {
  if (!strncmp(key_name, "ORDER", 5))
    *ier = ARKStepSetOrder(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "DENSE_ORDER", 11))
    *ier = ARKStepSetDenseOrder(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "LINEAR", 6))
    *ier = ARKStepSetLinear(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "NONLINEAR", 9))
    *ier = ARKStepSetNonlinear(ARK_arkodemem);
  else if (!strncmp(key_name, "EXPLICIT", 8))
    *ier = ARKStepSetExplicit(ARK_arkodemem);
  else if (!strncmp(key_name, "IMPLICIT", 8))
    *ier = ARKStepSetImplicit(ARK_arkodemem);
  else if (!strncmp(key_name, "IMEX", 4))
    *ier = ARKStepSetImEx(ARK_arkodemem);
  else if (!strncmp(key_name, "IRK_TABLE_NUM", 13))
    *ier = ARKStepSetTableNum(ARK_arkodemem, (int) *ival, -1);
  else if (!strncmp(key_name, "ERK_TABLE_NUM", 13))
    *ier = ARKStepSetTableNum(ARK_arkodemem, -1, (int) *ival);
  else if (!strncmp(key_name, "ARK_TABLE_NUM", 13))
    *ier = ARKStepSetTableNum(ARK_arkodemem, (int) ival[0], (int) ival[1]);
  else if (!strncmp(key_name, "MAX_NSTEPS", 10))
    *ier = ARKStepSetMaxNumSteps(ARK_arkodemem, (long int) *ival);
  else if (!strncmp(key_name, "HNIL_WARNS", 10))
    *ier = ARKStepSetMaxHnilWarns(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "PREDICT_METHOD", 14))
    *ier = ARKStepSetPredictorMethod(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "MAX_ERRFAIL", 11))
    *ier = ARKStepSetMaxErrTestFails(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "MAX_CONVFAIL", 12))
    *ier = ARKStepSetMaxConvFails(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "MAX_NITERS", 10))
    *ier = ARKStepSetMaxNonlinIters(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "ADAPT_SMALL_NEF", 15))
    *ier = ARKStepSetSmallNumEFails(ARK_arkodemem, (int) *ival);
  else if (!strncmp(key_name, "LSETUP_MSBP", 11))
    *ier = ARKStepSetMaxStepsBetweenLSet(ARK_arkodemem, (int) *ival);
  else {
    *ier = -99;
    fprintf(stderr, "FARKSETIIN: Unrecognized key.\n\n");
  }
  return;
}

/*=============================================================*/

/* Fortran interface to C "set" routines having real
   arguments; see farkode.h for further details */
void FARK_SETRIN(char key_name[], realtype *rval, int *ier) {
  if (!strncmp(key_name, "INIT_STEP", 9))
    *ier = ARKStepSetInitStep(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "MAX_STEP", 8))
    *ier = ARKStepSetMaxStep(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "MIN_STEP", 8))
    *ier = ARKStepSetMinStep(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "STOP_TIME", 9))
    *ier = ARKStepSetStopTime(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "NLCONV_COEF", 11))
    *ier = ARKStepSetNonlinConvCoef(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_CFL", 9))
    *ier = ARKStepSetCFLFraction(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_SAFETY", 12))
    *ier = ARKStepSetSafetyFactor(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_BIAS", 10))
    *ier = ARKStepSetErrorBias(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_GROWTH", 12))
    *ier = ARKStepSetMaxGrowth(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_BOUNDS", 12))
    *ier = ARKStepSetFixedStepBounds(ARK_arkodemem, rval[0], rval[1]);
  else if (!strncmp(key_name, "ADAPT_ETAMX1", 12))
    *ier = ARKStepSetMaxFirstGrowth(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_ETAMXF", 12))
    *ier = ARKStepSetMaxEFailGrowth(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "ADAPT_ETACF", 11))
    *ier = ARKStepSetMaxCFailGrowth(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "NONLIN_CRDOWN", 11))
    *ier = ARKStepSetNonlinCRDown(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "NONLIN_RDIV", 9))
    *ier = ARKStepSetNonlinRDiv(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "LSETUP_DGMAX", 12))
    *ier = ARKStepSetDeltaGammaMax(ARK_arkodemem, *rval);
  else if (!strncmp(key_name, "FIXED_STEP", 10))
    *ier = ARKStepSetFixedStep(ARK_arkodemem, *rval);
  else {
    *ier = -99;
    fprintf(stderr, "FARKSETRIN: Unrecognized key: %s\n\n",key_name);
  }
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetAdaptivityMethod;
   see farkode.h for further details */
void FARK_SETADAPTMETHOD(int *imethod, int *idefault, int *ipq,
                         realtype *params, int *ier) {

  *ier = ARKStepSetAdaptivityMethod(ARK_arkodemem, *imethod,
                                    *idefault, *ipq, params);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetTables; see
   farkode.h for further details */
void FARK_SETERKTABLE(int *s, int *q, int *p, realtype *c, realtype *A,
                      realtype *b, realtype *b2, int *ier) {
  ARKodeButcherTable Be;
  Be = ARKodeButcherTable_Create(*s, *q, *p, c, A, b, b2);
  *ier = ARKStepSetTables(ARK_arkodemem, *q, *p, NULL, Be);
  ARKodeButcherTable_Free(Be);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetTables; see
   farkode.h for further details */
void FARK_SETIRKTABLE(int *s, int *q, int *p, realtype *c, realtype *A,
                      realtype *b, realtype *b2, int *ier) {
  ARKodeButcherTable Bi;
  Bi = ARKodeButcherTable_Create(*s, *q, *p, c, A, b, b2);
  *ier = ARKStepSetTables(ARK_arkodemem, *q, *p, Bi, NULL);
  ARKodeButcherTable_Free(Bi);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetTables; see
   farkode.h for further details */
void FARK_SETARKTABLES(int *s, int *q, int *p, realtype *ci,
                       realtype *ce, realtype *Ai, realtype *Ae,
                       realtype *bi, realtype *be, realtype *b2i,
                       realtype *b2e, int *ier) {
  ARKodeButcherTable Bi, Be;
  Bi = ARKodeButcherTable_Create(*s, *q, *p, ci, Ai, bi, b2i);
  Be = ARKodeButcherTable_Create(*s, *q, *p, ce, Ae, be, b2e);
  *ier = ARKStepSetTables(ARK_arkodemem, *q, *p, Bi, Be);
  ARKodeButcherTable_Free(Bi);
  ARKodeButcherTable_Free(Be);
  return;
}

/*=============================================================*/

/* Fortran interface routine to set residual tolerance
   scalar/array; functions as an all-in-one interface to the C
   routines ARKStepResStolerance and ARKStepResVtolerance;
   see farkode.h for further details */
void FARK_SETRESTOLERANCE(int *itol, realtype *atol, int *ier) {

  N_Vector Vatol;
  realtype abstol;

  *ier = 0;

  /* Set tolerance, based on itol argument */
  abstol = ABSTOL;
  switch (*itol) {
  case 1:
    if (*atol > ZERO)  abstol = *atol;
    *ier = ARKStepResStolerance(ARK_arkodemem, abstol);
    break;
  case 2:
    Vatol = NULL;
    Vatol = N_VCloneEmpty(F2C_ARKODE_vec);
    if (Vatol == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    if (N_VMin(Vatol) <= ZERO)  N_VConst(abstol, Vatol);
    *ier = ARKStepResVtolerance(ARK_arkodemem, Vatol);
    N_VDestroy(Vatol);
    break;
  }

  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetDiagnostics; see
   farkode.h for further details */
void FARK_SETDIAGNOSTICS(char fname[], int *flen, int *ier) {
  char *filename=NULL;
  FILE *DFID=NULL;
  int i;

  /* copy fname into array of specified length */
  filename = (char *) malloc((*flen)*sizeof(char));
  for (i=0; i<*flen; i++)  filename[i] = fname[i];

  /* open diagnostics output file */
  DFID = fopen(filename,"w");
  if (DFID == NULL) {
    *ier = 1;
    return;
  }
  *ier = ARKStepSetDiagnostics(ARK_arkodemem, DFID);
  free(filename);
  return;
}

/*=============================================================*/

/* Fortran routine to close diagnostics output file; see farkode.h
   for further details */
void FARK_STOPDIAGNOSTICS(int *ier) {
  ARKodeMem ark_mem;
  if (ARK_arkodemem == NULL) {
    *ier = 1;
    return;
  }
  ark_mem = (ARKodeMem) ARK_arkodemem;

  if (ark_mem->diagfp == NULL) {
    *ier = 1;
    return;
  }
  *ier = fclose(ark_mem->diagfp);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetNonlinearSolver */
void FARK_NLSINIT(int *ier) {
  if ( (ARK_arkodemem == NULL) || (F2C_ARKODE_nonlinsol == NULL) ) {
    *ier = -1;
    return;
  }
  *ier = ARKStepSetNonlinearSolver(ARK_arkodemem, F2C_ARKODE_nonlinsol);
  return;
}

/*=============================================================*/

/* ---DEPRECATED---
   Fortran interface to C routine ARKStepSetLinearSolver; see
   farkode.h for further details */
void FARK_DLSINIT(int *ier)
{ FARK_LSINIT(ier); }

/* ---DEPRECATED---
   Fortran interface to C routine ARKStepSetMassLinearSolver; see
   farkode.h for further details */
void FARK_DLSMASSINIT(int *time_dep, int *ier)
{ FARK_LSMASSINIT(time_dep, ier); }

/*=============================================================*/

/* ---DEPRECATED---
   Fortran interface to C routine ARKStepSetLinearSolver; see
   farkode.h for further details */
void FARK_SPILSINIT(int *ier)
{ FARK_LSINIT(ier); }

/* ---DEPRECATED---
   Fortran interface to C routine ARKStepSetMassLinearSolver; see
   farkode.h for further details */
void FARK_SPILSMASSINIT(int *time_dep, int *ier)
{ FARK_LSMASSINIT(time_dep, ier); }

/*=============================================================*/

/* ---DEPRECATED---
   Fortran interfaces to C "set" routines for the ARKStep linear
   solver; see farkode.h for further details */
void FARK_SPILSSETEPSLIN(realtype *eplifac, int *ier)
{ FARK_LSSETEPSLIN(eplifac, ier); }

void FARK_SPILSSETMASSEPSLIN(realtype *eplifac, int *ier)
{ FARK_LSSETMASSEPSLIN(eplifac, ier); }

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetLinearSolver; see
   farkode.h for further details */
void FARK_LSINIT(int *ier) {
  if ( (ARK_arkodemem == NULL) || (F2C_ARKODE_linsol == NULL) ) {
    *ier = -1;
    return;
  }
  *ier = ARKStepSetLinearSolver(ARK_arkodemem, F2C_ARKODE_linsol,
                                F2C_ARKODE_matrix);
  ARK_ls = SUNTRUE;
  return;
}

/* Fortran interface to C routine ARKStepSetMassLinearSolver; see
   farkode.h for further details */
void FARK_LSMASSINIT(int *time_dep, int *ier) {
  if ( (ARK_arkodemem == NULL) || (F2C_ARKODE_mass_sol == NULL) ) {
    *ier = -1;
    return;
  }
  *ier = ARKStepSetMassLinearSolver(ARK_arkodemem,
                                    F2C_ARKODE_mass_sol,
                                    F2C_ARKODE_mass_matrix, 
                                    *time_dep);
  ARK_mass_ls = SUNTRUE;
  return;
}

/*=============================================================*/

/* Fortran interfaces to C "set" routines for the ARKStep linear
   solver; see farkode.h for further details */
void FARK_LSSETEPSLIN(realtype *eplifac, int *ier)
{ *ier = ARKStepSetEpsLin(ARK_arkodemem, *eplifac); }

void FARK_LSSETMASSEPSLIN(realtype *eplifac, int *ier)
{ *ier = ARKStepSetMassEpsLin(ARK_arkodemem, *eplifac); }

/*=============================================================*/

/* Fortran interface to C routine ARKStepEvolve (the main integrator)
   and optional output routines ARKStepGet*; see farkode.h for
   further details */
void FARK_ARKODE(realtype *tout, realtype *t, realtype *y,
                 int *itask, int *ier) {

  /* attach user solution array to solver memory */
  N_VSetArrayPointer(y, F2C_ARKODE_vec);

  /* call ARKStepEvolve solver */
  *ier = ARKStepEvolve(ARK_arkodemem, *tout, F2C_ARKODE_vec, t, *itask);

  /* detach user solution array from solver memory */
  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);

  /* Load optional outputs in iout & rout */
  ARKStepGetWorkSpace(ARK_arkodemem,
                      &ARK_iout[0],          /* LENRW   */
                      &ARK_iout[1]);         /* LENIW   */
  ARKStepGetStepStats(ARK_arkodemem,
                      &ARK_iout[2],    /* NST     */
                      &ARK_rout[0],    /* H0U     */
                      &ARK_rout[1],    /* HU      */
                      &ARK_rout[2],    /* HCUR    */
                      &ARK_rout[3]);   /* TCUR    */
  ARKStepGetTimestepperStats(ARK_arkodemem,
                             &ARK_iout[3],    /* NST_STB */
                             &ARK_iout[4],    /* NST_ACC */
                             &ARK_iout[5],    /* NST_ATT */
                             &ARK_iout[6],    /* NFE     */
                             &ARK_iout[7],    /* NFI     */
                             &ARK_iout[8],    /* NSETUPS */
                             &ARK_iout[9]);   /* NETF    */
  ARKStepGetTolScaleFactor(ARK_arkodemem,
                           &ARK_rout[4]);    /* TOLSFAC */
  ARKStepGetNonlinSolvStats(ARK_arkodemem,
                            &ARK_iout[10],   /* NNI     */
                            &ARK_iout[11]);  /* NCFN    */

  /* If root finding is on, load those outputs as well */
  if (ARK_nrtfn != 0)
    ARKStepGetNumGEvals(ARK_arkodemem, &ARK_iout[12]);  /* NGE */

  /* Attach linear solver outputs */
  if (ARK_ls) {
    ARKStepGetLinWorkSpace(ARK_arkodemem, &ARK_iout[13], &ARK_iout[14]); /* LENRWLS, LENIWLS */
    ARKStepGetLastLinFlag(ARK_arkodemem,     &ARK_iout[15]);            /* LSTF  */
    ARKStepGetNumLinRhsEvals(ARK_arkodemem,  &ARK_iout[16]);            /* NFELS */
    ARKStepGetNumJacEvals(ARK_arkodemem,     &ARK_iout[17]);            /* NJE   */
    ARKStepGetNumJTSetupEvals(ARK_arkodemem, &ARK_iout[18]);            /* NJTS  */
    ARKStepGetNumJtimesEvals(ARK_arkodemem,  &ARK_iout[19]);            /* NJTV  */
    ARKStepGetNumPrecEvals(ARK_arkodemem,    &ARK_iout[20]);            /* NPE   */
    ARKStepGetNumPrecSolves(ARK_arkodemem,   &ARK_iout[21]);            /* NPS   */
    ARKStepGetNumLinIters(ARK_arkodemem,     &ARK_iout[22]);            /* NLI   */
    ARKStepGetNumLinConvFails(ARK_arkodemem, &ARK_iout[23]);            /* NCFL  */
  }
  
  /* Attach mass matrix linear solver outputs */
  if(ARK_mass_ls) {
    ARKStepGetMassWorkSpace(ARK_arkodemem, &ARK_iout[24], &ARK_iout[25]); /* LENRWMS, LENIWMS */
    ARKStepGetLastMassFlag(ARK_arkodemem,      &ARK_iout[26]);            /* LSTMF  */
    ARKStepGetNumMassSetups(ARK_arkodemem,     &ARK_iout[27]);            /* NMSET  */
    ARKStepGetNumMassSolves(ARK_arkodemem,     &ARK_iout[28]);            /* NMSOL  */
    ARKStepGetNumMTSetups(ARK_arkodemem,       &ARK_iout[29]);            /* NMTSET */
    ARKStepGetNumMassMult(ARK_arkodemem,       &ARK_iout[30]);            /* NMMUL  */
    ARKStepGetNumMassPrecEvals(ARK_arkodemem,  &ARK_iout[31]);            /* NMPE   */
    ARKStepGetNumMassPrecSolves(ARK_arkodemem, &ARK_iout[32]);            /* NMPS   */
    ARKStepGetNumMassIters(ARK_arkodemem,      &ARK_iout[33]);            /* NMLI   */
    ARKStepGetNumMassConvFails(ARK_arkodemem,  &ARK_iout[34]);            /* NMCFL  */
  }
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepGetDky; see farkode.h
   for further details */
void FARK_DKY(realtype *t, int *k, realtype *dky, int *ier) {

  /* store pointer existing F2C_ARKODE_vec data array */
  realtype *f2c_data = N_VGetArrayPointer(F2C_ARKODE_vec);

  /* attach output data array to F2C_ARKODE_vec */
  N_VSetArrayPointer(dky, F2C_ARKODE_vec);

  /* call ARKStepGetDky */
  *ier = 0;
  *ier = ARKStepGetDky(ARK_arkodemem, *t, *k, F2C_ARKODE_vec);

  /* reattach F2C_ARKODE_vec to previous data array */
  N_VSetArrayPointer(f2c_data, F2C_ARKODE_vec);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepGetErrWeights; see
   farkode.h for further details */
void FARK_GETERRWEIGHTS(realtype *eweight, int *ier) {

  /* store pointer existing F2C_ARKODE_vec data array */
  realtype *f2c_data = N_VGetArrayPointer(F2C_ARKODE_vec);

  /* attach output data array to F2C_ARKODE_vec */
  N_VSetArrayPointer(eweight, F2C_ARKODE_vec);

  /* call ARKStepGetErrWeights */
  *ier = 0;
  *ier = ARKStepGetErrWeights(ARK_arkodemem, F2C_ARKODE_vec);

  /* reattach F2C_ARKODE_vec to previous data array */
  N_VSetArrayPointer(f2c_data, F2C_ARKODE_vec);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepGetResWeights; see
   farkode.h for further details */
void FARK_GETRESWEIGHTS(realtype *rweight, int *ier) {

  /* store pointer existing F2C_ARKODE_vec data array */
  realtype *f2c_data = N_VGetArrayPointer(F2C_ARKODE_vec);

  /* attach output data array to F2C_ARKODE_vec */
  N_VSetArrayPointer(rweight, F2C_ARKODE_vec);

  /* call ARKStepGetResWeights */
  *ier = 0;
  *ier = ARKStepGetResWeights(ARK_arkodemem, F2C_ARKODE_vec);

  /* reattach F2C_ARKODE_vec to previous data array */
  N_VSetArrayPointer(f2c_data, F2C_ARKODE_vec);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepGetEstLocalErrors; see
   farkode.h for further details */
void FARK_GETESTLOCALERR(realtype *ele, int *ier) {

  /* store pointer existing F2C_ARKODE_vec data array */
  realtype *f2c_data = N_VGetArrayPointer(F2C_ARKODE_vec);

  /* attach output data array to F2C_ARKODE_vec */
  N_VSetArrayPointer(ele, F2C_ARKODE_vec);

  /* call ARKStepGetEstLocalErrors */
  *ier = 0;
  *ier = ARKStepGetEstLocalErrors(ARK_arkodemem, F2C_ARKODE_vec);

  /* reattach F2C_ARKODE_vec to previous data array */
  N_VSetArrayPointer(f2c_data, F2C_ARKODE_vec);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepFree; see farkode.h for
   further details */
void FARK_FREE() {

  ARKodeMem ark_mem;
  ark_mem = (ARKodeMem) ARK_arkodemem;

  /* free user_data structure */
  if (ark_mem->user_data)
    free(ark_mem->user_data);
  ark_mem->user_data = NULL;

  /* free main integrator memory structure (internally
     frees time step module, rootfinding, interpolation structures) */
  ARKStepFree(&ARK_arkodemem);

  /* free interface vector / matrices / linear solvers */
  N_VSetArrayPointer(NULL, F2C_ARKODE_vec);
  N_VDestroy(F2C_ARKODE_vec);
  if (F2C_ARKODE_matrix)
    SUNMatDestroy(F2C_ARKODE_matrix);
  if (F2C_ARKODE_mass_matrix)
    SUNMatDestroy(F2C_ARKODE_mass_matrix);
  if (F2C_ARKODE_linsol)
    SUNLinSolFree(F2C_ARKODE_linsol);
  if (F2C_ARKODE_mass_sol)
    SUNLinSolFree(F2C_ARKODE_mass_sol);
  return;
}

/*=============================================================*/

/* Fortran interface to C routineARKStepWriteParameters; see
   farkode.h for further details */
void FARK_WRITEPARAMETERS(int *ier) {
  *ier += ARKStepWriteParameters(ARK_arkodemem, stdout);
  return;
}

/*=============================================================*/

/* C interface to user-supplied FORTRAN function FARKEFUN; see
   farkode.h for further details */
int FARKfe(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  int ier;
  realtype *ydata, *dydata;
  FARKUserData ARK_userdata;
  ydata  = N_VGetArrayPointer(y);
  dydata = N_VGetArrayPointer(ydot);
  ARK_userdata = (FARKUserData) user_data;

  FARK_EXP_FUN(&t, ydata, dydata, ARK_userdata->ipar,
               ARK_userdata->rpar, &ier);
  return(ier);
}

/*=============================================================*/

/* C interface to user-supplied FORTRAN function FARKIFUN; see
   farkode.h for further details */
int FARKfi(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  int ier;
  realtype *ydata, *dydata;
  FARKUserData ARK_userdata;
  ydata  = N_VGetArrayPointer(y);
  dydata = N_VGetArrayPointer(ydot);
  ARK_userdata = (FARKUserData) user_data;

  FARK_IMP_FUN(&t, ydata, dydata, ARK_userdata->ipar,
               ARK_userdata->rpar, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
