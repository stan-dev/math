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
 * The C functions FARKPSet and FARKPSol are to interface between 
 * the ARKLS module and the user-supplied preconditioner 
 * setup/solve routines FARKPSET and FARKPSOL. Note the use of 
 * the generic names FARK_PSET and FARK_PSOL in the code below.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_arkstep.h>

/*=============================================================*/

/* Prototype of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_PSET(realtype *T, realtype *Y, realtype *FY,
                        booleantype *JOK, booleantype *JCUR,
                        realtype *GAMMA, realtype *H,
                        long int *IPAR, realtype *RPAR, int *IER);
  extern void FARK_PSOL(realtype *T, realtype *Y, realtype *FY,
                        realtype *R, realtype *Z, 
                        realtype *GAMMA, realtype *DELTA,
                        int *LR, long int *IPAR, realtype *RPAR, 
                        int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* ---DEPRECATED---
   Fortran interface to C routine ARKStepSetPreconditioner; see 
   farkode.h for further details */
void FARK_SPILSSETPREC(int *flag, int *ier)
{ FARK_LSSETPREC(flag, ier); }

/* Fortran interface to C routine ARKStepSetPreconditioner; see 
   farkode.h for further details */
void FARK_LSSETPREC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = ARKStepSetPreconditioner(ARK_arkodemem, NULL, NULL);
  } else {
    *ier = ARKStepSetPreconditioner(ARK_arkodemem, 
                                    FARKPSet, FARKPSol);
  }
  return;
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKPSET; see 
   farkode.h for further details */
int FARKPSet(realtype t, N_Vector y, N_Vector fy, 
             booleantype jok, booleantype *jcurPtr, 
             realtype gamma, void *user_data)
{
  int ier = 0;
  realtype *ydata, *fydata;
  realtype h;
  FARKUserData ARK_userdata;

  ARKStepGetLastStep(ARK_arkodemem, &h);
  ydata  = N_VGetArrayPointer(y);
  fydata = N_VGetArrayPointer(fy);
  ARK_userdata = (FARKUserData) user_data;

  FARK_PSET(&t, ydata, fydata, &jok, jcurPtr, &gamma, &h,
            ARK_userdata->ipar, ARK_userdata->rpar, &ier);
  return(ier);
}


/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKPSOL; see 
   farkode.h for further details */
int FARKPSol(realtype t, N_Vector y, N_Vector fy, N_Vector r, 
             N_Vector z, realtype gamma, realtype delta,
             int lr, void *user_data)
{
  int ier = 0;
  realtype *ydata, *fydata, *rdata, *zdata;
  FARKUserData ARK_userdata;

  ydata  = N_VGetArrayPointer(y);
  fydata = N_VGetArrayPointer(fy);
  rdata  = N_VGetArrayPointer(r);
  zdata  = N_VGetArrayPointer(z);
  ARK_userdata = (FARKUserData) user_data;

  FARK_PSOL(&t, ydata, fydata, rdata, zdata, &gamma, &delta, &lr, 
            ARK_userdata->ipar, ARK_userdata->rpar, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
