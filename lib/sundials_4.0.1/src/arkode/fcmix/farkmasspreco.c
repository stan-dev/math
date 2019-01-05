/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 * The C function FARKPSet is to interface between the ARKLSMASS 
 * module and the user-supplied mass matrix preconditioner  
 * setup/solve routines FARKPSET and FARKPSOL. Note the use of 
 * the generic names FARK_PSET and  FARK_PSOL in the code below.
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

  extern void FARK_MASSPSET(realtype *T, long int *IPAR, 
                            realtype *RPAR, int *IER);
  extern void FARK_MASSPSOL(realtype *T, realtype *R, realtype *Z, 
                            realtype *DELTA, int *LR, long int *IPAR, 
                            realtype *RPAR, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* ---DEPRECATED---
   Fortran interface to C routine ARKStepSetMassPreconditioner; see 
   farkode.h for further details */
void FARK_SPILSSETMASSPREC(int *flag, int *ier)
{ FARK_LSSETMASSPREC(flag, ier); }

/* Fortran interface to C routine ARKStepSetMassPreconditioner; see 
   farkode.h for further details */
void FARK_LSSETMASSPREC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = ARKStepSetMassPreconditioner(ARK_arkodemem, NULL, NULL);
  } else {
    *ier = ARKStepSetMassPreconditioner(ARK_arkodemem, 
                                        FARKMassPSet, FARKMassPSol);
  }
  return;
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKMASSPSET; see 
   farkode.h for further details */
int FARKMassPSet(realtype t, void *user_data)
{
  int ier = 0;
  FARKUserData ARK_userdata;
  ARK_userdata = (FARKUserData) user_data;
  FARK_MASSPSET(&t, ARK_userdata->ipar, ARK_userdata->rpar, &ier);
  return(ier);
}


/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKMASSPSOL; see 
   farkode.h for further details */
int FARKMassPSol(realtype t, N_Vector r, N_Vector z, realtype delta,
                 int lr, void *user_data)
{
  int ier = 0;
  realtype *rdata, *zdata;
  FARKUserData ARK_userdata;

  rdata  = N_VGetArrayPointer(r);
  zdata  = N_VGetArrayPointer(z);
  ARK_userdata = (FARKUserData) user_data;

  FARK_MASSPSOL(&t, rdata, zdata, &delta, &lr, ARK_userdata->ipar, 
                ARK_userdata->rpar, &ier);
  return(ier);
}


/*===============================================================
   EOF
===============================================================*/
