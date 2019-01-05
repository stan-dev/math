/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
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
 * Fortran/C interface routines for ARKODE, for the case of a 
 * user-supplied explicit stability routine.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_EXPSTAB(realtype *Y, realtype *T, realtype *HSTAB, 
                           long int *IPAR, realtype *RPAR, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetStabilityFn; see 
   farkode.h for further information */
void FARK_EXPSTABSET(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = ARKStepSetStabilityFn(ARK_arkodemem, NULL, NULL);
  } else {
    *ier = ARKStepSetStabilityFn(ARK_arkodemem, FARKExpStab, 
                                 ARK_arkodemem);
  }
  return;
}

/*=============================================================*/

/* C interface to user-supplied fortran routine FARKEXPSTAB; see 
   farkode.h for further information */
int FARKExpStab(N_Vector y, realtype t, realtype *hstab, void *udata)
{
  int ier = 0;
  realtype *ydata;
  FARKUserData ARK_userdata;

  ydata = N_VGetArrayPointer(y);
  ARK_userdata = (FARKUserData) udata;

  FARK_EXPSTAB(ydata, &t, hstab, ARK_userdata->ipar, 
               ARK_userdata->rpar, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
