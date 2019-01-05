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
 * user-supplied error weight calculation routine.
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

  extern void FARK_EWT(realtype *Y, realtype *EWT,
                       long int *IPAR, realtype *RPAR,
                       int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKStepWFtolerances; see
   farkode.h for further information */
void FARK_EWTSET(int *flag, int *ier)
{
  if (*flag != 0) {
    *ier = ARKStepWFtolerances(ARK_arkodemem, FARKEwt);
  }
  return;
}

/*=============================================================*/

/* C interface to user-supplied fortran routine FARKEWT; see
   farkode.h for further information */
int FARKEwt(N_Vector y, N_Vector ewt, void *user_data)
{
  int ier = 0;
  realtype *ydata, *ewtdata;
  FARKUserData ARK_userdata;

  ydata  = N_VGetArrayPointer(y);
  ewtdata = N_VGetArrayPointer(ewt);
  ARK_userdata = (FARKUserData) user_data;

  FARK_EWT(ydata, ewtdata, ARK_userdata->ipar,
           ARK_userdata->rpar, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
