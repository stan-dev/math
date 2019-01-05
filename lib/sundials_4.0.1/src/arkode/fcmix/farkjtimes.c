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
 * The C functions FARKJTSetup and FARKJtimes are to interface 
 * between the ARKLS module and the user-supplied Jacobian-vector 
 * product routines FARKJTSETUP and FARKJTIMES. Note use of the
 * generic names FARK_JTSETUP and FARK_JTIMES in the code below.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_arkstep.h>

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_JTSETUP(realtype *T, realtype *Y, realtype *FY, 
                           realtype *H, long int *IPAR, 
                           realtype *RPAR, int *IER);

  extern void FARK_JTIMES(realtype *V, realtype *JV, realtype *T, 
                          realtype *Y, realtype *FY, realtype *H,
                          long int *IPAR, realtype *RPAR,
                          realtype *WRK, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* ---DEPRECATED---
   Fortran interface to C routine ARKStepSetJacTimes; see 
   farkode.h for further information */
void FARK_SPILSSETJAC(int *flag, int *ier)
{ FARK_LSSETJAC(flag,ier); }

/* Fortran interface to C routine ARKStepSetJacTimes; see 
   farkode.h for further information */
void FARK_LSSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = ARKStepSetJacTimes(ARK_arkodemem, NULL, NULL);
  } else {
    *ier = ARKStepSetJacTimes(ARK_arkodemem, FARKJTSetup, FARKJtimes);
  }
  return;
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKJTSETUP; see
   farkode.h for further information */
int FARKJTSetup(realtype t, N_Vector y, N_Vector fy, void *user_data)
{
  realtype *ydata, *fydata;
  realtype h;
  FARKUserData ARK_userdata;
  int ier = 0;
  
  /* Initialize all pointers to NULL */
  ydata = fydata = NULL;
  
  ARKStepGetLastStep(ARK_arkodemem, &h);
  ydata  = N_VGetArrayPointer(y);
  fydata = N_VGetArrayPointer(fy);
  ARK_userdata = (FARKUserData) user_data;
 
  FARK_JTSETUP(&t, ydata, fydata, &h, ARK_userdata->ipar, 
              ARK_userdata->rpar, &ier);
  return(ier);
}

/* C interface to user-supplied Fortran routine FARKJTIMES; see
   farkode.h for further information */
int FARKJtimes(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
               N_Vector fy, void *user_data, N_Vector work)
{
  realtype *vdata, *Jvdata, *ydata, *fydata, *wkdata;
  realtype h;
  FARKUserData ARK_userdata;
  int ier = 0;
  
  /* Initialize all pointers to NULL */
  vdata = Jvdata = ydata = fydata = wkdata = NULL;

  ARKStepGetLastStep(ARK_arkodemem, &h);

  vdata  = N_VGetArrayPointer(v);
  Jvdata = N_VGetArrayPointer(Jv);
  ydata  = N_VGetArrayPointer(y);
  fydata = N_VGetArrayPointer(fy);
  wkdata = N_VGetArrayPointer(work);

  ARK_userdata = (FARKUserData) user_data;
 
  FARK_JTIMES(vdata, Jvdata, &t, ydata, fydata, &h, ARK_userdata->ipar, 
              ARK_userdata->rpar, wkdata, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
