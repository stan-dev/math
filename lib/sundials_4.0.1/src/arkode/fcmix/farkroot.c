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
 * The FARKROOT module contains the routines necessary to use
 * the rootfinding feature of the ARKODE module and to interface
 * with the user-supplied Fortran subroutine.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "farkroot.h"
#include "arkode_impl.h"

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FARK_ROOTFN(realtype *T, realtype *Y,
                          realtype *G, long int *IPAR,
                          realtype *RPAR, int *ier);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKStepRootInit; see farkroot.h
   for further information. */
void FARK_ROOTINIT(int *nrtfn, int *ier)
{
  *ier = ARKStepRootInit(ARK_arkodemem, *nrtfn,
                         (ARKRootFn) FARKrootfunc);
  ARK_nrtfn = *nrtfn;
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepGetRootInfo; see
   farkroot.h for further information. */
void FARK_ROOTINFO(int *nrtfn, int *info, int *ier)
{
  *ier = ARKStepGetRootInfo(ARK_arkodemem, info);
  return;
}

/*=============================================================*/

/* Fortran interface to C routine ARKStepRootInit, used to free
   existing memory resources; see farkroot.h for further
   information. */
void FARK_ROOTFREE(void)
{
  ARKStepRootInit(ARK_arkodemem, 0, NULL);
  return;
}

/*=============================================================*/

/* C interface to user-supplied routine FARKROOTFN; see
   farkroot.h for further information. */
int FARKrootfunc(realtype t, N_Vector y,
                 realtype *gout, void *user_data)
{
  int ier;
  realtype *ydata;
  FARKUserData ARK_userdata;

  ydata = N_VGetArrayPointer(y);
  ARK_userdata = (FARKUserData) user_data;
  FARK_ROOTFN(&t, ydata, gout, ARK_userdata->ipar,
              ARK_userdata->rpar, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
