/*
 * -----------------------------------------------------------------
 * $Revision: 4294 $
 * $Date: 2014-12-15 13:18:40 -0800 (Mon, 15 Dec 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * The FCVROOT module contains the routines necessary to use
 * the rootfinding feature of the CVODE module and to interface
 * with the user-supplied Fortran subroutine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global variables */
#include "fcvroot.h"    /* prototypes of interfaces to CVODE                 */
#include "cvode_impl.h" /* definition of CVodeMem type                       */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_ROOTFN(realtype *, realtype*, realtype*,  /* T, Y, G    */
                         long int*, realtype*,              /* IPAR, RPAR */
                         int *ier);                         /* IER        */
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_ROOTINIT(int *nrtfn, int *ier)
{
  *ier = CVodeRootInit(CV_cvodemem, *nrtfn, (CVRootFn) FCVrootfunc);
  CV_nrtfn = *nrtfn;

  return;
}

/***************************************************************************/

void FCV_ROOTINFO(int *nrtfn, int *info, int *ier)
{
  *ier = CVodeGetRootInfo(CV_cvodemem, info);
  return; 
}

/***************************************************************************/

void FCV_ROOTFREE(void)
{
  CVodeRootInit(CV_cvodemem, 0, NULL);

  return;
}

/***************************************************************************/

int FCVrootfunc(realtype t, N_Vector y, realtype *gout, void *user_data)
{
  int ier;
  realtype *ydata;
  FCVUserData CV_userdata;

  ydata = N_VGetArrayPointer(y);

  CV_userdata = (FCVUserData) user_data;

  FCV_ROOTFN(&t, ydata, gout, CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}

