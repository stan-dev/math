/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
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
  extern void FCV_ROOTFN(realtype *T, realtype *Y, realtype *G,
                         long int *IPAR, realtype *RPAR,
                         int *ier);
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

