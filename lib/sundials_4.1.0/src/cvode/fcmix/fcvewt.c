/*
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
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
 * Fortran/C interface routines for CVODE, for the case of a 
 * user-supplied error weight calculation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"           /* actual fn. names, prototypes and global vars.  */
#include "cvode_impl.h"       /* definition of CVodeMem type                    */

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_EWT(realtype *Y, realtype *EWT,
                      long int *IPAR, realtype *RPAR,
                      int *IER);
#ifdef __cplusplus
}
#endif

/***************************************************************************/

/* 
 * User-callable function to interface to CVodeSetEwtFn.
 */

void FCV_EWTSET(int *flag, int *ier)
{
  if (*flag != 0) {
    *ier = CVodeWFtolerances(CV_cvodemem, FCVEwtSet);
  }
}

/***************************************************************************/

/* 
 * C function to interface between CVODE and a Fortran subroutine FCVEWT.
 */

int FCVEwtSet(N_Vector y, N_Vector ewt, void *user_data)
{
  int ier = 0;
  realtype *ydata, *ewtdata;
  FCVUserData CV_userdata;

  ydata  = N_VGetArrayPointer(y);
  ewtdata = N_VGetArrayPointer(ewt);

  CV_userdata = (FCVUserData) user_data;

  FCV_EWT(ydata, ewtdata, CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}
