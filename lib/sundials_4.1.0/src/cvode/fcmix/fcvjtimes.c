/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *    Alan C. Hindmarsh, Radu Serban and Aaron Collier @ LLNL
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
 * The C functions FCVJTSetup and FCVJtimes are to interface 
 * between the CVLS module and the user-supplied 
 * Jacobian-vector product routines FCVJTSETUP and FCVJTIMES. 
 * Note the use of the generic names FCV_JTSETUP and FCV_JTIMES 
 * in the code below.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_ls.h>

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FCV_JTSETUP(realtype *T, realtype *Y, realtype *FY, 
                          realtype *H, long int *IPAR, 
                          realtype *RPAR, int *IER);

  extern void FCV_JTIMES(realtype *V, realtype *JV, realtype *T, 
			  realtype *Y, realtype *FY, realtype *H,
			  long int *IPAR, realtype *RPAR,
			  realtype *WRK, int *IER);

#ifdef __cplusplus
}
#endif

/***************************************************************************/

/* ---DEPRECATED--- */
void FCV_SPILLSSETJAC(int *flag, int *ier)
{ FCV_LSSETJAC(flag, ier); }


void FCV_LSSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVodeSetJacTimes(CV_cvodemem, NULL, NULL);
  } else {
    *ier = CVodeSetJacTimes(CV_cvodemem, FCVJTSetup, FCVJtimes);
  }
}

/***************************************************************************/

/* C function  FCVJTSetup to interface between CVODE and user-supplied
   Fortran routine FCVJTSETUP for preparing a Jacobian * vector product.
   Addresses of t, y, fy and h are passed to FCVJTSETUP,
   using the routine N_VGetArrayPointer from NVECTOR.
   A return flag ier from FCVJTSETUP is returned by FCVJTSetup. */

int FCVJTSetup(realtype t, N_Vector y, N_Vector fy, void *user_data)
{
  realtype *ydata, *fydata;
  realtype h;
  FCVUserData CV_userdata;
  int ier = 0;
  
  CVodeGetLastStep(CV_cvodemem, &h);
  ydata  = N_VGetArrayPointer(y);
  fydata = N_VGetArrayPointer(fy);
  CV_userdata = (FCVUserData) user_data;
 
  FCV_JTSETUP(&t, ydata, fydata, &h, CV_userdata->ipar, 
	      CV_userdata->rpar, &ier);
  return(ier);
}

/* C function  FCVJtimes to interface between CVODE and  user-supplied
   Fortran routine FCVJTIMES for Jacobian * vector product.
   Addresses of v, Jv, t, y, fy, h, and work are passed to FCVJTIMES,
   using the routine N_VGetArrayPointer from NVECTOR.
   A return flag ier from FCVJTIMES is returned by FCVJtimes. */

int FCVJtimes(N_Vector v, N_Vector Jv, realtype t, 
              N_Vector y, N_Vector fy,
              void *user_data, N_Vector work)
{
  realtype *vdata, *Jvdata, *ydata, *fydata, *wkdata;
  realtype h;
  FCVUserData CV_userdata;

  int ier = 0;
  
  CVodeGetLastStep(CV_cvodemem, &h);

  vdata   = N_VGetArrayPointer(v);
  Jvdata  = N_VGetArrayPointer(Jv);
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  wkdata  = N_VGetArrayPointer(work);

  CV_userdata = (FCVUserData) user_data;
 
  FCV_JTIMES (vdata, Jvdata, &t, ydata, fydata, &h, 
              CV_userdata->ipar, CV_userdata->rpar, wkdata, &ier);

  return(ier);
}
