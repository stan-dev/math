/*
 * -----------------------------------------------------------------
 * $Revision: 4294 $
 * $Date: 2014-12-15 13:18:40 -0800 (Mon, 15 Dec 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
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
 * Fortran/C interface routines for CVODE/CVLAPACK, for the case
 * of a user-supplied band Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_lapack.h>

/***************************************************************************/

/* Prototype of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FCV_BJAC(long int*, long int*, long int*, long int*,    /* N,MU,ML,EBAND */
                       realtype*, realtype*, realtype*,  /* T, Y, FY         */
                       realtype*,                        /* LBJAC            */
                       realtype*,                        /* H                */
                       long int*, realtype*,             /* IPAR, RPAR       */
                       realtype*, realtype*, realtype*,  /* V1, V2, V3       */
                       int*);                            /* IER              */

#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_LAPACKBANDSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVDlsSetBandJacFn(CV_cvodemem, NULL);
  } else {
    *ier = CVDlsSetBandJacFn(CV_cvodemem, FCVLapackBandJac);
  }
}

/***************************************************************************/

/* The C function FCVLapackBandJac interfaces between CVODE and a 
 * Fortran subroutine FCVBJAC for the solution of a linear system using
 * Lapack with band Jacobian approximation.
 * Addresses of arguments are passed to FCVBJAC, using the macro 
 * BAND_COL and the routine N_VGetArrayPointer from NVECTOR.
 * The address passed for J is that of the element in column 0 with row 
 * index -mupper.  An extended bandwith equal to (J->smu) + mlower + 1 is
 * passed as the column dimension of the corresponding array.
 * Auxiliary data is assumed to be communicated by Common. 
 */

int FCVLapackBandJac(long int N, long int mupper, long int mlower,
                     realtype t, N_Vector y, N_Vector fy, 
                     DlsMat J, void *user_data,
                     N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;
  long int eband;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  eband = (J->s_mu) + mlower + 1;
  jacdata = BAND_COL(J,0) - mupper;

  CV_userdata = (FCVUserData) user_data;

  FCV_BJAC(&N, &mupper, &mlower, &eband, &t, ydata, fydata, jacdata, &h,
           CV_userdata->ipar, CV_userdata->rpar, v1data, v2data, v3data, &ier);

  return(ier);
}
