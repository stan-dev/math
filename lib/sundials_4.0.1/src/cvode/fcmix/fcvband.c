/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE/CVLS, for the case of 
 * a user-supplied Jacobian approximation routine.                
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_ls.h>
#include <sunmatrix/sunmatrix_band.h>

/******************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_BJAC(long int *N, long int *MU, long int *ML,
                       long int *EBAND, realtype *T, realtype *Y,
                       realtype *FY, realtype *BJAC, realtype *H,
                       long int *IPAR, realtype *RPAR, realtype *V1, 
                       realtype *V2, realtype *V3, int *IER);
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_BANDSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVodeSetJacFn(CV_cvodemem, NULL);
  } else {
    *ier = CVodeSetJacFn(CV_cvodemem, FCVBandJac);
  }
}

/***************************************************************************/

/* C function CVBandJac interfaces between CVODE and a Fortran subroutine
   FCVBJAC for solution of a linear system with band Jacobian approximation.
   Addresses of arguments are passed to FCVBJAC, using the accessor routines
   from the SUNBandMatrix and N_Vector modules.
   The address passed for J is that of the element in column 0 with row 
   index -mupper.  An extended bandwith equal to (J->smu) + mlower + 1 is
   passed as the column dimension of the corresponding array. */

int FCVBandJac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void *user_data, N_Vector vtemp1, N_Vector vtemp2,
               N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;
  long int N, mupper, mlower, smu, eband;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  N = SUNBandMatrix_Columns(J);
  mupper = SUNBandMatrix_UpperBandwidth(J);
  mlower = SUNBandMatrix_LowerBandwidth(J);
  smu = SUNBandMatrix_StoredUpperBandwidth(J);
  eband = smu + mlower + 1;
  jacdata = SUNBandMatrix_Column(J,0) - mupper;

  CV_userdata = (FCVUserData) user_data;

  FCV_BJAC(&N, &mupper, &mlower, &eband, &t, ydata, fydata,
           jacdata, &h, CV_userdata->ipar, CV_userdata->rpar,
           v1data, v2data, v3data, &ier);

  return(ier);
}
