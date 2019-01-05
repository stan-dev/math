/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *     Alan C. Hindmarsh, Radu Serban and Aaron Collier @ LLNL
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
 * The C functions FCVPSet and FCVPSol are to interface between the 
 * CVLS module and the user-supplied preconditioner setup/solve
 * routines FCVPSET and FCVPSOL. Note the use of the generic names
 * FCV_PSET and FCV_PSOL below.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_ls.h>

/*********************************************************************/

/* Prototype of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FCV_PSET(realtype *T, realtype *Y, realtype *FY,
			booleantype *JOK, booleantype *JCUR,
			realtype *GAMMA, realtype *H,
			long int *IPAR, realtype *RPAR, int *IER);

  extern void FCV_PSOL(realtype *T, realtype *Y, realtype *FY,
			realtype *R, realtype *Z, 
			realtype *GAMMA, realtype *DELTA,
			int *LR, long int *IPAR, realtype *RPAR, 
                        int *IER);

#ifdef __cplusplus
}
#endif

/***************************************************************************/

/* ---DEPRECATED--- */
void FCV_SPILSSETPREC(int *flag, int *ier)
{ FCV_LSSETPREC(flag, ier); }

void FCV_LSSETPREC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVodeSetPreconditioner(CV_cvodemem, NULL, NULL);
  } else {
    *ier = CVodeSetPreconditioner(CV_cvodemem, FCVPSet, FCVPSol);
  }
}

/***************************************************************************/

/* C function FCVPSet to interface between CVODE and a Fortran subroutine
   FCVPSET for setup of a Krylov preconditioner.
   Addresses of t, y, fy, jok, gamma, h, vtemp1, vtemp2, vtemp3, and the 
   address jcurPtr are passed to FCVPSET, using the routine
   N_VGetArrayPointer from NVECTOR.
   A return flag ier from FCVPSET is returned by FCVPSet.
   Auxiliary data is assumed to be communicated by common blocks. */

int FCVPSet(realtype t, N_Vector y, N_Vector fy, booleantype jok,
            booleantype *jcurPtr, realtype gamma,
            void *user_data)
{
  int ier = 0;
  realtype *ydata, *fydata;
  realtype h;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);

  CV_userdata = (FCVUserData) user_data;

  FCV_PSET(&t, ydata, fydata, &jok, jcurPtr, &gamma, &h,
           CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}

/***************************************************************************/

/* C function FCVPSol to interface between CVODE and a Fortran subroutine
   FCVPSOL for solution of a Krylov preconditioner.
   Addresses of t, y, fy, gamma, delta, lr, vtemp, r, and z are
   passed to FCVPSOL, using the routine N_VGetArrayPointer from NVECTOR.
   A return flag ier from FCVPSOL is returned by FCVPSol.
   Auxiliary data is assumed to be communicated by Common blocks. */

int FCVPSol(realtype t, N_Vector y, N_Vector fy, 
            N_Vector r, N_Vector z,
            realtype gamma, realtype delta,
            int lr, void *user_data)
{
  int ier = 0;
  realtype *ydata, *fydata, *rdata, *zdata;
  FCVUserData CV_userdata;

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  rdata   = N_VGetArrayPointer(r);
  zdata   = N_VGetArrayPointer(z);

  CV_userdata = (FCVUserData) user_data;

  FCV_PSOL(&t, ydata, fydata, rdata, zdata, &gamma, &delta, &lr, 
           CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}
