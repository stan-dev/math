/*
 * -----------------------------------------------------------------
 * $Revision: 4294 $
 * $Date: 2014-12-15 13:18:40 -0800 (Mon, 15 Dec 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh, Radu Serban and
 *                Aaron Collier @ LLNL
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
 * The C function FCVPSet is to interface between the CVSP*
 * module and the user-supplied preconditioner setup routine FCVPSET.
 * Note the use of the generic name FCV_PSET below.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_spils.h>

/*********************************************************************/

/* Prototype of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FCV_PSET(realtype*, realtype*, realtype*,  /* T, Y, FY */
                       booleantype*, booleantype*,       /* JOK, JCUR */
                       realtype*, realtype*,             /* GAMMA, H */
                       long int*, realtype*,             /* IPAR, RPAR */
                       realtype*, realtype*, realtype*,  /* W1, W2, W3 */
                       int*);                            /* IER */

  extern void FCV_PSOL(realtype*, realtype*, realtype*,  /* T, Y, FY */
                       realtype*, realtype*,             /* R, Z */
                       realtype*, realtype*,             /* GAMMA, DELTA */
                       int*,                             /* LR */
                       long int*, realtype*,             /* IPAR, RPAR */
                       realtype*,                        /* WRK */
                       int*);                            /* IER */

#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_SPILSSETPREC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVSpilsSetPreconditioner(CV_cvodemem, NULL, NULL);
  } else {
    *ier = CVSpilsSetPreconditioner(CV_cvodemem, FCVPSet, FCVPSol);
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
            void *user_data,
            N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier = 0;
  realtype *ydata, *fydata, *v1data, *v2data, *v3data;
  realtype h;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  CV_userdata = (FCVUserData) user_data;

  FCV_PSET(&t, ydata, fydata, &jok, jcurPtr, &gamma, &h,
           CV_userdata->ipar, CV_userdata->rpar,
           v1data, v2data, v3data, &ier);

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
            int lr, void *user_data, N_Vector vtemp)
{
  int ier = 0;
  realtype *ydata, *fydata, *vtdata, *rdata, *zdata;
  FCVUserData CV_userdata;

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  vtdata  = N_VGetArrayPointer(vtemp);
  rdata   = N_VGetArrayPointer(r);
  zdata   = N_VGetArrayPointer(z);

  CV_userdata = (FCVUserData) user_data;

  FCV_PSOL(&t, ydata, fydata, rdata, zdata, &gamma, &delta, &lr, 
           CV_userdata->ipar, CV_userdata->rpar, vtdata, &ier);

  return(ier);
}
