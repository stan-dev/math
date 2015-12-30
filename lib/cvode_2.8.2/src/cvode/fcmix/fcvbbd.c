/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
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
 * This module contains the routines necessary to interface with the
 * CVBBDPRE module and user-supplied Fortran routines.
 * The routines here call the generically named routines and provide
 * a standard interface to the C code of the CVBBDPRE package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"              /* actual function names, prototypes, global vars.*/
#include "fcvbbd.h"              /* prototypes of interfaces to CVBBDPRE           */

#include <cvode/cvode_bbdpre.h>  /* prototypes of CVBBDPRE functions and macros    */
#include <cvode/cvode_sptfqmr.h> /* prototypes of CVSPTFQMR interface routines     */
#include <cvode/cvode_spbcgs.h>  /* prototypes of CVSPBCG interface routines       */
#include <cvode/cvode_spgmr.h>   /* prototypes of CVSPGMR interface routines       */

/***************************************************************************/

/* Prototypes of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FCV_GLOCFN(long int*,                        /* NLOC          */
                         realtype*, realtype*, realtype*,  /* T, YLOC, GLOC */
                         long int*, realtype*,             /* IPAR, RPAR    */
                         int *ier);                        /* IER           */

  extern void FCV_COMMFN(long int*,                        /* NLOC          */
                         realtype*, realtype*,             /* T, Y          */
                         long int*, realtype*,             /* IPAR, RPAR    */
                         int *ier);                        /* IER           */

#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_BBDINIT(long int *Nloc, long int *mudq, long int *mldq, long int *mu, long int *ml, 
                 realtype* dqrely, int *ier)
{

  /* 
     First call CVBBDPrecInit to initialize CVBBDPRE module:
     Nloc       is the local vector size
     mudq,mldq  are the half-bandwidths for computing preconditioner blocks
     mu, ml     are the half-bandwidths of the retained preconditioner blocks
     dqrely     is the difference quotient relative increment factor
     FCVgloc    is a pointer to the CVLocalFn function
     FCVcfn     is a pointer to the CVCommFn function 
  */

  *ier = CVBBDPrecInit(CV_cvodemem, *Nloc, *mudq, *mldq, *mu, *ml, 
                       *dqrely, FCVgloc, FCVcfn);

  return; 
}

/***************************************************************************/

void FCV_BBDREINIT(long int *Nloc, long int *mudq, long int *mldq, realtype* dqrely, int *ier)
{
  /* 
     First call CVReInitBBD to re-initialize CVBBDPRE module:
     mudq,mldq   are the half-bandwidths for computing preconditioner blocks
     dqrely      is the difference quotient relative increment factor
     FCVgloc     is a pointer to the CVLocalFn function
     FCVcfn      is a pointer to the CVCommFn function 
  */

  *ier = CVBBDPrecReInit(CV_cvodemem, *mudq, *mldq, *dqrely);
}

/***************************************************************************/

/* C function FCVgloc to interface between CVBBDPRE module and a Fortran 
   subroutine FCVLOCFN. */

int FCVgloc(long int Nloc, realtype t, N_Vector yloc, N_Vector gloc, void *user_data)
{
  int ier;
  realtype *yloc_data, *gloc_data;
  FCVUserData CV_userdata;

  yloc_data = N_VGetArrayPointer(yloc);
  gloc_data = N_VGetArrayPointer(gloc);

  CV_userdata = (FCVUserData) user_data;

  FCV_GLOCFN(&Nloc, &t, yloc_data, gloc_data, 
             CV_userdata->ipar, CV_userdata->rpar,
             &ier);
  return(ier);
}

/***************************************************************************/

/* C function FCVcfn to interface between CVBBDPRE module and a Fortran 
   subroutine FCVCOMMF. */

int FCVcfn(long int Nloc, realtype t, N_Vector y, void *user_data)
{
  int ier;
  realtype *yloc;
  FCVUserData CV_userdata;

  yloc = N_VGetArrayPointer(y);

  CV_userdata = (FCVUserData) user_data;

  FCV_COMMFN(&Nloc, &t, yloc, CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}

/***************************************************************************/

/* C function FCVBBDOPT to access optional outputs from CVBBD_Data */

void FCV_BBDOPT(long int *lenrwbbd, long int *leniwbbd, long int *ngebbd)
{
  CVBBDPrecGetWorkSpace(CV_cvodemem, lenrwbbd, leniwbbd);
  CVBBDPrecGetNumGfnEvals(CV_cvodemem, ngebbd);
}
