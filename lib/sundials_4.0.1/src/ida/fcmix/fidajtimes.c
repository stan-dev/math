/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Aaron Collier and Radu Serban @ LLNL
 *-----------------------------------------------------------------
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
 *-----------------------------------------------------------------
 * The C functions FIDAJTSetup and FIDAJtimes are to interface 
 * between the IDALS module and the user-supplied 
 * Jacobian-vector product routines FIDAJTSETUP and FIDAJTIMES. 
 * Note the use of the generic names FIDA_JTSETUP and FIDA_JTIMES 
 * below.
 *-----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"     /* actual fn. names, prototypes and global vars.*/
#include "ida_impl.h" /* definition of IDAMem type                    */

#include <ida/ida_ls.h>

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_JTSETUP(realtype *T, realtype *Y, realtype *YP,
                           realtype *R, realtype *CJ, realtype *EWT, 
                           realtype *H, long int *IPAR,
                           realtype *RPAR, int *IER);

  extern void FIDA_JTIMES(realtype *T, realtype *Y, realtype *YP,
                          realtype *R, realtype *V, realtype *FJV,
                          realtype *CJ, realtype *EWT, realtype *H,
                          long int *IPAR, realtype *RPAR,
                          realtype *WK1, realtype *WK2, int *IER);

#ifdef __cplusplus
}
#endif

/*************************************************/

/*** DEPRECATED ***/
void FIDA_SPILSSETJAC(int *flag, int *ier)
{ FIDA_LSSETJAC(flag, ier); }


/* Fortran interface to C routine IDASetJacTimes; see 
   fida.h for further information */
void FIDA_LSSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = IDASetJacTimes(IDA_idamem, NULL, NULL);
  } else {
    if (F2C_IDA_ewtvec == NULL) {
      F2C_IDA_ewtvec = N_VClone(F2C_IDA_vec);
      if (F2C_IDA_ewtvec == NULL) {
        *ier = -1;
        return;
      }
    }
    *ier = IDASetJacTimes(IDA_idamem, FIDAJTSetup, FIDAJtimes);
  }
  return;
}

/*************************************************/

/* C interface to user-supplied Fortran routine FIDAJTSETUP; see
   fida.h for further information */
int FIDAJTSetup(realtype t, N_Vector y, N_Vector yp, 
                N_Vector r, realtype cj, void *user_data)
{
  realtype *ydata, *ypdata, *rdata, *ewtdata;
  realtype h;
  FIDAUserData IDA_userdata;
  int ier = 0;
  
  /* Initialize all pointers to NULL */
  ydata = ypdata = rdata = ewtdata = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;
  
  IDAGetErrWeights(IDA_idamem, F2C_IDA_ewtvec);
  IDAGetLastStep(IDA_idamem, &h);
  ydata   = N_VGetArrayPointer(y);
  ypdata  = N_VGetArrayPointer(yp);
  rdata   = N_VGetArrayPointer(r);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);

  IDA_userdata = (FIDAUserData) user_data;
 
  /* Call user-supplied routine */
  FIDA_JTSETUP(&t, ydata, ypdata, rdata, &cj, ewtdata, &h,
               IDA_userdata->ipar, IDA_userdata->rpar, &ier);
  return(ier);
}

/* C interface to user-supplied Fortran routine FIDAJTIMES; see
   fida.h for further information */
int FIDAJtimes(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
	       N_Vector v, N_Vector Jv, realtype c_j,
               void *user_data, N_Vector vtemp1, N_Vector vtemp2)
{
  realtype *yy_data, *yp_data, *rr_data, *vdata, *Jvdata, *ewtdata;
  realtype *v1data, *v2data;
  realtype h;
  FIDAUserData IDA_userdata;
  int ier;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = vdata = Jvdata = ewtdata = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  IDAGetErrWeights(IDA_idamem, F2C_IDA_ewtvec);
  IDAGetLastStep(IDA_idamem, &h);

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);
  vdata   = N_VGetArrayPointer(v);
  Jvdata  = N_VGetArrayPointer(Jv);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_JTIMES(&t, yy_data, yp_data, rr_data, vdata, Jvdata,
	      &c_j, ewtdata, &h, 
              IDA_userdata->ipar, IDA_userdata->rpar,
              v1data, v2data, &ier);

  return(ier);
}
