/* -----------------------------------------------------------------
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
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
 * The C function FIDAPSet is to interface between the IDALS
 * module and the user-supplied preconditioner setup routine FIDAPSET.
 * Note the use of the generic name FIDA_PSET below.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"     /* actual fn. names, prototypes and global vars.*/
#include "ida_impl.h" /* definition of IDAMem type                    */

#include <ida/ida_ls.h>

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_PSET(realtype* t,  realtype* yy,   realtype* yp,
                        realtype* rr, realtype* c_j,  realtype* ewt,
                        realtype* h,  long int* ipar, realtype* rpar,
                        int* ier);
  
  extern void FIDA_PSOL(realtype* t,    realtype* yy,    realtype* yp,
                        realtype* rr,   realtype* r,     realtype* z,
                        realtype* c_j,  realtype* delta, realtype* ewt,
                        long int* ipar, realtype* rpar,  int* ier);

#ifdef __cplusplus
}
#endif

/*************************************************/

/*** DEPRECATED ***/
void FIDA_SPILSSETPREC(int *flag, int *ier)
{ FIDA_LSSETPREC(flag, ier); }

void FIDA_LSSETPREC(int *flag, int *ier)
{
  *ier = 0;

  if (*flag == 0) {

    *ier = IDASetPreconditioner(IDA_idamem, NULL, NULL);

  } else {

    if (F2C_IDA_ewtvec == NULL) {
      F2C_IDA_ewtvec = N_VClone(F2C_IDA_vec);
      if (F2C_IDA_ewtvec == NULL) {
        *ier = -1;
        return;
      }
    }

    *ier = IDASetPreconditioner(IDA_idamem, FIDAPSet, FIDAPSol);
  }

  return;
}

/*************************************************/

int FIDAPSet(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
	     realtype c_j, void *user_data)
{
  realtype *yy_data, *yp_data, *rr_data, *ewtdata;
  realtype h;
  int ier;
  FIDAUserData IDA_userdata;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = ewtdata = NULL;

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

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_PSET(&t, yy_data, yp_data, rr_data, &c_j, ewtdata, &h,
            IDA_userdata->ipar, IDA_userdata->rpar, &ier);

  return(ier);
}

/*************************************************/

int FIDAPSol(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
	     N_Vector rvec, N_Vector zvec,
	     realtype c_j, realtype delta, void *user_data)
{
  realtype *yy_data, *yp_data, *rr_data, *ewtdata, *rdata, *zdata;
  int ier;
  FIDAUserData IDA_userdata;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = ewtdata = zdata = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  IDAGetErrWeights(IDA_idamem, F2C_IDA_ewtvec);

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);
  rdata   = N_VGetArrayPointer(rvec);
  zdata   = N_VGetArrayPointer(zvec);

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_PSOL(&t, yy_data, yp_data, rr_data, rdata, zdata,
	    &c_j, &delta, ewtdata, 
            IDA_userdata->ipar, IDA_userdata->rpar, &ier);

  return(ier);
}
