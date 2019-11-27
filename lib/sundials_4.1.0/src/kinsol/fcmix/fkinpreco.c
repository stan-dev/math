/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
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
 * This file contains the interfaces between KINSOL and the
 * user-supplied Fortran routines FK_PSET and FK_PSOL.
 *
 * The C function FKINPSet is used to interface between KINSOL and
 * the Fortran user-supplied preconditioner setup routine.
 *
 * The C function FKINPSol is used to interface between KINSOL and
 * the Fortran user-supplied preconditioner solve routine.
 *
 * Note: The use of the generic names FK_PSET and FK_PSOL below.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"
#include "kinsol_impl.h"

#include <kinsol/kinsol_ls.h>

/*------------------------------------------------------------------
  prototype of the user-supplied fortran routine
  ------------------------------------------------------------------*/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_PSET(realtype* uudata,     realtype* uscaledata,
                    realtype* fvaldata,   realtype* fscaledata,
                    int* ier);

extern void FK_PSOL(realtype* uudata,   realtype* uscaledata,
                    realtype* fvaldata, realtype* fscaledata,
                    realtype* vvdata,   int* ier);

#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------
  Function : FKIN_LSSETPREC
  ------------------------------------------------------------------*/
void FKIN_LSSETPREC(int *flag, int *ier)
{
  if ((*flag) == 0) {
    *ier = KINSetPreconditioner(KIN_kinmem, NULL, NULL);
  } else {
    *ier = KINSetPreconditioner(KIN_kinmem, FKINPSet, FKINPSol);
  }

  return;
}

/*------------------------------------------------------------------
  Function : FKIN_SPILSSETPREC -- DEPRECATED
  ------------------------------------------------------------------*/
void FKIN_SPILSSETPREC(int *flag, int *ier)
{ FKIN_LSSETPREC(flag,ier); }

/*------------------------------------------------------------------
  Function : FKINPSet
  ------------------------------------------------------------------
  C function FKINPSet is used to interface between FK_PSET and
  the user-supplied Fortran preconditioner setup routine.
  ------------------------------------------------------------------*/
int FKINPSet(N_Vector uu, N_Vector uscale,
             N_Vector fval, N_Vector fscale,
             void *user_data)
{
  realtype *udata, *uscaledata, *fdata, *fscaledata;
  int ier;

  /* Initialize all pointers to NULL */
  udata = uscaledata = fdata = fscaledata = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  udata      = N_VGetArrayPointer(uu);
  uscaledata = N_VGetArrayPointer(uscale);
  fdata      = N_VGetArrayPointer(fval);
  fscaledata = N_VGetArrayPointer(fscale);

  /* Call user-supplied routine */
  FK_PSET(udata, uscaledata, fdata, fscaledata, &ier);

  return(ier);
}

/*------------------------------------------------------------------
  Function : FKINPSol
  ------------------------------------------------------------------
  C function FKINPSol is used to interface between FK_PSOL and
  the user-supplied Fortran preconditioner solve routine.
  ------------------------------------------------------------------*/
int FKINPSol(N_Vector uu, N_Vector uscale, 
             N_Vector fval, N_Vector fscale, 
             N_Vector vv, void *user_data)
{
  realtype *udata, *uscaledata, *fdata, *fscaledata, *vvdata;
  int ier;

  /* Initialize all pointers to NULL */
  udata = uscaledata = fdata = fscaledata = vvdata = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  udata      = N_VGetArrayPointer(uu);
  uscaledata = N_VGetArrayPointer(uscale);
  fdata      = N_VGetArrayPointer(fval);
  fscaledata = N_VGetArrayPointer(fscale);
  vvdata     = N_VGetArrayPointer(vv);

  /* Call user-supplied routine */
  FK_PSOL(udata, uscaledata, fdata, fscaledata, vvdata, &ier);

  return(ier);
}
