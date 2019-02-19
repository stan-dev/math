/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 *                David J. Gardner @ LLNL
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
 * This module contains the routines necessary to interface with
 * the KINBBDPRE module and user-supplied Fortran routines. Generic
 * names are used (e.g. FK_COMMFN). The routines here call the
 * generically named routines and provide a standard interface to
 * the C code of the KINBBDPRE package.
 * ----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"               /* standard interfaces and global variables     */
#include "fkinbbd.h"               /* prototypes of interfaces to KINBBDPRE        */

#include <kinsol/kinsol_bbdpre.h>  /* prototypes of KINBBDPRE functions and macros */

/*
 * ----------------------------------------------------------------
 * private constants
 * ----------------------------------------------------------------
 */

#define ZERO RCONST(0.0)

/*
 * ----------------------------------------------------------------
 * prototypes of the user-supplied fortran routines
 * ----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_LOCFN(long int* NLOC, realtype* ULOC, realtype* GLOC, int* IER);
extern void FK_COMMFN(long int* NLOC, realtype* ULOC, int* IER);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDINIT
 * ----------------------------------------------------------------
 */

void FKIN_BBDINIT(long int *nlocal, long int *mudq, long int *mldq,
		  long int *mu, long int *ml, int *ier)
{
  *ier = KINBBDPrecInit(KIN_kinmem, *nlocal, *mudq, *mldq, *mu, *ml, ZERO,
                        (KINBBDLocalFn) FKINgloc, (KINBBDCommFn) FKINgcomm);

  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINgloc
 * ----------------------------------------------------------------
 * C function FKINgloc is the interface between the KINBBDPRE
 * module and the Fortran subroutine FK_LOCFN.
 * ----------------------------------------------------------------
 */

int FKINgloc(long int Nloc, N_Vector uu, N_Vector gval, void *user_data)
{
  realtype *uloc, *gloc;
  int ier;

  /* Initialize all pointers to NULL */
  uloc = gloc = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  uloc = N_VGetArrayPointer(uu);
  gloc = N_VGetArrayPointer(gval);

  /* Call user-supplied routine */
  FK_LOCFN(&Nloc, uloc, gloc, &ier);

  return(ier);
}

/*
 * ----------------------------------------------------------------
 * Function : FKINgcomm
 * ----------------------------------------------------------------
 * C function FKINgcomm is the interface between the KINBBDPRE
 * module and the Fortran subroutine FK_COMMFN.
 * ----------------------------------------------------------------
 */

int FKINgcomm(long int Nloc, N_Vector uu, void *user_data)
{
  realtype *uloc;
  int ier;

  /* Initialize all pointers to NULL */
  uloc = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  uloc = N_VGetArrayPointer(uu);
  
  /* Call user-supplied routine */
  FK_COMMFN(&Nloc, uloc, &ier);

  return(ier);
}

/*
 * ----------------------------------------------------------------
 * Function : FKIN_BBDOPT
 * ----------------------------------------------------------------
 * C function FKIN_BBDOPT is used to access optional outputs 
 * realated to the BBD preconditioner.
 * ----------------------------------------------------------------
 */

void FKIN_BBDOPT(long int *lenrpw, long int *lenipw, long int *nge)
{
  KINBBDPrecGetWorkSpace(KIN_kinmem, lenrpw, lenipw);
  KINBBDPrecGetNumGfnEvals(KIN_kinmem, nge);

  return;
}

