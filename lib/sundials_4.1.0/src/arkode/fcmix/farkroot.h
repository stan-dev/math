/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the Fortran interface include file for the rootfinding
 * feature of ARKODE.
 *--------------------------------------------------------------*/

/*===============================================================
  FARKROOT Interface Package
  
  The FARKROOT interface package allows programs written in 
  FORTRAN to use the rootfinding features of the ARKODE solver 
  module.  We refer the reader to the main ARKode documentation 
  (PDF and HTML) for usage information.
  ===============================================================*/

#ifndef _FARKROOT_H
#define _FARKROOT_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */
#if defined(SUNDIALS_F77_FUNC)

#define FARK_ROOTINIT SUNDIALS_F77_FUNC(farkrootinit, FARKROOTINIT)
#define FARK_ROOTINFO SUNDIALS_F77_FUNC(farkrootinfo, FARKROOTINFO)
#define FARK_ROOTFREE SUNDIALS_F77_FUNC(farkrootfree, FARKROOTFREE)
#define FARK_ROOTFN   SUNDIALS_F77_FUNC(farkrootfn,   FARKROOTFN)

#else

#define FARK_ROOTINIT farkrootinit_
#define FARK_ROOTINFO farkrootinfo_
#define FARK_ROOTFREE farkrootfree_
#define FARK_ROOTFN   farkrootfn_

#endif

/* Prototypes of exported function */
void FARK_ROOTINIT(int *nrtfn, int *ier);
void FARK_ROOTINFO(int *nrtfn, int *info, int *ier);
void FARK_ROOTFREE(void);

/* Prototype of function called by ARKODE module */
int FARKrootfunc(realtype t, N_Vector y, 
                 realtype *gout, void *user_data);

#ifdef __cplusplus
}
#endif

#endif

/*===============================================================
  EOF
  ===============================================================*/
