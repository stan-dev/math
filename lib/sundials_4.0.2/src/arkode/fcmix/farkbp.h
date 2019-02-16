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
 * This is the Fortran interface include file for the BAND
 * preconditioner (ARKBANDPRE).
 *--------------------------------------------------------------*/

/*===============================================================
  FARKBP Interface Package
  
  The FARKBP Interface Package is a package of C functions which,
  together with the FARKODE Interface Package, support the use of
  the ARKODE solver and serial, OpenMP or PThreads vector module
  with the ARKBANDPRE preconditioner module, for the solution of
  ODE systems in a mixed Fortran/C setting.  We refer the reader to 
  the main ARKode documentation PDF and HTML) for information on 
  usage of the FARKBBD interfce.
  ===============================================================*/

#ifndef _FARKBP_H
#define _FARKBP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* header files  */
/* Definitions of interface function names */
#if defined(SUNDIALS_F77_FUNC)

#define FARK_BPINIT    SUNDIALS_F77_FUNC(farkbpinit, FARKBPINIT)
#define FARK_BPOPT     SUNDIALS_F77_FUNC(farkbpopt,  FARKBPOPT)

#else

#define FARK_BPINIT    farkbpinit_
#define FARK_BPOPT     farkbpopt_

#endif

/* Prototypes of exported function */
void FARK_BPINIT(long int *N,
                 long int *mu,
                 long int *ml,
                 int *ier);
void FARK_BPOPT(long int *lenrwbp,
                long int *leniwbp,
                long int *nfebp);

#ifdef __cplusplus
}
#endif

#endif

/*===============================================================
  EOF
  ===============================================================*/
