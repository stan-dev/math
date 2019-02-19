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
 * This is the Fortran interface include file for the BBD
 * preconditioner (ARKBBDPRE)
 *--------------------------------------------------------------*/

/*===============================================================
  FARKBBD Interface Package

  The FARKBBD Interface Package is a package of C functions which,
  together with the FARKODE Interface Package, support the use of
  the ARKODE solver and MPI-parallel N_Vector module, along with
  the ARKBBDPRE preconditioner module, for the solution of ODE
  systems in a mixed Fortran/C setting.  We refer the reader to 
  the main ARKode documentation PDF and HTML) for information on 
  usage of the FARKBBD interfce.
  ===============================================================*/

#ifndef _FARKBBD_H
#define _FARKBBD_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* header files  */
/* Definitions of interface function names */
#if defined(SUNDIALS_F77_FUNC)

#define FARK_BBDINIT    SUNDIALS_F77_FUNC(farkbbdinit,   FARKBBDINIT)
#define FARK_BBDREINIT  SUNDIALS_F77_FUNC(farkbbdreinit, FARKBBDREINIT)
#define FARK_BBDOPT     SUNDIALS_F77_FUNC(farkbbdopt,    FARKBBDOPT)
#define FARK_GLOCFN     SUNDIALS_F77_FUNC(farkglocfn,    FARKGLOCFN)
#define FARK_COMMFN     SUNDIALS_F77_FUNC(farkcommfn,    FARKCOMMFN)

#else

#define FARK_BBDINIT    farkbbdinit_
#define FARK_BBDREINIT  farkbbdreinit_
#define FARK_BBDOPT     farkbbdopt_
#define FARK_GLOCFN     farkglocfn_
#define FARK_COMMFN     farkcommfn_

#endif

/* Prototypes of exported functions */
void FARK_BBDINIT(long int *Nloc, long int *mudq,
                  long int *mldq, long int *mu,
                  long int *ml, realtype* dqrely, int *ier);
void FARK_BBDREINIT(long int *mudq, long int *mldq,
                    realtype* dqrely, int *ier);
void FARK_BBDOPT(long int *lenrwbbd, long int *leniwbbd,
                 long int *ngebbd);

/* Prototypes: Functions Called by the ARKBBDPRE Module */
int FARKgloc(long int Nloc, realtype t, N_Vector yloc,
             N_Vector gloc, void *user_data);
int FARKcfn(long int Nloc, realtype t, N_Vector y, 
            void *user_data);

#ifdef __cplusplus
}
#endif

#endif

/*===============================================================
  EOF
  ===============================================================*/
