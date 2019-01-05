/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * Implementation header file for the ARKBANDPRE module.
 *--------------------------------------------------------------*/

#ifndef _ARKBANDPRE_IMPL_H
#define _ARKBANDPRE_IMPL_H

#include <arkode/arkode_bandpre.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
 Type: ARKBandPrecData
---------------------------------------------------------------*/

typedef struct ARKBandPrecDataRec {

  /* Data set by user in ARKBandPrecInit */
  sunindextype N;
  sunindextype ml, mu;

  /* Data set by ARKBandPrecSetup */
  SUNMatrix savedJ;
  SUNMatrix savedP;
  SUNLinearSolver LS;
  N_Vector tmp1;
  N_Vector tmp2;

  /* Rhs calls */
  long int nfeBP;

  /* Pointer to arkode_mem */
  void *arkode_mem;

} *ARKBandPrecData;


/*---------------------------------------------------------------
 ARKBANDPRE error messages
---------------------------------------------------------------*/

#define MSG_BP_MEM_NULL       "Integrator memory is NULL."
#define MSG_BP_LMEM_NULL      "Linear solver memory is NULL. The SPILS interface must be attached."
#define MSG_BP_MEM_FAIL       "A memory request failed."
#define MSG_BP_BAD_NVECTOR    "A required vector operation is not implemented."
#define MSG_BP_SUNMAT_FAIL    "An error arose from a SUNBandMatrix routine."
#define MSG_BP_SUNLS_FAIL     "An error arose from a SUNBandLinearSolver routine."
#define MSG_BP_PMEM_NULL      "Band preconditioner memory is NULL. ARKBandPrecInit must be called."
#define MSG_BP_RHSFUNC_FAILED "The right-hand side routine failed in an unrecoverable manner."


#ifdef __cplusplus
}
#endif

#endif
