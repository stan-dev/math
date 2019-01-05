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
 * Implementation header file for the ARKBBDPRE module.
 *--------------------------------------------------------------*/

#ifndef _ARKBBDPRE_IMPL_H
#define _ARKBBDPRE_IMPL_H

#include <arkode/arkode_bbdpre.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*---------------------------------------------------------------
 Type: ARKBBDPrecData
---------------------------------------------------------------*/
typedef struct ARKBBDPrecDataRec {

  /* passed by user to ARKBBDPrecAlloc and used by PrecSetup/PrecSolve */
  sunindextype mudq, mldq, mukeep, mlkeep;
  realtype dqrely;
  ARKLocalFn gloc;
  ARKCommFn cfn;

  /* set by ARKBBDPrecSetup and used by ARKBBDPrecSolve */
  SUNMatrix savedJ;
  SUNMatrix savedP;
  SUNLinearSolver LS;
  N_Vector tmp1;
  N_Vector tmp2;
  N_Vector tmp3;
  N_Vector zlocal;
  N_Vector rlocal;

  /* set by ARKBBDPrecAlloc and used by ARKBBDPrecSetup */
  sunindextype n_local;

  /* available for optional output */
  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* pointer to arkode_mem */
  void *arkode_mem;

} *ARKBBDPrecData;


/*---------------------------------------------------------------
 ARKBBDPRE error messages
---------------------------------------------------------------*/

#define MSG_BBD_MEM_NULL    "Integrator memory is NULL."
#define MSG_BBD_LMEM_NULL   "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSG_BBD_MEM_FAIL    "A memory request failed."
#define MSG_BBD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSG_BBD_SUNMAT_FAIL "An error arose from a SUNBandMatrix routine."
#define MSG_BBD_SUNLS_FAIL  "An error arose from a SUNBandLinearSolver routine."
#define MSG_BBD_PMEM_NULL   "BBD peconditioner memory is NULL. ARKBBDPrecInit must be called."
#define MSG_BBD_FUNC_FAILED "The gloc or cfn routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
