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
 * Implementation header file for ARKode's temporal interpolation
 * utility.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_INTERP_IMPL_H
#define _ARKODE_INTERP_IMPL_H

#include <stdarg.h>
#include <arkode/arkode.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARKode temporal interpolation constants
  ===============================================================*/

#define QDENSE_DEF  3       /* default dense output order */
#define ARK_INTERP_LRW  2   /* real workspace size */
#define ARK_INTERP_LIW  5   /* int/ptr workspace size */


/*===============================================================
  ARKode Temporal Interpolation Data Structure
  ===============================================================*/

/*---------------------------------------------------------------
 Types : struct ARKodeInterpMemRec, ARKodeInterpMem
-----------------------------------------------------------------
 The type ARKodeInterpMem is type pointer to struct
 ARKodeInterpMemRec.  This structure contains fields to
 perform temporal interpolation.
---------------------------------------------------------------*/
typedef struct ARKodeInterpMemRec {

  N_Vector fold;    /* f(t,y) at beginning of last successful step */
  N_Vector fnew;    /* f(t,y) at end of last successful step       */
  N_Vector yold;    /* y at beginning of last successful step      */
  N_Vector ynew;    /* y at end of last successful step            */
  N_Vector fa;      /* f(t,y) used in higher-order interpolation   */
  N_Vector fb;      /* f(t,y) used in higher-order interpolation   */
  realtype told;    /* t at beginning of last successful step      */
  realtype tnew;    /* t at end of last successful step            */
  realtype t_fa;    /* t when fa was last evaluated                */
  realtype t_fb;    /* t when fb was last evaluated                */
  realtype h;       /* last successful step size                   */
  int order;        /* interpolation order                         */

} *ARKodeInterpMem;


/*===============================================================
  ARKode Temporal Interpolation Routines
===============================================================*/

ARKodeInterpMem arkInterpCreate(void* arkode_mem);
int arkInterpResize(void* arkode_mem, ARKodeInterpMem interp_mem,
                    ARKVecResizeFn resize, void *resize_data,
                    sunindextype lrw_diff, sunindextype liw_diff,
                    N_Vector tmpl);
void arkInterpFree(ARKodeInterpMem *interp_mem);
void arkPrintInterpMem(ARKodeInterpMem interp_mem, FILE *outfile);
int arkInterpInit(void* arkode_mem, ARKodeInterpMem interp_mem,
                  realtype tnew);
int arkInterpUpdate(void* arkode_mem, ARKodeInterpMem interp_mem,
                    realtype tnew, booleantype forceRHS);
int arkInterpEvaluate(void* arkode_mem, ARKodeInterpMem interp_mem,
                      realtype tau, int d, int order, N_Vector yout);


#ifdef __cplusplus
}
#endif

#endif
