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
 * Header file for Butcher table structures in ARKode.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_BUTCHER_H
#define _ARKODE_BUTCHER_H

#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*---------------------------------------------------------------
 Types : struct ARKodeButcherTableMem, ARKodeButcherTable
---------------------------------------------------------------*/
typedef struct ARKodeButcherTableMem {

  int q;           /* method order of accuracy       */
  int p;           /* embedding order of accuracy    */
  int stages;      /* number of stages               */
  realtype **A;    /* Butcher table coefficients     */
  realtype *c;     /* canopy node coefficients       */
  realtype *b;     /* root node coefficients         */
  realtype *d;     /* embedding coefficients         */

} *ARKodeButcherTable;


/* Utility routines to allocate/free/output Butcher table structures */
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_Alloc(int stages,
                                                            booleantype embedded);
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_Create(int s, int q,
                                                             int p,
                                                             realtype *c,
                                                             realtype *A,
                                                             realtype *b,
                                                             realtype *d);
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_Copy(ARKodeButcherTable B);
SUNDIALS_EXPORT void ARKodeButcherTable_Space(ARKodeButcherTable B,
                                              sunindextype *liw,
                                              sunindextype *lrw);
SUNDIALS_EXPORT void ARKodeButcherTable_Free(ARKodeButcherTable B);
SUNDIALS_EXPORT void ARKodeButcherTable_Write(ARKodeButcherTable B, FILE *outfile);

SUNDIALS_EXPORT int ARKodeButcherTable_CheckOrder(ARKodeButcherTable B, int *q,
                                                  int *p, FILE *outfile);
SUNDIALS_EXPORT int ARKodeButcherTable_CheckARKOrder(ARKodeButcherTable B1,
                                                     ARKodeButcherTable B2,
                                                     int *q, int *p, FILE *outfile);

#ifdef __cplusplus
}
#endif

#endif
