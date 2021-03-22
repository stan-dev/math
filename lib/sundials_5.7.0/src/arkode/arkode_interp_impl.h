/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * Implementation header file for ARKode's temporal interpolation
 * utilities.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_INTERP_IMPL_H
#define _ARKODE_INTERP_IMPL_H

#include <stdarg.h>
#include <arkode/arkode.h>
#include "arkode_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARKode temporal interpolation constants
  ===============================================================*/

/* Numeric constants */
#define FOURTH RCONST(0.25)
#define THREE  RCONST(3.0)
#define SIX    RCONST(6.0)
#define TWELVE RCONST(12.0)



/*===============================================================
  ARKode Hermite Temporal Interpolation Data Structure
  ===============================================================*/

/* Hermite interpolation structure */

struct _ARKInterpContent_Hermite {
  int      degree;  /* maximum interpolant degree to use           */
  N_Vector fold;    /* f(t,y) at beginning of last successful step */
  N_Vector fnew;    /* f(t,y) at end of last successful step       */
  N_Vector yold;    /* y at beginning of last successful step      */
  N_Vector ynew;    /* y at end of last successful step            */
  N_Vector fa;      /* f(t,y) used in higher-order interpolation   */
  N_Vector fb;      /* f(t,y) used in higher-order interpolation   */
  realtype told;    /* t at beginning of last successful step      */
  realtype tnew;    /* t at end of last successful step            */
  realtype h;       /* last successful step size                   */
};

typedef struct _ARKInterpContent_Hermite *ARKInterpContent_Hermite;


/* Hermite structure accessor macros */

#define HINT_CONTENT(I)     ( (ARKInterpContent_Hermite)(I->content) )
#define HINT_DEGREE(I)      ( HINT_CONTENT(I)->degree )
#define HINT_FOLD(I)        ( HINT_CONTENT(I)->fold )
#define HINT_FNEW(I)        ( HINT_CONTENT(I)->fnew )
#define HINT_YOLD(I)        ( HINT_CONTENT(I)->yold )
#define HINT_YNEW(I)        ( HINT_CONTENT(I)->ynew )
#define HINT_FA(I)          ( HINT_CONTENT(I)->fa )
#define HINT_FB(I)          ( HINT_CONTENT(I)->fb )
#define HINT_TOLD(I)        ( HINT_CONTENT(I)->told )
#define HINT_TNEW(I)        ( HINT_CONTENT(I)->tnew )
#define HINT_H(I)           ( HINT_CONTENT(I)->h )


/* Hermite structure operations */

ARKInterp arkInterpCreate_Hermite(void* arkode_mem, int degree);

int arkInterpResize_Hermite(void* arkode_mem, ARKInterp interp,
                            ARKVecResizeFn resize, void *resize_data,
                            sunindextype lrw_diff, sunindextype liw_diff,
                            N_Vector tmpl);
void arkInterpFree_Hermite(void* arkode_mem, ARKInterp interp);
void arkInterpPrintMem_Hermite(ARKInterp interp, FILE *outfile);
int arkInterpSetDegree_Hermite(void *arkode_mem, ARKInterp interp, int degree);
int arkInterpInit_Hermite(void* arkode_mem, ARKInterp interp,
                          realtype tnew);
int arkInterpUpdate_Hermite(void* arkode_mem, ARKInterp interp, realtype tnew);
int arkInterpEvaluate_Hermite(void* arkode_mem, ARKInterp interp,
                              realtype tau, int d, int order, N_Vector yout);





/*===============================================================
  ARKode Lagrange Temporal Interpolation Data Structure
  ===============================================================*/

/* Lagrange interpolation structure */

struct _ARKInterpContent_Lagrange {
  int       nmax;      /* number of previous solutions to use      */
  int       nmaxalloc; /* vectors allocated for previous solutions */
  N_Vector *yhist;     /* previous solution vectors                */
  realtype *thist;     /* 't' values associated with yhist         */
  int       nhist;     /* number of 'active' vectors in yhist      */
  realtype  tround;    /* unit roundoff for 't' values             */
};

typedef struct _ARKInterpContent_Lagrange *ARKInterpContent_Lagrange;


/* Lagrange structure accessor macros */

#define LINT_CONTENT(I)     ( (ARKInterpContent_Lagrange)(I->content) )
#define LINT_NMAX(I)        ( LINT_CONTENT(I)->nmax )
#define LINT_NMAXALLOC(I)   ( LINT_CONTENT(I)->nmaxalloc )
#define LINT_YHIST(I)       ( LINT_CONTENT(I)->yhist )
#define LINT_THIST(I)       ( LINT_CONTENT(I)->thist )
#define LINT_YJ(I,j)        ( (LINT_YHIST(I))[j] )
#define LINT_TJ(I,j)        ( (LINT_THIST(I))[j] )
#define LINT_NHIST(I)       ( LINT_CONTENT(I)->nhist )
#define LINT_TROUND(I)      ( LINT_CONTENT(I)->tround )


/* Lagrange structure operations */

ARKInterp arkInterpCreate_Lagrange(void* arkode_mem, int degree);

int arkInterpResize_Lagrange(void* arkode_mem, ARKInterp interp,
                             ARKVecResizeFn resize, void *resize_data,
                             sunindextype lrw_diff, sunindextype liw_diff,
                             N_Vector tmpl);
void arkInterpFree_Lagrange(void* arkode_mem, ARKInterp interp);
void arkInterpPrintMem_Lagrange(ARKInterp interp, FILE *outfile);
int arkInterpSetDegree_Lagrange(void *arkode_mem, ARKInterp interp, int degree);
int arkInterpInit_Lagrange(void* arkode_mem, ARKInterp interp,
                           realtype tnew);
int arkInterpUpdate_Lagrange(void* arkode_mem, ARKInterp interp, realtype tnew);
int arkInterpEvaluate_Lagrange(void* arkode_mem, ARKInterp interp,
                               realtype tau, int d, int order, N_Vector yout);


/* Lagrange structure utility routines */
realtype LBasis(ARKInterp interp, int idx, realtype t);
realtype LBasisD(ARKInterp interp, int idx, realtype t);
realtype LBasisD2(ARKInterp interp, int idx, realtype t);
realtype LBasisD3(ARKInterp interp, int idx, realtype t);

#ifdef __cplusplus
}
#endif

#endif
