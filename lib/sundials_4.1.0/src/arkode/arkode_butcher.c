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
 * This is the implementation file for Butcher table structure
 * for the ARKode infrastructure.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

/* tolerance for checking order conditions */
#define TOL  (SUNRsqrt(UNIT_ROUNDOFF))

/* Private utility functions for checking method order */
static int __mv(realtype **A, realtype *x, int s, realtype *b);
static int __vv(realtype *x, realtype *y, int s, realtype *z);
static int __vp(realtype *x, int l, int s, realtype *z);
static int __dot(realtype *x, realtype *y, int s, realtype *d);
static booleantype __rowsum(realtype **A, realtype *c, int s);
static booleantype __order1(realtype *b, int s);
static booleantype __order2(realtype *b, realtype *c, int s);
static booleantype __order3a(realtype *b, realtype *c1, realtype *c2, int s);
static booleantype __order3b(realtype *b, realtype **A, realtype *c, int s);
static booleantype __order4a(realtype *b, realtype *c1, realtype *c2, realtype *c3, int s);
static booleantype __order4b(realtype *b, realtype *c1, realtype **A, realtype *c2, int s);
static booleantype __order4c(realtype *b, realtype **A, realtype *c1, realtype *c2, int s);
static booleantype __order4d(realtype *b, realtype **A1, realtype **A2, realtype *c, int s);
static booleantype __order5a(realtype *b, realtype *c1, realtype *c2, realtype *c3, realtype *c4, int s);
static booleantype __order5b(realtype *b, realtype *c1, realtype *c2, realtype **A, realtype *c3, int s);
static booleantype __order5c(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype *c2, int s);
static booleantype __order5d(realtype *b, realtype *c1, realtype **A, realtype *c2, realtype *c3, int s);
static booleantype __order5e(realtype *b, realtype **A, realtype *c1, realtype *c2, realtype *c3, int s);
static booleantype __order5f(realtype *b, realtype *c1, realtype **A1, realtype **A2, realtype *c2, int s);
static booleantype __order5g(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype *c2, int s);
static booleantype __order5h(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype *c2, int s);
static booleantype __order5i(realtype *b, realtype **A1, realtype **A2, realtype **A3, realtype *c, int s);
static booleantype __order6a(realtype *b, realtype *c1, realtype *c2, realtype *c3, realtype *c4, realtype *c5, int s);
static booleantype __order6b(realtype *b, realtype *c1, realtype *c2, realtype *c3, realtype **A, realtype *c4, int s);
static booleantype __order6c(realtype *b, realtype *c1, realtype **A1, realtype *c2, realtype **A2, realtype *c3, int s);
static booleantype __order6d(realtype *b, realtype *c1, realtype *c2, realtype **A, realtype *c3, realtype *c4, int s);
static booleantype __order6e(realtype *b, realtype *c1, realtype *c2, realtype **A1, realtype **A2, realtype *c3, int s);
static booleantype __order6f(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype **A3, realtype *c2, int s);
static booleantype __order6g(realtype *b, realtype *c1, realtype **A, realtype *c2, realtype *c3, realtype *c4, int s);
static booleantype __order6h(realtype *b, realtype *c1, realtype **A1, realtype *c2, realtype **A2, realtype *c3, int s);
static booleantype __order6i(realtype *b, realtype *c1, realtype **A1, realtype **A2, realtype *c2, realtype *c3, int s);
static booleantype __order6j(realtype *b, realtype *c1, realtype **A1, realtype **A2, realtype **A3, realtype *c2, int s);
static booleantype __order6k(realtype *b, realtype **A, realtype *c1, realtype *c2, realtype *c3, realtype *c4, int s);
static booleantype __order6l(realtype *b, realtype **A1, realtype *c1, realtype *c2, realtype **A2, realtype *c3, int s);
static booleantype __order6m(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype **A3, realtype *c2, int s);
static booleantype __order6n(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype *c2, realtype *c3, int s);
static booleantype __order6o(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype **A3, realtype *c2, int s);
static booleantype __order6p(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype *c2, realtype *c3, int s);
static booleantype __order6q(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype **A3, realtype *c2, int s);
static booleantype __order6r(realtype *b, realtype **A1, realtype **A2, realtype **A3, realtype *c1, realtype *c2, int s);
static booleantype __order6s(realtype *b, realtype **A1, realtype **A2, realtype **A3, realtype **A4, realtype *c, int s);
static int __ButcherSimplifyingAssumptions(realtype **A, realtype *b, realtype *c, int s);


/*---------------------------------------------------------------
  Routine to allocate an empty Butcher table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_Alloc(int stages, booleantype embedded)
{
  int i;
  ARKodeButcherTable B;

  /* Check for legal 'stages' value */
  if (stages < 1) return(NULL);

  /* Allocate Butcher table structure */
  B = NULL;
  B = (ARKodeButcherTable) malloc(sizeof(struct ARKodeButcherTableMem));
  if (B == NULL) return(NULL);

  /* initialize pointers in B structure to NULL */
  B->A = NULL;
  B->b = NULL;
  B->c = NULL;
  B->d = NULL;

  /* set stages into B structure */
  B->stages = stages;

  /*
   * Allocate fields within Butcher table structure
   */

  /* allocate rows of A */
  B->A = (realtype **) calloc( stages, sizeof(realtype*) );
  if (B->A == NULL) { ARKodeButcherTable_Free(B); return(NULL); }

  /* initialize each row of A to NULL */
  for (i=0; i<stages; i++)
    B->A[i] = NULL;

  /* allocate columns of A */
  for (i=0; i<stages; i++) {
    B->A[i] = (realtype *) calloc( stages, sizeof(realtype) );
    if (B->A[i] == NULL) { ARKodeButcherTable_Free(B); return(NULL); }
  }

  B->b = (realtype *) calloc( stages, sizeof(realtype) );
  if (B->b == NULL) { ARKodeButcherTable_Free(B); return(NULL); }

  B->c = (realtype *) calloc( stages, sizeof(realtype) );
  if (B->c == NULL) { ARKodeButcherTable_Free(B); return(NULL); }

  if (embedded) {
    B->d = (realtype *) calloc( stages, sizeof(realtype) );
    if (B->d == NULL) { ARKodeButcherTable_Free(B); return(NULL); }
  }

  /* initialize order parameters */
  B->q = 0;
  B->p = 0;

  return(B);
}


/*---------------------------------------------------------------
  Routine to allocate and fill a Butcher table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_Create(int s, int q, int p, realtype *c,
                                             realtype *A, realtype *b,
                                             realtype *d)
{
  int i, j;
  ARKodeButcherTable B;
  booleantype embedded;

  /* Check for legal number of stages */
  if (s < 1) return(NULL);

  /* Does the table have an embedding? */
  embedded = (d != NULL) ? SUNTRUE : SUNFALSE;

  /* Allocate Butcher table structure */
  B = ARKodeButcherTable_Alloc(s, embedded);
  if (B == NULL) return(NULL);

  /* set the relevant parameters */
  B->stages = s;
  B->q = q;
  B->p = p;

  for (i=0; i<s; i++) {
    B->c[i] = c[i];
    B->b[i] = b[i];
    for (j=0; j<s; j++) {
      B->A[i][j] = A[i*s + j];
    }
  }

  if (embedded)
    for (i=0; i<s; i++)
      B->d[i] = d[i];

  return(B);
}


/*---------------------------------------------------------------
  Routine to copy a Butcher table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_Copy(ARKodeButcherTable B)
{
  int i, j, s;
  ARKodeButcherTable Bcopy;
  booleantype embedded;

  /* Check for legal input */
  if (B == NULL) return(NULL);

  /* Get the number of stages */
  s = B->stages;

  /* Does the table have an embedding? */
  embedded = (B->d != NULL) ? SUNTRUE : SUNFALSE;

  /* Allocate Butcher table structure */
  Bcopy = ARKodeButcherTable_Alloc(s, embedded);
  if (Bcopy == NULL) return(NULL);

  /* set the relevant parameters */
  Bcopy->stages = B->stages;
  Bcopy->q = B->q;
  Bcopy->p = B->p;

  /* Copy Butcher table */
  for (i=0; i<s; i++) {
    Bcopy->c[i] = B->c[i];
    Bcopy->b[i] = B->b[i];
    for (j=0; j<s; j++) {
      Bcopy->A[i][j] = B->A[i][j];
    }
  }

  if (embedded)
    for (i=0; i<s; i++)
      Bcopy->d[i] = B->d[i];

  return(Bcopy);
}


/*---------------------------------------------------------------
  Routine to query the Butcher table structure workspace size
  ---------------------------------------------------------------*/
void ARKodeButcherTable_Space(ARKodeButcherTable B, sunindextype *liw,
                              sunindextype *lrw)
{
  /* initialize outputs and return if B is not allocated */
  *liw = 0;  *lrw = 0;
  if (B == NULL)  return;

  /* fill outputs based on B */
  *liw = 3;
  if (B->d != NULL) {
    *lrw = B->stages * (B->stages + 3);
  } else {
    *lrw = B->stages * (B->stages + 2);
  }
}


/*---------------------------------------------------------------
  Routine to free a Butcher table structure
  ---------------------------------------------------------------*/
void ARKodeButcherTable_Free(ARKodeButcherTable B)
{
  int i;

  /* Free each field within Butcher table structure, and then
     free structure itself */
  if (B != NULL) {
    if (B->d != NULL)  free(B->d);
    if (B->c != NULL)  free(B->c);
    if (B->b != NULL)  free(B->b);
    if (B->A != NULL) {
      for (i=0; i<B->stages; i++)
        if (B->A[i] != NULL)  free(B->A[i]);
      free(B->A);
    }

    free(B);
  }
}


/*---------------------------------------------------------------
  Routine to print a Butcher table structure
  ---------------------------------------------------------------*/
void ARKodeButcherTable_Write(ARKodeButcherTable B, FILE *outfile)
{
  int i, j;

  /* check for vaild table */
  if (B == NULL) return;
  if (B->A == NULL) return;
  for (i=0; i<B->stages; i++)
    if (B->A[i] == NULL) return;
  if (B->c == NULL) return;
  if (B->b == NULL) return;

  fprintf(outfile, "  A = \n");
  for (i=0; i<B->stages; i++) {
    fprintf(outfile, "      ");
    for (j=0; j<B->stages; j++)
      fprintf(outfile, "%"RSYM"  ", B->A[i][j]);
    fprintf(outfile, "\n");
  }

  fprintf(outfile, "  c = ");
  for (i=0; i<B->stages; i++)
    fprintf(outfile, "%"RSYM"  ", B->c[i]);
  fprintf(outfile, "\n");

  fprintf(outfile, "  b = ");
  for (i=0; i<B->stages; i++)
    fprintf(outfile, "%"RSYM"  ", B->b[i]);
  fprintf(outfile, "\n");

  if (B->d != NULL) {
    fprintf(outfile, "  d = ");
    for (i=0; i<B->stages; i++)
      fprintf(outfile, "%"RSYM"  ", B->d[i]);
    fprintf(outfile, "\n");
  }
}


/*---------------------------------------------------------------
  Routine to determine the analytical order of accuracy for a
  specified Butcher table.  We check the analytical [necessary]
  order conditions up through order 6.  After that, we revert to
  the [sufficient] Butcher simplifying assumptions.

  Inputs:
     B: Butcher table to check
     outfile: file pointer to print results; if NULL then no
        outputs are printed

  Outputs:
     q: measured order of accuracy for method
     p: measured order of accuracy for embedding [0 if not present]

  Return values:
     0 (success): internal {q,p} values match analytical order
     1 (warning): internal {q,p} values are lower than analytical
        order, or method achieves maximum order possible with this
        routine and internal {q,p} are higher.
    -1 (failure): internal p and q values are higher than analytical
         order
    -2 (failure): NULL-valued B (or critical contents)

  Note: for embedded methods, if the return flags for p and q would
  differ, failure takes precedence over warning, which takes
  precedence over success.
  ---------------------------------------------------------------*/
int ARKodeButcherTable_CheckOrder(ARKodeButcherTable B, int *q, int *p, FILE *outfile)
{
  /* local variables */
  int q_SA, p_SA, i, s;
  realtype **A, *b, *c, *d;
  booleantype alltrue;
  (*q) = (*p) = 0;

  /* verify non-NULL Butcher table structure and contents */
  if (B == NULL)  return(-2);
  if (B->stages < 1)  return(-2);
  if (B->A == NULL)  return(-2);
  for (i=0; i<B->stages; i++)
    if (B->A[i] == NULL)  return(-2);
  if (B->c == NULL)  return(-2);
  if (B->b == NULL)  return(-2);

  /* set shortcuts for Butcher table components */
  A = B->A;
  b = B->b;
  c = B->c;
  d = B->d;
  s = B->stages;

  /* check method order */
  if (outfile)  fprintf(outfile,"ARKodeButcherTable_CheckOrder:\n");

  /*    row sum condition */
  if (__rowsum(A, c, s)) {
    (*q) = 0;
  } else {
    (*q) = -1;
    if (outfile)  fprintf(outfile,"  method fails row sum condition\n");
  }
  /*    order 1 condition */
  if ((*q) == 0) {
    if (__order1(b, s)) {
      (*q) = 1;
    } else {
      if (outfile)  fprintf(outfile,"  method fails order 1 condition\n");
    }
  }
  /*    order 2 condition */
  if ((*q) == 1) {
    if (__order2(b, c, s)) {
      (*q) = 2;
    } else {
      if (outfile)  fprintf(outfile,"  method fails order 2 condition\n");
    }
  }
  /*    order 3 conditions */
  if ((*q) == 2) {
    alltrue = SUNTRUE;
    if (!__order3a(b, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 3 condition A\n");
    }
    if (!__order3b(b, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 3 condition B\n");
    }
    if (alltrue)  (*q) = 3;
  }
  /*    order 4 conditions */
  if ((*q) == 3) {
    alltrue = SUNTRUE;
    if (!__order4a(b, c, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 4 condition A\n");
    }
    if (!__order4b(b, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 4 condition B\n");
    }
    if (!__order4c(b, A, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 4 condition C\n");
    }
    if (!__order4d(b, A, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 4 condition D\n");
    }
    if (alltrue)  (*q) = 4;
  }
  /*    order 5 conditions */
  if ((*q) == 4) {
    alltrue = SUNTRUE;
    if (!__order5a(b, c, c, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 5 condition A\n");
    }
    if (!__order5b(b, c, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 5 condition B\n");
    }
    if (!__order5c(b, A, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 5 condition C\n");
    }
    if (!__order5d(b, c, A, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 5 condition D\n");
    }
    if (!__order5e(b, A, c, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 5 condition E\n");
    }
    if (!__order5f(b, c, A, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 5 condition F\n");
    }
    if (!__order5g(b, A, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 5 condition G\n");
    }
    if (!__order5h(b, A, A, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 5 condition H\n");
    }
    if (!__order5i(b, A, A, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 5 condition I\n");
    }
    if (alltrue)  (*q) = 5;
  }
  /*    order 6 conditions */
  if ((*q) == 5) {
    alltrue = SUNTRUE;
    if (!__order6a(b, c, c, c, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition A\n");
    }
    if (!__order6b(b, c, c, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition B\n");
    }
    if (!__order6c(b, c, A, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition C\n");
    }
    if (!__order6d(b, c, c, A, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition D\n");
    }
    if (!__order6e(b, c, c, A, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition E\n");
    }
    if (!__order6f(b, A, A, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition F\n");
    }
    if (!__order6g(b, c, A, c, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition G\n");
    }
    if (!__order6h(b, c, A, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition H\n");
    }
    if (!__order6i(b, c, A, A, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition I\n");
    }
    if (!__order6j(b, c, A, A, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition J\n");
    }
    if (!__order6k(b, A, c, c, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition K\n");
    }
    if (!__order6l(b, A, c, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition L\n");
    }
    if (!__order6m(b, A, A, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition M\n");
    }
    if (!__order6n(b, A, c, A, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition N\n");
    }
    if (!__order6o(b, A, c, A, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition O\n");
    }
    if (!__order6p(b, A, A, c, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition P\n");
    }
    if (!__order6q(b, A, A, c, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition Q\n");
    }
    if (!__order6r(b, A, A, A, c, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition R\n");
    }
    if (!__order6s(b, A, A, A, A, c, s)) {
      alltrue = SUNFALSE;
      if (outfile)  fprintf(outfile,"  method fails order 6 condition S\n");
    }
    if (alltrue)  (*q) = 6;
  }
  /*    higher order conditions (via simplifying assumptions) */
  if ((*q) == 6) {
    if (outfile)  fprintf(outfile,"  method order >= 6; reverting to simplifying assumptions\n");
    q_SA = __ButcherSimplifyingAssumptions(A, b, c, s);
    (*q) = SUNMAX((*q), q_SA);
    if (outfile)  fprintf(outfile,"  method order = %i\n", (*q));
  }

  /* check embedding order */
  if (d) {
    if (outfile)  fprintf(outfile,"\n");
    b = d;

    /*    row sum condition */
    if (__rowsum(A, c, s)) {
      (*p) = 0;
    } else {
      (*p) = -1;
      if (outfile)  fprintf(outfile,"  embedding fails row sum condition\n");
    }
    /*    order 1 condition */
    if ((*p) == 0) {
      if (__order1(b, s)) {
        (*p) = 1;
      } else {
        if (outfile)  fprintf(outfile,"  embedding fails order 1 condition\n");
      }
    }
    /*    order 2 condition */
    if ((*p) == 1) {
      if (__order2(b, c, s)) {
        (*p) = 2;
      } else {
        if (outfile)  fprintf(outfile,"  embedding fails order 2 condition\n");
      }
    }
    /*    order 3 conditions */
    if ((*p) == 2) {
      alltrue = SUNTRUE;
      if (!__order3a(b, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 3 condition A\n");
      }
      if (!__order3b(b, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 3 condition B\n");
      }
      if (alltrue)  (*p) = 3;
    }
    /*    order 4 conditions */
    if ((*p) == 3) {
      alltrue = SUNTRUE;
      if (!__order4a(b, c, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 4 condition A\n");
      }
      if (!__order4b(b, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 4 condition B\n");
      }
      if (!__order4c(b, A, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 4 condition C\n");
      }
      if (!__order4d(b, A, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 4 condition D\n");
      }
      if (alltrue)  (*p) = 4;
    }
    /*    order 5 conditions */
    if ((*p) == 4) {
      alltrue = SUNTRUE;
      if (!__order5a(b, c, c, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 5 condition A\n");
      }
      if (!__order5b(b, c, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 5 condition B\n");
      }
      if (!__order5c(b, A, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 5 condition C\n");
      }
      if (!__order5d(b, c, A, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 5 condition D\n");
      }
      if (!__order5e(b, A, c, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 5 condition E\n");
      }
      if (!__order5f(b, c, A, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 5 condition F\n");
      }
      if (!__order5g(b, A, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 5 condition G\n");
      }
      if (!__order5h(b, A, A, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 5 condition H\n");
      }
      if (!__order5i(b, A, A, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 5 condition I\n");
      }
      if (alltrue)  (*p) = 5;
    }
    /*    order 6 conditions */
    if ((*p) == 5) {
      alltrue = SUNTRUE;
      if (!__order6a(b, c, c, c, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition A\n");
      }
      if (!__order6b(b, c, c, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition B\n");
      }
      if (!__order6c(b, c, A, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition C\n");
      }
      if (!__order6d(b, c, c, A, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition D\n");
      }
      if (!__order6e(b, c, c, A, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition E\n");
      }
      if (!__order6f(b, A, A, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition F\n");
      }
      if (!__order6g(b, c, A, c, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition G\n");
      }
      if (!__order6h(b, c, A, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition H\n");
      }
      if (!__order6i(b, c, A, A, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition I\n");
      }
      if (!__order6j(b, c, A, A, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition J\n");
      }
      if (!__order6k(b, A, c, c, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition K\n");
      }
      if (!__order6l(b, A, c, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition L\n");
      }
      if (!__order6m(b, A, A, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition M\n");
      }
      if (!__order6n(b, A, c, A, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition N\n");
      }
      if (!__order6o(b, A, c, A, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition O\n");
      }
      if (!__order6p(b, A, A, c, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition P\n");
      }
      if (!__order6q(b, A, A, c, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition Q\n");
      }
      if (!__order6r(b, A, A, A, c, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition R\n");
      }
      if (!__order6s(b, A, A, A, A, c, s)) {
        alltrue = SUNFALSE;
        if (outfile)  fprintf(outfile,"  embedding fails order 6 condition S\n");
      }
      if (alltrue)  (*p) = 6;
    }
    /*    higher order conditions (via simplifying assumptions) */
    if ((*p) == 6) {
      if (outfile)  fprintf(outfile,"  embedding order >= 6; reverting to simplifying assumptions\n");
      p_SA = __ButcherSimplifyingAssumptions(A, b, c, s);
      (*p) = SUNMAX((*p), p_SA);
      if (outfile)  fprintf(outfile,"  embedding order = %i\n", (*p));
    }
  }

  /* compare results against stored values and return */

  /*    check failure modes first */
  if (((*q) < B->q) && ((*q) < 6))  return(-1);
  if (d)
    if (((*p) < B->p) && ((*p) < 6))  return(-1);

  /*    check warning modes */
  if ((*q) > B->q)  return(1);
  if (d)
    if ((*p) > B->p)  return(1);
  if (((*q) < B->q) && ((*q) >= 6)) return(1);
  if (d)
    if (((*p) < B->p) && ((*p) >= 6)) return(1);

  /*    return success */
  return(0);
}


/*---------------------------------------------------------------
  Routine to determine the analytical order of accuracy for a
  specified pair of Butcher tables in an ARK pair.  We check the
  analytical order conditions up through order 6.

  Inputs:
     B1, B2: Butcher tables to check
     outfile: file pointer to print results; if NULL then no
        outputs are printed

  Outputs:
     q: measured order of accuracy for method
     p: measured order of accuracy for embedding [0 if not present]

  Return values:
     0 (success): completed checks
     1 (warning): internal {q,p} values are lower than analytical
        order, or method achieves maximum order possible with this
        routine and internal {q,p} are higher.
    -1 (failure): NULL-valued B1, B2 (or critical contents)

  Note: for embedded methods, if the return flags for p and q would
  differ, warning takes precedence over success.
  ---------------------------------------------------------------*/
int ARKodeButcherTable_CheckARKOrder(ARKodeButcherTable B1,
                                     ARKodeButcherTable B2,
                                     int *q, int *p, FILE *outfile)
{
  /* local variables */
  int i, j, k, l, m, n, s;
  booleantype alltrue;
  realtype **A[2], *b[2], *c[2], *d[2];
  (*q) = (*p) = 0;

  /* verify non-NULL Butcher table structure and contents */
  if (B1 == NULL)  return(-1);
  if (B1->stages < 1)  return(-1);
  if (B1->A == NULL)  return(-1);
  for (i=0; i<B1->stages; i++)
    if (B1->A[i] == NULL)  return(-1);
  if (B1->c == NULL)  return(-1);
  if (B1->b == NULL)  return(-1);
  if (B2 == NULL)  return(-1);
  if (B2->stages < 1)  return(-1);
  if (B2->A == NULL)  return(-1);
  for (i=0; i<B2->stages; i++)
    if (B2->A[i] == NULL)  return(-1);
  if (B2->c == NULL)  return(-1);
  if (B2->b == NULL)  return(-1);
  if (B1->stages != B2->stages)  return(-1);

  /* set shortcuts for Butcher table components */
  A[0] = B1->A;
  b[0] = B1->b;
  c[0] = B1->c;
  d[0] = B1->d;
  A[1] = B2->A;
  b[1] = B2->b;
  c[1] = B2->c;
  d[1] = B1->d;
  s = B1->stages;

  /* check method order */
  if (outfile)  fprintf(outfile,"ARKodeButcherTable_CheckARKOrder:\n");

  /*    row sum conditions */
  if (__rowsum(A[0], c[0], s) && __rowsum(A[1], c[1], s)) {
    (*q) = 0;
  } else {
    (*q) = -1;
    if (outfile)  fprintf(outfile,"  method fails row sum conditions\n");
  }
  /*    order 1 conditions */
  if ((*q) == 0) {
    if (__order1(b[0], s) && __order1(b[1], s)) {
      (*q) = 1;
    } else {
      if (outfile)  fprintf(outfile,"  method fails order 1 conditions\n");
    }
  }
  /*    order 2 conditions */
  if ((*q) == 1) {
    alltrue = SUNTRUE;
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        alltrue = (alltrue && __order2(b[i], c[j], s));
    if (alltrue) {
      (*q) = 2;
    } else {
      if (outfile)  fprintf(outfile,"  method fails order 2 conditions\n");
    }
  }
  /*    order 3 conditions */
  if ((*q) == 2) {
    alltrue = SUNTRUE;
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          alltrue = (alltrue && __order3a(b[i], c[j], c[k], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 3 conditions A\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          alltrue = (alltrue && __order3b(b[i], A[j], c[k], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 3 conditions B\n");
    if (alltrue)  (*q) = 3;
  }
  /*    order 4 conditions */
  if ((*q) == 3) {
    alltrue = SUNTRUE;
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            alltrue = (alltrue && __order4a(b[i], c[j], c[k], c[l], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 4 conditions A\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            alltrue = (alltrue && __order4b(b[i], c[j], A[k], c[l], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 4 conditions B\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            alltrue = (alltrue && __order4c(b[i], A[j], c[k], c[l], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 4 conditions C\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            alltrue = (alltrue && __order4d(b[i], A[j], A[k], c[l], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 4 conditions D\n");
    if (alltrue)  (*q) = 4;
  }
  /*    order 5 conditions */
  if ((*q) == 4) {
    alltrue = SUNTRUE;
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              alltrue = (alltrue && __order5a(b[i], c[j], c[k], c[l], c[m], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 5 conditions A\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              alltrue = (alltrue && __order5b(b[i], c[j], c[k], A[l], c[m], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 5 conditions B\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              alltrue = (alltrue && __order5c(b[i], A[j], c[k], A[l], c[m], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 5 conditions C\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              alltrue = (alltrue && __order5d(b[i], c[j], A[k], c[l], c[m], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 5 conditions D\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              alltrue = (alltrue && __order5e(b[i], A[j], c[k], c[l], c[m], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 5 conditions E\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              alltrue = (alltrue && __order5f(b[i], c[j], A[k], A[l], c[m], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 5 conditions F\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              alltrue = (alltrue && __order5g(b[i], A[j], c[k], A[l], c[m], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 5 conditions G\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              alltrue = (alltrue && __order5h(b[i], A[j], A[k], c[l], c[m], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 5 conditions H\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              alltrue = (alltrue && __order5i(b[i], A[j], A[k], A[l], c[m], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 5 conditions I\n");
    if (alltrue)  (*q) = 5;
  }
  /*    order 6 conditions */
  if ((*q) == 5) {
    alltrue = SUNTRUE;
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6a(b[i], c[j], c[k], c[l], c[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions A\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6b(b[i], c[j], c[k], c[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions B\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6c(b[i], c[j], A[k], c[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions C\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6d(b[i], c[j], c[k], A[l], c[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions D\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6e(b[i], c[j], c[k], A[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions E\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6f(b[i], A[j], A[k], c[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions F\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6g(b[i], c[j], A[k], c[l], c[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions G\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6h(b[i], c[j], A[k], c[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions H\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6i(b[i], c[j], A[k], A[l], c[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions I\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6j(b[i], c[j], A[k], A[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions J\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6k(b[i], A[j], c[k], c[l], c[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions K\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6l(b[i], A[j], c[k], c[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions L\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6m(b[i], A[j], A[k], c[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions M\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6n(b[i], A[j], c[k], A[l], c[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions N\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6o(b[i], A[j], c[k], A[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions O\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6p(b[i], A[j], A[k], c[l], c[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions P\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6q(b[i], A[j], A[k], c[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions Q\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6r(b[i], A[j], A[k], A[l], c[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions R\n");
    for (i=0; i<2; i++)
      for (j=0; j<2; j++)
        for (k=0; k<2; k++)
          for (l=0; l<2; l++)
            for (m=0; m<2; m++)
              for (n=0; n<2; n++)
                alltrue = (alltrue && __order6s(b[i], A[j], A[k], A[l], A[m], c[n], s));
    if ( (!alltrue) && outfile)  fprintf(outfile,"  method fails order 6 conditions S\n");
    if (alltrue)  (*q) = 6;
  }

  /* check embedding order */
  if (d[0] && d[1]) {
    if (outfile)  fprintf(outfile,"\n");

    /*    row sum conditions */
    if (__rowsum(A[0], c[0], s) && __rowsum(A[1], c[1], s)) {
      (*p) = 0;
    } else {
      (*p) = -1;
      if (outfile)  fprintf(outfile,"  embedding fails row sum conditions\n");
    }
    /*    order 1 conditions */
    if ((*p) == 0) {
      if (__order1(d[0], s) && __order1(d[1], s)) {
        (*p) = 1;
      } else {
        if (outfile)  fprintf(outfile,"  embedding fails order 1 conditions\n");
      }
    }
    /*    order 2 conditions */
    if ((*p) == 1) {
      alltrue = SUNTRUE;
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          alltrue = (alltrue && __order2(d[i], c[j], s));
      if (alltrue) {
        (*p) = 2;
      } else {
        if (outfile)  fprintf(outfile,"  embedding fails order 2 conditions\n");
      }
    }
    /*    order 3 conditions */
    if ((*p) == 2) {
      alltrue = SUNTRUE;
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            alltrue = (alltrue && __order3a(d[i], c[j], c[k], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 3 conditions A\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            alltrue = (alltrue && __order3b(d[i], A[j], c[k], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 3 conditions B\n");
      if (alltrue)  (*p) = 3;
    }
    /*    order 4 conditions */
    if ((*p) == 3) {
      alltrue = SUNTRUE;
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              alltrue = (alltrue && __order4a(d[i], c[j], c[k], c[l], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 4 conditions A\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              alltrue = (alltrue && __order4b(d[i], c[j], A[k], c[l], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 4 conditions B\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              alltrue = (alltrue && __order4c(d[i], A[j], c[k], c[l], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 4 conditions C\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              alltrue = (alltrue && __order4d(d[i], A[j], A[k], c[l], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 4 conditions D\n");
      if (alltrue)  (*p) = 4;
    }
    /*    order 5 conditions */
    if ((*p) == 4) {
      alltrue = SUNTRUE;
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                alltrue = (alltrue && __order5a(d[i], c[j], c[k], c[l], c[m], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 5 conditions A\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                alltrue = (alltrue && __order5b(d[i], c[j], c[k], A[l], c[m], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 5 conditions B\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                alltrue = (alltrue && __order5c(d[i], A[j], c[k], A[l], c[m], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 5 conditions C\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                alltrue = (alltrue && __order5d(d[i], c[j], A[k], c[l], c[m], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 5 conditions D\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                alltrue = (alltrue && __order5e(d[i], A[j], c[k], c[l], c[m], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 5 conditions E\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                alltrue = (alltrue && __order5f(d[i], c[j], A[k], A[l], c[m], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 5 conditions F\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                alltrue = (alltrue && __order5g(d[i], A[j], c[k], A[l], c[m], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 5 conditions G\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                alltrue = (alltrue && __order5h(d[i], A[j], A[k], c[l], c[m], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 5 conditions H\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                alltrue = (alltrue && __order5i(d[i], A[j], A[k], A[l], c[m], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 5 conditions I\n");
      if (alltrue)  (*p) = 5;
    }
    /*    order 6 conditions */
    if ((*p) == 5) {
      alltrue = SUNTRUE;
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6a(d[i], c[j], c[k], c[l], c[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions A\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6b(d[i], c[j], c[k], c[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions B\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6c(d[i], c[j], A[k], c[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions C\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6d(d[i], c[j], c[k], A[l], c[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions D\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6e(d[i], c[j], c[k], A[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions E\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6f(d[i], A[j], A[k], c[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions F\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6g(d[i], c[j], A[k], c[l], c[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions G\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6h(d[i], c[j], A[k], c[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions H\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6i(d[i], c[j], A[k], A[l], c[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions I\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6j(d[i], c[j], A[k], A[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions J\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6k(d[i], A[j], c[k], c[l], c[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions K\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6l(d[i], A[j], c[k], c[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions L\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6m(d[i], A[j], A[k], c[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions M\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6n(d[i], A[j], c[k], A[l], c[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions N\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6o(d[i], A[j], c[k], A[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions O\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6p(d[i], A[j], A[k], c[l], c[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions P\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6q(d[i], A[j], A[k], c[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions Q\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6r(d[i], A[j], A[k], A[l], c[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions R\n");
      for (i=0; i<2; i++)
        for (j=0; j<2; j++)
          for (k=0; k<2; k++)
            for (l=0; l<2; l++)
              for (m=0; m<2; m++)
                for (n=0; n<2; n++)
                  alltrue = (alltrue && __order6s(d[i], A[j], A[k], A[l], A[m], c[n], s));
      if ( (!alltrue) && outfile)  fprintf(outfile,"  embedding fails order 6 conditions S\n");
      if (alltrue)  (*p) = 6;
    }
  }

  /* compare results against stored values and return */

  /*    check warning modes */
  if ((*q) > B1->q)  return(1);
  if ((*q) > B2->q)  return(1);
  if (d[0] && d[1]) {
    if ((*p) > B1->p)  return(1);
    if ((*p) > B2->p)  return(1);
  }
  if (((*q) < B1->q) && ((*q) == 6)) return(1);
  if (((*q) < B2->q) && ((*q) == 6)) return(1);
  if (d[0] && d[1]) {
    if (((*p) < B1->p) && ((*p) == 6)) return(1);
    if (((*p) < B2->p) && ((*p) == 6)) return(1);
  }

  /*    return success */
  return(0);
}


/*---------------------------------------------------------------
  Private utility routines for checking method order
  ---------------------------------------------------------------*/

/*---------------------------------------------------------------
  Utility routine to compute small dense matrix-vector product
       b = A*x
  Here A is (s x s), x and b are (s x 1).  Returns 0 on success,
  nonzero on failure.
  ---------------------------------------------------------------*/
static int __mv(realtype **A, realtype *x, int s, realtype *b)
{
  int i, j;
  if ((A == NULL) || (x == NULL) || (b == NULL) || (s < 1))
    return(1);
  for (i=0; i<s; i++)  b[i] = RCONST(0.0);
  for (i=0; i<s; i++)
    for (j=0; j<s; j++)
      b[i] += A[i][j]*x[j];
  return(0);
}


/*---------------------------------------------------------------
  Utility routine to compute small vector .* vector product
       z = x.*y   [Matlab notation]
  Here all vectors are (s x 1).   Returns 0 on success,
  nonzero on failure.
  ---------------------------------------------------------------*/
static int __vv(realtype *x, realtype *y, int s, realtype *z)
{
  int i;
  if ((x == NULL) || (y == NULL) || (z == NULL) || (s < 1))
    return(1);
  for (i=0; i<s; i++)
    z[i] = x[i]*y[i];
  return(0);
}


/*---------------------------------------------------------------
  Utility routine to compute small vector .^ int
       z = x.^l   [Matlab notation]
  Here all vectors are (s x 1).   Returns 0 on success,
  nonzero on failure.
  ---------------------------------------------------------------*/
static int __vp(realtype *x, int l, int s, realtype *z)
{
  int i;
  if ((x == NULL) || (z == NULL) || (s < 1) || (s < 0))
    return(1);
  for (i=0; i<s; i++)
    z[i] = SUNRpowerI(x[i],l);
  return(0);
}


/*---------------------------------------------------------------
  Utility routine to compute small vector dot product:
       d = dot(x,y)
  Here x and y are (s x 1), and d is scalar.   Returns 0 on success,
  nonzero on failure.
  ---------------------------------------------------------------*/
static int __dot(realtype *x, realtype *y, int s, realtype *d)
{
  int i;
  if ((x == NULL) || (y == NULL) || (d == NULL) || (s < 1))
    return(1);
  (*d) = RCONST(0.0);
  for (i=0; i<s; i++)
    (*d) += x[i]*y[i];
  return(0);
}


/*---------------------------------------------------------------
  Utility routines to check specific order conditions.  Each
  returns SUNTRUE on success, SUNFALSE on failure.
     Order 0:  __rowsum
     Order 1:  __order1
     Order 2:  __order2
     Order 3:  __order3a and __order3b
     Order 4:  __order4a through __order4d
     Order 5:  __order5a through __order5i
     Order 6:  __order6a through __order6s
  ---------------------------------------------------------------*/

/* c(i) = sum(A(i,:)) */
static booleantype __rowsum(realtype **A, realtype *c, int s)
{
  int i, j;
  realtype rsum;
  for (i=0; i<s; i++) {
    rsum = RCONST(0.0);
    for (j=0; j<s; j++)
      rsum += A[i][j];
    if (SUNRabs(rsum - c[i]) > TOL)
      return(SUNFALSE);
  }
  return(SUNTRUE);
}

/* b'*e = 1 */
static booleantype __order1(realtype *b, int s)
{
  int i;
  realtype err = RCONST(1.0);
  for (i=0; i<s; i++)
    err -= b[i];
  return (SUNRabs(err) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*c = 1/2 */
static booleantype __order2(realtype *b, realtype *c, int s)
{
  realtype bc;
  if (__dot(b,c,s,&bc))  return(SUNFALSE);
  return (SUNRabs(bc - RCONST(0.5)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(c1.*c2) = 1/3 */
static booleantype __order3a(realtype *b, realtype *c1, realtype *c2, int s)
{
  realtype bcc;
  realtype *tmp = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp)) { free(tmp); return(SUNFALSE); }
  if (__dot(b,tmp,s,&bcc))  return(SUNFALSE);
  free(tmp);
  return (SUNRabs(bcc - RCONST(1.0)/RCONST(3.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(A*c) = 1/6 */
static booleantype __order3b(realtype *b, realtype **A, realtype *c, int s)
{
  realtype bAc;
  realtype *tmp = calloc( s, sizeof(realtype) );
  if (__mv(A,c,s,tmp)) { free(tmp); return(SUNFALSE); }
  if (__dot(b,tmp,s,&bAc))  return(SUNFALSE);
  free(tmp);
  return (SUNRabs(bAc - RCONST(1.0)/RCONST(6.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(c1.*c2.*c3) = 1/4 */
static booleantype __order4a(realtype *b, realtype *c1, realtype *c2, realtype *c3, int s)
{
  realtype bccc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c3,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bccc))  return(SUNFALSE);
  free(tmp1);  free(tmp2);
  return (SUNRabs(bccc - RCONST(0.25)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* (b.*c1)'*(A*c2) = 1/8 */
static booleantype __order4b(realtype *b, realtype *c1, realtype **A, realtype *c2, int s)
{
  realtype bcAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(b,c1,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A,c2,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(tmp1,tmp2,s,&bcAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bcAc - RCONST(0.125)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A*(c1.*c2) = 1/12 */
static booleantype __order4c(realtype *b, realtype **A, realtype *c1, realtype *c2, int s)
{
  realtype bAcc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAcc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAcc - RCONST(1.0)/RCONST(12.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*A2*c = 1/24 */
static booleantype __order4d(realtype *b, realtype **A1, realtype **A2, realtype *c, int s)
{
  realtype bAAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A2,c,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAAc - RCONST(1.0)/RCONST(24.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(c1.*c2.*c3.*c4) = 1/5 */
static booleantype __order5a(realtype *b, realtype *c1, realtype *c2, realtype *c3, realtype *c4, int s)
{
  realtype bcccc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c3,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c4,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp1,s,&bcccc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bcccc - RCONST(0.2)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* (b.*c1.*c2)'*(A*c3) = 1/10 */
static booleantype __order5b(realtype *b, realtype *c1, realtype *c2, realtype **A, realtype *c3, int s)
{
  realtype bccAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(b,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A,c3,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(tmp1,tmp2,s,&bccAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bccAc - RCONST(0.1)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*((A1*c1).*(A2*c2)) = 1/20 */
static booleantype __order5c(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype *c2, int s)
{
  realtype bAcAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  realtype *tmp3 = calloc( s, sizeof(realtype) );
  if (__mv(A1,c1,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__mv(A2,c2,s,tmp2)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__vv(tmp1,tmp2,s,tmp3)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__dot(b,tmp3,s,&bAcAc))  return(SUNFALSE);
  free(tmp1); free(tmp2); free(tmp3);
  return (SUNRabs(bAcAc - RCONST(0.05)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* (b.*c1)'*A*(c2.*c3) = 1/15 */
static booleantype __order5d(realtype *b, realtype *c1, realtype **A, realtype *c2, realtype *c3, int s)
{
  realtype bcAcc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c2,c3,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(b,c1,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(tmp1,tmp2,s,&bcAcc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bcAcc - RCONST(1.0)/RCONST(15.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A*(c1.*c2.*c3) = 1/20 */
static booleantype __order5e(realtype *b, realtype **A, realtype *c1, realtype *c2, realtype *c3, int s)
{
  realtype bAccc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c3,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp1,s,&bAccc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAccc - RCONST(0.05)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* (b.*c1)'*A1*A2*c2 = 1/30 */
static booleantype __order5f(realtype *b, realtype *c1, realtype **A1, realtype **A2, realtype *c2, int s)
{
  realtype bcAAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A2,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(b,c1,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(tmp1,tmp2,s,&bcAAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bcAAc - RCONST(1.0)/RCONST(30.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*(c1.*(A2*c2)) = 1/40 */
static booleantype __order5g(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype *c2, int s)
{
  realtype bAcAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A2,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp1,s,&bAcAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAcAc - RCONST(1.0)/RCONST(40.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*A2*(c1.*c2) = 1/60 */
static booleantype __order5h(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype *c2, int s)
{
  realtype bAAcc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp1,s,&bAAcc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAAcc - RCONST(1.0)/RCONST(60.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*A2*A3*c = 1/120 */
static booleantype __order5i(realtype *b, realtype **A1, realtype **A2, realtype **A3, realtype *c, int s)
{
  realtype bAAAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A3,c,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp1,s,&bAAAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAAAc - RCONST(1.0)/RCONST(120.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(c1.*c2.*c3.*c4.*c5) = 1/6 */
static booleantype __order6a(realtype *b, realtype *c1, realtype *c2, realtype *c3, realtype *c4, realtype *c5, int s)
{
  realtype bccccc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c3,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c4,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c5,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bccccc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bccccc - RCONST(1.0)/RCONST(6.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* (b.*c1.*c2.*c3)'*(A*c4) = 1/12 */
static booleantype __order6b(realtype *b, realtype *c1, realtype *c2, realtype *c3, realtype **A, realtype *c4, int s)
{
  realtype bcccAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(b,c1,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c3,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A,c4,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(tmp1,tmp2,s,&bcccAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bcccAc - RCONST(1.0)/RCONST(12.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(c1.*(A1*c2).*(A2*c3)) = 1/24 */
static booleantype __order6c(realtype *b, realtype *c1, realtype **A1, realtype *c2, realtype **A2, realtype *c3, int s)
{
  realtype bcAc2;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  realtype *tmp3 = calloc( s, sizeof(realtype) );
  if (__mv(A2,c3,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__mv(A1,c2,s,tmp2)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__vv(tmp1,tmp2,s,tmp3)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__vv(c1,tmp3,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__dot(b,tmp1,s,&bcAc2))  return(SUNFALSE);
  free(tmp1); free(tmp2); free(tmp3);
  return (SUNRabs(bcAc2 - RCONST(1.0)/RCONST(24.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* (b.*c1.*c2)'*A*(c3.*c4) = 1/18 */
static booleantype __order6d(realtype *b, realtype *c1, realtype *c2, realtype **A, realtype *c3, realtype *c4, int s)
{
  realtype bccAcc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  realtype *tmp3 = calloc( s, sizeof(realtype) );
  if (__vv(c3,c4,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__mv(A,tmp1,s,tmp2)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__vv(b,tmp1,s,tmp3)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__dot(tmp2,tmp3,s,&bccAcc))  return(SUNFALSE);
  free(tmp1); free(tmp2);  free(tmp3);
  return (SUNRabs(bccAcc - RCONST(1.0)/RCONST(18.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* (b.*(c1.*c2))'*A1*A2*c3 = 1/36 */
static booleantype __order6e(realtype *b, realtype *c1, realtype *c2, realtype **A1, realtype **A2, realtype *c3, int s)
{
  realtype bccAAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  realtype *tmp3 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__vv(b,tmp1,s,tmp2)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__mv(A2,c3,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp3)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__dot(tmp2,tmp3,s,&bccAAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);  free(tmp3);
  return (SUNRabs(bccAAc - RCONST(1.0)/RCONST(36.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*((A1*A2*c1).*(A3*c2)) = 1/72 */
static booleantype __order6f(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype **A3, realtype *c2, int s)
{
  realtype bAAcAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  realtype *tmp3 = calloc( s, sizeof(realtype) );
  if (__mv(A2,c1,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__mv(A3,c2,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__vv(tmp1,tmp2,s,tmp3)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__dot(b,tmp3,s,&bAAcAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);  free(tmp3);
  return (SUNRabs(bAAcAc - RCONST(1.0)/RCONST(72.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(c1.*(A*(c2.*c3.*c4))) = 1/24 */
static booleantype __order6g(realtype *b, realtype *c1, realtype **A, realtype *c2, realtype *c3, realtype *c4, int s)
{
  realtype bcAccc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c2,c3,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c4,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bcAccc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bcAccc - RCONST(1.0)/RCONST(24.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(c1.*(A1*(c2.*(A2*c3)))) = 1/48 */
static booleantype __order6h(realtype *b, realtype *c1, realtype **A1, realtype *c2, realtype **A2, realtype *c3, int s)
{
  realtype bcAcAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A2,c3,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bcAcAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bcAcAc - RCONST(1.0)/RCONST(48.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(c1.*(A1*A2*(c2.*c3))) = 1/72 */
static booleantype __order6i(realtype *b, realtype *c1, realtype **A1, realtype **A2, realtype *c2, realtype *c3, int s)
{
  realtype bcAAcc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c2,c3,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bcAAcc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bcAAcc - RCONST(1.0)/RCONST(72.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*(c1.*(A1*A2*A3*c2)) = 1/144 */
static booleantype __order6j(realtype *b, realtype *c1, realtype **A1, realtype **A2, realtype **A3, realtype *c2, int s)
{
  realtype bcAAAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A3,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bcAAAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bcAAAc - RCONST(1.0)/RCONST(144.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A*(c1.*c2.*c3.*c4) = 1/30 */
static booleantype __order6k(realtype *b, realtype **A, realtype *c1, realtype *c2, realtype *c3, realtype *c4, int s)
{
  realtype bAcccc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c3,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c4,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAcccc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAcccc - RCONST(1.0)/RCONST(30.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*(c1.*c2.*(A2*c3)) = 1/60 */
static booleantype __order6l(realtype *b, realtype **A1, realtype *c1, realtype *c2, realtype **A2, realtype *c3, int s)
{
  realtype bAccAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A2,c3,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c1,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAccAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAccAc - RCONST(1.0)/RCONST(60.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*((A2*c1).*(A3*c2)) = 1/120 */
static booleantype __order6m(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype **A3, realtype *c2, int s)
{
  realtype bAAcAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  realtype *tmp3 = calloc( s, sizeof(realtype) );
  if (__mv(A3,c2,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__mv(A2,c1,s,tmp2)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__vv(tmp1,tmp2,s,tmp3)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp3,s,tmp1)) { free(tmp1); free(tmp2); free(tmp3); return(SUNFALSE); }
  if (__dot(b,tmp1,s,&bAAcAc))  return(SUNFALSE);
  free(tmp1); free(tmp2); free(tmp3);
  return (SUNRabs(bAAcAc - RCONST(1.0)/RCONST(120.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*(c1.*(A2*(c2.*c3))) = 1/90 */
static booleantype __order6n(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype *c2, realtype *c3, int s)
{
  realtype bAcAcc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c2,c3,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c1,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAcAcc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAcAcc - RCONST(1.0)/RCONST(90.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*(c1.*(A2*A3*c2)) = 1/180 */
static booleantype __order6o(realtype *b, realtype **A1, realtype *c1, realtype **A2, realtype **A3, realtype *c2, int s)
{
  realtype bAcAAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A3,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c1,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAcAAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAcAAc - RCONST(1.0)/RCONST(180.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*A2*(c1.*c2.*c3) = 1/120 */
static booleantype __order6p(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype *c2, realtype *c3, int s)
{
  realtype bAAccc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c3,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAAccc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAAccc - RCONST(1.0)/RCONST(120.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*A2*(c1.*(A3*c2)) = 1/240 */
static booleantype __order6q(realtype *b, realtype **A1, realtype **A2, realtype *c1, realtype **A3, realtype *c2, int s)
{
  realtype bAAcAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A3,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__vv(c1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAAcAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAAcAc - RCONST(1.0)/RCONST(240.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*A2*A3*(c1.*c2) = 1/360 */
static booleantype __order6r(realtype *b, realtype **A1, realtype **A2, realtype **A3, realtype *c1, realtype *c2, int s)
{
  realtype bAAAcc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__vv(c1,c2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A3,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAAAcc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAAAcc - RCONST(1.0)/RCONST(360.0)) > TOL) ? SUNFALSE : SUNTRUE;
}

/* b'*A1*A2*A3*A4*c = 1/720 */
static booleantype __order6s(realtype *b, realtype **A1, realtype **A2, realtype **A3, realtype **A4, realtype *c, int s)
{
  realtype bAAAAc;
  realtype *tmp1 = calloc( s, sizeof(realtype) );
  realtype *tmp2 = calloc( s, sizeof(realtype) );
  if (__mv(A4,c,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A2,tmp2,s,tmp1)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__mv(A1,tmp1,s,tmp2)) { free(tmp1); free(tmp2); return(SUNFALSE); }
  if (__dot(b,tmp2,s,&bAAAAc))  return(SUNFALSE);
  free(tmp1); free(tmp2);
  return (SUNRabs(bAAAAc - RCONST(1.0)/RCONST(720.0)) > TOL) ? SUNFALSE : SUNTRUE;
}


/*---------------------------------------------------------------
  Utility routine to check Butcher's simplifying assumptions.
  Returns the maximum predicted order.
  ---------------------------------------------------------------*/
static int __ButcherSimplifyingAssumptions(realtype **A, realtype *b, realtype *c, int s)
{
  int P, Q, R, i, j, k, q;
  realtype RHS, LHS;
  booleantype alltrue;
  realtype *tmp = calloc( s, sizeof(realtype) );

  /* B(P) */
  P = 0;
  for (i=1; i<1000; i++) {
    if (__vp(c,i-1,s,tmp)) { free(tmp); return(0); }
    if (__dot(b,tmp,s,&LHS)) { free(tmp); return(0); }
    RHS = RCONST(1.0)/i;
    if (SUNRabs(RHS-LHS) > TOL)
      break;
    P++;
  }

  /* C(Q) */
  Q = 0;
  for (k=1; k<1000; k++) {
    alltrue = SUNTRUE;
    for (i=0; i<s; i++) {
      if (__vp(c,k-1,s,tmp)) { free(tmp); return(0); }
      if (__dot(A[i],tmp,s,&LHS)) { free(tmp); return(0); }
      RHS = SUNRpowerI(c[i],k) / k;
      if (SUNRabs(RHS-LHS) > TOL) {
        alltrue = SUNFALSE;
        break;
      }
    }
    if (alltrue) {
      Q++;
    } else {
      break;
    }
  }

  /* D(R) */
  R = 0;
  for (k=1; k<1000; k++) {
    alltrue = SUNTRUE;
    for (j=0; j<s; j++) {
      LHS = RCONST(0.0);
      for (i=0; i<s; i++)
        LHS += A[i][j]*b[i]*SUNRpowerI(c[i],k-1);
      RHS = b[j]/k*(RCONST(1.0)-SUNRpowerI(c[j],k));
      if (SUNRabs(RHS-LHS) > TOL) {
        alltrue = SUNFALSE;
        break;
      }
    }
    if (alltrue) {
      R++;
    } else {
      break;
    }
  }

  /* determine q, clean up and return */
  q = 0;
  for (i=1; i<=P; i++) {
    if ((q > Q+R+1) || (q > 2*Q+2))
      break;
    q++;
  }
  free(tmp);
  return(q);
}

/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
