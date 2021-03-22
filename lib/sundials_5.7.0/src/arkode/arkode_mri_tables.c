/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation file for ARKode's MRIStepCoupling tables.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif

/* constants */
#define ZERO   RCONST(0.0)
#define ONE    RCONST(1.0)


/*---------------------------------------------------------------
  Returns MRIStepCoupling table structure for pre-set MRI methods.

  Input:  imeth -- integer key for the desired method (see below)

  Allowed 'method' names and properties are listed in the table
  below.

  The 'type' column denotes whether the method is explicit (E),
  or solve-decoupled implicit (ID).

  The 'QP' column denotes whether the coefficients of the method
  are known precisely enough for use in quad precision (128-bit)
  calculations.

     imeth                       order   type    QP
    ------------------------------------------------
     MIS_KW3                     3       E       Y
     MRI_GARK_ERK45a             4       E       Y?
     MRI_GARK_IRK21a             2       ID      Y
     MRI_GARK_ESDIRK34a          4       ID      N?
    ------------------------------------------------

  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_LoadTable(int imethod)
{

  MRIStepCoupling C;
  ARKodeButcherTable B;
  realtype beta;
  B = NULL;
  C = NULL;

  /* fill in coefficients based on method name */
  switch(imethod) {

  case(MIS_KW3):              /* Schlegel et al., JCAM 226:345-357, 2009 */
    B = ARKodeButcherTable_LoadERK(KNOTH_WOLKE_3_3);
    C = MRIStepCoupling_MIStoMRI(B, 3, 0);
    ARKodeButcherTable_Free(B);
    break;

  case(MRI_GARK_ERK45a):      /* A. Sandu, SINUM 57:2300-2327, 2019 */
    C = MRIStepCoupling_Alloc(2,6);
    C->q = 4;
    C->p = 0;
    C->c[1] = RCONST(0.2);
    C->c[2] = RCONST(0.4);
    C->c[3] = RCONST(0.6);
    C->c[4] = RCONST(0.8);
    C->c[5] = ONE;

    C->G[0][1][0] = RCONST(0.2);
    C->G[0][2][0] = -RCONST(53.0)/RCONST(16.0);
    C->G[0][2][1] = RCONST(281.0)/RCONST(80.0);
    C->G[0][3][0] = -RCONST(36562993.0)/RCONST(71394880.0);
    C->G[0][3][1] = RCONST(34903117.0)/RCONST(17848720.0);
    C->G[0][3][2] = -RCONST(88770499.0)/RCONST(71394880.0);
    C->G[0][4][0] = -RCONST(7631593.0)/RCONST(71394880.0);
    C->G[0][4][1] = -RCONST(166232021.0)/RCONST(35697440.0);
    C->G[0][4][2] = RCONST(6068517.0)/RCONST(1519040.0);
    C->G[0][4][3] = RCONST(8644289.0)/RCONST(8924360.0);
    C->G[0][5][0] = RCONST(277061.0)/RCONST(303808.0);
    C->G[0][5][1] = -RCONST(209323.0)/RCONST(1139280.0);
    C->G[0][5][2] = -RCONST(1360217.0)/RCONST(1139280.0);
    C->G[0][5][3] = -RCONST(148789.0)/RCONST(56964.0);
    C->G[0][5][4] = RCONST(147889.0)/RCONST(45120.0);

    C->G[1][2][0] = RCONST(503.0)/RCONST(80.0);
    C->G[1][2][1] = -RCONST(503.0)/RCONST(80.0);
    C->G[1][3][0] = -RCONST(1365537.0)/RCONST(35697440.0);
    C->G[1][3][1] = RCONST(4963773.0)/RCONST(7139488.0);
    C->G[1][3][2] = -RCONST(1465833.0)/RCONST(2231090.0);
    C->G[1][4][0] = RCONST(66974357.0)/RCONST(35697440.0);
    C->G[1][4][1] = RCONST(21445367.0)/RCONST(7139488.0);
    C->G[1][4][2] = -RCONST(3.0);
    C->G[1][4][3] = -RCONST(8388609.0)/RCONST(4462180.0);
    C->G[1][5][0] = -RCONST(18227.0)/RCONST(7520.0);
    C->G[1][5][1] = TWO;
    C->G[1][5][2] = ONE;
    C->G[1][5][3] = RCONST(5.0);
    C->G[1][5][4] = -RCONST(41933.0)/RCONST(7520.0);
    break;

  case(MRI_GARK_IRK21a):   /* A. Sandu, SINUM 57:2300-2327, 2019 */
    B = ARKodeButcherTable_Alloc(3, SUNFALSE);
    B->A[1][0] = ONE;
    B->A[2][0] = RCONST(0.5);
    B->A[2][2] = RCONST(0.5);
    B->b[0] = RCONST(0.5);
    B->b[2] = RCONST(0.5);
    B->c[1] = ONE;
    B->c[2] = ONE;
    B->q=2;
    C = MRIStepCoupling_MIStoMRI(B, 2, 0);
    ARKodeButcherTable_Free(B);
    break;

  case(MRI_GARK_ESDIRK34a):   /* A. Sandu, SINUM 57:2300-2327, 2019 */
    C = MRIStepCoupling_Alloc(1,7);
    beta = RCONST(0.435866521508458999416019);
    C->c[1] = ONE/RCONST(3.0);
    C->c[2] = ONE/RCONST(3.0);
    C->c[3] = TWO/RCONST(3.0);
    C->c[4] = TWO/RCONST(3.0);
    C->c[5] = ONE;
    C->c[6] = ONE;
    C->G[0][1][0] = ONE/RCONST(3.0);
    C->G[0][2][0] = -beta;
    C->G[0][2][2] = beta;
    C->G[0][3][0] = (RCONST(3.0)-RCONST(10.0)*beta)/(RCONST(24.0)*beta-RCONST(6.0));
    C->G[0][3][2] = (RCONST(5.0)-RCONST(18.0)*beta)/(RCONST(6.0)-RCONST(24.0)*beta);
    C->G[0][4][0] = (-RCONST(24.0)*beta*beta+RCONST(6.0)*beta+ONE)/(RCONST(6.0)-RCONST(24.0)*beta);
    C->G[0][4][2] = (-RCONST(48.0)*beta*beta+RCONST(12.0)*beta+ONE)/(RCONST(24.0)*beta-RCONST(6.0));
    C->G[0][4][4] = beta;
    C->G[0][5][0] = (RCONST(3.0)-RCONST(16.0)*beta)/(RCONST(12.0)-RCONST(48.0)*beta);
    C->G[0][5][2] = (RCONST(48.0)*beta*beta-RCONST(21.0)*beta+TWO)/(RCONST(12.0)*beta-RCONST(3.0));
    C->G[0][5][4] = (RCONST(3.0)-RCONST(16.0)*beta)/RCONST(4.0);
    C->G[0][6][0] = -beta;
    C->G[0][6][6] = beta;
    C->q = 3;
    C->p = 0;
    break;

  default:

    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode",
                    "MRIStepCoupling_LoadTable",
                    "Unknown coupling table");
    return(NULL);

  }

  return(C);

}



/*---------------------------------------------------------------
  Routine to allocate an empty MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages)
{
  int i, j;
  MRIStepCoupling MRIC;

  /* Check for legal 'nmat' and 'stages' values */
  if (nmat < 1)   return(NULL);
  if (stages < 1) return(NULL);

  /* Allocate coupling structure */
  MRIC = NULL;
  MRIC = (MRIStepCoupling) malloc(sizeof(struct MRIStepCouplingMem));
  if (MRIC == NULL) return(NULL);

  /* initialize pointers in G structure to NULL */
  MRIC->G = NULL;
  MRIC->c = NULL;

  /* set nmat and stages into G structure */
  MRIC->nmat   = nmat;
  MRIC->stages = stages;

  /* initialize values for method/embedding orders */
  MRIC->q = 0;
  MRIC->p = 0;

  /*
   * Allocate fields within coupling structure
   */

  /* allocate matrices in G */
  MRIC->G = NULL;
  MRIC->G = (realtype ***) calloc( nmat, sizeof(realtype**) );
  if (MRIC->G == NULL) { MRIStepCoupling_Free(MRIC); return(NULL); }

  /* allocate rows of each matrix in G */
  for (i=0; i<nmat; i++) {
    MRIC->G[i] = NULL;
    MRIC->G[i] = (realtype **) calloc( stages, sizeof(realtype*) );
    if (MRIC->G[i] == NULL) { MRIStepCoupling_Free(MRIC); return(NULL); }
  }

  /* allocate columns of each matrix in G */
  for (i=0; i<nmat; i++)
    for (j=0; j<stages; j++) {
      MRIC->G[i][j] = NULL;
      MRIC->G[i][j] = (realtype *) calloc( stages, sizeof(realtype) );
      if (MRIC->G[i][j] == NULL) { MRIStepCoupling_Free(MRIC); return(NULL); }
    }

  /* allocate abcissae */
  MRIC->c = NULL;
  MRIC->c = (realtype *) calloc( stages, sizeof(realtype) );
  if (MRIC->c == NULL) { MRIStepCoupling_Free(MRIC); return(NULL); }

  return(MRIC);
}


/*---------------------------------------------------------------
  Routine to allocate and fill a MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Create(int nmat, int stages,
                                       int q, int p, realtype *G,
                                       realtype *c)
{
  int i, j, k;
  MRIStepCoupling MRIC;

  /* Check for legal 'nmat' and 'stages' values */
  if (nmat < 1)   return(NULL);
  if (stages < 1) return(NULL);

  /* Check that G and c are non-NULL */
  if ((G == NULL) || (c == NULL)) return(NULL);

  /* Allocate MRIStepCoupling structure */
  MRIC = MRIStepCoupling_Alloc(nmat, stages);
  if (MRIC == NULL) return(NULL);

  /* Copy the inputs for method/embedding order */
  MRIC->q = q;
  MRIC->p = p;

  /* Copy the input for G into MRIC; note that G must be a 1D array
     of length nmat*stages*stages, stored in C (row-major) order */
  for (k=0; k<nmat; k++)
    for (i=0; i<stages; i++)
      for (j=0; j<stages; j++)
        MRIC->G[k][i][j] = G[stages*(stages*k + i) + j];

  /* Copy the input for c into MRIC */
  for (i=0; i<stages; i++)
    MRIC->c[i] = c[i];

  return(MRIC);
}


/*---------------------------------------------------------------
  Construct the MRI 'Gamma' matrix for an MIS method based on
  a given 'slow' Butcher table.
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B,
                                         int q, int p)
{
  int i, j, stages;
  booleantype padding;
  realtype Asum;
  MRIStepCoupling MRIC;
  const realtype tol = RCONST(100.0)*UNIT_ROUNDOFF;

  /* Check that input table is non-NULL */
  if (B == NULL)  return(NULL);

  /* Check that the input table is valid */
  /*     First stage is just old solution */
  Asum = SUNRabs(B->c[0]);
  for (j=0; j<B->stages; j++)  Asum += SUNRabs(B->A[0][j]);
  if (Asum > tol)  return(NULL);
  /*     Last stage exceeds 1 */
  if (B->c[B->stages-1] > ONE+tol)  return(NULL);
  /*     All stages are sorted */
  for (j=1; j<B->stages; j++)
    if ((B->c[j] - B->c[j-1]) < -tol)  return(NULL);
  /*     Each stage at most diagonally implicit */
  Asum = ZERO;
  for (i=0; i<B->stages; i++)
    for (j=i+1; j<B->stages; j++)
      Asum += SUNRabs(B->A[i][j]);
  if (Asum > tol)  return(NULL);

  /* determine whether the table needs padding */
  padding = SUNFALSE;
  /*     Last stage time should equal 1 */
  if (SUNRabs(B->c[B->stages-1] - ONE) > tol) {
    padding = SUNTRUE;
  }
  /*     Last row of A should equal b */
  for (j=0; j<B->stages; j++) {
    if (SUNRabs(B->A[B->stages-1][j] - B->b[j]) > tol)
      padding = SUNTRUE;
  }
  stages = (padding) ? B->stages+1 : B->stages;

  /*  construct MRIC structure */
  MRIC = MRIStepCoupling_Alloc(1, stages);
  if (MRIC == NULL) return(NULL);

  /* Copy method/embedding orders */
  MRIC->q = q;
  MRIC->p = p;

  /* Copy c from Butcher table into MRIC, and
     adjust for padding if needed */
  for (i=0; i<B->stages; i++)  MRIC->c[i] = B->c[i];
  if (padding)  MRIC->c[stages-1] = ONE;

  /* Construct G:
     First row is identically zero
     Remaining rows = A(2:end,:) - A(1:end-1,:)
     Padded row = b(:) - A(end,:) */
  for (i=0; i<stages; i++)
    for (j=0; j<stages; j++)
      MRIC->G[0][i][j] = ZERO;
  for (i=1; i<B->stages; i++)
    for (j=0; j<B->stages; j++)
      MRIC->G[0][i][j] = B->A[i][j] - B->A[i-1][j];
  if (padding)
    for (j=0; j<B->stages; j++)
      MRIC->G[0][stages-1][j] = B->b[j] - B->A[B->stages-1][j];

  return(MRIC);
}

/*---------------------------------------------------------------
  Routine to copy a MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling MRIC)
{
  int i, j, k, nmat, stages;
  MRIStepCoupling MRICcopy;

  /* Check for legal input */
  if (MRIC == NULL) return(NULL);

  /* Get the number of coupling matrices and stages */
  nmat = MRIC->nmat;
  stages = MRIC->stages;

  /* Allocate MRIStepCoupling structure */
  MRICcopy = MRIStepCoupling_Alloc(nmat, stages);
  if (MRICcopy == NULL) return(NULL);

  /* Copy MRIStepCoupling matrices */
  for (k=0; k<nmat; k++)
    for (i=0; i<stages; i++)
      for (j=0; j<stages; j++)
        MRICcopy->G[k][i][j] = MRIC->G[k][i][j];

  /* Copy MRIStepCoupling abcissae */
  for (i=0; i<stages; i++)
    MRICcopy->c[i] = MRIC->c[i];

  /* Copy MRIStepCoupling method/embedding orders */
  MRICcopy->q = MRIC->q;
  MRICcopy->p = MRIC->p;

  return(MRICcopy);
}


/*---------------------------------------------------------------
  Routine to query the MRIStepCoupling structure workspace size
  ---------------------------------------------------------------*/
void MRIStepCoupling_Space(MRIStepCoupling MRIC, sunindextype *liw,
                           sunindextype *lrw)
{
  /* initialize outputs and return if MRIC is not allocated */
  *liw = 0;  *lrw = 0;
  if (MRIC == NULL)  return;

  /* fill outputs based on MRIC */
  *liw = 4;
  *lrw = MRIC->stages*(1 + MRIC->nmat*MRIC->stages);
}


/*---------------------------------------------------------------
  Routine to free a MRIStepCoupling structure
  ---------------------------------------------------------------*/
void MRIStepCoupling_Free(MRIStepCoupling MRIC)
{
  int k, i;

  /* Free each field within MRIStepCoupling structure, and then
     free structure itself */
  if (MRIC != NULL) {
    if (MRIC->c != NULL)  free(MRIC->c);
    if (MRIC->G != NULL) {
      for (k=0; k<MRIC->nmat; k++)
        if (MRIC->G[k] != NULL) {
          for (i=0; i<MRIC->stages; i++)
            if (MRIC->G[k][i] != NULL) {
              free(MRIC->G[k][i]);
              MRIC->G[k][i] = NULL;
            }
          free(MRIC->G[k]);
          MRIC->G[k] = NULL;
        }
      free(MRIC->G);
    }
    free(MRIC);
  }
}


/*---------------------------------------------------------------
  Routine to print a MRIStepCoupling structure
  ---------------------------------------------------------------*/
void MRIStepCoupling_Write(MRIStepCoupling MRIC, FILE *outfile)
{
  int i, j, k;

  /* check for vaild coupling structure */
  if (MRIC == NULL) return;
  if (MRIC->G == NULL) return;
  for (i=0; i<MRIC->nmat; i++) {
    if (MRIC->G[i] == NULL) return;
    for (j=0; j<MRIC->stages; j++)
      if (MRIC->G[i][j] == NULL) return;
  }
  if (MRIC->c == NULL) return;

  fprintf(outfile, "  nmat = %i\n", MRIC->nmat);
  fprintf(outfile, "  stages = %i\n", MRIC->stages);
  fprintf(outfile, "  method order (q) = %i\n", MRIC->q);
  fprintf(outfile, "  embedding order (p) = %i\n", MRIC->p);

  fprintf(outfile, "  c = ");
  for (i=0; i<MRIC->stages; i++)
    fprintf(outfile, "%"RSYM"  ", MRIC->c[i]);
  fprintf(outfile, "\n");

  for (k=0; k<MRIC->nmat; k++) {
    fprintf(outfile, "  G[%i] = \n",k);
    for (i=0; i<MRIC->stages; i++) {
      fprintf(outfile, "      ");
      for (j=0; j<MRIC->stages; j++)
        fprintf(outfile, "%"RSYM"  ", MRIC->G[k][i][j]);
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
  }

}


/*===============================================================
  EOF
  ===============================================================*/
