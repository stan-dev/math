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
 * This is the implementation file for ARKode's time step
 * adaptivity utilities.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif



/*---------------------------------------------------------------
  arkAdaptInit:

  This routine creates and sets default values in an
  ARKodeHAdaptMem structure.  This returns a non-NULL structure
  if no errors occurred, or a NULL value otherwise.
  ---------------------------------------------------------------*/
ARKodeHAdaptMem arkAdaptInit()
{
  ARKodeHAdaptMem hadapt_mem;

  /* allocate structure */
  hadapt_mem = (ARKodeHAdaptMem) malloc(sizeof(struct ARKodeHAdaptMemRec));
  if (hadapt_mem == NULL)  return(NULL);

  /* initialize default values */
  memset(hadapt_mem, 0, sizeof(struct ARKodeHAdaptMemRec));
  hadapt_mem->etamx1      = ETAMX1;     /* max change on first step */
  hadapt_mem->etamxf      = ETAMXF;     /* max change on error-failed step */
  hadapt_mem->small_nef   = SMALL_NEF;  /* num error fails before ETAMXF enforced */
  hadapt_mem->etacf       = ETACF;      /* max change on convergence failure */
  hadapt_mem->HAdapt      = NULL;       /* step adaptivity fn */
  hadapt_mem->HAdapt_data = NULL;       /* step adaptivity data */
  hadapt_mem->imethod     = 0;          /* PID controller */
  hadapt_mem->cfl         = CFLFAC;     /* explicit stability factor */
  hadapt_mem->safety      = SAFETY;     /* step adaptivity safety factor  */
  hadapt_mem->bias        = BIAS;       /* step adaptivity error bias */
  hadapt_mem->growth      = GROWTH;     /* step adaptivity growth factor */
  hadapt_mem->lbound      = HFIXED_LB;  /* step adaptivity no-change lower bound */
  hadapt_mem->ubound      = HFIXED_UB;  /* step adaptivity no-change upper bound */
  hadapt_mem->k1          = AD0_K1;     /* step adaptivity parameter */
  hadapt_mem->k2          = AD0_K2;     /* step adaptivity parameter */
  hadapt_mem->k3          = AD0_K3;     /* step adaptivity parameter */
  hadapt_mem->ehist[0]    = ONE;
  hadapt_mem->ehist[1]    = ONE;
  hadapt_mem->ehist[2]    = ONE;
  hadapt_mem->hhist[0]    = ZERO;
  hadapt_mem->hhist[1]    = ZERO;
  hadapt_mem->hhist[2]    = ZERO;
  hadapt_mem->nst_acc     = 0;
  hadapt_mem->nst_exp     = 0;

  hadapt_mem->expstab     = arkExpStab;
  hadapt_mem->estab_data  = NULL;
  return(hadapt_mem);
}


/*---------------------------------------------------------------
  arkPrintAdaptMem

  This routine outputs the time step adaptivity memory structure
  to a specified file pointer.
  ---------------------------------------------------------------*/
void arkPrintAdaptMem(ARKodeHAdaptMem hadapt_mem, FILE *outfile)
{
  if (hadapt_mem != NULL) {
    fprintf(outfile, "ark_hadapt: etamax = %"RSYM"\n", hadapt_mem->etamax);
    fprintf(outfile, "ark_hadapt: etamx1 = %"RSYM"\n", hadapt_mem->etamx1);
    fprintf(outfile, "ark_hadapt: etamxf = %"RSYM"\n", hadapt_mem->etamxf);
    fprintf(outfile, "ark_hadapt: small_nef = %i\n", hadapt_mem->small_nef);
    fprintf(outfile, "ark_hadapt: etacf = %"RSYM"\n", hadapt_mem->etacf);
    fprintf(outfile, "ark_hadapt: imethod = %i\n", hadapt_mem->imethod);
    fprintf(outfile, "ark_hadapt: ehist =  %"RSYM"  %"RSYM"  %"RSYM"\n",
            hadapt_mem->ehist[0],
            hadapt_mem->ehist[1],
            hadapt_mem->ehist[2]);
    fprintf(outfile, "ark_hadapt: hhist =  %"RSYM"  %"RSYM"  %"RSYM"\n",
            hadapt_mem->hhist[0],
            hadapt_mem->hhist[1],
            hadapt_mem->hhist[2]);
    fprintf(outfile, "ark_hadapt: cfl = %"RSYM"\n", hadapt_mem->cfl);
    fprintf(outfile, "ark_hadapt: safety = %"RSYM"\n", hadapt_mem->safety);
    fprintf(outfile, "ark_hadapt: bias = %"RSYM"\n", hadapt_mem->bias);
    fprintf(outfile, "ark_hadapt: growth = %"RSYM"\n", hadapt_mem->growth);
    fprintf(outfile, "ark_hadapt: lbound = %"RSYM"\n", hadapt_mem->lbound);
    fprintf(outfile, "ark_hadapt: ubound = %"RSYM"\n", hadapt_mem->ubound);
    fprintf(outfile, "ark_hadapt: k1 = %"RSYM"\n", hadapt_mem->k1);
    fprintf(outfile, "ark_hadapt: k2 = %"RSYM"\n", hadapt_mem->k2);
    fprintf(outfile, "ark_hadapt: k3 = %"RSYM"\n", hadapt_mem->k3);
    fprintf(outfile, "ark_hadapt: nst_acc = %li\n", hadapt_mem->nst_acc);
    fprintf(outfile, "ark_hadapt: nst_exp = %li\n", hadapt_mem->nst_exp);
    if (hadapt_mem->expstab == arkExpStab) {
      fprintf(outfile, "  ark_hadapt: Default explicit stability function\n");
    } else {
      fprintf(outfile, "  ark_hadapt: User provided explicit stability function\n");
      fprintf(outfile, "  ark_hadapt: stability function data pointer = %p\n",
              hadapt_mem->estab_data);
    }
  }
}




/*---------------------------------------------------------------
  arkAdapt is the time step adaptivity wrapper function.  This
  computes and sets the value of ark_eta inside of the ARKodeMem
  data structure.
  ---------------------------------------------------------------*/
int arkAdapt(void* arkode_mem, ARKodeHAdaptMem hadapt_mem,
             N_Vector ycur, realtype tcur, realtype hcur,
             int q, int p, booleantype pq, long int nst)
{
  int ier, k;
  realtype h_acc, h_cfl, int_dir;
  ARKodeMem ark_mem;
  if (arkode_mem == NULL) {
    arkProcessError(NULL, ARK_MEM_NULL, "ARKode",
                    "arkAdapt", MSG_ARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Set k as either p or q, based on pq flag */
  k = (pq) ? q : p;
  
  /* Call algorithm-specific error adaptivity method */
  switch (hadapt_mem->imethod) {
  case(0):    /* PID controller */
    ier = arkAdaptPID(hadapt_mem, k, hcur, &h_acc);
    break;
  case(1):    /* PI controller */
    ier = arkAdaptPI(hadapt_mem, k, hcur, &h_acc);
    break;
  case(2):    /* I controller */
    ier = arkAdaptI(hadapt_mem, k, hcur, &h_acc);
    break;
  case(3):    /* explicit Gustafsson controller */
    ier = arkAdaptExpGus(hadapt_mem, k, nst, hcur, &h_acc);
    break;
  case(4):    /* implicit Gustafsson controller */
    ier = arkAdaptImpGus(hadapt_mem, k, nst, hcur, &h_acc);
    break;
  case(5):    /* imex Gustafsson controller */
    ier = arkAdaptImExGus(hadapt_mem, k, nst, hcur, &h_acc);
    break;
  case(-1):   /* user-supplied controller */
    ier = hadapt_mem->HAdapt(ycur, tcur,
                             hadapt_mem->hhist[0],
                             hadapt_mem->hhist[1],
                             hadapt_mem->hhist[2],
                             hadapt_mem->ehist[0],
                             hadapt_mem->ehist[1],
                             hadapt_mem->ehist[2],
                             q, p, &h_acc, hadapt_mem->HAdapt_data);
    break;
  default:
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkAdapt",
                    "Illegal imethod.");
    return (ARK_ILL_INPUT);
  }
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkAdapt",
                    "Error in accuracy-based adaptivity function.");
    return (ARK_ILL_INPUT);
  }

  /* determine direction of integration */
  int_dir = hcur / SUNRabs(hcur);

  /* Call explicit stability function */
  ier = hadapt_mem->expstab(ycur, tcur, &h_cfl, hadapt_mem->estab_data);
  if (ier != ARK_SUCCESS) {
    arkProcessError(ark_mem, ARK_ILL_INPUT, "ARKode", "arkAdapt",
                    "Error in explicit stability function.");
    return (ARK_ILL_INPUT);
  }
  if (h_cfl <= 0.0)  h_cfl = RCONST(1.0e30) * SUNRabs(hcur);

  /* Solver diagnostics reporting */
  if (ark_mem->report)
    fprintf(ark_mem->diagfp, "ARKadapt  adapt  %"RSYM"  %"RSYM"  %"RSYM"  %"RSYM"  %"RSYM"  %"RSYM"  %"RSYM"  %"RSYM"  ",
            hadapt_mem->ehist[0], hadapt_mem->ehist[1],
            hadapt_mem->ehist[2], hadapt_mem->hhist[0],
            hadapt_mem->hhist[1], hadapt_mem->hhist[2], h_acc, h_cfl);

  /* enforce safety factors */
  h_acc *= hadapt_mem->safety;
  h_cfl *= hadapt_mem->cfl * int_dir;

  /* enforce maximum bound on time step growth */
  h_acc = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(hadapt_mem->etamax*hcur));

  /* enforce minimum bound time step reduction */
  h_acc = int_dir * SUNMAX(SUNRabs(h_acc), SUNRabs(ETAMIN*hcur));

  /* Solver diagnostics reporting */
  if (ark_mem->report)
    fprintf(ark_mem->diagfp, "%"RSYM"  %"RSYM"  ", h_acc, h_cfl);

  /* increment the relevant step counter, set desired step */
  if (SUNRabs(h_acc) < SUNRabs(h_cfl))
    hadapt_mem->nst_acc++;
  else
    hadapt_mem->nst_exp++;
  h_acc = int_dir * SUNMIN(SUNRabs(h_acc), SUNRabs(h_cfl));

  /* enforce adaptivity bounds to retain Jacobian/preconditioner accuracy */
  if ( (SUNRabs(h_acc) > SUNRabs(hcur*hadapt_mem->lbound*ONEMSM)) &&
       (SUNRabs(h_acc) < SUNRabs(hcur*hadapt_mem->ubound*ONEPSM)) )
    h_acc = hcur;

  /* set basic value of ark_eta */
  ark_mem->eta = h_acc / hcur;

  /* enforce minimum time step size */
  ark_mem->eta = SUNMAX(ark_mem->eta,
                        ark_mem->hmin / SUNRabs(hcur));

  /* enforce maximum time step size */
  ark_mem->eta /= SUNMAX(ONE, SUNRabs(hcur) *
                         ark_mem->hmax_inv*ark_mem->eta);

  /* Solver diagnostics reporting */
  if (ark_mem->report)
    fprintf(ark_mem->diagfp, "%"RSYM"\n", ark_mem->eta);

  return(ier);
}


/*---------------------------------------------------------------
  arkAdaptPID implements a PID time step control algorithm.
  ---------------------------------------------------------------*/
int arkAdaptPID(ARKodeHAdaptMem hadapt_mem, int k, realtype hcur,
                realtype *hnew)
{
  realtype k1, k2, k3, e1, e2, e3, h_acc;

  /* set usable time-step adaptivity parameters */
  k1 = -hadapt_mem->k1 / k;
  k2 =  hadapt_mem->k2 / k;
  k3 = -hadapt_mem->k3 / k;
  e1 = SUNMAX(hadapt_mem->ehist[0], TINY);
  e2 = SUNMAX(hadapt_mem->ehist[1], TINY);
  e3 = SUNMAX(hadapt_mem->ehist[2], TINY);

  /* compute estimated optimal time step size, set into output */
  h_acc = hcur * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2) * SUNRpowerR(e3,k3);
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkAdaptPI implements a PI time step control algorithm.
  ---------------------------------------------------------------*/
int arkAdaptPI(ARKodeHAdaptMem hadapt_mem, int k, realtype hcur,
               realtype *hnew)
{
  realtype k1, k2, e1, e2, h_acc;

  /* set usable time-step adaptivity parameters */
  k1 = -hadapt_mem->k1 / k;
  k2 =  hadapt_mem->k2 / k;
  e1 = SUNMAX(hadapt_mem->ehist[0], TINY);
  e2 = SUNMAX(hadapt_mem->ehist[1], TINY);

  /* compute estimated optimal time step size, set into output */
  h_acc = hcur * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkAdaptI implements an I time step control algorithm.
  ---------------------------------------------------------------*/
int arkAdaptI(ARKodeHAdaptMem hadapt_mem, int k, realtype hcur,
              realtype *hnew)
{
  realtype k1, e1, h_acc;

  /* set usable time-step adaptivity parameters */
  k1 = -hadapt_mem->k1 / k;
  e1 = SUNMAX(hadapt_mem->ehist[0], TINY);

  /* compute estimated optimal time step size, set into output */
  h_acc = hcur * SUNRpowerR(e1,k1);
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkAdaptExpGus implements the explicit Gustafsson time step
  control algorithm.
  ---------------------------------------------------------------*/
int arkAdaptExpGus(ARKodeHAdaptMem hadapt_mem, int k, long int nst,
                   realtype hcur, realtype *hnew)
{
  realtype k1, k2, e1, e2, h_acc;

  /* modified method for first step */
  if (nst < 2) {

    k1 = -ONE / k;
    e1 = SUNMAX(hadapt_mem->ehist[0], TINY);
    h_acc = hcur * SUNRpowerR(e1,k1);

  /* general estimate */
  } else {

    k1 = -hadapt_mem->k1 / k;
    k2 = -hadapt_mem->k2 / k;
    e1 = SUNMAX(hadapt_mem->ehist[0], TINY);
    e2 = e1 / SUNMAX(hadapt_mem->ehist[1], TINY);
    h_acc = hcur * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);

  }
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkAdaptImpGus implements the implicit Gustafsson time step
  control algorithm.
  ---------------------------------------------------------------*/
int arkAdaptImpGus(ARKodeHAdaptMem hadapt_mem, int k, long int nst,
                   realtype hcur, realtype *hnew)
{
  realtype k1, k2, e1, e2, hrat, h_acc;

  /* modified method for first step */
  if (nst < 2) {

    k1 = -ONE / k;
    e1 = SUNMAX(hadapt_mem->ehist[0], TINY);
    h_acc = hcur * SUNRpowerR(e1,k1);

  /* general estimate */
  } else {

    k1 = -hadapt_mem->k1 / k;
    k2 = -hadapt_mem->k2 / k;
    e1 = SUNMAX(hadapt_mem->ehist[0], TINY);
    e2 = e1 / SUNMAX(hadapt_mem->ehist[1], TINY);
    hrat = hcur / hadapt_mem->hhist[1];
    h_acc = hcur * hrat * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2);

  }
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  arkAdaptImExGus implements a combination implicit/explicit
  Gustafsson time step control algorithm.
  ---------------------------------------------------------------*/
int arkAdaptImExGus(ARKodeHAdaptMem hadapt_mem, int k, long int nst,
                    realtype hcur, realtype *hnew)
{
  realtype k1, k2, k3, e1, e2, hrat, h_acc;

  /* modified method for first step */
  if (nst < 2) {

    k1 = -ONE / k;
    e1 = SUNMAX(hadapt_mem->ehist[0], TINY);
    h_acc = hcur * SUNRpowerR(e1,k1);

  /* general estimate */
  } else {

    k1 = -hadapt_mem->k1 / k;
    k2 = -hadapt_mem->k2 / k;
    k3 = -hadapt_mem->k3 / k;
    e1 = SUNMAX(hadapt_mem->ehist[0], TINY);
    e2 = e1 / SUNMAX(hadapt_mem->ehist[1], TINY);
    hrat = hcur / hadapt_mem->hhist[1];
    /* implicit estimate */
    h_acc = hcur * hrat * SUNRpowerR(e1,k3) * SUNRpowerR(e2,k3);
    /* explicit estimate */
    h_acc = SUNMIN(h_acc, hcur * SUNRpowerR(e1,k1) * SUNRpowerR(e2,k2));

  }
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*===============================================================
  EOF
  ===============================================================*/
