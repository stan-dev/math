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
 * Implementation header file for ARKode's time step adaptivity 
 * utilities.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_ADAPT_IMPL_H
#define _ARKODE_ADAPT_IMPL_H

#include <stdarg.h>
#include <arkode/arkode.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  ARKode Time Step Adaptivity Data Structure
===============================================================*/

/* size constants for the adaptivity memory structure */
#define ARK_ADAPT_LRW  19
#define ARK_ADAPT_LIW  8    /* includes functin/data pointers */
  
/*---------------------------------------------------------------
 Types : struct ARKodeHAdaptMemRec, ARKodeHAdaptMem
-----------------------------------------------------------------
 The type ARKodeHAdaptMem is type pointer to struct 
 ARKodeHAdaptMemRec.  This structure contains fields to 
 keep track of temporal adaptivity.
---------------------------------------------------------------*/
typedef struct ARKodeHAdaptMemRec {

  realtype     etamax;      /* eta <= etamax                              */
  realtype     etamx1;      /* max step size change on first step         */
  realtype     etamxf;      /* h reduction factor on multiple error fails */
  int          small_nef;   /* bound to determine 'multiple' above        */
  realtype     etacf;       /* h reduction factor on nonlinear conv fail  */
  ARKAdaptFn   HAdapt;      /* function to set the new time step size     */
  void        *HAdapt_data; /* user pointer passed to hadapt              */
  realtype     ehist[3];    /* error history for time adaptivity          */
  realtype     hhist[3];    /* step history for time adaptivity           */
  int          imethod;     /* step adaptivity method to use:
                               -1 -> User-specified function above
                                0 -> PID controller
                                1 -> PI controller
                                2 -> I controller
                                3 -> explicit Gustafsson controller
                                4 -> implicit Gustafsson controller
                                5 -> imex Gustafsson controller           */
  realtype     cfl;         /* cfl safety factor                          */
  realtype     safety;      /* accuracy safety factor on h                */
  realtype     bias;        /* accuracy safety factor on LTE              */
  realtype     growth;      /* maximum step growth safety factor          */
  realtype     lbound;      /* eta lower bound to leave h unchanged       */
  realtype     ubound;      /* eta upper bound to leave h unchanged       */
  realtype     k1;          /* method-specific adaptivity parameters      */
  realtype     k2;
  realtype     k3;

  ARKExpStabFn expstab;     /* step stability function                    */
  void        *estab_data;  /* user pointer passed to expstab             */

  long int     nst_acc;     /* num accuracy-limited internal steps        */
  long int     nst_exp;     /* num stability-limited internal steps       */
  
} *ARKodeHAdaptMem;


/*===============================================================
  ARKode Time Step Adaptivity Routines
===============================================================*/

ARKodeHAdaptMem arkAdaptInit();
void arkPrintAdaptMem(ARKodeHAdaptMem hadapt_mem, FILE *outfile);
int arkAdapt(void* arkode_mem, ARKodeHAdaptMem hadapt_mem,
             N_Vector ycur, realtype tcur, realtype hcur,
             int q, int p, booleantype pq, long int nst);
int arkAdaptPID(ARKodeHAdaptMem hadapt_mem, int k, 
                realtype hcur, realtype *hnew);
int arkAdaptPI(ARKodeHAdaptMem hadapt_mem, int k, 
               realtype hcur, realtype *hnew);
int arkAdaptI(ARKodeHAdaptMem hadapt_mem, int k, 
              realtype hcur, realtype *hnew);
int arkAdaptExpGus(ARKodeHAdaptMem hadapt_mem, int k, 
                   long int nst, realtype hcur, realtype *hnew);
int arkAdaptImpGus(ARKodeHAdaptMem hadapt_mem, int k, 
                   long int nst, realtype hcur, realtype *hnew);
int arkAdaptImExGus(ARKodeHAdaptMem hadapt_mem, int k, 
                    long int nst, realtype hcur, realtype *hnew);

  
#ifdef __cplusplus
}
#endif

#endif
