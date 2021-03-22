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
  ARKode Time Step Adaptivity Private Constants
  ===============================================================*/

/* size constants for the adaptivity memory structure */
#define ARK_ADAPT_LRW  19
#define ARK_ADAPT_LIW  8    /* includes function/data pointers */

/* Time step controller default values */
#define CFLFAC    RCONST(0.5)
#define SAFETY    RCONST(0.96)  /* CVODE uses 1.0  */
#define BIAS      RCONST(1.5)   /* CVODE uses 6.0  */
#define GROWTH    RCONST(20.0)  /* CVODE uses 10.0 */
#define HFIXED_LB RCONST(1.0)   /* CVODE uses 1.0  */
#define HFIXED_UB RCONST(1.5)   /* CVODE uses 1.5  */
#define AD0_K1    RCONST(0.58)  /* PID controller constants */
#define AD0_K2    RCONST(0.21)
#define AD0_K3    RCONST(0.1)
#define AD1_K1    RCONST(0.8)   /* PI controller constants */
#define AD1_K2    RCONST(0.31)
#define AD2_K1    RCONST(1.0)   /* I controller constants */
#define AD3_K1    RCONST(0.367) /* explicit Gustafsson controller */
#define AD3_K2    RCONST(0.268)
#define AD4_K1    RCONST(0.98)  /* implicit Gustafsson controller */
#define AD4_K2    RCONST(0.95)
#define AD5_K1    RCONST(0.367) /* imex Gustafsson controller */
#define AD5_K2    RCONST(0.268)
#define AD5_K3    RCONST(0.95)

#define ETAMX1    RCONST(10000.0)  /* maximum step size change on first step */
#define ETAMXF    RCONST(0.3)      /* step size reduction factor on multiple error
                                      test failures (multiple implies >= SMALL_NEF) */
#define ETAMIN    RCONST(0.1)      /* smallest allowable step size reduction factor
                                      on an error test failure */
#define ETACF     RCONST(0.25)     /* step size reduction factor on nonlinear
                                      convergence failure */
#define SMALL_NEF 2                /* if an error failure occurs and SMALL_NEF <= nef,
                                      then reset  eta = MIN(eta, ETAMXF) */


/*===============================================================
  ARKode Time Step Adaptivity Data Structure
  ===============================================================*/

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
  realtype     etamin;      /* eta >= etamin on error test fail           */
  int          small_nef;   /* bound to determine 'multiple' above        */
  realtype     etacf;       /* h reduction factor on nonlinear conv fail  */
  ARKAdaptFn   HAdapt;      /* function to set the new time step size     */
  void        *HAdapt_data; /* user pointer passed to hadapt              */
  realtype     ehist[2];    /* error history for time adaptivity          */
  realtype     hhist[2];    /* step history for time adaptivity           */
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
  int q;                    /* method order                               */
  int p;                    /* embedding order                            */
  booleantype pq;           /* choice of using p (0) vs q (1)             */

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
             realtype ecur, long int nst);
int arkAdaptPID(ARKodeHAdaptMem hadapt_mem, int k,
                realtype hcur, realtype ecur, realtype *hnew);
int arkAdaptPI(ARKodeHAdaptMem hadapt_mem, int k,
               realtype hcur, realtype ecur, realtype *hnew);
int arkAdaptI(ARKodeHAdaptMem hadapt_mem, int k,
              realtype hcur, realtype ecur, realtype *hnew);
int arkAdaptExpGus(ARKodeHAdaptMem hadapt_mem, int k, long int nst,
                   realtype hcur, realtype ecur, realtype *hnew);
int arkAdaptImpGus(ARKodeHAdaptMem hadapt_mem, int k, long int nst,
                   realtype hcur, realtype ecur, realtype *hnew);
int arkAdaptImExGus(ARKodeHAdaptMem hadapt_mem, int k, long int nst,
                    realtype hcur, realtype ecur, realtype *hnew);


#ifdef __cplusplus
}
#endif

#endif
