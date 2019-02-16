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
 * Implementation header file for ARKode's root-finding (in time) 
 * utility.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_ROOT_IMPL_H
#define _ARKODE_ROOT_IMPL_H

#include <stdarg.h>
#include <arkode/arkode.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  ARKode Root-finding constants
  ===============================================================*/

#define ARK_ROOT_LRW   5
#define ARK_ROOT_LIW  12   /* int, ptr, etc */
  
/*===============================================================
  ARKode Root-finding Data Structure
  ===============================================================*/

/*---------------------------------------------------------------
  Types : struct ARKodeRootMemRec, ARKodeRootMem
  -----------------------------------------------------------------
  The type ARKodeRootMem is type pointer to struct 
  ARKodeRootMemRec.  This structure contains data pertaining to 
  the use of root-finding capabilities in ARKode.
  ---------------------------------------------------------------*/
typedef struct ARKodeRootMemRec {

  ARKRootFn    gfun;        /* function g for roots sought                  */
  int          nrtfn;       /* number of components of g                    */
  int         *iroots;      /* array for root information                   */
  int         *rootdir;     /* array specifying direction of zero-crossing  */
  realtype     tlo;         /* nearest endpoint of interval in root search  */
  realtype     thi;         /* farthest endpoint of interval in root search */
  realtype     trout;       /* t value returned by rootfinding routine      */
  realtype    *glo;         /* saved array of g values at t = tlo           */
  realtype    *ghi;         /* saved array of g values at t = thi           */
  realtype    *grout;       /* array of g values at t = trout               */
  realtype     toutc;       /* copy of tout (if NORMAL mode)                */
  realtype     ttol;        /* tolerance on root location                   */
  int          taskc;       /* copy of parameter itask                      */
  int          irfnd;       /* flag showing whether last step had a root    */
  long int     nge;         /* counter for g evaluations                    */
  booleantype *gactive;     /* array with active/inactive event functions   */
  int          mxgnull;     /* num. warning messages about possible g==0    */

} *ARKodeRootMem;


/*===============================================================
  ARKode Root-finding Routines
===============================================================*/

int arkRootFree(void* arkode_mem);
int arkPrintRootMem(void* arkode_mem, FILE *outfile);
int arkRootCheck1(void* arkode_mem);
int arkRootCheck2(void* arkode_mem);
int arkRootCheck3(void* arkode_mem);
int arkRootfind(void* arkode_mem);

  
#ifdef __cplusplus
}
#endif

#endif
