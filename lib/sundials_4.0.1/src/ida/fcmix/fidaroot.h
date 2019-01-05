/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Alan C. Hindmarsh @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the Fortran interface include file for the rootfinding
 * feature of IDA.
 * -----------------------------------------------------------------
 */

/*
 * ==============================================================================
 *
 *                   FIDAROOT Interface Package
 *
 * The FIDAROOT interface package allows programs written in FORTRAN to
 * use the rootfinding feature of the IDA solver module.
 *
 * The user-callable functions constituting the FIDAROOT package are the
 * following: FIDAROOTINIT, FIDAROOTINFO, and FIDAROOTFREE. The corresponding
 * IDA subroutine called by each interface function is given below.
 *
 *   ------------------      ---------------------
 *  | FIDAROOT routine |    | IDA function called |
 *   ------------------      ---------------------
 *      FIDAROOTINIT     ->     IDARootInit
 *      FIDAROOTINFO     ->     IDAGetRootInfo
 *      FIDAROOTFREE     ->     IDARootInit
 *
 * FIDAROOTFN is a user-supplied subroutine defining the functions whose
 * roots are sought.
 *
 * ==============================================================================
 *
 *                     Usage of the FIDAROOT Interface Package
 *
 * 1. In order to use the rootfinding feature of the IDA package the user must
 * define the following subroutine:
 *
 *   SUBROUTINE FIDAROOTFN (T, Y, YP, G, IPAR, RPAR, IER)
 *   DIMENSION Y(*), YP(*), G(*)
 *
 * The arguments are:
 *   T  = independent variable value t  [input]
 *   Y  = dependent variable vector y  [input]
 *   YP = dependent variable derivative vector y'  [input]
 *   G  = function values g(t,y,y')  [output]
 *   IPAR, RPAR = user (long int and realtype) data [input/output]
 *   IER = return flag (set on 0 if successful, non-zero if an error occurred)
 *
 * 2. After calling FIDAMALLOC but prior to calling FIDASOLVE, the user must
 * allocate and initialize memory for the FIDAROOT module by making the
 * following call:
 *
 *   CALL FIDAROOTINIT (NRTFN, IER)
 *
 * The arguments are:
 *   NRTFN = total number of root functions  [input]
 *   IER   = return completion flag (0 = success, -1 = IDA memory NULL and
 *           -14 = memory allocation error)  [output]
 *
 * 3. After calling FIDA, to see whether a root was found, test the FIDA
 * return flag IER.  The value IER = 2 means one or more roots were found.
 *
 * 4. If a root was found, and if NRTFN > 1, then to determine which root
 * functions G(*) were found to have a root, make the following call:
 *     CALL FIDAROOTINFO (NRTFN, INFO, IER)
 * The arguments are:
 *   NRTFN = total number of root functions  [input]
 *   INFO  = integer array of length NRTFN, with values 0 or 1 [output]
 *           For i = 1,...,NRTFN, G(i) was found to have a root if INFO(i) = 1.
 *   IER   = completion flag (0 = success,  negative = failure)
 *
 * 5. The total number of calls made to the root function (FIDAROOTFN),
 * NGE, can be obtained from IOUT(12).
 *
 * If the FIDA/IDA memory block is reinitialized to solve a different
 * problem via a call to FIDAREINIT, then the counter variable NGE is cleared
 * (reset to zero).
 *
 * 6. To free the memory resources allocated by a prior call to FIDAROOTINIT,
 * make the following call:
 *   CALL FIDAROOTFREE
 * See the IDA documentation for additional information.
 *
 * ==============================================================================
 */

#ifndef _FIDAROOT_H
#define _FIDAROOT_H

/* header files */
#include <sundials/sundials_nvector.h> /* definition of type N_Vector          */
#include <sundials/sundials_types.h>   /* definition of SUNDIALS type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */

#if defined(SUNDIALS_F77_FUNC)

#define FIDA_ROOTINIT SUNDIALS_F77_FUNC(fidarootinit, FIDAROOTINIT)
#define FIDA_ROOTINFO SUNDIALS_F77_FUNC(fidarootinfo, FIDAROOTINFO)
#define FIDA_ROOTFREE SUNDIALS_F77_FUNC(fidarootfree, FIDAROOTFREE)
#define FIDA_ROOTFN   SUNDIALS_F77_FUNC(fidarootfn, FIDAROOTFN)

#else

#define FIDA_ROOTINIT fidarootinit_
#define FIDA_ROOTINFO fidarootinfo_
#define FIDA_ROOTFREE fidarootfree_
#define FIDA_ROOTFN   fidarootfn_

#endif

/* Prototypes of exported function */

void FIDA_ROOTINIT(int *nrtfn, int *ier);
void FIDA_ROOTINFO(int *nrtfn, int *info, int *ier);
void FIDA_ROOTFREE(void);

/* Prototype of function called by IDA module */

int FIDArootfunc(realtype t, N_Vector y, N_Vector yp, realtype *gout,
                 void *user_data);

#ifdef __cplusplus
}
#endif


#endif
