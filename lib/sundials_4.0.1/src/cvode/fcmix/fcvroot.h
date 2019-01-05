/*
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
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
 * feature of CVODE.
 * -----------------------------------------------------------------
 */

/*
 * ==============================================================================
 *
 *                   FCVROOT Interface Package
 *
 * The FCVROOT interface package allows programs written in FORTRAN to
 * use the rootfinding feature of the CVODE solver module.
 *
 * The user-callable functions constituting the FCVROOT package are the
 * following: FCVROOTINIT, FCVROOTINFO, and FCVROOTFREE. The corresponding
 * CVODE subroutine called by each interface function is given below.
 *
 *   -----------------      -----------------------
 *  | FCVROOT routine |    | CVODE function called |
 *   -----------------      -----------------------
 *      FCVROOTINIT     ->     CVodeRootInit
 *      FCVROOTINFO     ->     CVodeGetRootInfo
 *      FCVROOTFREE     ->     CVodeRootInit
 *
 * FCVROOTFN is a user-supplied subroutine defining the functions whose
 * roots are sought.
 *
 * ==============================================================================
 *
 *                     Usage of the FCVROOT Interface Package
 *
 * 1. In order to use the rootfinding feature of the CVODE package the user must
 * define the following subroutine:
 *
 *   SUBROUTINE FCVROOTFN (T, Y, G, IPAR, RPAR, IER)
 *   DIMENSION Y(*), G(*), IPAR(*), RPAR(*)
 *
 * The arguments are:
 *   T = independent variable value t  [input]
 *   Y = dependent variable vector y  [input]
 *   G = function values g(t,y)  [output]
 *   IPAR, RPAR = user (long int and realtype) data [input/output]
 *   IER = return flag (0 for success, a non-zero value if an error occurred.)
 *
 * 2. After calling FCVMALLOC but prior to calling FCVODE, the user must
 * allocate and initialize memory for the FCVROOT module by making the
 * following call:
 *
 *   CALL FCVROOTINIT (NRTFN, IER)
 *
 * The arguments are:
 *   NRTFN = total number of root functions  [input]
 *   IER   = return completion flag (0 = success, -1 = CVODE memory NULL and
 *           -11  memory allocation error)  [output]
 *
 * 3. After calling FCVODE, to see whether a root was found, test the FCVODE
 * return flag IER.  The value IER = 2 means one or more roots were found.
 *
 * 4. If a root was found, and if NRTFN > 1, then to determine which root
 * functions G(*) were found to have a root, make the following call:
 *     CALL FCVROOTINFO (NRTFN, INFO, IER)
 * The arguments are:
 *   NRTFN = total number of root functions  [input]
 *   INFO  = integer array of length NRTFN, with values 0 or 1 [output]
 *           For i = 1,...,NRTFN, G(i) was found to have a root if INFO(i) = 1.
 *   IER   = completion flag (0 = success,  negative = failure)
 *
 * 5. The total number of calls made to the root function (FCVROOTFN), NGE,
 * can be obtained from IOUT(12).
 *
 * If the FCVODE/CVODE memory block is reinitialized to solve a different
 * problem via a call to FCVREINIT, then the counter variable NGE is cleared
 * (reset to zero).
 *
 * 6. To free the memory resources allocated by a prior call to FCVROOTINIT make
 * the following call:
 *   CALL FCVROOTFREE
 * See the CVODE documentation for additional information.
 *
 * ==============================================================================
 */

#ifndef _FCVROOT_H
#define _FCVROOT_H

/* header files */
#include <sundials/sundials_nvector.h> /* definition of type N_Vector          */
#include <sundials/sundials_types.h>   /* definition of SUNDIALS type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */

#if defined(SUNDIALS_F77_FUNC)

#define FCV_ROOTINIT SUNDIALS_F77_FUNC(fcvrootinit, FCVROOTINIT)
#define FCV_ROOTINFO SUNDIALS_F77_FUNC(fcvrootinfo, FCVROOTINFO)
#define FCV_ROOTFREE SUNDIALS_F77_FUNC(fcvrootfree, FCVROOTFREE)
#define FCV_ROOTFN   SUNDIALS_F77_FUNC(fcvrootfn, FCVROOTFN)

#else

#define FCV_ROOTINIT fcvrootinit_
#define FCV_ROOTINFO fcvrootinfo_
#define FCV_ROOTFREE fcvrootfree_
#define FCV_ROOTFN   fcvrootfn_

#endif

/* Prototypes of exported function */

void FCV_ROOTINIT(int *nrtfn, int *ier);
void FCV_ROOTINFO(int *nrtfn, int *info, int *ier);
void FCV_ROOTFREE(void);

/* Prototype of function called by CVODE module */

int FCVrootfunc(realtype t, N_Vector y, realtype *gout, void *user_data);

#ifdef __cplusplus
}
#endif


#endif
