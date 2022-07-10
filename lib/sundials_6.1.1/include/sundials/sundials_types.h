/* -----------------------------------------------------------------
 * Programmer(s): Scott Cohen, Alan Hindmarsh, Radu Serban,
 *                Aaron Collier, and Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This header file exports three types: sunrealtype, sunindextype,
 * and sunbooleantype, as well as the constants SUNTRUE and SUNFALSE.
 *
 * Users should include the header file sundials_types.h in every
 * program file and use the exported name sunrealtype instead of
 * float, double or long double.
 *
 * The constants SUNDIALS_SINGLE_PRECISION, SUNDIALS_DOUBLE_PRECISION
 * and SUNDIALS_LONG_DOUBLE_PRECISION indicate the underlying data
 * type of sunrealtype.
 *
 * The legal types for sunrealtype are float, double and long double.
 *
 * The constants SUNDIALS_INT64_T and SUNDIALS_INT32_T indicate
 * the underlying data type of sunindextype -- the integer data type
 * used for vector and matrix indices.
 *
 * Data types are set at the configuration stage.
 *
 * The macro SUN_RCONST gives the user a convenient way to define
 * real-valued literal constants. To use the constant 1.0, for example,
 * the user should write the following:
 *
 *   #define ONE SUN_RCONST(1.0)
 *
 * If sunrealtype is defined as a double, then SUN_RCONST(1.0) expands
 * to 1.0. If sunrealtype is defined as a float, then SUN_RCONST(1.0)
 * expands to 1.0F. If sunrealtype is defined as a long double,
 * then SUN_RCONST(1.0) expands to 1.0L. There is never a need to
 * explicitly cast 1.0 to (sunrealtype). The macro can be used for
 * literal constants only. It cannot be used for expressions.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_TYPES_H
#define _SUNDIALS_TYPES_H

#include <sundials/sundials_config.h>

#include <float.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 *------------------------------------------------------------------
 * Macro _SUNDIALS_STRUCT_
 * The _SUNDIALS_STRUCT_ macro is defined as a `struct` unless
 * generating the SWIG interfaces - in that case it is defined as
 * nothing. This is needed to work around a bug in SWIG which prevents
 * it from properly parsing our generic module structures.
 *------------------------------------------------------------------
 */
#ifdef SWIG
#define _SUNDIALS_STRUCT_
#else
#define _SUNDIALS_STRUCT_ struct
#endif

/*
 *------------------------------------------------------------------
 * Type sunrealtype
 * Macro SUN_RCONST
 * Constants SUN_SMALL_REAL, SUN_BIG_REAL, and SUN_UNIT_ROUNDOFF
 *------------------------------------------------------------------
 */

#if defined(SUNDIALS_SINGLE_PRECISION)

/* deprecated */
typedef float realtype;
# define RCONST(x) x##F
# define BIG_REAL FLT_MAX
# define SMALL_REAL FLT_MIN
# define UNIT_ROUNDOFF FLT_EPSILON

typedef float sunrealtype;
# define SUN_RCONST(x) x##F
# define SUN_BIG_REAL FLT_MAX
# define SUN_SMALL_REAL FLT_MIN
# define SUN_UNIT_ROUNDOFF FLT_EPSILON

#elif defined(SUNDIALS_DOUBLE_PRECISION)

/* deprecated */
typedef double realtype;
# define RCONST(x) x
# define BIG_REAL DBL_MAX
# define SMALL_REAL DBL_MIN
# define UNIT_ROUNDOFF DBL_EPSILON

typedef double sunrealtype;
# define SUN_RCONST(x) x
# define SUN_BIG_REAL DBL_MAX
# define SUN_SMALL_REAL DBL_MIN
# define SUN_UNIT_ROUNDOFF DBL_EPSILON

#elif defined(SUNDIALS_EXTENDED_PRECISION)

/* deprecated */
typedef long double realtype;
# define RCONST(x) x##L
# define BIG_REAL LDBL_MAX
# define SMALL_REAL LDBL_MIN
# define UNIT_ROUNDOFF LDBL_EPSILON

typedef long double sunrealtype;
# define SUN_RCONST(x) x##L
# define SUN_BIG_REAL LDBL_MAX
# define SUN_SMALL_REAL LDBL_MIN
# define SUN_UNIT_ROUNDOFF LDBL_EPSILON

#endif


/*
 *------------------------------------------------------------------
 * Type : sunindextype
 *------------------------------------------------------------------
 * Defines integer type to be used for vector and matrix indices.
 * User can build sundials to use 32- or 64-bit signed integers.
 * If compiler does not support portable data types, the SUNDIALS
 * CMake build system tries to find a type of the desired size.
 *------------------------------------------------------------------
 */

typedef SUNDIALS_INDEX_TYPE sunindextype;

/*
 *------------------------------------------------------------------
 * Type : sunbooleantype
 *------------------------------------------------------------------
 * Constants : SUNFALSE and SUNTRUE
 *------------------------------------------------------------------
 * ANSI C does not have a built-in boolean data type. Below is the
 * definition for a new type called sunbooleantype. The advantage of
 * using the name sunbooleantype (instead of int) is an increase in
 * code readability. It also allows the programmer to make a
 * distinction between int and boolean data. Variables of type
 * sunbooleantype are intended to have only the two values SUNFALSE and
 * SUNTRUE which are defined below to be equal to 0 and 1,
 * respectively.
 *------------------------------------------------------------------
 */

/* deprecated */
#ifndef booleantype
#define booleantype int
#endif

#ifndef sunbooleantype
#define sunbooleantype int
#endif

#ifndef SUNFALSE
#define SUNFALSE 0
#endif

#ifndef SUNTRUE
#define SUNTRUE 1
#endif

#ifdef __cplusplus
}
#endif

#endif  /* _SUNDIALS_TYPES_H */
