/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Aaron Collier @ LLNL
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
 * This is the implementation file for a simple C-language math
 * library.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sundials/sundials_math.h>

static booleantype sunIsInf(realtype a)
{
#if (__STDC_VERSION__ >= 199901L)
  return(isinf(a));
#else
  return(a < -BIG_REAL || a > BIG_REAL);
#endif
}

static booleantype sunIsNaN(realtype a)
{
#if ( __STDC_VERSION__ >= 199901L)
  return(isnan(a));
#else
  /* Most compilers/platforms follow NaN != a,
   * but since C89 does not require this, it is
   * possible some platforms might not follow it.
   */
  return(a != a);
#endif
}

realtype SUNRpowerI(realtype base, int exponent)
{
  int i, expt;
  realtype prod;

  prod = RCONST(1.0);
  expt = abs(exponent);
  for(i = 1; i <= expt; i++) prod *= base;
  if (exponent < 0) prod = RCONST(1.0)/prod;
  return(prod);
}

realtype SUNRpowerR(realtype base, realtype exponent)
{
  if (base <= RCONST(0.0)) return(RCONST(0.0));

#if defined(SUNDIALS_USE_GENERIC_MATH)
  return((realtype) pow((double) base, (double) exponent));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  return(pow(base, exponent));
#elif defined(SUNDIALS_SINGLE_PRECISION)
  return(powf(base, exponent));
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  return(powl(base, exponent));
#endif
}

booleantype SUNRCompare(realtype a, realtype b)
{
  return(SUNRCompareTol(a, b, 10*UNIT_ROUNDOFF));
}

booleantype SUNRCompareTol(realtype a, realtype b, realtype tol)
{
  realtype diff;
  realtype norm;

  /* If a and b are exactly equal.
   * This also covers the case where a and b are both inf under IEEE 754.
   */
  if (a == b) return(SUNFALSE);

  /* If a or b are NaN */
  if (sunIsNaN(a) || sunIsNaN(b)) return(SUNTRUE);

  /* If one of a or b are Inf (since we handled both being inf above) */
  if (sunIsInf(a) || sunIsInf(b)) return(SUNTRUE);

  diff = SUNRabs(a - b);
  norm = SUNMIN(SUNRabs(a + b), BIG_REAL);

  /* When |a + b| is very small (less than 10*UNIT_ROUNDOFF) or zero, we use an
   * absolute difference:
   *    |a - b| >= 10*UNIT_ROUNDOFF
   * Otherwise we use a relative difference:
   *    |a - b| < tol * |a + b|
   * The choice to use |a + b| over max(a, b)
   * is arbitrary, as is the choice to use
   * 10*UNIT_ROUNDOFF.
   */
  return(diff >= SUNMAX(10*UNIT_ROUNDOFF, tol*norm));
}
