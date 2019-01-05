/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Aaron Collier @ LLNL
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
 * This is the implementation file for a simple C-language math
 * library.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

realtype SUNRpowerI(realtype base, int exponent)
{
  int i, expt;
  realtype prod;

  prod = ONE;
  expt = abs(exponent);
  for(i = 1; i <= expt; i++) prod *= base;
  if (exponent < 0) prod = ONE/prod;
  return(prod);
}

realtype SUNRpowerR(realtype base, realtype exponent)
{
  if (base <= ZERO) return(ZERO);

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
