#ifndef STAN_MATH_REV_SCAL_FUN_LDEXP_HPP
#define STAN_MATH_REV_SCAL_FUN_LDEXP_HPP

#include <stan/math/rev/core.hpp>
#include <cmath>

namespace stan {
  namespace math {
      
      /**
       * Returns the product of a (the significand) and 
       * 2 to power b (the exponent).
       *
       * @param a the significand
       * @param b an integer that is the exponent
       * @return product of a times 2 to the power b
       */
      
      inline var ldexp(var a, int b) {
        return a * std::pow(2.0, b);
      }
  }
}
#endif
