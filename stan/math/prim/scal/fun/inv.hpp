#ifndef STAN_MATH_PRIM_SCAL_FUN_INV_HPP
#define STAN_MATH_PRIM_SCAL_FUN_INV_HPP

#include <boost/math/tools/promotion.hpp>

namespace stan {
  namespace math {

    /**
     * Return the inverse of the specified value.  The inverse of x is
     * defined as 1/x.
     *
     * @param x Value to invert.
     * @return Inverse of value.
     */
    inline double inv(double x) {
      return 1.0 / x;
    }
  }
}

#endif
