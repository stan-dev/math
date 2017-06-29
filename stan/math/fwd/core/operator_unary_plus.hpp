#ifndef STAN_MATH_FWD_CORE_OPERATOR_UNARY_PLUS_HPP
#define STAN_MATH_FWD_CORE_OPERATOR_UNARY_PLUS_HPP

#include <stan/math/fwd/core/fvar.hpp>

namespace stan {
  namespace math {

    /**
     * Returns the argument.  The unary operator+ exists as a
     * complement to the unary operator-.
     */
    template <typename T>
    inline fvar<T> operator+(const fvar<T>& x) {
      return x;
    }
  }
}
#endif
