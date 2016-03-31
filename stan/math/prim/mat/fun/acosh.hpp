#ifndef STAN_MATH_PRIM_MAT_FUN_ACOSH_HPP
#define STAN_MATH_PRIM_MAT_FUN_ACOSH_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Structure to wrap acosh() so it can be vectorized.
     * @param x Argument variable >= 1.
     * @tparam T Argument type.
     * @return Inverse hyperbolic cosine of variable in radians. 
     */
    struct acosh_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::acosh;
        return acosh(x);
      }
    };

    /**
     * Vectorized version of acosh().
     * @param x Container of variables >= 1.
     * @tparam T Container type.
     * @return Inverse hyperbolic cosine of each variable 
     *         in the container, in radians. 
     */
    template <typename T>
    inline typename apply_scalar_unary<acosh_fun, T>::return_t
    acosh(const T& x) {
      return apply_scalar_unary<acosh_fun, T>::apply(x);
    }

  }
}

#endif
