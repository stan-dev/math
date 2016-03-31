#ifndef STAN_MATH_PRIM_MAT_FUN_ATANH_HPP
#define STAN_MATH_PRIM_MAT_FUN_ATANH_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Structure to wrap atan() so it can be vectorized.
     * @param x Variable in range [-1, 1].
     * @tparam T Variable type.
     * @return Inverse hyperbolic tangent of x in radians. 
     */
    struct atanh_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::atanh;
        return atanh(x);
      }
    };

    /**
     * Vectorized version of atan().
     * @param x Container of variables in range [-1, 1].
     * @tparam T Container type.
     * @return Inverse hyberbolic tangent of each value in x, in radians. 
     */
    template <typename T>
    inline typename apply_scalar_unary<atanh_fun, T>::return_t
    atanh(const T& x) {
      return apply_scalar_unary<atanh_fun, T>::apply(x);
    }

  }
}

#endif
