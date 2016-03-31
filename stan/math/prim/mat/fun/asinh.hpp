#ifndef STAN_MATH_PRIM_MAT_FUN_ASINH_HPP
#define STAN_MATH_PRIM_MAT_FUN_ASINH_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Structure to wrap asinh() so it can be vectorized.
     * @param x Variable.
     * @tparam T Variable type.
     * @return Inverse hyperbolic sine of x in radians. 
     */
    struct asinh_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::asinh;
        return asinh(x);
      }
    };

    /**
     * Vectorized version of acos().
     * @param x Container.
     * @tparam T Container type.
     * @return Inverse hyperbolic sine of each value in x, in radians. 
     */
    template <typename T>
    inline typename apply_scalar_unary<asinh_fun, T>::return_t
    asinh(const T& x) {
      return apply_scalar_unary<asinh_fun, T>::apply(x);
    }

  }
}

#endif
