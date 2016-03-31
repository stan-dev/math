#ifndef STAN_MATH_PRIM_MAT_FUN_LOG2_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG2_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Structure to wrap log2() so that it can be vectorized.
     * @param x Variable.
     * @tparam T Variable type.
     * @return Log base-2 of x.
     */
    struct log2_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::log2;
        return log2(x);
      }
    };

    /**
     * Vectorized version of log2().
     * @param x Container.
     * @tparam T Container type.
     * @return Log base-2 of each value in x.
     */
    template <typename T>
    inline typename apply_scalar_unary<log2_fun, T>::return_t
    log2(const T& x) {
      return apply_scalar_unary<log2_fun, T>::apply(x);
    }

  }
}

#endif
