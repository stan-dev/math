#ifndef STAN_MATH_PRIM_MAT_FUN_LOG_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Structure to wrap log() so that it can be vectorized.
     * @param x Variable.
     * @tparam T Variable type.
     * @return Natural log of x.
     */
    struct log_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::log;
        return log(x);
      }
    };

    /**
     * Vectorized version of log().
     * @param x Container.
     * @tparam T Container type.
     * @return Natural log applied to each value in x.
     */
    template <typename T>
    inline typename apply_scalar_unary<log_fun, T>::return_t
    log(const T& x) {
      return apply_scalar_unary<log_fun, T>::apply(x);
    }

  }
}

#endif
