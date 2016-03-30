#ifndef STAN_MATH_PRIM_MAT_FUN_LOG_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of log().
     */
    struct log_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::log;
        return log(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<log_fun, T>::return_t
    log(const T& x) {
      return apply_scalar_unary<log_fun, T>::apply(x);
    }

  }
}

#endif
