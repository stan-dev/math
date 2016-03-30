#ifndef STAN_MATH_PRIM_MAT_FUN_LOG2_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG2_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of log2().
     */
    struct log2_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::log2;
        return log2(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<log2_fun, T>::return_t
    log2(const T& x) {
      return apply_scalar_unary<log2_fun, T>::apply(x);
    }

  }
}

#endif
