#ifndef STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of log1p().
     */
    struct log1p_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::log1p;
        return log1p(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<log1p_fun, T>::return_t
    log1p(const T& x) {
      return apply_scalar_unary<log1p_fun, T>::apply(x);
    }

  }
}

#endif
