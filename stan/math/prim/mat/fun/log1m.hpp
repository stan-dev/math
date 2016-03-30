#ifndef STAN_MATH_PRIM_MAT_FUN_LOG1M_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG1M_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of log1m().
     */
    struct log1m_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::log1m;
        return log1m(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<log1m_fun, T>::return_t
    log1m(const T& x) {
      return apply_scalar_unary<log1m_fun, T>::apply(x);
    }

  }
}

#endif
