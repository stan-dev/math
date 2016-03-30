#ifndef STAN_MATH_PRIM_MAT_FUN_LOG1M_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG1M_EXP_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/log1m_exp.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of log1m_exp().
     */
    struct log1m_exp_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::log1m_exp;
        return log1m_exp(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<log1m_exp_fun, T>::return_t
    log1m_exp(const T& x) {
      return apply_scalar_unary<log1m_exp_fun, T>::apply(x);
    }

  }
}

#endif
