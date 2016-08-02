#ifndef STAN_MATH_PRIM_MAT_FUN_LOG1P_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG1P_EXP_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/log1p_exp.hpp>

namespace stan {
  namespace math {

    /**
     * Structure to wrap log1m_exp() so that it can be vectorized.
     * @param x Variable.
     * @tparam T Variable type.
     * @return Natural log of (1 + exp(x)).
     */
    struct log1p_exp_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::log1p_exp;
        return log1p_exp(x);
      }
    };

    /**
     * Vectorized version of log1m_exp().
     * @param x Container.
     * @tparam T Container type.
     * @return Natural log of (1 + exp()) applied to each value in x.
     */
    template <typename T>
    inline typename apply_scalar_unary<log1p_exp_fun, T>::return_t
    log1p_exp(const T& x) {
      return apply_scalar_unary<log1p_exp_fun, T>::apply(x);
    }

  }
}

#endif
