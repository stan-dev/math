#ifndef STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>

namespace stan {
  namespace math {

    /**
     * Structure to wrap log1p() so it can be vectorized.
     * @param x Variable.
     * @tparam T Variable type.
     * @return Natural log of (1 + x).
     */
    struct log1p_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::log1p;
        return log1p(x);
      }
    };

    /**
     * Vectorized version of log1m_exp().
     * @param x Container.
     * @tparam T Container type.
     * @return Natural log of one plus each value in x.
     */
    template <typename T>
    inline typename apply_scalar_unary<log1p_fun, T>::return_t
    log1p(const T& x) {
      return apply_scalar_unary<log1p_fun, T>::apply(x);
    }

  }
}

#endif
