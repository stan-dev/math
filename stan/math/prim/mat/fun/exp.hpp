#ifndef STAN_MATH_PRIM_MAT_FUN_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_EXP_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Structure to wrap exp() so that it can be vectorized.
     * @param x Variable.
     * @tparam T Variable type.
     * @return Natural exponential of x. 
     */
    struct exp_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::exp;
        return exp(x);
      }
    };

    /**
     * Vectorized version of exp().
     * @param x Container.
     * @tparam T Container type.
     * @return Natural exponential applied to each value in x. 
     */
    template <typename T>
    inline typename apply_scalar_unary<exp_fun, T>::return_t
    exp(const T& x) {
      return apply_scalar_unary<exp_fun, T>::apply(x);
    }

  }
}

#endif
