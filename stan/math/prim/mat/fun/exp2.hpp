#ifndef STAN_MATH_PRIM_MAT_FUN_EXP2_HPP
#define STAN_MATH_PRIM_MAT_FUN_EXP2_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/exp2.hpp>

namespace stan {
  namespace math {

    /**
     * Structure to wrap exp2() so that it can be vectorized.
     * @param x Variable.
     * @tparam T Variable type.
     * @return Base-2 exponential of x. 
     */
    struct exp2_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::exp2;
        return exp2(x);
      }
    };

    /**
     * Vectorized version of exp2().
     * @param x Container.
     * @tparam T Container type.
     * @return Base-2 exponential of each value in x. 
     */
    template <typename T>
    inline typename apply_scalar_unary<exp2_fun, T>::return_t
    exp2(const T& x) {
      return apply_scalar_unary<exp2_fun, T>::apply(x);
    }

  }
}

#endif
