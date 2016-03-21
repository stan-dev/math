#ifndef STAN_MATH_PRIM_MAT_FUN_EXP2_HPP
#define STAN_MATH_PRIM_MAT_FUN_EXP2_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/exp2.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of exp2().
     */
    struct exp2_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::exp2;
        return exp2(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<exp2_fun, T>::return_t
    exp2(const T& x) {
      return apply_scalar_unary<exp2_fun, T>::apply(x);
    }

  }
}

#endif
