#ifndef STAN_MATH_PRIM_MAT_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_MAT_FUN_LGAMMA_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of lgamma().
     */
    struct lgamma_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::lgamma;
        return lgamma(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<lgamma_fun, T>::return_t
    lgamma(const T& x) {
      return apply_scalar_unary<lgamma_fun, T>::apply(x);
    }

  }
}

#endif
