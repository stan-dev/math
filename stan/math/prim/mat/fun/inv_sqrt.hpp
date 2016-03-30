#ifndef STAN_MATH_PRIM_MAT_FUN_INV_SQRT_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_SQRT_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/inv_sqrt.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of inv_sqrt().
     */
    struct inv_sqrt_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::inv_sqrt;
        return inv_sqrt(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<inv_sqrt_fun, T>::return_t
    inv_sqrt(const T& x) {
      return apply_scalar_unary<inv_sqrt_fun, T>::apply(x);
    }

  }
}

#endif
