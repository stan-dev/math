#ifndef STAN_MATH_PRIM_MAT_FUN_INV_SQUARE_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_SQUARE_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/inv_square.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of inv_square().
     */
    struct inv_square_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::inv_square;
        return inv_square(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<inv_square_fun, T>::return_t
    inv_square(const T& x) {
      return apply_scalar_unary<inv_square_fun, T>::apply(x);
    }

  }
}

#endif
