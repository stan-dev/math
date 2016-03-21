#ifndef STAN_MATH_PRIM_MAT_FUN_SQUARE_HPP
#define STAN_MATH_PRIM_MAT_FUN_SQUARE_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/square.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of square().
     */
    struct square_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::square;
        return square(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<square_fun, T>::return_t
    square(const T& x) {
      return apply_scalar_unary<square_fun, T>::apply(x);
    }

  }
}

#endif
