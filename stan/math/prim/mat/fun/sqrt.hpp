#ifndef STAN_MATH_PRIM_MAT_FUN_SQRT_HPP
#define STAN_MATH_PRIM_MAT_FUN_SQRT_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of sqrt().
     */
    struct sqrt_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::sqrt;
        return sqrt(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<sqrt_fun, T>::return_t
    sqrt(const T& x) {
      return apply_scalar_unary<sqrt_fun, T>::apply(x);
    }

  }
}

#endif
