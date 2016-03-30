#ifndef STAN_MATH_PRIM_MAT_FUN_TAN_HPP
#define STAN_MATH_PRIM_MAT_FUN_TAN_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of tan().
     */
    struct tan_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::tan;
        return tan(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<tan_fun, T>::return_t
    tan(const T& x) {
      return apply_scalar_unary<tan_fun, T>::apply(x);
    }

  }
}

#endif
