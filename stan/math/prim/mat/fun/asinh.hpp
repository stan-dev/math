#ifndef STAN_MATH_PRIM_MAT_FUN_ASINH_HPP
#define STAN_MATH_PRIM_MAT_FUN_ASINH_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of asinh().
     */
    struct asinh_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::asinh;
        return asinh(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<asinh_fun, T>::return_t
    asinh(const T& x) {
      return apply_scalar_unary<asinh_fun, T>::apply(x);
    }

  }
}

#endif
