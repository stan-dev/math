#ifndef STAN_MATH_PRIM_MAT_FUN_ACOS_HPP
#define STAN_MATH_PRIM_MAT_FUN_ACOS_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of acos().
     */
    struct acos_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::acos;
        return acos(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<acos_fun, T>::return_t
    acos(const T& x) {
      return apply_scalar_unary<acos_fun, T>::apply(x);
    }

  }
}

#endif
