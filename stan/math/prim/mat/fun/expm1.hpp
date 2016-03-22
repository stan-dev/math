#ifndef STAN_MATH_PRIM_MAT_FUN_EXPM1_HPP
#define STAN_MATH_PRIM_MAT_FUN_EXPM1_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of expm1().
     */
    struct expm1_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using boost::math::expm1;
        return expm1(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<expm1_fun, T>::return_t
    expm1(const T& x) {
      return apply_scalar_unary<expm1_fun, T>::apply(x);
    }

  }
}

#endif
