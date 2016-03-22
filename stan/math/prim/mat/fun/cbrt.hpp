#ifndef STAN_MATH_PRIM_MAT_FUN_CBRT_HPP
#define STAN_MATH_PRIM_MAT_FUN_CBRT_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of cbrt().
     */
    struct cbrt_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using boost::math::cbrt;
        return cbrt(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<cbrt_fun, T>::return_t
    cbrt(const T& x) {
      return apply_scalar_unary<cbrt_fun, T>::apply(x);
    }

  }
}

#endif
