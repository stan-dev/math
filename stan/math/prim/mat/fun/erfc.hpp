#ifndef STAN_MATH_PRIM_MAT_FUN_ERFC_HPP
#define STAN_MATH_PRIM_MAT_FUN_ERFC_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <boost/math/special_functions/erf.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of erfc().
     */
    struct erfc_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using boost::math::erfc;
        return erfc(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<erfc_fun, T>::return_t
    erfc(const T& x) {
      return apply_scalar_unary<erfc_fun, T>::apply(x);
    }

  }
}

#endif
