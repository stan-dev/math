#ifndef STAN_MATH_PRIM_MAT_FUN_DIGAMMA_HPP
#define STAN_MATH_PRIM_MAT_FUN_DIGAMMA_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <boost/math/special_functions/digamma.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of digamma().
     */
    struct digamma_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using boost::math::digamma;
        return digamma(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<digamma_fun, T>::return_t
    digamma(const T& x) {
      return apply_scalar_unary<digamma_fun, T>::apply(x);
    }

  }
}

#endif
