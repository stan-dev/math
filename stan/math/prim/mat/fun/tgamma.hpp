#ifndef STAN_MATH_PRIM_MAT_FUN_TGAMMA_HPP
#define STAN_MATH_PRIM_MAT_FUN_TGAMMA_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of tgamma().
     */
    struct tgamma_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using boost::math::tgamma;
        return tgamma(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<tgamma_fun, T>::return_t
    tgamma(const T& x) {
      return apply_scalar_unary<tgamma_fun, T>::apply(x);
    }

  }
}

#endif
