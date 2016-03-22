#ifndef STAN_MATH_PRIM_MAT_FUN_ERF_HPP
#define STAN_MATH_PRIM_MAT_FUN_ERF_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <boost/math/special_functions/erf.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of erf().
     */
    struct erf_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using boost::math::erf;
        return erf(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<erf_fun, T>::return_t
    erf(const T& x) {
      return apply_scalar_unary<erf_fun, T>::apply(x);
    }

  }
}

#endif
