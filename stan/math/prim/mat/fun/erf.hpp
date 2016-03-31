#ifndef STAN_MATH_PRIM_MAT_FUN_ERF_HPP
#define STAN_MATH_PRIM_MAT_FUN_ERF_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <boost/math/special_functions/erf.hpp>

namespace stan {
  namespace math {

    /**
     * Structure to wrap erf() so it can be vectorized.
     * @param x Variable.
     * @tparam T Variable type.
     * @return Error function of x. 
     */
    struct erf_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using boost::math::erf;
        return erf(x);
      }
    };

    /**
     * Vectorized version of erf().
     * @param x Container.
     * @tparam T Container type.
     * @return Error function applied to each value in x. 
     */
    template <typename T>
    inline typename apply_scalar_unary<erf_fun, T>::return_t
    erf(const T& x) {
      return apply_scalar_unary<erf_fun, T>::apply(x);
    }

  }
}

#endif
