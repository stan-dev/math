#ifndef STAN_MATH_PRIM_MAT_FUN_INV_LOGIT_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_LOGIT_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of inv_logit().
     */
    struct inv_logit_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::inv_logit;
        return inv_logit(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<inv_logit_fun, T>::return_t
    inv_logit(const T& x) {
      return apply_scalar_unary<inv_logit_fun, T>::apply(x);
    }

  }
}

#endif
