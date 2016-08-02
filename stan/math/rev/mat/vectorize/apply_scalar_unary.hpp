#ifndef STAN_MATH_REV_MAT_VECTORIZE_APPLY_SCALAR_UNARY_HPP
#define STAN_MATH_REV_MAT_VECTORIZE_APPLY_SCALAR_UNARY_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/rev/core/var.hpp>

namespace stan {

  namespace math {

    /**
     * Template specialization to var for vectorizing a unary scalar
     * function.  This is a base scalar specialization.  It applies
     * the function specified by the template parameter to the
     * argument.  
     *
     * @tparam F Type of function to apply.
     */
    template <typename F>
    struct apply_scalar_unary<F, stan::math::var> {
      /**
       * Function return type, which is <code>var</code>.
       */
      typedef stan::math::var return_t;

      /**
       * Apply the function specified by F to the specified argument.  
       * 
       * @param x Argument variable.
       * @return Function applied to the variable.
       */
      static inline return_t apply(const stan::math::var& x) {
        return F::fun(x);
      }
    };

  }
}
#endif
