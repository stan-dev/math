#ifndef STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG1P_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
// #include <math.h>

namespace stan {
  namespace math {

    /**
     * Structure to wrap log1p() so it can be vectorized.
     */
    struct log1p_fun {
      /**
       * Return the natural log of one plus the specified argument.
       *
       * @param x Argument.
       * @return Natural log of one plus the argument.
       * @tparam T Argument type.
       */
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::log1p;
        return log1p(x);
      }
    };

    /**
     * Return the elementwise application of <code>log1p()</code> to
     * specified argument container.  The return type promotes the
     * underlying scalar argument type to double if it is an integer,
     * and otherwise is the argument type.
     *
     * @tparam T Container type.
     * @param x Container.
     * @return Elementwise log1p of members of container.
     */
    template <typename T>
    inline typename 
    boost::disable_if<boost::is_arithmetic<T>,
                      typename apply_scalar_unary<log1p_fun, T>::return_t>::type
    log1p(const T& x) {
      return apply_scalar_unary<log1p_fun, T>::apply(x);
    }

  }
}

#endif
