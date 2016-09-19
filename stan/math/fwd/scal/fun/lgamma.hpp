#ifndef STAN_MATH_FWD_SCAL_FUN_LGAMMA_HPP
#define STAN_MATH_FWD_SCAL_FUN_LGAMMA_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

namespace stan {
  namespace math {

    /**
     * Return the natural logarithm of the gamma function applied to
     * the specified argument.
     *
     * @tparam T Scalar type of autodiff variable.
     * @param x Argument.
     * @return Log gamma function of argument.
     */
    template <typename T>
    inline fvar<T> lgamma(const fvar<T>& x) {
      return fvar<T>(lgamma(x.val_),
                     x.d_ * boost::math::digamma(x.val_));
    }

  }
}
#endif
