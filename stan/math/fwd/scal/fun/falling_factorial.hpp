#ifndef STAN_MATH_FWD_SCAL_FUN_FALLING_FACTORIAL_HPP
#define STAN_MATH_FWD_SCAL_FUN_FALLING_FACTORIAL_HPP

#include <stan/math/fwd/core.hpp>

#include <stan/math/prim/scal/fun/falling_factorial.hpp>
#include <boost/math/special_functions/digamma.hpp>

namespace stan {
  namespace math {

    template<typename T>
    inline fvar<T>
    falling_factorial(const fvar<T>& x, int n) {
      using boost::math::digamma;

      T falling_fact(falling_factorial(x.val_, n));
      return fvar<T>(falling_fact,
                     falling_fact
                     * (digamma(x.val_ + 1) - digamma(x.val_ - n + 1))
                     * x.d_);
    }
  }
}
#endif
