#ifndef STAN_MATH_FWD_SCAL_FUN_INV_PHI_HPP
#define STAN_MATH_FWD_SCAL_FUN_INV_PHI_HPP

#include <stan/math/fwd/core.hpp>

#include <stan/math/prim/scal/fun/inv_Phi.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>

namespace stan {

  namespace math {

    template <typename T>
    inline fvar<T> inv_Phi(const fvar<T>& p) {
      using stan::math::inv_Phi;
      using std::exp;
      using std::log;
      T xv = inv_Phi(p.val_);
      return fvar<T>(xv,
                     p.d_ / exp(-0.5 * xv * xv) * SQRT_2_TIMES_SQRT_PI);
    }
  }
}
#endif
