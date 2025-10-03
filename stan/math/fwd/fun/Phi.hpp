#ifndef STAN_MATH_FWD_FUN_PHI_HPP
#define STAN_MATH_FWD_FUN_PHI_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/pow.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Phi.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> Phi(const fvar<T>& x) {
  T xv = x.val_;
  return fvar<T>(Phi(xv), x.d_ * exp(xv * xv / -2.0) * INV_SQRT_TWO_PI);
}

}  // namespace math
}  // namespace stan
#endif
