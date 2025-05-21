#ifndef STAN_MATH_FWD_FUN_INV_SQRT_HPP
#define STAN_MATH_FWD_FUN_INV_SQRT_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/sqrt.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> inv_sqrt(const fvar<T>& x) {
  T sqrt_x(sqrt(x.val_));
  return fvar<T>(1.0 / sqrt_x, -0.5 * x.d_ / (x.val_ * sqrt_x));
}
}  // namespace math
}  // namespace stan
#endif
