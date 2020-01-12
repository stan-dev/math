#ifndef STAN_MATH_FWD_FUN_SQRT_HPP
#define STAN_MATH_FWD_FUN_SQRT_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> sqrt(const fvar<T>& x) {
  using std::sqrt;
  return fvar<T>(sqrt(x.val_), 0.5 * x.d_ * inv_sqrt(x.val_));
}
}  // namespace math
}  // namespace stan
#endif
