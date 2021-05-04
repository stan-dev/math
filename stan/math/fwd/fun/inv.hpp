#ifndef STAN_MATH_FWD_FUN_INV_HPP
#define STAN_MATH_FWD_FUN_INV_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> inv(const fvar<T>& x) {
  return fvar<T>(1 / x.val_, -x.d_ / square(x.val_));
}
}  // namespace math
}  // namespace stan
#endif
