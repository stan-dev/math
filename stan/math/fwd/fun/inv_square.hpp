#ifndef STAN_MATH_FWD_FUN_INV_SQUARE_HPP
#define STAN_MATH_FWD_FUN_INV_SQUARE_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/square.hpp>
#include <stan/math/prim/fun/inv_square.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> inv_square(const fvar<T>& x) {
  return fvar<T>(inv_square(x.val_), -2 * ((x.d_ / x.val_) / x.val_) / x.val_);
}
}  // namespace math
}  // namespace stan
#endif
