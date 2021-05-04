#ifndef STAN_MATH_FWD_FUN_SQUARE_HPP
#define STAN_MATH_FWD_FUN_SQUARE_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> square(const fvar<T>& x) {
  return fvar<T>(square(x.val_), x.d_ * 2 * x.val_);
}
}  // namespace math
}  // namespace stan
#endif
