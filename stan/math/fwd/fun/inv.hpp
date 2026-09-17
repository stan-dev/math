#ifndef STAN_MATH_FWD_FUN_INV_HPP
#define STAN_MATH_FWD_FUN_INV_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/square.hpp>
#include <stan/math/prim/fun/inv.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> inv(const fvar<T>& x) {
  return 1.0 / x;
}
}  // namespace math
}  // namespace stan
#endif
