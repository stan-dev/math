#ifndef STAN_MATH_FWD_FUN_FLOOR_HPP
#define STAN_MATH_FWD_FUN_FLOOR_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/floor.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> floor(const fvar<T>& x) {
  return fvar<T>(floor(x.val_), 0);
}

}  // namespace math
}  // namespace stan
#endif
