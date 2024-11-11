#ifndef STAN_MATH_FWD_FUN_CEIL_HPP
#define STAN_MATH_FWD_FUN_CEIL_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/ceil.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> ceil(const fvar<T>& x) {
  return fvar<T>(ceil(x.val_), 0);
}

}  // namespace math
}  // namespace stan
#endif
