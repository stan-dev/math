#ifndef STAN_MATH_FWD_FUN_EXPM1_HPP
#define STAN_MATH_FWD_FUN_EXPM1_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> expm1(const fvar<T>& x) {
  return fvar<T>(expm1(x.val_), x.d_ * exp(x.val_));
}

}  // namespace math
}  // namespace stan
#endif
