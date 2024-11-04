#ifndef STAN_MATH_FWD_FUN_TAN_HPP
#define STAN_MATH_FWD_FUN_TAN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/tan.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> tan(const fvar<T>& x) {
  using std::cos;
  using std::tan;
  return fvar<T>(tan(x.val_), x.d_ / (cos(x.val_) * cos(x.val_)));
}


}  // namespace math
}  // namespace stan
#endif
