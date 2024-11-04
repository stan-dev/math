#ifndef STAN_MATH_FWD_FUN_SIN_HPP
#define STAN_MATH_FWD_FUN_SIN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/sin.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> sin(const fvar<T>& x) {
  using std::cos;
  using std::sin;
  return fvar<T>(sin(x.val_), x.d_ * cos(x.val_));
}


}  // namespace math
}  // namespace stan
#endif
