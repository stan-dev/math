#ifndef STAN_MATH_FWD_FUN_ERF_HPP
#define STAN_MATH_FWD_FUN_ERF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> erf(const fvar<T>& x) {
  using std::exp;
  return fvar<T>(erf(x.val_), x.d_ * exp(-square(x.val_)) * TWO_OVER_SQRT_PI);
}

}  // namespace math
}  // namespace stan
#endif
