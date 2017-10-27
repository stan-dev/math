#ifndef STAN_MATH_FWD_SCAL_FUN_EXP_HPP
#define STAN_MATH_FWD_SCAL_FUN_EXP_HPP

#include <cmath>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {
template <typename T>
inline fvar<T> exp(const fvar<T>& x) {
  using std::exp;
  return fvar<T>(exp(x.val_), x.d_ * exp(x.val_));
}

}  // namespace math
}  // namespace stan
#endif
