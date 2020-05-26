#ifndef STAN_MATH_FWD_FUN_LAMBERT_W_HPP
#define STAN_MATH_FWD_FUN_LAMBERT_W_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/fun/lambert_w.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> lambert_w0(const fvar<T>& x) {
  return fvar<T>(lambert_w0(x.val_),  (x.d_ / (x.val_ + exp(lambert_w0(x.val_)))));
}


}  // namespace math
}  // namespace stan
#endif
