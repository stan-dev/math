#ifndef STAN_MATH_FWD_FUN_TANH_HPP
#define STAN_MATH_FWD_FUN_TANH_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/prim/fun/tanh.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> tanh(const fvar<T>& x) {
  using std::tanh;
  T u = tanh(x.val_);
  return fvar<T>(u, x.d_ * (1 - u * u));
}


}  // namespace math
}  // namespace stan
#endif
