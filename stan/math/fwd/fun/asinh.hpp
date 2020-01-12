#ifndef STAN_MATH_FWD_FUN_ASINH_HPP
#define STAN_MATH_FWD_FUN_ASINH_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/asinh.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> asinh(const fvar<T>& x) {
  using std::sqrt;
  return fvar<T>(asinh(x.val_), x.d_ / sqrt(square(x.val_) + 1));
}

}  // namespace math
}  // namespace stan
#endif
