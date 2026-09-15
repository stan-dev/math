#ifndef STAN_MATH_FWD_FUN_ERFCX_HPP
#define STAN_MATH_FWD_FUN_ERFCX_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfcx.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the scaled complementary error function of the argument.
 *
 * The derivative `2 * x * erfcx(x) - 2 / sqrt(pi)` reuses the value, so the
 * tangent is as accurate as the value itself.
 *
 * @tparam T inner type of the fvar
 * @param x argument
 * @return scaled complementary error function of the argument
 */
template <typename T>
inline fvar<T> erfcx(const fvar<T>& x) {
  T v = erfcx(x.val_);
  return fvar<T>(v, x.d_ * (2.0 * x.val_ * v - TWO_OVER_SQRT_PI));
}

}  // namespace math
}  // namespace stan
#endif
