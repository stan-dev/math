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
 * The derivative comes from `internal::erfcx_derivative`. Below 4 it is
 * `2 * x * erfcx(x) - 2 / sqrt(pi)`, which reuses the value. At 4 and above
 * that difference cancels, so the tail rational supplies the derivative
 * directly.
 *
 * @tparam T inner type of the fvar
 * @param x argument
 * @return scaled complementary error function of the argument
 */
template <typename T>
inline fvar<T> erfcx(const fvar<T>& x) {
  T v = erfcx(x.val_);
  return fvar<T>(v, x.d_ * internal::erfcx_derivative(x.val_, v));
}

}  // namespace math
}  // namespace stan
#endif
