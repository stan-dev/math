#ifndef STAN_MATH_FWD_FUN_CBRT_HPP
#define STAN_MATH_FWD_FUN_CBRT_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/cbrt.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {

/**
 * Return cube root of specified argument.
 *
 * @tparam T inner type of the fvar
 * @param x Argument.
 * @return Cube root of argument.
 */
template <typename T>
inline fvar<T> cbrt(const fvar<T>& x) {
  return fvar<T>(cbrt(x.val_), x.d_ / (3 * square(cbrt(x.val_))));
}

}  // namespace math
}  // namespace stan
#endif
