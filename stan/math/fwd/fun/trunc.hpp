#ifndef STAN_MATH_FWD_FUN_TRUNC_HPP
#define STAN_MATH_FWD_FUN_TRUNC_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/trunc.hpp>

namespace stan {
namespace math {

/**
 * Return the nearest integral value that is not larger in
 * magnitude than the specified argument.
 *
 * @tparam T inner type of the fvar
 * @param[in] x Argument.
 * @return The truncated argument.
 */
template <typename T>
inline fvar<T> trunc(const fvar<T>& x) {
  return fvar<T>(trunc(x.val_), 0);
}

}  // namespace math
}  // namespace stan
#endif
