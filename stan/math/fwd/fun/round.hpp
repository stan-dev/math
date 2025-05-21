#ifndef STAN_MATH_FWD_FUN_ROUND_HPP
#define STAN_MATH_FWD_FUN_ROUND_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/is_nan.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/round.hpp>

namespace stan {
namespace math {

/**
 * Return the closest integer to the specified argument, with
 * halfway cases rounded away from zero.
 *
 * The derivative is always zero.
 *
 * @tparam T inner type of the fvar
 * @param x Argument.
 * @return The rounded value of the argument.
 */
template <typename T>
inline fvar<T> round(const fvar<T>& x) {
  return fvar<T>(round(x.val_), is_nan(x.val_) ? NOT_A_NUMBER : 0.0);
}

}  // namespace math
}  // namespace stan
#endif
