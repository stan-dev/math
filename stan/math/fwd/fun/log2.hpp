#ifndef STAN_MATH_FWD_FUN_LOG2_HPP
#define STAN_MATH_FWD_FUN_LOG2_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log2.hpp>

namespace stan {
namespace math {

/**
 * Return the base two logarithm of the specified argument.
 *
 * @tparam T inner type of the fvar
 * @param x argument
 * @return base two logarithm of argument
 */
template <typename T>
inline fvar<T> log2(const fvar<T>& x) {
  if (x.val_ < 0.0) {
    return fvar<T>(NOT_A_NUMBER, NOT_A_NUMBER);
  } else {
    return fvar<T>(log2(x.val_), x.d_ / (x.val_ * LOG_TWO));
  }
}

}  // namespace math
}  // namespace stan
#endif
