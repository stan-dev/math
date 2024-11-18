#ifndef STAN_MATH_FWD_FUN_DIGAMMA_HPP
#define STAN_MATH_FWD_FUN_DIGAMMA_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/trigamma.hpp>

#include <stan/math/prim/fun/digamma.hpp>

namespace stan {
namespace math {

/**
 * Return the derivative of the log gamma function at the
 * specified argument.
 *
 * @tparam T inner type of the fvar
 * @param[in] x argument
 * @return derivative of the log gamma function at the specified
 * argument
 */
template <typename T>
inline fvar<T> digamma(const fvar<T>& x) {
  return fvar<T>(digamma(x.val_), x.d_ * trigamma(x.val_));
}

}  // namespace math
}  // namespace stan
#endif
