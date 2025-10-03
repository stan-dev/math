#ifndef STAN_MATH_FWD_FUN_TGAMMA_HPP
#define STAN_MATH_FWD_FUN_TGAMMA_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/digamma.hpp>
#include <stan/math/prim/fun/tgamma.hpp>

namespace stan {
namespace math {

/**
 * Return the result of applying the gamma function to the
 * specified argument.
 *
 * @tparam T inner type of the fvar
 * @param x Argument.
 * @return Gamma function applied to argument.
 */
template <typename T>
inline fvar<T> tgamma(const fvar<T>& x) {
  T u = tgamma(x.val_);
  return fvar<T>(u, x.d_ * u * digamma(x.val_));
}

}  // namespace math
}  // namespace stan
#endif
