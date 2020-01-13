#ifndef STAN_MATH_FWD_FUN_INV_LOGIT_HPP
#define STAN_MATH_FWD_FUN_INV_LOGIT_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the inverse logit function applied to the argument.
 *
 * @tparam T inner type of the fvar
 * @param x argument
 * @return inverse logit of argument
 */
template <typename T>
inline fvar<T> inv_logit(const fvar<T>& x) {
  using std::exp;
  using std::pow;
  return fvar<T>(inv_logit(x.val_),
                 x.d_ * inv_logit(x.val_) * (1 - inv_logit(x.val_)));
}

}  // namespace math
}  // namespace stan
#endif
