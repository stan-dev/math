#ifndef STAN_MATH_FWD_FUN_LOG_COSH_HPP
#define STAN_MATH_FWD_FUN_LOG_COSH_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/tanh.hpp>
#include <stan/math/prim/fun/log_cosh.hpp>
#include <stan/math/prim/fun/tanh.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the hyperbolic cosine of the argument.
 *
 * @tparam T inner type of the forward-mode autodiff variable
 * @param x argument
 * @return `log(cosh(x))`, with derivative `tanh(x)`
 */
template <typename T>
inline fvar<T> log_cosh(const fvar<T>& x) {
  return fvar<T>(log_cosh(x.val_), x.d_ * tanh(x.val_));
}

}  // namespace math
}  // namespace stan
#endif
