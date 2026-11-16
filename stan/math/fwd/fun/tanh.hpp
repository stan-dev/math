#ifndef STAN_MATH_FWD_FUN_TANH_HPP
#define STAN_MATH_FWD_FUN_TANH_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/fwd/fun/value_of_rec.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/tanh.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the hyperbolic tangent of the argument. For `|x| >= 0.5` the
 * derivative `1 - tanh(x)^2` is computed as `4 e / (1 + e)^2` with
 * `e = exp(-2 |x|)`, which neither cancels nor overflows.
 *
 * @tparam T inner type of the forward-mode autodiff variable
 * @param[in] x argument
 * @return hyperbolic tangent of the argument
 */
template <typename T>
inline fvar<T> tanh(const fvar<T>& x) {
  const T tanh_x = tanh(x.val_);
  const double x_val = value_of_rec(x.val_);
  if (fabs(x_val) < 0.5) {
    return fvar<T>(tanh_x, x.d_ * (1 - tanh_x * tanh_x));
  }
  // exp(-2 |x|), with the sign of x fixed in this branch
  const T exp_neg_2abs_x = exp(x_val > 0 ? -2 * x.val_ : 2 * x.val_);
  const T one_plus_exp = 1 + exp_neg_2abs_x;
  return fvar<T>(tanh_x,
                 x.d_ * (4 * exp_neg_2abs_x / (one_plus_exp * one_plus_exp)));
}

/**
 * Return the hyperbolic tangent of the complex argument.
 *
 * @tparam T autodiff value type
 * @param[in] z argument
 * @return hyperbolic tangent of the argument
 */
template <typename T>
inline std::complex<fvar<T>> tanh(const std::complex<fvar<T>>& z) {
  return stan::math::internal::complex_tanh(z);
}

}  // namespace math
}  // namespace stan
#endif
