#ifndef STAN_MATH_FWD_FUN_OWENS_T_HPP
#define STAN_MATH_FWD_FUN_OWENS_T_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/owens_t.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return Owen's T function applied to the specified
 * arguments.
 *
 * @tparam T inner type of the fvar
 * @param x1 First argument.
 * @param x2 Second argument.
 * @return Owen's T function applied to the specified arguments.
 */
template <typename T>
inline fvar<T> owens_t(const fvar<T>& x1, const fvar<T>& x2) {
  using std::exp;

  T neg_x1_sq_div_2 = -square(x1.val_) * 0.5;
  T one_p_x2_sq = 1.0 + square(x2.val_);
  return fvar<T>(owens_t(x1.val_, x2.val_),
                 -x1.d_
                         * (erf(x2.val_ * x1.val_ * INV_SQRT_TWO)
                            * exp(neg_x1_sq_div_2) * INV_SQRT_TWO_PI * 0.5)
                     + x2.d_ * exp(neg_x1_sq_div_2 * one_p_x2_sq)
                           / (one_p_x2_sq * TWO_PI));
}

/**
 * Return Owen's T function applied to the specified arguments.
 *
 * @tparam T inner type of the fvar
 * @param x1 First argument.
 * @param x2 Second argument.
 * @return Owen's T function applied to the specified arguments.
 */
template <typename T>
inline fvar<T> owens_t(double x1, const fvar<T>& x2) {
  using std::exp;

  T neg_x1_sq_div_2 = -square(x1) * 0.5;
  T one_p_x2_sq = 1.0 + square(x2.val_);
  return fvar<T>(
      owens_t(x1, x2.val_),
      x2.d_ * exp(neg_x1_sq_div_2 * one_p_x2_sq) / (one_p_x2_sq * TWO_PI));
}

/**
 * Return Owen's T function applied to the specified arguments.
 *
 * @tparam T inner type of the fvar
 * @param x1 First argument.
 * @param x2 Second argument.
 * @return Owen's T function applied to the specified arguments.
 */
template <typename T>
inline fvar<T> owens_t(const fvar<T>& x1, double x2) {
  using std::exp;

  T neg_x1_sq_div_2 = -square(x1.val_) * 0.5;
  return fvar<T>(owens_t(x1.val_, x2),
                 -x1.d_
                     * (erf(x2 * x1.val_ * INV_SQRT_TWO) * exp(neg_x1_sq_div_2)
                        * INV_SQRT_TWO_PI * 0.5));
}

}  // namespace math
}  // namespace stan
#endif
