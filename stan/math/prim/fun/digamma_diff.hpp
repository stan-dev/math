#ifndef STAN_MATH_PRIM_FUN_DIGAMMA_DIFF_HPP
#define STAN_MATH_PRIM_FUN_DIGAMMA_DIFF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {
/**
 * Gradient code that needs digamma(x) - digamma(x + d) can use the plain
 * difference for x below this value: its absolute error is a few eps times
 * |digamma|, and x times that error (the error of a gradient in log(x))
 * stays at the rounding level. From this value on it must use digamma_diff.
 * The derivatives of lbeta use this rule.
 */
constexpr double digamma_diff_min_x = 10.0;
}  // namespace internal

/**
 * Return the difference of the digamma function at two arguments that
 * differ by a nonnegative offset,
 *
   \f[
   \mbox{digamma\_diff}(x, d) = \Psi(x + d) - \Psi(x).
   \f]
 *
 * The plain difference of two `digamma` calls loses all accuracy when `x`
 * is large and `d` is not: both values are close to \f$\log x\f$ and their
 * difference is close to \f$d / x\f$, so the absolute error is about
 * \f$\epsilon \log x\f$. This function has a relative error of a few ulp
 * for all \f$x > 0\f$ and \f$d \ge 0\f$. Densities with shape or
 * dispersion parameters use it for gradients such as
 * \f$\Psi(\alpha + n) - \Psi(\alpha)\f$.
 *
 * Method. Two cases have a cheaper form of the same accuracy.
 *
 * - If `d` is not an autodiff type and is an integer from 0 to 8 (a count),
 *   \f$\Psi(x + d) - \Psi(x) = \sum_{j=0}^{d-1} 1 / (x + j)\f$. Every term
 *   is positive, the terms are added from the smallest, and for \f$d = 1\f$
 *   the result is exactly \f$1/x\f$.
 * - If \f$x < 10\f$ and \f$d \ge 10\f$, the plain difference
 *   \f$\Psi(x + d) - \Psi(x)\f$ is used. The result is at least
 *   \f$\Psi(20) - \Psi(10) \approx 0.72\f$ and \f$\Psi(x) \le \Psi(10)
 *   \approx 2.25\f$, so the rounding errors of the two values grow by a
 *   factor of at most about 7 (nothing cancels where \f$\Psi(x) < 0\f$).
 *
 * Otherwise, for \f$x < 10\f$ the recurrence
 * \f$\Psi(y + 1) = \Psi(y) + 1/y\f$ shifts both arguments up by the same
 * integer \f$J\f$:
 *
   \f[
   \Psi(x + d) - \Psi(x) = \sum_{j=0}^{J-1} \frac{d}{(x + j)(x + j + d)}
     + \Psi(y + d) - \Psi(y), \quad y = x + J \ge 10.
   \f]
 *
 * Every term of the sum is positive. A term with \f$d \ge x + j\f$ is
 * split into \f$1/(x + j) - 1/(x + j + d)\f$, which does not cancel, and
 * the parts \f$1/(x + j)\f$ are added last, so that for small \f$x\f$ the
 * dominant \f$1/x\f$ keeps its correct rounding; for \f$d = 1\f$ the result
 * is exactly \f$1/x\f$. For \f$y \ge 10\f$ the asymptotic
 * expansion \f$\Psi(y) = \log y - 1/(2y) - \sum_i B_{2i} / (2i\,y^{2i})\f$
 * gives
 *
   \f[
   \Psi(y + d) - \Psi(y) = \log\left(1 + \frac{d}{y}\right)
     + \frac{d}{2y(y + d)}
     + \sum_{i=1}^{8} \frac{B_{2i}}{2i} \left(u^i - v^i\right),
   \f]
 *
 * with \f$u = y^{-2}\f$ and \f$v = (y + d)^{-2}\f$. The differences
 * \f$u^i - v^i = (u - v) h_i\f$, \f$h_1 = 1\f$,
 * \f$h_{i+1} = u h_i + v^i\f$, and
 * \f$u - v = u \frac{d}{y + d} \left(1 + \frac{y}{y + d}\right)\f$ are
 * formed without cancellation. The truncation error is below
 * \f$|B_{18}| y^{-18} \approx 6 \times 10^{-17}\f$ relative.
 *
 * @tparam T1 type of the first argument
 * @tparam T2 type of the second argument
 * @param x first argument, positive
 * @param d offset, nonnegative
 * @return \f$\Psi(x + d) - \Psi(x)\f$
 * @throw std::domain_error if `x` is not positive or `d` is negative
 */
template <typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> digamma_diff(const T1& x, const T2& d) {
  using T_ret = return_type_t<T1, T2>;
  static constexpr const char* function = "digamma_diff";
  if (is_any_nan(x, d)) {
    return NOT_A_NUMBER;
  }
  check_positive(function, "first argument", x);
  check_nonnegative(function, "second argument", d);
  // d = 0 needs no special case: every term below is then exactly 0, and
  // the derivative with respect to d stays correct
  if (is_inf(d)) {
    return INFTY;
  }
  if (is_inf(x)) {
    return T_ret(0.0);
  }

  // A count d from 0 to 8: sum of the positive terms 1 / (x + j), from the
  // smallest (exactly 1 / x for d = 1)
  if constexpr (std::is_arithmetic<T2>::value) {
    if (d <= 8 && d == std::floor(d)) {
      T_ret sum(0.0);
      for (int j = static_cast<int>(d) - 1; j >= 0; --j) {
        sum += inv(x + j);
      }
      return sum;
    }
  }
  // x < 10 and d >= 10: the plain difference does not cancel much
  if (value_of_rec(x) < 10.0 && value_of_rec(d) >= 10.0) {
    return digamma(x + d) - digamma(x);
  }

  // B_{2i} / (2i), i = 1..8
  static constexpr double coeffs[]
      = {1.0 / 12.0,  -1.0 / 120.0,     1.0 / 252.0, -1.0 / 240.0,
         1.0 / 132.0, -691.0 / 32760.0, 1.0 / 12.0,  -3617.0 / 8160.0};
  static constexpr int n_coeffs = 8;
  static constexpr double shift_to = 10.0;

  // Each shift term is d / (y (y + d)) = 1 / y - 1 / (y + d). For d < y it
  // is formed as (d / (y + d)) / y, which neither cancels nor overflows. For
  // d >= y the two parts do not cancel either, and the parts 1 / y are
  // summed separately and added last: for small x the result is dominated
  // by 1 / x, which then keeps its correct rounding (for d = 1 the result is
  // exactly 1 / x, as for digamma(x + 1) - digamma(x)).
  T_ret inv_sum(0.0);
  T_ret shift_sum(0.0);
  T_ret y = x;
  while (value_of_rec(y) < shift_to) {
    if (value_of_rec(d) >= value_of_rec(y)) {
      inv_sum += inv(y);
      shift_sum -= inv(y + d);
    } else {
      shift_sum += (d / (y + d)) / y;
    }
    y += 1.0;
  }

  const T_ret y_plus_d = y + d;
  const T_ret d_frac = d / y_plus_d;
  // square(inv(.)) underflows to 0 where inv(square(.)) would overflow
  const T_ret u = square(inv(y));
  const T_ret v = square(inv(y_plus_d));
  const T_ret u_minus_v = u * d_frac * (1.0 + y / y_plus_d);
  T_ret h(1.0);
  T_ret v_pow = v;
  T_ret series = coeffs[0];
  for (int i = 1; i < n_coeffs; ++i) {
    h = u * h + v_pow;
    v_pow *= v;
    series += coeffs[i] * h;
  }
  return inv_sum
         + (shift_sum + log1p(d / y) + 0.5 * d_frac / y + u_minus_v * series);
}

/**
 * Enables the vectorized application of the digamma_diff function, when the
 * first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return digamma_diff function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto digamma_diff(T1&& a, T2&& b) {
  return apply_scalar_binary(
      [](auto&& c, auto&& d) {
        return digamma_diff(std::forward<decltype(c)>(c),
                            std::forward<decltype(d)>(d));
      },
      std::forward<T1>(a), std::forward<T2>(b));
}

}  // namespace math
}  // namespace stan

#endif
