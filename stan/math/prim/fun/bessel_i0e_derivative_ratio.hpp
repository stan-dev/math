#ifndef STAN_MATH_PRIM_FUN_BESSEL_I0E_DERIVATIVE_RATIO_HPP
#define STAN_MATH_PRIM_FUN_BESSEL_I0E_DERIVATIVE_RATIO_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/log_modified_bessel_first_kind.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <boost/math/tools/rational.hpp>
#include <utility>

namespace stan {
namespace math {

namespace internal {

/**
 * Implementation of `bessel_i0e_derivative_ratio` without the argument
 * checks.
 *
 * For an arithmetic argument, \f$\rho_0(z) = (I_0(z) - I_1(z)) / I_0(z)\f$
 * is formed for \f$z < 25\f$ from Boost's approximations of \f$I_0\f$ and
 * \f$I_1\f$ (the coefficients of `log_modified_bessel_first_kind`): there
 * the subtraction loses less than two digits, and the continued fraction
 * would need up to 59 terms. Otherwise \f$\rho_{n+1}\f$ is computed as
 * \f$(2n + 3) / ((n + 2 + 2z) t)\f$ with
 * \f$t = 1 - c_{n+2} / (1 - c_{n+3} / (1 - \cdots))\f$ and
 * \f$c_m = (2m + 1) z / ((m + 2z) (m + 1 + 2z))\f$, by the three-term
 * recurrence of the convergents of \f$t\f$, until two successive convergents
 * agree to machine precision. Then \f$\rho_n\f$ follows from
 * \f$\rho_n = (2n + 1) / (n + 1 + 2z - z \rho_{n+1})\f$, a step that damps
 * the rounding error of the convergents. All denominators are divided by
 * \f$1 + z\f$ so that nothing overflows.
 *
 * For an autodiff argument the derivative is
 * \f$\rho_n (\rho_n - \rho_{n+1})\f$, from values of the same function.
 *
 * @tparam T type of the argument
 * @param n order of the derivative in the denominator
 * @param z argument
 * @return \f$\rho_n(z)\f$
 */
template <typename T>
inline return_type_t<T> bessel_i0e_derivative_ratio_impl(int n, const T& z) {
  if constexpr (is_autodiff_v<T>) {
    const auto z_val = value_of(z);
    const auto rho = bessel_i0e_derivative_ratio_impl(n, z_val);
    auto ops_partials = make_partials_propagator(z);
    partials<0>(ops_partials)[0]
        = rho * (rho - bessel_i0e_derivative_ratio_impl(n + 1, z_val));
    return ops_partials.build(rho);
  } else {
    const double x = z;
    if (is_inf(x)) {
      return 0.0;
    }
    if (n == 0 && x < 25.0) {
      using boost::math::tools::evaluate_polynomial;
      if (x < 7.75) {
        const double q = 0.25 * x * x;
        const double i0_minus_1 = q * evaluate_polynomial(bessel_i0_small_p, q);
        const double i1
            = 0.5 * x
              * (1.0
                 + q * (0.5 + q * evaluate_polynomial(bessel_i1_small_p, q)));
        return (1.0 + i0_minus_1 - i1) / (1.0 + i0_minus_1);
      }
      const double t = 1.0 / x;
      const double p0 = evaluate_polynomial(bessel_i0_medium_p, t);
      return (p0 - evaluate_polynomial(bessel_i1_medium_p, t)) / p0;
    }
    const double s = 1.0 / (1.0 + x);
    const double w = x * s;
    double inv_lower = 1.0 / ((n + 2.0) * s + 2.0 * w);
    const double head = (2.0 * n + 3.0) * s * inv_lower;
    double p_prev = 1.0;
    double p = 1.0;
    double q_prev = 0.0;
    double q = 1.0;
    double c_product = 1.0;
    int m = n + 2;
    // successive convergents p / q differ relatively by c_product / (p q_prev)
    do {
      const double inv_upper = 1.0 / ((m + 1.0) * s + 2.0 * w);
      const double c = (2.0 * m + 1.0) * w * s * inv_lower * inv_upper;
      const double p_next = p - c * p_prev;
      const double q_next = q - c * q_prev;
      p_prev = p;
      p = p_next;
      q_prev = q;
      q = q_next;
      inv_lower = inv_upper;
      c_product *= c;
      ++m;
    } while (c_product > EPSILON * p * q_prev);
    const double rho_next = head * q / p;
    return (2.0 * n + 1.0) * s / ((n + 1.0) * s + 2.0 * w - w * rho_next);
  }
}

}  // namespace internal

/**
 * Return \f$\rho_n(z) = -g^{(n+1)}(z) / g^{(n)}(z)\f$, the negated ratio of
 * successive derivatives of the exponentially scaled modified Bessel function
 * of the first kind of order zero, \f$g(z) = e^{-z} I_0(z)\f$. In particular
 * \f$\rho_0(z) = 1 - I_1(z) / I_0(z)\f$, without loss of relative precision
 * for large \f$z\f$, where it is close to \f$1 / (2z)\f$. It is the negative
 * derivative of \f$\log g\f$ (`log_bessel_i0e`).
 *
 * Since \f$g(z) = \frac{1}{\pi} \int_0^\pi e^{-z (1 - \cos\theta)} d\theta\f$,
 * the derivatives of \f$g\f$ alternate in sign and
 * \f$0 < \rho_n < 2\f$. Differentiating \f$z g'' + (2z + 1) g' + g = 0\f$
 * \f$n\f$ times gives the continued fraction
 * \f[
 *   \rho_n = \frac{2n + 1}{n + 1 + 2z - z \rho_{n+1}},
 * \f]
 * whose denominators are at least \f$n + 1\f$, and by the definition of
 * \f$\rho_n\f$,
 * \f[
 *   \rho_n' = \rho_n (\rho_n - \rho_{n+1}).
 * \f]
 * The second relation gives the derivatives of every order from values of
 * the same function, so autodiff types do not differentiate through the
 * iteration that evaluates the first. The relative error is below
 * \f$3 \cdot 10^{-14}\f$ for \f$\rho_0\f$ with \f$z < 25\f$, which is
 * formed from approximations of \f$I_0\f$ and \f$I_1\f$, below
 * \f$3 \cdot 10^{-15}\f$ otherwise for \f$n \le 4\f$, and below
 * \f$10^{-13}\f$ up to the largest order; beyond \f$n\f$ of about 2000 the
 * stopping rule of the iteration would underflow, hence the limit on
 * \f$n\f$.
 *
 * @tparam T type of the argument
 * @param n order of the derivative in the denominator, from 0 to 1000
 * @param z non-negative argument
 * @return \f$\rho_n(z)\f$
 * @throw std::domain_error if `n` is not in [0, 1000] or `z` is NaN or
 * negative
 */
template <typename T, require_stan_scalar_t<T>* = nullptr,
          require_not_complex_t<T>* = nullptr>
inline return_type_t<T> bessel_i0e_derivative_ratio(int n, T&& z) {
  static constexpr const char* function = "bessel_i0e_derivative_ratio";
  check_bounded(function, "n", n, 0, 1000);
  check_nonnegative(function, "z", z);
  return internal::bessel_i0e_derivative_ratio_impl(n, z);
}

/**
 * Return the elementwise `bessel_i0e_derivative_ratio()` of a container.
 *
 * @tparam T type of the container
 * @param n order of the derivative in the denominator, from 0 to 1000
 * @param z container of non-negative arguments
 * @return \f$\rho_n\f$ of each element
 * @throw std::domain_error if `n` is not in [0, 1000] or an element of `z`
 * is NaN or negative
 */
template <typename T, require_container_t<T>* = nullptr>
inline auto bessel_i0e_derivative_ratio(int n, T&& z) {
  return apply_scalar_binary(
      [](auto&& a, auto&& b) {
        return bessel_i0e_derivative_ratio(a, std::forward<decltype(b)>(b));
      },
      n, std::forward<T>(z));
}

}  // namespace math
}  // namespace stan

#endif
