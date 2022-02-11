#ifndef STAN_MATH_PRIM_FUN_RISING_FACTORIAL_HPP
#define STAN_MATH_PRIM_FUN_RISING_FACTORIAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <boost/math/special_functions/factorials.hpp>

namespace stan {
namespace math {

/**
 * Return the rising factorial function evaluated
 * at the inputs.
 *
 * @tparam T type of the first argument
 * @param x first argument
 * @param n second argument
 * @return Result of rising factorial function.
 * @throw std::domain_error if x is NaN
 * @throw std::domain_error if n is negative
 *
 \f[
 \mbox{rising\_factorial}(x, n) =
 \begin{cases}
 \textrm{error} & \mbox{if } x \leq 0\\
 x^{(n)} & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq \infty \\[6pt]
 \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
 \end{cases}
 \f]

 \f[
 \frac{\partial\, \mbox{rising\_factorial}(x, n)}{\partial x} =
 \begin{cases}
 \textrm{error} & \mbox{if } x \leq 0\\
 \frac{\partial\, x^{(n)}}{\partial x} & \mbox{if } x > 0 \textrm{ and } -\infty
 \leq n \leq \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n =
 \textrm{NaN} \end{cases} \f]

 \f[
 \frac{\partial\, \mbox{rising\_factorial}(x, n)}{\partial n} =
 \begin{cases}
 \textrm{error} & \mbox{if } x \leq 0\\
 \frac{\partial\, x^{(n)}}{\partial n} & \mbox{if } x > 0 \textrm{ and } -\infty
 \leq n \leq \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n =
 \textrm{NaN} \end{cases} \f]

 \f[
 x^{(n)}=\frac{\Gamma(x+n)}{\Gamma(x)}
 \f]

 \f[
 \frac{\partial \, x^{(n)}}{\partial x} = x^{(n)}(\Psi(x+n)-\Psi(x))
 \f]

 \f[
 \frac{\partial \, x^{(n)}}{\partial n} = (x)_n\Psi(x+n)
 \f]
 *
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline return_type_t<T> rising_factorial(const T& x, int n) {
  static const char* function = "rising_factorial";
  check_not_nan(function, "first argument", x);
  check_nonnegative(function, "second argument", n);
  return boost::math::rising_factorial(x, n, boost_policy_t<>());
}

/**
 * Enables the vectorised application of the rising_factorial
 * function, when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return rising_factorial function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto rising_factorial(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return rising_factorial(c, d);
  });
}

}  // namespace math
}  // namespace stan
#endif
