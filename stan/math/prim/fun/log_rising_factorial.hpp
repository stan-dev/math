#ifndef STAN_MATH_PRIM_FUN_LOG_RISING_FACTORIAL_HPP
#define STAN_MATH_PRIM_FUN_LOG_RISING_FACTORIAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the rising factorial from the
 * first specified argument to the second.
 *
   \f[
   \mbox{log\_rising\_factorial}(x, n) =
   \begin{cases}
     \textrm{error} & \mbox{if } x \leq 0\\
     \ln x^{(n)} & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_rising\_factorial}(x, n)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } x \leq 0\\
     \Psi(x+n) - \Psi(x) & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq
 \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_rising\_factorial}(x, n)}{\partial n} =
   \begin{cases}
     \textrm{error} & \mbox{if } x \leq 0\\
     \Psi(x+n) & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T1 type of first argument x
 * @tparam T2 type of second argument n
 * @param[in] x first argument
 * @param[in] n second argument
 * @return natural logarithm of the rising factorial from x to n
 * @throw std::domain_error if the first argument is not positive
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> log_rising_factorial(const T1& x, const T2& n) {
  if (is_any_nan(x, n)) {
    return NOT_A_NUMBER;
  }
  static const char* function = "log_rising_factorial";
  check_positive(function, "first argument", x);
  return lgamma(x + n) - lgamma(x);
}

/**
 * Enables the vectorized application of the log_rising_factorial
 * function, when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return log_rising_factorial function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto log_rising_factorial(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return log_rising_factorial(c, d);
  });
}

}  // namespace math
}  // namespace stan

#endif
