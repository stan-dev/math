#ifndef STAN_MATH_PRIM_FUN_LOG_FALLING_FACTORIAL_HPP
#define STAN_MATH_PRIM_FUN_LOG_FALLING_FACTORIAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 *
 * Return the natural log of the falling factorial of the
 * specified arguments.
 *
   \f[
   \mbox{log\_falling\_factorial}(x, n) =
   \begin{cases}
     \textrm{error} & \mbox{if } x \leq 0\\
     \ln (x)_n & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_falling\_factorial}(x, n)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } x \leq 0\\
     \Psi(x) & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_falling\_factorial}(x, n)}{\partial n} =
   \begin{cases}
     \textrm{error} & \mbox{if } x \leq 0\\
     -\Psi(n) & \mbox{if } x > 0 \textrm{ and } -\infty \leq n \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } n = \textrm{NaN}
   \end{cases}
   \f]

 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param[in] x First argument
 * @param[in] n Second argument
 * @return log of falling factorial of arguments
 * @throw std::domain_error if the first argument is not
 * positive
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> log_falling_factorial(const T1 x, const T2 n) {
  if (is_any_nan(x, n)) {
    return NOT_A_NUMBER;
  }
  static const char* function = "log_falling_factorial";
  check_positive(function, "first argument", x);
  return lgamma(x + 1) - lgamma(x - n + 1);
}

/**
 * Enables the vectorized application of the log_falling_factorial function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return log_falling_factorial function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto log_falling_factorial(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return log_falling_factorial(c, d);
  });
}

}  // namespace math
}  // namespace stan
#endif
