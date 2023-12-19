#ifndef STAN_MATH_PRIM_FUN_LOG_DIFF_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG_DIFF_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 * The natural logarithm of the difference of the natural exponentiation
 * of x and the natural exponentiation of y
 *
 * This function is only defined for x >= y
 *
   \f[
   \mbox{log\_diff\_exp}(x, y) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < y\\
     \ln(\exp(x)-\exp(y)) & \mbox{if } x \geq y \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_diff\_exp}(x, y)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x \leq y\\
     \frac{\exp(x)}{\exp(x)-\exp(y)} & \mbox{if } x > y \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_diff\_exp}(x, y)}{\partial y} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x \leq y\\
     -\frac{\exp(y)}{\exp(x)-\exp(y)} & \mbox{if } x > y \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T1 type of the first argument
 * @tparam T2 type of the second argument
 * @param x first argument
 * @param y second argument
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> log_diff_exp(const T1 x, const T2 y) {
  if (x <= y) {
    return (x < INFTY && x == y) ? NEGATIVE_INFTY : NOT_A_NUMBER;
  }
  return x + log1m_exp(y - x);
}

/**
 * Enables the vectorized application of the log_diff_exp function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return log_diff_exp function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto log_diff_exp(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return log_diff_exp(c, d); });
}

}  // namespace math
}  // namespace stan

#endif
