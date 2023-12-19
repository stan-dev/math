#ifndef STAN_MATH_PRIM_FUN_LOG_INV_LOGIT_DIFF_HPP
#define STAN_MATH_PRIM_FUN_LOG_INV_LOGIT_DIFF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 * Returns the natural logarithm of the difference of the
 * inverse logits of the specified arguments.
 *
   \f[
     \mathrm{log\_inv\_logit\_diff}(x,y) =
      \ln\left(\frac{1}{1+\exp(-x)}-\frac{1}{1+\exp(-y)}\right)
   \f]

   \f[
    \frac{\partial }{\partial x} = -\frac{e^x}{e^y-e^x}-\frac{e^x}{e^x+1}
   \f]

   \f[
    \frac{\partial }{\partial x} = -\frac{e^y}{e^x-e^y}-\frac{e^y}{e^y+1}
   \f]
 *
 * @tparam T1 type of x argument
 * @tparam T2 type of y argument
 * @param x first argument
 * @param y second argument
 * @return Result of log difference of inverse logits of arguments.
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> log_inv_logit_diff(const T1& x, const T2& y) {
  if (is_inf(x) && x >= 0) {
    return -log1p_exp(y);
  }
  return x - log1p_exp(x) + log1m_exp(y - x) - log1p_exp(y);
}

/**
 * Enables the vectorized application of the log_inv_logit_diff function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return log_inv_logit_diff function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto log_inv_logit_diff(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [&](const auto& c, const auto& d) {
    return log_inv_logit_diff(c, d);
  });
}

}  // namespace math
}  // namespace stan

#endif
