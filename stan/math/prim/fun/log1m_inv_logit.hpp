#ifndef STAN_MATH_PRIM_FUN_LOG1M_INV_LOGIT_HPP
#define STAN_MATH_PRIM_FUN_LOG1M_INV_LOGIT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the natural logarithm of 1 minus the inverse logit
 * of the specified argument.
 *
   \f[
   \mbox{log1m\_inv\_logit}(x) =
   \begin{cases}
     -\ln(\exp(x)+1) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log1m\_inv\_logit}(x)}{\partial x} =
   \begin{cases}
     -\frac{\exp(x)}{\exp(x)+1} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param u argument
 * @return log of one minus the inverse logit of the argument
 */
inline double log1m_inv_logit(double u) {
  using std::exp;
  if (u > 0.0) {
    return -u - log1p(exp(-u));  // prevent underflow
  }
  return -log1p(exp(u));
}

/**
 * Return the natural logarithm of one minus the inverse logit of
 * the specified argument.
 *
 * @param u argument
 * @return log of one minus the inverse logit of the argument
 */
inline double log1m_inv_logit(int u) {
  return log1m_inv_logit(static_cast<double>(u));
}

/**
 * Structure to wrap log1m_inv_logit() so it can be vectorized.
 */
struct log1m_inv_logit_fun {
  /**
   * Return the natural logarithm of one minus the inverse logit
   * of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return natural log of one minus inverse logit of argument
   */
  template <typename T>
  static inline T fun(const T& x) {
    return log1m_inv_logit(x);
  }
};

/**
 * Return the elementwise application of
 * <code>log1m_inv_logit()</code> to specified argument container.
 * The return type promotes the underlying scalar argument type to
 * double if it is an integer, and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise log1m_inv_logit of members of container.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto log1m_inv_logit(const T& x) {
  return apply_scalar_unary<log1m_inv_logit_fun, T>::apply(x);
}

/**
 * Version of log1m_inv_logit() that accepts Eigen Matrix/Array objects
 * or expressions. The input is cast to double (a no-op for objects already
 * of type double) as Eigen's log1p() does not take integer inputs.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return log(1 - inv_logit(x)) of each variable in the container.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto log1m_inv_logit(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    const auto& v_array = v.derived().template cast<double>().array().eval();
    return (v_array > 0.0).select(-v_array - (-v_array).exp().log1p(),
                                  -(v_array).exp().log1p());
  });
}

}  // namespace math
}  // namespace stan

#endif
