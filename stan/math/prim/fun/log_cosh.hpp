#ifndef STAN_MATH_PRIM_FUN_LOG_COSH_HPP
#define STAN_MATH_PRIM_FUN_LOG_COSH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the hyperbolic cosine of the argument,
 * from `expm1` for `|x| < 1` and from `exp(-|x|)` above, so that neither
 * branch subtracts two nearly equal numbers or overflows. The derivative is
 * `tanh(x)`.
 *
 * @tparam T arithmetic type
 * @param x argument
 * @return `log(cosh(x))`
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline double log_cosh(T x) {
  const double abs_x = fabs(static_cast<double>(x));
  if (abs_x < 1.0) {
    const double expm1_abs_x = expm1(abs_x);
    return log1p(expm1_abs_x * expm1_abs_x / (2.0 * (1.0 + expm1_abs_x)));
  }
  const double exp_neg_abs_x = exp(-abs_x);
  return abs_x + log1p(0.5 * (exp_neg_abs_x * exp_neg_abs_x - 1.0));
}

/**
 * Structure to wrap `log_cosh()` so it can be vectorized.
 */
struct log_cosh_fun {
  /**
   * Return the natural logarithm of the hyperbolic cosine of the argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return `log(cosh(x))`
   */
  template <typename T>
  static inline auto fun(const T& x) {
    return log_cosh(x);
  }
};

/**
 * Return the elementwise `log_cosh()` of the specified container.
 *
 * @tparam T type of container
 * @param x container
 * @return elementwise `log(cosh(x))`
 */
template <typename T, require_container_t<T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr,
          require_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
inline auto log_cosh(T&& x) {
  return apply_scalar_unary<log_cosh_fun, T>::apply(std::forward<T>(x));
}

}  // namespace math
}  // namespace stan

#endif
