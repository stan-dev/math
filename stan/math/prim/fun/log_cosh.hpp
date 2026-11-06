#ifndef STAN_MATH_PRIM_FUN_LOG_COSH_HPP
#define STAN_MATH_PRIM_FUN_LOG_COSH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {

namespace internal {
/**
 * Return `log(cosh(x))` and its derivative `tanh(x)`.
 *
 * For `|x| < 1`, with `u = expm1(|x|)`, `cosh(x) - 1 = u^2 / (2 (1 + u))`
 * and `tanh(|x|) = u (u + 2) / (u (u + 2) + 2)`. For `|x| >= 1`, with
 * `t = expm1(-2 |x|)`, `log(cosh(x)) = |x| + log1p(t / 2)` and
 * `tanh(|x|) = -t / (2 + t)`. No branch subtracts two nearly equal numbers
 * or overflows, so both results have a small relative error for every `x`.
 *
 * @param x argument
 * @return pair of `log(cosh(x))` and `tanh(x)`
 */
inline std::pair<double, double> log_cosh_and_tanh(double x) {
  const double a = std::fabs(x);
  double value;
  double tanh_a;
  if (a < 1.0) {
    const double u = std::expm1(a);
    const double w = u * (u + 2.0);
    value = std::log1p(u * u / (2.0 * (1.0 + u)));
    tanh_a = w / (w + 2.0);
  } else {
    const double t = std::expm1(-2.0 * a);
    value = a + std::log1p(0.5 * t);
    tanh_a = -t / (2.0 + t);
  }
  return {value, std::copysign(tanh_a, x)};
}
}  // namespace internal

/**
 * Return the natural logarithm of the hyperbolic cosine of the argument.
 *
 * The result does not overflow for large `|x|` and has a small relative
 * error for small `|x|`. The derivative is `tanh(x)`.
 *
 * @tparam T arithmetic type
 * @param x argument
 * @return `log(cosh(x))`
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline double log_cosh(T x) {
  return internal::log_cosh_and_tanh(static_cast<double>(x)).first;
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
inline auto log_cosh(const T& x) {
  return apply_scalar_unary<log_cosh_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
