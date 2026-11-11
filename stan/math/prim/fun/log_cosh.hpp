#ifndef STAN_MATH_PRIM_FUN_LOG_COSH_HPP
#define STAN_MATH_PRIM_FUN_LOG_COSH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/copysign.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <utility>

namespace stan {
namespace math {

namespace internal {
/**
 * `log(cosh(x))`, `tanh(x)` and `sech(x) = 1 / cosh(x)` of one argument.
 */
struct log_cosh_parts {
  double log_cosh;
  double tanh;
  double sech;
};

/**
 * Return `log(cosh(x))`, its derivative `tanh(x)`, and `sech(x)` from one
 * call of `expm1` (`|x| < 1`) or `exp` (`|x| >= 1`). Neither branch
 * subtracts two nearly equal numbers or overflows.
 *
 * @param x argument
 * @return `log(cosh(x))`, `tanh(x)` and `sech(x)`
 */
inline log_cosh_parts log_cosh_tanh_sech(double x) {
  const double a = fabs(x);
  if (a < 1.0) {
    const double u = expm1(a);
    const double d = u * u / (2.0 * (1.0 + u));
    const double w = u * (u + 2.0);
    return {log1p(d), copysign(w / (w + 2.0), x), 1.0 / (1.0 + d)};
  }
  const double e = exp(-a);
  const double e2 = e * e;
  return {a + log1p(0.5 * (e2 - 1.0)), copysign((1.0 - e2) / (1.0 + e2), x),
          2.0 * e / (1.0 + e2)};
}
}  // namespace internal

/**
 * Return the natural logarithm of the hyperbolic cosine of the argument,
 * without overflow for large `|x|`. The derivative is `tanh(x)`.
 *
 * @tparam T arithmetic type
 * @param x argument
 * @return `log(cosh(x))`
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline double log_cosh(T x) {
  return internal::log_cosh_tanh_sech(static_cast<double>(x)).log_cosh;
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
