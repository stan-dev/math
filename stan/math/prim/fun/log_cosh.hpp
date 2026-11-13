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

namespace internal {
/**
 * Return `sech(x) = 1 / cosh(x)` as `2 e / (1 + e^2)` with `e = exp(-|x|)`,
 * which underflows to 0 instead of overflowing for large `|x|`.
 *
 * @param x argument
 * @return `sech(x)`
 */
inline double sech(double x) {
  const double e = exp(-fabs(x));
  return 2.0 * e / (1.0 + e * e);
}
}  // namespace internal

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
  const double a = fabs(static_cast<double>(x));
  if (a < 1.0) {
    const double u = expm1(a);
    return log1p(u * u / (2.0 * (1.0 + u)));
  }
  const double e = exp(-a);
  return a + log1p(0.5 * (e * e - 1.0));
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
