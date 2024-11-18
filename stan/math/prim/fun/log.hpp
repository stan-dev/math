#ifndef STAN_MATH_PRIM_FUN_LOG_HPP
#define STAN_MATH_PRIM_FUN_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the arithmetic argument.
 *
 * @tparam V `Arithmetic` argument
 * @param[in] x argument
 * @return natural logarithm of the argument
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto log(const T x) {
  return std::log(x);
}

/**
 * Return the natural logarithm of the complex argument.
 *
 * @tparam V `complex<Arithmetic>` argument
 * @param[in] x argument
 * @return natural logarithm of the argument
 */
template <typename T, require_complex_t<T>* = nullptr>
inline auto log(const T x) {
  return std::log(x);
}

/**
 * Structure to wrap `log()` so that it can be vectorized.
 */
struct log_fun {
  /**
   * Return natural log of specified argument.
   *
   * @tparam T type of argument
   * @param[in] x argument
   * @return Natural log of x.
   */
  template <typename T>
  static inline auto fun(const T& x) {
    return log(x);
  }
};

/**
 * Return the elementwise natural log of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 * The return type is the same as the argument type.
 *
 * @tparam Container type of container
 * @param[in] x container
 * @return Elementwise application of natural log to the argument.
 */
template <typename Container, require_ad_container_t<Container>* = nullptr>
inline auto log(const Container& x) {
  return apply_scalar_unary<log_fun, Container>::apply(x);
}

/**
 * Version of `log()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Natural log of each variable in the container.
 */
template <typename Container,
          require_container_bt<std::is_arithmetic, Container>* = nullptr>
inline auto log(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().log(); });
}

namespace internal {
/**
 * Return the natural logarithm of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return natural logarithm of the argument
 */
template <typename V>
inline std::complex<V> complex_log(const std::complex<V>& z) {
  if ((is_nan(z.real()) && is_inf(z.imag()))
      || (is_inf(z.real()) && is_nan(z.imag()))) {
    return {INFTY, NOT_A_NUMBER};
  }
  V r = sqrt(norm(z));
  V theta = arg(z);
  return {log(r), theta};
}
}  // namespace internal

}  // namespace math
}  // namespace stan
#endif
