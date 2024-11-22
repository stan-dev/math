#ifndef STAN_MATH_PRIM_FUN_SINH_HPP
#define STAN_MATH_PRIM_FUN_SINH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the hyperbolic sine of the arithmetic argument.
 *
 * @tparam V An arithmetic argument
 * @param[in] x argument
 * @return hyperbolic sine of the argument
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto sinh(const T x) {
  return std::sinh(x);
}

/**
 * Return the hyperbolic sine of the complex argument.
 *
 * @tparam V `complex<Arithmetic>` argument
 * @param[in] x argument
 * @return hyperbolic sine of the argument
 */
template <typename T, require_complex_bt<std::is_arithmetic, T>* = nullptr>
inline auto sinh(const T x) {
  return std::sinh(x);
}

/**
 * Structure to wrap sinh() so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic sine of x.
 */
struct sinh_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return sinh(x);
  }
};

/**
 * Vectorized version of sinh().
 *
 * @tparam Container type of container
 * @param x container
 * @return Hyperbolic sine of each variable in x.
 */
template <typename Container, require_ad_container_t<Container>* = nullptr>
inline auto sinh(const Container& x) {
  return apply_scalar_unary<sinh_fun, Container>::apply(x);
}

/**
 * Version of sinh() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Hyperbolic sine of each variable in x.
 */
template <typename Container,
          require_container_bt<std::is_arithmetic, Container>* = nullptr>
inline auto sinh(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().sinh(); });
}

namespace internal {
/**
 * Return the hyperbolic sine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic sine of the argument
 */
template <typename V>
inline std::complex<V> complex_sinh(const std::complex<V>& z) {
  return 0.5 * (exp(z) - exp(-z));
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
