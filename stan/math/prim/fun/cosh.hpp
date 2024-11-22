#ifndef STAN_MATH_PRIM_FUN_COSH_HPP
#define STAN_MATH_PRIM_FUN_COSH_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the hyperbolic cosine of the arithmetic argument.
 *
 * @tparam V An arithmetic type
 * @param[in] x argument
 * @return hyperbolic cosine of the argument
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto cosh(const T x) {
  return std::cosh(x);
}

/**
 * Return the hyperbolic cosine of the complex argument.
 *
 * @tparam V `complex<Arithmetic>` type of argument
 * @param[in] x argument
 * @return hyperbolic cosine of the argument
 */
template <typename T, require_complex_bt<std::is_arithmetic, T>* = nullptr>
inline auto cosh(const T x) {
  return std::cosh(x);
}

/**
 * Structure to wrap `cosh()` so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic cosine of x.
 */
struct cosh_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return cosh(x);
  }
};

/**
 * Returns the elementwise `cosh()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x angles in radians
 * @return Hyberbolic cosine of x.
 */
template <typename Container, require_ad_container_t<Container>* = nullptr>
inline auto cosh(const Container& x) {
  return apply_scalar_unary<cosh_fun, Container>::apply(x);
}

/**
 * Version of `cosh()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Hyberbolic cosine of x.
 */
template <typename Container,
          require_container_bt<std::is_arithmetic, Container>* = nullptr>
inline auto cosh(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().cosh(); });
}

namespace internal {
/**
 * Return the hyperbolic cosine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic cosine of the argument
 */
template <typename V>
inline std::complex<V> complex_cosh(const std::complex<V>& z) {
  return 0.5 * (exp(z) + exp(-z));
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
