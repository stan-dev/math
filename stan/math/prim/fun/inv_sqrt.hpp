#ifndef STAN_MATH_PRIM_FUN_INV_SQRT_HPP
#define STAN_MATH_PRIM_FUN_INV_SQRT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto inv_sqrt(const T x) {
  return inv(std::sqrt(x));
}

template <typename T, require_complex_t<T>* = nullptr>
inline auto inv_sqrt(const T x) {
  return inv(sqrt(x));
}

/**
 * Structure to wrap `1 / sqrt(x)` so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return inverse square root of x.
 */
struct inv_sqrt_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return inv_sqrt(x);
  }
};

/**
 * Return the elementwise `1 / sqrt(x)}` of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x container
 * @return inverse square root of each value in x.
 */
template <typename Container, require_ad_container_t<Container>* = nullptr>
inline auto inv_sqrt(const Container& x) {
  return apply_scalar_unary<inv_sqrt_fun, Container>::apply(x);
}

/**
 * Version of `inv_sqrt()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return inverse square root each variable in the container.
 */
template <typename Container, require_not_var_matrix_t<Container>* = nullptr,
          require_container_bt<std::is_arithmetic, Container>* = nullptr>
inline auto inv_sqrt(const Container& x) {
// Eigen 3.4.0 has precision issues on ARM64 with vectorised rsqrt
// Resolved in current master branch, below can be removed on next release
#ifdef __aarch64__
  return apply_scalar_unary<inv_sqrt_fun, Container>::apply(x);
#else
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().rsqrt(); });
#endif
}

}  // namespace math
}  // namespace stan

#endif
