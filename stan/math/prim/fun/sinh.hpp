#ifndef STAN_MATH_PRIM_FUN_SINH_HPP
#define STAN_MATH_PRIM_FUN_SINH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap sinh() so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic sine of x.
 */
struct sinh_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sinh;
    return sinh(x);
  }
};

/**
 * Vectorized version of sinh().
 *
 * @tparam T type of container
 * @param x container
 * @return Hyperbolic sine of each variable in x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto sinh(const T& x) {
  return apply_scalar_unary<sinh_fun, T>::apply(x);
}

/**
 * Version of sinh() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Hyperbolic sine of each variable in x.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto sinh(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    return v.derived().array().sinh();
  });
}

}  // namespace math
}  // namespace stan

#endif
