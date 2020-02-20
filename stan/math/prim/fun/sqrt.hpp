#ifndef STAN_MATH_PRIM_FUN_SQRT_HPP
#define STAN_MATH_PRIM_FUN_SQRT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the square root of the specified argument.  This
 * version is required to disambiguate <code>sqrt(int)</code>.
 *
 * @param[in] x Argument.
 * @return Natural exponential of argument.
 */
inline double sqrt(int x) { return std::sqrt(x); }

/**
 * Structure to wrap sqrt() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Square root of x.
 */
struct sqrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sqrt;
    return sqrt(x);
  }
};

/**
 * Vectorized version of sqrt().
 *
 * @tparam T type of container
 * @param x container
 * @return Square root of each value in x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto sqrt(const T& x) {
  return apply_scalar_unary<sqrt_fun, T>::apply(x);
}

/**
 * Version of sqrt() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Square root of each value in x.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto sqrt(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [&](const auto& v) { return v.derived().array().sqrt(); });
}

}  // namespace math
}  // namespace stan

#endif
