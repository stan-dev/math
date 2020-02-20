#ifndef STAN_MATH_PRIM_FUN_ATAN_HPP
#define STAN_MATH_PRIM_FUN_ATAN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap atan() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Arctan of x in radians.
 */
struct atan_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::atan;
    return atan(x);
  }
};

/**
 * Vectorized version of atan().
 *
 * @tparam T type of container
 * @param x container
 * @return Arctan of each value in x, in radians.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto atan(const T& x) {
  return apply_scalar_unary<atan_fun, T>::apply(x);
}

/**
 * Version of atan() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Elementwise atan of members of container.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto atan(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [&](const auto& v) { return v.derived().array().atan(); });
}

}  // namespace math
}  // namespace stan

#endif
