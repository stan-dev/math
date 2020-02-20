#ifndef STAN_MATH_PRIM_FUN_TAN_HPP
#define STAN_MATH_PRIM_FUN_TAN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap tan() so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Tangent of x.
 */
struct tan_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::tan;
    return tan(x);
  }
};

/**
 * Vectorized version of tan().
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Tangent of each value in x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto tan(const T& x) {
  return apply_scalar_unary<tan_fun, T>::apply(x);
}

/**
 * Version of tan() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Tangent of each value in x.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto tan(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    return v.derived().array().tan();
  });
}

}  // namespace math
}  // namespace stan

#endif
