#ifndef STAN_MATH_PRIM_FUN_FLOOR_HPP
#define STAN_MATH_PRIM_FUN_FLOOR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap floor() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Greatest integer <= x.
 */
struct floor_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::floor;
    return floor(x);
  }
};

/**
 * Vectorized version of floor().
 *
 * @tparam T type of container
 * @param x container
 * @return Greatest integer <= each value in x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto floor(const T& x) {
  return apply_scalar_unary<floor_fun, T>::apply(x);
}

/**
 * Version of floor() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Greatest integer <= each value in x.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto floor(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    return v.derived().array().floor();
  });
}

}  // namespace math
}  // namespace stan

#endif
