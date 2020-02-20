#ifndef STAN_MATH_PRIM_FUN_SIN_HPP
#define STAN_MATH_PRIM_FUN_SIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap sin() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Sine of x.
 */
struct sin_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sin;
    return sin(x);
  }
};

/**
 * Vectorized version of sin().
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Sine of each value in x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto sin(const T& x) {
  return apply_scalar_unary<sin_fun, T>::apply(x);
}

/**
 * Version of sin() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Sine of each value in x.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto sin(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    return v.derived().array().sin();
  });
}

}  // namespace math
}  // namespace stan

#endif
