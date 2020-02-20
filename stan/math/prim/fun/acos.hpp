#ifndef STAN_MATH_PRIM_FUN_ACOS_HPP
#define STAN_MATH_PRIM_FUN_ACOS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap acos() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Arc cosine of variable in radians.
 */
struct acos_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::acos;
    return acos(x);
  }
};

/**
 * Vectorized version of acos().
 *
 * @tparam T type of container
 * @param x container
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto acos(const T& x) {
  return apply_scalar_unary<acos_fun, T>::apply(x);
}

/**
 * Version of acos() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto acos(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    return v.derived().array().acos();
  });
}

}  // namespace math
}  // namespace stan

#endif
