#ifndef STAN_MATH_PRIM_FUN_ASIN_HPP
#define STAN_MATH_PRIM_FUN_ASIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap asin() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x argument
 * @return Arcsine of x in radians.
 */
struct asin_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::asin;
    return asin(x);
  }
};

/**
 * Vectorized version of asin().
 *
 * @tparam T type of container
 * @param x container
 * @return Arcsine of each variable in the container, in radians.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto asin(const T& x) {
  return apply_scalar_unary<asin_fun, T>::apply(x);
}

/**
 * Version of asin() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Arcsine of each variable in the container, in radians.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto asin(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    return v.derived().array().asin();
  });
}

}  // namespace math
}  // namespace stan

#endif
