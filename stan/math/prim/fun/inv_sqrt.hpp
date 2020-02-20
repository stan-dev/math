#ifndef STAN_MATH_PRIM_FUN_INV_SQRT_HPP
#define STAN_MATH_PRIM_FUN_INV_SQRT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

inline double inv_sqrt(double x) {
  using std::sqrt;
  return inv(sqrt(x));
}

/**
 * Structure to wrap inv_sqrt() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 / sqrt of x.
 */
struct inv_sqrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_sqrt(x);
  }
};

/**
 * Vectorized version of inv_sqrt().
 *
 * @tparam T type of container
 * @param x container
 * @return 1 / sqrt of each value in x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto inv_sqrt(const T& x) {
  return apply_scalar_unary<inv_sqrt_fun, T>::apply(x);
}

/**
 * Version of inv_sqrt() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto inv_sqrt(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    return v.derived().array().rsqrt();
  });
}

}  // namespace math
}  // namespace stan

#endif
