#ifndef STAN_MATH_PRIM_FUN_INV_HPP
#define STAN_MATH_PRIM_FUN_INV_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

inline double inv(double x) { return 1.0 / x; }

/**
 * Structure to wrap inv() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 / x.
 */
struct inv_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv(x);
  }
};

/**
 * Vectorized version of inv().
 *
 * @tparam T type of container
 * @param x container
 * @return 1 divided by each value in x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto inv(const T& x) {
  return apply_scalar_unary<inv_fun, T>::apply(x);
}

/**
 * Version of inv() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @param x Matrix or matrix expression
 * @return 1 divided by each value in x.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto inv(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [&](const auto& v) { return v.derived().array().inverse(); });
}

}  // namespace math
}  // namespace stan

#endif
