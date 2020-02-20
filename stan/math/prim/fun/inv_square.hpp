#ifndef STAN_MATH_PRIM_FUN_INV_SQUARE_HPP
#define STAN_MATH_PRIM_FUN_INV_SQUARE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/inv_square.hpp>

namespace stan {
namespace math {

inline double inv_square(double x) { return inv(square(x)); }

/**
 * Vectorized version of inv_square().
 *
 * @tparam T type of container
 * @param x container
 * @return 1 / the square of each value in x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto inv_square(const T& x) {
  return inv(square(x));
}

/**
 * Version of inv_square() that accepts Eigen Matrix/Array objects or
 * expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return 1 / the square of each value in x.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto inv_square(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [&](const auto& v) { return v.derived().array().square().inverse(); });
}

}  // namespace math
}  // namespace stan

#endif
