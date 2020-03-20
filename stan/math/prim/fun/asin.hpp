#ifndef STAN_MATH_PRIM_FUN_ASIN_HPP
#define STAN_MATH_PRIM_FUN_ASIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the arcsine of the arithmetic input.
 *
 * @tparam T type of arithmetic variable
 * @param x variable
 * @return Arcsine of variable in radians.
 */
template <typename T, require_arithmetic_t<T>...>
auto asin(const T& x) {
  using std::asin;
  return asin(x);
}

/**
 * Version of `asin()` that accepts std::vectors, Eigen Matrix/Array objects,
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Arcsine of each variable in the container, in radians.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto asin(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().asin(); });
}

}  // namespace math
}  // namespace stan

#endif
