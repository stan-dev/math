#ifndef STAN_MATH_PRIM_FUN_FLOOR_HPP
#define STAN_MATH_PRIM_FUN_FLOOR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap `floor()` so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Greatest integer <= x.
 */
template <typename T, require_arithmetic_t<T>...>
auto floor(const T& x) {
  using std::floor;
  return floor(x);
}

/**
 * Version of `floor()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Greatest integer <= each value in x.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto floor(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().floor(); });
}

}  // namespace math
}  // namespace stan

#endif
