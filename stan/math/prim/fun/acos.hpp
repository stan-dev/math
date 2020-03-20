#ifndef STAN_MATH_PRIM_FUN_ACOS_HPP
#define STAN_MATH_PRIM_FUN_ACOS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the arc cosine of the arithmetic input.
 *
 * @tparam T type of arithmetic variable
 * @param x variable
 * @return Arc cosine of variable in radians.
 */
template <typename T, require_arithmetic_t<T>...>
auto acos(const T& x) {
  using std::acos;
  return acos(x);
}

/**
 * Version of `acos()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto acos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().acos(); });
}

}  // namespace math
}  // namespace stan

#endif
