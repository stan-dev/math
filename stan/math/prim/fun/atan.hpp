#ifndef STAN_MATH_PRIM_FUN_ATAN_HPP
#define STAN_MATH_PRIM_FUN_ATAN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the arc tangent of the arithmetic input.
 *
 * @tparam T type of arithmetic variable
 * @param x variable
 * @return Arc tangent of variable.
 */
template <typename T, require_arithmetic_t<T>...>
auto atan(const T& x) {
  using std::atan;
  return atan(x);
}

/**
 * Version of atan() that accepts std::vectors, Eigen Matrix/Array objects,
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Elementwise atan of members of container.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto atan(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().atan(); });
}

}  // namespace math
}  // namespace stan

#endif
