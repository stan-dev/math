#ifndef STAN_MATH_PRIM_FUN_FABS_HPP
#define STAN_MATH_PRIM_FUN_FABS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the absolute value of the arithmetic input.
 *
 * @tparam T type of arithmetic variable
 * @param x variable
 * @return Absolute value of variable.
 */
template <typename T, require_arithmetic_t<T>...>
auto fabs(const T& x) {
  using std::fabs;
  return fabs(x);
}

/**
 * Version of `fabs()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Absolute value of each value in x.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto fabs(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().abs(); });
}

}  // namespace math
}  // namespace stan

#endif
