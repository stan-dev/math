#ifndef STAN_MATH_PRIM_FUN_COS_HPP
#define STAN_MATH_PRIM_FUN_COS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the cosine of the arithmetic input.
 *
 * @tparam T type of arithmetic variable
 * @param x variable
 * @return Cosine of variable in radians.
 */
template <typename T, require_arithmetic_t<T>...>
auto cos(const T& x) {
  using std::cos;
  return cos(x);
}

/**
 * Version of `cos()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Cosine of each value in x.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto cos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().cos(); });
}

}  // namespace math
}  // namespace stan

#endif
