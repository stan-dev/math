#ifndef STAN_MATH_PRIM_FUN_CEIL_HPP
#define STAN_MATH_PRIM_FUN_CEIL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap `ceil()` so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Least integer >= x.
 */
template <typename T, require_arithmetic_t<T>...>
auto ceil(const T& x) {
  using std::ceil;
  return ceil(x);
}

/**
 * Version of `ceil()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Least integer >= each value in x.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto ceil(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().ceil(); });
}

}  // namespace math
}  // namespace stan

#endif
