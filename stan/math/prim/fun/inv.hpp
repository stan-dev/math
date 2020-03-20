#ifndef STAN_MATH_PRIM_FUN_INV_HPP
#define STAN_MATH_PRIM_FUN_INV_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the inverse of the arithmetic variable.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 / x.
 */
template <typename T, require_arithmetic_t<T>...>
auto inv(const T& x) {
  return 1.0 / x;
}

/**
 * Version of \c inv() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return 1 divided by each value in x.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto inv(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().inverse(); });
}

}  // namespace math
}  // namespace stan

#endif
