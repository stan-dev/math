
#ifndef STAN_MATH_VECTORIZED_FUN_ABS_HPP
#define STAN_MATH_VECTORIZED_FUN_ABS_HPP
#include <stan/math/prim/fun/abs.hpp>
namespace stan {
namespace math {

/**
 * Version of `abs()` that accepts std::vectors, Eigen Matrix/Array objects
 * or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x argument
 * @return Absolute value of each variable in the container.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto abs(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [&](const auto& v) { return v.array().abs(); });
}


} // namespace math
} // namespace stan
#endif 

