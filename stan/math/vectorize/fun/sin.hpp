
#ifndef STAN_MATH_VECTORIZED_FUN_SIN_HPP
#define STAN_MATH_VECTORIZED_FUN_SIN_HPP
#include <stan/math/prim/fun/sin.hpp>
namespace stan {
namespace math {

/**
 * Version of sin() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Sine of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto sin(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [&](const auto& v) { return v.array().sin(); });
}


} // namespace math
} // namespace stan
#endif 

