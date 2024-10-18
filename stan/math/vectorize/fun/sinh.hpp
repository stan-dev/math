
#ifndef STAN_MATH_VECTORIZED_FUN_SINH_HPP
#define STAN_MATH_VECTORIZED_FUN_SINH_HPP
#include <stan/math/prim/fun/sinh.hpp>
namespace stan {
namespace math {

/**
 * Version of sinh() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Hyperbolic sine of each variable in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto sinh(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().sinh(); });
}


} // namespace math
} // namespace stan
#endif 

