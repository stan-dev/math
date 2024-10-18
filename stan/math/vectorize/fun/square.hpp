
#ifndef STAN_MATH_VECTORIZED_FUN_SQUARE_HPP
#define STAN_MATH_VECTORIZED_FUN_SQUARE_HPP
#include <stan/math/prim/fun/square.hpp>
namespace stan {
namespace math {

/**
 * Version of square() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Each value in x squared.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto square(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().square(); });
}


} // namespace math
} // namespace stan
#endif 

