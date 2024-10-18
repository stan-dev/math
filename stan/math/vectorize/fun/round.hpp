
#ifndef STAN_MATH_VECTORIZED_FUN_ROUND_HPP
#define STAN_MATH_VECTORIZED_FUN_ROUND_HPP
#include <stan/math/prim/fun/round.hpp>
namespace stan {
namespace math {

/**
 * Version of `round()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Rounded value of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto round(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().round(); });
}


} // namespace math
} // namespace stan
#endif 

