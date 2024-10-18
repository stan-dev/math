
#ifndef STAN_MATH_VECTORIZED_FUN_TAN_HPP
#define STAN_MATH_VECTORIZED_FUN_TAN_HPP
#include <stan/math/prim/fun/tan.hpp>
namespace stan {
namespace math {

/**
 * Version of `tan()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Tangent of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto tan(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().tan(); });
}


} // namespace math
} // namespace stan
#endif 

