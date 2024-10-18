
#ifndef STAN_MATH_VECTORIZED_FUN_ASIN_HPP
#define STAN_MATH_VECTORIZED_FUN_ASIN_HPP
#include <stan/math/prim/fun/asin.hpp>
namespace stan {
namespace math {

/**
 * Version of `asin()` that accepts std::vectors, Eigen Matrix/Array objects,
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Arcsine of each variable in the container, in radians.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto asin(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().asin(); });
}


} // namespace math
} // namespace stan
#endif 

