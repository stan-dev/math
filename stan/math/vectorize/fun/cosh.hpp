
#ifndef STAN_MATH_VECTORIZED_FUN_COSH_HPP
#define STAN_MATH_VECTORIZED_FUN_COSH_HPP
#include <stan/math/prim/fun/cosh.hpp>
namespace stan {
namespace math {

/**
 * Version of `cosh()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Hyberbolic cosine of x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto cosh(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().cosh(); });
}


} // namespace math
} // namespace stan
#endif 

