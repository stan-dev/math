
#ifndef STAN_MATH_VECTORIZED_FUN_TANH_HPP
#define STAN_MATH_VECTORIZED_FUN_TANH_HPP
#include <stan/math/prim/fun/tanh.hpp>
namespace stan {
namespace math {

/**
 * Version of `tanh()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Hyperbolic tangent of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto tanh(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().tanh(); });
}


} // namespace math
} // namespace stan
#endif 

