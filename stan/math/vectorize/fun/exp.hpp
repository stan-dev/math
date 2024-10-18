
#ifndef STAN_MATH_VECTORIZED_FUN_EXP_HPP
#define STAN_MATH_VECTORIZED_FUN_EXP_HPP
#include <stan/math/prim/fun/exp.hpp>
namespace stan {
namespace math {

/**
 * Version of `exp()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Elementwise application of exponentiation to the argument.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto exp(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().exp(); });
}


} // namespace math
} // namespace stan
#endif 

