
#ifndef STAN_MATH_VECTORIZED_FUN_SQRT_HPP
#define STAN_MATH_VECTORIZED_FUN_SQRT_HPP
#include <stan/math/prim/fun/sqrt.hpp>
namespace stan {
namespace math {

/**
 * Version of `sqrt()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Square root of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr>
inline auto sqrt(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().sqrt(); });
}


} // namespace math
} // namespace stan
#endif 

