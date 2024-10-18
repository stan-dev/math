
#ifndef STAN_MATH_VECTORIZED_FUN_MINUS_HPP
#define STAN_MATH_VECTORIZED_FUN_MINUS_HPP
#include <stan/math/prim/fun/minus.hpp>
namespace stan {
namespace math {

/**
 * Return the negation of the each element of a vector
 *
 * @tparam T Type of container.
 * @param x Container.
 * @return Container where each element is negated.
 */
template <typename T>
inline auto minus(const std::vector<T>& x) {
  return apply_vector_unary<std::vector<T>>::apply(
      x, [](const auto& v) { return -v; });
}


} // namespace math
} // namespace stan
#endif 

