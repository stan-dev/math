
#ifndef STAN_MATH_VECTORIZED_FUN_ACOS_HPP
#define STAN_MATH_VECTORIZED_FUN_ACOS_HPP
#include <stan/math/prim/fun/acos.hpp>
namespace stan {
namespace math {

/**
 * Version of `acos()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x argument
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto acos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().acos(); });
}


} // namespace math
} // namespace stan
#endif 

