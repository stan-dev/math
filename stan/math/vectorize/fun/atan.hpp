
#ifndef STAN_MATH_VECTORIZED_FUN_ATAN_HPP
#define STAN_MATH_VECTORIZED_FUN_ATAN_HPP
#include <stan/math/prim/fun/atan.hpp>
namespace stan {
namespace math {

/**
 * Version of atan() that accepts std::vectors, Eigen Matrix/Array objects,
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Elementwise atan of members of container.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto atan(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().atan(); });
}


} // namespace math
} // namespace stan
#endif 

