
#ifndef STAN_MATH_VECTORIZED_FUN_COS_HPP
#define STAN_MATH_VECTORIZED_FUN_COS_HPP
#include <stan/math/prim/fun/cos.hpp>
namespace stan {
namespace math {

/**
 * Version of `cos()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Cosine of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto cos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [&](const auto& v) { return v.array().cos(); });
}


} // namespace math
} // namespace stan
#endif 

