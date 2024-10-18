
#ifndef STAN_MATH_VECTORIZED_FUN_LOG10_HPP
#define STAN_MATH_VECTORIZED_FUN_LOG10_HPP
#include <stan/math/prim/fun/log10.hpp>
namespace stan {
namespace math {

/**
 * Version of log10() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Log base-10 of each variable in the container.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto log10(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().log10(); });
}


} // namespace math
} // namespace stan
#endif 

