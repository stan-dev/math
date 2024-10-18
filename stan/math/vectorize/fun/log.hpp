
#ifndef STAN_MATH_VECTORIZED_FUN_LOG_HPP
#define STAN_MATH_VECTORIZED_FUN_LOG_HPP
#include <stan/math/prim/fun/log.hpp>
namespace stan {
namespace math {

/**
 * Version of `log()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Natural log of each variable in the container.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto log(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().log(); });
}


} // namespace math
} // namespace stan
#endif 

