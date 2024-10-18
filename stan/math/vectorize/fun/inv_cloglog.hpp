
#ifndef STAN_MATH_VECTORIZED_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_VECTORIZED_FUN_INV_CLOGLOG_HPP
#include <stan/math/prim/fun/inv_cloglog.hpp>
namespace stan {
namespace math {

/**
 * Version of inv_cloglog() that accepts std::vectors, Eigen Matrix/Array
 * objects or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return 1 - exp(-exp()) applied to each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto inv_cloglog(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return 1 - (-v.array().exp()).exp(); });
}


} // namespace math
} // namespace stan
#endif 

