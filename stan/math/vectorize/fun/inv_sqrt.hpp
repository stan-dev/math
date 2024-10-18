
#ifndef STAN_MATH_VECTORIZED_FUN_INV_SQRT_HPP
#define STAN_MATH_VECTORIZED_FUN_INV_SQRT_HPP
#include <stan/math/prim/fun/inv_sqrt.hpp>
namespace stan {
namespace math {

/**
 * Version of `inv_sqrt()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return inverse square root each variable in the container.
 */
template <typename Container, require_not_var_matrix_t<Container>* = nullptr,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto inv_sqrt(const Container& x) {
// Eigen 3.4.0 has precision issues on ARM64 with vectorised rsqrt
// Resolved in current master branch, below can be removed on next release
#ifdef __aarch64__
  return apply_scalar_unary<inv_sqrt_fun, Container>::apply(x);
#else
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().rsqrt(); });
#endif
}


} // namespace math
} // namespace stan
#endif 

