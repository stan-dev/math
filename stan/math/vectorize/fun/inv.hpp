
#ifndef STAN_MATH_VECTORIZED_FUN_INV_HPP
#define STAN_MATH_VECTORIZED_FUN_INV_HPP
#include <stan/math/prim/fun/inv.hpp>
namespace stan {
namespace math {

/**
 * Version of \c inv() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return 1 divided by each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr>
inline auto inv(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().inverse(); });
}


} // namespace math
} // namespace stan
#endif 

