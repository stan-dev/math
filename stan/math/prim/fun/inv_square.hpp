#ifndef STAN_MATH_PRIM_FUN_INV_SQUARE_HPP
#define STAN_MATH_PRIM_FUN_INV_SQUARE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Returns `1 / square(x)`.
 *
 * @tparam Container type of container
 * @param x container
 * @return `1 / square(x)` of each value in x.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr>
inline auto inv_square(const Container& x) {
  return inv(square(x));
}

/**
 * Version of inv_square() that accepts Eigen Matrix/Array objects or
 * expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return 1 / the square of each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto inv_square(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().square().inverse(); });
}

}  // namespace math
}  // namespace stan

#endif
