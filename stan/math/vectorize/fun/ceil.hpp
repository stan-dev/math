
#ifndef STAN_MATH_VECTORIZED_FUN_CEIL_HPP
#define STAN_MATH_VECTORIZED_FUN_CEIL_HPP
#include <stan/math/prim/fun/ceil.hpp>
namespace stan {
namespace math {

/**
 * Version of `ceil()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Least integer >= each value in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr>
inline auto ceil(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().ceil(); });
}


} // namespace math
} // namespace stan
#endif 

