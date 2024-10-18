
#ifndef STAN_MATH_VECTORIZED_FUN_TO_INT_HPP
#define STAN_MATH_VECTORIZED_FUN_TO_INT_HPP
#include <stan/math/prim/fun/to_int.hpp>
namespace stan {
namespace math {

/**
 * Returns the elementwise `to_int()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x argument
 * @return Integer value of each variable in the container.
 */
template <typename Container,
          require_std_vector_st<std::is_arithmetic, Container>* = nullptr>
inline auto to_int(const Container& x) {
  return apply_scalar_unary<to_int_fun, Container>::apply(x);
}


} // namespace math
} // namespace stan
#endif 

