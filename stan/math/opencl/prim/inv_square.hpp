#ifndef STAN_MATH_OPENCL_PRIM_FUN_INV_SQUARE_HPP
#define STAN_MATH_OPENCL_PRIM_FUN_INV_SQUARE_HPP

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise `1 / square(x)` of the specified kernel generator expression.
 *
 * @param x input kernel generator expression
 * @return inverse square of each value in x.
 */
template <typename T_a,
          typename = require_all_valid_expressions_and_none_scalar_t<T_a>>
inline auto inv_square(T_a&& a) {  // NOLINT
  return elewise_division(1.0, elewise_multiplication(std::forward<T_a>(a), std::forward<T_a>(a)));
}
}  // namespace math
}  // namespace stan

#endif
