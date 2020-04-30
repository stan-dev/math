#ifndef STAN_MATH_OPENCL_PRIM_FUN_INV_SQUARE_HPP
#define STAN_MATH_OPENCL_PRIM_FUN_INV_SQUARE_HPP

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise 1.0 / x of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @param x Input kernel generator expression.
 * @return elementwise 1.0 / x of the specified argument.
 */
template <typename T_a,
          typename = require_all_valid_expressions_and_none_scalar_t<T_a>>
inline auto inv(T_a&& a) {  // NOLINT
  return elewise_division(1.0, a);
}
}  // namespace math
}  // namespace stan

#endif