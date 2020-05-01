#ifndef STAN_MATH_OPENCL_PRIM_FUN_INV_SQRT_HPP
#define STAN_MATH_OPENCL_PRIM_FUN_INV_SQRT_HPP

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise `1 / sqrt(x)` of the specified kernel generator
 * expression.
 *
 * @param x input kernel generator expression
 * @return inverse square root of each value in x.
 */
template <typename T_x,
          typename = require_all_valid_expressions_and_none_scalar_t<T_x>>
inline auto inv_sqrt(T_x&& x) {  // NOLINT
  return rsqrt(std::forward<T_x>(x));
}
}  // namespace math
}  // namespace stan

#endif
