#ifndef STAN_MATH_OPENCL_PRIM_FUN_INV_SQUARE_HPP
#define STAN_MATH_OPENCL_PRIM_FUN_INV_SQUARE_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise 1.0 / x of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T_x type of input kernel generator expression x
 * @param x Input kernel generator expression.
 * @return elementwise 1.0 / x of the specified argument.
 */
template <typename T_x,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_x>>
inline auto inv(T_x&& x) {  // NOLINT
  return elt_divide(1.0, std::forward<T_x>(x));
}
}  // namespace math
}  // namespace stan

#endif
#endif
