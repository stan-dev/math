#ifndef STAN_MATH_OPENCL_PRIM_COLS_HPP
#define STAN_MATH_OPENCL_PRIM_COLS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Returns the number of columns in the specified kernel generator
 * expression.
 *
 * @tparam T_x type of input kernel generator expression x
 * @param x input kernel generator expression.
 *
 * @return number of columns in x
 */
template <typename T_x,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_x>>
inline int cols(const T_x& x) {
  return x.cols();
}
}  // namespace math
}  // namespace stan

#endif
#endif
