#ifndef STAN_MATH_OPENCL_PRIM_ROWS_HPP
#define STAN_MATH_OPENCL_PRIM_ROWS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <cstdint>

namespace stan {
namespace math {
/** \ingroup opencl
 * Returns the number of rows in the specified kernel generator
 * expression.
 *
 * @tparam T_x type of input kernel generator expression x
 * @param x input kernel generator expression.
 *
 * @return number of rows in x
 */

template <typename T_x,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr>
inline int64_t rows(const T_x& x) {
  return x.rows();
}

}  // namespace math
}  // namespace stan

#endif
#endif
