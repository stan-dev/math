#ifndef STAN_MATH_OPENCL_PRIM_COLS_HPP
#define STAN_MATH_OPENCL_PRIM_COLS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Returns the number of columns in matrix_cl.
 *
 * @tparam T_a type of input kernel generator expression a
 * @tparam T type of elements in the input matrix
 * @param x the input matrix_cl
 *
 * @return number of columns in x
 */
template <typename T_x,
          typename = require_all_valid_expressions_and_none_scalar_t<T_x>>
inline int cols(const matrix_cl<T>& x) {
  return x.cols();
}

}  // namespace math
}  // namespace stan

#endif
#endif
