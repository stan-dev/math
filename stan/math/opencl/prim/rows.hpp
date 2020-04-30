#ifndef STAN_MATH_OPENCL_PRIM_COLS_HPP
#define STAN_MATH_OPENCL_PRIM_COLS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Returns the number of rows in matrix_cl.
 *
 * @tparam T type of elements in the input matrix
 * @param x the input matrix_cl
 *
 * @return number of rows in x
 */
template <typename T>
inline int rows(const matrix_cl<T>& x) {
  return x.rows();
}

}  // namespace math
}  // namespace stan

#endif
#endif
