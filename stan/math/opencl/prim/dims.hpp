#ifndef STAN_MATH_OPENCL_PRIM_DIMS_HPP
#define STAN_MATH_OPENCL_PRIM_DIMS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Returns a vector of matrix_cl dimensions.
 *
 * @tparam T type of elements in the input matrix
 * @param x the input 1x1 matrix_cl
 * @param m number of rows in the results row_vector
 *
 * @return matrix_cl with replicated value from the input matrix
 *
 * @throw <code>domain_error</code> if the
 * requested dimensions are negative
 * @throw <code>invalid_argument</code> if input
 * element is not a matrix_cl of size 1
 *
 */
template <typename T>
inline std::vector<int> dims(const matrix_cl<T>& x) {
  return {x.rows(), x.cols()};
}

}  // namespace math
}  // namespace stan

#endif
#endif
