#ifndef STAN_MATH_OPENCL_PRIM_DIMS_HPP
#define STAN_MATH_OPENCL_PRIM_DIMS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>
#include <vector>

namespace stan {
namespace math {
/** \ingroup opencl
 * Returns a vector of matrix_cl dimensions.
 *
 * @tparam T type of elements in the input matrix
 * @param x the input matrix_cl
 *
 * @return std::vector of matrix_cl dimensions
 */
template <typename T>
inline std::vector<int> dims(const matrix_cl<T>& x) {
  return {x.rows(), x.cols()};
}

}  // namespace math
}  // namespace stan

#endif
#endif
