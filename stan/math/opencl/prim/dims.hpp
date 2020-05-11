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
 * @tparam T_x type of input kernel generator expression a
 * @param x the input matrix_cl
 *
 * @return std::vector of the dimensions of the input kernel generato expression
 */
template <typename T_x,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_x>>
inline std::vector<int> dims(const T_x& x) {
  return {x.rows(), x.cols()};
}
}  // namespace math
}  // namespace stan

#endif
#endif
