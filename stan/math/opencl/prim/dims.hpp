#ifndef STAN_MATH_OPENCL_PRIM_DIMS_HPP
#define STAN_MATH_OPENCL_PRIM_DIMS_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/fun/dims.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <vector>

namespace stan {
namespace math {
/** \ingroup opencl
 * matrix_cl overload of the dims helper function in prim/fun/dims.hpp.
 * Pushes the rows and columns to the result vector argument.
 *
 * @tparam T_x type of input kernel generator expression a
 * @param[in] x the input matrix_cl
 * @param[out] result the output vector of dimensions
 */
template <typename T_x,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr>
inline void dims(const T_x& x, std::vector<int>& result) {
  result.push_back(x.rows());
  result.push_back(x.cols());
}
}  // namespace math
}  // namespace stan

#endif
#endif
