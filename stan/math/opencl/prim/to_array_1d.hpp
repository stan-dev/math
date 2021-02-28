#ifndef STAN_MATH_OPENCL_PRIM_TO_ARRAY_1D_HPP
#define STAN_MATH_OPENCL_PRIM_TO_ARRAY_1D_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/to_matrix.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns input matrix reshaped into a vector.
 *
 * @tparam T_x type of the matrix
 *
 * @param x matrix
 * @return the vector representation of the input
 */
template <typename T_x,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr>
inline auto to_array_1d(T_x&& x) {
  return to_matrix(std::forward<T_x>(x), x.size(), 1);
}

}  // namespace math
}  // namespace stan
#endif
#endif
