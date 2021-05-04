#ifndef STAN_MATH_OPENCL_PRIM_TO_ROW_VECTOR_HPP
#define STAN_MATH_OPENCL_PRIM_TO_ROW_VECTOR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/to_matrix.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns input matrix reshaped into a row vector.
 *
 * @tparam T_x type of the matrix
 *
 * @param x matrix
 * @return the row vector representation of the input
 */
template <typename T_x,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr>
inline auto to_row_vector(T_x&& x) {
  return to_matrix(std::forward<T_x>(x), 1, x.size());
}

}  // namespace math
}  // namespace stan
#endif
#endif
