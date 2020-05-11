#ifndef STAN_MATH_OPENCL_SUB_BLOCK_HPP
#define STAN_MATH_OPENCL_SUB_BLOCK_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator/block.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <CL/cl2.hpp>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/** \ingroup matrix_cl_group
 * Write the contents of A into
 * `this` starting at the top left of `this`
 * @param A input matrix
 * @param A_i the offset row in A
 * @param A_j the offset column in A
 * @param this_i the offset row for the matrix to be subset into
 * @param this_j the offset col for the matrix to be subset into
 * @param nrows the number of rows in the submatrix
 * @param ncols the number of columns in the submatrix
 */
template <typename T>
inline void matrix_cl<T, require_arithmetic_t<T>>::sub_block(
    const matrix_cl<T, require_arithmetic_t<T>>& A, size_t A_i, size_t A_j,
    size_t this_i, size_t this_j, size_t nrows, size_t ncols) {
  block(*this, this_i, this_j, nrows, ncols) = block(A, A_i, A_j, nrows, ncols);
}

}  // namespace math
}  // namespace stan

#endif
#endif
