#ifndef STAN_MATH_GPU_KERNELS_SCALAR_MUL_DIAGONAL_HPP
#define STAN_MATH_GPU_KERNELS_SCALAR_MUL_DIAGONAL_HPP
#ifdef STAN_OPENCL

#include <stan/math/gpu/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
const char *scalar_mul_diagonal_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Multiplication of the matrix A diagonal with a scalar
     *
     * @param[in, out] A matrix A
     * @param[in] scalar the value with which to multiply the diagonal of A
     * @param[in] rows the number of rows in A
     * @param[in] cols the number of columns in A
     */
    __kernel void scalar_mul_diagonal(
        __global read_write double *A, const read_only double scalar,
        const read_only unsigned int rows, const read_only unsigned int min_dim) {
      int i = get_global_id(0);
      if (i < min_dim) {
        A(i, i) *= scalar;
      }
    }
    // \cond
);
// \endcond

/**
 * See the docs for \link kernels/scalar_mul_diagonal.hpp add() \endlink
 */
const global_range_kernel<cl::Buffer, double, int, int> scalar_mul_diagonal(
    "scalar_mul_diagonal", scalar_mul_diagonal_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
