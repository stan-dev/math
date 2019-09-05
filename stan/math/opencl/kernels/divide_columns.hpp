#ifndef STAN_MATH_OPENCL_KERNELS_DIVIDE_COLUMNS_HPP
#define STAN_MATH_OPENCL_KERNELS_DIVIDE_COLUMNS_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {
// \cond
static const std::string divide_columns_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Takes vector A and divides columns vector in A element-wise by the values
     * in vec
     * @param[out] A  Matrix to be divided elementwise
     * @param[in] vec  Vector to divide A by
     * @param vec_size Size of elementwise divisor.
     * @note Code is a <code>const char*</code> held in
     * <code>divide_columns_vec_kernel_code.</code>
     */
    __kernel void divide_columns_vec(__global double *A, __global double *vec,
                                     int vec_size) {
      const int i = get_global_id(0);
      A[i] /= vec[i % vec_size];
    }
    // \cond
);
// \endcond
/**
 * See the docs for \link kernels/divide_columns.hpp divide_columns_vec()
 * \endlink
 */
const kernel_cl<out_buffer, in_buffer, int> divide_columns_vec(
    "divide_columns_vec", {indexing_helpers, divide_columns_kernel_code});

// \cond
static const std::string divide_column_scalar_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Performs element-wise division on \c A
     * @param[out] A  Matrix to be divided elementwise
     * @param divisor  element to divide A by elementwise
     * @note Code is a <code>const char*</code> held in
     * <code>divide_column_scalar_kernel_code.</code>
     */
    __kernel void divide_columns_scalar(__global double *A, double divisor) {
      const int i = get_global_id(0);
      A[i] /= divisor;
    }
    // \cond
);
// \endcond
/**
 * See the docs for \link kernels/add.hpp add() \endlink
 */
const kernel_cl<out_buffer, double> divide_columns_scalar(
    "divide_columns_scalar",
    {indexing_helpers, divide_column_scalar_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
