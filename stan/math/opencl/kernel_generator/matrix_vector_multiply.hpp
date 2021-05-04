#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_MATRIX_VECTOR_MULTIPLY_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_MATRIX_VECTOR_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/opencl/kernel_generator/binary_operation.hpp>
#include <stan/math/opencl/kernel_generator/transpose.hpp>
#include <stan/math/opencl/kernel_generator/rowwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/broadcast.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Multiplies a matrix and a vector on an OpenCL device.
 * @tparam T_matrix type of the matrix or matrix expression
 * @tparam T_vector type of the vector or vector expression
 * @param matrix the matrix
 * @param vector the vector
 * @return matrix-vector multiplication expression
 * @throw invalid_argument vector has more than column or matrix has different
 * number of rows than vector columns.
 */
template <typename T_matrix, typename T_vector,
          typename = require_all_kernel_expressions_t<T_matrix, T_vector>>
inline auto matrix_vector_multiply(T_matrix&& matrix, T_vector&& vector) {
  return rowwise_sum(elt_multiply(
      std::forward<T_matrix>(matrix),
      colwise_broadcast(transpose(std::forward<T_vector>(vector)))));
}

}  // namespace math
}  // namespace stan

#endif
#endif
