#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_MATRIX_VECTOR_MULTIPLY_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_MATRIX_VECTOR_MULTIPLY_HPP

#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <stan/math/opencl/kernel_generator/binary_operation.hpp>
#include <stan/math/opencl/kernel_generator/transpose.hpp>
#include <stan/math/opencl/kernel_generator/rowwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/broadcast.hpp>

namespace stan {
namespace math {

template <typename T_matrix, typename T_vector,
          typename = require_all_valid_expressions_t<T_matrix, T_vector>>
inline auto matrix_vector_multiply(T_matrix&& matrix, T_vector&& vector) {
  return rowwise_sum(elewise_multiplication(
      std::forward<T_matrix>(matrix),
      colwise_broadcast(transpose(std::forward<T_vector>(vector)))));
}

}  // namespace math
}  // namespace stan

#endif
