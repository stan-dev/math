#ifndef STAN_MATH_OPENCL_PRIM_IDENTITY_MATRIX_HPP
#define STAN_MATH_OPENCL_PRIM_IDENTITY_MATRIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return a square identity matrix
 *
 * @param K size of the matrix
 * @return An identity matrix of size K.
 * @throw std::domain_error if K is negative.
 */
template <typename T_x, require_matrix_cl_t<T_x>* = nullptr>
inline auto identity_matrix(int K) {
  using T_val = value_type_t<T_x>;
  check_nonnegative("identity_matrix(OpenCL)", "size", K);
  return select(row_index(K, K) == col_index(), T_val(1), T_val(0));
}
}  // namespace math
}  // namespace stan
#endif
#endif
