#ifndef STAN_MATH_OPENCL_PRIM_FUN_TCROSSPROD_HPP
#define STAN_MATH_OPENCL_PRIM_FUN_TCROSSPROD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {
/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 * 
 * @tparam T type of elements in A
 * @param A input matrix
 * @return A * transpose(A)
 */
template <typename T_A,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_A>>
inline auto tcrossprod(T_A&& A) {
  const matrix_cl<double>& A_ref = A;
  return A_ref * transpose(A_ref);
}
}  // namespace math
}  // namespace stan
#endif
#endif
