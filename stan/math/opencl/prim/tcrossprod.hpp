#ifndef STAN_MATH_OPENCL_PRIM_FUN_TCROSSPROD_HPP
#define STAN_MATH_OPENCL_PRIM_FUN_TCROSSPROD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
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
inline matrix_cl<typename std::decay_t<T_A>::Scalar> tcrossprod(T_A&& A) {
  const matrix_cl<typename std::decay_t<T_A>::Scalar>& A_eval
      = std::forward<T_A>(A);
  return multiply_transpose(A_eval);
}
}  // namespace math
}  // namespace stan
#endif
#endif
