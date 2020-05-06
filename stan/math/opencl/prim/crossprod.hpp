#ifndef STAN_MATH_OPENCL_PRIM_FUN_CROSSPROD_HPP
#define STAN_MATH_OPENCL_PRIM_FUN_CROSSPROD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {
/**
 * Returns the result of pre-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T type of elements in A
 * @param A input matrix
 * @return transpose(A) * A
 */
template <typename T_A,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_A>>
inline auto crossprod(T_A&& A) {
  const T_A& A_ref = A;
  return transpose(A_ref) * A_ref;
}
}  // namespace math
}  // namespace stan
#endif
#endif
