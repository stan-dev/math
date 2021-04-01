#ifndef STAN_MATH_OPENCL_REV_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_REV_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/crossprod.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Transposes a matrix.
 * @param M input matrix
 * @return transposed matrix
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline auto transpose(const var_value<T>& M) {
  return M.transpose();
}

}  // namespace math
}  // namespace stan

#endif
#endif
