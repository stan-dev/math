#ifndef STAN_MATH_OPENCL_REV_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_REV_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/crossprod.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of pre-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T Type of the matrix
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
inline var_value<matrix_cl<double>> transpose(
    const var_value<matrix_cl<double>>& M) {
  var_value<matrix_cl<double>> res = transpose(M.val());

  reverse_pass_callback([M, res]() mutable {
    M.adj() = transpose(res.adj());
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
