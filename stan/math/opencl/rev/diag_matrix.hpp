#ifndef STAN_MATH_OPENCL_REV_DIAG_MATRIX_HPP
#define STAN_MATH_OPENCL_REV_DIAG_MATRIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/diag_matrix.hpp>
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
inline var_value<matrix_cl<double>> diag_matrix(
    const var_value<matrix_cl<double>>& M) {
  var_value<matrix_cl<double>> res = diag_matrix(M.val());

  reverse_pass_callback([M, res]() mutable {
    M.adj() = M.adj() + diagonal(res.adj());
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
