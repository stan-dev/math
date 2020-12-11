#ifndef STAN_MATH_OPENCL_REV_DIAGONAL_HPP
#define STAN_MATH_OPENCL_REV_DIAGONAL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/diag_matrix.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return a column vector of the diagonal elements of the
 * specified matrix.  The matrix is not required to be square.
 *
 * @param m Specified matrix.
 * @return Diagonal of the matrix.
 */
inline var_value<matrix_cl<double>> diagonal(
    const var_value<matrix_cl<double>>& M) {
  var_value<matrix_cl<double>> res = diagonal(M.val());

  reverse_pass_callback([M, res]() mutable {
    M.adj() = M.adj() + diag_matrix(res.adj());
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
#endif
