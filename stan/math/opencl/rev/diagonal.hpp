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
 * @param M Specified matrix.
 * @return Diagonal of the matrix.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> diagonal(const var_value<T>& M) {
  return make_callback_var(diagonal(M.val()),
                           [M](vari_value<matrix_cl<double>>& res) mutable {
                             diagonal(M.adj()) += res.adj();
                           });
}

}  // namespace math
}  // namespace stan

#endif
#endif
