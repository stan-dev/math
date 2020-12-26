#ifndef STAN_MATH_OPENCL_PRIM_DIAG_POST_MULTIPLY_HPP
#define STAN_MATH_OPENCL_PRIM_DIAG_POST_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/prim/diag_matrix.hpp>

namespace stan {
namespace math {

/**
 * Return the product of the matrix and the diagonal matrix
 * formed from the vector or row_vector.
 *
 * @tparam T1 type of input kernel generator expression for the
 * matrix
 * @tparam T2 type of input kernel generator expression for the
 * vector/row_vector
 * @param m1 input kernel generator expression for the matrix
 * @param m2 input kernel generator expression for the vector/row_vector
 *
 * @return product of the matrix and the diagonal matrix formed from the
 * vector or row_vector.
 */
template <
    typename T1, typename T2,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T1, T2>* = nullptr>
auto diag_post_multiply(const T1& m1, const T2& m2) {
  check_vector("diag_post_multiply (OpenCL)", "m2", m2);
  check_size_match("diag_post_multiply (OpenCL)", "m2.size()", m2.size(),
                   "m1.cols()", m1.cols());
  return m1 * diag_matrix(m2);
}

}  // namespace math
}  // namespace stan
#endif
