#ifndef STAN_MATH_OPENCL_PRIM_DIAG_PRE_MULTIPLY_HPP
#define STAN_MATH_OPENCL_PRIM_DIAG_PRE_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/prim/diag_matrix.hpp>

namespace stan {
namespace math {

/**
 * Return the product of the diagonal matrix formed from the
 * vector or row_vector and a matrix.
 *
 * @tparam T1 type of input kernel generator expression for the
 * vector/row_vector
 * @tparam T2 type of input kernel generator expression for the
 * matrix
 * @param m1 input kernel generator expression for the vector/row_vector
 * @param m2 input kernel generator expression for the matrix
 *
 * @return product of the diagonal matrix formed from the vector or row_vector
 * and a matrix.
 */
template <
    typename T1, typename T2,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T1, T2>* = nullptr>
auto diag_pre_multiply(const T1& m1, const T2& m2) {
  check_vector("diag_pre_multiply (OpenCL)", "m1", m1);
  check_size_match("diag_pre_multiply (OpenCL)", "m1.size()", m1.size(),
                   "m2.rows()", m2.rows());
  return diag_matrix(m1) * m2;
}

}  // namespace math
}  // namespace stan
#endif
