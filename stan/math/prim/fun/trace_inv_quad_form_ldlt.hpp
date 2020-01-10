#ifndef STAN_MATH_PRIM_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_PRIM_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>
#include <stan/math/prim/fun/mdivide_left_ldlt.hpp>
#include <stan/math/prim/fun/trace.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/fun/multiply.hpp>

namespace stan {
namespace math {

/*
 * Compute the trace of an inverse quadratic form.  I.E., this computes
 *       trace(B^T A^-1 B)
 * where the LDLT_factor of A is provided.
 *
 * @tparam T1 type of elements in the first matrix and the LDLT_factor
 * @tparam T2 type of elements in the second matrix
 * @tparam R1 number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam C1 number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam R2 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C2 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam R3 number of rows in the third matrix, can be Eigen::Dynamic
 * @tparam C3 number of columns in the third matrix, can be Eigen::Dynamic
 *
 */
template <typename T1, typename T2, int R2, int C2, int R3, int C3,
          typename = require_all_not_var_t<T1, T2>>
inline return_type_t<T1, T2> trace_inv_quad_form_ldlt(
    const LDLT_factor<T1, R2, C2> &A, const Eigen::Matrix<T2, R3, C3> &B) {
  if (A.rows() == 0 && B.rows() == 0) {
    return 0;
  }

  check_multiplicable("trace_inv_quad_form_ldlt", "A", A, "B", B);

  return trace(multiply(transpose(B), mdivide_left_ldlt(A, B)));
}

}  // namespace math
}  // namespace stan

#endif
