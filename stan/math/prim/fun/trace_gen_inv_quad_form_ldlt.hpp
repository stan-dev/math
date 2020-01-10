#ifndef STAN_MATH_PRIM_FUN_TRACE_GEN_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_PRIM_FUN_TRACE_GEN_INV_QUAD_FORM_LDLT_HPP

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
 *       trace(D B^T A^-1 B)
 * where D is a square matrix and the LDLT_factor of A is provided.
 *
 * @tparam T1 type of elements in the first matrix
 * @tparam T2 type of elements in the LDLT_factor
 * @tparam T3 type of elements in the third matrix
 * @tparam R1 number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam C1 number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam R2 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C2 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam R3 number of rows in the third matrix, can be Eigen::Dynamic
 * @tparam C3 number of columns in the third matrix, can be Eigen::Dynamic
 *
 * @param D multiplier
 * @param A LDLT_factor
 * @param B inner term in quadratic form
 * @return trace(D * B^T * A^-1 * B)
 * @throw std::domain_error if D is not square
 * @throw std::domain_error if A cannot be multiplied by B or B cannot
 * be multiplied by D.
 */
template <typename T1, typename T2, typename T3, int R1, int C1, int R2, int C2,
          int R3, int C3, typename = require_all_not_var_t<T1, T2, T3>>
inline return_type_t<T1, T2, T3> trace_gen_inv_quad_form_ldlt(
    const Eigen::Matrix<T1, R1, C1> &D, const LDLT_factor<T2, R2, C2> &A,
    const Eigen::Matrix<T3, R3, C3> &B) {
  check_square("trace_gen_inv_quad_form_ldlt", "D", D);
  if (D.size() == 0 && A.cols() == 0 && B.rows() == 0 && B.cols() == 0) {
    return 0;
  }

  check_multiplicable("trace_gen_inv_quad_form_ldlt", "A", A, "B", B);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "B", B, "D", D);

  return trace(multiply(multiply(D, transpose(B)), mdivide_left_ldlt(A, B)));
}

}  // namespace math
}  // namespace stan

#endif
