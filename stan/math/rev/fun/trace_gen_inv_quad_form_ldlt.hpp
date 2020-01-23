#ifndef STAN_MATH_REV_FUN_TRACE_GEN_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_REV_FUN_TRACE_GEN_INV_QUAD_FORM_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/trace_inv_quad_form_ldlt.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Compute the trace of an inverse quadratic form premultiplied by a
 * square matrix. This computes
 *       trace(D B^T A^-1 B)
 * where D is a square matrix and the LDLT_factor of A is provided.
 *
 * @tparam T1 type of elements in the first matrix
 * @tparam R1 number of rows, can be Eigen::Dynamic
 * @tparam C1 number of columns, can be Eigen::Dynamic
 * @tparam T2 type of elements in the LDLT_factor
 * @tparam R2 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C2 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam T3 type of elements in the second matrix
 * @tparam R3 number of rows, can be Eigen::Dynamic
 * @tparam C3 number of columns, can be Eigen::Dynamic
 *
 * @param D a square matrix
 * @param A an LDLT_factor
 * @param B a matrix
 * @return The trace of the inverse quadratic form.
 */
template <typename T1, int R1, int C1, typename T2, int R2, int C2, typename T3,
          int R3, int C3, require_any_var_t<T1, T2, T3>...>
inline var trace_gen_inv_quad_form_ldlt(const Eigen::Matrix<T1, R1, C1> &D,
                                        const LDLT_factor<T2, R2, C2> &A,
                                        const Eigen::Matrix<T3, R3, C3> &B) {
  check_square("trace_gen_inv_quad_form_ldlt", "D", D);
  if (D.size() == 0 && A.cols() == 0 && B.rows() == 0 && B.cols() == 0) {
    return 0;
  }

  check_multiplicable("trace_gen_inv_quad_form_ldlt", "A", A, "B", B);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "B", B, "D", D);

  internal::trace_inv_quad_form_ldlt_impl<T2, R2, C2, T3, R3, C3> *_impl
      = new internal::trace_inv_quad_form_ldlt_impl<T2, R2, C2, T3, R3, C3>(
          D, A, B);

  return var(
      new internal::trace_inv_quad_form_ldlt_vari<T2, R2, C2, T3, R3, C3>(
          _impl));
}

}  // namespace math
}  // namespace stan
#endif
