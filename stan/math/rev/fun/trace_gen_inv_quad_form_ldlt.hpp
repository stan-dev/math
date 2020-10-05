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
 * @tparam EigMat1 type of the first matrix
 * @tparam T2 type of elements in the LDLT_factor
 * @tparam R2 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C2 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam EigMat3 type of the second matrix
 *
 * @param D a square matrix
 * @param A an LDLT_factor
 * @param B a matrix
 * @return The trace of the inverse quadratic form.
 */
template <typename EigMat1, typename T2, int R2, int C2, typename EigMat3,
          require_all_eigen_t<EigMat1, EigMat3>* = nullptr,
          require_any_st_var<EigMat1, T2, EigMat3>* = nullptr>
inline var trace_gen_inv_quad_form_ldlt(const EigMat1& D,
                                        const LDLT_factor<T2, R2, C2>& A,
                                        const EigMat3& B) {
  using T3 = value_type_t<EigMat3>;
  constexpr int R3 = EigMat3::RowsAtCompileTime;
  constexpr int C3 = EigMat3::ColsAtCompileTime;
  check_square("trace_gen_inv_quad_form_ldlt", "D", D);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "A", A, "B", B);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "B", B, "D", D);
  if (D.size() == 0 || A.cols() == 0) {
    return 0;
  }

  auto* _impl
      = new internal::trace_inv_quad_form_ldlt_impl<T2, R2, C2, T3, R3, C3>(
          D, A, B);

  return var(
      new internal::trace_inv_quad_form_ldlt_vari<T2, R2, C2, T3, R3, C3>(
          _impl));
}

}  // namespace math
}  // namespace stan
#endif
