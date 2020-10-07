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

/**
 * Compute the trace of an inverse quadratic form.  I.E., this computes
 *       trace(B^T A^-1 B)
 * where the LDLT_factor of A is provided.
 *
 * @tparam T type of elements in the first matrix and the LDLT_factor
 * @tparam EigMat2 type of the second matrix
 * @tparam R number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C number of columns in the LDLT_factor, can be Eigen::Dynamic
 *
 */
template <typename T, typename EigMat2, int R, int C,
          typename = require_all_not_st_var<T, EigMat2>>
inline return_type_t<T, EigMat2> trace_inv_quad_form_ldlt(
    const LDLT_factor<T, R, C>& A, const EigMat2& B) {
  check_multiplicable("trace_inv_quad_form_ldlt", "A", A, "B", B);
  if (A.cols() == 0) {
    return 0;
  }

  return B.cwiseProduct(mdivide_left_ldlt(A, B)).sum();
}

}  // namespace math
}  // namespace stan

#endif
