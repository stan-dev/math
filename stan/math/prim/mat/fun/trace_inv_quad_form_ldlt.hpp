#ifndef STAN_MATH_PRIM_MAT_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_PRIM_MAT_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/LDLT_factor.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/fun/mdivide_left_ldlt.hpp>
#include <stan/math/prim/mat/fun/trace.hpp>
#include <stan/math/prim/mat/fun/transpose.hpp>
#include <stan/math/prim/mat/fun/multiply.hpp>

namespace stan {
namespace math {

/*
 * Compute the trace of an inverse quadratic form.  I.E., this computes
 *       trace(B^T A^-1 B)
 * where the LDLT_factor of A is provided.
 */
template <typename T1, typename T2, int R2, int C2, int R3, int C3,
          typename = require_any_not_var_t<T1, T2>>
inline return_type_t<T1, T2> trace_inv_quad_form_ldlt(
    const LDLT_factor<T1, R2, C2> &A, const Eigen::Matrix<T2, R3, C3> &B) {
  check_multiplicable("trace_inv_quad_form_ldlt", "A", A, "B", B);

  return trace(multiply(transpose(B), mdivide_left_ldlt(A, B)));
}

}  // namespace math
}  // namespace stan
#endif
