#ifndef STAN_MATH_PRIM_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP
#define STAN_MATH_PRIM_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {
/**
 * Compute trace(D B^T A B).
 **/
template <typename TD, int RD, int CD, typename TA, int RA, int CA, typename TB,
          int RB, int CB, typename = require_any_not_var_t<TD, TA, TB>>
inline return_type_t<TD, TA, TB> trace_gen_quad_form(
    const Eigen::Matrix<TD, RD, CD> &D, const Eigen::Matrix<TA, RA, CA> &A,
    const Eigen::Matrix<TB, RB, CB> &B) {
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);
  return (D * B.transpose() * A * B).trace();
}

}  // namespace math
}  // namespace stan
#endif
