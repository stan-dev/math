#ifndef STAN_MATH_PRIM_MAT_FUN_TRACE_QUAD_FORM_HPP
#define STAN_MATH_PRIM_MAT_FUN_TRACE_QUAD_FORM_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {
/**
 * Compute trace(B^T A B).
 **/
template <typename T1, typename T2, typename = enable_if_all_eigen<T1, T2>,
          typename = enable_if_all_scalar_arithmetic<T1, T2>>
inline double trace_quad_form(const T1 &A, const T2 &B) {
  check_square("trace_quad_form", "A", A);
  check_multiplicable("trace_quad_form", "A", A, "B", B);
  return (B.transpose() * A * B).trace();
}

}  // namespace math
}  // namespace stan
#endif
