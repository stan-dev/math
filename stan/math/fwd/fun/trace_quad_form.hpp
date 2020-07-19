#ifndef STAN_MATH_FWD_FUN_TRACE_QUAD_FORM_HPP
#define STAN_MATH_FWD_FUN_TRACE_QUAD_FORM_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/fun/trace.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_any_vt_fvar<EigMat1, EigMat2>* = nullptr>
inline return_type_t<EigMat1, EigMat2> trace_quad_form(const EigMat1& A,
                                                       const EigMat2& B) {
  check_square("trace_quad_form", "A", A);
  check_multiplicable("trace_quad_form", "A", A, "B", B);
  const auto& B_ref = to_ref(B);
  return B_ref.cwiseProduct(multiply(A, B_ref)).sum();
}

}  // namespace math
}  // namespace stan

#endif
