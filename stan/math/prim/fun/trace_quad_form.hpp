#ifndef STAN_MATH_PRIM_FUN_TRACE_QUAD_FORM_HPP
#define STAN_MATH_PRIM_FUN_TRACE_QUAD_FORM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Compute trace(B^T A B).
 *
 * @tparam EigMat1 type of the first matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param A matrix
 * @param B matrix
 * @return The trace of B^T A B
 * @throw std::domain_error if A is not square
 * @throw std::domain_error if A cannot be multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_vt<std::is_arithmetic, EigMat1, EigMat2>* = nullptr>
inline auto trace_quad_form(EigMat1&& A, EigMat2&& B) {
  check_square("trace_quad_form", "A", A);
  check_multiplicable("trace_quad_form", "A", A, "B", B);
  decltype(auto) B_ref = to_ref(std::forward<EigMat2>(B));
  return make_holder(
      [](auto&& A_, auto&& B_) { return B_.cwiseProduct(A_ * B_).sum(); },
      std::forward<EigMat1>(A), std::forward<decltype(B_ref)>(B_ref));
}

}  // namespace math
}  // namespace stan

#endif
