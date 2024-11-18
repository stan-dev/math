#ifndef STAN_MATH_REV_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_REV_FUN_QUAD_FORM_SYM_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/quad_form.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/quad_form_sym.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * Symmetry of the resulting matrix is guaranteed.
 *
 * @tparam EigMat1 type of the first (symmetric) matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param A symmetric matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix of size Cb.
 * If \c B is a column vector returns a scalar.
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_any_vt_var<EigMat1, EigMat2>* = nullptr>
inline auto quad_form_sym(const EigMat1& A, const EigMat2& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  const auto& A_ref = to_ref(A);
  check_symmetric("quad_form_sym", "A", A_ref);
  return quad_form(A_ref, B, true);
}

}  // namespace math
}  // namespace stan
#endif
