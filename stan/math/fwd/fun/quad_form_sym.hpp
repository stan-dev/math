#ifndef STAN_MATH_FWD_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_FWD_FUN_QUAD_FORM_SYM_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

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
 * @return The quadratic form, which is a symmetric matrix of size CB.
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_not_eigen_col_vector_t<EigMat2>* = nullptr,
          require_any_vt_fvar<EigMat1, EigMat2>* = nullptr>
inline promote_scalar_t<return_type_t<EigMat1, EigMat2>, EigMat2> quad_form_sym(
    const EigMat1& A, const EigMat2& B) {
  using T_ret = return_type_t<EigMat1, EigMat2>;
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  check_symmetric("quad_form_sym", "A", A);
  const auto& B_ref = to_ref(B);
  promote_scalar_t<T_ret, EigMat2> ret(
      multiply(B_ref.transpose(), multiply(A, B_ref)));
  return T_ret(0.5) * (ret + ret.transpose());
}

/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * @tparam EigMat type of the (symmetric) matrix
 * @tparam ColVec type of the vector
 *
 * @param A symmetric matrix
 * @param B vector
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename EigMat, typename ColVec, require_eigen_t<EigMat>* = nullptr,
          require_eigen_col_vector_t<ColVec>* = nullptr,
          require_any_vt_fvar<EigMat, ColVec>* = nullptr>
inline return_type_t<EigMat, ColVec> quad_form_sym(const EigMat& A,
                                                   const ColVec& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  check_symmetric("quad_form_sym", "A", A);
  const auto& B_ref = to_ref(B);
  return dot_product(B_ref, multiply(A, B_ref));
}

}  // namespace math
}  // namespace stan

#endif
