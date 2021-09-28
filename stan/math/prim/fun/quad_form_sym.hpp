#ifndef STAN_MATH_PRIM_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_PRIM_FUN_QUAD_FORM_SYM_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
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
 * @return The quadratic form, which is a symmetric matrix.
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_not_eigen_col_vector_t<EigMat2>* = nullptr,
          require_vt_same<EigMat1, EigMat2>* = nullptr,
          require_all_vt_arithmetic<EigMat1, EigMat2>* = nullptr>
inline plain_type_t<EigMat2> quad_form_sym(const EigMat1& A, const EigMat2& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  const auto& A_ref = to_ref(A);
  const auto& B_ref = to_ref(B);
  check_symmetric("quad_form_sym", "A", A_ref);
  return make_holder(
      [](const auto& ret) { return 0.5 * (ret + ret.transpose()); },
      (B_ref.transpose() * A_ref * B_ref).eval());
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
          require_vt_same<EigMat, ColVec>* = nullptr,
          require_all_vt_arithmetic<EigMat, ColVec>* = nullptr>
inline value_type_t<EigMat> quad_form_sym(const EigMat& A, const ColVec& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  const auto& A_ref = to_ref(A);
  const auto& B_ref = to_ref(B);
  check_symmetric("quad_form_sym", "A", A_ref);
  return B_ref.dot(A_ref * B_ref);
}

}  // namespace math
}  // namespace stan

#endif
