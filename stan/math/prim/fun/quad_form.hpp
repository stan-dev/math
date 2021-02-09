#ifndef STAN_MATH_PRIM_FUN_QUAD_FORM_HPP
#define STAN_MATH_PRIM_FUN_QUAD_FORM_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * @tparam EigMat1 type of the first (square) matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param A square matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix.
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_not_eigen_col_vector_t<EigMat2>* = nullptr,
          require_vt_same<EigMat1, EigMat2>* = nullptr,
          require_all_vt_arithmetic<EigMat1, EigMat2>* = nullptr>
inline auto quad_form(const EigMat1& A, const EigMat2& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);
  return make_holder(
      [](const auto& b, const auto& a) { return b.transpose() * a * b; },
      to_ref(B), to_ref(A));
}

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * @tparam EigMat type of the matrix
 * @tparam ColVec type of the vector
 *
 * @param A square matrix
 * @param B vector
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename EigMat, typename ColVec, require_eigen_t<EigMat>* = nullptr,
          require_eigen_col_vector_t<ColVec>* = nullptr,
          require_vt_same<EigMat, ColVec>* = nullptr,
          require_all_vt_arithmetic<EigMat, ColVec>* = nullptr>
inline value_type_t<EigMat> quad_form(const EigMat& A, const ColVec& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);
  const auto& B_ref = to_ref(B);
  return B_ref.dot(A * B_ref);
}

}  // namespace math
}  // namespace stan

#endif
