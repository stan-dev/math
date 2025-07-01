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
inline plain_type_t<EigMat2> quad_form_sym(EigMat1&& A, EigMat2&& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  auto&& A_ref = to_ref(std::forward<EigMat1>(A));
  check_symmetric("quad_form_sym", "A", A_ref);
  return make_holder([](auto&& a, auto&& b) {
    auto ret = (a.transpose() * a * b).eval();
    return 0.5 * (ret + ret.transpose());
  }, std::forward<decltype(A_ref)>(A_ref), to_ref(std::forward<EigMat2>(B)));
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
inline value_type_t<EigMat> quad_form_sym(EigMat&& A, ColVec&& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  auto&& A_ref = to_ref(std::forward<EigMat>(A));
  check_symmetric("quad_form_sym", "A", A_ref);
  return make_holder([](auto&& a, auto&& b) {
    return b.dot(a * b);
  }, std::forward<decltype(A_ref)>(A_ref), to_ref(std::forward<ColVec>(B)));
}

}  // namespace math
}  // namespace stan

#endif
