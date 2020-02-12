#ifndef STAN_MATH_PRIM_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_PRIM_FUN_QUAD_FORM_SYM_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * Symmetry of the resulting matrix is guaranteed.
 *
 * @tparam RA number of rows in the symmetric matrix, can be Eigen::Dynamic
 * @tparam CA number of columns in the symmetric matrix, can be Eigen::Dynamic
 * @tparam RB number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam CB number of columns in the second matrix, can be Eigen::Dynamic
 * @tparam TA type of elements
 *
 * @param A symmetric matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix of size CB.
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <int RA, int CA, int RB, int CB, typename T>
inline Eigen::Matrix<T, CB, CB> quad_form_sym(
    const Eigen::Matrix<T, RA, CA>& A, const Eigen::Matrix<T, RB, CB>& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  check_symmetric("quad_form_sym", "A", A);
  Eigen::Matrix<T, CB, CB> ret(B.transpose() * A * B);
  return T(0.5) * (ret + ret.transpose());
}

/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * @tparam RA number of rows in the symmetric matrix, can be Eigen::Dynamic
 * @tparam CA number of columns in the symmetric matrix, can be Eigen::Dynamic
 * @tparam RB number of rows in the vector, can be Eigen::Dynamic
 * @tparam T type of elements
 *
 * @param A symmetric matrix
 * @param B vector
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <int RA, int CA, int RB, typename T>
inline T quad_form_sym(const Eigen::Matrix<T, RA, CA>& A,
                       const Eigen::Matrix<T, RB, 1>& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  check_symmetric("quad_form_sym", "A", A);
  return B.dot(A * B);
}

}  // namespace math
}  // namespace stan

#endif
