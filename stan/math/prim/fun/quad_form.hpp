#ifndef STAN_MATH_PRIM_FUN_QUAD_FORM_HPP
#define STAN_MATH_PRIM_FUN_QUAD_FORM_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * @tparam RA number of rows in the square matrix, can be Eigen::Dynamic
 * @tparam CA number of columns in the square matrix, can be Eigen::Dynamic
 * @tparam RB number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam CB number of columns in the second matrix, can be Eigen::Dynamic
 * @tparam T type of elements
 *
 * @param A square matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix of size CB.
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <int RA, int CA, int RB, int CB, typename T>
inline Eigen::Matrix<T, CB, CB> quad_form(const Eigen::Matrix<T, RA, CA>& A,
                                          const Eigen::Matrix<T, RB, CB>& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);
  return B.transpose() * A * B;
}

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * @tparam RA number of rows in the square matrix, can be Eigen::Dynamic
 * @tparam CA number of columns in the square matrix, can be Eigen::Dynamic
 * @tparam RB number of rows in the vector, can be Eigen::Dynamic
 * @tparam T type of elements
 *
 * @param A square matrix
 * @param B vector
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <int RA, int CA, int RB, typename T>
inline T quad_form(const Eigen::Matrix<T, RA, CA>& A,
                   const Eigen::Matrix<T, RB, 1>& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);
  return B.dot(A * B);
}

}  // namespace math
}  // namespace stan

#endif
