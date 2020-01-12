#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_TRI_LOW_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_TRI_LOW_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>

namespace stan {
namespace math {

/**
 * Return the result of left dividing the second argument by the
 * first argument.  Calling <code>mdivide_left_tri_low(A,
 * b)</code> with divisor <code>A</code> and dividend
 * <code>b</code> is more arithmetically stable than calling
 * <code>inv(A) * b</code>.
 *
 * @tparam T1 type of elements in the divisor matrix
 * @tparam T2 type of elements in the dividend matrix
 * @tparam R1 number rows in the divisor matrix, can be Eigen::Dynamic
 * @tparam C1 number columns in the divisor matrix, can be Eigen::Dynamic
 * @tparam R2 number rows in the dividend matrix, can be Eigen::Dynamic
 * @tparam C2 number columns in the dividend matrix, can be Eigen::Dynamic
 *
 * @param A divisor, an invertible square matrix
 * @param b dividend, a matrix or vector with the same number of
 *   rows as the divisor has columns
 * @return left division of b by A
 * @throws std::invalid_argument if the divisor is not square or
 *   the dividend does not have the same number of rows as the
 *   divisor has columns.
 */
template <typename T1, typename T2, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<return_type_t<T1, T2>, R1, C2> mdivide_left_tri_low(
    const Eigen::Matrix<T1, R1, C1> &A, const Eigen::Matrix<T2, R2, C2> &b) {
  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  return mdivide_left_tri<Eigen::Lower>(A, b);
}

template <typename T, int R1, int C1>
inline Eigen::Matrix<T, R1, C1> mdivide_left_tri_low(
    const Eigen::Matrix<T, R1, C1> &A) {
  check_square("mdivide_left_tri_low", "A", A);
  return mdivide_left_tri<Eigen::Lower>(A);
}

}  // namespace math
}  // namespace stan

#endif
