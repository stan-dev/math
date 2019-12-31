#ifndef STAN_MATH_PRIM_MAT_FUN_MDIVIDE_RIGHT_SPD_HPP
#define STAN_MATH_PRIM_MAT_FUN_MDIVIDE_RIGHT_SPD_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/mdivide_left_spd.hpp>
#include <stan/math/prim/mat/fun/transpose.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system xA=b where A is symmetric
 * positive definite.
 *
 * @tparam T1 type of elements in the right-hand side matrix or vector
 * @tparam T2 type of elements in the second matrix
 * @tparam R1 number of rows in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam C1 number of columns in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam R2 number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam C2 number of columns in the second matrix, can be Eigen::Dynamic
 *
 * @param b right-hand side matrix or vector
 * @param A matrix
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename T1, typename T2, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<return_type_t<T1, T2>, R1, C2> mdivide_right_spd(
    const Eigen::Matrix<T1, R1, C1> &b, const Eigen::Matrix<T2, R2, C2> &A) {
  check_multiplicable("mdivide_right_spd", "b", b, "A", A);
  check_pos_definite("mdivide_right_spd", "A", A);
  // FIXME: After allowing for general MatrixBase in mdivide_left_spd,
  //        change to b.transpose()
  return mdivide_left_spd(A, transpose(b)).transpose();
}

}  // namespace math
}  // namespace stan

#endif
