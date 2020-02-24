#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_TRI_LOW_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_TRI_LOW_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_right_tri.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system x tri(A) = b when tri(A) is a
 * lower triangular view of the matrix A.
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
 * @return x = b * tri(A)^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename T1, typename T2, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<return_type_t<T1, T2>, R1, C2> mdivide_right_tri_low(
    const Eigen::Matrix<T1, R1, C1> &b, const Eigen::Matrix<T2, R2, C2> &A) {
  return mdivide_right_tri<Eigen::Lower>(
      Eigen::Matrix<return_type_t<T1, T2>, R1, C1>(b),
      Eigen::Matrix<return_type_t<T1, T2>, R2, C2>(A));
}

}  // namespace math
}  // namespace stan

#endif
