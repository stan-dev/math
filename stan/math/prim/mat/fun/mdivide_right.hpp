#ifndef STAN_MATH_PRIM_MAT_FUN_MDIVIDE_RIGHT_HPP
#define STAN_MATH_PRIM_MAT_FUN_MDIVIDE_RIGHT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/promote_common.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system xA=b.
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
 * @param A Matrix.
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename T1, typename T2, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<return_type_t<T1, T2>, R1, C2> mdivide_right(
    const Eigen::Matrix<T1, R1, C1> &b, const Eigen::Matrix<T2, R2, C2> &A) {
  check_square("mdivide_right", "A", A);
  check_multiplicable("mdivide_right", "b", b, "A", A);
  return promote_common<Eigen::Matrix<T1, R2, C2>, Eigen::Matrix<T2, R2, C2> >(
             A)
      .transpose()
      .lu()
      .solve(
          promote_common<Eigen::Matrix<T1, R1, C1>, Eigen::Matrix<T2, R1, C1> >(
              b)
              .transpose())
      .transpose();
}

}  // namespace math
}  // namespace stan

#endif
