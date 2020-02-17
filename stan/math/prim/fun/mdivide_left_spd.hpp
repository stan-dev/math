#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_SPD_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_SPD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b where A is symmetric positive
 * definite.
 *
 * @tparam T1 type of elements in the first matrix
 * @tparam T2 type of elements in the right-hand side matrix or vector
 * @tparam R1 number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam C1 number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam R2 number of rows in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam C2 number of columns in the right-hand side matrix, can be
 *         Eigen::Dynamic
 *
 * @param A Matrix.
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename T1, typename T2, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<return_type_t<T1, T2>, R1, C2> mdivide_left_spd(
    const Eigen::Matrix<T1, R1, C1> &A, const Eigen::Matrix<T2, R2, C2> &b) {
  static const char *function = "mdivide_left_spd";
  check_multiplicable(function, "A", A, "b", b);
  check_positive(function, "rows", A.rows());
  check_symmetric(function, "A", A);
  check_not_nan(function, "A", A);

  auto llt = Eigen::Matrix<return_type_t<T1, T2>, R1, C1>(A).llt();
  check_pos_definite(function, "A", llt);
  return llt.solve(Eigen::Matrix<return_type_t<T1, T2>, R2, C2>(b));
}

}  // namespace math
}  // namespace stan

#endif
