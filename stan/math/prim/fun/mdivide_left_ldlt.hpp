#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_LDLT_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_LDLT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>
#include <stan/math/prim/fun/promote_common.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam T1 type of elements in the LDLT_factor
 * @tparam T2 type of elements in right-hand side matrix or vector
 * @tparam R1 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C1 number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam R2 number of rows in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam C2 number of columns in the right-hand side matrix, can be
 *         Eigen::Dynamic
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <int R1, int C1, int R2, int C2, typename T1, typename T2>
inline Eigen::Matrix<return_type_t<T1, T2>, R1, C2> mdivide_left_ldlt(
    const LDLT_factor<T1, R1, C1> &A, const Eigen::Matrix<T2, R2, C2> &b) {
  if (A.cols() == 0 && b.rows() == 0) {
    return Eigen::Matrix<return_type_t<T1, T2>, R1, C2>(0, b.cols());
  }

  check_multiplicable("mdivide_left_ldlt", "A", A, "b", b);

  return A.solve(
      promote_common<Eigen::Matrix<T1, R2, C2>, Eigen::Matrix<T2, R2, C2> >(b));
}

}  // namespace math
}  // namespace stan

#endif
