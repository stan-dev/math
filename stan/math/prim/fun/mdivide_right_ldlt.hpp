#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_LDLT_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_LDLT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>
#include <stan/math/prim/fun/mdivide_left_ldlt.hpp>
#include <stan/math/prim/fun/transpose.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system xA=b given an LDLT_factor of A
 *
 * @tparam T1 type of elements in right-hand side matrix or vector
 * @tparam T2 type of elements in the LDLT_factor
 * @tparam R1 number of rows in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam C1 number of columns in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam R2 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C2 number of columns in the LDLT_factor, can be Eigen::Dynamic
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */

template <typename T1, typename T2, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<return_type_t<T1, T2>, R1, C2> mdivide_right_ldlt(
    const Eigen::Matrix<T1, R1, C1> &b, const LDLT_factor<T2, R2, C2> &A) {
  if (b.cols() == 0 && A.rows() == 0) {
    return Eigen::Matrix<return_type_t<T1, T2>, R1, C2>(b.rows(), 0);
  }

  check_multiplicable("mdivide_right_ldlt", "b", b, "A", A);

  return transpose(mdivide_left_ldlt(A, transpose(b)));
}

template <int R1, int C1, int R2, int C2>
inline Eigen::Matrix<double, R1, C2> mdivide_right_ldlt(
    const Eigen::Matrix<double, R1, C1> &b,
    const LDLT_factor<double, R2, C2> &A) {
  if (b.cols() == 0 && A.rows() == 0) {
    return Eigen::Matrix<double, R1, C2>(b.rows(), 0);
  }

  check_multiplicable("mdivide_right_ldlt", "b", b, "A", A);

  return A.solveRight(b);
}

}  // namespace math
}  // namespace stan

#endif
