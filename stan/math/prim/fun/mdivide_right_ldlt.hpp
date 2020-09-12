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
 * @tparam EigMat type of the matrix or vector
 * @tparam T type of elements in the LDLT_factor
 * @tparam R number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C number of columns in the LDLT_factor, can be Eigen::Dynamic
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <typename EigMat, typename T, int R, int C,
          require_eigen_t<EigMat>* = nullptr,
          require_any_not_arithmetic_t<value_type_t<EigMat>, T>* = nullptr>
inline Eigen::Matrix<return_type_t<EigMat, T>, EigMat::RowsAtCompileTime, C>
mdivide_right_ldlt(const EigMat& b, const LDLT_factor<T, R, C>& A) {
  check_multiplicable("mdivide_right_ldlt", "b", b, "A", A);
  if (A.rows() == 0) {
    return {b.rows(), 0};
  }

  return transpose(mdivide_left_ldlt(A, transpose(b)));
}

template <typename EigMat, typename T, int R, int C,
          require_eigen_t<EigMat>* = nullptr,
          require_all_arithmetic_t<value_type_t<EigMat>, T>* = nullptr>
inline Eigen::Matrix<T, EigMat::RowsAtCompileTime, C> mdivide_right_ldlt(
    const EigMat& b, const LDLT_factor<T, R, C>& A) {
  check_multiplicable("mdivide_right_ldlt", "b", b, "A", A);
  if (A.rows() == 0) {
    return {b.rows(), 0};
  }

  return A.solveRight(b);
}

}  // namespace math
}  // namespace stan

#endif
