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
 * @tparam EigMat type of the right hand side
 * @tparam T type of matrix in the LDLT_factor
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <typename EigMat, typename T,
          require_all_matrix_t<EigMat, T>* = nullptr,
          require_any_not_st_arithmetic<EigMat, T>* = nullptr>
inline auto mdivide_right_ldlt(const EigMat& b, LDLT_factor<T>& A) {
  check_multiplicable("mdivide_right_ldlt", "b", b, "A", A.matrix());
  check_ldlt_factor("mdivide_right_ldlt", "A", A);

  return mdivide_left_ldlt(A, b.transpose()).transpose().eval();
}

/**
 * Returns the solution of the system xA=b given an LDLT_factor of A
 *
 * Overload for arithmetic types
 *
 * @tparam EigMat type of the right hand side
 * @tparam T type of matrix in the LDLT_factor
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <typename EigMat, typename T,
          require_all_matrix_t<EigMat, T>* = nullptr,
          require_all_st_arithmetic<EigMat, T>* = nullptr>
inline Eigen::Matrix<double, EigMat::RowsAtCompileTime, T::ColsAtCompileTime>
mdivide_right_ldlt(const EigMat& b, LDLT_factor<T>& A) {
  check_multiplicable("mdivide_right_ldlt", "b", b, "A", A.matrix());
  check_ldlt_factor("mdivide_right_ldlt", "A", A);

  if (A.matrix().rows() == 0) {
    return {b.rows(), 0};
  }

  return A.ldlt().solve(b.transpose()).transpose();
}

}  // namespace math
}  // namespace stan

#endif
