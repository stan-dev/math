#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_LDLT_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_LDLT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam T type of elements in the LDLT_factor
 * @tparam R number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C number of columns in the LDLT_factor, can be Eigen::Dynamic
 * @tparam EigMat type of the matrix
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <int R, int C, typename T, typename EigMat,
          require_eigen_t<EigMat>* = nullptr,
          require_all_not_st_var<T, EigMat>* = nullptr,
          require_any_not_t<std::is_arithmetic<T>,
                            is_fvar<value_type_t<EigMat>>>* = nullptr>
inline Eigen::Matrix<return_type_t<T, EigMat>, R, EigMat::ColsAtCompileTime>
mdivide_left_ldlt(const LDLT_factor<T, R, C>& A, const EigMat& b) {
  check_multiplicable("mdivide_left_ldlt", "A", A, "b", b);
  if (A.cols() == 0) {
    return {0, b.cols()};
  }

  return A.solve(
      Eigen::Matrix<return_type_t<T, EigMat>, EigMat::RowsAtCompileTime,
                    EigMat::ColsAtCompileTime>(b));
}

}  // namespace math
}  // namespace stan

#endif
