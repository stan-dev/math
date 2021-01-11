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
 * @tparam T type of matrix in the LDLT_factor
 * @tparam EigMat type of the right hand side
 *
 * @param A LDLT_factor
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if rows of b don't match the size of A.
 */
template <typename T, typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_all_not_st_var<T, EigMat>* = nullptr,
          require_any_not_t<std::is_arithmetic<value_type_t<T>>,
                            is_fvar<value_type_t<EigMat>>>* = nullptr>
inline Eigen::Matrix<return_type_t<T, EigMat>, Eigen::Dynamic,
                     EigMat::ColsAtCompileTime>
mdivide_left_ldlt(LDLT_factor<T>& A, const EigMat& b) {
  check_multiplicable("mdivide_left_ldlt", "A", A.matrix(), "b", b);

  if (A.matrix().cols() == 0) {
    return {0, b.cols()};
  }

  return A.ldlt().solve(
      Eigen::Matrix<return_type_t<T, EigMat>, EigMat::RowsAtCompileTime,
                    EigMat::ColsAtCompileTime>(b));
}

}  // namespace math
}  // namespace stan

#endif
