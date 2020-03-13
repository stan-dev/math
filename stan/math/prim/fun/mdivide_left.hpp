#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b.
 *
 * @tparam T1 type of the first matrix
 * @tparam T2 type of the right-hand side matrix or vector
 *
 * @param A Matrix.
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename T1, typename T2,
          require_all_eigen_vt<std::is_arithmetic, T1, T2>* = nullptr>
inline Eigen::Matrix<return_type_t<T1, T2>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left(const T1& A, const T2& b) {
  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  return Eigen::Matrix<return_type_t<T1, T2>, T1::RowsAtCompileTime,
                       T1::ColsAtCompileTime>(A)
      .lu()
      .solve(Eigen::Matrix<return_type_t<T1, T2>, T2::RowsAtCompileTime,
                           T2::ColsAtCompileTime>(b));
}

}  // namespace math
}  // namespace stan

#endif
