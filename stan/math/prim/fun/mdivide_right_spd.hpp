#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_SPD_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_SPD_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_left_spd.hpp>
#include <stan/math/prim/fun/transpose.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system xA=b where A is symmetric
 * positive definite.
 *
 * @tparam EigMat1 type of the right-hand side matrix or vector
 * @tparam EigMat2 type of the second matrix
 *
 * @param b right-hand side matrix or vector
 * @param A matrix
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr>
inline Eigen::Matrix<return_type_t<EigMat1, EigMat2>,
                     EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime>
mdivide_right_spd(const EigMat1& b, const EigMat2& A) {
  static const char* function = "mdivide_right_spd";
  check_multiplicable(function, "b", b, "A", A);
  check_symmetric(function, "A", A);
  check_not_nan(function, "A", A);
  if (A.size() == 0) {
    return {b.rows(), 0};
  }

  return mdivide_left_spd(A, b.transpose()).transpose();
}

}  // namespace math
}  // namespace stan

#endif
