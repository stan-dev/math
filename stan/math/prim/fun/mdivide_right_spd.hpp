#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_SPD_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_SPD_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_left_spd.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

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
                     std::decay_t<EigMat1>::RowsAtCompileTime, std::decay_t<EigMat2>::ColsAtCompileTime>
mdivide_right_spd(EigMat1&& b, EigMat2&& A) {
  static constexpr const char* function = "mdivide_right_spd";
  check_multiplicable(function, "b", b, "A", A);
  auto&& A_ref = to_ref(std::forward<EigMat2>(A));
  check_symmetric(function, "A", A_ref);
  check_not_nan(function, "A", A_ref);
  if (A.size() == 0) {
    return {b.rows(), 0};
  }

  return mdivide_left_spd(std::forward<decltype(A_ref)>(A_ref), std::forward<EigMat1>(b).transpose()).transpose();
}

}  // namespace math
}  // namespace stan

#endif
