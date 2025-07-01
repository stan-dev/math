#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_SPD_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_SPD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b where A is symmetric positive
 * definite.
 *
 * @tparam EigMat1 type of the first matrix
 * @tparam EigMat2 type of the right-hand side matrix or vector
 *
 * @param A Matrix.
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_all_not_vt_var<EigMat1, EigMat2>* = nullptr>
inline Eigen::Matrix<return_type_t<EigMat1, EigMat2>,
                     std::decay_t<EigMat1>::RowsAtCompileTime, std::decay_t<EigMat2>::ColsAtCompileTime>
mdivide_left_spd(EigMat1&& A, EigMat2&& b) {
  static constexpr const char* function = "mdivide_left_spd";
  check_multiplicable(function, "A", A, "b", b);
  auto&& A_ref = to_ref(std::forward<EigMat1>(A));
  check_symmetric(function, "A", A_ref);
  check_not_nan(function, "A", A_ref);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  auto llt
      = Eigen::Matrix<return_type_t<EigMat1, EigMat2>,
                      std::decay_t<EigMat1>::RowsAtCompileTime, std::decay_t<EigMat1>::ColsAtCompileTime>(
            std::forward<decltype(A_ref)>(A_ref)).llt();
  check_pos_definite(function, "A", llt);
  return std::move(llt).solve(
      Eigen::Matrix<return_type_t<EigMat1, EigMat2>, std::decay_t<EigMat2>::RowsAtCompileTime,
                    std::decay_t<EigMat2>::ColsAtCompileTime>(std::forward<EigMat2>(b)));
}

}  // namespace math
}  // namespace stan

#endif
