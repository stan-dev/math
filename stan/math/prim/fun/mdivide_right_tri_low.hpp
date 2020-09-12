#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_TRI_LOW_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_TRI_LOW_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_right_tri.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system x tri(A) = b when tri(A) is a
 * lower triangular view of the matrix A.
 *
 * @tparam EigMat1 type of the right-hand side matrix or vector
 * @tparam EigMat2 type of the second matrix
 *
 * @param b right-hand side matrix or vector
 * @param A matrix
 * @return x = b * tri(A)^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_all_not_vt_fvar<EigMat1, EigMat2>* = nullptr>
inline Eigen::Matrix<return_type_t<EigMat1, EigMat2>,
                     EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime>
mdivide_right_tri_low(const EigMat1& b, const EigMat2& A) {
  return mdivide_right_tri<Eigen::Lower>(b, A);
}

}  // namespace math
}  // namespace stan

#endif
