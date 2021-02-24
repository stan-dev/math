#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_TRI_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_TRI_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system xA=b when A is triangular
 *
 * @tparam TriView Specifies whether A is upper (Eigen::Upper)
 * or lower triangular (Eigen::Lower).
 * @tparam EigMat1 type of the right-hand side matrix or vector
 * @tparam EigMat2 type of the triangular matrix
 *
 * @param A Triangular matrix.  Specify upper or lower with TriView
 * being Eigen::Upper or Eigen::Lower.
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <Eigen::UpLoType TriView, typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr>
inline Eigen::Matrix<return_type_t<EigMat1, EigMat2>,
                     EigMat1::RowsAtCompileTime, EigMat2::ColsAtCompileTime>
mdivide_right_tri(const EigMat1& b, const EigMat2& A) {
  check_square("mdivide_right_tri", "A", A);
  check_multiplicable("mdivide_right_tri", "b", b, "A", A);
  if (TriView != Eigen::Lower && TriView != Eigen::Upper) {
    throw_domain_error("mdivide_right_tri",
                       "triangular view must be Eigen::Lower or Eigen::Upper",
                       "", "");
  }
  if (A.rows() == 0) {
    return {b.rows(), 0};
  }

  return Eigen::Matrix<return_type_t<EigMat1, EigMat2>,
                       EigMat2::RowsAtCompileTime, EigMat2::ColsAtCompileTime>(
             A)
      .template triangularView<TriView>()
      .transpose()
      .solve(
          Eigen::Matrix<return_type_t<EigMat1, EigMat2>,
                        EigMat1::RowsAtCompileTime, EigMat1::ColsAtCompileTime>(
              b)
              .transpose())
      .transpose();
}

}  // namespace math
}  // namespace stan

#endif
