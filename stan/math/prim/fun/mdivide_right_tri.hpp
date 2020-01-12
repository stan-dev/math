#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_TRI_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_RIGHT_TRI_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>
#include <stan/math/prim/fun/promote_common.hpp>
#include <stan/math/prim/fun/transpose.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system xA=b when A is triangular
 *
 * @tparam TriView Specifies whether A is upper (Eigen::Upper)
 * or lower triangular (Eigen::Lower).
 * @tparam T1 type of elements in the right-hand side matrix or vector
 * @tparam T2 type of elements in the triangular matrix
 * @tparam R1 number of rows in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam C1 number of columns in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam R2 number of rows in the triangular matrix, can be Eigen::Dynamic
 * @tparam C2 number of columns in the triangular matrix, can be Eigen::Dynamic
 *
 * @param A Triangular matrix.  Specify upper or lower with TriView
 * being Eigen::Upper or Eigen::Lower.
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <Eigen::UpLoType TriView, typename T1, typename T2, int R1, int C1,
          int R2, int C2>
inline Eigen::Matrix<return_type_t<T1, T2>, R1, C2> mdivide_right_tri(
    const Eigen::Matrix<T1, R1, C1> &b, const Eigen::Matrix<T2, R2, C2> &A) {
  check_square("mdivide_right_tri", "A", A);
  check_multiplicable("mdivide_right_tri", "b", b, "A", A);
  if (TriView != Eigen::Lower && TriView != Eigen::Upper) {
    throw_domain_error("mdivide_left_tri",
                       "triangular view must be Eigen::Lower or Eigen::Upper",
                       "", "");
  }
  return promote_common<Eigen::Matrix<T1, R2, C2>, Eigen::Matrix<T2, R2, C2> >(
             A)
      .template triangularView<TriView>()
      .transpose()
      .solve(
          promote_common<Eigen::Matrix<T1, R1, C1>, Eigen::Matrix<T2, R1, C1> >(
              b)
              .transpose())
      .transpose();
}

/**
 * Returns the solution of the system xA=b when A is triangular
 * and A and b are matrices of doubles.
 *
 * @tparam TriView Specifies whether A is upper (Eigen::Upper)
 * or lower triangular (Eigen::Lower).
 * @tparam R1 number of rows in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam C1 number of columns in the right-hand side matrix, can be
 *         Eigen::Dynamic
 * @tparam R2 number of rows in the triangular matrix, can be Eigen::Dynamic
 * @tparam C2 number of columns in the triangular matrix, can be Eigen::Dynamic
 *
 * @param A Triangular matrix.  Specify upper or lower with TriView
 * being Eigen::Upper or Eigen::Lower.
 * @param b Right hand side matrix or vector.
 * @return x = b A^-1, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<double, R1, C2> mdivide_right_tri(
    const Eigen::Matrix<double, R1, C1> &b,
    const Eigen::Matrix<double, R2, C2> &A) {
  check_square("mdivide_right_tri", "A", A);
  check_multiplicable("mdivide_right_tri", "b", b, "A", A);
#ifdef STAN_OPENCL
  if (A.rows()
      >= opencl_context.tuning_opts().tri_inverse_size_worth_transfer) {
    matrix_cl<double> A_cl(A, from_eigen_uplo_type(TriView));
    matrix_cl<double> b_cl(b);
    matrix_cl<double> A_inv_cl = tri_inverse(A_cl);
    matrix_cl<double> C_cl = b_cl * A_inv_cl;
    return from_matrix_cl(C_cl);
  } else {
#endif
    return A.template triangularView<TriView>()
        .transpose()
        .solve(b.transpose())
        .transpose();
#ifdef STAN_OPENCL
  }
#endif
}

}  // namespace math
}  // namespace stan

#endif
