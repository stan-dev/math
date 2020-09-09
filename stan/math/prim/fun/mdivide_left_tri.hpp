#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_TRI_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_TRI_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl.hpp>
#endif

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b when A is triangular.
 *
 * @tparam TriView Specifies whether A is upper (Eigen::Upper)
 * or lower triangular (Eigen::Lower).
 * @tparam T1 type of the triangular matrix
 * @tparam T2 type of the right-hand side matrix or vector
 *
 * @param A Triangular matrix.
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_all_eigen_t<T1, T2> * = nullptr,
          require_any_not_vt_same<double, T1, T2> * = nullptr,
          require_all_not_eigen_vt<is_var, T1, T2> * = nullptr>
inline Eigen::Matrix<return_type_t<T1, T2>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri(const T1 &A, const T2 &b) {
  using T_return = return_type_t<T1, T2>;
  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "b", b);
  if (A.rows() == 0) {
    return {0, b.cols()};
  }

  return A.template cast<T_return>()
      .eval()
      .template triangularView<TriView>()
      .solve(b.template cast<T_return>().eval());
}

/**
 * Returns the solution of the system Ax=b when A is triangular and b=I.
 *
 * @tparam T type of the matrix
 *
 * @param A Triangular matrix.
 * @return x = A^-1 .
 * @throws std::domain_error if A is not square
 */
template <Eigen::UpLoType TriView, typename T, require_eigen_t<T> * = nullptr,
          require_not_vt_same<double, T> * = nullptr>
inline plain_type_t<T> mdivide_left_tri(const T &A) {
  check_square("mdivide_left_tri", "A", A);
  if (A.rows() == 0) {
    return {};
  }

  int n = A.rows();
  plain_type_t<T> b = plain_type_t<T>::Identity(n, n);
  A.template triangularView<TriView>().solveInPlace(b);
  return b;
}

/**
 * Returns the solution of the system Ax=b when A is triangular
 * and A and b are matrices of doubles.
 *
 * @tparam TriView Specifies whether A is upper (Eigen::Upper)
 * or lower triangular (Eigen::Lower).
 * @tparam T1 type of the triangular matrix
 * @tparam T2 type of the right-hand side matrix or vector
 *
 * @param A Triangular matrix.
 * @param b Right hand side matrix or vector.
 * @return x = A^-1 b, solution of the linear system.
 * @throws std::domain_error if A is not square or the rows of b don't
 * match the size of A.
 */
template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_all_eigen_t<T1, T2> * = nullptr,
          require_all_vt_same<double, T1, T2> * = nullptr>
inline Eigen::Matrix<double, T1::RowsAtCompileTime, T2::ColsAtCompileTime>
mdivide_left_tri(const T1 &A, const T2 &b) {
  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "b", b);
  if (A.rows() == 0) {
    return {0, b.cols()};
  }

#ifdef STAN_OPENCL
  if (A.rows()
      >= opencl_context.tuning_opts().tri_inverse_size_worth_transfer) {
    matrix_cl<double> A_cl(A, from_eigen_uplo_type(TriView));
    matrix_cl<double> b_cl(b);
    matrix_cl<double> A_inv_cl = tri_inverse(A_cl);
    matrix_cl<double> C_cl = A_inv_cl * b_cl;
    return from_matrix_cl(C_cl);
  } else {
#endif
    return to_ref(A).template triangularView<TriView>().solve(b);
#ifdef STAN_OPENCL
  }
#endif
}

/**
 * Returns the solution of the system Ax=b when A is triangular, b=I and
 * both are matrices of doubles.
 *
 * @tparam TriView Specifies whether A is upper (Eigen::Upper)
 * or lower triangular (Eigen::Lower).
 * @tparam T type of the matrix
 *
 * @param A Triangular matrix.
 * @return x = A^-1 .
 * @throws std::domain_error if A is not square
 */
template <Eigen::UpLoType TriView, typename T, require_eigen_t<T> * = nullptr,
          require_vt_same<double, T> * = nullptr>
inline plain_type_t<T> mdivide_left_tri(const T &A) {
  check_square("mdivide_left_tri", "A", A);
  if (A.rows() == 0) {
    return {};
  }

  const int n = A.rows();
#ifdef STAN_OPENCL
  if (A.rows()
      >= opencl_context.tuning_opts().tri_inverse_size_worth_transfer) {
    matrix_cl<double> A_cl(A, from_eigen_uplo_type(TriView));
    A_cl = tri_inverse(A_cl);
    return from_matrix_cl(A_cl);
  } else {
#endif
    plain_type_t<T> b = plain_type_t<T>::Identity(n, n);
    A.template triangularView<TriView>().solveInPlace(b);
    return b;
#ifdef STAN_OPENCL
  }
#endif
}

}  // namespace math
}  // namespace stan

#endif
