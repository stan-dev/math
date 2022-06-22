#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_TRI_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_TRI_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

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
          require_all_not_eigen_vt<is_var, T1, T2> * = nullptr>
inline auto mdivide_left_tri(const T1 &A, const T2 &b) {
  using T_return = return_type_t<T1, T2>;
  using ret_type = Eigen::Matrix<T_return, Eigen::Dynamic, Eigen::Dynamic>;
  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "b", b);
  if (A.rows() == 0) {
    return ret_type(0, b.cols());
  }

  return ret_type(A)
      .template triangularView<TriView>()
      .solve(ret_type(b))
      .eval();
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
template <Eigen::UpLoType TriView, typename T, require_eigen_t<T> * = nullptr>
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

}  // namespace math
}  // namespace stan

#endif
