#ifndef STAN_MATH_PRIM_FUN_APPEND_ROW_HPP
#define STAN_MATH_PRIM_FUN_APPEND_ROW_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the result of stacking the rows of the first argument
 * matrix on top of the second argument matrix.
 *
 * Given input types result in following outputs:
 * (matrix, matrix) -> matrix,
 * (matrix, row_vector) -> matrix,
 * (row_vector, matrix) -> matrix,
 * (row_vector, row_vector) -> matrix,
 * (vector, vector) -> vector.
 *
 * @tparam T1 type of the first matrix
 * @tparam T2 type of the second matrix
 *
 * @param A First matrix.
 * @param B Second matrix.
 * @return Result of stacking first matrix on top of second.
 */
template <typename T1, typename T2, require_all_eigen_t<T1, T2>* = nullptr>
inline auto append_row(const T1& A, const T2& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T_return = return_type_t<T1, T2>;
  constexpr int col_type
      = (T1::ColsAtCompileTime == 1 && T2::ColsAtCompileTime == 1)
            ? 1
            : Eigen::Dynamic;

  int Arows = A.rows();
  int Brows = B.rows();
  int Acols = A.cols();
  int Bcols = B.cols();
  check_size_match("append_row", "columns of A", Acols, "columns of B", Bcols);

  Matrix<T_return, Dynamic, col_type> result(Arows + Brows, Acols);
  result.topRows(Arows) = A.template cast<T_return>();
  result.bottomRows(Brows) = B.template cast<T_return>();
  return result;
}

/**
 * Return the result of stacking an scalar on top of the
 * a vector, with the result being a vector.
 *
 * This function applies to (scalar, vector) and returns a vector.
 *
 * @tparam Scal type of the scalar
 * @tparam ColVec type of the vector
 *
 * @param A scalar.
 * @param B vector.
 * @return Result of stacking the scalar on top of the vector.
 */
template <typename Scal, typename ColVec,
          require_stan_scalar_t<Scal>* = nullptr,
          require_t<is_eigen_col_vector<ColVec>>* = nullptr>
inline Eigen::Matrix<return_type_t<Scal, ColVec>, Eigen::Dynamic, 1> append_row(
    const Scal& A, const ColVec& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T_return = return_type_t<Scal, ColVec>;

  Matrix<T_return, Dynamic, 1> result(B.size() + 1);
  result << A, B.template cast<T_return>();
  return result;
}

/**
 * Return the result of stacking a vector on top of the
 * an scalar, with the result being a vector.
 *
 * This function applies to (vector, scalar) and returns a vector.
 *
 * @tparam ColVec type of the vector
 * @tparam Scal type of the scalar
 *
 * @param A vector.
 * @param B scalar.
 * @return Result of stacking the vector on top of the scalar.
 */
template <typename ColVec, typename Scal,
          require_t<is_eigen_col_vector<ColVec>>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
inline Eigen::Matrix<return_type_t<ColVec, Scal>, Eigen::Dynamic, 1> append_row(
    const ColVec& A, const Scal& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T_return = return_type_t<ColVec, Scal>;

  Matrix<T_return, Dynamic, 1> result(A.size() + 1);
  result << A.template cast<T_return>(), B;
  return result;
}

}  // namespace math
}  // namespace stan

#endif
