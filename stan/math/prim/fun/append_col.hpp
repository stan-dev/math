#ifndef STAN_MATH_PRIM_FUN_APPEND_COL_HPP
#define STAN_MATH_PRIM_FUN_APPEND_COL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the result of appending the second argument matrix after the
 * first argument matrix, that is, putting them side by side, with
 * the first matrix followed by the second matrix.
 *
 * Given input types result in following outputs:
 * (matrix, matrix) -> matrix,
 * (matrix, vector) -> matrix,
 * (vector, matrix) -> matrix,
 * (vector, vector) -> matrix,
 * (row vector, row vector) -> row_vector.
 *
 * @tparam EigMat1 type of the first matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param A First matrix.
 * @param B Second matrix.
 * @return Result of appending the first matrix followed by the
 * second matrix side by side.
 */
template <typename EigMat1, typename EigMat2,
          typename = require_all_eigen_t<EigMat1, EigMat2>>
inline auto append_col(const EigMat1& A, const EigMat2& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T_return = return_type_t<EigMat1, EigMat2>;
  constexpr int row_type
      = (EigMat1::RowsAtCompileTime == 1 && EigMat2::RowsAtCompileTime == 1)
            ? 1
            : Eigen::Dynamic;

  int Arows = A.rows();
  int Brows = B.rows();
  int Acols = A.cols();
  int Bcols = B.cols();
  check_size_match("append_col", "rows of A", Arows, "rows of B", Brows);

  Matrix<T_return, row_type, Dynamic> result(Arows, Acols + Bcols);
  result.leftCols(Acols) = A.template cast<T_return>();
  result.rightCols(Bcols) = B.template cast<T_return>();
  return result;
}

/**
 * Return the result of stacking an scalar on top of the
 * a row vector, with the result being a row vector.
 *
 * This function applies to (scalar, row vector) and returns a
 * row vector.
 *
 * @tparam Scal type of the scalar
 * @tparam RowVec type of the row vector
 *
 * @param a scalar.
 * @param B row vector.
 * @return Result of stacking the scalar on top of the row vector.
 */
template <typename Scal, typename RowVec,
          require_stan_scalar_t<Scal>* = nullptr,
          require_t<is_eigen_row_vector<RowVec>>* = nullptr>
inline Eigen::Matrix<return_type_t<Scal, RowVec>, 1, Eigen::Dynamic> append_col(
    const Scal& a, const RowVec& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T_return = return_type_t<Scal, RowVec>;

  Matrix<T_return, 1, Dynamic> result(B.size() + 1);
  result << a, B.template cast<T_return>();
  return result;
}

/**
 * Return the result of stacking a row vector on top of the
 * an scalar, with the result being a row vector.
 *
 * This function applies to (row vector, scalar) and returns a
 * row vector.
 *
 * @tparam RowVec type of the row vector
 * @tparam Scal type of the scalar
 *
 * @param A row vector.
 * @param b scalar.
 * @return Result of stacking the row vector on top of the scalar.
 */
template <typename RowVec, typename Scal,
          require_t<is_eigen_row_vector<RowVec>>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
inline Eigen::Matrix<return_type_t<RowVec, Scal>, 1, Eigen::Dynamic> append_col(
    const RowVec& A, const Scal& b) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T_return = return_type_t<RowVec, Scal>;

  Matrix<T_return, 1, Dynamic> result(A.size() + 1);
  result << A.template cast<T_return>(), b;
  return result;
}

}  // namespace math
}  // namespace stan

#endif
