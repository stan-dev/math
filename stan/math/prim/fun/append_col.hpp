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
 * The inputs can be
 * (matrix, matrix),
 * (matrix, vector),
 * (vector, matrix), or
 * (vector, vector)
 * and the output is always a matrix.
 *
 * @tparam T1 type of elements in the first matrix
 * @tparam T2 type of elements in the second matrix
 * @tparam R1 number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam C1 number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam R2 number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam C2 number of columns in the second matrix, can be Eigen::Dynamic
 *
 * @param A First matrix.
 * @param B Second matrix.
 * @return Result of appending the first matrix followed by the
 * second matrix side by side.
 */
template <typename T1, typename T2, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<return_type_t<T1, T2>, Eigen::Dynamic, Eigen::Dynamic>
append_col(const Eigen::Matrix<T1, R1, C1>& A,
           const Eigen::Matrix<T2, R2, C2>& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  int Arows = A.rows();
  int Brows = B.rows();
  int Acols = A.cols();
  int Bcols = B.cols();
  check_size_match("append_col", "rows of A", Arows, "rows of B", Brows);

  Matrix<return_type_t<T1, T2>, Dynamic, Dynamic> result(Arows, Acols + Bcols);
  for (int j = 0; j < Acols; j++) {
    for (int i = 0; i < Arows; i++) {
      result(i, j) = A(i, j);
    }
  }

  for (int j = Acols, k = 0; k < Bcols; j++, k++) {
    for (int i = 0; i < Arows; i++) {
      result(i, j) = B(i, k);
    }
  }
  return result;
}

/**
 * Return the result of concatenating the first row vector followed
 * by the second row vector side by side, with the result being a
 * row vector.
 *
 * This function applies to (row_vector, row_vector) and returns a
 * row_vector.
 *
 * @tparam T1 type of elements in the first row vector
 * @tparam T2 type of elements in the second row vector
 * @tparam C1 number of columns in the first row vector, can be Eigen::Dynamic
 * @tparam C2 number of columns in the second row vector, can be Eigen::Dynamic
 *
 * @param A First vector.
 * @param B Second vector
 * @return Result of appending the second row vector to the right
 * of the first row vector.
 */
template <typename T1, typename T2, int C1, int C2>
inline Eigen::Matrix<return_type_t<T1, T2>, 1, Eigen::Dynamic> append_col(
    const Eigen::Matrix<T1, 1, C1>& A, const Eigen::Matrix<T2, 1, C2>& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  int Asize = A.size();
  int Bsize = B.size();
  Matrix<return_type_t<T1, T2>, 1, Dynamic> result(Asize + Bsize);
  for (int i = 0; i < Asize; i++) {
    result(i) = A(i);
  }
  for (int i = 0, j = Asize; i < Bsize; i++, j++) {
    result(j) = B(i);
  }
  return result;
}

/**
 * Return the result of appending the second argument matrix after the
 * first argument matrix, that is, putting them side by side, with
 * the first matrix followed by the second matrix.   This is an
 * overloaded template function for the case when both matrices
 * have the same type.
 *
 * The inputs can be
 * (matrix, matrix),
 * (matrix, vector),
 * (vector, matrix), or
 * (vector, vector),
 * and the output is always a matrix.
 *
 * @tparam T type of elements in both matrices
 * @tparam R1 number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam C1 number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam R2 number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam C2 number of columns in the second matrix, can be Eigen::Dynamic
 *
 * @param A First matrix.
 * @param B Second matrix.
 * @return Result of appending the first matrix followed by the
 * second matrix side by side.
 */
template <typename T, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> append_col(
    const Eigen::Matrix<T, R1, C1>& A, const Eigen::Matrix<T, R2, C2>& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  check_size_match("append_col", "rows of A", A.rows(), "rows of B", B.rows());

  Matrix<T, Dynamic, Dynamic> result(A.rows(), A.cols() + B.cols());
  result << A, B;
  return result;
}

/**
 * Return the result of concatenating the first row vector followed
 * by the second row vector side by side, with the result being a
 * row vector.
 *
 * This function applies to (row_vector, row_vector) and returns a
 * row_vector.
 *
 * @tparam T type of elements in both vectors
 * @tparam C1 number of columns in the first row vector, can be Eigen::Dynamic
 * @tparam C2 number of columns in the second row vector, can be Eigen::Dynamic
 *
 * @param A First vector.
 * @param B Second vector
 * @return Result of appending the second row vector to the right
 * of the first row vector.
 */
template <typename T, int C1, int C2>
inline Eigen::Matrix<T, 1, Eigen::Dynamic> append_col(
    const Eigen::Matrix<T, 1, C1>& A, const Eigen::Matrix<T, 1, C2>& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  Matrix<T, 1, Dynamic> result(A.size() + B.size());
  result << A, B;
  return result;
}

/**
 * Return the result of stacking an scalar on top of the
 * a row vector, with the result being a row vector.
 *
 * This function applies to (scalar, row vector) and returns a
 * row vector.
 *
 * @tparam T1 type of the scalar
 * @tparam T2 type of elements in the row vector
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param A scalar.
 * @param B row vector.
 * @return Result of stacking the scalar on top of the row vector.
 */
template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<return_type_t<T1, T2>, 1, Eigen::Dynamic> append_col(
    const T1& A, const Eigen::Matrix<T2, R, C>& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using return_type = return_type_t<T1, T2>;

  Matrix<return_type, 1, Dynamic> result(B.size() + 1);
  result << A, B.template cast<return_type>();
  return result;
}

/**
 * Return the result of stacking a row vector on top of the
 * an scalar, with the result being a row vector.
 *
 * This function applies to (row vector, scalar) and returns a
 * row vector.
 *
 * @tparam T1 type of elements in the row vector
 * @tparam T2 type of the scalar
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param A row vector.
 * @param B scalar.
 * @return Result of stacking the row vector on top of the scalar.
 */
template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<return_type_t<T1, T2>, 1, Eigen::Dynamic> append_col(
    const Eigen::Matrix<T1, R, C>& A, const T2& B) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using return_type = return_type_t<T1, T2>;

  Matrix<return_type, 1, Dynamic> result(A.size() + 1);
  result << A.template cast<return_type>(), B;
  return result;
}

}  // namespace math
}  // namespace stan

#endif
