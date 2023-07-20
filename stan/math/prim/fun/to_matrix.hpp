#ifndef STAN_MATH_PRIM_FUN_TO_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_TO_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a matrix with dynamic dimensions constructed from an
 * Eigen matrix.
 *
 * @tparam EigMat type of the matrix
 *
 * @param x matrix
 * @return the matrix representation of the input
 */
template <typename EigMat, require_eigen_dense_dynamic_t<EigMat>* = nullptr>
inline EigMat to_matrix(EigMat&& x) {
  return std::forward<EigMat>(x);
}

/**
 * Returns a matrix with dynamic dimensions constructed from an
 * Eigen row or column vector.
 * The runtime dimensions will be the same as the input.
 *
 * @tparam EigMat type of the vector/row vector
 *
 * @param matrix input vector/row vector
 * @return the matrix representation of the input
 */
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
inline auto to_matrix(EigVec&& matrix) {
  return Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>(
      std::forward<EigVec>(matrix));
}

/**
 * Returns a matrix representation of a standard vector of Eigen
 * row vectors with the same dimensions and indexing order.
 *
 * @tparam T type of the elements in the vector
 * @param x Eigen vector of vectors of scalar values
 * @return the matrix representation of the input
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> to_matrix(
    const std::vector<Eigen::Matrix<T, 1, Eigen::Dynamic>>& x) {
  int rows = x.size();
  if (rows == 0) {
    return {};
  }
  int cols = x[0].size();
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(rows, cols);
  for (int i = 0, ij = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++, ij++) {
      result.coeffRef(ij) = x[j][i];
    }
  }
  return result;
}

/**
 * Returns a matrix representation of the standard vector of
 * standard vectors with the same dimensions and indexing order.
 *
 * @tparam T type of elements in the vector
 * @param x vector of vectors of scalar values
 * @return the matrix representation of the input
 */
template <typename T>
inline Eigen::Matrix<return_type_t<T, double>, Eigen::Dynamic, Eigen::Dynamic>
to_matrix(const std::vector<std::vector<T>>& x) {
  size_t rows = x.size();
  if (rows == 0) {
    return {};
  }
  size_t cols = x[0].size();
  Eigen::Matrix<return_type_t<T, double>, Eigen::Dynamic, Eigen::Dynamic>
      result(rows, cols);
  for (size_t i = 0, ij = 0; i < cols; i++) {
    for (size_t j = 0; j < rows; j++, ij++) {
      result.coeffRef(ij) = x[j][i];
    }
  }
  return result;
}

/**
 * Returns a matrix representation of the vector in column-major
 * order with the specified number of rows and columns.
 *
 * @tparam EigMat type of the matrix
 *
 * @param x matrix
 * @param m rows
 * @param n columns
 * @return Reshaped inputted matrix
 * @throw <code>std::invalid_argument</code> if the sizes
 * do not match
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>
to_matrix(EigMat&& x, int m, int n) {
  static constexpr const char* function = "to_matrix(matrix)";
  check_size_match(function, "rows * columns", m * n, "vector size", x.size());
  Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic> y
      = std::forward<EigMat>(x);
  y.resize(m, n);
  return y;
}

/**
 * Returns a matrix representation of the vector in column-major
 * order with the specified number of rows and columns.
 *
 * @tparam T type of elements in the vector
 * @param x vector of values
 * @param m rows
 * @param n columns
 * @return the matrix representation of the input
 * @throw <code>std::invalid_argument</code>
 * if the sizes do not match
 */
template <typename T>
inline auto to_matrix(const std::vector<T>& x, int m, int n) {
  static constexpr const char* function = "to_matrix(array)";
  check_size_match(function, "rows * columns", m * n, "vector size", x.size());
  return Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(
      &x[0], m, n);
}

/**
 * Returns a matrix representation of the vector in column-major
 * order with the specified number of rows and columns.
 *
 * @param x vector of values
 * @param m rows
 * @param n columns
 * @return the matrix representation of the input
 * @throw <code>std::invalid_argument</code>
 * if the sizes do not match
 */
inline Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> to_matrix(
    const std::vector<int>& x, int m, int n) {
  static constexpr const char* function = "to_matrix(array)";
  int x_size = x.size();
  check_size_match(function, "rows * columns", m * n, "vector size", x_size);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> result(m, n);
  for (int i = 0; i < x_size; i++) {
    result.coeffRef(i) = x[i];
  }
  return result;
}

/**
 * Returns a matrix representation of the vector in column-major
 * order with the specified number of rows and columns.
 *
 * @tparam EigMat type of the matrix
 *
 * @param x matrix
 * @param m rows
 * @param n columns
 * @param col_major column-major indicator:
 * if 1, output matrix is transversed in column-major order,
 * if 0, output matrix is transversed in row-major order,
 * otherwise function throws std::invalid_argument
 * @return Reshaped inputted matrix
 * @throw <code>std::invalid_argument</code>
 * if the sizes do not match
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>
to_matrix(EigMat&& x, int m, int n, bool col_major) {
  if (col_major) {
    return to_matrix(std::forward<EigMat>(x), m, n);
  } else {
    Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic> res
        = to_matrix(std::forward<EigMat>(x), n, m);
    res.transposeInPlace();
    return res;
  }
}

/**
 * Returns a matrix representation of the vector in column-major
 * order with the specified number of rows and columns.
 *
 * @tparam T type of elements in the vector
 * @param x vector of values
 * @param m rows
 * @param n columns
 * @param col_major column-major indicator:
 * if 1, output matrix is transversed in column-major order,
 * if 0, output matrix is transversed in row-major order,
 * otherwise function throws std::invalid_argument
 * @return the matrix representation of the input
 * @throw <code>std::invalid_argument</code>
 * if the sizes do not match
 */
template <typename T>
inline Eigen::Matrix<return_type_t<T, double>, Eigen::Dynamic, Eigen::Dynamic>
to_matrix(const std::vector<T>& x, int m, int n, bool col_major) {
  if (col_major) {
    return to_matrix(x, m, n);
  }
  check_size_match("to_matrix", "rows * columns", m * n, "matrix size",
                   x.size());
  Eigen::Matrix<return_type_t<T, double>, Eigen::Dynamic, Eigen::Dynamic>
      result(m, n);
  for (int i = 0, ij = 0; i < m; i++) {
    for (int j = 0; j < n; j++, ij++) {
      result.coeffRef(i, j) = x[ij];
    }
  }
  return result;
}

}  // namespace math
}  // namespace stan

#endif
