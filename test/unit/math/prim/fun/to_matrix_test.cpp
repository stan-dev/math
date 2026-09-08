#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>
#include <stdexcept>

template <typename T, int R, int C>
inline Eigen::Matrix<T, R, C> row_major_to_column_major(
    const Eigen::Matrix<T, R, C>& x) {
  int rows = x.rows();
  int cols = x.cols();
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(rows, cols);
  for (int i = 0, ij = 0; i < rows; i++)
    for (int j = 0; j < cols; j++, ij++)
      result(ij) = x(i, j);
  return result;
}

// [T] -> Matrix
inline void test_to_matrix_array_answers(int m, int n) {
  using stan::math::to_matrix;
  std::vector<double> vec(m * n);
  std::vector<int> vec_int(m * n);
  for (int i = 0; i < m * n; ++i) {
    vec[i] = i;
    vec_int[i] = i;
  }
  Eigen::MatrixXd a(m, n);
  for (int i = 0; i < m * n; ++i)
    a(i) = i;
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(vec, m, n));
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(vec, m, n, 1));
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(vec, m, n, -1));
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(vec, m, n, 2));
  EXPECT_MATRIX_FLOAT_EQ(a, row_major_to_column_major(to_matrix(vec, m, n, 0)));
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(vec_int, m, n));
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(vec_int, m, n, 1));
  EXPECT_MATRIX_FLOAT_EQ(
      a, row_major_to_column_major(to_matrix(vec_int, m, n, 0)));
}

TEST(ToMatrixArray, answers) {
  test_to_matrix_array_answers(0, 0);
  test_to_matrix_array_answers(3, 2);
  test_to_matrix_array_answers(3, 0);
  test_to_matrix_array_answers(0, 3);
}

TEST(ToMatrixArray, exceptions) {
  using stan::math::to_matrix;
  std::vector<double> vec(3);
  EXPECT_THROW(to_matrix(vec, 2, 2), std::invalid_argument);
  EXPECT_THROW(to_matrix(vec, 1, 2), std::invalid_argument);
  EXPECT_NO_THROW(to_matrix(vec, 1, 3));
}

// Matrix -> Matrix
inline void test_to_matrix_matrix_answers(int m, int n) {
  using stan::math::to_matrix;
  Eigen::MatrixXd a(m, n);
  for (int i = 0; i < m * n; ++i)
    a(i) = i;
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(a));
}

TEST(ToMatrixMatrix, answers) {
  test_to_matrix_matrix_answers(0, 0);
  test_to_matrix_matrix_answers(3, 2);
  test_to_matrix_matrix_answers(3, 0);
  test_to_matrix_matrix_answers(0, 3);
}

// Matrix -> Matrix (with reshape)
inline void test_to_matrix_matrix_reshape_answers(int m1, int n1, int m2,
                                                  int n2) {
  using stan::math::to_matrix;
  Eigen::MatrixXd a(m1, n1);
  Eigen::MatrixXd b(m2, n2);
  for (int i = 0; i < m1 * n1; ++i) {
    a(i) = static_cast<double>(i) / 1.26;
    b(i) = static_cast<double>(i) / 1.26;
  }
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(b, m1, n1));
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(b, m1, n1, 1));
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(b, m1, n1, -1));
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(b, m1, n1, 2));
  EXPECT_MATRIX_FLOAT_EQ(a, row_major_to_column_major(to_matrix(b, m1, n1, 0)));

  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2));
  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2, 1));
  EXPECT_MATRIX_FLOAT_EQ(b, row_major_to_column_major(to_matrix(a, m2, n2, 0)));

  if (n1 != 0) {
    EXPECT_THROW(to_matrix(a, m1 + 1, n1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m1 + 1, n1, 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m1 + 1, n1, 0), std::invalid_argument);
  }
  if (m1 != 0) {
    EXPECT_THROW(to_matrix(a, m1, n1 + 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m1, n1 + 1, 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m1, n1 + 1, 0), std::invalid_argument);
  }
  if (n2 != 0) {
    EXPECT_THROW(to_matrix(a, m2 + 1, n2), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2 + 1, n2, 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2 + 1, n2, 0), std::invalid_argument);
  }
  if (m2 != 0) {
    EXPECT_THROW(to_matrix(a, m2, n2 + 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2, n2 + 1, 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2, n2 + 1, 0), std::invalid_argument);
  }
}

TEST(ToMatrixMatrixReshape, answers) {
  test_to_matrix_matrix_reshape_answers(0, 0, 0, 0);
  test_to_matrix_matrix_reshape_answers(3, 2, 2, 3);
  test_to_matrix_matrix_reshape_answers(3, 2, 6, 1);
  test_to_matrix_matrix_reshape_answers(3, 0, 0, 3);
  test_to_matrix_matrix_reshape_answers(8, 2, 4, 4);
}

// Vector -> Matrix
inline void test_to_vector_matrix_answers(int m, int m2, int n2) {
  using stan::math::to_matrix;
  Eigen::VectorXd a(m);
  Eigen::MatrixXd b(m2, n2);
  Eigen::MatrixXd c(m, 1);
  for (int i = 0; i < m2 * n2; ++i) {
    a(i) = static_cast<double>(i) / 1.26;
    b(i) = static_cast<double>(i) / 1.26;
    c(i) = static_cast<double>(i) / 1.26;
  }
  // without reshape
  EXPECT_MATRIX_FLOAT_EQ(c, to_matrix(a));

  // with reshape
  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2));
  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2, 1));
  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2, -1));
  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2, 2));
  EXPECT_MATRIX_FLOAT_EQ(b, row_major_to_column_major(to_matrix(a, m2, n2, 0)));

  if (n2 != 0) {
    EXPECT_THROW(to_matrix(a, m2 + 1, n2), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2 + 1, n2, 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2 + 1, n2, 0), std::invalid_argument);
  }
  if (m2 != 0) {
    EXPECT_THROW(to_matrix(a, m2, n2 + 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2, n2 + 1, 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2, n2 + 1, 0), std::invalid_argument);
  }
}

TEST(ToMatrixVector, answers) {
  test_to_vector_matrix_answers(0, 0, 0);
  test_to_vector_matrix_answers(6, 2, 3);
  test_to_vector_matrix_answers(18, 6, 3);
  test_to_vector_matrix_answers(0, 0, 3);
  test_to_vector_matrix_answers(8, 1, 8);
}

// RowVector -> Matrix
inline void test_to_row_vector_matrix_answers(int n, int m2, int n2) {
  using stan::math::to_matrix;
  Eigen::RowVectorXd a(n);
  Eigen::MatrixXd b(m2, n2);
  Eigen::MatrixXd c(1, n);
  for (int i = 0; i < m2 * n2; ++i) {
    a(i) = static_cast<double>(i) / 1.26;
    b(i) = static_cast<double>(i) / 1.26;
    c(i) = static_cast<double>(i) / 1.26;
  }
  // without reshape
  EXPECT_MATRIX_FLOAT_EQ(c, to_matrix(a));

  // with reshape
  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2));
  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2, 1));
  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2, -1));
  EXPECT_MATRIX_FLOAT_EQ(b, to_matrix(a, m2, n2, 2));
  EXPECT_MATRIX_FLOAT_EQ(b, row_major_to_column_major(to_matrix(a, m2, n2, 0)));

  if (n2 != 0) {
    EXPECT_THROW(to_matrix(a, m2 + 1, n2), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2 + 1, n2, 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2 + 1, n2, 0), std::invalid_argument);
  }
  if (m2 != 0) {
    EXPECT_THROW(to_matrix(a, m2, n2 + 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2, n2 + 1, 1), std::invalid_argument);
    EXPECT_THROW(to_matrix(a, m2, n2 + 1, 0), std::invalid_argument);
  }
}

TEST(ToMatrixRowVector, answers) {
  test_to_row_vector_matrix_answers(0, 0, 0);
  test_to_row_vector_matrix_answers(6, 2, 3);
  test_to_row_vector_matrix_answers(18, 6, 3);
  test_to_row_vector_matrix_answers(0, 3, 0);
  test_to_row_vector_matrix_answers(8, 1, 8);
}

// [[T]] -> Matrix
inline void test_to_matrix_2darray_answers(int m, int n) {
  using stan::math::to_matrix;
  std::vector<std::vector<double>> vec(m, std::vector<double>(n));
  std::vector<std::vector<int>> vec_int(m, std::vector<int>(n));
  // Any vec (0, C) will become (0, 0)
  if (m == 0)
    n = 0;
  Eigen::MatrixXd a(m, n);

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      vec[i][j] = i * j;
      vec_int[i][j] = i * j;
      a(i, j) = i * j;
    }
  }
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(vec));
  EXPECT_MATRIX_FLOAT_EQ(a, to_matrix(vec_int));
}

TEST(ToMatrix2dArray, answers) {
  test_to_matrix_2darray_answers(0, 0);
  test_to_matrix_2darray_answers(3, 2);
  test_to_matrix_2darray_answers(3, 0);
  test_to_matrix_2darray_answers(0, 3);
}

TEST(ToMatrixVectorArray, answers) {
  using stan::math::to_matrix;
  using stan::math::to_vector_array;
  using vector_d = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  std::vector<vector_d> vecs(2, vector_d(3));
  vecs[0] << 1.1, 2.2, 3.3;
  vecs[1] << 4.4, 5.5, 6.6;

  Eigen::MatrixXd expected(3, 2);
  expected << 1.1, 4.4, 2.2, 5.5, 3.3, 6.6;

  Eigen::MatrixXd result = to_matrix(vecs);
  EXPECT_MATRIX_FLOAT_EQ(expected, result);

  std::vector<vector_d> round_trip = to_vector_array(result);
  ASSERT_EQ(vecs.size(), round_trip.size());
  for (size_t i = 0; i < vecs.size(); ++i) {
    EXPECT_MATRIX_FLOAT_EQ(vecs[i], round_trip[i]);
  }
}

TEST(ToMatrixVectorArray, empty) {
  using stan::math::to_matrix;
  using stan::math::to_vector_array;
  using vector_d = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  std::vector<vector_d> empty;
  Eigen::MatrixXd empty_result = to_matrix(empty);
  EXPECT_EQ(0, empty_result.rows());
  EXPECT_EQ(0, empty_result.cols());

  std::vector<vector_d> zero_rows(3, vector_d(0));
  Eigen::MatrixXd zero_rows_result = to_matrix(zero_rows);
  EXPECT_EQ(0, zero_rows_result.rows());
  EXPECT_EQ(3, zero_rows_result.cols());

  std::vector<vector_d> round_trip = to_vector_array(zero_rows_result);
  ASSERT_EQ(zero_rows.size(), round_trip.size());
  for (size_t i = 0; i < zero_rows.size(); ++i) {
    EXPECT_MATRIX_FLOAT_EQ(zero_rows[i], round_trip[i]);
  }
}

TEST(ToMatrixVectorArray, preservesScalarType) {
  using stan::math::to_matrix;
  using vector_i = Eigen::Matrix<int, Eigen::Dynamic, 1>;

  std::vector<vector_i> vecs(2, vector_i(2));
  vecs[0] << 1, 2;
  vecs[1] << 3, 4;

  auto result = to_matrix(vecs);
  static_assert(
      std::is_same<decltype(result),
                   Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>::value,
      "to_matrix should preserve the vector scalar type");
  EXPECT_EQ(1, result(0, 0));
  EXPECT_EQ(2, result(1, 0));
  EXPECT_EQ(3, result(0, 1));
  EXPECT_EQ(4, result(1, 1));
}

TEST(ToMatrixRowVectorArray, answers) {
  using stan::math::to_matrix;
  using stan::math::to_row_vector_array;
  using row_vector_d = Eigen::Matrix<double, 1, Eigen::Dynamic>;

  std::vector<row_vector_d> row_vecs(3, row_vector_d(2));
  row_vecs[0] << 1.1, 4.4;
  row_vecs[1] << 2.2, 5.5;
  row_vecs[2] << 3.3, 6.6;

  Eigen::MatrixXd expected(3, 2);
  expected << 1.1, 4.4, 2.2, 5.5, 3.3, 6.6;

  Eigen::MatrixXd result = to_matrix(row_vecs);
  EXPECT_MATRIX_FLOAT_EQ(expected, result);

  std::vector<row_vector_d> round_trip = to_row_vector_array(result);
  ASSERT_EQ(row_vecs.size(), round_trip.size());
  for (size_t i = 0; i < row_vecs.size(); ++i) {
    EXPECT_MATRIX_FLOAT_EQ(row_vecs[i], round_trip[i]);
  }
}

TEST(ToMatrixRowVectorArray, empty) {
  using stan::math::to_matrix;
  using stan::math::to_row_vector_array;
  using row_vector_d = Eigen::Matrix<double, 1, Eigen::Dynamic>;

  std::vector<row_vector_d> empty;
  Eigen::MatrixXd empty_result = to_matrix(empty);
  EXPECT_EQ(0, empty_result.rows());
  EXPECT_EQ(0, empty_result.cols());

  std::vector<row_vector_d> zero_cols(3, row_vector_d(0));
  Eigen::MatrixXd zero_cols_result = to_matrix(zero_cols);
  EXPECT_EQ(3, zero_cols_result.rows());
  EXPECT_EQ(0, zero_cols_result.cols());

  std::vector<row_vector_d> round_trip = to_row_vector_array(zero_cols_result);
  ASSERT_EQ(zero_cols.size(), round_trip.size());
  for (size_t i = 0; i < zero_cols.size(); ++i) {
    EXPECT_MATRIX_FLOAT_EQ(zero_cols[i], round_trip[i]);
  }
}
