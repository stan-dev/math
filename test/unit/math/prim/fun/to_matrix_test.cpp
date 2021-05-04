#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
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
void test_to_matrix_array_answers(int m, int n) {
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
void test_to_matrix_matrix_answers(int m, int n) {
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
void test_to_matrix_matrix_reshape_answers(int m1, int n1, int m2, int n2) {
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
void test_to_vector_matrix_answers(int m, int m2, int n2) {
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
void test_to_row_vector_matrix_answers(int n, int m2, int n2) {
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
void test_to_matrix_2darray_answers(int m, int n) {
  using stan::math::to_matrix;
  std::vector<std::vector<double> > vec(m, std::vector<double>(n));
  std::vector<std::vector<int> > vec_int(m, std::vector<int>(n));
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
