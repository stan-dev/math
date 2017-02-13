#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>

using stan::math::to_matrix;

// [T] -> Matrix
void test_to_matrix_array_answers(int m, int n) {
  std::vector<double> vec(m * n);
  for (int i = 0; i < m * n; ++i)
    vec[i] = i;
  Eigen::MatrixXd a(m, n);
  for (int i = 0; i < m * n; ++i)
    a(i) = i;
  expect_matrix_eq(a, to_matrix(vec, m, n));
}

TEST(ToMatrixArray, answers) {
  test_to_matrix_array_answers(0, 0);
  test_to_matrix_array_answers(3, 2);
  test_to_matrix_array_answers(3, 0);
  test_to_matrix_array_answers(0, 3);
}

TEST(ToMatrixArray, exceptions) {
  std::vector<double> vec(3);
  EXPECT_THROW(to_matrix(vec, 2, 2), std::invalid_argument);
  EXPECT_THROW(to_matrix(vec, 1, 2), std::invalid_argument);
  EXPECT_NO_THROW(to_matrix(vec, 1, 3));
}

// Matrix -> Matrix
void test_to_matrix_matrix_answers(int m, int n) {
  Eigen::MatrixXd a(m, n);
  for (int i = 0; i < m * n; ++i)
    a(i) = i;
  expect_matrix_eq(a, to_matrix(a));
}

TEST(ToMatrixMatrix, answers) {
  test_to_matrix_matrix_answers(0, 0);
  test_to_matrix_matrix_answers(3, 2);
  test_to_matrix_matrix_answers(3, 0);
  test_to_matrix_matrix_answers(0, 3);
}

// [[T]] -> Matrix
void test_to_matrix_2darray_answers(int m, int n) {
  std::vector<std::vector<double> > vec(m, std::vector<double>(n));
  if (m == 0) n = 0; // Any vec (0, C) will become (0, 0)
  Eigen::MatrixXd a(m, n);

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      vec[i][j] = i * j;
      a(i, j) = i * j;
    }
  }
  expect_matrix_eq(a, to_matrix(vec));
}

TEST(ToMatrix2dArray, answers) {
  test_to_matrix_2darray_answers(0, 0);
  test_to_matrix_2darray_answers(3, 2);
  test_to_matrix_2darray_answers(3, 0);
  test_to_matrix_2darray_answers(0, 3);
}

// [[int]] -> Matrix
void test_to_matrix_2d_int_array_answers(int m, int n) {
  std::vector<std::vector<int> > vec(m, std::vector<int>(n));
  if (m == 0) n = 0; // Any vec (0, C) will become (0, 0)
  Eigen::MatrixXd a(m, n);

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      vec[i][j] = i * j;
      a(i, j) = i * j;
    }
  }
  expect_matrix_eq(a, to_matrix(vec));
}

TEST(ToMatrix2dIntArray, answers) {
  test_to_matrix_2d_int_array_answers(0, 0);
  test_to_matrix_2d_int_array_answers(3, 2);
  test_to_matrix_2d_int_array_answers(3, 0);
  test_to_matrix_2d_int_array_answers(0, 3);
}
