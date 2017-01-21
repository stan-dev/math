#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATRIX_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATRIX_ERROR_HPP

#include <test/unit/math/prim/mat/vectorize/expect_binary_scalar_matrix_err_throw.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_matrix_scalar_err_throw.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_scalar_std_vector_matrix_err_throw.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_std_vector_matrix_scalar_err_throw.hpp>
#include <test/unit/math/prim/mat/vectorize/build_prim_binary_matrix.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>

template <typename F, typename T>
void expect_binary_matrix_size_error() {
  using std::vector;
  using Eigen::MatrixXd;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

  matrix_t badsize_tm1(4, 5);
  matrix_t badsize_tm2(4, 9);
  matrix_t badsize_tm3(7, 5);
  matrix_t badsize_tm4(3, 6);
  matrix_t badsize_tm5(5, 4);
  MatrixXd badsize_dm1(6, 13);
  MatrixXd badsize_dm2(4, 8);
  MatrixXd badsize_dm3(12, 5);
  MatrixXd badsize_dm4(5, 4);
  T test_val = F::valid_inputs1()[0];

  badsize_tm1 = build_prim_binary_matrix(test_val, badsize_tm1);
  badsize_tm2 = build_prim_binary_matrix(test_val, badsize_tm2);
  badsize_tm3 = build_prim_binary_matrix(test_val, badsize_tm3);
  badsize_tm4 = build_prim_binary_matrix(test_val, badsize_tm4);
  badsize_tm5 = build_prim_binary_matrix(test_val, badsize_tm5);
  badsize_dm1 = build_prim_binary_matrix(F::valid_inputs1()[0],
  badsize_dm1);
  badsize_dm2 = build_prim_binary_matrix(F::valid_inputs1()[0],
  badsize_dm2);
  badsize_dm3 = build_prim_binary_matrix(F::valid_inputs1()[0],
  badsize_dm3);
  badsize_dm4 = build_prim_binary_matrix(F::valid_inputs1()[0],
  badsize_dm4);

  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_dm1, badsize_tm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm2),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_dm2, badsize_tm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm3),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_dm3, badsize_tm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm4),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_dm4, badsize_tm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm2),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm2, badsize_tm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm3),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm3, badsize_tm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm4),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm4, badsize_tm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm5),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm5, badsize_tm1),
  std::invalid_argument);
}

template <typename F, typename T>
void expect_binary_matrix_value_error() {
  using std::vector;
  using Eigen::MatrixXd;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;


  vector<double> invalid_inputs1 = F::invalid_inputs1();
  if (invalid_inputs1.size() == 0) return;
  vector<double> invalid_inputs2 = F::invalid_inputs2();
  vector<int> int_invalid_inputs1 = F::int_invalid_inputs1();
  vector<int> int_invalid_inputs2 = F::int_invalid_inputs2();
  vector<T> y1(invalid_inputs1.begin(), invalid_inputs1.end());
  vector<T> y2(invalid_inputs2.begin(), invalid_inputs2.end());
  MatrixXd a1(3, invalid_inputs1.size());
  MatrixXd a2(3, invalid_inputs2.size());
  matrix_t b1(3, invalid_inputs1.size());
  matrix_t b2(3, invalid_inputs2.size());
  for (int i = 0; i < a1.rows(); ++i) {
    for (int j = 0; j < a1.cols(); ++j) {
      a1(i, j) = invalid_inputs1[j];
    }
  }
  for (int i = 0; i < a2.rows(); ++i) {
    for (int j = 0; j < a2.cols(); ++j) {
      a2(i, j) = invalid_inputs2[j];
    }
  }
  for (int i = 0; i < b1.rows(); ++i) {
    for (int j = 0; j < b1.cols(); ++j) {
      b1(i, j) = invalid_inputs1[j];
    }
  }
  for (int i = 0; i < b2.rows(); ++i) {
    for (int j = 0; j < b2.cols(); ++j) {
      b2(i, j) = invalid_inputs2[j];
    }
  }
  //scalar, matrix
  //int, T
  expect_binary_scalar_matrix_err_throw<F, T>(int_invalid_inputs1, b2);
  expect_binary_matrix_scalar_err_throw<F, T>(b1, int_invalid_inputs2);
  //double, T 
  expect_binary_scalar_matrix_err_throw<F, T>(y1, a2);
  expect_binary_matrix_scalar_err_throw<F, T>(a1, y2);
  expect_binary_scalar_matrix_err_throw<F, T>(invalid_inputs1, b2);
  expect_binary_matrix_scalar_err_throw<F, T>(b1, invalid_inputs2);
  //T, T
  expect_binary_scalar_matrix_err_throw<F, T>(y1, b2);
  expect_binary_matrix_scalar_err_throw<F, T>(b1, y2);
  
  //matrix, matrix
  EXPECT_THROW(F::template apply<matrix_t>(a1, b2), std::domain_error);
  EXPECT_THROW(F::template apply<matrix_t>(b1, a2), std::domain_error);
  EXPECT_THROW(F::template apply<matrix_t>(b1, b2), std::domain_error);
  EXPECT_THROW(F::template apply<matrix_t>(b1.block(1, 1, 1, 1), 
  a2.block(1, 1, 1, 1)), std::domain_error);
  EXPECT_THROW(F::template apply<matrix_t>(a1.block(1, 1, 1, 1), 
  b2.block(1, 1, 1, 1)), std::domain_error);
  EXPECT_THROW(F::template apply<matrix_t>(b1.block(1, 1, 1, 1), 
  b2.block(1, 1, 1, 1)), std::domain_error);

  //scalar, vector<matrix>
  //int, T
  expect_binary_scalar_std_vector_matrix_err_throw<F, T>(
  int_invalid_inputs1, b2);
  expect_binary_std_vector_matrix_scalar_err_throw<F, T>(
  b1, int_invalid_inputs2);
  //double, T
  expect_binary_scalar_std_vector_matrix_err_throw<F, T>(y1, a2);
  expect_binary_std_vector_matrix_scalar_err_throw<F, T>(a1, y2);
  expect_binary_scalar_std_vector_matrix_err_throw<F, T>(
  invalid_inputs1, b2);
  expect_binary_std_vector_matrix_scalar_err_throw<F, T>(
  b1, invalid_inputs2);
  //T, T
  expect_binary_scalar_std_vector_matrix_err_throw<F, T>(y1, b2);
  expect_binary_std_vector_matrix_scalar_err_throw<F, T>(b1, y2);
  //vector<matrix>, vector<matrix>
  vector<MatrixXd> d1;
  d1.push_back(a1);
  d1.push_back(a1);
  vector<MatrixXd> d2;
  d2.push_back(a2);
  d2.push_back(a2);
  vector<matrix_t> e1;
  e1.push_back(b1);
  e1.push_back(b1);
  vector<matrix_t> e2;
  e2.push_back(b2);
  e2.push_back(b2);
  EXPECT_THROW(F::template apply<vector<matrix_t> >(e1, d2), 
  std::domain_error);
  EXPECT_THROW(F::template apply<vector<matrix_t> >(d1, e2), 
  std::domain_error);
  EXPECT_THROW(F::template apply<vector<matrix_t> >(e1, e2), 
  std::domain_error);
}


template <typename F, typename T>
void expect_binary_matrix_error() {
  expect_binary_matrix_size_error<F, T>();
  expect_binary_matrix_value_error<F, T>();
}

#endif
