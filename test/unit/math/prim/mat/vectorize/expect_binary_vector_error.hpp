#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_VECTOR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_VECTOR_ERROR_HPP

#include <test/unit/math/prim/mat/vectorize/expect_binary_scalar_matrix_err_throw.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_scalar_std_vector_matrix_err_throw.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>

template <typename F, typename T, typename VT>
void expect_binary_vector_size_error() {
  using std::vector;
  using Eigen::VectorXd;
  typedef Eigen::Matrix<T, VT::RowsAtCompileTime, VT::ColsAtCompileTime> 
  vector_t;
  vector_t badsize_tm1(4);
  vector_t badsize_tm2(9);
  VT badsize_dm(11);

  EXPECT_THROW(F::template apply<vector_t>(badsize_tm1, badsize_dm),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<vector_t>(badsize_dm, badsize_tm1),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<vector_t>(badsize_tm1, badsize_tm2),
  std::invalid_argument);
  EXPECT_THROW(F::template apply<vector_t>(badsize_tm2, badsize_tm1),
  std::invalid_argument);
}

template <typename F, typename T, typename VT>
void expect_binary_vector_value_error() {
  using std::vector;
  using Eigen::VectorXd;
  typedef Eigen::Matrix<T, VT::RowsAtCompileTime, VT::ColsAtCompileTime> 
  vector_t;

  vector<double> invalid_inputs1 = F::invalid_inputs1();
  if (invalid_inputs1.size() == 0) return;
  vector<double> invalid_inputs2 = F::invalid_inputs2();
  vector<int> int_invalid_inputs1 = F::int_invalid_inputs1();
  vector<int> int_invalid_inputs2 = F::int_invalid_inputs2();
  vector<T> y1(invalid_inputs1.begin(), invalid_inputs1.end());
  vector<T> y2(invalid_inputs2.begin(), invalid_inputs2.end());
  VT a1(invalid_inputs2.size());
  VT a2(invalid_inputs2.size());
  vector_t b1(invalid_inputs1.size());
  vector_t b2(invalid_inputs2.size());
  for (int i = 0; i < a1.size(); ++i) {
      a1(i) = invalid_inputs1[i];
  }
  for (int i = 0; i < a2.size(); ++i) {
      a2(i) = invalid_inputs2[i];
  }
  for (int i = 0; i < b1.size(); ++i) {
      b1(i) = invalid_inputs1[i];
  }
  for (int i = 0; i < b2.size(); ++i) {
      b2(i) = invalid_inputs2[i];
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
  EXPECT_THROW(F::template apply<vector_t>(a1, b2), std::domain_error);
  EXPECT_THROW(F::template apply<vector_t>(b1, a2), std::domain_error);
  EXPECT_THROW(F::template apply<vector_t>(b1, b2), std::domain_error);

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
  vector<VT> d1;
  d1.push_back(a1);
  d1.push_back(a1);
  vector<VT> d2;
  d2.push_back(a2);
  d2.push_back(a2);
  vector<vector_t> e1;
  e1.push_back(b1);
  e1.push_back(b1);
  vector<vector_t> e2;
  e2.push_back(b2);
  e2.push_back(b2);
  EXPECT_THROW(F::template apply<vector<vector_t> >(e1, d2), 
  std::domain_error);
  EXPECT_THROW(F::template apply<vector<vector_t> >(d1, e2), 
  std::domain_error);
  EXPECT_THROW(F::template apply<vector<vector_t> >(e1, e2), 
  std::domain_error);
}

template <typename F, typename T, typename VT>
void expect_binary_vector_error() {
  expect_binary_vector_size_error<F, T, VT>();
  expect_binary_vector_value_error<F, T, VT>();
}
#endif
