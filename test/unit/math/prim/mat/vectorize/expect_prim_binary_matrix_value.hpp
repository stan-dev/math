#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_matrix_scalar_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_matrix_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_std_vector_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_std_vector_matrix_scalar_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_std_vector_matrix_std_vector_matrix_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F>
void expect_prim_binary_matrix_value() {
  using Eigen::MatrixXd;
  using std::vector;

  vector<int> int_valid_inputs1 = F::int_valid_inputs1();
  vector<int> int_valid_inputs2 = F::int_valid_inputs2();
  vector<double> valid_inputs1 = F::valid_inputs1();
  vector<double> valid_inputs2 = F::valid_inputs2();

  //Tests int, matrix
  for (size_t i = 0; i < int_valid_inputs1.size(); ++i) {
    MatrixXd a1 = MatrixXd::Constant(5, 3, valid_inputs1[i]);
    MatrixXd a2 = MatrixXd::Constant(5, 3, valid_inputs2[i]);
    expect_prim_binary_scalar_matrix_eq<F>(int_valid_inputs1[i], a2);
    expect_prim_binary_matrix_scalar_eq<F>(a1, int_valid_inputs2[i]);
  }

  //Tests double, matrix
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    MatrixXd a1 = MatrixXd::Constant(5, 3, valid_inputs1[i]);
    MatrixXd a2 = MatrixXd::Constant(5, 3, valid_inputs2[i]);
    expect_prim_binary_scalar_matrix_eq<F>(valid_inputs1[i], a2);
    expect_prim_binary_matrix_scalar_eq<F>(a1, valid_inputs2[i]);
  }

  //Tests matrix, matrix
  MatrixXd a1(valid_inputs1.size(), 3);
  for (int i = 0; i < a1.size(); i++) {
    a1(i) =  valid_inputs1[(i % valid_inputs1.size())];
  }
  MatrixXd a2(valid_inputs2.size(), 3);
  for (int i = 0; i < a2.size(); i++) {
    a2(i) =  valid_inputs1[(i % valid_inputs2.size())];
  }

  expect_prim_binary_matrix_matrix_eq<F>(a1, a2);

  MatrixXd ab1 = a1.block(1, 1, 1, 1);
  MatrixXd ab2 = a2.block(1, 1, 1, 1);
  expect_prim_binary_matrix_matrix_eq<F>(ab1, ab2);

  //Tests int, vector<matrix>
  for (size_t i = 0; i < int_valid_inputs1.size(); ++i) {
    MatrixXd b1 = MatrixXd::Constant(5, 3, valid_inputs1[i]);
    MatrixXd b2 = MatrixXd::Constant(5, 3, valid_inputs2[i]);
    vector<MatrixXd> d1;
    d1.push_back(b1);
    d1.push_back(b1);
    vector<MatrixXd> d2;
    d2.push_back(b2);
    d2.push_back(b2);
    expect_prim_binary_scalar_std_vector_matrix_eq<F>(
      int_valid_inputs1[i], d2);
    expect_prim_binary_std_vector_matrix_scalar_eq<F>(
      d1, int_valid_inputs2[i]);
  }

  //Tests double, vector<matrix>
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    MatrixXd b1 = MatrixXd::Constant(5, 3, valid_inputs1[i]);
    MatrixXd b2 = MatrixXd::Constant(5, 3, valid_inputs2[i]);
    vector<MatrixXd> d1;
    d1.push_back(b1);
    d1.push_back(b1);
    vector<MatrixXd> d2;
    d2.push_back(b2);
    d2.push_back(b2);
    expect_prim_binary_scalar_std_vector_matrix_eq<F>(valid_inputs1[i], d2);
    expect_prim_binary_std_vector_matrix_scalar_eq<F>(d1, valid_inputs2[i]);
  }

  vector<MatrixXd> d1;
  d1.push_back(a1);
  d1.push_back(a1);
  vector<MatrixXd> d2;
  d2.push_back(a2);
  d2.push_back(a2);
  expect_prim_binary_std_vector_matrix_std_vector_matrix_eq<F>(d1, d2);
}

#endif
