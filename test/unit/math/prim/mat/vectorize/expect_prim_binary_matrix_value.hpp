#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_matrix_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_std_vector_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_std_vector_matrix_std_vector_matrix_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F>
void expect_prim_binary_matrix_value() {
  using stan::math::is_nan;
  using Eigen::MatrixXd;
  using std::vector;

  vector<int> int_valid_inputs = F::int_valid_inputs();
  vector<double> valid_inputs = F::valid_inputs();
  MatrixXd a(valid_inputs.size(), 3);
  for (int i = 0; i < a.size(); i++) {
    a(i) =  valid_inputs[(i % valid_inputs.size())];
  }

  //Tests int, matrix
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    expect_prim_binary_scalar_matrix_eq<F>(int_valid_inputs[i], a);
  }

  //Tests double, matrix
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    expect_prim_binary_scalar_matrix_eq<F>(valid_inputs[i], a);
  }

  //Tests matrix, matrix
  expect_prim_binary_matrix_matrix_eq<F>(a);

  MatrixXd ab = a.block(1, 1, 1, 1);
  expect_prim_binary_matrix_matrix_eq<F>(ab);

  vector<MatrixXd> d;
  d.push_back(a);
  d.push_back(a);
  //Tests int, vector<matrix>
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    expect_prim_binary_scalar_std_vector_matrix_eq<F>(
    int_valid_inputs[i], d, a.size());
  }

  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    expect_prim_binary_scalar_std_vector_matrix_eq<F>(
    valid_inputs[i], d, a.size());
  }

  expect_prim_binary_std_vector_matrix_std_vector_matrix_eq<F>(d, a.size());
  
}

#endif
