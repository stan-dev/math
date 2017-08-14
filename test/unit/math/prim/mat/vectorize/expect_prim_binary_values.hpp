#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_VALUES_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_VALUES_HPP

#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_value.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_std_vector_value.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_matrix_value.hpp>
#include <Eigen/Dense>

template <typename F>
void expect_prim_binary_values() {
  expect_prim_binary_scalar_value<F>();
  expect_prim_binary_std_vector_value<F>();
  Eigen::MatrixXd model_matrix;
  Eigen::VectorXd model_vector;
  Eigen::RowVectorXd model_row_vector;
  expect_prim_binary_matrix_value<F>(model_matrix);
  expect_prim_binary_matrix_value<F>(model_vector);
  expect_prim_binary_matrix_value<F>(model_row_vector);
}

#endif
