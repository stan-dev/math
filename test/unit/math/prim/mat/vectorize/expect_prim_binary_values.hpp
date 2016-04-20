#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_VALUES_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_VALUES_HPP

#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_value.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_std_vector_value.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_matrix_value.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_vector_value.hpp>
#include <Eigen/Dense>

template <typename F>
void expect_prim_binary_values() {
  expect_prim_binary_scalar_value<F>();
  expect_prim_binary_std_vector_value<F>();
  expect_prim_binary_matrix_value<F>();
  expect_prim_binary_vector_value<F, Eigen::VectorXd>();
  expect_prim_binary_vector_value<F, Eigen::RowVectorXd>();
}

#endif
