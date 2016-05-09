#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_VALUES_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_VALUES_HPP

#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_scalar_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_std_vector_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_matrix_value.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_vector_value.hpp>
#include <Eigen/Dense>

//Also will test derivatives
template <typename F>
void expect_rev_binary_values() {
  expect_rev_binary_scalar_value<F>();
  expect_rev_binary_std_vector_value<F>();
  expect_rev_binary_matrix_value<F>();
  expect_rev_binary_vector_value<F, Eigen::VectorXd>();
  expect_rev_binary_vector_value<F, Eigen::RowVectorXd>();
}

#endif
