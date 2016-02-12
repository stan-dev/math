#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_PRIM_EXPECT_ERRORS_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_PRIM_EXPECT_ERRORS_HPP

#include <stdexcept>
#include <test/unit/math/prim/mat/vectorize/expect_errors.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

// CALL THIS TO TEST EVERYTHING
// see: apply_scalar_unary_test.cpp for an example
template <typename F>
void expect_errors() {
  expect_int_scalar_error<F>();
  expect_int_std_vectors_error<F>();

  expect_scalar_error<F, double>();
  expect_std_vectors_error<F, double>();
  expect_matrix_error<F, double>();
  expect_vector_error<F, double>();
  expect_row_vector_error<F, double>();
}

#endif
