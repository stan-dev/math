#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_ERRORS_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_ERRORS_HPP

#include <test/unit/math/prim/mat/vectorize/expect_binary_int_scalar_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_scalar_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_int_std_vector_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_std_vector_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_matrix_error.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_binary_vector_error.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>

/**
 * Tests that the function defined statically in the template class
 * throws errors when supplied with illegal values.  The illegal
 * values are themselves defined statically in the template class.
 *
 * Tests integer, scalar, standard vector, Eigen matrix, Eigen vector,
 * and Eigen row vecor cases, and nested standard vectors of Eigen
 * types. 
 *
 * @tparam F Test class used to define test case.
 */
template <typename F>
void expect_prim_binary_errors() {
  expect_binary_int_scalar_error<F>();
  expect_binary_int_std_vector_error<F>();
  expect_binary_scalar_error<F, double>();
  expect_binary_std_vector_error<F, double>();
  expect_binary_matrix_error<F, double>();
  expect_vector_error<F, double, Eigen::VectorXd>();
  expect_vector_error<F, double, Eigen::RowVectorXd>();
}

#endif
