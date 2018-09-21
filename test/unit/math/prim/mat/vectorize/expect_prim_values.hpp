#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_VALUES_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_VALUES_HPP

#include <test/unit/math/prim/mat/vectorize/expect_prim_scalar_value.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_std_vector_value.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_matrix_value.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_vector_value.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_row_vector_value.hpp>

// CALL THIS TO TEST EVERYTHING
// see: apply_scalar_unary_test.cpp for an example
template <typename F>
void expect_prim_values() {
  expect_prim_scalar_value<F>();
  expect_prim_std_vector_value<F>();
  expect_prim_matrix_value<F>();
  expect_prim_vector_value<F>();
  expect_prim_row_vector_value<F>();
}

#endif
