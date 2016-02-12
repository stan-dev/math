#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_REV_EXPECT_ERRORS_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_REV_EXPECT_ERRORS_HPP

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_errors.hpp>

// CALL THIS TO TEST EVERYTHING
// see: apply_scalar_unary_test.cpp for an example
template <typename F>
void expect_errors() {
  using stan::math::var;

  expect_scalar_error<F, var>();
  expect_std_vectors_error<F, var>();
  expect_matrix_error<F, var>();
  expect_vector_error<F, var>();
  expect_row_vector_error<F, var>();
}

#endif
